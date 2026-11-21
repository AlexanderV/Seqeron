using System.Globalization;
using System.IO.Compression;
using System.Text.Json;
using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// A5-9 (B09 F173): the <c>-x</c> presets of minimap2 2.26 (<c>options.c mm_set_opt</c>: <see cref="MinimapPresets"/>), <c>-F</c> / <c>-G</c> and the PAF text (<see cref="Minimap2Chaining.FormatPaf"/>) against the installed program:
/// 220 random cases per preset + 40 cases with options after the preset (A6-7) (TestData/Alignment/minimap2_preset_oracle.jsonl.gz, scripts/alignment/generate_minimap2_preset_oracle.py), the WHOLE stdout of <c>minimap2 -x PRESET [-F|-G] ref.fa query.fa</c> (columns 1 - 12 and the
/// tags tp cm s1 s2 dv rl of every line) compared as text. Primary and secondary hits, both strands, all-vs-all self hits, short reads and spliced queries are in the sets.
/// </summary>
[TestFixture]
public class Minimap2Presets_Tests
{
    private sealed record Case(string Preset, string[] Extra, string[] Names, string[] Targets, string QName, string Query, string Paf);

    private static readonly Lazy<List<Case>> Oracle = new(() =>
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.minimap2_preset_oracle.jsonl.gz";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name)
            ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var gz = new GZipStream(stream, CompressionMode.Decompress);
        using var reader = new StreamReader(gz);
        var cases = new List<Case>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            using var doc = JsonDocument.Parse(line);
            var r = doc.RootElement;
            string[] Strings(string p) => [.. r.GetProperty(p).EnumerateArray().Select(e => e.GetString()!)];
            cases.Add(new Case(r.GetProperty("preset").GetString()!, Strings("extra"), Strings("names"), Strings("targets"), r.GetProperty("qname").GetString()!, r.GetProperty("query").GetString()!, r.GetProperty("paf").GetString()!));
        }

        return cases;
    });

    /// <summary>The options after <c>-x PRESET</c> on the command line of the case: <see cref="MinimapPresets.Apply(string?, MinimapOverrides)"/> (the program applies the preset first, every option after it wins).</summary>
    private static MinimapMapOptions Options(Case c)
    {
        var o = new MinimapOverrides();
        string[] x = c.Extra;
        int i = 0;
        while (i < x.Length)
        {
            string option = x[i++];
            string Value() => x[i++];
            int Int() => int.Parse(Value(), CultureInfo.InvariantCulture);
            double Dbl() => double.Parse(Value(), CultureInfo.InvariantCulture);
            switch (option)
            {
                case "-k": o = o with { KmerSize = Int() }; break;
                case "-w": o = o with { WindowSize = Int() }; break;
                case "-g": o = o with { MaxGap = Int() }; break;
                case "-n": o = o with { MinCount = Int() }; break;
                case "-m": o = o with { MinChainScore = Int() }; break;
                case "-N": o = o with { BestN = Int() }; break;
                case "-e": o = o with { OccDist = Int() }; break;
                case "-F": o = o with { MaxFragmentLength = Int() }; break;
                case "-G": o = o with { MaxIntronLength = Int() }; break;
                case "--max-chain-skip": o = o with { MaxChainSkip = Int() }; break;
                case "-p": o = o with { PrimaryRatio = Dbl() }; break;
                case "-M": o = o with { MaskLevel = Dbl() }; break;
                case "-f": o = o with { MidOccFraction = Dbl() }; break;
                case "-r":
                    var b = Value().Split(',');
                    o = o with { Bandwidth = int.Parse(b[0], CultureInfo.InvariantCulture), BandwidthLong = b.Length > 1 ? int.Parse(b[1], CultureInfo.InvariantCulture) : o.BandwidthLong };
                    break;
                case "-U":
                    var u = Value().Split(',');
                    o = o with { MinMidOcc = int.Parse(u[0], CultureInfo.InvariantCulture), MaxMidOcc = u.Length > 1 ? int.Parse(u[1], CultureInfo.InvariantCulture) : o.MaxMidOcc };
                    break;
                case "--rmq=yes": o = o with { UseRmq = true }; break;
                case "--rmq=no": o = o with { UseRmq = false }; break;
                case "--secondary=yes": o = o with { OutputSecondary = true }; break;
                case "--secondary=no": o = o with { OutputSecondary = false }; break;
                case "-H": o = o with { HomopolymerCompressed = true }; break;
                case "--no-long-join": o = o with { NoLongJoin = true }; break;
                default: throw new InvalidOperationException($"option {option} of the oracle is not handled by the test");
            }
        }

        return MinimapPresets.Apply(c.Preset, o);
    }

    private static bool HasOverride(Case c) => c.Extra.Length > 0 && c.Extra[0] is not ("-F" or "-G");          // the cases of A5-9 carry one -F / -G, the cases of A6-7 the options of the generator's OVERRIDES table

    /// <summary>The rows of a preset: all cases, the cases with <c>-F</c> / <c>-G</c> (after the preset: sr honours -F, the splice presets -G, the others ignore -G), the cases with the option overrides of A6-7, the PAF lines.</summary>
    private static void Check(string preset, int cases, int withFragmentOrIntron, int withOverrides, int pafLines)
    {
        var all = Oracle.Value.Where(c => c.Preset == preset).ToList();
        Assert.That(all.Count, Is.EqualTo(cases), preset);
        var failures = new List<string>();
        int lines = 0;
        foreach (var c in all)
        {
            var result = Minimap2Chaining.Map(c.Targets, c.Query, Options(c), c.Names, c.QName);
            string got = Minimap2Chaining.FormatPaf(result, c.QName, c.Query.Length, c.Names, [.. c.Targets.Select(t => t.Length)]);
            lines += c.Paf.Count(ch => ch == '\n');
            if (got != c.Paf)
                failures.Add($"{c.QName} {string.Join(' ', c.Extra)} ({c.Targets.Length} targets, {c.Query.Length} bp)\n   expected {c.Paf.Replace("\n", " || ")}\n   actual   {got.Replace("\n", " || ")}");
        }

        Assert.That(failures, Is.Empty, $"{failures.Count} of {all.Count}:\n" + string.Join("\n", failures.Take(3)));
        Assert.That(lines, Is.EqualTo(pafLines), "PAF lines in the set");
        Assert.That(all.Count(c => c.Extra.Contains("-F") || c.Extra.Contains("-G")), Is.EqualTo(withFragmentOrIntron), "cases with -F / -G");
        Assert.That(all.Count(HasOverride), Is.EqualTo(withOverrides), "cases with option overrides after the preset");
    }

    [Test] public void MapOnt_EqualsMinimap2_OnAllCases() => Check("map-ont", 260, 38, 40, 260);
    [Test] public void MapPb_HomopolymerCompressed_EqualsMinimap2_OnAllCases() => Check("map-pb", 260, 35, 40, 240);
    [Test] public void MapHifi_EqualsMinimap2_OnAllCases() => Check("map-hifi", 260, 54, 40, 261);
    [Test] public void AvaOnt_AllChainsNoSelfNoDual_EqualsMinimap2_OnAllCases() => Check("ava-ont", 260, 54, 40, 507);
    [Test] public void AvaPb_AllChainsNoSelfNoDual_EqualsMinimap2_OnAllCases() => Check("ava-pb", 260, 42, 40, 381);
    [Test] public void Asm5_Rmq_EqualsMinimap2_OnAllCases() => Check("asm5", 260, 49, 40, 263);
    [Test] public void Asm10_Rmq_EqualsMinimap2_OnAllCases() => Check("asm10", 260, 47, 40, 260);
    [Test] public void Asm20_Rmq_EqualsMinimap2_OnAllCases() => Check("asm20", 260, 54, 40, 257);
    [Test] public void Sr_ShortReads_WithAndWithoutFragmentLength_EqualsMinimap2_OnAllCases() => Check("sr", 260, 63, 40, 258);
    [Test] public void Splice_CdnaChaining_WithAndWithoutMaxIntronLength_EqualsMinimap2_OnAllCases() => Check("splice", 260, 49, 40, 292);
    [Test] public void SpliceHq_EqualsMinimap2_OnAllCases() => Check("splice:hq", 260, 45, 40, 270);

    [Test]
    public void Presets_AreTheOptionsC_Tables()
    {
        // options.c mm_set_opt (v2.26)
        var ont = MinimapPresets.Apply("map-ont");
        Assert.That((ont.Seed.KmerSize, ont.Seed.WindowSize, ont.Bandwidth, ont.MaxGap, ont.MinChainScore), Is.EqualTo((15, 10, 500, 5000, 40)));
        var hifi = MinimapPresets.Apply("map-hifi");
        Assert.That((hifi.Seed.KmerSize, hifi.Seed.WindowSize, hifi.MaxGap, hifi.Seed.MinMidOcc, hifi.Seed.MaxMidOcc, hifi.Scoring.Match, hifi.Scoring.Mismatch, hifi.Scoring.GapOpen, hifi.Scoring.GapOpen2), Is.EqualTo((19, 19, 10000, 50, 500, 1, 4, 6, 26)));
        Assert.That(MinimapPresets.Apply("map10k"), Is.EqualTo(MinimapPresets.Apply("map-pb")));
        Assert.That(MinimapPresets.Apply("map-ccs"), Is.EqualTo(MinimapPresets.Apply("map-hifi")));
        var pb = MinimapPresets.Apply("map-pb");
        Assert.That((pb.Seed.HomopolymerCompressed, pb.Seed.KmerSize, pb.Seed.WindowSize), Is.EqualTo((true, 19, 10)));
        var ava = MinimapPresets.Apply("ava-ont");
        Assert.That((ava.AllChains, ava.NoLongJoin, ava.Seed.SkipSelf, ava.Seed.SkipDual, ava.MinChainScore, ava.PrimaryRatio, ava.Bandwidth, ava.BandwidthLong, ava.Seed.OccDist, ava.Seed.WindowSize), Is.EqualTo((true, true, true, true, 100, 0f, 2000, 2000, 0, 5)));
        var avaPb = MinimapPresets.Apply("ava-pb");
        Assert.That((avaPb.Seed.HomopolymerCompressed, avaPb.Seed.KmerSize, avaPb.BandwidthLong), Is.EqualTo((true, 19, 500)));
        var asm5 = MinimapPresets.Apply("asm5");
        Assert.That((asm5.UseRmq, asm5.BestN, asm5.Bandwidth, asm5.BandwidthLong, asm5.Scoring.Mismatch, asm5.Scoring.GapOpen, asm5.Scoring.GapOpen2, asm5.Scoring.GapExtend, asm5.Scoring.ZDrop), Is.EqualTo((true, 50, 1000, 100000, 19, 39, 81, 3, 200)));
        Assert.That(MinimapPresets.Apply("asm20").Seed.WindowSize, Is.EqualTo(10));
        Assert.That(MinimapPresets.Apply("asm10").Scoring.GapOpen2, Is.EqualTo(41));
        var sr = MinimapPresets.Apply("sr");
        Assert.That((sr.ShortRead, sr.OutputSecondary, sr.Seed.KmerSize, sr.Seed.WindowSize, sr.MaxFragmentLength, sr.MaxGap, sr.MinCount, sr.MinChainScore, sr.PrimaryRatio, sr.BestN, sr.Seed.MidOcc, sr.Seed.MaxOcc, sr.Scoring.EndBonus),
            Is.EqualTo((true, false, 21, 11, 800, 100, 2, 25, 0.5f, 20, 1000, 5000, 10)));
        Assert.That(MinimapPresets.Apply("short"), Is.EqualTo(sr));
        var splice = MinimapPresets.Apply("splice");
        Assert.That((splice.Splice, splice.MaxGap, splice.MaxGapRef, splice.Bandwidth, splice.Scoring.JunctionBonus, splice.Scoring.NonCanonical, splice.Scoring.GapOpen2, splice.Scoring.GapExtend2), Is.EqualTo((true, 2000, 200000, 200000, 9, 9, 32, 0)));
        var hq = MinimapPresets.Apply("splice:hq");
        Assert.That((hq.Scoring.JunctionBonus, hq.Scoring.Mismatch, hq.Scoring.GapOpen, hq.Scoring.GapOpen2), Is.EqualTo((5, 4, 6, 24)));
        Assert.That(MinimapPresets.Apply(null), Is.EqualTo(new MinimapMapOptions()));
        Assert.Throws<ArgumentException>(() => MinimapPresets.Apply("lr:hq"));        // not a preset of 2.26
    }

    [Test]
    public void WithMaxIntronLength_OnlyInSpliceMode()
    {
        var ont = MinimapPresets.WithMaxIntronLength(MinimapPresets.Apply("map-ont"), 7000);
        Assert.That(ont, Is.EqualTo(MinimapPresets.Apply("map-ont")));
        var spliced = MinimapPresets.WithMaxIntronLength(MinimapPresets.Apply("splice"), 7000);
        Assert.That((spliced.MaxGapRef, spliced.Bandwidth, spliced.BandwidthLong), Is.EqualTo((7000, 7000, 7000)));
    }

    [Test]
    public void FormatPaf_PrintsNothingForNoHit_UnlessPafNoHit()
    {
        var result = Minimap2Chaining.Map(["ACGTACGTAGCTAGCTAGCATCGATCGATGCATGCATGCA"], "TTTTTTTTTTTTTTTTTTTT", new MinimapMapOptions());
        Assert.That(Minimap2Chaining.FormatPaf(result, "q", 20, ["t"], [40]), Is.Empty);
        Assert.That(Minimap2Chaining.FormatPaf(result, "q", 20, ["t"], [40], printUnmapped: true), Is.EqualTo($"q\t20\t0\t0\t*\t*\t0\t0\t0\t0\t0\t0\trl:i:{result.RepeatLength}\n"));
    }
}
