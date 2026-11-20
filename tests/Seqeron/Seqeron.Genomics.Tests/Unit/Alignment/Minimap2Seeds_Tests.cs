using System.Text.Json;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 (audit A4-7, B09 F165): the seed stage of minimap2 2.26 beyond the defaults, <see cref="Minimap2Chaining.CollectSeeds"/>: the frequent-minimizer rules (<c>mid_occ</c>, <c>mm_seed_select</c>,
/// <c>mm_seed_mz_flt</c>, the repeat length <c>rl</c>), <c>-f</c> / <c>-U</c>, <c>-H</c>, the hidden <c>-T</c> (sdust), <c>-X</c> and <c>--for-only</c> / <c>--rev-only</c>.
/// Oracle (TestData/Alignment/minimap2_seed_oracle.jsonl.gz, scripts/alignment/generate_minimap2_seed_oracle.py): per case <c>minimap2 [options] --print-seeds</c>: the anchor array, <c>RS</c> (repeat length) and the
/// <c>mid_occ = N</c> line of the log, in six groups of >= 300 cases.
/// </summary>
[TestFixture]
public class Minimap2Seeds_Tests
{
    private sealed record Case(string Group, MinimapSeedOptions Options, string[] Names, string QName, string[] Targets, string Query, int MidOcc, int RepLen, MinimapAnchor[] Seeds);

    private static readonly Lazy<List<Case>> Oracle = new(ReadOracle);

    private static List<Case> ReadOracle()
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.minimap2_seed_oracle.jsonl.gz";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name)
            ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var gzip = new System.IO.Compression.GZipStream(stream, System.IO.Compression.CompressionMode.Decompress);
        using var reader = new StreamReader(gzip);
        var cases = new List<Case>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length == 0) continue;
            using var doc = JsonDocument.Parse(line);
            var r = doc.RootElement;
            var o = r.GetProperty("opts");
            var options = new MinimapSeedOptions
            {
                KmerSize = o.GetProperty("k").GetInt32(),
                WindowSize = o.GetProperty("w").GetInt32(),
                HomopolymerCompressed = o.GetProperty("H").GetBoolean(),
                SdustThreshold = o.GetProperty("T").GetInt32(),
                // -f x[,y]: below 1 the fraction of repeated minimizers to filter, otherwise the occurrence limit itself (main.c); the part after the comma is max_occ
                MidOccFraction = o.GetProperty("f").ValueKind == JsonValueKind.Null || o.GetProperty("f").GetDouble() >= 1.0 ? 2e-4f : (float)o.GetProperty("f").GetDouble(),
                MidOcc = o.GetProperty("f").ValueKind != JsonValueKind.Null && o.GetProperty("f").GetDouble() >= 1.0 ? (int)(o.GetProperty("f").GetDouble() + .499) : 0,
                MaxOcc = o.GetProperty("fmax").ValueKind == JsonValueKind.Null ? 0 : o.GetProperty("fmax").GetInt32(),
                OccDist = o.GetProperty("e").ValueKind == JsonValueKind.Null ? 500 : o.GetProperty("e").GetInt32(),
                MinMidOcc = o.GetProperty("U").ValueKind == JsonValueKind.Null ? 10 : o.GetProperty("U")[0].GetInt32(),
                MaxMidOcc = o.GetProperty("U").ValueKind == JsonValueKind.Null ? 1_000_000 : o.GetProperty("U")[1].GetInt32(),
                QueryOccurrenceFraction = o.GetProperty("q_occ_frac").ValueKind == JsonValueKind.Null ? 0.01f : (float)o.GetProperty("q_occ_frac").GetDouble(),
                SkipSelf = o.GetProperty("X").GetBoolean(),
                SkipDual = o.GetProperty("X").GetBoolean(),
                ForwardOnly = o.GetProperty("for_only").GetBoolean(),
                ReverseOnly = o.GetProperty("rev_only").GetBoolean(),
            };
            var seeds = r.GetProperty("seeds").EnumerateArray()
                .Select(s => new MinimapAnchor(s[0].GetInt32(), s[1].GetInt32() == 1, s[2].GetInt32(), s[3].GetInt32(), s[4].GetInt32())).ToArray();
            cases.Add(new Case(r.GetProperty("group").GetString()!, options, r.GetProperty("names").EnumerateArray().Select(e => e.GetString()!).ToArray(), r.GetProperty("qname").GetString()!,
                r.GetProperty("targets").EnumerateArray().Select(e => e.GetString()!).ToArray(), r.GetProperty("query").GetString()!, r.GetProperty("mid_occ").GetInt32(), r.GetProperty("rep_len").GetInt32(), seeds));
        }

        return cases;
    }

    [Test]
    public void CollectSeeds_EqualsMinimap2_AllGroups()
    {
        var cases = Oracle.Value;
        var perGroup = new Dictionary<string, int>();
        var failures = new List<string>();
        int repeats = 0, defaultMid = 0, customMid = 0, selected = 0, dropped = 0;
        foreach (var c in cases)
        {
            perGroup[c.Group] = perGroup.GetValueOrDefault(c.Group) + 1;
            var actual = Minimap2Chaining.CollectSeeds(c.Targets, c.Query, c.Options, c.Names, c.QName);
            if (actual.MidOcc != c.MidOcc || actual.RepeatLength != c.RepLen || !actual.Anchors.SequenceEqual(c.Seeds))
                failures.Add($"{c.Group} k={c.Options.KmerSize} w={c.Options.WindowSize} H={c.Options.HomopolymerCompressed} T={c.Options.SdustThreshold} f={c.Options.MidOccFraction} U={c.Options.MinMidOcc},{c.Options.MaxMidOcc} X={c.Options.SkipSelf}: " +
                    $"mid_occ {actual.MidOcc} vs {c.MidOcc}, rl {actual.RepeatLength} vs {c.RepLen}, seeds {actual.Anchors.Count} vs {c.Seeds.Length} (first difference {FirstDifference(actual.Anchors, c.Seeds)})");
            if (c.RepLen > 0) repeats++;
            if (c.MidOcc == 10) defaultMid++; else customMid++;
            if (c.Group == "freq" && c.RepLen > 0 && c.Seeds.Length > 0) selected++;
            if (c.Group == "ava" && c.Seeds.Length > 0) dropped++;
        }

        Assert.That(failures, Is.Empty, $"{failures.Count} of {cases.Count}:\n" + string.Join("\n", failures.Take(4)));
        foreach (var g in new[] { "freq", "fu", "hpc", "dust", "ava", "strand" })
            Assert.That(perGroup[g], Is.GreaterThanOrEqualTo(300), g);
        Assert.That(repeats, Is.GreaterThan(150));          // cases with frequent minimizers (rl > 0)
        Assert.That(selected, Is.GreaterThan(80));          // ... that still have seeds (mm_seed_select kept some / other minimizers)
        Assert.That(customMid, Is.GreaterThan(200));        // mid_occ above its floor of 10
        Assert.That(defaultMid, Is.GreaterThan(200));
        Assert.That(dropped, Is.GreaterThan(100));
    }

    private static string FirstDifference(IReadOnlyList<MinimapAnchor> actual, MinimapAnchor[] expected)
    {
        for (int i = 0; i < Math.Min(actual.Count, expected.Length); i++)
            if (actual[i] != expected[i])
                return $"#{i} {actual[i]} vs {expected[i]}";
        return "none (length)";
    }

    [Test]
    public void FrequentMinimizers_AreHandledLikeTheProgram_NotRefused()
    {
        // a 60 kb reference with 3 000 poly-A: the poly-A minimizer is the one frequent key (mid_occ = 10 + ... ); the old port refused this query, the program drops the minimizer except for one per 500 query bases of the streak
        var rng = new Random(5);
        string reference = new string(Enumerable.Range(0, 60_000).Select(_ => "ACGT"[rng.Next(4)]).ToArray()) + new string('A', 3000);
        var seeds = Minimap2Chaining.CollectSeeds([reference], new string('A', 24) + reference[..200]);
        Assert.That(seeds.MidOcc, Is.GreaterThanOrEqualTo(10));
        Assert.That(seeds.RepeatLength, Is.GreaterThan(0));
        Assert.That(Minimap2Chaining.CollectAnchors([reference], new string('A', 24) + reference[..200]), Is.EqualTo(seeds.Anchors));
        Assert.Throws<ArgumentOutOfRangeException>(() => Minimap2Chaining.CollectSeeds(["ACGT"], "ACGT", new MinimapSeedOptions { KmerSize = 29 }));
        Assert.Throws<ArgumentOutOfRangeException>(() => Minimap2Chaining.CollectSeeds(["ACGT"], "ACGT", new MinimapSeedOptions { WindowSize = 256 }));
        Assert.That(Minimap2Chaining.CollectSeeds(["ACGT"], "ACGT").Anchors, Is.Empty);
    }
}
