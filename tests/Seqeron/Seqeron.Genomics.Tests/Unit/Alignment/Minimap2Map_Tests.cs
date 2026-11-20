using System.Globalization;
using System.Text.Json;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 (audit A4-7, B09 F166): <see cref="Minimap2Chaining.Map"/> = minimap2 2.26 up to the end of the chaining stage: the primary / secondary marking (<c>mm_set_parent</c>, <c>mm_select_sub</c>, <c>-N</c> <c>-p</c>
/// <c>--mask-level</c> <c>--mask-len</c> <c>--hard-mask-level</c>), the long-join re-chaining (<c>-r b,B</c>, <c>mg_lchain_rmq</c>), <c>--rmq</c>, the divergence estimate (<c>dv</c>), the mapping quality and the repeat length <c>rl</c>.
/// Oracle (TestData/Alignment/minimap2_map_oracle.jsonl.gz, scripts/alignment/generate_minimap2_map_oracle.py): the PAF lines (without <c>-c</c>) of the installed minimap2 2.26 in four groups of >= 300 cases: columns 3 - 12 and the tags
/// <c>tp cm s1 s2 dv rl</c> of every line.
/// </summary>
[TestFixture]
public class Minimap2Map_Tests
{
    private sealed record Case(string Group, MinimapMapOptions Options, string[] Names, string QName, string[] Targets, string Query, object?[][] Hits);

    private static readonly Lazy<List<Case>> Oracle = new(ReadOracle);

    private static T? Opt<T>(JsonElement o, string key, Func<JsonElement, T> read) where T : struct => o.GetProperty(key).ValueKind == JsonValueKind.Null ? null : read(o.GetProperty(key));

    private static List<Case> ReadOracle()
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.minimap2_map_oracle.jsonl.gz";
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
            var defaults = new MinimapMapOptions();
            int? bw = null, bwLong = null;
            if (o.GetProperty("r").ValueKind != JsonValueKind.Null)
            {
                bw = o.GetProperty("r")[0].GetInt32();
                bwLong = o.GetProperty("r")[1].GetInt32();
            }

            var options = new MinimapMapOptions
            {
                Seed = new MinimapSeedOptions { KmerSize = o.GetProperty("k").GetInt32(), WindowSize = o.GetProperty("w").GetInt32() },
                MaxGap = Opt(o, "g", e => e.GetInt32()) ?? defaults.MaxGap,
                Bandwidth = bw ?? defaults.Bandwidth,
                BandwidthLong = bwLong ?? defaults.BandwidthLong,
                MinCount = Opt(o, "n", e => e.GetInt32()) ?? defaults.MinCount,
                MinChainScore = Opt(o, "m", e => e.GetInt32()) ?? defaults.MinChainScore,
                // main.c: "-N 0" is changed to the previous value of -N with --secondary=no
                BestN = Opt(o, "N", e => e.GetInt32()) is > 0 and var bestN ? bestN : defaults.BestN,
                OutputSecondary = Opt(o, "N", e => e.GetInt32()) != 0,
                PrimaryRatio = Opt(o, "p", e => (float)e.GetDouble()) ?? defaults.PrimaryRatio,
                MaskLevel = Opt(o, "mask_level", e => (float)e.GetDouble()) ?? defaults.MaskLevel,
                MaskLength = Opt(o, "mask_len", e => e.GetInt32()) ?? defaults.MaskLength,
                HardMaskLevel = o.GetProperty("hard_mask").GetBoolean(),
                UseRmq = o.GetProperty("rmq").GetBoolean(),
                ChainGapScale = Opt(o, "chain_gap_scale", e => (float)e.GetDouble()) ?? defaults.ChainGapScale,
                ChainSkipScale = Opt(o, "chain_skip_scale", e => (float)e.GetDouble()) ?? defaults.ChainSkipScale,
                MaxChainSkip = Opt(o, "max_chain_skip", e => e.GetInt32()) ?? defaults.MaxChainSkip,
                MaxChainIterations = Opt(o, "max_chain_iter", e => e.GetInt32()) ?? defaults.MaxChainIterations,
            };
            var hits = r.GetProperty("hits").EnumerateArray().Select(h => h.EnumerateArray().Select(e => e.ValueKind switch
            {
                JsonValueKind.Null => null,
                JsonValueKind.String => (object)e.GetString()!,
                _ => e.GetInt32(),
            }).ToArray()).ToArray();
            cases.Add(new Case(r.GetProperty("group").GetString()!, options, r.GetProperty("names").EnumerateArray().Select(e => e.GetString()!).ToArray(), r.GetProperty("qname").GetString()!,
                r.GetProperty("targets").EnumerateArray().Select(e => e.GetString()!).ToArray(), r.GetProperty("query").GetString()!, hits));
        }

        return cases;
    }

    // [qs, qe, strand, tname, ts, te, mlen, blen, mapq, tp, cm, s1, s2, dv, rl] as the PAF prints them
    private static object?[][] Paf(Case c, MinimapMapResult result) => result.Hits.Select(h => new object?[]
    {
        h.QueryStart, h.QueryEnd, h.Reverse ? "-" : "+", c.Names[h.TargetId], h.TargetStart, h.TargetEnd, h.MatchLength, h.BlockLength, h.MapQ, h.IsPrimary ? "P" : "S", h.AnchorCount, h.Score,
        h.IsPrimary ? h.SecondaryScore : null,
        h.Divergence is >= 0f and <= 1f ? (h.Divergence == 0f ? "0" : ((double)h.Divergence).ToString("0.0000", CultureInfo.InvariantCulture)) : null,
        result.RepeatLength,
    }).ToArray();

    private static string Show(object?[][] hits) => string.Join(" | ", hits.Select(h => string.Join(",", h.Select(x => x?.ToString() ?? "-"))));

    private static (int Cases, int Hits, int Secondary, int Multi, int Reverse) CheckGroup(string group, List<string> failures)
    {
        int cases = 0, hits = 0, secondary = 0, multi = 0, reverse = 0;
        foreach (var c in Oracle.Value.Where(c => c.Group == group))
        {
            cases++;
            MinimapMapResult actual;
            try
            {
                actual = Minimap2Chaining.Map(c.Targets, c.Query, c.Options, c.Names, c.QName);
            }
            catch (Exception ex) when (ex is NotSupportedException)
            {
                failures.Add($"{group}: {ex.Message}");
                continue;
            }

            var mine = Paf(c, actual);
            bool equal = mine.Length == c.Hits.Length && mine.Zip(c.Hits).All(p => p.First.SequenceEqual(p.Second));
            if (!equal)
                failures.Add($"{group} k={c.Options.Seed.KmerSize} w={c.Options.Seed.WindowSize} N={c.Options.BestN} p={c.Options.PrimaryRatio} r={c.Options.Bandwidth},{c.Options.BandwidthLong} rmq={c.Options.UseRmq}\n expected {Show(c.Hits)}\n actual   {Show(mine)}");
            hits += c.Hits.Length;
            secondary += c.Hits.Count(h => (string)h[9]! == "S");
            if (c.Hits.Length > 1) multi++;
            if (c.Hits.Any(h => (string)h[2]! == "-")) reverse++;
        }

        return (cases, hits, secondary, multi, reverse);
    }

    [Test]
    public void Map_SecondaryChainMarking_EqualsMinimap2()
    {
        var failures = new List<string>();
        var (cases, _, secondary, multi, reverse) = CheckGroup("sec", failures);
        Assert.That(failures, Is.Empty, $"{failures.Count} of {cases}:\n" + string.Join("\n", failures.Take(3)));
        Assert.That(cases, Is.GreaterThanOrEqualTo(300));
        Assert.That(secondary, Is.GreaterThan(150));
        Assert.That(multi, Is.GreaterThan(200));
        Assert.That(reverse, Is.GreaterThan(60));
    }

    [Test]
    public void Map_ChainingOptions_EqualMinimap2()
    {
        var failures = new List<string>();
        var (cases, hits, _, multi, reverse) = CheckGroup("misc", failures);
        Assert.That(failures, Is.Empty, $"{failures.Count} of {cases}:\n" + string.Join("\n", failures.Take(3)));
        Assert.That(cases, Is.GreaterThanOrEqualTo(300));
        Assert.That(hits, Is.GreaterThan(400));
        Assert.That(multi, Is.GreaterThan(100));
        Assert.That(reverse, Is.GreaterThan(30));
    }

    [Test]
    public void Map_LongJoinRechaining_EqualsMinimap2()
    {
        var failures = new List<string>();
        var (cases, hits, _, multi, _) = CheckGroup("lj", failures);
        Assert.That(failures, Is.Empty, $"{failures.Count} of {cases}:\n" + string.Join("\n", failures.Take(3)));
        Assert.That(cases, Is.GreaterThanOrEqualTo(300));
        Assert.That(hits, Is.GreaterThan(300));
        Assert.That(multi, Is.GreaterThan(100));
    }

    [Test]
    public void Map_RmqChaining_EqualsMinimap2()
    {
        var failures = new List<string>();
        var (cases, hits, _, multi, _) = CheckGroup("rmq", failures);
        Assert.That(failures, Is.Empty, $"{failures.Count} of {cases}:\n" + string.Join("\n", failures.Take(3)));
        Assert.That(cases, Is.GreaterThanOrEqualTo(300));
        Assert.That(hits, Is.GreaterThan(300));
        Assert.That(multi, Is.GreaterThan(80));
    }
}
