using System.Globalization;
using System.Text.Json;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 / GENOMIC-TANDEM-001 (audit A4-4, B09 F156): the option boundaries of the minimap2 chaining, mgaps clustering and etandem ports against the installed programs
/// (<c>option_boundary_oracle.json</c>, scripts/alignment/generate_option_boundary_oracle.py): every option value that the program accepts is accepted and gives the program's result (negative or zero bandwidth, gap limit,
/// counts, scores, separations; NaN / infinite / huge scales and fractions through the C float-to-integer conversion), the values it rejects or dies on are refused with a clear exception.
/// </summary>
[TestFixture]
public class OptionBoundaries_Oracle_Tests
{
    private static readonly Lazy<JsonDocument> Oracle = new(() =>
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.option_boundary_oracle.json";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name)
            ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        return JsonDocument.Parse(stream);
    });

    private static float F(string v) => float.Parse(v, NumberStyles.Float, CultureInfo.InvariantCulture);
    private static double D(string v) => double.Parse(v, NumberStyles.Float, CultureInfo.InvariantCulture);
    private static int I(string v) => int.Parse(v, CultureInfo.InvariantCulture);

    private static MinimapChainOptions ChainOptions(string property, string value) => property switch
    {
        "Bandwidth" => new MinimapChainOptions { Bandwidth = I(value) },
        "MaxTargetGap" => new MinimapChainOptions { MaxTargetGap = I(value), MaxQueryGap = I(value) },
        "MinCount" => new MinimapChainOptions { MinCount = I(value) },
        "MinScore" => new MinimapChainOptions { MinScore = I(value) },
        "ChainGapScale" => new MinimapChainOptions { ChainGapScale = F(value) },
        "ChainSkipScale" => new MinimapChainOptions { ChainSkipScale = F(value) },
        "MaxChainSkip" => new MinimapChainOptions { MaxChainSkip = I(value) },
        "MaxChainIterations" => new MinimapChainOptions { MaxChainIterations = I(value) },
        _ => throw new InvalidOperationException(property),
    };

    private static NucmerClusterOptions ClusterOptions(string property, string value)
    {
        var baseline = new NucmerClusterOptions { MinClusterLength = 1, MaxSeparation = 90, DiagonalDifference = 5, DiagonalFraction = 0.12 };
        return property switch
        {
            "MaxSeparation" => baseline with { MaxSeparation = long.Parse(value, CultureInfo.InvariantCulture) },
            "MinClusterLength" => baseline with { MinClusterLength = long.Parse(value, CultureInfo.InvariantCulture) },
            "DiagonalDifference" => baseline with { DiagonalDifference = I(value) },
            "DiagonalFraction" => baseline with { DiagonalFraction = D(value) },
            "UseExtent" => baseline with { UseExtent = bool.Parse(value) },
            _ => throw new InvalidOperationException(property),
        };
    }

    [Test]
    public void MinimapChainOptions_EveryValueTheProgramAccepts_GivesTheProgramsChains()
    {
        var root = Oracle.Value.RootElement.GetProperty("minimap2");
        string reference = root.GetProperty("reference").GetString()!, query = root.GetProperty("query").GetString()!;
        var failures = new List<string>();
        int rows = 0, empty = 0, nonFinite = 0;
        foreach (var row in root.GetProperty("rows").EnumerateArray())
        {
            string property = row.GetProperty("property").GetString()!, value = row.GetProperty("value").GetString()!;
            Assert.That(row.GetProperty("rc").GetInt32(), Is.EqualTo(0), $"{property} {value}: the program runs");
            var expected = row.GetProperty("chains").EnumerateArray().Select(c => (c[0].GetInt32(), c[1].GetInt32())).ToList();
            var actual = Minimap2Chaining.MapChains([reference], query, 15, 10, ChainOptions(property, value))
                .Select(c => (c.Score, c.Anchors.Count)).OrderByDescending(c => c.Score).ThenByDescending(c => c.Item2).ToList();
            if (!expected.SequenceEqual(actual))
                failures.Add($"{property} = {value}: expected [{string.Join(", ", expected)}] actual [{string.Join(", ", actual)}]");
            rows++;
            if (expected.Count == 0) empty++;
            if (value is "NaN" or "Infinity" or "-Infinity" or "1e30") nonFinite++;
        }

        Assert.That(failures, Is.Empty, string.Join("\n", failures));
        Assert.That(rows, Is.EqualTo(30));
        Assert.That(empty, Is.EqualTo(1), "a bandwidth of -5 leaves no chain");
        Assert.That(nonFinite, Is.GreaterThanOrEqualTo(6), "NaN / infinite / huge scales go through the C float-to-int conversion");
    }

    [Test]
    public void MinimapChainOptions_KmerSizeOutside1To28_IsRefused_LikeMmCheckOptAndMmSketch()
    {
        var anchors = new[] { new MinimapAnchor(0, false, 100, 100, 15) };
        foreach (int k in new[] { 0, -3, 29, 40, int.MinValue, int.MaxValue })
        {
            var ex = Assert.Throws<ArgumentOutOfRangeException>(() => Minimap2Chaining.Chain(anchors, new MinimapChainOptions { KmerSize = k }), $"k = {k}")!;
            Assert.That(ex.Message, Does.Contain("-k and -w must be positive").And.Contain("k <= 28"));
            Assert.That(ex.ParamName, Is.EqualTo("options"));
        }

        Assert.DoesNotThrow(() => Minimap2Chaining.Chain(anchors, new MinimapChainOptions { KmerSize = 1 }));
        Assert.DoesNotThrow(() => Minimap2Chaining.Chain(anchors, new MinimapChainOptions { KmerSize = 28 }));
    }

    [Test]
    public void CollectAnchors_QueryOccurrenceFraction_NaNNegativeInfinite_SwitchTheFilterOff_LikeTheProgram()
    {
        // minimap2 2.26: --q-occ-frac nan / -1 / inf give the chains of the default (the comparisons q_occ_frac > 0 and count > n x q_occ_frac are false)
        var root = Oracle.Value.RootElement.GetProperty("minimap2");
        string reference = root.GetProperty("reference").GetString()!, query = root.GetProperty("query").GetString()!;
        var dflt = Minimap2Chaining.CollectAnchors([reference], query);
        foreach (float f in new[] { float.NaN, -1f, 0f, float.PositiveInfinity })
            Assert.That(Minimap2Chaining.CollectAnchors([reference], query, queryOccurrenceFraction: f), Is.EqualTo(dflt), $"fraction {f}");
    }

    [Test]
    public void NucmerClusterOptions_EveryValueMgapsAccepts_GivesMgapsClusters()
    {
        var root = Oracle.Value.RootElement.GetProperty("mgaps");
        var matches = root.GetProperty("matches").EnumerateArray().Select(m => (m[0].GetInt64(), m[1].GetInt64(), m[2].GetInt64())).ToList();
        var failures = new List<string>();
        int rows = 0, nonFinite = 0;
        foreach (var row in root.GetProperty("rows").EnumerateArray())
        {
            string property = row.GetProperty("property").GetString()!, value = row.GetProperty("value").GetString()!;
            Assert.That(row.GetProperty("rc").GetInt32(), Is.EqualTo(0), $"{property} {value}: mgaps runs");
            string expected = string.Join(" # ", row.GetProperty("clusters").EnumerateArray().Select(c => string.Join(" ", c.EnumerateArray().Select(m => $"{m[0]}:{m[1]}:{m[2]}"))));
            string actual = string.Join(" # ", NucmerClusters.Cluster(matches, ClusterOptions(property, value)).Select(c => string.Join(" ", c.Matches.Select(m => $"{m.Start1}:{m.Start2}:{m.Length}"))));
            if (expected != actual)
                failures.Add($"{property} = {value}: expected '{expected}' actual '{actual}'");
            rows++;
            if (value is "NaN" or "Infinity" or "1e30") nonFinite++;
        }

        Assert.That(failures, Is.Empty, string.Join("\n", failures));
        Assert.That(rows, Is.EqualTo(17));
        Assert.That(nonFinite, Is.EqualTo(3), "NaN / infinite / beyond-64-bit products become LONG_MIN (no fraction added), exactly as mgaps computes them");
    }
}
