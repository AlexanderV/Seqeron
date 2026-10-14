// CODON-RARE-001 — MCP wrappers codon_minmax_profile / find_rare_codon_clusters (FIN-B02 WP5, A1-19).
// The tools delegate to CodonOptimizer.CalculateMinMaxProfile / FindRareCodonClusters; the values
// below are the reference-locked values of
// tests/Seqeron/Seqeron.Genomics.Tests/Unit/MolTools/CodonOptimizer_RareCodonClusters_Tests.cs
// (Clarke & Clark 2008 %MinMax; Chartier et al. 2012 Sherlocc rule), E. coli K12 preset.
using NUnit.Framework;
using Seqeron.Genomics.MolTools;
using Seqeron.Mcp.MolTools.Models;
using Seqeron.Mcp.MolTools.Tools;

namespace Seqeron.Mcp.MolTools.Tests;

[TestFixture]
public class RareCodonClusterToolsTests
{
    private static CodonUsageTableInput EColi() => new(Preset: "EColiK12");

    [Test]
    public void CodonMinMaxProfile_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => MolToolsTools.codon_minmax_profile("AGAAGAAGA", EColi(), 3));
        Assert.Throws<ArgumentException>(() => MolToolsTools.codon_minmax_profile("", EColi()));
        Assert.Throws<ArgumentException>(() => MolToolsTools.codon_minmax_profile(null!, EColi()));
        Assert.Throws<ArgumentNullException>(() => MolToolsTools.codon_minmax_profile("AGAAGAAGA", null!));
        Assert.Throws<ArgumentOutOfRangeException>(() => MolToolsTools.codon_minmax_profile("AGAAGAAGA", EColi(), 0));
    }

    // MM4 of the unit tests: 5 codons, window 3 -> starts 0,1,2; window 2 (3x AGA) = -86.3636..%.
    [Test]
    public void CodonMinMaxProfile_MatchesLockedApiValues()
    {
        var windows = MolToolsTools.codon_minmax_profile("CGCCGCAGAAGAAGA", EColi(), window_size: 3).Windows;
        var api = CodonOptimizer.CalculateMinMaxProfile("CGCCGCAGAAGAAGA", CodonOptimizer.EColiK12, 3);

        Assert.Multiple(() =>
        {
            Assert.That(windows.Select(w => w.WindowStartCodon), Is.EqualTo(new[] { 0, 1, 2 }));
            Assert.That(windows[2].PercentMinMax, Is.EqualTo(-86.36363636363637).Within(1e-10));
            Assert.That(windows.Select(w => w.PercentMinMax), Is.EqualTo(api.Select(w => w.PercentMinMax)));
        });

        // MM2 (+100 for an all-most-common window) and MM3 (mixed CUG+AGA = 36.4705..%).
        Assert.That(MolToolsTools.codon_minmax_profile("CGCCGCCGC", EColi(), 3).Windows[0].PercentMinMax,
            Is.EqualTo(100.0).Within(1e-10));
        Assert.That(MolToolsTools.codon_minmax_profile("CUGAGA", EColi(), 2).Windows[0].PercentMinMax,
            Is.EqualTo(36.470588235294116).Within(1e-10));
    }

    [Test]
    public void CodonMinMaxProfile_DefaultWindowIs18()
    {
        // 2 codons < default 18 -> empty, as the C# default.
        Assert.That(MolToolsTools.codon_minmax_profile("AGAAGA", EColi()).Windows, Is.Empty);
    }

    [Test]
    public void FindRareCodonClusters_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => MolToolsTools.find_rare_codon_clusters("AGAAGAAGA", EColi()));
        Assert.Throws<ArgumentException>(() => MolToolsTools.find_rare_codon_clusters("", EColi()));
        Assert.Throws<ArgumentException>(() => MolToolsTools.find_rare_codon_clusters(null!, EColi()));
        Assert.Throws<ArgumentNullException>(() => MolToolsTools.find_rare_codon_clusters("AGA", null!));
        Assert.Throws<ArgumentOutOfRangeException>(() => MolToolsTools.find_rare_codon_clusters("AGA", EColi(), window_size: 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => MolToolsTools.find_rare_codon_clusters("AGA", EColi(), min_rare_codons: 0));
    }

    // C-two-runs of the unit tests: 7xAGA + 8xCGC + 7xAGA (Sherlocc defaults 0.15 / 7 / 4)
    // -> clusters 0..9 (7 rare) and 12..21 (7 rare).
    [Test]
    public void FindRareCodonClusters_MatchesLockedApiValues()
    {
        string run = string.Concat(Enumerable.Repeat("AGA", 7));
        string sequence = run + string.Concat(Enumerable.Repeat("CGC", 8)) + run;

        var clusters = MolToolsTools.find_rare_codon_clusters(sequence, EColi()).Clusters;
        var api = CodonOptimizer.FindRareCodonClusters(sequence, CodonOptimizer.EColiK12);

        Assert.Multiple(() =>
        {
            Assert.That(clusters, Is.EqualTo(new[]
            {
                new RareCodonClusterDto(0, 9, 7),
                new RareCodonClusterDto(12, 21, 7),
            }));
            Assert.That(clusters.Select(c => (c.StartCodon, c.EndCodon, c.RareCount)),
                Is.EqualTo(api.Select(c => (c.StartCodon, c.EndCodon, c.RareCount))));
        });

        // Custom 3/3 rule on 3x AGA (unit test: one cluster, 3 rare).
        var custom = MolToolsTools.find_rare_codon_clusters("AGAAGAAGA", EColi(), window_size: 3, min_rare_codons: 3).Clusters;
        Assert.That(custom, Is.EqualTo(new[] { new RareCodonClusterDto(0, 2, 3) }));
    }
}
