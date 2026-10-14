// CODON-RARE-001 — MCP tool find_conserved_rare_codon_clusters (FIN-B02 WP15, audit A2-3; B02 F49).
// Delegates to CodonOptimizer.FindConservedRareCodonClusters (F40, Sherlocc 1.0 stage 3). The expected
// values are the sherlocc.pl-locked cases A/B/C of
// tests/Seqeron/Seqeron.Genomics.Tests/Unit/MolTools/CodonOptimizer_SherloccConservedClusters_Tests.cs
// (sherlocc.pl, github.com/mtthchrtr/sherlocc commit b3b4d15, run with perl 5.38 on these alignments).
using NUnit.Framework;
using Seqeron.Genomics.MolTools;
using Seqeron.Mcp.MolTools.Models;
using Seqeron.Mcp.MolTools.Tools;

namespace Seqeron.Mcp.MolTools.Tests;

[TestFixture]
public class FindConservedRareCodonClustersToolTests
{
    private static readonly string[] AlignmentA =
    {
        "GAAAAAGATGAAAAAGATGAACGGCGCCGACCGCTCGGGGCGGAAAAAGATGAAAAAGAT",
        "GAAAAA---GAAAAGGATGAACGCCGGCCGCGAGGGCTCGCAGAAAAAGATGAA---GAT",
        "GAGAAAGATGAAAAAGACGAACGACGGCGCCTCCCGGCGGGGGAAAATGATGAAAAAGAT",
    };

    private static readonly string[] AlignmentB =
    {
        "GAAAAAGATCGGCGCCGA---CCGCTCGGGGCGCGG---CGACGCGAAAAAGAT",
        "GAAAAAGATCGACGGCGC---CTCCCGGCGGGGCGC---CGGCGAGAAAAAGAT",
    };

    private static readonly string[] AlignmentC =
    {
        "GAAAAACGGCGCCGACCGCTCGGGGCGCGGCGAGAAAAAGATGAAAAAGATGAAAAACGACGGCCGCGCGCGCTCGGGCGACGCGAAAAAGATGAAAAAGAT",
        "GAGAAGCGTTCGCCGACGGCGCGTTCGACGCCGGAGAAGCAGGAGAAGCAGGAGAAGTCGCGTACG---CCGGCGCGTTCGACGGAGAAGCAGGAGAAGCAG",
        "GAAAAACGACGGCCGCGCGCGCTCGGGCGACGCGAAAAAGATGAAAAAGATGAAAAACGGCGCCGACCGCTCGGGGCGCGGCGAGAAAAAGATGAAAAAGAT",
    };

    private const double Tol = 1e-12;
    private static CodonUsageTableInput Preset(string id) => new(Preset: id);

    private static void AssertEqualsApi(ConservedRareCodonClustersResult mcp, CodonOptimizer.SherloccResult api)
    {
        Assert.Multiple(() =>
        {
            Assert.That(mcp.SequenceCount, Is.EqualTo(api.SequenceCount));
            Assert.That(mcp.AlignmentLength, Is.EqualTo(api.AlignmentLength));
            Assert.That(mcp.ColumnAverages, Is.EqualTo(api.ColumnAverages));
            Assert.That(mcp.ColumnUsageCounts, Is.EqualTo(api.ColumnUsageCounts));
            Assert.That(mcp.PauseColumns, Is.EqualTo(api.PauseColumns));
            Assert.That(mcp.Clusters, Is.EqualTo(api.Clusters.Select(c =>
                new ConservedRareCodonClusterDto(c.StartColumn, c.EndColumn, c.PauseCount, c.MeanUsage, c.MiddlePoint))));
            Assert.That(mcp.TotalClusterLength, Is.EqualTo(api.TotalClusterLength));
            Assert.That(mcp.ClusterFraction, Is.EqualTo(api.ClusterFraction));
            Assert.That(mcp.MeanClusterUsage, Is.EqualTo(api.MeanClusterUsage));
        });
    }

    [Test]
    public void AlignmentA_MatchesSherlocc_AndCSharp()
    {
        var r = MolToolsTools.find_conserved_rare_codon_clusters(AlignmentA, 15, Preset("Yeast"));

        // Perl: total_nb_of_pauses=4; clusters=9-12|9.83928571428571.
        Assert.Multiple(() =>
        {
            Assert.That(r.SequenceCount, Is.EqualTo(3));
            Assert.That(r.AlignmentLength, Is.EqualTo(20));
            Assert.That(r.ColumnUsageCounts, Is.EqualTo(new[] { 11, 14, 17, 20, 20, 20, 21, 21, 21, 21, 21, 21, 21, 21, 21, 20, 20, 17, 14, 11 }));
            Assert.That(r.ColumnAverages[0], Is.Null);
            Assert.That(r.ColumnAverages[8]!.Value, Is.EqualTo(13.657142857142858).Within(Tol));
            Assert.That(r.PauseColumns, Is.EqualTo(new[] { 8, 9, 10, 11 }));
            Assert.That(r.Clusters, Has.Count.EqualTo(1));
            Assert.That((r.Clusters[0].StartColumn + 1, r.Clusters[0].EndColumn + 1), Is.EqualTo((9, 12)));
            Assert.That(r.Clusters[0].PauseCount, Is.EqualTo(4));
            Assert.That(r.Clusters[0].MeanUsage, Is.EqualTo(9.839285714285714).Within(Tol));
            Assert.That(r.Clusters[0].MiddlePoint, Is.EqualTo(10));
            Assert.That(r.TotalClusterLength, Is.EqualTo(4));
            Assert.That(r.ClusterFraction, Is.EqualTo(0.2).Within(Tol));
            Assert.That(r.MeanClusterUsage!.Value, Is.EqualTo(9.839285714285714).Within(Tol));
        });
        AssertEqualsApi(r, CodonOptimizer.FindConservedRareCodonClusters(AlignmentA, CodonOptimizer.Yeast, 15));
    }

    [Test]
    public void AlignmentB_NonRepresentativeColumnInsideCluster_MatchesSherlocc_AndCSharp()
    {
        var r = MolToolsTools.find_conserved_rare_codon_clusters(AlignmentB, 15, Preset("Yeast"));

        // Perl: pauses 5,6,7,8,10,11,12 (0-based; column 9 "NA"); clusters=6-13|5.265625.
        Assert.Multiple(() =>
        {
            Assert.That(r.ColumnUsageCounts[9], Is.EqualTo(10));
            Assert.That(r.ColumnAverages[9], Is.Null);
            Assert.That(r.PauseColumns, Is.EqualTo(new[] { 5, 6, 7, 8, 10, 11, 12 }));
            Assert.That(r.Clusters, Has.Count.EqualTo(1));
            Assert.That((r.Clusters[0].StartColumn + 1, r.Clusters[0].EndColumn + 1), Is.EqualTo((6, 13)));
            Assert.That(r.Clusters[0].PauseCount, Is.EqualTo(7));
            Assert.That(r.Clusters[0].MeanUsage, Is.EqualTo(5.265625).Within(Tol));
            Assert.That(r.Clusters[0].MiddlePoint, Is.EqualTo(9));
            Assert.That(r.ClusterFraction, Is.EqualTo(8.0 / 18).Within(Tol));
        });
        AssertEqualsApi(r, CodonOptimizer.FindConservedRareCodonClusters(AlignmentB, CodonOptimizer.Yeast, 15));
    }

    [Test]
    public void AlignmentC_PerSequenceTables_MatchesSherlocc_AndCSharp()
    {
        var r = MolToolsTools.find_conserved_rare_codon_clusters(
            AlignmentC, 16, target_organisms_per_sequence: new[] { Preset("Yeast"), Preset("Human"), Preset("Yeast") });

        // Perl: clusters=4-10|9.20962052101093,21-27|9.30926559372343; clusters_freq_avrg=9.25944305736718.
        Assert.Multiple(() =>
        {
            Assert.That(r.PauseColumns, Is.EqualTo(new[] { 3, 4, 5, 6, 7, 8, 9, 20, 21, 22, 23, 24, 25, 26 }));
            Assert.That(r.Clusters.Select(c => (c.StartColumn + 1, c.EndColumn + 1)), Is.EqualTo(new[] { (4, 10), (21, 27) }));
            Assert.That(r.Clusters[0].MeanUsage, Is.EqualTo(9.209620521010933).Within(Tol));
            Assert.That(r.Clusters[1].MeanUsage, Is.EqualTo(9.309265593723435).Within(Tol));
            Assert.That(r.Clusters.Select(c => c.MiddlePoint), Is.EqualTo(new[] { 6, 23 }));
            Assert.That(r.MeanClusterUsage!.Value, Is.EqualTo(9.25944305736718).Within(1e-13));
            Assert.That(r.TotalClusterLength, Is.EqualTo(14));
            Assert.That(r.ClusterFraction, Is.EqualTo(14.0 / 34).Within(Tol));
        });
        AssertEqualsApi(r, CodonOptimizer.FindConservedRareCodonClusters(
            AlignmentC, new[] { CodonOptimizer.Yeast, CodonOptimizer.Human, CodonOptimizer.Yeast }, 16));
    }

    [Test]
    public void InlineTable_ValuesAreUsedAsPerThousand()
    {
        // The Yeast preset's per-thousand column supplied inline (DNA keys) reproduces the preset result.
        var inline = new CodonUsageTableInput(
            OrganismName: "yeast per thousand",
            CodonFrequencies: CodonOptimizer.Yeast.PerThousand!.ToDictionary(kv => kv.Key.Replace('U', 'T'), kv => kv.Value));
        var viaInline = MolToolsTools.find_conserved_rare_codon_clusters(AlignmentA, 15, inline);
        var viaPreset = MolToolsTools.find_conserved_rare_codon_clusters(AlignmentA, 15, Preset("Yeast"));

        Assert.Multiple(() =>
        {
            Assert.That(viaInline.ColumnAverages, Is.EqualTo(viaPreset.ColumnAverages));
            Assert.That(viaInline.Clusters, Is.EqualTo(viaPreset.Clusters));
            Assert.That(viaInline.MeanClusterUsage, Is.EqualTo(viaPreset.MeanClusterUsage));
        });
    }

    [Test]
    public void UseCodonFrequencies_ScoresFractions_EColiDefaultScoresPerThousand()
    {
        string[] aln = { string.Concat(Enumerable.Repeat("AGA", 7)), string.Concat(Enumerable.Repeat("AGG", 7)) };
        // SRC RO1 (2026-10-14): was Assert.Throws (EColiK12 had no per-thousand column). The preset now
        // carries the Kazusa 316407 /1000 column (AGA 2.0, AGG 1.1) → column 3 averages 1.55, pauses 2–4.
        var perThousand = MolToolsTools.find_conserved_rare_codon_clusters(aln, 15, Preset("EColiK12"));
        Assert.Multiple(() =>
        {
            Assert.That(perThousand.ColumnAverages[3]!.Value, Is.EqualTo(1.55).Within(Tol));
            Assert.That(perThousand.PauseColumns, Is.EqualTo(new[] { 2, 3, 4 }));
            Assert.That(perThousand.Clusters, Is.Empty);
        });
        AssertEqualsApi(perThousand, CodonOptimizer.FindConservedRareCodonClusters(aln, CodonOptimizer.EColiK12, 15));

        // E. coli fractions: AGA 0.04, AGG 0.02 → column 3 averages 0.03 (C# lock).
        var r = MolToolsTools.find_conserved_rare_codon_clusters(aln, 0.15, Preset("EColiK12"), use_codon_frequencies: true);
        Assert.Multiple(() =>
        {
            Assert.That(r.ColumnAverages[3]!.Value, Is.EqualTo(0.03).Within(Tol));
            Assert.That(r.PauseColumns, Is.EqualTo(new[] { 2, 3, 4 }));
            Assert.That(r.Clusters, Is.Empty);
        });
        AssertEqualsApi(r, CodonOptimizer.FindConservedRareCodonClusters(
            aln, CodonOptimizer.EColiK12, 0.15, new CodonOptimizer.SherloccOptions { UseCodonFrequencies = true }));
    }

    [Test]
    public void Options_AreForwarded()
    {
        var options = new CodonOptimizer.SherloccOptions { AveragingWindow = 5, ClusterWindow = 5, MinPauses = 3, MinCoverage = 0.6 };
        var r = MolToolsTools.find_conserved_rare_codon_clusters(
            AlignmentA, 15, Preset("Yeast"), averaging_window: 5, cluster_window: 5, min_pauses: 3, min_coverage: 0.6);
        AssertEqualsApi(r, CodonOptimizer.FindConservedRareCodonClusters(AlignmentA, CodonOptimizer.Yeast, 15, options));
    }

    [Test]
    public void InvalidInput_Throws()
    {
        var y = Preset("Yeast");
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentException>(() => MolToolsTools.find_conserved_rare_codon_clusters(null!, 15, y));
            Assert.Throws<ArgumentException>(() => MolToolsTools.find_conserved_rare_codon_clusters(Array.Empty<string>(), 15, y));
            Assert.Throws<ArgumentException>(() => MolToolsTools.find_conserved_rare_codon_clusters(AlignmentA, 15));
            Assert.Throws<ArgumentException>(() => MolToolsTools.find_conserved_rare_codon_clusters(
                AlignmentA, 15, y, new[] { y, y, y }));
            Assert.Throws<ArgumentException>(() => MolToolsTools.find_conserved_rare_codon_clusters(
                AlignmentA, 15, target_organisms_per_sequence: new[] { y }));
            Assert.Throws<ArgumentException>(() => MolToolsTools.find_conserved_rare_codon_clusters(new[] { "AAAA" }, 15, y));
            Assert.Throws<ArgumentOutOfRangeException>(() => MolToolsTools.find_conserved_rare_codon_clusters(AlignmentA, double.NaN, y));
            Assert.Throws<ArgumentOutOfRangeException>(() => MolToolsTools.find_conserved_rare_codon_clusters(AlignmentA, 15, y, averaging_window: 6));
            Assert.Throws<ArgumentOutOfRangeException>(() => MolToolsTools.find_conserved_rare_codon_clusters(AlignmentA, 15, y, min_coverage: 1.5));
            Assert.Throws<ArgumentOutOfRangeException>(() => MolToolsTools.find_conserved_rare_codon_clusters(AlignmentA, 15, y, min_pauses: 0));
            Assert.Throws<ArgumentException>(() => MolToolsTools.find_conserved_rare_codon_clusters(AlignmentA, 15, Preset("Martian")));
        });
    }
}
