// CODON-RARE-001 — review 2026-10 (FIN-B02 WP8, audit A1-2): Sherlocc conserved rare-codon clusters
// over a codon alignment (family averaging, 75 % representativeness filter, cluster table).
// Evidence: docs/Evidence/CODON-RARE-001-Evidence.md; TestSpec: tests/TestSpecs/CODON-RARE-001.md (S1–S12)
// Source / oracle: Sherlocc 1.0 sherlocc.pl (github.com/mtthchrtr/sherlocc, commit b3b4d15), stages 2+3 run
//   with perl 5.38 on the three alignments below (stage-1 records written in the l. 369–374 format, the
//   two stage-1 directory resets l. 15/17 removed — stage 1 needs live UniProt/ENA queries); files/cut
//   rows = the Yeast / Human preset per-thousand values (taxid 4932 / 9606). Perl printed cluster
//   averages with %.15g and column averages with %.2f (HTML "Average" row); the unrounded column
//   averages come from an independent Python port of l. 772–803 that reproduces every Perl value.

namespace Seqeron.Genomics.Tests.Unit.MolTools;

[TestFixture]
public class CodonOptimizer_SherloccConservedClusters_Tests
{
    // A: yeast, 3 rows × 20 columns, rare Arg/Pro/Leu/Gly/Ala stretch in columns 7–13; two gaps.
    private static readonly string[] AlignmentA =
    {
        "GAAAAAGATGAAAAAGATGAACGGCGCCGACCGCTCGGGGCGGAAAAAGATGAAAAAGAT",
        "GAAAAA---GAAAAGGATGAACGCCGGCCGCGAGGGCTCGCAGAAAAAGATGAA---GAT",
        "GAGAAAGATGAAAAAGACGAACGACGGCGCCTCCCGGCGGGGGAAAATGATGAAAAAGAT",
    };

    // B: yeast, 2 rows × 18 columns; columns 6 and 12 gapped in both rows, so only column 9
    // (window 6..12: 10 valued cells < 0.75 × 2 × 7 = 10.5) is non-representative, inside the cluster.
    private static readonly string[] AlignmentB =
    {
        "GAAAAAGATCGGCGCCGA---CCGCTCGGGGCGCGG---CGACGCGAAAAAGAT",
        "GAAAAAGATCGACGGCGC---CTCCCGGCGGGGCGC---CGGCGAGAAAAAGAT",
    };

    // C: per-sequence tables yeast / human / yeast, 3 rows × 34 columns, two separate rare stretches.
    private static readonly string[] AlignmentC =
    {
        "GAAAAACGGCGCCGACCGCTCGGGGCGCGGCGAGAAAAAGATGAAAAAGATGAAAAACGACGGCCGCGCGCGCTCGGGCGACGCGAAAAAGATGAAAAAGAT",
        "GAGAAGCGTTCGCCGACGGCGCGTTCGACGCCGGAGAAGCAGGAGAAGCAGGAGAAGTCGCGTACG---CCGGCGCGTTCGACGGAGAAGCAGGAGAAGCAG",
        "GAAAAACGACGGCCGCGCGCGCTCGGGCGACGCGAAAAAGATGAAAAAGATGAAAAACGGCGCCGACCGCTCGGGGCGCGGCGAGAAAAAGATGAAAAAGAT",
    };

    private const double Tol = 1e-12;

    [Test]
    public void AlignmentA_ColumnAverages_MatchSherlocc()
    {
        var r = CodonOptimizer.FindConservedRareCodonClusters(AlignmentA, CodonOptimizer.Yeast, 15);

        double?[] expected =
        {
            null, null, 38.71176470588236, 39.74500000000001, 34.59, 28.605000000000008, 24.180952380952384,
            18.31904761904762, 13.657142857142858, 9.95238095238095, 4.79047619047619, 10.957142857142857,
            16.361904761904764, 21.214285714285715, 27.076190476190476, 31.785000000000004, 36.54500000000001,
            41.32352941176471, null, null,
        };
        Assert.Multiple(() =>
        {
            Assert.That(r.SequenceCount, Is.EqualTo(3));
            Assert.That(r.AlignmentLength, Is.EqualTo(20));
            Assert.That(r.ColumnUsageCounts, Is.EqualTo(new[] { 11, 14, 17, 20, 20, 20, 21, 21, 21, 21, 21, 21, 21, 21, 21, 20, 20, 17, 14, 11 }));
            for (int i = 0; i < expected.Length; i++)
            {
                if (expected[i] is null)
                    Assert.That(r.ColumnAverages[i], Is.Null, $"column {i} (Sherlocc NA)");
                else
                    Assert.That(r.ColumnAverages[i]!.Value, Is.EqualTo(expected[i]!.Value).Within(Tol), $"column {i}");
            }

            // Perl HTML "Average" row, 2 dp.
            Assert.That(Math.Round(r.ColumnAverages[8]!.Value, 2), Is.EqualTo(13.66));
            Assert.That(Math.Round(r.ColumnAverages[11]!.Value, 2), Is.EqualTo(10.96));
        });
    }

    [Test]
    public void AlignmentA_OneMergedCluster_MatchesSherloccResults()
    {
        var r = CodonOptimizer.FindConservedRareCodonClusters(AlignmentA, CodonOptimizer.Yeast, 15);

        // Perl: total_nb_of_pauses=4; clusters=9-12|9.83928571428571 (four qualifying windows merged).
        Assert.Multiple(() =>
        {
            Assert.That(r.PauseColumns, Is.EqualTo(new[] { 8, 9, 10, 11 }));
            Assert.That(r.Clusters, Has.Count.EqualTo(1));
            var c = r.Clusters[0];
            Assert.That((c.StartColumn + 1, c.EndColumn + 1), Is.EqualTo((9, 12)));
            Assert.That(c.PauseCount, Is.EqualTo(4));
            Assert.That(c.MeanUsage, Is.EqualTo(9.839285714285714).Within(Tol));
            Assert.That(c.MiddlePoint, Is.EqualTo(10));
            Assert.That(r.TotalClusterLength, Is.EqualTo(4));
            Assert.That(r.ClusterFraction, Is.EqualTo(0.2).Within(Tol));
            Assert.That(r.MeanClusterUsage!.Value, Is.EqualTo(9.839285714285714).Within(Tol));
        });
    }

    [Test]
    public void AlignmentB_NonRepresentativeColumnInsideCluster_CountsAsZero()
    {
        var r = CodonOptimizer.FindConservedRareCodonClusters(AlignmentB, CodonOptimizer.Yeast, 15);

        // Perl: pauses 5,6,7,8,10,11,12 (0-based; column 9 "NA"); clusters=6-13|5.265625.
        Assert.Multiple(() =>
        {
            Assert.That(r.ColumnUsageCounts[9], Is.EqualTo(10));
            Assert.That(r.ColumnAverages[9], Is.Null);
            Assert.That(r.ColumnAverages[2]!.Value, Is.EqualTo(22.066666666666666).Within(Tol));
            Assert.That(r.ColumnAverages[12]!.Value, Is.EqualTo(10.85).Within(Tol));
            Assert.That(r.ColumnAverages[15], Is.Null);
            Assert.That(r.PauseColumns, Is.EqualTo(new[] { 5, 6, 7, 8, 10, 11, 12 }));
            Assert.That(r.Clusters, Has.Count.EqualTo(1));
            var c = r.Clusters[0];
            Assert.That((c.StartColumn + 1, c.EndColumn + 1), Is.EqualTo((6, 13)));
            Assert.That(c.PauseCount, Is.EqualTo(7));
            // Column 9 contributes 0 (Perl undefined value); excluding it would give 6.017857….
            Assert.That(c.MeanUsage, Is.EqualTo(5.265625).Within(Tol));
            Assert.That(c.MiddlePoint, Is.EqualTo(9));
            Assert.That(r.ClusterFraction, Is.EqualTo(8.0 / 18).Within(Tol));
        });
    }

    [Test]
    public void AlignmentC_PerSequenceTables_TwoClusters_MatchSherlocc()
    {
        var tables = new[] { CodonOptimizer.Yeast, CodonOptimizer.Human, CodonOptimizer.Yeast };
        var r = CodonOptimizer.FindConservedRareCodonClusters(AlignmentC, tables, 16);

        // Perl: total_nb_of_pauses=14; clusters=4-10|9.20962052101093,21-27|9.30926559372343;
        // clusters_freq_avrg=9.25944305736718.
        Assert.Multiple(() =>
        {
            Assert.That(r.PauseColumns, Is.EqualTo(new[] { 3, 4, 5, 6, 7, 8, 9, 20, 21, 22, 23, 24, 25, 26 }));
            Assert.That(r.ColumnAverages[2]!.Value, Is.EqualTo(16.31023785563948).Within(Tol));
            Assert.That(r.ColumnAverages[20]!.Value, Is.EqualTo(15.25921407007553).Within(Tol));
            Assert.That(r.ColumnAverages[27]!.Value, Is.EqualTo(19.62654223857159).Within(Tol));
            Assert.That(r.Clusters.Select(c => (c.StartColumn + 1, c.EndColumn + 1)), Is.EqualTo(new[] { (4, 10), (21, 27) }));
            Assert.That(r.Clusters[0].MeanUsage, Is.EqualTo(9.209620521010933).Within(Tol));
            Assert.That(r.Clusters[1].MeanUsage, Is.EqualTo(9.309265593723435).Within(Tol));
            Assert.That(r.Clusters.Select(c => c.MiddlePoint), Is.EqualTo(new[] { 6, 23 }));
            Assert.That(r.Clusters.Select(c => c.PauseCount), Is.EqualTo(new[] { 7, 7 }));
            Assert.That(r.MeanClusterUsage!.Value, Is.EqualTo(9.25944305736718).Within(1e-13));
            Assert.That(r.TotalClusterLength, Is.EqualTo(14));
            Assert.That(r.ClusterFraction, Is.EqualTo(14.0 / 34).Within(Tol));
        });
    }

    [Test]
    public void Threshold_IsInclusive_PerSherlocc()
    {
        // Column 8 average 13.657142857142858: a pause at exactly that threshold (l. 803 `<=`),
        // not just below it — then only 3 pauses remain and no window reaches 4.
        double col8 = CodonOptimizer.FindConservedRareCodonClusters(AlignmentA, CodonOptimizer.Yeast, 15).ColumnAverages[8]!.Value;
        var at = CodonOptimizer.FindConservedRareCodonClusters(AlignmentA, CodonOptimizer.Yeast, col8);
        var below = CodonOptimizer.FindConservedRareCodonClusters(AlignmentA, CodonOptimizer.Yeast, Math.BitDecrement(col8));
        Assert.Multiple(() =>
        {
            Assert.That(at.PauseColumns, Is.EqualTo(new[] { 8, 9, 10, 11 }));
            Assert.That(at.Clusters, Has.Count.EqualTo(1));
            Assert.That(below.PauseColumns, Is.EqualTo(new[] { 9, 10, 11 }));
            Assert.That(below.Clusters, Is.Empty);
            Assert.That(below.MeanClusterUsage, Is.Null);
            Assert.That(below.ClusterFraction, Is.EqualTo(0.0));
        });
    }

    [Test]
    public void SingleTableOverload_EqualsPerSequenceOverload_AndAcceptsRnaLowerCase()
    {
        var single = CodonOptimizer.FindConservedRareCodonClusters(AlignmentA, CodonOptimizer.Yeast, 15);
        var perSeq = CodonOptimizer.FindConservedRareCodonClusters(
            AlignmentA, new[] { CodonOptimizer.Yeast, CodonOptimizer.Yeast, CodonOptimizer.Yeast }, 15);
        var rna = CodonOptimizer.FindConservedRareCodonClusters(
            AlignmentA.Select(s => s.Replace('T', 'U').ToLowerInvariant()).ToArray(), CodonOptimizer.Yeast, 15);
        var dotGaps = CodonOptimizer.FindConservedRareCodonClusters(
            AlignmentA.Select(s => s.Replace('-', '.')).ToArray(), CodonOptimizer.Yeast, 15);
        foreach (var other in new[] { perSeq, rna, dotGaps })
        {
            Assert.That(other.ColumnAverages, Is.EqualTo(single.ColumnAverages));
            Assert.That(other.Clusters, Is.EqualTo(single.Clusters));
        }
    }

    [Test]
    public void Defaults_AreSherloccConstants()
    {
        var o = new CodonOptimizer.SherloccOptions();
        Assert.Multiple(() =>
        {
            Assert.That(o.AveragingWindow, Is.EqualTo(7));
            Assert.That(o.ClusterWindow, Is.EqualTo(7));
            Assert.That(o.MinPauses, Is.EqualTo(4));
            Assert.That(o.MinCoverage, Is.EqualTo(0.75));
            Assert.That(o.UseCodonFrequencies, Is.False);
        });
    }

    [Test]
    public void TableWithoutPerThousand_Rejected_UnlessCodonFrequenciesRequested()
    {
        string[] aln = { string.Concat(Enumerable.Repeat("AGA", 7)), string.Concat(Enumerable.Repeat("AGG", 7)) };
        Assert.Throws<ArgumentException>(() => CodonOptimizer.FindConservedRareCodonClusters(aln, CodonOptimizer.EColiK12, 15));

        // E. coli fractions: AGA 0.04, AGG 0.02 → column 3 (the only full window) averages 0.03.
        var r = CodonOptimizer.FindConservedRareCodonClusters(
            aln, CodonOptimizer.EColiK12, 0.15, new CodonOptimizer.SherloccOptions { UseCodonFrequencies = true });
        Assert.Multiple(() =>
        {
            Assert.That(r.ColumnAverages[3]!.Value, Is.EqualTo(0.03).Within(Tol));
            Assert.That(r.PauseColumns, Is.EqualTo(new[] { 2, 3, 4 }));
            Assert.That(r.Clusters, Is.Empty); // only 3 representative pauses in the single 7-column window
        });
    }

    [Test]
    public void InvalidInput_Throws()
    {
        var y = CodonOptimizer.Yeast;
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentNullException>(() => CodonOptimizer.FindConservedRareCodonClusters(null!, y, 15));
            Assert.Throws<ArgumentException>(() => CodonOptimizer.FindConservedRareCodonClusters(Array.Empty<string>(), y, 15));
            Assert.Throws<ArgumentException>(() => CodonOptimizer.FindConservedRareCodonClusters(new[] { "AAAA" }, y, 15));
            Assert.Throws<ArgumentException>(() => CodonOptimizer.FindConservedRareCodonClusters(new[] { "AAA", "AAAAAA" }, y, 15));
            Assert.Throws<ArgumentException>(() => CodonOptimizer.FindConservedRareCodonClusters(new[] { "AAA", "AAA" }, new[] { y }, 15));
            Assert.Throws<ArgumentOutOfRangeException>(() => CodonOptimizer.FindConservedRareCodonClusters(new[] { "AAA" }, y, double.NaN));
            Assert.Throws<ArgumentOutOfRangeException>(() => CodonOptimizer.FindConservedRareCodonClusters(
                new[] { "AAA" }, y, 15, new CodonOptimizer.SherloccOptions { AveragingWindow = 6 }));
            Assert.Throws<ArgumentOutOfRangeException>(() => CodonOptimizer.FindConservedRareCodonClusters(
                new[] { "AAA" }, y, 15, new CodonOptimizer.SherloccOptions { MinCoverage = 1.5 }));
            Assert.Throws<ArgumentOutOfRangeException>(() => CodonOptimizer.FindConservedRareCodonClusters(
                new[] { "AAA" }, y, 15, new CodonOptimizer.SherloccOptions { MinPauses = 0 }));
        });
    }

    [Test]
    public void ShortAlignment_NoClusterWindow_NoClusters()
    {
        // Fewer columns than the cluster window: Perl's window loop (l. 838) does not run.
        var r = CodonOptimizer.FindConservedRareCodonClusters(new[] { "CGGCGGCGGCGG" }, CodonOptimizer.Yeast, 15);
        Assert.Multiple(() =>
        {
            Assert.That(r.AlignmentLength, Is.EqualTo(4));
            Assert.That(r.ColumnAverages, Is.All.Null); // ≤ 4 of 7 cells < 5.25
            Assert.That(r.Clusters, Is.Empty);
        });
    }
}
