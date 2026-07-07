// ASSEMBLY-STATS-001 — Assembly Statistics (N50 / L50 / Nx / Lx / auN, gaps)
// Evidence: docs/Evidence/ASSEMBLY-STATS-001-Evidence.md
// TestSpec: tests/TestSpecs/ASSEMBLY-STATS-001.md
// Source: Miller JR, Koren S, Sutton G (2010). Genomics 95(6):315-327, §1.2;
//         Wikipedia "N50, L50, and related statistics" (worked example);
//         QUAST quast_libs/N50.py (cumulative >= threshold% convention);
//         Li H (2020), auN: a new metric to measure assembly contiguity (auN = Sum(l^2)/Sum(l)).

namespace Seqeron.Genomics.Tests.Unit.Chromosome;

[TestFixture]
public class GenomeAssemblyAnalyzer_AssemblyStatistics_Tests
{
    // Wikipedia worked example "Assembly A": lengths 80,70,50,40,30,20 (kbp), total 290.
    // Source: Wikipedia N50 article — N50 = 70, L50 = 2.
    private static readonly int[] AssemblyA = { 80, 70, 50, 40, 30, 20 };

    // Wikipedia "Assembly B": Assembly A plus 10 and 5, total 305 — N50 = 50, L50 = 3.
    private static readonly int[] AssemblyB = { 80, 70, 50, 40, 30, 20, 10, 5 };

    private static List<(string Id, string Sequence)> AsContigs(IEnumerable<int> lengths) =>
        lengths.Select((len, i) => ($"contig{i + 1}", new string('G', len))).ToList();

    #region CalculateNx (core overload)

    // M1 — Assembly A {80,70,50,40,30,20}, total 290: 80+70=150 >= 145 (50% of 290) => N50=70, L50=2.
    // Source: Wikipedia worked example; Miller 2010 §1.2.
    [Test]
    public void CalculateNx_AssemblyAThreshold50_ReturnsN50_70_L50_2()
    {
        var sorted = AssemblyA.OrderByDescending(l => l).ToList();
        long total = sorted.Sum(l => (long)l);

        var n50 = GenomeAssemblyAnalyzer.CalculateNx(sorted, total, 50);

        Assert.Multiple(() =>
        {
            Assert.That(n50.Nx, Is.EqualTo(70),
                "N50 of {80,70,50,40,30,20} is 70: 80+70=150 reaches >=50% of 290 (=145) at the second contig (Wikipedia; Miller 2010).");
            Assert.That(n50.Lx, Is.EqualTo(2),
                "L50 is the count of contigs in that prefix = 2 (Wikipedia).");
            Assert.That(n50.CumulativeLength, Is.EqualTo(150),
                "Cumulative length at the N50 contig is 80+70=150.");
        });
    }

    // M3 — Assembly A, threshold 90: 90% of 290 = 261; cumulative 80+70+50+40+30=270 >= 261 => N90=30, L90=5.
    // Source: Wikipedia N90 definition + worked-example lengths.
    [Test]
    public void CalculateNx_AssemblyAThreshold90_ReturnsN90_30_L90_5()
    {
        var sorted = AssemblyA.OrderByDescending(l => l).ToList();
        long total = sorted.Sum(l => (long)l);

        var n90 = GenomeAssemblyAnalyzer.CalculateNx(sorted, total, 90);

        Assert.Multiple(() =>
        {
            Assert.That(n90.Nx, Is.EqualTo(30),
                "N90 of {80,70,50,40,30,20} is 30: cumulative 270 first reaches >=90% of 290 (=261) at the fifth contig (Wikipedia).");
            Assert.That(n90.Lx, Is.EqualTo(5),
                "L90 is 5 contigs.");
        });
    }

    // M4 — Monotonicity: raising the threshold cannot increase Nx (N90 <= N50) nor decrease Lx (L90 >= L50).
    // Source: Wikipedia ("N90 is less than or equal to the N50").
    [Test]
    public void CalculateNx_HigherThreshold_NxNonIncreasingAndLxNonDecreasing()
    {
        var sorted = AssemblyA.OrderByDescending(l => l).ToList();
        long total = sorted.Sum(l => (long)l);

        var n50 = GenomeAssemblyAnalyzer.CalculateNx(sorted, total, 50);
        var n90 = GenomeAssemblyAnalyzer.CalculateNx(sorted, total, 90);

        Assert.Multiple(() =>
        {
            Assert.That(n90.Nx, Is.LessThanOrEqualTo(n50.Nx),
                "INV-03: N90 (30) must be <= N50 (70) because a larger threshold extends the prefix to shorter contigs.");
            Assert.That(n90.Lx, Is.GreaterThanOrEqualTo(n50.Lx),
                "INV-03: L90 (5) must be >= L50 (2).");
        });
    }

    // M5 — Inclusive boundary: lengths {50,50}, total 100; cumulative 50 equals exactly 50% => Nx selected at first contig.
    // Source: Miller 2010 "at least 50%"; QUAST N50.py (s <= limit). INV-05.
    [Test]
    public void CalculateNx_CumulativeExactlyAtThreshold_SelectsThatContig()
    {
        var sorted = new List<int> { 50, 50 };
        long total = 100;

        var n50 = GenomeAssemblyAnalyzer.CalculateNx(sorted, total, 50);

        Assert.Multiple(() =>
        {
            Assert.That(n50.Nx, Is.EqualTo(50),
                "INV-05: cumulative 50 equals exactly 50% of 100; the inclusive '>= threshold' boundary selects the first contig.");
            Assert.That(n50.Lx, Is.EqualTo(1),
                "Only one contig is needed to reach exactly 50%.");
        });
    }

    // S1 — Single contig is the whole assembly: N50 = its length, L50 = 1.
    [Test]
    public void CalculateNx_SingleContig_ReturnsThatLengthAndCountOne()
    {
        var n50 = GenomeAssemblyAnalyzer.CalculateNx(new List<int> { 100 }, 100, 50);

        Assert.Multiple(() =>
        {
            Assert.That(n50.Nx, Is.EqualTo(100), "A single contig of length 100 is the entire assembly, so N50 = 100.");
            Assert.That(n50.Lx, Is.EqualTo(1), "L50 = 1 for a single contig.");
        });
    }

    // C2 — Empty input: documented edge — returns Nx=Lx=0 without throwing.
    // Source: Evidence ASSUMPTION 1 (QUAST returns None; repository returns zeros).
    [Test]
    public void CalculateNx_EmptyInput_ReturnsZeros()
    {
        var n50 = GenomeAssemblyAnalyzer.CalculateNx(new List<int>(), 0, 50);

        Assert.Multiple(() =>
        {
            Assert.That(n50.Nx, Is.EqualTo(0), "Empty assembly has no defined N50; the repository returns 0.");
            Assert.That(n50.Lx, Is.EqualTo(0), "Empty assembly returns L50 = 0.");
        });
    }

    #endregion

    #region CalculateNx (lengths, threshold) and CalculateN50 (delegates)

    // M8 — 2-arg overload sorts internally; shuffled input still yields N50=70, L50=2.
    // Source: delegates to the core overload (M1).
    [Test]
    public void CalculateNx_UnsortedLengthsThreshold50_ReturnsN50_70_L50_2()
    {
        var shuffled = new[] { 20, 80, 50, 30, 70, 40 }; // same multiset as Assembly A

        var n50 = GenomeAssemblyAnalyzer.CalculateNx(shuffled, 50);

        Assert.Multiple(() =>
        {
            Assert.That(n50.Nx, Is.EqualTo(70), "The 2-arg overload sorts descending internally, so order does not matter: N50 = 70.");
            Assert.That(n50.Lx, Is.EqualTo(2), "L50 = 2.");
        });
    }

    // M7 — CalculateN50 delegate returns the N50 length for unsorted input.
    // Source: Wikipedia worked example (Assembly A).
    [Test]
    public void CalculateN50_UnsortedAssemblyA_Returns70()
    {
        var shuffled = new[] { 30, 50, 80, 20, 70, 40 };

        int n50 = GenomeAssemblyAnalyzer.CalculateN50(shuffled);

        Assert.That(n50, Is.EqualTo(70),
            "CalculateN50 delegates to CalculateNx(.,50).Nx; Assembly A's N50 is 70 (Wikipedia).");
    }

    // M2 — Assembly B {…,10,5}, total 305: 80+70+50=200 >= 152.5 => N50=50, L50=3.
    // Source: Wikipedia worked example (Assembly B).
    [Test]
    public void CalculateNx_AssemblyB_ReturnsN50_50_L50_3()
    {
        var n50 = GenomeAssemblyAnalyzer.CalculateNx(AssemblyB, 50);

        Assert.Multiple(() =>
        {
            Assert.That(n50.Nx, Is.EqualTo(50),
                "Assembly B (total 305): 80+70+50=200 first reaches >=50% of 305 (=152.5) at the third contig => N50 = 50 (Wikipedia).");
            Assert.That(n50.Lx, Is.EqualTo(3), "L50 = 3 (Wikipedia).");
        });
    }

    #endregion

    #region CalculateAuN

    // M6 — auN = Sum(l^2)/Sum(l) = (100^2+80^2+60^2+40^2+20^2)/300 = 22000/300 = 73.3333...
    // Source: Li H (2020); QUAST au_metric.
    [Test]
    public void CalculateAuN_KnownLengths_EqualsSumOfSquaresOverTotal()
    {
        var lengths = new[] { 100, 80, 60, 40, 20 };

        double auN = GenomeAssemblyAnalyzer.CalculateAuN(lengths);

        Assert.That(auN, Is.EqualTo(22000.0 / 300.0).Within(1e-10),
            "auN = Sum(l^2)/Sum(l) = 22000/300 = 73.3333... (Li 2020; QUAST au_metric).");
    }

    // M6b — auN for Assembly A: (6400+4900+2500+1600+900+400)/290 = 16700/290.
    // Source: Li H (2020); QUAST au_metric.
    [Test]
    public void CalculateAuN_AssemblyA_Equals16700Over290()
    {
        double auN = GenomeAssemblyAnalyzer.CalculateAuN(AssemblyA);

        Assert.That(auN, Is.EqualTo(16700.0 / 290.0).Within(1e-10),
            "auN of Assembly A = (80^2+70^2+50^2+40^2+30^2+20^2)/290 = 16700/290 (Li 2020).");
    }

    // C3 — Empty input returns 0 (documented edge).
    [Test]
    public void CalculateAuN_EmptyInput_ReturnsZero()
    {
        double auN = GenomeAssemblyAnalyzer.CalculateAuN(new List<int>());

        Assert.That(auN, Is.EqualTo(0.0).Within(1e-10), "auN of an empty assembly is 0 (documented edge).");
    }

    #endregion

    #region CalculateStatistics

    // M9 — Full aggregation over Assembly A as contigs: N50/L50/N90/L90/largest/smallest/total.
    // Source: Wikipedia worked example; lengths realised as all-G contigs.
    [Test]
    public void CalculateStatistics_AssemblyA_AggregatesContiguityMetrics()
    {
        var stats = GenomeAssemblyAnalyzer.CalculateStatistics(AsContigs(AssemblyA));

        Assert.Multiple(() =>
        {
            Assert.That(stats.TotalSequences, Is.EqualTo(6), "Six contigs in Assembly A.");
            Assert.That(stats.TotalLength, Is.EqualTo(290), "Total assembly length is 290.");
            Assert.That(stats.N50, Is.EqualTo(70), "N50 = 70 (Wikipedia).");
            Assert.That(stats.L50, Is.EqualTo(2), "L50 = 2 (Wikipedia).");
            Assert.That(stats.N90, Is.EqualTo(30), "N90 = 30 (cumulative 270 >= 261).");
            Assert.That(stats.L90, Is.EqualTo(5), "L90 = 5.");
            Assert.That(stats.LargestContig, Is.EqualTo(80), "Largest contig is 80.");
            Assert.That(stats.SmallestContig, Is.EqualTo(20), "Smallest contig is 20.");
        });
    }

    // S4 — GC content over an all-GC contig is 1.0 exactly.
    [Test]
    public void CalculateStatistics_AllGcContig_GcContentIsOne()
    {
        var stats = GenomeAssemblyAnalyzer.CalculateStatistics(
            new List<(string, string)> { ("seq1", "GCGCGCGC") });

        Assert.That(stats.GcContent, Is.EqualTo(1.0).Within(1e-10),
            "An all-G/C contig has GC fraction 1.0 over its non-N bases.");
    }

    // S5 — Gap aggregation: one 4-N gap in a length-12 contig => GapPercentage = 100*4/12.
    [Test]
    public void CalculateStatistics_WithGap_AggregatesGapStats()
    {
        var stats = GenomeAssemblyAnalyzer.CalculateStatistics(
            new List<(string, string)> { ("seq1", "ACGTNNNNACGT") });

        Assert.Multiple(() =>
        {
            Assert.That(stats.TotalGaps, Is.EqualTo(1), "One N-run gap.");
            Assert.That(stats.TotalGapLength, Is.EqualTo(4), "The gap is 4 bases long.");
            Assert.That(stats.GapPercentage, Is.EqualTo(100.0 * 4 / 12).Within(1e-10),
                "GapPercentage = 100 * gapLength / totalLength = 100*4/12.");
        });
    }

    // C1 — Empty input returns the all-zero statistics record without throwing.
    // Source: Evidence ASSUMPTION 1.
    [Test]
    public void CalculateStatistics_EmptyInput_ReturnsAllZeros()
    {
        var stats = GenomeAssemblyAnalyzer.CalculateStatistics(new List<(string, string)>());

        Assert.Multiple(() =>
        {
            Assert.That(stats.TotalSequences, Is.EqualTo(0), "No sequences.");
            Assert.That(stats.TotalLength, Is.EqualTo(0), "Zero total length.");
            Assert.That(stats.N50, Is.EqualTo(0), "No N50 for an empty assembly.");
            Assert.That(stats.L50, Is.EqualTo(0), "No L50 for an empty assembly.");
        });
    }

    // C4 — All-N contig: whole sequence is a gap (TotalGapLength = length, GapPercentage = 100).
    [Test]
    public void CalculateStatistics_AllNContig_IsEntirelyGap()
    {
        var stats = GenomeAssemblyAnalyzer.CalculateStatistics(
            new List<(string, string)> { ("seq1", new string('N', 100)) });

        Assert.Multiple(() =>
        {
            Assert.That(stats.TotalLength, Is.EqualTo(100), "Length counts all bases including N.");
            Assert.That(stats.TotalGapLength, Is.EqualTo(100), "All 100 bases are gap.");
            Assert.That(stats.GapPercentage, Is.EqualTo(100.0).Within(1e-10), "GapPercentage is 100%.");
        });
    }

    #endregion

    #region FindGaps

    // M10 — Single interior gap "ACGTNNNNACGT": one gap at 0-based inclusive [4,7], length 4.
    // Source: INV-06.
    [Test]
    public void FindGaps_SingleInteriorGap_ReportsExactCoordinates()
    {
        var gaps = GenomeAssemblyAnalyzer.FindGaps(
            new List<(string, string)> { ("seq1", "ACGTNNNNACGT") }).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(gaps, Has.Count.EqualTo(1), "Exactly one maximal N-run.");
            Assert.That(gaps[0].Start, Is.EqualTo(4), "Gap starts at 0-based index 4.");
            Assert.That(gaps[0].End, Is.EqualTo(7), "Gap ends at inclusive index 7.");
            Assert.That(gaps[0].Length, Is.EqualTo(4), "Length = End - Start + 1 = 4 (INV-06).");
        });
    }

    // M11 — Leading gap "NNNNACGT": gap at [0,3], length 4.
    [Test]
    public void FindGaps_LeadingGap_ReportsFromIndexZero()
    {
        var gaps = GenomeAssemblyAnalyzer.FindGaps(
            new List<(string, string)> { ("seq1", "NNNNACGT") }).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(gaps, Has.Count.EqualTo(1), "One leading N-run.");
            Assert.That(gaps[0].Start, Is.EqualTo(0), "Leading gap starts at index 0.");
            Assert.That(gaps[0].End, Is.EqualTo(3), "Leading gap ends at index 3.");
            Assert.That(gaps[0].Length, Is.EqualTo(4), "Length 4.");
        });
    }

    // M12 — Trailing gap "ACGTNNNN": gap at [4,7], length 4.
    [Test]
    public void FindGaps_TrailingGap_ReportsToLastIndex()
    {
        var gaps = GenomeAssemblyAnalyzer.FindGaps(
            new List<(string, string)> { ("seq1", "ACGTNNNN") }).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(gaps, Has.Count.EqualTo(1), "One trailing N-run.");
            Assert.That(gaps[0].Start, Is.EqualTo(4), "Trailing gap starts at index 4.");
            Assert.That(gaps[0].End, Is.EqualTo(7), "Trailing gap ends at the last index 7.");
            Assert.That(gaps[0].Length, Is.EqualTo(4), "Length 4.");
        });
    }

    // M13 — minGapLength filter: "ACGTNNACGTNNNNNNACGT" with minGapLength=5 keeps only the 6-N gap at [10,15].
    // Source: INV-06; documented boundary handling.
    [Test]
    public void FindGaps_MinGapLengthFilter_KeepsOnlyLongGap()
    {
        var gaps = GenomeAssemblyAnalyzer.FindGaps(
            new List<(string, string)> { ("seq1", "ACGTNNACGTNNNNNNACGT") },
            minGapLength: 5).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(gaps, Has.Count.EqualTo(1), "The 2-N gap is below minGapLength=5 and is filtered out.");
            Assert.That(gaps[0].Start, Is.EqualTo(10), "Surviving gap starts at index 10.");
            Assert.That(gaps[0].End, Is.EqualTo(15), "Surviving gap ends at index 15.");
            Assert.That(gaps[0].Length, Is.EqualTo(6), "Surviving gap length is 6.");
        });
    }

    // S2 — No N runs => empty result.
    [Test]
    public void FindGaps_NoGaps_ReturnsEmpty()
    {
        var gaps = GenomeAssemblyAnalyzer.FindGaps(
            new List<(string, string)> { ("seq1", "ACGTACGT") }).ToList();

        Assert.That(gaps, Is.Empty, "A gap-free contig yields no gaps.");
    }

    // S3 — Two separated gaps "ACGTNNNACGTNNNNNNACGT": gap1 [4,6] len3, gap2 [11,16] len6.
    [Test]
    public void FindGaps_TwoSeparatedGaps_ReportsBothExactly()
    {
        var gaps = GenomeAssemblyAnalyzer.FindGaps(
            new List<(string, string)> { ("seq1", "ACGTNNNACGTNNNNNNACGT") }).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(gaps, Has.Count.EqualTo(2), "Two separated N-runs.");
            Assert.That(gaps[0].Start, Is.EqualTo(4), "First gap starts at 4.");
            Assert.That(gaps[0].Length, Is.EqualTo(3), "First gap length 3.");
            Assert.That(gaps[1].Start, Is.EqualTo(11), "Second gap starts at 11.");
            Assert.That(gaps[1].Length, Is.EqualTo(6), "Second gap length 6.");
        });
    }

    #endregion
}
