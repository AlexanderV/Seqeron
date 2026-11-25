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

    // A3-1 — auNG (Li 2020; QUAST au_metric(numlist, reference_length)): Σ l² / G with the genome size G instead of
    // Σ l. Oracle: tools/review/oracles/B10/wp4_oracle.py (QUAST au_metric verbatim + exact Fraction), 7 seeds x 1500
    // random cases (genome sizes 1 .. 3e9, lengths to 2^31) agree within 1e-9 relative with QUAST's float value
    // and exactly with the correctly rounded quotient.
    [TestCase(new[] { 100, 80, 60, 40, 20 }, 1000L, 22.0)]                       // 22000 / 1000
    [TestCase(new[] { 80, 70, 50, 40, 30, 20 }, 300L, 55.666666666666664)]       // 16700 / 300 (Assembly A, G > assembly)
    [TestCase(new[] { 5 }, 100L, 0.25)]
    [TestCase(new[] { 2000000000, 1999999999, 1500000001, 3, 1 }, 3000000000L, 3416666666.3333335)]
    [TestCase(new[] { 2000000000, 1999999999, 1500000001, 3, 1 }, 100L, 1.0249999999e+17)]     // Σl² > 2^63, exact
    [TestCase(new[] { 94906267, 94906267, 1 }, 200000000L, 90071995.15875289)]
    public void CalculateAuNG_GenomeSize_EqualsSumOfSquaresOverGenomeSize(int[] lengths, long genomeSize, double expected)
    {
        Assert.That(GenomeAssemblyAnalyzer.CalculateAuN(lengths, genomeSize), Is.EqualTo(expected).Within(1).Ulps);
    }

    // A3-1 — with G = Σ l the auNG is the auN (QUAST: the denominator is the assembly length when no reference is given).
    [Test]
    public void CalculateAuNG_GenomeSizeEqualToAssemblyLength_IsAuN()
    {
        var lengths = new[] { 100, 80, 60, 40, 20 };

        Assert.That(GenomeAssemblyAnalyzer.CalculateAuN(lengths, 300), Is.EqualTo(GenomeAssemblyAnalyzer.CalculateAuN(lengths)));
    }

    [Test]
    public void CalculateAuNG_EmptyLengths_ReturnsZero_AndInvalidInputsThrow()
    {
        Assert.Multiple(() =>
        {
            Assert.That(GenomeAssemblyAnalyzer.CalculateAuN(new List<int>(), 1000), Is.EqualTo(0.0));
            Assert.Throws<ArgumentNullException>(() => GenomeAssemblyAnalyzer.CalculateAuN(null!, 1000));
            Assert.Throws<ArgumentException>(() => GenomeAssemblyAnalyzer.CalculateAuN(new[] { 3, -1 }, 1000));
            Assert.Throws<ArgumentOutOfRangeException>(() => GenomeAssemblyAnalyzer.CalculateAuN(new[] { 3 }, 0));
            Assert.Throws<ArgumentOutOfRangeException>(() => GenomeAssemblyAnalyzer.CalculateAuN(new[] { 3 }, -5));
        });
    }

    // C3 — Empty input returns 0 (documented edge).
    [Test]
    public void CalculateAuN_EmptyInput_ReturnsZero()
    {
        double auN = GenomeAssemblyAnalyzer.CalculateAuN(new List<int>());

        Assert.That(auN, Is.EqualTo(0.0).Within(1e-10), "auN of an empty assembly is 0 (documented edge).");
    }

    #endregion

    #region QUAST N50.py oracle (B10 ASSEMBLY-STATS-001 review)

    // M14 — Expected values computed by running QUAST quast_libs/N50.py (ablab/quast master,
    // N50_and_L50 / NG50_and_LG50) on each length set; see Evidence "QUAST oracle table".
    [TestCase(new[] { 4, 3, 2 }, 50, 3, 2, TestName = "QuastN50_432_Odd9_N50is3")]
    [TestCase(new[] { 4, 3, 2 }, 90, 2, 3, TestName = "QuastN90_432")]
    [TestCase(new[] { 2, 1, 1, 1 }, 50, 1, 2, TestName = "QuastN50_2111_N50is1")]
    [TestCase(new[] { 2, 1, 1, 1 }, 90, 1, 4, TestName = "QuastN90_2111")]
    [TestCase(new[] { 5, 5, 5, 5 }, 50, 5, 2, TestName = "QuastN50_Ties")]
    [TestCase(new[] { 10, 9, 8, 7, 6, 5, 4, 3, 2, 1 }, 50, 7, 4, TestName = "QuastN50_OddTotal55")]
    [TestCase(new[] { 10, 9, 8, 7, 6, 5, 4, 3, 2, 1 }, 90, 3, 8, TestName = "QuastN90_OddTotal55")]
    [TestCase(new[] { 10, 9, 8, 7, 6, 5, 4, 3, 2, 1 }, 0, 10, 1, TestName = "QuastN0_IsLargest")]
    [TestCase(new[] { 10, 9, 8, 7, 6, 5, 4, 3, 2, 1 }, 100, 1, 10, TestName = "QuastN100_IsSmallest")]
    [TestCase(new[] { 80, 70, 50, 40, 30, 20, 10, 5 }, 90, 20, 6, TestName = "QuastN90_AssemblyB")]
    [TestCase(new[] { 80, 70, 50, 40, 30, 20, 10, 5 }, 100, 5, 8, TestName = "QuastN100_AssemblyB")]
    [TestCase(new[] { 3, 0, 0 }, 100, 3, 1, TestName = "QuastN100_TrailingZeroLengths")]
    [TestCase(new[] { 0, 0 }, 50, 0, 1, TestName = "QuastN50_AllZeroLengths")]
    public void CalculateNx_MatchesQuastN50Py(int[] lengths, int x, int expectedNx, int expectedLx)
    {
        var nx = GenomeAssemblyAnalyzer.CalculateNx(lengths, x);

        Assert.Multiple(() =>
        {
            Assert.That(nx.Nx, Is.EqualTo(expectedNx), "Nx = QUAST N50_and_L50(lengths, x)[0].");
            Assert.That(nx.Lx, Is.EqualTo(expectedLx), "Lx = QUAST N50_and_L50(lengths, x)[1].");
        });
    }

    // M15 — NGx/LGx: the core overload with a genome size (QUAST NG50_and_LG50). Assembly A with
    // reference 400: NG50 = 50, LG50 = 3; NG90 is never reached (QUAST None) => 0/0.
    // {4,3,2} with reference 7: NG50 = 4/1, NG90 = 3/2.
    [TestCase(new[] { 80, 70, 50, 40, 30, 20 }, 400L, 50, 50, 3)]
    [TestCase(new[] { 80, 70, 50, 40, 30, 20 }, 400L, 90, 0, 0)]
    [TestCase(new[] { 80, 70, 50, 40, 30, 20 }, 1000L, 50, 0, 0)]
    [TestCase(new[] { 4, 3, 2 }, 7L, 50, 4, 1)]
    [TestCase(new[] { 4, 3, 2 }, 7L, 90, 3, 2)]
    public void CalculateNx_GenomeSizeReference_MatchesQuastNG50(
        int[] lengthsDescending, long reference, int x, int expectedNg, int expectedLg)
    {
        var ng = GenomeAssemblyAnalyzer.CalculateNx(lengthsDescending, reference, x);

        Assert.Multiple(() =>
        {
            Assert.That(ng.Nx, Is.EqualTo(expectedNg), "NGx = QUAST NG50_and_LG50 (0 where QUAST gives None).");
            Assert.That(ng.Lx, Is.EqualTo(expectedLg), "LGx = QUAST NG50_and_LG50 (0 where QUAST gives None).");
        });
    }

    // M16 — Totals above int.MaxValue are summed exactly (QUAST: {2e9,2e9,1} -> N50 2e9/L50 2, N100 1/3).
    [Test]
    public void CalculateNx_TotalAboveInt32_IsExact()
    {
        int[] lengths = { 2_000_000_000, 2_000_000_000, 1 };

        var n50 = GenomeAssemblyAnalyzer.CalculateNx(lengths, 50);
        var n100 = GenomeAssemblyAnalyzer.CalculateNx(lengths, 100);

        Assert.Multiple(() =>
        {
            Assert.That((n50.Nx, n50.Lx), Is.EqualTo((2_000_000_000, 2)));
            Assert.That((n100.Nx, n100.Lx, n100.CumulativeLength), Is.EqualTo((1, 3, 4_000_000_001L)));
        });
    }

    // M17 — Out-of-contract input is rejected: QUAST asserts 0 <= x <= 100; the core overload needs
    // non-negative lengths in non-increasing order (an unsorted list would give a wrong Nx silently).
    [Test]
    public void CalculateNx_InvalidInput_Throws()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => GenomeAssemblyAnalyzer.CalculateNx(new[] { 3, 2 }, 5, -1));
            Assert.Throws<ArgumentOutOfRangeException>(() => GenomeAssemblyAnalyzer.CalculateNx(new[] { 3, 2 }, 5, 101));
            Assert.Throws<ArgumentOutOfRangeException>(() => GenomeAssemblyAnalyzer.CalculateNx(new[] { 3, 2 }, -5, 50));
            Assert.Throws<ArgumentException>(() => GenomeAssemblyAnalyzer.CalculateNx(new[] { 2, 3 }, 5, 50));
            Assert.Throws<ArgumentException>(() => GenomeAssemblyAnalyzer.CalculateNx(new[] { 3, -1 }, 2, 50));
            // A3-3: the curve validates eagerly (no enumeration needed to observe the error).
            Assert.Throws<ArgumentOutOfRangeException>(() => GenomeAssemblyAnalyzer.CalculateNxCurve(new[] { 3, 2 }, 150));
        });
    }

    // A3-3 — null / invalid inputs are rejected eagerly with ArgumentNullException / ArgumentException
    // (no NullReferenceException, no lazy iterator failure).
    [Test]
    public void NullAndInvalidInputs_ThrowEagerly()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentNullException>(() => GenomeAssemblyAnalyzer.CalculateStatistics(null!));
            Assert.Throws<ArgumentException>(() => GenomeAssemblyAnalyzer.CalculateStatistics(
                new[] { ("a", "ACGT"), ("b", (string)null!) }));
            Assert.Throws<ArgumentNullException>(() => GenomeAssemblyAnalyzer.CalculateNx((IEnumerable<int>)null!, 50));
            Assert.Throws<ArgumentNullException>(() => GenomeAssemblyAnalyzer.CalculateNx((IReadOnlyList<int>)null!, 5, 50));
            Assert.Throws<ArgumentNullException>(() => GenomeAssemblyAnalyzer.CalculateN50(null!));
            Assert.Throws<ArgumentNullException>(() => GenomeAssemblyAnalyzer.CalculateAuN(null!));
            Assert.Throws<ArgumentException>(() => GenomeAssemblyAnalyzer.CalculateAuN(new[] { 3, -1 }));
            Assert.Throws<ArgumentNullException>(() => GenomeAssemblyAnalyzer.CalculateLengthDistribution(null!));
            Assert.Throws<ArgumentNullException>(() => GenomeAssemblyAnalyzer.CalculateLengthDistribution(new[] { 1 }, null!));
        });
    }

    // A3-4 — auN accumulates Σl² and Σl exactly (QUAST au_metric uses Python integers). With
    // lengths whose squares exceed 2^53 a double accumulation loses the low bits; oracle values
    // from tools/review/oracles/B10/stats_oracle.py (correctly rounded Fraction(Σl², Σl)).
    [TestCase(new[] { 2000000000, 1999999999, 1500000001, 3, 1 }, 1863636362.0991735)]   // Σ=5 500 000 004, Σl²≈1.0e19 > 2^63
    [TestCase(new[] { 94906267, 94906267, 1 }, 94906266.50000001)]                       // Σl² = 18 014 399 031 750 579 > 2^53
    public void CalculateAuN_BeyondDoublePrecision_MatchesExactQuastValue(int[] lengths, double expected)
    {
        Assert.That(GenomeAssemblyAnalyzer.CalculateAuN(lengths), Is.EqualTo(expected).Within(1).Ulps);
    }

    // A3-4 — 50 × int.MaxValue + 1: Σl² > 2^63 (would overflow a 64-bit sum); oracle 2147483646.98.
    [Test]
    public void CalculateAuN_SumOfSquaresAbove2Pow63_NoOverflow()
    {
        var lengths = Enumerable.Repeat(int.MaxValue, 50).Append(1).ToArray();
        Assert.That(GenomeAssemblyAnalyzer.CalculateAuN(lengths), Is.EqualTo(2147483646.98).Within(1).Ulps);
    }

    #region CalculateNxCurve (A3-6)

    // Oracle (QUAST NG50_and_LG50 with reference = Σl): Assembly B {80,70,50,40,30,20,10,5}.
    private static readonly (int T, int Nx, int Lx)[] AssemblyBDecileCurve =
    {
        (10, 80, 1), (20, 80, 1), (30, 70, 2), (40, 70, 2), (50, 50, 3), (60, 50, 3), (70, 40, 4), (80, 30, 5), (90, 20, 6),
    };

    [Test]
    public void CalculateNxCurve_DefaultThresholds_AreDecilesInAscendingOrder_MatchQuastOracle()
    {
        var curve = GenomeAssemblyAnalyzer.CalculateNxCurve(AssemblyB.Reverse()).ToList();

        Assert.That(curve.Select(c => (c.Threshold, c.Nx, c.Lx)), Is.EqualTo(AssemblyBDecileCurve),
            "default thresholds 10..90 in ascending order, each point = QUAST NG50_and_LG50 (input order irrelevant)");
    }

    [Test]
    public void CalculateNxCurve_Oracle_500_400_300_200_100()
    {
        var curve = GenomeAssemblyAnalyzer.CalculateNxCurve(new[] { 100, 300, 500, 200, 400 }).ToList();

        Assert.That(curve.Select(c => (c.Threshold, c.Nx, c.Lx)), Is.EqualTo(new[]
        {
            (10, 500, 1), (20, 500, 1), (30, 500, 1), (40, 400, 2), (50, 400, 2), (60, 400, 2),
            (70, 300, 3), (80, 300, 3), (90, 200, 4),
        }));
    }

    [Test]
    public void CalculateNxCurve_ExplicitThresholds_SortedAscending_DuplicatesKept()
    {
        var curve = GenomeAssemblyAnalyzer.CalculateNxCurve(AssemblyB, 90, 50, 50, 0, 100).ToList();

        Assert.That(curve.Select(c => c.Threshold), Is.EqualTo(new[] { 0, 50, 50, 90, 100 }));
        Assert.That(curve[1], Is.EqualTo(curve[2]), "a repeated threshold yields an identical repeated entry");
        Assert.That(curve[1].Nx, Is.EqualTo(50));
        Assert.That(curve[0].Nx, Is.EqualTo(80), "N0: the first contig already reaches 0% (QUAST: s <= limit at the first length)");
        Assert.That(curve[4].Nx, Is.EqualTo(5), "N100 is the shortest contig (all 305 bases needed)");
        Assert.That(curve[4].Lx, Is.EqualTo(8));
    }

    [Test]
    public void CalculateNxCurve_PointsEqualCalculateNx()
    {
        foreach (var point in GenomeAssemblyAnalyzer.CalculateNxCurve(AssemblyA, 25, 50, 75))
            Assert.That(point, Is.EqualTo(GenomeAssemblyAnalyzer.CalculateNx(AssemblyA, point.Threshold)));
    }

    [Test]
    public void CalculateNxCurve_EmptyLengths_ZerosAtEveryThreshold()
    {
        var curve = GenomeAssemblyAnalyzer.CalculateNxCurve(Array.Empty<int>()).ToList();
        Assert.That(curve, Has.Count.EqualTo(9));
        Assert.That(curve.All(c => c.Nx == 0 && c.Lx == 0));
    }

    [Test]
    public void CalculateNxCurve_NullAndInvalid_ThrowAtTheCall()
    {
        // No .ToList()/enumeration: validation is eager (A3-3).
        Assert.Throws<ArgumentNullException>(() => GenomeAssemblyAnalyzer.CalculateNxCurve(null!));
        Assert.Throws<ArgumentNullException>(() => GenomeAssemblyAnalyzer.CalculateNxCurve(AssemblyA, null!));
        Assert.Throws<ArgumentOutOfRangeException>(() => GenomeAssemblyAnalyzer.CalculateNxCurve(AssemblyA, 50, -1));
        Assert.Throws<ArgumentOutOfRangeException>(() => GenomeAssemblyAnalyzer.CalculateNxCurve(AssemblyA, 101));
        Assert.Throws<ArgumentException>(() => GenomeAssemblyAnalyzer.CalculateNxCurve(new[] { 3, -2 }));
    }

    #endregion

    // M18 — auN = QUAST au_metric: Assembly B 16825/305 = 55.16393442622951; {4,3,2} 29/9.
    [TestCase(new[] { 80, 70, 50, 40, 30, 20, 10, 5 }, 55.16393442622951)]
    [TestCase(new[] { 4, 3, 2 }, 3.2222222222222223)]
    [TestCase(new[] { 2, 1, 1, 1 }, 1.4)]
    public void CalculateAuN_MatchesQuastAuMetric(int[] lengths, double expected)
    {
        Assert.That(GenomeAssemblyAnalyzer.CalculateAuN(lengths), Is.EqualTo(expected).Within(1e-12));
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

    // M19 — GC% is QUAST GC_content: pooled (G+C)/(A+C+G+T); correct_seq turns IUPAC codes into N
    // and upper-cases, so ambiguity codes and N are in neither count. QUAST values:
    // ["ACGTRYSW"] -> 50.0 %; ["GGCCnnAT","acgtNNNN"] -> 60.0 %; ["GCGCRRRR"] -> 100.0 %.
    [TestCase(new[] { "ACGTRYSW" }, 0.5)]
    [TestCase(new[] { "GGCCnnAT", "acgtNNNN" }, 0.6)]
    [TestCase(new[] { "GCGCRRRR" }, 1.0)]
    [TestCase(new[] { "ATATATAT", "GGGG" }, 4.0 / 12.0)]
    public void CalculateStatistics_GcContent_MatchesQuastGcContent(string[] sequences, double expected)
    {
        var stats = GenomeAssemblyAnalyzer.CalculateStatistics(
            sequences.Select((q, i) => ($"c{i}", q)));

        Assert.That(stats.GcContent, Is.EqualTo(expected).Within(1e-12),
            "GC fraction = (G+C)/(A+C+G+T), IUPAC codes and N excluded (QUAST GC_content).");
    }

    // A3-2 — QUAST correct_seq reads every non-ACGT symbol as N: IUPAC codes (and, leniently, U, '-', ...) are gap
    // symbols, not bases, in the gap fields and in the GC denominator (len - count('N')). nonAcgtAsGap = false keeps
    // the previous reading (only N/n is a gap; U/IUPAC are bases of the length, outside the A/C/G/T denominator).
    // Oracle (wp4_oracle.stats_expected, QUAST correct_seq + split_by_ns + GC_content formula):
    // "ACGU" -> GC 2/3, gaps 1 of length 1; "ACRYGT" -> gaps 1 of length 2, total without gaps 4.
    [Test]
    public void CalculateStatistics_NonAcgtAsGap_FollowsQuastCorrectSeq()
    {
        var u = GenomeAssemblyAnalyzer.CalculateStatistics(new[] { ("a", "ACGU") });
        var uBase = GenomeAssemblyAnalyzer.CalculateStatistics(new[] { ("a", "ACGU") }, nonAcgtAsGap: false);
        var iupac = GenomeAssemblyAnalyzer.CalculateStatistics(new[] { ("a", "ACRYGT") });
        var iupacBase = GenomeAssemblyAnalyzer.CalculateStatistics(new[] { ("a", "ACRYGT") }, nonAcgtAsGap: false);

        Assert.Multiple(() =>
        {
            Assert.That(u.GcContent, Is.EqualTo(2.0 / 3.0).Within(1e-12), "U is an N in QUAST: 2 GC of 3 ACGT");
            Assert.That((u.TotalGaps, u.TotalGapLength), Is.EqualTo((1, 1L)));
            Assert.That(uBase.GcContent, Is.EqualTo(0.5).Within(1e-12), "previous reading: U counted as a valid base");
            Assert.That(uBase.TotalGaps, Is.EqualTo(0));
            Assert.That((iupac.TotalGaps, iupac.TotalGapLength, iupac.TotalLengthNoGaps), Is.EqualTo((1, 2L, 4L)));
            Assert.That(iupac.GapPercentage, Is.EqualTo(100.0 * 2 / 6).Within(1e-12));
            Assert.That((iupacBase.TotalGaps, iupacBase.TotalGapLength), Is.EqualTo((0, 0L)));
            Assert.That(iupac.GcContent, Is.EqualTo(0.5).Within(1e-12), "(G + C) / (A + C + G + T) = 2 / 4 in both readings");
            Assert.That(iupacBase.GcContent, Is.EqualTo(0.5).Within(1e-12));
        });
    }

    // M20 — Median contig length is the sample median: {4,3,2,1} -> 2.5 (numpy/R median), odd {5,3,1} -> 3.
    [TestCase(new[] { 4, 3, 2, 1 }, 2.5)]
    [TestCase(new[] { 5, 3, 1 }, 3.0)]
    public void CalculateStatistics_MedianLength_IsSampleMedian(int[] lengths, double expected)
    {
        var stats = GenomeAssemblyAnalyzer.CalculateStatistics(AsContigs(lengths));

        Assert.That(stats.MedianLength, Is.EqualTo(expected));
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

    // C5 — Length histogram bins are half-open [previous, bound); bounds must be strictly increasing.
    [Test]
    public void CalculateLengthDistribution_BinsAndValidation()
    {
        var dist = GenomeAssemblyAnalyzer.CalculateLengthDistribution(new[] { 50, 100, 499, 500, 2000 }, 100, 500);

        Assert.Multiple(() =>
        {
            Assert.That(dist["<100"], Is.EqualTo(1));
            Assert.That(dist["<500"], Is.EqualTo(2));
            Assert.That(dist[">=500"], Is.EqualTo(2));
            Assert.Throws<ArgumentException>(() => GenomeAssemblyAnalyzer.CalculateLengthDistribution(new[] { 50 }, 500, 100));
            Assert.Throws<ArgumentException>(() => GenomeAssemblyAnalyzer.CalculateLengthDistribution(new[] { 50 }, 100, 100));
        });
    }
}
