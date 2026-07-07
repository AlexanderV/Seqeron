using static Seqeron.Genomics.Oncology.OncologyAnalyzer;

namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the Oncology INTRATUMOUR-HETEROGENEITY scoring unit — ONCO-HETERO-001.
/// The unit under test is the deterministic statistic suite implemented in
/// src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs:
/// <see cref="OncologyAnalyzer.CalculateITH(IReadOnlyList{double})"/> — the MATH score over a VAF
/// distribution — and the aggregate
/// <see cref="OncologyAnalyzer.AnalyzeHeterogeneity(IReadOnlyList{double}, IReadOnlyList{double}, int)"/>
/// (MATH + Shannon diversity + subclone count + subclonal fraction).
///
/// This file is scoped to the MATH/heterogeneity SCORING contract. CCF clustering itself is ONCO-CCF-001
/// (OncologyCcfFuzzTests.cs) and clonal/subclonal CLASSIFICATION of a single CCF is ONCO-CLONAL-001
/// (OncologyClonalityFuzzTests.cs); both are out of scope here.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate / boundary / extreme inputs to a unit and asserts the code NEVER fails in an
/// undisciplined way: no hang, no nonsense output, no *unhandled* runtime fault. Every input must resolve to
/// EITHER a well-defined, theory-correct value OR a *documented, intentional* outcome — here, an
/// ArgumentException for an empty list, a non-finite / out-of-[0,1] value, or a ZERO MEDIAN (the MATH ratio
/// divides by the median and is undefined). The headline contract for this BE row is the median-of-zero
/// hazard: with every VAF = 0 the median denominator is 0, so a naïve `100·MAD/median` would leak +∞ / NaN;
/// the contract REJECTS it with ArgumentException instead. No Inf or NaN ever leaks from CalculateITH /
/// AnalyzeHeterogeneity for any accepted input.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: ONCO-HETERO-001 — intratumour heterogeneity scoring (Oncology)
/// Checklist: docs/checklists/03_FUZZING.md, row 116.
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — граничні значення: 0, -1, MaxInt, empty.
///     Targets (checklist row 116): "single VAF, all-equal VAFs, VAF=0".
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// The targets map onto the documented contract (Tumor_Heterogeneity_Analysis.md) as follows:
///   • single VAF     ⇔ one mutation ⇒ median = that value, every abs-deviation = 0 ⇒ MAD = 0 ⇒ MATH = 0.
///                       No variance-of-one crash, no NaN; a singleton is the trivially homogeneous case
///                       (§6.1 "Single VAF", INV-02).
///   • all-equal VAFs ⇔ every fᵢ = the median ⇒ abs-deviations all 0 ⇒ MAD = 0 ⇒ MATH = exactly 0. The
///                       canonical homogeneous tumour (§6.1 "All VAFs identical", INV-02:
///                       MATH = 0 ⇔ MAD = 0).
///   • VAF = 0        ⇔ all-zero VAFs ⇒ median = 0 ⇒ the MATH `100·MAD/median` denominator is 0. This is a
///                       DivideByZero / +∞ / NaN hazard. The contract REJECTS a zero median with
///                       ArgumentException (§3.3, §6.1 "Median VAF = 0"); NO Inf/NaN leaks. A SINGLE non-zero
///                       VAF among zeros (median still 0 when zeros dominate) is exercised too.
///
/// ───────────────────────────────────────────────────────────────────────────
/// The documented contract under test (Tumor_Heterogeneity_Analysis.md)
/// ───────────────────────────────────────────────────────────────────────────
///   • MATH = 100 · MAD / median(f),  MAD = 1.4826 · median(|fᵢ − median(f)|)      (§2.2 Core Model, §4.1)
///   • MAD consistency constant 1.4826 = 1/Φ⁻¹(3/4); percent scale 100              (§4.2)
///   • Even-count median = mean of the two central order statistics (R convention)  (§4.2, §5.4)
///   • INV-01: MATH ≥ 0 (MAD ≥ 0, median > 0 enforced)                              (§2.4)
///   • INV-02: MATH = 0 ⇔ MAD = 0 (all VAFs equal the median)                       (§2.4)
///   • Single VAF ⇒ MATH = 0; all-equal ⇒ MATH = 0                                  (§6.1)
///   • Median = 0 ⇒ ArgumentException (division by zero)                            (§3.3, §6.1)
///   • empty / non-finite / out-of-[0,1] ⇒ ArgumentException                        (§3.3)
///   • null list ⇒ ArgumentNullException                                            (§3.3)
///   • Worked example: CalculateITH({0.10,0.20,0.30,0.40,0.50}) == 49.42            (§7.1)
///   • AnalyzeHeterogeneity: ShannonDiversity ≥ 0; SubcloneCount ∈ [1, k];
///     SubclonalFraction ∈ [0, 1] (CCF < 0.95)                                      (§2.4 INV-03..05, §3.2)
///   • Inputs are not mutated (median helper clones before sorting)                 (§5.2)
///
/// No source bug was found; no test was weakened. The median-of-zero guard already throws ArgumentException.
/// All randomness is LOCALLY seeded (new Random(seed)); no shared static Rng.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public sealed class OncologyHeterogeneityFuzzTests
{
    private const double Tolerance = 1e-9;
    private const double MadConsistencyConstant = 1.4826; // 1/Φ⁻¹(3/4); §4.2
    private const double MathPercentScale = 100.0;        // §4.2

    // ── Well-formed-result assertion helpers ─────────────────────────────────
    // Pin the documented numeric contract on EVERY accepted result so a fuzz
    // test cannot rubber-stamp a malformed value (NaN/Inf MATH, negative MATH,
    // negative Shannon, out-of-range subclone count / subclonal fraction) green.

    private static void AssertWellFormedMath(double math)
    {
        double.IsNaN(math).Should().BeFalse("MATH must never be NaN (zero-median is rejected, not divided)");
        double.IsInfinity(math).Should().BeFalse("MATH must never be ±∞ (the median denominator is guarded > 0)");
        math.Should().BeGreaterThanOrEqualTo(0.0, "MATH ≥ 0: MAD ≥ 0 and median > 0 (INV-01)");
    }

    private static void AssertWellFormedResult(HeterogeneityResult r, int n, int clusterCount)
    {
        AssertWellFormedMath(r.MathScore);

        double.IsNaN(r.ShannonDiversity).Should().BeFalse("Shannon H must never be NaN");
        double.IsInfinity(r.ShannonDiversity).Should().BeFalse("Shannon H must never be ±∞");
        r.ShannonDiversity.Should().BeGreaterThanOrEqualTo(-Tolerance, "H = −Σ p ln p ≥ 0 (INV-03)");

        r.SubcloneCount.Should().BeInRange(1, clusterCount,
            "occupied clusters ∈ [1, k] (INV-05)");

        r.SubclonalFraction.Should().BeInRange(0.0, 1.0,
            "subclonal fraction = #(CCF < 0.95)/n ∈ [0, 1] (INV-05)");
    }

    // Independent reference implementation of the documented MATH closed form,
    // used to cross-check CalculateITH on random in-contract inputs (mutation
    // testing of the formula: any drift in the 1.4826 / 100 constants, the MAD
    // definition, or the even-count median rule is caught).
    private static double ReferenceMath(IReadOnlyList<double> f)
    {
        double median = ReferenceMedian(f);
        double[] dev = f.Select(x => Math.Abs(x - median)).ToArray();
        double mad = ReferenceMedian(dev);
        return MathPercentScale * MadConsistencyConstant * mad / median;
    }

    private static double ReferenceMedian(IReadOnlyList<double> values)
    {
        double[] s = values.ToArray();
        Array.Sort(s);
        int n = s.Length;
        int mid = n / 2;
        return (n % 2 == 1) ? s[mid] : (s[mid - 1] + s[mid]) / 2.0;
    }

    #region ONCO-HETERO-001 — positive sanity (worked example, tight vs wide spread, formula cross-check)

    [Test]
    [CancelAfter(5_000)]
    public void CalculateITH_DocWorkedExample_Yields4942()
    {
        // Tumor_Heterogeneity_Analysis.md §7.1 (odd count):
        //   {0.10,0.20,0.30,0.40,0.50}: median = 0.30; abs-dev {0.20,0.10,0,0.10,0.20};
        //   raw MAD = 0.10; scaled MAD = 1.4826·0.10 = 0.14826; MATH = 100·0.14826/0.30 = 49.42.
        double math = CalculateITH(new[] { 0.10, 0.20, 0.30, 0.40, 0.50 });

        math.Should().BeApproximately(49.42, 1e-2, "hand-computed worked example (§7.1)");
        math.Should().BeApproximately(100.0 * 1.4826 * 0.10 / 0.30, Tolerance,
            "exact = 100·1.4826·rawMAD/median (§2.2)");
        AssertWellFormedMath(math);
    }

    [Test]
    [CancelAfter(5_000)]
    public void CalculateITH_WideSpread_ScoresHigherThanTightCluster()
    {
        // Heterogeneity monotonicity: a tightly clustered VAF set (homogeneous tumour) must score lower
        // than a widely dispersed one (heterogeneous tumour), both with the same central tendency (§2.1).
        double tight = CalculateITH(new[] { 0.48, 0.49, 0.50, 0.51, 0.52 });
        double wide = CalculateITH(new[] { 0.10, 0.30, 0.50, 0.70, 0.90 });

        tight.Should().BeGreaterThanOrEqualTo(0.0);
        wide.Should().BeGreaterThan(tight,
            "a wider, more dispersed VAF distribution is more heterogeneous ⇒ higher MATH (§2.1)");
        AssertWellFormedMath(tight);
        AssertWellFormedMath(wide);
    }

    [Test]
    [CancelAfter(15_000)]
    public void CalculateITH_RandomInContractInputs_MatchReferenceFormulaExactly()
    {
        // Mutation testing of the formula over random spread sets with a strictly positive median.
        var rng = new Random(116_001);
        for (int i = 0; i < 5_000; i++)
        {
            int n = rng.Next(1, 30);
            var f = new double[n];
            for (int k = 0; k < n; k++)
            {
                // (0, 1]: floor at a small positive value so the median is guaranteed > 0.
                f[k] = rng.NextDouble() * 0.9 + 0.05;
            }

            double math = CalculateITH(f);
            double expected = ReferenceMath(f);

            math.Should().BeApproximately(expected, 1e-9,
                "MATH must equal 100·1.4826·median(|f−median|)/median (§2.2)");
            AssertWellFormedMath(math);
        }
    }

    #endregion

    #region ONCO-HETERO-001 — BE: single VAF (no variance-of-one crash; MATH = 0)

    [Test]
    [CancelAfter(5_000)]
    public void CalculateITH_SingleVaf_YieldsExactlyZero_NoCrash()
    {
        // §6.1 "Single VAF": median = the value, every abs-deviation = 0 ⇒ MAD = 0 ⇒ MATH = 0.
        // A singleton must NOT raise a variance-of-one / divide-by-(n−1) fault.
        double math = CalculateITH(new[] { 0.37 });

        math.Should().Be(0.0, "a single mutation is trivially homogeneous ⇒ MATH = 0 (§6.1, INV-02)");
        AssertWellFormedMath(math);
    }

    [Test]
    [CancelAfter(5_000)]
    public void CalculateITH_SingleVaf_AcrossManyValues_AlwaysZero()
    {
        // Property: for ANY single positive VAF the MATH score is exactly 0 (centre = the value, spread 0).
        var rng = new Random(116_002);
        for (int i = 0; i < 1_000; i++)
        {
            double v = rng.NextDouble() * 0.999 + 0.001; // (0, 1]
            double math = CalculateITH(new[] { v });
            math.Should().Be(0.0, "singleton ⇒ MAD = 0 ⇒ MATH = 0 regardless of the value (§6.1)");
        }
    }

    #endregion

    #region ONCO-HETERO-001 — BE: all-equal VAFs (canonical homogeneous case; MATH = 0)

    [Test]
    [CancelAfter(5_000)]
    public void CalculateITH_AllEqualVafs_YieldsExactlyZero()
    {
        // §6.1 "All VAFs identical" / INV-02: every fᵢ = median ⇒ MAD = 0 ⇒ MATH = exactly 0.
        double math = CalculateITH(new[] { 0.42, 0.42, 0.42, 0.42, 0.42, 0.42 });

        math.Should().Be(0.0, "identical VAFs ⇒ zero dispersion ⇒ MATH = 0 (INV-02)");
        AssertWellFormedMath(math);
    }

    [Test]
    [CancelAfter(10_000)]
    public void CalculateITH_AllEqualVafs_AnyValueAnyCount_AlwaysZero()
    {
        // Property: a constant VAF vector of any length / value (positive) has MATH = 0 exactly.
        // Covers BOTH even and odd counts (the even-count median is the mean of the two central values,
        // which still equals the common value, so MAD stays 0).
        var rng = new Random(116_003);
        for (int i = 0; i < 1_000; i++)
        {
            double v = rng.NextDouble() * 0.999 + 0.001; // (0, 1]
            int n = rng.Next(2, 25);
            var f = Enumerable.Repeat(v, n).ToArray();

            double math = CalculateITH(f);
            math.Should().Be(0.0,
                "MATH = 0 ⇔ MAD = 0 (all VAFs equal the median) for any n, any value (INV-02)");
        }
    }

    #endregion

    #region ONCO-HETERO-001 — BE: VAF = 0 / zero median (headline DivideByZero hazard — guarded, no Inf/NaN)

    [Test]
    [CancelAfter(5_000)]
    public void CalculateITH_AllZeroVafs_ThrowsArgumentException_NoInfinityLeak()
    {
        // §6.1 "Median VAF = 0" / §3.3: every VAF = 0 ⇒ median = 0 ⇒ the 100·MAD/median ratio divides by 0.
        // The contract REJECTS it with ArgumentException rather than leaking +∞ / NaN (the headline BE case).
        Action act = () => CalculateITH(new[] { 0.0, 0.0, 0.0, 0.0 });

        act.Should().Throw<ArgumentException>(
                "a zero median makes the MATH ratio undefined (division by zero) (§3.3, §6.1)")
            .Which.ParamName.Should().Be("ccfDistribution");
    }

    [Test]
    [CancelAfter(5_000)]
    public void CalculateITH_SingleZeroVaf_ThrowsArgumentException()
    {
        // A single VAF = 0: median = 0 ⇒ rejected (the singleton MAD-is-0 path must NOT mask the zero median
        // and return 0; the documented behaviour is the ArgumentException guard).
        Action act = () => CalculateITH(new[] { 0.0 });

        act.Should().Throw<ArgumentException>("median = 0 ⇒ undefined MATH ratio (§3.3, §6.1)");
    }

    [Test]
    [CancelAfter(5_000)]
    public void CalculateITH_ZeroDominatedMedianZero_ThrowsArgumentException()
    {
        // Zeros dominate so the median is still 0 even with a few non-zero VAFs present
        // ({0,0,0,0.4,0.5}: sorted central value = 0). The zero-median guard still fires.
        Action act = () => CalculateITH(new[] { 0.0, 0.0, 0.0, 0.4, 0.5 });

        act.Should().Throw<ArgumentException>("median order-statistic is 0 ⇒ undefined ratio (§3.3)")
            .Which.ParamName.Should().Be("ccfDistribution");
    }

    [Test]
    [CancelAfter(10_000)]
    public void CalculateITH_ZeroMedian_RandomZeroMajority_AlwaysThrows_NeverInfNaN()
    {
        // Property over random zero-majority vectors (more than half are 0 ⇒ the central order statistic is 0,
        // or, for even counts where both central statistics are 0, the mean is 0): ALWAYS ArgumentException,
        // NEVER a leaked Inf/NaN/finite-nonsense result.
        var rng = new Random(116_004);
        for (int i = 0; i < 1_000; i++)
        {
            int n = rng.Next(2, 21);
            int zeros = n / 2 + 1 + rng.Next(0, n / 2 + 1); // strictly > n/2
            zeros = Math.Min(zeros, n);
            var f = new double[n];
            for (int k = zeros; k < n; k++)
            {
                f[k] = rng.NextDouble() * 0.9 + 0.1; // remaining non-zero in [0.1, 1.0)
            }

            // Shuffle so zeros are not contiguous (median must still be 0 by majority).
            for (int k = n - 1; k > 0; k--)
            {
                int j = rng.Next(k + 1);
                (f[k], f[j]) = (f[j], f[k]);
            }

            double captured = double.NaN;
            Action act = () => captured = CalculateITH(f);

            act.Should().Throw<ArgumentException>(
                "a zero majority forces median = 0 ⇒ ArgumentException, never an Inf/NaN leak (§3.3)");
            _ = captured; // never assigned a finite nonsense value
        }
    }

    #endregion

    #region ONCO-HETERO-001 — BE: invalid inputs (empty / null / non-finite / out-of-range)

    [Test]
    public void CalculateITH_Null_ThrowsArgumentNull()
    {
        Action act = () => CalculateITH(null!);
        act.Should().Throw<ArgumentNullException>("null list is rejected (§3.3)");
    }

    [Test]
    public void CalculateITH_Empty_ThrowsArgumentException()
    {
        Action act = () => CalculateITH(Array.Empty<double>());
        act.Should().Throw<ArgumentException>("at least one allele fraction is required (§3.3)")
            .Which.ParamName.Should().Be("ccfDistribution");
    }

    [TestCase(double.NaN)]
    [TestCase(double.PositiveInfinity)]
    [TestCase(double.NegativeInfinity)]
    [TestCase(-0.01)]
    [TestCase(1.01)]
    public void CalculateITH_OutOfRangeOrNonFinite_ThrowsArgumentException(double bad)
    {
        // Non-finite or out-of-[0,1] values are rejected BEFORE any division (§3.3); a NaN/Inf input must
        // never flow into the median/MAD arithmetic and leak out.
        Action act = () => CalculateITH(new[] { 0.3, bad, 0.5 });
        act.Should().Throw<ArgumentException>("VAF must be finite ∈ [0, 1] (§3.3)")
            .Which.ParamName.Should().Be("ccfDistribution");
    }

    [Test]
    [CancelAfter(5_000)]
    public void CalculateITH_DoesNotMutateInput()
    {
        // §5.2: the median helper clones before sorting; the caller's array order/contents are preserved.
        var input = new[] { 0.5, 0.1, 0.9, 0.3 };
        var copy = (double[])input.Clone();

        _ = CalculateITH(input);

        input.Should().Equal(copy, "CalculateITH must not mutate the caller's array (§5.2)");
    }

    #endregion

    #region ONCO-HETERO-001 — BE: AnalyzeHeterogeneity aggregate (Shannon / subclones, same VAF guards)

    [Test]
    [CancelAfter(5_000)]
    public void AnalyzeHeterogeneity_SingleMutation_HomogeneousAllZeroExceptMath()
    {
        // One mutation (single VAF, single CCF, k = 1): MATH = 0 (singleton), one occupied clone ⇒ H = 0,
        // subclone count = 1; CCF = 1.0 (clonal) ⇒ subclonal fraction = 0.
        var r = AnalyzeHeterogeneity(new[] { 0.5 }, new[] { 1.0 }, clusterCount: 1);

        r.MathScore.Should().Be(0.0, "singleton VAF ⇒ MATH = 0 (§6.1)");
        r.ShannonDiversity.Should().BeApproximately(0.0, Tolerance, "one occupied clone ⇒ H = 0 (INV-03)");
        r.SubcloneCount.Should().Be(1, "a single mutation occupies one cluster (INV-05)");
        r.SubclonalFraction.Should().Be(0.0, "CCF = 1.0 ≥ 0.95 ⇒ clonal ⇒ subclonal fraction 0");
        AssertWellFormedResult(r, n: 1, clusterCount: 1);
    }

    [Test]
    [CancelAfter(5_000)]
    public void AnalyzeHeterogeneity_AllZeroVafs_ThrowsArgumentException_NoInfNaN()
    {
        // The aggregate computes MATH first; all-zero VAFs ⇒ zero median ⇒ ArgumentException propagates,
        // and NO HeterogeneityResult with an Inf/NaN MathScore is ever returned.
        Action act = () => AnalyzeHeterogeneity(
            new[] { 0.0, 0.0, 0.0 }, new[] { 0.2, 0.5, 0.9 }, clusterCount: 2);

        act.Should().Throw<ArgumentException>("zero VAF median ⇒ undefined MATH ⇒ rejected (§3.3, §6.1)")
            .Which.ParamName.Should().Be("ccfDistribution");
    }

    [Test]
    [CancelAfter(15_000)]
    public void AnalyzeHeterogeneity_RandomInContract_AlwaysWellFormed()
    {
        // Fuzz the aggregate over random in-contract inputs (positive VAF median, CCF ∈ [0,1], k ∈ [1,n]):
        // it must ALWAYS return a well-formed result — finite non-negative MATH, H ≥ 0,
        // subclone count ∈ [1,k], subclonal fraction ∈ [0,1] — never a crash or NaN/Inf.
        var rng = new Random(116_005);
        for (int i = 0; i < 3_000; i++)
        {
            int n = rng.Next(1, 30);
            var vaf = new double[n];
            var ccf = new double[n];
            for (int k = 0; k < n; k++)
            {
                vaf[k] = rng.NextDouble() * 0.9 + 0.05; // (0,1], guarantees median > 0
                ccf[k] = rng.NextDouble();              // [0,1)
            }

            int clusterCount = rng.Next(1, n + 1); // [1, n]

            var r = AnalyzeHeterogeneity(vaf, ccf, clusterCount);
            AssertWellFormedResult(r, n, clusterCount);
        }
    }

    [Test]
    [CancelAfter(5_000)]
    public void AnalyzeHeterogeneity_MismatchedLengths_ThrowsArgumentException()
    {
        Action act = () => AnalyzeHeterogeneity(new[] { 0.3, 0.5 }, new[] { 0.4 }, clusterCount: 1);
        act.Should().Throw<ArgumentException>("VAF and CCF lists must have equal length (§3.3)");
    }

    [Test]
    public void AnalyzeHeterogeneity_NullLists_ThrowArgumentNull()
    {
        Action a1 = () => AnalyzeHeterogeneity(null!, new[] { 0.4 }, 1);
        Action a2 = () => AnalyzeHeterogeneity(new[] { 0.4 }, null!, 1);
        a1.Should().Throw<ArgumentNullException>();
        a2.Should().Throw<ArgumentNullException>();
    }

    #endregion
}
