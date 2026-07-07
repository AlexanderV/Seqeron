using static Seqeron.Genomics.Oncology.OncologyAnalyzer;

namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the Oncology cancer-cell-fraction (CCF) POINT-ESTIMATE computation — ONCO-CCF-001.
/// The unit under test is the deterministic closed-form estimator implemented in
/// src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs:
/// <see cref="OncologyAnalyzer.EstimateCcf(double, double, int, int)"/> —
/// <c>EstimateCcf(double vaf, double purity, int tumorCopyNumber, int multiplicity)</c> — which converts a
/// somatic variant's VAF, the sample purity ρ, the local tumor copy number N_T and the integer mutation
/// multiplicity m into the fraction of cancer cells carrying the mutation.
///
/// This file is scoped to the CCF COMPUTATION itself. Classification of clonality from a GIVEN CCF is
/// ONCO-CLONAL-001 (OncologyClonalityFuzzTests.cs) and is out of scope here.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate / boundary / extreme inputs to a unit and asserts the code NEVER fails in an
/// undisciplined way: no hang, no nonsense output, and no *unhandled* runtime fault. Every input must resolve to
/// EITHER a well-defined, theory-correct value OR a *documented, intentional* outcome — here, an
/// ArgumentOutOfRangeException for vaf ∉ [0,1] / purity ∉ (0,1] / tumorCopyNumber &lt; 1, or an ArgumentException
/// for multiplicity ∉ [1, tumorCopyNumber]. There is NO Inf or NaN leaking out of the division by ρ·m: the
/// formula's only denominator, ρ·m, is structurally protected by the (0,1] purity guard and the [1, N_T]
/// multiplicity guard, so a finite, defined CCF is ALWAYS returned for in-contract inputs.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: ONCO-CCF-001 — cancer cell fraction point estimation (Oncology)
/// Checklist: docs/checklists/03_FUZZING.md, row 115.
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — граничні значення: 0, -1, MaxInt, empty.
///     Targets (checklist row 115): "purity=0, CN=0, VAF&gt;purity".
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// The targets map onto the documented contract (Cancer_Cell_Fraction_Estimation.md) as follows:
///   • purity = 0   ⇔ ρ in the denominator (ρ·m) would force a divide-by-zero / +∞. The contract restricts
///                     purity to (0, 1]; purity = 0 (and any ρ ≤ 0) is REJECTED with ArgumentOutOfRangeException,
///                     NOT a leaked Inf/NaN (§3.1, §3.3; INV-CCF-01).
///   • CN = 0       ⇔ a local tumor copy number of 0 (homozygous deletion at the locus) makes the model degenerate.
///                     The contract requires N_T ≥ 1; tumorCopyNumber &lt; 1 is REJECTED with
///                     ArgumentOutOfRangeException, NOT a divide-by-zero or a 0/0 NaN (§3.1, §3.3).
///   • VAF &gt; purity ⇔ a VAF higher than purity alone could explain drives the RAW formula above 1.0. The
///                     reported CCF is CAPPED to 1.0 (a mutation in all cancer cells has CCF = 1) while the
///                     uncapped value is exposed as RawCcf — no &gt;1 nonsense leaks through the reported field
///                     (§3.2, §4.A step 3, §6.1; INV-CCF-01).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The documented contract under test (Cancer_Cell_Fraction_Estimation.md)
/// ───────────────────────────────────────────────────────────────────────────
///   • CCF formula: RawCcf = VAF·(ρ·N_T + 2(1−ρ)) / (ρ·m)                       (§2.A Core Model, §4.A, §5.3.A)
///   • Reported Ccf = min(1, RawCcf), Ccf ∈ [0, 1]                              (§3.2, §4.A step 3; INV-CCF-01)
///   • RawCcf may exceed 1 (uncapped); Ccf never does                          (§3.2, §6.1 "raw CCF > 1")
///   • CCF = 0 ⇔ VAF = 0 (numerator is VAF·positive constant)                  (INV-CCF-03, §6.1 "VAF = 0")
///   • CCF strictly increases with VAF, other inputs fixed (positive slope)    (INV-CCF-02)
///   • vaf ∉ [0,1] (or NaN) ⇒ ArgumentOutOfRangeException                      (§3.3)
///   • purity ∉ (0,1] (≤ 0, &gt; 1, or NaN) ⇒ ArgumentOutOfRangeException        (§3.3)
///   • tumorCopyNumber &lt; 1 ⇒ ArgumentOutOfRangeException                      (§3.3)
///   • multiplicity ∉ [1, tumorCopyNumber] ⇒ ArgumentException                 (§3.3)
///   • Worked example: EstimateCcf(0.20, 0.80, 2, 1).Ccf == 0.5                 (§7.1)
///
/// No source bug was found; no test was weakened.
/// All randomness is LOCALLY seeded (new Random(seed)); no shared static Rng.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public sealed class OncologyCcfFuzzTests
{
    private const double Tolerance = 1e-9;

    // ── Well-formed-result assertion helper ──────────────────────────────────
    // Pin the documented numeric contract on EVERY returned estimate so a fuzz
    // test cannot rubber-stamp a malformed result (NaN/Inf CCF, reported CCF
    // outside [0, 1], or a reported value diverging from the documented cap of
    // the raw formula) green.
    private static void AssertWellFormedEstimate(CcfEstimate est)
    {
        double.IsNaN(est.Ccf).Should().BeFalse("reported CCF must never be NaN (no 0/0 from a guarded denominator)");
        double.IsInfinity(est.Ccf).Should().BeFalse("reported CCF must never be ±∞ (ρ·m denominator is guarded)");
        double.IsNaN(est.RawCcf).Should().BeFalse("raw CCF must never be NaN for in-contract inputs");
        double.IsInfinity(est.RawCcf).Should().BeFalse("raw CCF must never be ±∞ for in-contract inputs");

        est.Ccf.Should().BeGreaterThanOrEqualTo(0.0, "reported CCF ∈ [0, 1] (INV-CCF-01)");
        est.Ccf.Should().BeLessThanOrEqualTo(1.0, "reported CCF is capped at 1 (INV-CCF-01, §4.A step 3)");
        est.RawCcf.Should().BeGreaterThanOrEqualTo(0.0, "raw CCF is VAF·(positive constant)/(positive denominator)");
        est.Ccf.Should().Be(Math.Min(1.0, est.RawCcf),
            "reported CCF is exactly min(1, raw) — capping is the ONLY transform (§4.A step 3)");
    }

    // Reference implementation of the documented closed form, used to cross-check
    // EstimateCcf on random in-contract inputs (mutation testing of the formula).
    private static double ReferenceRawCcf(double vaf, double purity, int tumorCopyNumber, int multiplicity)
        => vaf * (purity * tumorCopyNumber + 2.0 * (1.0 - purity)) / (purity * multiplicity);

    #region ONCO-CCF-001 — positive sanity (worked example, clonal ≈ 1, subclonal < 1)

    [Test]
    [CancelAfter(5_000)]
    public void EstimateCcf_DocWorkedExample_YieldsExactlyHalf()
    {
        // Cancer_Cell_Fraction_Estimation.md §7.1:
        //   EstimateCcf(vaf: 0.20, purity: 0.80, tumorCopyNumber: 2, multiplicity: 1) ⇒ Ccf == 0.5
        //   0.20·(0.8·2 + 2·0.2)/(0.8·1) = 0.20·2.0/0.8 = 0.5
        var est = EstimateCcf(vaf: 0.20, purity: 0.80, tumorCopyNumber: 2, multiplicity: 1);

        est.RawCcf.Should().BeApproximately(0.5, Tolerance, "hand-computed worked example (§7.1)");
        est.Ccf.Should().BeApproximately(0.5, Tolerance, "raw is ≤ 1 so reported equals raw");
        AssertWellFormedEstimate(est);
    }

    [Test]
    [CancelAfter(5_000)]
    public void EstimateCcf_ClonalConfiguration_YieldsCcfApproximatelyOne()
    {
        // A clonal heterozygous SNV in a pure, diploid region: VAF ≈ ½ at ρ = 1, N_T = 2, m = 1.
        //   0.5·(1·2 + 2·0)/(1·1) = 0.5·2 = 1.0  ⇒ present in all cancer cells.
        var est = EstimateCcf(vaf: 0.5, purity: 1.0, tumorCopyNumber: 2, multiplicity: 1);

        est.Ccf.Should().BeApproximately(1.0, Tolerance, "a clonal mutation has CCF = 1 (INV-CCF-01)");
        est.RawCcf.Should().BeApproximately(1.0, Tolerance, "raw equals 1.0, no capping needed");
        AssertWellFormedEstimate(est);
    }

    [Test]
    [CancelAfter(5_000)]
    public void EstimateCcf_SubclonalConfiguration_YieldsCcfBelowOne()
    {
        // Half the clonal VAF in the same pure diploid context ⇒ half the CCF, strictly subclonal.
        //   0.25·(1·2 + 0)/(1·1) = 0.5
        var est = EstimateCcf(vaf: 0.25, purity: 1.0, tumorCopyNumber: 2, multiplicity: 1);

        est.Ccf.Should().BeApproximately(0.5, Tolerance, "half the clonal VAF ⇒ half the CCF");
        est.Ccf.Should().BeLessThan(1.0, "a subclonal mutation has CCF < 1");
        AssertWellFormedEstimate(est);
    }

    [Test]
    [CancelAfter(5_000)]
    public void EstimateCcf_RandomInContractInputs_MatchReferenceFormulaExactly()
    {
        // Mutation testing of the formula: any drift from VAF·(ρ·N_T + 2(1−ρ))/(ρ·m) is caught.
        var rng = new Random(115_001);
        for (int i = 0; i < 5_000; i++)
        {
            double vaf = rng.NextDouble();                 // [0, 1)
            double purity = rng.NextDouble() * 0.999 + 0.001; // (0, 1]
            int cn = rng.Next(1, 9);                        // [1, 8]
            int m = rng.Next(1, cn + 1);                    // [1, N_T]

            var est = EstimateCcf(vaf, purity, cn, m);

            double expectedRaw = ReferenceRawCcf(vaf, purity, cn, m);
            est.RawCcf.Should().BeApproximately(expectedRaw, 1e-9,
                "RawCcf must equal the documented closed form (§2.A Core Model)");
            est.Ccf.Should().BeApproximately(Math.Min(1.0, expectedRaw), 1e-9,
                "reported CCF is exactly min(1, raw) (§4.A step 3)");
            AssertWellFormedEstimate(est);
        }
    }

    #endregion

    #region ONCO-CCF-001 — BE: purity = 0 (and ρ ≤ 0 / ρ > 1 / NaN ⇒ no Inf/NaN leak)

    [Test]
    public void EstimateCcf_PurityZero_ThrowsArgumentOutOfRange_NoInfinityLeak()
    {
        // ρ = 0 would make the ρ·m denominator zero ⇒ a leaked +∞ / NaN. The (0,1] guard rejects it.
        Action act = () => EstimateCcf(vaf: 0.3, purity: 0.0, tumorCopyNumber: 2, multiplicity: 1);

        act.Should().Throw<ArgumentOutOfRangeException>("purity ∈ (0, 1]; the CCF formula divides by purity (§3.3)")
            .Which.ParamName.Should().Be("purity");
    }

    [TestCase(-0.0)]
    [TestCase(-1.0)]
    [TestCase(-1e-12)]
    [TestCase(double.MinValue)]
    public void EstimateCcf_NonPositivePurity_ThrowsArgumentOutOfRange(double purity)
    {
        Action act = () => EstimateCcf(vaf: 0.3, purity: purity, tumorCopyNumber: 2, multiplicity: 1);

        act.Should().Throw<ArgumentOutOfRangeException>("purity must be strictly positive (§3.3)")
            .Which.ParamName.Should().Be("purity");
    }

    [TestCase(1.0000001)]
    [TestCase(2.0)]
    [TestCase(double.MaxValue)]
    [TestCase(double.PositiveInfinity)]
    [TestCase(double.NaN)]
    public void EstimateCcf_PurityAboveOneOrNaN_ThrowsArgumentOutOfRange(double purity)
    {
        Action act = () => EstimateCcf(vaf: 0.3, purity: purity, tumorCopyNumber: 2, multiplicity: 1);

        act.Should().Throw<ArgumentOutOfRangeException>("purity ∉ (0, 1] (or NaN) is rejected (§3.3)")
            .Which.ParamName.Should().Be("purity");
    }

    [Test]
    [CancelAfter(5_000)]
    public void EstimateCcf_PurityJustAboveZero_StaysFinite_NoOverflowToInfinity()
    {
        // The smallest in-contract purity exercises the ρ·m denominator at its minimum without leaking ∞:
        // RawCcf grows large but must remain a finite double, and reported CCF is capped to 1.
        var est = EstimateCcf(vaf: 1.0, purity: 1e-9, tumorCopyNumber: 1, multiplicity: 1);

        AssertWellFormedEstimate(est);
        est.RawCcf.Should().BeGreaterThan(1.0, "a tiny purity inflates the raw estimate well past 1");
        est.Ccf.Should().Be(1.0, "the inflated raw value is capped to 1 (INV-CCF-01)");
    }

    #endregion

    #region ONCO-CCF-001 — BE: CN = 0 (homozygous deletion ⇒ no divide-by-zero / NaN, rejected)

    [Test]
    public void EstimateCcf_CopyNumberZero_ThrowsArgumentOutOfRange_NoDegenerateResult()
    {
        // N_T = 0 is a homozygous deletion at the locus — the VAF→CCF model degenerates. The N_T ≥ 1 guard
        // rejects it rather than emitting a degenerate (possibly 0/0 or biased) CCF.
        Action act = () => EstimateCcf(vaf: 0.3, purity: 0.8, tumorCopyNumber: 0, multiplicity: 1);

        act.Should().Throw<ArgumentOutOfRangeException>("tumor copy number must be ≥ 1 (§3.3)")
            .Which.ParamName.Should().Be("tumorCopyNumber");
    }

    [TestCase(-1)]
    [TestCase(int.MinValue)]
    public void EstimateCcf_NegativeCopyNumber_ThrowsArgumentOutOfRange(int cn)
    {
        Action act = () => EstimateCcf(vaf: 0.3, purity: 0.8, tumorCopyNumber: cn, multiplicity: 1);

        act.Should().Throw<ArgumentOutOfRangeException>("a negative copy number is out of range (§3.3)")
            .Which.ParamName.Should().Be("tumorCopyNumber");
    }

    [Test]
    public void EstimateCcf_CopyNumberZero_RejectedEvenWhenMultiplicityWouldAlsoBeInvalid()
    {
        // When N_T = 0, [1, N_T] is empty so multiplicity can never be valid either; the copy-number guard
        // must fire first and deterministically (validation order is documented N_T then m, §4.A step 1).
        Action act = () => EstimateCcf(vaf: 0.3, purity: 0.8, tumorCopyNumber: 0, multiplicity: 1);

        act.Should().Throw<ArgumentOutOfRangeException>("copy number is validated before multiplicity (§3.3)");
    }

    [Test]
    [CancelAfter(5_000)]
    public void EstimateCcf_CopyNumberOne_MinimumValid_NoDivideByZero()
    {
        // N_T = 1 is the smallest legal copy number; the estimate must be finite and well-formed.
        var est = EstimateCcf(vaf: 0.2, purity: 0.8, tumorCopyNumber: 1, multiplicity: 1);

        AssertWellFormedEstimate(est);
        est.RawCcf.Should().BeApproximately(ReferenceRawCcf(0.2, 0.8, 1, 1), Tolerance,
            "minimum copy number still follows the documented formula");
    }

    #endregion

    #region ONCO-CCF-001 — BE: VAF > purity (raw exceeds 1 ⇒ reported CCF capped, no >1 leak)

    [Test]
    [CancelAfter(5_000)]
    public void EstimateCcf_VafExceedingPurity_RawAboveOne_ReportedCappedToOne()
    {
        // VAF (0.9) far above what purity (0.5) alone could explain in a diploid region drives the raw value
        // past 1: 0.9·(0.5·2 + 2·0.5)/(0.5·1) = 0.9·2/0.5 = 3.6. Reported CCF must be CAPPED to 1.0.
        var est = EstimateCcf(vaf: 0.9, purity: 0.5, tumorCopyNumber: 2, multiplicity: 1);

        est.RawCcf.Should().BeApproximately(3.6, Tolerance, "uncapped raw value is exposed (§3.2)");
        est.RawCcf.Should().BeGreaterThan(1.0, "VAF > purity overshoots the [0,1] range");
        est.Ccf.Should().Be(1.0, "reported CCF is capped to 1 — no >1 nonsense leaks (INV-CCF-01, §6.1)");
        AssertWellFormedEstimate(est);
    }

    [Test]
    [CancelAfter(5_000)]
    public void EstimateCcf_VafEqualToPurity_DiploidSingleCopy_RawIsTwo_Capped()
    {
        // VAF == purity in a diploid single-copy locus: 0.6·(0.6·2 + 2·0.4)/(0.6·1) = 0.6·2/0.6 = 2.0 ⇒ capped.
        var est = EstimateCcf(vaf: 0.6, purity: 0.6, tumorCopyNumber: 2, multiplicity: 1);

        est.RawCcf.Should().BeApproximately(2.0, Tolerance);
        est.Ccf.Should().Be(1.0, "raw 2.0 is capped to 1");
        AssertWellFormedEstimate(est);
    }

    [Test]
    [CancelAfter(10_000)]
    public void EstimateCcf_RandomVafAbovePurity_NeverReportsAboveOne()
    {
        // Sweep many (VAF > purity) configurations; the reported CCF must NEVER exceed 1 even as raw does.
        var rng = new Random(115_002);
        bool sawRawAboveOne = false;
        for (int i = 0; i < 5_000; i++)
        {
            double purity = rng.NextDouble() * 0.5 + 0.05;   // (0.05, 0.55]
            double vaf = purity + rng.NextDouble() * (1.0 - purity); // strictly ≥ purity, ≤ 1
            int cn = rng.Next(1, 5);
            int m = rng.Next(1, cn + 1);

            var est = EstimateCcf(vaf, purity, cn, m);

            est.Ccf.Should().BeLessThanOrEqualTo(1.0, "reported CCF is always capped (INV-CCF-01)");
            AssertWellFormedEstimate(est);
            if (est.RawCcf > 1.0)
            {
                sawRawAboveOne = true;
            }
        }

        sawRawAboveOne.Should().BeTrue("the VAF>purity regime is expected to produce raw values above 1");
    }

    #endregion

    #region ONCO-CCF-001 — BE: VAF / multiplicity boundaries & documented invariants

    [TestCase(-1e-12)]
    [TestCase(-1.0)]
    [TestCase(1.0000001)]
    [TestCase(2.0)]
    [TestCase(double.NaN)]
    [TestCase(double.PositiveInfinity)]
    public void EstimateCcf_VafOutOfRangeOrNaN_ThrowsArgumentOutOfRange(double vaf)
    {
        Action act = () => EstimateCcf(vaf: vaf, purity: 0.8, tumorCopyNumber: 2, multiplicity: 1);

        act.Should().Throw<ArgumentOutOfRangeException>("vaf ∉ [0, 1] (or NaN) is rejected (§3.3)")
            .Which.ParamName.Should().Be("vaf");
    }

    [Test]
    [CancelAfter(5_000)]
    public void EstimateCcf_VafZero_YieldsCcfZero()
    {
        // INV-CCF-03 / §6.1: CCF = 0 ⇔ VAF = 0 (numerator is VAF·positive constant). No divide-by-zero.
        var est = EstimateCcf(vaf: 0.0, purity: 0.7, tumorCopyNumber: 3, multiplicity: 2);

        est.Ccf.Should().Be(0.0, "CCF = 0 ⇔ VAF = 0 (INV-CCF-03)");
        est.RawCcf.Should().Be(0.0, "raw numerator is 0 when VAF is 0");
        AssertWellFormedEstimate(est);
    }

    [TestCase(0)]   // m < 1
    [TestCase(-1)]
    [TestCase(3)]   // m > tumorCopyNumber (= 2)
    [TestCase(int.MaxValue)]
    [TestCase(int.MinValue)]
    public void EstimateCcf_MultiplicityOutOfRange_ThrowsArgumentException(int multiplicity)
    {
        Action act = () => EstimateCcf(vaf: 0.3, purity: 0.8, tumorCopyNumber: 2, multiplicity: multiplicity);

        act.Should().Throw<ArgumentException>("multiplicity ∉ [1, tumorCopyNumber] is rejected (§3.3)")
            .Which.ParamName.Should().Be("multiplicity");
    }

    [Test]
    [CancelAfter(5_000)]
    public void EstimateCcf_StrictlyIncreasesWithVaf_OtherInputsFixed()
    {
        // INV-CCF-02: CCF strictly increases with VAF (linear, positive slope). Verify on the RAW value so
        // the cap does not mask monotonicity once both samples are below 1.
        var rng = new Random(115_003);
        for (int i = 0; i < 2_000; i++)
        {
            double purity = rng.NextDouble() * 0.9 + 0.1;  // [0.1, 1.0)
            int cn = rng.Next(1, 6);
            int m = rng.Next(1, cn + 1);
            double vafLow = rng.NextDouble() * 0.4;        // [0, 0.4)
            double vafHigh = vafLow + (rng.NextDouble() * 0.4 + 1e-3); // strictly greater, ≤ ~0.8

            var low = EstimateCcf(vafLow, purity, cn, m);
            var high = EstimateCcf(vafHigh, purity, cn, m);

            high.RawCcf.Should().BeGreaterThan(low.RawCcf, "raw CCF is linear in VAF with positive slope (INV-CCF-02)");
            AssertWellFormedEstimate(low);
            AssertWellFormedEstimate(high);
        }
    }

    #endregion

    #region ONCO-CCF-001 — robustness: large / extreme in-contract sweeps (no hang, contract holds)

    [Test]
    [CancelAfter(15_000)]
    public void EstimateCcf_LargeRandomSweep_AlwaysWellFormed()
    {
        var rng = new Random(115_010);
        for (int i = 0; i < 50_000; i++)
        {
            double vaf = rng.NextDouble();
            double purity = rng.NextDouble() * 0.999 + 1e-3;
            int cn = rng.Next(1, 50);
            int m = rng.Next(1, cn + 1);

            var est = EstimateCcf(vaf, purity, cn, m);
            AssertWellFormedEstimate(est);
        }
    }

    [Test]
    [CancelAfter(5_000)]
    public void EstimateCcf_ExtremeCopyNumber_StaysFinite()
    {
        // A massively amplified locus (large N_T) inflates the numerator but the result must stay finite.
        var est = EstimateCcf(vaf: 1.0, purity: 1.0, tumorCopyNumber: int.MaxValue, multiplicity: 1);

        AssertWellFormedEstimate(est);
        est.RawCcf.Should().BeApproximately(ReferenceRawCcf(1.0, 1.0, int.MaxValue, 1), 1e-3,
            "extreme but in-contract copy number still follows the formula");
    }

    #endregion
}
