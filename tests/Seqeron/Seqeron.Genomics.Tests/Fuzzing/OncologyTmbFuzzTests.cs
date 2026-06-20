using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Oncology;
using static Seqeron.Genomics.Oncology.OncologyAnalyzer;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Fuzz tests for the Oncology tumor-mutational-burden area — ONCO-TMB-001.
/// The units under test are the deterministic arithmetic-ratio entry points
/// <see cref="OncologyAnalyzer.CalculateTMB(int, double)"/>, its
/// <see cref="OncologyAnalyzer.CalculateTMB(System.Collections.Generic.IEnumerable{OncologyAnalyzer.SomaticCall}, double)"/>
/// overload, and the cutoff classifier <see cref="OncologyAnalyzer.ClassifyTMB"/>,
/// implemented in src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate / boundary / extreme inputs to a unit and asserts
/// that the code NEVER fails in an undisciplined way: no hang, no nonsense
/// output, and no *unhandled* runtime fault (DivideByZero / Overflow / a NaN or
/// ±Infinity TMB leaking out). Every input must resolve to EITHER a well-defined,
/// theory-correct value OR a *documented, intentional* outcome (here, an
/// <see cref="ArgumentOutOfRangeException"/> for a non-positive / non-finite
/// panel size, a negative count, or a non-finite / negative TMB to classify).
/// For TMB the headline hazards are:
///   • panel size 0 → the denominator is 0; the contract is a documented
///     ArgumentOutOfRangeException ("TMB is undefined at 0 Mb"), NEVER a
///     DivideByZero, a +Infinity TMB, nor a NaN escaping to the caller (§3.3, §6.1);
///   • negative panel size → documented ArgumentOutOfRangeException, never a
///     silently NEGATIVE TMB (§3.3);
///   • huge counts (mutationCount near Int32.MaxValue) and/or a tiny panel →
///     the division MUST be performed in double, so NO integer overflow and the
///     result stays FINITE and non-negative (§2.2, INV-02);
///   • zero mutations → TMB = 0 (quotient of 0), classified TMB-Low, no throw
///     and no DivideByZero (§6.1).
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: ONCO-TMB-001 — Tumor mutational burden (Oncology)
/// Checklist: docs/checklists/03_FUZZING.md, row 92.
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — граничні значення: 0, -1, MaxInt, empty.
///     Targets (checklist row 92): "zero mutations, panel size 0, negative size,
///     huge counts".
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The documented contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// Tumor_Mutational_Burden.md (docs/algorithms/Oncology/Tumor_Mutational_Burden.md):
///   • TMB = mutationCount / targetRegionMb, in mut/Mb                    (§2.2, INV-01)
///   • TMB ≥ 0 for mutationCount ≥ 0, targetRegionMb > 0                  (INV-02)
///   • TMB is non-decreasing in count (fixed Mb) and non-increasing in Mb
///       (fixed count)                                                    (INV-03)
///   • mutationCount = 0 ⇒ TMB = 0                                        (§6.1)
///   • targetRegionMb is NaN / ±Infinity / ≤ 0 ⇒ ArgumentOutOfRangeException
///       ("TMB is undefined at 0 Mb")                                     (§3.3, §6.1)
///   • mutationCount < 0 ⇒ ArgumentOutOfRangeException                    (§3.3)
///   • Overload counts ONLY SomaticStatus.Somatic calls (Germline /
///       NotDetected excluded), then delegates                           (§4.1)
///   • Overload throws ArgumentNullException on a null collection         (§3.3)
///   • ClassifyTMB(tmb): High ⇔ tmb ≥ 10.0 (inclusive), else Low         (§2.2, INV-04)
///   • ClassifyTMB throws ArgumentOutOfRangeException for negative or
///       non-finite tmb                                                   (§3.3)
///
/// All randomness is LOCALLY seeded (new Random(seed)); no shared static Rng.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public sealed class OncologyTmbFuzzTests
{
    private const double Cutoff = TmbHighThreshold; // 10.0 mut/Mb, FDA inclusive

    // ── Well-formed-TMB assertion helper ─────────────────────────────────────
    // Pins the documented numeric contract on EVERY returned TMB: it must be
    // FINITE (no DivideByZero +Infinity, no NaN leaking through) and NON-negative
    // (INV-02). This is what stops a fuzz test from rubber-stamping an Infinity /
    // NaN / negative TMB green.
    private static void AssertWellFormedTmb(double tmb)
    {
        double.IsNaN(tmb).Should().BeFalse("TMB must never be NaN");
        double.IsInfinity(tmb).Should().BeFalse("TMB = count/Mb must be finite (Mb > 0)");
        tmb.Should().BeGreaterThanOrEqualTo(0.0, "TMB ≥ 0 for count ≥ 0, Mb > 0 (INV-02)");
    }

    private static SomaticCall Call(SomaticStatus status) =>
        new(
            new VariantObservation("chr1", 100, "A", "T", 30, 100, 0, 100),
            TumorVaf: 0.30,
            NormalVaf: 0.00,
            Status: status,
            SomaticScore: status == SomaticStatus.Somatic ? 0.30 : 0.0);

    #region ONCO-TMB-001 — Tumor mutational burden (positive sanity)

    // ── POSITIVE sanity: hand-computed TMB and the correct High/Low class ─────
    [Test]
    public void CalculateTMB_KnownCountAndPanel_HandComputedValueAndClass()
    {
        // 30 somatic mutations / 1.5 Mb = 20 mut/Mb; 20 ≥ 10 ⇒ TMB-High.
        double tmb = CalculateTMB(mutationCount: 30, targetRegionMb: 1.5);

        AssertWellFormedTmb(tmb);
        tmb.Should().BeApproximately(20.0, 1e-12, "30 / 1.5 = 20 mut/Mb (§2.2)");
        ClassifyTMB(tmb).Should().Be(TmbStatus.High, "20 ≥ 10 ⇒ TMB-High (INV-04)");
    }

    [Test]
    public void CalculateTMB_DocWorkedExample_ExactlyAtCutoffIsHigh()
    {
        // Docs §7.1: 11 mutations / 1.1 Mb = 10.0 mut/Mb = exactly the FDA cutoff.
        double tmb = CalculateTMB(mutationCount: 11, targetRegionMb: 1.1);

        AssertWellFormedTmb(tmb);
        tmb.Should().BeApproximately(10.0, 1e-9, "11 / 1.1 = 10.0 mut/Mb (§7.1)");
        ClassifyTMB(tmb).Should().Be(TmbStatus.High, "cutoff is inclusive ≥ 10 (INV-04)");
    }

    [Test]
    public void CalculateTMB_LowBurden_ClassifiesLow()
    {
        // 5 mutations / 1.1 Mb ≈ 4.55 mut/Mb < 10 ⇒ TMB-Low.
        double tmb = CalculateTMB(mutationCount: 5, targetRegionMb: 1.1);

        AssertWellFormedTmb(tmb);
        tmb.Should().BeApproximately(5.0 / 1.1, 1e-12);
        ClassifyTMB(tmb).Should().Be(TmbStatus.Low, "≈4.55 < 10 ⇒ TMB-Low (INV-04)");
    }

    #endregion

    #region ONCO-TMB-001 — BE: zero mutations (TMB = 0, no DivideByZero, TMB-Low)

    [Test]
    public void CalculateTMB_ZeroMutations_IsZeroAndLow_NoDivideByZero()
    {
        // §6.1: mutationCount = 0 ⇒ TMB = 0 (quotient of 0), classified TMB-Low.
        // The hazard guarded here is that a 0 numerator must NOT trip a guard or
        // produce a NaN — it is a perfectly valid, finite TMB of exactly 0.
        double tmb = CalculateTMB(mutationCount: 0, targetRegionMb: 1.1);

        AssertWellFormedTmb(tmb);
        tmb.Should().Be(0.0, "0 / Mb = 0 (§6.1)");
        ClassifyTMB(tmb).Should().Be(TmbStatus.Low, "0 < 10 ⇒ TMB-Low");
    }

    [Test]
    public void CalculateTMB_ZeroMutations_AcrossFuzzedPanelSizes_AlwaysZeroLow()
    {
        var rng = new Random(92_0001);
        for (int i = 0; i < 500; i++)
        {
            // Any positive, finite panel size — tiny to huge.
            double mb = Math.Pow(10, rng.NextDouble() * 8 - 4); // [1e-4, 1e4)
            double tmb = CalculateTMB(0, mb);

            AssertWellFormedTmb(tmb);
            tmb.Should().Be(0.0, "zero mutations ⇒ TMB = 0 for every Mb > 0");
            ClassifyTMB(tmb).Should().Be(TmbStatus.Low);
        }
    }

    [Test]
    public void CalculateTMB_Overload_NoSomaticCalls_IsZeroAndLow()
    {
        // The overload counts ONLY Somatic calls; an all-Germline/NotDetected set
        // (or an empty set) yields a count of 0 ⇒ TMB = 0, no DivideByZero.
        var calls = new[]
        {
            Call(SomaticStatus.Germline),
            Call(SomaticStatus.NotDetected),
            Call(SomaticStatus.Germline),
        };

        double tmb = CalculateTMB(calls, targetRegionMb: 1.1);
        AssertWellFormedTmb(tmb);
        tmb.Should().Be(0.0, "no Somatic calls ⇒ count 0 ⇒ TMB 0 (§4.1)");

        CalculateTMB(Array.Empty<SomaticCall>(), 1.1).Should().Be(0.0, "empty ⇒ 0");
    }

    #endregion

    #region ONCO-TMB-001 — BE: panel size 0 (documented throw, NEVER DivideByZero/Infinity)

    [Test]
    public void CalculateTMB_PanelSizeZero_Throws_NotDivideByZeroNorInfinity()
    {
        // §3.3 / §6.1: a 0-Mb denominator makes TMB undefined; the contract is a
        // documented ArgumentOutOfRangeException — NOT a DivideByZero, and NOT a
        // +Infinity TMB leaking out of double division.
        var act = () => CalculateTMB(mutationCount: 7, targetRegionMb: 0.0);

        act.Should().Throw<ArgumentOutOfRangeException>()
            .Which.ParamName.Should().Be("targetRegionMb");
    }

    [Test]
    public void CalculateTMB_PanelSizeZero_EvenWithZeroCount_Throws()
    {
        // The 0-Mb guard fires regardless of the numerator: 0/0 must not become NaN.
        var act = () => CalculateTMB(mutationCount: 0, targetRegionMb: 0.0);
        act.Should().Throw<ArgumentOutOfRangeException>()
            .Which.ParamName.Should().Be("targetRegionMb");
    }

    [Test]
    public void CalculateTMB_NegativeZeroPanel_Throws()
    {
        // -0.0 is ≤ 0 and must be rejected like +0.0 (no signed-zero loophole).
        var act = () => CalculateTMB(5, -0.0);
        act.Should().Throw<ArgumentOutOfRangeException>();
    }

    [Test]
    public void CalculateTMB_NonFinitePanel_Throws()
    {
        // NaN / ±Infinity panel sizes are rejected (§3.3) — they must not produce
        // a NaN or 0 TMB.
        foreach (double bad in new[] { double.NaN, double.PositiveInfinity, double.NegativeInfinity })
        {
            var act = () => CalculateTMB(5, bad);
            act.Should().Throw<ArgumentOutOfRangeException>(
                "panel size {0} is not finite > 0 (§3.3)", bad);
        }
    }

    [Test]
    public void CalculateTMB_Overload_PanelSizeZero_Throws_NoDivideByZero()
    {
        // The overload delegates to the int/double core, so the 0-Mb guard must
        // still fire even when there ARE somatic calls to count.
        var calls = new[] { Call(SomaticStatus.Somatic), Call(SomaticStatus.Somatic) };
        var act = () => CalculateTMB(calls, targetRegionMb: 0.0);
        act.Should().Throw<ArgumentOutOfRangeException>()
            .Which.ParamName.Should().Be("targetRegionMb");
    }

    #endregion

    #region ONCO-TMB-001 — BE: negative panel size (documented throw, NEVER negative TMB)

    [Test]
    public void CalculateTMB_NegativePanelSize_Throws_NoNegativeTmb()
    {
        // §3.3: a negative megabase denominator is rejected; a silently NEGATIVE
        // TMB (count / -Mb) must never escape.
        var act = () => CalculateTMB(mutationCount: 12, targetRegionMb: -1.1);
        act.Should().Throw<ArgumentOutOfRangeException>()
            .Which.ParamName.Should().Be("targetRegionMb");
    }

    [Test]
    public void CalculateTMB_FuzzedNegativePanelSizes_AlwaysThrow()
    {
        var rng = new Random(92_0002);
        for (int i = 0; i < 500; i++)
        {
            double negMb = -Math.Pow(10, rng.NextDouble() * 8 - 4); // (-1e4, -1e-4]
            int count = rng.Next(0, int.MaxValue);

            var act = () => CalculateTMB(count, negMb);
            act.Should().Throw<ArgumentOutOfRangeException>(
                "negative Mb {0} must be rejected, not give a negative TMB", negMb);
        }
    }

    [Test]
    public void CalculateTMB_NegativeMutationCount_Throws()
    {
        // §3.3: a negative count is rejected (it has no biological meaning).
        var act = () => CalculateTMB(mutationCount: -1, targetRegionMb: 1.1);
        act.Should().Throw<ArgumentOutOfRangeException>()
            .Which.ParamName.Should().Be("mutationCount");

        var actMin = () => CalculateTMB(int.MinValue, 1.1);
        actMin.Should().Throw<ArgumentOutOfRangeException>()
            .Which.ParamName.Should().Be("mutationCount");
    }

    #endregion

    #region ONCO-TMB-001 — BE: huge counts (double division, NO integer overflow)

    [Test]
    public void CalculateTMB_MaxIntCount_TinyPanel_FiniteNonNegative_NoOverflow()
    {
        // BE headline: a near-Int32.MaxValue count over a tiny panel must be
        // computed in DOUBLE (int / double → double). There must be NO integer
        // overflow and the result must stay FINITE and exactly count/Mb.
        const int count = int.MaxValue;          // 2_147_483_647
        const double mb = 0.001;                  // 1 kb panel

        double tmb = CalculateTMB(count, mb);

        AssertWellFormedTmb(tmb);
        tmb.Should().BeApproximately((double)count / mb, Math.Abs(tmb) * 1e-12);
        tmb.Should().Be((double)count / mb, "division is performed in double, exact");
        ClassifyTMB(tmb).Should().Be(TmbStatus.High, "astronomically high TMB ⇒ High");
    }

    [Test]
    public void CalculateTMB_MaxIntCount_LargePanel_MatchesExactDoubleRatio()
    {
        // Even the largest sane panels keep the result finite and correct.
        const int count = int.MaxValue;
        const double mb = 50.0; // whole-exome scale

        double tmb = CalculateTMB(count, mb);

        AssertWellFormedTmb(tmb);
        tmb.Should().Be((double)count / mb);
        ClassifyTMB(tmb).Should().Be(TmbStatus.High);
    }

    [Test]
    public void CalculateTMB_FuzzedHugeCountsAndPanels_StayFiniteAndExact()
    {
        var rng = new Random(92_0003);
        for (int i = 0; i < 1000; i++)
        {
            // Bias toward huge counts but include small ones too.
            int count = rng.NextDouble() < 0.5
                ? rng.Next(int.MaxValue - 1_000_000, int.MaxValue) + (rng.NextDouble() < 0.5 ? 1 : 0)
                : rng.Next(0, 10_000);
            double mb = Math.Pow(10, rng.NextDouble() * 8 - 4); // [1e-4, 1e4)

            double tmb = CalculateTMB(count, mb);

            AssertWellFormedTmb(tmb);
            tmb.Should().Be(count / mb, "TMB = count/Mb computed exactly in double");

            // Monotonicity (INV-03): non-decreasing in count at fixed Mb.
            if (count < int.MaxValue)
            {
                CalculateTMB(count + 1, mb).Should()
                    .BeGreaterThanOrEqualTo(tmb, "TMB non-decreasing in count (INV-03)");
            }
        }
    }

    [Test]
    public void CalculateTMB_Overload_AllSomatic_LargeSet_CountsAndDividesFinite()
    {
        // The overload's Somatic count is an O(n) int sum; with a large set the
        // resulting TMB must still be finite and equal count/Mb. (n kept modest
        // to avoid a slow test while still exercising the counting path.)
        const int n = 5000;
        const double mb = 1.1;
        var calls = Enumerable.Range(0, n).Select(_ => Call(SomaticStatus.Somatic)).ToArray();

        double tmb = CalculateTMB(calls, mb);
        AssertWellFormedTmb(tmb);
        tmb.Should().Be(n / mb, "overload counts n Somatic calls then divides");
        ClassifyTMB(tmb).Should().Be(TmbStatus.High);
    }

    #endregion

    #region ONCO-TMB-001 — BE: ClassifyTMB boundary & malformed inputs

    [Test]
    public void ClassifyTMB_AtCutoff_IsHigh_JustBelow_IsLow()
    {
        // INV-04: the FDA cutoff is INCLUSIVE at 10.0.
        ClassifyTMB(Cutoff).Should().Be(TmbStatus.High, "10.0 ≥ 10 ⇒ High (inclusive)");
        ClassifyTMB(Math.BitDecrement(Cutoff)).Should()
            .Be(TmbStatus.Low, "the next double below 10 ⇒ Low");
        ClassifyTMB(0.0).Should().Be(TmbStatus.Low, "0 ⇒ Low");
    }

    [Test]
    public void ClassifyTMB_NegativeOrNonFinite_Throws()
    {
        // §3.3: a negative or non-finite TMB cannot be classified.
        foreach (double bad in new[]
                 {
                     -0.0001, -1.0, double.NaN,
                     double.PositiveInfinity, double.NegativeInfinity,
                 })
        {
            var act = () => ClassifyTMB(bad);
            act.Should().Throw<ArgumentOutOfRangeException>(
                "TMB {0} is not finite ≥ 0 (§3.3)", bad);
        }
    }

    [Test]
    [CancelAfter(30_000)]
    public void CalculateThenClassify_FuzzedValidInputs_NeverThrow_NeverInfinite()
    {
        // End-to-end BE sweep: every (count ≥ 0, Mb > 0) pair must produce a
        // finite TMB that classifies without throwing, and the class must agree
        // with the documented inclusive ≥ 10 cutoff.
        var rng = new Random(92_0004);
        for (int i = 0; i < 2000; i++)
        {
            int count = rng.Next(0, int.MaxValue);
            double mb = Math.Pow(10, rng.NextDouble() * 8 - 4); // [1e-4, 1e4) > 0

            double tmb = CalculateTMB(count, mb);
            AssertWellFormedTmb(tmb);

            TmbStatus status = ClassifyTMB(tmb);
            status.Should().Be(tmb >= Cutoff ? TmbStatus.High : TmbStatus.Low,
                "class must match the inclusive ≥ 10 cutoff (INV-04)");
        }
    }

    #endregion

    #region ONCO-TMB-001 — BE: overload null guard

    [Test]
    public void CalculateTMB_Overload_NullCalls_Throws()
    {
        var act = () => CalculateTMB((IEnumerable<SomaticCall>)null!, 1.1);
        act.Should().Throw<ArgumentNullException>("§3.3: null collection is rejected");
    }

    #endregion
}
