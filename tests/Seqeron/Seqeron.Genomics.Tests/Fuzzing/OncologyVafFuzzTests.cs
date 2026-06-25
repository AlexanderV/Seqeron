using System;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Oncology;
using static Seqeron.Genomics.Oncology.OncologyAnalyzer;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Fuzz tests for the Oncology variant-allele-frequency area — ONCO-VAF-001.
/// The unit under test is the empirical VAF computation
/// <see cref="OncologyAnalyzer.CalculateVAF"/> together with the two
/// VAF-specific public computations that consume it — the Wilson-score
/// confidence interval <see cref="OncologyAnalyzer.CalculateVAFConfidenceInterval"/>
/// and the purity/ploidy correction <see cref="OncologyAnalyzer.AdjustVAFForPurity"/> —
/// implemented in src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs.
///
/// This is a DISTINCT unit from ONCO-SOMATIC-001 (row 87,
/// OncologySomaticFuzzTests.cs). Row 87 exercised <see cref="OncologyAnalyzer.CalculateVAF"/>
/// only as a per-variant helper inside the tumor/normal classifier. Here the VAF
/// formula is tested as a UNIT in its own right: its boundary contract on raw
/// read-count pairs, the headline integer-width / huge-depth overflow concern,
/// and the dedicated VAF-derived public surface (Wilson CI, purity correction)
/// that the somatic file never touches.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate / boundary / malformed inputs to a unit and asserts
/// that the code NEVER fails in an undisciplined way: no hang, no nonsense
/// output, and no *unhandled* runtime exception (DivideByZero / Overflow / NaN).
/// Every input must resolve to EITHER a well-defined, theory-correct value OR a
/// *documented, intentional* outcome (here, an
/// <see cref="ArgumentOutOfRangeException"/> for malformed read counts).
/// For the empirical VAF the headline hazards are:
///   • a DivideByZeroException / NaN computing alt/total when total == 0
///     (the documented contract is VAF = 0 at zero coverage, §6.1 — NOT a throw,
///     NOT a NaN);
///   • a VAF that escapes [0, 1] — in particular alt &gt; total must NOT yield a
///     VAF &gt; 1 (malformed counts are a documented ArgumentOutOfRangeException,
///     §3.3, §6.1);
///   • a negative read count silently producing a negative or out-of-range VAF
///     (documented ArgumentOutOfRangeException, §3.3);
///   • INTEGER OVERFLOW on huge depth — counts near Int32.MaxValue must not wrap
///     in the alt≤total comparison nor in the ratio; the ratio is taken in
///     double width so VAF stays correct and in [0, 1] (the headline new
///     contract for this unit vs row 87).
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: ONCO-VAF-001 — Variant allele frequency (Oncology)
/// Checklist: docs/checklists/03_FUZZING.md, row 88.
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — граничні значення: 0, -1, MaxInt, empty.
///     Targets (checklist row 88): "ref=alt=0, alt>total, negative counts,
///     huge depth".
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
/// The API takes (altReads, totalReads) where totalReads = ref + alt; thus the
/// target "ref=alt=0" maps to total = 0 (zero coverage), "alt>total" is the
/// malformed read-support case, and "huge depth" is counts near Int32.MaxValue.
///
/// ───────────────────────────────────────────────────────────────────────────
/// The documented contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// Variant_Allele_Frequency.md (docs/algorithms/Oncology/Variant_Allele_Frequency.md):
///   • VAF = altReads / totalReads, in [0, 1]                            (§2.2, INV-01)
///   • totalReads == 0 ⇒ VAF = 0 (no coverage ⇒ allele absent)           (§3.3, §6.1)
///   • altReads &gt; totalReads ⇒ ArgumentOutOfRangeException (VAF &gt; 1
///       artifact); negative counts ⇒ ArgumentOutOfRangeException        (§3.3, §6.1)
///   • CalculateVAFConfidenceInterval: Wilson score interval, bounded in
///       [0, 1] with lower ≤ center ≤ upper (INV-02, INV-03); throws for
///       totalReads == 0 (undefined with no trials) and for confidence
///       ≠ 0.95 / outside (0, 1)                                         (§2.2, §3.3, §6.1)
///   • AdjustVAFForPurity = vaf·(2(1−π) + π·ploidy)/π; throws for
///       vaf ∉ [0, 1], purity ∉ (0, 1], or ploidy ≤ 0;
///       AdjustVAFForPurity(π/2, π, 2) = 1 (diploid het round-trip)      (§2.2, §2.4 INV-04, §3.3)
///   • Worked example: CalculateVAF(25, 100) = 0.25                      (§7.1)
///
/// All randomness is LOCALLY seeded (new Random(seed)); no shared static Rng.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public sealed class OncologyVafFuzzTests
{
    // ── Well-formed-VAF assertion helper ─────────────────────────────────────
    // Pins the documented numeric contract on EVERY accepted VAF: finite, never
    // NaN/Infinity (no DivideByZero leak), and inside [0, 1] (INV-01). This is
    // what stops a fuzz test from rubber-stamping a NaN / out-of-range value.
    private static void AssertWellFormedVaf(double vaf)
    {
        double.IsNaN(vaf).Should().BeFalse("VAF must never be NaN (no 0/0)");
        double.IsInfinity(vaf).Should().BeFalse("VAF must be finite");
        vaf.Should().BeInRange(0.0, 1.0, "VAF = alt/total ∈ [0, 1] (§2.2, INV-01)");
    }

    #region ONCO-VAF-001 — Variant allele frequency (positive sanity)

    [Test]
    public void CalculateVAF_KnownCounts_MatchesHandComputedValue()
    {
        // Docs §7.1 worked example: 25/100 = 0.25, exactly.
        CalculateVAF(25, 100).Should().BeApproximately(0.25, 1e-12);

        // ref=70, alt=30 ⇒ total = 100, VAF = 0.30 (hand-computed).
        CalculateVAF(30, 100).Should().BeApproximately(0.30, 1e-12);
    }

    [Test]
    public void CalculateVAF_AltZero_ReturnsZero()
    {
        // No alt support over real coverage ⇒ VAF = 0 (lower boundary).
        CalculateVAF(0, 100).Should().Be(0.0);
    }

    [Test]
    public void CalculateVAF_RefZero_AltPositive_ReturnsOne()
    {
        // ref = 0, alt > 0 ⇒ total = alt ⇒ VAF = 1.0 (upper boundary, alt == total).
        CalculateVAF(100, 100).Should().Be(1.0);
        CalculateVAF(1, 1).Should().Be(1.0);
    }

    #endregion

    #region ONCO-VAF-001 — BE: ref=alt=0 (zero coverage ⇒ VAF=0, no DivideByZero)

    [Test]
    public void CalculateVAF_RefAndAltZero_ReturnsZero_NoDivideByZero()
    {
        // ref = alt = 0 ⇒ totalReads = 0. The hazard is alt/total = 0/0 →
        // DivideByZero / NaN; the documented contract is VAF = 0 (§3.3, §6.1).
        var act = () => CalculateVAF(0, 0);

        act.Should().NotThrow();
        double vaf = act();
        AssertWellFormedVaf(vaf);
        vaf.Should().Be(0.0, "no coverage ⇒ allele absent ⇒ VAF 0 (§6.1)");
    }

    [Test]
    public void CalculateVAFConfidenceInterval_ZeroCoverage_ThrowsArgumentOutOfRange()
    {
        // A Wilson interval is undefined with zero trials ⇒ documented throw
        // (§3.3, §6.1) — NOT a degenerate [0, 0] or a NaN-bearing interval.
        var act = () => CalculateVAFConfidenceInterval(0, 0);

        act.Should().Throw<ArgumentOutOfRangeException>();
    }

    #endregion

    #region ONCO-VAF-001 — BE: alt > total (malformed, VAF must never exceed 1)

    [Test]
    public void CalculateVAF_AltExceedsTotal_ThrowsArgumentOutOfRange_NeverReturnsAboveOne()
    {
        // Malformed counts: alt (50) > total (10). Documented contract (§3.3):
        // ArgumentOutOfRangeException — the result must NEVER be a VAF of 5.0.
        var act = () => CalculateVAF(50, 10);

        act.Should().Throw<ArgumentOutOfRangeException>();
    }

    [Test]
    public void CalculateVAF_AltExceedsTotalByOne_ThrowsArgumentOutOfRange()
    {
        // Tightest boundary: alt = total + 1 must still be rejected, not yield a
        // VAF just above 1.
        var act = () => CalculateVAF(101, 100);

        act.Should().Throw<ArgumentOutOfRangeException>();
    }

    [Test]
    public void CalculateVAFConfidenceInterval_AltExceedsTotal_ThrowsArgumentOutOfRange()
    {
        var act = () => CalculateVAFConfidenceInterval(60, 50);

        act.Should().Throw<ArgumentOutOfRangeException>();
    }

    #endregion

    #region ONCO-VAF-001 — BE: negative counts (documented reject, no negative VAF)

    [Test]
    public void CalculateVAF_NegativeAlt_ThrowsArgumentOutOfRange()
    {
        // -1 boundary: a negative alt count must be rejected, not produce a
        // negative VAF (§3.3).
        var act = () => CalculateVAF(-1, 100);

        act.Should().Throw<ArgumentOutOfRangeException>();
    }

    [Test]
    public void CalculateVAF_NegativeTotal_ThrowsArgumentOutOfRange()
    {
        var act = () => CalculateVAF(0, -5);

        act.Should().Throw<ArgumentOutOfRangeException>();
    }

    [Test]
    public void CalculateVAF_BothNegative_ThrowsArgumentOutOfRange()
    {
        // alt < 0 is checked first; the call must throw regardless.
        var act = () => CalculateVAF(-3, -10);

        act.Should().Throw<ArgumentOutOfRangeException>();
    }

    [Test]
    public void CalculateVAF_IntMinValueCounts_ThrowsArgumentOutOfRange()
    {
        // Extreme negative boundary (Int32.MinValue) must be rejected cleanly,
        // never wrap to a positive value.
        var actAlt = () => CalculateVAF(int.MinValue, 100);
        var actTotal = () => CalculateVAF(0, int.MinValue);

        actAlt.Should().Throw<ArgumentOutOfRangeException>();
        actTotal.Should().Throw<ArgumentOutOfRangeException>();
    }

    #endregion

    #region ONCO-VAF-001 — BE: huge depth (no integer overflow; VAF stays in [0,1])

    [Test]
    public void CalculateVAF_MaxIntAltAndTotal_ReturnsExactlyOne_NoOverflow()
    {
        // Huge depth at the alt == total boundary: Int32.MaxValue / Int32.MaxValue.
        // The alt > total comparison must not overflow (no alt+total sum), and the
        // ratio must be taken in double width ⇒ exactly 1.0, not a wrapped value.
        double vaf = CalculateVAF(int.MaxValue, int.MaxValue);

        AssertWellFormedVaf(vaf);
        vaf.Should().Be(1.0);
    }

    [Test]
    public void CalculateVAF_HugeTotal_ZeroAlt_ReturnsZero_NoOverflow()
    {
        double vaf = CalculateVAF(0, int.MaxValue);

        AssertWellFormedVaf(vaf);
        vaf.Should().Be(0.0);
    }

    [Test]
    public void CalculateVAF_HugeDepth_HalfAlt_ApproximatesOneHalf_NoOverflow()
    {
        // alt = total/2 at huge depth ⇒ ~0.5 with no integer wrap. If alt+total
        // were summed in int32, total ≈ Int32.MaxValue would overflow on any
        // ref+alt reconstruction; the ratio alt/total stays correct because it is
        // computed directly in double.
        int total = int.MaxValue;
        int alt = total / 2;

        double vaf = CalculateVAF(alt, total);

        AssertWellFormedVaf(vaf);
        vaf.Should().BeApproximately(0.5, 1e-6);
    }

    [Test]
    public void CalculateVAF_HugeDepth_AltJustBelowTotal_StaysInRange_NoOverflow()
    {
        // alt = total - 1 at Int32.MaxValue: VAF must be just below 1, never wrap
        // above it.
        int total = int.MaxValue;
        double vaf = CalculateVAF(total - 1, total);

        AssertWellFormedVaf(vaf);
        vaf.Should().BeLessThan(1.0);
        vaf.Should().BeApproximately(1.0, 1e-6);
    }

    [Test]
    public void CalculateVAF_HugeAltExceedsSlightlySmallerTotal_ThrowsArgumentOutOfRange()
    {
        // alt = Int32.MaxValue, total = Int32.MaxValue - 1: alt > total even at
        // the int ceiling — the comparison must hold (no overflow) and reject.
        var act = () => CalculateVAF(int.MaxValue, int.MaxValue - 1);

        act.Should().Throw<ArgumentOutOfRangeException>();
    }

    [Test]
    public void CalculateVAFConfidenceInterval_HugeDepth_WellFormedInterval_NoOverflow()
    {
        // Wilson interval at huge n: n = (double)totalReads must not overflow, the
        // bounds must stay in [0, 1] with lower ≤ vaf-ish center ≤ upper, and the
        // margin must shrink toward 0 as n → ∞ (INV-02, INV-03).
        int total = int.MaxValue;
        int alt = total / 4;

        var ci = CalculateVAFConfidenceInterval(alt, total);

        AssertWellFormedVaf(ci.Vaf);
        ci.Vaf.Should().BeApproximately(0.25, 1e-6);
        ci.Lower.Should().BeInRange(0.0, 1.0);
        ci.Upper.Should().BeInRange(0.0, 1.0);
        ci.Lower.Should().BeLessThanOrEqualTo(ci.Upper);
        // Vast n ⇒ negligible margin: the interval collapses onto the point.
        (ci.Upper - ci.Lower).Should().BeLessThan(1e-3);
    }

    #endregion

    #region ONCO-VAF-001 — BE: broad random fuzz over raw count pairs

    [Test]
    [CancelAfter(30000)]
    public void CalculateVAF_RandomCountPairs_AlwaysGuardedOrInRange()
    {
        // Fuzz: random non-negative count pairs spanning small to huge depths
        // (including the Int32.MaxValue neighbourhood). When alt ≤ total the VAF
        // must be a finite value in [0, 1] (and exactly 0 at zero coverage);
        // when alt > total the call must THROW — never a nonsense VAF > 1.
        var rng = new Random(8801);
        for (int i = 0; i < 20000; i++)
        {
            // Mix ordinary depths with values near Int32.MaxValue to probe overflow.
            int total = rng.Next(2) == 0 ? rng.Next(0, 5000)
                                         : int.MaxValue - rng.Next(0, 5000);
            int alt = rng.Next(2) == 0 ? rng.Next(0, 5000)
                                       : int.MaxValue - rng.Next(0, 5000);

            if (alt <= total)
            {
                double vaf = CalculateVAF(alt, total);
                AssertWellFormedVaf(vaf);
                if (total == 0) vaf.Should().Be(0.0);
                if (alt == total && total > 0) vaf.Should().Be(1.0);
            }
            else
            {
                var act = () => CalculateVAF(alt, total);
                act.Should().Throw<ArgumentOutOfRangeException>();
            }
        }
    }

    [Test]
    [CancelAfter(30000)]
    public void CalculateVAFConfidenceInterval_RandomWellFormed_BoundsAlwaysInRange()
    {
        // INV-02 / INV-03 fuzz: for any covered (total > 0, alt ≤ total) locus the
        // Wilson bounds stay within [0, 1] and ordered lower ≤ upper, and the
        // point VAF matches alt/total.
        var rng = new Random(192837);
        for (int i = 0; i < 8000; i++)
        {
            int total = rng.Next(1, 200000);
            int alt = rng.Next(0, total + 1);

            var ci = CalculateVAFConfidenceInterval(alt, total);

            AssertWellFormedVaf(ci.Vaf);
            ci.Vaf.Should().BeApproximately((double)alt / total, 1e-12);
            ci.Lower.Should().BeInRange(0.0, 1.0, "Wilson interval ⊆ [0, 1] (INV-03)");
            ci.Upper.Should().BeInRange(0.0, 1.0, "Wilson interval ⊆ [0, 1] (INV-03)");
            ci.Lower.Should().BeLessThanOrEqualTo(ci.Upper, "lower ≤ upper (INV-02)");
            ci.Confidence.Should().Be(0.95);
        }
    }

    #endregion

    #region ONCO-VAF-001 — Wilson CI: documented degenerate boundaries (p̂ = 0, p̂ = 1)

    [Test]
    public void CalculateVAFConfidenceInterval_KnownExample_MatchesDocumentedBounds()
    {
        // Docs §7.1 numerical walk-through: 25/100 at 95% ⇒ [0.1754509, 0.3430465].
        var ci = CalculateVAFConfidenceInterval(25, 100);

        ci.Vaf.Should().BeApproximately(0.25, 1e-9);
        ci.Lower.Should().BeApproximately(0.1754509, 1e-6);
        ci.Upper.Should().BeApproximately(0.3430465, 1e-6);
    }

    [Test]
    public void CalculateVAFConfidenceInterval_PHatZero_LowerZero_UpperPositive()
    {
        // p̂ = 0: Wilson lower clamps to 0, upper is strictly positive with
        // non-zero width (§2.5, §6.1).
        var ci = CalculateVAFConfidenceInterval(0, 100);

        ci.Vaf.Should().Be(0.0);
        ci.Lower.Should().Be(0.0);
        ci.Upper.Should().BeGreaterThan(0.0);
        ci.Upper.Should().BeLessThanOrEqualTo(1.0);
    }

    [Test]
    public void CalculateVAFConfidenceInterval_PHatOne_UpperOne_LowerBelowOne()
    {
        // p̂ = 1: Wilson upper clamps to 1, lower is strictly below 1 (§6.1).
        var ci = CalculateVAFConfidenceInterval(100, 100);

        ci.Vaf.Should().Be(1.0);
        // Upper is clamped to ≤ 1 and sits at the boundary (within fp drift, §5.2/§6.1).
        ci.Upper.Should().BeLessThanOrEqualTo(1.0);
        ci.Upper.Should().BeApproximately(1.0, 1e-12);
        ci.Lower.Should().BeLessThan(1.0);
        ci.Lower.Should().BeGreaterThanOrEqualTo(0.0);
    }

    [Test]
    public void CalculateVAFConfidenceInterval_UnsupportedConfidence_ThrowsArgumentOutOfRange()
    {
        // Only 0.95 (z = 1.96) is supported; other levels / out-of-(0,1) ⇒ throw
        // (§3.3, §5.3 "Intentionally simplified").
        var actOther = () => CalculateVAFConfidenceInterval(25, 100, confidence: 0.99);
        var actHigh = () => CalculateVAFConfidenceInterval(25, 100, confidence: 1.5);
        var actLow = () => CalculateVAFConfidenceInterval(25, 100, confidence: 0.0);

        actOther.Should().Throw<ArgumentOutOfRangeException>();
        actHigh.Should().Throw<ArgumentOutOfRangeException>();
        actLow.Should().Throw<ArgumentOutOfRangeException>();
    }

    #endregion

    #region ONCO-VAF-001 — Purity correction: INV-04 round-trip and validation

    [Test]
    public void AdjustVAFForPurity_DiploidHetRoundTrip_ReturnsOne()
    {
        // INV-04 / §7.1: AdjustVAFForPurity(π/2, π, 2) = 1 for a clonal het SNV in
        // a diploid segment. Observed VAF 0.4 at purity 0.8 ⇒ 1.0.
        AdjustVAFForPurity(0.40, 0.80, 2).Should().BeApproximately(1.0, 1e-12);
    }

    [Test]
    public void AdjustVAFForPurity_VafZero_ReturnsZero()
    {
        // Lower boundary: 0 observed mutant fraction ⇒ 0 corrected fraction.
        AdjustVAFForPurity(0.0, 0.80, 2).Should().Be(0.0);
    }

    [Test]
    public void AdjustVAFForPurity_InvalidArguments_ThrowArgumentOutOfRange()
    {
        // §3.3: vaf ∉ [0, 1], purity ∉ (0, 1] (divides by purity), ploidy ≤ 0.
        var actVafHigh = () => AdjustVAFForPurity(1.5, 0.80, 2);
        var actVafNeg = () => AdjustVAFForPurity(-0.1, 0.80, 2);
        var actPurityZero = () => AdjustVAFForPurity(0.4, 0.0, 2);   // no divide-by-zero leak
        var actPurityHigh = () => AdjustVAFForPurity(0.4, 1.5, 2);
        var actPloidyZero = () => AdjustVAFForPurity(0.4, 0.80, 0);
        var actPloidyNeg = () => AdjustVAFForPurity(0.4, 0.80, -2);

        actVafHigh.Should().Throw<ArgumentOutOfRangeException>();
        actVafNeg.Should().Throw<ArgumentOutOfRangeException>();
        actPurityZero.Should().Throw<ArgumentOutOfRangeException>();
        actPurityHigh.Should().Throw<ArgumentOutOfRangeException>();
        actPloidyZero.Should().Throw<ArgumentOutOfRangeException>();
        actPloidyNeg.Should().Throw<ArgumentOutOfRangeException>();
    }

    [Test]
    public void AdjustVAFForPurity_NaNArguments_ThrowArgumentOutOfRange()
    {
        // NaN inputs must be rejected, not propagate a NaN result.
        var actVaf = () => AdjustVAFForPurity(double.NaN, 0.80, 2);
        var actPurity = () => AdjustVAFForPurity(0.4, double.NaN, 2);
        var actPloidy = () => AdjustVAFForPurity(0.4, 0.80, double.NaN);

        actVaf.Should().Throw<ArgumentOutOfRangeException>();
        actPurity.Should().Throw<ArgumentOutOfRangeException>();
        actPloidy.Should().Throw<ArgumentOutOfRangeException>();
    }

    [Test]
    [CancelAfter(20000)]
    public void AdjustVAFForPurity_RandomValidInputs_FiniteAndNonNegative()
    {
        // Fuzz the valid domain: vaf ∈ [0, 1], purity ∈ (0, 1], ploidy > 0. The
        // corrected mutant fraction must always be finite and non-negative
        // (ratio of non-negatives over a positive purity), never NaN/Infinity.
        var rng = new Random(556677);
        for (int i = 0; i < 10000; i++)
        {
            double vaf = rng.NextDouble();
            double purity = rng.NextDouble() * (1.0 - 1e-6) + 1e-6; // (0, 1]
            double ploidy = rng.NextDouble() * 10.0 + 1e-6;          // > 0

            double adjusted = AdjustVAFForPurity(vaf, purity, ploidy);

            double.IsNaN(adjusted).Should().BeFalse();
            double.IsInfinity(adjusted).Should().BeFalse();
            adjusted.Should().BeGreaterThanOrEqualTo(0.0);
        }
    }

    #endregion
}
