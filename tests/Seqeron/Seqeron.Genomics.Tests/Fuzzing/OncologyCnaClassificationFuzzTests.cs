using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Oncology;
using static Seqeron.Genomics.Oncology.OncologyAnalyzer;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Fuzz tests for the Oncology copy-number-alteration classification area —
/// ONCO-CNA-001. The unit under test is the log2-ratio → copy-number layer in
/// src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs:
///   • <see cref="OncologyAnalyzer.Log2RatioToCopyNumber"/> — continuous absolute
///     copy number n = ploidy·2^log2 (CNVkit <c>_log2_ratio_to_absolute_pure</c>);
///   • <see cref="OncologyAnalyzer.CallCopyNumber"/> — CNVkit hard-threshold
///     integer copy number (<c>absolute_threshold</c>);
///   • <see cref="OncologyAnalyzer.ClassifyCopyNumber"/> — full
///     <see cref="OncologyAnalyzer.CopyNumberCall"/> (absolute + integer + state);
///   • <see cref="OncologyAnalyzer.ClassifyCopyNumbers"/> — per-region batch
///     (length/order preserving).
///
/// This is the foundational ratio/classification member of the CNA family
/// (rows 103–105). ONCO-CNA-002 (row 104, focal amplification) and ONCO-CNA-003
/// (row 105, homozygous deletion) build on these primitives; this file scopes
/// strictly to CNA-001 — the log2(tumor/normal) → CNA-state map. Segmentation
/// (CBS) is explicitly NOT part of this unit (docs §5.3 "Not implemented" defers
/// segmentation to StructuralVariantAnalyzer.SegmentCopyNumber / SV-CNV-001), so
/// the "single bin" target here is exercised against the per-region batch
/// classifier — a length-1 profile must classify without crash and preserve
/// length 1, the degenerate analogue of a single-bin segmentation input.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate / boundary / malformed inputs to a unit and asserts
/// that the code NEVER fails in an undisciplined way: no hang, no nonsense
/// output, no *unhandled* runtime exception (DivideByZero / Overflow / a NaN or
/// ±∞ silently leaking into the downstream classification). Every input must
/// resolve to EITHER a well-defined, theory-correct value OR a *documented,
/// intentional* outcome (an ArgumentException for malformed thresholds,
/// ArgumentOutOfRangeException for non-positive ploidy, ArgumentNullException for
/// a null batch). The headline hazards for the log2→CN map are:
///   • log2 = +∞ / a huge finite log2 → 2^log2 = +∞: the integer copy number
///     (the downstream call) must NOT become a wrapped/garbage value — it stays
///     the saturated Math.Ceiling result (Int32.MaxValue) and the state stays the
///     mathematically-correct Amplification; no ±∞ leaks into the INTEGER call or
///     the STATE that segmentation/calling consumes;
///   • log2 = −∞ → 2^log2 = 0: must give the lowest finite state (DeepDeletion,
///     CN 0, absolute 0) — no NaN, no negative CN;
///   • NaN log2 → the DOCUMENTED no-call: Neutral, CN = rounded ploidy, absolute
///     = ploidy (NOT a NaN propagating into a segment mean / CN) (§3.3, §6.1);
///   • a single-bin batch → a length-1 result (no DivideByZero / variance-of-one
///     crash that a naive segmentation would hit on length 1), and an empty batch
///     → empty result (§6.1).
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: ONCO-CNA-001 — Copy-number alteration classification (Oncology)
/// Checklist: docs/checklists/03_FUZZING.md, row 103.
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — граничні значення: 0, -1, MaxInt, empty.
///     Targets (checklist row 103): "log2=±∞, NaN ratio, single bin".
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The documented contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// Copy_Number_Alteration_Classification.md
/// (docs/algorithms/Oncology/Copy_Number_Alteration_Classification.md):
///   • absolute n = ploidy·2^log2; v=0⇒2, v=1⇒4, v=−1⇒1                  (§2.2, INV-01)
///   • hard-threshold integer CN = index of first ascending cutoff with
///     log2 ≤ cutoff; above the last cutoff CN = ⌈ploidy·2^log2⌉ (≥4)      (§2.2, §4)
///   • default tumor cutoffs (−1.1, −0.25, 0.2, 0.7) partition the log2 axis
///     into states 0/1/2/3/≥4                                             (§2.2, §4.2)
///   • CN→state: 0→DeepDeletion, 1→Loss, 2→Neutral, 3→Gain, ≥4→Amplification (INV-04)
///   • boundary comparison is inclusive (log2 ≤ cutoff ⇒ lower state)      (§3.3, §6.1)
///   • NaN log2 ⇒ no-call ⇒ Neutral, CN = rounded ploidy, absolute=ploidy  (§3.3, §6.1)
///   • integer CN is non-decreasing in log2 (INV-02) and ≥ 0 for finite log2 (INV-03)
///   • ClassifyCopyNumbers: output length = input length, element i ↔ input i (INV-05)
///   • thresholds must be exactly 4 strictly ascending non-NaN values ⇒ else
///     ArgumentException; ploidy > 0 ⇒ else ArgumentOutOfRangeException;
///     batch not null ⇒ else ArgumentNullException                        (§3.3)
///   • worked example: ClassifyCopyNumber(1.0) ⇒ abs 4, CN 4, Amplification (§7.1)
///
/// All randomness is LOCALLY seeded (new Random(seed)); no shared static Rng.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public sealed class OncologyCnaClassificationFuzzTests
{
    // ── Well-formed-call assertion helper ────────────────────────────────────
    // Pins the documented downstream contract on EVERY accepted classification:
    // the INTEGER copy number (what calling/segmentation consume) is finite and
    // ≥ 0, and the STATE is a defined enum member — i.e. no NaN/±∞ ever leaks
    // into the integer call or the state, even when the continuous absolute CN is
    // an IEEE ±∞/0 by construction (2^±∞). This is what stops a fuzz test from
    // rubber-stamping a wrapped integer or an out-of-band state.
    private static void AssertWellFormedCall(in CopyNumberCall call)
    {
        call.IntegerCopyNumber.Should().BeGreaterThanOrEqualTo(
            0, "integer CN is ≥ 0 (DeepDeletion is the floor; INV-03)");
        Enum.IsDefined(typeof(CopyNumberState), call.State).Should().BeTrue(
            "state must be one of the five defined CNA states (INV-04)");
        // The integer call must equal the state's class (the CN↔state map, INV-04).
        StateForCn(call.IntegerCopyNumber).Should().Be(
            call.State, "state must be derived from the integer CN (INV-04)");
    }

    // Documented CN→state map (INV-04), mirrored here so the test owns the
    // expected mapping independently of the implementation's private helper.
    private static CopyNumberState StateForCn(int cn) => cn switch
    {
        <= 0 => CopyNumberState.DeepDeletion, // CN 0 (negative CN is impossible for finite log2)
        1 => CopyNumberState.Loss,
        2 => CopyNumberState.Neutral,
        3 => CopyNumberState.Gain,
        _ => CopyNumberState.Amplification    // CN ≥ 4
    };

    #region ONCO-CNA-001 — Positive sanity (documented gain / loss / neutral values)

    [Test]
    public void ClassifyCopyNumber_DocumentedAmplificationExample_MatchesExactValues()
    {
        // Docs §7.1 worked example: tumor = 2× normal ⇒ log2(2) = +1.0 ⇒
        // absolute n = 2·2^1 = 4, and 1.0 exceeds every default cutoff so
        // CN = ⌈4⌉ = 4 ⇒ Amplification (NOT a single-copy Gain — log2 1.0 > 0.7).
        var call = ClassifyCopyNumber(1.0);

        call.AbsoluteCopyNumber.Should().BeApproximately(4.0, 1e-12);
        call.IntegerCopyNumber.Should().Be(4);
        call.State.Should().Be(CopyNumberState.Amplification);
        AssertWellFormedCall(call);
    }

    [Test]
    public void ClassifyCopyNumber_EqualCoverage_IsNeutralWithTwoCopies()
    {
        // Equal tumor/normal coverage ⇒ log2(1) = 0 ⇒ neutral diploid:
        // absolute = 2·2^0 = 2, 0 ≤ 0.2 (cutoff index 2) ⇒ CN 2 ⇒ Neutral.
        var call = ClassifyCopyNumber(0.0);

        call.AbsoluteCopyNumber.Should().BeApproximately(2.0, 1e-12);
        call.IntegerCopyNumber.Should().Be(2);
        call.State.Should().Be(CopyNumberState.Neutral);
        AssertWellFormedCall(call);
    }

    [Test]
    public void ClassifyCopyNumber_HalfCoverage_IsSingleCopyLoss()
    {
        // tumor = 0.5× normal ⇒ log2(0.5) = −1.0 ⇒ absolute = 2·2^-1 = 1,
        // and −1.1 < −1.0 ≤ −0.25 (cutoff index 1) ⇒ CN 1 ⇒ Loss.
        var call = ClassifyCopyNumber(-1.0);

        call.AbsoluteCopyNumber.Should().BeApproximately(1.0, 1e-12);
        call.IntegerCopyNumber.Should().Be(1);
        call.State.Should().Be(CopyNumberState.Loss);
        AssertWellFormedCall(call);
    }

    [Test]
    public void ClassifyCopyNumber_SingleCopyGain_IsGainNotAmplification()
    {
        // A true single-copy gain (CN 3) is log2(3/2) = 0.585, which sits in
        // 0.2 < log2 ≤ 0.7 (cutoff index 3) ⇒ CN 3 ⇒ Gain (distinct from the
        // log2 1.0 amplification above — guards the upper-bin boundary).
        double log2Gain = Math.Log2(3.0 / 2.0); // ≈ 0.585

        var call = ClassifyCopyNumber(log2Gain);

        call.IntegerCopyNumber.Should().Be(3);
        call.State.Should().Be(CopyNumberState.Gain);
        AssertWellFormedCall(call);
    }

    [Test]
    public void ClassifyCopyNumber_DeepDeletion_IsCnZero()
    {
        // log2 ≤ −1.1 (cutoff index 0) ⇒ CN 0 ⇒ DeepDeletion.
        var call = ClassifyCopyNumber(-2.0);

        call.IntegerCopyNumber.Should().Be(0);
        call.State.Should().Be(CopyNumberState.DeepDeletion);
        AssertWellFormedCall(call);
    }

    #endregion

    #region ONCO-CNA-001 — BE: log2 = +∞ / huge log2 (no garbage CN, sensible state)

    [Test]
    public void CallCopyNumber_PositiveInfinityLog2_SaturatesNonNegative_NoIntWrap()
    {
        // log2 = +∞ ⇒ 2^log2 = +∞. The HAZARD is the (int)Math.Ceiling(+∞) cast
        // wrapping to Int32.MinValue (a negative CN) and mis-classifying as a
        // low state. The disciplined contract: the integer call must remain
        // ≥ 0 (here saturated to Int32.MaxValue) ⇒ Amplification — never a
        // negative/garbage CN leaking into the call.
        var act = () => CallCopyNumber(double.PositiveInfinity);

        act.Should().NotThrow();
        int cn = act();
        cn.Should().BeGreaterThanOrEqualTo(AmplificationFloorForTest,
            "+∞ log2 is an unbounded gain ⇒ at least the amplification class, never wrapped negative");
    }

    [Test]
    public void ClassifyCopyNumber_PositiveInfinityLog2_IsAmplification_NoGarbageCall()
    {
        // The full call for +∞: the continuous absolute CN is the IEEE value +∞
        // by construction (2^+∞), but the INTEGER CN and STATE — the products
        // calling/segmentation consume — must stay disciplined: a non-negative
        // saturated CN and the mathematically-correct Amplification state.
        var call = ClassifyCopyNumber(double.PositiveInfinity);

        AssertWellFormedCall(call);
        call.State.Should().Be(CopyNumberState.Amplification);
        double.IsNaN(call.IntegerCopyNumber).Should().BeFalse();
    }

    [Test]
    public void ClassifyCopyNumber_HugeFiniteLog2_OverflowsToAmplification_NoNegativeCn()
    {
        // A huge FINITE log2 (1024) makes 2·2^1024 overflow to +∞ in double.
        // Even so the integer call must saturate non-negative (Int32.MaxValue)
        // and classify as Amplification — INV-03 (CN ≥ 0 for finite log2) holds.
        var call = ClassifyCopyNumber(1024.0);

        AssertWellFormedCall(call);
        call.IntegerCopyNumber.Should().BeGreaterThanOrEqualTo(AmplificationFloorForTest);
        call.State.Should().Be(CopyNumberState.Amplification);
    }

    #endregion

    #region ONCO-CNA-001 — BE: log2 = −∞ (lowest finite state, no NaN)

    [Test]
    public void ClassifyCopyNumber_NegativeInfinityLog2_IsDeepDeletion_AbsoluteZero()
    {
        // log2 = −∞ ⇒ 2^log2 = 0 ⇒ absolute = 0 (finite, not NaN), and −∞ ≤ −1.1
        // (first cutoff) ⇒ CN 0 ⇒ DeepDeletion. No −∞ leaks: 2^−∞ underflows to a
        // clean 0, the floor of the state ladder.
        var call = ClassifyCopyNumber(double.NegativeInfinity);

        call.AbsoluteCopyNumber.Should().Be(0.0);
        double.IsNaN(call.AbsoluteCopyNumber).Should().BeFalse();
        double.IsInfinity(call.AbsoluteCopyNumber).Should().BeFalse();
        call.IntegerCopyNumber.Should().Be(0);
        call.State.Should().Be(CopyNumberState.DeepDeletion);
        AssertWellFormedCall(call);
    }

    #endregion

    #region ONCO-CNA-001 — BE: NaN ratio (documented Neutral no-call, no NaN leak)

    [Test]
    public void ClassifyCopyNumber_NaNLog2_IsNeutralNoCall_AbsoluteEqualsPloidy()
    {
        // Docs §3.3/§6.1: a NaN log2 is a no-call ⇒ Neutral, CN = rounded ploidy
        // (= 2 diploid), absolute = ploidy. The hazard is a NaN propagating into
        // a downstream segment mean / CN; the contract pins a clean Neutral 2.
        var call = ClassifyCopyNumber(double.NaN);

        call.IntegerCopyNumber.Should().Be(2);
        call.State.Should().Be(CopyNumberState.Neutral);
        call.AbsoluteCopyNumber.Should().Be(2.0);
        double.IsNaN(call.AbsoluteCopyNumber).Should().BeFalse(
            "the no-call absolute CN is the reference ploidy, never a propagated NaN");
        AssertWellFormedCall(call);
    }

    [Test]
    public void CallCopyNumber_NaNLog2_ReturnsRoundedPloidy_NoComparisonAgainstNaN()
    {
        // A naive cutoff scan (log2 ≤ cutoff) against NaN is always false and
        // would fall through to ⌈ploidy·2^NaN⌉ = ⌈NaN⌉ = a wrapped int. The
        // documented guard short-circuits NaN to the rounded ploidy BEFORE the
        // scan. Verify with a non-default ploidy that the no-call tracks ploidy.
        CallCopyNumber(double.NaN).Should().Be(2);              // default diploid
        CallCopyNumber(double.NaN, ploidy: 3.0).Should().Be(3); // rounded ploidy
        CallCopyNumber(double.NaN, ploidy: 4.4).Should().Be(4); // round(4.4) = 4
    }

    [Test]
    public void ClassifyCopyNumbers_NaNAmongRealRatios_DoesNotPoisonNeighbours()
    {
        // A NaN no-call in the middle of a batch must be classified independently
        // as Neutral and must NOT corrupt the real calls around it (per-element
        // map, INV-05) — no NaN bleeding into a segment-mean-style aggregate.
        var input = new[] { -2.0, double.NaN, 1.0 };

        var calls = ClassifyCopyNumbers(input);

        calls.Should().HaveCount(3);
        calls[0].State.Should().Be(CopyNumberState.DeepDeletion);
        calls[1].State.Should().Be(CopyNumberState.Neutral); // the NaN no-call
        calls[1].IntegerCopyNumber.Should().Be(2);
        calls[2].State.Should().Be(CopyNumberState.Amplification);
        foreach (var c in calls) AssertWellFormedCall(c);
    }

    #endregion

    #region ONCO-CNA-001 — BE: single bin / empty batch (length-preserving, no crash)

    [Test]
    public void ClassifyCopyNumbers_SingleBin_ReturnsLengthOne_NoSegmentationCrash()
    {
        // "Single bin": a profile of exactly one log2 ratio. A naive segmenter
        // would risk a DivideByZero / variance-of-one crash; the classification
        // layer must return exactly one call, classified correctly (INV-05).
        var calls = ClassifyCopyNumbers(new[] { 1.0 });

        calls.Should().HaveCount(1);
        calls[0].IntegerCopyNumber.Should().Be(4);
        calls[0].State.Should().Be(CopyNumberState.Amplification);
        AssertWellFormedCall(calls[0]);
    }

    [Test]
    public void ClassifyCopyNumbers_SingleNeutralBin_IsNeutral()
    {
        // Single bin at the neutral baseline ⇒ a single Neutral call.
        var calls = ClassifyCopyNumbers(new[] { 0.0 });

        calls.Should().ContainSingle();
        calls[0].State.Should().Be(CopyNumberState.Neutral);
        calls[0].AbsoluteCopyNumber.Should().BeApproximately(2.0, 1e-12);
    }

    [Test]
    public void ClassifyCopyNumbers_EmptyBatch_ReturnsEmpty_NoThrow()
    {
        // Empty profile ⇒ empty result (§6.1), the zero-bin boundary — no throw.
        var calls = ClassifyCopyNumbers(Array.Empty<double>());

        calls.Should().BeEmpty();
    }

    [Test]
    public void ClassifyCopyNumbers_NullBatch_ThrowsArgumentNull()
    {
        // Documented null guard (§3.3) — a null enumerable is rejected cleanly,
        // not a NullReference deep in the loop.
        var act = () => ClassifyCopyNumbers(null!);

        act.Should().Throw<ArgumentNullException>();
    }

    [Test]
    public void ClassifyCopyNumbers_TwoSegmentProfile_SplitsAtTheBoundary()
    {
        // A clean two-segment profile: a run of neutral bins (log2 0) followed by
        // a run of amplified bins (log2 1.0). The per-bin classification must flip
        // exactly at the boundary index — the analogue of a clean segment break.
        var input = new[] { 0.0, 0.0, 0.0, 1.0, 1.0 };

        var calls = ClassifyCopyNumbers(input);

        calls.Should().HaveCount(5);
        calls.Take(3).Should().OnlyContain(c => c.State == CopyNumberState.Neutral);
        calls.Skip(3).Should().OnlyContain(c => c.State == CopyNumberState.Amplification);
    }

    #endregion

    #region ONCO-CNA-001 — BE: threshold / ploidy validation (documented throws)

    [Test]
    public void CallCopyNumber_NonAscendingThresholds_ThrowsArgument()
    {
        // §3.3: thresholds must be strictly ascending — a flat/descending list
        // does not partition the log2 axis ⇒ ArgumentException.
        var descending = new[] { 0.7, 0.2, -0.25, -1.1 };
        var act = () => CallCopyNumber(0.0, descending);

        act.Should().Throw<ArgumentException>();
    }

    [Test]
    public void CallCopyNumber_WrongThresholdCount_ThrowsArgument()
    {
        // Exactly four cutoffs are required (five states) — a 3- or 5-element list
        // is rejected (§3.3).
        var three = new[] { -1.0, 0.0, 1.0 };
        var five = new[] { -2.0, -1.0, 0.0, 1.0, 2.0 };

        ((Action)(() => CallCopyNumber(0.0, three))).Should().Throw<ArgumentException>();
        ((Action)(() => CallCopyNumber(0.0, five))).Should().Throw<ArgumentException>();
    }

    [Test]
    public void CallCopyNumber_NaNInThresholds_ThrowsArgument()
    {
        // A NaN cutoff cannot order the axis ⇒ ArgumentException (§3.3) — not a
        // silent mis-call against an unorderable boundary.
        var withNaN = new[] { -1.1, double.NaN, 0.2, 0.7 };
        var act = () => CallCopyNumber(0.0, withNaN);

        act.Should().Throw<ArgumentException>();
    }

    [Test]
    public void Log2RatioToCopyNumber_NonPositiveOrNaNPloidy_ThrowsArgumentOutOfRange()
    {
        // §3.3: ploidy must be positive (and finite-defined) — 0, negative, and
        // NaN are rejected (a non-positive ploidy would invert/zero the scale).
        ((Action)(() => Log2RatioToCopyNumber(0.0, 0.0))).Should().Throw<ArgumentOutOfRangeException>();
        ((Action)(() => Log2RatioToCopyNumber(0.0, -2.0))).Should().Throw<ArgumentOutOfRangeException>();
        ((Action)(() => Log2RatioToCopyNumber(0.0, double.NaN))).Should().Throw<ArgumentOutOfRangeException>();
    }

    [Test]
    public void ClassifyCopyNumber_NonPositivePloidy_ThrowsArgumentOutOfRange()
    {
        var act = () => ClassifyCopyNumber(0.0, ploidy: -1.0);

        act.Should().Throw<ArgumentOutOfRangeException>();
    }

    #endregion

    #region ONCO-CNA-001 — Boundary: inclusive cutoff edges (lower state of the bin)

    [Test]
    public void CallCopyNumber_ExactlyOnCutoffs_AssignsLowerStateOfBin()
    {
        // §3.3/§6.1: comparison is inclusive (log2 ≤ cutoff), so a value exactly
        // on cutoff index i is assigned CN i (the lower bin). Pin each default
        // cutoff exactly.
        CallCopyNumber(-1.1).Should().Be(0);  // ≤ −1.1 ⇒ DeepDeletion
        CallCopyNumber(-0.25).Should().Be(1); // ≤ −0.25 ⇒ Loss
        CallCopyNumber(0.2).Should().Be(2);   // ≤ 0.2  ⇒ Neutral
        CallCopyNumber(0.7).Should().Be(3);   // ≤ 0.7  ⇒ Gain
    }

    #endregion

    #region ONCO-CNA-001 — BE: broad random fuzz (monotone, finite, non-negative)

    [Test]
    [CancelAfter(30000)]
    public void ClassifyCopyNumber_RandomFiniteLog2_AlwaysWellFormed(
        [Values(20260103, 777, 424242)] int seed)
    {
        // Fuzz finite log2 ratios across a wide range (including extreme but
        // finite magnitudes). Every call must be well-formed: CN ≥ 0, a defined
        // state, state derived from CN. For finite log2 the absolute CN must be
        // finite and non-negative — no NaN/±∞ leaking from 2^log2 within the
        // finite domain, and no negative CN (INV-03).
        var rng = new Random(seed);
        for (int i = 0; i < 20000; i++)
        {
            // Span ordinary ratios and large-but-finite magnitudes (±50).
            double log2 = (rng.NextDouble() * 2.0 - 1.0) * 50.0;

            var call = ClassifyCopyNumber(log2);

            AssertWellFormedCall(call);
            double.IsNaN(call.AbsoluteCopyNumber).Should().BeFalse(
                "finite log2 ⇒ finite-or-overflow absolute CN, never NaN");
            call.AbsoluteCopyNumber.Should().BeGreaterThanOrEqualTo(
                0.0, "absolute CN = ploidy·2^log2 ≥ 0");
        }
    }

    [Test]
    [CancelAfter(30000)]
    public void CallCopyNumber_IntegerCnIsNonDecreasingInLog2(
        [Values(13, 99, 2026)] int seed)
    {
        // INV-02: integer CN is non-decreasing in log2 (ascending cutoffs + a
        // monotone ceiling else-branch). Fuzz ordered pairs and assert the
        // monotonicity holds — a sorting/comparison bug would surface here.
        var rng = new Random(seed);
        for (int i = 0; i < 20000; i++)
        {
            double a = (rng.NextDouble() * 2.0 - 1.0) * 10.0;
            double b = (rng.NextDouble() * 2.0 - 1.0) * 10.0;
            double lo = Math.Min(a, b);
            double hi = Math.Max(a, b);

            int cnLo = CallCopyNumber(lo);
            int cnHi = CallCopyNumber(hi);

            cnHi.Should().BeGreaterThanOrEqualTo(
                cnLo, "integer CN is non-decreasing in log2 (INV-02)");
        }
    }

    [Test]
    [CancelAfter(30000)]
    public void ClassifyCopyNumbers_RandomBatch_PreservesLengthAndOrder()
    {
        // INV-05: the batch is a length/order-preserving per-element map. Fuzz
        // batches of mixed length (including 0 and 1) with NaN no-calls sprinkled
        // in, and assert element-wise agreement with the single-region classifier.
        var rng = new Random(31337);
        for (int t = 0; t < 2000; t++)
        {
            int n = rng.Next(0, 8); // includes empty and single-bin
            var input = new double[n];
            for (int i = 0; i < n; i++)
            {
                input[i] = rng.Next(10) == 0
                    ? double.NaN // occasional no-call
                    : (rng.NextDouble() * 2.0 - 1.0) * 5.0;
            }

            var calls = ClassifyCopyNumbers(input);

            calls.Should().HaveCount(n);
            for (int i = 0; i < n; i++)
            {
                var single = ClassifyCopyNumber(input[i]);
                calls[i].IntegerCopyNumber.Should().Be(single.IntegerCopyNumber);
                calls[i].State.Should().Be(single.State);
                AssertWellFormedCall(calls[i]);
            }
        }
    }

    #endregion

    // CN floor of the amplification class (CNVkit AMP(4) ≥ +0.7), mirrored locally
    // so the test owns the expected threshold independently of the source.
    private const int AmplificationFloorForTest = 4;
}
