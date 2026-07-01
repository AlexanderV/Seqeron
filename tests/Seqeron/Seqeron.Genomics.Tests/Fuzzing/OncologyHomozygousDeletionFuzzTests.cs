using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Oncology;
using static Seqeron.Genomics.Oncology.OncologyAnalyzer;

namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the Oncology homozygous (deep) deletion calling area —
/// ONCO-CNA-003. The unit under test is the CN-0 / Deep-Deletion filter in
/// src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs:
///   • <see cref="OncologyAnalyzer.IsHomozygousDeletion"/> — single-segment
///     predicate (integer copy number == 0);
///   • <see cref="OncologyAnalyzer.DetectHomozygousDeletions"/> — order-preserving
///     filter selecting the CN-0 segments from a segment stream;
///   • <see cref="OncologyAnalyzer.IdentifyDeletedTumorSuppressors"/> —
///     arm → tumour-suppressor panel map over the homozygous deletions.
///
/// This is the homozygous-deletion member of the CNA family (rows 103–105). It is
/// defined purely on top of the ONCO-CNA-001 integer copy-number call
/// (<see cref="OncologyAnalyzer.CallCopyNumber"/>, CNVkit <c>absolute_threshold</c>)
/// and the ONCO-CNA-002 <see cref="OncologyAnalyzer.CopyNumberArmSegment"/> /
/// <c>ValidateArmSegment</c> (docs §5.2). No new copy-number threshold is
/// introduced: a segment is a homozygous deletion <b>iff</b> its hard-threshold
/// integer copy number is 0 — i.e. its mean log2 ratio is ≤ the first (deepest)
/// cutoff, −1.1 by default (docs §2.2, INV-01). Segmentation itself is upstream
/// (SV-CNV-001), so the BE "single-bin deletion" target is exercised as a width-1
/// arm segment whose log2 still classifies to CN 0 — the degenerate analogue of a
/// one-bin homozygous-deletion call.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate / boundary / malformed inputs to a unit and asserts
/// that the code NEVER fails in an undisciplined way: no hang, no nonsense output,
/// no *unhandled* runtime exception (a DivideByZero on a width-1 / zero-arm-length
/// segment, a NaN/±∞ leaking through the integer-CN comparison as a false or
/// missed deep-deletion call, or an off-by-one at the CN-0 / CN-1 boundary). Every
/// input must resolve to EITHER a well-defined, theory-correct outcome OR a
/// *documented, intentional* one (ArgumentNullException for a null stream,
/// ArgumentException for ArmLength ≤ 0 / End ≤ Start, ArgumentOutOfRangeException
/// for non-positive ploidy). The headline hazards for the CN-0 rule are:
///   • the CN-0 / CN-1 boundary: log2 EXACTLY −1.1 is a homozygous deletion (the
///     CNVkit cutoff is INCLUSIVE, "≤ each threshold in sequence"), and log2 just
///     ABOVE −1.1 is a single-copy loss (CN 1), NOT homozygous — no off-by-one
///     (§6.1);
///   • a true copy-number-0 region has log2 = log2(0/2) = −∞; −∞ ≤ −1.1 must call
///     CN 0 (DeepDeletion) WITHOUT a crash or a wrapped value — log2(0) = −∞ must
///     not leak unguarded into the decision as anything but the deepest state
///     (§6.1 "log2 exactly −1.1", extended to the limiting CN-0 value);
///   • a single-bin (width-1) deep-deletion segment must classify without a
///     DivideByZero / length crash and, if its log2 ≤ −1.1, be reported as a
///     homozygous deletion (ArmLength &gt; 0 is validated up front);
///   • a NaN log2 is the documented CNVkit no-call ⇒ neutral reference CN
///     (rounded ploidy = 2) ⇒ NOT reported — a NaN must never slip through the
///     `== 0` test as a false deep deletion (§3.3, §6.1);
///   • a single-copy (heterozygous) loss (CN 1, e.g. log2 −0.5) is NEVER reported
///     (INV-02) — the filter does not over-call shallow losses as deep ones.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: ONCO-CNA-003 — Homozygous (Deep) Deletion Detection (Oncology)
/// Checklist: docs/checklists/03_FUZZING.md, row 105.
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — граничні значення: 0, -1, MaxInt, empty.
///     Targets (checklist row 105): "CN exactly 0, near-0, single-bin deletion".
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The documented contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// Homozygous_Deletion_Detection.md
/// (docs/algorithms/Oncology/Homozygous_Deletion_Detection.md):
///   • Homozygous deletion ⇔ integer copy number == 0 (DeepDeletion)      (§2.2, INV-01)
///   • integer CN == 0 ⇔ log2 ≤ first cutoff (−1.1, INCLUSIVE)            (§2.2, §6.1)
///   • log2 EXACTLY −1.1 ⇒ reported (CN 0)                                (§6.1)
///   • log2 just above −1.1 ⇒ NOT reported (CN 1, single-copy loss)       (§6.1, INV-02)
///   • single-copy (heterozygous) loss (CN 1) is never reported          (INV-02)
///   • NaN log2 ⇒ neutral no-call (CN = rounded ploidy) ⇒ NOT reported    (§3.3, §6.1)
///   • output is a subset of the input in input order (a filter)         (INV-03)
///   • a tumour suppressor is reported only for an arm carrying a hom-del (INV-04)
///   • arm matching is case-insensitive, panel order, each gene once     (INV-04)
///   • null segments/deletions ⇒ ArgumentNullException                   (§3.3, §6.1)
///   • empty input ⇒ empty result                                        (§6.1)
///   • ArmLength ≤ 0 or End ≤ Start ⇒ ArgumentException                   (§3.3)
///   • non-positive ploidy ⇒ ArgumentOutOfRangeException                 (§3.3)
///   • worked example: 9p [0,1000) log2 −2.0 ⇒ CN 0 ⇒ CDKN2A; 3p log2
///     −0.5 ⇒ CN 1 ⇒ not a homozygous deletion                          (§7.1)
///
/// All randomness is LOCALLY seeded (new Random(seed)); no shared static Rng.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public sealed class OncologyHomozygousDeletionFuzzTests
{
    // Fixed reference arm length (1 Mb). Only Log2Ratio drives the CN-0 decision;
    // Start/End/ArmLength matter only for the upstream validation guards.
    private const long ArmLen = 1_000_000L;

    // The deepest (first) CNVkit cutoff: log2 ≤ this ⇒ integer CN 0 ⇒ homozygous
    // deletion. Default is −1.1 (docs §2.2). Read from the source default list so
    // the test tracks the documented constant rather than a literal copy.
    private static readonly double DeepDeletionLog2Cutoff = DefaultCopyNumberThresholds[0];

    private static CopyNumberArmSegment Seg(
        string arm, double log2, long start = 0, long end = 100_000, long armLength = ArmLen) =>
        new(arm, start, end, armLength, log2);

    // ── Well-formed-result assertion helper ──────────────────────────────────
    // Pins the documented filter contract on EVERY accepted detection: the result
    // is a SUBSET of the input in input order (INV-03), and — the load-bearing
    // theory check — every reported segment genuinely has integer copy number 0
    // (log2 ≤ −1.1). This is what stops a fuzz test from rubber-stamping a
    // single-copy loss / neutral / gain segment as a homozygous deletion.
    private static void AssertWellFormedDetection(
        IReadOnlyList<CopyNumberArmSegment> result,
        IReadOnlyList<CopyNumberArmSegment> input)
    {
        result.Should().BeSubsetOf(input, "the detector is a filter, not a constructor (INV-03)");

        // Order-preserving subset: the result is the input in order with some elements dropped.
        int j = 0;
        foreach (var seg in input)
        {
            if (j < result.Count && result[j].Equals(seg)) j++;
        }
        j.Should().Be(result.Count, "reported segments appear in input order (INV-03)");

        foreach (var seg in result)
        {
            CallCopyNumber(seg.Log2Ratio).Should().Be(
                0, "every reported segment has integer copy number 0 (DeepDeletion, INV-01)");
        }
    }

    #region ONCO-CNA-003 — Positive sanity (documented deep / shallow / neutral)

    [Test]
    public void DetectHomozygousDeletions_DocumentedWorkedExample_MatchesDocSpec()
    {
        // Docs §7.1: a 9p segment at log2 −2.0 (≤ −1.1 ⇒ CN 0) IS a homozygous
        // deletion; a 3p segment at log2 −0.5 (CN 1, shallow) is NOT. Only the 9p
        // is reported, and it maps to CDKN2A (9p21.3).
        var input = new[]
        {
            Seg("9p", -2.0, end: 1_000, armLength: 40_000_000), // CN 0 ⇒ homozygous
            Seg("3p", -0.5, end: 1_000, armLength: 90_000_000), // CN 1 ⇒ shallow loss
        };

        var hom = DetectHomozygousDeletions(input);

        hom.Should().ContainSingle();
        hom[0].Arm.Should().Be("9p");
        IdentifyDeletedTumorSuppressors(hom).Should().Equal("CDKN2A");
        AssertWellFormedDetection(hom, input);
    }

    [Test]
    public void DetectHomozygousDeletions_DeepLoss_IsCalledHomozygous()
    {
        // POSITIVE sanity: a clearly deep loss (log2 −3.0 ≪ −1.1) is the canonical
        // homozygous deletion — reported, integer CN 0.
        var input = new[] { Seg("9p", -3.0) };

        DetectHomozygousDeletions(input).Should().ContainSingle();
        IsHomozygousDeletion(input[0]).Should().BeTrue();
        CallCopyNumber(input[0].Log2Ratio).Should().Be(0);
    }

    [Test]
    public void DetectHomozygousDeletions_NeutralSegment_IsNotHomozygous()
    {
        // POSITIVE sanity: a copy-neutral segment (log2 0.0 ⇒ CN 2) is NOT a
        // deletion at all — never reported. Guards against a false deep-del call.
        var input = new[] { Seg("9p", 0.0) };

        DetectHomozygousDeletions(input).Should().BeEmpty("CN-neutral is not a deletion");
        IsHomozygousDeletion(input[0]).Should().BeFalse();
        CallCopyNumber(input[0].Log2Ratio).Should().Be(2);
    }

    [Test]
    public void DetectHomozygousDeletions_SingleCopyLoss_IsNotHomozygous()
    {
        // POSITIVE sanity / INV-02: a single-copy (heterozygous) loss (log2 −0.5 ⇒
        // CN 1, cBioPortal "−1" shallow) is NOT homozygous — one allele remains.
        var input = new[] { Seg("9p", -0.5) };

        DetectHomozygousDeletions(input).Should().BeEmpty("single-copy loss is shallow, not deep (INV-02)");
        IsHomozygousDeletion(input[0]).Should().BeFalse();
        CallCopyNumber(input[0].Log2Ratio).Should().Be(1);
    }

    #endregion

    #region ONCO-CNA-003 — BE: CN exactly 0 (true zero, log2 = −∞)

    [Test]
    public void IsHomozygousDeletion_NegativeInfinityLog2_IsHomozygous_TrueCopyNumberZero()
    {
        // CN EXACTLY 0: a region with literally zero copies has log2 = log2(0/2) =
        // −∞. −∞ ≤ −1.1 ⇒ CN 0 ⇒ homozygous deletion. The hazard is log2(0) = −∞
        // leaking unguarded into the decision as anything but the deepest state.
        var seg = Seg("9p", double.NegativeInfinity);

        var act = () => IsHomozygousDeletion(seg);

        act.Should().NotThrow();
        act().Should().BeTrue("CN exactly 0 (log2 −∞) is the canonical homozygous deletion");
        CallCopyNumber(double.NegativeInfinity).Should().Be(0);
        Log2RatioToCopyNumber(double.NegativeInfinity).Should().Be(0.0, "2^−∞ = 0 absolute copies");
    }

    [Test]
    public void DetectHomozygousDeletions_NegativeInfinityLog2_IsReported_NoCrash()
    {
        // The CN-0 limiting value through the full filter: reported, no crash.
        var input = new[] { Seg("10q", double.NegativeInfinity) };

        var act = () => DetectHomozygousDeletions(input);

        act.Should().NotThrow();
        act().Should().ContainSingle();
        IdentifyDeletedTumorSuppressors(act()).Should().Equal("PTEN");
    }

    [Test]
    public void DetectHomozygousDeletions_VeryDeepFiniteLog2_IsReported()
    {
        // A very deep but finite loss (log2 −50) is still CN 0 — far below −1.1.
        var input = new[] { Seg("9p", -50.0) };

        DetectHomozygousDeletions(input).Should().ContainSingle();
    }

    #endregion

    #region ONCO-CNA-003 — BE: near-0 / CN-0↔CN-1 boundary (inclusive cutoff, no off-by-one)

    [Test]
    public void IsHomozygousDeletion_Log2ExactlyAtDeepCutoff_IsHomozygous_InclusiveBoundary()
    {
        // BOUNDARY: log2 EXACTLY −1.1 is a homozygous deletion — the CNVkit cutoff
        // is INCLUSIVE ("≤ each threshold in sequence", §6.1). Guards the off-by-one
        // that would mis-classify the boundary as CN 1.
        var seg = Seg("9p", DeepDeletionLog2Cutoff);

        IsHomozygousDeletion(seg).Should().BeTrue("log2 exactly at the deep cutoff is CN 0 (inclusive, §6.1)");
        CallCopyNumber(DeepDeletionLog2Cutoff).Should().Be(0);
    }

    [Test]
    public void IsHomozygousDeletion_Log2JustAboveDeepCutoff_IsNotHomozygous_Cn1()
    {
        // NEAR-0 (the other side of the boundary): log2 just ABOVE −1.1 is a
        // single-copy loss (CN 1), NOT homozygous. Pairs with the exactly-at-cutoff
        // case to pin the boundary direction — no off-by-one.
        double justAbove = Math.BitIncrement(DeepDeletionLog2Cutoff);
        var seg = Seg("9p", justAbove);

        IsHomozygousDeletion(seg).Should().BeFalse("log2 just above the deep cutoff is CN 1 (shallow), not CN 0");
        justAbove.Should().BeGreaterThan(DeepDeletionLog2Cutoff);
        CallCopyNumber(justAbove).Should().Be(1);
    }

    [Test]
    public void IsHomozygousDeletion_Log2JustBelowDeepCutoff_IsHomozygous()
    {
        // Just INSIDE the deep side: log2 just below −1.1 is still CN 0 ⇒ reported.
        double justBelow = Math.BitDecrement(DeepDeletionLog2Cutoff);
        var seg = Seg("9p", justBelow);

        IsHomozygousDeletion(seg).Should().BeTrue("log2 just below the deep cutoff is CN 0");
        justBelow.Should().BeLessThan(DeepDeletionLog2Cutoff);
        CallCopyNumber(justBelow).Should().Be(0);
    }

    [Test]
    public void DetectHomozygousDeletions_BoundaryTrio_KeepsOnlyTheCn0Ones()
    {
        // The three boundary points through the filter: exactly −1.1 (CN 0, kept),
        // just above (CN 1, dropped), just below (CN 0, kept) — in input order.
        var atCutoff = Seg("9p", DeepDeletionLog2Cutoff);
        var justAbove = Seg("10q", Math.BitIncrement(DeepDeletionLog2Cutoff));
        var justBelow = Seg("17p", Math.BitDecrement(DeepDeletionLog2Cutoff));
        var input = new[] { atCutoff, justAbove, justBelow };

        var hom = DetectHomozygousDeletions(input);

        hom.Should().Equal(atCutoff, justBelow);
        AssertWellFormedDetection(hom, input);
    }

    #endregion

    #region ONCO-CNA-003 — BE: single-bin deletion (width-1 segment, no length crash)

    [Test]
    public void DetectHomozygousDeletions_SingleBinWidthOne_DeepLoss_IsReported_NoCrash()
    {
        // SINGLE-BIN DELETION: a width-1 segment (End = Start + 1, Length = 1) whose
        // log2 ≤ −1.1. The hazard is a DivideByZero / length-0 crash; the decision
        // depends only on log2, so a width-1 deep-loss segment is reported. No throw.
        var input = new[] { Seg("9p", -2.0, start: 500_000, end: 500_001) };

        var act = () => DetectHomozygousDeletions(input);

        act.Should().NotThrow();
        act().Should().ContainSingle();
        input[0].Length.Should().Be(1);
    }

    [Test]
    public void IsHomozygousDeletion_SingleBinWidthOne_DoesNotThrow()
    {
        // The single-segment predicate on a width-1 segment: must evaluate cleanly.
        var seg = Seg("9p", -2.0, start: 0, end: 1);

        var act = () => IsHomozygousDeletion(seg);

        act.Should().NotThrow();
        act().Should().BeTrue("width-1, deep-loss ⇒ homozygous deletion");
        seg.Length.Should().Be(1);
    }

    [Test]
    public void IsHomozygousDeletion_SingleBinOnArmLengthOne_NoDivideByZero()
    {
        // Degenerate arm of length 1 with a width-1 segment: the CN-0 decision never
        // touches ArmLength (it is purely a log2 comparison), so there is no
        // DivideByZero even though ArmFraction = 1/1. log2 −2.0 ⇒ still homozygous.
        var seg = Seg("9p", -2.0, start: 0, end: 1, armLength: 1);

        var act = () => IsHomozygousDeletion(seg);

        act.Should().NotThrow();
        act().Should().BeTrue();
    }

    #endregion

    #region ONCO-CNA-003 — BE: NaN log2 (documented no-call, no false deep-del)

    [Test]
    public void IsHomozygousDeletion_NaNLog2_IsNotHomozygous_DocumentedNoCall()
    {
        // NaN log2 is the CNVkit no-call ⇒ neutral reference CN (rounded ploidy 2)
        // ⇒ NOT a homozygous deletion (§3.3/§6.1). The hazard is a NaN slipping
        // through the `== 0` test as a false deep deletion.
        var seg = Seg("9p", double.NaN);

        IsHomozygousDeletion(seg).Should().BeFalse("NaN log2 is a neutral no-call, not CN 0 (§3.3)");
        CallCopyNumber(double.NaN).Should().Be(2, "no-call returns rounded ploidy");
    }

    [Test]
    public void DetectHomozygousDeletions_NaNLog2_IsNotReported()
    {
        var input = new[] { Seg("9p", double.NaN) };

        DetectHomozygousDeletions(input).Should().BeEmpty();
    }

    [Test]
    public void IsHomozygousDeletion_PositiveInfinityLog2_IsNotHomozygous()
    {
        // log2 = +∞ is an unbounded amplification, the opposite extreme — never a
        // homozygous deletion. +∞ ≤ −1.1 is false ⇒ not CN 0.
        var seg = Seg("9p", double.PositiveInfinity);

        IsHomozygousDeletion(seg).Should().BeFalse("a +∞ amplification is not a deep deletion");
    }

    #endregion

    #region ONCO-CNA-003 — BE: empty / null / malformed segment (documented guards)

    [Test]
    public void DetectHomozygousDeletions_EmptyInput_ReturnsEmpty_NoThrow()
    {
        // Empty stream ⇒ empty result (§6.1) — the zero-segment boundary.
        DetectHomozygousDeletions(Array.Empty<CopyNumberArmSegment>()).Should().BeEmpty();
    }

    [Test]
    public void DetectHomozygousDeletions_NullSegments_ThrowsArgumentNull()
    {
        // Documented null guard (§3.3) — not a NullReference deep in the loop.
        var act = () => DetectHomozygousDeletions(null!);

        act.Should().Throw<ArgumentNullException>();
    }

    [Test]
    public void IdentifyDeletedTumorSuppressors_NullDeletions_ThrowsArgumentNull()
    {
        var act = () => IdentifyDeletedTumorSuppressors(null!);

        act.Should().Throw<ArgumentNullException>();
    }

    [Test]
    public void DetectHomozygousDeletions_NonPositiveArmLength_ThrowsArgument()
    {
        // ArmLength ≤ 0 is rejected by ValidateArmSegment (§3.3) up front.
        var zeroArm = new[] { Seg("9p", -2.0, armLength: 0) };
        var negArm = new[] { Seg("9p", -2.0, armLength: -1) };

        ((Action)(() => DetectHomozygousDeletions(zeroArm))).Should().Throw<ArgumentException>();
        ((Action)(() => DetectHomozygousDeletions(negArm))).Should().Throw<ArgumentException>();
    }

    [Test]
    public void DetectHomozygousDeletions_EndNotAfterStart_ThrowsArgument()
    {
        // End ≤ Start (zero/negative-width segment) ⇒ ArgumentException (§3.3).
        var zeroWidth = new[] { Seg("9p", -2.0, start: 100_000, end: 100_000) };
        var inverted = new[] { Seg("9p", -2.0, start: 100_000, end: 50_000) };

        ((Action)(() => DetectHomozygousDeletions(zeroWidth))).Should().Throw<ArgumentException>();
        ((Action)(() => DetectHomozygousDeletions(inverted))).Should().Throw<ArgumentException>();
    }

    [Test]
    public void DetectHomozygousDeletions_NonPositivePloidy_ThrowsArgumentOutOfRange()
    {
        // Non-positive ploidy ⇒ ArgumentOutOfRangeException (§3.3, via CallCopyNumber).
        var input = new[] { Seg("9p", -2.0) };

        ((Action)(() => DetectHomozygousDeletions(input, ploidy: 0.0)))
            .Should().Throw<ArgumentOutOfRangeException>();
        ((Action)(() => DetectHomozygousDeletions(input, ploidy: -2.0)))
            .Should().Throw<ArgumentOutOfRangeException>();
    }

    [Test]
    public void DetectHomozygousDeletions_MalformedThresholds_ThrowsArgument()
    {
        // thresholds not four strictly ascending values ⇒ ArgumentException (§3.3).
        var input = new[] { Seg("9p", -2.0) };
        var tooFew = new[] { -1.1, -0.25, 0.2 };
        var nonAscending = new[] { -1.1, -0.25, 0.2, 0.1 };

        ((Action)(() => DetectHomozygousDeletions(input, thresholds: tooFew)))
            .Should().Throw<ArgumentException>();
        ((Action)(() => DetectHomozygousDeletions(input, thresholds: nonAscending)))
            .Should().Throw<ArgumentException>();
    }

    #endregion

    #region ONCO-CNA-003 — BE: tumour-suppressor mapping (case-insensitive, panel order, subset)

    [Test]
    public void IdentifyDeletedTumorSuppressors_EmptyDeletions_ReturnsEmpty()
    {
        IdentifyDeletedTumorSuppressors(Array.Empty<CopyNumberArmSegment>()).Should().BeEmpty();
    }

    [Test]
    public void IdentifyDeletedTumorSuppressors_ArmCaseInsensitive_MatchesPanel()
    {
        // INV-04: arm matching is Ordinal-ignore-case. A "9P" hom-del still flags CDKN2A.
        var dels = new[] { Seg("9P", -2.0) };

        IdentifyDeletedTumorSuppressors(dels).Should().Equal("CDKN2A");
    }

    [Test]
    public void IdentifyDeletedTumorSuppressors_SharedArm_ReportsBothPanelGenesInOrder()
    {
        // 13q carries both RB1 and BRCA2; a single 13q hom-del reports both, in
        // panel order (RB1 before BRCA2), each once (INV-04, distinct).
        var dels = new[] { Seg("13q", -2.0) };

        IdentifyDeletedTumorSuppressors(dels).Should().Equal("RB1", "BRCA2");
    }

    [Test]
    public void IdentifyDeletedTumorSuppressors_NonPanelArm_ReportsNothing()
    {
        // A homozygous deletion on an arm outside the six-gene panel maps to nothing.
        var dels = new[] { Seg("1p", -2.0) };

        IdentifyDeletedTumorSuppressors(dels).Should().BeEmpty();
    }

    [Test]
    public void IdentifyDeletedTumorSuppressors_DuplicateArm_ReportsGeneOnce()
    {
        // Two hom-dels on the same arm still report each panel gene exactly once (INV-04).
        var dels = new[] { Seg("17p", -2.0), Seg("17p", -3.0) };

        IdentifyDeletedTumorSuppressors(dels).Should().Equal("TP53");
    }

    #endregion

    #region ONCO-CNA-003 — BE: broad random fuzz (well-formed filter, predicate agreement)

    [Test]
    [CancelAfter(30000)]
    public void DetectHomozygousDeletions_RandomSegments_AlwaysWellFormedFilter(
        [Values(20260105, 777, 313131)] int seed)
    {
        // Fuzz mixed streams of arm segments whose log2 spans the deep-deletion
        // boundary (and the deep/shallow/neutral/gain spectrum), including width-1
        // and whole-arm spans plus the ±∞ extremes. Every result must be an
        // order-preserving subset whose members all have integer CN 0 — no shallow/
        // neutral leak, no crash.
        var rng = new Random(seed);
        var arms = new[] { "9p", "10q", "17p", "13q", "17q", "1p" };

        for (int t = 0; t < 5000; t++)
        {
            int n = rng.Next(0, 8);
            var input = new List<CopyNumberArmSegment>(n);
            for (int i = 0; i < n; i++)
            {
                string arm = arms[rng.Next(arms.Length)];
                long len = (long)(rng.NextDouble() * ArmLen) + 1; // 1 bp .. whole arm
                long start = rng.Next(0, 50_000);
                // Concentrate around the −1.1 boundary, with occasional extremes.
                double log2 = rng.Next(0, 20) switch
                {
                    0 => double.NegativeInfinity,
                    1 => double.PositiveInfinity,
                    2 => double.NaN,
                    _ => -1.1 + (rng.NextDouble() * 2.0 - 1.0) * 0.5, // ±0.5 around −1.1
                };
                input.Add(Seg(arm, log2, start: start, end: start + len));
            }

            var hom = DetectHomozygousDeletions(input);

            AssertWellFormedDetection(hom, input);

            // Cross-check the predicate: a segment is in the result iff the predicate says so.
            var expected = input.Where(s => IsHomozygousDeletion(s)).ToList();
            hom.Should().Equal(expected, "filter agrees with the single-segment predicate (INV-01)");
        }
    }

    [Test]
    [CancelAfter(30000)]
    public void IsHomozygousDeletion_EquivalentToIntegerCopyNumberZero(
        [Values(11, 2024, 99999)] int seed)
    {
        // Property: homozygous ⇔ CallCopyNumber(log2) == 0 ⇔ log2 ≤ −1.1 (for finite
        // log2). Fuzz the log2 axis densely around the cutoff and assert the
        // predicate equals BOTH the integer-CN call and the direct ≤ comparison —
        // a single off-by-one would break this.
        var rng = new Random(seed);

        for (int i = 0; i < 20000; i++)
        {
            double log2 = -1.1 + (rng.NextDouble() * 2.0 - 1.0) * 1.0; // ±1.0 around −1.1 (finite)
            var seg = Seg("9p", log2);

            bool predicate = IsHomozygousDeletion(seg);
            predicate.Should().Be(CallCopyNumber(log2) == 0, "predicate ⇔ integer CN 0 (INV-01)");
            predicate.Should().Be(log2 <= DeepDeletionLog2Cutoff,
                "for finite log2, CN 0 ⇔ log2 ≤ first cutoff (inclusive, §2.2)");
        }
    }

    #endregion
}
