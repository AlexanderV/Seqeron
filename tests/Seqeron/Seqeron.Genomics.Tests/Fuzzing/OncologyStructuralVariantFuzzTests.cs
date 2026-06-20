using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Oncology;
using static Seqeron.Genomics.Oncology.OncologyAnalyzer;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Fuzz tests for the Oncology somatic complex-rearrangement (chromothripsis)
/// structural-variant layer — ONCO-SV-001. The unit under test is the
/// chromothripsis-inference screen in
/// src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs:
///   • <see cref="OncologyAnalyzer.CountCopyNumberStateOscillations"/> — counts
///     adjacent per-segment copy-number state transitions (the "oscillating CN
///     changes" first-pass quantity);
///   • <see cref="OncologyAnalyzer.TestBreakpointClustering"/> — coefficient-of-
///     variation clustering test of inter-breakpoint gaps against the exponential
///     random-breakpoint null (criterion A);
///   • <see cref="OncologyAnalyzer.ClassifyComplexRearrangement"/> — full
///     Chromothripsis / NotComplex call with the confidence tier (criterion B +
///     the Cortés-Ciriano SV-burden floor).
///
/// SCOPE. This file is scoped strictly to ONCO-SV-001 — the SV/chromothripsis
/// PATTERN layer that consumes already-segmented per-region copy numbers, a
/// clustered SV burden and breakpoint coordinates. The generic per-event SV typing
/// (DEL/DUP/INV/TRA from read orientation) lives in StructuralVariantAnalyzer and
/// the fusion-breakpoint analysis (rows 100–102, ONCO-FUSION-*) is a SEPARATE unit;
/// neither is exercised here (docs §1, §2.5, §5.2).
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate / boundary inputs to a unit and asserts that the code
/// NEVER fails in an undisciplined way: no hang, no nonsense output, no *unhandled*
/// runtime exception (DivideByZero on an empty/single breakpoint set, Overflow /
/// non-termination on an enormous breakpoint count, a NaN CV silently leaking into
/// the IsClustered flag). Every input must resolve to EITHER a well-defined,
/// theory-correct value OR a documented, intentional outcome (ArgumentNullException
/// for null lists). The headline hazards for THIS unit are:
///   • zero breakpoints / zero segments → no DivideByZero on the empty gap set or
///     the empty distinct-state set; clustering returns IsClustered = false (no CV
///     from 0 gaps), classification returns NotComplex (§3.3, §6.1, INV-01);
///   • a single breakpoint / single segment → no variance-of-one or
///     divide-by-(n−1) crash on a 1-element set; with < 3 breakpoints clustering is
///     undefined (< 2 gaps) ⇒ IsClustered = false; a single segment has 0
///     oscillations ⇒ NotComplex (§3.3, §6.1);
///   • genome-wide shattering — very many clustered breakpoints with an oscillating
///     CN profile — must terminate (no non-termination), not overflow on a large
///     breakpoint/SV count, and be called Chromothripsis with High confidence per
///     the documented gate (INV-02, INV-03).
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: ONCO-SV-001 — Somatic complex-rearrangement classification (Oncology)
/// Checklist: docs/checklists/03_FUZZING.md, row 119.
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — граничні значення: 0, -1, MaxInt, empty.
///     Targets (checklist row 119): "zero breakpoints, single breakpoint,
///     genome-wide shattering".
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The documented contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// Complex_Rearrangement_Classification.md
/// (docs/algorithms/Oncology/Complex_Rearrangement_Classification.md):
///   • oscillation = adjacent CN-state transition; for n segments count ∈ [0, n−1]
///     and < 2 segments ⇒ 0                                          (§2.2, INV-01)
///   • oscillating-segment count = transitions + 1 when > 0, else 0  (§3.2, §4.1)
///   • Chromothripsis ⇔ distinct states ∈ [2, 3] AND oscillations ≥ 10 AND
///     SV burden ≥ 6                                                 (§4.1, INV-02)
///   • a monotone / >3-distinct-state profile is NEVER Chromothripsis (INV-04)
///   • confidence = High iff oscillating segments ≥ 7; Low iff in [4, 6]; None
///     iff < 4                                                       (§3.2, INV-03)
///   • thresholds: MinOscillatingCopyNumberChanges 10, MaxChromothripsis…States 3,
///     MinChromothripsisSvBurden 6, High…Segments 7, Low…Segments 4 (§4.2)
///   • clustering: CV = sd/mean of inter-breakpoint gaps; IsClustered ⇔ CV > 1;
///     < 3 breakpoints ⇒ CV undefined ⇒ IsClustered = false; all-coincident
///     breakpoints (mean ≤ 0) ⇒ not assessable                  (§2.2, §3.3, INV-05)
///   • regular-spacing breakpoints ⇒ CV ≈ 0 (not clustered)         (INV-05)
///   • null SegmentCopyNumbers / breakpointPositions ⇒ ArgumentNullException (§3.3)
///   • worked example: profile 2,1,2,1,2,1,2,1,2,1,2 (10 oscillations, 2 states)
///     with SV burden 12 ⇒ Chromothripsis, High, OscillationCount 10            (§7.1)
///
/// All randomness is LOCALLY seeded (new Random(seed)); no shared static Rng.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public sealed class OncologyStructuralVariantFuzzTests
{
    // Documented thresholds (§4.2), mirrored locally so the test owns the expected
    // gate independently of the source constants.
    private const int MinOscillations = 10;   // first-pass screen
    private const int MaxStates = 3;          // two-or-three-state hallmark
    private const int MinSvBurden = 6;        // focal <6-SV exclusion
    private const int HighSegments = 7;       // ≥7 adjacent oscillating segments
    private const int LowSegments = 4;        // [4, 6] adjacent oscillating segments

    // ── Well-formed-result assertion helper ──────────────────────────────────
    // Pins the documented structural invariants on EVERY accepted classification,
    // regardless of input: the oscillation count is in [0, n−1] (INV-01), the
    // oscillating-segment count is the documented transitions+1/0 derivation, a
    // Chromothripsis call implies the full conjunction gate (INV-02/INV-04), and
    // the confidence tier matches the oscillating-segment thresholds (INV-03).
    // This is what stops a fuzz test from rubber-stamping a green call.
    private static void AssertWellFormedResult(in ComplexRearrangementResult r, int segmentCount)
    {
        r.OscillationCount.Should().BeInRange(
            0, Math.Max(0, segmentCount - 1),
            "oscillation count is one-per-adjacent-pair ⇒ ∈ [0, n−1] (INV-01)");

        int expectedSegments = r.OscillationCount > 0 ? r.OscillationCount + 1 : 0;
        r.OscillatingSegmentCount.Should().Be(
            expectedSegments, "oscillating segments = transitions+1 when >0 else 0 (§4.1)");

        Enum.IsDefined(typeof(ComplexRearrangementType), r.Type).Should().BeTrue();
        Enum.IsDefined(typeof(ChromothripsisConfidence), r.Confidence).Should().BeTrue();

        // INV-03: confidence tier is a pure function of the oscillating-segment count.
        ChromothripsisConfidence expectedConf =
            r.OscillatingSegmentCount >= HighSegments ? ChromothripsisConfidence.High
            : r.OscillatingSegmentCount >= LowSegments ? ChromothripsisConfidence.Low
            : ChromothripsisConfidence.None;
        r.Confidence.Should().Be(expectedConf, "confidence is derived from oscillating segments (INV-03)");

        // INV-02 / INV-04: a Chromothripsis call implies the full conjunction gate.
        if (r.Type == ComplexRearrangementType.Chromothripsis)
        {
            r.DistinctStateCount.Should().BeInRange(2, MaxStates,
                "Chromothripsis ⇒ two-/three-state hallmark (INV-02/INV-04)");
            r.OscillationCount.Should().BeGreaterThanOrEqualTo(MinOscillations,
                "Chromothripsis ⇒ ≥10 oscillations first-pass screen (INV-02)");
            r.StructuralVariantCount.Should().BeGreaterThanOrEqualTo(MinSvBurden,
                "Chromothripsis ⇒ ≥6 clustered SV burden (INV-02)");
        }
    }

    // Builds an alternating two-state profile of the requested segment length:
    // 2,1,2,1,… — n segments give n−1 oscillations and exactly 2 distinct states.
    private static int[] AlternatingProfile(int segments)
    {
        var profile = new int[segments];
        for (int i = 0; i < segments; i++)
        {
            profile[i] = i % 2 == 0 ? 2 : 1;
        }

        return profile;
    }

    #region ONCO-SV-001 — Positive sanity (documented chromothripsis / non-call)

    [Test]
    public void ClassifyComplexRearrangement_DocumentedWorkedExample_IsHighConfidenceChromothripsis()
    {
        // Docs §7.1 worked example: profile 2,1,2,1,2,1,2,1,2,1,2 (11 segments) ⇒
        // 10 transitions, 2 distinct states {1,2}, 11 oscillating segments ≥ 7 ⇒
        // High; 10 ≥ 10 and SV burden 12 ≥ 6 and states ∈ [2,3] ⇒ Chromothripsis.
        var input = new ComplexRearrangementInput(
            new[] { 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2 },
            StructuralVariantCount: 12);

        var result = ClassifyComplexRearrangement(input);

        result.Type.Should().Be(ComplexRearrangementType.Chromothripsis);
        result.Confidence.Should().Be(ChromothripsisConfidence.High);
        result.OscillationCount.Should().Be(10);
        result.OscillatingSegmentCount.Should().Be(11);
        result.DistinctStateCount.Should().Be(2);
        result.StructuralVariantCount.Should().Be(12);
        AssertWellFormedResult(result, segmentCount: 11);
    }

    [Test]
    public void ClassifyComplexRearrangement_MonotoneAscendingProfile_IsNotComplex()
    {
        // A progressive amplification profile (1,2,3,4,…) oscillates on every step
        // BUT has many distinct states (> 3) — the two-state hallmark fails, so it
        // is NotComplex regardless of length (INV-04, §6.1). The distinguishing
        // criterion B at work: clustered transitions are not enough.
        var profile = Enumerable.Range(1, 12).ToArray(); // 1..12 ⇒ 11 transitions, 12 states
        var input = new ComplexRearrangementInput(profile, StructuralVariantCount: 50);

        var result = ClassifyComplexRearrangement(input);

        result.Type.Should().Be(ComplexRearrangementType.NotComplex);
        result.OscillationCount.Should().Be(11);
        result.DistinctStateCount.Should().Be(12);
        AssertWellFormedResult(result, profile.Length);
    }

    [Test]
    public void TestBreakpointClustering_TightClusterPlusOutliers_IsClustered()
    {
        // Criterion A: a tight cluster of breakpoints plus a couple of distant
        // outliers gives many short gaps + a few long gaps ⇒ over-dispersion ⇒
        // CV > 1 ⇒ IsClustered = true (INV-05). Hand-built clustered case.
        var positions = new long[]
        {
            100, 101, 102, 103, 104, 105, 106, 107, 108, 109, // tight cluster
            5_000_000, 50_000_000                              // distant outliers
        };

        var result = TestBreakpointClustering(positions);

        result.BreakpointCount.Should().Be(positions.Length);
        result.CoefficientOfVariation.Should().BeGreaterThan(1.0);
        result.IsClustered.Should().BeTrue("over-dispersed gaps (CV>1) flag clustering (INV-05)");
    }

    [Test]
    public void TestBreakpointClustering_EvenlySpacedBreakpoints_IsNotClustered()
    {
        // INV-05: perfectly regular spacing ⇒ all gaps equal ⇒ CV = 0 (far below
        // the exponential-null CV of 1) ⇒ NOT clustered. The negative control for
        // the clustering test — uniform breakpoints must not false-positive.
        var positions = Enumerable.Range(0, 20).Select(i => (long)(i * 1000)).ToArray();

        var result = TestBreakpointClustering(positions);

        result.CoefficientOfVariation.Should().BeApproximately(0.0, 1e-12);
        result.IsClustered.Should().BeFalse("uniform gaps give CV≈0, not clustered (INV-05)");
    }

    #endregion

    #region ONCO-SV-001 — BE: zero breakpoints / zero segments (no DivideByZero)

    [Test]
    public void TestBreakpointClustering_EmptyBreakpointSet_NotClustered_NoDivideByZero()
    {
        // "Zero breakpoints": an empty set has NO gaps. A naive CV = sd/mean would
        // divide by zero (mean of an empty gap set); the contract short-circuits
        // < 3 breakpoints to a defined not-assessable result (§3.3, §6.1).
        var act = () => TestBreakpointClustering(Array.Empty<long>());

        act.Should().NotThrow();
        var result = act();
        result.BreakpointCount.Should().Be(0);
        result.MeanGap.Should().Be(0.0);
        result.CoefficientOfVariation.Should().Be(0.0);
        result.IsClustered.Should().BeFalse();
    }

    [Test]
    public void ClassifyComplexRearrangement_EmptySegments_IsNotComplex_NoCrash()
    {
        // "Zero segments": an empty CN profile has 0 oscillations and 0 distinct
        // states. The distinct-state set is empty (the source guards Count==0 to
        // avoid Distinct().Count() oddities), so the call is NotComplex with None
        // confidence — no DivideByZero / empty-aggregate crash (§6.1, INV-01).
        var input = new ComplexRearrangementInput(Array.Empty<int>(), StructuralVariantCount: 100);

        var result = ClassifyComplexRearrangement(input);

        result.Type.Should().Be(ComplexRearrangementType.NotComplex);
        result.Confidence.Should().Be(ChromothripsisConfidence.None);
        result.OscillationCount.Should().Be(0);
        result.OscillatingSegmentCount.Should().Be(0);
        result.DistinctStateCount.Should().Be(0);
        AssertWellFormedResult(result, segmentCount: 0);
    }

    [Test]
    public void CountCopyNumberStateOscillations_EmptyProfile_IsZero()
    {
        // Boundary: 0 segments ⇒ 0 transitions (the loop body never runs).
        CountCopyNumberStateOscillations(Array.Empty<int>()).Should().Be(0);
    }

    [Test]
    public void ClassifyComplexRearrangement_ZeroSvBurden_IsNotComplex_EvenWithManyOscillations()
    {
        // Boundary SV burden = 0: even a perfect 10-oscillation two-state profile
        // is NotComplex when the clustered SV burden is 0 (< 6 focal exclusion,
        // INV-02). Decouples the CN-oscillation gate from the SV-burden gate.
        var input = new ComplexRearrangementInput(AlternatingProfile(11), StructuralVariantCount: 0);

        var result = ClassifyComplexRearrangement(input);

        result.OscillationCount.Should().Be(10);
        result.DistinctStateCount.Should().Be(2);
        result.Type.Should().Be(ComplexRearrangementType.NotComplex);
        AssertWellFormedResult(result, segmentCount: 11);
    }

    #endregion

    #region ONCO-SV-001 — BE: single breakpoint / single segment (no variance-of-one crash)

    [Test]
    public void TestBreakpointClustering_SingleBreakpoint_NotClustered_NoVarianceCrash()
    {
        // "Single breakpoint": one position defines ZERO gaps — variance/CV are
        // undefined on a 0-element gap set. A naive variance/(n−1) would crash or
        // produce NaN; the < 3 guard returns a defined not-assessable result (§3.3).
        var act = () => TestBreakpointClustering(new long[] { 42 });

        act.Should().NotThrow();
        var result = act();
        result.BreakpointCount.Should().Be(1);
        result.CoefficientOfVariation.Should().Be(0.0);
        double.IsNaN(result.CoefficientOfVariation).Should().BeFalse("no NaN CV leaks from a 1-element set");
        result.IsClustered.Should().BeFalse();
    }

    [Test]
    public void TestBreakpointClustering_TwoBreakpoints_NotClustered_OneGapCannotDefineCv()
    {
        // Two breakpoints ⇒ exactly ONE gap. A CV needs ≥ 2 gaps (≥ 3 breakpoints)
        // to be meaningful; the documented < 3 guard returns IsClustered = false
        // rather than a degenerate CV = 0 over a single observation (§3.3, §6.1).
        var result = TestBreakpointClustering(new long[] { 100, 500 });

        result.BreakpointCount.Should().Be(2);
        result.IsClustered.Should().BeFalse("< 3 breakpoints ⇒ < 2 gaps ⇒ clustering undefined");
    }

    [Test]
    public void ClassifyComplexRearrangement_SingleSegment_IsNotComplex_NoCrash()
    {
        // "Single segment": one CN value ⇒ 0 transitions (no adjacent pair), 1
        // distinct state. distinct < 2 fails the hallmark ⇒ NotComplex with None
        // confidence — no DivideByZero on a 1-element profile (§6.1, INV-01).
        var input = new ComplexRearrangementInput(new[] { 2 }, StructuralVariantCount: 100);

        var result = ClassifyComplexRearrangement(input);

        result.Type.Should().Be(ComplexRearrangementType.NotComplex);
        result.Confidence.Should().Be(ChromothripsisConfidence.None);
        result.OscillationCount.Should().Be(0);
        result.OscillatingSegmentCount.Should().Be(0);
        result.DistinctStateCount.Should().Be(1);
        AssertWellFormedResult(result, segmentCount: 1);
    }

    [Test]
    public void CountCopyNumberStateOscillations_SingleSegment_IsZero()
    {
        // Boundary: 1 segment ⇒ 0 transitions (no adjacent pair to compare).
        CountCopyNumberStateOscillations(new[] { 7 }).Should().Be(0);
    }

    [Test]
    public void TestBreakpointClustering_AllCoincidentBreakpoints_NotAssessable_NoDivideByMeanZero()
    {
        // Degenerate: ≥ 3 breakpoints all at the SAME coordinate ⇒ every gap is 0
        // ⇒ mean gap 0. CV = sd/mean would divide by zero; the documented mean ≤ 0
        // guard returns a clean not-assessable result instead of NaN/∞ (§3.3).
        var result = TestBreakpointClustering(new long[] { 1000, 1000, 1000, 1000 });

        result.MeanGap.Should().Be(0.0);
        result.CoefficientOfVariation.Should().Be(0.0);
        double.IsNaN(result.CoefficientOfVariation).Should().BeFalse();
        result.IsClustered.Should().BeFalse();
    }

    #endregion

    #region ONCO-SV-001 — BE: genome-wide shattering (large counts terminate, no overflow)

    [Test]
    [CancelAfter(30000)]
    public void ClassifyComplexRearrangement_GenomeWideShattering_IsHighConfidenceChromothripsis_NoOverflow()
    {
        // "Genome-wide shattering": a very long two-state oscillating profile with a
        // large clustered SV burden — the canonical chromothripsis signature. Must
        // TERMINATE (CancelAfter) and call Chromothripsis / High without overflow on
        // the large counts. 50_000 segments ⇒ 49_999 oscillations.
        var profile = AlternatingProfile(50_000);
        var input = new ComplexRearrangementInput(profile, StructuralVariantCount: 10_000);

        var result = ClassifyComplexRearrangement(input);

        result.Type.Should().Be(ComplexRearrangementType.Chromothripsis);
        result.Confidence.Should().Be(ChromothripsisConfidence.High);
        result.OscillationCount.Should().Be(49_999);
        result.OscillatingSegmentCount.Should().Be(50_000);
        result.DistinctStateCount.Should().Be(2);
        AssertWellFormedResult(result, profile.Length);
    }

    [Test]
    [CancelAfter(30000)]
    public void ClassifyComplexRearrangement_MaxIntSvBurden_DoesNotOverflowGate()
    {
        // Boundary MaxInt SV burden: Int32.MaxValue ≥ 6 must satisfy the SV-burden
        // gate without overflow in the comparison; the echoed count is preserved.
        var input = new ComplexRearrangementInput(AlternatingProfile(11), StructuralVariantCount: int.MaxValue);

        var result = ClassifyComplexRearrangement(input);

        result.StructuralVariantCount.Should().Be(int.MaxValue);
        result.Type.Should().Be(ComplexRearrangementType.Chromothripsis);
        AssertWellFormedResult(result, segmentCount: 11);
    }

    [Test]
    [CancelAfter(30000)]
    public void TestBreakpointClustering_HugeShatteredBreakpointSet_Terminates_NoOverflow()
    {
        // Genome-wide shattering at the clustering layer: tens of thousands of
        // breakpoints with extreme coordinate magnitudes (up to genome scale). The
        // sort + gap math must terminate and not overflow the long subtraction or
        // produce a NaN/∞ CV. Tightly clustered + a far outlier ⇒ CV > 1.
        var rng = new Random(20260620);
        var positions = new List<long>(40_001);
        for (int i = 0; i < 40_000; i++)
        {
            positions.Add(rng.NextInt64(0, 200_000)); // dense cluster near origin
        }

        positions.Add(3_000_000_000L); // a far outlier beyond Int32 range (no overflow)

        var result = TestBreakpointClustering(positions);

        result.BreakpointCount.Should().Be(40_001);
        double.IsNaN(result.CoefficientOfVariation).Should().BeFalse();
        double.IsInfinity(result.CoefficientOfVariation).Should().BeFalse();
        result.CoefficientOfVariation.Should().BeGreaterThan(1.0);
        result.IsClustered.Should().BeTrue();
    }

    #endregion

    #region ONCO-SV-001 — Boundary: exact thresholds of the gate (off-by-one edges)

    [Test]
    public void ClassifyComplexRearrangement_ExactlyTenOscillations_MeetsScreen()
    {
        // Exactly MinOscillations (10): the ≥ comparison is inclusive, so a profile
        // with precisely 10 transitions (11 alternating segments) passes the
        // first-pass screen ⇒ Chromothripsis with sufficient SV burden (INV-02).
        var input = new ComplexRearrangementInput(AlternatingProfile(11), StructuralVariantCount: 6);

        var result = ClassifyComplexRearrangement(input);

        result.OscillationCount.Should().Be(10);
        result.Type.Should().Be(ComplexRearrangementType.Chromothripsis);
    }

    [Test]
    public void ClassifyComplexRearrangement_NineOscillations_BelowScreen_IsNotComplex()
    {
        // One below the screen (9 oscillations, 10 alternating segments) ⇒ fails
        // the ≥10 first-pass screen ⇒ NotComplex even with ample SV burden, but
        // High confidence (10 oscillating segments ≥ 7). Confirms the screen and
        // confidence tier are independent gates (INV-02 vs INV-03).
        var input = new ComplexRearrangementInput(AlternatingProfile(10), StructuralVariantCount: 100);

        var result = ClassifyComplexRearrangement(input);

        result.OscillationCount.Should().Be(9);
        result.Type.Should().Be(ComplexRearrangementType.NotComplex);
        result.Confidence.Should().Be(ChromothripsisConfidence.High);
    }

    [Test]
    public void ClassifyComplexRearrangement_FiveSvBurden_BelowFloor_IsNotComplex()
    {
        // One below the SV floor (5 < 6): a 10-oscillation two-state profile is
        // NotComplex when the clustered SV burden is 5 (focal <6-SV exclusion,
        // INV-02). The complementary boundary to the exact-10-oscillation test.
        var input = new ComplexRearrangementInput(AlternatingProfile(11), StructuralVariantCount: 5);

        var result = ClassifyComplexRearrangement(input);

        result.Type.Should().Be(ComplexRearrangementType.NotComplex);
    }

    [Test]
    public void ClassifyComplexRearrangement_ThreeStateProfile_StillEligible_FourStateNot()
    {
        // The hallmark cap is ≤ 3 distinct states. A 3-state oscillating profile
        // (2,1,3,2,1,3,…) with ≥10 oscillations + SV burden is STILL Chromothripsis
        // (states ∈ [2,3]); adding a 4th state breaks the hallmark ⇒ NotComplex
        // (INV-02/INV-04 distinct-state boundary).
        var threeState = new[] { 2, 1, 3, 2, 1, 3, 2, 1, 3, 2, 1, 3 }; // 11 transitions, 3 states
        var fourState = new[] { 2, 1, 3, 4, 2, 1, 3, 4, 2, 1, 3, 4 };  // 11 transitions, 4 states

        var threeResult = ClassifyComplexRearrangement(
            new ComplexRearrangementInput(threeState, StructuralVariantCount: 8));
        var fourResult = ClassifyComplexRearrangement(
            new ComplexRearrangementInput(fourState, StructuralVariantCount: 8));

        threeResult.DistinctStateCount.Should().Be(3);
        threeResult.Type.Should().Be(ComplexRearrangementType.Chromothripsis);
        fourResult.DistinctStateCount.Should().Be(4);
        fourResult.Type.Should().Be(ComplexRearrangementType.NotComplex);
    }

    [Test]
    public void ClassifyComplexRearrangement_ConfidenceTierBoundaries_HighLowNone()
    {
        // INV-03 confidence boundaries by oscillating-segment count (transitions+1):
        //   3 oscillations ⇒ 4 segments ⇒ Low (≥4);
        //   6 oscillations ⇒ 7 segments ⇒ High (≥7);
        //   2 oscillations ⇒ 3 segments ⇒ None (<4).
        var low = ClassifyComplexRearrangement(
            new ComplexRearrangementInput(AlternatingProfile(4), StructuralVariantCount: 0));
        var high = ClassifyComplexRearrangement(
            new ComplexRearrangementInput(AlternatingProfile(7), StructuralVariantCount: 0));
        var none = ClassifyComplexRearrangement(
            new ComplexRearrangementInput(AlternatingProfile(3), StructuralVariantCount: 0));

        low.OscillatingSegmentCount.Should().Be(4);
        low.Confidence.Should().Be(ChromothripsisConfidence.Low);
        high.OscillatingSegmentCount.Should().Be(7);
        high.Confidence.Should().Be(ChromothripsisConfidence.High);
        none.OscillatingSegmentCount.Should().Be(3);
        none.Confidence.Should().Be(ChromothripsisConfidence.None);
    }

    #endregion

    #region ONCO-SV-001 — Null guards (documented ArgumentNullException)

    [Test]
    public void CountCopyNumberStateOscillations_NullProfile_ThrowsArgumentNull()
    {
        var act = () => CountCopyNumberStateOscillations(null!);

        act.Should().Throw<ArgumentNullException>();
    }

    [Test]
    public void TestBreakpointClustering_NullPositions_ThrowsArgumentNull()
    {
        var act = () => TestBreakpointClustering(null!);

        act.Should().Throw<ArgumentNullException>();
    }

    [Test]
    public void ClassifyComplexRearrangement_NullSegments_ThrowsArgumentNull()
    {
        var act = () => ClassifyComplexRearrangement(
            new ComplexRearrangementInput(null!, StructuralVariantCount: 6));

        act.Should().Throw<ArgumentNullException>();
    }

    #endregion

    #region ONCO-SV-001 — BE: broad random fuzz (always well-formed, terminates)

    [Test]
    [CancelAfter(30000)]
    public void ClassifyComplexRearrangement_RandomProfiles_AlwaysWellFormed(
        [Values(20260620, 1234, 98765)] int seed)
    {
        // Fuzz random CN profiles of mixed length (incl. 0 and 1) over a small CN
        // alphabet, with random SV burdens (incl. 0 and MaxInt). Every result must
        // satisfy the structural invariants (INV-01..04) and the call must never
        // throw or hang — the contract holds across the whole input space.
        var rng = new Random(seed);
        for (int t = 0; t < 20_000; t++)
        {
            int n = rng.Next(0, 30); // includes empty (0) and single (1)
            var profile = new int[n];
            int alphabet = rng.Next(1, 6); // 1..5 distinct CN states available
            for (int i = 0; i < n; i++)
            {
                profile[i] = rng.Next(alphabet);
            }

            int svBurden = rng.Next(10) == 0 ? int.MaxValue : rng.Next(0, 30);

            var result = ClassifyComplexRearrangement(
                new ComplexRearrangementInput(profile, svBurden));

            AssertWellFormedResult(result, n);
            result.StructuralVariantCount.Should().Be(svBurden, "SV burden is echoed unchanged");
        }
    }

    [Test]
    [CancelAfter(30000)]
    public void TestBreakpointClustering_RandomBreakpointSets_AlwaysWellFormed(
        [Values(7, 555, 31337)] int seed)
    {
        // Fuzz random breakpoint sets of mixed size (incl. 0, 1, 2 — the sub-3
        // not-assessable cases) over a wide coordinate range. The CV must never be
        // NaN/∞ and IsClustered must be the exact CV > 1 predicate; sets with < 3
        // breakpoints are always not-clustered (§3.3, INV-05).
        var rng = new Random(seed);
        for (int t = 0; t < 5_000; t++)
        {
            int n = rng.Next(0, 12); // includes 0, 1, 2 and small clustered sets
            var positions = new long[n];
            for (int i = 0; i < n; i++)
            {
                positions[i] = rng.NextInt64(0, 10_000_000);
            }

            var result = TestBreakpointClustering(positions);

            result.BreakpointCount.Should().Be(n);
            double.IsNaN(result.CoefficientOfVariation).Should().BeFalse();
            double.IsInfinity(result.CoefficientOfVariation).Should().BeFalse();
            result.MeanGap.Should().BeGreaterThanOrEqualTo(0.0);

            if (n < 3)
            {
                result.IsClustered.Should().BeFalse("< 3 breakpoints ⇒ clustering undefined (§3.3)");
            }
            else if (result.MeanGap > 0.0)
            {
                result.IsClustered.Should().Be(
                    result.CoefficientOfVariation > 1.0,
                    "IsClustered is exactly the CV>1 predicate (INV-05)");
            }
        }
    }

    #endregion
}
