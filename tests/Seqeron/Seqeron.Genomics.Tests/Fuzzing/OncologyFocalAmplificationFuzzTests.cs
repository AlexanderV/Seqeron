using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Oncology;
using static Seqeron.Genomics.Oncology.OncologyAnalyzer;

namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the Oncology focal-amplification calling area — ONCO-CNA-002.
/// The unit under test is the GISTIC2 focal/broad length-rule filter in
/// src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs:
///   • <see cref="OncologyAnalyzer.IsFocalAmplification"/> — single-segment predicate
///     (amplified ∧ focal);
///   • <see cref="OncologyAnalyzer.DetectFocalAmplifications"/> — order-preserving
///     filter selecting the focal amplifications from a segment stream;
///   • <see cref="OncologyAnalyzer.IdentifyAmplifiedOncogenes"/> — arm→oncogene panel
///     map over focal amplifications.
///
/// This is the focal-vs-broad member of the CNA family (rows 103–105). ONCO-CNA-001
/// (row 103, log2→copy-number classification) supplies the per-bin amplitude state;
/// THIS unit operates one level up, on already-segmented arm-anchored intervals
/// (<see cref="OncologyAnalyzer.CopyNumberArmSegment"/>), and decides — purely by
/// LENGTH relative to the chromosome arm, plus an amplitude gate — whether an
/// amplified segment is a *focal* (narrow, therapeutically actionable) event or a
/// *broad / arm-level* (genome-wide-scale) event. Segmentation itself is upstream
/// (docs §5.2 — StructuralVariantAnalyzer.SegmentCopyNumber, SV-CNV-001), so the
/// "single-bin focal" target here is exercised as a width-1 arm segment, the
/// degenerate analogue of a one-bin focal call.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate / boundary / malformed inputs to a unit and asserts
/// that the code NEVER fails in an undisciplined way: no hang, no nonsense output,
/// no *unhandled* runtime exception (a DivideByZero on a width-1 / zero-arm-length
/// segment, a NaN/±∞ leaking through the amplitude or length test, a false focal
/// call on a genome-wide event, or an off-by-one at the strict length/amplitude
/// cutoff). Every input must resolve to EITHER a well-defined, theory-correct
/// outcome OR a *documented, intentional* one (ArgumentNullException for a null
/// stream, ArgumentException for ArmLength ≤ 0 or End ≤ Start). The headline
/// hazards for the GISTIC2 length rule are:
///   • a genome-wide / arm-level amplification (length ≥ 98% of the arm, any
///     amplitude) must NOT be reported as focal — the broad event is filtered out,
///     never a false focal call (§6.1, INV-01);
///   • a single-bin (width-1) focal segment must classify without a
///     DivideByZero/length crash and, if it clears t_amp, be reported as focal
///     (ArmLength &gt; 0 is validated; the width-1 fraction 1/ArmLength is &lt; 0.98);
///   • the strict boundaries: a segment occupying EXACTLY 0.98 of the arm is
///     arm-level (focal test is `&lt; 0.98`, not `≤`), and a segment at log2 EXACTLY
///     0.1 is NOT amplified (amplitude test is `&gt; 0.1`, not `≥`) — no off-by-one;
///   • NaN log2 / NaN-producing arm fractions must not slip through an `&gt;`/`&lt;`
///     comparison as a false positive (NaN comparisons are false ⇒ not focal).
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: ONCO-CNA-002 — Focal Amplification Detection (Oncology)
/// Checklist: docs/checklists/03_FUZZING.md, row 104.
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — граничні значення: 0, -1, MaxInt, empty.
///     Targets (checklist row 104): "genome-wide amp, single-bin focal, threshold edge".
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The documented contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// Focal_Amplification_Detection.md
/// (docs/algorithms/Oncology/Focal_Amplification_Detection.md):
///   • Amplified ⇔ r &gt; t_amp, t_amp = 0.1 (strict)                     (§2.2, INV-02)
///   • Focal ⇔ (L / A) &lt; c, c = broad_len_cutoff = 0.98 (strict)         (§2.2, INV-01)
///   • Focal amplification ⇔ Amplified ∧ Focal                          (§2.2)
///   • segment exactly 0.98 of the arm ⇒ arm-level (NOT focal)          (§3.3, §6.1)
///   • log2 ≤ t_amp ⇒ not amplified                                     (§6.1)
///   • output is a subset of the input in input order (a filter)        (INV-03)
///   • an oncogene is reported only for an arm carrying a focal amp      (INV-04)
///   • arm matching is case-insensitive (Ordinal-ignore-case)           (§3.3)
///   • null segments/amplifications ⇒ ArgumentNullException             (§3.3, §6.1)
///   • empty input ⇒ empty result                                       (§6.1)
///   • ArmLength ≤ 0 or End ≤ Start ⇒ ArgumentException                 (§3.3, §6.1)
///   • worked example: 17q [100k,600k) of 1M arm, log2 1.0 ⇒ focal;
///     8q [0,990k) of 1M arm, log2 1.5 ⇒ arm-level                      (§7.1)
///
/// All randomness is LOCALLY seeded (new Random(seed)); no shared static Rng.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public sealed class OncologyFocalAmplificationFuzzTests
{
    // Fixed reference arm length so a segment's [Start,End) span maps directly to a
    // fraction of the arm (e.g. [0, 500_000) of a 1_000_000 arm ⇒ 0.50).
    private const long ArmLen = 1_000_000L;

    private static CopyNumberArmSegment Seg(
        string arm, long start, long end, double log2, long armLength = ArmLen) =>
        new(arm, start, end, armLength, log2);

    // ── Well-formed-result assertion helper ──────────────────────────────────
    // Pins the documented filter contract on EVERY accepted detection: the result
    // is a SUBSET of the input in input order (INV-03), and — the load-bearing
    // theory check — every reported segment genuinely satisfies BOTH the strict
    // amplitude gate (log2 > t_amp) AND the strict focal-length gate
    // (ArmFraction < broad_len_cutoff). This is what stops a fuzz test from
    // rubber-stamping a broad/arm-level segment as focal.
    private static void AssertWellFormedDetection(
        IReadOnlyList<CopyNumberArmSegment> result,
        IReadOnlyList<CopyNumberArmSegment> input,
        FocalAmplificationThresholds cutoffs)
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
            seg.Log2Ratio.Should().BeGreaterThan(
                cutoffs.AmplificationLog2Threshold,
                "every reported amplification strictly exceeds t_amp (INV-02)");
            seg.ArmFraction.Should().BeLessThan(
                cutoffs.BroadLengthCutoff,
                "every reported amplification is strictly below the broad-length cutoff (INV-01)");
        }
    }

    #region ONCO-CNA-002 — Positive sanity (documented focal / broad / sub-amplitude)

    [Test]
    public void DetectFocalAmplifications_DocumentedWorkedExample_MatchesDocSpec()
    {
        // Docs §7.1: a 17q segment [100k,600k) of a 1M arm (fraction 0.50 < 0.98)
        // at log2 1.0 (> 0.1) is a focal amplification; an 8q segment [0,990k)
        // (fraction 0.99 ≥ 0.98) at log2 1.5 is arm-level. Only the 17q is focal.
        var input = new[]
        {
            Seg("17q", 100_000, 600_000, 1.0), // focal amp
            Seg("8q", 0, 990_000, 1.5),         // arm-level (broad)
        };

        var focal = DetectFocalAmplifications(input);

        focal.Should().ContainSingle();
        focal[0].Arm.Should().Be("17q");
        IdentifyAmplifiedOncogenes(focal).Should().Equal("ERBB2");
        AssertWellFormedDetection(focal, input, FocalAmplificationThresholds.Default);
    }

    [Test]
    public void DetectFocalAmplifications_NarrowHighAmp_IsCalledFocal()
    {
        // POSITIVE sanity: a narrow (0.10 of arm), high-amplitude (log2 2.0) segment
        // is the canonical focal amplification — amplified AND focal.
        var input = new[] { Seg("7p", 0, 100_000, 2.0) };

        var focal = DetectFocalAmplifications(input);

        focal.Should().ContainSingle();
        IsFocalAmplification(input[0], FocalAmplificationThresholds.Default).Should().BeTrue();
    }

    [Test]
    public void DetectFocalAmplifications_GenomeWideHighAmp_IsNotFocal()
    {
        // POSITIVE sanity (the BE "genome-wide amp" target): a segment spanning the
        // ENTIRE arm at very high amplitude (log2 5.0) is broad/arm-level, NOT focal
        // (fraction 1.0 ≥ 0.98) — amplitude alone never makes a genome-wide event focal.
        var input = new[] { Seg("8q", 0, ArmLen, 5.0) };

        var focal = DetectFocalAmplifications(input);

        focal.Should().BeEmpty("a whole-arm span is arm-level regardless of amplitude (INV-01)");
        IsFocalAmplification(input[0], FocalAmplificationThresholds.Default).Should().BeFalse();
    }

    [Test]
    public void DetectFocalAmplifications_SubAmplitudeNarrowSegment_IsNotFocal()
    {
        // POSITIVE sanity: a narrow segment (0.10 of arm, focal by length) whose
        // amplitude is at/below t_amp is NOT a focal amplification (the amplitude
        // gate rejects low-level segments). log2 0.05 < 0.1.
        var input = new[] { Seg("17q", 0, 100_000, 0.05) };

        var focal = DetectFocalAmplifications(input);

        focal.Should().BeEmpty("sub-t_amp amplitude is not amplified (INV-02)");
        IsFocalAmplification(input[0], FocalAmplificationThresholds.Default).Should().BeFalse();
    }

    #endregion

    #region ONCO-CNA-002 — BE: genome-wide / arm-level amp (no false focal call)

    [Test]
    public void DetectFocalAmplifications_ExactlyBroadCutoff_IsArmLevel_StrictLessThan()
    {
        // THRESHOLD EDGE (length): a segment occupying EXACTLY 0.98 of the arm is
        // arm-level — the focal test is strict `< 0.98`, so fraction == 0.98 fails.
        // [0, 980_000) of 1_000_000 = exactly 0.98. Guards an off-by-one (`≤`) bug.
        var input = new[] { Seg("17q", 0, 980_000, 1.0) };

        var focal = DetectFocalAmplifications(input);

        focal.Should().BeEmpty("fraction exactly at broad_len_cutoff is arm-level (§3.3/§6.1)");
        input[0].ArmFraction.Should().Be(0.98);
    }

    [Test]
    public void DetectFocalAmplifications_JustBelowBroadCutoff_IsFocal()
    {
        // Just inside the focal side: 0.979 of the arm (< 0.98) at log2 1.0 IS focal.
        // Pairs with the exactly-0.98 case to pin the boundary direction.
        var input = new[] { Seg("17q", 0, 979_000, 1.0) };

        var focal = DetectFocalAmplifications(input);

        focal.Should().ContainSingle();
        input[0].ArmFraction.Should().BeLessThan(0.98);
    }

    [Test]
    public void DetectFocalAmplifications_OverWholeArm_IsArmLevel()
    {
        // A segment LONGER than its arm (fraction > 1.0; an over-segmented call) is
        // unambiguously arm-level — no crash, simply not focal.
        var input = new[] { Seg("8q", 0, 2_000_000, 3.0) };

        var focal = DetectFocalAmplifications(input);

        focal.Should().BeEmpty("fraction > 1.0 is broad (INV-01)");
    }

    [Test]
    public void DetectFocalAmplifications_GenomeWideMixedWithFocal_KeepsOnlyFocal()
    {
        // A genome-wide arm-level amp interleaved with a true focal amp: only the
        // focal one survives the filter, in input order (no false focal on the broad).
        var input = new[]
        {
            Seg("8q", 0, ArmLen, 4.0),          // genome-wide ⇒ broad
            Seg("17q", 100_000, 300_000, 1.0),  // focal
            Seg("11q", 0, 990_000, 2.0),        // 0.99 ⇒ broad
        };

        var focal = DetectFocalAmplifications(input);

        focal.Should().ContainSingle();
        focal[0].Arm.Should().Be("17q");
        AssertWellFormedDetection(focal, input, FocalAmplificationThresholds.Default);
    }

    #endregion

    #region ONCO-CNA-002 — BE: single-bin focal (width-1 segment, no length crash)

    [Test]
    public void DetectFocalAmplifications_SingleBinWidthOne_IsFocal_NoLengthCrash()
    {
        // SINGLE-BIN FOCAL: a width-1 segment (End = Start + 1, Length = 1). The
        // hazard is a DivideByZero / length-0 crash; ArmFraction = 1/ArmLength is a
        // clean tiny positive fraction (< 0.98), so a width-1 amplified segment is
        // focal. No throw, reported as focal.
        var input = new[] { Seg("17q", 500_000, 500_001, 1.0) };

        var act = () => DetectFocalAmplifications(input);

        act.Should().NotThrow();
        var focal = act();
        focal.Should().ContainSingle();
        input[0].Length.Should().Be(1);
        input[0].ArmFraction.Should().BeApproximately(1.0 / ArmLen, 1e-15);
    }

    [Test]
    public void IsFocalAmplification_SingleBinWidthOne_DoesNotThrow()
    {
        // The single-segment predicate on a width-1 segment: must evaluate cleanly.
        var seg = Seg("7p", 0, 1, 1.0);

        var act = () => IsFocalAmplification(seg, FocalAmplificationThresholds.Default);

        act.Should().NotThrow();
        act().Should().BeTrue("width-1, high-amplitude ⇒ focal");
    }

    [Test]
    public void DetectFocalAmplifications_SingleBinOnArmLengthOne_IsArmLevel_NoCrash()
    {
        // Degenerate arm of length 1 with a width-1 segment ⇒ fraction = 1/1 = 1.0,
        // i.e. the segment IS the whole arm ⇒ arm-level. Still no DivideByZero (arm
        // length is the validated positive divisor).
        var input = new[] { Seg("8q", 0, 1, 2.0, armLength: 1) };

        var act = () => DetectFocalAmplifications(input);

        act.Should().NotThrow();
        act().Should().BeEmpty("fraction 1/1 = 1.0 ≥ 0.98 ⇒ arm-level");
    }

    #endregion

    #region ONCO-CNA-002 — BE: amplitude threshold edge (strict >, no off-by-one)

    [Test]
    public void DetectFocalAmplifications_Log2ExactlyAtTAmp_IsNotAmplified_StrictGreaterThan()
    {
        // THRESHOLD EDGE (amplitude): log2 EXACTLY 0.1 is NOT amplified — the test
        // is strict `> 0.1`. Narrow segment so length is not the disqualifier.
        var input = new[] { Seg("17q", 0, 100_000, DefaultAmplificationLog2Threshold) };

        var focal = DetectFocalAmplifications(input);

        focal.Should().BeEmpty("log2 exactly at t_amp is not amplified (strict >, §6.1)");
    }

    [Test]
    public void DetectFocalAmplifications_Log2JustAboveTAmp_IsAmplified()
    {
        // Just above the amplitude cutoff (0.1 + tiny epsilon) on a focal-length
        // segment IS focal — pairs with the exactly-0.1 case to pin the direction.
        double justAbove = Math.BitIncrement(DefaultAmplificationLog2Threshold);
        var input = new[] { Seg("17q", 0, 100_000, justAbove) };

        var focal = DetectFocalAmplifications(input);

        focal.Should().ContainSingle();
        justAbove.Should().BeGreaterThan(DefaultAmplificationLog2Threshold);
    }

    [Test]
    public void DetectFocalAmplifications_NegativeLog2_IsNotAmplified()
    {
        // A focal-length DELETION-direction segment (log2 < 0) is not an
        // amplification at all (deletions are ONCO-CNA-003 scope) — never reported.
        var input = new[] { Seg("17q", 0, 100_000, -1.5) };

        DetectFocalAmplifications(input).Should().BeEmpty();
    }

    #endregion

    #region ONCO-CNA-002 — BE: NaN / ±∞ log2 (no false focal via comparison)

    [Test]
    public void IsFocalAmplification_NaNLog2_IsNotFocal_NoComparisonLeak()
    {
        // A NaN log2 makes `NaN > t_amp` false ⇒ not amplified ⇒ not focal. The
        // hazard is a NaN slipping through an inverted comparison as a false focal.
        var seg = Seg("17q", 0, 100_000, double.NaN);

        IsFocalAmplification(seg, FocalAmplificationThresholds.Default)
            .Should().BeFalse("NaN amplitude cannot exceed t_amp (no false focal)");
    }

    [Test]
    public void DetectFocalAmplifications_PositiveInfinityLog2_FocalLength_IsFocal()
    {
        // log2 = +∞ is unboundedly amplified (+∞ > 0.1) and, on a focal-length
        // segment, IS focal — the amplitude gate admits it without overflow/crash.
        var input = new[] { Seg("17q", 0, 100_000, double.PositiveInfinity) };

        var focal = DetectFocalAmplifications(input);

        focal.Should().ContainSingle();
    }

    [Test]
    public void DetectFocalAmplifications_PositiveInfinityLog2_WholeArm_IsArmLevel()
    {
        // Even +∞ amplitude does NOT rescue a genome-wide span: fraction 1.0 ≥ 0.98
        // ⇒ arm-level. Amplitude and length are independent gates (both must pass).
        var input = new[] { Seg("8q", 0, ArmLen, double.PositiveInfinity) };

        DetectFocalAmplifications(input).Should().BeEmpty();
    }

    [Test]
    public void IsFocalAmplification_NaNArmFraction_IsNotFocal()
    {
        // Construct a NaN arm fraction: Length 0 over... not possible (End>Start
        // validated). Instead a huge length such that the ratio is finite; this test
        // pins that a non-finite fraction (only reachable via overflow) never passes
        // the strict `< 0.98` focal test. long.MaxValue length over a tiny arm ⇒ a
        // huge finite fraction (>> 0.98) ⇒ not focal.
        var seg = Seg("17q", 0, long.MaxValue, 1.0, armLength: 1);

        IsFocalAmplification(seg, FocalAmplificationThresholds.Default)
            .Should().BeFalse("an astronomically large fraction is broad, never focal");
    }

    #endregion

    #region ONCO-CNA-002 — BE: empty / null / malformed segment (documented guards)

    [Test]
    public void DetectFocalAmplifications_EmptyInput_ReturnsEmpty_NoThrow()
    {
        // Empty stream ⇒ empty result (§6.1) — the zero-segment boundary.
        DetectFocalAmplifications(Array.Empty<CopyNumberArmSegment>()).Should().BeEmpty();
    }

    [Test]
    public void DetectFocalAmplifications_NullSegments_ThrowsArgumentNull()
    {
        // Documented null guard (§3.3) — not a NullReference deep in the loop.
        var act = () => DetectFocalAmplifications(null!);

        act.Should().Throw<ArgumentNullException>();
    }

    [Test]
    public void IdentifyAmplifiedOncogenes_NullAmplifications_ThrowsArgumentNull()
    {
        var act = () => IdentifyAmplifiedOncogenes(null!);

        act.Should().Throw<ArgumentNullException>();
    }

    [Test]
    public void DetectFocalAmplifications_NonPositiveArmLength_ThrowsArgument()
    {
        // ArmLength ≤ 0 is the divisor hazard for ArmFraction; the validator rejects
        // it up front (§3.3/§6.1) rather than producing a NaN/∞ fraction.
        var zeroArm = new[] { Seg("17q", 0, 100_000, 1.0, armLength: 0) };
        var negArm = new[] { Seg("17q", 0, 100_000, 1.0, armLength: -1) };

        ((Action)(() => DetectFocalAmplifications(zeroArm))).Should().Throw<ArgumentException>();
        ((Action)(() => DetectFocalAmplifications(negArm))).Should().Throw<ArgumentException>();
    }

    [Test]
    public void DetectFocalAmplifications_EndNotAfterStart_ThrowsArgument()
    {
        // End ≤ Start (zero/negative-width segment) ⇒ ArgumentException (§3.3/§6.1).
        var zeroWidth = new[] { Seg("17q", 100_000, 100_000, 1.0) };
        var inverted = new[] { Seg("17q", 100_000, 50_000, 1.0) };

        ((Action)(() => DetectFocalAmplifications(zeroWidth))).Should().Throw<ArgumentException>();
        ((Action)(() => DetectFocalAmplifications(inverted))).Should().Throw<ArgumentException>();
    }

    #endregion

    #region ONCO-CNA-002 — BE: oncogene mapping (case-insensitive, panel order, subset)

    [Test]
    public void IdentifyAmplifiedOncogenes_EmptyAmplifications_ReturnsEmpty()
    {
        IdentifyAmplifiedOncogenes(Array.Empty<CopyNumberArmSegment>()).Should().BeEmpty();
    }

    [Test]
    public void IdentifyAmplifiedOncogenes_ArmCaseInsensitive_MatchesPanel()
    {
        // §3.3: arm matching is Ordinal-ignore-case. An "17Q" focal amp still flags ERBB2.
        var amps = new[] { Seg("17Q", 0, 100_000, 1.0) };

        IdentifyAmplifiedOncogenes(amps).Should().Equal("ERBB2");
    }

    [Test]
    public void IdentifyAmplifiedOncogenes_SharedArm_ReportsBothPanelGenesInOrder()
    {
        // 12q carries both MDM2 and CDK4; a single 12q focal amp reports both, in
        // panel order (MDM2 before CDK4), each once (INV-04, distinct).
        var amps = new[] { Seg("12q", 0, 100_000, 1.0) };

        IdentifyAmplifiedOncogenes(amps).Should().Equal("MDM2", "CDK4");
    }

    [Test]
    public void IdentifyAmplifiedOncogenes_NonPanelArm_ReportsNothing()
    {
        // An amplified arm outside the six-gene panel maps to no oncogene.
        var amps = new[] { Seg("1p", 0, 100_000, 1.0) };

        IdentifyAmplifiedOncogenes(amps).Should().BeEmpty();
    }

    #endregion

    #region ONCO-CNA-002 — BE: broad random fuzz (well-formed, gates independent)

    [Test]
    [CancelAfter(30000)]
    public void DetectFocalAmplifications_RandomSegments_AlwaysWellFormedFilter(
        [Values(20260104, 555, 909090)] int seed)
    {
        // Fuzz mixed streams of arm segments spanning the focal/broad boundary and
        // the amplitude boundary, including width-1 and whole-arm spans. Every
        // result must be an order-preserving subset whose members all pass BOTH
        // strict gates — no broad/sub-amplitude leak, no crash.
        var rng = new Random(seed);
        var arms = new[] { "17q", "8q", "7p", "11q", "12q", "1p" };
        var cutoffs = FocalAmplificationThresholds.Default;

        for (int t = 0; t < 5000; t++)
        {
            int n = rng.Next(0, 8);
            var input = new List<CopyNumberArmSegment>(n);
            for (int i = 0; i < n; i++)
            {
                string arm = arms[rng.Next(arms.Length)];
                // Length from 1 bp (single-bin) up to slightly over the whole arm.
                long len = (long)(rng.NextDouble() * (ArmLen + 50_000)) + 1;
                long start = rng.Next(0, 50_000);
                double log2 = (rng.NextDouble() * 2.0 - 1.0) * 4.0; // span ±4 around 0
                input.Add(Seg(arm, start, start + len, log2));
            }

            var focal = DetectFocalAmplifications(input);

            AssertWellFormedDetection(focal, input, cutoffs);

            // Cross-check the predicate: a segment is in the result iff the predicate says so.
            var expected = input.Where(s => IsFocalAmplification(s, cutoffs)).ToList();
            focal.Should().Equal(expected, "filter agrees with the single-segment predicate");
        }
    }

    [Test]
    [CancelAfter(30000)]
    public void IsFocalAmplification_GatesAreIndependentAndStrict(
        [Values(7, 4242, 31337)] int seed)
    {
        // Property: focal ⇔ (log2 > t_amp) ∧ (fraction < cutoff). Fuzz the two axes
        // independently and assert the predicate equals the conjunction of the two
        // strict comparisons — a single combined off-by-one would break this.
        var rng = new Random(seed);
        var cutoffs = FocalAmplificationThresholds.Default;

        for (int i = 0; i < 20000; i++)
        {
            double log2 = (rng.NextDouble() * 2.0 - 1.0) * 3.0;
            long len = (long)(rng.NextDouble() * (ArmLen + 20_000)) + 1;
            var seg = Seg("17q", 0, len, log2);

            bool amplified = log2 > cutoffs.AmplificationLog2Threshold;
            bool focal = (double)len / ArmLen < cutoffs.BroadLengthCutoff;

            IsFocalAmplification(seg, cutoffs).Should().Be(amplified && focal);
        }
    }

    #endregion
}
