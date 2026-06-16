// ONCO-SV-001 — Somatic Complex Rearrangement Classification (Chromothripsis Inference)
// Evidence: docs/Evidence/ONCO-SV-001-Evidence.md
// TestSpec: tests/TestSpecs/ONCO-SV-001.md
// Source: Korbel JO, Campbell PJ (2013). Cell 152(6):1226-1236. doi:10.1016/j.cell.2013.02.023
//         Cortés-Ciriano I et al. (2020). Nat Genet 52:331-341. doi:10.1038/s41588-019-0576-7
//         Magrangeas F et al. (2011). Blood 118(3):675-678 (>=10 oscillating-CN-change screen).

using System;
using System.Collections.Generic;
using NUnit.Framework;
using Seqeron.Genomics.Oncology;
using Type = Seqeron.Genomics.Oncology.OncologyAnalyzer.ComplexRearrangementType;
using Confidence = Seqeron.Genomics.Oncology.OncologyAnalyzer.ChromothripsisConfidence;
using Input = Seqeron.Genomics.Oncology.OncologyAnalyzer.ComplexRearrangementInput;

namespace Seqeron.Genomics.Tests;

[TestFixture]
public class OncologyAnalyzer_ClassifyComplexRearrangement_Tests
{
    private const double Tolerance = 1e-10;

    #region CountCopyNumberStateOscillations

    // M1 — Two-state alternating profile 2,1,...,2 (11 segments) has 10 adjacent state transitions
    //      (the >=10 first-pass oscillating-CN-change screen; Magrangeas 2011 / Korbel & Campbell 2013).
    [Test]
    public void CountCopyNumberStateOscillations_TwoStateAlternating_Returns10()
    {
        int[] states = { 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2 };
        int count = OncologyAnalyzer.CountCopyNumberStateOscillations(states);
        Assert.That(count, Is.EqualTo(10),
            "11-segment 2,1,...,2 profile has 10 adjacent CN-state transitions (oscillating CN changes).");
    }

    // M2 — All-equal profile has zero transitions (INV-1 lower bound).
    [Test]
    public void CountCopyNumberStateOscillations_AllEqual_ReturnsZero()
    {
        int[] states = { 2, 2, 2, 2 };
        int count = OncologyAnalyzer.CountCopyNumberStateOscillations(states);
        Assert.That(count, Is.EqualTo(0),
            "A constant copy-number profile has no state transitions.");
    }

    // M3 — Monotone rising profile: every adjacent pair differs, so n-1 transitions (here 4).
    [Test]
    public void CountCopyNumberStateOscillations_MonotoneRising_ReturnsNMinus1()
    {
        int[] states = { 2, 3, 4, 5, 6 };
        int count = OncologyAnalyzer.CountCopyNumberStateOscillations(states);
        Assert.That(count, Is.EqualTo(4),
            "Strictly increasing 5-segment profile has 4 transitions (each adjacent pair differs).");
    }

    // C2 — Empty and single-segment profiles have zero transitions (INV-1 boundary).
    [Test]
    public void CountCopyNumberStateOscillations_EmptyOrSingle_ReturnsZero()
    {
        Assert.Multiple(() =>
        {
            Assert.That(OncologyAnalyzer.CountCopyNumberStateOscillations(Array.Empty<int>()), Is.EqualTo(0),
                "Empty profile has no adjacent pairs, so 0 transitions.");
            Assert.That(OncologyAnalyzer.CountCopyNumberStateOscillations(new[] { 2 }), Is.EqualTo(0),
                "Single-segment profile has no adjacent pairs, so 0 transitions.");
        });
    }

    // C1 — Null profile is rejected.
    [Test]
    public void CountCopyNumberStateOscillations_Null_Throws()
    {
        Assert.Throws<ArgumentNullException>(
            () => OncologyAnalyzer.CountCopyNumberStateOscillations(null!),
            "Null segment list is invalid input.");
    }

    #endregion

    #region ClassifyComplexRearrangement

    // M4 — Canonical two-state oscillation: 11 segments (10 oscillations), 2 states, 12 SVs ->
    //      Chromothripsis (Korbel & Campbell criterion B + >=10 screen + >=6 SV burden).
    [Test]
    public void ClassifyComplexRearrangement_TwoStateTenOscillations_IsChromothripsis()
    {
        var input = new Input(new[] { 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2 }, 12);
        var r = OncologyAnalyzer.ClassifyComplexRearrangement(input);
        Assert.Multiple(() =>
        {
            Assert.That(r.Type, Is.EqualTo(Type.Chromothripsis),
                "Two-state profile with 10 oscillations and 12 SVs meets all hallmark gates.");
            Assert.That(r.OscillationCount, Is.EqualTo(10), "10 adjacent CN-state transitions.");
            Assert.That(r.DistinctStateCount, Is.EqualTo(2), "Profile uses exactly states {1,2}.");
        });
    }

    // M5 — Below the >=10 first-pass screen: 6 segments -> 5 oscillations -> NotComplex.
    [Test]
    public void ClassifyComplexRearrangement_FiveOscillations_NotComplex()
    {
        var input = new Input(new[] { 2, 1, 2, 1, 2, 1 }, 12);
        var r = OncologyAnalyzer.ClassifyComplexRearrangement(input);
        Assert.Multiple(() =>
        {
            Assert.That(r.OscillationCount, Is.EqualTo(5), "Six-segment alternation has 5 transitions.");
            Assert.That(r.Type, Is.EqualTo(Type.NotComplex),
                "5 oscillations is below the >=10 first-pass chromothripsis screen (Magrangeas 2011).");
        });
    }

    // M6 — Progressive amplification: monotone rising 11 segments (10 oscillations) but 11 distinct
    //      states -> NotComplex (two-state hallmark excludes many-state ascending profiles; INV-4).
    [Test]
    public void ClassifyComplexRearrangement_MonotoneManyStates_NotComplex()
    {
        var input = new Input(new[] { 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 }, 12);
        var r = OncologyAnalyzer.ClassifyComplexRearrangement(input);
        Assert.Multiple(() =>
        {
            Assert.That(r.OscillationCount, Is.EqualTo(10), "Each adjacent pair differs: 10 transitions.");
            Assert.That(r.DistinctStateCount, Is.EqualTo(11), "All 11 copy numbers are distinct.");
            Assert.That(r.Type, Is.EqualTo(Type.NotComplex),
                "More than 3 distinct states is progressive amplification, not the two-state chromothripsis hallmark.");
        });
    }

    // M7 — Valid two-state 10-oscillation profile but only 5 SVs -> excluded by the <6 SV focal floor.
    [Test]
    public void ClassifyComplexRearrangement_SvBurdenBelowSix_NotComplex()
    {
        var input = new Input(new[] { 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2 }, 5);
        var r = OncologyAnalyzer.ClassifyComplexRearrangement(input);
        Assert.That(r.Type, Is.EqualTo(Type.NotComplex),
            "Events with fewer than 6 SVs are excluded (Cortés-Ciriano 2020), even with 10 oscillations.");
    }

    // M8 — High confidence: 11 oscillating segments (>=7) -> High (Cortés-Ciriano 2020 segment threshold).
    [Test]
    public void ClassifyComplexRearrangement_ElevenOscillatingSegments_HighConfidence()
    {
        var input = new Input(new[] { 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2 }, 12);
        var r = OncologyAnalyzer.ClassifyComplexRearrangement(input);
        Assert.Multiple(() =>
        {
            Assert.That(r.OscillatingSegmentCount, Is.EqualTo(11),
                "10 transitions span 11 oscillating segments.");
            Assert.That(r.Confidence, Is.EqualTo(Confidence.High),
                ">=7 adjacent oscillating segments is high-confidence (Cortés-Ciriano 2020).");
        });
    }

    // S1 — Low-confidence tier: 6 oscillating segments (in [4,6]) but only 5 oscillations -> Low + NotComplex.
    [Test]
    public void ClassifyComplexRearrangement_SixOscillatingSegments_LowConfidenceNotComplex()
    {
        var input = new Input(new[] { 2, 1, 2, 1, 2, 1 }, 12); // 5 transitions -> 6 oscillating segments
        var r = OncologyAnalyzer.ClassifyComplexRearrangement(input);
        Assert.Multiple(() =>
        {
            Assert.That(r.OscillatingSegmentCount, Is.EqualTo(6), "5 transitions span 6 oscillating segments.");
            Assert.That(r.Confidence, Is.EqualTo(Confidence.Low),
                "4-6 adjacent oscillating segments is the low-confidence tier (Cortés-Ciriano 2020).");
            Assert.That(r.Type, Is.EqualTo(Type.NotComplex),
                "Low-confidence signal: still below the >=10 oscillation screen.");
        });
    }

    // S2 — Three-state oscillation {1,2,3} with 10 oscillations and 12 SVs -> Chromothripsis
    //      (the "two-or-three state" tolerance of criterion B; <=3 distinct states is allowed).
    [Test]
    public void ClassifyComplexRearrangement_ThreeStateOscillation_IsChromothripsis()
    {
        // 2,1,2,3,2,1,2,3,2,1,2 -> states {1,2,3}=3 distinct, 10 transitions.
        var input = new Input(new[] { 2, 1, 2, 3, 2, 1, 2, 3, 2, 1, 2 }, 12);
        var r = OncologyAnalyzer.ClassifyComplexRearrangement(input);
        Assert.Multiple(() =>
        {
            Assert.That(r.DistinctStateCount, Is.EqualTo(3), "Profile uses states {1,2,3}.");
            Assert.That(r.OscillationCount, Is.EqualTo(10), "10 adjacent transitions.");
            Assert.That(r.Type, Is.EqualTo(Type.Chromothripsis),
                "Up to 3 distinct oscillating states is permitted (Korbel & Campbell two-or-three hallmark).");
        });
    }

    // C4 — Null segment list inside the input is rejected.
    [Test]
    public void ClassifyComplexRearrangement_NullSegments_Throws()
    {
        var input = new Input(null!, 12);
        Assert.Throws<ArgumentNullException>(
            () => OncologyAnalyzer.ClassifyComplexRearrangement(input),
            "Null segment list is invalid input.");
    }

    #endregion

    #region TestBreakpointClustering

    // M9 — Regularly spaced breakpoints: equal gaps -> CV = 0 -> not clustered (exponential null CV=1).
    [Test]
    public void TestBreakpointClustering_RegularSpacing_NotClustered()
    {
        long[] positions = { 100, 200, 300, 400, 500 };
        var r = OncologyAnalyzer.TestBreakpointClustering(positions);
        Assert.Multiple(() =>
        {
            Assert.That(r.CoefficientOfVariation, Is.EqualTo(0.0).Within(Tolerance),
                "Equal inter-breakpoint gaps have zero variance, so CV = 0.");
            Assert.That(r.IsClustered, Is.False,
                "CV = 0 < 1 (exponential null): not over-dispersed, so not clustered.");
        });
    }

    // M10 — Tight cluster plus a far outlier: over-dispersed gaps -> CV > 1 -> clustered.
    [Test]
    public void TestBreakpointClustering_TightClusterPlusOutlier_Clustered()
    {
        // gaps {1,1,1,997}: mean = 1000/4 = 250; population var = ((1-250)^2*3 + (997-250)^2)/4
        //   = (62001*3 + 558009)/4 = 744012/4 = 186003; sd = sqrt(186003); CV = sqrt(186003)/250.
        const double expectedMean = 250.0;
        double expectedCv = Math.Sqrt(186003.0) / 250.0; // = 1.7251226043386019
        long[] positions = { 0, 1, 2, 3, 1000 };
        var r = OncologyAnalyzer.TestBreakpointClustering(positions);
        Assert.Multiple(() =>
        {
            Assert.That(r.MeanGap, Is.EqualTo(expectedMean).Within(Tolerance), "Gaps {1,1,1,997} have mean 250.");
            Assert.That(r.CoefficientOfVariation, Is.EqualTo(expectedCv).Within(Tolerance),
                "sqrt(186003)/250 = 1.7251226043386019 (population sd / mean of inter-breakpoint gaps).");
            Assert.That(r.IsClustered, Is.True,
                "CV > 1 (exponential null CV = 1): over-dispersed, so clustered (Korbel & Campbell criterion A).");
        });
    }

    // S3 — Exact CV on a small known gap set {1,3}: mean 2, var ((1-2)^2+(3-2)^2)/2 = 1, sd 1, CV = 0.5.
    [Test]
    public void TestBreakpointClustering_KnownGaps_ExactCoefficientOfVariation()
    {
        long[] positions = { 0, 1, 4 }; // gaps {1,3}
        var r = OncologyAnalyzer.TestBreakpointClustering(positions);
        Assert.Multiple(() =>
        {
            Assert.That(r.MeanGap, Is.EqualTo(2.0).Within(Tolerance), "Gaps {1,3} have mean 2.");
            Assert.That(r.CoefficientOfVariation, Is.EqualTo(0.5).Within(Tolerance),
                "sd = 1, mean = 2 -> CV = 0.5 (population standard deviation).");
            Assert.That(r.IsClustered, Is.False, "CV = 0.5 < 1: not clustered.");
        });
    }

    // C3 — Null, empty, and fewer-than-three breakpoint sets: clustering undefined / rejected.
    [Test]
    public void TestBreakpointClustering_NullOrTooFew_HandledExplicitly()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentNullException>(
                () => OncologyAnalyzer.TestBreakpointClustering(null!),
                "Null breakpoint list is invalid input.");
            Assert.That(OncologyAnalyzer.TestBreakpointClustering(Array.Empty<long>()).IsClustered, Is.False,
                "No breakpoints: clustering cannot be assessed.");
            Assert.That(OncologyAnalyzer.TestBreakpointClustering(new long[] { 5, 9 }).IsClustered, Is.False,
                "Two breakpoints give a single gap; CV is undefined, so not clustered.");
        });
    }

    #endregion
}
