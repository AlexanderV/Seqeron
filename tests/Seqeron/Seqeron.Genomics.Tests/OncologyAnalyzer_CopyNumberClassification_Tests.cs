// ONCO-CNA-001 — Copy-Number Alteration Classification
// Evidence: docs/Evidence/ONCO-CNA-001-Evidence.md
// TestSpec: tests/TestSpecs/ONCO-CNA-001.md
// Source: CNVkit cnvlib/call.py — absolute_threshold, _log2_ratio_to_absolute_pure, do_call
//         (default thresholds (-1.1, -0.25, 0.2, 0.7); n = ploidy*2^log2).
//         Mermel CH et al. (2011). GISTIC2.0. Genome Biology 12(4):R41.
//         https://doi.org/10.1186/gb-2011-12-4-r41

using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using Seqeron.Genomics.Oncology;
using CopyNumberState = Seqeron.Genomics.Oncology.OncologyAnalyzer.CopyNumberState;

namespace Seqeron.Genomics.Tests;

[TestFixture]
public class OncologyAnalyzer_CopyNumberClassification_Tests
{
    private const double Tolerance = 1e-10;

    #region Log2RatioToCopyNumber

    // M11 — n = ploidy*2^log2 (CNVkit _log2_ratio_to_absolute_pure). Diploid: 0->2, 1->4, -1->1.
    [Test]
    public void Log2RatioToCopyNumber_DiploidReference_MatchesFormula()
    {
        Assert.Multiple(() =>
        {
            Assert.That(OncologyAnalyzer.Log2RatioToCopyNumber(0.0), Is.EqualTo(2.0).Within(Tolerance),
                "log2 0 with diploid ploidy 2 gives n = 2*2^0 = 2.0 (neutral diploid).");
            Assert.That(OncologyAnalyzer.Log2RatioToCopyNumber(1.0), Is.EqualTo(4.0).Within(Tolerance),
                "log2 1 gives n = 2*2^1 = 4.0.");
            Assert.That(OncologyAnalyzer.Log2RatioToCopyNumber(-1.0), Is.EqualTo(1.0).Within(Tolerance),
                "log2 -1 gives n = 2*2^-1 = 1.0 (single-copy loss).");
            // log2(3/2) = 0.5849625... -> n = 2*(3/2) = 3.0 exactly.
            Assert.That(OncologyAnalyzer.Log2RatioToCopyNumber(Math.Log2(1.5)), Is.EqualTo(3.0).Within(Tolerance),
                "log2(3/2) gives n = 2*2^log2(3/2) = 3.0 (single-copy gain).");
        });
    }

    // E4 — n = ploidy*2^log2 requires ploidy > 0.
    [Test]
    public void Log2RatioToCopyNumber_NonPositivePloidy_Throws()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.Log2RatioToCopyNumber(0.0, 0.0),
                "Ploidy 0 is invalid (n = ploidy*2^log2).");
            Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.Log2RatioToCopyNumber(0.0, -1.0),
                "Negative ploidy is invalid.");
        });
    }

    #endregion

    #region CallCopyNumber

    // M1,M2,M3,M4,M5,M6 — integer CN per CNVkit absolute_threshold default cutoffs.
    [TestCase(-2.0, 0, TestName = "CallCopyNumber_DeepDeletionRange_Returns0")]
    [TestCase(-1.0, 1, TestName = "CallCopyNumber_LossRange_Returns1")]
    [TestCase(0.0, 2, TestName = "CallCopyNumber_NeutralRange_Returns2")]
    [TestCase(0.5849625007211562, 3, TestName = "CallCopyNumber_GainRange_Returns3")] // log2(3/2)
    [TestCase(1.0, 4, TestName = "CallCopyNumber_AmplificationRange_Returns4")]       // ceil(2*2^1)=4
    [TestCase(2.0, 8, TestName = "CallCopyNumber_HighAmplification_Returns8")]        // ceil(2*2^2)=8
    public void CallCopyNumber_DefaultThresholds_ReturnsExpectedInteger(double log2, int expected)
    {
        int cn = OncologyAnalyzer.CallCopyNumber(log2);
        Assert.That(cn, Is.EqualTo(expected),
            $"log2 {log2} maps to integer CN {expected} under CNVkit default cutoffs (-1.1,-0.25,0.2,0.7).");
    }

    // M7,M8,M9,M10 — boundary inclusivity: CNVkit uses 'log2 <= thresh', lower state wins on the edge.
    [TestCase(-1.1, 0, TestName = "CallCopyNumber_BoundaryDeepDeletion_Returns0")]
    [TestCase(-0.25, 1, TestName = "CallCopyNumber_BoundaryLoss_Returns1")]
    [TestCase(0.2, 2, TestName = "CallCopyNumber_BoundaryNeutral_Returns2")]
    [TestCase(0.7, 3, TestName = "CallCopyNumber_BoundaryGain_Returns3")]
    public void CallCopyNumber_OnThresholdBoundary_AssignsLowerState(double log2, int expected)
    {
        int cn = OncologyAnalyzer.CallCopyNumber(log2);
        Assert.That(cn, Is.EqualTo(expected),
            $"log2 exactly {log2} is on a cutoff; CNVkit 'log2 <= thresh' assigns CN {expected}.");
    }

    // M12 — NaN log2 is a no-call -> neutral reference copy number (rounded diploid = 2).
    [Test]
    public void CallCopyNumber_NaNLog2_ReturnsNeutralCopyNumber()
    {
        int cn = OncologyAnalyzer.CallCopyNumber(double.NaN);
        Assert.That(cn, Is.EqualTo(2),
            "CNVkit replaces a NaN log2 with the neutral reference copy number (diploid = 2).");
    }

    // S3 — above the last cutoff, CNVkit uses ceil(2*2^log2), NOT round. log2 0.8 -> 2*2^0.8=3.482 -> 4.
    [Test]
    public void CallCopyNumber_AboveLastThreshold_UsesCeilingNotRound()
    {
        int cn = OncologyAnalyzer.CallCopyNumber(0.8);
        // round(3.482) would be 3 (= Gain); ceil = 4 (= Amplification). A wrong (round) impl fails here.
        Assert.That(cn, Is.EqualTo(4),
            "log2 0.8 -> 2*2^0.8 = 3.482; CNVkit rounds UP (ceil) to CN 4, not 3.");
    }

    // S1 — custom (germline-tuned) thresholds override the defaults.
    [Test]
    public void CallCopyNumber_CustomThresholds_UsesProvidedCutoffs()
    {
        var custom = new[] { -0.4, -0.1, 0.1, 0.4 };
        // log2 -0.2: under custom cutoffs it is in (-0.4,-0.1] -> CN 1 (Loss); under defaults it is in
        // (-0.25,0.2] -> CN 2 (Neutral). The custom thresholds therefore change the call.
        int cnCustom = OncologyAnalyzer.CallCopyNumber(-0.2, custom);
        int cnDefault = OncologyAnalyzer.CallCopyNumber(-0.2);
        Assert.Multiple(() =>
        {
            Assert.That(cnCustom, Is.EqualTo(1),
                "Under custom cutoffs (-0.4,-0.1,0.1,0.4), log2 -0.2 is in (-0.4,-0.1] -> CN 1 (Loss).");
            Assert.That(cnDefault, Is.EqualTo(2),
                "Under default cutoffs, log2 -0.2 is in (-0.25,0.2] -> CN 2 (Neutral); the custom set changes the call.");
        });
    }

    #endregion

    #region ClassifyCopyNumber

    // M1–M6 — full classification: state + integer CN + absolute CN.
    [TestCase(-2.0, CopyNumberState.DeepDeletion, 0, TestName = "ClassifyCopyNumber_DeepDeletion")]
    [TestCase(-1.0, CopyNumberState.Loss, 1, TestName = "ClassifyCopyNumber_Loss")]
    [TestCase(0.0, CopyNumberState.Neutral, 2, TestName = "ClassifyCopyNumber_Neutral")]
    [TestCase(0.5849625007211562, CopyNumberState.Gain, 3, TestName = "ClassifyCopyNumber_Gain")]
    [TestCase(1.0, CopyNumberState.Amplification, 4, TestName = "ClassifyCopyNumber_Amplification")]
    public void ClassifyCopyNumber_DefaultThresholds_ReturnsExpectedState(
        double log2, CopyNumberState expectedState, int expectedCn)
    {
        var call = OncologyAnalyzer.ClassifyCopyNumber(log2);
        Assert.Multiple(() =>
        {
            Assert.That(call.State, Is.EqualTo(expectedState),
                $"log2 {log2} classifies as {expectedState} per CNVkit DEL/LOSS/neutral/GAIN/AMP mapping.");
            Assert.That(call.IntegerCopyNumber, Is.EqualTo(expectedCn),
                $"log2 {log2} integer CN is {expectedCn}.");
            Assert.That(call.Log2Ratio, Is.EqualTo(log2).Within(Tolerance),
                "The input log2 ratio is echoed in the call.");
        });
    }

    // M11 — absolute copy number carried on the call equals n = 2*2^log2.
    [Test]
    public void ClassifyCopyNumber_CarriesContinuousAbsoluteCopyNumber()
    {
        var amp = OncologyAnalyzer.ClassifyCopyNumber(1.0);
        var neutral = OncologyAnalyzer.ClassifyCopyNumber(0.0);
        Assert.Multiple(() =>
        {
            Assert.That(amp.AbsoluteCopyNumber, Is.EqualTo(4.0).Within(Tolerance),
                "log2 1 -> absolute CN 2*2^1 = 4.0.");
            Assert.That(neutral.AbsoluteCopyNumber, Is.EqualTo(2.0).Within(Tolerance),
                "log2 0 -> absolute CN 2*2^0 = 2.0.");
        });
    }

    // M12 — NaN classifies as Neutral with absolute CN = ploidy (no-call).
    [Test]
    public void ClassifyCopyNumber_NaN_IsNeutralNoCall()
    {
        var call = OncologyAnalyzer.ClassifyCopyNumber(double.NaN);
        Assert.Multiple(() =>
        {
            Assert.That(call.State, Is.EqualTo(CopyNumberState.Neutral),
                "NaN log2 is a no-call -> Neutral (CNVkit replaces with neutral reference CN).");
            Assert.That(call.IntegerCopyNumber, Is.EqualTo(2), "No-call integer CN is the diploid baseline 2.");
            Assert.That(call.AbsoluteCopyNumber, Is.EqualTo(2.0).Within(Tolerance),
                "No-call absolute CN is the reference ploidy 2.0.");
        });
    }

    #endregion

    #region ClassifyCopyNumbers (batch)

    // M13 — batch preserves order and length; per-element classification.
    [Test]
    public void ClassifyCopyNumbers_Batch_PreservesOrderAndLength()
    {
        var input = new[] { -2.0, 0.0, 1.0 };
        var calls = OncologyAnalyzer.ClassifyCopyNumbers(input);
        Assert.Multiple(() =>
        {
            Assert.That(calls, Has.Count.EqualTo(3), "One call per input log2 ratio (INV-05).");
            Assert.That(calls[0].State, Is.EqualTo(CopyNumberState.DeepDeletion),
                "Index 0 (log2 -2.0) -> DeepDeletion.");
            Assert.That(calls[1].State, Is.EqualTo(CopyNumberState.Neutral),
                "Index 1 (log2 0.0) -> Neutral.");
            Assert.That(calls[2].State, Is.EqualTo(CopyNumberState.Amplification),
                "Index 2 (log2 1.0) -> Amplification; order preserved.");
        });
    }

    // C1 — empty batch -> empty result.
    [Test]
    public void ClassifyCopyNumbers_EmptyInput_ReturnsEmpty()
    {
        var calls = OncologyAnalyzer.ClassifyCopyNumbers(Array.Empty<double>());
        Assert.That(calls, Is.Empty, "An empty batch produces no calls.");
    }

    // E5 — null batch throws.
    [Test]
    public void ClassifyCopyNumbers_NullInput_Throws()
    {
        Assert.Throws<ArgumentNullException>(() => OncologyAnalyzer.ClassifyCopyNumbers(null!), "A null log2-ratio sequence is rejected.");
    }

    // S2 / INV-02 — integer CN is monotonically non-decreasing in log2 ratio.
    [Test]
    public void ClassifyCopyNumbers_AscendingLog2_ProducesNonDecreasingCopyNumber()
    {
        var ascending = new[] { -3.0, -1.5, -1.0, -0.3, 0.0, 0.1, 0.5, 0.7, 1.0, 2.5 };
        var calls = OncologyAnalyzer.ClassifyCopyNumbers(ascending);
        for (int i = 1; i < calls.Count; i++)
        {
            Assert.That(calls[i].IntegerCopyNumber, Is.GreaterThanOrEqualTo(calls[i - 1].IntegerCopyNumber),
                $"CN must not decrease as log2 increases (INV-02): index {i - 1}->{i}.");
        }
    }

    #endregion

    #region Validation

    // E1 — null thresholds is treated as default (not an error) for single call.
    [Test]
    public void CallCopyNumber_NullThresholds_UsesDefaults()
    {
        int cn = OncologyAnalyzer.CallCopyNumber(0.0, thresholds: null);
        Assert.That(cn, Is.EqualTo(2), "Null thresholds fall back to the documented defaults -> Neutral CN 2.");
    }

    // E2 — wrong threshold count throws.
    [Test]
    public void CallCopyNumber_WrongThresholdCount_Throws()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentException>(() => OncologyAnalyzer.CallCopyNumber(0.0, new[] { -1.0, 0.0, 1.0 }),
                "Three thresholds cannot define five states.");
            Assert.Throws<ArgumentException>(() => OncologyAnalyzer.CallCopyNumber(0.0, new[] { -1.0, 0.0, 1.0, 2.0, 3.0 }),
                "Five thresholds is also invalid; exactly four required.");
        });
    }

    // E3 — non-ascending thresholds throw.
    [Test]
    public void CallCopyNumber_NonAscendingThresholds_Throws()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentException>(() => OncologyAnalyzer.CallCopyNumber(0.0, new[] { -1.0, 0.5, 0.2, 0.7 }),
                "Thresholds out of ascending order are rejected.");
            Assert.Throws<ArgumentException>(() => OncologyAnalyzer.CallCopyNumber(0.0, new[] { -1.0, -1.0, 0.2, 0.7 }),
                "Equal adjacent thresholds (not strictly ascending) are rejected.");
            Assert.Throws<ArgumentException>(() => OncologyAnalyzer.CallCopyNumber(0.0, new[] { -1.0, double.NaN, 0.2, 0.7 }),
                "NaN thresholds are rejected.");
        });
    }

    // E4 — non-positive ploidy throws on classify path.
    [Test]
    public void ClassifyCopyNumber_NonPositivePloidy_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.ClassifyCopyNumber(0.0, ploidy: 0.0),
            "Ploidy must be positive.");
    }

    #endregion
}
