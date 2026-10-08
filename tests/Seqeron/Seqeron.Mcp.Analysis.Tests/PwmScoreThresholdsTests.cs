using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// pwm_score_thresholds — delegation to <c>PositionWeightMatrix.ScoreDistribution</c> and its
/// ThresholdFpr / ThresholdFnr / ThresholdBalanced / ThresholdPatser. Expected values: Biopython 1.88
/// Bio.motifs.thresholds on the Wikipedia PWM, locked in MotifFinder_PwmStrandsAndThresholds_Tests.
/// </summary>
[TestFixture]
public class PwmScoreThresholdsTests
{
    [Test]
    public void PwmScoreThresholds_Schema_ValidatesCorrectly()
    {
        var pwm = ScanWithPwmBothStrandsTests.WikipediaPwm();
        Assert.DoesNotThrow(() => AnalysisTools.PwmScoreThresholds(pwm));
        Assert.Throws<ArgumentException>(() => AnalysisTools.PwmScoreThresholds(null!));
        Assert.Throws<ArgumentException>(() => AnalysisTools.PwmScoreThresholds(pwm, new[] { 0.5, 0.5 }));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.PwmScoreThresholds(pwm, precision: 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.PwmScoreThresholds(pwm, fpr: 1.5));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.PwmScoreThresholds(pwm, fnr: -0.1));
        var infinite = new PwmInput(new[]
        {
            new[] { double.NegativeInfinity }, new[] { 1.0 }, new[] { 1.0 }, new[] { 1.0 }
        }, 1);
        Assert.Throws<ArgumentException>(() => AnalysisTools.PwmScoreThresholds(infinite));
    }

    [Test]
    public void PwmScoreThresholds_WikipediaDefaults_EqualBiopython()
    {
        var r = AnalysisTools.PwmScoreThresholds(ScanWithPwmBothStrandsTests.WikipediaPwm());
        Assert.Multiple(() =>
        {
            Assert.That(r.MinScore, Is.EqualTo(-14.88138790352414).Within(1e-12));
            Assert.That(r.Step, Is.EqualTo(0.002954652535685415).Within(1e-15));
            Assert.That(r.PointCount, Is.EqualTo(9000));
            Assert.That(r.MeanScore, Is.EqualTo(5.522241422369563).Within(1e-12));
            Assert.That(r.ThresholdFpr, Is.EqualTo(4.028388324862519).Within(1e-9));
            Assert.That(r.ThresholdFnr, Is.EqualTo(1.2303323735684302).Within(1e-9));
            Assert.That(r.ThresholdBalanced, Is.EqualTo(0.1430202404361971).Within(1e-9));
            Assert.That(r.BalancedFalsePositiveRate, Is.EqualTo(0.06385040283203125).Within(1e-12));
            Assert.That(r.ThresholdPatser, Is.EqualTo(2.5924271925194056).Within(1e-9));
        });
    }

    [Test]
    public void PwmScoreThresholds_BackgroundAndPrecision_EqualBiopython()
    {
        var r = AnalysisTools.PwmScoreThresholds(ScanWithPwmBothStrandsTests.WikipediaPwm(),
            new[] { 0.3, 0.2, 0.2, 0.3 }, precision: 100, fpr: 0.01, fnr: 0.05, rateProportion: 2.0);
        Assert.Multiple(() =>
        {
            Assert.That(r.ThresholdFpr, Is.EqualTo(4.254351868562161).Within(1e-9));
            Assert.That(r.ThresholdFnr, Is.EqualTo(-0.4482487864018605).Within(1e-9));
            Assert.That(r.ThresholdBalanced, Is.EqualTo(3.3079165166197164).Within(1e-9));
            Assert.That(r.ThresholdPatser, Is.EqualTo(3.4262209356125233).Within(1e-9));
            Assert.That(r.MeanScore, Is.EqualTo(5.975137019521005).Within(1e-12));
        });
    }
}
