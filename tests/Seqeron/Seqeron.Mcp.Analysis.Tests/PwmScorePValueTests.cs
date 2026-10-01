using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// pwm_score_pvalue — delegation to <c>MotifFinder.PwmScorePValue</c> / <c>PwmScoreThresholdForPValue</c>.
/// Expected values: exhaustive enumeration of all 4^9 words of the Wikipedia PWM (and TFM-Pvalue for the
/// non-tied thresholds), locked in MotifFinder_PwmPValue_Tests.
/// </summary>
[TestFixture]
public class PwmScorePValueTests
{
    [Test]
    public void PwmScorePValue_Schema_ValidatesCorrectly()
    {
        var pwm = ScanWithPwmBothStrandsTests.WikipediaPwm();
        Assert.DoesNotThrow(() => AnalysisTools.PwmScorePValue(pwm, score: 0.0));
        Assert.DoesNotThrow(() => AnalysisTools.PwmScorePValue(pwm, pValue: 1e-3));
        Assert.Throws<ArgumentException>(() => AnalysisTools.PwmScorePValue(null!, score: 0.0));
        Assert.Throws<ArgumentException>(() => AnalysisTools.PwmScorePValue(pwm));
        Assert.Throws<ArgumentException>(() => AnalysisTools.PwmScorePValue(pwm, 0.0, 0.01));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.PwmScorePValue(pwm, score: double.NaN));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.PwmScorePValue(pwm, pValue: 1.5));
        Assert.Throws<ArgumentException>(() => AnalysisTools.PwmScorePValue(pwm, score: 0.0, background: new[] { 0.5, 0.5 }));
        var infinite = new PwmInput(new[]
        {
            new[] { double.NegativeInfinity }, new[] { 1.0 }, new[] { 1.0 }, new[] { 1.0 }
        }, 1);
        Assert.Throws<ArgumentException>(() => AnalysisTools.PwmScorePValue(infinite, score: 0.0));
    }

    [Test]
    public void PwmScorePValue_Binding_EqualsExhaustiveEnumeration()
    {
        var pwm = ScanWithPwmBothStrandsTests.WikipediaPwm();
        var site = AnalysisTools.PwmScorePValue(pwm, score: 4.77915994208994);
        Assert.Multiple(() =>
        {
            Assert.That(site.PValue, Is.EqualTo(0.006160736083984375).Within(1e-15));
            Assert.That(site.IsExact, Is.True);
            Assert.That(site.ScoreAboveMaximum, Is.False);
        });

        var t = AnalysisTools.PwmScorePValue(pwm, pValue: 1e-3);
        Assert.That(t.Score, Is.EqualTo(7.150252343998389).Within(1e-12));
        Assert.That(t.PValue, Is.EqualTo(0.000980377197265625).Within(1e-18));

        var bg = AnalysisTools.PwmScorePValue(pwm, pValue: 1e-4, background: new[] { 0.3, 0.2, 0.2, 0.3 });
        Assert.That(bg.Score, Is.EqualTo(9.567600003964435).Within(1e-12));
        Assert.That(bg.PValue, Is.EqualTo(9.5256e-05).Within(1e-15));

        var none = AnalysisTools.PwmScorePValue(pwm, pValue: 1e-6);
        Assert.That(none.ScoreAboveMaximum, Is.True);
        Assert.That(none.Score, Is.Null);
        Assert.That(none.PValue, Is.Zero);
    }
}
