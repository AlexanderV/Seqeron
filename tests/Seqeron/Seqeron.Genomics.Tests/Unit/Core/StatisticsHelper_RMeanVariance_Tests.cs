// StatisticsHelper.RMean / RSampleVariance — R mean.default (refined by the mean residual) and var(x), the shared
// canonical behind OneSampleTTestPValue / WelchTTestPValue and Battenberg's t-test / bootstrap sd (ONCO-ASCAT-001,
// FIN-B24 duplication sweep). Reference values: R 4.3.3, printed with %.17g.

using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Core;

[TestFixture]
public class StatisticsHelper_RMeanVariance_Tests
{
    [Test]
    public void RMean_MatchesR_IncludingTheResidualRefinement()
    {
        // Plain left-to-right sum/n gives 0.52750000000000008; R mean (and the refinement) gives 0.52749999999999997.
        double[] x = { 0.59, 0.13, 0.92, 0.47 };
        Assert.Multiple(() =>
        {
            Assert.That(StatisticsHelper.RMean(x), Is.EqualTo(0.52749999999999997));
            Assert.That(StatisticsHelper.RMean(new double[] { 1, 2, 3, 4 }), Is.EqualTo(2.5));
            Assert.That(StatisticsHelper.RMean(new[] { 0.61, 0.64, 0.66, 0.59, 0.70 }), Is.EqualTo(0.64000000000000001));
        });
    }

    [Test]
    public void RSampleVariance_MatchesR()
    {
        Assert.Multiple(() =>
        {
            Assert.That(StatisticsHelper.RSampleVariance(new[] { 0.59, 0.13, 0.92, 0.47 }), Is.EqualTo(0.10642500000000001));
            Assert.That(StatisticsHelper.RSampleVariance(new double[] { 1, 2, 3, 4 }), Is.EqualTo(1.6666666666666667));
            Assert.That(StatisticsHelper.RSampleVariance(new[] { 0.61, 0.64, 0.66, 0.59, 0.70 }), Is.EqualTo(0.0018499999999999999));
        });
    }

    [Test]
    public void RMean_RSampleVariance_EdgeCases_AsR()
    {
        Assert.Multiple(() =>
        {
            Assert.That(StatisticsHelper.RMean(Array.Empty<double>()), Is.NaN);              // mean(numeric(0)) = NaN
            Assert.That(StatisticsHelper.RMean(new[] { 1.0, double.NaN }), Is.NaN);          // NA propagates
            Assert.That(StatisticsHelper.RMean(new[] { 3.5 }), Is.EqualTo(3.5));
            Assert.That(StatisticsHelper.RSampleVariance(new[] { 3.5 }), Is.NaN);              // var(3.5) = NA
            Assert.That(StatisticsHelper.RSampleVariance(new[] { 1.0, double.NaN, 2.0 }), Is.NaN);
            Assert.Throws<ArgumentNullException>(() => StatisticsHelper.RMean(null!));
            Assert.Throws<ArgumentNullException>(() => StatisticsHelper.RSampleVariance(null!));
        });
    }
}
