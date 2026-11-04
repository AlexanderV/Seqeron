// StatisticsHelper.SampleQuantileType7 — R quantile(x, p, type = 7, names = FALSE), added for Battenberg's filterMarkS4
// (ONCO-ASCAT-001, B24 F40). Reference values: R 4.3.3, printed with %.17g.

using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Core;

[TestFixture]
public class StatisticsHelper_QuantileType7_Tests
{
    private static readonly double[] Sample = { 0.31, 0.12, 0.97, 0.55, 0.12, 0.44, 0.08 };

    [TestCase(0.0, 0.080000000000000002)]
    [TestCase(0.05, 0.091999999999999998)]
    [TestCase(0.25, 0.12)]
    [TestCase(0.5, 0.31)]
    [TestCase(0.85, 0.59199999999999986)]
    [TestCase(0.88, 0.66760000000000019)]
    [TestCase(1.0, 0.96999999999999997)]
    public void SampleQuantileType7_MatchesR(double p, double expected)
    {
        Assert.That(StatisticsHelper.SampleQuantileType7(Sample, p), Is.EqualTo(expected).Within(1e-15));
    }

    [Test]
    public void SampleQuantileType7_DocExample_And_Guards()
    {
        Assert.Multiple(() =>
        {
            Assert.That(StatisticsHelper.SampleQuantileType7(new double[] { 3, 1, 4, 1, 5 }, 0.85), Is.EqualTo(4.4000000000000004));
            Assert.That(StatisticsHelper.SampleQuantileType7(new[] { 2.5 }, 0.3), Is.EqualTo(2.5));
            Assert.Throws<ArgumentNullException>(() => StatisticsHelper.SampleQuantileType7(null!, 0.5));
            Assert.Throws<ArgumentException>(() => StatisticsHelper.SampleQuantileType7(Array.Empty<double>(), 0.5));
            Assert.Throws<ArgumentException>(() => StatisticsHelper.SampleQuantileType7(new[] { 1.0, double.NaN }, 0.5));
            Assert.Throws<ArgumentOutOfRangeException>(() => StatisticsHelper.SampleQuantileType7(new[] { 1.0 }, 1.5));
        });
    }
}
