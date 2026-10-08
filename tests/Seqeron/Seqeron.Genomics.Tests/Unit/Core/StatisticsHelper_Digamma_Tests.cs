// StatisticsHelper.Digamma — shared ψ(x) added for the Grassberger (2003) k-mer entropy estimator (SEQ-COMPLEX-KMER-001,
// B04 F54). Reference values: mpmath.digamma at 40 digits (= scipy.special.digamma to the last printed digit).

using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Core;

[TestFixture]
public class StatisticsHelper_Digamma_Tests
{
    [TestCase(1.0, -0.5772156649015329)]
    [TestCase(0.5, -1.9635100260214235)]
    [TestCase(2.0, 0.42278433509846713)]
    [TestCase(3.7, 1.1671535393615113)]
    [TestCase(10.0, 2.251752589066721)]
    [TestCase(10.5, 2.3030010342976865)]
    [TestCase(123.456, 4.811829323828985)]
    [TestCase(1e6, 13.815510057964191)]
    [TestCase(1e-3, -1000.5755719318103)]
    public void Digamma_MatchesMpmath(double x, double expected)
    {
        Assert.That(StatisticsHelper.Digamma(x), Is.EqualTo(expected).Within(2e-15 * Math.Max(1, Math.Abs(expected))));
    }

    [Test]
    public void Digamma_SatisfiesRecurrence([Values(0.25, 1.5, 7.3, 9.99, 42.0)] double x)
    {
        Assert.That(StatisticsHelper.Digamma(x + 1), Is.EqualTo(StatisticsHelper.Digamma(x) + 1 / x).Within(1e-14));
    }

    [TestCase(0.0)]
    [TestCase(-1.0)]
    [TestCase(double.NaN)]
    [TestCase(double.PositiveInfinity)]
    public void Digamma_OutsideDomain_Throws(double x)
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => StatisticsHelper.Digamma(x));
    }
}
