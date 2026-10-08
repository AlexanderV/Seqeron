// StatisticsHelper.LogPoissonRangeProbability / LogNegativeBinomialRangeProbability — shared helpers added for
// MOTIF-DISCOVER-001 (RSAT oligo-analysis -calibN / -calib1: RSAT::stats::sum_of_poisson / sum_of_negbin2 over the
// right tail [occ, n]).
// References: mpmath (50 digits) Σ_{x=a..b} of the exact pmf (Poisson: e^-λ λ^x / x!; negative binomial in RSAT's
// parameterisation p = v/m − 1, k = m/p, q = 1 + p: Γ(k+x)/(Γ(k) x!) p^x q^-(k+x)); scipy.stats.poisson / nbinom(k, 1/q)
// cdf differences agree (e.g. poisson.cdf(61, 1.5) − poisson.cdf(5, 1.5) = 0.004455980775247892).

namespace Seqeron.Genomics.Tests.Unit.Core;

[TestFixture]
public class StatisticsHelper_PoissonNegBinRange_Tests
{
    private static readonly object[] PoissonCases =
    {
        new object[] { 6L, 61L, 1.5, 0.0044559807752478491891, -5.4135080907183799402 },
        new object[] { 4L, 61L, 0.8, 0.0090798578001539854022, -4.701696747267559156 },
        new object[] { 2L, 100L, 10.0, 0.99950060077261266663, -0.00049952396871370105852 },
        new object[] { 40L, 1000L, 10.0, 7.3416363145604714222e-13, -27.940044460002594768 },
        new object[] { 150L, 10000L, 3.0, 3.2894660604843631623e-193, -443.20819768789034093 },
        new object[] { 0L, 5L, 10.0, 0.067085962879031782286, -2.7017804539230561654 },
        new object[] { 3L, 7L, 10.0, 0.2174512508861873643, -1.5257805874317547531 },
        new object[] { 12L, 12L, 10.0, 0.094780330091767650746, -2.3561933797333379413 },
    };

    private static readonly object[] NegBinCases =
    {
        new object[] { 6L, 61L, 0.9, 1.44, 0.006714235225703917789, -5.003525345931676168 },
        new object[] { 5L, 61L, 2.0, 2.25, 0.063244578406644241171, -2.7607458719394102222 },
        new object[] { 1L, 50L, 3.0, 30.0, 0.53531171692864821826, -0.62490605329508067143 },
        new object[] { 60L, 1000L, 3.0, 4.0, 5.8202515584905271664e-28, -62.711039119915686328 },
        new object[] { 2L, 4L, 5.0, 6.0, 0.40219682609332264036, -0.91081369303891271025 },
    };

    [TestCaseSource(nameof(PoissonCases))]
    public void LogPoissonRange_EqualsMpmath(long from, long to, double lambda, double expected, double expectedLog)
    {
        double log = StatisticsHelper.LogPoissonRangeProbability(from, to, lambda);
        Assert.Multiple(() =>
        {
            Assert.That(log, Is.EqualTo(expectedLog).Within(1e-12 * Math.Max(1, Math.Abs(expectedLog))));
            Assert.That(Math.Exp(log), Is.EqualTo(expected).Within(1e-10).Percent);
        });
    }

    [TestCaseSource(nameof(NegBinCases))]
    public void LogNegativeBinomialRange_EqualsMpmath(long from, long to, double mean, double variance, double expected, double expectedLog)
    {
        double log = StatisticsHelper.LogNegativeBinomialRangeProbability(from, to, mean, variance);
        Assert.Multiple(() =>
        {
            Assert.That(log, Is.EqualTo(expectedLog).Within(1e-12 * Math.Max(1, Math.Abs(expectedLog))));
            Assert.That(Math.Exp(log), Is.EqualTo(expected).Within(1e-10).Percent);
        });
    }

    [Test]
    public void RangeFromZero_IsAtMostOne_AndWholeSupportIsOne()
    {
        // mpmath: P(0 <= X <= 61) for Poisson(1.5) = 1 − 3.0e-66; NB(m 2, v 2.25): ln = −4.2685856536351676701e-45.
        Assert.Multiple(() =>
        {
            Assert.That(StatisticsHelper.LogPoissonRangeProbability(0, 61, 1.5), Is.LessThanOrEqualTo(0.0).And.GreaterThan(-1e-15));
            Assert.That(StatisticsHelper.LogNegativeBinomialRangeProbability(0, 61, 2.0, 2.25), Is.LessThanOrEqualTo(0.0).And.GreaterThan(-1e-15));
            Assert.That(StatisticsHelper.LogPoissonRangeProbability(0, long.MaxValue, 7.0), Is.EqualTo(0.0));
            Assert.That(StatisticsHelper.LogPoissonRangeProbability(5, 4, 7.0), Is.EqualTo(double.NegativeInfinity));
        });
    }

    [Test]
    public void InvalidParameters_Throw()
    {
        Assert.Multiple(() =>
        {
            Assert.That(() => StatisticsHelper.LogPoissonRangeProbability(1, 10, 0.0), NUnit.Framework.Throws.TypeOf<ArgumentOutOfRangeException>());
            Assert.That(() => StatisticsHelper.LogPoissonRangeProbability(1, 10, double.NaN), NUnit.Framework.Throws.TypeOf<ArgumentOutOfRangeException>());
            Assert.That(() => StatisticsHelper.LogNegativeBinomialRangeProbability(1, 10, 2.0, 2.0), NUnit.Framework.Throws.TypeOf<ArgumentOutOfRangeException>());
            Assert.That(() => StatisticsHelper.LogNegativeBinomialRangeProbability(1, 10, 0.0, 2.0), NUnit.Framework.Throws.TypeOf<ArgumentOutOfRangeException>());
        });
    }
}
