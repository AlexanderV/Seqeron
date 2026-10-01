// StatisticsHelper.BinomialUpperTail / LogBinomialUpperTail / Log1P / ExpM1 — shared helpers added for
// MOTIF-DISCOVER-001 / MOTIF-SHARED-001 (RSAT oligo-analysis occ_P / ms_P, RSAT::stats::sum_of_binomials).
// References: mpmath exact sums (cross-checked with scipy.stats.binom.sf(k - 1, n, p) / binom.logsf); for p below
// the double range, mpmath (50 digits)  ln Σ_{x=k..n} C(n,x) p^x (1-p)^(n-x).

namespace Seqeron.Genomics.Tests.Unit.Core;

[TestFixture]
public class StatisticsHelper_BinomialUpperTail_Tests
{
    // (successes, trials, p, P(X >= k), ln P(X >= k)): exact values from mpmath (60 digits, p taken as the exact double);
    // scipy.stats.binom.sf / logsf agree to <= 5e-14 relative (e.g. sf(5, 60, 1/256) = 1.484494210307282e-07).
    private static readonly object[] ReferenceCases =
    {
        new object[] { 6L, 60L, 1.0 / 256, 1.4844942103072811137e-07, -15.723021535840387106 },
        new object[] { 2L, 60L, 1.0 / 256, 0.02325031310108469454, -3.7614366803261269081 },
        new object[] { 1100L, 1_000_000L, 0.001, 0.00095746697002854186641, -6.9512193335582177704 },
        new object[] { 250L, 1000L, 0.3, 0.9998014526737669589, -0.00019856703936278755562 },
        new object[] { 200L, 1000L, 0.01, 2.2214826738393917669e-188, -432.08782263882422573 },
        new object[] { 10L, 10L, 0.5, 0.0009765625, -6.9314718055994530942 },
        new object[] { 1L, 1L, 0.9, 0.9000000000000000222, -0.10536051565782627656 },
    };

    [TestCaseSource(nameof(ReferenceCases))]
    public void BinomialUpperTail_EqualsExactSurvivalFunction(long k, long n, double p, double sf, double logSf)
    {
        Assert.Multiple(() =>
        {
            Assert.That(StatisticsHelper.BinomialUpperTail(k, n, p), Is.EqualTo(sf).Within(1e-10).Percent,
                "P(X >= k) must equal scipy.stats.binom.sf(k - 1, n, p).");
            Assert.That(StatisticsHelper.LogBinomialUpperTail(k, n, Math.Log(p)), Is.EqualTo(logSf).Within(1e-10).Percent,
                "ln P(X >= k) must equal scipy.stats.binom.logsf(k - 1, n, p).");
        });
    }

    [Test]
    public void LogBinomialUpperTail_ProbabilityBelowDoubleRange_EqualsExactLogSum()
    {
        // p = e^-800 (not representable); mpmath: ln Σ_{x>=2} C(1000,x) p^x (1-p)^(1000-x) = -1586.8786371229292547.
        double logTail = StatisticsHelper.LogBinomialUpperTail(2, 1000, -800.0);

        Assert.That(logTail, Is.EqualTo(-1586.8786371229292547).Within(1e-10).Percent);
    }

    [Test]
    public void BinomialUpperTail_TrivialBounds()
    {
        Assert.Multiple(() =>
        {
            Assert.That(StatisticsHelper.BinomialUpperTail(0, 10, 0.3), Is.EqualTo(1.0), "P(X >= 0) = 1.");
            Assert.That(StatisticsHelper.BinomialUpperTail(11, 10, 0.3), Is.EqualTo(0.0), "P(X >= n+1) = 0.");
            Assert.That(StatisticsHelper.BinomialUpperTail(1, 10, 0.0), Is.EqualTo(0.0), "p = 0: X = 0 surely.");
            Assert.That(StatisticsHelper.BinomialUpperTail(10, 10, 1.0), Is.EqualTo(1.0), "p = 1: X = n surely.");
            Assert.That(StatisticsHelper.BinomialUpperTail(5, 0, 0.5), Is.EqualTo(0.0), "no trials, no successes.");
        });
    }

    [Test]
    public void BinomialUpperTail_InvalidArguments_Throw()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => StatisticsHelper.BinomialUpperTail(1, 10, -0.1));
            Assert.Throws<ArgumentOutOfRangeException>(() => StatisticsHelper.BinomialUpperTail(1, 10, 1.5));
            Assert.Throws<ArgumentOutOfRangeException>(() => StatisticsHelper.BinomialUpperTail(1, 10, double.NaN));
            Assert.Throws<ArgumentOutOfRangeException>(() => StatisticsHelper.BinomialUpperTail(1, -1, 0.5));
            Assert.Throws<ArgumentOutOfRangeException>(() => StatisticsHelper.LogBinomialUpperTail(1, 10, 0.1));
        });
    }

    [Test]
    public void BinomialUpperTail_MonotoneNonIncreasingInSuccesses_AndSumsWithPmf()
    {
        // P(X >= k) - P(X >= k+1) = P(X = k) = C(n,k) p^k q^(n-k) (exact definition), and the tail never increases in k.
        const long n = 40;
        const double p = 0.17;
        double previous = 1.0;
        for (long k = 0; k <= n + 1; k++)
        {
            double tail = StatisticsHelper.BinomialUpperTail(k, n, p);
            Assert.That(tail, Is.LessThanOrEqualTo(previous * (1 + 1e-14)), $"tail must not increase at k = {k}");
            if (k <= n)
            {
                double pmf = Math.Exp(LogChoose(n, k) + k * Math.Log(p) + (n - k) * Math.Log(1 - p));
                double diff = tail - StatisticsHelper.BinomialUpperTail(k + 1, n, p);
                Assert.That(diff, Is.EqualTo(pmf).Within(1e-12), $"tail difference at k = {k} must be the pmf");
            }
            previous = tail;
        }
    }

    [Test]
    public void Log1P_ExpM1_AccurateForTinyArguments()
    {
        // numpy.log1p(1e-20) = 1e-20, numpy.expm1(1e-20) = 1e-20 (double.LogP1 / double.ExpM1 return 0 here);
        // numpy.log1p(-0.5) = -0.6931471805599453, numpy.expm1(1.0) = 1.718281828459045.
        Assert.Multiple(() =>
        {
            Assert.That(StatisticsHelper.Log1P(1e-20), Is.EqualTo(1e-20));
            Assert.That(StatisticsHelper.ExpM1(1e-20), Is.EqualTo(1e-20));
            Assert.That(StatisticsHelper.Log1P(-0.5), Is.EqualTo(-0.6931471805599453).Within(1e-15));
            Assert.That(StatisticsHelper.ExpM1(1.0), Is.EqualTo(1.718281828459045).Within(1e-15));
            Assert.That(StatisticsHelper.Log1P(-1.0), Is.EqualTo(double.NegativeInfinity));
            Assert.That(StatisticsHelper.ExpM1(-800), Is.EqualTo(-1.0));
        });
    }

    [Test]
    public void LogAddExp_MatchesNumpyLogaddexp()
    {
        // numpy 2.x np.logaddexp: (-1000, -1001) → -999.6867383124818; (ln ¼, ln ½) → ln ¾ = -0.2876820724517809;
        // (-inf, -inf) → -inf; (-inf, -3.5) → -3.5; (700, 710) → 710.0000453988993 (no overflow).
        Assert.Multiple(() =>
        {
            Assert.That(StatisticsHelper.LogAddExp(-1000.0, -1001.0), Is.EqualTo(-999.6867383124818).Within(1e-12));
            Assert.That(StatisticsHelper.LogAddExp(Math.Log(0.25), Math.Log(0.5)), Is.EqualTo(-0.2876820724517809).Within(1e-15));
            Assert.That(StatisticsHelper.LogAddExp(double.NegativeInfinity, double.NegativeInfinity), Is.EqualTo(double.NegativeInfinity));
            Assert.That(StatisticsHelper.LogAddExp(double.NegativeInfinity, -3.5), Is.EqualTo(-3.5));
            Assert.That(StatisticsHelper.LogAddExp(-3.5, double.NegativeInfinity), Is.EqualTo(-3.5));
            Assert.That(StatisticsHelper.LogAddExp(700.0, 710.0), Is.EqualTo(710.0000453988993).Within(1e-12));
        });
    }

    // numpy 2.4.6 np.logaddexp (npy_logaddexp: x == y → x + ln 2, handling same-sign infinities), printed with repr.
    private static readonly object[] NumpyLogAddExpSpecialCases =
    {
        new object[] { double.PositiveInfinity, double.PositiveInfinity, double.PositiveInfinity },
        new object[] { double.NegativeInfinity, double.NegativeInfinity, double.NegativeInfinity },
        new object[] { double.PositiveInfinity, double.NegativeInfinity, double.PositiveInfinity },
        new object[] { double.NegativeInfinity, double.PositiveInfinity, double.PositiveInfinity },
        new object[] { double.PositiveInfinity, 1.0, double.PositiveInfinity },
        new object[] { 1.0, double.PositiveInfinity, double.PositiveInfinity },
        new object[] { double.NaN, 1.0, double.NaN },
        new object[] { 1.0, double.NaN, double.NaN },
        new object[] { double.NaN, double.NegativeInfinity, double.NaN },
        new object[] { double.NegativeInfinity, double.NaN, double.NaN },
        new object[] { double.PositiveInfinity, double.NaN, double.NaN },
        new object[] { double.NaN, double.PositiveInfinity, double.NaN },
        new object[] { 1.0, 1.0, 1.6931471805599454 },
        new object[] { -3.5, -3.5, -2.8068528194400546 },
        new object[] { 1e308, 1e308, 1e308 },
    };

    [TestCaseSource(nameof(NumpyLogAddExpSpecialCases))]
    public void LogAddExp_SpecialValues_EqualNumpy(double a, double b, double expected)
        => Assert.That(StatisticsHelper.LogAddExp(a, b), Is.EqualTo(expected));

    [Test]
    public void LogAddExp_NearEqualFinite_MatchesNumpy()
    {
        // numpy 2.4.6 np.logaddexp, repr-printed.
        Assert.Multiple(() =>
        {
            Assert.That(StatisticsHelper.LogAddExp(1.0, 1.0 + Math.Pow(2, -52)), Is.EqualTo(1.6931471805599454).Within(1e-15));
            Assert.That(StatisticsHelper.LogAddExp(-700.0, -700.0000000001), Is.EqualTo(-699.3068528194901).Within(1e-12));
            Assert.That(StatisticsHelper.LogAddExp(0.5, 0.25), Is.EqualTo(1.0759394198788437).Within(1e-15));
        });
    }

    [Test]
    public void LogAddExp_EqualFinite_BitIdenticalToMaxPlusLogSum()
    {
        // a == b branch (numpy a + ln 2) equals the general formula max + ln(e^0 + e^0) bit for bit for finite a.
        var rng = new Random(7);
        for (int i = 0; i < 10000; i++)
        {
            double a = (rng.NextDouble() - 0.5) * Math.Pow(10, rng.Next(-20, 300));
            double general = a + Math.Log(Math.Exp(a - a) + Math.Exp(a - a));
            Assert.That(BitConverter.DoubleToInt64Bits(StatisticsHelper.LogAddExp(a, a)), Is.EqualTo(BitConverter.DoubleToInt64Bits(general)));
        }
    }

    private static double LogChoose(long n, long k)
    {
        double s = 0;
        for (long i = 1; i <= k; i++) s += Math.Log(n - k + i) - Math.Log(i);
        return s;
    }
}
