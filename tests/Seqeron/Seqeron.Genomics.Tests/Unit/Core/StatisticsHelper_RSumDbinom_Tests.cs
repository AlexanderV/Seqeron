// StatisticsHelper.ExtendedPrecisionSum — R sum() / colSums() / matprod "internal" 80-bit long-double accumulation — and
// StatisticsHelper.BinomialLogDensity — R dbinom(x, n, p, log = TRUE) (nmath dbinom_raw), added for the BMix port behind
// CNAqc's mixture peaks (ONCO-PURITY-001, B24 F63). Reference values: R 4.3.3 on x86-64 printed with %.17g
// (Evidence ONCO-PURITY-001 § "BMix and n_bootstrap").

using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Core;

[TestFixture]
public class StatisticsHelper_RSumDbinom_Tests
{
    // R: sum(c(1, 2^-53, 2^-53)) = 1.0000000000000002 (the halves survive in long double; a double loop gives 1);
    // sum(c(1, 2^-53, 2^-64)) = 1 (x87 double rounding: 64-bit tie to even, then 53-bit tie to even — the correctly
    // rounded sum would be 1 + 2^-52); sum(c(1e16, 1, 1)) = 10000000000000002 (double loop: 1e16).
    [TestCase(new[] { 1.0, 1.1102230246251565e-16, 1.1102230246251565e-16 }, 1.0000000000000002)]
    [TestCase(new[] { 1.0, 1.1102230246251565e-16, 5.4210108624275222e-20 }, 1.0)]
    [TestCase(new[] { 1e16, 1.0, 1.0 }, 10000000000000002.0)]
    [TestCase(new[] { 0.1, 0.2, 0.3 }, 0.59999999999999998)]
    [TestCase(new[] { -1.0, 1e-30 }, -1.0)]
    [TestCase(new[] { 3.0, -3.0 }, 0.0)]
    [TestCase(new[] { 4.9406564584124654e-324, 4.9406564584124654e-324 }, 9.8813129168249309e-324)]
    public void ExtendedPrecisionSum_MatchesRSum(double[] values, double expected)
    {
        Assert.That(StatisticsHelper.ExtendedPrecisionSum(values), Is.EqualTo(expected));
    }

    // R: sum(c(1e308, 1e308)) = Inf; sum(c(1, NaN)) = NaN; sum(c(1, NaN, 2), na.rm = TRUE) = 3; sum(c(Inf, -Inf)) = NaN.
    [Test]
    public void ExtendedPrecisionSum_SpecialValues_MatchR()
    {
        Assert.That(StatisticsHelper.ExtendedPrecisionSum(new[] { 1e308, 1e308 }), Is.EqualTo(double.PositiveInfinity));
        Assert.That(StatisticsHelper.ExtendedPrecisionSum(new[] { 1.0, double.NaN }), Is.NaN);
        Assert.That(StatisticsHelper.ExtendedPrecisionSum(new[] { 1.0, double.NaN, 2.0 }, skipNaN: true), Is.EqualTo(3.0));
        Assert.That(StatisticsHelper.ExtendedPrecisionSum(new[] { double.PositiveInfinity, double.NegativeInfinity }), Is.NaN);
        Assert.That(StatisticsHelper.ExtendedPrecisionSum(ReadOnlySpan<double>.Empty), Is.EqualTo(0.0));
    }

    // R dbinom(..., log = TRUE): interior (saddle point), x = 0 with p < 0.1 (bd0 branch) and p ≥ 0.1 (n·log q), x = n
    // with q < 0.1 and q ≥ 0.1, large n. ≤ 2 ulp-level differences come only from log1p (C library vs managed).
    [TestCase(3, 81, 0.2, -10.879346567225376)]
    [TestCase(0, 50, 0.05, -2.5646647193775265)]
    [TestCase(0, 50, 0.5, -34.657359027997266)]
    [TestCase(50, 50, 0.95, -2.5646647193775287)]
    [TestCase(50, 50, 0.5, -34.657359027997266)]
    [TestCase(120, 400, 0.3, -3.1351306610944456)]
    [TestCase(7000, 20000, 0.35, -5.1303939375288108)]
    public void BinomialLogDensity_MatchesRDbinom(double x, double n, double p, double expected)
    {
        Assert.That(StatisticsHelper.BinomialLogDensity(x, n, p), Is.EqualTo(expected).Within(4e-16 * Math.Abs(expected)));
    }

    // R: dbinom(5, 10, 0) = -Inf, dbinom(0, 10, 0) = 0, dbinom(10, 10, 1) = 0, dbinom(11, 10, .5) = -Inf,
    // dbinom(2.5, 10, .5) = -Inf (non-integer x, warning), dbinom(2, 10, 1.5) = NaN, dbinom(0, 0, .3) = 0.
    [Test]
    public void BinomialLogDensity_Boundaries_MatchR()
    {
        Assert.That(StatisticsHelper.BinomialLogDensity(5, 10, 0), Is.EqualTo(double.NegativeInfinity));
        Assert.That(StatisticsHelper.BinomialLogDensity(0, 10, 0), Is.EqualTo(0.0));
        Assert.That(StatisticsHelper.BinomialLogDensity(10, 10, 1), Is.EqualTo(0.0));
        Assert.That(StatisticsHelper.BinomialLogDensity(11, 10, 0.5), Is.EqualTo(double.NegativeInfinity));
        Assert.That(StatisticsHelper.BinomialLogDensity(2.5, 10, 0.5), Is.EqualTo(double.NegativeInfinity));
        Assert.That(StatisticsHelper.BinomialLogDensity(2, 10, 1.5), Is.NaN);
        Assert.That(StatisticsHelper.BinomialLogDensity(0, 0, 0.3), Is.EqualTo(0.0));
        Assert.That(StatisticsHelper.BinomialLogDensity(double.NaN, 10, 0.3), Is.NaN);
    }

    // ── ExtendedPrecisionMean / ExtendedPrecisionVariance (B24 F66): R mean() (summary.c real_mean) and var() (cov.c),
    // long-double two-pass. R 4.3.3 x86-64, sprintf("%.17g"). Rows chosen where the double-accumulated forms differ:
    // c(0.28, 0.01, 0.37, 0.96): R 0.40499999999999997, Σx/n and RMean 0.405; the 6-value rows: RMean is 1 ulp off.
    // Scratch corpus (3004 random vectors, n 1–1000): mean 3004/3004, var and bw.nrd0 2851/2851 bit-identical
    // (Σx/n 2137, RMean 2719, RSampleVariance 1455, old bw.nrd0 2360).
    [TestCase(new[] { 0.28, 0.01, 0.37, 0.96 }, 0.40499999999999997)]
    [TestCase(new[] { 0.25645005441148117, 0.235205005518964, 0.04509703448406843, -0.17313677835742297, 0.015755779681899895, 0.3420512029727458 }, 0.12023704978528939)]
    [TestCase(new[] { 2372146.1335775964, 370148.2528286602, 1221543.1720335376, 2220276.004596779, -1565453.9652750245, 3387.17783466971 }, 770341.1292660364)]
    [TestCase(new[] { 0.1, 0.2, 0.3 }, 0.20000000000000001)]
    [TestCase(new[] { 1e308, 1e308, -1e308 }, 3.3333333333333332e+307)]
    [TestCase(new[] { 1e308, 1e308 }, 1e308)] // Σx overflows a double: R's Σ(x/n) fallback
    public void ExtendedPrecisionMean_MatchesRMean(double[] values, double expected)
    {
        Assert.That(StatisticsHelper.ExtendedPrecisionMean(values), Is.EqualTo(expected));
    }

    // R: mean(c(0.28, NA, 0.01, 0.37, NaN, 0.96), na.rm = TRUE) = 0.40499999999999997; mean(numeric(0)) = NaN;
    // mean(c(1, NA)) = NA; mean(c(Inf, 1)) = Inf; mean(c(Inf, -Inf)) = NaN.
    [Test]
    public void ExtendedPrecisionMean_SpecialValues_MatchR()
    {
        Assert.That(StatisticsHelper.ExtendedPrecisionMean(new[] { 0.28, double.NaN, 0.01, 0.37, double.NaN, 0.96 }, skipNaN: true),
            Is.EqualTo(0.40499999999999997));
        Assert.That(StatisticsHelper.ExtendedPrecisionMean(ReadOnlySpan<double>.Empty), Is.NaN);
        Assert.That(StatisticsHelper.ExtendedPrecisionMean(new[] { double.NaN }, skipNaN: true), Is.NaN);
        Assert.That(StatisticsHelper.ExtendedPrecisionMean(new[] { 1.0, double.NaN }), Is.NaN);
        Assert.That(StatisticsHelper.ExtendedPrecisionMean(new[] { double.PositiveInfinity, 1.0 }), Is.EqualTo(double.PositiveInfinity));
        Assert.That(StatisticsHelper.ExtendedPrecisionMean(new[] { double.PositiveInfinity, double.NegativeInfinity }), Is.NaN);
    }

    // R var(): c(0.52, 0.05, 0.08, 0.41, 0.58) = 0.061469999999999997 (RSampleVariance 0.061470000000000004), so
    // sd = 0.24793144213673263 and bw.nrd0 = 0.16172610052228839; c(1, 2, 3, 4) = 1.6666666666666667; one value → NA;
    // c(1, Inf) → NaN.
    [Test]
    public void ExtendedPrecisionVariance_MatchesRVar()
    {
        double[] y = { 0.52, 0.05, 0.08, 0.41, 0.58 };
        Assert.That(StatisticsHelper.ExtendedPrecisionVariance(y), Is.EqualTo(0.061469999999999997));
        Assert.That(Math.Sqrt(StatisticsHelper.ExtendedPrecisionVariance(y)), Is.EqualTo(0.24793144213673263));
        Assert.That(StatisticsHelper.BandwidthNrd0(y), Is.EqualTo(0.16172610052228839));
        Assert.That(StatisticsHelper.ExtendedPrecisionVariance(new[] { 1.0, 2, 3, 4 }), Is.EqualTo(1.6666666666666667));
        Assert.That(StatisticsHelper.ExtendedPrecisionVariance(new[] { 0.3 }), Is.NaN);
        Assert.That(StatisticsHelper.ExtendedPrecisionVariance(new[] { 1.0, double.PositiveInfinity }), Is.NaN);
    }
}
