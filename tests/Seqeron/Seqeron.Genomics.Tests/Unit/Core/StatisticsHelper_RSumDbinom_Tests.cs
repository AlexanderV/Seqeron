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
}
