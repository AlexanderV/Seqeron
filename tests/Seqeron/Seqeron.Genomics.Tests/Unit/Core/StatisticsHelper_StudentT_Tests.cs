// StatisticsHelper.RegularizedIncompleteBeta / StudentTCdf / OneSampleTTestPValue — shared Student-t machinery added for
// Battenberg's per-segment SNP-BAF t-test (ONCO-ASCAT-001, B24 F39). Reference values: R 4.3.3 pbeta / pt /
// t.test(x, mu)$p.value printed with %.17g (Evidence ONCO-ASCAT-001 §F39).

using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Core;

[TestFixture]
public class StatisticsHelper_StudentT_Tests
{
    private const double RelTol = 1e-13;

    private static void AssertRel(double actual, double expected, string message = "")
        => Assert.That(actual, Is.EqualTo(expected).Within(RelTol * Math.Max(Math.Abs(expected), 1e-300)), message);

    [TestCase(0.3, 2.0, 3.0, 0.34829999999999989)]
    [TestCase(0.5, 0.5, 0.5, 0.49999999999999956)]
    [TestCase(0.9, 10.0, 0.5, 0.15164090963470997)]
    [TestCase(0.999, 500.0, 0.5, 0.31731044730971736)]
    [TestCase(0.2, 1e-3, 5.0, 0.99978344564642763)]
    [TestCase(0.7, 25.0, 30.0, 0.99991848204532263)]
    [TestCase(0.01, 0.5, 200.0, 0.95490589286292027)]
    [TestCase(0.5, 1e4, 1e4, 0.5)]
    [TestCase(0.97, 150.0, 4.5, 0.41508270595853397)]
    public void RegularizedIncompleteBeta_MatchesRPbeta(double x, double a, double b, double expected)
    {
        AssertRel(StatisticsHelper.RegularizedIncompleteBeta(x, a, b), expected);
    }

    [Test]
    public void RegularizedIncompleteBeta_Endpoints()
    {
        Assert.Multiple(() =>
        {
            Assert.That(StatisticsHelper.RegularizedIncompleteBeta(0.0, 2.0, 3.0), Is.EqualTo(0.0));
            Assert.That(StatisticsHelper.RegularizedIncompleteBeta(1.0, 2.0, 3.0), Is.EqualTo(1.0));
        });
    }

    [TestCase(-2.5, 4.0, 0.033383272405994063)]
    [TestCase(2.5, 4.0, 0.96661672759400596)]
    [TestCase(0.0, 7.0, 0.5)]
    [TestCase(-0.3, 1.0, 0.40722642092225769)]
    [TestCase(1.96, 1000.0, 0.97486340752212564)]
    [TestCase(-12.0, 3.0, 0.00062250790039466801)]
    [TestCase(-40.0, 25.0, 1.7349901388795852e-24)]
    [TestCase(-3.1, 2.5, 0.033894189162460162)]
    [TestCase(0.75, 9999.0, 0.77336382453143315)]
    [TestCase(-6.0, 200000.0, 9.8827494396739874e-10)]
    [TestCase(-1e-8, 10.0, 0.49999999610891616)]
    public void StudentTCdf_MatchesRPt(double t, double df, double expected)
    {
        AssertRel(StatisticsHelper.StudentTCdf(t, df), expected);
    }

    [Test]
    public void StudentTCdf_Infinities()
    {
        Assert.Multiple(() =>
        {
            Assert.That(StatisticsHelper.StudentTCdf(double.NegativeInfinity, 3.0), Is.EqualTo(0.0));
            Assert.That(StatisticsHelper.StudentTCdf(double.PositiveInfinity, 3.0), Is.EqualTo(1.0));
        });
    }

    [TestCase(new[] { 0.61, 0.64, 0.66, 0.59, 0.70 }, 0.6, 0.10608282922154902)]
    [TestCase(new[] { 0.61, 0.64 }, 0.6, 0.34404173924526132)]
    [TestCase(new[] { 0.71, 0.69, 0.74, 0.73, 0.70, 0.72, 0.75, 0.68 }, 0.6666666666666666, 0.00083230863131458646)]
    [TestCase(new[] { -1.2, 3.4, 0.5, 2.2, -0.7, 1.1 }, 0.0, 0.26838941143845474)]
    [TestCase(new[] { 0.5, 0.5000001 }, 0.2, 1.0610327765486823e-07)]
    public void OneSampleTTestPValue_MatchesRTTest(double[] values, double mu, double expected)
    {
        Assert.That(StatisticsHelper.OneSampleTTestPValue(values, mu), Is.EqualTo(expected).Within(1e-12 * expected));
    }

    // R t.test stops with an error for n < 2 and for "data are essentially constant"; the helper returns NaN there.
    [Test]
    public void OneSampleTTestPValue_RErrorCases_ReturnNaN()
    {
        Assert.Multiple(() =>
        {
            Assert.That(StatisticsHelper.OneSampleTTestPValue(new[] { 0.6 }, 0.5), Is.NaN, "n = 1.");
            Assert.That(StatisticsHelper.OneSampleTTestPValue(Array.Empty<double>(), 0.5), Is.NaN, "n = 0.");
            Assert.That(StatisticsHelper.OneSampleTTestPValue(new[] { 0.7, 0.7, 0.7 }, 0.5), Is.NaN, "Constant data.");
        });
    }

    [Test]
    public void InvalidArguments_Throw()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => StatisticsHelper.RegularizedIncompleteBeta(1.5, 1, 1));
            Assert.Throws<ArgumentOutOfRangeException>(() => StatisticsHelper.RegularizedIncompleteBeta(0.5, 0, 1));
            Assert.Throws<ArgumentOutOfRangeException>(() => StatisticsHelper.RegularizedIncompleteBeta(0.5, 1, double.NaN));
            Assert.Throws<ArgumentOutOfRangeException>(() => StatisticsHelper.StudentTCdf(1.0, 0.0));
            Assert.Throws<ArgumentOutOfRangeException>(() => StatisticsHelper.StudentTCdf(double.NaN, 3.0));
            Assert.Throws<ArgumentNullException>(() => StatisticsHelper.OneSampleTTestPValue(null!, 0.0));
            Assert.Throws<ArgumentException>(() => StatisticsHelper.OneSampleTTestPValue(new[] { 0.1, double.NaN }, 0.0));
        });
    }
}
