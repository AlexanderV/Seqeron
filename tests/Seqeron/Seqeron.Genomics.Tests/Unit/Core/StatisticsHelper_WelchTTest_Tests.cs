// StatisticsHelper.WelchTTestPValue — R t.test(x, y)$p.value (var.equal = FALSE), added for Battenberg merge_segments'
// logR / BAF tests between adjacent segments (ONCO-ASCAT-001, B24 F60). Reference values: R 4.3.3 t.test printed with
// %.17g (Evidence ONCO-ASCAT-001 §F60).

using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Core;

[TestFixture]
public class StatisticsHelper_WelchTTest_Tests
{
    private static void AssertRel(double actual, double expected, string message = "")
        => Assert.That(actual, Is.EqualTo(expected).Within(1e-13 * Math.Abs(expected)), message);

    [Test]
    public void WelchTTestPValue_MatchesRTTest()
    {
        Assert.Multiple(() =>
        {
            // ν = 6.9010638100365718 (non-integer Welch–Satterthwaite df).
            AssertRel(StatisticsHelper.WelchTTestPValue(new[] { 0.61, 0.64, 0.66, 0.59, 0.70 }, new[] { 0.55, 0.58, 0.62, 0.57 }),
                0.042862029945299411, "small samples");
            AssertRel(StatisticsHelper.WelchTTestPValue(new[] { 1.0, 2.0 }, new[] { 1.5, 2.5, 3.5 }), 0.28502284001427436, "n = 2");
            AssertRel(StatisticsHelper.WelchTTestPValue(new[] { -3.2, 1.1, 0.4, 2.2, -0.7, 0.05 }, new[] { 10.1, 9.8, 10.4 }),
                2.2976696502632257e-05, "far apart");
            var x = Enumerable.Range(1, 200).Select(i => Math.Sin(i)).ToArray();
            var y = Enumerable.Range(1, 150).Select(i => Math.Cos(i) * 3 + 0.1).ToArray();
            AssertRel(StatisticsHelper.WelchTTestPValue(x, y), 0.64314653287604995, "n = 200 / 150");
            AssertRel(StatisticsHelper.WelchTTestPValue(new[] { 1e6 + 1, 1e6 + 2, 1e6 + 3 }, new[] { 1e6 + 1.5, 1e6 + 2.5 }), 1.0,
                "equal means, large offset");
        });
    }

    [Test]
    public void WelchTTestPValue_RErrorCases_AreNaN()
    {
        Assert.Multiple(() =>
        {
            // R: "data are essentially constant".
            Assert.That(StatisticsHelper.WelchTTestPValue(new[] { 2.0, 2.0, 2.0 }, new[] { 2.0, 2.0 }), Is.NaN);
            // R: "not enough 'x' / 'y' observations".
            Assert.That(StatisticsHelper.WelchTTestPValue(new[] { 1.0 }, new[] { 2.0, 3.0 }), Is.NaN);
            Assert.That(StatisticsHelper.WelchTTestPValue(new[] { 1.0, 2.0 }, Array.Empty<double>()), Is.NaN);
        });
    }

    [Test]
    public void WelchTTestPValue_InvalidInput_Throws()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentNullException>(() => StatisticsHelper.WelchTTestPValue(null!, new[] { 1.0, 2.0 }));
            Assert.Throws<ArgumentNullException>(() => StatisticsHelper.WelchTTestPValue(new[] { 1.0, 2.0 }, null!));
            Assert.Throws<ArgumentException>(() => StatisticsHelper.WelchTTestPValue(new[] { 1.0, double.NaN }, new[] { 1.0, 2.0 }));
            Assert.Throws<ArgumentException>(() => StatisticsHelper.WelchTTestPValue(new[] { 1.0, 2.0 }, new[] { 1.0, double.PositiveInfinity }));
        });
    }
}
