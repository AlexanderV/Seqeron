using SuffixTree.Esa.Bench;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.Bench;

/// <summary>The arithmetic of the paired A/B protocol (the runner is <c>apps/SuffixTree.Esa.Bench/AbRunner.cs</c>; this file pins <c>PairedStats.cs</c> on synthetic data).</summary>
[TestFixture]
public class PairedStatsTests
{
    [Test]
    public void MedianAndMad()
    {
        Assert.That(PairedStats.Median([3, 1, 2]), Is.EqualTo(2));
        Assert.That(PairedStats.Median([4, 1, 2, 3]), Is.EqualTo(2.5));
        Assert.That(PairedStats.Mad([1, 2, 3, 4, 100]), Is.EqualTo(1)); // |x - 3| = 2 1 0 1 97
        Assert.Throws<ArgumentException>(() => PairedStats.Median([]));
    }

    [Test]
    public void Ratios_AreBOverA_PairByPair()
    {
        Assert.That(PairedStats.Ratios([1.0, 2.0, 4.0], [0.5, 2.0, 8.0]), Is.EqualTo(new[] { 0.5, 1.0, 2.0 }));
        Assert.Throws<ArgumentException>(() => PairedStats.Ratios([1.0], [1.0, 2.0]));
    }

    [Test]
    public void NoiseBand_IsTheOffsetPlusTwiceTheStandardErrorOfTheMedian()
    {
        Assert.That(PairedStats.NoiseBand([1, 1, 1, 1, 1, 1, 1]), Is.EqualTo((1.0, 1.0)));
        var (lowBias, highBias) = PairedStats.NoiseBand([1.02, 1.02, 1.02, 1.02, 1.02]);
        Assert.That(lowBias, Is.EqualTo(0.98).Within(1e-12));
        Assert.That(highBias, Is.EqualTo(1.02).Within(1e-12));
        // median 1.0, MAD 0.01, sigma 0.014826, se = 1.2533 * sigma / sqrt(7), half width 2 se
        var (low, high) = PairedStats.NoiseBand([0.98, 1.0, 1.02, 0.99, 1.01, 1.0, 1.0]);
        double h = 2 * 1.2533 * 1.4826 * 0.01 / Math.Sqrt(7);
        Assert.That(low, Is.EqualTo(1 - h).Within(1e-9));
        Assert.That(high, Is.EqualTo(1 + h).Within(1e-9));
    }

    [Test]
    public void MedianInterval_UsesOrderStatisticsWithAtLeast95PercentCoverage()
    {
        Assert.That(PairedStats.MedianInterval([1, 2, 3, 4, 5, 6, 7, 8, 9]), Is.EqualTo((2.0, 8.0)));
        Assert.That(PairedStats.MedianInterval([5, 4, 3, 2, 1]), Is.EqualTo((1.0, 5.0)), "fewer than 6 values: the full range");
        Assert.That(PairedStats.MedianInterval([1, 2, 3, 4, 5, 6, 7]), Is.EqualTo((1.0, 7.0)));
        Assert.That(PairedStats.MedianInterval(Enumerable.Range(1, 11).Select(i => (double)i).ToArray()), Is.EqualTo((2.0, 10.0)));
    }

    [Test]
    public void Decide_NeedsAtLeast3Percent_APositionOutsideTheNoiseBand_AndAnIntervalBelowOne()
    {
        (double, double) band = (0.985, 1.015), below = (0.90, 0.99), reaching = (0.90, 1.01);
        Assert.That(PairedStats.Decide(0.94, band, below).Accept, Is.True);
        Assert.That(PairedStats.Decide(0.97, band, below), Is.EqualTo((true, true, "median paired ratio 0.9700 (B/A), acceptance bar <= 0.97, A/A noise band [0.9850, 1.0150], 95 % interval of the median [0.9000, 0.9900]: at least 3 % faster, outside the noise band, interval below 1; BORDERLINE (within 1 percentage point of the bar): repeat with more pairs")), "exactly 3 %");
        Assert.That(PairedStats.Decide(0.975, band, below).Accept, Is.False, "less than 3 %");
        Assert.That(PairedStats.Decide(0.96, (0.90, 1.10), below).Reason, Does.Contain("inside the noise band"));
        Assert.That(PairedStats.Decide(0.94, band, reaching).Reason, Does.Contain("reaches 1.0"));
        Assert.That(PairedStats.Decide(0.94, band, reaching).Accept, Is.False);
        Assert.That(PairedStats.Decide(1.04, band, (1.0, 1.1)).Accept, Is.False, "slower");
        Assert.That(PairedStats.Decide(0.90, band, below).Borderline, Is.False);
        Assert.That(PairedStats.Decide(0.9726, band, below).Borderline, Is.True);
    }

    [Test]
    public void TheRule_OnSimulatedSessions_NeverAcceptsNoChange_AndAlmostAlwaysAcceptsARealFivePercent()
    {
        // per-pair ratio noise 2.4 % (two runs of 1.7 % each, as measured in the S1a review; the live A/A session of 1o.0 saw up to 3.5 %, which the 15-pair sessions of the iterations take into account), lognormal; 7 A/A pairs and 9 A/B pairs per session
        var rng = new SplitMix64(0xAB1234UL);
        double Gauss()
        {
            double u1 = (((rng.Next() >> 11) + 1) / 9007199254740992.0), u2 = (rng.Next() >> 11) / 9007199254740992.0;
            return Math.Sqrt(-2 * Math.Log(u1)) * Math.Cos(2 * Math.PI * u2);
        }
        List<double> Sample(int n, double trueRatio) => [.. Enumerable.Range(0, n).Select(_ => trueRatio * Math.Exp(0.024 * Gauss()))];

        static bool Accepts(List<double> ratios, (double Low, double High) band) => PairedStats.Decide(PairedStats.Median(ratios), band, PairedStats.MedianInterval(ratios)).Accept;
        int falseAccepts = 0, hits = 0, mid = 0;
        const int sessions = 2000;
        for (int s = 0; s < sessions; s++)
        {
            var band = PairedStats.NoiseBand(Sample(7, 1.0));
            if (Accepts(Sample(9, 1.0), band)) falseAccepts++;
            if (Accepts(Sample(9, 0.95), band)) hits++;
            if (Accepts(Sample(9, 0.90), band)) mid++;
        }
        TestContext.Out.WriteLine($"no change: {falseAccepts}/{sessions} accepted; true -5 %: {hits}/{sessions}; true -10 %: {mid}/{sessions}");
        Assert.That(falseAccepts, Is.LessThanOrEqualTo(sessions / 100), "false acceptance of an unchanged kernel");
        Assert.That(mid, Is.GreaterThanOrEqualTo(sessions * 99 / 100), "a 10 % gain must be accepted");
        Assert.That(hits, Is.GreaterThanOrEqualTo(sessions * 85 / 100), "a real 5 % gain must be accepted in most sessions");
    }
}
