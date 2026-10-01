using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>tandem_repeat_bernoulli_statistics</c> MCP tool. Expected counts = instrumented TRF 4.10.0
/// get_statistics on the same tracts, locked in RepeatFinder_ApproximateTandemRepeats_Tests (B1–B6). NOT the wrapper's output.
/// </summary>
[TestFixture]
public class TandemRepeatBernoulliStatisticsTests
{
    [Test]
    public void TandemRepeatBernoulliStatistics_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.TandemRepeatBernoulliStatistics("CACACACACA", 2));
        Assert.Throws<ArgumentException>(() => AnalysisTools.TandemRepeatBernoulliStatistics("", 2));
        Assert.Throws<ArgumentException>(() => AnalysisTools.TandemRepeatBernoulliStatistics(null!, 2));
        Assert.Throws<ArgumentException>(() => AnalysisTools.TandemRepeatBernoulliStatistics("CAG", 3));
    }

    [Test]
    public void TandemRepeatBernoulliStatistics_Binding_OneSubstitution()
    {
        // CAG x6 with copy 4 = TAG: TRF m=13 mm=2 ind=0.
        var s = AnalysisTools.TandemRepeatBernoulliStatistics("CAGCAGCAGTAGCAGCAG", 3);
        Assert.Multiple(() =>
        {
            Assert.That((s.Period, s.AdjacentCopyPairs, s.BernoulliTrials, s.Matches, s.Mismatches, s.Indels),
                Is.EqualTo((3, 5, 15, 13, 2, 0)));
            Assert.That(s.MatchProbability, Is.EqualTo(13.0 / 15.0).Within(1e-12));
            Assert.That(s.ExpectedMatches, Is.EqualTo(13.0).Within(1e-12));
            Assert.That(s.MeetsExpectedMatchProbability, Is.True);
        });
        Assert.That(AnalysisTools.TandemRepeatBernoulliStatistics("CAGCAGCAGTAGCAGCAG", 3, 0.90).MeetsExpectedMatchProbability, Is.False);
    }

    [Test]
    public void TandemRepeatBernoulliStatistics_DeletionTract_MatchesTrf()
    {
        // TRF m=25 mm=0 ind=2, 9 adjacent pairs.
        var s = AnalysisTools.TandemRepeatBernoulliStatistics("CAGCAGCAGCAGCAGAGCAGCAGCAGCAG", 3);
        Assert.That((s.Matches, s.Mismatches, s.Indels, s.BernoulliTrials, s.AdjacentCopyPairs), Is.EqualTo((25, 0, 2, 27, 9)));
        Assert.That(s.IndelProbability, Is.EqualTo(2.0 / 27).Within(1e-12));
    }
}
