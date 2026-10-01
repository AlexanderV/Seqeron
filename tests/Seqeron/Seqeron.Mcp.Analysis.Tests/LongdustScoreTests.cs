using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>longdust_score</c> MCP tool. Expected = S_L(x) from lh3/longdust 1.4-r97's own f() table
/// (ld_cal_f / ld_cal_f2), locked in SequenceComplexity_CalculateDustScore_Tests. NOT the wrapper's output.
/// </summary>
[TestFixture]
public class LongdustScoreTests
{
    [Test]
    public void LongdustScore_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.LongdustScore("ACGTACGT"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.LongdustScore(""));
        Assert.Throws<ArgumentException>(() => AnalysisTools.LongdustScore(null!));
        Assert.Throws<ArgumentException>(() => AnalysisTools.LongdustScore("XYZ"));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.LongdustScore("ACGT", k: 15));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.LongdustScore("ACGT", gcContent: 0.0));
    }

    [TestCase("AAAAAAAAAA", 3, null, 10.263913947178825)]
    [TestCase("AAAAAAAAAA", 3, 0.41, 10.230970927904716)]
    [TestCase("ACACACACACACACACACACACACACACAC", 7, null, 39.962247232206003)]
    [TestCase("ACGTNACGTACGTACGT", 4, null, 4.9941414253104837)]
    public void LongdustScore_Binding_MatchesLongdustF(string seq, int k, double? gc, double expected)
    {
        var r = AnalysisTools.LongdustScore(seq, k, gc);
        Assert.That(r.Score, Is.EqualTo(expected).Within(1e-12));
        Assert.That(r.KmerPositions, Is.EqualTo(seq.Length - k + 1));
    }

    [Test]
    public void LongdustScore_NoKmerPosition_IsZero()
    {
        var r = AnalysisTools.LongdustScore("ACGTAC");
        Assert.That((r.Score, r.KmerPositions), Is.EqualTo((0.0, 0)));
    }
}
