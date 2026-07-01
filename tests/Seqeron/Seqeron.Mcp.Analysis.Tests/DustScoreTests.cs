using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>dust_score</c> MCP tool.
/// Expected values computed by hand from SequenceComplexity.CalculateDustScore:
/// score = sum_t c_t*(c_t-1)/2 / (L - wordSize + 1). NOT from the wrapper's output.
/// </summary>
[TestFixture]
public class DustScoreTests
{
    [Test]
    public void DustScore_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.DustScore("AAAAA", 3));
        Assert.Throws<ArgumentException>(() => AnalysisTools.DustScore("", 3));
        Assert.Throws<ArgumentException>(() => AnalysisTools.DustScore(null!, 3));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.DustScore("AAAAA", 0));
    }

    [Test]
    public void DustScore_Binding_InvokesSuccessfully()
    {
        Assert.Multiple(() =>
        {
            // "AAAAA" wordSize 3: words AAA,AAA,AAA (count 3), wordCount = 3.
            //   sum = 3*2/2 = 3; score = 3/3 = 1.0.
            Assert.That(AnalysisTools.DustScore("AAAAA", 3).Score, Is.EqualTo(1.0).Within(1e-12));

            // "ACGT" wordSize 3: words ACG,CGT (all distinct) -> sum 0 -> score 0.
            Assert.That(AnalysisTools.DustScore("ACGT", 3).Score, Is.EqualTo(0.0).Within(1e-12));

            // "ACGTACGT" wordSize 3: ACG=2,CGT=2,GTA=1,TAC=1, wordCount = 6.
            //   sum = 1 + 1 = 2; score = 2/6 = 0.3333...
            Assert.That(AnalysisTools.DustScore("ACGTACGT", 3).Score, Is.EqualTo(2.0 / 6.0).Within(1e-12));
        });
    }
}
