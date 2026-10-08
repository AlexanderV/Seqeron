using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>dust_score</c> MCP tool.
/// Expected values computed by hand from SequenceComplexity.CalculateDustScore:
/// score = sum_t c_t*(c_t-1)/2 / (l - 1), l = L - wordSize + 1 words (Morgulis 2006;
/// NCBI symdust / lh3/sdust normaliser). NOT from the wrapper's output.
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
        // B04 F34: DUST is defined for triplets only (Morgulis 2006; symdust; lh3/sdust SD_WLEN = 3).
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.DustScore("AAAAA", 2));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.DustScore("AAAAA", 4));
    }

    [Test]
    public void DustScore_Binding_InvokesSuccessfully()
    {
        Assert.Multiple(() =>
        {
            // "AAAAA" wordSize 3: words AAA,AAA,AAA (count 3), l = 3.
            //   sum = 3*2/2 = 3; score = 3/(l-1) = 1.5.
            Assert.That(AnalysisTools.DustScore("AAAAA", 3).Score, Is.EqualTo(1.5).Within(1e-12));

            // "ACGT" wordSize 3: words ACG,CGT (all distinct) -> sum 0 -> score 0.
            Assert.That(AnalysisTools.DustScore("ACGT", 3).Score, Is.EqualTo(0.0).Within(1e-12));

            // "ACGTACGT" wordSize 3: ACG=2,CGT=2,GTA=1,TAC=1, l = 6.
            //   sum = 1 + 1 = 2; score = 2/(l-1) = 0.4.
            Assert.That(AnalysisTools.DustScore("ACGTACGT", 3).Score, Is.EqualTo(2.0 / 5.0).Within(1e-12));
        });
    }
}
