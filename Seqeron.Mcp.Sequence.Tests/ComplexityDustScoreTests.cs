using NUnit.Framework;
using Seqeron.Mcp.Sequence.Tools;

namespace Seqeron.Mcp.Sequence.Tests;

[TestFixture]
public class ComplexityDustScoreTests
{
    [Test]
    public void ComplexityDustScore_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => SequenceTools.ComplexityDustScore("ATGCGATCGATCG", 3));
        Assert.Throws<ArgumentException>(() => SequenceTools.ComplexityDustScore("", 3));
        Assert.Throws<ArgumentException>(() => SequenceTools.ComplexityDustScore(null!, 3));
        Assert.Throws<ArgumentException>(() => SequenceTools.ComplexityDustScore("ATGC", 0)); // wordSize < 1
    }

    [Test]
    public void ComplexityDustScore_Binding_InvokesSuccessfully()
    {
        var result = SequenceTools.ComplexityDustScore("ATGCGATCGATCG", 3);
        Assert.That(result.DustScore, Is.GreaterThanOrEqualTo(0));
        Assert.That(result.WordSize, Is.EqualTo(3));

        // Low complexity (repetitive) sequence has higher DUST score
        var lowComplexity = SequenceTools.ComplexityDustScore("AAAAAAAAAAAAAAAAAA", 3);
        var highComplexity = SequenceTools.ComplexityDustScore("ATGCATGCATGCATGCAT", 3);
        Assert.That(lowComplexity.DustScore, Is.GreaterThan(highComplexity.DustScore));
    }
}
