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
        // B04 F34: DUST is defined for triplets only (Morgulis 2006; symdust; lh3/sdust SD_WLEN = 3).
        Assert.Throws<ArgumentException>(() => SequenceTools.ComplexityDustScore("ATGC", 2));
        Assert.Throws<ArgumentException>(() => SequenceTools.ComplexityDustScore("ATGC", 4));
    }

    [Test]
    public void ComplexityDustScore_DocumentedExamples()
    {
        // Σ c(c−1)/2 / (ℓ − 1): ATGCGATCGATCG → CGA, GAT, ATC, TCG twice: 4 / 10; A×12 → 45 / 9.
        Assert.That(SequenceTools.ComplexityDustScore("ATGCGATCGATCG", 3).DustScore, Is.EqualTo(0.4).Within(1e-12));
        Assert.That(SequenceTools.ComplexityDustScore("AAAAAAAAAAAA", 3).DustScore, Is.EqualTo(5.0).Within(1e-12));
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
