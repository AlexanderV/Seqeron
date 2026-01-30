using NUnit.Framework;
using Seqeron.Mcp.Sequence.Tools;

namespace Seqeron.Mcp.Sequence.Tests;

[TestFixture]
public class ComplexityCompressionRatioTests
{
    [Test]
    public void ComplexityCompressionRatio_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => SequenceTools.ComplexityCompressionRatio("ATGCGATCGATCG"));
        Assert.Throws<ArgumentException>(() => SequenceTools.ComplexityCompressionRatio(""));
        Assert.Throws<ArgumentException>(() => SequenceTools.ComplexityCompressionRatio(null!));
    }

    [Test]
    public void ComplexityCompressionRatio_Binding_InvokesSuccessfully()
    {
        var result = SequenceTools.ComplexityCompressionRatio("ATGCGATCGATCG");
        Assert.That(result.CompressionRatio, Is.GreaterThan(0));
        Assert.That(result.CompressionRatio, Is.LessThanOrEqualTo(1));

        // Repetitive sequence has lower compression ratio (fewer unique substrings)
        var lowComplexity = SequenceTools.ComplexityCompressionRatio("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
        var highComplexity = SequenceTools.ComplexityCompressionRatio("ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC");
        Assert.That(highComplexity.CompressionRatio, Is.GreaterThan(lowComplexity.CompressionRatio));
    }
}
