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
        // CompressionRatio is the normalized Lempel-Ziv complexity (SEQ-COMPLEX-COMPRESS-001):
        // c / (n / log_b n). It is positive for any non-trivial sequence and may exceed 1 for
        // finite sequences (the random-sequence asymptote is 1; cf. Hu et al. 2006), so no <=1 bound.
        Assert.That(result.CompressionRatio, Is.GreaterThan(0));

        // A periodic sequence has lower compression-based complexity (fewer Lempel-Ziv
        // components) than a diverse one. Both compared sequences use the SAME length (40)
        // and the SAME 4-letter alphabet so the log_b(n) normalization factor is identical,
        // isolating repetitiveness as the only difference.
        var lowComplexity = SequenceTools.ComplexityCompressionRatio("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT");
        var highComplexity = SequenceTools.ComplexityCompressionRatio("ACAGGTCATGCTAGGCATTCGATCGATGCCATGTCAGCTA");
        Assert.That(highComplexity.CompressionRatio, Is.GreaterThan(lowComplexity.CompressionRatio));
    }
}
