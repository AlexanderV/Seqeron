using NUnit.Framework;
using Seqeron.Mcp.Sequence.Tools;

namespace Seqeron.Mcp.Sequence.Tests;

[TestFixture]
public class ComplexityShannonTests
{
    [Test]
    public void ComplexityShannon_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => SequenceTools.ComplexityShannon("ATGCGATCGATCG"));
        Assert.Throws<ArgumentException>(() => SequenceTools.ComplexityShannon(""));
        Assert.Throws<ArgumentException>(() => SequenceTools.ComplexityShannon(null!));
    }

    [Test]
    public void ComplexityShannon_Binding_InvokesSuccessfully()
    {
        var result = SequenceTools.ComplexityShannon("ATGCGATCGATCG");
        Assert.That(result.Entropy, Is.GreaterThanOrEqualTo(0));
        Assert.That(result.Entropy, Is.LessThanOrEqualTo(2)); // Max for DNA is 2 bits

        // Low complexity (repetitive) sequence has lower entropy
        var lowComplexity = SequenceTools.ComplexityShannon("AAAAAAAAAAAAAAAA");
        var highComplexity = SequenceTools.ComplexityShannon("ATGCATGCATGCATGC");
        Assert.That(highComplexity.Entropy, Is.GreaterThan(lowComplexity.Entropy));
    }
}
