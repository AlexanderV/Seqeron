using NUnit.Framework;
using Seqeron.Mcp.Sequence.Tools;

namespace Seqeron.Mcp.Sequence.Tests;

[TestFixture]
public class ComplexityKmerEntropyTests
{
    [Test]
    public void ComplexityKmerEntropy_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => SequenceTools.ComplexityKmerEntropy("ATGCGATCGATCG", 2));
        Assert.Throws<ArgumentException>(() => SequenceTools.ComplexityKmerEntropy("", 2));
        Assert.Throws<ArgumentException>(() => SequenceTools.ComplexityKmerEntropy(null!, 2));
        Assert.Throws<ArgumentException>(() => SequenceTools.ComplexityKmerEntropy("ATGC", 0)); // k < 1
        Assert.Throws<ArgumentException>(() => SequenceTools.ComplexityKmerEntropy("ZZZCGA", 2)); // invalid DNA
    }

    [Test]
    public void ComplexityKmerEntropy_Binding_InvokesSuccessfully()
    {
        var result = SequenceTools.ComplexityKmerEntropy("ATGCGATCGATCG", 2);
        Assert.That(result.Entropy, Is.GreaterThanOrEqualTo(0));
        Assert.That(result.K, Is.EqualTo(2));

        // Low complexity (repetitive) sequence has lower entropy
        var lowComplexity = SequenceTools.ComplexityKmerEntropy("AAAAAAAAAAAAAAAA", 2);
        var highComplexity = SequenceTools.ComplexityKmerEntropy("ATGCATGCATGCATGC", 2);
        Assert.That(highComplexity.Entropy, Is.GreaterThan(lowComplexity.Entropy));
    }
}
