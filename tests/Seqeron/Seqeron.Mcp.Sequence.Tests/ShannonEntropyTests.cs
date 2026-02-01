using NUnit.Framework;
using Seqeron.Mcp.Sequence.Tools;

namespace Seqeron.Mcp.Sequence.Tests;

[TestFixture]
public class ShannonEntropyTests
{
    [Test]
    public void ShannonEntropy_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => SequenceTools.ShannonEntropy("ATGCGATCGATCG"));
        Assert.Throws<ArgumentException>(() => SequenceTools.ShannonEntropy(""));
        Assert.Throws<ArgumentException>(() => SequenceTools.ShannonEntropy(null!));
    }

    [Test]
    public void ShannonEntropy_Binding_InvokesSuccessfully()
    {
        var result = SequenceTools.ShannonEntropy("ATGCGATCGATCG");
        Assert.That(result.Entropy, Is.TypeOf<double>());
        Assert.That(result.Entropy, Is.GreaterThanOrEqualTo(0));

        // Low complexity sequence (homopolymer) has lower entropy
        var lowComplexity = SequenceTools.ShannonEntropy("AAAAAAAAAA");
        var highComplexity = SequenceTools.ShannonEntropy("ATGCATGCAT");
        Assert.That(highComplexity.Entropy, Is.GreaterThan(lowComplexity.Entropy));
    }
}
