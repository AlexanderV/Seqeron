using NUnit.Framework;
using SuffixTree.Mcp.Sequence.Tools;

namespace SuffixTree.Mcp.Sequence.Tests;

[TestFixture]
public class KmerEntropyTests
{
    [Test]
    public void KmerEntropy_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => SequenceTools.KmerEntropy("ATGCGATCGATCG", 2));
        Assert.Throws<ArgumentException>(() => SequenceTools.KmerEntropy("", 2));
        Assert.Throws<ArgumentException>(() => SequenceTools.KmerEntropy(null!, 2));
        Assert.Throws<ArgumentException>(() => SequenceTools.KmerEntropy("ATGC", 0)); // k < 1
    }

    [Test]
    public void KmerEntropy_Binding_InvokesSuccessfully()
    {
        var result = SequenceTools.KmerEntropy("ATGCGATCGATCG", 2);
        Assert.That(result.Entropy, Is.GreaterThanOrEqualTo(0));
        Assert.That(result.K, Is.EqualTo(2));

        // Low complexity sequence has lower entropy
        var lowComplexity = SequenceTools.KmerEntropy("AAAAAAAAAA", 2);
        var highComplexity = SequenceTools.KmerEntropy("ATGCATGCAT", 2);
        Assert.That(highComplexity.Entropy, Is.GreaterThan(lowComplexity.Entropy));
    }
}
