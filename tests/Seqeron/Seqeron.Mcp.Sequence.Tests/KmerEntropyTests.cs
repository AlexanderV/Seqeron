using NUnit.Framework;
using Seqeron.Mcp.Sequence.Tools;

namespace Seqeron.Mcp.Sequence.Tests;

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

    [Test]
    public void KmerEntropy_DocExample_MatchesScipy()
    {
        // docs/mcp/tools/sequence/kmer_entropy.md Example 1: 2-mers AT 3, TG 2, GC 2, CA 2;
        // scipy.stats.entropy([3, 2, 2, 2], base=2) = 1.974937501201927 (the old doc value 3.17 was wrong).
        Assert.That(SequenceTools.KmerEntropy("ATGCATGCAT", 2).Entropy, Is.EqualTo(1.974937501201927).Within(1e-12));
        Assert.That(SequenceTools.KmerEntropy("AAAAAAAAAA", 2).Entropy, Is.EqualTo(0.0));
    }
}
