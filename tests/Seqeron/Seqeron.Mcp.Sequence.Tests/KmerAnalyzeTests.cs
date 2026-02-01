using NUnit.Framework;
using Seqeron.Mcp.Sequence.Tools;

namespace Seqeron.Mcp.Sequence.Tests;

[TestFixture]
public class KmerAnalyzeTests
{
    [Test]
    public void KmerAnalyze_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => SequenceTools.KmerAnalyze("ATGCGATCGATCG", 3));
        Assert.Throws<ArgumentException>(() => SequenceTools.KmerAnalyze("", 3));
        Assert.Throws<ArgumentException>(() => SequenceTools.KmerAnalyze(null!, 3));
        Assert.Throws<ArgumentException>(() => SequenceTools.KmerAnalyze("ATGC", 0)); // k < 1
    }

    [Test]
    public void KmerAnalyze_Binding_InvokesSuccessfully()
    {
        var result = SequenceTools.KmerAnalyze("ATGCATGCATGC", 3);
        Assert.That(result.K, Is.EqualTo(3));
        Assert.That(result.TotalKmers, Is.GreaterThan(0));
        Assert.That(result.UniqueKmers, Is.GreaterThan(0));
        Assert.That(result.MaxCount, Is.GreaterThanOrEqualTo(result.MinCount));
        Assert.That(result.AverageCount, Is.GreaterThan(0));
        Assert.That(result.Entropy, Is.GreaterThanOrEqualTo(0));
    }
}
