using NUnit.Framework;
using Seqeron.Mcp.Sequence.Tools;

namespace Seqeron.Mcp.Sequence.Tests;

[TestFixture]
public class KmerCountTests
{
    [Test]
    public void KmerCount_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => SequenceTools.KmerCount("ATGCGATCGATCG", 3));
        Assert.Throws<ArgumentException>(() => SequenceTools.KmerCount("", 3));
        Assert.Throws<ArgumentException>(() => SequenceTools.KmerCount(null!, 3));
        Assert.Throws<ArgumentException>(() => SequenceTools.KmerCount("ATGC", 0)); // k < 1
    }

    [Test]
    public void KmerCount_Binding_InvokesSuccessfully()
    {
        var result = SequenceTools.KmerCount("ATGCATGC", 3);
        Assert.That(result.K, Is.EqualTo(3));
        Assert.That(result.Counts, Is.Not.Null);
        Assert.That(result.UniqueKmers, Is.GreaterThan(0));
        Assert.That(result.TotalKmers, Is.GreaterThan(0));

        // ATG appears twice
        Assert.That(result.Counts.ContainsKey("ATG"));
        Assert.That(result.Counts["ATG"], Is.EqualTo(2));
    }
}
