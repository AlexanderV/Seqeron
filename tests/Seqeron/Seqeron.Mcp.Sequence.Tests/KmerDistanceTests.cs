using NUnit.Framework;
using Seqeron.Mcp.Sequence.Tools;

namespace Seqeron.Mcp.Sequence.Tests;

[TestFixture]
public class KmerDistanceTests
{
    [Test]
    public void KmerDistance_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => SequenceTools.KmerDistance("ATGCGATCG", "ATGCGATCG", 3));
        Assert.Throws<ArgumentException>(() => SequenceTools.KmerDistance("", "ATGC", 3));
        Assert.Throws<ArgumentException>(() => SequenceTools.KmerDistance("ATGC", "", 3));
        Assert.Throws<ArgumentException>(() => SequenceTools.KmerDistance(null!, "ATGC", 3));
        Assert.Throws<ArgumentException>(() => SequenceTools.KmerDistance("ATGC", null!, 3));
        Assert.Throws<ArgumentException>(() => SequenceTools.KmerDistance("ATGC", "ATGC", 0)); // k < 1
    }

    [Test]
    public void KmerDistance_Binding_InvokesSuccessfully()
    {
        // Identical sequences should have distance 0
        var identical = SequenceTools.KmerDistance("ATGCATGC", "ATGCATGC", 3);
        Assert.That(identical.Distance, Is.EqualTo(0).Within(0.0001));
        Assert.That(identical.K, Is.EqualTo(3));

        // Different sequences should have distance > 0
        var different = SequenceTools.KmerDistance("ATGCATGC", "CCCCCCCC", 3);
        Assert.That(different.Distance, Is.GreaterThan(0));
    }
}
