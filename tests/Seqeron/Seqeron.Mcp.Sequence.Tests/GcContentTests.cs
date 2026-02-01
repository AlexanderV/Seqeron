using NUnit.Framework;
using Seqeron.Mcp.Sequence.Tools;

namespace Seqeron.Mcp.Sequence.Tests;

[TestFixture]
public class GcContentTests
{
    [Test]
    public void GcContent_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => SequenceTools.GcContent("ATGCGATCGATCG"));
        Assert.Throws<ArgumentException>(() => SequenceTools.GcContent(""));
        Assert.Throws<ArgumentException>(() => SequenceTools.GcContent(null!));
    }

    [Test]
    public void GcContent_Binding_InvokesSuccessfully()
    {
        var result = SequenceTools.GcContent("ATGCGATCGATCG");
        Assert.That(result.GcContent, Is.InRange(0, 100));
        Assert.That(result.GcCount, Is.GreaterThanOrEqualTo(0));
        Assert.That(result.TotalCount, Is.EqualTo(13));

        // 50% GC content
        var fiftyPercent = SequenceTools.GcContent("ATGC");
        Assert.That(fiftyPercent.GcContent, Is.EqualTo(50));
        Assert.That(fiftyPercent.GcCount, Is.EqualTo(2));

        // 100% GC content
        var hundredPercent = SequenceTools.GcContent("GCGC");
        Assert.That(hundredPercent.GcContent, Is.EqualTo(100));
    }
}
