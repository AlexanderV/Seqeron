using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;

namespace Seqeron.Mcp.Alignment.Tests;

[TestFixture]
public class FindAllOverlapsTests
{
    [Test]
    public void FindAllOverlaps_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AlignmentTools.FindAllOverlaps(new[] { "AAAAA", "AAAAA" }, 5, 0.9));
        Assert.That(AlignmentTools.FindAllOverlaps(Array.Empty<string>(), 5, 0.9).Items, Is.Empty);
    }

    [Test]
    public void FindAllOverlaps_Binding_InvokesSuccessfully()
    {
        // read0 suffix CCCCC == read1 prefix CCCCC -> a single edge 0 -> 1, length 5,
        // starting at pos 5 in read0 and pos 0 in read1. No overlap the other way.
        var r = AlignmentTools.FindAllOverlaps(new[] { "AAAAACCCCC", "CCCCCGGGGG" }, 5, 0.9);
        Assert.Multiple(() =>
        {
            Assert.That(r.Items, Has.Length.EqualTo(1));
            Assert.That(r.Items[0].ReadIndex1, Is.EqualTo(0));
            Assert.That(r.Items[0].ReadIndex2, Is.EqualTo(1));
            Assert.That(r.Items[0].OverlapLength, Is.EqualTo(5));
            Assert.That(r.Items[0].Position1, Is.EqualTo(5));
            Assert.That(r.Items[0].Position2, Is.EqualTo(0));
        });
    }
}
