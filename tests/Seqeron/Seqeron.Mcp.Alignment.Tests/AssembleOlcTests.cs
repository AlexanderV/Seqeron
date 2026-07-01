using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;

namespace Seqeron.Mcp.Alignment.Tests;

[TestFixture]
public class AssembleOlcTests
{
    [Test]
    public void AssembleOlc_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AlignmentTools.AssembleOlc(new[] { "AAAAACCCCC", "CCCCCGGGGG" }, 5, 0.9, 4, 5));
        var empty = AlignmentTools.AssembleOlc(Array.Empty<string>(), 5, 0.9, 4, 5);
        Assert.That(empty.Contigs, Is.Empty);
        Assert.That(empty.TotalLength, Is.EqualTo(0));
    }

    [Test]
    public void AssembleOlc_Binding_InvokesSuccessfully()
    {
        // Two reads with a 5-bp suffix/prefix overlap (CCCCC) merge into one 15-bp contig
        // via overlap-layout-consensus: AAAAACCCCC + GGGGG.
        var r = AlignmentTools.AssembleOlc(new[] { "AAAAACCCCC", "CCCCCGGGGG" }, 5, 0.9, 4, 5);
        Assert.Multiple(() =>
        {
            Assert.That(r.Contigs, Has.Length.EqualTo(1));
            Assert.That(r.Contigs[0], Is.EqualTo("AAAAACCCCCGGGGG"));
            Assert.That(r.TotalReads, Is.EqualTo(2));
            Assert.That(r.LongestContig, Is.EqualTo(15));
            Assert.That(r.TotalLength, Is.EqualTo(15));
            Assert.That(r.N50, Is.EqualTo(15));
        });
    }
}
