using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;

namespace Seqeron.Mcp.Alignment.Tests;

[TestFixture]
public class AssembleDeBruijnTests
{
    [Test]
    public void AssembleDeBruijn_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AlignmentTools.AssembleDeBruijn(new[] { "AAABBB", "AABBBC" }, 20, 0.9, 3, 3));
        // Empty read set yields an empty assembly (no throw).
        var empty = AlignmentTools.AssembleDeBruijn(Array.Empty<string>(), 20, 0.9, 3, 3);
        Assert.That(empty.Contigs, Is.Empty);
        Assert.That(empty.TotalLength, Is.EqualTo(0));
    }

    [Test]
    public void AssembleDeBruijn_Binding_InvokesSuccessfully()
    {
        // k=3 de Bruijn nodes are 2-mers; the Eulerian walk over reads AAABBB / AABBBC
        // spells the superstring "AAABBBBBBC" (10 bp). Contigs shorter than 3 bp are dropped.
        var r = AlignmentTools.AssembleDeBruijn(new[] { "AAABBB", "AABBBC" }, 20, 0.9, 3, 3);
        Assert.Multiple(() =>
        {
            Assert.That(r.Contigs, Has.Length.EqualTo(1));
            Assert.That(r.Contigs[0], Is.EqualTo("AAABBBBBBC"));
            Assert.That(r.TotalReads, Is.EqualTo(2));
            Assert.That(r.LongestContig, Is.EqualTo(10));
            Assert.That(r.TotalLength, Is.EqualTo(10));
            Assert.That(r.N50, Is.EqualTo(10));
        });
    }
}
