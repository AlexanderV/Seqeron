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
        Assert.Throws<ArgumentOutOfRangeException>(() => AlignmentTools.AssembleDeBruijn(new[] { "ACGT" }, 20, 0.9, 1, 1));
    }

    [Test]
    public void AssembleDeBruijn_Binding_InvokesSuccessfully()
    {
        // Overlapping reads share the 4-mers GGCG/GCGT, counted once: one non-branching path.
        var r = AlignmentTools.AssembleDeBruijn(new[] { "ATGGCGT", "GGCGTGCA" }, 20, 0.9, 4, 4);
        Assert.Multiple(() =>
        {
            Assert.That(r.Contigs, Is.EqualTo(new[] { "ATGGCGTGCA" }));
            Assert.That(r.TotalReads, Is.EqualTo(2));
            Assert.That(r.LongestContig, Is.EqualTo(10));
            Assert.That(r.TotalLength, Is.EqualTo(10));
            Assert.That(r.N50, Is.EqualTo(10));
        });
    }

    [Test]
    public void AssembleDeBruijn_BranchingGraph_EmitsUnitigs()
    {
        // k=3: AA and BB carry self-loops (branch nodes), so the contigs are the maximal
        // non-branching paths AAA, AABB, BBB, BBC (Compeau & Pevzner BA3M; Python oracle).
        var r = AlignmentTools.AssembleDeBruijn(new[] { "AAABBB", "AABBBC" }, 20, 0.9, 3, 3);
        Assert.Multiple(() =>
        {
            Assert.That(r.Contigs, Is.EqualTo(new[] { "AAA", "AABB", "BBB", "BBC" }));
            Assert.That(r.TotalLength, Is.EqualTo(13));
            Assert.That(r.LongestContig, Is.EqualTo(4));
        });
    }

    [Test]
    public void AssembleDeBruijn_ReverseComplement_JoinsStrands()
    {
        var r = AlignmentTools.AssembleDeBruijn(new[] { "GATTACAGGT", "TAAGGACCTG" }, 20, 0.9, 5, 1, includeReverseComplement: true);
        Assert.That(r.Contigs, Is.EqualTo(new[] { "GATTACAGGTCCTTA" }));
    }
}
