using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;

namespace Seqeron.Mcp.Alignment.Tests;

[TestFixture]
public class FindOverlapTests
{
    [Test]
    public void FindOverlap_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AlignmentTools.FindOverlap("AAAAACCCCC", "CCCCCGGGGG", 5, 0.9));
        // No qualifying overlap -> null result (not an exception).
        var none = AlignmentTools.FindOverlap("AAAAA", "TTTTT", 5, 0.9);
        Assert.That(none.Overlap, Is.Null);
    }

    [Test]
    public void FindOverlap_Binding_InvokesSuccessfully()
    {
        // Longest suffix(seq1)/prefix(seq2) match with identity >= 0.9: CCCCC, length 5,
        // starting at pos 5 in seq1 and pos 0 in seq2.
        var r = AlignmentTools.FindOverlap("AAAAACCCCC", "CCCCCGGGGG", 5, 0.9);
        Assert.Multiple(() =>
        {
            Assert.That(r.Overlap, Is.Not.Null);
            Assert.That(r.Overlap!.Length, Is.EqualTo(5));
            Assert.That(r.Overlap.Position1, Is.EqualTo(5));
            Assert.That(r.Overlap.Position2, Is.EqualTo(0));
        });
    }
}
