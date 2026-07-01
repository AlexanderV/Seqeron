using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;

namespace Seqeron.Mcp.Alignment.Tests;

[TestFixture]
public class MergeContigsTests
{
    [Test]
    public void MergeContigs_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AlignmentTools.MergeContigs("AAAAACCCCC", "CCCCCGGGGG", 5));
        // A non-positive overlap falls back to plain concatenation (no throw).
        Assert.DoesNotThrow(() => AlignmentTools.MergeContigs("AAA", "GGG", 0));
    }

    [Test]
    public void MergeContigs_Binding_InvokesSuccessfully()
    {
        // Collapsing a 5-bp suffix/prefix overlap: |c1|+|c2|-overlap = 10+10-5 = 15.
        var merged = AlignmentTools.MergeContigs("AAAAACCCCC", "CCCCCGGGGG", 5);
        Assert.That(merged.Merged, Is.EqualTo("AAAAACCCCCGGGGG"));

        // overlapLength 0 -> concatenation.
        var concat = AlignmentTools.MergeContigs("AAA", "GGG", 0);
        Assert.That(concat.Merged, Is.EqualTo("AAAGGG"));
    }
}
