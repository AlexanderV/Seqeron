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
        Assert.DoesNotThrow(() => AlignmentTools.MergeContigs("AAA", "GGG", 0));
        // An overlap outside [0, shorter contig length] is not a suffix/prefix overlap.
        Assert.Throws<ArgumentOutOfRangeException>(() => AlignmentTools.MergeContigs("AAA", "GGG", -1));
        Assert.Throws<ArgumentOutOfRangeException>(() => AlignmentTools.MergeContigs("AAA", "GGGG", 4));
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

    [Test]
    public void MergeContigs_ResolveOverlapConflicts_UsesEmbossMergerRule()
    {
        // Overlap column G (contig1, near its end) vs T (contig2, interior): EMBOSS merger window-5
        // quality 14/11 vs 18/11 -> T. Default keeps contig1's G.
        Assert.Multiple(() =>
        {
            Assert.That(AlignmentTools.MergeContigs("AAAACCCCGGGG", "CCCCGTGGTTTT", 8).Merged,
                Is.EqualTo("AAAACCCCGGGGTTTT"));
            Assert.That(AlignmentTools.MergeContigs("AAAACCCCGGGG", "CCCCGTGGTTTT", 8, resolveOverlapConflicts: true).Merged,
                Is.EqualTo("AAAACCCCGTGGTTTT"));
        });
    }
}
