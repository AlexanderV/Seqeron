using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;

namespace Seqeron.Mcp.Alignment.Tests;

[TestFixture]
public class ScaffoldContigsTests
{
    [Test]
    public void ScaffoldContigs_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AlignmentTools.ScaffoldContigs(
            new[] { "AAA", "GGG" }, new[] { new ContigLinkInput(0, 1, 3) }, "N"));
        // gapCharacter must be exactly one character.
        Assert.Throws<ArgumentException>(() => AlignmentTools.ScaffoldContigs(
            new[] { "AAA", "GGG" }, new[] { new ContigLinkInput(0, 1, 3) }, "NN"));
        Assert.Throws<ArgumentException>(() => AlignmentTools.ScaffoldContigs(
            new[] { "AAA", "GGG" }, new[] { new ContigLinkInput(0, 1, 3) }, ""));
    }

    [Test]
    public void ScaffoldContigs_Binding_InvokesSuccessfully()
    {
        // Link 0 -> 1 with gapSize 3 places contig[1] after contig[0] separated by 3 'N's:
        // "AAA" + "NNN" + "GGG".
        var r = AlignmentTools.ScaffoldContigs(
            new[] { "AAA", "GGG" }, new[] { new ContigLinkInput(0, 1, 3) }, "N");
        Assert.Multiple(() =>
        {
            Assert.That(r.Scaffolds, Has.Length.EqualTo(1));
            Assert.That(r.Scaffolds[0], Is.EqualTo("AAANNNGGG"));
        });
    }

    [Test]
    public void ScaffoldContigs_OrientationWeightAndOverlap_DelegateToGreedyScaffolder()
    {
        // Heavier 0→2 link wins over 0→1 (Huson et al. 2002 greedy); contig 2 reverse-complemented.
        var weighted = AlignmentTools.ScaffoldContigs(
            new[] { "AAAA", "CCCC", "GGTT" },
            new[] { new ContigLinkInput(0, 1, 2), new ContigLinkInput(0, 2, 3, Contig2Reversed: true, Weight: 4) });
        // A link into a lower-index contig is joined; a −3 agreeing overlap is collapsed.
        var overlap = AlignmentTools.ScaffoldContigs(
            new[] { "TACGGA", "ACGTAC" }, new[] { new ContigLinkInput(1, 0, -3) });
        Assert.Multiple(() =>
        {
            Assert.That(weighted.Scaffolds, Is.EqualTo(new[] { "AAAANNNAACC", "CCCC" }));
            Assert.That(overlap.Scaffolds, Is.EqualTo(new[] { "ACGTACGGA" }));
            Assert.Throws<ArgumentOutOfRangeException>(() => AlignmentTools.ScaffoldContigs(
                new[] { "AA", "CC" }, new[] { new ContigLinkInput(0, 1, 1, Weight: 0) }));
        });
    }
}
