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
}
