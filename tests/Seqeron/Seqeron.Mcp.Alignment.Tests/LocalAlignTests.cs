using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;

namespace Seqeron.Mcp.Alignment.Tests;

[TestFixture]
public class LocalAlignTests
{
    [Test]
    public void LocalAlign_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AlignmentTools.LocalAlign("ACGT", "ACGT"));
        Assert.Throws<ArgumentException>(() => AlignmentTools.LocalAlign("", "ACGT"));
        Assert.Throws<ArgumentException>(() => AlignmentTools.LocalAlign("ACGT", null!));
    }

    [Test]
    public void LocalAlign_Binding_InvokesSuccessfully()
    {
        // Smith-Waterman finds the best-scoring local substring. The shared core "ACGT"
        // scores 4 and is reported without the flanking, non-matching regions.
        var r = AlignmentTools.LocalAlign("TTTACGTTTT", "GGGACGTGGG");
        Assert.Multiple(() =>
        {
            Assert.That(r.Score, Is.EqualTo(4));
            Assert.That(r.AlignedSequence1, Is.EqualTo("ACGT"));
            Assert.That(r.AlignedSequence2, Is.EqualTo("ACGT"));
            Assert.That(r.AlignmentType, Is.EqualTo("Local"));
            Assert.That(r.StartPosition1, Is.EqualTo(3));
            Assert.That(r.EndPosition1, Is.EqualTo(6));
        });
    }
}
