using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;

namespace Seqeron.Mcp.Alignment.Tests;

[TestFixture]
public class GlobalAlignTests
{
    [Test]
    public void GlobalAlign_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AlignmentTools.GlobalAlign("ACGT", "ACGT"));
        Assert.Throws<ArgumentException>(() => AlignmentTools.GlobalAlign("", "ACGT"));
        Assert.Throws<ArgumentException>(() => AlignmentTools.GlobalAlign("ACGT", null!));
    }

    [Test]
    public void GlobalAlign_Binding_InvokesSuccessfully()
    {
        // Needleman-Wunsch, default scoring (match +1). Identical 7-mers align end-to-end
        // with score = 7 (7 matches, no gaps).
        var same = AlignmentTools.GlobalAlign("GATTACA", "GATTACA");
        Assert.Multiple(() =>
        {
            Assert.That(same.Score, Is.EqualTo(7));
            Assert.That(same.AlignedSequence1, Is.EqualTo("GATTACA"));
            Assert.That(same.AlignedSequence2, Is.EqualTo("GATTACA"));
            Assert.That(same.AlignmentType, Is.EqualTo("Global"));
        });

        // One deletion: GATTACA vs GATACA. Best global alignment inserts one gap in seq2,
        // giving 6 matches and 1 gap: score = 6*1 + 1*(-1) = 5.
        var gap = AlignmentTools.GlobalAlign("GATTACA", "GATACA");
        Assert.Multiple(() =>
        {
            Assert.That(gap.Score, Is.EqualTo(5));
            Assert.That(gap.AlignedSequence1, Is.EqualTo("GATTACA"));
            Assert.That(gap.AlignedSequence2, Is.EqualTo("GA-TACA"));
        });
    }
}
