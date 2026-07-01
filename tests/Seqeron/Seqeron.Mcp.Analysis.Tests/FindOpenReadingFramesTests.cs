using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>find_open_reading_frames</c> MCP tool.
/// Expected values from the ATG-to-first-in-frame-stop definition (length divisible
/// by 3, frame 1-3, codonCount = length/3), NOT the wrapper output.
/// </summary>
[TestFixture]
public class FindOpenReadingFramesTests
{
    [Test]
    public void FindOpenReadingFrames_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.FindOpenReadingFrames("ATGTAA", 6));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindOpenReadingFrames("", 6));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindOpenReadingFrames(null!, 6));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindOpenReadingFrames("XYZ", 6));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.FindOpenReadingFrames("ATGTAA", 0));
    }

    [Test]
    public void FindOpenReadingFrames_Binding_InvokesSuccessfully()
    {
        // Minimal forward ORF ATG-TAA (RC = TTACAT has no ATG).
        var minimal = AnalysisTools.FindOpenReadingFrames("ATGTAA", 6).Items;
        Assert.Multiple(() =>
        {
            Assert.That(minimal, Has.Length.EqualTo(1));
            Assert.That(minimal[0].Sequence, Is.EqualTo("ATGTAA"));
            Assert.That(minimal[0].Position, Is.EqualTo(0));
            Assert.That(minimal[0].Frame, Is.EqualTo(1));
            Assert.That(minimal[0].IsReverseComplement, Is.False);
            Assert.That(minimal[0].Length, Is.EqualTo(6));
            Assert.That(minimal[0].CodonCount, Is.EqualTo(2));
        });

        // ATG-AAA-TAG (RC = CTATTTCAT has no ATG).
        var extra = AnalysisTools.FindOpenReadingFrames("ATGAAATAG", 9).Items;
        Assert.Multiple(() =>
        {
            Assert.That(extra, Has.Length.EqualTo(1));
            Assert.That(extra[0].Sequence, Is.EqualTo("ATGAAATAG"));
            Assert.That(extra[0].Length, Is.EqualTo(9));
            Assert.That(extra[0].CodonCount, Is.EqualTo(3));
        });
    }
}
