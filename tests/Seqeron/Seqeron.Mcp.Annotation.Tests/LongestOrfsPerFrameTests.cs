using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class LongestOrfsPerFrameTests
{
    [Test]
    public void LongestOrfsPerFrame_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.LongestOrfsPerFrame("ATGAAAAAAAAATAA"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.LongestOrfsPerFrame(""));
        Assert.Throws<ArgumentException>(() => AnnotationTools.LongestOrfsPerFrame(null!));
    }

    [Test]
    public void LongestOrfsPerFrame_ForwardOnly_ReturnsThreeFrames()
    {
        // Mirrors GenomeAnnotator_ORF_Tests.FindLongestOrfsPerFrame_ForwardOnly_Returns3Keys.
        var result = AnnotationTools.LongestOrfsPerFrame("ATGAAAAAAAAATAA", searchBothStrands: false);

        Assert.That(result.Frames, Has.Count.EqualTo(3));
        Assert.That(result.Frames.Select(f => f.Frame), Is.EqualTo(new[] { 1, 2, 3 }));

        var frame1 = result.Frames.Single(f => f.Frame == 1).Orf;
        Assert.Multiple(() =>
        {
            Assert.That(frame1, Is.Not.Null);
            Assert.That(frame1!.Start, Is.EqualTo(0));
            Assert.That(frame1.End, Is.EqualTo(15));
            Assert.That(frame1.ProteinSequence, Is.EqualTo("MKKK*"));
            Assert.That(frame1.IsReverseComplement, Is.False);
            // Frames with no ORF surface as a zero-value sentinel OrfDto (FirstOrDefault over a
            // value type): End == 0 and no protein, never a positive-length hit.
            var frame2 = result.Frames.Single(f => f.Frame == 2).Orf;
            var frame3 = result.Frames.Single(f => f.Frame == 3).Orf;
            Assert.That(frame2!.End, Is.EqualTo(0));
            Assert.That(frame2.ProteinSequence, Is.Null.Or.Empty);
            Assert.That(frame3!.End, Is.EqualTo(0));
            Assert.That(frame3.ProteinSequence, Is.Null.Or.Empty);
        });
    }

    [Test]
    public void LongestOrfsPerFrame_BothStrands_ReturnsSixFrames()
    {
        // Mirrors GenomeAnnotator_ORF_Tests.FindLongestOrfsPerFrame_BothStrands_Returns6Keys.
        var result = AnnotationTools.LongestOrfsPerFrame("ATGAAAAAAAAATAA", searchBothStrands: true);

        Assert.That(result.Frames, Has.Count.EqualTo(6));
        Assert.That(result.Frames.Select(f => f.Frame), Is.EqualTo(new[] { 1, 2, 3, -1, -2, -3 }));
        Assert.That(result.Frames.Single(f => f.Frame == 1).Orf!.ProteinSequence, Is.EqualTo("MKKK*"));
    }

    [Test]
    public void LongestOrfsPerFrame_ReturnsLongestWhenMultipleInFrame()
    {
        // Mirrors FindLongestOrfsPerFrame_ReturnsLongestPerFrame: short ORF + GGG + long ORF (both frame 1).
        // short = ATG + 10x AAA + TAA; long = ATG + 50x AAA + TAA.
        string shortOrf = "ATG" + string.Concat(Enumerable.Repeat("AAA", 10)) + "TAA";
        string longOrf = "ATG" + string.Concat(Enumerable.Repeat("AAA", 50)) + "TAA";
        var result = AnnotationTools.LongestOrfsPerFrame(shortOrf + "GGG" + longOrf, searchBothStrands: false);

        var frame1 = result.Frames.Single(f => f.Frame == 1).Orf;
        Assert.That(frame1, Is.Not.Null);
        Assert.That(frame1!.ProteinSequence.TrimEnd('*'), Has.Length.EqualTo(51));
    }
}
