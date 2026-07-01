using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class FindSplitReadsTests
{
    // 30S70M: 30 nt left soft-clip + 70 nt match. Sequence = 30 'A' + 70 'C' = 100 nt.
    private static readonly string Seq30S70M = new string('A', 30) + new string('C', 70);

    private static List<AlignmentInputDto> Alignments() =>
        new()
        {
            new AlignmentInputDto("split", "chr1", 1000, "30S70M", Seq30S70M),
            // 90M10S: right clip of 10 nt is below the default minClipLength (20) -> not a split read.
            new AlignmentInputDto("short-clip", "chr1", 2000, "90M10S", new string('C', 90) + new string('T', 10)),
        };

    [Test]
    public void FindSplitReads_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.FindSplitReads(Alignments()));
        Assert.Throws<ArgumentNullException>(() => AnnotationTools.FindSplitReads(null!));
    }

    [Test]
    public void FindSplitReads_Binding_InvokesSuccessfully()
    {
        // StructuralVariantAnalyzer.FindSplitReads: soft-clips >= minClipLength (20) become split reads.
        // 30S70M -> left clip 30 nt at PrimaryPosition 1000; clipped sequence = first 30 bases.
        var result = AnnotationTools.FindSplitReads(Alignments());

        Assert.That(result.Reads, Has.Count.EqualTo(1));
        Assert.Multiple(() =>
        {
            Assert.That(result.Reads[0].ReadId, Is.EqualTo("split"));
            Assert.That(result.Reads[0].PrimaryPosition, Is.EqualTo(1000));
            Assert.That(result.Reads[0].ClipLength, Is.EqualTo(30));
            Assert.That(result.Reads[0].ClippedSequence, Is.EqualTo(new string('A', 30)));
        });
    }

    [Test]
    public void FindSplitReads_NoClips_ReturnsEmpty()
    {
        // A fully-matched read (100M) has no soft-clips.
        var aligns = new List<AlignmentInputDto>
        {
            new AlignmentInputDto("aligned", "chr1", 500, "100M", new string('G', 100)),
        };
        var result = AnnotationTools.FindSplitReads(aligns);
        Assert.That(result.Reads, Is.Empty);
    }
}
