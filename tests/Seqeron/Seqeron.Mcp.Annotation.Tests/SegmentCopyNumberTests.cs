using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class SegmentCopyNumberTests
{
    private static List<CopyNumberProbeDto> Probes(double logRatio) =>
        Enumerable.Range(1, 5)
            .Select(i => new CopyNumberProbeDto("chr1", i * 100, logRatio, 0.5))
            .ToList();

    [Test]
    public void SegmentCopyNumber_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.SegmentCopyNumber(Probes(0.0)));
        Assert.Throws<ArgumentNullException>(() => AnnotationTools.SegmentCopyNumber(null!));
    }

    [Test]
    public void SegmentCopyNumber_Binding_InvokesSuccessfully()
    {
        // StructuralVariantAnalyzer.SegmentCopyNumber: >= minProbes (5) probes with a stable log-ratio
        // form one segment. CopyNumber = round(2 * 2^meanLogR). meanLogR 0 -> CN 2 (diploid).
        var result = AnnotationTools.SegmentCopyNumber(Probes(0.0));

        Assert.That(result.Segments, Has.Count.EqualTo(1));
        Assert.Multiple(() =>
        {
            Assert.That(result.Segments[0].LogRatio, Is.EqualTo(0.0).Within(1e-10));
            Assert.That(result.Segments[0].CopyNumber, Is.EqualTo(2));
            Assert.That(result.Segments[0].ProbeCount, Is.EqualTo(5));
            Assert.That(result.Segments[0].Start, Is.EqualTo(100));
            Assert.That(result.Segments[0].End, Is.EqualTo(500));
        });
    }

    [Test]
    public void SegmentCopyNumber_HalfDepth_CopyNumberOne()
    {
        // meanLogR -1 -> round(2 * 2^-1) = round(1) = 1 (single-copy deletion).
        var result = AnnotationTools.SegmentCopyNumber(Probes(-1.0));
        Assert.That(result.Segments, Has.Count.EqualTo(1));
        Assert.That(result.Segments[0].CopyNumber, Is.EqualTo(1));
    }

    [Test]
    public void SegmentCopyNumber_TooFewProbes_ReturnsEmpty()
    {
        // Fewer than minProbes (5) -> no segment emitted.
        var few = Enumerable.Range(1, 3)
            .Select(i => new CopyNumberProbeDto("chr1", i * 100, 0.0, 0.5)).ToList();
        Assert.That(AnnotationTools.SegmentCopyNumber(few).Segments, Is.Empty);
    }
}
