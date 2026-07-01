using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class IdentifyCnvsTests
{
    private static CopyNumberSegmentDto Seg(string chr, int start, int end, double logR, int cn) =>
        new(chr, start, end, logR, cn, 0.5, 10);

    private static List<CopyNumberSegmentDto> Segments() =>
        new()
        {
            Seg("chr1", 1000, 20000, -1.0, 1),   // CN 1, len 19000 -> Deletion
            Seg("chr1", 30000, 60000, 1.0, 4),   // CN 4, len 30000 -> Duplication
            Seg("chr1", 70000, 90000, 0.0, 2),   // CN 2 (baseline) -> skipped
            Seg("chr1", 95000, 96000, -1.0, 1),  // CN 1 but len 1000 < 10000 -> skipped
        };

    [Test]
    public void IdentifyCnvs_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.IdentifyCnvs(Segments()));
        Assert.Throws<ArgumentNullException>(() => AnnotationTools.IdentifyCnvs(null!));
    }

    [Test]
    public void IdentifyCnvs_Binding_InvokesSuccessfully()
    {
        // StructuralVariantAnalyzer.IdentifyCNVs: non-baseline segments >= minLength (10000).
        // CN < normal -> Deletion, CN > normal -> Duplication. Quality = |logRatio| * 50.
        var result = AnnotationTools.IdentifyCnvs(Segments());

        Assert.That(result.Variants, Has.Count.EqualTo(2));
        Assert.Multiple(() =>
        {
            Assert.That(result.Variants[0].Type, Is.EqualTo("Deletion"));
            Assert.That(result.Variants[0].Length, Is.EqualTo(19000));
            Assert.That(result.Variants[0].Quality, Is.EqualTo(50.0).Within(1e-9));
            Assert.That(result.Variants[1].Type, Is.EqualTo("Duplication"));
            Assert.That(result.Variants[1].Length, Is.EqualTo(30000));
        });
    }

    [Test]
    public void IdentifyCnvs_AllBaselineOrShort_ReturnsEmpty()
    {
        var segs = new List<CopyNumberSegmentDto>
        {
            Seg("chr1", 0, 50000, 0.0, 2),       // baseline
            Seg("chr1", 60000, 61000, -1.0, 1),  // short
        };
        Assert.That(AnnotationTools.IdentifyCnvs(segs).Variants, Is.Empty);
    }
}
