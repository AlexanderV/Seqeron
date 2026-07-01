using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class AnnotateRegulatoryElementsTests
{
    private static AnnotatorVariantDto Snv() =>
        new("chr1", 100, "A", "G", "SNV", null, null);

    private static RegulatoryRegionInputDto Region(string chrom, int start, int end, string type) =>
        new(chrom, start, end, type, null, null, new List<string>());

    [Test]
    public void AnnotateRegulatoryElements_Schema_ValidatesCorrectly()
    {
        var regions = new List<RegulatoryRegionInputDto> { Region("chr1", 50, 150, "promoter") };
        Assert.DoesNotThrow(() => AnnotationTools.AnnotateRegulatoryElements(Snv(), regions));
        Assert.Throws<ArgumentNullException>(() => AnnotationTools.AnnotateRegulatoryElements(null!, regions));
        Assert.Throws<ArgumentNullException>(() => AnnotationTools.AnnotateRegulatoryElements(Snv(), null!));
    }

    [Test]
    public void AnnotateRegulatoryElements_Binding_InvokesSuccessfully()
    {
        // Variant at chr1:100 (varStart=varEnd=100). Overlap rule: varEnd>=Start && varStart<=End
        // and matching chromosome (VariantAnnotator.AnnotateRegulatoryElements).
        var regions = new List<RegulatoryRegionInputDto>
        {
            Region("chr1", 50, 150, "promoter"),  // overlaps -> included
            Region("chr1", 200, 300, "enhancer"), // no overlap -> excluded
            Region("chr2", 50, 150, "silencer"),  // wrong chromosome -> excluded
        };

        var result = AnnotationTools.AnnotateRegulatoryElements(Snv(), regions);

        Assert.That(result.Annotations, Has.Count.EqualTo(1));
        Assert.Multiple(() =>
        {
            Assert.That(result.Annotations[0].FeatureType, Is.EqualTo("promoter"));
            Assert.That(result.Annotations[0].Chromosome, Is.EqualTo("chr1"));
            Assert.That(result.Annotations[0].Start, Is.EqualTo(50));
            Assert.That(result.Annotations[0].End, Is.EqualTo(150));
        });
    }

    [Test]
    public void AnnotateRegulatoryElements_NoOverlap_ReturnsEmpty()
    {
        var regions = new List<RegulatoryRegionInputDto> { Region("chr1", 200, 300, "enhancer") };
        var result = AnnotationTools.AnnotateRegulatoryElements(Snv(), regions);
        Assert.That(result.Annotations, Is.Empty);
    }
}
