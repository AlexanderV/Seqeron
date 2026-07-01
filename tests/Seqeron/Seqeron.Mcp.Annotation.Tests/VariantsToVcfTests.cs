using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class VariantsToVcfTests
{
    [Test]
    public void VariantsToVcf_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.VariantsToVcf(new[] { new VariantDto(2, "G", "T", "SNP", 2) }));
        Assert.Throws<ArgumentNullException>(() => AnnotationTools.VariantsToVcf(null!));
    }

    [Test]
    public void VariantsToVcf_EmitsHeaderAndRecord()
    {
        // SNP G>T at 0-based position 2 -> VCF POS 3.
        var result = AnnotationTools.VariantsToVcf(new[] { new VariantDto(2, "G", "T", "SNP", 2) }, "chr1", "SAMPLE");

        Assert.Multiple(() =>
        {
            Assert.That(result.Lines[0], Is.EqualTo("##fileformat=VCFv4.2"));
            Assert.That(result.Lines[1], Is.EqualTo("##source=Seqeron.Genomics"));
            Assert.That(result.Lines[2], Is.EqualTo("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE"));
            Assert.That(result.Lines[3], Is.EqualTo("chr1\t3\t.\tG\tT\t.\tPASS\t.\tGT\t0/1"));
        });
    }

    [Test]
    public void VariantsToVcf_GapAllelesBecomeDot()
    {
        // Insertion: reference gap '-' becomes '.' in REF; POS = position + 1.
        var result = AnnotationTools.VariantsToVcf(new[] { new VariantDto(3, "-", "T", "Insertion", 3) }, "chrX", "S1");
        Assert.That(result.Lines[3], Is.EqualTo("chrX\t4\t.\t.\tT\t.\tPASS\t.\tGT\t0/1"));
    }

    [Test]
    public void VariantsToVcf_CustomSampleName_AppearsInHeader()
    {
        var result = AnnotationTools.VariantsToVcf(Array.Empty<VariantDto>(), "chr1", "TUMOR");
        Assert.That(result.Lines[2], Does.EndWith("\tTUMOR"));
        Assert.That(result.Lines, Has.Count.EqualTo(3)); // header only, no records
    }
}
