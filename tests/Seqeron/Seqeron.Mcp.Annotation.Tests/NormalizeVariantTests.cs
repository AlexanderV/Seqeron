using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class NormalizeVariantTests
{
    [Test]
    public void NormalizeVariant_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.NormalizeVariant("chr1", 100, "ACG", "ATG"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.NormalizeVariant("", 100, "ACG", "ATG"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.NormalizeVariant("chr1", 100, "", "ATG"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.NormalizeVariant("chr1", 100, "ACG", ""));
    }

    [Test]
    public void NormalizeVariant_CommonSuffix_Trimmed()
    {
        // Mirrors VariantAnnotatorTests.NormalizeVariant_CommonSuffix_Trimmed: ACG/ATG@100 -> C/T@101 SNV.
        var result = AnnotationTools.NormalizeVariant("chr1", 100, "ACG", "ATG");

        Assert.Multiple(() =>
        {
            Assert.That(result.Variant.Position, Is.EqualTo(101));
            Assert.That(result.Variant.Reference, Is.EqualTo("C"));
            Assert.That(result.Variant.Alternate, Is.EqualTo("T"));
            Assert.That(result.Variant.Type, Is.EqualTo("SNV"));
            Assert.That(result.Variant.Chromosome, Is.EqualTo("chr1"));
        });
    }

    [Test]
    public void NormalizeVariant_PreservesChromosome()
    {
        var result = AnnotationTools.NormalizeVariant("chrX", 100, "A", "T");
        Assert.Multiple(() =>
        {
            Assert.That(result.Variant.Chromosome, Is.EqualTo("chrX"));
            Assert.That(result.Variant.Position, Is.EqualTo(100));
            Assert.That(result.Variant.Type, Is.EqualTo("SNV"));
        });
    }
}
