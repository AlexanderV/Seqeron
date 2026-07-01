using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class ParseVcfVariantTests
{
    [Test]
    public void ParseVcfVariant_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.ParseVcfVariant("chr1", 100, ".", "A", "G"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.ParseVcfVariant("", 100, ".", "A", "G"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.ParseVcfVariant("chr1", 100, ".", "", "G"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.ParseVcfVariant("chr1", 100, ".", "A", ""));
    }

    [Test]
    public void ParseVcfVariant_Binding_InvokesSuccessfully()
    {
        // VariantAnnotator.ParseVcfVariant preserves fields verbatim and classifies via ClassifyVariant.
        var result = AnnotationTools.ParseVcfVariant("chr1", 100, "rs123", "A", "G", 42.0);

        Assert.Multiple(() =>
        {
            Assert.That(result.Variant.Chromosome, Is.EqualTo("chr1"));
            Assert.That(result.Variant.Position, Is.EqualTo(100));
            Assert.That(result.Variant.Id, Is.EqualTo("rs123"));
            Assert.That(result.Variant.Reference, Is.EqualTo("A"));
            Assert.That(result.Variant.Alternate, Is.EqualTo("G"));
            Assert.That(result.Variant.Quality, Is.EqualTo(42.0));
            Assert.That(result.Variant.Type, Is.EqualTo("SNV"));
        });
    }

    [Test]
    public void ParseVcfVariant_ClassifiesVariantType()
    {
        // ClassifyVariant rules (VariantAnnotator.cs#L196):
        //   1bp/1bp -> SNV; ref len1 & alt starts-with ref -> Insertion;
        //   alt len1 & ref starts-with alt -> Deletion; equal length >1 -> MNV.
        Assert.Multiple(() =>
        {
            Assert.That(AnnotationTools.ParseVcfVariant("chr1", 1, ".", "A", "G").Variant.Type,
                Is.EqualTo("SNV"));
            Assert.That(AnnotationTools.ParseVcfVariant("chr1", 1, ".", "A", "AT").Variant.Type,
                Is.EqualTo("Insertion"));
            Assert.That(AnnotationTools.ParseVcfVariant("chr1", 1, ".", "AT", "A").Variant.Type,
                Is.EqualTo("Deletion"));
            Assert.That(AnnotationTools.ParseVcfVariant("chr1", 1, ".", "AT", "GC").Variant.Type,
                Is.EqualTo("MNV"));
        });
    }

    [Test]
    public void ParseVcfVariant_NullQuality_IsPreserved()
    {
        var result = AnnotationTools.ParseVcfVariant("chr1", 100, ".", "A", "G");
        Assert.That(result.Variant.Quality, Is.Null);
    }
}
