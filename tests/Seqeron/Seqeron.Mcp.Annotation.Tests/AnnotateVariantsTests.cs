using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class AnnotateVariantsTests
{
    [Test]
    public void AnnotateVariants_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.AnnotateVariants("ATGCAT", "ATTCAT"));

        Assert.Throws<ArgumentException>(() => AnnotationTools.AnnotateVariants("", "ATTCAT"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.AnnotateVariants(null!, "ATTCAT"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.AnnotateVariants("ATGCAT", ""));
        Assert.Throws<ArgumentException>(() => AnnotationTools.AnnotateVariants("ATGCAT", null!));
    }

    [Test]
    public void AnnotateVariants_CodingMissense_ReturnsMissenseTransversion()
    {
        // Reference codon ATG (Met) -> ATT (Ile) at CDS position 2 (G>T).
        // G(purine)->T(pyrimidine) is a transversion; M!=I and not stop -> Missense.
        var result = AnnotationTools.AnnotateVariants("ATGCAT", "ATTCAT", isCodingSequence: true);

        Assert.That(result.Annotated, Has.Count.EqualTo(1));
        var a = result.Annotated[0];
        Assert.Multiple(() =>
        {
            Assert.That(a.Variant.Position, Is.EqualTo(2));
            Assert.That(a.Variant.ReferenceAllele, Is.EqualTo("G"));
            Assert.That(a.Variant.AlternateAllele, Is.EqualTo("T"));
            Assert.That(a.Variant.Type, Is.EqualTo("SNP"));
            Assert.That(a.Effect, Is.EqualTo("Missense"));
            Assert.That(a.MutationType, Is.EqualTo("Transversion"));
        });
    }

    [Test]
    public void AnnotateVariants_NonCoding_EffectIsUnknown()
    {
        // Same SNP, but not treated as coding -> effect is Unknown, mutation type still classified.
        var result = AnnotationTools.AnnotateVariants("ATGCAT", "ATTCAT", isCodingSequence: false);

        Assert.That(result.Annotated, Has.Count.EqualTo(1));
        var a = result.Annotated[0];
        Assert.Multiple(() =>
        {
            Assert.That(a.Effect, Is.EqualTo("Unknown"));
            Assert.That(a.MutationType, Is.EqualTo("Transversion"));
        });
    }

    [Test]
    public void AnnotateVariants_IdenticalSequences_ReturnsNoVariants()
    {
        var result = AnnotationTools.AnnotateVariants("ATGCAT", "ATGCAT");
        Assert.That(result.Annotated, Is.Empty);
    }
}
