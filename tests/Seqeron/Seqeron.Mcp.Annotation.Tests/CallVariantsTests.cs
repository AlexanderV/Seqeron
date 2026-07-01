using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class CallVariantsTests
{
    [Test]
    public void CallVariants_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.CallVariants("ATGC", "ATTC"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.CallVariants("", "ATTC"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.CallVariants(null!, "ATTC"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.CallVariants("ATGC", ""));
        Assert.Throws<ArgumentException>(() => AnnotationTools.CallVariants("ATGC", null!));
    }

    [Test]
    public void CallVariants_SingleSnp_DetectsGtoT()
    {
        // Mirrors VariantCallerTests.CallVariants_SingleSnp_DetectsIt: ATGC -> ATTC, G>T at position 2.
        var result = AnnotationTools.CallVariants("ATGC", "ATTC");

        Assert.That(result.Variants, Has.Count.EqualTo(1));
        var v = result.Variants[0];
        Assert.Multiple(() =>
        {
            Assert.That(v.Type, Is.EqualTo("SNP"));
            Assert.That(v.ReferenceAllele, Is.EqualTo("G"));
            Assert.That(v.AlternateAllele, Is.EqualTo("T"));
            Assert.That(v.Position, Is.EqualTo(2));
        });
    }

    [Test]
    public void CallVariants_IdenticalSequences_NoVariants()
    {
        var result = AnnotationTools.CallVariants("ATGCATGC", "ATGCATGC");
        Assert.That(result.Variants, Is.Empty);
    }

    [Test]
    public void CallVariants_TwoSnps_DetectsBoth()
    {
        // Mirrors CallVariants_MultipleSnps_DetectsAll: AAAA -> TATA has two SNPs.
        var result = AnnotationTools.CallVariants("AAAA", "TATA");
        Assert.That(result.Variants.Count(v => v.Type == "SNP"), Is.EqualTo(2));
    }
}
