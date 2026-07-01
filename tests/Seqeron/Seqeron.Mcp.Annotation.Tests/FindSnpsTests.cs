using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class FindSnpsTests
{
    [Test]
    public void FindSnps_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.FindSnps("ATGC", "ATTC"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.FindSnps("", "ATTC"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.FindSnps(null!, "ATTC"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.FindSnps("ATGC", ""));
        Assert.Throws<ArgumentException>(() => AnnotationTools.FindSnps("ATGC", null!));
    }

    [Test]
    public void FindSnps_SingleSubstitution_ReturnsSnpOnly()
    {
        // Mirrors VariantCaller_FindSnps_Tests.FindSnps_SubstitutionOnlyInput_ReturnsSnpsOnly.
        var result = AnnotationTools.FindSnps("ATGCATGC", "ATGAATGC");

        Assert.That(result.Variants, Has.Count.EqualTo(1));
        var v = result.Variants[0];
        Assert.Multiple(() =>
        {
            Assert.That(v.Type, Is.EqualTo("SNP"));
            Assert.That(v.Position, Is.EqualTo(3));
            Assert.That(v.ReferenceAllele, Is.EqualTo("C"));
            Assert.That(v.AlternateAllele, Is.EqualTo("A"));
        });
    }

    [Test]
    public void FindSnps_IdenticalSequences_ReturnsNoSnps()
    {
        var result = AnnotationTools.FindSnps("ATGCATGC", "ATGCATGC");
        Assert.That(result.Variants, Is.Empty);
    }

    [Test]
    public void FindSnps_ReturnsSnpsOnly_ExcludesIndels()
    {
        // Query has an inserted base; FindSnps must exclude the resulting indel.
        var result = AnnotationTools.FindSnps("ATGCATGC", "ATGCAATGC");
        Assert.That(result.Variants.All(v => v.Type == "SNP"), Is.True);
    }
}
