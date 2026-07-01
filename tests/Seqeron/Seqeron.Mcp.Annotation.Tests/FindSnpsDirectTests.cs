using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class FindSnpsDirectTests
{
    [Test]
    public void FindSnpsDirect_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.FindSnpsDirect("ATGC", "ATTC"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.FindSnpsDirect("", "ATTC"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.FindSnpsDirect(null!, "ATTC"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.FindSnpsDirect("ATGC", ""));
    }

    [Test]
    public void FindSnpsDirect_SingleSubstitution_ReturnsExactSnp()
    {
        // Mirrors VariantCaller_FindSnps_Tests.FindSnpsDirect_SingleSubstitution_ReturnsExactSnp.
        var result = AnnotationTools.FindSnpsDirect("ATGC", "ATTC");

        Assert.That(result.Variants, Has.Count.EqualTo(1));
        var v = result.Variants[0];
        Assert.Multiple(() =>
        {
            Assert.That(v.Type, Is.EqualTo("SNP"));
            Assert.That(v.Position, Is.EqualTo(2));
            Assert.That(v.ReferenceAllele, Is.EqualTo("G"));
            Assert.That(v.AlternateAllele, Is.EqualTo("T"));
            Assert.That(v.QueryPosition, Is.EqualTo(2));
        });
    }

    [Test]
    public void FindSnpsDirect_MultipleSubstitutions_ExactPositions()
    {
        // Mirrors FindSnpsDirect_MultipleSubstitutions_ReturnsSnpsAtExactPositions: AAAA vs TGTA.
        var result = AnnotationTools.FindSnpsDirect("AAAA", "TGTA");

        Assert.That(result.Variants, Has.Count.EqualTo(3));
        Assert.Multiple(() =>
        {
            Assert.That(result.Variants.Select(v => v.Position), Is.EqualTo(new[] { 0, 1, 2 }));
            Assert.That(result.Variants.Select(v => v.AlternateAllele), Is.EqualTo(new[] { "T", "G", "T" }));
        });
    }

    [Test]
    public void FindSnpsDirect_UnequalLengths_ComparesCommonPrefixOnly()
    {
        // Mirrors FindSnpsDirect_UnequalLengths_ComparesCommonPrefixOnly: ATGCAA vs ATTC -> 1 SNP at pos 2.
        var result = AnnotationTools.FindSnpsDirect("ATGCAA", "ATTC");

        Assert.That(result.Variants, Has.Count.EqualTo(1));
        Assert.Multiple(() =>
        {
            Assert.That(result.Variants[0].Position, Is.EqualTo(2));
            Assert.That(result.Variants[0].ReferenceAllele, Is.EqualTo("G"));
            Assert.That(result.Variants[0].AlternateAllele, Is.EqualTo("T"));
        });
    }
}
