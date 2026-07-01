using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class FindInsertionsTests
{
    [Test]
    public void FindInsertions_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.FindInsertions("ATGCAT", "ATGTCAT"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.FindInsertions("", "ATGTCAT"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.FindInsertions("ATGCAT", ""));
    }

    [Test]
    public void FindInsertions_SingleInsertedBase_ReturnsExactInsertion()
    {
        // Mirrors VariantCaller_FindIndels_Tests.FindInsertions_SingleInsertedBase_ReturnsExactInsertion.
        var result = AnnotationTools.FindInsertions("ATGCAT", "ATGTCAT");

        Assert.That(result.Variants, Has.Count.EqualTo(1));
        var v = result.Variants[0];
        Assert.Multiple(() =>
        {
            Assert.That(v.Type, Is.EqualTo("Insertion"));
            Assert.That(v.ReferenceAllele, Is.EqualTo("-"));
            Assert.That(v.AlternateAllele, Is.EqualTo("T"));
            Assert.That(v.Position, Is.EqualTo(3));
        });
    }

    [Test]
    public void FindInsertions_SubstitutionOnlyInput_ReturnsEmpty()
    {
        // Mirrors FindInsertions_SubstitutionOnlyInput_ReturnsEmpty.
        var result = AnnotationTools.FindInsertions("ATGCATGC", "ATGAATGC");
        Assert.That(result.Variants, Is.Empty);
    }

    [Test]
    public void FindInsertions_ReturnsInsertionsOnly()
    {
        var result = AnnotationTools.FindInsertions("ATGCATGC", "ATGTCATGG");
        Assert.That(result.Variants.All(v => v.Type == "Insertion"), Is.True);
    }
}
