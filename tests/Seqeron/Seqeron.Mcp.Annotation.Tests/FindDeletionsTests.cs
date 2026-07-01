using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class FindDeletionsTests
{
    [Test]
    public void FindDeletions_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.FindDeletions("ATGTCAT", "ATGCAT"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.FindDeletions("", "ATGCAT"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.FindDeletions("ATGTCAT", ""));
    }

    [Test]
    public void FindDeletions_SingleDeletedBase_ReturnsExactDeletion()
    {
        // Mirrors VariantCaller_FindIndels_Tests.FindDeletions_SingleDeletedBase_ReturnsExactDeletion.
        var result = AnnotationTools.FindDeletions("ATGTCAT", "ATGCAT");

        Assert.That(result.Variants, Has.Count.EqualTo(1));
        var v = result.Variants[0];
        Assert.Multiple(() =>
        {
            Assert.That(v.Type, Is.EqualTo("Deletion"));
            Assert.That(v.ReferenceAllele, Is.EqualTo("T"));
            Assert.That(v.AlternateAllele, Is.EqualTo("-"));
            Assert.That(v.Position, Is.EqualTo(3));
        });
    }

    [Test]
    public void FindDeletions_SubstitutionOnlyInput_ReturnsEmpty()
    {
        var result = AnnotationTools.FindDeletions("ATGCATGC", "ATGAATGC");
        Assert.That(result.Variants, Is.Empty);
    }

    [Test]
    public void FindDeletions_ReturnsDeletionsOnly()
    {
        // Insertion input yields no deletions.
        var result = AnnotationTools.FindDeletions("ATGCAT", "ATGTCAT");
        Assert.That(result.Variants, Is.Empty);
    }
}
