using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class FindIndelsTests
{
    [Test]
    public void FindIndels_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.FindIndels("ATGCAT", "ATGTCAT"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.FindIndels("", "ATGTCAT"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.FindIndels("ATGCAT", ""));
    }

    [Test]
    public void FindIndels_Insertion_ReturnedAsIndel()
    {
        // Insertion of T (ATGCAT -> ATGTCAT) is reported by FindIndels.
        var result = AnnotationTools.FindIndels("ATGCAT", "ATGTCAT");

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
    public void FindIndels_Deletion_ReturnedAsIndel()
    {
        var result = AnnotationTools.FindIndels("ATGTCAT", "ATGCAT");

        Assert.That(result.Variants, Has.Count.EqualTo(1));
        Assert.Multiple(() =>
        {
            Assert.That(result.Variants[0].Type, Is.EqualTo("Deletion"));
            Assert.That(result.Variants[0].ReferenceAllele, Is.EqualTo("T"));
            Assert.That(result.Variants[0].AlternateAllele, Is.EqualTo("-"));
            Assert.That(result.Variants[0].Position, Is.EqualTo(3));
        });
    }

    [Test]
    public void FindIndels_ExcludesSnps()
    {
        // Pure substitution yields no indels.
        var result = AnnotationTools.FindIndels("ATGCATGC", "ATGAATGC");
        Assert.That(result.Variants, Is.Empty);
    }
}
