using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class FindMicrohomologyTests
{
    [Test]
    public void FindMicrohomology_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.FindMicrohomology("AAACGT", "CGTTTT"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.FindMicrohomology("", "CGT"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.FindMicrohomology("CGT", ""));
        Assert.Throws<ArgumentException>(() => AnnotationTools.FindMicrohomology(null!, "CGT"));
    }

    [Test]
    public void FindMicrohomology_Binding_InvokesSuccessfully()
    {
        // StructuralVariantAnalyzer.FindMicrohomology: longest suffix of leftFlank equal to a prefix
        // of rightFlank. "AAACGT" / "CGTTTT" -> "CGT" (length 3).
        var result = AnnotationTools.FindMicrohomology("AAACGT", "CGTTTT");
        Assert.Multiple(() =>
        {
            Assert.That(result.MicrohomologyLength, Is.EqualTo(3));
            Assert.That(result.Sequence, Is.EqualTo("CGT"));
        });
    }

    [Test]
    public void FindMicrohomology_NoOverlap_ReturnsZero()
    {
        var result = AnnotationTools.FindMicrohomology("AAAA", "TTTT");
        Assert.Multiple(() =>
        {
            Assert.That(result.MicrohomologyLength, Is.EqualTo(0));
            Assert.That(result.Sequence, Is.EqualTo(""));
        });
    }

    [Test]
    public void FindMicrohomology_CaseInsensitive()
    {
        // Lower-case input is upper-cased before comparison.
        var result = AnnotationTools.FindMicrohomology("aaacgt", "cgtttt");
        Assert.That(result.Sequence, Is.EqualTo("CGT"));
    }
}
