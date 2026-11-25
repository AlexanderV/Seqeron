using NUnit.Framework;
using Seqeron.Genomics.Annotation;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class CalculateMfeIndexTests
{
    [Test]
    public void CalculateMfeIndex_Schema_ValidatesCorrectly()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => AnnotationTools.CalculateMfeIndex(-35.13, 0, 48.6));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnnotationTools.CalculateMfeIndex(-35.13, -1, 48.6));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnnotationTools.CalculateMfeIndex(double.NaN, 72, 48.6));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnnotationTools.CalculateMfeIndex(double.PositiveInfinity, 72, 48.6));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnnotationTools.CalculateMfeIndex(-35.13, 72, 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnnotationTools.CalculateMfeIndex(-35.13, 72, double.NaN));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnnotationTools.CalculateMfeIndex(-35.13, 72, 100.5));
    }

    [Test]
    public void CalculateMfeIndex_Binding_MatchesZhang2006Values()
    {
        // Library MF9: dG -48.48, 57 nt, GC 25/57 -> 1.939200; hsa-mir-21 CL1: -35.13, 72 nt, GC 35/72 -> 1.0037142857142858.
        Assert.That(AnnotationTools.CalculateMfeIndex(-48.48, 57, 100.0 * 25 / 57).Mfei, Is.EqualTo(1.939200).Within(1e-6));
        var r = AnnotationTools.CalculateMfeIndex(-35.13, 72, 100.0 * 35 / 72);
        Assert.That(r.Mfei, Is.EqualTo(1.0037142857142858).Within(1e-12));
        Assert.That(r.Mfei, Is.EqualTo(MiRnaAnalyzer.CalculateMfeIndex(-35.13, 72, 100.0 * 35 / 72)));
        // hand: GC 100% -> AMFE/100 = (100*30/60)/100 = 0.5 (fuzz INV-CL-3)
        Assert.That(AnnotationTools.CalculateMfeIndex(-30.0, 60, 100.0).Mfei, Is.EqualTo(0.5).Within(1e-15));
        // magnitude: sign of dG is irrelevant
        Assert.That(AnnotationTools.CalculateMfeIndex(30.0, 60, 100.0).Mfei, Is.EqualTo(0.5).Within(1e-15));
    }
}
