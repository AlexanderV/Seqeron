using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class CpGObservedExpectedTests
{
    [Test]
    public void CpGObservedExpected_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.CpGObservedExpected("CGCG"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.CpGObservedExpected(""));
        Assert.Throws<ArgumentException>(() => AnnotationTools.CpGObservedExpected(null!));
    }

    [Test]
    public void CpGObservedExpected_Binding_InvokesSuccessfully()
    {
        // EpigeneticsAnalyzer.CalculateCpGObservedExpected: observed CpG / ((C * G) / length).
        // CGCGCG: length 6, C=3, G=3, observed CpG=3 -> expected=(3*3)/6=1.5 -> ratio=3/1.5=2.0.
        var result = AnnotationTools.CpGObservedExpected("CGCGCG");
        Assert.That(result.Ratio, Is.EqualTo(2.0).Within(1e-9));
    }

    [Test]
    public void CpGObservedExpected_NoCpG_ReturnsZero()
    {
        // No C or G (expected=0) -> the guarded branch returns 0.
        Assert.That(AnnotationTools.CpGObservedExpected("ATAT").Ratio, Is.EqualTo(0));
        // C and G present but no CG dinucleotide (GC only) -> observed=0 -> ratio 0.
        Assert.That(AnnotationTools.CpGObservedExpected("ATGCATGC").Ratio, Is.EqualTo(0).Within(1e-9));
    }
}
