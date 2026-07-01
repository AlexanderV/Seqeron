using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class SiteAccessibilityTests
{
    // Mirrors Seqeron.Genomics.Tests MiRnaAnalyzerMutationTests.CalculateSiteAccessibility_KnownWindow.
    // GAAAAUAAAC (len 10): window covers whole sequence.
    // Watson-Crick non-wobble pairs with j>=i+4: (G0,C9), (A1,U5) -> structureScore=2.
    // maxPairs=(10*6)/2=30; accessibility = 1 - (2/30)*10 = 1 - 20/30.
    private const string AccSeq = "GAAAAUAAAC";
    private static readonly double AccExpected = 1.0 - 20.0 / 30.0;

    [Test]
    public void SiteAccessibility_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.SiteAccessibility(AccSeq, 2, 7));
        Assert.Throws<ArgumentException>(() => AnnotationTools.SiteAccessibility("", 0, 0));
        // Out-of-range site indices are rejected by the wrapper guard.
        Assert.Throws<ArgumentOutOfRangeException>(() => AnnotationTools.SiteAccessibility("AAAA", -1, 2));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnnotationTools.SiteAccessibility("AAAA", 0, 100));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnnotationTools.SiteAccessibility("AAAA", 3, 1));
    }

    [Test]
    public void SiteAccessibility_Binding_InvokesSuccessfully()
    {
        // MiRnaAnalyzer.CalculateSiteAccessibility: accessibility = max(0, 1 - density*10).
        var result = AnnotationTools.SiteAccessibility(AccSeq, 2, 7);
        Assert.That(result.Accessibility, Is.EqualTo(AccExpected).Within(1e-9));
    }

    [Test]
    public void SiteAccessibility_SiteStartZero_StillComputes()
    {
        // siteStart == 0 is valid (guard is strict siteStart < 0); same window -> same value.
        var result = AnnotationTools.SiteAccessibility(AccSeq, 0, 7);
        Assert.That(result.Accessibility, Is.EqualTo(AccExpected).Within(1e-9));
    }
}
