using NUnit.Framework;
using Seqeron.Genomics.Analysis;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class SiteAccessibilityTests
{
    // MiRnaAnalyzer.CalculateSiteAccessibility = P(site entirely unpaired) = Z_open/Z (Turner 2004
    // McCaskill) via the canonical RnaSecondaryStructure.CalculateRegionUnpairedProbability; for a
    // sequence ≤ 80 nt the fold context is the whole sequence.
    private const string AccSeq = "GAAAAUAAAC";

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
        double expected = RnaSecondaryStructure.CalculateRegionUnpairedProbability(AccSeq, windowEnd: 7, windowLength: 6);
        var result = AnnotationTools.SiteAccessibility(AccSeq, 2, 7);
        Assert.That(result.Accessibility, Is.EqualTo(expected));
    }

    [Test]
    public void SiteAccessibility_UnstructuredSite_MatchesViennaRna()
    {
        // ViennaRNA 2.x (dangles=2) Z_c/Z for site 3..10 = 0.970119.
        var result = AnnotationTools.SiteAccessibility("AUGCUACCUCAAAAAAAAAAAAAAAAA", 3, 10);
        Assert.That(result.Accessibility, Is.EqualTo(0.970119).Within(0.01));
    }
}
