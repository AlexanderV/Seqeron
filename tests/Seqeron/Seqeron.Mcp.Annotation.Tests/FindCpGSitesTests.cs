using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class FindCpGSitesTests
{
    [Test]
    public void FindCpGSites_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.FindCpGSites("ACGT"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.FindCpGSites(""));
        Assert.Throws<ArgumentException>(() => AnnotationTools.FindCpGSites(null!));
    }

    [Test]
    public void FindCpGSites_Binding_InvokesSuccessfully()
    {
        // EpigeneticsAnalyzer.FindCpGSites returns 0-based start index of each CG dinucleotide.
        // ACGTACGT: C at index 1 -> G at 2; C at index 5 -> G at 6.
        var result = AnnotationTools.FindCpGSites("ACGTACGT");
        Assert.That(result.Positions, Is.EqualTo(new[] { 1, 5 }));
    }

    [Test]
    public void FindCpGSites_CaseInsensitive_And_NoSites()
    {
        // Lower case is upper-cased before scanning.
        Assert.That(AnnotationTools.FindCpGSites("acg").Positions, Is.EqualTo(new[] { 1 }));
        // No CG dinucleotide -> empty.
        Assert.That(AnnotationTools.FindCpGSites("AAATTT").Positions, Is.Empty);
    }
}
