using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class FindAcceptorSitesTests
{
    [Test]
    public void FindAcceptorSites_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.FindAcceptorSites("UUUUUUUUUUUUUUUUCAGGG", minScore: 0.1));
        Assert.Throws<ArgumentException>(() => AnnotationTools.FindAcceptorSites(""));
        Assert.Throws<ArgumentException>(() => AnnotationTools.FindAcceptorSites(null!));
    }

    [Test]
    public void FindAcceptorSites_Binding_InvokesSuccessfully()
    {
        // Mirrors SpliceSitePredictor_AcceptorSite_Tests M1: U×16 + CAG|G — Shapiro & Senapathy (1987)
        // primate 3' ss score: best 8 of 10 pyrimidine percentages and CAGG all maximal → 1.0.
        var result = AnnotationTools.FindAcceptorSites("UUUUUUUUUUUUUUUUCAGGG", minScore: 0.1);

        Assert.That(result.Sites, Has.Count.EqualTo(1));
        Assert.Multiple(() =>
        {
            Assert.That(result.Sites[0].Type, Is.EqualTo("Acceptor"));
            Assert.That(result.Sites[0].Position, Is.EqualTo(18));
            Assert.That(result.Sites[0].Score, Is.EqualTo(1.0));
        });
    }

    [Test]
    public void FindAcceptorSites_ShapiroSenapathy_PublishedTable5Example()
    {
        // S&S (1987) Table 5, human HLA intron 6: TTATTCTACTCCAG/G = 95.8 (primate).
        var site = AnnotationTools.FindAcceptorSites("TTATTCTACTCCAGG", minScore: 0.0).Sites.Single();

        Assert.Multiple(() =>
        {
            Assert.That(site.Position, Is.EqualTo(13));
            Assert.That(site.Motif, Is.EqualTo("UUAUUCUACUCCAGG"));
            Assert.That(site.Score, Is.EqualTo(0.9578125).Within(1e-15));
        });
    }

    [Test]
    public void FindAcceptorSites_MaxEntModel_DelegatesToScore3()
    {
        var site = AnnotationTools.FindAcceptorSites("TTCCAAACGAACTTTTGTAGGGA", minScore: 0.0, model: "maxent").Sites.Single();

        Assert.Multiple(() =>
        {
            Assert.That(site.Position, Is.EqualTo(19));
            Assert.That(site.Score, Is.EqualTo(2.8867730651152104));
        });
    }

    [Test]
    public void FindAcceptorSites_U12_BranchSiteAcceptor_ReportedWhenNonCanonical()
    {
        // U12 branch site UUCCUUAAC 8 nt upstream of the AC terminus → U12Acceptor 9/9.
        var sites = AnnotationTools.FindAcceptorSites("GUUCCUUAACGGGGGACGGG", minScore: 0.5, includeNonCanonical: true).Sites;

        Assert.That(sites.Single(s => s.Type == "U12Acceptor").Position, Is.EqualTo(16));
    }

    [Test]
    public void FindAcceptorSites_UnknownModelOrCategory_Throws()
    {
        Assert.Throws<ArgumentException>(() => AnnotationTools.FindAcceptorSites("UUUUUUUUUUUUUUUUCAGGG", model: "pwm"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.FindAcceptorSites("UUUUUUUUUUUUUUUUCAGGG", category: "fish"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.FindAcceptorSites("UUUUUUUUUUUUUUUUCAGGG", category: "3"));
    }

    [Test]
    public void FindAcceptorSites_NoAG_ReturnsEmpty()
    {
        // No AG dinucleotide (all pyrimidines) -> no acceptor sites.
        var result = AnnotationTools.FindAcceptorSites("CCCCCCCCCCCCCCCCCCCC", minScore: 0.1);
        Assert.That(result.Sites, Is.Empty);
    }
}
