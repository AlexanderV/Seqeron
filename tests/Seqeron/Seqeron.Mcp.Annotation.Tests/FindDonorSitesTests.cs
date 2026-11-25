using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class FindDonorSitesTests
{
    [Test]
    public void FindDonorSites_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.FindDonorSites("CAGGUAAGU", minScore: 0.3));
        Assert.Throws<ArgumentException>(() => AnnotationTools.FindDonorSites(""));
        Assert.Throws<ArgumentException>(() => AnnotationTools.FindDonorSites(null!));
    }

    [Test]
    public void FindDonorSites_Binding_InvokesSuccessfully()
    {
        // Mirrors Seqeron.Genomics.Tests SpliceSitePredictor_DonorSite_Tests: CAG|GUAAGU carries the
        // highest Shapiro & Senapathy (1987) primate percentage at every scored position -2..+6
        // -> one donor site at position 3, score/confidence 1.0.
        var result = AnnotationTools.FindDonorSites("CAGGUAAGU", minScore: 0.3);

        Assert.That(result.Sites, Has.Count.EqualTo(1));
        Assert.Multiple(() =>
        {
            Assert.That(result.Sites[0].Type, Is.EqualTo("Donor"));
            Assert.That(result.Sites[0].Position, Is.EqualTo(3));
            Assert.That(result.Sites[0].Score, Is.EqualTo(1.0).Within(1e-10));
            Assert.That(result.Sites[0].Confidence, Is.EqualTo(1.0).Within(1e-10));
        });
    }

    [Test]
    public void FindDonorSites_NoGU_ReturnsEmpty()
    {
        // No GU/GT dinucleotide -> no donor sites.
        var result = AnnotationTools.FindDonorSites("AAAAAAAAA", minScore: 0.3);
        Assert.That(result.Sites, Is.Empty);
    }

    [Test]
    public void FindDonorSites_ShapiroSenapathy_PublishedWorkedExample()
    {
        // Shapiro & Senapathy (1987) Methods: primate AAGTGAGT -> 100(509-47)/(595-47) = 84.3.
        var result = AnnotationTools.FindDonorSites("AAGTGAGT", minScore: 0.0);
        Assert.That(result.Sites.Single().Score, Is.EqualTo(462.0 / 548.0).Within(1e-15));

        // Table 5: invertebrate (C. elegans myosin intron 3) AC/GTAAGT -> 87.8.
        var inv = AnnotationTools.FindDonorSites("ACGTAAGT", minScore: 0.0, category: "invertebrate");
        Assert.That(Math.Round(inv.Sites.Single().Score * 100, 1), Is.EqualTo(87.8));
    }

    [Test]
    public void FindDonorSites_MaxEntModel_DelegatesToScore5()
    {
        // maxentpy score5('CAGGTAAGT') = 10.858313101356437 bits.
        var result = AnnotationTools.FindDonorSites("CAGGTAAGT", minScore: 0.0, model: "maxent");
        Assert.That(result.Sites.Single().Score, Is.EqualTo(10.858313101356437).Within(1e-12));
    }

    [Test]
    public void FindDonorSites_U12_GTAG_ReportedWhenNonCanonical()
    {
        // /GTATCCTT matches /RTATCCTY at 8/8; the U2 S&S score (376/548 = 0.686) is below 0.9.
        var result = AnnotationTools.FindDonorSites("CAGGTATCCTTA", minScore: 0.9, includeNonCanonical: true);
        Assert.That(result.Sites.Select(s => s.Type), Is.EquivalentTo(new[] { "U12Donor" }));
    }

    [Test]
    public void FindDonorSites_UnknownModelOrCategory_Throws()
    {
        Assert.Throws<ArgumentException>(() => AnnotationTools.FindDonorSites("CAGGTAAGT", model: "pwm"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.FindDonorSites("CAGGTAAGT", category: "fungi"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.FindDonorSites("CAGGTAAGT", category: "99"));
    }
}
