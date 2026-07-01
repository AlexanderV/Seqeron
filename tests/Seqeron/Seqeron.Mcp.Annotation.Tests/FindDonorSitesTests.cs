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
        // Mirrors Seqeron.Genomics.Tests SpliceSitePredictor_DonorSite_Tests:
        // perfect MAG|GURAGU consensus CAG|GUAAGU -> one donor site at position 3, score/confidence 1.0.
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
}
