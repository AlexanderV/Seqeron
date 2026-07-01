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
        // Mirrors Seqeron.Genomics.Tests SpliceSitePredictor_AcceptorSite_Tests:
        // strong polypyrimidine tract + CAG|G consensus -> one acceptor site at position 18, score ~0.8393.
        var result = AnnotationTools.FindAcceptorSites("UUUUUUUUUUUUUUUUCAGGG", minScore: 0.1);

        Assert.That(result.Sites, Has.Count.EqualTo(1));
        Assert.Multiple(() =>
        {
            Assert.That(result.Sites[0].Type, Is.EqualTo("Acceptor"));
            Assert.That(result.Sites[0].Position, Is.EqualTo(18));
            Assert.That(result.Sites[0].Score, Is.EqualTo(0.8393).Within(0.001));
        });
    }

    [Test]
    public void FindAcceptorSites_NoAG_ReturnsEmpty()
    {
        // No AG dinucleotide (all pyrimidines) -> no acceptor sites.
        var result = AnnotationTools.FindAcceptorSites("CCCCCCCCCCCCCCCCCCCC", minScore: 0.1);
        Assert.That(result.Sites, Is.Empty);
    }
}
