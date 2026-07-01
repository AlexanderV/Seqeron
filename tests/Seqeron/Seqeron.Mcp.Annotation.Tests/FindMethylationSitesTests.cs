using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class FindMethylationSitesTests
{
    [Test]
    public void FindMethylationSites_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.FindMethylationSites("CG"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.FindMethylationSites(""));
        Assert.Throws<ArgumentException>(() => AnnotationTools.FindMethylationSites(null!));
    }

    [Test]
    public void FindMethylationSites_Binding_InvokesSuccessfully()
    {
        // EpigeneticsAnalyzer.FindMethylationSites classifies each cytosine into CpG/CHG/CHH context
        // and reports 0-based position + up-to-3-base context window; level/coverage are 0 (no reads).
        // CAGCTTCG: idx0 CHG "CAG"; idx3 CHH "CTT"; idx6 CpG "CG" (2-base terminal window).
        var result = AnnotationTools.FindMethylationSites("CAGCTTCG");

        Assert.That(result.Sites, Has.Count.EqualTo(3));
        Assert.Multiple(() =>
        {
            Assert.That(result.Sites[0].Position, Is.EqualTo(0));
            Assert.That(result.Sites[0].Type, Is.EqualTo("CHG"));
            Assert.That(result.Sites[0].Context, Is.EqualTo("CAG"));

            Assert.That(result.Sites[1].Position, Is.EqualTo(3));
            Assert.That(result.Sites[1].Type, Is.EqualTo("CHH"));
            Assert.That(result.Sites[1].Context, Is.EqualTo("CTT"));

            Assert.That(result.Sites[2].Position, Is.EqualTo(6));
            Assert.That(result.Sites[2].Type, Is.EqualTo("CpG"));
            Assert.That(result.Sites[2].Context, Is.EqualTo("CG"));

            // Sequence-only input carries no bisulfite evidence.
            Assert.That(result.Sites[0].MethylationLevel, Is.EqualTo(0));
            Assert.That(result.Sites[0].Coverage, Is.EqualTo(0));
        });
    }

    [Test]
    public void FindMethylationSites_NoCytosineContext_ReturnsEmpty()
    {
        // No classifiable cytosine (all A/T) -> no sites.
        Assert.That(AnnotationTools.FindMethylationSites("ATATAT").Sites, Is.Empty);
    }
}
