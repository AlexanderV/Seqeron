using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class FindBranchPointsTests
{
    // Sequence with a UACUAAC branch consensus; DNA T is normalised to U internally.
    private const string Seq = "GGGGGUACUAACGGGGGUUUUUCAGG";

    [Test]
    public void FindBranchPoints_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.FindBranchPoints(Seq, minScore: 0.0));
        Assert.Throws<ArgumentException>(() => AnnotationTools.FindBranchPoints(""));
        Assert.Throws<ArgumentException>(() => AnnotationTools.FindBranchPoints(null!));
        Assert.Throws<ArgumentException>(() => AnnotationTools.FindBranchPoints(Seq, model: "mouse"));
    }

    [Test]
    public void FindBranchPoints_Binding_InvokesSuccessfully()
    {
        // Human yUnAy (Gao et al. 2008): CUAAC → branch A at 10, score 1.0; UUCAG → A at 23,
        // (0.790 + 0.746 + 0.923) / 3.210.
        var result = AnnotationTools.FindBranchPoints(Seq, minScore: 0.75);

        Assert.That(result.Sites.Select(s => (s.Position, s.Motif)), Is.EqualTo(new[] { (10, "CUAAC"), (23, "UUCAG") }));
        Assert.Multiple(() =>
        {
            Assert.That(result.Sites[0].Score, Is.EqualTo(1.0));
            Assert.That(result.Sites[1].Score, Is.EqualTo((0.790 + 0.746 + 0.923) / 3.210).Within(1e-15));
            Assert.That(result.Sites.All(s => s.Type == "Branch"), Is.True);
        });
    }

    [Test]
    public void FindBranchPoints_YeastModel_UACUAAC()
    {
        var site = AnnotationTools.FindBranchPoints(Seq, minScore: 0.75, model: "yeast").Sites.Single();

        Assert.That((site.Position, site.Motif, site.Score), Is.EqualTo((10, "UACUAAC", 1.0)));
    }

    [Test]
    public void FindBranchPoints_SearchRange_RestrictsResults()
    {
        // Motif starts restricted to [0, 6]: only the CUAAC at start 7 is excluded → no 1.0 site.
        var all = AnnotationTools.FindBranchPoints(Seq, 0, -1, 0.3);
        var restricted = AnnotationTools.FindBranchPoints(Seq, 0, 6, 0.3);
        Assert.That(restricted.Sites.Count, Is.LessThan(all.Sites.Count));
        Assert.That(restricted.Sites.All(s => s.Position - 3 <= 6), Is.True);
    }

    [Test]
    public void FindBranchPoints_TooShort_ReturnsEmpty()
    {
        // Sequences shorter than the 5-nt yUnAy motif have no placement.
        var result = AnnotationTools.FindBranchPoints("AUCG", minScore: 0.0);
        Assert.That(result.Sites, Is.Empty);
    }
}
