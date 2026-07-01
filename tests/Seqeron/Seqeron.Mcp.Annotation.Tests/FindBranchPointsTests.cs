using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class FindBranchPointsTests
{
    // Sequence with a UACUAAC-like branch consensus; DNA T is normalised to U internally.
    private const string Seq = "GGGGGUACUAACGGGGGUUUUUCAGG";

    [Test]
    public void FindBranchPoints_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.FindBranchPoints(Seq, minScore: 0.0));
        Assert.Throws<ArgumentException>(() => AnnotationTools.FindBranchPoints(""));
        Assert.Throws<ArgumentException>(() => AnnotationTools.FindBranchPoints(null!));
    }

    [Test]
    public void FindBranchPoints_Binding_InvokesSuccessfully()
    {
        // SpliceSitePredictor.FindBranchPoints: window i scored by the YNYURAC-like PWM;
        // reported Position = i + 5 (branch adenosine), Motif = the 7-mer at window start (i = Position-5),
        // Type = Branch. All reported sites must satisfy these invariants and score >= minScore.
        var result = AnnotationTools.FindBranchPoints(Seq, minScore: 0.0);
        var normalized = Seq.ToUpperInvariant().Replace('T', 'U');

        Assert.That(result.Sites, Is.Not.Empty);
        Assert.Multiple(() =>
        {
            foreach (var s in result.Sites)
            {
                Assert.That(s.Type, Is.EqualTo("Branch"));
                Assert.That(s.Score, Is.GreaterThanOrEqualTo(0.0));
                int windowStart = s.Position - 5;
                Assert.That(s.Motif, Is.EqualTo(normalized.Substring(windowStart, 7)));
            }
        });
    }

    [Test]
    public void FindBranchPoints_SearchRange_RestrictsResults()
    {
        // Restricting the search window cannot yield more sites than an unrestricted scan.
        var all = AnnotationTools.FindBranchPoints(Seq, 0, -1, 0.3);
        var restricted = AnnotationTools.FindBranchPoints(Seq, 0, 6, 0.3);
        Assert.That(restricted.Sites.Count, Is.LessThanOrEqualTo(all.Sites.Count));
    }

    [Test]
    public void FindBranchPoints_TooShort_ReturnsEmpty()
    {
        // Sequences shorter than 7 nt cannot contain a 7-mer branch window.
        var result = AnnotationTools.FindBranchPoints("AUCGA", minScore: 0.0);
        Assert.That(result.Sites, Is.Empty);
    }
}
