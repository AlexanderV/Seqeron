using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class FindRetainedIntronCandidatesTests
{
    private const string Exon1 = "AUGCCCAAAGGGCCCUUUAAAGGGCCCUUUAAAGC";
    private const string Donor = "GUAAGU";
    private const string IntronBody = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
    private const string Ppt = "UUUUUUUUUUUUUU";
    private const string Acceptor = "CAG";
    private const string Exon2 = "GCCUUUAAAGGGCCCUUUAAAGGGCCCUUUAAAGC";
    private const string TwoExon = Exon1 + Donor + IntronBody + Ppt + Acceptor + Exon2;

    [Test]
    public void FindRetainedIntronCandidates_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.FindRetainedIntronCandidates(TwoExon, minScore: 0.2));
        Assert.Throws<ArgumentException>(() => AnnotationTools.FindRetainedIntronCandidates(""));
        Assert.Throws<ArgumentException>(() => AnnotationTools.FindRetainedIntronCandidates(null!));
    }

    [Test]
    public void FindRetainedIntronCandidates_Binding_InvokesSuccessfully()
    {
        // SpliceSitePredictor.FindRetainedIntronCandidates = PredictIntrons(seq, 60, 500, minScore)
        // filtered to short (Length < 500) and moderate-score (Score < 0.8) introns.
        var result = AnnotationTools.FindRetainedIntronCandidates(TwoExon, minScore: 0.2);

        Assert.That(result.Introns, Is.Not.Empty);
        Assert.Multiple(() =>
        {
            Assert.That(result.Introns.All(i => i.Length < 500), Is.True);
            Assert.That(result.Introns.All(i => i.Score < 0.8), Is.True);
        });
    }

    [Test]
    public void FindRetainedIntronCandidates_NoIntron_ReturnsEmpty()
    {
        // No GU donor -> no introns -> no retained-intron candidates.
        var result = AnnotationTools.FindRetainedIntronCandidates(
            "AACCAACCAACCAACCAACCAACCAACCAACCAACCAACCAACCAACCAA", minScore: 0.2);
        Assert.That(result.Introns, Is.Empty);
    }
}
