using NUnit.Framework;
using Seqeron.Genomics.Annotation;
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
        // Delegates to SpliceSitePredictor.FindRetainedIntronCandidates(seq, criteria, minScore).
        var criteria = new SpliceSitePredictor.RetainedIntronCriteria(
            RequireOpenReadingFrame: false, RequireLengthMultipleOf3: false);
        var expected = SpliceSitePredictor.FindRetainedIntronCandidates(TwoExon, criteria, 0.2).ToList();
        var result = AnnotationTools.FindRetainedIntronCandidates(TwoExon, minScore: 0.2,
            requireOpenReadingFrame: false, requireLengthMultipleOf3: false);

        Assert.That(result.Introns, Is.Not.Empty);
        Assert.That(result.Introns.Select(i => (i.Start, i.End)),
            Is.EqualTo(expected.Select(i => (i.Start, i.End))));
    }

    [Test]
    public void FindRetainedIntronCandidates_OptionalFilters_Apply()
    {
        var result = AnnotationTools.FindRetainedIntronCandidates(TwoExon, minScore: 0.2,
            maxIntronLength: 10, requireOpenReadingFrame: false, requireLengthMultipleOf3: false);
        Assert.That(result.Introns, Is.Empty);
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
