using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class PredictIntronsTests
{
    // Mirrors Seqeron.Genomics.Tests SpliceSitePredictor_GeneStructure_Tests two-exon design:
    // exon1(35) + GUAAGU donor + 60xA + 14xU PPT + CAG acceptor + exon2(35).
    private const string Exon1 = "AUGCCCAAAGGGCCCUUUAAAGGGCCCUUUAAAGC";
    private const string Donor = "GUAAGU";
    private const string IntronBody = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"; // 60 A
    private const string Ppt = "UUUUUUUUUUUUUU";
    private const string Acceptor = "CAG";
    private const string Exon2 = "GCCUUUAAAGGGCCCUUUAAAGGGCCCUUUAAAGC";
    private const string TwoExon = Exon1 + Donor + IntronBody + Ppt + Acceptor + Exon2;

    [Test]
    public void PredictIntrons_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.PredictIntrons(TwoExon, minIntronLength: 60, minScore: 0.2));
        Assert.Throws<ArgumentException>(() => AnnotationTools.PredictIntrons(""));
        Assert.Throws<ArgumentException>(() => AnnotationTools.PredictIntrons(null!));
    }

    [Test]
    public void PredictIntrons_Binding_InvokesSuccessfully()
    {
        // GU-AG intron detected; every returned intron respects the length window;
        // GU-donor introns are classified U2 (major spliceosome, Burge 1999).
        var result = AnnotationTools.PredictIntrons(TwoExon, minIntronLength: 60, minScore: 0.2);

        Assert.That(result.Introns, Is.Not.Empty);
        Assert.Multiple(() =>
        {
            Assert.That(result.Introns.All(i => i.Length >= 60), Is.True);
            var guIntrons = result.Introns.Where(i => i.Sequence.StartsWith("GU")).ToList();
            Assert.That(guIntrons, Is.Not.Empty);
            Assert.That(guIntrons.All(i => i.Type == "U2"), Is.True);
            Assert.That(guIntrons.All(i => i.Sequence.EndsWith("AG")), Is.True);
        });
    }

    [Test]
    public void PredictIntrons_MinLengthFiltersAll()
    {
        // Raising the minimum above all candidate lengths yields no introns.
        var result = AnnotationTools.PredictIntrons(TwoExon, minIntronLength: 200, minScore: 0.2);
        Assert.That(result.Introns, Is.Empty);
    }
}
