using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class PredictGeneStructureTests
{
    private const string Exon1 = "AUGCCCAAAGGGCCCUUUAAAGGGCCCUUUAAAGC"; // 35
    private const string Donor = "GUAAGU";
    private const string IntronBody = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"; // 60
    private const string Ppt = "UUUUUUUUUUUUUU"; // 14
    private const string Acceptor = "CAG";
    private const string Exon2 = "GCCUUUAAAGGGCCCUUUAAAGGGCCCUUUAAAGC"; // 35
    private const string TwoExon = Exon1 + Donor + IntronBody + Ppt + Acceptor + Exon2;
    private static readonly int IntronPartLength = (Donor + IntronBody + Ppt + Acceptor).Length; // 83

    [Test]
    public void PredictGeneStructure_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.PredictGeneStructure(TwoExon, minExonLength: 5, minScore: 0.2));
        Assert.Throws<ArgumentException>(() => AnnotationTools.PredictGeneStructure(""));
        Assert.Throws<ArgumentException>(() => AnnotationTools.PredictGeneStructure(null!));
    }

    [Test]
    public void PredictGeneStructure_Binding_InvokesSuccessfully()
    {
        // Mirrors Seqeron.Genomics.Tests: one GT-AG intron -> 1 intron + 2 exons; intron starts after
        // exon1 (index 35) with length 83; exon+intron lengths cover the whole sequence.
        var result = AnnotationTools.PredictGeneStructure(TwoExon, minExonLength: 5, minScore: 0.2);

        Assert.Multiple(() =>
        {
            Assert.That(result.Introns, Has.Count.EqualTo(1));
            Assert.That(result.Exons, Has.Count.EqualTo(2));
            Assert.That(result.Introns[0].Start, Is.EqualTo(Exon1.Length));
            Assert.That(result.Introns[0].Length, Is.EqualTo(IntronPartLength));
            Assert.That(result.Exons[0].Type, Is.EqualTo("Initial"));
            Assert.That(result.Exons[1].Type, Is.EqualTo("Terminal"));

            int covered = result.Exons.Sum(e => e.Length) + result.Introns.Sum(i => i.Length);
            Assert.That(covered, Is.EqualTo(TwoExon.Length));
        });
    }

    [Test]
    public void PredictGeneStructure_SingleExon_NoIntrons()
    {
        // A 50 nt sequence with no GU dinucleotide is a single exon with no introns.
        const string singleExon = "AACCAACCAACCAACCAACCAACCAACCAACCAACCAACCAACCAACCAA";
        var result = AnnotationTools.PredictGeneStructure(singleExon, minExonLength: 5, minScore: 0.2);
        Assert.Multiple(() =>
        {
            Assert.That(result.Introns, Is.Empty);
            Assert.That(result.Exons, Has.Count.EqualTo(1));
            Assert.That(result.Exons[0].Type, Is.EqualTo("Single"));
        });
    }
}
