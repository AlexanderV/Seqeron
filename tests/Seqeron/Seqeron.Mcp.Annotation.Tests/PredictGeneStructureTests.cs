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

    // Ensembl ENST00000335295 (HBB), gene strand, GRCh38 11:5225464-5227071:-1; CDS starts at index 50.
    private const string Hbb =
        "ACATTTGCTTCTGACACAACTGTGTTCACTAGCAACCTCAAACAGACACCATGGTGCATCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAGGTTGGTATCAAGGTTACAAGACAGGTTTAAGGAGACCAATAGAAACTGGGCATGTGGAGACAGAGAAGACTCTTGGGTTTCTGATAGGCACTGACTCTCTCTGCCTATTGGTCTATTTTCCCACCCTTAGGCTGCTGGTGGTCTACCCTTGGACCCAGAGGTTCTTTGAGTCCTTTGGGGATCTGTCCACTCCTGATGCTGTTATGGGCAACCCTAAGGTGAAGGCTCATGGCAAGAAAGTGCTCGGTGCCTTTAGTGATGGCCTGGCTCACCTGGACAACCTCAAGGGCACCTTTGCCACACTGAGTGAGCTGCACTGTGACAAGCTGCACGTGGATCCTGAGAACTTCAGGGTGAGTCTATGGGACGCTTGATGTTTTCTTTCCCCTTCTTTTCTATGGTTAAGTTCATGTCATAGGAAGGGGATAAGTAACAGGGTACAGTTTAGAATGGGAAACAGACGAATGATTGCATCAGTGTGGAAGTCTCAGGATCGTTTTAGTTTCTTTTATTTGCTGTTCATAACAATTGTTTTCTTTTGTTTAATTCTTGCTTTCTTTTTTTTTCTTCTCCGCAATTTTTACTATTATACTTAATGCCTTAACATTGTGTATAACAAAAGGAAATATCTCTGAGATACATTAAGTAACTTAAAAAAAAACTTTACACAGTCTGCCTAGTACATTACTATTTGGAATATATGTGTGCTTATTTGCATATTCATAATCTCCCTACTTTATTTTCTTTTATTTTTAATTGATACATAATCATTATACATATTTATGGGTTAAAGTGTAATGTTTTAATATGTGTACACATATTGACCAAATCAGGGTAATTTTGCATTTGTAATTTTAAAAAATGCTTTCTTCTTTTAATATACTTTTTTGTTTATCTTATTTCTAATACTTTCCCTAATCTCTTTCTTTCAGGGCAATAATGATACAATGTATCATGCCTCTTTGCACCATTCTAAAGAATAACAGTGATAATTTCTGGGTTAAGGCAATAGCAATATCTCTGCATATAAATATTTCTGCATATAAATTGTAACTGATGTAAGAGGTTTCATATTGCTAATAGCAGCTACAATCCAGCTACCATTCTGCTTTTATTTTATGGTTGGGATAAGGCTGGATTATTCTGAGTCCAAGCTAGGCCCTTTTGCTAATCATGTTCATACCTCTTATCTTCCTCCCACAGCTCCTGGGCAACGTGCTGGTCTGTGTGCTGGCCCATCACTTTGGCAAAGAATTCACCCCACCAGTGCAGGCTGCCTATCAGAAAGTGGTGGCTGGTGTGGCTAATGCCCTGGCCCACAAGTATCACTAAGCTCGCTTTCTTGCTGTCCAATTTCTATTAAAGGTTCCTTTGTTCCCTAAGTCCAACTACTAAACTGGGGGATATTATGAAGGGCCTTGAGCATCTGGATTCTGCCTAATAAAAAACATTTATTTTCATTGCAA";

    [Test]
    public void PredictGeneStructure_SpliceDp_RecoversHbbIntronsAndOrf()
    {
        var r = AnnotationTools.PredictGeneStructure(Hbb, model: "splice_dp", startPosition: 50);

        Assert.Multiple(() =>
        {
            Assert.That(r.Introns.Select(i => (i.Start, i.End)), Is.EqualTo(new[] { (142, 271), (495, 1344) }));
            Assert.That(r.Exons, Has.Count.EqualTo(3));
            Assert.That(r.CdsStart, Is.EqualTo(50));
            Assert.That(r.CodingSequence, Has.Length.EqualTo(444));
            Assert.That(r.CodingSequence, Does.StartWith("AUGGUGCAUCUGACUCCUGAGGAG"));
            Assert.That(r.IntronScore, Is.GreaterThan(0));
            Assert.That(r.OverallScore, Is.EqualTo(r.IntronScore).Within(1e-9));
            Assert.That(r.Exons.Sum(e => e.Length) + r.Introns.Sum(i => i.Length), Is.EqualTo(Hbb.Length));
        });
    }

    [Test]
    public void PredictGeneStructure_SpliceDp_DefaultStartIsFirstAug()
    {
        var r = AnnotationTools.PredictGeneStructure(Hbb, model: "splice_dp");
        Assert.That(r.CdsStart, Is.EqualTo(50));
    }

    [Test]
    public void PredictGeneStructure_SpliceDp_NoOpenReadingFrame_Throws() =>
        Assert.Throws<ArgumentException>(() => AnnotationTools.PredictGeneStructure("CCAUGGCCGCCGCCGCC", model: "splice_dp"));

    [Test]
    public void PredictGeneStructure_UnknownModel_Throws() =>
        Assert.Throws<ArgumentException>(() => AnnotationTools.PredictGeneStructure(TwoExon, model: "viterbi"));

    [Test]
    public void PredictGeneStructure_GreedyDefault_HasNoCdsFields()
    {
        var r = AnnotationTools.PredictGeneStructure(TwoExon, minExonLength: 5, minScore: 0.2);
        Assert.Multiple(() =>
        {
            Assert.That(r.CdsStart, Is.Null);
            Assert.That(r.CodingSequence, Is.Null);
        });
    }
}
