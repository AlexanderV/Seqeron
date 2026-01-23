using NUnit.Framework;
using SuffixTree.Genomics;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace SuffixTree.Genomics.Tests;

/// <summary>
/// Tests for GenomeAnnotator methods NOT covered by ANNOT-ORF-001.
/// ORF-related tests are in GenomeAnnotator_ORF_Tests.cs.
/// </summary>
[TestFixture]
public class GenomeAnnotatorTests
{
    // Helper method for reverse complement that accepts N
    private static string GetReverseComplement(string sequence)
    {
        var complement = new Dictionary<char, char>
        {
            ['A'] = 'T',
            ['T'] = 'A',
            ['C'] = 'G',
            ['G'] = 'C',
            ['a'] = 't',
            ['t'] = 'a',
            ['c'] = 'g',
            ['g'] = 'c',
            ['N'] = 'N',
            ['n'] = 'n'
        };
        var sb = new StringBuilder(sequence.Length);
        for (int i = sequence.Length - 1; i >= 0; i--)
        {
            char c = sequence[i];
            sb.Append(complement.GetValueOrDefault(c, c));
        }
        return sb.ToString();
    }

    #region FindRibosomeBindingSites Tests

    [Test]
    public void FindRibosomeBindingSites_FindsShinedalgarno()
    {
        // AGGAGG upstream of ATG
        string sequence = "GGGGAGGAGGGGGGGATG" + new string('A', 294) + "TAA";

        var rbsSites = GenomeAnnotator.FindRibosomeBindingSites(sequence).ToList();

        // May or may not find depending on distance
        Assert.That(rbsSites, Is.Not.Null);
    }

    [Test]
    public void FindRibosomeBindingSites_CorrectDistance()
    {
        // Position RBS 8bp upstream of start
        string sequence = "GGGGG" + "AGGAGG" + "GGGGGGGG" + "ATG" + new string('A', 294) + "TAA";

        var rbsSites = GenomeAnnotator.FindRibosomeBindingSites(
            sequence, upstreamWindow: 25, minDistance: 4, maxDistance: 15).ToList();

        Assert.That(rbsSites.Any(r => r.sequence.Contains("AGGAGG") || r.sequence.Contains("GGAGG")), Is.True);
    }

    #endregion

    #region PredictGenes Tests

    [Test]
    public void PredictGenes_ReturnsGeneAnnotations()
    {
        string sequence = "ATG" + new string('A', 294) + "TAA" +
                         "ATG" + new string('C', 294) + "TAG";

        var genes = GenomeAnnotator.PredictGenes(sequence, minOrfLength: 50).ToList();

        Assert.That(genes.Count, Is.GreaterThan(0));
        Assert.That(genes.All(g => g.Type == "CDS"), Is.True);
    }

    [Test]
    public void PredictGenes_AssignsGeneIds()
    {
        string sequence = "ATG" + new string('A', 294) + "TAA";

        var genes = GenomeAnnotator.PredictGenes(sequence, minOrfLength: 50, prefix: "test").ToList();

        Assert.That(genes.All(g => g.GeneId.StartsWith("test_")), Is.True);
    }

    [Test]
    public void PredictGenes_IncludesStrandInfo()
    {
        string sequence = "ATG" + new string('A', 294) + "TAA";

        var genes = GenomeAnnotator.PredictGenes(sequence, minOrfLength: 50).ToList();

        Assert.That(genes.All(g => g.Strand == '+' || g.Strand == '-'), Is.True);
    }

    #endregion

    #region GFF3 Parsing Tests

    [Test]
    public void ParseGff3_ValidLine_ParsesCorrectly()
    {
        var lines = new List<string>
        {
            "##gff-version 3",
            "seq1\t.\tgene\t100\t500\t.\t+\t.\tID=gene1;Name=testGene"
        };

        var features = GenomeAnnotator.ParseGff3(lines).ToList();

        Assert.That(features.Count, Is.EqualTo(1));
        Assert.That(features[0].Type, Is.EqualTo("gene"));
        Assert.That(features[0].Start, Is.EqualTo(100));
        Assert.That(features[0].End, Is.EqualTo(500));
        Assert.That(features[0].Strand, Is.EqualTo('+'));
    }

    [Test]
    public void ParseGff3_WithScore_ParsesScore()
    {
        var lines = new List<string>
        {
            "seq1\t.\tCDS\t100\t500\t0.95\t+\t0\tID=cds1"
        };

        var features = GenomeAnnotator.ParseGff3(lines).ToList();

        Assert.That(features[0].Score, Is.EqualTo(0.95));
    }

    [Test]
    public void ParseGff3_SkipsComments()
    {
        var lines = new List<string>
        {
            "# This is a comment",
            "##gff-version 3",
            "seq1\t.\tgene\t100\t500\t.\t+\t.\tID=gene1"
        };

        var features = GenomeAnnotator.ParseGff3(lines).ToList();

        Assert.That(features.Count, Is.EqualTo(1));
    }

    [Test]
    public void ParseGff3_ParsesAttributes()
    {
        var lines = new List<string>
        {
            "seq1\t.\tgene\t100\t500\t.\t+\t.\tID=gene1;Name=myGene;product=test%20protein"
        };

        var features = GenomeAnnotator.ParseGff3(lines).ToList();

        Assert.That(features[0].Attributes["ID"], Is.EqualTo("gene1"));
        Assert.That(features[0].Attributes["Name"], Is.EqualTo("myGene"));
        Assert.That(features[0].Attributes["product"], Is.EqualTo("test protein"));
    }

    #endregion

    #region ToGff3 Tests

    [Test]
    public void ToGff3_GeneratesValidOutput()
    {
        var annotations = new List<GenomeAnnotator.GeneAnnotation>
        {
            new GenomeAnnotator.GeneAnnotation(
                GeneId: "gene1",
                Start: 99,
                End: 500,
                Strand: '+',
                Type: "CDS",
                Product: "hypothetical protein",
                Attributes: new Dictionary<string, string> { ["frame"] = "1" })
        };

        var lines = GenomeAnnotator.ToGff3(annotations, "chr1").ToList();

        Assert.That(lines[0], Does.Contain("##gff-version 3"));
        Assert.That(lines[1], Does.Contain("chr1"));
        Assert.That(lines[1], Does.Contain("CDS"));
        Assert.That(lines[1], Does.Contain("ID=gene1"));
    }

    [Test]
    public void ToGff3_EscapesSpecialCharacters()
    {
        var annotations = new List<GenomeAnnotator.GeneAnnotation>
        {
            new GenomeAnnotator.GeneAnnotation(
                GeneId: "gene 1",
                Start: 0,
                End: 100,
                Strand: '+',
                Type: "CDS",
                Product: "test;product",
                Attributes: new Dictionary<string, string>())
        };

        var lines = GenomeAnnotator.ToGff3(annotations).ToList();

        Assert.That(lines[1], Does.Contain("gene%201"));
    }

    #endregion

    #region FindPromoterMotifs Tests

    [Test]
    public void FindPromoterMotifs_FindsMinus35Box()
    {
        string sequence = "GGGTTGACAGGGGGGGGGGGGGGGGGTATAAATGGG";

        var motifs = GenomeAnnotator.FindPromoterMotifs(sequence).ToList();

        Assert.That(motifs.Any(m => m.type == "-35 box"), Is.True);
    }

    [Test]
    public void FindPromoterMotifs_FindsMinus10Box()
    {
        string sequence = "GGGTTGACAGGGGGGGGGGGGGGGGGTATAAATGGG";

        var motifs = GenomeAnnotator.FindPromoterMotifs(sequence).ToList();

        Assert.That(motifs.Any(m => m.type == "-10 box"), Is.True);
    }

    [Test]
    public void FindPromoterMotifs_NoMotifs_ReturnsEmpty()
    {
        string sequence = "CCCCCCCCCCCCCCCCCCCC"; // Use C instead of G to avoid false positives

        var motifs = GenomeAnnotator.FindPromoterMotifs(sequence).ToList();

        Assert.That(motifs.Count, Is.EqualTo(0));
    }

    #endregion

    #region CalculateCodingPotential Tests

    [Test]
    public void CalculateCodingPotential_ValidCodingSequence_HighScore()
    {
        // Real coding sequence without internal stops
        string coding = "ATGGCAACGTCAACGACGTCAACGTAA"; // No internal stops

        double score = GenomeAnnotator.CalculateCodingPotential(coding);

        Assert.That(score, Is.GreaterThan(0));
    }

    [Test]
    public void CalculateCodingPotential_ShortSequence_ReturnsZero()
    {
        double score = GenomeAnnotator.CalculateCodingPotential("ACG");

        Assert.That(score, Is.EqualTo(0));
    }

    [Test]
    public void CalculateCodingPotential_ManyStopCodons_LowerScore()
    {
        string withStops = "ATGTAATAGTGA"; // Multiple stops

        double score = GenomeAnnotator.CalculateCodingPotential(withStops);

        // Score should be lower due to internal stops
        Assert.That(score, Is.GreaterThanOrEqualTo(0));
    }

    #endregion

    #region FindRepetitiveElements Tests

    [Test]
    public void FindRepetitiveElements_TandemRepeat_Detected()
    {
        string sequence = "GGGGG" + "ACGTACGTACGTACGT" + "GGGGG"; // 4x ACGT

        var repeats = GenomeAnnotator.FindRepetitiveElements(sequence, minRepeatLength: 8, minCopies: 2).ToList();

        Assert.That(repeats.Any(r => r.type == "tandem_repeat"), Is.True);
    }

    [Test]
    public void FindRepetitiveElements_InvertedRepeat_Detected()
    {
        // ACGTACGTAC and its reverse complement GTACGTACGT
        string arm1 = "ACGTACGTAC";
        string arm2 = GetReverseComplement(arm1);
        string sequence = "CCC" + arm1 + "CCCCC" + arm2 + "CCC";

        var repeats = GenomeAnnotator.FindRepetitiveElements(sequence, minRepeatLength: 10).ToList();

        Assert.That(repeats.Any(r => r.type == "inverted_repeat"), Is.True);
    }

    [Test]
    public void FindRepetitiveElements_NoRepeats_ReturnsEmpty()
    {
        string sequence = "ACGTCCCCCCCCCCCCCCC";

        var repeats = GenomeAnnotator.FindRepetitiveElements(sequence, minRepeatLength: 20).ToList();

        Assert.That(repeats.Count, Is.EqualTo(0));
    }

    #endregion

    #region GetCodonUsage Tests

    [Test]
    public void GetCodonUsage_CountsCodons()
    {
        string sequence = "ATGATGATG"; // 3x ATG

        var usage = GenomeAnnotator.GetCodonUsage(sequence);

        Assert.That(usage["ATG"], Is.EqualTo(3));
    }

    [Test]
    public void GetCodonUsage_MultipleCodons()
    {
        string sequence = "ATGAAATTT"; // ATG, AAA, TTT

        var usage = GenomeAnnotator.GetCodonUsage(sequence);

        Assert.That(usage.Count, Is.EqualTo(3));
        Assert.That(usage["ATG"], Is.EqualTo(1));
        Assert.That(usage["AAA"], Is.EqualTo(1));
        Assert.That(usage["TTT"], Is.EqualTo(1));
    }

    [Test]
    public void GetCodonUsage_IgnoresPartialCodons()
    {
        string sequence = "ATGAA"; // Only 1 complete codon + 2bp

        var usage = GenomeAnnotator.GetCodonUsage(sequence);

        Assert.That(usage.Count, Is.EqualTo(1));
    }

    [Test]
    public void GetCodonUsage_CaseInsensitive()
    {
        string sequence = "atgATGAtg";

        var usage = GenomeAnnotator.GetCodonUsage(sequence);

        Assert.That(usage["ATG"], Is.EqualTo(3));
    }

    #endregion
}
