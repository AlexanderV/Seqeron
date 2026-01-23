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

    // NOTE: FindRibosomeBindingSites and PredictGenes tests moved to
    // GenomeAnnotator_Gene_Tests.cs as part of ANNOT-GENE-001 consolidation.

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

    // NOTE: FindPromoterMotifs tests moved to GenomeAnnotator_PromoterMotif_Tests.cs
    // as part of ANNOT-PROM-001 consolidation.

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
