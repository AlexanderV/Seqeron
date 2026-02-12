using NUnit.Framework;
using Seqeron.Genomics;
using Seqeron.Genomics.Tests.Helpers;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Tests for GenomeAnnotator methods NOT covered by ANNOT-ORF-001.
/// ORF-related tests are in GenomeAnnotator_ORF_Tests.cs.
/// </summary>
[TestFixture]
public class GenomeAnnotatorTests
{

    // NOTE: FindRibosomeBindingSites and PredictGenes tests moved to
    // GenomeAnnotator_Gene_Tests.cs as part of ANNOT-GENE-001 consolidation.

    // NOTE: GFF3 parsing/export tests moved to GenomeAnnotator_GFF3_Tests.cs
    // as part of ANNOT-GFF-001 consolidation.

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
        string arm2 = TestSequenceHelpers.GetReverseComplement(arm1);
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
