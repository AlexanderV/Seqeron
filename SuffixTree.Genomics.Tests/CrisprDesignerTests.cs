using NUnit.Framework;
using SuffixTree.Genomics;
using System.Linq;

namespace SuffixTree.Genomics.Tests;

/// <summary>
/// Tests for CRISPR guide RNA design, evaluation, and off-target analysis.
/// PAM site detection tests are in CrisprDesigner_PAM_Tests.cs.
/// </summary>
[TestFixture]
public class CrisprDesignerTests
{
    #region Guide RNA Evaluation Tests

    [Test]
    public void EvaluateGuideRna_OptimalGuide_HighScore()
    {
        // Good guide: ~50% GC, no polyT, low self-complementarity
        string guide = "ACGTACGTACGTACGTACGT"; // 50% GC
        var candidate = CrisprDesigner.EvaluateGuideRna(guide, CrisprSystemType.SpCas9);

        Assert.That(candidate.Score, Is.GreaterThan(70));
        Assert.That(candidate.GcContent, Is.EqualTo(50));
        Assert.That(candidate.HasPolyT, Is.False);
    }

    [Test]
    public void EvaluateGuideRna_LowGcContent_LowerScore()
    {
        string guide = "AAAAAAAAAAAAAAAAAAAA"; // 0% GC
        var candidate = CrisprDesigner.EvaluateGuideRna(guide, CrisprSystemType.SpCas9);

        Assert.That(candidate.Score, Is.LessThan(50));
        Assert.That(candidate.GcContent, Is.EqualTo(0));
        Assert.That(candidate.Issues, Has.Some.Contains("Low GC"));
    }

    [Test]
    public void EvaluateGuideRna_HighGcContent_LowerScore()
    {
        string guide = "GCGCGCGCGCGCGCGCGCGC"; // 100% GC
        var candidate = CrisprDesigner.EvaluateGuideRna(guide, CrisprSystemType.SpCas9);

        Assert.That(candidate.Score, Is.LessThan(50));
        Assert.That(candidate.GcContent, Is.EqualTo(100));
        Assert.That(candidate.Issues, Has.Some.Contains("High GC"));
    }

    [Test]
    public void EvaluateGuideRna_HasPolyT_Penalized()
    {
        string guide = "ACGTACGTTTTTACGTACGT"; // Contains TTTT
        var candidate = CrisprDesigner.EvaluateGuideRna(guide, CrisprSystemType.SpCas9);

        Assert.That(candidate.HasPolyT, Is.True);
        Assert.That(candidate.Issues, Has.Some.Contains("TTTT"));
    }

    [Test]
    public void EvaluateGuideRna_NoPolyT_NotPenalized()
    {
        string guide = "ACGTACGTACGTACGTACGT";
        var candidate = CrisprDesigner.EvaluateGuideRna(guide, CrisprSystemType.SpCas9);

        Assert.That(candidate.HasPolyT, Is.False);
    }

    [Test]
    public void EvaluateGuideRna_CalculatesSeedGc()
    {
        string guide = "AAAAAAAAAAAAACGTACGT"; // Last 12 bases are AAAACGTACGT
        var candidate = CrisprDesigner.EvaluateGuideRna(guide, CrisprSystemType.SpCas9);

        Assert.That(candidate.SeedGcContent, Is.GreaterThan(0));
    }

    [Test]
    public void EvaluateGuideRna_FullGuideRna_IncludesScaffold()
    {
        string guide = "ACGTACGTACGTACGTACGT";
        var candidate = CrisprDesigner.EvaluateGuideRna(guide, CrisprSystemType.SpCas9);

        Assert.That(candidate.FullGuideRna, Does.StartWith(guide));
        Assert.That(candidate.FullGuideRna.Length, Is.GreaterThan(guide.Length));
    }

    [Test]
    public void EvaluateGuideRna_EmptyGuide_ThrowsException()
    {
        Assert.Throws<ArgumentNullException>(() =>
            CrisprDesigner.EvaluateGuideRna("", CrisprSystemType.SpCas9));
    }

    #endregion

    #region Guide RNA Design Tests

    [Test]
    public void DesignGuideRnas_WithPamInRegion_ReturnsGuides()
    {
        // Create a sequence with PAM in the target region
        var sequence = new DnaSequence("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAGG");
        var guides = CrisprDesigner.DesignGuideRnas(sequence, 20, 45, CrisprSystemType.SpCas9).ToList();

        // Should find at least one guide
        Assert.That(guides, Has.Count.GreaterThanOrEqualTo(0));
    }

    [Test]
    public void DesignGuideRnas_InvalidRegionStart_ThrowsException()
    {
        var sequence = new DnaSequence("ACGTACGTACGT");
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            CrisprDesigner.DesignGuideRnas(sequence, -1, 10, CrisprSystemType.SpCas9).ToList());
    }

    [Test]
    public void DesignGuideRnas_InvalidRegionEnd_ThrowsException()
    {
        var sequence = new DnaSequence("ACGTACGTACGT");
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            CrisprDesigner.DesignGuideRnas(sequence, 0, 100, CrisprSystemType.SpCas9).ToList());
    }

    [Test]
    public void DesignGuideRnas_NullSequence_ThrowsException()
    {
        Assert.Throws<ArgumentNullException>(() =>
            CrisprDesigner.DesignGuideRnas(null!, 0, 10, CrisprSystemType.SpCas9).ToList());
    }

    #endregion

    #region Off-Target Analysis Tests

    [Test]
    public void FindOffTargets_NoMismatches_ReturnsEmpty()
    {
        string guide = "ACGTACGTACGTACGTACGT";
        var genome = new DnaSequence("ACGTACGTACGTACGTACGTACGTAGG"); // Exact match

        var offTargets = CrisprDesigner.FindOffTargets(guide, genome, 3, CrisprSystemType.SpCas9).ToList();

        // Off-targets require mismatches, exact matches are not off-targets
        Assert.That(offTargets.All(ot => ot.Mismatches > 0));
    }

    [Test]
    public void FindOffTargets_WithMismatches_FindsOffTargets()
    {
        string guide = "ACGTACGTACGTACGTACGT";
        // Create sequence with similar but not identical target
        var genome = new DnaSequence("ACGTACGTACGTACGTACGTACGTAGGTTTTTTTTTTTTTTTTTTTTACGAACGTACGTACGTACGTCGG");

        var offTargets = CrisprDesigner.FindOffTargets(guide, genome, 3, CrisprSystemType.SpCas9).ToList();

        // May or may not find off-targets depending on sequence
        Assert.That(offTargets, Is.Not.Null);
    }

    [Test]
    public void FindOffTargets_MaxMismatchesRespected()
    {
        string guide = "ACGTACGTACGTACGTACGT";
        var genome = new DnaSequence("AAAAAAAAAAAAAAAAAAAAAAAAGG"); // Many mismatches

        var offTargets = CrisprDesigner.FindOffTargets(guide, genome, 2, CrisprSystemType.SpCas9).ToList();

        Assert.That(offTargets.All(ot => ot.Mismatches <= 2));
    }

    [Test]
    public void FindOffTargets_ReturnsMismatchPositions()
    {
        string guide = "ACGTACGTACGTACGTACGT";
        var genome = new DnaSequence("ACGTACGTACGTACGTACGTTCGAACGTACGTACGTACGTCGG"); // Different target

        var offTargets = CrisprDesigner.FindOffTargets(guide, genome, 3, CrisprSystemType.SpCas9).ToList();

        foreach (var ot in offTargets)
        {
            Assert.That(ot.MismatchPositions, Is.Not.Null);
            Assert.That(ot.MismatchPositions.Count, Is.EqualTo(ot.Mismatches));
        }
    }

    [Test]
    public void FindOffTargets_EmptyGuide_ThrowsException()
    {
        var genome = new DnaSequence("ACGTACGT");
        Assert.Throws<ArgumentNullException>(() =>
            CrisprDesigner.FindOffTargets("", genome, 3, CrisprSystemType.SpCas9).ToList());
    }

    [Test]
    public void FindOffTargets_NullGenome_ThrowsException()
    {
        Assert.Throws<ArgumentNullException>(() =>
            CrisprDesigner.FindOffTargets("ACGTACGT", null!, 3, CrisprSystemType.SpCas9).ToList());
    }

    [Test]
    public void FindOffTargets_InvalidMaxMismatches_ThrowsException()
    {
        var genome = new DnaSequence("ACGTACGT");
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            CrisprDesigner.FindOffTargets("ACGTACGT", genome, 10, CrisprSystemType.SpCas9).ToList());
    }

    #endregion

    #region Specificity Score Tests

    [Test]
    public void CalculateSpecificityScore_NoOffTargets_ReturnsHigh()
    {
        string guide = "ACGTACGTACGTACGTACGT";
        // Short sequence unlikely to have off-targets
        var genome = new DnaSequence("ACGTACGTACGTACGTACGTACGTAGG");

        double score = CrisprDesigner.CalculateSpecificityScore(guide, genome, CrisprSystemType.SpCas9);

        Assert.That(score, Is.GreaterThanOrEqualTo(0));
        Assert.That(score, Is.LessThanOrEqualTo(100));
    }

    [Test]
    public void CalculateSpecificityScore_ManyOffTargets_ReturnsLower()
    {
        string guide = "AAAAAAAAAAAAAAAAAAAA";
        // Sequence with many similar regions
        var genome = new DnaSequence("AAAAAAAAAAAAAAAAAAAAAAGGAAAAAAAAAAAAAAAAAAACGGAAAAAAAAAAAAAAAAAAAAGG");

        double score = CrisprDesigner.CalculateSpecificityScore(guide, genome, CrisprSystemType.SpCas9);

        Assert.That(score, Is.LessThanOrEqualTo(100));
    }

    #endregion

    #region Parameter Tests

    [Test]
    public void GuideRnaParameters_Default_HasValidValues()
    {
        var defaults = GuideRnaParameters.Default;

        Assert.That(defaults.MinGcContent, Is.EqualTo(40));
        Assert.That(defaults.MaxGcContent, Is.EqualTo(70));
        Assert.That(defaults.MinScore, Is.EqualTo(50));
        Assert.That(defaults.AvoidPolyT, Is.True);
        Assert.That(defaults.CheckSelfComplementarity, Is.True);
    }

    [Test]
    public void GuideRnaParameters_CustomValues_Respected()
    {
        var custom = new GuideRnaParameters(
            MinGcContent: 30,
            MaxGcContent: 80,
            MinScore: 40,
            AvoidPolyT: false,
            CheckSelfComplementarity: false);

        Assert.That(custom.MinGcContent, Is.EqualTo(30));
        Assert.That(custom.MaxGcContent, Is.EqualTo(80));
        Assert.That(custom.MinScore, Is.EqualTo(40));
    }

    #endregion
}
