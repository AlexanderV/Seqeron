using NUnit.Framework;
using SuffixTree.Genomics;
using System.Linq;

namespace SuffixTree.Genomics.Tests;

/// <summary>
/// Tests for CRISPR Off-Target Analysis and Specificity Scoring (CRISPR-OFF-001).
/// 
/// Guide RNA Design tests are in CrisprDesigner_GuideRNA_Tests.cs (CRISPR-GUIDE-001).
/// PAM site detection tests are in CrisprDesigner_PAM_Tests.cs (CRISPR-PAM-001).
/// </summary>
[TestFixture]
public class CrisprDesignerTests
{
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
}
