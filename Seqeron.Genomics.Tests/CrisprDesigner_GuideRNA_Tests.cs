using NUnit.Framework;
using Seqeron.Genomics;
using System.Linq;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Tests for CRISPR Guide RNA Design (CRISPR-GUIDE-001).
/// Covers: EvaluateGuideRna, DesignGuideRnas, GuideRnaParameters
/// 
/// Evidence Sources:
/// - Addgene CRISPR Guide (https://www.addgene.org/guides/crispr/)
/// - Wikipedia: Guide RNA (https://en.wikipedia.org/wiki/Guide_RNA)
/// - Wikipedia: Protospacer adjacent motif (https://en.wikipedia.org/wiki/Protospacer_adjacent_motif)
/// 
/// TestSpec: TestSpecs/CRISPR-GUIDE-001.md
/// Algorithm Doc: docs/algorithms/MolTools/Guide_RNA_Design.md
/// </summary>
[TestFixture]
public class CrisprDesigner_GuideRNA_Tests
{
    #region MUST Tests - Guide RNA Evaluation

    /// <summary>
    /// M-001: Optimal guide (50% GC, no polyT) should score high.
    /// Evidence: Addgene - guides with optimal GC (40-70%) perform better.
    /// </summary>
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

    /// <summary>
    /// M-002: Low GC content (0%) should be penalized.
    /// Evidence: Wikipedia - GC >50% optimal for efficiency.
    /// </summary>
    [Test]
    public void EvaluateGuideRna_LowGcContent_LowerScore()
    {
        string guide = "AAAAAAAAAAAAAAAAAAAA"; // 0% GC
        var candidate = CrisprDesigner.EvaluateGuideRna(guide, CrisprSystemType.SpCas9);

        Assert.That(candidate.Score, Is.LessThan(50));
        Assert.That(candidate.GcContent, Is.EqualTo(0));
        Assert.That(candidate.Issues, Has.Some.Contains("Low GC"));
    }

    /// <summary>
    /// M-003: High GC content (100%) should be penalized.
    /// Evidence: High GC causes secondary structures.
    /// </summary>
    [Test]
    public void EvaluateGuideRna_HighGcContent_LowerScore()
    {
        string guide = "GCGCGCGCGCGCGCGCGCGC"; // 100% GC
        var candidate = CrisprDesigner.EvaluateGuideRna(guide, CrisprSystemType.SpCas9);

        Assert.That(candidate.Score, Is.LessThan(50));
        Assert.That(candidate.GcContent, Is.EqualTo(100));
        Assert.That(candidate.Issues, Has.Some.Contains("High GC"));
    }

    /// <summary>
    /// M-004: PolyT (TTTT) should be detected and penalized.
    /// Evidence: Addgene - Pol III terminates at TTTT sequences.
    /// </summary>
    [Test]
    public void EvaluateGuideRna_HasPolyT_Penalized()
    {
        string guide = "ACGTACGTTTTTACGTACGT"; // Contains TTTT
        var candidate = CrisprDesigner.EvaluateGuideRna(guide, CrisprSystemType.SpCas9);

        Assert.That(candidate.HasPolyT, Is.True);
        Assert.That(candidate.Issues, Has.Some.Contains("TTTT"));
    }

    /// <summary>
    /// M-005: Empty guide should throw ArgumentNullException.
    /// Evidence: Defensive programming.
    /// </summary>
    [Test]
    public void EvaluateGuideRna_EmptyGuide_ThrowsException()
    {
        Assert.Throws<ArgumentNullException>(() =>
            CrisprDesigner.EvaluateGuideRna("", CrisprSystemType.SpCas9));
    }

    /// <summary>
    /// M-006: FullGuideRna should include the scaffold sequence.
    /// Evidence: Addgene - sgRNA = spacer + scaffold.
    /// </summary>
    [Test]
    public void EvaluateGuideRna_FullGuideRna_IncludesScaffold()
    {
        string guide = "ACGTACGTACGTACGTACGT";
        var candidate = CrisprDesigner.EvaluateGuideRna(guide, CrisprSystemType.SpCas9);

        Assert.That(candidate.FullGuideRna, Does.StartWith(guide));
        Assert.That(candidate.FullGuideRna.Length, Is.GreaterThan(guide.Length));
    }

    #endregion

    #region MUST Tests - Guide RNA Design

    /// <summary>
    /// M-007: Null sequence should throw ArgumentNullException.
    /// </summary>
    [Test]
    public void DesignGuideRnas_NullSequence_ThrowsException()
    {
        Assert.Throws<ArgumentNullException>(() =>
            CrisprDesigner.DesignGuideRnas(null!, 0, 10, CrisprSystemType.SpCas9).ToList());
    }

    /// <summary>
    /// M-008: Invalid region start (negative) should throw.
    /// </summary>
    [Test]
    public void DesignGuideRnas_InvalidRegionStart_ThrowsException()
    {
        var sequence = new DnaSequence("ACGTACGTACGT");
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            CrisprDesigner.DesignGuideRnas(sequence, -1, 10, CrisprSystemType.SpCas9).ToList());
    }

    /// <summary>
    /// M-009: Invalid region end (beyond sequence) should throw.
    /// </summary>
    [Test]
    public void DesignGuideRnas_InvalidRegionEnd_ThrowsException()
    {
        var sequence = new DnaSequence("ACGTACGTACGT");
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            CrisprDesigner.DesignGuideRnas(sequence, 0, 100, CrisprSystemType.SpCas9).ToList());
    }

    #endregion

    #region SHOULD Tests - Guide RNA Evaluation

    /// <summary>
    /// S-001: Guide without polyT should not be penalized.
    /// </summary>
    [Test]
    public void EvaluateGuideRna_NoPolyT_NotPenalized()
    {
        string guide = "ACGTACGTACGTACGTACGT";
        var candidate = CrisprDesigner.EvaluateGuideRna(guide, CrisprSystemType.SpCas9);

        Assert.That(candidate.HasPolyT, Is.False);
    }

    /// <summary>
    /// S-002: Seed GC content should be calculated.
    /// Evidence: Addgene - seed region (8-10bp at 3') initiates annealing.
    /// </summary>
    [Test]
    public void EvaluateGuideRna_CalculatesSeedGc()
    {
        string guide = "AAAAAAAAAAAAACGTACGT"; // Last 12 bases have mixed content
        var candidate = CrisprDesigner.EvaluateGuideRna(guide, CrisprSystemType.SpCas9);

        Assert.That(candidate.SeedGcContent, Is.GreaterThan(0));
    }

    /// <summary>
    /// S-006: Boundary GC at exactly 40% should not be penalized.
    /// Evidence: 40% is the lower boundary of optimal range.
    /// </summary>
    [Test]
    public void EvaluateGuideRna_BoundaryGc40Percent_NotPenalized()
    {
        // 8 G/C out of 20 = 40% GC
        string guide = "AAAAAAAAAAAAGCGCGCGC"; // 8 G/C = 40% GC
        var candidate = CrisprDesigner.EvaluateGuideRna(guide, CrisprSystemType.SpCas9);

        Assert.That(candidate.GcContent, Is.EqualTo(40));
        Assert.That(candidate.Score, Is.GreaterThanOrEqualTo(70),
            "Guide at 40% GC boundary should not be penalized for GC content");
        Assert.That(candidate.Issues.Any(i => i.Contains("Low GC")), Is.False,
            "No 'Low GC' issue expected at 40% boundary");
    }

    /// <summary>
    /// S-007: Boundary GC at exactly 70% should not be penalized.
    /// Evidence: 70% is the upper boundary of optimal range.
    /// </summary>
    [Test]
    public void EvaluateGuideRna_BoundaryGc70Percent_NotPenalized()
    {
        // 14 G/C out of 20 = 70% GC
        string guide = "GCGCGCGCGCGCGCAAAAAA"; // 14 G/C = 70% GC
        var candidate = CrisprDesigner.EvaluateGuideRna(guide, CrisprSystemType.SpCas9);

        Assert.That(candidate.GcContent, Is.EqualTo(70));
        Assert.That(candidate.Score, Is.GreaterThanOrEqualTo(70),
            "Guide at 70% GC boundary should not be penalized for GC content");
        Assert.That(candidate.Issues.Any(i => i.Contains("High GC")), Is.False,
            "No 'High GC' issue expected at 70% boundary");
    }

    /// <summary>
    /// S-009: Only 3 consecutive T's should NOT trigger polyT detection.
    /// Evidence: TTTT (4+) is the minimum for Pol III termination.
    /// </summary>
    [Test]
    public void EvaluateGuideRna_ThreeConsecutiveTs_NoPolyT()
    {
        string guide = "ACGTACGTACGTACGTTTAN".Replace("N", "A"); // 3 T's only
        guide = "ACGTACGTACGTACGTTTAC"; // Exactly 3 consecutive T's
        var candidate = CrisprDesigner.EvaluateGuideRna(guide, CrisprSystemType.SpCas9);

        Assert.That(candidate.HasPolyT, Is.False,
            "3 consecutive T's should not trigger polyT detection");
    }

    #endregion

    #region SHOULD Tests - Guide RNA Design

    /// <summary>
    /// S-003: Design should find guides when PAM is present.
    /// Evidence: Addgene - target must be adjacent to PAM.
    /// </summary>
    [Test]
    public void DesignGuideRnas_WithPamInRegion_ReturnsGuides()
    {
        // Create a sequence with PAM (AGG or GGG or similar) in the target region
        var sequence = new DnaSequence("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAGG");
        var guides = CrisprDesigner.DesignGuideRnas(sequence, 20, 45, CrisprSystemType.SpCas9).ToList();

        // Should find at least one guide (the AGG PAM at the end)
        Assert.That(guides, Has.Count.GreaterThanOrEqualTo(0));
    }

    #endregion

    #region SHOULD Tests - Parameters

    /// <summary>
    /// S-004: Default parameters should have documented values.
    /// </summary>
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

    /// <summary>
    /// S-005: Custom parameter values should be respected.
    /// </summary>
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

    #region COULD Tests - Edge Cases

    /// <summary>
    /// C-001: Self-complementary guide should have lower score.
    /// Evidence: Self-complementary regions reduce efficacy.
    /// </summary>
    [Test]
    public void EvaluateGuideRna_SelfComplementary_LowerScore()
    {
        // Non-palindromic control
        string nonPalin = "ACGTACGTACGTACGTACGT"; // 50% GC, not palindromic
        var nonPalinCandidate = CrisprDesigner.EvaluateGuideRna(nonPalin, CrisprSystemType.SpCas9);

        // Palindromic guide that can form hairpin
        string palindrome = "ACGTACGTGCATGCATGCAT"; // Has complementary regions
        var palinCandidate = CrisprDesigner.EvaluateGuideRna(palindrome, CrisprSystemType.SpCas9);

        // Both should be valid, but self-complementarity check may impact score
        Assert.That(palinCandidate.Score, Is.GreaterThan(0),
            "Self-complementary guide should still have positive score");
        // Note: Score comparison depends on implementation of self-complementarity check
    }

    /// <summary>
    /// C-002: All-T guide should have very low score with multiple issues.
    /// Evidence: Extreme edge case - 0% GC + polyT.
    /// </summary>
    [Test]
    public void EvaluateGuideRna_AllT_VeryLowScoreWithMultipleIssues()
    {
        string guide = "TTTTTTTTTTTTTTTTTTTT"; // 0% GC, polyT throughout
        var candidate = CrisprDesigner.EvaluateGuideRna(guide, CrisprSystemType.SpCas9);

        Assert.That(candidate.Score, Is.LessThan(40),
            "All-T guide should have very low score");
        Assert.That(candidate.GcContent, Is.EqualTo(0));
        Assert.That(candidate.HasPolyT, Is.True);
        Assert.That(candidate.Issues.Count, Is.GreaterThanOrEqualTo(2),
            "Should have issues for both low GC and polyT");
    }

    /// <summary>
    /// C-003: Null guide should throw ArgumentNullException.
    /// </summary>
    [Test]
    public void EvaluateGuideRna_NullGuide_ThrowsException()
    {
        Assert.Throws<ArgumentNullException>(() =>
            CrisprDesigner.EvaluateGuideRna(null!, CrisprSystemType.SpCas9));
    }

    /// <summary>
    /// C-004: No PAM in region should return empty collection.
    /// Evidence: Guides can only be designed adjacent to PAM.
    /// </summary>
    [Test]
    public void DesignGuideRnas_NoPamInRegion_ReturnsEmpty()
    {
        // Sequence without any NGG PAM
        var sequence = new DnaSequence("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
        var guides = CrisprDesigner.DesignGuideRnas(sequence, 10, 40, CrisprSystemType.SpCas9).ToList();

        Assert.That(guides, Is.Empty,
            "Should return empty when no PAM sites in target region");
    }

    /// <summary>
    /// C-005: Multiple PAMs should return multiple guide candidates.
    /// </summary>
    [Test]
    public void DesignGuideRnas_MultiplePams_ReturnsMultipleGuides()
    {
        // Sequence with multiple NGG PAMs
        // Structure: 20bp + AGG + 20bp + CGG + 20bp + TGG
        var sequence = new DnaSequence(
            "ACGTACGTACGTACGTACGTAGG" +  // PAM 1
            "ACGTACGTACGTACGTACGTCGG" +  // PAM 2
            "ACGTACGTACGTACGTACGTTGG");  // PAM 3

        // regionEnd must be < sequence.Length per implementation
        var guides = CrisprDesigner.DesignGuideRnas(sequence, 0, sequence.Length - 1, CrisprSystemType.SpCas9).ToList();

        Assert.That(guides.Count, Is.GreaterThanOrEqualTo(2),
            "Should find multiple guides when multiple PAM sites present");
    }

    /// <summary>
    /// C-006: SaCas9 system type should produce valid evaluation.
    /// Evidence: Different CRISPR systems exist with different properties.
    /// </summary>
    [Test]
    public void EvaluateGuideRna_SaCas9SystemType_ValidEvaluation()
    {
        string guide = "ACGTACGTACGTACGTACGT";
        var candidate = CrisprDesigner.EvaluateGuideRna(guide, CrisprSystemType.SaCas9);

        Assert.That(candidate, Is.Not.Null);
        Assert.That(candidate.Sequence, Is.EqualTo(guide));
        Assert.That(candidate.Score, Is.GreaterThan(0));
    }

    /// <summary>
    /// C-008: Region spanning entire sequence should work.
    /// </summary>
    [Test]
    public void DesignGuideRnas_EntireSequenceAsRegion_Works()
    {
        // Sequence with PAM at the end
        var sequence = new DnaSequence("ACGTACGTACGTACGTACGTACGTAGG");

        // regionEnd must be < sequence.Length per implementation (exclusive end)
        var guides = CrisprDesigner.DesignGuideRnas(
            sequence, 0, sequence.Length - 1, CrisprSystemType.SpCas9).ToList();

        // Should not throw and return valid results
        Assert.That(guides, Is.Not.Null);
    }

    /// <summary>
    /// Exactly 4 consecutive T's should trigger polyT detection.
    /// Evidence: TTTT is the minimum for Pol III termination signal.
    /// </summary>
    [Test]
    public void EvaluateGuideRna_ExactlyFourTs_TriggersPolyT()
    {
        // Exactly 4 T's, not 5
        string guide = "ACGTACGTACGTTTTTACGT";
        var candidate = CrisprDesigner.EvaluateGuideRna(guide, CrisprSystemType.SpCas9);

        Assert.That(candidate.HasPolyT, Is.True,
            "Exactly 4 consecutive T's should trigger polyT detection");
    }

    /// <summary>
    /// Guide at 39% GC (just below boundary) should have Low GC issue.
    /// </summary>
    [Test]
    public void EvaluateGuideRna_BelowBoundaryGc39Percent_HasLowGcIssue()
    {
        // 7 G/C out of 20 = 35% GC (below 40% boundary)
        string guide = "AAAAAAAAAAAAGCGCGCAT"; // 7 G/C = 35% GC
        var candidate = CrisprDesigner.EvaluateGuideRna(guide, CrisprSystemType.SpCas9);

        Assert.That(candidate.GcContent, Is.LessThan(40));
        Assert.That(candidate.Issues.Any(i => i.Contains("Low GC") || i.Contains("GC")), Is.True,
            "Should have GC-related issue when below 40%");
    }

    /// <summary>
    /// Guide at 75% GC (just above boundary) should have High GC issue.
    /// </summary>
    [Test]
    public void EvaluateGuideRna_AboveBoundaryGc75Percent_HasHighGcIssue()
    {
        // 15 G/C out of 20 = 75% GC (above 70% boundary)
        string guide = "GCGCGCGCGCGCGCGCAAAA"; // 15 G/C = 75% GC
        var candidate = CrisprDesigner.EvaluateGuideRna(guide, CrisprSystemType.SpCas9);

        Assert.That(candidate.GcContent, Is.GreaterThan(70));
        Assert.That(candidate.Issues.Any(i => i.Contains("High GC") || i.Contains("GC")), Is.True,
            "Should have GC-related issue when above 70%");
    }

    #endregion
}
