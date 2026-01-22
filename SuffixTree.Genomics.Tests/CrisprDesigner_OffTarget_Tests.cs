using NUnit.Framework;
using SuffixTree.Genomics;
using System;
using System.Linq;

namespace SuffixTree.Genomics.Tests;

/// <summary>
/// Tests for CRISPR Off-Target Analysis (CRISPR-OFF-001).
/// 
/// Evidence Sources:
/// - Wikipedia: Off-target genome editing
/// - Hsu et al. (2013) Nature Biotechnology - DNA targeting specificity of RNA-guided Cas9 nucleases
/// - Fu et al. (2013) Nature Biotechnology - High-frequency off-target mutagenesis
/// 
/// Key Evidence:
/// - Off-target sites have sequence similarity to guide with 1-5 mismatches
/// - Seed region (PAM-proximal 10-12bp) is critical for specificity
/// - PAM is required at off-target sites for cleavage
/// </summary>
[TestFixture]
public class CrisprDesigner_OffTarget_Tests
{
    #region Input Validation Tests (M-001 to M-003)

    /// <summary>
    /// M-001: Empty guide should throw ArgumentNullException.
    /// Evidence: Defensive programming - null/empty input undefined.
    /// </summary>
    [Test]
    public void FindOffTargets_EmptyGuide_ThrowsArgumentNullException()
    {
        var genome = new DnaSequence("ACGTACGTACGTACGTACGTACGTAGG");

        Assert.Throws<ArgumentNullException>(() =>
            CrisprDesigner.FindOffTargets("", genome, 3, CrisprSystemType.SpCas9).ToList());
    }

    /// <summary>
    /// M-002: Null genome should throw ArgumentNullException.
    /// Evidence: Defensive programming.
    /// </summary>
    [Test]
    public void FindOffTargets_NullGenome_ThrowsArgumentNullException()
    {
        Assert.Throws<ArgumentNullException>(() =>
            CrisprDesigner.FindOffTargets("ACGTACGTACGTACGTACGT", null!, 3, CrisprSystemType.SpCas9).ToList());
    }

    /// <summary>
    /// M-003: MaxMismatches > 5 should throw.
    /// Evidence: Hsu et al. (2013) - practical limit is 5 mismatches for detectable off-targets.
    /// </summary>
    [TestCase(-1)]
    [TestCase(6)]
    [TestCase(10)]
    public void FindOffTargets_InvalidMaxMismatches_ThrowsArgumentOutOfRangeException(int maxMismatches)
    {
        var genome = new DnaSequence("ACGTACGTACGTACGTACGTACGTAGG");

        Assert.Throws<ArgumentOutOfRangeException>(() =>
            CrisprDesigner.FindOffTargets("ACGTACGTACGTACGTACGT", genome, maxMismatches, CrisprSystemType.SpCas9).ToList());
    }

    #endregion

    #region Core Off-Target Detection Tests (M-004 to M-009)

    /// <summary>
    /// M-004: Exact matches are on-targets, not off-targets.
    /// Evidence: Hsu et al. (2013) - off-targets are sites with mismatches.
    /// </summary>
    [Test]
    public void FindOffTargets_ExactMatch_NotReturnedAsOffTarget()
    {
        // Guide sequence with exact match in genome (including PAM)
        string guide = "ACGTACGTACGTACGTACGT"; // 20bp guide
        // Genome: exact target + NGG PAM
        var genome = new DnaSequence("ACGTACGTACGTACGTACGTAGG");

        var offTargets = CrisprDesigner.FindOffTargets(guide, genome, 3, CrisprSystemType.SpCas9).ToList();

        // All returned off-targets must have mismatches > 0
        Assert.That(offTargets.All(ot => ot.Mismatches > 0), Is.True,
            "Exact matches should not be returned as off-targets");
    }

    /// <summary>
    /// M-005: Single mismatch within maxMismatches returns off-target.
    /// Evidence: Hsu et al. (2013) - single mismatches are tolerated, especially in PAM-distal region.
    /// </summary>
    [Test]
    public void FindOffTargets_SingleMismatch_ReturnsOffTarget()
    {
        string guide = "ACGTACGTACGTACGTACGT"; // 20bp guide
        // Genome with 1 mismatch at position 0 (T instead of A) + NGG PAM
        var genome = new DnaSequence("TCGTACGTACGTACGTACGTAGG");

        var offTargets = CrisprDesigner.FindOffTargets(guide, genome, 3, CrisprSystemType.SpCas9).ToList();

        Assert.That(offTargets, Has.Count.GreaterThanOrEqualTo(1));
        var offTarget = offTargets.First();
        Assert.That(offTarget.Mismatches, Is.EqualTo(1));
    }

    /// <summary>
    /// M-006: MaxMismatches is respected - no off-targets with more mismatches returned.
    /// Evidence: Algorithm specification.
    /// </summary>
    [Test]
    public void FindOffTargets_MaxMismatchesRespected_AllResultsWithinLimit()
    {
        string guide = "ACGTACGTACGTACGTACGT";
        // Genome with multiple sites having different mismatch counts
        // Site with many mismatches should not be returned when maxMismatches=2
        var genome = new DnaSequence("AAAAAAAAAAAAAAAAAAAAAGG");

        var offTargets = CrisprDesigner.FindOffTargets(guide, genome, 2, CrisprSystemType.SpCas9).ToList();

        Assert.That(offTargets.All(ot => ot.Mismatches <= 2), Is.True,
            "All off-targets should have mismatches <= maxMismatches");
    }

    /// <summary>
    /// M-007: MismatchPositions count equals Mismatches count.
    /// Evidence: Hsu et al. (2013) - position of mismatches affects activity, must be tracked.
    /// </summary>
    [Test]
    public void FindOffTargets_MismatchPositions_CountMatchesMismatches()
    {
        string guide = "ACGTACGTACGTACGTACGT";
        // Genome with 2 mismatches at known positions + NGG PAM
        // Positions 0 and 1 are different (TT instead of AC)
        var genome = new DnaSequence("TTGTACGTACGTACGTACGTAGG");

        var offTargets = CrisprDesigner.FindOffTargets(guide, genome, 3, CrisprSystemType.SpCas9).ToList();

        Assert.That(offTargets, Has.Count.GreaterThanOrEqualTo(1));

        foreach (var ot in offTargets)
        {
            Assert.Multiple(() =>
            {
                Assert.That(ot.MismatchPositions, Is.Not.Null);
                Assert.That(ot.MismatchPositions.Count, Is.EqualTo(ot.Mismatches),
                    $"MismatchPositions.Count should equal Mismatches for off-target at position {ot.Position}");
            });
        }
    }

    /// <summary>
    /// M-007b: MismatchPositions contains correct positions.
    /// Evidence: Positions must be accurately reported for scoring.
    /// </summary>
    [Test]
    public void FindOffTargets_MismatchPositions_ContainsCorrectPositions()
    {
        string guide = "ACGTACGTACGTACGTACGT";
        // Single mismatch at position 0 (T instead of A)
        var genome = new DnaSequence("TCGTACGTACGTACGTACGTAGG");

        var offTargets = CrisprDesigner.FindOffTargets(guide, genome, 3, CrisprSystemType.SpCas9).ToList();

        Assert.That(offTargets, Has.Count.GreaterThanOrEqualTo(1));
        var offTarget = offTargets.First(ot => ot.Mismatches == 1);
        Assert.That(offTarget.MismatchPositions, Does.Contain(0),
            "Mismatch at position 0 should be recorded");
    }

    /// <summary>
    /// M-008: Off-targets require PAM at the site.
    /// Evidence: Wikipedia, Hsu et al. (2013) - PAM is required for Cas9 cleavage.
    /// </summary>
    [Test]
    public void FindOffTargets_NoPam_NoOffTargetReturned()
    {
        string guide = "ACGTACGTACGTACGTACGT";
        // Genome with similar sequence but NO valid PAM (ending in TTT instead of NGG)
        var genome = new DnaSequence("TCGTACGTACGTACGTACGTTTT");

        var offTargets = CrisprDesigner.FindOffTargets(guide, genome, 3, CrisprSystemType.SpCas9).ToList();

        Assert.That(offTargets, Is.Empty,
            "No off-targets should be found when PAM is not present");
    }

    /// <summary>
    /// M-009: Off-targets found on reverse strand.
    /// Evidence: CRISPR can target either strand.
    /// </summary>
    [Test]
    public void FindOffTargets_ReverseStrand_ReturnsOffTarget()
    {
        // For reverse strand off-target:
        // Guide: AAAAAAAAAAAAAAAAAAAA (20 A's)
        // Off-target on reverse strand would be: guide's revcomp = TTTTTTTTTTTTTTTTTTTT
        // With 1 mismatch: ATTTTTTTTTTTTTTTTTTT + CCG (PAM for reverse is CCN on forward)
        // So genome: CCG + ATTTTTTTTTTTTTTTTTTT = CCGATTTTTTTTTTTTTTTTTTT (23bp)

        string testGuide = "AAAAAAAAAAAAAAAAAAAA"; // 20 A's
        var genome = new DnaSequence("CCGATTTTTTTTTTTTTTTTTTT");

        var offTargets = CrisprDesigner.FindOffTargets(testGuide, genome, 3, CrisprSystemType.SpCas9).ToList();

        // Should find at least one off-target on reverse strand
        Assert.That(offTargets.Any(ot => !ot.IsForwardStrand), Is.True,
            "Should find off-target on reverse strand");
    }

    #endregion

    #region Specificity Score Tests (M-010 to M-012)

    /// <summary>
    /// M-010: SpecificityScore returns value in valid range.
    /// Evidence: Score should be normalized percentage.
    /// </summary>
    [Test]
    public void CalculateSpecificityScore_ReturnsValueInValidRange()
    {
        string guide = "ACGTACGTACGTACGTACGT";
        var genome = new DnaSequence("ACGTACGTACGTACGTACGTACGTAGG");

        double score = CrisprDesigner.CalculateSpecificityScore(guide, genome, CrisprSystemType.SpCas9);

        Assert.Multiple(() =>
        {
            Assert.That(score, Is.GreaterThanOrEqualTo(0), "Score should be >= 0");
            Assert.That(score, Is.LessThanOrEqualTo(100), "Score should be <= 100");
        });
    }

    /// <summary>
    /// M-011: No off-targets returns maximum specificity score.
    /// Evidence: No off-targets = highest specificity.
    /// </summary>
    [Test]
    public void CalculateSpecificityScore_NoOffTargets_Returns100()
    {
        string guide = "ACGTACGTACGTACGTACGT";
        // Very short genome with only exact match (not off-target)
        var genome = new DnaSequence("ACGTACGTACGTACGTACGTAGG");

        double score = CrisprDesigner.CalculateSpecificityScore(guide, genome, CrisprSystemType.SpCas9);

        Assert.That(score, Is.EqualTo(100),
            "Score should be 100 when no off-targets exist");
    }

    /// <summary>
    /// M-012: Off-targets reduce specificity score.
    /// Evidence: More off-targets = lower specificity.
    /// </summary>
    [Test]
    public void CalculateSpecificityScore_WithOffTargets_ScoreReducedFromMaximum()
    {
        string guide = "ACGTACGTACGTACGTACGT";
        // Genome with off-target site (single mismatch)
        var genome = new DnaSequence("TCGTACGTACGTACGTACGTAGG");

        double score = CrisprDesigner.CalculateSpecificityScore(guide, genome, CrisprSystemType.SpCas9);

        Assert.That(score, Is.LessThan(100),
            "Score should be less than 100 when off-targets exist");
    }

    #endregion

    #region Position-Dependent Scoring Tests (S-001, S-004)

    /// <summary>
    /// S-001: Seed region mismatches receive higher penalty score.
    /// Evidence: Hsu et al. (2013) - PAM-proximal mismatches less tolerated.
    /// Implementation: Seed mismatches score 5 points vs 2 for distal.
    /// </summary>
    [Test]
    public void FindOffTargets_SeedMismatch_HigherOffTargetScore()
    {
        // For SpCas9, seed is last 12bp (positions 8-19 of 20bp guide)
        string guide = "ACGTACGTACGTACGTACGT";

        // Off-target with mismatch at position 0 (PAM-distal)
        var genomeDistal = new DnaSequence("TCGTACGTACGTACGTACGTAGG");

        // Off-target with mismatch at position 19 (PAM-proximal/seed - last position)
        var genomeSeed = new DnaSequence("ACGTACGTACGTACGTACGAAGG");

        var offTargetsDistal = CrisprDesigner.FindOffTargets(guide, genomeDistal, 3, CrisprSystemType.SpCas9).ToList();
        var offTargetsSeed = CrisprDesigner.FindOffTargets(guide, genomeSeed, 3, CrisprSystemType.SpCas9).ToList();

        Assert.That(offTargetsDistal, Has.Count.GreaterThanOrEqualTo(1), "Should find distal off-target");
        Assert.That(offTargetsSeed, Has.Count.GreaterThanOrEqualTo(1), "Should find seed off-target");

        var distalScore = offTargetsDistal.First().OffTargetScore;
        var seedScore = offTargetsSeed.First().OffTargetScore;

        Assert.That(seedScore, Is.GreaterThan(distalScore),
            $"Seed mismatch score ({seedScore}) should be greater than distal mismatch score ({distalScore})");
    }

    #endregion

    #region Multiple Mismatches Tests (S-002)

    /// <summary>
    /// S-002: Multiple mismatches are correctly counted and reported.
    /// Evidence: Hsu et al. (2013) - aggregate effect of multiple mismatches.
    /// </summary>
    [Test]
    public void FindOffTargets_MultipleMismatches_AllReported()
    {
        string guide = "ACGTACGTACGTACGTACGT";
        // 3 mismatches at positions 0, 1, 2 (TTT instead of ACG)
        var genome = new DnaSequence("TTTTACGTACGTACGTACGTAGG");

        var offTargets = CrisprDesigner.FindOffTargets(guide, genome, 3, CrisprSystemType.SpCas9).ToList();

        Assert.That(offTargets, Has.Count.GreaterThanOrEqualTo(1));
        var offTarget = offTargets.First(ot => ot.Mismatches == 3);

        Assert.Multiple(() =>
        {
            Assert.That(offTarget.Mismatches, Is.EqualTo(3));
            Assert.That(offTarget.MismatchPositions.Count, Is.EqualTo(3));
            Assert.That(offTarget.MismatchPositions, Does.Contain(0));
            Assert.That(offTarget.MismatchPositions, Does.Contain(1));
            Assert.That(offTarget.MismatchPositions, Does.Contain(2));
        });
    }

    #endregion

    #region Different CRISPR Systems (S-003)

    /// <summary>
    /// S-003: Cas12a system uses correct PAM and guide length.
    /// Evidence: Cas12a uses TTTV PAM before target, 23bp guide.
    /// </summary>
    [Test]
    public void FindOffTargets_Cas12a_UsesCorrectParameters()
    {
        // Cas12a: PAM (TTTA/C/G) is BEFORE the target, guide is 23bp
        string guide = "ACGTACGTACGTACGTACGTACG"; // 23bp guide

        // Genome: TTTA (PAM) + 23bp with 1 mismatch
        var genome = new DnaSequence("TTTATCGTACGTACGTACGTACGTACG");

        var offTargets = CrisprDesigner.FindOffTargets(guide, genome, 3, CrisprSystemType.Cas12a).ToList();

        Assert.That(offTargets, Has.Count.GreaterThanOrEqualTo(1),
            "Should find off-target with Cas12a system");
        Assert.That(offTargets.First().Mismatches, Is.EqualTo(1));
    }

    #endregion

    #region Edge Cases

    /// <summary>
    /// Edge case: Empty genome returns empty results (no crash).
    /// ASSUMPTION A-001: Empty genome handled gracefully.
    /// </summary>
    [Test]
    public void FindOffTargets_EmptyGenome_ReturnsEmpty()
    {
        string guide = "ACGTACGTACGTACGTACGT";
        var genome = new DnaSequence("");

        var offTargets = CrisprDesigner.FindOffTargets(guide, genome, 3, CrisprSystemType.SpCas9).ToList();

        Assert.That(offTargets, Is.Empty);
    }

    /// <summary>
    /// Edge case: Genome shorter than guide + PAM returns empty.
    /// </summary>
    [Test]
    public void FindOffTargets_GenomeTooShort_ReturnsEmpty()
    {
        string guide = "ACGTACGTACGTACGTACGT"; // 20bp
        var genome = new DnaSequence("ACGT"); // Only 4bp

        var offTargets = CrisprDesigner.FindOffTargets(guide, genome, 3, CrisprSystemType.SpCas9).ToList();

        Assert.That(offTargets, Is.Empty);
    }

    /// <summary>
    /// Edge case: MaxMismatches = 0 should return empty (0 mismatches = exact match = on-target).
    /// </summary>
    [Test]
    public void FindOffTargets_MaxMismatchesZero_ReturnsEmpty()
    {
        string guide = "ACGTACGTACGTACGTACGT";
        // Genome with 1 mismatch
        var genome = new DnaSequence("TCGTACGTACGTACGTACGTAGG");

        var offTargets = CrisprDesigner.FindOffTargets(guide, genome, 0, CrisprSystemType.SpCas9).ToList();

        Assert.That(offTargets, Is.Empty,
            "With maxMismatches=0, no off-targets should be found (exact match is not off-target)");
    }

    #endregion

    #region Invariant Tests

    /// <summary>
    /// Invariant: All off-targets have OffTargetScore >= 0.
    /// </summary>
    [Test]
    public void FindOffTargets_OffTargetScore_IsNonNegative()
    {
        string guide = "ACGTACGTACGTACGTACGT";
        var genome = new DnaSequence("TTTTACGTACGTACGTACGTAGG");

        var offTargets = CrisprDesigner.FindOffTargets(guide, genome, 5, CrisprSystemType.SpCas9).ToList();

        Assert.That(offTargets.All(ot => ot.OffTargetScore >= 0), Is.True,
            "All off-target scores should be non-negative");
    }

    /// <summary>
    /// Invariant: Results are deterministic (same input = same output).
    /// </summary>
    [Test]
    public void FindOffTargets_SameInput_DeterministicOutput()
    {
        string guide = "ACGTACGTACGTACGTACGT";
        var genome = new DnaSequence("TTGTACGTACGTACGTACGTAGG");

        var offTargets1 = CrisprDesigner.FindOffTargets(guide, genome, 3, CrisprSystemType.SpCas9).ToList();
        var offTargets2 = CrisprDesigner.FindOffTargets(guide, genome, 3, CrisprSystemType.SpCas9).ToList();

        Assert.That(offTargets1.Count, Is.EqualTo(offTargets2.Count));
        for (int i = 0; i < offTargets1.Count; i++)
        {
            Assert.Multiple(() =>
            {
                Assert.That(offTargets1[i].Position, Is.EqualTo(offTargets2[i].Position));
                Assert.That(offTargets1[i].Mismatches, Is.EqualTo(offTargets2[i].Mismatches));
                Assert.That(offTargets1[i].OffTargetScore, Is.EqualTo(offTargets2[i].OffTargetScore));
            });
        }
    }

    #endregion
}
