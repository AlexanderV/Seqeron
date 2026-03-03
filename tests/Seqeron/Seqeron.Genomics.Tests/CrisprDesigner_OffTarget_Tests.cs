using NUnit.Framework;
using Seqeron.Genomics;
using System;
using System.Linq;

namespace Seqeron.Genomics.Tests;

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
/// - Seed region (PAM-proximal 8-10bp per Addgene) is critical for specificity
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

    /// <summary>
    /// M-003b: Guide length must match system's expected length.
    /// Evidence: SpCas9 uses 20bp guide (Hsu 2013); Cas12a uses 23bp.
    /// Mismatched guide length is an error, not silently ignored.
    /// </summary>
    [TestCase("ACGTACGTACGT", CrisprSystemType.SpCas9, Description = "12bp guide for SpCas9 (expects 20bp)")]
    [TestCase("ACGTACGTACGTACGTACGTACGT", CrisprSystemType.SpCas9, Description = "24bp guide for SpCas9 (expects 20bp)")]
    [TestCase("ACGTACGTACGTACGTACGT", CrisprSystemType.Cas12a, Description = "20bp guide for Cas12a (expects 23bp)")]
    public void FindOffTargets_GuideLengthMismatch_ThrowsArgumentException(string guide, CrisprSystemType system)
    {
        var genome = new DnaSequence("ACGTACGTACGTACGTACGTACGTACGTACGTAGG");

        Assert.Throws<ArgumentException>(() =>
            CrisprDesigner.FindOffTargets(guide, genome, 3, system).ToList());
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
        // Genome: exact target + NGG PAM — only one PAM site, exact match
        var genome = new DnaSequence("ACGTACGTACGTACGTACGTAGG");

        var offTargets = CrisprDesigner.FindOffTargets(guide, genome, 3, CrisprSystemType.SpCas9).ToList();

        Assert.That(offTargets, Is.Empty,
            "Exact match is on-target, not off-target — collection must be empty");
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

        Assert.That(offTargets, Has.Count.EqualTo(1), "Exactly one off-target site");
        var ot = offTargets[0];
        Assert.Multiple(() =>
        {
            Assert.That(ot.Mismatches, Is.EqualTo(1));
            Assert.That(ot.MismatchPositions, Is.EqualTo(new[] { 0 }), "Mismatch at position 0");
            Assert.That(ot.IsForwardStrand, Is.True);
            Assert.That(ot.OffTargetScore, Is.EqualTo(2.0),
                "Position 0 is distal (outside seed region 8-19), penalty = 2");
        });
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

        Assert.That(offTargets, Has.Count.EqualTo(1), "Exactly one off-target site");
        var ot = offTargets[0];
        Assert.Multiple(() =>
        {
            Assert.That(ot.Mismatches, Is.EqualTo(2));
            Assert.That(ot.MismatchPositions, Is.EqualTo(new[] { 0, 1 }),
                "Mismatches at positions 0 (A→T) and 1 (C→T)");
            Assert.That(ot.OffTargetScore, Is.EqualTo(4.0),
                "Both positions are distal (outside seed 8-19), 2 × 2 = 4");
        });
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

        Assert.That(offTargets, Has.Count.EqualTo(1), "Exactly one off-target site");
        Assert.That(offTargets[0].MismatchPositions, Is.EqualTo(new[] { 0 }),
            "Single mismatch at position 0 (A→T)");
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
        // Reverse-strand off-target construction:
        // Guide: 20× A. RevComp target on forward strand = CCN(PAM) + near-match.
        // Genome: CCG + ATTTTTTTTTTTTTTTTTTT = 23bp.
        // RevComp: AAAAAAAAAAAAAAAAAAAAT + CGG → target = AAAAAAAAAAAAAAAAAAAAT,
        //   mismatch at position 19 (A vs T), which is in seed region (8-19).
        string testGuide = "AAAAAAAAAAAAAAAAAAAA"; // 20 A's
        var genome = new DnaSequence("CCGATTTTTTTTTTTTTTTTTTT");

        var offTargets = CrisprDesigner.FindOffTargets(testGuide, genome, 3, CrisprSystemType.SpCas9).ToList();

        Assert.That(offTargets, Has.Count.EqualTo(1), "Exactly one reverse-strand off-target");
        var ot = offTargets[0];
        Assert.Multiple(() =>
        {
            Assert.That(ot.IsForwardStrand, Is.False, "Off-target is on reverse strand");
            Assert.That(ot.Mismatches, Is.EqualTo(1));
            Assert.That(ot.MismatchPositions, Is.EqualTo(new[] { 19 }),
                "Mismatch at position 19 (last position)");
            Assert.That(ot.OffTargetScore, Is.EqualTo(5.0),
                "Position 19 is in seed region (8-19), penalty = 5");
        });
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
    /// Genome has 1 off-target with distal mismatch (score=2), so SpecificityScore = 100 − 2 = 98.
    /// </summary>
    [Test]
    public void CalculateSpecificityScore_WithOffTargets_ScoreReducedFromMaximum()
    {
        string guide = "ACGTACGTACGTACGTACGT";
        // Genome with off-target site (single mismatch at position 0, distal)
        var genome = new DnaSequence("TCGTACGTACGTACGTACGTAGG");

        double score = CrisprDesigner.CalculateSpecificityScore(guide, genome, CrisprSystemType.SpCas9);

        Assert.That(score, Is.EqualTo(98.0),
            "1 off-target with distal mismatch (penalty=2): 100 − 2 = 98");
    }

    #endregion

    #region Position-Dependent Scoring Tests (S-001, S-004)

    /// <summary>
    /// S-001: Seed region mismatches receive higher penalty score.
    /// Evidence: Hsu et al. (2013) - PAM-proximal mismatches less tolerated.
    /// Implementation: Seed mismatches score 5 points vs 2 for distal.
    /// Seed region = last 12bp (positions 8-19) for SpCas9 (PamAfterTarget=true).
    /// </summary>
    [Test]
    public void FindOffTargets_SeedMismatch_HigherOffTargetScore()
    {
        string guide = "ACGTACGTACGTACGTACGT";

        // Off-target with mismatch at position 0 (PAM-distal, outside seed 8-19)
        var genomeDistal = new DnaSequence("TCGTACGTACGTACGTACGTAGG");

        // Off-target with mismatch at position 19 (PAM-proximal/seed, last position)
        var genomeSeed = new DnaSequence("ACGTACGTACGTACGTACGAAGG");

        var offTargetsDistal = CrisprDesigner.FindOffTargets(guide, genomeDistal, 3, CrisprSystemType.SpCas9).ToList();
        var offTargetsSeed = CrisprDesigner.FindOffTargets(guide, genomeSeed, 3, CrisprSystemType.SpCas9).ToList();

        Assert.That(offTargetsDistal, Has.Count.EqualTo(1), "Exactly one distal off-target");
        Assert.That(offTargetsSeed, Has.Count.EqualTo(1), "Exactly one seed off-target");

        var distalScore = offTargetsDistal[0].OffTargetScore;
        var seedScore = offTargetsSeed[0].OffTargetScore;

        Assert.Multiple(() =>
        {
            Assert.That(distalScore, Is.EqualTo(2.0), "Distal mismatch penalty = 2");
            Assert.That(seedScore, Is.EqualTo(5.0), "Seed mismatch penalty = 5");
            Assert.That(seedScore, Is.GreaterThan(distalScore),
                "Seed mismatch score should exceed distal mismatch score");
        });
    }

    #endregion

    /// <summary>
    /// S-004: CalculateSpecificityScore penalizes seed mismatches more than distal.
    /// Evidence: Hsu et al. (2013) — seed (PAM-proximal) mismatches contribute higher penalty.
    /// Distal mismatch (pos 0, penalty=2) → score 98. Seed mismatch (pos 19, penalty=5) → score 95.
    /// </summary>
    [Test]
    public void CalculateSpecificityScore_SeedMismatch_LowerThanDistal()
    {
        string guide = "ACGTACGTACGTACGTACGT";

        // Genome with 1 distal mismatch (position 0): penalty=2 → score=98
        var genomeDistal = new DnaSequence("TCGTACGTACGTACGTACGTAGG");

        // Genome with 1 seed mismatch (position 19): penalty=5 → score=95
        var genomeSeed = new DnaSequence("ACGTACGTACGTACGTACGAAGG");

        double scoreDistal = CrisprDesigner.CalculateSpecificityScore(guide, genomeDistal, CrisprSystemType.SpCas9);
        double scoreSeed = CrisprDesigner.CalculateSpecificityScore(guide, genomeSeed, CrisprSystemType.SpCas9);

        Assert.Multiple(() =>
        {
            Assert.That(scoreDistal, Is.EqualTo(98.0), "Distal mismatch penalty=2: 100 − 2 = 98");
            Assert.That(scoreSeed, Is.EqualTo(95.0), "Seed mismatch penalty=5: 100 − 5 = 95");
            Assert.That(scoreSeed, Is.LessThan(scoreDistal),
                "Seed mismatch should reduce SpecificityScore more than distal mismatch");
        });
    }

#endregion

    #region Multiple Mismatches Tests (S-002)

    /// <summary>
    /// S-002: Multiple mismatches are correctly counted and reported.
    /// Evidence: Hsu et al. (2013) - aggregate effect of multiple mismatches.
    /// 3 distal mismatches at positions 0, 1, 2 → score = 3 × 2 = 6.
    /// </summary>
    [Test]
    public void FindOffTargets_MultipleMismatches_AllReported()
    {
        string guide = "ACGTACGTACGTACGTACGT";
        // 3 mismatches at positions 0, 1, 2 (TTT instead of ACG)
        var genome = new DnaSequence("TTTTACGTACGTACGTACGTAGG");

        var offTargets = CrisprDesigner.FindOffTargets(guide, genome, 3, CrisprSystemType.SpCas9).ToList();

        Assert.That(offTargets, Has.Count.EqualTo(1), "Exactly one off-target site");
        var ot = offTargets[0];
        Assert.Multiple(() =>
        {
            Assert.That(ot.Mismatches, Is.EqualTo(3));
            Assert.That(ot.MismatchPositions, Is.EqualTo(new[] { 0, 1, 2 }),
                "Mismatches at positions 0 (A→T), 1 (C→T), 2 (G→T)");
            Assert.That(ot.OffTargetScore, Is.EqualTo(6.0),
                "All 3 mismatches are distal (outside seed 8-19): 3 × 2 = 6");
        });
    }

    #endregion

    #region Different CRISPR Systems (S-003)

    /// <summary>
    /// S-003: Cas12a system uses correct PAM and guide length.
    /// Evidence: Cas12a uses TTTV PAM before target, 23bp guide.
    /// Seed region for Cas12a (PamAfterTarget=false): first 12bp (positions 0-11).
    /// Mismatch at position 0 is in seed → score = 5.
    /// </summary>
    [Test]
    public void FindOffTargets_Cas12a_UsesCorrectParameters()
    {
        // Cas12a: PAM (TTTA/C/G) is BEFORE the target, guide is 23bp
        string guide = "ACGTACGTACGTACGTACGTACG"; // 23bp guide

        // Genome: TTTA (PAM) + 23bp with 1 mismatch at position 0 (T instead of A)
        var genome = new DnaSequence("TTTATCGTACGTACGTACGTACGTACG");

        var offTargets = CrisprDesigner.FindOffTargets(guide, genome, 3, CrisprSystemType.Cas12a).ToList();

        Assert.That(offTargets, Has.Count.EqualTo(1), "Exactly one off-target with Cas12a");
        var ot = offTargets[0];
        Assert.Multiple(() =>
        {
            Assert.That(ot.Mismatches, Is.EqualTo(1));
            Assert.That(ot.MismatchPositions, Is.EqualTo(new[] { 0 }));
            Assert.That(ot.OffTargetScore, Is.EqualTo(5.0),
                "Cas12a seed is positions 0-11; position 0 is seed → penalty = 5");
        });
    }

    #endregion

    #region Edge Cases

    /// <summary>
    /// Edge case: Empty genome returns empty results (no crash).
    /// Natural behavior: no PAM sites found → no off-targets.
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
