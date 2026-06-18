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

    #region Seed vs Distal Specificity Tests (S-004)

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

    #region MIT / Hsu 2013 Score Tests (CRISPR-OFF-001, C7)

    // Source of the model + weights (re-grounded this session):
    //   Hsu et al. Nat Biotechnol 31:827-832 (2013), PMID 23873081; scoring scheme per
    //   crispr.mit.edu, transcribed in CRISPOR's calcHitScore / calcMitGuideScore:
    //   https://github.com/maximilianh/crisporWebsite/blob/master/crispor.py
    // Published W vector (index 0 = PAM-distal .. 19 = PAM-proximal):
    //   [0,0,0.014,0,0,0.395,0.317,0,0.389,0.079,0.445,0.508,0.613,0.851,0.732,0.828,0.615,0.804,0.685,0.583]
    // All expected values below come from an INDEPENDENT Python run of the reference formula,
    // not from this C# code. Reference guide used: "GACGCATAAAGATGAGACGC".

    private const string MitGuide = "GACGCATAAAGATGAGACGC"; // 20 nt; index 5 = 'A', 15 = 'G', 19 = 'C'

    /// <summary>
    /// MIT-001 (boundary): an exact match scores exactly 100 (on-target, no mismatches).
    /// </summary>
    [Test]
    public void CalculateMitHitScore_PerfectMatch_Returns100()
    {
        double score = CrisprDesigner.CalculateMitHitScore(MitGuide, MitGuide);
        Assert.That(score, Is.EqualTo(100.0).Within(1e-9));
    }

    /// <summary>
    /// MIT-002: single mismatch at PAM-distal position 0 (W[0]=0) → 100*(1-0)=100.
    /// Verifies the zero-weight positions contribute no penalty (only nmm term, =1 for 1 mm).
    /// </summary>
    [Test]
    public void CalculateMitHitScore_SingleMismatchPos0_ZeroWeight_Returns100()
    {
        string ot = "T" + MitGuide.Substring(1); // pos 0: G→T
        double score = CrisprDesigner.CalculateMitHitScore(MitGuide, ot);
        Assert.That(score, Is.EqualTo(100.0).Within(1e-9));
    }

    /// <summary>
    /// MIT-003 (known penalty): single mismatch at position 5 (W[5]=0.395) → 100*(1-0.395)=60.5.
    /// </summary>
    [Test]
    public void CalculateMitHitScore_SingleMismatchPos5_Returns60Point5()
    {
        var chars = MitGuide.ToCharArray();
        chars[5] = 'T'; // guide[5]='A' → 'T'
        double score = CrisprDesigner.CalculateMitHitScore(MitGuide, new string(chars));
        Assert.That(score, Is.EqualTo(60.5).Within(1e-9));
    }

    /// <summary>
    /// MIT-004 (known penalty, PAM-proximal): single mismatch at position 19 (W[19]=0.583)
    /// → 100*(1-0.583)=41.7. PAM-proximal mismatch is penalized far more than the distal pos-0.
    /// </summary>
    [Test]
    public void CalculateMitHitScore_SingleMismatchPos19_Returns41Point7()
    {
        var chars = MitGuide.ToCharArray();
        chars[19] = 'A'; // guide[19]='C' → 'A'
        double score = CrisprDesigner.CalculateMitHitScore(MitGuide, new string(chars));
        Assert.That(score, Is.EqualTo(41.7).Within(1e-9));
    }

    /// <summary>
    /// MIT-005 (two mismatches → all three terms active): positions 5 and 15
    /// (W[5]=0.395, W[15]=0.828). Independent Python computation:
    ///   score1 = (1-0.395)*(1-0.828) = 0.10406
    ///   meanDist = 15-5 = 10; score2 = 1/(((19-10)/19)*4 + 1) = 0.34545454545...
    ///   score3 = 1/(2^2) = 0.25
    ///   hitScore = 0.10406 * 0.345454.. * 0.25 * 100 = 0.8987
    /// </summary>
    [Test]
    public void CalculateMitHitScore_TwoMismatches_AllTermsActive()
    {
        var chars = MitGuide.ToCharArray();
        chars[5] = 'T';  // 'A' → 'T'
        chars[15] = 'A'; // 'G' → 'A'
        double score = CrisprDesigner.CalculateMitHitScore(MitGuide, new string(chars));
        Assert.That(score, Is.EqualTo(0.8987).Within(1e-9));
    }

    /// <summary>
    /// MIT-005b (W-vector ORIENTATION guard): pins index 0 = PAM-distal (5') .. index 19 =
    /// PAM-proximal (3'). A single mismatch at the maximum-weight position 13 (W[13]=0.851)
    /// must score 100*(1-0.851)=14.9. If the published W vector were transcribed REVERSED
    /// (index 0 = PAM-proximal), position 13 would map to W=0.317 and this would score 68.3,
    /// failing this test. Conversely a PAM-distal mismatch at the zero-weight position 0 must
    /// score 100. Asserting BOTH locks the orientation: the PAM-proximal/seed end is the
    /// high-penalty end, exactly as in CRISPOR calcHitScore (hitScoreM indexed 5'→3').
    /// Independent Python: 100*(1-0.851)=14.9; 100*(1-0)=100.
    /// </summary>
    [Test]
    public void CalculateMitHitScore_WeightOrientation_PamProximalIsHighPenalty()
    {
        // Max-weight PAM-proximal position 13 (W=0.851) — reversed vector would give 68.3.
        var proximal = MitGuide.ToCharArray();
        proximal[13] = proximal[13] == 'A' ? 'C' : 'A'; // ensure a real mismatch
        double proximalScore = CrisprDesigner.CalculateMitHitScore(MitGuide, new string(proximal));

        // Zero-weight PAM-distal position 0 (W=0) — reversed vector would give 41.7.
        var distal = MitGuide.ToCharArray();
        distal[0] = distal[0] == 'A' ? 'C' : 'A';
        double distalScore = CrisprDesigner.CalculateMitHitScore(MitGuide, new string(distal));

        Assert.Multiple(() =>
        {
            Assert.That(proximalScore, Is.EqualTo(14.9).Within(1e-9),
                "Max-weight PAM-proximal pos 13 (W=0.851) → 14.9; reversed W would give 68.3");
            Assert.That(distalScore, Is.EqualTo(100.0).Within(1e-9),
                "Zero-weight PAM-distal pos 0 (W=0) → 100; reversed W would give 41.7");
            Assert.That(proximalScore, Is.LessThan(distalScore),
                "PAM-proximal (seed) mismatch must be penalised far more than PAM-distal");
        });
    }

    /// <summary>
    /// MIT-006 (aggregate, boundary): empty off-target set → specificity 100.
    /// </summary>
    [Test]
    public void CalculateMitSpecificityScore_NoHits_Returns100()
    {
        double score = CrisprDesigner.CalculateMitSpecificityScore(System.Array.Empty<double>());
        Assert.That(score, Is.EqualTo(100.0).Within(1e-12));
    }

    /// <summary>
    /// MIT-007 (aggregate, known): one single-hit score of 60.5
    /// → 100/(100+60.5)*100 = 62.305295950155... (independent Python).
    /// </summary>
    [Test]
    public void CalculateMitSpecificityScore_SingleHit_KnownAggregate()
    {
        double score = CrisprDesigner.CalculateMitSpecificityScore(new[] { 60.5 });
        Assert.That(score, Is.EqualTo(62.30529595015576).Within(1e-9));
    }

    /// <summary>
    /// MIT-008 (aggregate, additivity): two hits 60.5 and 41.7
    /// → 100/(100+102.2)*100 = 49.455984174... (independent Python).
    /// </summary>
    [Test]
    public void CalculateMitSpecificityScore_TwoHits_KnownAggregate()
    {
        double score = CrisprDesigner.CalculateMitSpecificityScore(new[] { 60.5, 41.7 });
        // sum = 102.2; 100/(100+102.2)*100 = 49.45598417... (independent Python).
        Assert.That(score, Is.EqualTo(100.0 / (100.0 + 102.2) * 100.0).Within(1e-12));
        Assert.That(score, Is.EqualTo(49.455984174).Within(1e-6));
    }

    /// <summary>
    /// MIT-009: genome-scanning overload — a guide with one single-mismatch off-target.
    /// The off-target protospacer differs at exactly one position; the aggregate specificity
    /// equals 100/(100 + hitScore)*100 with the MIT/Hsu single-hit score of that one site.
    /// </summary>
    [Test]
    public void CalculateMitSpecificityScore_Genome_SingleOffTarget_MatchesFormula()
    {
        // Guide "GACGCATAAAGATGAGACGC"; off-target = same with pos 0 (PAM-distal, W=0) changed,
        // plus an NGG PAM, so the off-target protospacer differs by 1 mismatch at a zero-weight
        // position → single-hit score 100 → aggregate 100/(100+100)*100 = 50.
        string offProto = "TACGCATAAAGATGAGACGC"; // pos0 G→T (W[0]=0)
        var genome = new DnaSequence(offProto + "AGG");

        double score = CrisprDesigner.CalculateMitSpecificityScore(
            MitGuide, genome, 3, CrisprSystemType.SpCas9);

        Assert.That(score, Is.EqualTo(50.0).Within(1e-9));
    }

    /// <summary>
    /// MIT-010: genome with no off-targets → specificity 100.
    /// </summary>
    [Test]
    public void CalculateMitSpecificityScore_Genome_NoOffTargets_Returns100()
    {
        var genome = new DnaSequence(MitGuide + "AGG"); // only the exact on-target site
        double score = CrisprDesigner.CalculateMitSpecificityScore(
            MitGuide, genome, 3, CrisprSystemType.SpCas9);
        Assert.That(score, Is.EqualTo(100.0).Within(1e-12));
    }

    /// <summary>
    /// MIT-011: input-length validation — non-20-nt sequences throw.
    /// </summary>
    [Test]
    public void CalculateMitHitScore_WrongLength_Throws()
    {
        Assert.Throws<ArgumentException>(() =>
            CrisprDesigner.CalculateMitHitScore("ACGT", "ACGT"));
        Assert.Throws<ArgumentNullException>(() =>
            CrisprDesigner.CalculateMitHitScore(null!, MitGuide));
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
