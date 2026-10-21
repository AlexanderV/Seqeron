namespace Seqeron.Genomics.Tests.Unit.MolTools;

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
        string ot = string.Concat("T", MitGuide.AsSpan(1)); // pos 0: G→T
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

    #region F9 — MIT aggregate on non-20-nt systems (the reference's own length fudge)

    // Defect F9: CalculateMitSpecificityScore(string, DnaSequence, int, CrisprSystemType) filtered its
    // hits with `.Where(ot => ot.Sequence.Length == 20 && guide.Length == 20)`, so for EVERY system
    // whose guide is not 20 nt (SaCas9 21, Cas12a/AsCas12a 23, LbCas12a 24) all hits were dropped and
    // the aggregate returned a perfect 100 despite real off-targets.
    //
    // The reference handles the lengths itself and picks the model per system (annotateOfftargets,
    // crispor.py L1934-1946):
    //   SaCas9 PAM  -> calcSaHitScore = 100 x SaCas9Scorer.calcScore (Tycko 2018 trained model)
    //   otherwise   -> calcHitScore, whose own length fudge is len > 20 -> last 20 nt ("S. aureus
    //                  requires 21bp long guides. We fudge by using only last 20bp") and len == 19 ->
    //                  left-pad both with "A" (crispor.py L2047-2068)
    // F14 routes SaCas9 through the Tycko model, as the reference does, so the SaCas9 expectations
    // below are its values (independent Python run of predictSingle.py's semantics against the
    // downloaded rstan_mult_21_coeffs.txt):
    //   calcSaHitScore("ACGTACGTACGTACGTACGTA", "ACGTACGTACGTACGTACGTC") = 4.2613097364785695
    //       -> calcMitGuideScore 95.91285612347565
    //   calcSaHitScore("ACGTACGTACGTACGTACGTA", "ACGTCCGTACGTACGTACGTC") = 2.5134348817970387
    //       -> calcMitGuideScore 97.54818977171612
    //   calcSaHitScore("ACGTACGTACGTACGTACGTA", "TCGTACGTACGTACGTACGTA") = 100.0 (the 5'-most position
    //       contributes nothing in either model) -> calcMitGuideScore 50.0

    private const string SaGuide21 = "ACGTACGTACGTACGTACGTA";     // SaCas9: 21-nt guide, NNGRRT PAM
    private const string SaPam = "AAGAAT";                        // NNGRRT

    /// <summary>
    /// F9-001 (regression): a 21-nt SaCas9 guide with ONE real off-target must NOT score 100.
    /// The single off-target protospacer differs at index 20 (the PAM-proximal base), so the
    /// reference single-hit score is calcSaHitScore = 4.2613097364785695 and the aggregate is
    /// 100/(100+4.2613097364785695)*100 = 95.91285612347565. Before the fix this returned exactly
    /// 100.0 for EVERY SaCas9 guide, however many off-targets the genome held.
    /// </summary>
    [Test]
    public void CalculateMitSpecificityScore_SaCas9_21ntGuide_RealOffTarget_IsNot100()
    {
        var genome = new DnaSequence("ACGTACGTACGTACGTACGTC" + SaPam);

        var offTargets = CrisprDesigner.FindOffTargets(SaGuide21, genome, 3, CrisprSystemType.SaCas9).ToList();
        double score = CrisprDesigner.CalculateMitSpecificityScore(
            SaGuide21, genome, 3, CrisprSystemType.SaCas9);

        Assert.Multiple(() =>
        {
            // The hit set is asserted independently of the score, so the expected aggregate below
            // belongs to exactly this one off-target.
            Assert.That(offTargets, Has.Count.EqualTo(1));
            Assert.That(offTargets[0].Sequence, Is.EqualTo("ACGTACGTACGTACGTACGTC"));
            Assert.That(offTargets[0].Mismatches, Is.EqualTo(1));
            Assert.That(offTargets[0].MismatchPositions, Is.EqualTo(new[] { 20 }));

            Assert.That(score, Is.Not.EqualTo(100.0), "F9: a real off-target must lower the aggregate");
            Assert.That(score, Is.EqualTo(95.91285612347565).Within(1e-9));
        });
    }

    /// <summary>
    /// F9-002: the same 21-nt SaCas9 path with TWO mismatches — the Tycko model multiplies the two
    /// positions' coefficients. Reference single-hit score 2.5134348817970387
    /// → aggregate 97.54818977171612.
    /// </summary>
    [Test]
    public void CalculateMitSpecificityScore_SaCas9_21ntGuide_TwoMismatches_MatchesReference()
    {
        var genome = new DnaSequence("ACGTCCGTACGTACGTACGTC" + SaPam);

        var offTargets = CrisprDesigner.FindOffTargets(SaGuide21, genome, 3, CrisprSystemType.SaCas9).ToList();
        double score = CrisprDesigner.CalculateMitSpecificityScore(
            SaGuide21, genome, 3, CrisprSystemType.SaCas9);

        Assert.Multiple(() =>
        {
            Assert.That(offTargets, Has.Count.EqualTo(1));
            Assert.That(offTargets[0].MismatchPositions, Is.EqualTo(new[] { 4, 20 }));
            Assert.That(score, Is.EqualTo(97.54818977171612).Within(1e-9));
        });
    }

    /// <summary>
    /// F9-003 (locked from the reference): a mismatch at the 5'-most position contributes nothing —
    /// the Tycko model skips index 0 outright ("Mismatches at the 5' position are not truly determined
    /// by the model") and calcHitScore's last-20 window does not even see it. The hit therefore scores
    /// a perfect 100 and the aggregate is 100/(100+100)*100 = 50, not 100. This distinguishes
    /// "scored, and it is a strong hit" from F9's "never scored at all".
    /// </summary>
    [Test]
    public void CalculateMitSpecificityScore_SaCas9_MismatchAtFivePrimeEnd_ScoredAsPerfectHit()
    {
        var genome = new DnaSequence("TCGTACGTACGTACGTACGTA" + SaPam);

        var offTargets = CrisprDesigner.FindOffTargets(SaGuide21, genome, 3, CrisprSystemType.SaCas9).ToList();
        double score = CrisprDesigner.CalculateMitSpecificityScore(
            SaGuide21, genome, 3, CrisprSystemType.SaCas9);

        Assert.Multiple(() =>
        {
            Assert.That(offTargets, Has.Count.EqualTo(1));
            Assert.That(offTargets[0].MismatchPositions, Is.EqualTo(new[] { 0 }));
            Assert.That(score, Is.EqualTo(50.0).Within(1e-9));
        });
    }

    /// <summary>
    /// F9-004: Cas12a-class (PAM-before-target) systems must THROW rather than report a specificity.
    /// The reference reports no MIT/CFD off-target score for Cpf1 PAMs at all: in annotateOfftargets
    /// (crispor.py L1939-1942) every hit gets mitScore = cfdScore = 0.0 when pamIsCpf1(pam), and the
    /// guide aggregates are the sentinels guideScore = guideCfdScore = -1 (L1986-1988). The test first
    /// proves there ARE real off-targets, so the throw cannot be mistaken for "nothing to score".
    /// </summary>
    [Test]
    public void CalculateMitSpecificityScore_Cas12aClass_WithRealOffTargets_Throws()
    {
        const string guide23 = "ACGTACGTACGTACGTACGTACG";
        var genome = new DnaSequence("TTTA" + "ACGTACGTACGTACGTACGTACC"); // TTTV PAM, 1 mismatch at index 22

        var offTargets = CrisprDesigner.FindOffTargets(guide23, genome, 3, CrisprSystemType.Cas12a).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(offTargets, Is.Not.Empty, "the genome must contain a real Cas12a off-target");
            Assert.Throws<ArgumentException>(() => CrisprDesigner.CalculateMitSpecificityScore(
                guide23, genome, 3, CrisprSystemType.Cas12a));
            Assert.Throws<ArgumentException>(() => CrisprDesigner.CalculateCfdSpecificityScore(
                guide23, genome, 3, CrisprSystemType.Cas12a));
        });
    }

    /// <summary>
    /// F9-005: every PAM-before-target system is refused, including CasX (TTCN), which the
    /// reference's pamIsCpf1 list (crispor.py L889-891) does NOT contain — CRISPOR would route it
    /// through calcHitScore, but its own PAM stripping there takes the 3' end, the wrong end for a
    /// PAM-before system, and the Hsu weight vector is indexed PAM-proximal-last. Refusing is the
    /// deliberate divergence (recorded in docs/Evidence/CRISPR-OFF-001-Evidence.md).
    /// </summary>
    [TestCase(CrisprSystemType.Cas12a, 23)]
    [TestCase(CrisprSystemType.AsCas12a, 23)]
    [TestCase(CrisprSystemType.LbCas12a, 24)]
    [TestCase(CrisprSystemType.CasX, 20)]
    public void CalculateMitSpecificityScore_PamBeforeTargetSystems_AlwaysThrow(
        CrisprSystemType system, int guideLength)
    {
        string guide = new string('A', guideLength);
        var genome = new DnaSequence("TTTAACGTACGTACGTACGTACGTACGT");

        Assert.Throws<ArgumentException>(() =>
            CrisprDesigner.CalculateMitSpecificityScore(guide, genome, 3, system));
    }

    /// <summary>
    /// F9-006 (honesty guard): the unsourced heuristic and the published MIT aggregate must be
    /// visibly different numbers on the same input, so neither can be mistaken for the other.
    /// Heuristic: one PAM-distal mismatch → 100 − 2 = 98. MIT: single-hit 100 → 50.
    /// </summary>
    [Test]
    public void CalculateSpecificityScore_HeuristicAndMitAggregate_AreDistinct()
    {
        const string guide = "GACGCATAAAGATGAGACGC";
        var genome = new DnaSequence("TACGCATAAAGATGAGACGC" + "AGG");

        double heuristic = CrisprDesigner.CalculateSpecificityScore(guide, genome, CrisprSystemType.SpCas9);
        double mit = CrisprDesigner.CalculateMitSpecificityScore(guide, genome, 4, CrisprSystemType.SpCas9);

        Assert.Multiple(() =>
        {
            Assert.That(heuristic, Is.EqualTo(98.0).Within(1e-12));
            Assert.That(mit, Is.EqualTo(50.0).Within(1e-9));
        });
    }

    /// <summary>
    /// F39-001 (B08 WP2 / D13): the unsourced heuristic's primary name is
    /// <c>CalculateHeuristicSpecificityScore</c> — the word "heuristic" is in the signature, so the
    /// declaration is structural and not doc-only — and the historical
    /// <c>CalculateSpecificityScore</c> is a pure forwarding alias kept for backward compatibility.
    /// This test is what keeps the alias from drifting away from the method it forwards to; it is
    /// also why the alias is NOT <c>[Obsolete]</c> (the repo builds warnings-as-errors and the MCP
    /// tool name <c>crispr_specificity_score</c> stays unchanged, a published catalog name).
    /// Values are the same ones F9-006 locks (heuristic: one PAM-distal mismatch → 100 − 2 = 98).
    /// </summary>
    [Test]
    public void CalculateHeuristicSpecificityScore_LegacyAlias_ForwardsIdentically()
    {
        const string guide = "GACGCATAAAGATGAGACGC";
        var genome = new DnaSequence("TACGCATAAAGATGAGACGC" + "AGG");

        double primary = CrisprDesigner.CalculateHeuristicSpecificityScore(guide, genome, CrisprSystemType.SpCas9);
        double alias = CrisprDesigner.CalculateSpecificityScore(guide, genome, CrisprSystemType.SpCas9);

        Assert.Multiple(() =>
        {
            Assert.That(primary, Is.EqualTo(98.0).Within(1e-12));
            Assert.That(alias, Is.EqualTo(primary), "the alias must forward with no behaviour of its own");
        });
    }

    #endregion

    #region F14 — SaCas9 off-target model (Tycko et al. 2018)

    // Source: Tycko J, Wainberg M, Marinov GK, et al., "Pairwise library screen systematically
    // interrogates Staphylococcus aureus Cas9 specificity in human cells", Nat Commun 9:2962 (2018),
    // doi:10.1038/s41467-018-05391-2.
    // Reference implementation: CRISPOR calcSaHitScore (crispor.py L2025-2043) -> SaCas9Scorer in
    // bin/src/pairwise-library-screen/predictSingle.py + models/rstan_mult_21_coeffs.txt, both
    // downloaded this session from raw.githubusercontent.com/maximilianh/crisporWebsite/master/.
    // The three expected values below are the reference script's OWN asserted values (its __main__
    // block), reproduced by an independent Python run against the downloaded coefficient table:
    //   calcScore(identical)                                                      = 1.0
    //   calcScore("GCAACCACAAACCCACGAGGG", "ACAAACACATACCCACAAGGA") = 0.041293256150042444
    //   calcScore("GGGTGAGTGAGTGTGTGCGTG", "GAGCGAGTGGGTGTGTGCGTG") = 0.13441983202750601
    // calcSaHitScore returns 100 x calcScore ("to be compatible with the MIT score, has to be in the
    // range 0-100 for the MIT aggregate guide specificity score").

    /// <summary>SA-001: a perfect match scores exactly 100 (the reference's own first assert).</summary>
    [Test]
    public void CalculateSaCas9HitScore_PerfectMatch_Returns100()
    {
        Assert.That(
            CrisprDesigner.CalculateSaCas9HitScore("GCAACCACAAACCCACGAGGG", "GCAACCACAAACCCACGAGGG"),
            Is.EqualTo(100.0).Within(1e-9));
    }

    /// <summary>SA-002 / SA-003: the reference script's two published multi-mismatch oracles.</summary>
    [TestCase("GCAACCACAAACCCACGAGGG", "ACAAACACATACCCACAAGGA", 4.1293256150042446)]
    [TestCase("GGGTGAGTGAGTGTGTGCGTG", "GAGCGAGTGGGTGTGTGCGTG", 13.4419832027506)]
    public void CalculateSaCas9HitScore_ReferenceOracles_Match(string guide, string offTarget, double expected)
    {
        Assert.That(CrisprDesigner.CalculateSaCas9HitScore(guide, offTarget),
            Is.EqualTo(expected).Within(1e-9));
    }

    /// <summary>
    /// SA-004 (coefficient ORIENTATION guard): the shipped table indexes Position from the 3'/
    /// PAM-proximal end and SaCas9Scorer remaps it to a 1-based index from the 5' end
    /// (guideLenDict[guide_len - Position + 1]). The SAME mismatch type "dT:rC" (guide C, off-target A)
    /// therefore carries 0.546734106002513 at 5'-adjacent index 1 and 0.105693733152788 at
    /// PAM-proximal index 20 → 54.6734106002513 vs 10.569373315278801. A table transcribed without the
    /// remap would swap these two values, so asserting both pins the orientation.
    /// </summary>
    [Test]
    public void CalculateSaCas9HitScore_CoefficientOrientation_NotReversed()
    {
        const string guide = "ACGTACGTACGTACGTACGTC";

        Assert.Multiple(() =>
        {
            Assert.That(CrisprDesigner.CalculateSaCas9HitScore(guide, "AAGTACGTACGTACGTACGTC"),
                Is.EqualTo(54.6734106002513).Within(1e-9), "dT:rC at 5'-adjacent index 1");
            Assert.That(CrisprDesigner.CalculateSaCas9HitScore(guide, "ACGTACGTACGTACGTACGTA"),
                Is.EqualTo(10.569373315278801).Within(1e-9), "dT:rC at PAM-proximal index 20");
        });
    }

    /// <summary>
    /// SA-005: "Mismatches at the 5' position are not truly determined by the model" — index 0 is
    /// skipped outright, so a guide differing only there scores a perfect 100.
    /// </summary>
    [Test]
    public void CalculateSaCas9HitScore_MismatchAtFivePrimeMostPosition_Ignored()
    {
        Assert.That(
            CrisprDesigner.CalculateSaCas9HitScore("ACGTACGTACGTACGTACGTC", "CCGTACGTACGTACGTACGTC"),
            Is.EqualTo(100.0).Within(1e-9));
    }

    /// <summary>
    /// SA-006: only the 21-nt coefficient table is transcribed (21 nt is this library's SaCas9 guide
    /// length), and the model supports neither gaps/bulges nor non-A/C/G/T bases — all refused rather
    /// than approximated.
    /// </summary>
    [Test]
    public void CalculateSaCas9HitScore_InvalidInputs_Throw()
    {
        const string guide = "ACGTACGTACGTACGTACGTC";

        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentNullException>(() => CrisprDesigner.CalculateSaCas9HitScore(null!, guide));
            Assert.Throws<ArgumentNullException>(() => CrisprDesigner.CalculateSaCas9HitScore(guide, null!));
            Assert.Throws<ArgumentNullException>(() => CrisprDesigner.CalculateSaCas9HitScore("", guide));
            // 20 nt and 22 nt tables exist in the reference but are not transcribed here.
            Assert.Throws<ArgumentException>(() => CrisprDesigner.CalculateSaCas9HitScore(
                "ACGTACGTACGTACGTACGT", "ACGTACGTACGTACGTACGT"));
            Assert.Throws<ArgumentException>(() => CrisprDesigner.CalculateSaCas9HitScore(guide, "ACGTACGTACGTACGTACGT"));
            // Gaps/bulges: "This model does not support gaps/bulges."
            Assert.Throws<ArgumentException>(() => CrisprDesigner.CalculateSaCas9HitScore(
                guide, "A-GTACGTACGTACGTACGTC"));
            Assert.Throws<ArgumentException>(() => CrisprDesigner.CalculateSaCas9HitScore(
                guide, "ANGTACGTACGTACGTACGTC"));
        });
    }

    /// <summary>SA-007: lower-case input is folded, like every other scoring entry point here.</summary>
    [Test]
    public void CalculateSaCas9HitScore_LowerCaseInput_TreatedAsUpperCase()
    {
        Assert.That(
            CrisprDesigner.CalculateSaCas9HitScore("gcaaccacaaacccacgaggg", "acaaacacatacccacaagga"),
            Is.EqualTo(4.1293256150042446).Within(1e-9));
    }

    #endregion

    #region Mismatch-count / position / penalty consistency and enumeration order

    /// <summary>
    /// Locks INV-03 plus the penalty's agreement with it for a mixed-case input: Mismatches,
    /// MismatchPositions.Count and OffTargetScore are all derived from ONE case-folding rule, so they
    /// cannot disagree. (DnaSequence itself upper-cases and validates its input, so the genome side is
    /// already normalised — see DnaSequence.cs ValidateSequence; the guide side is folded in the scan.)
    /// </summary>
    [Test]
    public void FindOffTargets_MixedCaseInput_CountPositionsAndPenaltyAgree()
    {
        var offTargets = CrisprDesigner.FindOffTargets(
            "gacgcataaagatgagacgc",                               // lower-case guide
            new DnaSequence("tacgcataaagatgagacgcagg"),            // lower-case genome
            3, CrisprSystemType.SpCas9).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(offTargets, Has.Count.EqualTo(1));
            Assert.That(offTargets[0].Mismatches, Is.EqualTo(1));
            Assert.That(offTargets[0].MismatchPositions.Count, Is.EqualTo(offTargets[0].Mismatches));
            Assert.That(offTargets[0].MismatchPositions, Is.EqualTo(new[] { 0 }));
            // index 0 is PAM-distal for a 20-nt Cas9 guide (seed = last 12 nt) → penalty 2, not 5.
            Assert.That(offTargets[0].OffTargetScore, Is.EqualTo(2.0).Within(1e-12));
        });
    }

    /// <summary>
    /// A DnaSequence can never carry a non-A/C/G/T base (DnaSequence.ValidateSequence →
    /// SequenceExtensions.IndexOfInvalidDna), which is why the reference's "[^ATCG] → cfdScore = -1"
    /// exclusion inside the CFD genome overload is defensive and unreachable. Locked here so the
    /// Evidence doc's reachability claim stays true if DnaSequence's alphabet ever widens.
    /// </summary>
    [Test]
    public void DnaSequence_RejectsAmbiguityCodes_SoNonAcgtOffTargetsCannotReachTheScan()
    {
        Assert.Throws<ArgumentException>(() => new DnaSequence("ACGTACGTNCGTACGTACGTAGG"));
    }

    /// <summary>
    /// Determinism has a documented shape: all forward-strand hits in ascending position order, then
    /// all reverse-strand hits. Locked with a genome carrying one off-target per strand.
    /// </summary>
    [Test]
    public void FindOffTargets_EnumerationOrder_ForwardSitesThenReverseSites()
    {
        const string guide = "GACGCATAAAGATGAGACGC";
        // forward off-target (1 mm at index 0) + NGG, then the reverse-complement of another
        // off-target + NGG, so one hit comes from each strand.
        var genome = new DnaSequence("TACGCATAAAGATGAGACGCAGG" + "CCTGCGTCACATCTTTATGCGTA");

        var offTargets = CrisprDesigner.FindOffTargets(guide, genome, 4, CrisprSystemType.SpCas9).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(offTargets.Select(ot => ot.IsForwardStrand),
                Is.Ordered.Descending, "forward-strand hits must precede reverse-strand hits");
            Assert.That(offTargets.Where(ot => ot.IsForwardStrand).Select(ot => ot.Position),
                Is.Ordered.Ascending);
            Assert.That(offTargets.Any(ot => ot.IsForwardStrand), Is.True);
            Assert.That(offTargets.Any(ot => !ot.IsForwardStrand), Is.True);
        });
    }

    /// <summary>
    /// Dedup lock: <c>CountMismatches</c> delegates equal-length inputs to the canonical
    /// <see cref="SequenceExtensions.HammingDistance(ReadOnlySpan{char}, ReadOnlySpan{char})"/>, so the
    /// reported <c>Mismatches</c> must equal that primitive on the very same pair — including for
    /// mixed-case input, since the canonical folds ASCII case and the replaced private loop did not.
    /// </summary>
    [Test]
    public void FindOffTargets_Mismatches_EqualCanonicalHammingDistance()
    {
        const string guide = "gacgcataaagatgagacgc";
        var genome = new DnaSequence("tacgcataaagatgagacgcagg" + "CCTGCGTCACATCTTTATGCGTA");

        var offTargets = CrisprDesigner.FindOffTargets(guide, genome, 4, CrisprSystemType.SpCas9).ToList();

        Assert.That(offTargets, Is.Not.Empty);
        Assert.Multiple(() =>
        {
            foreach (var ot in offTargets)
            {
                Assert.That(
                    ot.Mismatches,
                    Is.EqualTo(guide.AsSpan().HammingDistance(ot.Sequence.AsSpan())),
                    $"reported Mismatches must equal the canonical HammingDistance for {ot.Sequence}");
                Assert.That(ot.MismatchPositions.Count, Is.EqualTo(ot.Mismatches));
            }
        });
    }

    /// <summary>
    /// Dedup lock, the other direction: <c>CfdComplement</c> is deliberately NOT the canonical
    /// <see cref="SequenceExtensions.GetComplementBase(char)"/>. The canonical is IUPAC-complete and
    /// maps <c>U → A</c>, whereas CFD must refuse anything outside A/C/G/T — the reference
    /// <c>calcCfdScore</c> screens both sequences with <c>re.search('[^ATCG]')</c> and returns the
    /// sentinel <c>-1</c> (crispor.py L2175-2185). This pins both halves of that argument, so the
    /// non-delegation cannot be "cleaned up" into a silently wrong score.
    /// </summary>
    [Test]
    public void CalculateCfdScore_NonAcgtBases_RefusedAlthoughCanonicalComplementAcceptsThem()
    {
        const string guide = "GACGCATAAAGATGAGACGC";

        Assert.Multiple(() =>
        {
            // The canonical complement would happily map these — which is exactly why CFD cannot use it.
            Assert.That(SequenceExtensions.GetComplementBase('U'), Is.EqualTo('A'));
            Assert.That(SequenceExtensions.GetComplementBase('N'), Is.EqualTo('N'));
            Assert.That(SequenceExtensions.GetComplementBase('R'), Is.EqualTo('Y'));

            // CFD refuses them instead of scoring them.
            Assert.Throws<ArgumentException>(() => CrisprDesigner.CalculateCfdScore(
                guide, "UACGCATAAAGATGAGACGC", "AGG"));
            Assert.Throws<ArgumentException>(() => CrisprDesigner.CalculateCfdScore(
                guide, "NACGCATAAAGATGAGACGC", "AGG"));
            Assert.Throws<ArgumentException>(() => CrisprDesigner.CalculateCfdScore(
                "RACGCATAAAGATGAGACGC", "TACGCATAAAGATGAGACGC", "AGG"));
            Assert.Throws<ArgumentException>(() => CrisprDesigner.CalculateCfdScore(
                guide, "TACGCATAAAGATGAGACGC", "ANG"));
        });
    }

    /// <summary>
    /// Recorded divergence (Evidence 10): every exact match is dropped genome-wide, because this API is
    /// given a bare reference sequence and no on-target locus. CRISPOR's skip
    /// (<c>annotateOfftargets</c>, crispor.py L1920-1929) is scoped to the caller-supplied on-target
    /// interval — or, with none, to the first 0-mismatch hit it meets — and scores any FURTHER perfect
    /// match as a full <c>calcHitScore = 100</c> off-target, which would give an aggregate of 50 here.
    /// This pins the documented behaviour so the divergence stays visible rather than drifting.
    /// </summary>
    [Test]
    public void FindOffTargets_TwoPerfectMatches_BothDroppedGenomeWide()
    {
        const string guide = "GACGCATAAAGATGAGACGC";
        // The guide's own protospacer + NGG, twice: two 0-mismatch hits and nothing else.
        var genome = new DnaSequence(guide + "AGG" + "TTTT" + guide + "AGG");

        var offTargets = CrisprDesigner.FindOffTargets(guide, genome, 4, CrisprSystemType.SpCas9).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(offTargets, Is.Empty, "both exact matches are treated as on-targets");
            Assert.That(CrisprDesigner.CalculateMitSpecificityScore(guide, genome, 4), Is.EqualTo(100.0));
            Assert.That(CrisprDesigner.CalculateCfdSpecificityScore(guide, genome, 4), Is.EqualTo(100.0));
            // CRISPOR would keep the second perfect match and report 100/(100+100)*100 = 50.
        });
    }

    #endregion
    #region F50 — CRISPOR's alternative-PAM x0.2 MIT downweight (crispor.py L1947-1958)

    // Reference rule (annotateOfftargets, crispor.py L1947-1958), applied AFTER the per-hit MIT score:
    //     mitScore = calcHitScore(guideNoPam, otSeqNoPam)
    //     if pam=="NGG" and otSeq[-2:]!="GG":
    //         mitScore = mitScore * 0.2
    //     cfdScore = calcCfdScore(guideSeq, otSeq)        # NOT downweighted
    // and the guide-level value stays calcMitGuideScore(sum(mitOtScores)).
    // Alternative PAMs come from the reference's own table (offtargetPams, L243-253):
    // NGG -> [NAG, NGA], NGA -> [NGG], NNGRRT -> [NNGRRN], TTTV -> [TTTN].
    //
    // Expected values below come from a verbatim Python copy of crispor.py's calcHitScore
    // (L2045-2096), calcMitGuideScore (L2098-2104, unrounded) and that downweight, run this session:
    //   calcHitScore(guide, AltPamOffTarget)  = 60.5                 (1 mismatch at index 5)
    //   calcHitScore(guide, ProximalOffTarget)= 31.499999999999993   (1 mismatch at index 18)
    //   calcHitScore(guide, TwoMmOffTarget)   = 2.6227599999999995   (mismatches at 2 and 17)
    //   sum[60.5, 31.4999...]                 -> 52.083333333333336  (both main PAM)
    //   sum[60.5, 0.2 x 31.4999...]           -> 59.95203836930455   (second hit on an alt PAM)
    //   sum[60.5, 0.2 x 31.4999..., 0.2 x 2.62276] -> 59.7640924805823

    private const string AltPamGuide = "GACGCATAAAGATGAGACGC";
    private const string AltPamOffTarget = "GACGCTTAAAGATGAGACGC";      // 1 mismatch, index 5
    private const string ProximalOffTarget = "GACGCATAAAGATGAGACTC";    // 1 mismatch, index 18
    private const string TwoMmOffTarget = "GATGCATAAAGATGAGAAGC";       // mismatches at 2 and 17

    /// <summary>
    /// A genome holding exactly two off-targets of <see cref="AltPamGuide"/>: one on the main NGG PAM
    /// and one on an alternative NAG PAM. It contains no "CC", so the reverse strand carries no NGG/NAG
    /// site at all, and the hit sets are asserted in the tests rather than assumed.
    /// </summary>
    private static DnaSequence AltPamGenome() =>
        new(AltPamOffTarget + "AGG" + "TTTT" + ProximalOffTarget + "TAG");

    private static PamAnnotatedOffTarget Hit(string protospacer, string pam) =>
        new(new OffTargetSite(0, protospacer, 1, new[] { 0 }, true, 0.0), pam);

    /// <summary>
    /// F50-001: the same two hits aggregate to 52.083333333333336 when both carry a main NGG PAM and to
    /// 59.95203836930455 when the second carries an NAG PAM — i.e. the second hit's contribution is
    /// exactly one fifth. Reference values from the Python copy of calcHitScore/calcMitGuideScore.
    /// </summary>
    [Test]
    public void CalculateMitSpecificityScore_Annotated_AlternativePamHit_ContributesOneFifth()
    {
        var bothMain = new[] { Hit(AltPamOffTarget, "AGG"), Hit(ProximalOffTarget, "CGG") };
        var secondIsAlt = new[] { Hit(AltPamOffTarget, "AGG"), Hit(ProximalOffTarget, "TAG") };

        Assert.Multiple(() =>
        {
            Assert.That(CrisprDesigner.CalculateMitSpecificityScore(AltPamGuide, bothMain),
                Is.EqualTo(52.083333333333336).Within(1e-9));
            Assert.That(CrisprDesigner.CalculateMitSpecificityScore(AltPamGuide, secondIsAlt),
                Is.EqualTo(59.95203836930455).Within(1e-9));
        });
    }

    /// <summary>
    /// F50-002: the reference tests the hit's last two PAM bases (<c>otSeq[-2:]!="GG"</c>), not
    /// membership of the alternative-PAM table, so NGA is downweighted exactly like NAG, and an
    /// arbitrary non-*GG PAM is too. A *GG PAM other than the canonical AGG is NOT downweighted.
    /// </summary>
    [Test]
    public void CalculateMitSpecificityScore_Annotated_DownweightIsOnTheLastTwoPamBases()
    {
        double nag = CrisprDesigner.CalculateMitSpecificityScore(
            AltPamGuide, new[] { Hit(AltPamOffTarget, "AGG"), Hit(ProximalOffTarget, "TAG") });
        double nga = CrisprDesigner.CalculateMitSpecificityScore(
            AltPamGuide, new[] { Hit(AltPamOffTarget, "AGG"), Hit(ProximalOffTarget, "CGA") });
        double tgg = CrisprDesigner.CalculateMitSpecificityScore(
            AltPamGuide, new[] { Hit(AltPamOffTarget, "GGG"), Hit(ProximalOffTarget, "TGG") });

        Assert.Multiple(() =>
        {
            Assert.That(nga, Is.EqualTo(nag).Within(1e-12), "NGA and NAG are both non-*GG PAMs");
            Assert.That(nag, Is.EqualTo(59.95203836930455).Within(1e-9));
            Assert.That(tgg, Is.EqualTo(52.083333333333336).Within(1e-9), "every *GG PAM is a main-PAM hit");
        });
    }

    /// <summary>
    /// F50-003: three hits, two of them on alternative PAMs (one NAG, one hand-supplied NGA — a PAM no
    /// <see cref="CrisprSystemType"/> in this library scans) aggregate to 59.7640924805823.
    /// </summary>
    [Test]
    public void CalculateMitSpecificityScore_Annotated_ThreeHitsTwoAlternative_MatchesReference()
    {
        var hits = new[]
        {
            Hit(AltPamOffTarget, "AGG"),
            Hit(ProximalOffTarget, "TAG"),
            Hit(TwoMmOffTarget, "CGA"),
        };

        Assert.That(CrisprDesigner.CalculateMitSpecificityScore(AltPamGuide, hits),
            Is.EqualTo(59.7640924805823).Within(1e-9));
    }

    /// <summary>
    /// F50-004: the reference gates the downweight on the SEARCHED PAM pattern being literally "NGG"
    /// (<c>if pam=="NGG" and ...</c>), so an SpCas9-NAG guide gets no downweight on any hit — not even
    /// on a hit whose own PAM is NAG. Both hit sets must therefore give the undiscounted
    /// 52.083333333333336.
    /// </summary>
    [Test]
    public void CalculateMitSpecificityScore_Annotated_NagSystem_NeverDownweights()
    {
        var allNag = new[] { Hit(AltPamOffTarget, "TAG"), Hit(ProximalOffTarget, "CAG") };
        var mixed = new[] { Hit(AltPamOffTarget, "TAG"), Hit(ProximalOffTarget, "AGG") };

        Assert.Multiple(() =>
        {
            Assert.That(
                CrisprDesigner.CalculateMitSpecificityScore(AltPamGuide, allNag, CrisprSystemType.SpCas9_NAG),
                Is.EqualTo(52.083333333333336).Within(1e-9));
            Assert.That(
                CrisprDesigner.CalculateMitSpecificityScore(AltPamGuide, mixed, CrisprSystemType.SpCas9_NAG),
                Is.EqualTo(52.083333333333336).Within(1e-9));
        });
    }

    /// <summary>
    /// F50-005: the annotated overload keeps F14's per-system model routing — a SaCas9 guide is scored
    /// with the Tycko 2018 model, and SaCas9's main PAM is not NGG, so nothing is downweighted. The
    /// expected value is the one F9-001 locks for the same single off-target (calcSaHitScore
    /// 4.2613097364785695 -> 95.91285612347565).
    /// </summary>
    [Test]
    public void CalculateMitSpecificityScore_Annotated_SaCas9_UsesTyckoModelAndNoDownweight()
    {
        var hits = new[]
        {
            new PamAnnotatedOffTarget(
                new OffTargetSite(0, "ACGTACGTACGTACGTACGTC", 1, new[] { 20 }, true, 0.0), SaPam),
        };

        Assert.That(
            CrisprDesigner.CalculateMitSpecificityScore(SaGuide21, hits, CrisprSystemType.SaCas9),
            Is.EqualTo(95.91285612347565).Within(1e-9));
    }

    /// <summary>
    /// F50-006: the alternative-PAM table is the reference's own <c>offtargetPams</c> entry for each
    /// system's main PAM (crispor.py L243-253); a PAM the table does not list yields an empty list.
    /// </summary>
    [Test]
    public void GetAlternativeOffTargetPams_MatchesReferenceTable()
    {
        Assert.Multiple(() =>
        {
            Assert.That(CrisprDesigner.GetAlternativeOffTargetPams(CrisprSystemType.SpCas9),
                Is.EqualTo(new[] { "NAG", "NGA" }));
            Assert.That(CrisprDesigner.GetAlternativeOffTargetPams(CrisprSystemType.SaCas9),
                Is.EqualTo(new[] { "NNGRRN" }));
            Assert.That(CrisprDesigner.GetAlternativeOffTargetPams(CrisprSystemType.Cas12a),
                Is.EqualTo(new[] { "TTTN" }));
            // NAG and TTCN are not keys of the reference's table (only values / absent).
            Assert.That(CrisprDesigner.GetAlternativeOffTargetPams(CrisprSystemType.SpCas9_NAG), Is.Empty);
            Assert.That(CrisprDesigner.GetAlternativeOffTargetPams(CrisprSystemType.CasX), Is.Empty);
        });
    }

    /// <summary>
    /// F50-007: the genome-level alternative-PAM overload folds the NAG scan into the NGG scan with the
    /// 0.2 weight. The two scans' hit sets are asserted first, so the expected aggregate belongs to
    /// exactly those two hits: 60.5 (NGG) + 0.2 x 31.4999... (NAG) -> 59.95203836930455. The
    /// single-system overload is unchanged on the same genome: 100/(100+60.5)*100 = 62.30529595015576.
    /// </summary>
    [Test]
    public void CalculateMitSpecificityScoreWithAlternativePams_FoldsInTheNagScan()
    {
        var genome = AltPamGenome();

        var mainHits = CrisprDesigner.FindOffTargetsWithPam(AltPamGuide, genome, 4).ToList();
        var nagHits = CrisprDesigner
            .FindOffTargetsWithPam(AltPamGuide, genome, 4, CrisprSystemType.SpCas9_NAG).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(mainHits, Has.Count.EqualTo(1));
            Assert.That(mainHits[0].Site.Sequence, Is.EqualTo(AltPamOffTarget));
            Assert.That(mainHits[0].Pam, Is.EqualTo("AGG"));
            Assert.That(nagHits, Has.Count.EqualTo(1));
            Assert.That(nagHits[0].Site.Sequence, Is.EqualTo(ProximalOffTarget));
            Assert.That(nagHits[0].Pam, Is.EqualTo("TAG"));

            Assert.That(
                CrisprDesigner.CalculateMitSpecificityScoreWithAlternativePams(AltPamGuide, genome, 4),
                Is.EqualTo(59.95203836930455).Within(1e-9));
            Assert.That(
                CrisprDesigner.CalculateMitSpecificityScore(AltPamGuide, genome, 4),
                Is.EqualTo(62.30529595015576).Within(1e-9),
                "the single-system overload must keep its behaviour");
        });
    }

    /// <summary>
    /// F50-008: with no off-target at all both overloads still return 100, and for a system whose
    /// alternative PAM no <see cref="CrisprSystemType"/> provides (SaCas9 -> NNGRRN) the alternative-PAM
    /// overload reduces to the single-system one.
    /// </summary>
    [Test]
    public void CalculateMitSpecificityScoreWithAlternativePams_NoScannableAlternative_EqualsSingleSystem()
    {
        var saGenome = new DnaSequence("ACGTACGTACGTACGTACGTC" + SaPam);
        var empty = new DnaSequence(string.Concat(Enumerable.Repeat("T", 60)));

        Assert.Multiple(() =>
        {
            Assert.That(
                CrisprDesigner.CalculateMitSpecificityScoreWithAlternativePams(
                    SaGuide21, saGenome, 3, CrisprSystemType.SaCas9),
                Is.EqualTo(CrisprDesigner.CalculateMitSpecificityScore(
                    SaGuide21, saGenome, 3, CrisprSystemType.SaCas9)).Within(1e-12));
            Assert.That(
                CrisprDesigner.CalculateMitSpecificityScoreWithAlternativePams(AltPamGuide, empty, 4),
                Is.EqualTo(100.0));
        });
    }

    /// <summary>
    /// F50-009: <see cref="CrisprDesigner.FindOffTargetsWithPam"/> reports exactly the sites
    /// <see cref="CrisprDesigner.FindOffTargets"/> reports, in the same order, and repeats its guards.
    /// </summary>
    [Test]
    public void FindOffTargetsWithPam_MatchesFindOffTargets_AndValidatesItsArguments()
    {
        var genome = AltPamGenome();
        var plain = CrisprDesigner.FindOffTargets(AltPamGuide, genome, 4).ToList();
        var annotated = CrisprDesigner.FindOffTargetsWithPam(AltPamGuide, genome, 4).ToList();

        Assert.Multiple(() =>
        {
            // Compared field by field: OffTargetSite's record equality is reference equality on its
            // MismatchPositions list, so two equal-valued sites from two scans are never Equal.
            Assert.That(annotated, Has.Count.EqualTo(plain.Count));
            Assert.That(annotated.Select(h => h.Site.Position), Is.EqualTo(plain.Select(s => s.Position)));
            Assert.That(annotated.Select(h => h.Site.Sequence), Is.EqualTo(plain.Select(s => s.Sequence)));
            Assert.That(annotated.Select(h => h.Site.Mismatches), Is.EqualTo(plain.Select(s => s.Mismatches)));
            Assert.That(annotated.Select(h => h.Site.MismatchPositions),
                Is.EqualTo(plain.Select(s => s.MismatchPositions)));
            Assert.That(annotated.Select(h => h.Site.IsForwardStrand),
                Is.EqualTo(plain.Select(s => s.IsForwardStrand)));
            Assert.That(annotated.Select(h => h.Site.OffTargetScore),
                Is.EqualTo(plain.Select(s => s.OffTargetScore)));
            Assert.That(annotated.Select(h => h.Pam), Is.All.Not.Empty);

            Assert.Throws<ArgumentNullException>(() =>
                CrisprDesigner.FindOffTargetsWithPam("", genome, 4).ToList());
            Assert.Throws<ArgumentNullException>(() =>
                CrisprDesigner.FindOffTargetsWithPam(AltPamGuide, null!, 4).ToList());
            Assert.Throws<ArgumentOutOfRangeException>(() =>
                CrisprDesigner.FindOffTargetsWithPam(AltPamGuide, genome, 6).ToList());
            Assert.Throws<ArgumentOutOfRangeException>(() =>
                CrisprDesigner.FindOffTargetsWithPam(AltPamGuide, genome, -1).ToList());
            Assert.Throws<ArgumentException>(() =>
                CrisprDesigner.FindOffTargetsWithPam("ACGT", genome, 4).ToList());
        });
    }

    /// <summary>
    /// F50-010: the annotated aggregate validates its inputs — the guide, the list, its elements, each
    /// element's PAM (the rule is defined on it) and the protospacer length — and refuses
    /// PAM-before-target systems exactly as the other MIT entry points do.
    /// </summary>
    [Test]
    public void CalculateMitSpecificityScore_Annotated_InvalidInputs_Throw()
    {
        var ok = new[] { Hit(AltPamOffTarget, "AGG") };

        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentNullException>(() =>
                CrisprDesigner.CalculateMitSpecificityScore("", ok));
            Assert.Throws<ArgumentNullException>(() =>
                CrisprDesigner.CalculateMitSpecificityScore(
                    AltPamGuide, (IEnumerable<PamAnnotatedOffTarget>)null!));
            Assert.Throws<ArgumentException>(() =>
                CrisprDesigner.CalculateMitSpecificityScore("ACGT", ok));
            Assert.Throws<ArgumentException>(() =>
                CrisprDesigner.CalculateMitSpecificityScore(
                    AltPamGuide, new PamAnnotatedOffTarget[] { null! }));
            Assert.Throws<ArgumentException>(() =>
                CrisprDesigner.CalculateMitSpecificityScore(AltPamGuide, new[] { Hit(AltPamOffTarget, "") }));
            Assert.Throws<ArgumentException>(() =>
                CrisprDesigner.CalculateMitSpecificityScore(AltPamGuide, new[] { Hit("ACGTACGT", "AGG") }));
            Assert.Throws<ArgumentException>(() =>
                CrisprDesigner.CalculateMitSpecificityScore(AltPamGuide, ok, CrisprSystemType.Cas12a));
            Assert.Throws<ArgumentException>(() =>
                CrisprDesigner.CalculateMitSpecificityScoreWithAlternativePams(
                    AltPamGuide, AltPamGenome(), 4, CrisprSystemType.CasX));
        });
    }

    #endregion

    #region B08 round 2 E1 — NGA / NNGRRN alternative-PAM scans and ALTPAMMINSCORE

    // Oracle: crispor.py (maximilianh/crisporWebsite, master) calcHitScore / hitScoreM / calcMitGuideScore
    // copied verbatim into Python and run on the constructed genomes, with an independent brute-force
    // two-strand PAM scan confirming exactly which hits exist:
    //   calcHitScore("GACGCATAAAGATGAGACGC", "GACGCTTAAAGATGAGACGC") = 60.5   (1 mismatch, idx 5)
    //   calcHitScore(guide, "GACGCATAAAGATGAGACTC") = 31.499999999999993       (1 mismatch, idx 18)
    //   calcHitScore(guide, "GTCGCTTCAAGATGAGACTC") = 0.3128679435483871       (4 mismatches)
    // crispor.py L1955: an alternative-PAM (non-GG) hit of an NGG guide is multiplied by 0.2.

    private const string NgaPam = "CGA";
    private const string FourMismatchOffTarget = "GTCGCTTCAAGATGAGACTC";   // 4 mismatches, calcHitScore 0.3128679435483871

    /// <summary>
    /// E1 oracle: one 1-mismatch hit carrying an NGA PAM. Before the fix the SpCas9 alt-PAM aggregate
    /// returned 100.0 (the NGA scan was silently dropped); CRISPOR gives 100/(100 + 0.2 x 60.5) x 100.
    /// </summary>
    [Test]
    public void CalculateMitSpecificityScoreWithAlternativePams_NgaHit_IsScannedAndDownweighted()
    {
        var genome = new DnaSequence(AltPamOffTarget + NgaPam + "TTTTTT");

        var ngaHits = CrisprDesigner.FindOffTargetsWithPam(AltPamGuide, genome, 4, CrisprSystemType.SpCas9_NGA).ToList();
        double hitScore = CrisprDesigner.CalculateMitHitScore(AltPamGuide, AltPamOffTarget);

        Assert.Multiple(() =>
        {
            Assert.That(ngaHits, Has.Count.EqualTo(1));
            Assert.That(ngaHits[0].Pam, Is.EqualTo(NgaPam));
            Assert.That(CrisprDesigner.FindOffTargetsWithPam(AltPamGuide, genome, 4), Is.Empty,
                "no NGG site: the single-system aggregate cannot see this hit");

            // library vs the independent Python calcHitScore
            Assert.That(hitScore, Is.EqualTo(60.5).Within(1e-9));
            Assert.That(100.0 / (100.0 + 0.2 * hitScore) * 100.0, Is.EqualTo(89.20606601248885).Within(1e-9));

            Assert.That(
                CrisprDesigner.CalculateMitSpecificityScoreWithAlternativePams(AltPamGuide, genome, 4),
                Is.EqualTo(89.20606601248885).Within(1e-9));
            Assert.That(
                CrisprDesigner.CalculateMitSpecificityScore(AltPamGuide, genome, 4),
                Is.EqualTo(100.0), "the single-system overload is unchanged");
        });
    }

    /// <summary>
    /// E1: the SaCas9 NNGRRT -> NNGRRN alternative. crispor.py gates the 0.2 factor on pam=="NGG"
    /// (L1955), so an NNGRRN hit of a SaCas9 guide enters unweighted: calcSaHitScore 4.2613097364785695
    /// (the value F9-001 locks) -> 100/(100+4.2613097364785695) x 100 = 95.91285612347565.
    /// A site that matches both NNGRRT and NNGRRN is ONE off-target (reference: one search), so adding
    /// an NNGRRT twin gives 100/(100+2 x 4.2613...) x 100 = 92.14668839146395 and not three hits' worth.
    /// </summary>
    [Test]
    public void CalculateMitSpecificityScoreWithAlternativePams_SaCas9_NnggrrnHit_IsScannedUnweighted_AndNotDoubleCounted()
    {
        var altOnly = new DnaSequence("ACGTACGTACGTACGTACGTC" + "AAGAAC" + "TTTTTT");
        var both = new DnaSequence("ACGTACGTACGTACGTACGTC" + "AAGAAC" + "TTTT" + "ACGTACGTACGTACGTACGTC" + SaPam);

        var relaxedHits = CrisprDesigner.FindOffTargetsWithPam(SaGuide21, altOnly, 3, CrisprSystemType.SaCas9_NNGRRN).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(relaxedHits, Has.Count.EqualTo(1));
            Assert.That(relaxedHits[0].Pam, Is.EqualTo("AAGAAC"));
            Assert.That(CrisprDesigner.FindOffTargetsWithPam(SaGuide21, altOnly, 3, CrisprSystemType.SaCas9), Is.Empty);

            Assert.That(
                CrisprDesigner.CalculateMitSpecificityScoreWithAlternativePams(
                    SaGuide21, altOnly, 3, CrisprSystemType.SaCas9),
                Is.EqualTo(95.91285612347565).Within(1e-9));
            Assert.That(
                CrisprDesigner.CalculateMitSpecificityScore(SaGuide21, altOnly, 3, CrisprSystemType.SaCas9),
                Is.EqualTo(100.0), "the single-system overload is unchanged");

            // NNGRRT hit + NNGRRN-only hit: both are counted once each.
            Assert.That(
                CrisprDesigner.CalculateMitSpecificityScoreWithAlternativePams(
                    SaGuide21, both, 3, CrisprSystemType.SaCas9),
                Is.EqualTo(92.14668839146395).Within(1e-9));
            Assert.That(
                CrisprDesigner.CalculateMitSpecificityScore(SaGuide21, both, 3, CrisprSystemType.SaCas9),
                Is.EqualTo(95.91285612347565).Within(1e-9));
        });
    }

    /// <summary>
    /// E1: the two appended systems carry the reference's PAM / guide length / orientation, and the
    /// alternative-PAM table now resolves a scannable system for every alternative the reference lists
    /// for NGG and NNGRRT.
    /// </summary>
    [Test]
    public void NewAlternativePamSystems_HaveReferenceMetadata_AndExistingEnumValuesAreUnchanged()
    {
        var nga = CrisprDesigner.GetSystem(CrisprSystemType.SpCas9_NGA);
        var sa = CrisprDesigner.GetSystem(CrisprSystemType.SaCas9_NNGRRN);

        Assert.Multiple(() =>
        {
            Assert.That((nga.PamSequence, nga.GuideLength, nga.PamAfterTarget), Is.EqualTo(("NGA", 20, true)));
            Assert.That((sa.PamSequence, sa.GuideLength, sa.PamAfterTarget), Is.EqualTo(("NNGRRN", 21, true)));
            // appended at the end: pre-existing numeric values are unchanged
            Assert.That((int)CrisprSystemType.SpCas9, Is.EqualTo(0));
            Assert.That((int)CrisprSystemType.CasX, Is.EqualTo(6));
            Assert.That((int)CrisprSystemType.SpCas9_NGA, Is.EqualTo(7));
            Assert.That((int)CrisprSystemType.SaCas9_NNGRRN, Is.EqualTo(8));
            // NGA is a key of the reference's table too (NGA -> NGG)
            Assert.That(CrisprDesigner.GetAlternativeOffTargetPams(CrisprSystemType.SpCas9_NGA),
                Is.EqualTo(new[] { "NGG" }));
            Assert.That(CrisprDesigner.GetAlternativeOffTargetPams(CrisprSystemType.SaCas9_NNGRRN), Is.Empty);
        });
    }

    /// <summary>
    /// E1/F67: <c>minAlternativePamHitScore</c> is crispor.py's ALTPAMMINSCORE (L4108-4122, default 1.0
    /// there and here since G6). Genome = NGG hit (60.5), NAG hit with 4 mismatches (0.3128679435483871) and NGA
    /// hit (31.499999999999993). Aggregates (Python, crispor formulae):
    /// default (= 1.0) and explicit 1.0 drop only the NAG hit -> 59.95203836930455;
    /// 40 drops both alternative hits but never the main NGG hit -> 62.30529595015576.
    /// </summary>
    [Test]
    public void CalculateMitSpecificityScoreWithAlternativePams_MinAlternativePamHitScore_FiltersAlternativeHitsOnly()
    {
        var genome = new DnaSequence(
            AltPamOffTarget + "AGG" + "TTTT" + FourMismatchOffTarget + "TAG" + "TTTT" + ProximalOffTarget + NgaPam);

        Assert.Multiple(() =>
        {
            Assert.That(CrisprDesigner.CalculateMitHitScore(AltPamGuide, FourMismatchOffTarget),
                Is.EqualTo(0.3128679435483871).Within(1e-12));
            Assert.That(CrisprDesigner.CalculateMitHitScore(AltPamGuide, ProximalOffTarget),
                Is.EqualTo(31.499999999999993).Within(1e-12));

            // G6: the default is CRISPOR's own ALTPAMMINSCORE = 1.0 (it was 0 = no filtering); explicit 0 keeps all hits.
            Assert.That(CrisprDesigner.CalculateMitSpecificityScoreWithAlternativePams(AltPamGuide, genome, 4),
                Is.EqualTo(59.95203836930455).Within(1e-9), "default == explicit 1.0");
            Assert.That(
                CrisprDesigner.CalculateMitSpecificityScoreWithAlternativePams(
                    AltPamGuide, genome, 4, CrisprSystemType.SpCas9, 0.0),
                Is.EqualTo(59.92955631050283).Within(1e-9), "explicit 0 = no filtering");
            Assert.That(
                CrisprDesigner.CalculateMitSpecificityScoreWithAlternativePams(
                    AltPamGuide, genome, 4, CrisprSystemType.SpCas9, 1.0),
                Is.EqualTo(59.95203836930455).Within(1e-9));
            Assert.That(
                CrisprDesigner.CalculateMitSpecificityScoreWithAlternativePams(
                    AltPamGuide, genome, 4, CrisprSystemType.SpCas9, 40.0),
                Is.EqualTo(62.30529595015576).Within(1e-9));
            Assert.That(
                CrisprDesigner.CalculateMitSpecificityScore(AltPamGuide, genome, 4),
                Is.EqualTo(62.30529595015576).Within(1e-9));

            Assert.Throws<ArgumentOutOfRangeException>(() =>
                CrisprDesigner.CalculateMitSpecificityScoreWithAlternativePams(
                    AltPamGuide, genome, 4, CrisprSystemType.SpCas9, -1.0));
            Assert.Throws<ArgumentOutOfRangeException>(() =>
                CrisprDesigner.CalculateMitSpecificityScoreWithAlternativePams(
                    AltPamGuide, genome, 4, CrisprSystemType.SpCas9, double.NaN));
        });
    }

    /// <summary>
    /// E1/F67: reference quirk kept on purpose — crispor.py tests <c>genomePamSeq in altPamSeqs</c> BEFORE
    /// the main list, and the expansion of NNGRRN contains every NNGRRT PAM, so for SaCas9 the cutoff also
    /// applies to main-PAM hits. The single SaCas9 hit scores 4.2613097364785695: kept at cutoff 4,
    /// dropped (aggregate 100) at cutoff 5.
    /// </summary>
    [Test]
    public void CalculateMitSpecificityScoreWithAlternativePams_SaCas9_CutoffAppliesToMainPamHitsToo()
    {
        var genome = new DnaSequence("ACGTACGTACGTACGTACGTC" + SaPam + "TTTTTT");

        Assert.Multiple(() =>
        {
            Assert.That(
                CrisprDesigner.CalculateMitSpecificityScoreWithAlternativePams(
                    SaGuide21, genome, 3, CrisprSystemType.SaCas9, 4.0),
                Is.EqualTo(95.91285612347565).Within(1e-9));
            Assert.That(
                CrisprDesigner.CalculateMitSpecificityScoreWithAlternativePams(
                    SaGuide21, genome, 3, CrisprSystemType.SaCas9, 5.0),
                Is.EqualTo(100.0));
        });
    }

    #endregion

    #region B08 round 3 G1 — CFD aggregate over the alternative PAMs (crispor.py annotateOfftargets L1905-1990)

    // Oracle (B08 round 3 G1, F80): crispor.py L2043-2185 (hitScoreM, calcHitScore, calcMitGuideScore,
    // calcCfdGuideScore, calc_cfd, calcCfdScore) exec'd verbatim in Python with the pickled CFD tables
    // (mismatch_score.pkl / pam_scores.pkl, pam_scores["AG"] = 0.25925925899999996; the file's own doctests
    // reproduce), plus an independent two-strand brute-force PAM scan reproducing the L4104-4122 pre-filter
    // (alternative-PAM hit kept iff calcHitScore >= ALTPAMMINSCORE) and annotateOfftargets (cfdScore =
    // calcCfdScore(guideSeq, otSeq), NO 0.2 factor; guide score = 100/(100 + sum(cfd)) x 100, unrounded):
    //   NAG-only genome  : one hit (TAG, 1 mm, calcHitScore 60.5), CFD 0.18518518492592592 -> 99.81515711670932
    //   NGA-only genome  : one hit (CGA, calcHitScore 60.5),        CFD 0.049603174265873025 -> 99.95042141828442
    //   mixed genome     : AGG hit (60.5, CFD 0.7142857140000001), TAG 4-mm hit (calcHitScore 0.3128679435483871,
    //                      CFD 0.03848003846945646), CGA hit (31.499999999999993, CFD 0.04629629602314815)
    //                      cutoff 1.0 (default) -> 99.24515917350747 ; 0 -> 99.20727233740708 ;
    //                      0.3 -> 99.20727233740708 (the NAG hit's MIT score 0.3129 passes, its CFD 0.0385 would not) ;
    //                      40 -> 99.29078014212564 (= main-PAM-only aggregate)
    //   NAG-only genome, ParkinsonNormalized (100 / sum of raw CFD) -> 540.000000756

    private const string NagPam = "TAG";

    [Test]
    public void CalculateCfdSpecificityScoreWithAlternativePams_NagOnlyGuide_IsNoLongerPerfect()
    {
        var genome = new DnaSequence(AltPamOffTarget + NagPam + "TTTTTT");

        Assert.Multiple(() =>
        {
            Assert.That(CrisprDesigner.CalculateCfdScore(AltPamGuide, AltPamOffTarget, NagPam),
                Is.EqualTo(0.18518518492592592).Within(1e-12), "per-hit CFD vs the Python calcCfdScore");
            Assert.That(CrisprDesigner.CalculateCfdSpecificityScore(AltPamGuide, genome, 4),
                Is.EqualTo(100.0), "main-PAM-only overload cannot see a NAG site (unchanged)");
            Assert.That(CrisprDesigner.CalculateCfdSpecificityScoreWithAlternativePams(AltPamGuide, genome, 4),
                Is.EqualTo(99.81515711670932).Within(1e-9));
            Assert.That(
                CrisprDesigner.CalculateCfdSpecificityScoreWithAlternativePams(
                    AltPamGuide, genome, 4, CrisprSystemType.SpCas9, CfdAggregation.ParkinsonNormalized),
                Is.EqualTo(540.000000756).Within(1e-6));
            // no 0.2 factor on the CFD path: the value is the plain sum of the raw CFD scores
            Assert.That(
                CrisprDesigner.CalculateCfdSpecificityScore(new[] { 0.18518518492592592 }),
                Is.EqualTo(CrisprDesigner.CalculateCfdSpecificityScoreWithAlternativePams(AltPamGuide, genome, 4))
                    .Within(1e-12));
        });
    }

    [Test]
    public void CalculateCfdSpecificityScoreWithAlternativePams_NgaOnlyGuide_IsScanned()
    {
        var genome = new DnaSequence(AltPamOffTarget + NgaPam + "TTTTTT");

        Assert.Multiple(() =>
        {
            Assert.That(CrisprDesigner.CalculateCfdSpecificityScore(AltPamGuide, genome, 4), Is.EqualTo(100.0));
            Assert.That(CrisprDesigner.CalculateCfdSpecificityScoreWithAlternativePams(AltPamGuide, genome, 4),
                Is.EqualTo(99.95042141828442).Within(1e-9));
        });
    }

    /// <summary>
    /// Mixed genome (NGG 60.5 + NAG 4-mismatch 0.3129 + NGA 31.5, the F67 genome): the cutoff acts on the
    /// MIT hit score (as L4104-4122 does), never on the CFD value, and never on a main-PAM hit.
    /// </summary>
    [Test]
    public void CalculateCfdSpecificityScoreWithAlternativePams_MixedGenome_AppliesTheMitHitScoreCutoff()
    {
        var genome = new DnaSequence(
            AltPamOffTarget + "AGG" + "TTTT" + FourMismatchOffTarget + NagPam + "TTTT" + ProximalOffTarget + NgaPam);

        Assert.Multiple(() =>
        {
            Assert.That(CrisprDesigner.CalculateCfdSpecificityScoreWithAlternativePams(AltPamGuide, genome, 4),
                Is.EqualTo(99.24515917350747).Within(1e-9), "default == ALTPAMMINSCORE 1.0: the 0.3129 NAG hit is dropped");
            Assert.That(CrisprDesigner.CalculateCfdSpecificityScoreWithAlternativePams(
                    AltPamGuide, genome, 4, CrisprSystemType.SpCas9, CfdAggregation.MitStyle, 0.0),
                Is.EqualTo(99.20727233740708).Within(1e-9));
            Assert.That(CrisprDesigner.CalculateCfdSpecificityScoreWithAlternativePams(
                    AltPamGuide, genome, 4, CrisprSystemType.SpCas9, CfdAggregation.MitStyle, 0.3),
                Is.EqualTo(99.20727233740708).Within(1e-9), "0.3 is a MIT-scale cutoff: NAG hit (MIT 0.3129, CFD 0.0385) is kept");
            Assert.That(CrisprDesigner.CalculateCfdSpecificityScoreWithAlternativePams(
                    AltPamGuide, genome, 4, CrisprSystemType.SpCas9, CfdAggregation.MitStyle, 40.0),
                Is.EqualTo(99.29078014212564).Within(1e-9));
            Assert.That(CrisprDesigner.CalculateCfdSpecificityScore(AltPamGuide, genome, 4),
                Is.EqualTo(99.29078014212564).Within(1e-9), "main-PAM-only aggregate");

            Assert.Throws<ArgumentOutOfRangeException>(() =>
                CrisprDesigner.CalculateCfdSpecificityScoreWithAlternativePams(
                    AltPamGuide, genome, 4, CrisprSystemType.SpCas9, CfdAggregation.MitStyle, -1.0));
            Assert.Throws<ArgumentOutOfRangeException>(() =>
                CrisprDesigner.CalculateCfdSpecificityScoreWithAlternativePams(
                    AltPamGuide, genome, 4, CrisprSystemType.SpCas9, CfdAggregation.MitStyle, double.NaN));
        });
    }

    [Test]
    public void CalculateCfdSpecificityScoreWithAlternativePams_KeepsTheCfdGuards()
    {
        var genome = new DnaSequence(AltPamOffTarget + NagPam + "TTTTTT");

        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentNullException>(() =>
                CrisprDesigner.CalculateCfdSpecificityScoreWithAlternativePams("", genome));
            Assert.Throws<ArgumentNullException>(() =>
                CrisprDesigner.CalculateCfdSpecificityScoreWithAlternativePams(AltPamGuide, null!));
            Assert.Throws<ArgumentOutOfRangeException>(() =>
                CrisprDesigner.CalculateCfdSpecificityScoreWithAlternativePams(AltPamGuide, genome, 6));
            Assert.Throws<ArgumentException>(() =>
                CrisprDesigner.CalculateCfdSpecificityScoreWithAlternativePams(SaGuide21, genome, 3, CrisprSystemType.SaCas9),
                "21-nt SaCas9: CFD refused as in the main-PAM overload");
            Assert.Throws<ArgumentException>(() =>
                CrisprDesigner.CalculateCfdSpecificityScoreWithAlternativePams(
                    AltPamGuide, genome, 3, CrisprSystemType.Cas12a));
            Assert.That(CrisprDesigner.CalculateCfdSpecificityScoreWithAlternativePams(
                    AltPamGuide.Replace('T', 'U'), genome, 4),
                Is.EqualTo(99.81515711670932).Within(1e-9), "RNA spelling == DNA spelling (F77 rule)");
        });
    }

    #endregion
}
