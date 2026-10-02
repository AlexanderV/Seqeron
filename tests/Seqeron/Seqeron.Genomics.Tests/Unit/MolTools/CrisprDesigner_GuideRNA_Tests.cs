namespace Seqeron.Genomics.Tests.Unit.MolTools;

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
    // Standard SpCas9 scaffold (76 nt)
    private const string SpCas9Scaffold =
        "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC";

    #region MUST Tests - Guide RNA Evaluation

    /// <summary>
    /// M-001: Optimal guide (50% GC, no polyT, no penalties) → perfect score.
    /// Evidence: Addgene — guides with optimal GC (40-70%) perform better.
    /// Input: "ACGTACGTACGTACGTACGT" — 50% GC, no polyT, selfComp=0.15 (below 0.3 threshold).
    /// Score: 100 (base) − 0 = 100.
    /// </summary>
    [Test]
    public void EvaluateGuideRna_OptimalGuide_HighScore()
    {
        string guide = "ACGTACGTACGTACGTACGT"; // 50% GC, no polyT
        var candidate = CrisprDesigner.EvaluateGuideRna(guide, CrisprSystemType.SpCas9);

        Assert.That(candidate.Score, Is.EqualTo(100));
        Assert.That(candidate.GcContent, Is.EqualTo(50));
        Assert.That(candidate.SeedGcContent, Is.EqualTo(50));
        Assert.That(candidate.HasPolyT, Is.False);
        Assert.That(candidate.SelfComplementarityScore, Is.EqualTo(0.15));
        Assert.That(candidate.Issues, Is.Empty);
    }

    /// <summary>
    /// M-002: 0% GC → strong penalty: (40-0)×2 = 80, plus seed GC 0% → −5.
    /// Evidence: Wikipedia — "GC content of sgRNA should optimally be over 50%."
    /// Score: 100 − 80 − 5 = 15.
    /// </summary>
    [Test]
    public void EvaluateGuideRna_LowGcContent_LowerScore()
    {
        string guide = "AAAAAAAAAAAAAAAAAAAA"; // 0% GC
        var candidate = CrisprDesigner.EvaluateGuideRna(guide, CrisprSystemType.SpCas9);

        Assert.That(candidate.Score, Is.EqualTo(15));
        Assert.That(candidate.GcContent, Is.EqualTo(0));
        Assert.That(candidate.SeedGcContent, Is.EqualTo(0));
        Assert.That(candidate.HasPolyT, Is.False);
        Assert.That(candidate.Issues, Has.Count.EqualTo(2));
        Assert.That(candidate.Issues, Has.Some.Contains("Low GC"));
        Assert.That(candidate.Issues, Has.Some.Contains("Suboptimal seed region GC"));
    }

    /// <summary>
    /// M-003: 100% GC → penalty: (100-70)×2 = 60, plus seed GC 100% → −5.
    /// Score: 100 − 60 − 5 = 35.
    /// </summary>
    [Test]
    public void EvaluateGuideRna_HighGcContent_LowerScore()
    {
        string guide = "GCGCGCGCGCGCGCGCGCGC"; // 100% GC
        var candidate = CrisprDesigner.EvaluateGuideRna(guide, CrisprSystemType.SpCas9);

        Assert.That(candidate.Score, Is.EqualTo(35));
        Assert.That(candidate.GcContent, Is.EqualTo(100));
        Assert.That(candidate.SeedGcContent, Is.EqualTo(100));
        Assert.That(candidate.HasPolyT, Is.False);
        Assert.That(candidate.Issues, Has.Count.EqualTo(2));
        Assert.That(candidate.Issues, Has.Some.Contains("High GC"));
        Assert.That(candidate.Issues, Has.Some.Contains("Suboptimal seed region GC"));
    }

    /// <summary>
    /// M-004: PolyT (TTTT) detected → −20 penalty.
    /// Evidence: Addgene — "RNA polymerase III terminates at poly-T sequences."
    /// Input: "ACGTACGTTTTTACGTACGT" — 40% GC (no GC penalty), polyT present.
    /// Score: 100 − 20 = 80.
    /// </summary>
    [Test]
    public void EvaluateGuideRna_HasPolyT_Penalized()
    {
        string guide = "ACGTACGTTTTTACGTACGT"; // 40% GC, contains 5 consecutive T's
        var candidate = CrisprDesigner.EvaluateGuideRna(guide, CrisprSystemType.SpCas9);

        Assert.That(candidate.Score, Is.EqualTo(80));
        Assert.That(candidate.GcContent, Is.EqualTo(40));
        Assert.That(candidate.HasPolyT, Is.True);
        Assert.That(candidate.Issues, Has.Count.EqualTo(1));
        Assert.That(candidate.Issues[0], Does.Contain("TTTT"));
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
    /// M-006: FullGuideRna = spacer + scaffold (76 nt).
    /// Evidence: Addgene — "sgRNA composed of a scaffold sequence necessary for Cas-binding and a spacer."
    /// Total length: 20 + 76 = 96.
    /// </summary>
    [Test]
    public void EvaluateGuideRna_FullGuideRna_IncludesScaffold()
    {
        string guide = "ACGTACGTACGTACGTACGT"; // 20bp
        var candidate = CrisprDesigner.EvaluateGuideRna(guide, CrisprSystemType.SpCas9);

        Assert.That(candidate.FullGuideRna, Is.EqualTo(guide + SpCas9Scaffold));
        Assert.That(candidate.FullGuideRna.Length, Is.EqualTo(96));
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
    /// S-001: Restriction site penalty (−5).
    /// Evidence: Common restriction sites interfere with cloning.
    /// Input: "ACGTGAATTCACGTACGTAC" — 45% GC, contains EcoRI site (GAATTC).
    /// Score: 100 − 5 = 95.
    /// </summary>
    [Test]
    public void EvaluateGuideRna_RestrictionSite_Penalized()
    {
        string guide = "ACGTGAATTCACGTACGTAC"; // Contains GAATTC (EcoRI)
        var candidate = CrisprDesigner.EvaluateGuideRna(guide, CrisprSystemType.SpCas9);

        Assert.That(candidate.Score, Is.EqualTo(95));
        Assert.That(candidate.GcContent, Is.EqualTo(45));
        Assert.That(candidate.Issues, Has.Count.EqualTo(1));
        Assert.That(candidate.Issues[0], Does.Contain("restriction site"));
    }

    /// <summary>
    /// S-002: Seed GC content is calculated from last 10bp.
    /// Evidence: Addgene — seed region (8-10bp at 3') initiates annealing; implementation uses 10bp upper bound.
    /// Input: "AAAAAAAAAAAAACGTACGT" — overall GC 20%, seed (last 10) GC 40%.
    /// Score: 100 − (40−20)×2 = 60 (Low GC penalty only; seed in [30,80]).
    /// </summary>
    [Test]
    public void EvaluateGuideRna_CalculatesSeedGc()
    {
        string guide = "AAAAAAAAAAAAACGTACGT"; // Overall GC=20%, seed last 10 = "AACGTACGT" → 40% GC
        var candidate = CrisprDesigner.EvaluateGuideRna(guide, CrisprSystemType.SpCas9);

        Assert.That(candidate.SeedGcContent, Is.EqualTo(40));
        Assert.That(candidate.GcContent, Is.EqualTo(20));
        Assert.That(candidate.Score, Is.EqualTo(60));
        Assert.That(candidate.Issues, Has.Count.EqualTo(1));
        Assert.That(candidate.Issues[0], Does.Contain("Low GC"));
    }

    /// <summary>
    /// S-006: Boundary GC at exactly 40% → no GC penalty.
    /// Evidence: 40% is MinGcContent default — lower boundary inclusive.
    /// Input: "AAAAAAAAAAAAGCGCGCGC" — 8 G/C of 20 = 40%, seed GC = 80%.
    /// Score: 100 (no penalties — seed GC 80% is ≤ 80 threshold).
    /// </summary>
    [Test]
    public void EvaluateGuideRna_BoundaryGc40Percent_NotPenalized()
    {
        string guide = "AAAAAAAAAAAAGCGCGCGC"; // 8 G/C = 40% GC
        var candidate = CrisprDesigner.EvaluateGuideRna(guide, CrisprSystemType.SpCas9);

        Assert.That(candidate.GcContent, Is.EqualTo(40));
        Assert.That(candidate.SeedGcContent, Is.EqualTo(80));
        Assert.That(candidate.Score, Is.EqualTo(100));
        Assert.That(candidate.Issues, Is.Empty);
    }

    /// <summary>
    /// S-007: Boundary GC at exactly 70% → no GC penalty.
    /// Evidence: 70% is MaxGcContent default — upper boundary inclusive.
    /// Input: "GCGCGCGCGCGCGCAAAAAA" — 14 G/C of 20 = 70%, seed GC = 40%.
    /// Score: 100 (no penalties).
    /// </summary>
    [Test]
    public void EvaluateGuideRna_BoundaryGc70Percent_NotPenalized()
    {
        string guide = "GCGCGCGCGCGCGCAAAAAA"; // 14 G/C = 70% GC
        var candidate = CrisprDesigner.EvaluateGuideRna(guide, CrisprSystemType.SpCas9);

        Assert.That(candidate.GcContent, Is.EqualTo(70));
        Assert.That(candidate.SeedGcContent, Is.EqualTo(40));
        Assert.That(candidate.Score, Is.EqualTo(100));
        Assert.That(candidate.Issues, Is.Empty);
    }

    /// <summary>
    /// S-008: Exactly 4 consecutive T's triggers polyT detection.
    /// Evidence: Addgene — TTTT is the minimum for Pol III termination.
    /// Input: "ACGTACGTACGATTTTACGT" — note 'A' at position 11 prevents merge with preceding T.
    /// Score: 100 − 20 = 80.
    /// </summary>
    [Test]
    public void EvaluateGuideRna_ExactlyFourTs_TriggersPolyT()
    {
        string guide = "ACGTACGTACGATTTTACGT"; // Exactly 4 consecutive T's (positions 12-15)
        var candidate = CrisprDesigner.EvaluateGuideRna(guide, CrisprSystemType.SpCas9);

        Assert.That(candidate.HasPolyT, Is.True);
        Assert.That(candidate.GcContent, Is.EqualTo(40));
        Assert.That(candidate.Score, Is.EqualTo(80));
        Assert.That(candidate.Issues, Has.Count.EqualTo(1));
        Assert.That(candidate.Issues[0], Does.Contain("TTTT"));
    }

    /// <summary>
    /// S-009: 3 consecutive T's should NOT trigger polyT detection.
    /// Evidence: TTTT (4+) is the minimum for Pol III termination.
    /// Input: "ACGTACGTACGTACGTTTAC" — 3 consecutive T's at positions 15-17.
    /// Score: 100 (no penalties — 45% GC in range, no polyT, seed 40%).
    /// The guide does carry the Graf et al. 2019 "TT-motif" (last four bases "TTAC" contain TT plus a
    /// C), which CRISPOR reports as a warning without deducting from any score, so it is listed as an
    /// issue while the score stays 100.
    /// </summary>
    [Test]
    public void EvaluateGuideRna_ThreeConsecutiveTs_NoPolyT()
    {
        string guide = "ACGTACGTACGTACGTTTAC"; // 3 consecutive T's (positions 15-17)
        var candidate = CrisprDesigner.EvaluateGuideRna(guide, CrisprSystemType.SpCas9);

        Assert.That(candidate.HasPolyT, Is.False);
        Assert.That(candidate.GcContent, Is.EqualTo(45));
        Assert.That(candidate.Score, Is.EqualTo(100));
        Assert.That(candidate.Issues, Has.Count.EqualTo(1));
        Assert.That(candidate.Issues[0], Does.Contain("Graf 2019 TT-motif"));
        Assert.That(candidate.GrafMotif, Is.EqualTo(GrafMotifType.TtMotif));
    }

    /// <summary>
    /// S-010: Suboptimal seed GC (low) → −5 penalty.
    /// Evidence: Seed region (last 10bp) GC outside 30-80% → penalty.
    /// Input: "GCGCGCGCGCAAAAAAAAAA" — overall GC 50% (no GC penalty), seed GC 0% (all A's) → −5.
    /// Score: 100 − 5 = 95.
    /// </summary>
    [Test]
    public void EvaluateGuideRna_SeedGcLow_Penalized()
    {
        string guide = "GCGCGCGCGCAAAAAAAAAA"; // Overall 50% GC, seed last 10 = all A → 0% GC
        var candidate = CrisprDesigner.EvaluateGuideRna(guide, CrisprSystemType.SpCas9);

        Assert.That(candidate.GcContent, Is.EqualTo(50));
        Assert.That(candidate.SeedGcContent, Is.EqualTo(0));
        Assert.That(candidate.Score, Is.EqualTo(95));
        Assert.That(candidate.Issues, Has.Count.EqualTo(1));
        Assert.That(candidate.Issues[0], Does.Contain("Suboptimal seed region GC"));
    }

    /// <summary>
    /// S-011: Suboptimal seed GC (high) → −5 penalty.
    /// Evidence: Seed region (last 10bp) GC outside 30-80% → penalty.
    /// Input: "AAAAAAAAAAGGGGGGGGGG" — overall GC 50%, seed GC 100% (all G) → −5.
    /// Score: 100 − 5 = 95.
    /// </summary>
    [Test]
    public void EvaluateGuideRna_SeedGcHigh_Penalized()
    {
        string guide = "AAAAAAAAAAGGGGGGGGGG"; // Overall 50% GC, seed last 10 = all G → 100% GC
        var candidate = CrisprDesigner.EvaluateGuideRna(guide, CrisprSystemType.SpCas9);

        Assert.That(candidate.GcContent, Is.EqualTo(50));
        Assert.That(candidate.SeedGcContent, Is.EqualTo(100));
        Assert.That(candidate.Score, Is.EqualTo(95));
        Assert.That(candidate.Issues, Has.Count.EqualTo(1));
        Assert.That(candidate.Issues[0], Does.Contain("Suboptimal seed region GC"));
    }

    #endregion

    #region SHOULD Tests - Guide RNA Design

    /// <summary>
    /// S-003: Design finds guides when PAM is present in target region.
    /// Evidence: Addgene — "target is present immediately adjacent to a PAM."
    /// Returns 1 guide at position 24 (forward strand) with Score 100.
    /// </summary>
    [Test]
    public void DesignGuideRnas_WithPamInRegion_ReturnsGuides()
    {
        var sequence = new DnaSequence("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAGG");
        var guides = CrisprDesigner.DesignGuideRnas(sequence, 20, 45, CrisprSystemType.SpCas9).ToList();

        Assert.That(guides, Has.Count.EqualTo(1));
        Assert.That(guides[0].Position, Is.EqualTo(24));
        Assert.That(guides[0].Score, Is.EqualTo(100));
        Assert.That(guides[0].IsForwardStrand, Is.True);
    }

    #endregion

    #region SHOULD Tests - Parameters

    /// <summary>
    /// S-004: Default parameters have documented values.
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
    /// S-005: Custom parameter values are preserved.
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
        Assert.That(custom.AvoidPolyT, Is.False);
        Assert.That(custom.CheckSelfComplementarity, Is.False);
    }

    #endregion

    #region COULD Tests - Edge Cases

    /// <summary>
    /// C-001: Self-complementarity > 0.3 triggers penalty.
    /// Evidence: Self-complementary regions form secondary structures reducing efficacy.
    /// Uses 8bp period-2 palindrome "GCGCGCGC" (selfComp=0.3125, > 0.3 threshold).
    /// Score: 100 − 60(GC) − 9.375(selfComp×30) − 5(seedGC) = 25.625.
    /// Control: "ACGTACGT" (8bp, selfComp=0.1875, below threshold) → Score 100.
    /// </summary>
    [Test]
    public void EvaluateGuideRna_SelfComplementary_PenaltyTriggered()
    {
        // High self-comp: 8bp period-2 palindrome exceeds 0.3 threshold
        string highSelfComp = "GCGCGCGC";
        var result = CrisprDesigner.EvaluateGuideRna(highSelfComp, CrisprSystemType.SpCas9);

        Assert.That(result.SelfComplementarityScore, Is.EqualTo(0.3125));
        Assert.That(result.Score, Is.EqualTo(25.625));
        Assert.That(result.Issues, Has.Some.Contains("self-complementarity"));

        // Control: same length, below threshold → no self-comp penalty
        string lowSelfComp = "ACGTACGT";
        var control = CrisprDesigner.EvaluateGuideRna(lowSelfComp, CrisprSystemType.SpCas9);

        Assert.That(control.SelfComplementarityScore, Is.EqualTo(0.1875));
        Assert.That(control.Score, Is.EqualTo(100));
        Assert.That(control.Issues.Any(i => i.Contains("self-complementarity")), Is.False);
    }

    /// <summary>
    /// C-002: All-T guide — maximal penalties: low GC + polyT + seed GC → clamped to 0.
    /// Score: 100 − 80(GC) − 20(polyT) − 5(seedGC) = −5 → clamped to 0.
    /// The guide also ends in TTT, i.e. the Graf et al. 2019 "TT-motif", which is reported as a
    /// fourth issue without a score deduction (CRISPOR warns only).
    /// </summary>
    [Test]
    public void EvaluateGuideRna_AllT_VeryLowScoreWithMultipleIssues()
    {
        string guide = "TTTTTTTTTTTTTTTTTTTT"; // 0% GC, polyT throughout
        var candidate = CrisprDesigner.EvaluateGuideRna(guide, CrisprSystemType.SpCas9);

        Assert.That(candidate.Score, Is.EqualTo(0));
        Assert.That(candidate.GcContent, Is.EqualTo(0));
        Assert.That(candidate.SeedGcContent, Is.EqualTo(0));
        Assert.That(candidate.HasPolyT, Is.True);
        Assert.That(candidate.Issues, Has.Count.EqualTo(4));
        Assert.That(candidate.Issues, Has.Some.Contains("Low GC"));
        Assert.That(candidate.Issues, Has.Some.Contains("TTTT"));
        Assert.That(candidate.Issues, Has.Some.Contains("Suboptimal seed region GC"));
        Assert.That(candidate.Issues, Has.Some.Contains("Graf 2019 TT-motif"));
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
    /// C-004: No PAM in region returns empty collection.
    /// Evidence: Guides can only be designed adjacent to PAM.
    /// </summary>
    [Test]
    public void DesignGuideRnas_NoPamInRegion_ReturnsEmpty()
    {
        var sequence = new DnaSequence("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
        var guides = CrisprDesigner.DesignGuideRnas(sequence, 10, 40, CrisprSystemType.SpCas9).ToList();

        Assert.That(guides, Is.Empty);
    }

    /// <summary>
    /// C-005: Multiple PAMs produce multiple guides.
    /// Input: 3 × (20bp + NGG) structure → 3 forward-strand guide candidates.
    /// All guides are "ACGTACGTACGTACGTACGT" scoring 100.
    /// </summary>
    [Test]
    public void DesignGuideRnas_MultiplePams_ReturnsMultipleGuides()
    {
        var sequence = new DnaSequence(
            "ACGTACGTACGTACGTACGTAGG" +  // PAM 1 (AGG)
            "ACGTACGTACGTACGTACGTCGG" +  // PAM 2 (CGG)
            "ACGTACGTACGTACGTACGTTGG");  // PAM 3 (TGG)

        var guides = CrisprDesigner.DesignGuideRnas(sequence, 0, sequence.Length - 1, CrisprSystemType.SpCas9).ToList();

        Assert.That(guides, Has.Count.EqualTo(3));
        Assert.That(guides.All(g => g.Sequence == "ACGTACGTACGTACGTACGT"), Is.True);
        Assert.That(guides.All(g => g.Score == 100), Is.True);
    }

    /// <summary>
    /// C-006: SaCas9 system type evaluates correctly.
    /// Evidence: SaCas9 (Staphylococcus aureus) — NNGRRT PAM, 21bp guide, PAM after target.
    /// Same 20bp input scores identically (no system-specific penalty in EvaluateGuideRna).
    /// </summary>
    [Test]
    public void EvaluateGuideRna_SaCas9SystemType_ValidEvaluation()
    {
        string guide = "ACGTACGTACGTACGTACGT";
        var candidate = CrisprDesigner.EvaluateGuideRna(guide, CrisprSystemType.SaCas9);

        Assert.That(candidate.Sequence, Is.EqualTo(guide));
        Assert.That(candidate.Score, Is.EqualTo(100));
        Assert.That(candidate.GcContent, Is.EqualTo(50));
        Assert.That(candidate.System.Name, Is.EqualTo("SaCas9"));
        Assert.That(candidate.System.GuideLength, Is.EqualTo(21));
        Assert.That(candidate.Issues, Is.Empty);
    }

    /// <summary>
    /// C-007: GC just below 40% boundary triggers Low GC issue.
    /// Input: "AAAAAAAAAAAAGCGCGCAT" — 6 G/C of 20 = 30% GC.
    /// Score: 100 − (40−30)×2 = 80.
    /// </summary>
    [Test]
    public void EvaluateGuideRna_BelowBoundaryGc_HasLowGcIssue()
    {
        string guide = "AAAAAAAAAAAAGCGCGCAT"; // 6 G/C = 30% GC
        var candidate = CrisprDesigner.EvaluateGuideRna(guide, CrisprSystemType.SpCas9);

        Assert.That(candidate.GcContent, Is.EqualTo(30));
        Assert.That(candidate.Score, Is.EqualTo(80));
        Assert.That(candidate.Issues, Has.Count.EqualTo(1));
        Assert.That(candidate.Issues[0], Does.Contain("Low GC"));
    }

    /// <summary>
    /// C-008: Region spanning entire sequence works without error.
    /// Returns exactly 1 guide at position 4.
    /// </summary>
    [Test]
    public void DesignGuideRnas_EntireSequenceAsRegion_Works()
    {
        var sequence = new DnaSequence("ACGTACGTACGTACGTACGTACGTAGG");

        var guides = CrisprDesigner.DesignGuideRnas(
            sequence, 0, sequence.Length - 1, CrisprSystemType.SpCas9).ToList();

        Assert.That(guides, Has.Count.EqualTo(1));
        Assert.That(guides[0].Position, Is.EqualTo(4));
        Assert.That(guides[0].Score, Is.EqualTo(100));
    }

    /// <summary>
    /// C-009: GC above 70% boundary triggers High GC issue.
    /// Input: "GCGCGCGCGCGCGCGCAAAA" — 16 G/C of 20 = 80% GC.
    /// Score: 100 − (80−70)×2 = 80.
    /// </summary>
    [Test]
    public void EvaluateGuideRna_AboveBoundaryGc_HasHighGcIssue()
    {
        string guide = "GCGCGCGCGCGCGCGCAAAA"; // 16 G/C = 80% GC
        var candidate = CrisprDesigner.EvaluateGuideRna(guide, CrisprSystemType.SpCas9);

        Assert.That(candidate.GcContent, Is.EqualTo(80));
        Assert.That(candidate.Score, Is.EqualTo(80));
        Assert.That(candidate.Issues, Has.Count.EqualTo(1));
        Assert.That(candidate.Issues[0], Does.Contain("High GC"));
    }

    /// <summary>
    /// C-010: DesignGuideRnas filters guides below MinScore.
    /// A guide scoring 100 is excluded when MinScore = 101.
    /// </summary>
    [Test]
    public void DesignGuideRnas_MinScoreFiltering_ExcludesLowScoreGuides()
    {
        var sequence = new DnaSequence("ACGTACGTACGTACGTACGTACGTAGG");

        // With default MinScore (50), the guide (score 100) is included
        var defaultGuides = CrisprDesigner.DesignGuideRnas(
            sequence, 0, sequence.Length - 1, CrisprSystemType.SpCas9).ToList();
        Assert.That(defaultGuides, Has.Count.EqualTo(1));

        // With MinScore > 100, no guide qualifies
        var strictParams = new GuideRnaParameters(
            MinGcContent: 40, MaxGcContent: 70, MinScore: 101,
            AvoidPolyT: true, CheckSelfComplementarity: true);
        var strictGuides = CrisprDesigner.DesignGuideRnas(
            sequence, 0, sequence.Length - 1, CrisprSystemType.SpCas9, strictParams).ToList();
        Assert.That(strictGuides, Is.Empty);
    }

    #endregion

    #region 2026-09 review — cut-site frame, ranking, system metadata, parameter flags

    private static readonly GuideRnaParameters NoFilter =
        new(MinGcContent: 0, MaxGcContent: 100, MinScore: 0, AvoidPolyT: false, CheckSelfComplementarity: false);

    /// <summary>
    /// S-012: a reverse-strand guide whose cleavage site is inside the region must be designed.
    /// The SpCas9 blunt cut sits 3 bp 5' of the PAM *on the PAM-bearing strand* (CRISPOR
    /// <c>crispor.py</c>: "the expected cleavage position located -3bp 5' of the PAM site"; its
    /// schematic draws the marker as <c>startFt = start - 3</c> on '+' and as <c>ftSeq + "---"</c>
    /// on '-'), i.e. forward coordinate <c>Position + pamLen + 2</c> on the minus strand.
    /// Repro of the defect: "CCAACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT" has exactly one PAM site,
    /// a reverse-strand TGG at forward position 0 whose protospacer starts at 3 and whose cut is at
    /// forward coordinate 5 — the old code evaluated 0 − 3 = −3 and therefore dropped the guide for
    /// EVERY region, including the whole sequence.
    /// </summary>
    [Test]
    public void DesignGuideRnas_ReverseStrandGuide_CutSiteInRegion_IsDesigned()
    {
        var sequence = new DnaSequence("CCAACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT");

        var site = CrisprDesigner.FindPamSites(sequence, CrisprSystemType.SpCas9).Single();
        Assert.That(site.IsForwardStrand, Is.False);
        Assert.That(site.Position, Is.EqualTo(0));
        Assert.That(site.TargetStart, Is.EqualTo(3));
        Assert.That(CrisprDesigner.GetCutSite(site), Is.EqualTo(5));

        var guides = CrisprDesigner.DesignGuideRnas(
            sequence, 0, sequence.Length - 1, CrisprSystemType.SpCas9, NoFilter).ToList();

        Assert.That(guides, Has.Count.EqualTo(1));
        Assert.That(guides[0].IsForwardStrand, Is.False);
        Assert.That(guides[0].Position, Is.EqualTo(3));
        Assert.That(guides[0].Sequence, Is.EqualTo("ACGTACGTACGTACGTACGT"));
    }

    /// <summary>
    /// S-013: the same reverse-strand guide is excluded when its cleavage site (forward coordinate 5)
    /// falls outside the requested region.
    /// </summary>
    [Test]
    public void DesignGuideRnas_ReverseStrandGuide_CutSiteOutsideRegion_IsExcluded()
    {
        var sequence = new DnaSequence("CCAACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT");

        Assert.That(CrisprDesigner.DesignGuideRnas(sequence, 6, 30, CrisprSystemType.SpCas9, NoFilter),
            Is.Empty, "cut at 5 is left of regionStart 6");
        Assert.That(CrisprDesigner.DesignGuideRnas(sequence, 0, 4, CrisprSystemType.SpCas9, NoFilter),
            Is.Empty, "cut at 5 is right of regionEnd 4");
        Assert.That(CrisprDesigner.DesignGuideRnas(sequence, 5, 5, CrisprSystemType.SpCas9, NoFilter),
            Has.Exactly(1).Items, "cut at 5 is the only base of the region");
    }

    /// <summary>
    /// S-014: <see cref="CrisprDesigner.GetCutSite"/> on both strands and both PAM geometries.
    /// Cas9: 3 bp 5' of the PAM on the PAM strand (CRISPOR, see S-012). Cas12a/Cpf1: the staggered
    /// cut falls "after the 18th base on the non-targeted strand which has the TTTV PAM motif"
    /// (Zetsche et al. 2015 Fig. 3, as described by CRISPOR), so the 19th protospacer base is
    /// returned: <c>Position + pamLen + 18</c> on '+' and <c>Position − 19</c> on '−'.
    /// </summary>
    [Test]
    public void GetCutSite_FollowsCrisporConventionOnBothStrands()
    {
        var spCas9 = CrisprDesigner.GetSystem(CrisprSystemType.SpCas9);
        var cas12a = CrisprDesigner.GetSystem(CrisprSystemType.Cas12a);

        var fwdCas9 = new PamSite(40, "AGG", new string('A', 20), 20, true, spCas9);
        var revCas9 = new PamSite(40, "AGG", new string('A', 20), 43, false, spCas9);
        Assert.That(CrisprDesigner.GetCutSite(fwdCas9), Is.EqualTo(37));
        Assert.That(CrisprDesigner.GetCutSite(revCas9), Is.EqualTo(45));

        var fwdCas12a = new PamSite(40, "TTTA", new string('A', 23), 44, true, cas12a);
        var revCas12a = new PamSite(40, "TTTA", new string('A', 23), 17, false, cas12a);
        Assert.That(CrisprDesigner.GetCutSite(fwdCas12a), Is.EqualTo(62));
        Assert.That(CrisprDesigner.GetCutSite(revCas12a), Is.EqualTo(21));
    }

    /// <summary>
    /// S-015: the designed list is ranked best-first, as CRISPOR sorts its guide table
    /// (<c>crispor.py mergeGuideInfo</c>: <c>guideData.sort(reverse=True, key=…)</c>).
    /// Input plants three PAMs whose protospacers score 100 (50% GC), 80 (30% GC) and 0 (all-T).
    /// </summary>
    [Test]
    public void DesignGuideRnas_Output_IsRankedByScoreDescending()
    {
        var sequence = new DnaSequence(
            "AAAAAAAAAAAAGCGCGCATAGG" +      // 30% GC -> 80
            "TTTTTTTTTTTTTTTTTTTTCGG" +      // all-T  -> 0
            "ACGTACGTACGTACGTACGTTGG");      // 50% GC -> 100

        var guides = CrisprDesigner.DesignGuideRnas(
            sequence, 0, sequence.Length - 1, CrisprSystemType.SpCas9,
            GuideRnaParameters.Default with { MinScore = 0 }).ToList();

        Assert.That(guides.Select(g => g.Score), Is.EqualTo(new[] { 100.0, 80.0, 0.0 }));
        Assert.That(guides[0].Sequence, Is.EqualTo("ACGTACGTACGTACGTACGT"));
        Assert.That(guides[^1].Sequence, Is.EqualTo("TTTTTTTTTTTTTTTTTTTT"));
    }

    /// <summary>
    /// S-016: designed candidates keep the metadata of the system they were designed with. A
    /// name-based remap previously replaced LbCas12a/AsCas12a/CasX/SpCas9-NAG with SpCas9, which also
    /// flipped the seed-region orientation (<c>PamAfterTarget</c>) and therefore changed the score:
    /// the 24-nt LbCas12a protospacer "CGTACGTACGTACGTACGTACGTA" was seeded from its last 10 bases
    /// ("GTACGTACGTA"[^10..] -> GC 40%) instead of its PAM-proximal first 10 ("CGTACGTACG" -> GC 60%).
    /// </summary>
    [Test]
    public void DesignGuideRnas_KeepsTheRequestedSystemMetadata()
    {
        var sequence = new DnaSequence("TTTACGTACGTACGTACGTACGTACGTACGTACGTACGT");

        var lb = CrisprDesigner.DesignGuideRnas(
            sequence, 0, sequence.Length - 1, CrisprSystemType.LbCas12a, NoFilter).ToList();

        Assert.That(lb, Has.Count.EqualTo(1));
        Assert.That(lb[0].System.Name, Is.EqualTo("LbCas12a"));
        Assert.That(lb[0].System.GuideLength, Is.EqualTo(24));
        Assert.That(lb[0].System.PamAfterTarget, Is.False);
        Assert.That(lb[0].Sequence, Has.Length.EqualTo(24));
        // Seed = the 10 PAM-proximal (5') bases "CGTACGTACG" -> 6 G/C of 10 = 60%.
        Assert.That(lb[0].SeedGcContent, Is.EqualTo(60));
    }

    /// <summary>
    /// S-017: <see cref="GuideRnaParameters.AvoidPolyT"/> = <c>false</c> disables the TTTT penalty and
    /// its issue while the measurement <see cref="GuideRnaCandidate.HasPolyT"/> is still reported.
    /// Repro of the defect: the flag was stored but never read, so M-004's guide still scored 80.
    /// </summary>
    [Test]
    public void EvaluateGuideRna_AvoidPolyTFalse_DisablesThePolyTPenalty()
    {
        const string guide = "ACGTACGTTTTTACGTACGT"; // 40% GC, contains TTTTT

        var on = CrisprDesigner.EvaluateGuideRna(guide, CrisprSystemType.SpCas9,
            GuideRnaParameters.Default with { AvoidPolyT = true });
        var off = CrisprDesigner.EvaluateGuideRna(guide, CrisprSystemType.SpCas9,
            GuideRnaParameters.Default with { AvoidPolyT = false });

        Assert.That(on.Score, Is.EqualTo(80));
        Assert.That(off.Score, Is.EqualTo(100));
        Assert.That(off.HasPolyT, Is.True);
        Assert.That(off.Issues, Is.Empty);
    }

    /// <summary>
    /// S-018: <see cref="GuideRnaParameters.CheckSelfComplementarity"/> = <c>false</c> disables the
    /// self-complementarity penalty and its issue while the measurement is still reported.
    /// </summary>
    [Test]
    public void EvaluateGuideRna_CheckSelfComplementarityFalse_DisablesThePenalty()
    {
        const string guide = "GCGCGCGC"; // selfComp 0.3125 (> 0.3), seed GC 100%
        var wide = new GuideRnaParameters(0, 100, 0, AvoidPolyT: true, CheckSelfComplementarity: true);

        var on = CrisprDesigner.EvaluateGuideRna(guide, CrisprSystemType.SpCas9, wide);
        var off = CrisprDesigner.EvaluateGuideRna(guide, CrisprSystemType.SpCas9,
            wide with { CheckSelfComplementarity = false });

        Assert.That(on.Score, Is.EqualTo(100 - 0.3125 * 30 - 5));
        Assert.That(off.Score, Is.EqualTo(95));
        Assert.That(off.SelfComplementarityScore, Is.EqualTo(0.3125));
        Assert.That(off.Issues, Has.Count.EqualTo(1));
        Assert.That(off.Issues[0], Does.Contain("Suboptimal seed region GC"));
    }

    #endregion

    #region 2026-09 review — Doench 2016 Rule Set 2 on designed guides

    // Oracle row from tests/.../TestData/Azimuth/nopos_oracle.csv (agrees == 1, i.e. the verified
    // reference prediction and Microsoft Azimuth's own regression fixture concur):
    //   guide30 = AAAAAAAAAAAAAAAGCTAACAGCAGGAGT, ref_score = 0.529037
    // Layout: 4 nt 5' flank "AAAA" | protospacer "AAAAAAAAAAAGCTAACAGC" | PAM "AGG" | 3' flank "AGT".
    private const string Oracle30Mer = "AAAAAAAAAAAAAAAGCTAACAGCAGGAGT";
    private const double Oracle30MerScore = 0.529037;

    /// <summary>
    /// C-011: a designed forward-strand guide carries the published 30-nt Rule Set 2 context and the
    /// Doench et al. 2016 / Azimuth on-target score for it. CRISPOR builds exactly this window —
    /// <c>crisporEffScores.calcAllScores</c> scores <c>trimSeqs(seqs, -24, 6)</c>, i.e. 4 nt 5' flank
    /// + 20-nt protospacer + 3-nt PAM + 3 nt 3' flank — and reports the resulting Azimuth value as
    /// the guide's efficiency score.
    /// </summary>
    [Test]
    public void DesignGuideRnas_ForwardGuide_CarriesRuleSet2ContextAndScore()
    {
        var guides = CrisprDesigner.DesignGuideRnas(
            new DnaSequence(Oracle30Mer), 0, Oracle30Mer.Length - 1, CrisprSystemType.SpCas9, NoFilter).ToList();

        var guide = guides.Single(g => g.Position == 4 && g.IsForwardStrand);
        Assert.That(guide.Sequence, Is.EqualTo("AAAAAAAAAAAGCTAACAGC"));
        Assert.That(guide.Context30Mer, Is.EqualTo(Oracle30Mer));
        Assert.That(guide.OnTargetScore, Is.Not.Null);
        Assert.That(guide.OnTargetScore!.Value, Is.EqualTo(Oracle30MerScore).Within(1e-5));
    }

    /// <summary>
    /// C-012: the same locus presented on the other strand (the whole window reverse-complemented)
    /// yields a reverse-strand guide whose 30-nt context — read on the protospacer strand, as the
    /// model requires — and whose Azimuth score are identical to C-011's.
    /// </summary>
    [Test]
    public void DesignGuideRnas_ReverseGuide_CarriesTheSameRuleSet2ContextAndScore()
    {
        string reverseInput = DnaSequence.GetReverseComplementString(Oracle30Mer);

        var guides = CrisprDesigner.DesignGuideRnas(
            new DnaSequence(reverseInput), 0, reverseInput.Length - 1, CrisprSystemType.SpCas9, NoFilter).ToList();

        var guide = guides.Single(g => !g.IsForwardStrand && g.Context30Mer is not null);
        Assert.That(guide.Position, Is.EqualTo(6));
        Assert.That(guide.Sequence, Is.EqualTo("AAAAAAAAAAAGCTAACAGC"));
        Assert.That(guide.Context30Mer, Is.EqualTo(Oracle30Mer));
        Assert.That(guide.OnTargetScore!.Value, Is.EqualTo(Oracle30MerScore).Within(1e-5));
    }

    /// <summary>
    /// C-013: the Rule Set 2 context is <c>null</c> (and so is the score) whenever the model's 30-nt
    /// window is unavailable: standalone evaluation, a truncated flank, or a system that is not a
    /// 20-nt NGG system. CRISPOR likewise scores only guides for which the 100-mer context exists and
    /// returns −1 for windows containing N.
    /// </summary>
    [Test]
    public void DesignGuideRnas_WithoutA30MerContext_LeavesTheOnTargetScoreNull()
    {
        // Standalone evaluation never has genomic context.
        var standalone = CrisprDesigner.EvaluateGuideRna("AAAAAAAAAAAGCTAACAGC", CrisprSystemType.SpCas9);
        Assert.That(standalone.Context30Mer, Is.Null);
        Assert.That(standalone.OnTargetScore, Is.Null);

        // Truncated 3' flank (one base short of the 30-mer).
        string truncated = Oracle30Mer[..^1];
        var truncatedGuide = CrisprDesigner.DesignGuideRnas(
                new DnaSequence(truncated), 0, truncated.Length - 1, CrisprSystemType.SpCas9, NoFilter)
            .Single(g => g.Position == 4 && g.IsForwardStrand);
        Assert.That(truncatedGuide.Context30Mer, Is.Null);
        Assert.That(truncatedGuide.OnTargetScore, Is.Null);

        // Truncated 5' flank.
        string shortFlank = Oracle30Mer[1..];
        var shortFlankGuide = CrisprDesigner.DesignGuideRnas(
                new DnaSequence(shortFlank), 0, shortFlank.Length - 1, CrisprSystemType.SpCas9, NoFilter)
            .Single(g => g.Position == 3 && g.IsForwardStrand);
        Assert.That(shortFlankGuide.Context30Mer, Is.Null);

        // Cas12a is not an NGG 20-nt system.
        var cas12a = CrisprDesigner.DesignGuideRnas(
            new DnaSequence("TTTACGTACGTACGTACGTACGTACGTACGTACGTACGT"), 0, 38,
            CrisprSystemType.Cas12a, NoFilter);
        Assert.That(cas12a.All(g => g.Context30Mer is null && g.OnTargetScore is null), Is.True);
    }

    /// <summary>
    /// C-014: <see cref="GuideRnaRanking.OnTargetRuleSet2"/> ranks by the published Doench 2016 score
    /// instead of the composition heuristic (CRISPOR offers the same choice through its
    /// <c>sortBy=fusi</c> column), and guides without a 30-nt context rank last.
    /// </summary>
    [Test]
    public void DesignGuideRnas_RankingByRuleSet2_OrdersByTheAzimuthScore()
    {
        // A leading protospacer whose PAM sits at offset 20 (no room for the 4-nt 5' flank, so it
        // cannot be scored), followed by two loci that do carry the full 30-nt window.
        var sequence = new DnaSequence(
            "ACGTACGTACGTACGTACGTAGG" +
            Oracle30Mer +
            "TTTTACGTACGTACGTACGTACGTAGGCCC");

        var guides = CrisprDesigner.DesignGuideRnas(
            sequence, 0, sequence.Length - 1, CrisprSystemType.SpCas9,
            NoFilter with { Ranking = GuideRnaRanking.OnTargetRuleSet2 }).ToList();

        var scored = guides.Where(g => g.OnTargetScore is not null).ToList();
        Assert.That(scored, Has.Count.GreaterThanOrEqualTo(2));
        Assert.That(scored.Select(g => g.OnTargetScore!.Value), Is.Ordered.Descending);
        Assert.That(guides.SkipWhile(g => g.OnTargetScore is not null).All(g => g.OnTargetScore is null),
            Is.True, "unscored guides must come after every scored guide");
    }

    #endregion

    #region 2026-09 review — Graf et al. 2019 inefficiency motifs

    /// <summary>
    /// C-015: <see cref="CrisprDesigner.GetGrafMotif"/> reproduces <c>getGrafType</c> from CRISPOR's
    /// <c>crisporEffScores.py</c> (Graf et al., Cell Reports 2019, 26:1098–1103). Every expectation
    /// below was produced by running that reference function, not by reading the C# output.
    /// </summary>
    [TestCase("ACGTACGTACGTACGTACGT", GrafMotifType.None)]
    [TestCase("ACGTACGTACGTACGTTTAC", GrafMotifType.TtMotif)]
    [TestCase("TTTTTTTTTTTTTTTTTTTT", GrafMotifType.TtMotif)]
    [TestCase("AAAAAAAAAAAAAAAAATTC", GrafMotifType.TtMotif)]
    [TestCase("AAAAAAAAAAAAAAAAATTT", GrafMotifType.TtMotif)]
    [TestCase("AAAAAAAAAAAAAAAACTCT", GrafMotifType.TtMotif)]
    [TestCase("AAAAAAAAAAAAAAAATCTC", GrafMotifType.TtMotif)]
    [TestCase("AAAAAAAAAAAAAAAACCCT", GrafMotifType.None)]
    [TestCase("AAAAAAAAAAAAAAAAATGC", GrafMotifType.None)]
    [TestCase("AAAAAAAAAAAAAAAAAGCC", GrafMotifType.GccMotif)]
    [TestCase("GAAAAAAAAAAAAAAAAGCC", GrafMotifType.GccMotif)]
    [TestCase("TAAAAAAAAAAAAAAAAGCC", GrafMotifType.GccMotif)]
    [TestCase("CAAAAAAAAAAAAAAAAGCC", GrafMotifType.GccMotif)]
    [TestCase("AAAAAAAAAAAAAAAAGCCT", GrafMotifType.GccMotif)]
    [TestCase("AAAAAAAAAAAAAAAAATTA", GrafMotifType.None)]
    [TestCase("AAAAAAAAAAAAAAAAGTTA", GrafMotifType.None)]
    [TestCase("AAAAAAAAAAAAAAAAACTT", GrafMotifType.TtMotif)]
    [TestCase("AAAAAAAAAAAAAAAACCCC", GrafMotifType.None)]
    [TestCase("ACGTACGTACGTACGTACGC", GrafMotifType.None)]
    [TestCase("ACGTACGTACGTACGTACTT", GrafMotifType.TtMotif)]
    public void GetGrafMotif_MatchesTheCrisporReference(string guide, GrafMotifType expected)
    {
        Assert.That(CrisprDesigner.GetGrafMotif(guide), Is.EqualTo(expected));
        Assert.That(CrisprDesigner.GetGrafMotif(guide.ToLowerInvariant()), Is.EqualTo(expected));
    }

    /// <summary>
    /// C-016: exhaustive agreement with the reference over every 3'-terminal 6-mer (4^6 = 4096 guides
    /// on a fixed 14-nt A prefix). The reference <c>getGrafType</c> classifies exactly 496 as the
    /// TT-motif, 64 as the GCC-motif and 3536 as neither.
    /// </summary>
    [Test]
    public void GetGrafMotif_ExhaustiveSuffixCensus_MatchesTheReferenceCounts()
    {
        const string bases = "ACGT";
        var census = new Dictionary<GrafMotifType, int>
        {
            [GrafMotifType.None] = 0, [GrafMotifType.TtMotif] = 0, [GrafMotifType.GccMotif] = 0,
        };

        for (int code = 0; code < 4096; code++)
        {
            var suffix = new char[6];
            int c = code;
            for (int i = 5; i >= 0; i--) { suffix[i] = bases[c % 4]; c /= 4; }
            census[CrisprDesigner.GetGrafMotif(new string('A', 14) + new string(suffix))]++;
        }

        Assert.That(census[GrafMotifType.TtMotif], Is.EqualTo(496));
        Assert.That(census[GrafMotifType.GccMotif], Is.EqualTo(64));
        Assert.That(census[GrafMotifType.None], Is.EqualTo(3536));
    }

    /// <summary>
    /// C-017: the Graf motif is surfaced as an issue for NGG systems only (as in CRISPOR, which gates
    /// the check on <c>pam == "NGG"</c>) and never deducts from the score — CRISPOR only prints a
    /// warning for it.
    /// </summary>
    [Test]
    public void EvaluateGuideRna_GrafMotif_IsReportedForNggSystemsWithoutAPenalty()
    {
        const string guide = "ACGTACGTACGTACGTACTT"; // TT-motif, 45% GC, no polyT

        var spCas9 = CrisprDesigner.EvaluateGuideRna(guide, CrisprSystemType.SpCas9);
        Assert.That(spCas9.GrafMotif, Is.EqualTo(GrafMotifType.TtMotif));
        Assert.That(spCas9.Score, Is.EqualTo(100), "CRISPOR warns about the Graf motifs, it does not score them");
        Assert.That(spCas9.Issues, Has.Some.Contains("Graf 2019 TT-motif"));

        var saCas9 = CrisprDesigner.EvaluateGuideRna(guide, CrisprSystemType.SaCas9);
        Assert.That(saCas9.GrafMotif, Is.EqualTo(GrafMotifType.None));
        Assert.That(saCas9.Issues, Is.Empty);
    }

    #endregion

    #region 2026-09 review — CHOPCHOP self-complementarity stems

    // Expectations produced by running CHOPCHOP's selfComp() (chopchop.py, STEM_LEN = 4) directly,
    // with no backbone and with the standard backbone "AGGCTAGTCCGT" (which CHOPCHOP reverse-
    // complements to "ACGGACTAGCCT" before matching).
    [TestCase("ACGTACGTACGTACGTACGT", 10, 10)]
    [TestCase("GCGCGCGCGCGCGCGCGCGC", 11, 11)]
    [TestCase("AAAAAAAAAAAAAAAAAAAA", 0, 0)]
    [TestCase("TTTTTTTTTTTTTTTTTTTT", 0, 0)]
    [TestCase("ACGTGAATTCACGTACGTAC", 3, 3)]
    [TestCase("GGGGCCCCGGGGCCCCGGGG", 8, 8)]
    [TestCase("AAAAAAAAAAAAGCGCGCGC", 0, 0)]
    [TestCase("ACGGACTAGCCTACGTACGT", 0, 8)]
    [TestCase("GCGCGCGC", 0, 0)]
    [TestCase("ACGTACGT", 0, 0)]
    [TestCase("ACGT", 0, 0)]
    [TestCase("ACG", 0, 0)]
    [TestCase("GGGGCCCC", 0, 0)]
    public void CountSelfComplementaryStems_MatchesTheChopchopReference(
        string guide, int expectedNoBackbone, int expectedStandardBackbone)
    {
        Assert.That(CrisprDesigner.CountSelfComplementaryStems(guide), Is.EqualTo(expectedNoBackbone));
        Assert.That(CrisprDesigner.CountSelfComplementaryStems(guide.ToLowerInvariant()),
            Is.EqualTo(expectedNoBackbone));
        Assert.That(
            CrisprDesigner.CountSelfComplementaryStems(
                guide, 4, new[] { CrisprDesigner.StandardSgRnaBackboneRegion }),
            Is.EqualTo(expectedStandardBackbone));
    }

    /// <summary>
    /// Exhaustive agreement with the reference over every 9-mer (4^9 = 262 144 guides): CHOPCHOP's
    /// selfComp sums to 704 stems without a backbone (704 guides carrying exactly one) and to 41 632
    /// with the standard backbone.
    /// </summary>
    [Test]
    public void CountSelfComplementaryStems_ExhaustiveNineMerCensus_MatchesTheReference()
    {
        const string bases = "ACGT";
        var backbone = new[] { CrisprDesigner.StandardSgRnaBackboneRegion };
        long sumNoBackbone = 0, sumBackbone = 0;
        int nonZero = 0;
        var buffer = new char[9];

        for (int code = 0; code < 262144; code++)
        {
            int c = code;
            for (int i = 8; i >= 0; i--) { buffer[i] = bases[c % 4]; c /= 4; }
            string guide = new string(buffer);

            int plain = CrisprDesigner.CountSelfComplementaryStems(guide);
            sumNoBackbone += plain;
            if (plain > 0) nonZero++;
            sumBackbone += CrisprDesigner.CountSelfComplementaryStems(guide, 4, backbone);
        }

        Assert.That(sumNoBackbone, Is.EqualTo(704));
        Assert.That(nonZero, Is.EqualTo(704));
        Assert.That(sumBackbone, Is.EqualTo(41632));
    }

    /// <summary>
    /// The stem count is reported on every candidate, and
    /// <see cref="GuideRnaParameters.MaxSelfComplementaryStems"/> turns it into CHOPCHOP's
    /// <c>filterSelfCompMax</c> filter. The default (null) must keep the pre-review behaviour, i.e.
    /// filter nothing — CHOPCHOP's own CLI default is -1 (no filter).
    /// </summary>
    [Test]
    public void DesignGuideRnas_MaxSelfComplementaryStems_FiltersLikeChopchop()
    {
        // The forward protospacer here is "ACGTACGTACGTACGTACGT" -> 10 stems (see the table above).
        var sequence = new DnaSequence("ACGTACGTACGTACGTACGTAGG");

        var unfiltered = CrisprDesigner.DesignGuideRnas(
            sequence, 0, sequence.Length - 1, CrisprSystemType.SpCas9, NoFilter).ToList();
        Assert.That(unfiltered, Has.Count.EqualTo(1));
        Assert.That(unfiltered[0].SelfComplementaryStems, Is.EqualTo(10));
        Assert.That(unfiltered[0].Issues, Is.Empty, "no filter is set, so no issue is raised");

        Assert.That(CrisprDesigner.DesignGuideRnas(sequence, 0, sequence.Length - 1,
            CrisprSystemType.SpCas9, NoFilter with { MaxSelfComplementaryStems = 10 }),
            Has.Exactly(1).Items, "the limit is inclusive");
        Assert.That(CrisprDesigner.DesignGuideRnas(sequence, 0, sequence.Length - 1,
            CrisprSystemType.SpCas9, NoFilter with { MaxSelfComplementaryStems = 9 }),
            Is.Empty);

        var flagged = CrisprDesigner.EvaluateGuideRna("ACGTACGTACGTACGTACGT", CrisprSystemType.SpCas9,
            GuideRnaParameters.Default with { MaxSelfComplementaryStems = 2 });
        Assert.That(flagged.SelfComplementaryStems, Is.EqualTo(10));
        Assert.That(flagged.Issues, Has.Some.Contains("self-complementary 4-bp stems"));
        Assert.That(flagged.Score, Is.EqualTo(100), "the limit filters, it does not deduct");
    }

    /// <summary>
    /// The scaffold regions are matched the way CHOPCHOP matches them: the caller supplies them on
    /// the guide's strand and the reference reverse-complements them internally
    /// (<c>args.backbone = [revComp(el) for el in …]</c>).
    /// </summary>
    [Test]
    public void CountSelfComplementaryStems_BackboneRegions_AreReverseComplementedLikeChopchop()
    {
        // "ACGGACTAGCCT" is revComp("AGGCTAGTCCGT"): as a guide it pairs against the standard
        // backbone in 8 windows but has no internal stem.
        const string guide = "ACGGACTAGCCTACGTACGT";

        Assert.That(CrisprDesigner.CountSelfComplementaryStems(guide), Is.EqualTo(0));
        Assert.That(
            CrisprDesigner.CountSelfComplementaryStems(guide, 4, new[] { "AGGCTAGTCCGT" }),
            Is.EqualTo(8));
        // Passing the already-reverse-complemented string is NOT what CHOPCHOP expects: the
        // orientation collapses the 8 stems to the single window that is its own match in both
        // orientations ("CTAG" occurs in AGGCTAGTCCGT and in ACGGACTAGCCT alike).
        Assert.That(
            CrisprDesigner.CountSelfComplementaryStems(guide, 4, new[] { "ACGGACTAGCCT" }),
            Is.EqualTo(1));
    }

    [Test]
    public void CountSelfComplementaryStems_InvalidArguments_Throw()
    {
        Assert.Throws<ArgumentNullException>(() => CrisprDesigner.CountSelfComplementaryStems(null!));
        Assert.Throws<ArgumentNullException>(() => CrisprDesigner.CountSelfComplementaryStems(""));
        Assert.Throws<ArgumentOutOfRangeException>(
            () => CrisprDesigner.CountSelfComplementaryStems("ACGTACGTACGT", 0));
    }

    #endregion
}
