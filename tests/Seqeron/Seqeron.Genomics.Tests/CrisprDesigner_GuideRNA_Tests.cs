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
    /// </summary>
    [Test]
    public void EvaluateGuideRna_ThreeConsecutiveTs_NoPolyT()
    {
        string guide = "ACGTACGTACGTACGTTTAC"; // 3 consecutive T's (positions 15-17)
        var candidate = CrisprDesigner.EvaluateGuideRna(guide, CrisprSystemType.SpCas9);

        Assert.That(candidate.HasPolyT, Is.False);
        Assert.That(candidate.GcContent, Is.EqualTo(45));
        Assert.That(candidate.Score, Is.EqualTo(100));
        Assert.That(candidate.Issues, Is.Empty);
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
        Assert.That(candidate.Issues, Has.Count.EqualTo(3));
        Assert.That(candidate.Issues, Has.Some.Contains("Low GC"));
        Assert.That(candidate.Issues, Has.Some.Contains("TTTT"));
        Assert.That(candidate.Issues, Has.Some.Contains("Suboptimal seed region GC"));
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
}
