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
    /// Input: "ACGGACTAGCCTACGTACGT" — 55% GC, no polyT, 0 CHOPCHOP self-complementary stems.
    /// Score: 100 (base) − 0 = 100.
    /// </summary>
    [Test]
    public void EvaluateGuideRna_OptimalGuide_HighScore()
    {
        string guide = "ACGGACTAGCCTACGTACGT"; // 55% GC, no polyT, 0 stems (CHOPCHOP oracle)
        var candidate = CrisprDesigner.EvaluateGuideRna(guide, CrisprSystemType.SpCas9);

        Assert.That(candidate.Score, Is.EqualTo(100));
        Assert.That(candidate.GcContent, Is.EqualTo(55).Within(1e-9));
        Assert.That(candidate.SeedGcContent, Is.EqualTo(50));
        Assert.That(candidate.HasPolyT, Is.False);
        Assert.That(candidate.SelfComplementaryStems, Is.EqualTo(0));
        Assert.That(candidate.Issues, Is.Empty);
    }

    /// <summary>
    /// M-002: 0% GC is outside CHOPCHOP's GC window → flat <c>SCORE['CRISPR_BAD_GC'] = 300</c> penalty
    /// (clamped: Score 0); the seed GC 0% is reported as an issue only.
    /// Evidence: Wikipedia — "GC content of sgRNA should optimally be over 50%."; chopchop.py calcGCContent.
    /// Score: clamp(100 − 300, 0, 100) = 0 (oracle: chopchop.py scoring lines, see TestSpec).
    /// </summary>
    [Test]
    public void EvaluateGuideRna_LowGcContent_LowerScore()
    {
        string guide = "AAAAAAAAAAAAAAAAAAAA"; // 0% GC
        var candidate = CrisprDesigner.EvaluateGuideRna(guide, CrisprSystemType.SpCas9);

        Assert.That(candidate.Score, Is.EqualTo(0));
        Assert.That(candidate.GcContent, Is.EqualTo(0));
        Assert.That(candidate.SeedGcContent, Is.EqualTo(0));
        Assert.That(candidate.HasPolyT, Is.False);
        Assert.That(candidate.Issues, Has.Count.EqualTo(2));
        Assert.That(candidate.Issues, Has.Some.Contains("Low GC"));
        Assert.That(candidate.Issues, Has.Some.Contains("Suboptimal seed region GC"));
    }

    /// <summary>
    /// M-003: 100% GC → CHOPCHOP's flat 300 penalty + 11 stems = 311 → Score 0; seed GC 100% is an issue only.
    /// </summary>
    [Test]
    public void EvaluateGuideRna_HighGcContent_LowerScore()
    {
        string guide = "GCGCGCGCGCGCGCGCGCGC"; // 100% GC, 11 CHOPCHOP stems
        var candidate = CrisprDesigner.EvaluateGuideRna(guide, CrisprSystemType.SpCas9);

        Assert.That(candidate.Score, Is.EqualTo(0)); // clamp(100 − (300 + 11), 0, 100)
        Assert.That(candidate.GcContent, Is.EqualTo(100));
        Assert.That(candidate.SeedGcContent, Is.EqualTo(100));
        Assert.That(candidate.HasPolyT, Is.False);
        Assert.That(candidate.Issues, Has.Count.EqualTo(3)); // + the stem issue
        Assert.That(candidate.Issues, Has.Some.Contains("High GC"));
        Assert.That(candidate.Issues, Has.Some.Contains("Suboptimal seed region GC"));
    }

    /// <summary>
    /// M-004: PolyT (TTTT) detected → reported as an issue; no score effect (CRISPOR warns only).
    /// Evidence: Addgene — "RNA polymerase III terminates at poly-T sequences."
    /// Input: "ACGTACGTTTTTACGTACGT" — 40% GC (no GC penalty), polyT present.
    /// Score: 100 − 5 (stems) = 95.
    /// </summary>
    [Test]
    public void EvaluateGuideRna_HasPolyT_Penalized()
    {
        string guide = "ACGTACGTTTTTACGTACGT"; // 40% GC, contains 5 consecutive T's, 5 CHOPCHOP stems
        var candidate = CrisprDesigner.EvaluateGuideRna(guide, CrisprSystemType.SpCas9);

        Assert.That(candidate.Score, Is.EqualTo(95)); // 100 − 5 (stems); TTTT is report-only
        Assert.That(candidate.GcContent, Is.EqualTo(40));
        Assert.That(candidate.HasPolyT, Is.True);
        Assert.That(candidate.Issues, Has.Count.EqualTo(2));
        Assert.That(candidate.Issues, Has.Some.Contains("TTTT"));
    }

    /// <summary>
    /// M-005: Empty guide should throw ArgumentNullException.
    /// Evidence: Defensive programming.
    /// </summary>
    [Test]
    public void EvaluateGuideRna_EmptyGuide_ThrowsArgumentException()
    {
        Assert.Throws<ArgumentException>(() =>
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
        Assert.That(candidate.FullGuideRna!.Length, Is.EqualTo(96));
    }

    // Round 4 H9/H41. Oracle: CRISPOR's constant primers (crispor.py L7576-7586): the SpCas9 tracrPrimer is the reverse
    // complement of scaffold + "TTTT"; the SaCas9 template is the scaffold + "TTTTTT" terminator.
    private const string CrisporSpCas9TracrPrimer =
        "AAAAGCACCGACTCGGTGCCACTTTTTCAAGTTGATAACGGACTAGCCTTATTTTAACTTGCTATTTCTAGCTCTAAAAC";
    private const string CrisporSaCas9TracrTemplate =
        "GTTTTAGTACTCTGGAAACAGAATCTACTAAAACAAGGCAAAATGCCGTGTTTATCTCGTCAACTTGTTGGCGAGATTTTTTT";

    [Test]
    public void SgRnaScaffolds_AreTheOnesCrisporDerivesItsConstantPrimersFrom()
    {
        Assert.Multiple(() =>
        {
            Assert.That(DnaSequence.GetReverseComplementString(CrisporSpCas9TracrPrimer),
                Is.EqualTo(CrisprDesigner.SpCas9SgRnaScaffold + "TTTT"));
            Assert.That(CrisprDesigner.SpCas9SgRnaScaffold, Has.Length.EqualTo(76));
            Assert.That(CrisporSaCas9TracrTemplate, Is.EqualTo(CrisprDesigner.SaCas9SgRnaScaffold + "TTTTTT"));
            Assert.That(CrisprDesigner.SaCas9SgRnaScaffold, Has.Length.EqualTo(77));
        });
    }

    /// <summary>
    /// H9: the scaffold follows the guide's own system. SpCas9 and its alternative-PAM variants share one protein,
    /// SaCas9 has a different tracr, and a Cas12a / CasX spacer has no scaffold in the reference tools (Cas12a's
    /// repeat is 5' of the spacer): <c>null</c>, not the SpCas9 scaffold.
    /// </summary>
    [TestCase(CrisprSystemType.SpCas9, 20, "SpCas9")]
    [TestCase(CrisprSystemType.SpCas9_NAG, 20, "SpCas9")]
    [TestCase(CrisprSystemType.SpCas9_NGA, 20, "SpCas9")]
    [TestCase(CrisprSystemType.SaCas9, 21, "SaCas9")]
    [TestCase(CrisprSystemType.SaCas9_NNGRRN, 21, "SaCas9")]
    [TestCase(CrisprSystemType.Cas12a, 23, null)]
    [TestCase(CrisprSystemType.AsCas12a, 23, null)]
    [TestCase(CrisprSystemType.LbCas12a, 23, null)]
    [TestCase(CrisprSystemType.CasX, 20, null)]
    public void FullGuideRna_FollowsTheGuidesOwnSystem(CrisprSystemType system, int length, string? scaffold)
    {
        string guide = new string('A', length - 8) + "CGTACGTA";

        string? full = CrisprDesigner.EvaluateGuideRna(guide, system).FullGuideRna;

        string? expected = scaffold switch
        {
            "SpCas9" => guide + CrisprDesigner.SpCas9SgRnaScaffold,
            "SaCas9" => guide + CrisprDesigner.SaCas9SgRnaScaffold,
            _ => null
        };
        Assert.That(full, Is.EqualTo(expected));
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
    /// S-001: Restriction site is reported as an issue (no score effect).
    /// Evidence: Common restriction sites interfere with cloning.
    /// Input: "ACGTGAATTCACGTACGTAC" — 45% GC, contains EcoRI site (GAATTC).
    /// Score: 100 − 3 (stems) = 97; the site is report-only.
    /// </summary>
    [Test]
    public void EvaluateGuideRna_RestrictionSite_Penalized()
    {
        string guide = "ACGTGAATTCACGTACGTAC"; // Contains GAATTC (EcoRI); 3 CHOPCHOP stems
        var candidate = CrisprDesigner.EvaluateGuideRna(guide, CrisprSystemType.SpCas9);

        Assert.That(candidate.Score, Is.EqualTo(97)); // 100 − 3 (stems); the site is report-only
        Assert.That(candidate.GcContent, Is.EqualTo(45));
        Assert.That(candidate.Issues, Has.Count.EqualTo(2));
        Assert.That(candidate.Issues[0], Does.Contain("restriction site"));
    }

    /// <summary>
    /// S-002: Seed GC content is calculated from last 10bp.
    /// Evidence: Addgene — seed region (8-10bp at 3') initiates annealing; implementation uses 10bp upper bound.
    /// Input: "AAAAAAAAAAAAACGTACGT" — overall GC 20%, seed (last 10) GC 40%.
    /// Score: 0 (GC 20% is outside the window: CHOPCHOP's flat 300 penalty; seed in [30,80], no seed issue).
    /// </summary>
    [Test]
    public void EvaluateGuideRna_CalculatesSeedGc()
    {
        string guide = "AAAAAAAAAAAAACGTACGT"; // Overall GC=20%, seed last 10 = "AACGTACGT" → 40% GC
        var candidate = CrisprDesigner.EvaluateGuideRna(guide, CrisprSystemType.SpCas9);

        Assert.That(candidate.SeedGcContent, Is.EqualTo(40));
        Assert.That(candidate.GcContent, Is.EqualTo(20));
        Assert.That(candidate.Score, Is.EqualTo(0));
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
        string guide = "GCGCGCGCGCGCGCAAAAAA"; // 14 G/C = 70% GC; 5 CHOPCHOP stems
        var candidate = CrisprDesigner.EvaluateGuideRna(guide, CrisprSystemType.SpCas9);

        Assert.That(candidate.GcContent, Is.EqualTo(70));
        Assert.That(candidate.SeedGcContent, Is.EqualTo(40));
        Assert.That(candidate.Score, Is.EqualTo(95)); // no GC deduction at the boundary; 100 − 5 (stems)
        Assert.That(candidate.Issues, Has.Count.EqualTo(1));
        Assert.That(candidate.Issues[0], Does.Contain("5 self-complementary 4-bp stems"));
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
        string guide = "ACGTACGTACGATTTTACGT"; // Exactly 4 consecutive T's (positions 12-15); 4 CHOPCHOP stems
        var candidate = CrisprDesigner.EvaluateGuideRna(guide, CrisprSystemType.SpCas9);

        Assert.That(candidate.HasPolyT, Is.True);
        Assert.That(candidate.GcContent, Is.EqualTo(40));
        Assert.That(candidate.Score, Is.EqualTo(96)); // 100 − 4 (stems); TTTT is report-only
        Assert.That(candidate.Issues, Has.Count.EqualTo(2));
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
        string guide = "ACGTACGTACGTACGTTTAC"; // 3 consecutive T's (positions 15-17); 6 CHOPCHOP stems
        var candidate = CrisprDesigner.EvaluateGuideRna(guide, CrisprSystemType.SpCas9);

        Assert.That(candidate.HasPolyT, Is.False);
        Assert.That(candidate.GcContent, Is.EqualTo(45));
        Assert.That(candidate.Score, Is.EqualTo(94)); // no polyT / Graf deduction; 100 − 6 (stems)
        Assert.That(candidate.Issues, Has.Count.EqualTo(2));
        Assert.That(candidate.Issues[0], Does.Contain("Graf 2019 TT-motif"));
        Assert.That(candidate.GrafMotif, Is.EqualTo(GrafMotifType.TtMotif));
    }

    /// <summary>
    /// S-010: Suboptimal seed GC (low) → issue only (no score effect).
    /// Evidence: Seed region (last 10bp) GC outside 30-80% → penalty.
    /// Input: "GCGCGCGCGCAAAAAAAAAA" — overall GC 50% (no GC penalty), seed GC 0% (all A's) → −5.
    /// Score: 100 − 1 (stem) = 99.
    /// </summary>
    [Test]
    public void EvaluateGuideRna_SeedGcLow_Penalized()
    {
        string guide = "GCGCGCGCGCAAAAAAAAAA"; // Overall 50% GC, seed last 10 = all A → 0% GC; 1 CHOPCHOP stem
        var candidate = CrisprDesigner.EvaluateGuideRna(guide, CrisprSystemType.SpCas9);

        Assert.That(candidate.GcContent, Is.EqualTo(50));
        Assert.That(candidate.SeedGcContent, Is.EqualTo(0));
        Assert.That(candidate.Score, Is.EqualTo(99)); // 100 − 1 (stem); seed GC is report-only
        Assert.That(candidate.Issues, Has.Count.EqualTo(2));
        Assert.That(candidate.Issues[0], Does.Contain("Suboptimal seed region GC"));
    }

    /// <summary>
    /// S-011: Suboptimal seed GC (high) → issue only (no score effect).
    /// Evidence: Seed region (last 10bp) GC outside 30-80% → penalty.
    /// Input: "AAAAAAAAAAGGGGGGGGGG" — overall GC 50%, seed GC 100% (all G) → −5.
    /// Score: 100 (0 stems).
    /// </summary>
    [Test]
    public void EvaluateGuideRna_SeedGcHigh_Penalized()
    {
        string guide = "AAAAAAAAAAGGGGGGGGGG"; // Overall 50% GC, seed last 10 = all G → 100% GC
        var candidate = CrisprDesigner.EvaluateGuideRna(guide, CrisprSystemType.SpCas9);

        Assert.That(candidate.GcContent, Is.EqualTo(50));
        Assert.That(candidate.SeedGcContent, Is.EqualTo(100));
        Assert.That(candidate.Score, Is.EqualTo(100));
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
        Assert.That(guides[0].Score, Is.EqualTo(90), "ACGT×5 carries 10 CHOPCHOP stems");
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
    /// C-001: each CHOPCHOP self-complementary 4-bp stem deducts 1 point (<c>SCORE['FOLDING'] = 1</c>).
    /// Evidence: chopchop.py <c>Guide.calcSelfComplementarity</c> ends in
    /// <c>self.score += self.folding * SCORE['FOLDING']</c> (E3). Oracle = the reference function run verbatim
    /// in Python: "GCGCGCGCGCGCGCGCGCGC" → 11 stems; "ACGGACTAGCCTACGTACGT" → 0 stems.
    /// Score: 100 − 60 (GC 100 %: (100−70)×2) − 11 (stems) − 5 (seed GC 100 % &gt; 80) = 24.
    /// Control: the 8-mer palindrome "GCGCGCGC" has 0 stems in the reference (the stem window must fit
    /// inside rvs[0:L−i], L = 8−4−1 = 3), so the old pair fraction's 0.3125 no longer penalises it.
    /// </summary>
    [Test]
    public void EvaluateGuideRna_SelfComplementary_PenaltyTriggered()
    {
        var result = CrisprDesigner.EvaluateGuideRna("GCGCGCGCGCGCGCGCGCGC", CrisprSystemType.SpCas9);

        Assert.That(result.SelfComplementaryStems, Is.EqualTo(11));
        // chopchop.py: 11 stems x SCORE['FOLDING'] (1) + SCORE['CRISPR_BAD_GC'] (300, GC 100%) = 311 -> clamp -> 0
        Assert.That(result.Score, Is.EqualTo(0));
        var gcTermOff = CrisprDesigner.EvaluateGuideRna("GCGCGCGCGCGCGCGCGCGC", CrisprSystemType.SpCas9,
            GuideRnaParameters.Default with { ScoreGcContent = false });
        Assert.That(gcTermOff.Score, Is.EqualTo(100 - 11 * CrisprDesigner.SelfComplementarityStemPenalty));
        Assert.That(gcTermOff.Score, Is.EqualTo(89), "CHOPCHOP default (GC term off): 100 - 11 stems");
        Assert.That(result.Issues, Has.Some.Contains("11 self-complementary 4-bp stems"));

        // The 8-mer is no SpCas9 guide (EvaluateGuideRna enforces the system's 20 nt, H17), so the stem count and the
        // window fraction are checked on the functions that take any length.
        Assert.That(CrisprDesigner.CountSelfComplementaryStems("GCGCGCGC"), Is.EqualTo(0));
        Assert.That(CrisprDesigner.SelfComplementarityWindowFraction(8, 0), Is.EqualTo(0.0), "0 stems / (8 - 4) windows (G13)");
        var control = CrisprDesigner.EvaluateGuideRna("ACGGACTAGCCTACGTACGT", CrisprSystemType.SpCas9);
        Assert.That(control.SelfComplementaryStems, Is.EqualTo(0));
        Assert.That(control.Score, Is.EqualTo(100), "55% GC, 0 stems: nothing is deducted");
        Assert.That(control.Issues.Any(i => i.Contains("self-complementary")), Is.False);
    }

    /// <summary>
    /// C-002: All-T guide — GC 0% is outside the window (CHOPCHOP flat 300 → clamped to 0); polyT and
    /// seed GC are reported as issues only.
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
        Assert.That(guides.All(g => g.Score == 90), Is.True, "ACGT×5 carries 10 CHOPCHOP stems");
    }

    /// <summary>
    /// C-006: SaCas9 system type evaluates correctly.
    /// Evidence: SaCas9 (Staphylococcus aureus) — NNGRRT PAM, 21bp guide, PAM after target.
    /// A 20-nt spacer is refused for SaCas9 (H17: the system's guide length is enforced); the 21-nt spacer below has 11
    /// CHOPCHOP stems and 10/21 GC (chopchop_sc.py, a line-by-line port of CHOPCHOP's selfComp / GC).
    /// </summary>
    [Test]
    public void EvaluateGuideRna_SaCas9SystemType_ValidEvaluation()
    {
        string guide = "ACGTACGTACGTACGTACGTA";
        var candidate = CrisprDesigner.EvaluateGuideRna(guide, CrisprSystemType.SaCas9);

        Assert.That(candidate.Sequence, Is.EqualTo(guide));
        Assert.That(candidate.Score, Is.EqualTo(89), "ACGT×5+A carries 11 CHOPCHOP stems");
        Assert.That(candidate.GcContent, Is.EqualTo(100.0 * 10 / 21).Within(1e-12));
        Assert.That(candidate.System.Name, Is.EqualTo("SaCas9"));
        Assert.That(candidate.System.GuideLength, Is.EqualTo(21));
        Assert.That(candidate.Issues, Has.Count.EqualTo(1));
        Assert.That(candidate.Issues[0], Does.Contain("11 self-complementary 4-bp stems"));
        Assert.Throws<ArgumentException>(() =>
            CrisprDesigner.EvaluateGuideRna("ACGTACGTACGTACGTACGT", CrisprSystemType.SaCas9));
    }

    /// <summary>
    /// C-007: GC just below 40% boundary triggers Low GC issue.
    /// Input: "AAAAAAAAAAAAGCGCGCAT" — 6 G/C of 20 = 30% GC.
    /// Score: 0 (outside the window → flat 300).
    /// </summary>
    [Test]
    public void EvaluateGuideRna_BelowBoundaryGc_HasLowGcIssue()
    {
        string guide = "AAAAAAAAAAAAGCGCGCAT"; // 6 G/C = 30% GC
        var candidate = CrisprDesigner.EvaluateGuideRna(guide, CrisprSystemType.SpCas9);

        Assert.That(candidate.GcContent, Is.EqualTo(30));
        Assert.That(candidate.Score, Is.EqualTo(0));
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
        Assert.That(guides[0].Score, Is.EqualTo(90), "ACGT×5 carries 10 CHOPCHOP stems");
    }

    /// <summary>
    /// C-009: GC above 70% boundary triggers High GC issue.
    /// Input: "GCGCGCGCGCGCGCGCAAAA" — 16 G/C of 20 = 80% GC.
    /// Score: 0 (outside the window → flat 300, plus 7 stems).
    /// </summary>
    [Test]
    public void EvaluateGuideRna_AboveBoundaryGc_HasHighGcIssue()
    {
        string guide = "GCGCGCGCGCGCGCGCAAAA"; // 16 G/C = 80% GC; 7 CHOPCHOP stems
        var candidate = CrisprDesigner.EvaluateGuideRna(guide, CrisprSystemType.SpCas9);

        Assert.That(candidate.GcContent, Is.EqualTo(80));
        Assert.That(candidate.Score, Is.EqualTo(0)); // clamp(100 − (300 + 7), 0, 100)
        Assert.That(candidate.Issues, Has.Count.EqualTo(2));
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
    /// Input plants four PAMs whose protospacers score 100 (50% GC, 0 stems), 95 (50% GC, 5 stems), 90 (10 stems)
    /// and 0 (all-T, GC outside the window). Round 4 H8: the GC window is an independent filter, so the all-T
    /// guide is no longer listed (it used to be, as a trailing score-0 row, with MinScore = 0). The stem counts are
    /// CHOPCHOP's own <c>selfComp</c> (a line-by-line port, scratch oracle chopchop_sc.py): 0, 5 and 10.
    /// </summary>
    [Test]
    public void DesignGuideRnas_Output_IsRankedByScoreDescending()
    {
        var sequence = new DnaSequence(
            "CGTAACGTGCTTGCAACTGAAGG" +      // 50% GC, 0 stems -> 100 (no CC: no reverse-strand PAM)
            "AGTCTCATCTCACTTGAGACGCGG" +     // 50% GC, 5 stems -> 95 (leading A: no NGG across the previous PAM)
            "TTTTTTTTTTTTTTTTTTTTCGG" +      // all-T  -> GC outside the window: filtered
            "ACGTACGTACGTACGTACGTTGG");      // 50% GC, 10 CHOPCHOP stems -> 90

        var guides = CrisprDesigner.DesignGuideRnas(
            sequence, 0, sequence.Length - 1, CrisprSystemType.SpCas9,
            GuideRnaParameters.Default with { MinScore = 0 }).ToList();

        Assert.That(guides.Select(g => g.Score).ToArray(), Is.EqualTo(new[] { 100.0, 95.0, 90.0 }));
        Assert.That(guides[0].Sequence, Is.EqualTo("CGTAACGTGCTTGCAACTGA"));
        Assert.That(guides[1].Sequence, Is.EqualTo("GTCTCATCTCACTTGAGACG"));
        Assert.That(guides.Any(g => g.Sequence == "TTTTTTTTTTTTTTTTTTTT"), Is.False,
            "a guide outside the GC window is excluded, not listed with score 0");
    }

    /// <summary>
    /// S-019: equal ranking keys are broken by ascending forward-strand <see cref="GuideRnaCandidate.Position"/>
    /// (then forward strand first), in both ranking modes — the deterministic tie-break documented on
    /// <see cref="CrisprDesigner.DesignGuideRnas"/> (B08 F5). The PAM scan yields forward-strand sites first,
    /// so without the tie-break the forward guide (position 23) would precede the reverse guide (position 3).
    /// Both protospacers are "ACGTACGTACGTACGTACGT" (equal composition score 100) and neither has a full
    /// 30-nt Rule Set 2 window (both on-target scores null), so the keys tie in both modes.
    /// </summary>
    [TestCase(GuideRnaRanking.Score)]
    [TestCase(GuideRnaRanking.OnTargetRuleSet2)]
    public void DesignGuideRnas_EqualRankingKeys_AreOrderedByPositionDeterministically(GuideRnaRanking ranking)
    {
        // Reverse-strand TGG (forward "CCA") at 0 -> protospacer at 3; forward AGG at 43 -> protospacer at 23.
        var sequence = new DnaSequence("CCA" + "ACGTACGTACGTACGTACGT" + "ACGTACGTACGTACGTACGT" + "AGG");

        var guides = CrisprDesigner.DesignGuideRnas(
            sequence, 0, sequence.Length - 1, CrisprSystemType.SpCas9, NoFilter with { Ranking = ranking }).ToList();

        Assert.That(guides, Has.Count.EqualTo(2));
        Assert.That(guides.Select(g => g.Score), Is.All.EqualTo(100.0));
        Assert.That(guides.Select(g => g.OnTargetScore), Is.All.Null);
        Assert.That(guides.Select(g => (g.Position, g.IsForwardStrand)),
            Is.EqualTo(new[] { (3, false), (23, true) }));
    }

    /// <summary>
    /// S-016: designed candidates keep the metadata of the system they were designed with. A
    /// name-based remap previously replaced LbCas12a/AsCas12a/CasX/SpCas9-NAG with SpCas9, which also
    /// flipped the seed-region orientation (<c>PamAfterTarget</c>) and therefore changed the score:
    /// the 23-nt LbCas12a protospacer "CGTACGTACGTACGTACGTACGT" was seeded from its last 10 bases
    /// ("GTACGTACGT" -> GC 50%) instead of its PAM-proximal first 10 ("CGTACGTACG" -> GC 60%).
    /// </summary>
    [Test]
    public void DesignGuideRnas_KeepsTheRequestedSystemMetadata()
    {
        var sequence = new DnaSequence("TTTACGTACGTACGTACGTACGTACGTACGTACGTACGT");

        var lb = CrisprDesigner.DesignGuideRnas(
            sequence, 0, sequence.Length - 1, CrisprSystemType.LbCas12a, NoFilter).ToList();

        Assert.That(lb, Has.Count.EqualTo(1));
        Assert.That(lb[0].System.Name, Is.EqualTo("LbCas12a"));
        Assert.That(lb[0].System.GuideLength, Is.EqualTo(23));
        Assert.That(lb[0].System.PamAfterTarget, Is.False);
        Assert.That(lb[0].Sequence, Has.Length.EqualTo(23));
        // Seed = the 10 PAM-proximal (5') bases "CGTACGTACG" -> 6 G/C of 10 = 60%.
        Assert.That(lb[0].SeedGcContent, Is.EqualTo(60));
    }

    /// <summary>
    /// S-017: <see cref="GuideRnaParameters.AvoidPolyT"/> = <c>false</c> disables the TTTT issue (TTTT never
    /// affected the score after B08 round 3 G5; CRISPOR only warns) while the measurement <see cref="GuideRnaCandidate.HasPolyT"/> is still reported.
    /// Repro of the defect: the flag was stored but never read, so M-004's guide still scored 80.
    /// </summary>
    [Test]
    public void EvaluateGuideRna_AvoidPolyTFalse_DisablesThePolyTIssue_ScoreUnchanged()
    {
        const string guide = "ACGTACGTTTTTACGTACGT"; // 40% GC, contains TTTTT

        var on = CrisprDesigner.EvaluateGuideRna(guide, CrisprSystemType.SpCas9,
            GuideRnaParameters.Default with { AvoidPolyT = true });
        var off = CrisprDesigner.EvaluateGuideRna(guide, CrisprSystemType.SpCas9,
            GuideRnaParameters.Default with { AvoidPolyT = false });

        Assert.That(on.Score, Is.EqualTo(95)); // 100 − 5 (stems); TTTT is report-only
        Assert.That(off.Score, Is.EqualTo(95)); // identical: the flag only controls the issue
        Assert.That(on.Issues, Has.Some.Contains("TTTT"));
        Assert.That(off.HasPolyT, Is.True);
        Assert.That(off.Issues, Has.None.Contains("TTTT"));
    }

    /// <summary>
    /// S-018: <see cref="GuideRnaParameters.CheckSelfComplementarity"/> = <c>false</c> disables the
    /// self-complementarity penalty and its issue while the measurement is still reported.
    /// </summary>
    [Test]
    public void EvaluateGuideRna_CheckSelfComplementarityFalse_DisablesThePenalty()
    {
        const string guide = "GCGCGCGCGCGCGCGCGCGC"; // 11 CHOPCHOP stems, seed GC 100%
        var wide = new GuideRnaParameters(0, 100, 0, AvoidPolyT: true, CheckSelfComplementarity: true);

        var on = CrisprDesigner.EvaluateGuideRna(guide, CrisprSystemType.SpCas9, wide);
        var off = CrisprDesigner.EvaluateGuideRna(guide, CrisprSystemType.SpCas9,
            wide with { CheckSelfComplementarity = false });

        Assert.That(on.Score, Is.EqualTo(100 - 11));
        Assert.That(off.Score, Is.EqualTo(100));
        Assert.That(off.SelfComplementaryStems, Is.EqualTo(11), "the measurement is still reported");
        Assert.That(off.SelfComplementarityScore, Is.EqualTo(11.0 / 16.0), "11 stems / (20 - 4) windows, still reported (G13)");
        Assert.That(off.Issues, Has.Count.EqualTo(1));
        Assert.That(off.Issues[0], Does.Contain("Suboptimal seed region GC"));
    }

    #endregion

    #region 2026-09 review — Doench 2016 Rule Set 2 on designed guides

    // Oracle row from tests/.../TestData/Azimuth/nopos_oracle.csv (the verified reference prediction
    // and Microsoft Azimuth's own regression fixture `truth nopos` concur on every row):
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
        Assert.That(spCas9.Score, Is.EqualTo(100 - 8), "CRISPOR warns about the Graf motifs, it does not score them; 8 CHOPCHOP stems");
        Assert.That(spCas9.Issues, Has.Some.Contains("Graf 2019 TT-motif"));

        var saCas9 = CrisprDesigner.EvaluateGuideRna("A" + guide, CrisprSystemType.SaCas9); // 21 nt, same 3' end
        Assert.That(saCas9.GrafMotif, Is.EqualTo(GrafMotifType.None));
        Assert.That(saCas9.Issues, Has.None.Contains("Graf"));
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
        Assert.That(flagged.Score, Is.EqualTo(90), "the limit does not deduct; only the 10 stems do (1 point each)");
        Assert.That(flagged.Issues, Has.Some.Contains("(limit 2)"));
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
        Assert.Throws<ArgumentException>(() => CrisprDesigner.CountSelfComplementaryStems(""));
        Assert.Throws<ArgumentOutOfRangeException>(
            () => CrisprDesigner.CountSelfComplementaryStems("ACGTACGTACGT", 0));
    }

    #endregion

    #region 2026-10 review (B08 round 2, E3) — CHOPCHOP stem penalty, raw-string alphabet, ASCII fold

    /// <summary>
    /// Oracle: <c>Guide.calcSelfComplementarity</c> of CHOPCHOP's chopchop.py (sha256 3c7a4344…), transcribed
    /// verbatim in Python and run under Biopython 1.88, with the default empty backbone and with the standard
    /// backbone "AGGCTAGTCCGT". The deduction per stem is <c>SCORE['FOLDING'] = 1</c>.
    /// </summary>
    [TestCase("ACGTACGTACGTACGTACGT", 10, 10)]
    [TestCase("GCGCGCGCGCGCGCGCGCGC", 11, 11)]
    [TestCase("ACGGACTAGCCTACGTACGT", 0, 8)]
    public void EvaluateGuideRna_StemPenalty_IsOnePointPerChopchopStem(string guide, int stems, int stemsWithBackbone)
    {
        var wide = new GuideRnaParameters(0, 100, 0, AvoidPolyT: false, CheckSelfComplementarity: true);
        var off = wide with { CheckSelfComplementarity = false };
        var withBackbone = wide with
        {
            SelfComplementarityBackboneRegions = new[] { CrisprDesigner.StandardSgRnaBackboneRegion },
        };

        var on = CrisprDesigner.EvaluateGuideRna(guide, CrisprSystemType.SpCas9, wide);
        var none = CrisprDesigner.EvaluateGuideRna(guide, CrisprSystemType.SpCas9, off);
        var bb = CrisprDesigner.EvaluateGuideRna(guide, CrisprSystemType.SpCas9, withBackbone);

        Assert.That(on.SelfComplementaryStems, Is.EqualTo(stems));
        Assert.That(bb.SelfComplementaryStems, Is.EqualTo(stemsWithBackbone));
        Assert.That(on.Score, Is.EqualTo(none.Score - stems));
        Assert.That(bb.Score, Is.EqualTo(none.Score - stemsWithBackbone), "the backbone stems are penalised too");
    }

    /// <summary>
    /// B08 round 3 G5 oracle. Expected values are the output of CHOPCHOP's own scoring lines
    /// (<c>chopchop.py</c> sha256 3c7a4344...0d5960: <c>calcGCContent</c> body, <c>selfComp</c>, <c>SCORE</c>,
    /// <c>GC_LOW = 40</c>, <c>GC_HIGH = 70</c>) run verbatim in Python 3 + Biopython 1.88 on each guide:
    /// <c>penalty = stems x 1 + (300 if GC &gt; 70 or GC &lt; 40)</c>; Score = clamp(100 - penalty, 0, 100).
    /// Columns: guide, GC %, stems, penalty with the GC term on, penalty with it off.
    /// </summary>
    [TestCase("AGTAGAATTAGACAAAATAA", 20.0, 0, 300, 0)]
    [TestCase("ATATTCTAAATGTCTTCAAC", 25.0, 0, 300, 0)]
    [TestCase("CTCACTCATCACGTATAAAA", 35.0, 0, 300, 0)]
    [TestCase("TCCGATCCGTCTTATAATAG", 40.0, 0, 0, 0)]       // lower edge inclusive
    [TestCase("AAATCATTGCCTCACTACGC", 45.0, 0, 0, 0)]
    [TestCase("ACGTACGTACGTACGTACGT", 50.0, 10, 10, 10)]
    [TestCase("GGACAAGGGGGTGATGACGA", 60.0, 0, 0, 0)]
    [TestCase("GATCCGCGCGAGACCCTGAG", 70.0, 0, 0, 0)]       // upper edge inclusive
    [TestCase("CTCGGGCCTGGCGGGGTAAG", 75.0, 0, 300, 0)]
    [TestCase("GTGCGATGCGCCCCCGACGC", 80.0, 0, 300, 0)]
    [TestCase("GCGCGCGCGCGCGCGCGCGC", 100.0, 11, 311, 11)]
    [TestCase("TTTTTTTTTTTTTTTTTTTT", 0.0, 0, 300, 0)]      // TTTT: report-only
    [TestCase("ACGTACGTTTTTACGTACGT", 40.0, 5, 5, 5)]       // TTTT inside the window: report-only
    [TestCase("AAAAAAAAAAGGGGGGGGGG", 50.0, 0, 0, 0)]       // seed GC 100 %: report-only
    [TestCase("GGGGGGGGGGAAAAAAAAAA", 50.0, 0, 0, 0)]       // seed GC 0 %: report-only
    [TestCase("GAATTCGACGTCAGCGATCG", 55.0, 0, 0, 0)]       // EcoRI site: report-only
    public void EvaluateGuideRna_Score_MatchesChopchopPenaltyOracle(
        string guide, double gc, int stems, int penaltyGcOn, int penaltyGcOff)
    {
        var on = CrisprDesigner.EvaluateGuideRna(guide, CrisprSystemType.SpCas9);
        var off = CrisprDesigner.EvaluateGuideRna(guide, CrisprSystemType.SpCas9,
            GuideRnaParameters.Default with { ScoreGcContent = false });

        Assert.That(on.GcContent, Is.EqualTo(gc).Within(1e-9));
        Assert.That(on.SelfComplementaryStems, Is.EqualTo(stems));
        Assert.That(on.Score, Is.EqualTo(Math.Clamp(100.0 - penaltyGcOn, 0, 100)));
        Assert.That(off.Score, Is.EqualTo(Math.Clamp(100.0 - penaltyGcOff, 0, 100)));
        Assert.That(on.Issues.Any(i => i.Contains("GC content")), Is.EqualTo(gc < 40 || gc > 70),
            "the GC window is reported as an issue whatever the score switch says");
        Assert.That(off.Issues.Any(i => i.Contains("GC content")), Is.EqualTo(gc < 40 || gc > 70));
    }

    /// <summary>
    /// The composition flags that CHOPCHOP does not score (TTTT, seed GC, cloning site) change the issue list
    /// only, never the score.
    /// </summary>
    [Test]
    public void EvaluateGuideRna_TtttSeedGcCloningSite_AreIssuesOnly()
    {
        var plain = CrisprDesigner.EvaluateGuideRna("AAATCATTGCCTCACTACGC", CrisprSystemType.SpCas9);
        var tttt = CrisprDesigner.EvaluateGuideRna("ACGTACGTTTTTACGTACGT", CrisprSystemType.SpCas9);
        var seed = CrisprDesigner.EvaluateGuideRna("AAAAAAAAAAGGGGGGGGGG", CrisprSystemType.SpCas9);
        var site = CrisprDesigner.EvaluateGuideRna("GAATTCGACGTCAGCGATCG", CrisprSystemType.SpCas9);

        Assert.That(plain.Score, Is.EqualTo(100));
        Assert.That(plain.Issues, Is.Empty);
        Assert.That(tttt.Issues, Has.Some.Contains("TTTT"));
        Assert.That(tttt.Score, Is.EqualTo(100 - tttt.SelfComplementaryStems));
        Assert.That(seed.Issues, Has.Some.Contains("Suboptimal seed region GC"));
        Assert.That(seed.Score, Is.EqualTo(100));
        Assert.That(site.Issues, Has.Some.Contains("restriction site"));
        Assert.That(site.Score, Is.EqualTo(100));
    }

    /// <summary>
    /// <see cref="GuideRnaParameters.ScoreGcContent"/> defaults to <c>true</c> on every construction path; a custom
    /// window moves the edges but the penalty stays CHOPCHOP's flat 300.
    /// </summary>
    [Test]
    public void GuideRnaParameters_ScoreGcContent_DefaultAndCustomWindow()
    {
        Assert.That(GuideRnaParameters.Default.ScoreGcContent, Is.True);
        Assert.That(new GuideRnaParameters(40, 70, 50, true, true).ScoreGcContent, Is.True);
        Assert.That(CrisprDesigner.BadGcContentPenalty, Is.EqualTo(300.0));

        // 80 % GC, 0 stems: outside the default window, inside a widened one.
        var wide = new GuideRnaParameters(30, 85, 50, true, true);
        Assert.That(CrisprDesigner.EvaluateGuideRna("GTGCGATGCGCCCCCGACGC", CrisprSystemType.SpCas9).Score, Is.EqualTo(0));
        Assert.That(CrisprDesigner.EvaluateGuideRna("GTGCGATGCGCCCCCGACGC", CrisprSystemType.SpCas9, wide).Score, Is.EqualTo(100));
    }

    /// <summary>
    /// G13 oracle: stems from <c>chopchop.py selfComp</c> (verbatim, Python 3 + Biopython 1.88) and the window count
    /// <c>len - 4</c> of its <c>range(0, len(fwd) - STEM_LEN)</c> loop; the expected fraction is stems / windows.
    /// </summary>
    [TestCase("ACGTACGTACGTACGTACGT", 10, 16, 0.625)]
    [TestCase("GCGCGCGCGCGCGCGCGCGC", 11, 16, 0.6875)]
    [TestCase("ACGGACTAGCCTACGTACGT", 0, 16, 0.0)]
    public void EvaluateGuideRna_SelfComplementarityScore_MatchesOracleStemFraction(
        string guide, int stems, int windows, double expected)
    {
        var c = CrisprDesigner.EvaluateGuideRna(guide, CrisprSystemType.SpCas9);

        Assert.That(c.SelfComplementaryStems, Is.EqualTo(stems));
        Assert.That(guide.Length - 4 > 0 ? guide.Length - 4 : 0, Is.EqualTo(windows));
        Assert.That(c.SelfComplementarityScore, Is.EqualTo(expected).Within(1e-12));
        Assert.That(c.SelfComplementarityScore, Is.InRange(0.0, 1.0));
    }

    /// <summary>The window fraction itself, for lengths no system's guide has (the same G13 oracle rows).</summary>
    [TestCase("ACGTACGTACGT", 2, 8, 0.25)]
    [TestCase("GCGCGCGC", 0, 4, 0.0)]
    [TestCase("ACGT", 0, 0, 0.0)]                      // no 4-bp window fits -> 0, never NaN
    public void SelfComplementarityWindowFraction_MatchesOracleForAnyLength(
        string guide, int stems, int windows, double expected)
    {
        Assert.That(CrisprDesigner.CountSelfComplementaryStems(guide), Is.EqualTo(stems));
        Assert.That(CrisprDesigner.SelfComplementarityWindowFraction(guide.Length, stems), Is.EqualTo(expected).Within(1e-12));
        Assert.That(guide.Length - 4 > 0 ? guide.Length - 4 : 0, Is.EqualTo(windows));
    }

    /// <summary>The backbone-aware stem count feeds the fraction too (oracle: `ACGGACTAGCCTACGTACGT` 8 stems / 16 windows).</summary>
    [Test]
    public void EvaluateGuideRna_SelfComplementarityScore_FollowsBackboneStems()
    {
        var c = CrisprDesigner.EvaluateGuideRna("ACGGACTAGCCTACGTACGT", CrisprSystemType.SpCas9,
            GuideRnaParameters.Default with
            {
                SelfComplementarityBackboneRegions = new[] { CrisprDesigner.StandardSgRnaBackboneRegion },
            });

        Assert.That(c.SelfComplementaryStems, Is.EqualTo(8));
        Assert.That(c.SelfComplementarityScore, Is.EqualTo(0.5));
    }

    /// <summary>
    /// B08 round 3 G13: <see cref="GuideRnaCandidate.SelfComplementarityScore"/> is CHOPCHOP's stem count over its
    /// own loop bound (<c>len - 4</c> windows), not the former unsourced pair fraction. The 8-mer palindrome
    /// <c>GCGCGCGC</c> had 0.3125 under the old fraction and has 0 stems in the reference.
    /// </summary>
    [Test]
    public void EvaluateGuideRna_SelfComplementarityScore_IsTheChopchopStemFraction_AndDoesNotEnterTheScore()
    {
        // 0 stems: nothing is deducted. (The 8-mer palindrome of the old test is no SpCas9 guide any more; its 0 stems and
        // 0.0 fraction are pinned by SelfComplementarityWindowFraction_MatchesOracleForAnyLength.)
        var c = CrisprDesigner.EvaluateGuideRna("ACGGACTAGCCTACGTACGT", CrisprSystemType.SpCas9,
            new GuideRnaParameters(0, 100, 0, AvoidPolyT: true, CheckSelfComplementarity: true));

        Assert.That(c.SelfComplementarityScore, Is.EqualTo(0.0));
        Assert.That(c.SelfComplementaryStems, Is.EqualTo(0));
        Assert.That(c.Score, Is.EqualTo(100), "0 stems: nothing is deducted (the seed-GC flag is report-only)");
    }

    [TestCase("XXXXXXXXXXXXXXXXXXXX")]
    [TestCase("ACGTACGTACGTACGTACGN")]   // IUPAC ambiguity codes are not unambiguous nucleotides
    [TestCase("ACGTACGTACGTACGTACG-")]
    [TestCase("ACGTACGTACGT ACGTACG")]
    [TestCase("ACGTACGTACGTACGTACGſ")]   // U+017F must not be folded onto IUPAC 'S' and must not pass
    public void RawStringGuideEntryPoints_RejectNonNucleotides(string guide)
    {
        var e1 = Assert.Throws<ArgumentException>(() => CrisprDesigner.EvaluateGuideRna(guide));
        Assert.That(e1!.ParamName, Is.EqualTo("guideSequence"));
        var e2 = Assert.Throws<ArgumentException>(() => CrisprDesigner.CountSelfComplementaryStems(guide));
        Assert.That(e2!.ParamName, Is.EqualTo("guideSequence"));
        Assert.Throws<ArgumentException>(
            () => CrisprDesigner.CountSelfComplementaryStems("ACGTACGTACGTACGTACGT", 4, new[] { guide }));
    }

    [Test]
    public void RawStringGuideEntryPoints_FoldCaseAndReadUAsT()
    {
        const string dna = "ACGGACTTTTCCTACGTACG";
        string rna = dna.Replace('T', 'U');

        var a = CrisprDesigner.EvaluateGuideRna(dna);
        var b = CrisprDesigner.EvaluateGuideRna(rna.ToLowerInvariant());

        Assert.That(b.Sequence, Is.EqualTo(dna));
        Assert.That(b.HasPolyT, Is.True, "UUUU is TTTT");
        Assert.That(b.Score, Is.EqualTo(a.Score));
        Assert.That(CrisprDesigner.CountSelfComplementaryStems(rna),
            Is.EqualTo(CrisprDesigner.CountSelfComplementaryStems(dna)));
        // the RNA scaffold is accepted as well
        Assert.That(
            CrisprDesigner.CountSelfComplementaryStems(dna, 4, new[] { "AGGCUAGUCCGU" }),
            Is.EqualTo(CrisprDesigner.CountSelfComplementaryStems(dna, 4, new[] { "AGGCTAGTCCGT" })));
    }

    /// <summary>
    /// The case fold of the raw-string scorers is ASCII-only (D20 / F49): <c>string.ToUpperInvariant</c> maps
    /// U+017F 'ſ' onto 'S', so a guide carrying 'ſ' used to compare equal to one carrying 'S'. Since G2 (B08
    /// round 3) neither 'S' nor 'ſ' is a nucleotide, so the hit score refuses both (consistency with F70/C-023);
    /// the fold is still locked: 'ſ' is never silently treated as 'S'.
    /// </summary>
    [Test]
    public void CalculateMitHitScore_LongS_IsRejectedNotFoldedOntoS()
    {
        const string guide = "ACGTACGTACGTACGTACGA";
        const string longS = "ACGTACGTACGTACGTACGſ";
        const string plainMismatch = "ACGTACGTACGTACGTACGC";

        Assert.That(CrisprDesigner.CalculateMitHitScore(guide, guide), Is.EqualTo(100.0), "control: identical");
        Assert.That(CrisprDesigner.CalculateMitHitScore(guide, plainMismatch), Is.LessThan(100.0), "control: mismatch");
        var e = Assert.Throws<ArgumentException>(() => CrisprDesigner.CalculateMitHitScore(guide, longS));
        Assert.That(e!.ParamName, Is.EqualTo("offTarget20"));
        Assert.Throws<ArgumentException>(() => CrisprDesigner.CalculateMitHitScore("ACGTACGTACGTACGTACGS", longS));
    }

    // ------------------------------------------------------------------------------------------
    // G2 (B08 round 3): every raw-string guide entry point shares F70's alphabet rule — ASCII fold,
    // U read as T, anything else is an ArgumentException naming the parameter. Before, the off-target
    // scan compared the U-spelling literally (the on-target became a 5-mismatch hit) and a guide of 'N'
    // (or any junk) mismatched everything, so the aggregates returned a false perfect 100.
    // ------------------------------------------------------------------------------------------

    private const string G2Guide = "ACGGACTTCCTACGTAGCAT";
    private const string G2Rna = "ACGGACUUCCUACGUAGCAU";

    /// <summary>Genome with the on-target (AGG), a 2-mismatch NGG site and a 4-mismatch NGG site.</summary>
    private static DnaSequence G2Genome() => new(
        "TTTT" + G2Guide + "AGG" +
        "CCCC" + "ACGCACTTCCTACGTAGCTT" + "TGG" +
        "GGGG" + "ACGGTCTTCCTACGTAGAAA" + "CGG" + "TTTT");

    [Test]
    public void FindOffTargets_RnaSpellingOfTheGuide_EqualsDnaSpelling_AndDropsTheOnTarget()
    {
        var genome = G2Genome();
        var dna = CrisprDesigner.FindOffTargets(G2Guide, genome, 4).ToList();
        var rna = CrisprDesigner.FindOffTargets(G2Rna, genome, 4).ToList();
        var lower = CrisprDesigner.FindOffTargets(G2Rna.ToLowerInvariant(), genome, 4).ToList();

        Assert.That(dna, Is.Not.Empty, "control: the fixture has off-targets");
        Assert.That(dna.All(h => h.Mismatches > 0), Is.True, "the on-target (0 mismatches) is never an off-target");
        Assert.That(Describe(rna), Is.EqualTo(Describe(dna)));
        Assert.That(Describe(lower), Is.EqualTo(Describe(dna)));
        Assert.That(
            Describe(CrisprDesigner.FindOffTargetsWithPam(G2Rna, genome, 4).Select(h => h.Site)),
            Is.EqualTo(Describe(CrisprDesigner.FindOffTargetsWithPam(G2Guide, genome, 4).Select(h => h.Site))));
        Assert.That(
            CrisprDesigner.FindOffTargetsWithPam(G2Rna, genome, 4).Select(h => h.Pam),
            Is.EqualTo(CrisprDesigner.FindOffTargetsWithPam(G2Guide, genome, 4).Select(h => h.Pam)));
    }

    // OffTargetSite carries a list, so record equality is by reference: compare by value.
    private static string[] Describe(IEnumerable<OffTargetSite> sites) => sites
        .Select(h => $"{h.Position}|{h.Sequence}|{h.Mismatches}|{string.Join(',', h.MismatchPositions)}|{h.IsForwardStrand}|{h.OffTargetScore}")
        .ToArray();

    [Test]
    public void SpecificityAggregates_RnaSpellingOfTheGuide_EqualsDnaSpelling()
    {
        var genome = G2Genome();

        Assert.That(CrisprDesigner.CalculateHeuristicSpecificityScore(G2Rna, genome),
            Is.EqualTo(CrisprDesigner.CalculateHeuristicSpecificityScore(G2Guide, genome)));
        Assert.That(CrisprDesigner.CalculateMitSpecificityScore(G2Rna, genome),
            Is.EqualTo(CrisprDesigner.CalculateMitSpecificityScore(G2Guide, genome)));
        Assert.That(CrisprDesigner.CalculateMitSpecificityScoreWithAlternativePams(G2Rna, genome),
            Is.EqualTo(CrisprDesigner.CalculateMitSpecificityScoreWithAlternativePams(G2Guide, genome)));
        Assert.That(CrisprDesigner.CalculateCfdSpecificityScore(G2Rna, genome),
            Is.EqualTo(CrisprDesigner.CalculateCfdSpecificityScore(G2Guide, genome)));
        Assert.That(CrisprDesigner.CalculateMitSpecificityScore(G2Guide, genome), Is.LessThan(100.0),
            "control: the fixture is not trivially perfect");

        var hits = CrisprDesigner.FindOffTargetsWithPam(G2Guide, genome, 4).ToList();
        Assert.That(CrisprDesigner.CalculateMitSpecificityScore(G2Rna, hits),
            Is.EqualTo(CrisprDesigner.CalculateMitSpecificityScore(G2Guide, hits)));
    }

    [TestCase("NNNNNNNNNNNNNNNNNNNN")]   // the false-perfect-100 probe
    [TestCase("XXXXXXXXXXXXXXXXXXXX")]
    [TestCase("ACGGACTTCCTACGTAGCA-")]
    [TestCase("ACGGACTTCCTACGTAGCAſ")]
    public void RawStringOffTargetEntryPoints_RejectNonNucleotides(string guide)
    {
        var genome = G2Genome();
        var hits = CrisprDesigner.FindOffTargetsWithPam(G2Guide, genome, 4).ToList();

        void Named(TestDelegate call)
        {
            var e = Assert.Throws<ArgumentException>(call);
            Assert.That(e!.ParamName, Is.EqualTo("guideSequence"));
        }

        Named(() => CrisprDesigner.FindOffTargets(guide, genome));
        Named(() => CrisprDesigner.FindOffTargetsWithPam(guide, genome));
        Named(() => CrisprDesigner.CalculateHeuristicSpecificityScore(guide, genome));
        Named(() => CrisprDesigner.CalculateSpecificityScore(guide, genome));
        Named(() => CrisprDesigner.CalculateMitSpecificityScore(guide, genome));
        Named(() => CrisprDesigner.CalculateMitSpecificityScore(guide, hits));
        Named(() => CrisprDesigner.CalculateMitSpecificityScoreWithAlternativePams(guide, genome));
        Named(() => CrisprDesigner.CalculateCfdSpecificityScore(guide, genome));
        Named(() => CrisprDesigner.GetGrafMotif(guide));
    }

    [Test]
    public void RawStringOffTargetEntryPoints_KeepTheirGuardOrder()
    {
        var genome = G2Genome();
        // null/empty first, then the genome, then the range, then the length, only then the alphabet.
        Assert.Throws<ArgumentException>(() => CrisprDesigner.FindOffTargets("", genome));
        Assert.Throws<ArgumentNullException>(() => CrisprDesigner.FindOffTargets("NNNN", null!));
        Assert.Throws<ArgumentOutOfRangeException>(() => CrisprDesigner.FindOffTargets("NNNN", genome, 6));
        var e = Assert.Throws<ArgumentException>(() => CrisprDesigner.FindOffTargets("NNNN", genome));
        Assert.That(e!.Message, Does.Contain("Guide length"), "the length guard precedes the alphabet guard");
    }

    [Test]
    public void GetGrafMotif_RnaSpellingAndCase_EqualDnaSpelling()
    {
        foreach (var guide in new[] { "ACGTACGTACGTACGTACTC", "ACGTACGTACGTACGTAGCC", "ACGTACGTACGTACGTATTC", "ACGTACGTACGTACGTCTTC" })
        {
            var expected = CrisprDesigner.GetGrafMotif(guide);
            Assert.That(CrisprDesigner.GetGrafMotif(guide.Replace('T', 'U')), Is.EqualTo(expected), guide);
            Assert.That(CrisprDesigner.GetGrafMotif(guide.Replace('T', 'U').ToLowerInvariant()), Is.EqualTo(expected), guide);
        }
        Assert.That(CrisprDesigner.GetGrafMotif("AAAAAAAAAAAAAAAAAUUC"), Is.EqualTo(GrafMotifType.TtMotif));
    }

    [Test]
    public void MitHitScore_RejectsNonNucleotides_AndReadsUAsT()
    {
        const string g = "ACGGACTTCCTACGTAGCAT";
        const string o = "ACGGACTTCCTACGTAGCTA";
        const string n20 = "NNNNNNNNNNNNNNNNNNNN";

        var e1 = Assert.Throws<ArgumentException>(() => CrisprDesigner.CalculateMitHitScore(n20, o));
        Assert.That(e1!.ParamName, Is.EqualTo("guide20"));
        var e2 = Assert.Throws<ArgumentException>(() => CrisprDesigner.CalculateMitHitScore(g, n20));
        Assert.That(e2!.ParamName, Is.EqualTo("offTarget20"));
        // CFD / SaCas9 hit scores already refused junk (and keep refusing U: CFD's reference screens [^ATCG]).
        Assert.Throws<ArgumentException>(() => CrisprDesigner.CalculateCfdScore(n20, o, "GG"));
        Assert.Throws<ArgumentException>(() => CrisprDesigner.CalculateCfdScore(g, n20, "GG"));
        Assert.Throws<ArgumentException>(() => CrisprDesigner.CalculateSaCas9HitScore(n20 + "N", o + "A"));

        Assert.That(CrisprDesigner.CalculateMitHitScore(g.Replace('T', 'U'), o.Replace('T', 'U')),
            Is.EqualTo(CrisprDesigner.CalculateMitHitScore(g, o)));
    }

    #endregion
    #region Round 4 H8 - the GC window is an independent filter (CHOPCHOP filterGCmin / filterGCmax)

    // Oracle: chopchop.py L1911-1914 (sha256 3c7a4344...): `gc = GC(dna[0:-len(PAM)]); if gc < filterGCmin or
    // gc > filterGCmax: add = False` - GC in percent of the protospacer, inclusive edges, independent of the
    // SCORE_GC term. Evaluated for a 40..70 window: the passing GC counts are 8..14 of 20 nt, 9..14 of 21 nt.

    private static string GuideWithGcCount(int length, int gcCount)
    {
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = i < gcCount ? (i % 2 == 0 ? 'G' : 'C') : (i % 2 == 0 ? 'A' : 'T');
        return new string(chars);
    }

    private static bool DesignedGuideListed(
        string guide, string pam, CrisprSystemType system, GuideRnaParameters parameters)
    {
        string text = "ACTGACTCAT" + guide + pam + "CATCACTGAT";
        return CrisprDesigner
            .DesignGuideRnas(new DnaSequence(text), 0, text.Length - 1, system, parameters)
            .Any(c => c.IsForwardStrand && c.Sequence == guide);
    }

    [Test]
    public void DesignGuideRnas_GcWindow_FiltersWhetherOrNotTheScoreTermIsOn_Spcas9([Range(0, 20)] int gcCount)
    {
        string guide = GuideWithGcCount(20, gcCount);
        bool expected = gcCount is >= 8 and <= 14;
        var noScoreTerm = GuideRnaParameters.Default with
        { ScoreGcContent = false, MinScore = 0, CheckSelfComplementarity = false };
        var scoreTerm = GuideRnaParameters.Default with { MinScore = 0, CheckSelfComplementarity = false };

        Assert.Multiple(() =>
        {
            Assert.That(DesignedGuideListed(guide, "AGG", CrisprSystemType.SpCas9, noScoreTerm), Is.EqualTo(expected),
                "ScoreGcContent = false must not switch the window off (a 100 %-GC guide used to pass)");
            Assert.That(DesignedGuideListed(guide, "AGG", CrisprSystemType.SpCas9, scoreTerm), Is.EqualTo(expected));
        });
    }

    [Test]
    public void DesignGuideRnas_GcWindow_FiltersFor21NtSaCas9([Range(0, 21)] int gcCount)
    {
        string guide = GuideWithGcCount(21, gcCount);
        var parameters = GuideRnaParameters.Default with
        { ScoreGcContent = false, MinScore = 0, CheckSelfComplementarity = false };

        Assert.That(DesignedGuideListed(guide, "CCGAAT", CrisprSystemType.SaCas9, parameters),
            Is.EqualTo(gcCount is >= 9 and <= 14), $"{gcCount}/21 = {100.0 * gcCount / 21:F2} %");
    }

    /// <summary>CHOPCHOP's filter defaults to 0..100 (no filter): the same window values switch the filter off.</summary>
    [Test]
    public void DesignGuideRnas_FullGcWindow_ListsEveryGuide_AndEvaluateOnlyReports()
    {
        string allGc = GuideWithGcCount(20, 20);
        var open = GuideRnaParameters.Default with
        { MinGcContent = 0, MaxGcContent = 100, ScoreGcContent = false, MinScore = 0, CheckSelfComplementarity = false };

        Assert.Multiple(() =>
        {
            Assert.That(DesignedGuideListed(allGc, "AGG", CrisprSystemType.SpCas9, open), Is.True);
            // EvaluateGuideRna evaluates the guide it is given: it reports the GC and does not drop it.
            var evaluated = CrisprDesigner.EvaluateGuideRna(allGc, CrisprSystemType.SpCas9,
                GuideRnaParameters.Default with { ScoreGcContent = false });
            Assert.That(evaluated.GcContent, Is.EqualTo(100.0));
            Assert.That(evaluated.Issues, Has.Some.Contains("High GC content"));
        });
    }

    #endregion
}
