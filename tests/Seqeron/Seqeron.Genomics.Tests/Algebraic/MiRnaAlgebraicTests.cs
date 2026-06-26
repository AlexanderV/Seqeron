using System.Linq;

namespace Seqeron.Genomics.Tests.Algebraic;

/// <summary>
/// Algebraic-law tests for the MiRNA area — TargetScan/context++ site scoring, the
/// Friedman/TargetScan PCT sigmoid, the pre-miRNA logistic classifier, and the
/// Drosha/Dicer cleavage rulers.
///
/// Algebraic testing pins the additive feature-coefficient dot product of the
/// context++ score, the published-logistic closed forms of PCT and the pre-miRNA
/// classifier, the integer measuring rules of the cleavage model, and the
/// determinism of all four.
/// — docs/checklists/06_ALGEBRAIC_TESTING.md §Description, rows 252, 253, 254, 255.
/// </summary>
[TestFixture]
[Category("Algebraic")]
[Category("MiRNA")]
public class MiRnaAlgebraicTests
{
    // ═══════════════════════════════════════════════════════════════════════
    // Unit: MIRNA-CONTEXT-001 — context++ target-site score (MiRNA), row 252.
    // ID/DIST — the context++ score is a feature-coefficient DOT PRODUCT: the total equals the
    //           intercept plus the sum of every per-feature contribution (Agarwal et al. 2015).
    // IDEMP — the score is a pure, deterministic function.
    //   — MiRnaAnalyzer.ScoreTargetSiteContextPlusPlus; TestSpec MIRNA-CONTEXT-001.
    // ═══════════════════════════════════════════════════════════════════════

    private static readonly MiRnaAnalyzer.MiRna Let7a =
        MiRnaAnalyzer.CreateMiRna("let-7a", "UGAGGUAGUAGGUUGUAUAGUU");
    private const string Let7aSeedRC = "CUACCUC";

    private static MiRnaAnalyzer.ContextPlusPlusScore ScoreLet7a8mer()
    {
        string mrna = "GGGGG" + Let7aSeedRC + "A" + "GGGGG"; // canonical 8mer fixture
        var site = MiRnaAnalyzer.FindTargetSites(mrna, Let7a, minScore: 0.0).Single();
        return MiRnaAnalyzer.ScoreTargetSiteContextPlusPlus(mrna, Let7a, site);
    }

    [Test]
    public void Context_Distributive_TotalIsInterceptPlusSumOfContributions()
    {
        var ctx = ScoreLet7a8mer();
        double sum = ctx.Intercept
                     + ctx.LocalAuContribution + ctx.SRna1Contribution + ctx.SRna8Contribution
                     + ctx.Site8Contribution + ctx.SaContribution + ctx.ThreePrimePairingContribution
                     + ctx.MinDistContribution + ctx.Len3UtrContribution + ctx.Off6mContribution
                     + ctx.SpsContribution + ctx.TaContribution + ctx.LenOrfContribution
                     + ctx.Orf8mContribution + ctx.PctContribution;
        ctx.ContextScorePartial.Should().BeApproximately(sum, 1e-12,
            "the context++ score is the dot product of feature values and coefficients (sum of contributions)");
    }

    [Test]
    public void Context_Identity_ReferenceWindowReproducesPerlScore()
    {
        // Hand-derived 8mer partial context++ score reproduced verbatim from the perl reference.
        ScoreLet7a8mer().ContextScorePartial.Should().BeApproximately(-0.7561913315126536, 1e-9);
    }

    [Test]
    public void Context_Idempotent_Deterministic()
    {
        ScoreLet7a8mer().ContextScorePartial.Should().Be(ScoreLet7a8mer().ContextScorePartial);
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: MIRNA-PCT-001 — probability of conserved targeting (MiRNA), row 253.
    // ID — PCT(Bls) = B0 + B1/(1 + e^(−B2·Bls + B3)): the published TargetScan logistic
    //      (targetscan_70_BL_PCT.pl). For (B0,B1,B2,B3)=(0,1,1,0), PCT(3) = 1/(1+e^−3).
    // IDEMP — the mapping is a pure, deterministic function.
    //   — MiRnaAnalyzer.PctFromBranchLength; Friedman et al. (2009); TargetScan.
    // ═══════════════════════════════════════════════════════════════════════

    [Test]
    public void Pct_Identity_BranchLengthMapsToPublishedLogistic()
    {
        var p = new MiRnaAnalyzer.PctSigmoidParameters(B0: 0.0, B1: 1.0, B2: 1.0, B3: 0.0);
        double pct = MiRnaAnalyzer.PctFromBranchLength(3.0, p);
        // 1/(1+e^-3) = 0.9525741268224334.
        pct.Should().BeApproximately(1.0 / (1.0 + System.Math.Exp(-3.0)), 1e-12);
        pct.Should().BeApproximately(0.952574126822433, 1e-12);
    }

    [Test]
    public void Pct_Identity_NegativeValueTruncatedToZero()
    {
        // targetscan_70_BL_PCT.pl: if ($pct < 0) { $pct = 0 }.
        var p = new MiRnaAnalyzer.PctSigmoidParameters(B0: -0.5, B1: 0.3, B2: 1.0, B3: 5.0);
        MiRnaAnalyzer.PctFromBranchLength(0.0, p).Should().Be(0.0);
    }

    [Test]
    public void Pct_Idempotent_Deterministic()
    {
        var p = new MiRnaAnalyzer.PctSigmoidParameters(0.0, 1.0, 1.0, 0.0);
        MiRnaAnalyzer.PctFromBranchLength(2.5, p).Should().Be(MiRnaAnalyzer.PctFromBranchLength(2.5, p));
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: MIRNA-CLASSIFY-001 — pre-miRNA logistic classifier (MiRNA), row 254.
    // ID — the natural probability is the LOGISTIC of the (standardised) feature vector, so it lies
    //      strictly in (0,1) and IsNatural ⇔ probability ≥ 0.5. A real precursor scores near 1.
    // IDEMP — the classifier is a pure, deterministic function.
    //   — MiRnaAnalyzer.ClassifyPreMiRna; logistic regression over MFE/AMFE/MFEI/GC/%paired.
    // ═══════════════════════════════════════════════════════════════════════

    private const string HsaMir21 =
        "UGUCGGGUAGCUUAUCAGACUGAUGUUGACUGUUGAAUCUCAUGGCAACACCAGUCGAUGGGCUGUCUGACA";

    [Test]
    public void Classify_Identity_ProbabilityIsLogisticAndDrivesCall()
    {
        var c = MiRnaAnalyzer.ClassifyPreMiRna(HsaMir21);
        c.Should().NotBeNull();
        // Logistic output is a probability in (0,1)…
        c!.Value.NaturalProbability.Should().BeInRange(0.0, 1.0);
        c.Value.NaturalProbability.Should().BeGreaterThan(0.0).And.BeLessThan(1.0);
        // …and the boolean call is exactly the 0.5 threshold of that logistic.
        c.Value.IsNatural.Should().Be(c.Value.NaturalProbability >= 0.5);
        // hsa-mir-21 is a genuine precursor ⇒ scores near 1.
        c.Value.IsNatural.Should().BeTrue();
        c.Value.NaturalProbability.Should().BeGreaterThan(0.95);
    }

    [Test]
    public void Classify_Idempotent_Deterministic()
    {
        var a = MiRnaAnalyzer.ClassifyPreMiRna(HsaMir21);
        var b = MiRnaAnalyzer.ClassifyPreMiRna(HsaMir21);
        a!.Value.NaturalProbability.Should().Be(b!.Value.NaturalProbability);
        a.Value.IsNatural.Should().Be(b.Value.IsNatural);
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: MIRNA-CLEAVAGE-001 — Drosha/Dicer cleavage rulers (MiRNA), row 255.
    // ID — the integer measuring rules place the mature 5p product: Drosha cuts ~11 bp from the
    //      basal junction (Han 2006) and Dicer's 5' counting rule fixes the mature at ~22 nt
    //      (Park 2011). With basalJunction 0 the mature is exactly pri[11 .. 11+22).
    // IDEMP — the prediction is a pure, deterministic function.
    //   — MiRnaAnalyzer.PredictDroshaDicerCleavage; TestSpec MIRNA-CLEAVAGE-001.
    // ═══════════════════════════════════════════════════════════════════════

    // 11-nt lower stem + the miR-21 stem region (≥ 11 + 22 + 2 nt).
    private const string CleavagePri =
        "CCCCCCCCCCC" + "UAGCUUAUCAGACUGAUGUUGACUGUUGAAUCUCAUGGCAACACCAGUCGAUGGGCUGU";

    [Test]
    public void Cleavage_Identity_ElevenBpAndTwentyTwoNtRulesYieldMature()
    {
        var cut = MiRnaAnalyzer.PredictDroshaDicerCleavage(CleavagePri, basalJunction: 0);
        cut.Should().NotBeNull();

        // 11-bp ruler: Drosha 5' cut at basalJunction + 11 = 11; the mature starts there.
        cut!.Value.DroshaCut5Prime.Should().Be(11);
        cut.Value.MatureStart.Should().Be(11);
        // 22-nt ruler: the mature length and inclusive span are exactly 22 nt.
        cut.Value.MatureSequence.Length.Should().Be(22);
        (cut.Value.MatureEnd - cut.Value.MatureStart + 1).Should().Be(22);
        // The mature sequence is exactly the measured window of the precursor.
        cut.Value.MatureSequence.Should().Be(CleavagePri.Substring(11, 22));
    }

    [Test]
    public void Cleavage_Idempotent_Deterministic()
    {
        var a = MiRnaAnalyzer.PredictDroshaDicerCleavage(CleavagePri, 0);
        var b = MiRnaAnalyzer.PredictDroshaDicerCleavage(CleavagePri, 0);
        a!.Value.MatureSequence.Should().Be(b!.Value.MatureSequence);
        a.Value.DroshaCut5Prime.Should().Be(b.Value.DroshaCut5Prime);
    }
}
