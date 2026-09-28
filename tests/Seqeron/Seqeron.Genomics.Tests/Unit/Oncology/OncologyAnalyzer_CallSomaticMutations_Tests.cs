// ONCO-SOMATIC-001 — Somatic Mutation Calling
// Evidence: docs/Evidence/ONCO-SOMATIC-001-Evidence.md
// TestSpec: tests/TestSpecs/ONCO-SOMATIC-001.md
// Source: Saunders CT et al. (2012). Bioinformatics 28(14):1811–1817. https://doi.org/10.1093/bioinformatics/bts271
//         Yan YH et al. (2021). Sci. Rep. 11:11640. https://doi.org/10.1038/s41598-021-91142-1

using VO = Seqeron.Genomics.Oncology.OncologyAnalyzer.VariantObservation;
using Status = Seqeron.Genomics.Oncology.OncologyAnalyzer.SomaticStatus;

namespace Seqeron.Genomics.Tests.Unit.Oncology;

[TestFixture]
public class OncologyAnalyzer_CallSomaticMutations_Tests
{
    private static VO Make(int tAlt, int tTot, int nAlt, int nTot)
        => new("chr1", 100, "A", "T", tAlt, tTot, nAlt, nTot);

    #region CallSomaticMutations / Classify Tests

    // M1 — Clear somatic: f_t=0.25 ≥ 0.05, f_n=0.00 ≤ 0.01 → Somatic; score = 0.25 − 0.00.
    // Evidence: Saunders 2012 — S={f_t≠f_n}, ref/ref normal.
    [Test]
    public void CallSomaticMutations_HighTumorVafZeroNormal_ClassifiesSomatic()
    {
        var calls = OncologyAnalyzer.CallSomaticMutations(new[] { Make(25, 100, 0, 100) });

        Assert.Multiple(() =>
        {
            Assert.That(calls[0].Status, Is.EqualTo(Status.Somatic), "Present in tumor, absent in normal is somatic");
            Assert.That(calls[0].TumorVaf, Is.EqualTo(0.25).Within(1e-10), "f_t = 25/100");
            Assert.That(calls[0].NormalVaf, Is.EqualTo(0.00).Within(1e-10), "f_n = 0/100");
            Assert.That(calls[0].SomaticScore, Is.EqualTo(0.25).Within(1e-10), "score = f_t - f_n = 0.25");
        });
    }

    // M2 — Germline het: f_t=0.48, f_n=0.50 > 0.01 → Germline; score 0.
    // Evidence: Mutect2 (Benjamin 2019) skips variants present in matched normal.
    [Test]
    public void CallSomaticMutations_PresentInTumorAndNormal_ClassifiesGermline()
    {
        var calls = OncologyAnalyzer.CallSomaticMutations(new[] { Make(48, 100, 50, 100) });

        Assert.Multiple(() =>
        {
            Assert.That(calls[0].Status, Is.EqualTo(Status.Germline), "Present in both tumor and normal is germline");
            Assert.That(calls[0].SomaticScore, Is.EqualTo(0.0).Within(1e-10), "Germline score is 0");
        });
    }

    // M3 — Sub-LoD tumor: f_t=0.02 < 0.05 → NotDetected.
    // Evidence: Yan 2021 — WES VAF LoD = 5%; ≤5% frequently errors.
    [Test]
    public void CallSomaticMutations_TumorVafBelowLimitOfDetection_ClassifiesNotDetected()
    {
        var calls = OncologyAnalyzer.CallSomaticMutations(new[] { Make(2, 100, 0, 100) });

        Assert.Multiple(() =>
        {
            Assert.That(calls[0].Status, Is.EqualTo(Status.NotDetected), "Tumor VAF below 5% LoD is not detected");
            Assert.That(calls[0].SomaticScore, Is.EqualTo(0.0).Within(1e-10), "NotDetected score is 0");
        });
    }

    // M4 — Tumor-only mode: normal total 0 ⇒ f_n=0 → Somatic; score 0.20.
    // Evidence: Mutect2 (Benjamin 2019) — "If we have no matched normal, ℓ_n = 1".
    [Test]
    public void CallSomaticMutations_TumorOnlyNoNormalCoverage_ClassifiesSomatic()
    {
        var calls = OncologyAnalyzer.CallSomaticMutations(new[] { Make(20, 100, 0, 0) });

        Assert.Multiple(() =>
        {
            Assert.That(calls[0].NormalVaf, Is.EqualTo(0.0).Within(1e-10), "Uncovered normal yields VAF 0");
            Assert.That(calls[0].Status, Is.EqualTo(Status.Somatic), "Tumor-only present variant is somatic");
            Assert.That(calls[0].SomaticScore, Is.EqualTo(0.20).Within(1e-10), "score = 0.20 - 0.00");
        });
    }

    // M5 — Tumor threshold boundary: f_t exactly 0.05 is inclusive (present).
    // Evidence: DefaultTumorVafThreshold = 0.05.
    [Test]
    public void CallSomaticMutations_TumorVafAtThreshold_IsPresentAndSomatic()
    {
        var calls = OncologyAnalyzer.CallSomaticMutations(new[] { Make(5, 100, 0, 100) });

        Assert.That(calls[0].Status, Is.EqualTo(Status.Somatic),
            "f_t == 0.05 meets the inclusive presence threshold");
    }

    // M6 — Normal threshold boundary: f_n exactly 0.01 is inclusive (absent) → Somatic; score 0.30-0.01.
    // Evidence: DefaultNormalVafThreshold = 0.01.
    [Test]
    public void CallSomaticMutations_NormalVafAtThreshold_IsAbsentAndSomatic()
    {
        var calls = OncologyAnalyzer.CallSomaticMutations(new[] { Make(30, 100, 1, 100) });

        Assert.Multiple(() =>
        {
            Assert.That(calls[0].Status, Is.EqualTo(Status.Somatic), "f_n == 0.01 meets the inclusive absence ceiling");
            Assert.That(calls[0].SomaticScore, Is.EqualTo(0.29).Within(1e-10), "score = 0.30 - 0.01");
        });
    }

    // M7 — CHIP-like normal: f_n=0.03 > 0.01 → Germline.
    // Evidence: Mutect2 (Benjamin 2019) — clonal-hematopoiesis contamination present in normal.
    [Test]
    public void CallSomaticMutations_LowLevelNormalContamination_ClassifiesGermline()
    {
        var calls = OncologyAnalyzer.CallSomaticMutations(new[] { Make(30, 100, 3, 100) });

        Assert.That(calls[0].Status, Is.EqualTo(Status.Germline),
            "Normal VAF above 1% ceiling is treated as present in normal (germline/contamination)");
    }

    // M10 — Order and count preserved across a mixed panel.
    // Evidence: INV-6 implementation contract.
    [Test]
    public void CallSomaticMutations_MixedPanel_PreservesOrderAndCount()
    {
        var input = new[]
        {
            Make(25, 100, 0, 100),  // Somatic
            Make(48, 100, 50, 100), // Germline
            Make(2, 100, 0, 100),   // NotDetected
        };

        var calls = OncologyAnalyzer.CallSomaticMutations(input);

        Assert.Multiple(() =>
        {
            Assert.That(calls, Has.Count.EqualTo(3), "One call per input variant");
            Assert.That(calls[0].Status, Is.EqualTo(Status.Somatic), "Order preserved: first is somatic");
            Assert.That(calls[1].Status, Is.EqualTo(Status.Germline), "Order preserved: second is germline");
            Assert.That(calls[2].Status, Is.EqualTo(Status.NotDetected), "Order preserved: third not detected");
        });
    }

    // S1 — Normal just above threshold: f_n=0.02 → Germline.
    [Test]
    public void CallSomaticMutations_NormalVafJustAboveThreshold_ClassifiesGermline()
    {
        var calls = OncologyAnalyzer.CallSomaticMutations(new[] { Make(30, 100, 2, 100) });

        Assert.That(calls[0].Status, Is.EqualTo(Status.Germline), "f_n = 0.02 exceeds the 0.01 absence ceiling");
    }

    // S2 — Custom thresholds reclassify a sub-default variant as somatic.
    [Test]
    public void CallSomaticMutations_LowerTumorThreshold_ReclassifiesAsSomatic()
    {
        var variant = new[] { Make(3, 100, 0, 100) }; // f_t = 0.03

        var defaultCalls = OncologyAnalyzer.CallSomaticMutations(variant);
        var customCalls = OncologyAnalyzer.CallSomaticMutations(variant, tumorVafThreshold: 0.02);

        Assert.Multiple(() =>
        {
            Assert.That(defaultCalls[0].Status, Is.EqualTo(Status.NotDetected), "0.03 < default 0.05");
            Assert.That(customCalls[0].Status, Is.EqualTo(Status.Somatic), "0.03 ≥ custom 0.02 → somatic");
        });
    }

    // C1 — Empty input yields an empty result.
    [Test]
    public void CallSomaticMutations_EmptyInput_ReturnsEmpty()
    {
        var calls = OncologyAnalyzer.CallSomaticMutations(System.Array.Empty<VO>());

        Assert.That(calls, Is.Empty, "No variants in, no calls out");
    }

    [Test]
    public void CallSomaticMutations_NullInput_Throws()
    {
        Assert.Throws<System.ArgumentNullException>(
            () => OncologyAnalyzer.CallSomaticMutations(null!),
            "Null variant collection is rejected");
    }

    [Test]
    public void CallSomaticMutations_ThresholdOutOfRange_Throws()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<System.ArgumentOutOfRangeException>(
                () => OncologyAnalyzer.CallSomaticMutations(System.Array.Empty<VO>(), tumorVafThreshold: 1.5),
                "Tumor threshold above 1 is rejected");
            Assert.Throws<System.ArgumentOutOfRangeException>(
                () => OncologyAnalyzer.CallSomaticMutations(System.Array.Empty<VO>(), normalVafThreshold: -0.1),
                "Negative normal threshold is rejected");
        });
    }

    [Test]
    public void Classify_AltExceedsTotal_Throws()
    {
        Assert.Throws<System.ArgumentOutOfRangeException>(
            () => OncologyAnalyzer.Classify(Make(120, 100, 0, 100)),
            "Alt reads cannot exceed total reads");
    }

    // Read counts must be non-negative (contract §3.3); the negative-read guard is part of
    // this unit's classification path (CalculateVaf), not only the standalone CalculateVAF.
    [Test]
    public void Classify_NegativeReadCounts_Throw()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<System.ArgumentOutOfRangeException>(
                () => OncologyAnalyzer.Classify(Make(-1, 100, 0, 100)),
                "Negative tumor alt reads are rejected");
            Assert.Throws<System.ArgumentOutOfRangeException>(
                () => OncologyAnalyzer.Classify(Make(10, 100, -2, 100)),
                "Negative normal alt reads are rejected");
            Assert.Throws<System.ArgumentOutOfRangeException>(
                () => OncologyAnalyzer.Classify(Make(10, -100, 0, 100)),
                "Negative tumor total reads are rejected");
        });
    }

    // Direct coverage of the public single-variant Classify overload: all three Stage-A branches.
    // f_t=0.25,f_n=0 → Somatic (Saunders 2012); f_t=0.48,f_n=0.50 → Germline (Benjamin 2019);
    // f_t=0.02<0.05 → NotDetected (Yan 2021). Threshold defaults τ_t=0.05, τ_n=0.01.
    [Test]
    public void Classify_SingleVariant_CoversAllThreeBranches()
    {
        Assert.Multiple(() =>
        {
            Assert.That(OncologyAnalyzer.Classify(Make(25, 100, 0, 100)).Status,
                Is.EqualTo(Status.Somatic), "f_t=0.25 ≥ 0.05, f_n=0.00 ≤ 0.01 → Somatic");
            Assert.That(OncologyAnalyzer.Classify(Make(48, 100, 50, 100)).Status,
                Is.EqualTo(Status.Germline), "f_t=0.48, f_n=0.50 > 0.01 → Germline");
            Assert.That(OncologyAnalyzer.Classify(Make(2, 100, 0, 100)).Status,
                Is.EqualTo(Status.NotDetected), "f_t=0.02 < 0.05 → NotDetected");
        });
    }

    #endregion

    #region FilterGermlineVariants Tests

    // M8 — Filter returns exactly the somatic subset, in order.
    // Evidence: Mutect2 (Benjamin 2019) — calls somatic variants only.
    [Test]
    public void FilterGermlineVariants_MixedPanel_ReturnsOnlySomatic()
    {
        var input = new[]
        {
            Make(25, 100, 0, 100),  // Somatic
            Make(48, 100, 50, 100), // Germline (removed)
            Make(2, 100, 0, 100),   // NotDetected (removed)
            Make(40, 100, 0, 100),  // Somatic
        };

        var somatic = OncologyAnalyzer.FilterGermlineVariants(input);

        Assert.Multiple(() =>
        {
            Assert.That(somatic, Has.Count.EqualTo(2), "Only the two somatic variants remain");
            Assert.That(somatic.All(c => c.Status == Status.Somatic), Is.True, "Every retained call is somatic");
            Assert.That(somatic[0].TumorVaf, Is.EqualTo(0.25).Within(1e-10), "First somatic kept in order");
            Assert.That(somatic[1].TumorVaf, Is.EqualTo(0.40).Within(1e-10), "Second somatic kept in order");
        });
    }

    [Test]
    public void FilterGermlineVariants_NullInput_Throws()
    {
        Assert.Throws<System.ArgumentNullException>(
            () => OncologyAnalyzer.FilterGermlineVariants(null!),
            "Null variant collection is rejected");
    }

    #endregion

    #region CalculateSomaticScore Tests

    // M9 — Score = f_t − f_n = 0.25 − 0.05 = 0.20.
    // Evidence: INV-3 separation score (ASSUMPTION, documented).
    [Test]
    public void CalculateSomaticScore_SeparationBetweenTumorAndNormal_ReturnsDifference()
    {
        double score = OncologyAnalyzer.CalculateSomaticScore(Make(25, 100, 5, 100));

        Assert.That(score, Is.EqualTo(0.20).Within(1e-10), "score = 0.25 - 0.05 = 0.20");
    }

    // S3 — Score is 0 when the normal carries the allele at or above the tumor level.
    [Test]
    public void CalculateSomaticScore_NormalAtOrAboveTumor_ReturnsZero()
    {
        double score = OncologyAnalyzer.CalculateSomaticScore(Make(10, 100, 20, 100));

        Assert.That(score, Is.EqualTo(0.0).Within(1e-10), "No somatic separation when f_n ≥ f_t");
    }

    #endregion

    #region Mutect2 somatic likelihoods model (CallSomaticMutationsMutect2)

    // Reference values: a per-read Python port of GATK SomaticLikelihoodsEngine.logEvidence,
    // SomaticGenotypingEngine.somaticLogOdds / diploidAltLogOdds, PairHMM tri-state likelihoods (ε/3) and the
    // Q45 global-mismapping cap (broadinstitute/gatk master, fetched 2026-09-28). The variational TLOD was
    // also checked against the exact flat-prior marginal likelihood (numerical integration): e.g.
    // 25 alt / 75 ref at Q30: variational 61.5425 vs exact 61.5429 log10 units.
    private const double LodTolerance = 1e-8;

    [TestCase(75, 25, 30, 61.54246801643632)]
    [TestCase(52, 48, 30, 135.9216837925928)]
    [TestCase(98, 2, 30, 1.2661263748027127)]
    [TestCase(97, 3, 30, 3.2295883512710644)]
    [TestCase(97, 3, 20, 0.3263769386302512)]
    [TestCase(95, 5, 30, 7.515738706497197)]
    [TestCase(80, 20, 40, 66.80951519855469)]
    [TestCase(80, 20, 50, 67.26788516897696)] // Q50: mismatch likelihood floored by the Q45 mismapping cap
    [TestCase(0, 10, 30, 33.72555096640962)]
    [TestCase(100, 0, 30, -2.0042635241874462)]
    public void CalculateMutect2TumorLog10Odds_MatchesGatkReference(int refReads, int altReads, int q, double expected)
    {
        double tlod = OncologyAnalyzer.CalculateMutect2TumorLog10Odds(refReads, altReads, q);

        Assert.That(tlod, Is.EqualTo(expected).Within(LodTolerance),
            $"TLOD({refReads} ref, {altReads} alt, Q{q}) must equal the GATK somatic-likelihoods reference");
    }

    [TestCase(100, 0, 30, 30.088511009736504)]
    [TestCase(50, 50, 30, -143.74582613754575)]
    [TestCase(97, 3, 30, 19.658450780899575)]
    [TestCase(99, 1, 30, 26.61182426679086)]
    [TestCase(7, 0, 30, 2.106195770681557)]
    [TestCase(8, 0, 30, 2.407080880778922)]
    [TestCase(80, 20, 50, -59.89837377162566)]
    public void CalculateMutect2NormalLog10Odds_MatchesGatkReference(int refReads, int altReads, int q, double expected)
    {
        double nlod = OncologyAnalyzer.CalculateMutect2NormalLog10Odds(refReads, altReads, q);

        Assert.That(nlod, Is.EqualTo(expected).Within(LodTolerance),
            $"NLOD({refReads} ref, {altReads} alt, Q{q}) = log10 P(hom-ref)/P(het)");
    }

    [Test]
    public void CalculateMutect2LogOdds_EmptyPileup_IsZero()
    {
        Assert.Multiple(() =>
        {
            Assert.That(OncologyAnalyzer.CalculateMutect2TumorLog10Odds(0, 0, 30), Is.EqualTo(0.0),
                "GATK: evidence 0 for both allele sets when there are no reads");
            Assert.That(OncologyAnalyzer.CalculateMutect2NormalLog10Odds(0, 0, 30), Is.EqualTo(0.0),
                "No normal reads: hom-ref and het equally likely");
        });
    }

    [Test]
    public void CallSomaticMutationsMutect2_ClassifiesByTlodAndNlod()
    {
        var input = new[]
        {
            Make(25, 100, 0, 100),  // TLOD 61.54 > 3, NLOD 30.09 > 2.2 → Somatic
            Make(48, 100, 50, 100), // TLOD 135.92 > 3, NLOD −143.75 ≤ 2.2 → Germline
            Make(2, 100, 0, 100),   // TLOD 1.27 ≤ 3 → NotDetected
            Make(30, 100, 3, 100),  // normal 3/100: NLOD 19.66 > 2.2 → Somatic (not a het germline)
        };

        var calls = OncologyAnalyzer.CallSomaticMutationsMutect2(input, baseQuality: 30);

        Assert.Multiple(() =>
        {
            Assert.That(calls, Has.Count.EqualTo(4), "One call per variant, in order");
            Assert.That(calls[0].Status, Is.EqualTo(Status.Somatic));
            Assert.That(calls[0].TumorLog10Odds, Is.EqualTo(61.54246801643632).Within(LodTolerance));
            Assert.That(calls[0].NormalLog10Odds, Is.EqualTo(30.088511009736504).Within(LodTolerance));
            Assert.That(calls[0].TumorVaf, Is.EqualTo(0.25).Within(1e-12));
            Assert.That(calls[1].Status, Is.EqualTo(Status.Germline));
            Assert.That(calls[2].Status, Is.EqualTo(Status.NotDetected));
            Assert.That(calls[3].Status, Is.EqualTo(Status.Somatic));
            Assert.That(calls[3].NormalLog10Odds, Is.EqualTo(19.658450780899575).Within(LodTolerance));
        });
    }

    // Base quality matters: 3/100 alt reads pass the TLOD 3.0 emission threshold at Q30 (3.23) but not at Q20 (0.33).
    [Test]
    public void CallSomaticMutationsMutect2_LowBaseQuality_RaisesDetectionBar()
    {
        var variant = new[] { Make(3, 100, 0, 100) };

        Assert.Multiple(() =>
        {
            Assert.That(OncologyAnalyzer.CallSomaticMutationsMutect2(variant, 30)[0].Status, Is.EqualTo(Status.Somatic));
            Assert.That(OncologyAnalyzer.CallSomaticMutationsMutect2(variant, 20)[0].Status, Is.EqualTo(Status.NotDetected));
        });
    }

    // NLOD per hom-ref Q30 read = log10(2(1−ε)/(1−ε+ε/3)) ≈ 0.3009 ⇒ 7 reads 2.106 ≤ 2.2 (germline not excluded),
    // 8 reads 2.407 > 2.2 (somatic): Mutect2 needs normal depth to rule out a germline het.
    [Test]
    public void CallSomaticMutationsMutect2_NormalDepthRequiredToExcludeGermline()
    {
        Assert.Multiple(() =>
        {
            Assert.That(OncologyAnalyzer.CallSomaticMutationsMutect2(new[] { Make(25, 100, 0, 7) }, 30)[0].Status,
                Is.EqualTo(Status.Germline), "7 normal ref reads: NLOD 2.106 ≤ 2.2");
            Assert.That(OncologyAnalyzer.CallSomaticMutationsMutect2(new[] { Make(25, 100, 0, 8) }, 30)[0].Status,
                Is.EqualTo(Status.Somatic), "8 normal ref reads: NLOD 2.407 > 2.2");
        });
    }

    // mutect.tex: "If we have no matched normal, ℓ_n = 1" — tumor-only mode skips the NLOD test. With a matched
    // normal of zero coverage NLOD = 0 ≤ 2.2 and Mutect2 skips the allele as possibly germline.
    [Test]
    public void CallSomaticMutationsMutect2_TumorOnlyVersusUncoveredNormal()
    {
        var variant = new[] { Make(20, 100, 0, 0) };

        var tumorOnly = OncologyAnalyzer.CallSomaticMutationsMutect2(variant, 30, hasMatchedNormal: false)[0];
        var uncovered = OncologyAnalyzer.CallSomaticMutationsMutect2(variant, 30, hasMatchedNormal: true)[0];

        Assert.Multiple(() =>
        {
            Assert.That(tumorOnly.Status, Is.EqualTo(Status.Somatic));
            Assert.That(double.IsNaN(tumorOnly.NormalLog10Odds), Is.True, "No NLOD without a matched normal");
            Assert.That(uncovered.NormalLog10Odds, Is.EqualTo(0.0));
            Assert.That(uncovered.Status, Is.EqualTo(Status.Germline));
        });
    }

    [Test]
    public void CallSomaticMutationsMutect2_InvalidInput_Throws()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<System.ArgumentNullException>(() => OncologyAnalyzer.CallSomaticMutationsMutect2(null!, 30));
            Assert.Throws<System.ArgumentOutOfRangeException>(
                () => OncologyAnalyzer.CallSomaticMutationsMutect2(System.Array.Empty<VO>(), 0), "Q must be ≥ 1");
            Assert.Throws<System.ArgumentOutOfRangeException>(
                () => OncologyAnalyzer.CallSomaticMutationsMutect2(System.Array.Empty<VO>(), 30, tumorLog10OddsThreshold: double.NaN));
            Assert.Throws<System.ArgumentOutOfRangeException>(
                () => OncologyAnalyzer.CallSomaticMutationsMutect2(new[] { Make(120, 100, 0, 100) }, 30), "alt > total");
            Assert.Throws<System.ArgumentOutOfRangeException>(
                () => OncologyAnalyzer.CalculateMutect2TumorLog10Odds(-1, 3, 30));
        });
    }

    #endregion
}
