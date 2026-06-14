// ONCO-SOMATIC-001 — Somatic Mutation Calling
// Evidence: docs/Evidence/ONCO-SOMATIC-001-Evidence.md
// TestSpec: tests/TestSpecs/ONCO-SOMATIC-001.md
// Source: Saunders CT et al. (2012). Bioinformatics 28(14):1811–1817. https://doi.org/10.1093/bioinformatics/bts271
//         Yan YH et al. (2021). Sci. Rep. 11:11640. https://doi.org/10.1038/s41598-021-91142-1

using VO = Seqeron.Genomics.Oncology.OncologyAnalyzer.VariantObservation;
using Status = Seqeron.Genomics.Oncology.OncologyAnalyzer.SomaticStatus;

namespace Seqeron.Genomics.Tests;

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
}
