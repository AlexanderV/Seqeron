// ONCO-TMB-001 — Tumor Mutational Burden
// Evidence: docs/Evidence/ONCO-TMB-001-Evidence.md
// TestSpec: tests/TestSpecs/ONCO-TMB-001.md
// Source: Chalmers ZR et al. (2017). Genome Medicine 9:34. https://doi.org/10.1186/s13073-017-0424-2
//         Marcus L et al. (2021). FDA Approval Summary: Pembrolizumab for TMB-High Solid Tumors.
//                 Clin Cancer Res 27(17):4685-4689. https://doi.org/10.1158/1078-0432.CCR-21-0327

namespace Seqeron.Genomics.Tests.Unit.Oncology;

[TestFixture]
public class OncologyAnalyzer_CalculateTMB_Tests
{
    #region CalculateTMB(int, double)

    // M1 — 315-gene FoundationOne panel: 11 somatic mutations / 1.1 Mb = 10.0 mut/Mb (Chalmers 2017).
    [Test]
    public void CalculateTMB_11over1_1Mb_Returns10()
    {
        double tmb = OncologyAnalyzer.CalculateTMB(11, 1.1);

        Assert.That(tmb, Is.EqualTo(10.0).Within(1e-10),
            "TMB = mutations/Mb = 11/1.1 = 10.0 mut/Mb (Chalmers 2017, 1.1 Mb panel)");
    }

    // M2 — whole-exome example: 300 mutations / 30 Mb = 10.0 mut/Mb (TMB = mut/Mb).
    [Test]
    public void CalculateTMB_300over30Mb_Returns10()
    {
        double tmb = OncologyAnalyzer.CalculateTMB(300, 30.0);

        Assert.That(tmb, Is.EqualTo(10.0).Within(1e-10),
            "TMB = 300/30 = 10.0 mut/Mb (definition: somatic mutations per megabase)");
    }

    // M3 — high count: 150 / 10 Mb = 15.0 mut/Mb.
    [Test]
    public void CalculateTMB_150over10Mb_Returns15()
    {
        double tmb = OncologyAnalyzer.CalculateTMB(150, 10.0);

        Assert.That(tmb, Is.EqualTo(15.0).Within(1e-10),
            "TMB = 150/10 = 15.0 mut/Mb");
    }

    // M4 — zero mutations over a non-zero region is TMB 0 (no somatic burden).
    [Test]
    public void CalculateTMB_ZeroMutations_ReturnsZero()
    {
        double tmb = OncologyAnalyzer.CalculateTMB(0, 10.0);

        Assert.That(tmb, Is.EqualTo(0.0).Within(1e-10),
            "0 mutations / 10 Mb = 0.0 mut/Mb");
    }

    // M5 — region = 0 Mb: TMB is undefined (division by zero) → throws.
    [Test]
    public void CalculateTMB_ZeroRegion_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.CalculateTMB(5, 0.0),
            "TMB is undefined when the megabase denominator is 0 (Registry edge case: division by zero)");
    }

    // S1 — small panel (<0.5 Mb): value is still mathematically defined; no exception (Chalmers 2017
    // documents instability, not an error). 2 / 0.3 = 6.6667 mut/Mb.
    [Test]
    public void CalculateTMB_SmallPanel_ComputesRatioWithoutThrowing()
    {
        double tmb = OncologyAnalyzer.CalculateTMB(2, 0.3);

        Assert.That(tmb, Is.EqualTo(2.0 / 0.3).Within(1e-10),
            "Sub-0.5-Mb panels still yield a defined ratio (2/0.3 ≈ 6.6667); instability is documented, not an error");
    }

    // S2 — negative mutation count is invalid.
    [Test]
    public void CalculateTMB_NegativeCount_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.CalculateTMB(-1, 1.0),
            "A negative mutation count is invalid");
    }

    // S3 — non-positive / non-finite region is invalid.
    [Test]
    public void CalculateTMB_InvalidRegion_Throws()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.CalculateTMB(5, -1.0),
                "A negative region size is invalid");
            Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.CalculateTMB(5, double.NaN),
                "A NaN region size is invalid");
            Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.CalculateTMB(5, double.PositiveInfinity),
                "An infinite region size is invalid");
        });
    }

    #endregion

    #region CalculateTMB(IEnumerable<SomaticCall>, double)

    // M10 — only Somatic-status calls are counted: 3 somatic + 1 germline + 1 not-detected / 1.0 Mb = 3.0.
    [Test]
    public void CalculateTMB_FromSomaticCalls_CountsOnlySomatic()
    {
        var calls = new[]
        {
            MakeCall(OncologyAnalyzer.SomaticStatus.Somatic),
            MakeCall(OncologyAnalyzer.SomaticStatus.Somatic),
            MakeCall(OncologyAnalyzer.SomaticStatus.Somatic),
            MakeCall(OncologyAnalyzer.SomaticStatus.Germline),
            MakeCall(OncologyAnalyzer.SomaticStatus.NotDetected),
        };

        double tmb = OncologyAnalyzer.CalculateTMB(calls, 1.0);

        Assert.That(tmb, Is.EqualTo(3.0).Within(1e-10),
            "Only the 3 Somatic calls count toward TMB (germline/not-detected excluded): 3/1.0 = 3.0 mut/Mb");
    }

    // Empty somatic-call collection: no somatic mutations → TMB = 0/Mb = 0.0 (no burden).
    [Test]
    public void CalculateTMB_FromEmptySomaticCalls_ReturnsZero()
    {
        double tmb = OncologyAnalyzer.CalculateTMB(
            System.Array.Empty<OncologyAnalyzer.SomaticCall>(), 1.1);

        Assert.That(tmb, Is.EqualTo(0.0).Within(1e-10),
            "No somatic calls → 0 mutations / 1.1 Mb = 0.0 mut/Mb");
    }

    // S4 — null collection is rejected.
    [Test]
    public void CalculateTMB_NullCalls_Throws()
    {
        Assert.Throws<ArgumentNullException>(
            () => OncologyAnalyzer.CalculateTMB((IEnumerable<OncologyAnalyzer.SomaticCall>)null!, 1.0),
            "A null somatic-call collection is invalid");
    }

    #endregion

    #region ClassifyTMB

    // M6 — below the FDA cutoff (9.9 < 10) is Low.
    [Test]
    public void ClassifyTMB_BelowCutoff_IsLow()
    {
        Assert.That(OncologyAnalyzer.ClassifyTMB(9.9), Is.EqualTo(OncologyAnalyzer.TmbStatus.Low),
            "9.9 mut/Mb is below the FDA TMB-High cutoff of 10 → Low");
    }

    // M7 — at the FDA cutoff (10.0) is High (inclusive boundary, Marcus 2021 "TMB ≥10 mut/Mb").
    [Test]
    public void ClassifyTMB_AtCutoff_IsHigh()
    {
        Assert.That(OncologyAnalyzer.ClassifyTMB(10.0), Is.EqualTo(OncologyAnalyzer.TmbStatus.High),
            "Exactly 10 mut/Mb is TMB-High (cutoff is inclusive: TMB ≥ 10, FDA pembrolizumab)");
    }

    // M8 — above the cutoff (15.0) is High.
    [Test]
    public void ClassifyTMB_AboveCutoff_IsHigh()
    {
        Assert.That(OncologyAnalyzer.ClassifyTMB(15.0), Is.EqualTo(OncologyAnalyzer.TmbStatus.High),
            "15 mut/Mb is above the FDA TMB-High cutoff → High");
    }

    // M9 — zero TMB is Low.
    [Test]
    public void ClassifyTMB_Zero_IsLow()
    {
        Assert.That(OncologyAnalyzer.ClassifyTMB(0.0), Is.EqualTo(OncologyAnalyzer.TmbStatus.Low),
            "0 mut/Mb is below the cutoff → Low");
    }

    // ClassifyTMB rejects invalid TMB values.
    [Test]
    public void ClassifyTMB_InvalidValue_Throws()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.ClassifyTMB(-0.1),
                "Negative TMB is invalid");
            Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.ClassifyTMB(double.NaN),
                "NaN TMB is invalid");
        });
    }

    #endregion

    #region Invariants (property-based)

    // C1 — INV-03: TMB is non-decreasing in count (fixed region) and non-increasing in region (fixed count).
    [Test]
    public void CalculateTMB_Monotonicity_HoldsOverSweep()
    {
        Assert.Multiple(() =>
        {
            // Non-decreasing in mutation count for a fixed region.
            double region = 5.0;
            double previous = double.NegativeInfinity;
            for (int count = 0; count <= 100; count++)
            {
                double tmb = OncologyAnalyzer.CalculateTMB(count, region);
                Assert.That(tmb, Is.GreaterThanOrEqualTo(previous),
                    $"TMB must not decrease as count rises (count={count}, region={region})");
                previous = tmb;
            }

            // Non-increasing in region size for a fixed count.
            int fixedCount = 100;
            previous = double.PositiveInfinity;
            for (double r = 1.0; r <= 50.0; r += 1.0)
            {
                double tmb = OncologyAnalyzer.CalculateTMB(fixedCount, r);
                Assert.That(tmb, Is.LessThanOrEqualTo(previous),
                    $"TMB must not increase as the region grows (count={fixedCount}, region={r})");
                previous = tmb;
            }
        });
    }

    // C2 — INV-04: classification flips from Low to High exactly at 10 mut/Mb.
    [TestCase(9.99, OncologyAnalyzer.TmbStatus.Low)]
    [TestCase(10.0, OncologyAnalyzer.TmbStatus.High)]
    [TestCase(10.01, OncologyAnalyzer.TmbStatus.High)]
    public void ClassifyTMB_BoundarySweep_FlipsOnlyAtTen(double tmb, OncologyAnalyzer.TmbStatus expected)
    {
        Assert.That(OncologyAnalyzer.ClassifyTMB(tmb), Is.EqualTo(expected),
            $"Classification flips only at the inclusive 10 mut/Mb cutoff (tmb={tmb})");
    }

    #endregion

    private static OncologyAnalyzer.SomaticCall MakeCall(OncologyAnalyzer.SomaticStatus status)
    {
        var variant = new OncologyAnalyzer.VariantObservation(
            "chr1", 100, "A", "T", 10, 20, 0, 20);
        return new OncologyAnalyzer.SomaticCall(variant, 0.5, 0.0, status, 0.5);
    }
}
