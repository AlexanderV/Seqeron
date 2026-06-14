// ONCO-VAF-001 — Variant Allele Frequency Analysis
// Evidence: docs/Evidence/ONCO-VAF-001-Evidence.md
// TestSpec: tests/TestSpecs/ONCO-VAF-001.md
// Source: Wilson E.B. (1927). JASA 22(158):209-212. https://doi.org/10.1080/01621459.1927.10502953
//         GATK Mutect2 FAQ (empirical allele fraction = alt/total). https://gatk.broadinstitute.org/hc/en-us/articles/360050722212-FAQ-for-Mutect2
//         Tarabichi et al. (2017). PMC5538405. CNAqc, Genome Biology (2024). https://doi.org/10.1186/s13059-024-03170-5

namespace Seqeron.Genomics.Tests;

[TestFixture]
public class OncologyAnalyzer_CalculateVAF_Tests
{
    #region CalculateVAF

    // M1/M2/M3/M4 — Empirical VAF = altReads / totalReads (GATK AD-based definition).
    [TestCase(25, 100, 0.25)]
    [TestCase(50, 100, 0.50)]
    [TestCase(10, 10, 1.00)]
    [TestCase(0, 10, 0.00)]
    public void CalculateVAF_ValidCounts_ReturnsAltOverTotal(int alt, int total, double expected)
    {
        double vaf = OncologyAnalyzer.CalculateVAF(alt, total);

        Assert.That(vaf, Is.EqualTo(expected).Within(1e-10),
            $"VAF must equal altReads/totalReads = {alt}/{total} = {expected}");
    }

    // M5 — Zero coverage: VAF defined as 0 (no reads => allele absent).
    [Test]
    public void CalculateVAF_ZeroCoverage_ReturnsZero()
    {
        Assert.That(OncologyAnalyzer.CalculateVAF(0, 0), Is.EqualTo(0.0).Within(1e-10),
            "An uncovered locus (totalReads=0) yields VAF 0, not a division error");
    }

    // M6 — altReads > totalReads (VAF > 1 alignment artifact) is invalid input.
    [Test]
    public void CalculateVAF_AltExceedsTotal_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.CalculateVAF(11, 10),
                "altReads cannot exceed totalReads (would give VAF > 1)");
    }

    // M7 — Negative counts are invalid.
    [Test]
    public void CalculateVAF_NegativeCounts_Throw()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.CalculateVAF(-1, 10),
                "Negative altReads is invalid");
            Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.CalculateVAF(0, -5),
                "Negative totalReads is invalid");
        });
    }

    #endregion

    #region CalculateVAFConfidenceInterval (Wilson score interval, 95%, z=1.96)

    // M8 — Wilson 95% CI for 25/100. Values computed from the source formula (Wilson 1927).
    [Test]
    public void CalculateVAFConfidenceInterval_25of100_MatchesWilsonFormula()
    {
        var ci = OncologyAnalyzer.CalculateVAFConfidenceInterval(25, 100);

        Assert.Multiple(() =>
        {
            Assert.That(ci.Vaf, Is.EqualTo(0.25).Within(1e-10), "Point estimate is 25/100");
            Assert.That(ci.Lower, Is.EqualTo(0.1754509400).Within(1e-9), "Wilson lower bound for p=0.25, n=100");
            Assert.That(ci.Upper, Is.EqualTo(0.3430464637).Within(1e-9), "Wilson upper bound for p=0.25, n=100");
            Assert.That(ci.Confidence, Is.EqualTo(0.95).Within(1e-10), "Reported confidence is 0.95");
        });
    }

    // M9 — Wilson 95% CI for 50/100. Center is exactly 0.5; bounds are symmetric about it.
    [Test]
    public void CalculateVAFConfidenceInterval_50of100_IsSymmetricAboutHalf()
    {
        var ci = OncologyAnalyzer.CalculateVAFConfidenceInterval(50, 100);

        Assert.Multiple(() =>
        {
            Assert.That(ci.Lower, Is.EqualTo(0.4038298286).Within(1e-9), "Wilson lower bound for p=0.5, n=100");
            Assert.That(ci.Upper, Is.EqualTo(0.5961701714).Within(1e-9), "Wilson upper bound for p=0.5, n=100");
            // Wilson center equals 0.5 when p=0.5 (symmetry); bounds are equidistant from 0.5.
            Assert.That((ci.Lower + ci.Upper) / 2.0, Is.EqualTo(0.5).Within(1e-9),
                "Midpoint of bounds is the Wilson center 0.5 at p=0.5");
        });
    }

    // M10 — Wilson 95% CI for small n (5/20).
    [Test]
    public void CalculateVAFConfidenceInterval_5of20_MatchesWilsonFormula()
    {
        var ci = OncologyAnalyzer.CalculateVAFConfidenceInterval(5, 20);

        Assert.Multiple(() =>
        {
            Assert.That(ci.Lower, Is.EqualTo(0.1118600528).Within(1e-9), "Wilson lower bound for p=0.25, n=20");
            Assert.That(ci.Upper, Is.EqualTo(0.4687050100).Within(1e-9), "Wilson upper bound for p=0.25, n=20");
        });
    }

    // M11 — No overshoot at p=0: lower bound is exactly 0, upper > 0 (non-zero width).
    [Test]
    public void CalculateVAFConfidenceInterval_ZeroAltReads_LowerBoundIsZeroNoOvershoot()
    {
        var ci = OncologyAnalyzer.CalculateVAFConfidenceInterval(0, 10);

        Assert.Multiple(() =>
        {
            Assert.That(ci.Lower, Is.EqualTo(0.0).Within(1e-12), "Wilson lower bound at p=0 is 0 (no overshoot)");
            Assert.That(ci.Upper, Is.EqualTo(0.2775401688).Within(1e-9), "Wilson upper bound for p=0, n=10");
        });
    }

    // M12 — No overshoot at p=1: upper bound is exactly 1, lower < 1 (non-zero width).
    [Test]
    public void CalculateVAFConfidenceInterval_AllAltReads_UpperBoundIsOneNoOvershoot()
    {
        var ci = OncologyAnalyzer.CalculateVAFConfidenceInterval(10, 10);

        Assert.Multiple(() =>
        {
            Assert.That(ci.Upper, Is.EqualTo(1.0).Within(1e-12), "Wilson upper bound at p=1 is 1 (no overshoot)");
            Assert.That(ci.Lower, Is.EqualTo(0.7224598312).Within(1e-9), "Wilson lower bound for p=1, n=10");
        });
    }

    // S1 — totalReads=0 has no defined interval.
    [Test]
    public void CalculateVAFConfidenceInterval_ZeroCoverage_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.CalculateVAFConfidenceInterval(0, 0),
                "An interval is undefined with zero trials");
    }

    // S1 — confidence outside (0,1) or unsupported level is rejected.
    [TestCase(0.0)]
    [TestCase(1.0)]
    [TestCase(1.5)]
    [TestCase(0.90)]
    public void CalculateVAFConfidenceInterval_UnsupportedConfidence_Throws(double confidence)
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.CalculateVAFConfidenceInterval(25, 100, confidence),
                "Only the source-cited 0.95 level (z=1.96) is supported");
    }

    #endregion

    #region AdjustVAFForPurity

    // M13 — Diploid clonal heterozygous: VAF 0.4 at purity 0.8 => 1.0 (Tarabichi 2017 worked example).
    [Test]
    public void AdjustVAFForPurity_DiploidHet_080_RecoversFullFraction()
    {
        double adjusted = OncologyAnalyzer.AdjustVAFForPurity(0.40, 0.80, 2.0);

        Assert.That(adjusted, Is.EqualTo(1.0).Within(1e-10),
            "Diploid het at 80% purity: expected VAF=0.4 => adjusted m*CCF=1.0");
    }

    // M14 — Diploid: VAF 0.2 at purity 0.5, ploidy 2 => 0.8 (CNAqc inversion).
    [Test]
    public void AdjustVAFForPurity_DiploidHalfPurity_ScalesByAverageCopies()
    {
        double adjusted = OncologyAnalyzer.AdjustVAFForPurity(0.20, 0.50, 2.0);

        Assert.That(adjusted, Is.EqualTo(0.8).Within(1e-10),
            "0.2 * (2*(1-0.5) + 0.5*2) / 0.5 = 0.2 * 2 / 0.5 = 0.8");
    }

    // M15 — Tetraploid segment: VAF 0.3 at purity 0.5, ploidy 4 => 1.8 (CNAqc inversion).
    [Test]
    public void AdjustVAFForPurity_TetraploidSegment_AccountsForCopyNumber()
    {
        double adjusted = OncologyAnalyzer.AdjustVAFForPurity(0.30, 0.50, 4.0);

        Assert.That(adjusted, Is.EqualTo(1.8).Within(1e-10),
            "0.3 * (2*(1-0.5) + 0.5*4) / 0.5 = 0.3 * 3 / 0.5 = 1.8");
    }

    // M16 — purity=0 makes the correction (division by purity) undefined.
    [Test]
    public void AdjustVAFForPurity_ZeroPurity_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.AdjustVAFForPurity(0.3, 0.0, 2.0),
                "Purity 0 (pure normal) divides by zero");
    }

    // S2/S3 — input validation for vaf and ploidy.
    [Test]
    public void AdjustVAFForPurity_InvalidInputs_Throw()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.AdjustVAFForPurity(-0.1, 0.8, 2.0),
                "VAF below 0 is invalid");
            Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.AdjustVAFForPurity(1.1, 0.8, 2.0),
                "VAF above 1 is invalid");
            Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.AdjustVAFForPurity(0.3, 1.5, 2.0),
                "Purity above 1 is invalid");
            Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.AdjustVAFForPurity(0.3, 0.8, 0.0),
                "Ploidy must be positive");
        });
    }

    // INV-04 — diploid het round-trip: AdjustVAFForPurity(pi/2, pi, 2) = 1 for any pi in (0,1].
    [TestCase(0.10)]
    [TestCase(0.37)]
    [TestCase(0.95)]
    [TestCase(1.00)]
    public void AdjustVAFForPurity_DiploidHetRoundTrip_IsOne(double purity)
    {
        double expectedVaf = purity / 2.0; // diploid het expected VAF = pi/2 (CNAqc/Tarabichi)
        double adjusted = OncologyAnalyzer.AdjustVAFForPurity(expectedVaf, purity, 2.0);

        Assert.That(adjusted, Is.EqualTo(1.0).Within(1e-10),
            "Diploid het expected VAF pi/2 must invert to mutant fraction 1.0");
    }

    #endregion

    #region Invariants (property-based)

    // C1 — INV-01/02/03: 0<=VAF<=1, lower<=center<=upper, bounds within [0,1] over a sweep.
    [Test]
    public void CalculateVAFConfidenceInterval_InvariantsHoldOverSweep()
    {
        Assert.Multiple(() =>
        {
            for (int total = 1; total <= 50; total++)
            {
                for (int alt = 0; alt <= total; alt++)
                {
                    var ci = OncologyAnalyzer.CalculateVAFConfidenceInterval(alt, total);
                    Assert.That(ci.Vaf, Is.InRange(0.0, 1.0), $"INV-01 VAF in [0,1] for {alt}/{total}");
                    Assert.That(ci.Lower, Is.LessThanOrEqualTo(ci.Vaf + 1e-12).And.GreaterThanOrEqualTo(-1e-12),
                        $"INV-02/03 lower<=vaf and >=0 for {alt}/{total}");
                    Assert.That(ci.Upper, Is.GreaterThanOrEqualTo(ci.Vaf - 1e-12).And.LessThanOrEqualTo(1.0 + 1e-12),
                        $"INV-02/03 upper>=vaf and <=1 for {alt}/{total}");
                }
            }
        });
    }

    #endregion
}
