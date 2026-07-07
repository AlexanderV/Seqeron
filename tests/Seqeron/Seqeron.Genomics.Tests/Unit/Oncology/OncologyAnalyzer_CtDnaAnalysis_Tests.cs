// ONCO-CTDNA-001 — ctDNA Analysis
// Evidence: docs/Evidence/ONCO-CTDNA-001-Evidence.md
// TestSpec: tests/TestSpecs/ONCO-CTDNA-001.md
// Source: Newman et al. (2014). Nat Med 20(5):548-554 (CAPP-Seq; detection range, mean-VAF reporters).
//         Avanzini et al. (2020). Sci Adv 6(50):eabc4308 / US Patent 11,085,084 B2 (Poisson p = 1 - e^(-ndk)).
//         Devonshire et al. (2014) PMC4182654 (3.3 pg/haploid genome); Alcaide et al. (2020) Sci Rep 10:12564.
//         Antonello et al. (2024). Genome Biology 25:38 (clonal het diploid: tumour fraction = 2 * VAF).
//         Pessoa et al. (2023) PMC10314661 (n=15000, d=0.001 => 15 mutant molecules).

namespace Seqeron.Genomics.Tests.Unit.Oncology;

[TestFixture]
public class OncologyAnalyzer_CtDnaAnalysis_Tests
{
    // Helper: a variant whose tumour (plasma) VAF = altReads/totalReads. Normal counts unused here.
    private static OncologyAnalyzer.VariantObservation Reporter(int altReads, int totalReads) =>
        new("1", 100, "A", "T", altReads, totalReads, 0, 100);

    #region CtDnaDetectionProbability

    // M1 — Patent US11085084: p = 1 - e^(-ndk). n=15000, d=0.001, k=1 => lambda=15 (Pessoa 2023: 15 molecules).
    [Test]
    public void DetectionProbability_WorkedExample_15Molecules()
    {
        // Arrange / Act
        double p = OncologyAnalyzer.CtDnaDetectionProbability(15000, 0.001, 1);

        // Assert — 1 - e^(-15) = 0.9999996939215850 (hand-derived from the cited Poisson formula).
        Assert.That(p, Is.EqualTo(1.0 - Math.Exp(-15.0)).Within(1e-12),
            "p must equal 1 - e^(-n*d*k) with lambda = 15000*0.001*1 = 15 (Patent US11085084).");
        Assert.That(p, Is.EqualTo(0.99999969409767953).Within(1e-12),
            "Numeric value of 1 - e^(-15) per the cited Poisson detection model.");
    }

    // M2 — single expected molecule: n=1000, d=0.001, k=1 => lambda=1 => p = 1 - e^(-1) = 0.6321205588.
    [Test]
    public void DetectionProbability_OneExpectedMolecule_Returns1MinusEInverse()
    {
        double p = OncologyAnalyzer.CtDnaDetectionProbability(1000, 0.001, 1);

        Assert.That(p, Is.EqualTo(0.6321205588285577).Within(1e-12),
            "lambda = 1 => p = 1 - e^(-1) = 0.63212... (Patent US11085084).");
    }

    // M3 — k reporters raise p: n=1000, d=0.001, k=10 => lambda=10 => p = 1 - e^(-10), strictly > M2.
    [Test]
    public void DetectionProbability_TenReporters_HigherThanSingle()
    {
        double pTen = OncologyAnalyzer.CtDnaDetectionProbability(1000, 0.001, 10);
        double pOne = OncologyAnalyzer.CtDnaDetectionProbability(1000, 0.001, 1);

        Assert.Multiple(() =>
        {
            Assert.That(pTen, Is.EqualTo(0.9999546000702375).Within(1e-12),
                "lambda = 10 => p = 1 - e^(-10) (Patent US11085084: p = 1 - e^(-ndk)).");
            Assert.That(pTen, Is.GreaterThan(pOne),
                "More independent reporters must increase the detection probability (INV-03).");
        });
    }

    // M4 — lambda = 0 => p = 0 (n=0 or d=0). INV-02.
    [Test]
    public void DetectionProbability_ZeroLambda_ReturnsZero()
    {
        Assert.Multiple(() =>
        {
            Assert.That(OncologyAnalyzer.CtDnaDetectionProbability(0, 0.5, 5), Is.EqualTo(0.0).Within(1e-15),
                "n = 0 => lambda = 0 => p = 1 - e^0 = 0 (INV-02).");
            Assert.That(OncologyAnalyzer.CtDnaDetectionProbability(10000, 0.0, 5), Is.EqualTo(0.0).Within(1e-15),
                "d = 0 => lambda = 0 => p = 0 (INV-02).");
        });
    }

    // M12 — allele fraction out of [0,1] throws.
    [Test]
    public void DetectionProbability_AlleleFractionOutOfRange_Throws()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.CtDnaDetectionProbability(100, -0.1, 1),
                "d < 0 is not a valid allele fraction.");
            Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.CtDnaDetectionProbability(100, 1.1, 1),
                "d > 1 is not a valid allele fraction.");
            Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.CtDnaDetectionProbability(100, double.NaN, 1),
                "d = NaN is invalid.");
        });
    }

    // S1 — negative genome equivalents throws.
    [Test]
    public void DetectionProbability_NegativeGenomeEquivalents_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.CtDnaDetectionProbability(-1, 0.01, 1),
            "Genome equivalents (n) cannot be negative.");
    }

    // S2 — reporter count < 1 throws.
    [Test]
    public void DetectionProbability_ReporterCountBelowOne_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.CtDnaDetectionProbability(100, 0.01, 0),
            "Reporter count (k) must be at least 1.");
    }

    // C1 — strictly increasing in d at fixed n,k (INV-03).
    [Test]
    public void DetectionProbability_IncreasingAlleleFraction_IsMonotone()
    {
        double pLow = OncologyAnalyzer.CtDnaDetectionProbability(1000, 0.001, 1);
        double pHigh = OncologyAnalyzer.CtDnaDetectionProbability(1000, 0.002, 1);

        Assert.That(pHigh, Is.GreaterThan(pLow),
            "p must strictly increase with mutant allele fraction d (INV-03).");
    }

    // C2 — bounded by 1 for large lambda (INV-01).
    [Test]
    public void DetectionProbability_LargeLambda_BoundedByOne()
    {
        double p = OncologyAnalyzer.CtDnaDetectionProbability(1_000_000, 0.5, 5);

        Assert.Multiple(() =>
        {
            Assert.That(p, Is.LessThanOrEqualTo(1.0), "p must never exceed 1 (INV-01).");
            Assert.That(p, Is.EqualTo(1.0).Within(1e-12), "For huge lambda, p -> 1 (INV-01).");
        });
    }

    #endregion

    #region IsCtDnaDetected

    // M5 — above LoD: lambda=15>=1 and p~1>=0.95 => detected.
    [Test]
    public void IsCtDnaDetected_AboveLod_ReturnsTrue()
    {
        Assert.That(OncologyAnalyzer.IsCtDnaDetected(15000, 0.001, 1, 0.95), Is.True,
            "lambda = 15 (>= 1) and p ~ 1 (>= 0.95) => detected (Patent US11085084; Newman 2014).");
    }

    // M6 — below detection: n=100, d=0.0001 => lambda=0.01 < 1 => not detected (p=0.00995 << 0.95).
    [Test]
    public void IsCtDnaDetected_BelowOneMolecule_ReturnsFalse()
    {
        Assert.That(OncologyAnalyzer.IsCtDnaDetected(100, 0.0001, 1, 0.95), Is.False,
            "lambda = 0.01 < 1 mutant molecule => not detectable regardless of probability threshold.");
    }

    // M5b — lambda >= 1 but probability below threshold => not detected (the second AND-condition).
    // n=1000, d=0.001, k=1 => lambda=1 (>= 1) but p = 1 - e^(-1) = 0.6321... < 0.95 default threshold.
    [Test]
    public void IsCtDnaDetected_LambdaAtLeastOneButProbabilityBelowThreshold_ReturnsFalse()
    {
        // lambda = 1 satisfies the >=1 molecule floor, yet p = 0.6321205588 < 0.95 (Patent US11085084 formula).
        Assert.That(OncologyAnalyzer.IsCtDnaDetected(1000, 0.001, 1), Is.False,
            "lambda = 1 (>= 1) but p = 1 - e^(-1) = 0.6321 < 0.95 default => not detected.");
        // The same inputs ARE detected if the caller lowers the threshold below the actual probability.
        Assert.That(OncologyAnalyzer.IsCtDnaDetected(1000, 0.001, 1, 0.5), Is.True,
            "Same lambda=1, p=0.6321 >= 0.5 caller threshold => detected.");
    }

    // S6 — minDetectionProbability outside (0, 1] throws (documented contract bound).
    [Test]
    public void IsCtDnaDetected_ThresholdOutOfRange_Throws()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentOutOfRangeException>(
                () => OncologyAnalyzer.IsCtDnaDetected(15000, 0.001, 1, 0.0),
                "Threshold must be > 0 (a 0 threshold is not a valid operating point).");
            Assert.Throws<ArgumentOutOfRangeException>(
                () => OncologyAnalyzer.IsCtDnaDetected(15000, 0.001, 1, 1.1),
                "Threshold must be <= 1 (a probability cannot exceed 1).");
            Assert.Throws<ArgumentOutOfRangeException>(
                () => OncologyAnalyzer.IsCtDnaDetected(15000, 0.001, 1, double.NaN),
                "Threshold NaN is invalid.");
        });
    }

    #endregion

    #region ExpectedMutantMolecules

    // Supports M1/M5/M6: lambda = n*d*k exactly (Pessoa 2023 worked example: 15000*0.001=15).
    [Test]
    public void ExpectedMutantMolecules_WorkedExample_Returns15()
    {
        Assert.That(OncologyAnalyzer.ExpectedMutantMolecules(15000, 0.001, 1), Is.EqualTo(15.0).Within(1e-12),
            "lambda = n*d*k = 15000*0.001*1 = 15 (Pessoa 2023).");
    }

    // lambda is multiplicative in k: 15000*0.001*10 = 150 (Patent US11085084: mean lambda = n*d*k).
    [Test]
    public void ExpectedMutantMolecules_TenReporters_Returns150()
    {
        Assert.That(OncologyAnalyzer.ExpectedMutantMolecules(15000, 0.001, 10), Is.EqualTo(150.0).Within(1e-12),
            "lambda = n*d*k = 15000*0.001*10 = 150.");
    }

    // Domain guards mirror CtDnaDetectionProbability (negative n, d out of [0,1], k < 1).
    [Test]
    public void ExpectedMutantMolecules_InvalidArguments_Throw()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.ExpectedMutantMolecules(-1, 0.01, 1),
                "Genome equivalents (n) cannot be negative.");
            Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.ExpectedMutantMolecules(100, 1.1, 1),
                "Mutant allele fraction (d) must be in [0, 1].");
            Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.ExpectedMutantMolecules(100, double.NaN, 1),
                "Mutant allele fraction NaN is invalid.");
            Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.ExpectedMutantMolecules(100, 0.01, 0),
                "Reporter count (k) must be at least 1.");
        });
    }

    #endregion

    #region CalculateTumorFraction

    // M7 — tumour fraction = 2 * mean VAF. Two clonal het SNVs VAF 0.10 & 0.20 => mean 0.15 => TF 0.30.
    [Test]
    public void CalculateTumorFraction_TwoClonalHetSnvs_ReturnsTwiceMeanVaf()
    {
        var variants = new[] { Reporter(10, 100), Reporter(20, 100) }; // VAF 0.10, 0.20

        double tf = OncologyAnalyzer.CalculateTumorFraction(variants);

        Assert.That(tf, Is.EqualTo(0.30).Within(1e-10),
            "TF = 2 * mean(0.10, 0.20) = 2 * 0.15 = 0.30 (Antonello 2024: v = pi/2 => TF = 2v).");
    }

    // M8 — TF clamped to 1.0. Use VAF 0.5 (max for het diploid) twice => 2*0.5 = 1.0 (boundary, not clamped beyond).
    //      To exercise the clamp, mix VAF 0.5 and 0.5: mean 0.5 => 1.0. (>0.5 per-variant is rejected, see S4.)
    [Test]
    public void CalculateTumorFraction_MaxVafs_ClampedToOne()
    {
        var variants = new[] { Reporter(50, 100), Reporter(50, 100) }; // VAF 0.50, 0.50 => 2*0.5 = 1.0

        double tf = OncologyAnalyzer.CalculateTumorFraction(variants);

        Assert.That(tf, Is.EqualTo(1.0).Within(1e-10),
            "2 * mean(0.5, 0.5) = 1.0; tumour fraction cannot exceed 1 (INV-04).");
    }

    // M13 — null variants throws.
    [Test]
    public void CalculateTumorFraction_Null_Throws()
    {
        Assert.Throws<ArgumentNullException>(
            () => OncologyAnalyzer.CalculateTumorFraction(null!),
            "Null variant collection is invalid.");
    }

    // M14 — empty variants throws (TF undefined).
    [Test]
    public void CalculateTumorFraction_Empty_Throws()
    {
        Assert.Throws<ArgumentException>(
            () => OncologyAnalyzer.CalculateTumorFraction(new List<OncologyAnalyzer.VariantObservation>()),
            "Tumour fraction is undefined for an empty variant set.");
    }

    // S4 — a per-variant VAF > 0.5 is impossible for a diploid het SNV => throws.
    [Test]
    public void CalculateTumorFraction_VafAboveHalf_Throws()
    {
        var variants = new[] { Reporter(60, 100) }; // VAF 0.60 > 0.5

        Assert.Throws<ArgumentOutOfRangeException>(
            () => OncologyAnalyzer.CalculateTumorFraction(variants),
            "A clonal heterozygous diploid SNV cannot have VAF > 0.5 (Antonello 2024).");
    }

    #endregion

    #region CalculateMeanVaf

    // M9 — mean of per-reporter VAFs: (0.05 + 0.30 + 0.01)/3 = 0.12.
    [Test]
    public void CalculateMeanVaf_ThreeReporters_ReturnsArithmeticMean()
    {
        var variants = new[] { Reporter(5, 100), Reporter(30, 100), Reporter(1, 100) };

        double mean = OncologyAnalyzer.CalculateMeanVaf(variants);

        Assert.That(mean, Is.EqualTo(0.12).Within(1e-10),
            "mean VAF = (0.05 + 0.30 + 0.01)/3 = 0.12 (Newman 2014: fraction across reporters).");
    }

    // S3 — null and empty throw.
    [Test]
    public void CalculateMeanVaf_NullAndEmpty_Throws()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentNullException>(() => OncologyAnalyzer.CalculateMeanVaf(null!),
                "Null variant collection is invalid.");
            Assert.Throws<ArgumentException>(
                () => OncologyAnalyzer.CalculateMeanVaf(new List<OncologyAnalyzer.VariantObservation>()),
                "Mean VAF is undefined for an empty variant set.");
        });
    }

    #endregion

    #region HaploidGenomeEquivalents

    // M10 — 1 ng => 1000/3.3 = 303.0303... haploid genome equivalents (Devonshire 2014; Alcaide 2020).
    [Test]
    public void HaploidGenomeEquivalents_OneNanogram_Returns303()
    {
        double ge = OncologyAnalyzer.HaploidGenomeEquivalents(1.0);

        Assert.That(ge, Is.EqualTo(1000.0 / 3.3).Within(1e-10),
            "1 ng / 3.3 pg-per-haploid-genome = 1000/3.3 = 303.0303... GE (Devonshire 2014).");
        Assert.That(ge, Is.EqualTo(303.030303030303).Within(1e-9),
            "Numeric value ~303 GE per ng (Alcaide 2020).");
    }

    // M11 — 0 ng => 0 GE (INV-05).
    [Test]
    public void HaploidGenomeEquivalents_Zero_ReturnsZero()
    {
        Assert.That(OncologyAnalyzer.HaploidGenomeEquivalents(0.0), Is.EqualTo(0.0).Within(1e-15),
            "0 ng yields 0 genome equivalents (INV-05).");
    }

    // S5 — negative mass throws.
    [Test]
    public void HaploidGenomeEquivalents_Negative_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(
            () => OncologyAnalyzer.HaploidGenomeEquivalents(-1.0),
            "cfDNA mass cannot be negative.");
    }

    #endregion
}
