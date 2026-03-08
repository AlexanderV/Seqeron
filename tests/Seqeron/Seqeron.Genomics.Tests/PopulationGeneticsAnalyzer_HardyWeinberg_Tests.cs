using NUnit.Framework;
using Seqeron.Genomics;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Tests for Hardy-Weinberg Equilibrium testing.
/// Test Unit: POP-HW-001
/// </summary>
/// <remarks>
/// Evidence Sources:
/// - Wikipedia: Hardy-Weinberg principle (https://en.wikipedia.org/wiki/Hardy-Weinberg_principle)
/// - Ford (1971): Ecological Genetics - Scarlet tiger moth dataset
/// - Hardy (1908), Weinberg (1908): Original formulation
/// </remarks>
[TestFixture]
public class PopulationGeneticsAnalyzer_HardyWeinberg_Tests
{
    private const double Tolerance = 0.01;
    private const double ChiSquareCriticalValue_Alpha05_Df1 = 3.841;

    #region Published Dataset Tests (Evidence-Based)

    /// <summary>
    /// HW-M01: Ford's Scarlet Tiger Moth data from Wikipedia.
    /// Source: Ford (1971) Ecological Genetics, cited in Wikipedia HWE article.
    /// Observed: AA=1469, Aa=138, aa=5, n=1612
    /// Expected χ² ≈ 0.83, population IS in HWE.
    /// </summary>
    [Test]
    public void TestHardyWeinberg_FordsMothData_IsInEquilibrium()
    {
        var result = PopulationGeneticsAnalyzer.TestHardyWeinberg(
            variantId: "FORD_MOTH",
            observedAA: 1469,
            observedAa: 138,
            observedaa: 5);

        Assert.Multiple(() =>
        {
            // Verify chi-square: exact = 0.8309 (Wikipedia rounds to 0.83)
            Assert.That(result.ChiSquare, Is.EqualTo(0.8309).Within(0.01),
                "Chi-square should be approximately 0.8309 for Ford's data");

            // Verify equilibrium decision
            Assert.That(result.InEquilibrium, Is.True,
                "Ford's moth data should be in Hardy-Weinberg equilibrium");

            // Verify expected counts (Wikipedia: E_AA≈1467.4, E_Aa≈141.2, E_aa≈3.4)
            Assert.That(result.ExpectedAA, Is.EqualTo(1467.40).Within(0.1),
                "Expected AA count should be approximately 1467.40");
            Assert.That(result.ExpectedAa, Is.EqualTo(141.21).Within(0.1),
                "Expected Aa count should be approximately 141.21");
            Assert.That(result.Expectedaa, Is.EqualTo(3.40).Within(0.1),
                "Expected aa count should be approximately 3.40");

            // P-value should be > 0.05 (well above threshold)
            Assert.That(result.PValue, Is.GreaterThan(0.05),
                "P-value should exceed 0.05 significance level");
        });
    }

    /// <summary>
    /// HW-M02: Perfect Hardy-Weinberg equilibrium with p=0.5.
    /// When observed genotype counts exactly match expected frequencies.
    /// </summary>
    [Test]
    public void TestHardyWeinberg_PerfectEquilibrium_ChiSquareNearZero()
    {
        // For p=0.5, n=100: E(AA)=25, E(Aa)=50, E(aa)=25
        var result = PopulationGeneticsAnalyzer.TestHardyWeinberg(
            variantId: "PERFECT_HWE",
            observedAA: 25,
            observedAa: 50,
            observedaa: 25);

        Assert.Multiple(() =>
        {
            Assert.That(result.ChiSquare, Is.EqualTo(0).Within(Tolerance),
                "Chi-square should be zero for perfect HWE");
            Assert.That(result.InEquilibrium, Is.True,
                "Perfect HWE should report equilibrium");
            Assert.That(result.ExpectedAA, Is.EqualTo(25).Within(Tolerance));
            Assert.That(result.ExpectedAa, Is.EqualTo(50).Within(Tolerance));
            Assert.That(result.Expectedaa, Is.EqualTo(25).Within(Tolerance));
        });
    }

    #endregion

    #region Deviation Detection Tests

    /// <summary>
    /// HW-M03: Excess heterozygotes pattern.
    /// Source: Wikipedia - deviation from HWE example.
    /// </summary>
    [Test]
    public void TestHardyWeinberg_ExcessHeterozygotes_DeviatesFromEquilibrium()
    {
        // Observed: AA=10, Aa=80, aa=10, n=100, p=0.5
        // Expected: AA=25, Aa=50, aa=25
        // χ² = (10-25)²/25 + (80-50)²/50 + (10-25)²/25 = 9 + 18 + 9 = 36
        var result = PopulationGeneticsAnalyzer.TestHardyWeinberg(
            variantId: "EXCESS_HET",
            observedAA: 10,
            observedAa: 80,
            observedaa: 10);

        Assert.Multiple(() =>
        {
            Assert.That(result.ChiSquare, Is.EqualTo(36.0).Within(0.01),
                "Chi-square should be exactly 36: (10-25)²/25 + (80-50)²/50 + (10-25)²/25 = 9+18+9");
            Assert.That(result.ChiSquare, Is.GreaterThan(ChiSquareCriticalValue_Alpha05_Df1),
                "Chi-square should exceed critical value for excess heterozygotes");
            Assert.That(result.InEquilibrium, Is.False,
                "Excess heterozygotes should deviate from HWE");
            Assert.That(result.PValue, Is.LessThan(0.05),
                "P-value should be less than 0.05");
        });
    }

    /// <summary>
    /// HW-M04: All heterozygotes - extreme deviation.
    /// This is biologically implausible under random mating.
    /// p = (0 + 100) / 200 = 0.5, q = 0.5
    /// E(AA) = 25, E(Aa) = 50, E(aa) = 25
    /// χ² = (0-25)²/25 + (100-50)²/50 + (0-25)²/25 = 25 + 50 + 25 = 100
    /// </summary>
    [Test]
    public void TestHardyWeinberg_AllHeterozygotes_SignificantDeviation()
    {
        var result = PopulationGeneticsAnalyzer.TestHardyWeinberg(
            variantId: "ALL_HET",
            observedAA: 0,
            observedAa: 100,
            observedaa: 0);

        Assert.Multiple(() =>
        {
            Assert.That(result.ChiSquare, Is.EqualTo(100.0).Within(0.01),
                "Chi-square should be exactly 100 for all-heterozygote case");
            Assert.That(result.InEquilibrium, Is.False,
                "All heterozygotes cannot be in HWE");
            Assert.That(result.PValue, Is.LessThan(0.05),
                "P-value should be far below 0.05");
            Assert.That(result.ExpectedAA, Is.EqualTo(25.0).Within(Tolerance),
                "E(AA) = 0.5² × 100 = 25");
            Assert.That(result.ExpectedAa, Is.EqualTo(50.0).Within(Tolerance),
                "E(Aa) = 2×0.5×0.5 × 100 = 50");
            Assert.That(result.Expectedaa, Is.EqualTo(25.0).Within(Tolerance),
                "E(aa) = 0.5² × 100 = 25");
        });
    }

    /// <summary>
    /// HW-M05: Deficit of heterozygotes - inbreeding pattern.
    /// p = (2×45 + 10) / 200 = 0.5, q = 0.5
    /// E(AA) = 25, E(Aa) = 50, E(aa) = 25
    /// χ² = (45-25)²/25 + (10-50)²/50 + (45-25)²/25 = 16 + 32 + 16 = 64
    /// </summary>
    [Test]
    public void TestHardyWeinberg_DeficitHeterozygotes_InbreedingPattern()
    {
        // p = 0.5, Expected: AA=25, Aa=50, aa=25
        // Observed: AA=45, Aa=10, aa=45 (excess homozygotes)
        var result = PopulationGeneticsAnalyzer.TestHardyWeinberg(
            variantId: "DEFICIT_HET",
            observedAA: 45,
            observedAa: 10,
            observedaa: 45);

        Assert.Multiple(() =>
        {
            Assert.That(result.ChiSquare, Is.EqualTo(64.0).Within(0.01),
                "Chi-square should be exactly 64 for deficit heterozygotes");
            Assert.That(result.InEquilibrium, Is.False,
                "Inbreeding pattern should deviate from HWE");
            Assert.That(result.PValue, Is.LessThan(0.05),
                "P-value should be far below 0.05");
            Assert.That(result.ExpectedAA, Is.EqualTo(25.0).Within(Tolerance),
                "E(AA) = 0.5² × 100 = 25");
            Assert.That(result.ExpectedAa, Is.EqualTo(50.0).Within(Tolerance),
                "E(Aa) = 2×0.5×0.5 × 100 = 50");
            Assert.That(result.Expectedaa, Is.EqualTo(25.0).Within(Tolerance),
                "E(aa) = 0.5² × 100 = 25");
        });
    }

    #endregion

    #region Expected Count Verification Tests

    /// <summary>
    /// HW-M06: Fixed major allele (all AA).
    /// p = 1.0, q = 0.0
    /// </summary>
    [Test]
    public void TestHardyWeinberg_FixedMajorAllele_ExpectedCountsCorrect()
    {
        var result = PopulationGeneticsAnalyzer.TestHardyWeinberg(
            variantId: "FIXED_MAJOR",
            observedAA: 100,
            observedAa: 0,
            observedaa: 0);

        Assert.Multiple(() =>
        {
            // p = 1.0: Expected = (100, 0, 0)
            Assert.That(result.ExpectedAA, Is.EqualTo(100).Within(Tolerance),
                "E(AA) should be 100 when p=1.0");
            Assert.That(result.ExpectedAa, Is.EqualTo(0).Within(Tolerance),
                "E(Aa) should be 0 when p=1.0");
            Assert.That(result.Expectedaa, Is.EqualTo(0).Within(Tolerance),
                "E(aa) should be 0 when p=1.0");
            Assert.That(result.ChiSquare, Is.EqualTo(0).Within(Tolerance),
                "Chi-square should be 0 when observed equals expected");
            Assert.That(result.InEquilibrium, Is.True,
                "Fixed allele should be in equilibrium");
        });
    }

    /// <summary>
    /// HW-M07: Fixed minor allele (all aa).
    /// p = 0.0, q = 1.0
    /// </summary>
    [Test]
    public void TestHardyWeinberg_FixedMinorAllele_ExpectedCountsCorrect()
    {
        var result = PopulationGeneticsAnalyzer.TestHardyWeinberg(
            variantId: "FIXED_MINOR",
            observedAA: 0,
            observedAa: 0,
            observedaa: 100);

        Assert.Multiple(() =>
        {
            Assert.That(result.ExpectedAA, Is.EqualTo(0).Within(Tolerance));
            Assert.That(result.ExpectedAa, Is.EqualTo(0).Within(Tolerance));
            Assert.That(result.Expectedaa, Is.EqualTo(100).Within(Tolerance));
            Assert.That(result.ChiSquare, Is.EqualTo(0).Within(Tolerance));
            Assert.That(result.InEquilibrium, Is.True);
        });
    }

    /// <summary>
    /// HW-M08: Verify expected counts follow HWE formulas.
    /// Using p=0.6 (from genotype counts 36, 48, 16).
    /// </summary>
    [Test]
    public void TestHardyWeinberg_VerifyExpectedCountFormulas()
    {
        // Observed: AA=36, Aa=48, aa=16, n=100
        // p = (2*36 + 48) / 200 = 120/200 = 0.6
        // q = 0.4
        // E(AA) = 0.36 * 100 = 36
        // E(Aa) = 0.48 * 100 = 48
        // E(aa) = 0.16 * 100 = 16
        var result = PopulationGeneticsAnalyzer.TestHardyWeinberg(
            variantId: "P60_Q40",
            observedAA: 36,
            observedAa: 48,
            observedaa: 16);

        const double p = 0.6;
        const double q = 0.4;
        const int n = 100;

        Assert.Multiple(() =>
        {
            Assert.That(result.ExpectedAA, Is.EqualTo(p * p * n).Within(Tolerance),
                "E(AA) should equal p² × n");
            Assert.That(result.ExpectedAa, Is.EqualTo(2 * p * q * n).Within(Tolerance),
                "E(Aa) should equal 2pq × n");
            Assert.That(result.Expectedaa, Is.EqualTo(q * q * n).Within(Tolerance),
                "E(aa) should equal q² × n");
        });
    }

    #endregion

    #region Edge Case Tests

    /// <summary>
    /// HW-M09: Zero samples edge case.
    /// </summary>
    [Test]
    public void TestHardyWeinberg_ZeroSamples_ReturnsEquilibrium()
    {
        var result = PopulationGeneticsAnalyzer.TestHardyWeinberg(
            variantId: "ZERO_SAMPLES",
            observedAA: 0,
            observedAa: 0,
            observedaa: 0);

        Assert.Multiple(() =>
        {
            Assert.That(result.InEquilibrium, Is.True,
                "Zero samples should return equilibrium (no evidence against)");
            Assert.That(result.PValue, Is.EqualTo(1),
                "P-value should be 1 for zero samples");
            Assert.That(result.ChiSquare, Is.EqualTo(0),
                "Chi-square should be 0 for zero samples");
        });
    }

    /// <summary>
    /// HW-M10: Single sample - minimum valid input.
    /// </summary>
    [Test]
    public void TestHardyWeinberg_SingleSample_ReturnsValidResult()
    {
        // AA=1, Aa=0, aa=0 → p=1.0, q=0.0
        // E(AA)=1, E(Aa)=0, E(aa)=0; Observed=Expected → χ²=0, PValue=1
        var result = PopulationGeneticsAnalyzer.TestHardyWeinberg(
            variantId: "SINGLE",
            observedAA: 1,
            observedAa: 0,
            observedaa: 0);

        Assert.Multiple(() =>
        {
            Assert.That(result.ObservedAA, Is.EqualTo(1));
            Assert.That(result.ExpectedAA, Is.EqualTo(1).Within(Tolerance));
            Assert.That(result.ExpectedAa, Is.EqualTo(0).Within(Tolerance));
            Assert.That(result.Expectedaa, Is.EqualTo(0).Within(Tolerance));
            Assert.That(result.ChiSquare, Is.EqualTo(0).Within(Tolerance),
                "Chi-square should be 0 when observed equals expected");
            Assert.That(result.PValue, Is.EqualTo(1).Within(Tolerance),
                "P-value should be 1 when χ²=0");
            Assert.That(result.InEquilibrium, Is.True,
                "Single monomorphic sample should be in equilibrium");
        });
    }

    /// <summary>
    /// HW-M11: VariantId is preserved in result.
    /// </summary>
    [Test]
    public void TestHardyWeinberg_VariantIdPreserved()
    {
        const string testId = "rs12345";

        var result = PopulationGeneticsAnalyzer.TestHardyWeinberg(
            variantId: testId,
            observedAA: 50,
            observedAa: 40,
            observedaa: 10);

        Assert.That(result.VariantId, Is.EqualTo(testId),
            "VariantId should be preserved in result");
    }

    #endregion

    #region Invariant Tests

    /// <summary>
    /// HW-M12: Chi-square is always non-negative.
    /// </summary>
    [TestCase(10, 80, 10)]
    [TestCase(25, 50, 25)]
    [TestCase(100, 0, 0)]
    [TestCase(0, 100, 0)]
    [TestCase(45, 10, 45)]
    [TestCase(1, 1, 1)]
    public void TestHardyWeinberg_ChiSquareNonNegative(int observedAA, int observedAa, int observedaa)
    {
        var result = PopulationGeneticsAnalyzer.TestHardyWeinberg("INV_TEST", observedAA, observedAa, observedaa);

        Assert.That(result.ChiSquare, Is.GreaterThanOrEqualTo(0),
            "Chi-square must always be non-negative");
    }

    /// <summary>
    /// HW-M13: P-value is in valid probability range [0, 1].
    /// </summary>
    [TestCase(10, 80, 10)]
    [TestCase(25, 50, 25)]
    [TestCase(100, 0, 0)]
    [TestCase(0, 100, 0)]
    [TestCase(45, 10, 45)]
    public void TestHardyWeinberg_PValueInValidRange(int observedAA, int observedAa, int observedaa)
    {
        var result = PopulationGeneticsAnalyzer.TestHardyWeinberg("PVAL_TEST", observedAA, observedAa, observedaa);

        Assert.That(result.PValue, Is.InRange(0, 1),
            "P-value must be in [0, 1]");
    }

    /// <summary>
    /// HW-M14: InEquilibrium is consistent with PValue and significance level.
    /// </summary>
    [TestCase(25, 50, 25, 0.05)]   // Perfect HWE: p=1.0, InEquilibrium=true
    [TestCase(10, 80, 10, 0.05)]   // Excess heterozygotes: p≈0, InEquilibrium=false
    [TestCase(10, 80, 10, 0.001)]  // Even stricter α, still fails
    [TestCase(46, 49, 5, 0.05)]    // Borderline: p≈0.075, InEquilibrium=true
    [TestCase(46, 49, 5, 0.10)]    // Borderline: p≈0.075, InEquilibrium=false
    public void TestHardyWeinberg_EquilibriumConsistentWithPValue(
        int observedAA, int observedAa, int observedaa, double significance)
    {
        var result = PopulationGeneticsAnalyzer.TestHardyWeinberg(
            "CONSIST_TEST", observedAA, observedAa, observedaa, significance);

        // Invariant: InEquilibrium should be true iff PValue >= significance
        bool expectedFromPValue = result.PValue >= significance;
        Assert.That(result.InEquilibrium, Is.EqualTo(expectedFromPValue),
            $"InEquilibrium should equal (PValue >= significance). PValue={result.PValue}, sig={significance}");
    }

    /// <summary>
    /// HW-M15: Expected counts sum to total sample size.
    /// </summary>
    [TestCase(36, 48, 16)]
    [TestCase(100, 100, 100)]
    [TestCase(1, 2, 3)]
    public void TestHardyWeinberg_ExpectedCountsSumToN(int observedAA, int observedAa, int observedaa)
    {
        int n = observedAA + observedAa + observedaa;
        var result = PopulationGeneticsAnalyzer.TestHardyWeinberg("SUM_TEST", observedAA, observedAa, observedaa);

        double expectedSum = result.ExpectedAA + result.ExpectedAa + result.Expectedaa;
        Assert.That(expectedSum, Is.EqualTo(n).Within(Tolerance),
            "Expected counts must sum to total sample size");
    }

    #endregion

    #region Significance Level Tests

    /// <summary>
    /// HW-M16: Custom significance level is respected.
    /// Uses borderline data (46,49,5) where χ²≈3.17, p-value≈0.075.
    /// At α=0.01: p-value (0.075) ≥ 0.01 → InEquilibrium = true
    /// At α=0.10: p-value (0.075) < 0.10 → InEquilibrium = false
    /// </summary>
    [Test]
    public void TestHardyWeinberg_CustomSignificanceLevel()
    {
        var resultAt01 = PopulationGeneticsAnalyzer.TestHardyWeinberg(
            "SIG_TEST", 46, 49, 5, significanceLevel: 0.01);
        var resultAt10 = PopulationGeneticsAnalyzer.TestHardyWeinberg(
            "SIG_TEST", 46, 49, 5, significanceLevel: 0.10);

        Assert.Multiple(() =>
        {
            // P-value should be identical for same data regardless of significance level
            Assert.That(resultAt01.PValue, Is.EqualTo(resultAt10.PValue).Within(Tolerance),
                "P-value should be same regardless of significance level");

            // P-value should be between 0.01 and 0.10 (verified from χ²≈3.17)
            Assert.That(resultAt01.PValue, Is.GreaterThan(0.01).And.LessThan(0.10),
                "P-value must fall between the two significance levels for this test to be meaningful");

            // At α=0.01: fail to reject H₀
            Assert.That(resultAt01.InEquilibrium, Is.True,
                "At α=0.01, borderline case should be in equilibrium");

            // At α=0.10: reject H₀
            Assert.That(resultAt10.InEquilibrium, Is.False,
                "At α=0.10, borderline case should NOT be in equilibrium");
        });
    }

    /// <summary>
    /// HW-M17: Default significance level is 0.05.
    /// Uses borderline data (46,49,5) where p-value≈0.075, which is >0.05.
    /// Both default and explicit α=0.05 should yield InEquilibrium=true.
    /// </summary>
    [Test]
    public void TestHardyWeinberg_DefaultSignificanceIs005()
    {
        // p-value ≈ 0.075 for this data → InEquilibrium=true at α=0.05
        var resultDefault = PopulationGeneticsAnalyzer.TestHardyWeinberg(
            "DEFAULT_SIG", 46, 49, 5);
        var resultExplicit = PopulationGeneticsAnalyzer.TestHardyWeinberg(
            "EXPLICIT_SIG", 46, 49, 5, significanceLevel: 0.05);

        Assert.Multiple(() =>
        {
            Assert.That(resultDefault.InEquilibrium, Is.EqualTo(resultExplicit.InEquilibrium),
                "Default significance should be 0.05");
            Assert.That(resultDefault.ChiSquare, Is.EqualTo(resultExplicit.ChiSquare).Within(Tolerance),
                "Chi-square should be identical");
            Assert.That(resultDefault.PValue, Is.EqualTo(resultExplicit.PValue).Within(Tolerance),
                "P-value should be identical");

            // Verify this is a non-trivial case: p-value is between 0.05 and 0.10
            Assert.That(resultDefault.PValue, Is.GreaterThan(0.05).And.LessThan(0.10),
                "Test data should produce a p-value between 0.05 and 0.10 to be meaningful");
            Assert.That(resultDefault.InEquilibrium, Is.True,
                "With default α=0.05 and p-value>0.05, should be in equilibrium");
        });
    }

    #endregion

    #region SHOULD Tests

    /// <summary>
    /// HW-S01: Large sample size for numerical stability.
    /// </summary>
    [Test]
    public void TestHardyWeinberg_LargeSample_NumericallyStable()
    {
        // n = 10000, perfect HWE with p = 0.5
        var result = PopulationGeneticsAnalyzer.TestHardyWeinberg(
            variantId: "LARGE_SAMPLE",
            observedAA: 2500,
            observedAa: 5000,
            observedaa: 2500);

        Assert.Multiple(() =>
        {
            Assert.That(result.ChiSquare, Is.EqualTo(0).Within(Tolerance));
            Assert.That(result.InEquilibrium, Is.True);
            Assert.That(double.IsFinite(result.PValue), Is.True,
                "P-value should be finite for large samples");
        });
    }

    /// <summary>
    /// HW-S02: Various allele frequencies — data chosen such that observed = expected.
    /// All are perfect HWE: genotype counts exactly equal p², 2pq, q² × n.
    /// </summary>
    [TestCase(81, 18, 1, 0.9)]   // p = 0.9: E=(81,18,1) = Observed
    [TestCase(49, 42, 9, 0.7)]   // p = 0.7: E=(49,42,9) = Observed
    [TestCase(9, 42, 49, 0.3)]   // p = 0.3: E=(9,42,49) = Observed
    [TestCase(1, 18, 81, 0.1)]   // p = 0.1: E=(1,18,81) = Observed
    public void TestHardyWeinberg_VariousAlleleFrequencies_CorrectExpected(
        int observedAA, int observedAa, int observedaa, double expectedP)
    {
        int n = observedAA + observedAa + observedaa;
        var result = PopulationGeneticsAnalyzer.TestHardyWeinberg("FREQ_TEST", observedAA, observedAa, observedaa);

        double q = 1 - expectedP;
        Assert.Multiple(() =>
        {
            Assert.That(result.ExpectedAA, Is.EqualTo(expectedP * expectedP * n).Within(Tolerance),
                $"E(AA) should equal {expectedP}² × {n}");
            Assert.That(result.ExpectedAa, Is.EqualTo(2 * expectedP * q * n).Within(Tolerance),
                $"E(Aa) should equal 2×{expectedP}×{q}×{n}");
            Assert.That(result.Expectedaa, Is.EqualTo(q * q * n).Within(Tolerance),
                $"E(aa) should equal {q}² × {n}");

            // These are perfect HWE cases (observed = expected), so χ² should be 0
            Assert.That(result.ChiSquare, Is.EqualTo(0).Within(Tolerance),
                "Chi-square should be 0 when observed counts exactly match HWE expectations");
            Assert.That(result.InEquilibrium, Is.True,
                "Perfect HWE data should report equilibrium");
        });
    }

    /// <summary>
    /// HW-S03: Borderline significance - decision boundary testing.
    /// Data: AA=46, Aa=49, aa=5, n=100, p=0.705
    /// χ² ≈ 3.17, p-value ≈ 0.075
    /// At α=0.05: p-value (0.075) ≥ 0.05 → InEquilibrium = true
    /// At α=0.10: p-value (0.075) &lt; 0.10 → InEquilibrium = false
    /// Source: Computed from HWE formulas per Wikipedia.
    /// </summary>
    [Test]
    public void TestHardyWeinberg_BorderlineSignificance_DecisionDependsOnAlpha()
    {
        // χ² ≈ 3.17, p-value ≈ 0.075 — between α=0.05 and α=0.10
        var resultAt05 = PopulationGeneticsAnalyzer.TestHardyWeinberg(
            variantId: "BORDERLINE",
            observedAA: 46,
            observedAa: 49,
            observedaa: 5,
            significanceLevel: 0.05);

        var resultAt10 = PopulationGeneticsAnalyzer.TestHardyWeinberg(
            variantId: "BORDERLINE",
            observedAA: 46,
            observedAa: 49,
            observedaa: 5,
            significanceLevel: 0.10);

        Assert.Multiple(() =>
        {
            // Chi-square should be between critical values for α=0.10 (2.706) and α=0.05 (3.841)
            Assert.That(resultAt05.ChiSquare, Is.GreaterThan(2.706).And.LessThan(3.841),
                "Chi-square should be between α=0.10 and α=0.05 critical values");

            // Same p-value regardless of significance level
            Assert.That(resultAt05.PValue, Is.EqualTo(resultAt10.PValue).Within(Tolerance),
                "P-value should be identical for same data");

            // P-value should be between 0.05 and 0.10
            Assert.That(resultAt05.PValue, Is.GreaterThan(0.05).And.LessThan(0.10),
                "P-value should fall between the two significance levels");

            // At α=0.05: fail to reject H₀ (in equilibrium)
            Assert.That(resultAt05.InEquilibrium, Is.True,
                "At α=0.05, borderline case should be in equilibrium");

            // At α=0.10: reject H₀ (not in equilibrium)
            Assert.That(resultAt10.InEquilibrium, Is.False,
                "At α=0.10, borderline case should NOT be in equilibrium");
        });
    }

    #endregion

    #region COULD Tests

    /// <summary>
    /// HW-C01: Symmetric genotypes produce the same χ² and p-value.
    /// Swapping AA↔aa relabels alleles (A↔a) but does not change the HWE test.
    /// This follows from the symmetry of p²+2pq+q² under p↔q.
    /// </summary>
    [TestCase(70, 25, 5)]
    [TestCase(10, 80, 10)]
    [TestCase(40, 45, 15)]
    public void TestHardyWeinberg_SymmetricGenotypes_SameChiSquareAndPValue(
        int observedAA, int observedAa, int observedaa)
    {
        var original = PopulationGeneticsAnalyzer.TestHardyWeinberg(
            "SYM_ORIG", observedAA, observedAa, observedaa);
        var swapped = PopulationGeneticsAnalyzer.TestHardyWeinberg(
            "SYM_SWAP", observedaa, observedAa, observedAA);

        Assert.Multiple(() =>
        {
            Assert.That(swapped.ChiSquare, Is.EqualTo(original.ChiSquare).Within(Tolerance),
                "Chi-square should be identical when AA and aa are swapped");
            Assert.That(swapped.PValue, Is.EqualTo(original.PValue).Within(Tolerance),
                "P-value should be identical when AA and aa are swapped");
            Assert.That(swapped.InEquilibrium, Is.EqualTo(original.InEquilibrium),
                "Equilibrium decision should be identical when AA and aa are swapped");

            // Expected counts should be swapped: E(AA)↔E(aa), E(Aa) unchanged
            Assert.That(swapped.ExpectedAA, Is.EqualTo(original.Expectedaa).Within(Tolerance),
                "Swapped E(AA) should equal original E(aa)");
            Assert.That(swapped.Expectedaa, Is.EqualTo(original.ExpectedAA).Within(Tolerance),
                "Swapped E(aa) should equal original E(AA)");
            Assert.That(swapped.ExpectedAa, Is.EqualTo(original.ExpectedAa).Within(Tolerance),
                "E(Aa) should be unchanged under allele relabeling");
        });
    }

    /// <summary>
    /// HW-C02: Extreme frequency (p=0.99) with large n.
    /// Rare variants should be handled correctly without numerical issues.
    /// n=10000, p=0.99: E(AA)=9801, E(Aa)=198, E(aa)=1.
    /// Observed = Expected → χ²=0.
    /// </summary>
    [Test]
    public void TestHardyWeinberg_ExtremeFrequency_RareVariantStable()
    {
        // p = (2*9801 + 198) / 20000 = 19800/20000 = 0.99
        // E(AA) = 0.99² × 10000 = 9801
        // E(Aa) = 2×0.99×0.01 × 10000 = 198
        // E(aa) = 0.01² × 10000 = 1
        var result = PopulationGeneticsAnalyzer.TestHardyWeinberg(
            variantId: "RARE_VARIANT",
            observedAA: 9801,
            observedAa: 198,
            observedaa: 1);

        Assert.Multiple(() =>
        {
            Assert.That(result.ExpectedAA, Is.EqualTo(9801.0).Within(Tolerance),
                "E(AA) = 0.99² × 10000 = 9801");
            Assert.That(result.ExpectedAa, Is.EqualTo(198.0).Within(Tolerance),
                "E(Aa) = 2×0.99×0.01 × 10000 = 198");
            Assert.That(result.Expectedaa, Is.EqualTo(1.0).Within(Tolerance),
                "E(aa) = 0.01² × 10000 = 1");
            Assert.That(result.ChiSquare, Is.EqualTo(0).Within(Tolerance),
                "Chi-square should be 0 for perfect HWE");
            Assert.That(result.InEquilibrium, Is.True,
                "Perfect HWE should report equilibrium");
            Assert.That(double.IsFinite(result.PValue), Is.True,
                "P-value should be finite even with extreme allele frequency");
        });
    }

    #endregion
}
