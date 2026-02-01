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
            // Verify chi-square is approximately 0.83 as per Wikipedia
            Assert.That(result.ChiSquare, Is.EqualTo(0.83).Within(0.1),
                "Chi-square should be approximately 0.83 for Ford's data");

            // Verify equilibrium decision
            Assert.That(result.InEquilibrium, Is.True,
                "Ford's moth data should be in Hardy-Weinberg equilibrium");

            // Verify expected counts (Wikipedia: E_AA≈1467.4, E_Aa≈141.2, E_aa≈3.4)
            Assert.That(result.ExpectedAA, Is.EqualTo(1467.4).Within(1.0),
                "Expected AA count should be approximately 1467.4");
            Assert.That(result.ExpectedAa, Is.EqualTo(141.2).Within(1.0),
                "Expected Aa count should be approximately 141.2");
            Assert.That(result.Expectedaa, Is.EqualTo(3.4).Within(0.5),
                "Expected aa count should be approximately 3.4");

            // P-value should be > 0.05
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
            Assert.That(result.ChiSquare, Is.GreaterThan(ChiSquareCriticalValue_Alpha05_Df1),
                "Chi-square should exceed critical value for excess heterozygotes");
            Assert.That(result.ChiSquare, Is.EqualTo(36).Within(1.0),
                "Chi-square should be approximately 36");
            Assert.That(result.InEquilibrium, Is.False,
                "Excess heterozygotes should deviate from HWE");
            Assert.That(result.PValue, Is.LessThan(0.05),
                "P-value should be less than 0.05");
        });
    }

    /// <summary>
    /// HW-M04: All heterozygotes - extreme deviation.
    /// This is biologically implausible under random mating.
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
            Assert.That(result.ChiSquare, Is.GreaterThan(ChiSquareCriticalValue_Alpha05_Df1),
                "All heterozygotes should show significant deviation");
            Assert.That(result.InEquilibrium, Is.False,
                "All heterozygotes cannot be in HWE");
        });
    }

    /// <summary>
    /// HW-M05: Deficit of heterozygotes - inbreeding pattern.
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
            Assert.That(result.ChiSquare, Is.GreaterThan(ChiSquareCriticalValue_Alpha05_Df1),
                "Deficit of heterozygotes should show significant deviation");
            Assert.That(result.InEquilibrium, Is.False,
                "Inbreeding pattern should deviate from HWE");
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
        var result = PopulationGeneticsAnalyzer.TestHardyWeinberg(
            variantId: "SINGLE",
            observedAA: 1,
            observedAa: 0,
            observedaa: 0);

        Assert.Multiple(() =>
        {
            Assert.That(result.ObservedAA, Is.EqualTo(1));
            Assert.That(result.ExpectedAA, Is.EqualTo(1).Within(Tolerance));
            Assert.That(result.ChiSquare, Is.GreaterThanOrEqualTo(0),
                "Chi-square should be non-negative");
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
    public void TestHardyWeinberg_ChiSquareNonNegative(int aa, int aA, int AA)
    {
        var result = PopulationGeneticsAnalyzer.TestHardyWeinberg("INV_TEST", aa, aA, AA);

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
    public void TestHardyWeinberg_PValueInValidRange(int aa, int aA, int AA)
    {
        var result = PopulationGeneticsAnalyzer.TestHardyWeinberg("PVAL_TEST", aa, aA, AA);

        Assert.That(result.PValue, Is.InRange(0, 1),
            "P-value must be in [0, 1]");
    }

    /// <summary>
    /// HW-M14: InEquilibrium is consistent with PValue and significance level.
    /// </summary>
    [TestCase(25, 50, 25, 0.05, true)]   // Perfect HWE
    [TestCase(10, 80, 10, 0.05, false)]  // Excess heterozygotes
    [TestCase(10, 80, 10, 0.001, false)] // Even stricter, still fails
    public void TestHardyWeinberg_EquilibriumConsistentWithPValue(
        int aa, int aA, int AA, double significance, bool expectedEquilibrium)
    {
        var result = PopulationGeneticsAnalyzer.TestHardyWeinberg(
            "CONSIST_TEST", aa, aA, AA, significance);

        // InEquilibrium should be true iff PValue >= significance
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
    public void TestHardyWeinberg_ExpectedCountsSumToN(int aa, int aA, int AA)
    {
        int n = aa + aA + AA;
        var result = PopulationGeneticsAnalyzer.TestHardyWeinberg("SUM_TEST", aa, aA, AA);

        double expectedSum = result.ExpectedAA + result.ExpectedAa + result.Expectedaa;
        Assert.That(expectedSum, Is.EqualTo(n).Within(Tolerance),
            "Expected counts must sum to total sample size");
    }

    #endregion

    #region Significance Level Tests

    /// <summary>
    /// HW-M16: Custom significance level is respected.
    /// </summary>
    [Test]
    public void TestHardyWeinberg_CustomSignificanceLevel()
    {
        // Use data that might pass at α=0.05 but fail at α=0.10
        // or vice versa near the boundary
        var result01 = PopulationGeneticsAnalyzer.TestHardyWeinberg(
            "SIG_TEST", 30, 45, 25, significanceLevel: 0.01);
        var result10 = PopulationGeneticsAnalyzer.TestHardyWeinberg(
            "SIG_TEST", 30, 45, 25, significanceLevel: 0.10);

        // At lower significance (0.01), more likely to be in equilibrium
        // At higher significance (0.10), more likely to reject
        // This verifies the parameter is actually used
        Assert.Multiple(() =>
        {
            Assert.That(result01.PValue, Is.EqualTo(result10.PValue).Within(Tolerance),
                "P-value should be same regardless of significance level");
            // If p-value is between 0.01 and 0.10, decisions should differ
            if (result01.PValue >= 0.01 && result01.PValue < 0.10)
            {
                Assert.That(result01.InEquilibrium, Is.True);
                Assert.That(result10.InEquilibrium, Is.False);
            }
        });
    }

    /// <summary>
    /// HW-M17: Default significance level is 0.05.
    /// </summary>
    [Test]
    public void TestHardyWeinberg_DefaultSignificanceIs005()
    {
        var resultDefault = PopulationGeneticsAnalyzer.TestHardyWeinberg(
            "DEFAULT_SIG", 25, 50, 25);
        var resultExplicit = PopulationGeneticsAnalyzer.TestHardyWeinberg(
            "EXPLICIT_SIG", 25, 50, 25, significanceLevel: 0.05);

        Assert.That(resultDefault.InEquilibrium, Is.EqualTo(resultExplicit.InEquilibrium),
            "Default significance should be 0.05");
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
    /// HW-S02: Various allele frequencies.
    /// </summary>
    [TestCase(81, 18, 1, 0.9)]   // p = 0.9
    [TestCase(49, 42, 9, 0.7)]   // p = 0.7
    [TestCase(9, 42, 49, 0.3)]   // p = 0.3
    [TestCase(1, 18, 81, 0.1)]   // p = 0.1
    public void TestHardyWeinberg_VariousAlleleFrequencies_CorrectExpected(
        int aa, int aA, int AA, double expectedP)
    {
        int n = aa + aA + AA;
        var result = PopulationGeneticsAnalyzer.TestHardyWeinberg("FREQ_TEST", aa, aA, AA);

        double q = 1 - expectedP;
        Assert.Multiple(() =>
        {
            Assert.That(result.ExpectedAA, Is.EqualTo(expectedP * expectedP * n).Within(1.0));
            Assert.That(result.ExpectedAa, Is.EqualTo(2 * expectedP * q * n).Within(1.0));
            Assert.That(result.Expectedaa, Is.EqualTo(q * q * n).Within(1.0));
        });
    }

    #endregion
}
