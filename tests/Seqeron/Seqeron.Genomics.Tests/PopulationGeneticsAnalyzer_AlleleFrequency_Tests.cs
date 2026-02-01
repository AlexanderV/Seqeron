using NUnit.Framework;
using Seqeron.Genomics;
using System.Collections.Generic;
using System.Linq;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Canonical tests for POP-FREQ-001: Allele Frequency Calculation.
/// 
/// Methods tested:
/// - CalculateAlleleFrequencies(int, int, int)
/// - CalculateMAF(IEnumerable&lt;int&gt;)
/// - FilterByMAF(IEnumerable&lt;Variant&gt;, double, double)
/// 
/// Evidence: Wikipedia (Allele frequency, Minor allele frequency, Genotype frequency),
///           Gillespie (2004) Population Genetics: A Concise Guide
/// </summary>
[TestFixture]
public class PopulationGeneticsAnalyzer_AlleleFrequency_Tests
{
    #region CalculateAlleleFrequencies Tests

    /// <summary>
    /// AF-M01: Allele frequencies must sum to 1.0 for non-empty populations.
    /// Source: Wikipedia - "p and q must sum to 1"
    /// </summary>
    [Test]
    [TestCase(49, 42, 9)]      // Wikipedia flower example
    [TestCase(25, 50, 25)]     // HWE balanced
    [TestCase(100, 0, 0)]      // All homozygous major
    [TestCase(0, 0, 100)]      // All homozygous minor
    [TestCase(0, 100, 0)]      // All heterozygous
    [TestCase(1, 0, 0)]        // Single sample
    [TestCase(10000, 5000, 2500)] // Large population
    public void CalculateAlleleFrequencies_NonEmptyPopulation_FrequenciesSumToOne(
        int homMaj, int het, int homMin)
    {
        var (major, minor) = PopulationGeneticsAnalyzer.CalculateAlleleFrequencies(
            homMaj, het, homMin);

        Assert.That(major + minor, Is.EqualTo(1.0).Within(1e-10),
            "Allele frequencies must sum to 1.0 (p + q = 1)");
    }

    /// <summary>
    /// AF-M02: Wikipedia flower example - canonical test case.
    /// Source: Wikipedia Genotype Frequency - "49 AA, 42 Aa, 9 aa → p=0.70, q=0.30"
    /// </summary>
    [Test]
    public void CalculateAlleleFrequencies_WikipediaFlowerExample_ReturnsCorrectFrequencies()
    {
        // Wikipedia example: 100 four-o'-clock plants
        // 49 red (AA), 42 pink (Aa), 9 white (aa)
        // Expected: p(A) = (2×49 + 42) / 200 = 140/200 = 0.70
        //           q(a) = (2×9 + 42) / 200 = 60/200 = 0.30
        var (major, minor) = PopulationGeneticsAnalyzer.CalculateAlleleFrequencies(
            homozygousMajor: 49,
            heterozygous: 42,
            homozygousMinor: 9);

        Assert.Multiple(() =>
        {
            Assert.That(major, Is.EqualTo(0.70).Within(0.001),
                "Major allele frequency should be 0.70");
            Assert.That(minor, Is.EqualTo(0.30).Within(0.001),
                "Minor allele frequency should be 0.30");
        });
    }

    /// <summary>
    /// AF-M03: All homozygous major yields (1.0, 0.0).
    /// Source: Derived from allele counting formula
    /// </summary>
    [Test]
    public void CalculateAlleleFrequencies_AllHomozygousMajor_ReturnsMajorFrequencyOne()
    {
        var (major, minor) = PopulationGeneticsAnalyzer.CalculateAlleleFrequencies(
            homozygousMajor: 100,
            heterozygous: 0,
            homozygousMinor: 0);

        Assert.Multiple(() =>
        {
            Assert.That(major, Is.EqualTo(1.0).Within(1e-10));
            Assert.That(minor, Is.EqualTo(0.0).Within(1e-10));
        });
    }

    /// <summary>
    /// AF-M04: All homozygous minor yields (0.0, 1.0).
    /// Source: Derived from allele counting formula
    /// </summary>
    [Test]
    public void CalculateAlleleFrequencies_AllHomozygousMinor_ReturnsMinorFrequencyOne()
    {
        var (major, minor) = PopulationGeneticsAnalyzer.CalculateAlleleFrequencies(
            homozygousMajor: 0,
            heterozygous: 0,
            homozygousMinor: 100);

        Assert.Multiple(() =>
        {
            Assert.That(major, Is.EqualTo(0.0).Within(1e-10));
            Assert.That(minor, Is.EqualTo(1.0).Within(1e-10));
        });
    }

    /// <summary>
    /// AF-M05: All heterozygous yields (0.5, 0.5).
    /// Source: Derived from allele counting formula (each het has one of each allele)
    /// </summary>
    [Test]
    public void CalculateAlleleFrequencies_AllHeterozygous_ReturnsEqualFrequencies()
    {
        var (major, minor) = PopulationGeneticsAnalyzer.CalculateAlleleFrequencies(
            homozygousMajor: 0,
            heterozygous: 100,
            homozygousMinor: 0);

        Assert.Multiple(() =>
        {
            Assert.That(major, Is.EqualTo(0.5).Within(1e-10));
            Assert.That(minor, Is.EqualTo(0.5).Within(1e-10));
        });
    }

    /// <summary>
    /// AF-M06: HWE balanced genotypes (25-50-25) yield (0.5, 0.5).
    /// Source: Hardy-Weinberg equilibrium expectation for p=q=0.5
    /// </summary>
    [Test]
    public void CalculateAlleleFrequencies_HweBalanced_ReturnsEqualFrequencies()
    {
        var (major, minor) = PopulationGeneticsAnalyzer.CalculateAlleleFrequencies(
            homozygousMajor: 25,
            heterozygous: 50,
            homozygousMinor: 25);

        Assert.Multiple(() =>
        {
            Assert.That(major, Is.EqualTo(0.5).Within(1e-10));
            Assert.That(minor, Is.EqualTo(0.5).Within(1e-10));
        });
    }

    /// <summary>
    /// AF-M07: Zero samples returns (0, 0) gracefully.
    /// Source: Edge case handling
    /// </summary>
    [Test]
    public void CalculateAlleleFrequencies_ZeroSamples_ReturnsZeroZero()
    {
        var (major, minor) = PopulationGeneticsAnalyzer.CalculateAlleleFrequencies(0, 0, 0);

        Assert.Multiple(() =>
        {
            Assert.That(major, Is.EqualTo(0));
            Assert.That(minor, Is.EqualTo(0));
        });
    }

    /// <summary>
    /// AF-M08/M09: Frequencies are in valid range [0, 1].
    /// Source: Range constraint from allele frequency definition
    /// </summary>
    [Test]
    [TestCase(49, 42, 9)]
    [TestCase(100, 0, 0)]
    [TestCase(0, 0, 100)]
    [TestCase(1, 1, 1)]
    [TestCase(7, 21, 72)]
    public void CalculateAlleleFrequencies_AnyValidInput_FrequenciesInValidRange(
        int homMaj, int het, int homMin)
    {
        var (major, minor) = PopulationGeneticsAnalyzer.CalculateAlleleFrequencies(
            homMaj, het, homMin);

        Assert.Multiple(() =>
        {
            Assert.That(major, Is.GreaterThanOrEqualTo(0), "Major frequency must be ≥ 0");
            Assert.That(major, Is.LessThanOrEqualTo(1), "Major frequency must be ≤ 1");
            Assert.That(minor, Is.GreaterThanOrEqualTo(0), "Minor frequency must be ≥ 0");
            Assert.That(minor, Is.LessThanOrEqualTo(1), "Minor frequency must be ≤ 1");
        });
    }

    /// <summary>
    /// AF-S01: Large population does not overflow.
    /// Source: Implementation robustness
    /// </summary>
    [Test]
    public void CalculateAlleleFrequencies_LargePopulation_NoOverflow()
    {
        var (major, minor) = PopulationGeneticsAnalyzer.CalculateAlleleFrequencies(
            homozygousMajor: 100_000,
            heterozygous: 50_000,
            homozygousMinor: 25_000);

        Assert.Multiple(() =>
        {
            Assert.That(major + minor, Is.EqualTo(1.0).Within(1e-10));
            Assert.That(double.IsNaN(major), Is.False);
            Assert.That(double.IsNaN(minor), Is.False);
        });
    }

    #endregion

    #region CalculateMAF Tests

    /// <summary>
    /// MAF-M01: Alt frequency &lt; 0.5 returns alt frequency as MAF.
    /// Source: MAF definition - "frequency of second most common allele"
    /// </summary>
    [Test]
    public void CalculateMAF_AltFrequencyBelowHalf_ReturnsAltFrequency()
    {
        // 10 samples: 5 hom ref (0), 4 het (1), 1 hom alt (2)
        // Alt alleles = 0×5 + 1×4 + 2×1 = 6
        // Total = 20, alt freq = 0.3, MAF = 0.3
        var genotypes = new[] { 0, 0, 0, 0, 0, 1, 1, 1, 1, 2 };

        double maf = PopulationGeneticsAnalyzer.CalculateMAF(genotypes);

        Assert.That(maf, Is.EqualTo(0.3).Within(0.001));
    }

    /// <summary>
    /// MAF-M02: Alt frequency > 0.5 returns 1 - alt frequency as MAF.
    /// Source: MAF = min(f, 1-f)
    /// </summary>
    [Test]
    public void CalculateMAF_AltFrequencyAboveHalf_ReturnsOneMinusAltFrequency()
    {
        // 10 samples: 1 hom ref (0), 4 het (1), 5 hom alt (2)
        // Alt alleles = 0×1 + 1×4 + 2×5 = 14
        // Total = 20, alt freq = 0.7, MAF = 0.3
        var genotypes = new[] { 0, 1, 1, 1, 1, 2, 2, 2, 2, 2 };

        double maf = PopulationGeneticsAnalyzer.CalculateMAF(genotypes);

        Assert.That(maf, Is.EqualTo(0.3).Within(0.001));
    }

    /// <summary>
    /// MAF-M03/M04: MAF is always in range [0, 0.5].
    /// Source: MAF invariant
    /// </summary>
    [Test]
    [TestCase(new[] { 0, 0, 0, 0, 0 })]           // All ref
    [TestCase(new[] { 2, 2, 2, 2, 2 })]           // All alt
    [TestCase(new[] { 1, 1, 1, 1, 1 })]           // All het
    [TestCase(new[] { 0, 1, 2 })]                 // Mixed
    [TestCase(new[] { 0, 0, 0, 0, 1, 1, 2 })]     // Low alt freq
    [TestCase(new[] { 2, 2, 2, 2, 1, 1, 0 })]     // High alt freq
    public void CalculateMAF_AnyGenotypes_MafInValidRange(int[] genotypes)
    {
        double maf = PopulationGeneticsAnalyzer.CalculateMAF(genotypes);

        Assert.Multiple(() =>
        {
            Assert.That(maf, Is.GreaterThanOrEqualTo(0), "MAF must be ≥ 0");
            Assert.That(maf, Is.LessThanOrEqualTo(0.5), "MAF must be ≤ 0.5");
        });
    }

    /// <summary>
    /// MAF-M05: Monomorphic reference (all 0) yields MAF = 0.
    /// Source: Fixed allele has no minor allele
    /// </summary>
    [Test]
    public void CalculateMAF_MonomorphicReference_ReturnsZero()
    {
        var genotypes = new[] { 0, 0, 0, 0, 0 };

        double maf = PopulationGeneticsAnalyzer.CalculateMAF(genotypes);

        Assert.That(maf, Is.EqualTo(0));
    }

    /// <summary>
    /// MAF-M06: Monomorphic alternate (all 2) yields MAF = 0.
    /// Source: Fixed allele has no minor allele (ref becomes minor at freq 0)
    /// </summary>
    [Test]
    public void CalculateMAF_MonomorphicAlternate_ReturnsZero()
    {
        var genotypes = new[] { 2, 2, 2, 2, 2 };

        double maf = PopulationGeneticsAnalyzer.CalculateMAF(genotypes);

        Assert.That(maf, Is.EqualTo(0));
    }

    /// <summary>
    /// MAF-M07: Empty genotypes returns 0.
    /// Source: Edge case handling
    /// </summary>
    [Test]
    public void CalculateMAF_EmptyGenotypes_ReturnsZero()
    {
        var genotypes = System.Array.Empty<int>();

        double maf = PopulationGeneticsAnalyzer.CalculateMAF(genotypes);

        Assert.That(maf, Is.EqualTo(0));
    }

    /// <summary>
    /// MAF-M08: Perfect 50/50 split yields MAF = 0.5.
    /// Source: Maximum MAF value at balanced polymorphism
    /// </summary>
    [Test]
    public void CalculateMAF_PerfectBalance_ReturnsHalf()
    {
        // 4 samples: 2 het gives alt_freq = 0.5
        var genotypes = new[] { 0, 1, 1, 2 };
        // Alt alleles = 0 + 1 + 1 + 2 = 4, total = 8, freq = 0.5

        double maf = PopulationGeneticsAnalyzer.CalculateMAF(genotypes);

        Assert.That(maf, Is.EqualTo(0.5).Within(0.001));
    }

    /// <summary>
    /// MAF-S01: Various genotype patterns compute correctly.
    /// Source: Algorithm verification
    /// </summary>
    [Test]
    [TestCase(new[] { 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 }, 0.05)]  // 1/20 = 0.05
    [TestCase(new[] { 0, 0, 0, 0, 0, 0, 0, 0, 1, 1 }, 0.10)]  // 2/20 = 0.10
    [TestCase(new[] { 0, 0, 0, 0, 0, 0, 0, 0, 1, 2 }, 0.15)]  // 3/20 = 0.15
    public void CalculateMAF_VariousPatterns_CalculatesCorrectly(int[] genotypes, double expectedMaf)
    {
        double maf = PopulationGeneticsAnalyzer.CalculateMAF(genotypes);

        Assert.That(maf, Is.EqualTo(expectedMaf).Within(0.001));
    }

    #endregion

    #region FilterByMAF Tests

    /// <summary>
    /// FLT-M01: Variants below minMAF are filtered out.
    /// Source: Filter logic
    /// </summary>
    [Test]
    public void FilterByMAF_VariantsBelowMinMaf_AreExcluded()
    {
        var variants = new List<PopulationGeneticsAnalyzer.Variant>
        {
            new("V1", "chr1", 100, "A", "G", 0.005, 100), // MAF = 0.005 < 0.01
            new("V2", "chr1", 200, "A", "G", 0.05, 100),  // MAF = 0.05 ≥ 0.01
        };

        var filtered = PopulationGeneticsAnalyzer.FilterByMAF(variants, minMAF: 0.01).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(filtered, Has.Count.EqualTo(1));
            Assert.That(filtered[0].Id, Is.EqualTo("V2"));
        });
    }

    /// <summary>
    /// FLT-M02: Variants at/above minMAF are included.
    /// Source: Filter logic - boundary case
    /// </summary>
    [Test]
    public void FilterByMAF_VariantsAtMinMaf_AreIncluded()
    {
        var variants = new List<PopulationGeneticsAnalyzer.Variant>
        {
            new("V1", "chr1", 100, "A", "G", 0.01, 100), // MAF = 0.01 = minMAF
        };

        var filtered = PopulationGeneticsAnalyzer.FilterByMAF(variants, minMAF: 0.01).ToList();

        Assert.That(filtered, Has.Count.EqualTo(1));
    }

    /// <summary>
    /// FLT-M03: Variants above maxMAF are filtered out.
    /// Source: Filter logic - high MAF filtering (for near-monomorphic)
    /// </summary>
    [Test]
    public void FilterByMAF_VariantsAboveMaxMaf_AreExcluded()
    {
        var variants = new List<PopulationGeneticsAnalyzer.Variant>
        {
            new("V1", "chr1", 100, "A", "G", 0.3, 100),  // MAF = 0.3 < 0.4
            new("V2", "chr1", 200, "A", "G", 0.45, 100), // MAF = 0.45 > 0.4
        };

        var filtered = PopulationGeneticsAnalyzer.FilterByMAF(variants, minMAF: 0.01, maxMAF: 0.4).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(filtered, Has.Count.EqualTo(1));
            Assert.That(filtered[0].Id, Is.EqualTo("V1"));
        });
    }

    /// <summary>
    /// FLT-M04: Variants at/below maxMAF are included.
    /// Source: Filter logic - boundary case
    /// </summary>
    [Test]
    public void FilterByMAF_VariantsAtMaxMaf_AreIncluded()
    {
        var variants = new List<PopulationGeneticsAnalyzer.Variant>
        {
            new("V1", "chr1", 100, "A", "G", 0.4, 100), // MAF = 0.4 = maxMAF
        };

        var filtered = PopulationGeneticsAnalyzer.FilterByMAF(variants, minMAF: 0.01, maxMAF: 0.4).ToList();

        Assert.That(filtered, Has.Count.EqualTo(1));
    }

    /// <summary>
    /// FLT-M05: Empty input yields empty output.
    /// Source: Edge case handling
    /// </summary>
    [Test]
    public void FilterByMAF_EmptyInput_ReturnsEmptyEnumerable()
    {
        var variants = Enumerable.Empty<PopulationGeneticsAnalyzer.Variant>();

        var filtered = PopulationGeneticsAnalyzer.FilterByMAF(variants, minMAF: 0.01).ToList();

        Assert.That(filtered, Is.Empty);
    }

    /// <summary>
    /// FLT-M06: All variants filtered yields empty output.
    /// Source: Edge case handling
    /// </summary>
    [Test]
    public void FilterByMAF_AllVariantsFiltered_ReturnsEmptyEnumerable()
    {
        var variants = new List<PopulationGeneticsAnalyzer.Variant>
        {
            new("V1", "chr1", 100, "A", "G", 0.001, 100), // MAF < 0.01
            new("V2", "chr1", 200, "A", "G", 0.002, 100), // MAF < 0.01
        };

        var filtered = PopulationGeneticsAnalyzer.FilterByMAF(variants, minMAF: 0.01).ToList();

        Assert.That(filtered, Is.Empty);
    }

    /// <summary>
    /// FLT-M07: All variants pass when all are in range.
    /// Source: Edge case handling
    /// </summary>
    [Test]
    public void FilterByMAF_AllVariantsPass_ReturnsAllVariants()
    {
        var variants = new List<PopulationGeneticsAnalyzer.Variant>
        {
            new("V1", "chr1", 100, "A", "G", 0.1, 100),
            new("V2", "chr1", 200, "A", "G", 0.2, 100),
            new("V3", "chr1", 300, "A", "G", 0.3, 100),
        };

        var filtered = PopulationGeneticsAnalyzer.FilterByMAF(variants, minMAF: 0.05, maxMAF: 0.4).ToList();

        Assert.That(filtered, Has.Count.EqualTo(3));
    }

    /// <summary>
    /// FLT-S01: Filter handles high allele frequencies correctly (MAF calculation).
    /// Source: MAF is min(AF, 1-AF)
    /// </summary>
    [Test]
    public void FilterByMAF_HighAlleleFrequency_CalculatesMafCorrectly()
    {
        // AF = 0.95 → MAF = 0.05
        var variants = new List<PopulationGeneticsAnalyzer.Variant>
        {
            new("V1", "chr1", 100, "A", "G", 0.95, 100), // MAF = 0.05
        };

        var filtered = PopulationGeneticsAnalyzer.FilterByMAF(variants, minMAF: 0.01, maxMAF: 0.1).ToList();

        Assert.That(filtered, Has.Count.EqualTo(1),
            "Variant with AF=0.95 should have MAF=0.05 and pass filter");
    }

    /// <summary>
    /// FLT-S02: Filter preserves order of input variants.
    /// Source: Implementation expectation
    /// </summary>
    [Test]
    public void FilterByMAF_MultipleVariants_PreservesOrder()
    {
        var variants = new List<PopulationGeneticsAnalyzer.Variant>
        {
            new("V1", "chr1", 100, "A", "G", 0.1, 100),
            new("V2", "chr1", 200, "A", "G", 0.2, 100),
            new("V3", "chr1", 300, "A", "G", 0.3, 100),
        };

        var filtered = PopulationGeneticsAnalyzer.FilterByMAF(variants, minMAF: 0.05).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(filtered[0].Id, Is.EqualTo("V1"));
            Assert.That(filtered[1].Id, Is.EqualTo("V2"));
            Assert.That(filtered[2].Id, Is.EqualTo("V3"));
        });
    }

    #endregion
}
