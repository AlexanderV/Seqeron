using NUnit.Framework;
using Seqeron.Genomics;
using System;
using System.Collections.Generic;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Tests for MetagenomicsAnalyzer.CalculateAlphaDiversity (META-ALPHA-001).
/// 
/// Evidence Sources:
/// - Wikipedia (Diversity Index, Alpha Diversity, Species Richness, Species Evenness)
/// - Shannon (1948): Shannon entropy formula
/// - Simpson (1949): Simpson concentration index
/// - Hill (1973): Unified diversity notation
/// - Chao (1984): Chao1 estimator
/// </summary>
[TestFixture]
public class MetagenomicsAnalyzer_AlphaDiversity_Tests
{
    // Mathematical constants for test assertions
    private const double Ln2 = 0.6931471805599453;  // ln(2)
    private const double Ln4 = 1.3862943611198906;  // ln(4)
    private const double Tolerance = 1e-10;

    #region Empty and Null Input Tests

    /// <summary>
    /// M1: Empty abundances should return all metrics = 0.
    /// Evidence: Standard robustness — no data implies no diversity.
    /// </summary>
    [Test]
    public void CalculateAlphaDiversity_EmptyAbundances_ReturnsAllZeros()
    {
        var abundances = new Dictionary<string, double>();

        var diversity = MetagenomicsAnalyzer.CalculateAlphaDiversity(abundances);

        Assert.Multiple(() =>
        {
            Assert.That(diversity.ObservedSpecies, Is.EqualTo(0), "ObservedSpecies");
            Assert.That(diversity.ShannonIndex, Is.EqualTo(0), "ShannonIndex");
            Assert.That(diversity.SimpsonIndex, Is.EqualTo(0), "SimpsonIndex");
            Assert.That(diversity.InverseSimpson, Is.EqualTo(0), "InverseSimpson");
            Assert.That(diversity.PielouEvenness, Is.EqualTo(0), "PielouEvenness");
            Assert.That(diversity.Chao1Estimate, Is.EqualTo(0), "Chao1Estimate");
        });
    }

    /// <summary>
    /// M2: Null abundances should return all metrics = 0.
    /// Evidence: Standard robustness — null safety.
    /// </summary>
    [Test]
    public void CalculateAlphaDiversity_NullAbundances_ReturnsAllZeros()
    {
        var diversity = MetagenomicsAnalyzer.CalculateAlphaDiversity(null!);

        Assert.Multiple(() =>
        {
            Assert.That(diversity.ObservedSpecies, Is.EqualTo(0), "ObservedSpecies");
            Assert.That(diversity.ShannonIndex, Is.EqualTo(0), "ShannonIndex");
            Assert.That(diversity.SimpsonIndex, Is.EqualTo(0), "SimpsonIndex");
            Assert.That(diversity.InverseSimpson, Is.EqualTo(0), "InverseSimpson");
            Assert.That(diversity.PielouEvenness, Is.EqualTo(0), "PielouEvenness");
            Assert.That(diversity.Chao1Estimate, Is.EqualTo(0), "Chao1Estimate");
        });
    }

    #endregion

    #region Single Species Tests

    /// <summary>
    /// M3-M7: Single species returns:
    /// - Shannon = 0 (no uncertainty: −1·ln(1) = 0)
    /// - Simpson = 1.0 (complete dominance: 1² = 1)
    /// - InverseSimpson = 1.0 (1/1 = 1)
    /// - Pielou = 0 (undefined when S=1, convention returns 0)
    /// - ObservedSpecies = 1
    /// 
    /// Evidence: Shannon (1948), Simpson (1949), Hill (1973)
    /// </summary>
    [Test]
    public void CalculateAlphaDiversity_SingleSpecies_ReturnsCorrectMetrics()
    {
        var abundances = new Dictionary<string, double>
        {
            { "OnlySpecies", 1.0 }
        };

        var diversity = MetagenomicsAnalyzer.CalculateAlphaDiversity(abundances);

        Assert.Multiple(() =>
        {
            Assert.That(diversity.ObservedSpecies, Is.EqualTo(1), "ObservedSpecies");
            Assert.That(diversity.ShannonIndex, Is.EqualTo(0).Within(Tolerance), "Shannon = −1·ln(1) = 0");
            Assert.That(diversity.SimpsonIndex, Is.EqualTo(1.0).Within(Tolerance), "Simpson = 1² = 1");
            Assert.That(diversity.InverseSimpson, Is.EqualTo(1.0).Within(Tolerance), "InverseSimpson = 1/1");
            Assert.That(diversity.PielouEvenness, Is.EqualTo(0).Within(Tolerance), "Pielou undefined for S=1");
            Assert.That(diversity.Chao1Estimate, Is.EqualTo(1), "Chao1 = ObservedSpecies");
        });
    }

    #endregion

    #region Even Distribution Tests

    /// <summary>
    /// M8-M10: Two equal species (0.5, 0.5) returns:
    /// - Shannon = ln(2) ≈ 0.693 (−2×(0.5×ln(0.5)) = ln(2))
    /// - Simpson = 0.5 (0.5² + 0.5² = 0.5)
    /// - InverseSimpson = 2.0 (1/0.5)
    /// - Pielou = 1.0 (perfect evenness: ln(2)/ln(2) = 1)
    /// 
    /// Evidence: Shannon (1948), Simpson (1949)
    /// </summary>
    [Test]
    public void CalculateAlphaDiversity_TwoEqualSpecies_ReturnsCorrectMetrics()
    {
        var abundances = new Dictionary<string, double>
        {
            { "Species1", 0.5 },
            { "Species2", 0.5 }
        };

        var diversity = MetagenomicsAnalyzer.CalculateAlphaDiversity(abundances);

        Assert.Multiple(() =>
        {
            Assert.That(diversity.ObservedSpecies, Is.EqualTo(2), "ObservedSpecies");
            Assert.That(diversity.ShannonIndex, Is.EqualTo(Ln2).Within(Tolerance), "Shannon = ln(2)");
            Assert.That(diversity.SimpsonIndex, Is.EqualTo(0.5).Within(Tolerance), "Simpson = 2×0.5² = 0.5");
            Assert.That(diversity.InverseSimpson, Is.EqualTo(2.0).Within(Tolerance), "InverseSimpson = 1/0.5 = 2");
            Assert.That(diversity.PielouEvenness, Is.EqualTo(1.0).Within(Tolerance), "Pielou = ln(2)/ln(2) = 1");
        });
    }

    /// <summary>
    /// M11-M13: Four equal species (0.25 each) returns:
    /// - Shannon = ln(4) ≈ 1.386
    /// - Simpson = 0.25 (4×0.25² = 0.25)
    /// - InverseSimpson = 4.0 (equals species count for even distribution)
    /// - Pielou = 1.0 (perfect evenness)
    /// 
    /// Evidence: Shannon (1948), Simpson (1949), Hill (1973)
    /// </summary>
    [Test]
    public void CalculateAlphaDiversity_FourEqualSpecies_ReturnsCorrectMetrics()
    {
        var abundances = new Dictionary<string, double>
        {
            { "Species1", 0.25 },
            { "Species2", 0.25 },
            { "Species3", 0.25 },
            { "Species4", 0.25 }
        };

        var diversity = MetagenomicsAnalyzer.CalculateAlphaDiversity(abundances);

        Assert.Multiple(() =>
        {
            Assert.That(diversity.ObservedSpecies, Is.EqualTo(4), "ObservedSpecies");
            Assert.That(diversity.ShannonIndex, Is.EqualTo(Ln4).Within(Tolerance), "Shannon = ln(4)");
            Assert.That(diversity.SimpsonIndex, Is.EqualTo(0.25).Within(Tolerance), "Simpson = 4×0.25² = 0.25");
            Assert.That(diversity.InverseSimpson, Is.EqualTo(4.0).Within(Tolerance), "InverseSimpson = S for even dist");
            Assert.That(diversity.PielouEvenness, Is.EqualTo(1.0).Within(Tolerance), "Pielou = 1.0 for even dist");
        });
    }

    /// <summary>
    /// S4: InverseSimpson equals species count for perfectly even distribution.
    /// Evidence: Hill (1973) — effective number of species.
    /// </summary>
    [TestCase(2)]
    [TestCase(5)]
    [TestCase(10)]
    [TestCase(100)]
    public void CalculateAlphaDiversity_EvenDistribution_InverseSimpsonEqualsSpeciesCount(int speciesCount)
    {
        var abundances = new Dictionary<string, double>();
        double proportion = 1.0 / speciesCount;
        for (int i = 0; i < speciesCount; i++)
        {
            abundances[$"Species{i}"] = proportion;
        }

        var diversity = MetagenomicsAnalyzer.CalculateAlphaDiversity(abundances);

        Assert.Multiple(() =>
        {
            Assert.That(diversity.ObservedSpecies, Is.EqualTo(speciesCount), "ObservedSpecies");
            Assert.That(diversity.InverseSimpson, Is.EqualTo(speciesCount).Within(1e-6),
                $"InverseSimpson should equal {speciesCount} for even distribution");
            Assert.That(diversity.PielouEvenness, Is.EqualTo(1.0).Within(1e-6), "Pielou = 1.0 for even dist");
        });
    }

    #endregion

    #region Uneven Distribution Tests

    /// <summary>
    /// S1: Highly uneven distribution produces low Shannon, high Simpson.
    /// Evidence: Diversity theory — dominance reduces diversity.
    /// </summary>
    [Test]
    public void CalculateAlphaDiversity_HighlyUneven_LowDiversityHighDominance()
    {
        var abundances = new Dictionary<string, double>
        {
            { "DominantSpecies", 0.99 },
            { "RareSpecies", 0.01 }
        };

        // Expected: H = −(0.99×ln(0.99) + 0.01×ln(0.01)) ≈ 0.056
        // Expected: λ = 0.99² + 0.01² ≈ 0.9802
        double expectedShannon = -(0.99 * Math.Log(0.99) + 0.01 * Math.Log(0.01));
        double expectedSimpson = 0.99 * 0.99 + 0.01 * 0.01;

        var diversity = MetagenomicsAnalyzer.CalculateAlphaDiversity(abundances);

        Assert.Multiple(() =>
        {
            Assert.That(diversity.ShannonIndex, Is.EqualTo(expectedShannon).Within(1e-6), "Shannon low for uneven");
            Assert.That(diversity.SimpsonIndex, Is.EqualTo(expectedSimpson).Within(1e-6), "Simpson high for uneven");
            Assert.That(diversity.ShannonIndex, Is.LessThan(0.1), "Shannon < 0.1 for extreme dominance");
            Assert.That(diversity.SimpsonIndex, Is.GreaterThan(0.9), "Simpson > 0.9 for extreme dominance");
            Assert.That(diversity.PielouEvenness, Is.LessThan(0.2), "Pielou < 0.2 for extreme unevenness");
        });
    }

    /// <summary>
    /// Uneven distribution (0.9, 0.05, 0.05) produces predictable values.
    /// </summary>
    [Test]
    public void CalculateAlphaDiversity_UnevenThreeSpecies_CalculatesCorrectly()
    {
        var abundances = new Dictionary<string, double>
        {
            { "DominantSpecies", 0.9 },
            { "RareSpecies1", 0.05 },
            { "RareSpecies2", 0.05 }
        };

        // Expected: H = −(0.9×ln(0.9) + 2×0.05×ln(0.05))
        double expectedShannon = -(0.9 * Math.Log(0.9) + 0.05 * Math.Log(0.05) + 0.05 * Math.Log(0.05));
        double expectedSimpson = 0.9 * 0.9 + 0.05 * 0.05 + 0.05 * 0.05;
        double expectedPielou = expectedShannon / Math.Log(3);

        var diversity = MetagenomicsAnalyzer.CalculateAlphaDiversity(abundances);

        Assert.Multiple(() =>
        {
            Assert.That(diversity.ObservedSpecies, Is.EqualTo(3), "ObservedSpecies");
            Assert.That(diversity.ShannonIndex, Is.EqualTo(expectedShannon).Within(1e-6), "Shannon");
            Assert.That(diversity.SimpsonIndex, Is.EqualTo(expectedSimpson).Within(1e-6), "Simpson");
            Assert.That(diversity.PielouEvenness, Is.EqualTo(expectedPielou).Within(1e-6), "Pielou");
        });
    }

    #endregion

    #region Edge Case Tests

    /// <summary>
    /// M14: Zero abundances are filtered out (ln(0) is undefined).
    /// Evidence: Implementation requirement — zero abundance means absent species.
    /// </summary>
    [Test]
    public void CalculateAlphaDiversity_ZeroAbundances_FilteredOut()
    {
        var abundances = new Dictionary<string, double>
        {
            { "Present1", 0.5 },
            { "Absent", 0.0 },
            { "Present2", 0.5 }
        };

        var diversity = MetagenomicsAnalyzer.CalculateAlphaDiversity(abundances);

        Assert.Multiple(() =>
        {
            Assert.That(diversity.ObservedSpecies, Is.EqualTo(2), "Zero abundance species not counted");
            Assert.That(diversity.ShannonIndex, Is.EqualTo(Ln2).Within(Tolerance), "Shannon as if 2 equal species");
            Assert.That(diversity.SimpsonIndex, Is.EqualTo(0.5).Within(Tolerance), "Simpson as if 2 equal species");
        });
    }

    /// <summary>
    /// M14: All zero abundances behave like empty input.
    /// </summary>
    [Test]
    public void CalculateAlphaDiversity_AllZeroAbundances_ReturnsZeros()
    {
        var abundances = new Dictionary<string, double>
        {
            { "Species1", 0.0 },
            { "Species2", 0.0 }
        };

        var diversity = MetagenomicsAnalyzer.CalculateAlphaDiversity(abundances);

        Assert.Multiple(() =>
        {
            Assert.That(diversity.ObservedSpecies, Is.EqualTo(0), "No species with positive abundance");
            Assert.That(diversity.ShannonIndex, Is.EqualTo(0), "Shannon = 0");
            Assert.That(diversity.SimpsonIndex, Is.EqualTo(0), "Simpson = 0");
        });
    }

    /// <summary>
    /// M15: Unnormalized abundances are normalized before calculation.
    /// Evidence: Implementation — abundances may not sum to 1.
    /// </summary>
    [Test]
    public void CalculateAlphaDiversity_UnnormalizedAbundances_Normalized()
    {
        // Abundances sum to 2.0, not 1.0
        var abundances = new Dictionary<string, double>
        {
            { "Species1", 1.0 },
            { "Species2", 1.0 }
        };

        var diversity = MetagenomicsAnalyzer.CalculateAlphaDiversity(abundances);

        // After normalization: (0.5, 0.5) → same as two equal species
        Assert.Multiple(() =>
        {
            Assert.That(diversity.ObservedSpecies, Is.EqualTo(2), "ObservedSpecies");
            Assert.That(diversity.ShannonIndex, Is.EqualTo(Ln2).Within(Tolerance), "Shannon = ln(2) after normalization");
            Assert.That(diversity.SimpsonIndex, Is.EqualTo(0.5).Within(Tolerance), "Simpson = 0.5 after normalization");
            Assert.That(diversity.PielouEvenness, Is.EqualTo(1.0).Within(Tolerance), "Pielou = 1.0 after normalization");
        });
    }

    /// <summary>
    /// C1: Large taxon count handled correctly (scalability).
    /// </summary>
    [Test]
    public void CalculateAlphaDiversity_LargeTaxonCount_HandledCorrectly()
    {
        const int taxonCount = 10000;
        var abundances = new Dictionary<string, double>();
        double proportion = 1.0 / taxonCount;
        for (int i = 0; i < taxonCount; i++)
        {
            abundances[$"Taxon{i}"] = proportion;
        }

        var diversity = MetagenomicsAnalyzer.CalculateAlphaDiversity(abundances);

        double expectedShannon = Math.Log(taxonCount);
        double expectedSimpson = 1.0 / taxonCount;

        Assert.Multiple(() =>
        {
            Assert.That(diversity.ObservedSpecies, Is.EqualTo(taxonCount), "ObservedSpecies");
            Assert.That(diversity.ShannonIndex, Is.EqualTo(expectedShannon).Within(1e-6), "Shannon = ln(S)");
            Assert.That(diversity.SimpsonIndex, Is.EqualTo(expectedSimpson).Within(1e-9), "Simpson = 1/S");
            Assert.That(diversity.InverseSimpson, Is.EqualTo(taxonCount).Within(0.001), "InverseSimpson = S");
            Assert.That(diversity.PielouEvenness, Is.EqualTo(1.0).Within(1e-6), "Pielou = 1.0");
        });
    }

    /// <summary>
    /// C2: Numerical stability with very small abundances.
    /// </summary>
    [Test]
    public void CalculateAlphaDiversity_VerySmallAbundances_NumericallyStable()
    {
        var abundances = new Dictionary<string, double>
        {
            { "Dominant", 0.999999 },
            { "VeryRare", 0.000001 }
        };

        var diversity = MetagenomicsAnalyzer.CalculateAlphaDiversity(abundances);

        Assert.Multiple(() =>
        {
            Assert.That(diversity.ObservedSpecies, Is.EqualTo(2), "Both species counted");
            Assert.That(diversity.ShannonIndex, Is.GreaterThan(0), "Shannon > 0");
            Assert.That(diversity.ShannonIndex, Is.LessThan(0.001), "Shannon very low for extreme dominance");
            Assert.That(double.IsFinite(diversity.ShannonIndex), Is.True, "Shannon is finite");
            Assert.That(double.IsFinite(diversity.SimpsonIndex), Is.True, "Simpson is finite");
            Assert.That(double.IsNaN(diversity.ShannonIndex), Is.False, "Shannon is not NaN");
        });
    }

    #endregion

    #region Invariant Tests

    /// <summary>
    /// M16-M19: All results satisfy theoretical bounds.
    /// Evidence: Shannon ≥ 0, Simpson ∈ [0,1], Pielou ∈ [0,1], Chao1 ≥ S.
    /// </summary>
    [TestCase(0.5, 0.5)]
    [TestCase(0.9, 0.1)]
    [TestCase(0.25, 0.25, 0.25, 0.25)]
    [TestCase(0.7, 0.2, 0.1)]
    [TestCase(0.99, 0.01)]
    [TestCase(1.0)]
    public void CalculateAlphaDiversity_VariousDistributions_SatisfyTheoreticalBounds(params double[] proportions)
    {
        var abundances = new Dictionary<string, double>();
        for (int i = 0; i < proportions.Length; i++)
        {
            abundances[$"Species{i}"] = proportions[i];
        }

        var diversity = MetagenomicsAnalyzer.CalculateAlphaDiversity(abundances);

        Assert.Multiple(() =>
        {
            // M16: Shannon ≥ 0
            Assert.That(diversity.ShannonIndex, Is.GreaterThanOrEqualTo(0), "Shannon ≥ 0");

            // M17: Simpson ∈ [0, 1]
            Assert.That(diversity.SimpsonIndex, Is.InRange(0, 1), "Simpson ∈ [0, 1]");

            // M18: Pielou ∈ [0, 1] (for S > 1) or = 0 (for S ≤ 1)
            Assert.That(diversity.PielouEvenness, Is.InRange(0, 1), "Pielou ∈ [0, 1]");

            // M19: Chao1 ≥ ObservedSpecies
            Assert.That(diversity.Chao1Estimate, Is.GreaterThanOrEqualTo(diversity.ObservedSpecies),
                "Chao1 ≥ ObservedSpecies");

            // InverseSimpson ≥ 1 when Simpson > 0
            if (diversity.SimpsonIndex > 0)
            {
                Assert.That(diversity.InverseSimpson, Is.GreaterThanOrEqualTo(1), "InverseSimpson ≥ 1");
            }

            // All values finite
            Assert.That(double.IsFinite(diversity.ShannonIndex), Is.True, "Shannon finite");
            Assert.That(double.IsFinite(diversity.SimpsonIndex), Is.True, "Simpson finite");
            Assert.That(double.IsFinite(diversity.InverseSimpson), Is.True, "InverseSimpson finite");
            Assert.That(double.IsFinite(diversity.PielouEvenness), Is.True, "Pielou finite");
        });
    }

    /// <summary>
    /// S2: Shannon increases with increasing richness for even distribution.
    /// Evidence: Shannon theory — more species = more uncertainty.
    /// </summary>
    [Test]
    public void CalculateAlphaDiversity_EvenDistribution_ShannonIncreasesWithRichness()
    {
        double previousShannon = 0;

        foreach (int speciesCount in new[] { 1, 2, 4, 8, 16 })
        {
            var abundances = new Dictionary<string, double>();
            double proportion = 1.0 / speciesCount;
            for (int i = 0; i < speciesCount; i++)
            {
                abundances[$"Species{i}"] = proportion;
            }

            var diversity = MetagenomicsAnalyzer.CalculateAlphaDiversity(abundances);

            if (speciesCount > 1)
            {
                Assert.That(diversity.ShannonIndex, Is.GreaterThan(previousShannon),
                    $"Shannon for S={speciesCount} should exceed S={speciesCount / 2}");
            }

            previousShannon = diversity.ShannonIndex;
        }
    }

    /// <summary>
    /// S3: Simpson decreases with increasing richness for even distribution.
    /// Evidence: Simpson theory — more species = less dominance.
    /// </summary>
    [Test]
    public void CalculateAlphaDiversity_EvenDistribution_SimpsonDecreasesWithRichness()
    {
        double previousSimpson = 1.0;

        foreach (int speciesCount in new[] { 1, 2, 4, 8, 16 })
        {
            var abundances = new Dictionary<string, double>();
            double proportion = 1.0 / speciesCount;
            for (int i = 0; i < speciesCount; i++)
            {
                abundances[$"Species{i}"] = proportion;
            }

            var diversity = MetagenomicsAnalyzer.CalculateAlphaDiversity(abundances);

            if (speciesCount > 1)
            {
                Assert.That(diversity.SimpsonIndex, Is.LessThan(previousSimpson),
                    $"Simpson for S={speciesCount} should be less than S={speciesCount / 2}");
            }

            previousSimpson = diversity.SimpsonIndex;
        }
    }

    #endregion
}
