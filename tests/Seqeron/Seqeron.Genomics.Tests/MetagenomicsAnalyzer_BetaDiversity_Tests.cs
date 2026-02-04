using NUnit.Framework;
using Seqeron.Genomics;
using System.Collections.Generic;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Tests for META-BETA-001: Beta Diversity Analysis
/// 
/// Tests the CalculateBetaDiversity method which computes Bray-Curtis dissimilarity
/// and Jaccard distance between two ecological samples.
/// 
/// Evidence Sources:
/// - Bray & Curtis (1957): Original Bray-Curtis dissimilarity formula
/// - Jaccard (1901): Original Jaccard index definition  
/// - Whittaker (1960): Beta diversity concept
/// </summary>
[TestFixture]
public class MetagenomicsAnalyzer_BetaDiversity_Tests
{
    #region Mathematical Correctness Tests

    [Test]
    public void CalculateBetaDiversity_IdenticalSamples_ReturnsZeroDistance()
    {
        // Arrange
        var sample = new Dictionary<string, double>
        {
            { "Species_A", 0.4 },
            { "Species_B", 0.6 }
        };

        // Act
        var result = MetagenomicsAnalyzer.CalculateBetaDiversity("Sample1", sample, "Sample2", sample);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(result.BrayCurtis, Is.EqualTo(0).Within(0.001), "Identical samples should have Bray-Curtis = 0");
            Assert.That(result.JaccardDistance, Is.EqualTo(0).Within(0.001), "Identical samples should have Jaccard distance = 0");
            Assert.That(result.Sample1, Is.EqualTo("Sample1"), "Sample1 name should be preserved");
            Assert.That(result.Sample2, Is.EqualTo("Sample2"), "Sample2 name should be preserved");
            Assert.That(result.SharedSpecies, Is.EqualTo(2), "Should identify all species as shared");
            Assert.That(result.UniqueToSample1, Is.EqualTo(0), "No species unique to sample1");
            Assert.That(result.UniqueToSample2, Is.EqualTo(0), "No species unique to sample2");
        });
    }

    [Test]
    public void CalculateBetaDiversity_CompletelyDisjointSamples_ReturnsMaximumDistance()
    {
        // Arrange
        var sample1 = new Dictionary<string, double> { { "Species_A", 1.0 } };
        var sample2 = new Dictionary<string, double> { { "Species_B", 1.0 } };

        // Act
        var result = MetagenomicsAnalyzer.CalculateBetaDiversity("S1", sample1, "S2", sample2);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(result.BrayCurtis, Is.EqualTo(1.0).Within(0.001), "Disjoint samples should have Bray-Curtis = 1");
            Assert.That(result.JaccardDistance, Is.EqualTo(1.0).Within(0.001), "Disjoint samples should have Jaccard distance = 1");
            Assert.That(result.SharedSpecies, Is.EqualTo(0), "No shared species");
            Assert.That(result.UniqueToSample1, Is.EqualTo(1), "One species unique to sample1");
            Assert.That(result.UniqueToSample2, Is.EqualTo(1), "One species unique to sample2");
        });
    }

    [Test]
    public void CalculateBetaDiversity_SymmetryProperty_ReturnsIdenticalResults()
    {
        // Arrange
        var sample1 = new Dictionary<string, double>
        {
            { "Species_A", 0.3 },
            { "Species_B", 0.7 }
        };
        var sample2 = new Dictionary<string, double>
        {
            { "Species_B", 0.4 },
            { "Species_C", 0.6 }
        };

        // Act
        var result1 = MetagenomicsAnalyzer.CalculateBetaDiversity("S1", sample1, "S2", sample2);
        var result2 = MetagenomicsAnalyzer.CalculateBetaDiversity("S2", sample2, "S1", sample1);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(result2.BrayCurtis, Is.EqualTo(result1.BrayCurtis).Within(0.001), "Bray-Curtis should be symmetric");
            Assert.That(result2.JaccardDistance, Is.EqualTo(result1.JaccardDistance).Within(0.001), "Jaccard distance should be symmetric");
            Assert.That(result2.SharedSpecies, Is.EqualTo(result1.SharedSpecies), "Shared species count should be symmetric");
            Assert.That(result2.UniqueToSample1, Is.EqualTo(result1.UniqueToSample2), "Unique species counts should swap");
            Assert.That(result2.UniqueToSample2, Is.EqualTo(result1.UniqueToSample1), "Unique species counts should swap");
        });
    }

    #endregion

    #region Edge Case Tests

    [Test]
    public void CalculateBetaDiversity_BothSamplesEmpty_HandlesGracefully()
    {
        // Arrange
        var emptySample = new Dictionary<string, double>();

        // Act & Assert - Should not throw
        var result = MetagenomicsAnalyzer.CalculateBetaDiversity("Empty1", emptySample, "Empty2", emptySample);

        Assert.Multiple(() =>
        {
            Assert.That(result.Sample1, Is.EqualTo("Empty1"), "Sample names preserved");
            Assert.That(result.Sample2, Is.EqualTo("Empty2"), "Sample names preserved");
            Assert.That(result.SharedSpecies, Is.EqualTo(0), "No species in empty samples");
            Assert.That(result.UniqueToSample1, Is.EqualTo(0), "No unique species in empty sample");
            Assert.That(result.UniqueToSample2, Is.EqualTo(0), "No unique species in empty sample");
        });
    }

    [Test]
    public void CalculateBetaDiversity_OneSampleEmpty_ReturnsMaximumDistance()
    {
        // Arrange
        var nonEmptySample = new Dictionary<string, double> { { "Species_A", 1.0 } };
        var emptySample = new Dictionary<string, double>();

        // Act
        var result = MetagenomicsAnalyzer.CalculateBetaDiversity("NonEmpty", nonEmptySample, "Empty", emptySample);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(result.BrayCurtis, Is.EqualTo(1.0).Within(0.001), "Empty vs non-empty should have maximum Bray-Curtis");
            Assert.That(result.JaccardDistance, Is.EqualTo(1.0).Within(0.001), "Empty vs non-empty should have maximum Jaccard distance");
            Assert.That(result.SharedSpecies, Is.EqualTo(0), "No shared species with empty sample");
            Assert.That(result.UniqueToSample1, Is.EqualTo(1), "All species unique to non-empty sample");
            Assert.That(result.UniqueToSample2, Is.EqualTo(0), "No unique species in empty sample");
        });
    }

    [Test]
    public void CalculateBetaDiversity_SameSingleSpeciesSameAbundance_ReturnsZeroDistance()
    {
        // Arrange - same species AND same abundance
        var sample1 = new Dictionary<string, double> { { "Species_X", 0.5 } };
        var sample2 = new Dictionary<string, double> { { "Species_X", 0.5 } };

        // Act
        var result = MetagenomicsAnalyzer.CalculateBetaDiversity("S1", sample1, "S2", sample2);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(result.BrayCurtis, Is.EqualTo(0).Within(0.001), "Same species + same abundance should have Bray-Curtis = 0");
            Assert.That(result.JaccardDistance, Is.EqualTo(0).Within(0.001), "Same species should have Jaccard distance = 0");
            Assert.That(result.SharedSpecies, Is.EqualTo(1), "One shared species");
            Assert.That(result.UniqueToSample1, Is.EqualTo(0), "No unique species");
            Assert.That(result.UniqueToSample2, Is.EqualTo(0), "No unique species");
        });
    }

    [Test]
    public void CalculateBetaDiversity_SameSingleSpeciesDifferentAbundance_BrayCurtisNonZero()
    {
        // Arrange - same species but different abundance
        // Bray-Curtis is abundance-sensitive: BC = 1 - 2*min(a1,a2)/(a1+a2)
        var sample1 = new Dictionary<string, double> { { "Species_X", 0.8 } };
        var sample2 = new Dictionary<string, double> { { "Species_X", 0.3 } };

        // Expected: BC = 1 - 2*0.3/1.1 ≈ 0.455
        // Jaccard should still be 0 (presence/absence only)

        // Act
        var result = MetagenomicsAnalyzer.CalculateBetaDiversity("S1", sample1, "S2", sample2);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(result.BrayCurtis, Is.GreaterThan(0), "Different abundances should give non-zero Bray-Curtis");
            Assert.That(result.BrayCurtis, Is.EqualTo(0.455).Within(0.01), "Calculated Bray-Curtis for 0.8 vs 0.3");
            Assert.That(result.JaccardDistance, Is.EqualTo(0).Within(0.001), "Same species should have Jaccard distance = 0");
            Assert.That(result.SharedSpecies, Is.EqualTo(1), "One shared species");
        });
    }

    [Test]
    public void CalculateBetaDiversity_DifferentSingleSpecies_ReturnsMaximumDistance()
    {
        // Arrange
        var sample1 = new Dictionary<string, double> { { "Species_X", 1.0 } };
        var sample2 = new Dictionary<string, double> { { "Species_Y", 1.0 } };

        // Act
        var result = MetagenomicsAnalyzer.CalculateBetaDiversity("S1", sample1, "S2", sample2);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(result.BrayCurtis, Is.EqualTo(1.0).Within(0.001), "Different species should have Bray-Curtis = 1");
            Assert.That(result.JaccardDistance, Is.EqualTo(1.0).Within(0.001), "Different species should have Jaccard distance = 1");
            Assert.That(result.SharedSpecies, Is.EqualTo(0), "No shared species");
            Assert.That(result.UniqueToSample1, Is.EqualTo(1), "One unique to sample1");
            Assert.That(result.UniqueToSample2, Is.EqualTo(1), "One unique to sample2");
        });
    }

    [Test]
    public void CalculateBetaDiversity_ZeroAbundanceHandling_TreatsAsAbsent()
    {
        // Arrange
        var sample1 = new Dictionary<string, double>
        {
            { "Species_A", 0.5 },
            { "Species_B", 0.0 },  // Zero abundance = absent
            { "Species_C", 0.5 }
        };
        var sample2 = new Dictionary<string, double>
        {
            { "Species_A", 0.3 },
            { "Species_B", 0.7 },  // Present in sample2
            { "Species_D", 0.0 }   // Zero abundance = absent
        };

        // Act
        var result = MetagenomicsAnalyzer.CalculateBetaDiversity("S1", sample1, "S2", sample2);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(result.SharedSpecies, Is.EqualTo(1), "Only Species_A should be counted as shared");
            Assert.That(result.UniqueToSample1, Is.EqualTo(1), "Species_C unique to sample1");
            Assert.That(result.UniqueToSample2, Is.EqualTo(1), "Species_B unique to sample2");
            // Species with 0 abundance should not be counted
        });
    }

    #endregion

    #region Range and Invariant Tests

    [Test]
    public void CalculateBetaDiversity_DistanceRangeConstraints_AlwaysValid()
    {
        // Arrange - Test multiple scenarios
        var scenarios = new[]
        {
            (new Dictionary<string, double> { { "A", 1 } }, new Dictionary<string, double> { { "A", 1 } }),
            (new Dictionary<string, double> { { "A", 1 } }, new Dictionary<string, double> { { "B", 1 } }),
            (new Dictionary<string, double> { { "A", 0.3 }, { "B", 0.7 } }, new Dictionary<string, double> { { "B", 0.4 }, { "C", 0.6 } }),
            (new Dictionary<string, double>(), new Dictionary<string, double> { { "A", 1 } })
        };

        foreach (var (sample1, sample2) in scenarios)
        {
            // Act
            var result = MetagenomicsAnalyzer.CalculateBetaDiversity("S1", sample1, "S2", sample2);

            // Assert
            Assert.Multiple(() =>
            {
                Assert.That(result.BrayCurtis, Is.InRange(0.0, 1.0), "Bray-Curtis must be in [0,1]");
                Assert.That(result.JaccardDistance, Is.InRange(0.0, 1.0), "Jaccard distance must be in [0,1]");
            });
        }
    }

    [Test]
    public void CalculateBetaDiversity_SpeciesCountAccuracy_MatchesExpected()
    {
        // Arrange
        var sample1 = new Dictionary<string, double>
        {
            { "Shared1", 0.2 },
            { "Shared2", 0.3 },
            { "Unique1", 0.5 }
        };
        var sample2 = new Dictionary<string, double>
        {
            { "Shared1", 0.4 },
            { "Shared2", 0.1 },
            { "Unique2A", 0.25 },
            { "Unique2B", 0.25 }
        };

        // Act
        var result = MetagenomicsAnalyzer.CalculateBetaDiversity("S1", sample1, "S2", sample2);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(result.SharedSpecies, Is.EqualTo(2), "Two species shared");
            Assert.That(result.UniqueToSample1, Is.EqualTo(1), "One species unique to sample1");
            Assert.That(result.UniqueToSample2, Is.EqualTo(2), "Two species unique to sample2");

            // Verify total species count integrity
            var totalUniqueSpecies = result.SharedSpecies + result.UniqueToSample1 + result.UniqueToSample2;
            Assert.That(totalUniqueSpecies, Is.EqualTo(5), "Total species count should match");
        });
    }

    #endregion

    #region Abundance vs Presence Tests

    [Test]
    public void CalculateBetaDiversity_BrayCurtisRespondsToAbundance_JaccardIgnoresIt()
    {
        // Arrange - Same species, different abundance distributions
        var sample1Balanced = new Dictionary<string, double>
        {
            { "Species_A", 0.5 },
            { "Species_B", 0.5 }
        };
        var sample1Skewed = new Dictionary<string, double>
        {
            { "Species_A", 0.9 },
            { "Species_B", 0.1 }
        };
        var sample2 = new Dictionary<string, double>
        {
            { "Species_A", 0.5 },
            { "Species_B", 0.5 }
        };

        // Act
        var resultBalanced = MetagenomicsAnalyzer.CalculateBetaDiversity("Balanced", sample1Balanced, "S2", sample2);
        var resultSkewed = MetagenomicsAnalyzer.CalculateBetaDiversity("Skewed", sample1Skewed, "S2", sample2);

        // Assert
        Assert.Multiple(() =>
        {
            // Jaccard should be identical (presence/absence)
            Assert.That(resultSkewed.JaccardDistance, Is.EqualTo(resultBalanced.JaccardDistance).Within(0.001),
                "Jaccard distance should ignore abundance differences");

            // Bray-Curtis should be different (abundance-sensitive)
            Assert.That(resultSkewed.BrayCurtis, Is.Not.EqualTo(resultBalanced.BrayCurtis),
                "Bray-Curtis should respond to abundance differences");

            // Species counts should be identical
            Assert.That(resultSkewed.SharedSpecies, Is.EqualTo(resultBalanced.SharedSpecies),
                "Shared species count should be identical");
        });
    }

    #endregion

    #region Sample Name Preservation

    [Test]
    public void CalculateBetaDiversity_ComplexSampleNames_PreservedCorrectly()
    {
        // Arrange
        var complexName1 = "Site-α_2024.01!@#";
        var complexName2 = "Control_Group_µM_Concentration";
        var sample1 = new Dictionary<string, double> { { "Species", 1.0 } };
        var sample2 = new Dictionary<string, double> { { "Species", 1.0 } };

        // Act
        var result = MetagenomicsAnalyzer.CalculateBetaDiversity(complexName1, sample1, complexName2, sample2);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(result.Sample1, Is.EqualTo(complexName1), "Complex sample name 1 should be preserved");
            Assert.That(result.Sample2, Is.EqualTo(complexName2), "Complex sample name 2 should be preserved");
        });
    }

    #endregion
}