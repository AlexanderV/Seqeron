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

    /// <summary>
    /// Verifies Bray-Curtis against the worked example from Wikipedia (Bray–Curtis dissimilarity).
    /// Source: https://en.wikipedia.org/wiki/Bray%E2%80%93Curtis_dissimilarity#Example
    /// Tank 1: Goldfish=6, Guppy=7, Rainbow fish=4 (S_j=17)
    /// Tank 2: Goldfish=10, Guppy=0, Rainbow fish=6 (S_k=16)
    /// C_jk = min(6,10)+min(7,0)+min(4,6) = 6+0+4 = 10
    /// BC = 1 - 2*10/(17+16) = 1 - 20/33 ≈ 0.3939
    /// Jaccard: shared={Goldfish,Rainbow}=2, unique1={Guppy}=1, unique2=0
    /// J_dist = 1 - 2/3 = 1/3 ≈ 0.3333
    /// </summary>
    [Test]
    public void CalculateBetaDiversity_WikipediaAquariumExample_MatchesPublishedResult()
    {
        // Arrange — exact data from Wikipedia Bray-Curtis dissimilarity article
        var tank1 = new Dictionary<string, double>
        {
            { "Goldfish", 6 },
            { "Guppy", 7 },
            { "Rainbow fish", 4 }
        };
        var tank2 = new Dictionary<string, double>
        {
            { "Goldfish", 10 },
            { "Guppy", 0 },
            { "Rainbow fish", 6 }
        };

        // Act
        var result = MetagenomicsAnalyzer.CalculateBetaDiversity("Tank 1", tank1, "Tank 2", tank2);

        // Assert — values computed directly from Wikipedia formulas
        Assert.Multiple(() =>
        {
            // BC = 1 - 2*10/33 = 1 - 20/33 = 13/33
            Assert.That(result.BrayCurtis, Is.EqualTo(13.0 / 33.0).Within(1e-10),
                "Bray-Curtis must match Wikipedia worked example: 1 - 20/33 ≈ 0.3939");

            // Jaccard: shared={Goldfish,Rainbow}=2, union={Goldfish,Guppy,Rainbow}=3
            // J_dist = 1 - 2/3 = 1/3
            Assert.That(result.JaccardDistance, Is.EqualTo(1.0 / 3.0).Within(1e-10),
                "Jaccard distance = 1 - |intersection|/|union| = 1 - 2/3 = 1/3");

            Assert.That(result.SharedSpecies, Is.EqualTo(2), "Goldfish and Rainbow fish shared");
            Assert.That(result.UniqueToSample1, Is.EqualTo(1), "Guppy unique to Tank 1");
            Assert.That(result.UniqueToSample2, Is.EqualTo(0), "No species unique to Tank 2");
        });
    }

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
            Assert.That(result.BrayCurtis, Is.EqualTo(0).Within(1e-10), "Identical samples should have Bray-Curtis = 0");
            Assert.That(result.JaccardDistance, Is.EqualTo(0).Within(1e-10), "Identical samples should have Jaccard distance = 0");
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
            Assert.That(result.BrayCurtis, Is.EqualTo(1.0).Within(1e-10), "Disjoint samples should have Bray-Curtis = 1");
            Assert.That(result.JaccardDistance, Is.EqualTo(1.0).Within(1e-10), "Disjoint samples should have Jaccard distance = 1");
            Assert.That(result.SharedSpecies, Is.EqualTo(0), "No shared species");
            Assert.That(result.UniqueToSample1, Is.EqualTo(1), "One species unique to sample1");
            Assert.That(result.UniqueToSample2, Is.EqualTo(1), "One species unique to sample2");
        });
    }

    [Test]
    public void CalculateBetaDiversity_SymmetryProperty_ReturnsIdenticalResults()
    {
        // Arrange
        // sample1={A:0.3, B:0.7}, sample2={B:0.4, C:0.6}
        // allSpecies={A, B, C}, shared={B}=1, unique1={A}=1, unique2={C}=1
        // BC: sumMin = min(0.3,0)+min(0.7,0.4)+min(0,0.6) = 0+0.4+0 = 0.4
        //     sumTotal = 0.3+1.1+0.6 = 2.0
        //     BC = 1 - 2*0.4/2.0 = 1 - 0.4 = 3/5
        // Jaccard: 1 - 1/3 = 2/3
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

        // Assert — first verify exact values, then symmetry
        Assert.Multiple(() =>
        {
            // Exact values for forward direction
            Assert.That(result1.BrayCurtis, Is.EqualTo(3.0 / 5.0).Within(1e-10),
                "BC = 1 - 2*0.4/2.0 = 3/5");
            Assert.That(result1.JaccardDistance, Is.EqualTo(2.0 / 3.0).Within(1e-10),
                "Jaccard = 1 - 1/3 = 2/3");
            Assert.That(result1.SharedSpecies, Is.EqualTo(1), "Species_B shared");
            Assert.That(result1.UniqueToSample1, Is.EqualTo(1), "Species_A unique to sample1");
            Assert.That(result1.UniqueToSample2, Is.EqualTo(1), "Species_C unique to sample2");
            Assert.That(result1.Sample1, Is.EqualTo("S1"), "Sample1 name preserved");
            Assert.That(result1.Sample2, Is.EqualTo("S2"), "Sample2 name preserved");

            // Symmetry: distances identical, species counts swap
            Assert.That(result2.BrayCurtis, Is.EqualTo(result1.BrayCurtis).Within(1e-10),
                "Bray-Curtis should be symmetric");
            Assert.That(result2.JaccardDistance, Is.EqualTo(result1.JaccardDistance).Within(1e-10),
                "Jaccard distance should be symmetric");
            Assert.That(result2.SharedSpecies, Is.EqualTo(result1.SharedSpecies),
                "Shared species count should be symmetric");
            Assert.That(result2.UniqueToSample1, Is.EqualTo(result1.UniqueToSample2),
                "Unique species counts should swap");
            Assert.That(result2.UniqueToSample2, Is.EqualTo(result1.UniqueToSample1),
                "Unique species counts should swap");
            Assert.That(result2.Sample1, Is.EqualTo("S2"), "Sample names follow argument order");
            Assert.That(result2.Sample2, Is.EqualTo("S1"), "Sample names follow argument order");
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

        // Convention: two identical empty inputs → zero distance (see Evidence Design Decisions)
        Assert.Multiple(() =>
        {
            Assert.That(result.BrayCurtis, Is.EqualTo(0.0).Within(1e-10),
                "Empty samples convention: BC = 0 (identical empty inputs)");
            Assert.That(result.JaccardDistance, Is.EqualTo(0.0).Within(1e-10),
                "Empty samples convention: Jaccard = 0 (identical empty inputs)");
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
            Assert.That(result.BrayCurtis, Is.EqualTo(1.0).Within(1e-10), "Empty vs non-empty should have maximum Bray-Curtis");
            Assert.That(result.JaccardDistance, Is.EqualTo(1.0).Within(1e-10), "Empty vs non-empty should have maximum Jaccard distance");
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
            Assert.That(result.BrayCurtis, Is.EqualTo(0).Within(1e-10), "Same species + same abundance should have Bray-Curtis = 0");
            Assert.That(result.JaccardDistance, Is.EqualTo(0).Within(1e-10), "Same species should have Jaccard distance = 0");
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

        // BC = 1 - 2*min(0.8,0.3)/(0.8+0.3) = 1 - 0.6/1.1 = 1 - 6/11 = 5/11
        // Jaccard: same species present in both → J_dist = 0

        // Act
        var result = MetagenomicsAnalyzer.CalculateBetaDiversity("S1", sample1, "S2", sample2);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(result.BrayCurtis, Is.EqualTo(5.0 / 11.0).Within(1e-10),
                "BC = 1 - 2*min(0.8,0.3)/(0.8+0.3) = 5/11");
            Assert.That(result.JaccardDistance, Is.EqualTo(0).Within(1e-10), "Same species should have Jaccard distance = 0");
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
            Assert.That(result.BrayCurtis, Is.EqualTo(1.0).Within(1e-10), "Different species should have Bray-Curtis = 1");
            Assert.That(result.JaccardDistance, Is.EqualTo(1.0).Within(1e-10), "Different species should have Jaccard distance = 1");
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
        // Active species: sample1={A:0.5, C:0.5}, sample2={A:0.3, B:0.7}
        // BC: sumMin = min(0.5,0.3)+min(0,0.7)+min(0.5,0)+min(0,0) = 0.3
        //     sumTotal = 0.8+0.7+0.5+0 = 2.0
        //     BC = 1 - 2*0.3/2.0 = 1 - 0.3 = 0.7
        // Jaccard: shared=1(A), unique1=1(C), unique2=1(B), J_dist = 1 - 1/3 = 2/3
        Assert.Multiple(() =>
        {
            Assert.That(result.SharedSpecies, Is.EqualTo(1), "Only Species_A should be counted as shared");
            Assert.That(result.UniqueToSample1, Is.EqualTo(1), "Species_C unique to sample1");
            Assert.That(result.UniqueToSample2, Is.EqualTo(1), "Species_B unique to sample2");
            Assert.That(result.BrayCurtis, Is.EqualTo(0.7).Within(1e-10),
                "BC = 1 - 2*0.3/2.0 = 0.7 (zero-abundance species excluded from sums)");
            Assert.That(result.JaccardDistance, Is.EqualTo(2.0 / 3.0).Within(1e-10),
                "Jaccard = 1 - 1/3 = 2/3 (zero-abundance species treated as absent)");
        });
    }

    #endregion

    #region Range and Invariant Tests

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
        // BC: sumMin=min(0.2,0.4)+min(0.3,0.1)+0+0+0 = 0.3, sumTotal=0.6+0.4+0.5+0.25+0.25 = 2.0
        //     BC = 1 - 2*0.3/2.0 = 0.7
        // Jaccard: shared=2, unique1=1, unique2=2, J_dist = 1 - 2/5 = 3/5
        Assert.Multiple(() =>
        {
            Assert.That(result.SharedSpecies, Is.EqualTo(2), "Two species shared");
            Assert.That(result.UniqueToSample1, Is.EqualTo(1), "One species unique to sample1");
            Assert.That(result.UniqueToSample2, Is.EqualTo(2), "Two species unique to sample2");

            var totalUniqueSpecies = result.SharedSpecies + result.UniqueToSample1 + result.UniqueToSample2;
            Assert.That(totalUniqueSpecies, Is.EqualTo(5), "Total species count should match");

            Assert.That(result.BrayCurtis, Is.EqualTo(0.7).Within(1e-10),
                "BC = 1 - 2*0.3/2.0 = 0.7");
            Assert.That(result.JaccardDistance, Is.EqualTo(3.0 / 5.0).Within(1e-10),
                "Jaccard = 1 - 2/5 = 3/5");
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
        // Balanced vs sample2: identical → BC=0
        // Skewed vs sample2: sumMin = min(0.9,0.5)+min(0.1,0.5) = 0.6
        //   sumTotal = 1.4+0.6 = 2.0, BC = 1 - 2*0.6/2.0 = 0.4
        Assert.Multiple(() =>
        {
            Assert.That(resultBalanced.BrayCurtis, Is.EqualTo(0.0).Within(1e-10),
                "Balanced identical to sample2 → BC = 0");
            Assert.That(resultSkewed.BrayCurtis, Is.EqualTo(0.4).Within(1e-10),
                "Skewed vs sample2: BC = 1 - 2*0.6/2.0 = 0.4");

            // Jaccard should be identical (presence/absence)
            Assert.That(resultSkewed.JaccardDistance, Is.EqualTo(resultBalanced.JaccardDistance).Within(1e-10),
                "Jaccard distance should ignore abundance differences");

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