using NUnit.Framework;
using Seqeron.Genomics;
using System.Collections.Generic;
using System.Linq;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Tests for POP-DIV-001: Diversity Statistics
/// Covers: CalculateNucleotideDiversity, CalculateWattersonTheta, 
///         CalculateTajimasD, CalculateDiversityStatistics
/// </summary>
/// <remarks>
/// Evidence: Wikipedia (Nucleotide diversity, Watterson estimator, Tajima's D),
///           Nei & Li (1979), Watterson (1975), Tajima (1989)
/// </remarks>
[TestFixture]
public class PopulationGeneticsAnalyzer_Diversity_Tests
{
    #region CalculateNucleotideDiversity Tests

    [Test]
    [Description("ND-M01: Identical sequences should have π = 0")]
    public void CalculateNucleotideDiversity_IdenticalSequences_ReturnsZero()
    {
        // Arrange - all sequences identical (no polymorphism)
        var sequences = new List<IReadOnlyList<char>>
        {
            "ACGT".ToList(),
            "ACGT".ToList(),
            "ACGT".ToList()
        };

        // Act
        double pi = PopulationGeneticsAnalyzer.CalculateNucleotideDiversity(sequences);

        // Assert
        Assert.That(pi, Is.EqualTo(0));
    }

    [Test]
    [Description("ND-M02: Two sequences differing at all positions should have π = 1.0")]
    public void CalculateNucleotideDiversity_AllDifferent_ReturnsOne()
    {
        // Arrange - two sequences that differ at every position
        var sequences = new List<IReadOnlyList<char>>
        {
            "AAAA".ToList(),
            "TTTT".ToList()
        };

        // Act
        double pi = PopulationGeneticsAnalyzer.CalculateNucleotideDiversity(sequences);

        // Assert - 4 differences / (1 comparison × 4 positions) = 1.0
        Assert.That(pi, Is.EqualTo(1.0));
    }

    [Test]
    [Description("ND-M03: Single sequence returns π = 0 (undefined for n < 2)")]
    public void CalculateNucleotideDiversity_SingleSequence_ReturnsZero()
    {
        // Arrange
        var sequences = new List<IReadOnlyList<char>>
        {
            "ACGT".ToList()
        };

        // Act
        double pi = PopulationGeneticsAnalyzer.CalculateNucleotideDiversity(sequences);

        // Assert
        Assert.That(pi, Is.EqualTo(0));
    }

    [Test]
    [Description("ND-M05: Empty input returns π = 0")]
    public void CalculateNucleotideDiversity_EmptyInput_ReturnsZero()
    {
        // Arrange
        var sequences = new List<IReadOnlyList<char>>();

        // Act
        double pi = PopulationGeneticsAnalyzer.CalculateNucleotideDiversity(sequences);

        // Assert
        Assert.That(pi, Is.EqualTo(0));
    }

    [Test]
    [Description("ND-M06: Nucleotide diversity is always non-negative")]
    public void CalculateNucleotideDiversity_VariousInputs_AlwaysNonNegative()
    {
        // Arrange - various test cases
        var testCases = new List<List<IReadOnlyList<char>>>
        {
            new() { "ACGT".ToList(), "ACGT".ToList() },
            new() { "AAAA".ToList(), "TTTT".ToList() },
            new() { "ACGT".ToList(), "ACGA".ToList(), "ACTT".ToList() },
            new() { "A".ToList() },
            new()
        };

        // Act & Assert
        foreach (var sequences in testCases)
        {
            double pi = PopulationGeneticsAnalyzer.CalculateNucleotideDiversity(sequences);
            Assert.That(pi, Is.GreaterThanOrEqualTo(0),
                $"π should be ≥ 0 for input with {sequences.Count} sequences");
        }
    }

    [Test]
    [Description("ND-S01: Partial polymorphism with 3+ sequences")]
    public void CalculateNucleotideDiversity_PartialPolymorphism_CalculatesCorrectly()
    {
        // Arrange - 3 sequences with some differences
        // Seq1: ACGT, Seq2: ACGA (1 diff), Seq3: ACGT (0 diff from Seq1, 1 from Seq2)
        // Pairwise: (1,2)=1, (1,3)=0, (2,3)=1 → total=2
        // π = 2 / (3 × 4) = 2/12 = 0.1667
        var sequences = new List<IReadOnlyList<char>>
        {
            "ACGT".ToList(),
            "ACGA".ToList(),
            "ACGT".ToList()
        };

        // Act
        double pi = PopulationGeneticsAnalyzer.CalculateNucleotideDiversity(sequences);

        // Assert
        Assert.That(pi, Is.EqualTo(2.0 / 12.0).Within(0.0001));
    }

    #endregion

    #region CalculateWattersonTheta Tests

    [Test]
    [Description("WT-M01: Wikipedia example - S=10, n=10, L=1000 → θ ≈ 0.00353")]
    public void CalculateWattersonTheta_WikipediaExample_CalculatesCorrectly()
    {
        // Arrange
        // a₁ = 1 + 1/2 + 1/3 + ... + 1/9 ≈ 2.8289...
        // θ = 10 / (2.8289 × 1000) ≈ 0.00353
        int segregatingSites = 10;
        int sampleSize = 10;
        int sequenceLength = 1000;

        // Act
        double theta = PopulationGeneticsAnalyzer.CalculateWattersonTheta(
            segregatingSites, sampleSize, sequenceLength);

        // Assert
        Assert.That(theta, Is.EqualTo(0.00353).Within(0.0005));
    }

    [Test]
    [Description("WT-M02: Zero segregating sites returns θ = 0")]
    public void CalculateWattersonTheta_ZeroSegregatingSites_ReturnsZero()
    {
        // Arrange
        int segregatingSites = 0;
        int sampleSize = 10;
        int sequenceLength = 1000;

        // Act
        double theta = PopulationGeneticsAnalyzer.CalculateWattersonTheta(
            segregatingSites, sampleSize, sequenceLength);

        // Assert
        Assert.That(theta, Is.EqualTo(0));
    }

    [Test]
    [Description("WT-M03: Sample size < 2 returns θ = 0")]
    public void CalculateWattersonTheta_SampleSizeLessThanTwo_ReturnsZero()
    {
        // Arrange & Act
        double theta1 = PopulationGeneticsAnalyzer.CalculateWattersonTheta(5, 1, 100);
        double theta0 = PopulationGeneticsAnalyzer.CalculateWattersonTheta(5, 0, 100);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(theta1, Is.EqualTo(0), "n=1 should return 0");
            Assert.That(theta0, Is.EqualTo(0), "n=0 should return 0");
        });
    }

    [Test]
    [Description("WT-M04: Sequence length ≤ 0 returns θ = 0")]
    public void CalculateWattersonTheta_InvalidSequenceLength_ReturnsZero()
    {
        // Arrange & Act
        double thetaZero = PopulationGeneticsAnalyzer.CalculateWattersonTheta(5, 10, 0);
        double thetaNeg = PopulationGeneticsAnalyzer.CalculateWattersonTheta(5, 10, -1);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(thetaZero, Is.EqualTo(0), "L=0 should return 0");
            Assert.That(thetaNeg, Is.EqualTo(0), "L<0 should return 0");
        });
    }

    [Test]
    [Description("WT-M05: Minimum valid sample size n=2 → θ = S / L")]
    public void CalculateWattersonTheta_MinimumSampleSize_CalculatesCorrectly()
    {
        // Arrange - for n=2, a₁ = 1
        // θ = S / (1 × L) = S / L
        int segregatingSites = 5;
        int sampleSize = 2;
        int sequenceLength = 100;

        // Act
        double theta = PopulationGeneticsAnalyzer.CalculateWattersonTheta(
            segregatingSites, sampleSize, sequenceLength);

        // Assert
        Assert.That(theta, Is.EqualTo(5.0 / 100.0).Within(0.0001));
    }

    [Test]
    [Description("WT-M06: Watterson's theta is always non-negative")]
    public void CalculateWattersonTheta_VariousInputs_AlwaysNonNegative()
    {
        // Arrange - various valid inputs
        var testCases = new[] { (5, 10, 100), (0, 5, 50), (10, 2, 200), (1, 3, 10) };

        // Act & Assert
        foreach (var (s, n, l) in testCases)
        {
            double theta = PopulationGeneticsAnalyzer.CalculateWattersonTheta(s, n, l);
            Assert.That(theta, Is.GreaterThanOrEqualTo(0),
                $"θ should be ≥ 0 for S={s}, n={n}, L={l}");
        }
    }

    [Test]
    [Description("WT-S01: Various sample sizes verify harmonic number calculation")]
    public void CalculateWattersonTheta_VariousSampleSizes_HarmonicNumberCorrect()
    {
        // Arrange - S=10, L=100, varying n
        // Expected a₁ values: n=3→1.5, n=5→2.083, n=10→2.829
        int s = 10;
        int l = 100;

        // Act
        double theta3 = PopulationGeneticsAnalyzer.CalculateWattersonTheta(s, 3, l);
        double theta5 = PopulationGeneticsAnalyzer.CalculateWattersonTheta(s, 5, l);
        double theta10 = PopulationGeneticsAnalyzer.CalculateWattersonTheta(s, 10, l);

        // Assert - θ decreases as n increases (a₁ increases)
        Assert.Multiple(() =>
        {
            Assert.That(theta3, Is.GreaterThan(theta5));
            Assert.That(theta5, Is.GreaterThan(theta10));
            Assert.That(theta3, Is.EqualTo(10.0 / (1.5 * 100)).Within(0.001)); // a₁(3) = 1.5
        });
    }

    #endregion

    #region CalculateTajimasD Tests

    [Test]
    [Description("TD-M01: When π ≈ θ (neutral evolution), D ≈ 0")]
    public void CalculateTajimasD_NeutralEvolution_NearZero()
    {
        // Arrange - equal diversity estimates suggest neutral evolution
        double nucleotideDiversity = 0.01;
        double wattersonTheta = 0.01;
        int segregatingSites = 100;
        int sampleSize = 50;

        // Act
        double d = PopulationGeneticsAnalyzer.CalculateTajimasD(
            nucleotideDiversity, wattersonTheta, segregatingSites, sampleSize);

        // Assert - D should be very close to 0
        Assert.That(d, Is.EqualTo(0).Within(0.5));
    }

    [Test]
    [Description("TD-M02: When π << θ (positive selection/expansion), D < 0")]
    public void CalculateTajimasD_PositiveSelection_Negative()
    {
        // Arrange - excess of rare variants (π << θ)
        double nucleotideDiversity = 0.001;
        double wattersonTheta = 0.01;
        int segregatingSites = 100;
        int sampleSize = 50;

        // Act
        double d = PopulationGeneticsAnalyzer.CalculateTajimasD(
            nucleotideDiversity, wattersonTheta, segregatingSites, sampleSize);

        // Assert
        Assert.That(d, Is.LessThan(0));
    }

    [Test]
    [Description("TD-M03: When π >> θ (balancing selection), D > 0")]
    public void CalculateTajimasD_BalancingSelection_Positive()
    {
        // Arrange - deficit of rare variants (π >> θ)
        double nucleotideDiversity = 0.02;
        double wattersonTheta = 0.005;
        int segregatingSites = 50;
        int sampleSize = 50;

        // Act
        double d = PopulationGeneticsAnalyzer.CalculateTajimasD(
            nucleotideDiversity, wattersonTheta, segregatingSites, sampleSize);

        // Assert
        Assert.That(d, Is.GreaterThan(0));
    }

    [Test]
    [Description("TD-M04: Zero segregating sites returns D = 0")]
    public void CalculateTajimasD_NoSegregatingSites_ReturnsZero()
    {
        // Arrange
        double nucleotideDiversity = 0;
        double wattersonTheta = 0;
        int segregatingSites = 0;
        int sampleSize = 50;

        // Act
        double d = PopulationGeneticsAnalyzer.CalculateTajimasD(
            nucleotideDiversity, wattersonTheta, segregatingSites, sampleSize);

        // Assert
        Assert.That(d, Is.EqualTo(0));
    }

    [Test]
    [Description("TD-M05: Sample size < 3 returns D = 0 (undefined)")]
    public void CalculateTajimasD_SampleSizeLessThanThree_ReturnsZero()
    {
        // Arrange & Act
        double d1 = PopulationGeneticsAnalyzer.CalculateTajimasD(0.01, 0.01, 10, 1);
        double d2 = PopulationGeneticsAnalyzer.CalculateTajimasD(0.01, 0.01, 10, 2);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(d1, Is.EqualTo(0), "n=1 should return 0");
            Assert.That(d2, Is.EqualTo(0), "n=2 should return 0");
        });
    }

    [Test]
    [Description("TD-M06: Minimum sample size n=3 produces valid D")]
    public void CalculateTajimasD_MinimumValidSampleSize_CalculatesValue()
    {
        // Arrange
        double nucleotideDiversity = 0.01;
        double wattersonTheta = 0.015;
        int segregatingSites = 10;
        int sampleSize = 3;

        // Act
        double d = PopulationGeneticsAnalyzer.CalculateTajimasD(
            nucleotideDiversity, wattersonTheta, segregatingSites, sampleSize);

        // Assert - should produce a valid (possibly non-zero) result
        Assert.That(double.IsNaN(d), Is.False);
        Assert.That(double.IsInfinity(d), Is.False);
    }

    #endregion

    #region CalculateDiversityStatistics Tests

    [Test]
    [Description("DS-M01: Returns all metrics with correct sample size")]
    public void CalculateDiversityStatistics_ReturnsAllMetrics()
    {
        // Arrange
        var sequences = new List<IReadOnlyList<char>>
        {
            "ACGTACGT".ToList(),
            "ACGTATGT".ToList(),
            "ACGTACGA".ToList()
        };

        // Act
        var stats = PopulationGeneticsAnalyzer.CalculateDiversityStatistics(sequences);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(stats.SampleSize, Is.EqualTo(3));
            Assert.That(stats.SegregratingSites, Is.GreaterThan(0));
            Assert.That(stats.NucleotideDiversity, Is.GreaterThanOrEqualTo(0));
            Assert.That(stats.WattersonTheta, Is.GreaterThanOrEqualTo(0));
            Assert.That(stats.HeterozygosityObserved, Is.GreaterThanOrEqualTo(0));
            Assert.That(stats.HeterozygosityExpected, Is.GreaterThanOrEqualTo(0));
        });
    }

    [Test]
    [Description("DS-M02: Single sequence returns all zeros")]
    public void CalculateDiversityStatistics_SingleSequence_ReturnsZeros()
    {
        // Arrange
        var sequences = new List<IReadOnlyList<char>>
        {
            "ACGT".ToList()
        };

        // Act
        var stats = PopulationGeneticsAnalyzer.CalculateDiversityStatistics(sequences);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(stats.NucleotideDiversity, Is.EqualTo(0));
            Assert.That(stats.WattersonTheta, Is.EqualTo(0));
            Assert.That(stats.TajimasD, Is.EqualTo(0));
            Assert.That(stats.SegregratingSites, Is.EqualTo(0));
            Assert.That(stats.SampleSize, Is.EqualTo(1));
        });
    }

    [Test]
    [Description("DS-M03: Empty input returns all zeros")]
    public void CalculateDiversityStatistics_EmptyInput_ReturnsZeros()
    {
        // Arrange
        var sequences = new List<IReadOnlyList<char>>();

        // Act
        var stats = PopulationGeneticsAnalyzer.CalculateDiversityStatistics(sequences);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(stats.NucleotideDiversity, Is.EqualTo(0));
            Assert.That(stats.WattersonTheta, Is.EqualTo(0));
            Assert.That(stats.TajimasD, Is.EqualTo(0));
            Assert.That(stats.SegregratingSites, Is.EqualTo(0));
            Assert.That(stats.SampleSize, Is.EqualTo(0));
        });
    }

    [Test]
    [Description("DS-M04: SampleSize matches input count")]
    public void CalculateDiversityStatistics_SampleSizeMatchesInputCount()
    {
        // Arrange
        var sequences = new List<IReadOnlyList<char>>
        {
            "ACGT".ToList(),
            "ACGA".ToList(),
            "TCGT".ToList(),
            "ACGT".ToList(),
            "ACCT".ToList()
        };

        // Act
        var stats = PopulationGeneticsAnalyzer.CalculateDiversityStatistics(sequences);

        // Assert
        Assert.That(stats.SampleSize, Is.EqualTo(5));
    }

    [Test]
    [Description("DS-M05: Segregating sites counted correctly")]
    public void CalculateDiversityStatistics_SegregatingSitesCountedCorrectly()
    {
        // Arrange - positions 0 and 3 are polymorphic
        // Pos 0: A,T (2 alleles), Pos 1-2: same, Pos 3: T,A (2 alleles)
        var sequences = new List<IReadOnlyList<char>>
        {
            "ACGT".ToList(),
            "TCGA".ToList()
        };

        // Act
        var stats = PopulationGeneticsAnalyzer.CalculateDiversityStatistics(sequences);

        // Assert - positions 0, 3 differ
        Assert.That(stats.SegregratingSites, Is.EqualTo(2));
    }

    [Test]
    [Description("DS-S01: Heterozygosity values are in valid range [0,1]")]
    public void CalculateDiversityStatistics_HeterozygosityInValidRange()
    {
        // Arrange
        var sequences = new List<IReadOnlyList<char>>
        {
            "ACGTACGT".ToList(),
            "TCGATCGA".ToList(),
            "GCTAGCTA".ToList(),
            "CGATCGAT".ToList()
        };

        // Act
        var stats = PopulationGeneticsAnalyzer.CalculateDiversityStatistics(sequences);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(stats.HeterozygosityObserved, Is.InRange(0, 1));
            Assert.That(stats.HeterozygosityExpected, Is.InRange(0, 1));
        });
    }

    [Test]
    [Description("DS-S02: Identical sequences have zero segregating sites")]
    public void CalculateDiversityStatistics_IdenticalSequences_ZeroSegregatingSites()
    {
        // Arrange
        var sequences = new List<IReadOnlyList<char>>
        {
            "ACGTACGT".ToList(),
            "ACGTACGT".ToList(),
            "ACGTACGT".ToList()
        };

        // Act
        var stats = PopulationGeneticsAnalyzer.CalculateDiversityStatistics(sequences);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(stats.SegregratingSites, Is.EqualTo(0));
            Assert.That(stats.NucleotideDiversity, Is.EqualTo(0));
            Assert.That(stats.WattersonTheta, Is.EqualTo(0));
        });
    }

    #endregion

    #region Integration Tests

    [Test]
    [Description("Comprehensive statistics for diverse population")]
    public void CalculateDiversityStatistics_DiversePopulation_AllMetricsReasonable()
    {
        // Arrange - moderately diverse population
        var sequences = new List<IReadOnlyList<char>>
        {
            "ACGTACGTACGT".ToList(),
            "ACGTACGTATGT".ToList(),
            "ACGTACGTACGA".ToList(),
            "ACGTACGTACGT".ToList(),
            "TCGTACGTACGT".ToList()
        };

        // Act
        var stats = PopulationGeneticsAnalyzer.CalculateDiversityStatistics(sequences);

        // Assert - sanity checks
        Assert.Multiple(() =>
        {
            Assert.That(stats.SampleSize, Is.EqualTo(5));
            Assert.That(stats.SegregratingSites, Is.GreaterThan(0));
            Assert.That(stats.NucleotideDiversity, Is.GreaterThan(0));
            Assert.That(stats.WattersonTheta, Is.GreaterThan(0));
            Assert.That(double.IsNaN(stats.TajimasD), Is.False);
            Assert.That(double.IsInfinity(stats.TajimasD), Is.False);
        });
    }

    #endregion
}
