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
    [Description("ND-M04: Wikipedia Tajima's D example dataset → k̂ = 2.0, π = 0.1")]
    public void CalculateNucleotideDiversity_WikipediaExample_CalculatesCorrectly()
    {
        // Arrange - Wikipedia Tajima's D example: 5 sequences, length 20
        // Source: https://en.wikipedia.org/wiki/Tajima%27s_D#Example
        // Person Y: 00000000000000000000
        // Person A: 00100000000010000010
        // Person B: 00000000000010000010
        // Person C: 00000010000000000010
        // Person D: 00000010000010000010
        // Using '0'/'1' characters to represent alleles
        var sequences = new List<IReadOnlyList<char>>
        {
            "00000000000000000000".ToList(), // Y
            "00100000000010000010".ToList(), // A
            "00000000000010000010".ToList(), // B
            "00000010000000000010".ToList(), // C
            "00000010000010000010".ToList(), // D
        };

        // Act
        double pi = PopulationGeneticsAnalyzer.CalculateNucleotideDiversity(sequences);

        // Assert
        // Pairwise differences (from Wikipedia):
        //   Y-A:3, Y-B:2, Y-C:2, Y-D:3, A-B:1, A-C:3, A-D:2, B-C:2, B-D:1, C-D:1
        //   Total = 20, Comparisons = C(5,2) = 10, k̂ = 20/10 = 2.0
        //   π = k̂ / L = 2.0 / 20 = 0.1
        Assert.That(pi, Is.EqualTo(0.1).Within(0.0001));
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

        // Assert — exact values from formula θ = S/(a₁×L):
        // a₁(3)  = 1 + 1/2 = 3/2          → θ = 10/(1.5×100)   = 0.06667
        // a₁(5)  = 1 + 1/2 + 1/3 + 1/4 = 25/12 → θ = 10/(2.08333×100) = 0.048
        // a₁(10) = Σ(1/i, i=1..9) ≈ 2.82897      → θ = 10/(2.82897×100) ≈ 0.03535
        Assert.Multiple(() =>
        {
            Assert.That(theta3, Is.EqualTo(10.0 / (1.5 * 100)).Within(0.0001), "a₁(3) = 1.5");
            Assert.That(theta5, Is.EqualTo(10.0 / (25.0 / 12.0 * 100)).Within(0.0001), "a₁(5) = 25/12");
            Assert.That(theta10, Is.EqualTo(0.03535).Within(0.0005), "a₁(10) ≈ 2.82897");
            Assert.That(theta3, Is.GreaterThan(theta5), "θ decreases as n increases");
            Assert.That(theta5, Is.GreaterThan(theta10), "θ decreases as n increases");
        });
    }

    #endregion

    #region CalculateTajimasD Tests

    [Test]
    [Description("TD-M01: When k\u0302 \u2248 S/a\u2081 (neutral evolution), D \u2248 0")]
    public void CalculateTajimasD_NeutralEvolution_NearZero()
    {
        // Arrange - when average pairwise differences match Watterson estimate
        // For n=50, a\u2081 \u2248 4.499, S=100 \u2192 S/a\u2081 \u2248 22.23
        // Set k\u0302 = S/a\u2081 to simulate neutral evolution
        double a1 = 0;
        for (int i = 1; i < 50; i++) a1 += 1.0 / i;
        double kHat = 100.0 / a1; // k\u0302 = S/a\u2081 (neutral expectation)
        int segregatingSites = 100;
        int sampleSize = 50;

        // Act
        double d = PopulationGeneticsAnalyzer.CalculateTajimasD(
            kHat, segregatingSites, sampleSize);

        // Assert - D should be exactly 0 when k\u0302 = S/a\u2081
        Assert.That(d, Is.EqualTo(0).Within(0.001));
    }

    [Test]
    [Description("TD-M02: When k\u0302 << S/a\u2081 (positive selection/expansion), D < 0")]
    public void CalculateTajimasD_PositiveSelection_Negative()
    {
        // Arrange - excess of rare variants (k\u0302 << S/a\u2081)
        // For n=50, a\u2081 \u2248 4.499, S=100 \u2192 S/a\u2081 \u2248 22.23
        // Set k\u0302 much lower than S/a\u2081
        double kHat = 5.0; // much less than S/a\u2081 \u2248 22.23
        int segregatingSites = 100;
        int sampleSize = 50;

        // Act
        double d = PopulationGeneticsAnalyzer.CalculateTajimasD(
            kHat, segregatingSites, sampleSize);

        // Assert
        Assert.That(d, Is.LessThan(0));
    }

    [Test]
    [Description("TD-M03: When k\u0302 >> S/a\u2081 (balancing selection), D > 0")]
    public void CalculateTajimasD_BalancingSelection_Positive()
    {
        // Arrange - deficit of rare variants (k\u0302 >> S/a\u2081)
        // For n=50, a\u2081 \u2248 4.499, S=50 \u2192 S/a\u2081 \u2248 11.11
        // Set k\u0302 much higher than S/a\u2081
        double kHat = 40.0; // much more than S/a\u2081 \u2248 11.11
        int segregatingSites = 50;
        int sampleSize = 50;

        // Act
        double d = PopulationGeneticsAnalyzer.CalculateTajimasD(
            kHat, segregatingSites, sampleSize);

        // Assert
        Assert.That(d, Is.GreaterThan(0));
    }

    [Test]
    [Description("TD-M04: Zero segregating sites returns D = 0")]
    public void CalculateTajimasD_NoSegregatingSites_ReturnsZero()
    {
        // Arrange
        double kHat = 0;
        int segregatingSites = 0;
        int sampleSize = 50;

        // Act
        double d = PopulationGeneticsAnalyzer.CalculateTajimasD(
            kHat, segregatingSites, sampleSize);

        // Assert
        Assert.That(d, Is.EqualTo(0));
    }

    [Test]
    [Description("TD-M05: Sample size < 3 returns D = 0 (undefined)")]
    public void CalculateTajimasD_SampleSizeLessThanThree_ReturnsZero()
    {
        // Arrange & Act
        double d1 = PopulationGeneticsAnalyzer.CalculateTajimasD(5.0, 10, 1);
        double d2 = PopulationGeneticsAnalyzer.CalculateTajimasD(5.0, 10, 2);

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
        double kHat = 5.0; // average pairwise differences
        int segregatingSites = 10;
        int sampleSize = 3;

        // Act
        double d = PopulationGeneticsAnalyzer.CalculateTajimasD(
            kHat, segregatingSites, sampleSize);

        // Assert - should produce a valid (possibly non-zero) result
        Assert.That(double.IsNaN(d), Is.False);
        Assert.That(double.IsInfinity(d), Is.False);
    }

    [Test]
    [Description("TD-C01: Wikipedia Tajima's D full example — exact numerical verification")]
    public void CalculateTajimasD_WikipediaExample_ExactValue()
    {
        // Source: https://en.wikipedia.org/wiki/Tajima%27s_D#Example
        // n=5, S=4, k̂=2.0
        //
        // a₁ = 1 + 1/2 + 1/3 + 1/4 = 25/12 ≈ 2.08333
        // a₂ = 1 + 1/4 + 1/9 + 1/16 = 205/144 ≈ 1.42361
        // S/a₁ = 4/2.08333 ≈ 1.92
        // d = k̂ − S/a₁ = 2.0 − 1.92 = 0.08
        //
        // b₁ = 6/12 = 0.5
        // b₂ = 66/180 ≈ 0.36667
        // c₁ = 0.5 − 1/2.08333 = 0.02
        // c₂ = 0.36667 − 7/10.41667 + 1.42361/4.34028 ≈ 0.02267
        // e₁ = 0.02/2.08333 ≈ 0.009600
        // e₂ = 0.02267/5.76389 ≈ 0.003933
        //
        // Var = e₁·S + e₂·S·(S−1) = 0.0384 + 0.04720 ≈ 0.08560
        // D = 0.08 / √0.08560 ≈ 0.2734
        double kHat = 2.0;
        int segregatingSites = 4;
        int sampleSize = 5;

        double d = PopulationGeneticsAnalyzer.CalculateTajimasD(kHat, segregatingSites, sampleSize);

        Assert.That(d, Is.EqualTo(0.273).Within(0.005),
            "Wikipedia example: D should be approximately 0.273");
    }

    [Test]
    [Description("TD-C02: Full end-to-end: Wikipedia sequences → CalculateDiversityStatistics → correct D")]
    public void CalculateDiversityStatistics_WikipediaExample_CorrectTajimasD()
    {
        // Source: https://en.wikipedia.org/wiki/Tajima%27s_D#Example
        // Same dataset as TD-C01 and ND-M04
        var sequences = new List<IReadOnlyList<char>>
        {
            "00000000000000000000".ToList(), // Y
            "00100000000010000010".ToList(), // A
            "00000000000010000010".ToList(), // B
            "00000010000000000010".ToList(), // C
            "00000010000010000010".ToList(), // D
        };

        var stats = PopulationGeneticsAnalyzer.CalculateDiversityStatistics(sequences);

        Assert.Multiple(() =>
        {
            Assert.That(stats.SampleSize, Is.EqualTo(5));
            Assert.That(stats.SegregratingSites, Is.EqualTo(4), "S = 4 (positions 3,7,13,19)");
            Assert.That(stats.NucleotideDiversity, Is.EqualTo(0.1).Within(0.001), "π = 2.0/20 = 0.1");
            Assert.That(stats.WattersonTheta, Is.EqualTo(0.096).Within(0.001), "θ_W = 4/(2.083×20) ≈ 0.096");
            Assert.That(stats.TajimasD, Is.EqualTo(0.273).Within(0.005), "D ≈ 0.273 per Wikipedia hand-calculation");
        });
    }

    #endregion

    #region CalculateDiversityStatistics Tests

    [Test]
    [Description("DS-M01: Returns all metrics with exact values")]
    public void CalculateDiversityStatistics_ReturnsAllMetrics()
    {
        // Arrange - 3 sequences, length 8
        // Seq1: ACGTACGT
        // Seq2: ACGTATGT (pos 5: C→T)
        // Seq3: ACGTACGA (pos 7: T→A)
        var sequences = new List<IReadOnlyList<char>>
        {
            "ACGTACGT".ToList(),
            "ACGTATGT".ToList(),
            "ACGTACGA".ToList()
        };

        // Act
        var stats = PopulationGeneticsAnalyzer.CalculateDiversityStatistics(sequences);

        // Assert — all values hand-computed from theory:
        // S = 2 (positions 5, 7)
        // Pairwise diffs: (1,2)=1, (1,3)=1, (2,3)=2, total=4, C(3,2)=3
        // k̂ = 4/3, π = 4/(3×8) = 1/6
        // a₁(3) = 1 + 1/2 = 3/2, θ_W = 2/(1.5×8) = 1/6
        // k̂ = S/a₁ = 4/3 = 4/3 → D = 0
        // H_exp = (4/9 + 4/9)/8 = 1/9, H_obs = 3/2 × 1/9 = 1/6
        Assert.Multiple(() =>
        {
            Assert.That(stats.SampleSize, Is.EqualTo(3));
            Assert.That(stats.SegregratingSites, Is.EqualTo(2), "S = 2 (positions 5, 7)");
            Assert.That(stats.NucleotideDiversity, Is.EqualTo(1.0 / 6.0).Within(0.0001), "π = 4/(3×8) = 1/6");
            Assert.That(stats.WattersonTheta, Is.EqualTo(1.0 / 6.0).Within(0.0001), "θ_W = 2/(1.5×8) = 1/6");
            Assert.That(stats.TajimasD, Is.EqualTo(0).Within(0.001), "k̂ = S/a₁ → D = 0");
            Assert.That(stats.HeterozygosityObserved, Is.EqualTo(1.0 / 6.0).Within(0.0001), "H_obs = n/(n-1) × H_exp = 1/6");
            Assert.That(stats.HeterozygosityExpected, Is.EqualTo(1.0 / 9.0).Within(0.0001), "H_exp = (8/9)/8 = 1/9");
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
}
