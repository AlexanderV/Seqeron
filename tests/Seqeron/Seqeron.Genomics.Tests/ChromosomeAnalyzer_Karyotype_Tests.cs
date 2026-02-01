using NUnit.Framework;
using Seqeron.Genomics;
using System.Collections.Generic;
using System.Linq;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Tests for ChromosomeAnalyzer karyotype analysis methods.
/// Test Unit: CHROM-KARYO-001
/// </summary>
[TestFixture]
public class ChromosomeAnalyzer_Karyotype_Tests
{
    #region AnalyzeKaryotype Tests

    [Test]
    [Description("Normal diploid karyotype with autosome pairs and sex chromosomes")]
    public void AnalyzeKaryotype_WithNormalDiploidSet_ReturnsCorrectCounts()
    {
        // Arrange - diploid set with 2 autosome pairs + XY sex chromosomes
        var chromosomes = new List<(string Name, long Length, bool IsSexChromosome)>
        {
            ("chr1_1", 248956422, false),
            ("chr1_2", 248956422, false),
            ("chr2_1", 242193529, false),
            ("chr2_2", 242193529, false),
            ("chrX", 156040895, true),
            ("chrY", 57227415, true)
        };

        // Act
        var result = ChromosomeAnalyzer.AnalyzeKaryotype(chromosomes);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(result.TotalChromosomes, Is.EqualTo(6));
            Assert.That(result.AutosomeCount, Is.EqualTo(4));
            Assert.That(result.SexChromosomes.Count, Is.EqualTo(2));
            Assert.That(result.HasAneuploidy, Is.False);
            Assert.That(result.Abnormalities, Is.Empty);
            Assert.That(result.PloidyLevel, Is.EqualTo(2));
        });
    }

    [Test]
    [Description("Trisomy (3 copies of a chromosome) is detected as aneuploidy")]
    public void AnalyzeKaryotype_WithTrisomy_DetectsAneuploidy()
    {
        // Arrange - Trisomy 21 (Down syndrome scenario)
        var chromosomes = new List<(string Name, long Length, bool IsSexChromosome)>
        {
            ("chr21_1", 46709983, false),
            ("chr21_2", 46709983, false),
            ("chr21_3", 46709983, false) // Extra copy = Trisomy
        };

        // Act
        var result = ChromosomeAnalyzer.AnalyzeKaryotype(chromosomes);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(result.HasAneuploidy, Is.True);
            Assert.That(result.Abnormalities, Has.Some.Contains("Trisomy"));
            Assert.That(result.TotalChromosomes, Is.EqualTo(3));
        });
    }

    [Test]
    [Description("Monosomy (1 copy of a chromosome) is detected as aneuploidy")]
    public void AnalyzeKaryotype_WithMonosomy_DetectsAneuploidy()
    {
        // Arrange - Monosomy (only one copy when diploid expected)
        var chromosomes = new List<(string Name, long Length, bool IsSexChromosome)>
        {
            ("chr1_1", 248956422, false) // Only one copy = Monosomy
        };

        // Act
        var result = ChromosomeAnalyzer.AnalyzeKaryotype(chromosomes);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(result.HasAneuploidy, Is.True);
            Assert.That(result.Abnormalities, Has.Some.Contains("Monosomy"));
        });
    }

    [Test]
    [Description("Empty input returns empty karyotype with no aneuploidy")]
    public void AnalyzeKaryotype_EmptyInput_ReturnsEmptyKaryotype()
    {
        // Act
        var result = ChromosomeAnalyzer.AnalyzeKaryotype(
            Enumerable.Empty<(string, long, bool)>());

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(result.TotalChromosomes, Is.EqualTo(0));
            Assert.That(result.AutosomeCount, Is.EqualTo(0));
            Assert.That(result.SexChromosomes, Is.Empty);
            Assert.That(result.HasAneuploidy, Is.False);
            Assert.That(result.Abnormalities, Is.Empty);
            Assert.That(result.TotalGenomeSize, Is.EqualTo(0));
        });
    }

    [Test]
    [Description("Sex chromosomes are correctly separated from autosomes")]
    public void AnalyzeKaryotype_CorrectlySeparatesSexChromosomes()
    {
        // Arrange
        var chromosomes = new List<(string Name, long Length, bool IsSexChromosome)>
        {
            ("chr1_1", 248956422, false),
            ("chr1_2", 248956422, false),
            ("chrX", 156040895, true),
            ("chrY", 57227415, true)
        };

        // Act
        var result = ChromosomeAnalyzer.AnalyzeKaryotype(chromosomes);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(result.SexChromosomes, Contains.Item("chrX"));
            Assert.That(result.SexChromosomes, Contains.Item("chrY"));
            Assert.That(result.AutosomeCount, Is.EqualTo(2));
        });
    }

    [Test]
    [Description("Total genome size is sum of all chromosome lengths")]
    public void AnalyzeKaryotype_CalculatesTotalGenomeSizeCorrectly()
    {
        // Arrange
        var chromosomes = new List<(string Name, long Length, bool IsSexChromosome)>
        {
            ("chr1_1", 100000, false),
            ("chr1_2", 100000, false),
            ("chr2_1", 50000, false),
            ("chr2_2", 50000, false)
        };
        long expectedTotalSize = 300000;

        // Act
        var result = ChromosomeAnalyzer.AnalyzeKaryotype(chromosomes);

        // Assert
        Assert.That(result.TotalGenomeSize, Is.EqualTo(expectedTotalSize));
    }

    [Test]
    [Description("Mean chromosome length is calculated correctly")]
    public void AnalyzeKaryotype_CalculatesMeanLengthCorrectly()
    {
        // Arrange
        var chromosomes = new List<(string Name, long Length, bool IsSexChromosome)>
        {
            ("chr1_1", 100000, false),
            ("chr1_2", 100000, false),
            ("chr2_1", 200000, false),
            ("chr2_2", 200000, false)
        };
        double expectedMean = 150000.0;

        // Act
        var result = ChromosomeAnalyzer.AnalyzeKaryotype(chromosomes);

        // Assert
        Assert.That(result.MeanChromosomeLength, Is.EqualTo(expectedMean));
    }

    [Test]
    [Description("Custom ploidy level (e.g., tetraploid) works correctly")]
    public void AnalyzeKaryotype_WithTetraploidExpectation_DetectsCorrectly()
    {
        // Arrange - 4 copies of each chromosome (tetraploid)
        var chromosomes = new List<(string Name, long Length, bool IsSexChromosome)>
        {
            ("chr1_1", 100000, false),
            ("chr1_2", 100000, false),
            ("chr1_3", 100000, false),
            ("chr1_4", 100000, false)
        };

        // Act - expect tetraploid (4 copies)
        var result = ChromosomeAnalyzer.AnalyzeKaryotype(chromosomes, expectedPloidyLevel: 4);

        // Assert
        Assert.That(result.HasAneuploidy, Is.False);
    }

    [Test]
    [Description("Multiple chromosome groups with mixed aneuploidy")]
    public void AnalyzeKaryotype_WithMultipleAneuploidies_DetectsAll()
    {
        // Arrange - Trisomy chr21, Monosomy chr22
        var chromosomes = new List<(string Name, long Length, bool IsSexChromosome)>
        {
            ("chr21_1", 46709983, false),
            ("chr21_2", 46709983, false),
            ("chr21_3", 46709983, false), // Trisomy
            ("chr22_1", 50818468, false)  // Monosomy
        };

        // Act
        var result = ChromosomeAnalyzer.AnalyzeKaryotype(chromosomes);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(result.HasAneuploidy, Is.True);
            Assert.That(result.Abnormalities.Count, Is.EqualTo(2));
            Assert.That(result.Abnormalities, Has.Some.Contains("Trisomy"));
            Assert.That(result.Abnormalities, Has.Some.Contains("Monosomy"));
        });
    }

    [Test]
    [Description("Invariant: TotalChromosomes = AutosomeCount + SexChromosomes.Count")]
    public void AnalyzeKaryotype_Invariant_TotalEqualsAutosomesPlusSexChromosomes()
    {
        // Arrange
        var chromosomes = new List<(string Name, long Length, bool IsSexChromosome)>
        {
            ("chr1_1", 248956422, false),
            ("chr1_2", 248956422, false),
            ("chr2_1", 242193529, false),
            ("chr2_2", 242193529, false),
            ("chr3_1", 198295559, false),
            ("chr3_2", 198295559, false),
            ("chrX", 156040895, true),
            ("chrY", 57227415, true)
        };

        // Act
        var result = ChromosomeAnalyzer.AnalyzeKaryotype(chromosomes);

        // Assert - Invariant check
        Assert.That(result.TotalChromosomes,
            Is.EqualTo(result.AutosomeCount + result.SexChromosomes.Count),
            "TotalChromosomes must equal AutosomeCount + SexChromosomes.Count");
    }

    [Test]
    [Description("Invariant: HasAneuploidy correlates with non-empty Abnormalities")]
    public void AnalyzeKaryotype_Invariant_AneuploidyCorrelatesWithAbnormalities()
    {
        // Arrange - Normal karyotype
        var normalChromosomes = new List<(string Name, long Length, bool IsSexChromosome)>
        {
            ("chr1_1", 100000, false),
            ("chr1_2", 100000, false)
        };

        // Arrange - Abnormal karyotype
        var abnormalChromosomes = new List<(string Name, long Length, bool IsSexChromosome)>
        {
            ("chr1_1", 100000, false),
            ("chr1_2", 100000, false),
            ("chr1_3", 100000, false) // Trisomy
        };

        // Act
        var normalResult = ChromosomeAnalyzer.AnalyzeKaryotype(normalChromosomes);
        var abnormalResult = ChromosomeAnalyzer.AnalyzeKaryotype(abnormalChromosomes);

        // Assert
        Assert.Multiple(() =>
        {
            // Normal: no aneuploidy, empty abnormalities
            Assert.That(normalResult.HasAneuploidy, Is.False);
            Assert.That(normalResult.Abnormalities, Is.Empty);

            // Abnormal: has aneuploidy, non-empty abnormalities
            Assert.That(abnormalResult.HasAneuploidy, Is.True);
            Assert.That(abnormalResult.Abnormalities, Is.Not.Empty);
        });
    }

    #endregion

    #region DetectPloidy Tests

    [Test]
    [Description("Diploid depth (ratio ≈ 1.0) returns ploidy 2")]
    public void DetectPloidy_WithDiploidDepth_ReturnsPloidy2()
    {
        // Arrange - Uniform depth at expected diploid level
        var depths = Enumerable.Repeat(1.0, 100);

        // Act
        var (ploidy, confidence) = ChromosomeAnalyzer.DetectPloidy(depths);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(ploidy, Is.EqualTo(2));
            Assert.That(confidence, Is.GreaterThan(0.9));
        });
    }

    [Test]
    [Description("Tetraploid depth (ratio ≈ 2.0) returns ploidy 4")]
    public void DetectPloidy_WithTetraploidDepth_ReturnsPloidy4()
    {
        // Arrange - Depth at 2x expected diploid level
        var depths = Enumerable.Repeat(2.0, 100);

        // Act
        var (ploidy, confidence) = ChromosomeAnalyzer.DetectPloidy(depths);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(ploidy, Is.EqualTo(4));
            Assert.That(confidence, Is.GreaterThan(0.9));
        });
    }

    [Test]
    [Description("Haploid depth (ratio ≈ 0.5) returns ploidy 1")]
    public void DetectPloidy_WithHaploidDepth_ReturnsPloidy1()
    {
        // Arrange - Depth at 0.5x expected diploid level
        var depths = Enumerable.Repeat(0.5, 100);

        // Act
        var (ploidy, confidence) = ChromosomeAnalyzer.DetectPloidy(depths);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(ploidy, Is.EqualTo(1));
            Assert.That(confidence, Is.GreaterThan(0.9));
        });
    }

    [Test]
    [Description("Triploid depth (ratio ≈ 1.5) returns ploidy 3")]
    public void DetectPloidy_WithTriploidDepth_ReturnsPloidy3()
    {
        // Arrange - Depth at 1.5x expected diploid level
        var depths = Enumerable.Repeat(1.5, 100);

        // Act
        var (ploidy, confidence) = ChromosomeAnalyzer.DetectPloidy(depths);

        // Assert
        Assert.That(ploidy, Is.EqualTo(3));
    }

    [Test]
    [Description("Empty input returns default diploid with zero confidence")]
    public void DetectPloidy_EmptyInput_ReturnsDefault()
    {
        // Act
        var (ploidy, confidence) = ChromosomeAnalyzer.DetectPloidy(
            Enumerable.Empty<double>());

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(ploidy, Is.EqualTo(2), "Default ploidy should be diploid (2)");
            Assert.That(confidence, Is.EqualTo(0), "Confidence should be 0 for empty input");
        });
    }

    [Test]
    [Description("High ploidy values are clamped to maximum of 8")]
    public void DetectPloidy_ExtremeHighDepth_ClampedToMax()
    {
        // Arrange - Very high depth (would suggest ploidy > 8)
        var depths = Enumerable.Repeat(10.0, 100); // 10x depth = would be 20-ploid

        // Act
        var (ploidy, _) = ChromosomeAnalyzer.DetectPloidy(depths);

        // Assert
        Assert.That(ploidy, Is.LessThanOrEqualTo(8), "Ploidy should be clamped to 8");
    }

    [Test]
    [Description("Very low ploidy values are clamped to minimum of 1")]
    public void DetectPloidy_ExtremeLowDepth_ClampedToMin()
    {
        // Arrange - Very low depth
        var depths = Enumerable.Repeat(0.1, 100);

        // Act
        var (ploidy, _) = ChromosomeAnalyzer.DetectPloidy(depths);

        // Assert
        Assert.That(ploidy, Is.GreaterThanOrEqualTo(1), "Ploidy should be at least 1");
    }

    [Test]
    [Description("Noisy data reduces confidence")]
    public void DetectPloidy_WithNoisyData_ReducedConfidence()
    {
        // Arrange - Mix of different depths (noisy)
        var depths = new[] { 0.8, 1.0, 1.2, 0.9, 1.1, 1.0, 0.85, 1.15 };

        // Arrange - Uniform data
        var uniformDepths = Enumerable.Repeat(1.0, 100);

        // Act
        var (_, noisyConfidence) = ChromosomeAnalyzer.DetectPloidy(depths);
        var (_, uniformConfidence) = ChromosomeAnalyzer.DetectPloidy(uniformDepths);

        // Assert - Both should be high confidence for diploid
        Assert.That(noisyConfidence, Is.GreaterThan(0.5));
        Assert.That(uniformConfidence, Is.GreaterThan(noisyConfidence).Or.EqualTo(noisyConfidence));
    }

    [Test]
    [Description("Invariant: PloidyLevel is always in range [1, 8]")]
    [TestCase(0.1)]
    [TestCase(0.5)]
    [TestCase(1.0)]
    [TestCase(2.0)]
    [TestCase(4.0)]
    [TestCase(10.0)]
    public void DetectPloidy_Invariant_PloidyAlwaysInValidRange(double depthValue)
    {
        // Arrange
        var depths = Enumerable.Repeat(depthValue, 100);

        // Act
        var (ploidy, _) = ChromosomeAnalyzer.DetectPloidy(depths);

        // Assert
        Assert.That(ploidy, Is.InRange(1, 8));
    }

    [Test]
    [Description("Invariant: Confidence is always in range [0, 1]")]
    [TestCase(0.1)]
    [TestCase(0.5)]
    [TestCase(1.0)]
    [TestCase(1.25)]
    [TestCase(2.0)]
    public void DetectPloidy_Invariant_ConfidenceAlwaysInValidRange(double depthValue)
    {
        // Arrange
        var depths = Enumerable.Repeat(depthValue, 100);

        // Act
        var (_, confidence) = ChromosomeAnalyzer.DetectPloidy(depths);

        // Assert
        Assert.That(confidence, Is.InRange(0, 1));
    }

    [Test]
    [Description("Single value input works correctly")]
    public void DetectPloidy_SingleValue_WorksCorrectly()
    {
        // Arrange
        var depths = new[] { 1.0 };

        // Act
        var (ploidy, confidence) = ChromosomeAnalyzer.DetectPloidy(depths);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(ploidy, Is.EqualTo(2));
            Assert.That(confidence, Is.GreaterThanOrEqualTo(0));
        });
    }

    [Test]
    [Description("Custom expected diploid depth works correctly")]
    public void DetectPloidy_WithCustomExpectedDepth_WorksCorrectly()
    {
        // Arrange - Depths at 60x when expecting 30x for diploid
        var depths = Enumerable.Repeat(60.0, 100);
        double expectedDiploidDepth = 30.0;

        // Act
        var (ploidy, confidence) = ChromosomeAnalyzer.DetectPloidy(depths, expectedDiploidDepth);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(ploidy, Is.EqualTo(4)); // 60/30 = 2.0 ratio → tetraploid
            Assert.That(confidence, Is.GreaterThan(0.9));
        });
    }

    #endregion
}
