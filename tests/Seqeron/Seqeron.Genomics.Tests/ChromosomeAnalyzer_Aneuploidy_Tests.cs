using NUnit.Framework;
using Seqeron.Genomics;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Tests for ChromosomeAnalyzer aneuploidy detection methods.
/// Test Unit: CHROM-ANEU-001
/// 
/// Evidence: Wikipedia (Aneuploidy, Copy Number Variation), Griffiths et al. (2000)
/// </summary>
[TestFixture]
public class ChromosomeAnalyzer_Aneuploidy_Tests
{
    private const double MedianDepth = 30.0;
    private const int DefaultBinSize = 1_000_000;

    #region DetectAneuploidy - Core Functionality

    [Test]
    [Description("Normal diploid depth returns copy number 2")]
    public void DetectAneuploidy_NormalDepth_ReturnsCopyNumber2()
    {
        // Arrange: depth = median (ratio = 1.0) → CN = 2
        var depthData = Enumerable.Range(0, 100)
            .Select(i => ("chr1", i * DefaultBinSize, MedianDepth));

        // Act
        var results = ChromosomeAnalyzer.DetectAneuploidy(
            depthData, medianDepth: MedianDepth, binSize: DefaultBinSize).ToList();

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(results, Is.Not.Empty);
            Assert.That(results.All(r => r.CopyNumber == 2), Is.True,
                "All bins should have copy number 2 for normal depth");
            Assert.That(results[0].LogRatio, Is.EqualTo(0.0).Within(1e-10),
                "log2(1.0) = 0.0");
            Assert.That(results[0].Confidence, Is.EqualTo(1.0).Within(1e-10),
                "Exact ratio match yields confidence 1.0");
        });
    }

    [Test]
    [Description("Trisomy detection: 1.5× depth → copy number 3 (Wikipedia: Down syndrome chr21)")]
    public void DetectAneuploidy_Trisomy_ReturnsCopyNumber3()
    {
        // Arrange: depth = 1.5× median → CN = 3
        var depthData = Enumerable.Range(0, 100)
            .Select(i => ("chr21", i * DefaultBinSize, MedianDepth * 1.5));

        // Act
        var results = ChromosomeAnalyzer.DetectAneuploidy(
            depthData, medianDepth: MedianDepth, binSize: DefaultBinSize).ToList();

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(results, Is.Not.Empty);
            Assert.That(results.All(r => r.CopyNumber == 3), Is.True,
                "1.5× depth should yield copy number 3 (trisomy)");
            Assert.That(results[0].LogRatio, Is.EqualTo(Math.Log2(1.5)).Within(1e-10),
                "log2(1.5) ≈ 0.585");
            Assert.That(results[0].Confidence, Is.EqualTo(1.0).Within(1e-10),
                "Exact ratio match yields confidence 1.0");
        });
    }

    [Test]
    [Description("Monosomy detection: 0.5× depth → copy number 1 (Wikipedia: Turner syndrome)")]
    public void DetectAneuploidy_Monosomy_ReturnsCopyNumber1()
    {
        // Arrange: depth = 0.5× median → CN = 1
        var depthData = Enumerable.Range(0, 100)
            .Select(i => ("chrX", i * DefaultBinSize, MedianDepth * 0.5));

        // Act
        var results = ChromosomeAnalyzer.DetectAneuploidy(
            depthData, medianDepth: MedianDepth, binSize: DefaultBinSize).ToList();

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(results, Is.Not.Empty);
            Assert.That(results.All(r => r.CopyNumber == 1), Is.True,
                "0.5× depth should yield copy number 1 (monosomy)");
            Assert.That(results[0].LogRatio, Is.EqualTo(-1.0).Within(1e-10),
                "log2(0.5) = -1.0");
            Assert.That(results[0].Confidence, Is.EqualTo(1.0).Within(1e-10),
                "Exact ratio match yields confidence 1.0");
        });
    }

    [Test]
    [Description("Nullisomy detection: 0× depth → copy number 0 (Wikipedia)")]
    public void DetectAneuploidy_Nullisomy_ReturnsCopyNumber0()
    {
        // Arrange: depth = 0 → CN = 0
        var depthData = Enumerable.Range(0, 10)
            .Select(i => ("chr1", i * DefaultBinSize, 0.0));

        // Act
        var results = ChromosomeAnalyzer.DetectAneuploidy(
            depthData, medianDepth: MedianDepth, binSize: DefaultBinSize).ToList();

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(results, Is.Not.Empty);
            Assert.That(results.All(r => r.CopyNumber == 0), Is.True,
                "Zero depth should yield copy number 0 (nullisomy)");
            Assert.That(results[0].LogRatio, Is.EqualTo(double.NegativeInfinity),
                "log2(0) = -∞");
            Assert.That(results[0].Confidence, Is.EqualTo(1.0).Within(1e-10),
                "Nullisomy boundary yields confidence 1.0");
        });
    }

    [Test]
    [Description("Tetrasomy detection: 2× depth → copy number 4 (Wikipedia)")]
    public void DetectAneuploidy_Tetrasomy_ReturnsCopyNumber4()
    {
        // Arrange: depth = 2× median → CN = 4
        var depthData = Enumerable.Range(0, 50)
            .Select(i => ("chr1", i * DefaultBinSize, MedianDepth * 2.0));

        // Act
        var results = ChromosomeAnalyzer.DetectAneuploidy(
            depthData, medianDepth: MedianDepth, binSize: DefaultBinSize).ToList();

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(results, Is.Not.Empty);
            Assert.That(results.All(r => r.CopyNumber == 4), Is.True,
                "2× depth should yield copy number 4 (tetrasomy)");
            Assert.That(results[0].LogRatio, Is.EqualTo(1.0).Within(1e-10),
                "log2(2.0) = 1.0");
            Assert.That(results[0].Confidence, Is.EqualTo(1.0).Within(1e-10),
                "Exact ratio match yields confidence 1.0");
        });
    }

    #endregion

    #region DetectAneuploidy - Edge Cases

    [Test]
    [Description("Empty input returns empty result")]
    public void DetectAneuploidy_EmptyInput_ReturnsEmpty()
    {
        // Act
        var results = ChromosomeAnalyzer.DetectAneuploidy(
            Enumerable.Empty<(string, int, double)>(), MedianDepth).ToList();

        // Assert
        Assert.That(results, Is.Empty);
    }

    [Test]
    [Description("Zero median depth returns empty to prevent division by zero")]
    public void DetectAneuploidy_ZeroMedianDepth_ReturnsEmpty()
    {
        // Arrange
        var depthData = new[] { ("chr1", 0, 30.0) };

        // Act
        var results = ChromosomeAnalyzer.DetectAneuploidy(
            depthData, medianDepth: 0.0).ToList();

        // Assert
        Assert.That(results, Is.Empty);
    }

    [Test]
    [Description("Negative median depth returns empty")]
    public void DetectAneuploidy_NegativeMedianDepth_ReturnsEmpty()
    {
        // Arrange
        var depthData = new[] { ("chr1", 0, 30.0) };

        // Act
        var results = ChromosomeAnalyzer.DetectAneuploidy(
            depthData, medianDepth: -10.0).ToList();

        // Assert
        Assert.That(results, Is.Empty);
    }

    [Test]
    [Description("Copy number clamped to maximum of 10")]
    public void DetectAneuploidy_VeryHighDepth_CopyNumberClampedTo10()
    {
        // Arrange: depth = 10× median would give CN = 20, but should clamp to 10
        var depthData = new[] { ("chr1", 0, MedianDepth * 10.0) };

        // Act
        var results = ChromosomeAnalyzer.DetectAneuploidy(
            depthData, medianDepth: MedianDepth).ToList();

        // Assert
        Assert.That(results.Single().CopyNumber, Is.EqualTo(10),
            "Copy number should be clamped to 10");
    }

    [Test]
    [Description("Copy number minimum is 0")]
    public void DetectAneuploidy_VeryLowDepth_CopyNumberMinimumIsZero()
    {
        // Arrange: very low but non-zero depth
        var depthData = new[] { ("chr1", 0, MedianDepth * 0.01) };

        // Act
        var results = ChromosomeAnalyzer.DetectAneuploidy(
            depthData, medianDepth: MedianDepth).ToList();

        // Assert
        Assert.That(results.Single().CopyNumber, Is.EqualTo(0),
            "Depth ratio 0.01 \u2192 round(0.02) = 0");
    }

    #endregion

    #region DetectAneuploidy - Grouping and Binning

    [Test]
    [Description("Multiple chromosomes are grouped correctly")]
    public void DetectAneuploidy_MultipleChromosomes_GroupsCorrectly()
    {
        // Arrange: chr1 normal, chr21 trisomy
        var depthData = new List<(string, int, double)>();
        for (int i = 0; i < 10; i++)
        {
            depthData.Add(("chr1", i * DefaultBinSize, MedianDepth));      // Normal
            depthData.Add(("chr21", i * DefaultBinSize, MedianDepth * 1.5)); // Trisomy
        }

        // Act
        var results = ChromosomeAnalyzer.DetectAneuploidy(
            depthData, medianDepth: MedianDepth, binSize: DefaultBinSize).ToList();

        var chr1Results = results.Where(r => r.Chromosome == "chr1").ToList();
        var chr21Results = results.Where(r => r.Chromosome == "chr21").ToList();

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(chr1Results.All(r => r.CopyNumber == 2), Is.True,
                "chr1 should have CN=2");
            Assert.That(chr21Results.All(r => r.CopyNumber == 3), Is.True,
                "chr21 should have CN=3");
        });
    }

    [Test]
    [Description("Binning aggregates multiple depth values by averaging")]
    public void DetectAneuploidy_BinAggregation_AveragesDepth()
    {
        // Arrange: Multiple points in same bin with varying depth, average = median
        var depthData = new[]
        {
            ("chr1", 100, 20.0),   // Same bin (position / binSize = 0)
            ("chr1", 200, 40.0),   // Same bin
            ("chr1", 300, 30.0)    // Same bin - average = 30 = median
        };

        // Act
        var results = ChromosomeAnalyzer.DetectAneuploidy(
            depthData, medianDepth: MedianDepth, binSize: DefaultBinSize).ToList();

        // Assert
        Assert.That(results.Count, Is.EqualTo(1), "Should produce single bin");
        Assert.That(results.Single().CopyNumber, Is.EqualTo(2),
            "Average depth of 30 should yield CN=2");
    }

    [Test]
    [Description("Single data point works correctly")]
    public void DetectAneuploidy_SingleDataPoint_WorksCorrectly()
    {
        // Arrange
        var depthData = new[] { ("chr1", 0, MedianDepth) };

        // Act
        var results = ChromosomeAnalyzer.DetectAneuploidy(
            depthData, medianDepth: MedianDepth).ToList();

        // Assert
        Assert.That(results.Count, Is.EqualTo(1));
        Assert.That(results.Single().CopyNumber, Is.EqualTo(2));
    }

    [Test]
    [Description("Larger bin size aggregates more data points, producing fewer output bins")]
    public void DetectAneuploidy_LargeBinSize_ReducesOutput()
    {
        // Arrange: 10 data points at standard 1 Mb spacing
        var depthData = Enumerable.Range(0, 10)
            .Select(i => ("chr1", i * DefaultBinSize, MedianDepth));

        // Act: Compare output counts with different bin sizes
        var smallBinResults = ChromosomeAnalyzer.DetectAneuploidy(
            depthData, MedianDepth, binSize: DefaultBinSize).ToList();
        var largeBinResults = ChromosomeAnalyzer.DetectAneuploidy(
            depthData, MedianDepth, binSize: 5 * DefaultBinSize).ToList();

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(smallBinResults, Has.Count.EqualTo(10),
                "1 Mb bins → 10 output bins");
            Assert.That(largeBinResults, Has.Count.EqualTo(2),
                "5 Mb bins → 2 output bins");
        });
    }

    #endregion

    #region DetectAneuploidy - Output Invariants

    [Test]
    [Description("Confidence equals 1.0 for exact copy-number boundary ratios")]
    public void DetectAneuploidy_Confidence_ExactForBoundaryRatios()
    {
        // Arrange: All depth values are exact CN boundaries
        var depthData = new[]
        {
            ("chr1", 0, 0.0),           // Nullisomy
            ("chr1", 1000000, 15.0),    // Monosomy
            ("chr1", 2000000, 30.0),    // Normal
            ("chr1", 3000000, 45.0),    // Trisomy
            ("chr1", 4000000, 60.0)     // Tetrasomy
        };

        // Act
        var results = ChromosomeAnalyzer.DetectAneuploidy(
            depthData, medianDepth: MedianDepth).ToList();

        // Assert: exact boundary ratios yield confidence = 1.0
        Assert.Multiple(() =>
        {
            foreach (var r in results)
            {
                Assert.That(r.Confidence, Is.EqualTo(1.0).Within(1e-10),
                    $"Exact CN boundary (CN={r.CopyNumber}) should give confidence \u2248 1.0");
            }
        });
    }

    [Test]
    [Description("Start is always less than End for each bin")]
    public void DetectAneuploidy_OutputPositions_StartLessThanEnd()
    {
        // Arrange
        var depthData = Enumerable.Range(0, 10)
            .Select(i => ("chr1", i * DefaultBinSize, MedianDepth));

        // Act
        var results = ChromosomeAnalyzer.DetectAneuploidy(
            depthData, medianDepth: MedianDepth, binSize: DefaultBinSize).ToList();

        // Assert
        Assert.That(results.All(r => r.Start < r.End), Is.True,
            "Start should be less than End for all bins");
    }

    [Test]
    [Description("LogRatio matches formula log2(depth/medianDepth)")]
    public void DetectAneuploidy_LogRatio_MatchesFormula()
    {
        // Arrange: Known depth ratios with exact log2 values
        var depthData = new[]
        {
            ("chr1", 0 * DefaultBinSize, MedianDepth * 0.5),  // log2(0.5) = -1.0
            ("chr1", 1 * DefaultBinSize, MedianDepth * 1.0),  // log2(1.0) =  0.0
            ("chr1", 2 * DefaultBinSize, MedianDepth * 2.0),  // log2(2.0) = +1.0
        };

        // Act
        var results = ChromosomeAnalyzer.DetectAneuploidy(
            depthData, medianDepth: MedianDepth, binSize: DefaultBinSize).ToList();

        // Assert: Verify logRatio = log2(depth / medianDepth)
        Assert.Multiple(() =>
        {
            Assert.That(results[0].LogRatio, Is.EqualTo(-1.0).Within(1e-10),
                "log2(0.5) should be -1.0");
            Assert.That(results[1].LogRatio, Is.EqualTo(0.0).Within(1e-10),
                "log2(1.0) should be 0.0");
            Assert.That(results[2].LogRatio, Is.EqualTo(1.0).Within(1e-10),
                "log2(2.0) should be +1.0");
        });
    }

    [Test]
    [Description("Output bins are ordered by position within each chromosome")]
    public void DetectAneuploidy_OutputOrdered_ByPosition()
    {
        // Arrange: Data in reverse order
        var depthData = new[]
        {
            ("chr1", 4 * DefaultBinSize, MedianDepth),
            ("chr1", 1 * DefaultBinSize, MedianDepth),
            ("chr1", 3 * DefaultBinSize, MedianDepth),
            ("chr1", 0 * DefaultBinSize, MedianDepth),
            ("chr1", 2 * DefaultBinSize, MedianDepth),
        };

        // Act
        var results = ChromosomeAnalyzer.DetectAneuploidy(
            depthData, medianDepth: MedianDepth, binSize: DefaultBinSize).ToList();

        // Assert: Bins should be ordered by Start position
        var starts = results.Select(r => r.Start).ToList();
        Assert.That(starts, Is.Ordered,
            "Output bins should be ordered by position");
    }

    #endregion

    #region IdentifyWholeChromosomeAneuploidy - Classifications

    [Test]
    [Description("Trisomy classification for CN=3 chromosome (Wikipedia: Down syndrome)")]
    public void IdentifyWholeChromosomeAneuploidy_Trisomy_IdentifiesTrisomy()
    {
        // Arrange
        var states = Enumerable.Range(0, 100)
            .Select(i => new ChromosomeAnalyzer.CopyNumberState(
                "chr21", i * DefaultBinSize, (i + 1) * DefaultBinSize - 1,
                CopyNumber: 3, LogRatio: 0.58, Confidence: 0.9));

        // Act
        var results = ChromosomeAnalyzer.IdentifyWholeChromosomeAneuploidy(states).ToList();

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(results.Count, Is.EqualTo(1));
            Assert.That(results[0].Chromosome, Is.EqualTo("chr21"));
            Assert.That(results[0].CopyNumber, Is.EqualTo(3));
            Assert.That(results[0].Type, Is.EqualTo("Trisomy"));
        });
    }

    [Test]
    [Description("Monosomy classification for CN=1 chromosome (Wikipedia: Turner syndrome)")]
    public void IdentifyWholeChromosomeAneuploidy_Monosomy_IdentifiesMonosomy()
    {
        // Arrange
        var states = Enumerable.Range(0, 100)
            .Select(i => new ChromosomeAnalyzer.CopyNumberState(
                "chrX", i * DefaultBinSize, (i + 1) * DefaultBinSize - 1,
                CopyNumber: 1, LogRatio: -1.0, Confidence: 0.9));

        // Act
        var results = ChromosomeAnalyzer.IdentifyWholeChromosomeAneuploidy(states).ToList();

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(results.Count, Is.EqualTo(1));
            Assert.That(results[0].Chromosome, Is.EqualTo("chrX"));
            Assert.That(results[0].CopyNumber, Is.EqualTo(1));
            Assert.That(results[0].Type, Is.EqualTo("Monosomy"));
        });
    }

    [Test]
    [Description("Nullisomy classification for CN=0 chromosome")]
    public void IdentifyWholeChromosomeAneuploidy_Nullisomy_IdentifiesNullisomy()
    {
        // Arrange
        var states = Enumerable.Range(0, 100)
            .Select(i => new ChromosomeAnalyzer.CopyNumberState(
                "chr1", i * DefaultBinSize, (i + 1) * DefaultBinSize - 1,
                CopyNumber: 0, LogRatio: double.NegativeInfinity, Confidence: 0.9));

        // Act
        var results = ChromosomeAnalyzer.IdentifyWholeChromosomeAneuploidy(states).ToList();

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(results.Count, Is.EqualTo(1));
            Assert.That(results[0].CopyNumber, Is.EqualTo(0));
            Assert.That(results[0].Type, Is.EqualTo("Nullisomy"));
        });
    }

    [Test]
    [Description("Tetrasomy classification for CN=4 chromosome")]
    public void IdentifyWholeChromosomeAneuploidy_Tetrasomy_IdentifiesTetrasomy()
    {
        // Arrange
        var states = Enumerable.Range(0, 100)
            .Select(i => new ChromosomeAnalyzer.CopyNumberState(
                "chr1", i * DefaultBinSize, (i + 1) * DefaultBinSize - 1,
                CopyNumber: 4, LogRatio: 1.0, Confidence: 0.9));

        // Act
        var results = ChromosomeAnalyzer.IdentifyWholeChromosomeAneuploidy(states).ToList();

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(results.Count, Is.EqualTo(1));
            Assert.That(results[0].CopyNumber, Is.EqualTo(4));
            Assert.That(results[0].Type, Is.EqualTo("Tetrasomy"));
        });
    }

    [Test]
    [Description("Normal chromosome (CN=2) is not flagged as aneuploidy")]
    public void IdentifyWholeChromosomeAneuploidy_Normal_ReturnsEmpty()
    {
        // Arrange
        var states = Enumerable.Range(0, 100)
            .Select(i => new ChromosomeAnalyzer.CopyNumberState(
                "chr1", i * DefaultBinSize, (i + 1) * DefaultBinSize - 1,
                CopyNumber: 2, LogRatio: 0.0, Confidence: 0.95));

        // Act
        var results = ChromosomeAnalyzer.IdentifyWholeChromosomeAneuploidy(states).ToList();

        // Assert
        Assert.That(results, Is.Empty,
            "Chromosome with CN=2 should not be identified as aneuploid");
    }

    [Test]
    [Description("Pentasomy classification for CN=5 chromosome (Wikipedia: Tetrasomy and pentasomy)")]
    public void IdentifyWholeChromosomeAneuploidy_Pentasomy_IdentifiesPentasomy()
    {
        // Arrange
        var states = Enumerable.Range(0, 100)
            .Select(i => new ChromosomeAnalyzer.CopyNumberState(
                "chr1", i * DefaultBinSize, (i + 1) * DefaultBinSize - 1,
                CopyNumber: 5, LogRatio: 1.32, Confidence: 0.9));

        // Act
        var results = ChromosomeAnalyzer.IdentifyWholeChromosomeAneuploidy(states).ToList();

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(results.Count, Is.EqualTo(1));
            Assert.That(results[0].CopyNumber, Is.EqualTo(5));
            Assert.That(results[0].Type, Is.EqualTo("Pentasomy"));
        });
    }

    [Test]
    [Description("High copy number (>5) formats as 'Copy number = N'")]
    public void IdentifyWholeChromosomeAneuploidy_HighCopyNumber_FormatsCorrectly()
    {
        // Arrange
        var states = Enumerable.Range(0, 100)
            .Select(i => new ChromosomeAnalyzer.CopyNumberState(
                "chr1", i * DefaultBinSize, (i + 1) * DefaultBinSize - 1,
                CopyNumber: 7, LogRatio: 1.8, Confidence: 0.9));

        // Act
        var results = ChromosomeAnalyzer.IdentifyWholeChromosomeAneuploidy(states).ToList();

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(results.Count, Is.EqualTo(1));
            Assert.That(results[0].CopyNumber, Is.EqualTo(7));
            Assert.That(results[0].Type, Is.EqualTo("Copy number = 7"));
        });
    }

    #endregion

    #region IdentifyWholeChromosomeAneuploidy - Threshold Behavior

    [Test]
    [Description("Mixed copy number below threshold does not trigger aneuploidy")]
    public void IdentifyWholeChromosomeAneuploidy_MixedCN_BelowThreshold_ReturnsEmpty()
    {
        // Arrange: 70% CN=3, 30% CN=2 (below 80% threshold)
        var states = new List<ChromosomeAnalyzer.CopyNumberState>();
        for (int i = 0; i < 70; i++)
        {
            states.Add(new ChromosomeAnalyzer.CopyNumberState(
                "chr21", i * DefaultBinSize, (i + 1) * DefaultBinSize - 1,
                CopyNumber: 3, LogRatio: 0.58, Confidence: 0.9));
        }
        for (int i = 70; i < 100; i++)
        {
            states.Add(new ChromosomeAnalyzer.CopyNumberState(
                "chr21", i * DefaultBinSize, (i + 1) * DefaultBinSize - 1,
                CopyNumber: 2, LogRatio: 0.0, Confidence: 0.9));
        }

        // Act
        var results = ChromosomeAnalyzer.IdentifyWholeChromosomeAneuploidy(states).ToList();

        // Assert
        Assert.That(results, Is.Empty,
            "Should not identify aneuploidy when dominant CN is below minFraction threshold");
    }

    [Test]
    [Description("Mixed copy number at threshold triggers aneuploidy")]
    public void IdentifyWholeChromosomeAneuploidy_MixedCN_AtThreshold_IdentifiesAneuploidy()
    {
        // Arrange: 80% CN=3, 20% CN=2 (exactly at threshold)
        var states = new List<ChromosomeAnalyzer.CopyNumberState>();
        for (int i = 0; i < 80; i++)
        {
            states.Add(new ChromosomeAnalyzer.CopyNumberState(
                "chr21", i * DefaultBinSize, (i + 1) * DefaultBinSize - 1,
                CopyNumber: 3, LogRatio: 0.58, Confidence: 0.9));
        }
        for (int i = 80; i < 100; i++)
        {
            states.Add(new ChromosomeAnalyzer.CopyNumberState(
                "chr21", i * DefaultBinSize, (i + 1) * DefaultBinSize - 1,
                CopyNumber: 2, LogRatio: 0.0, Confidence: 0.9));
        }

        // Act
        var results = ChromosomeAnalyzer.IdentifyWholeChromosomeAneuploidy(states).ToList();

        // Assert
        Assert.That(results.Count, Is.EqualTo(1),
            "Should identify aneuploidy when dominant CN is at or above minFraction threshold");
        Assert.That(results[0].Type, Is.EqualTo("Trisomy"));
    }

    [Test]
    [Description("Custom minFraction parameter works correctly")]
    public void IdentifyWholeChromosomeAneuploidy_CustomMinFraction_Works()
    {
        // Arrange: 60% CN=3, 40% CN=2
        var states = new List<ChromosomeAnalyzer.CopyNumberState>();
        for (int i = 0; i < 60; i++)
        {
            states.Add(new ChromosomeAnalyzer.CopyNumberState(
                "chr21", i * DefaultBinSize, (i + 1) * DefaultBinSize - 1,
                CopyNumber: 3, LogRatio: 0.58, Confidence: 0.9));
        }
        for (int i = 60; i < 100; i++)
        {
            states.Add(new ChromosomeAnalyzer.CopyNumberState(
                "chr21", i * DefaultBinSize, (i + 1) * DefaultBinSize - 1,
                CopyNumber: 2, LogRatio: 0.0, Confidence: 0.9));
        }

        // Act: With minFraction = 0.5, 60% should pass
        var results = ChromosomeAnalyzer.IdentifyWholeChromosomeAneuploidy(
            states, minFraction: 0.5).ToList();

        // Assert
        Assert.That(results.Count, Is.EqualTo(1),
            "Should identify aneuploidy with custom minFraction = 0.5");
    }

    [Test]
    [Description("Multiple chromosomes are handled independently")]
    public void IdentifyWholeChromosomeAneuploidy_MultipleChromosomes_IndependentClassification()
    {
        // Arrange: chr21 trisomy, chr1 normal
        var states = new List<ChromosomeAnalyzer.CopyNumberState>();
        for (int i = 0; i < 100; i++)
        {
            states.Add(new ChromosomeAnalyzer.CopyNumberState(
                "chr21", i * DefaultBinSize, (i + 1) * DefaultBinSize - 1,
                CopyNumber: 3, LogRatio: 0.58, Confidence: 0.9));
            states.Add(new ChromosomeAnalyzer.CopyNumberState(
                "chr1", i * DefaultBinSize, (i + 1) * DefaultBinSize - 1,
                CopyNumber: 2, LogRatio: 0.0, Confidence: 0.95));
        }

        // Act
        var results = ChromosomeAnalyzer.IdentifyWholeChromosomeAneuploidy(states).ToList();

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(results.Count, Is.EqualTo(1));
            Assert.That(results[0].Chromosome, Is.EqualTo("chr21"));
            Assert.That(results.Any(r => r.Chromosome == "chr1"), Is.False);
        });
    }

    #endregion

    #region IdentifyWholeChromosomeAneuploidy - Edge Cases

    [Test]
    [Description("Empty input returns empty result")]
    public void IdentifyWholeChromosomeAneuploidy_EmptyInput_ReturnsEmpty()
    {
        // Act
        var results = ChromosomeAnalyzer.IdentifyWholeChromosomeAneuploidy(
            Enumerable.Empty<ChromosomeAnalyzer.CopyNumberState>()).ToList();

        // Assert
        Assert.That(results, Is.Empty);
    }

    [Test]
    [Description("Single bin is classified correctly")]
    public void IdentifyWholeChromosomeAneuploidy_SingleBin_ClassifiedCorrectly()
    {
        // Arrange
        var states = new[]
        {
            new ChromosomeAnalyzer.CopyNumberState(
                "chr21", 0, DefaultBinSize - 1,
                CopyNumber: 3, LogRatio: 0.58, Confidence: 0.9)
        };

        // Act
        var results = ChromosomeAnalyzer.IdentifyWholeChromosomeAneuploidy(states).ToList();

        // Assert
        Assert.That(results.Count, Is.EqualTo(1));
        Assert.That(results[0].Type, Is.EqualTo("Trisomy"));
    }

    #endregion
}
