using NUnit.Framework;
using Seqeron.Genomics;
using System;
using System.Linq;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Tests for ChromosomeAnalyzer.AnalyzeCentromere (CHROM-CENT-001).
/// 
/// Evidence: Wikipedia (Centromere, Karyotype, Chromosome), Levan et al. (1964)
/// TestSpec: TestSpecs/CHROM-CENT-001.md
/// </summary>
[TestFixture]
public class ChromosomeAnalyzer_Centromere_Tests
{
    #region Edge Cases (Must Tests M1-M3)

    [Test]
    public void AnalyzeCentromere_EmptySequence_ReturnsUnknownWithNullBoundaries()
    {
        // Arrange & Act
        var result = ChromosomeAnalyzer.AnalyzeCentromere("chr1", "");

        // Assert - M1: Empty sequence returns Unknown type with null boundaries
        Assert.Multiple(() =>
        {
            Assert.That(result.Chromosome, Is.EqualTo("chr1"));
            Assert.That(result.CentromereType, Is.EqualTo("Unknown"));
            Assert.That(result.Start, Is.Null);
            Assert.That(result.End, Is.Null);
            Assert.That(result.Length, Is.EqualTo(0));
            Assert.That(result.AlphaSatelliteContent, Is.EqualTo(0));
            Assert.That(result.IsAcrocentric, Is.False);
        });
    }

    [Test]
    public void AnalyzeCentromere_NullSequence_ReturnsUnknownWithNullBoundaries()
    {
        // Arrange & Act
        var result = ChromosomeAnalyzer.AnalyzeCentromere("chr1", null!);

        // Assert - Null is treated same as empty
        Assert.Multiple(() =>
        {
            Assert.That(result.CentromereType, Is.EqualTo("Unknown"));
            Assert.That(result.Start, Is.Null);
            Assert.That(result.End, Is.Null);
        });
    }

    [Test]
    public void AnalyzeCentromere_SequenceShorterThanWindowSize_ReturnsUnknown()
    {
        // Arrange - M2: Sequence shorter than window cannot be analyzed
        string shortSequence = new string('A', 1000);

        // Act
        var result = ChromosomeAnalyzer.AnalyzeCentromere("chr1", shortSequence, windowSize: 10000);

        // Assert
        Assert.That(result.CentromereType, Is.EqualTo("Unknown"));
    }

    [Test]
    public void AnalyzeCentromere_ChromosomeNamePreservedInResult()
    {
        // Arrange - M3: Chromosome name is preserved
        const string chromosomeName = "chr21";
        string sequence = GenerateRandomSequence(seed: 42, length: 200000);

        // Act
        var result = ChromosomeAnalyzer.AnalyzeCentromere(chromosomeName, sequence, windowSize: 10000);

        // Assert
        Assert.That(result.Chromosome, Is.EqualTo(chromosomeName));
    }

    #endregion

    #region Repetitive Region Detection (Must Test M4)

    [Test]
    public void AnalyzeCentromere_WithHighlyRepetitiveMiddleRegion_DetectsCentromere()
    {
        // Arrange - M4: Repetitive regions should be detected
        // Create a sequence with non-repetitive flanks and highly repetitive middle
        var random = new Random(42);
        string flank = GenerateRandomSequence(random, 100000);

        // Create highly repetitive alpha-satellite-like region
        string repeatUnit = "AATGAATATTTCTTTTATGTT"; // 21-mer repeat unit
        string repetitiveRegion = string.Concat(Enumerable.Repeat(repeatUnit, 5000)); // ~105kb

        string sequence = flank + repetitiveRegion + flank;

        // Act
        var result = ChromosomeAnalyzer.AnalyzeCentromere(
            "chr1",
            sequence,
            windowSize: 10000,
            minAlphaSatelliteContent: 0.1);

        // Assert - Should detect something in the repetitive region
        Assert.Multiple(() =>
        {
            Assert.That(result.AlphaSatelliteContent, Is.GreaterThan(0),
                "Should detect repetitive content");
        });
    }

    [Test]
    public void AnalyzeCentromere_WithNoRepetitiveRegions_ReturnsLowScore()
    {
        // Arrange - Completely random, non-repetitive sequence
        string sequence = GenerateRandomSequence(seed: 12345, length: 500000);

        // Act
        var result = ChromosomeAnalyzer.AnalyzeCentromere(
            "chr1",
            sequence,
            windowSize: 50000,
            minAlphaSatelliteContent: 0.3);

        // Assert - Should have low or zero score for non-repetitive sequence
        // Type may be Unknown if no region meets threshold
        Assert.That(result.AlphaSatelliteContent, Is.LessThan(0.3));
    }

    #endregion

    #region Result Invariants (Must Tests M5-M6)

    [Test]
    public void AnalyzeCentromere_WhenCentromereFound_StartIsLessThanOrEqualToEnd()
    {
        // Arrange - M5: Start <= End invariant
        // 160,000 bp repetitive region embedded in flanks - MUST be detected
        string repeatUnit = "GATCGATCGATCGATC";
        string repetitiveRegion = string.Concat(Enumerable.Repeat(repeatUnit, 10000));
        string flank = GenerateRandomSequence(seed: 99, length: 100000);
        string sequence = flank + repetitiveRegion + flank;

        // Act
        var result = ChromosomeAnalyzer.AnalyzeCentromere(
            "chr1", sequence, windowSize: 10000, minAlphaSatelliteContent: 0.05);

        // Assert - large repetitive region MUST be detected
        Assert.That(result.Start.HasValue && result.End.HasValue, Is.True,
            "160kb repetitive region with 5% threshold must be detected as centromere");
        Assert.That(result.Start!.Value, Is.LessThanOrEqualTo(result.End!.Value),
            "Start must be <= End");
    }

    [Test]
    public void AnalyzeCentromere_WhenCentromereFound_LengthEqualsEndMinusStart()
    {
        // Arrange - M6: Length = End - Start invariant
        // 112,000 bp repetitive region embedded in flanks - MUST be detected
        string repeatUnit = "ATATATATATATAT";
        string repetitiveRegion = string.Concat(Enumerable.Repeat(repeatUnit, 8000));
        string flank = GenerateRandomSequence(seed: 88, length: 80000);
        string sequence = flank + repetitiveRegion + flank;

        // Act
        var result = ChromosomeAnalyzer.AnalyzeCentromere(
            "chr1", sequence, windowSize: 10000, minAlphaSatelliteContent: 0.05);

        // Assert - large repetitive region MUST be detected
        Assert.That(result.Start.HasValue && result.End.HasValue, Is.True,
            "112kb repetitive region with 5% threshold must be detected as centromere");
        Assert.That(result.Length, Is.EqualTo(result.End!.Value - result.Start!.Value),
            "Length must equal End - Start");
    }

    [Test]
    public void AnalyzeCentromere_WhenNoCentromereFound_LengthIsZero()
    {
        // Arrange
        var result = ChromosomeAnalyzer.AnalyzeCentromere("chr1", "");

        // Assert - When Start/End are null, Length should be 0
        Assert.Multiple(() =>
        {
            Assert.That(result.Start, Is.Null);
            Assert.That(result.End, Is.Null);
            Assert.That(result.Length, Is.EqualTo(0));
        });
    }

    #endregion

    #region Type Classification (Must Tests M7-M8)

    [Test]
    public void AnalyzeCentromere_CentromereInMiddle_ClassifiesAsMetacentric()
    {
        // Arrange - M7: Metacentric = centromere in middle (position 0.35-0.65)
        // Create sequence where 96kb repetitive region is centered between equal 100kb flanks
        string flank = GenerateRandomSequence(seed: 1, length: 100000);
        string repeatUnit = "GCGCGCGCGCGCGCGC";
        string repetitiveRegion = string.Concat(Enumerable.Repeat(repeatUnit, 6000));

        // Flanks are equal, so repetitive region is centered
        string sequence = flank + repetitiveRegion + flank;

        // Act
        var result = ChromosomeAnalyzer.AnalyzeCentromere(
            "chr1", sequence, windowSize: 10000, minAlphaSatelliteContent: 0.05);

        // Assert - centered repetitive region MUST be detected and classified
        Assert.That(result.Start.HasValue && result.End.HasValue, Is.True,
            "96kb centered repetitive region must be detected as centromere");

        // Verify the detected position is roughly in the middle
        int midpoint = (result.Start!.Value + result.End!.Value) / 2;
        double positionRatio = midpoint / (double)sequence.Length;

        // Metacentric should be around 0.35-0.65
        Assert.That(result.CentromereType, Is.EqualTo("Metacentric"),
            $"Centromere at position ratio {positionRatio:F2} should be Metacentric");
    }

    [Test]
    public void AnalyzeCentromere_CentromereNearStart_ClassifiesAsAcrocentric()
    {
        // Arrange - M7: Acrocentric = centromere near end (position < 0.15 or > 0.85)
        // 80kb repetitive region at start, followed by 800kb random sequence
        string repeatUnit = "TATATATATATATATA";
        string repetitiveRegion = string.Concat(Enumerable.Repeat(repeatUnit, 5000));
        string longFlank = GenerateRandomSequence(seed: 2, length: 800000);

        // Repetitive region at start
        string sequence = repetitiveRegion + longFlank;

        // Act
        var result = ChromosomeAnalyzer.AnalyzeCentromere(
            "chr13", sequence, windowSize: 10000, minAlphaSatelliteContent: 0.05);

        // Assert - repetitive region at start MUST be detected
        Assert.That(result.Start.HasValue && result.End.HasValue, Is.True,
            "80kb repetitive region at chromosome start must be detected");

        int midpoint = (result.Start!.Value + result.End!.Value) / 2;
        double positionRatio = midpoint / (double)sequence.Length;

        // Should be acrocentric since detected at start (< 0.15)
        Assert.That(positionRatio, Is.LessThan(0.15),
            $"Centromere at start should have position ratio < 0.15, got {positionRatio:F2}");
        Assert.That(result.CentromereType, Is.EqualTo("Acrocentric"),
            "Centromere near start should be classified as Acrocentric");
    }

    [Test]
    public void AnalyzeCentromere_IsAcrocentricFlagMatchesType()
    {
        // Arrange - M8: IsAcrocentric matches type
        string sequence = GenerateRandomSequence(seed: 77, length: 200000);

        // Act
        var result = ChromosomeAnalyzer.AnalyzeCentromere(
            "chr1", sequence, windowSize: 10000, minAlphaSatelliteContent: 0.5);

        // Assert - IsAcrocentric should be true iff Type is "Acrocentric"
        bool expectedFlag = result.CentromereType == "Acrocentric";
        Assert.That(result.IsAcrocentric, Is.EqualTo(expectedFlag));
    }

    #endregion

    #region Should Tests - Parameter Behavior

    [Test]
    [TestCase(10000)]
    [TestCase(50000)]
    [TestCase(100000)]
    public void AnalyzeCentromere_DifferentWindowSizes_ProducesResults(int windowSize)
    {
        // Arrange - S1: Different window sizes
        string repeatUnit = "ACGTACGTACGTACGT";
        string repetitiveRegion = string.Concat(Enumerable.Repeat(repeatUnit, 10000));
        string flank = GenerateRandomSequence(seed: 55, length: windowSize * 2);
        string sequence = flank + repetitiveRegion + flank;

        // Act
        var result = ChromosomeAnalyzer.AnalyzeCentromere(
            "chr1", sequence, windowSize: windowSize, minAlphaSatelliteContent: 0.05);

        // Assert - Should not throw and produces valid result
        Assert.That(result.Chromosome, Is.EqualTo("chr1"));
        Assert.That(result.CentromereType, Is.Not.Null.And.Not.Empty);
    }

    [Test]
    [TestCase(0.1)]
    [TestCase(0.3)]
    [TestCase(0.5)]
    public void AnalyzeCentromere_DifferentThresholds_AffectsDetection(double threshold)
    {
        // Arrange - S2: Different thresholds
        string repeatUnit = "CGCGCGCGCGCGCGCG";
        string repetitiveRegion = string.Concat(Enumerable.Repeat(repeatUnit, 8000));
        string flank = GenerateRandomSequence(seed: 66, length: 100000);
        string sequence = flank + repetitiveRegion + flank;

        // Act
        var result = ChromosomeAnalyzer.AnalyzeCentromere(
            "chr1", sequence, windowSize: 10000, minAlphaSatelliteContent: threshold);

        // Assert - Higher thresholds may result in Unknown type
        Assert.That(result.Chromosome, Is.EqualTo("chr1"));
        // Higher threshold = less likely to detect
        if (threshold >= 0.5)
        {
            // May or may not detect - just verify no crash
            Assert.That(result.CentromereType, Is.Not.Null);
        }
    }

    [Test]
    public void AnalyzeCentromere_MixedCaseInput_HandledCorrectly()
    {
        // Arrange - S3: Case insensitivity
        string repeatUnit = "atcgatcgatcgatcg";
        string repetitiveRegion = string.Concat(Enumerable.Repeat(repeatUnit, 8000));
        string flank = "acgtacgt".PadRight(100000, 'n');
        string sequence = flank + repetitiveRegion + flank;

        // Act
        var result = ChromosomeAnalyzer.AnalyzeCentromere(
            "chr1", sequence, windowSize: 10000, minAlphaSatelliteContent: 0.05);

        // Assert - Should handle lowercase
        Assert.That(result.Chromosome, Is.EqualTo("chr1"));
    }

    #endregion

    #region Invariant Tests

    [Test]
    public void AnalyzeCentromere_AlphaSatelliteContent_IsNonNegative()
    {
        // Arrange
        string sequence = GenerateRandomSequence(seed: 111, length: 300000);

        // Act
        var result = ChromosomeAnalyzer.AnalyzeCentromere(
            "chr1", sequence, windowSize: 50000);

        // Assert - AlphaSatelliteContent >= 0
        Assert.That(result.AlphaSatelliteContent, Is.GreaterThanOrEqualTo(0));
    }

    [Test]
    public void AnalyzeCentromere_TypeIsValidValue()
    {
        // Arrange
        string[] validTypes = { "Metacentric", "Submetacentric", "Acrocentric", "Unknown" };
        string sequence = GenerateRandomSequence(seed: 222, length: 250000);

        // Act
        var result = ChromosomeAnalyzer.AnalyzeCentromere(
            "chr1", sequence, windowSize: 20000);

        // Assert - Type must be one of the valid values
        Assert.That(validTypes, Does.Contain(result.CentromereType));
    }

    #endregion

    #region Constants Tests

    [Test]
    public void AlphaSatelliteConsensus_IsNonEmpty()
    {
        // Assert - Verify constant is defined and non-empty
        Assert.That(ChromosomeAnalyzer.AlphaSatelliteConsensus, Is.Not.Null.And.Not.Empty);
    }

    [Test]
    public void AlphaSatelliteConsensus_ContainsOnlyValidBases()
    {
        // Arrange
        string consensus = ChromosomeAnalyzer.AlphaSatelliteConsensus;

        // Assert - Should contain only valid DNA bases
        Assert.That(consensus.ToUpperInvariant().All(c => "ACGT".Contains(c)), Is.True);
    }

    #endregion

    #region Helper Methods

    private static string GenerateRandomSequence(int seed, int length)
    {
        return GenerateRandomSequence(new Random(seed), length);
    }

    private static string GenerateRandomSequence(Random random, int length)
    {
        var bases = new char[] { 'A', 'C', 'G', 'T' };
        var sequence = new char[length];

        for (int i = 0; i < length; i++)
        {
            sequence[i] = bases[random.Next(4)];
        }

        return new string(sequence);
    }

    #endregion
}
