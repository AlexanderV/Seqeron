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

    [Test]
    public void AnalyzeCentromere_AllSameBase_DetectsHighRepeatContent()
    {
        // Arrange - Edge case: degenerate sequence of uniform base
        // All-same-base has maximum k-mer repeat content and zero GC variability
        string sequence = new string('A', 300000);

        // Act
        var result = ChromosomeAnalyzer.AnalyzeCentromere(
            "chr1", sequence, windowSize: 10000, minAlphaSatelliteContent: 0.3);

        // Assert - Uniform sequence is maximally repetitive → must be detected
        Assert.Multiple(() =>
        {
            Assert.That(result.Start, Is.Not.Null,
                "All-same-base sequence is maximally repetitive and must be detected");
            Assert.That(result.End, Is.Not.Null);
            Assert.That(result.AlphaSatelliteContent, Is.GreaterThan(0.3),
                "Uniform base has maximum repeat content");
            Assert.That(result.CentromereType, Is.Not.EqualTo("Unknown"));
        });
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

        // Assert - Centromere must be detected in the repetitive region
        Assert.Multiple(() =>
        {
            Assert.That(result.Start, Is.Not.Null,
                "105kb repetitive region must produce non-null Start");
            Assert.That(result.End, Is.Not.Null,
                "105kb repetitive region must produce non-null End");
            Assert.That(result.CentromereType, Is.Not.EqualTo("Unknown"),
                "Detected centromere must have a classification");
            Assert.That(result.AlphaSatelliteContent, Is.GreaterThan(0.1),
                "Repeat content score must exceed the threshold");
        });
    }

    [Test]
    public void AnalyzeCentromere_WithNoRepetitiveRegions_ReturnsUnknown()
    {
        // Arrange - Completely random, non-repetitive sequence
        string sequence = GenerateRandomSequence(seed: 12345, length: 500000);

        // Act
        var result = ChromosomeAnalyzer.AnalyzeCentromere(
            "chr1",
            sequence,
            windowSize: 50000,
            minAlphaSatelliteContent: 0.3);

        // Assert - Non-repetitive sequence: no centromere detected
        Assert.Multiple(() =>
        {
            Assert.That(result.CentromereType, Is.EqualTo("Unknown"),
                "Random non-repetitive sequence should not be classified");
            Assert.That(result.Start, Is.Null);
            Assert.That(result.End, Is.Null);
            Assert.That(result.Length, Is.EqualTo(0));
            Assert.That(result.AlphaSatelliteContent, Is.LessThan(0.3),
                "Score must be below threshold for non-repetitive input");
        });
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

    #endregion

    #region Type Classification per Levan et al. (1964) — Must Tests M7-M8

    // Levan A, Fredga K, Sandberg AA (1964). "Nomenclature for centromeric position
    // on chromosomes". Hereditas. 52(2):201-220.
    //
    // Classification is based on arm ratio (q/p) where:
    //   p = short arm, q = long arm, ratio = q/p
    //
    //   ratio <= 1.7  → Metacentric
    //   ratio 1.7-3.0 → Submetacentric
    //   ratio 3.0-7.0 → Subtelocentric
    //   ratio >= 7.0  → Acrocentric
    //   p = 0         → Telocentric

    [Test]
    public void AnalyzeCentromere_CentromereInMiddle_ClassifiesAsMetacentric()
    {
        // Arrange - Metacentric: arm ratio <= 1.7 → equal flanks give ratio ~1.0
        string flank = GenerateRandomSequence(seed: 1, length: 100000);
        string repeatUnit = "GCGCGCGCGCGCGCGC";
        string repetitiveRegion = string.Concat(Enumerable.Repeat(repeatUnit, 6000));

        string sequence = flank + repetitiveRegion + flank;

        // Act
        var result = ChromosomeAnalyzer.AnalyzeCentromere(
            "chr1", sequence, windowSize: 10000, minAlphaSatelliteContent: 0.05);

        // Assert
        Assert.That(result.Start.HasValue && result.End.HasValue, Is.True,
            "96kb centered repetitive region must be detected as centromere");

        int midpoint = (result.Start!.Value + result.End!.Value) / 2;
        int pArm = Math.Min(midpoint, sequence.Length - midpoint);
        int qArm = Math.Max(midpoint, sequence.Length - midpoint);
        double armRatio = (double)qArm / pArm;

        Assert.That(armRatio, Is.LessThanOrEqualTo(1.7),
            $"Equal flanks should give arm ratio ~1.0, got {armRatio:F2}");
        Assert.That(result.CentromereType, Is.EqualTo("Metacentric"),
            $"Arm ratio {armRatio:F2} <= 1.7 → Metacentric per Levan (1964)");
    }

    [Test]
    public void AnalyzeCentromere_CentromereOffCenter_ClassifiesAsSubmetacentric()
    {
        // Arrange - Submetacentric: arm ratio 1.7-3.0
        // Target ratio ~2.0 → p:q = 1:2 → short flank : long flank = 1:2
        // repetitive region ~40kb in between
        string repeatUnit = "TGTGTGTGTGTGTGTG";
        string repetitiveRegion = string.Concat(Enumerable.Repeat(repeatUnit, 2500)); // 40kb
        string shortFlank = GenerateRandomSequence(seed: 10, length: 180000);
        string longFlank = GenerateRandomSequence(seed: 11, length: 380000);

        string sequence = shortFlank + repetitiveRegion + longFlank;
        // Total ~600kb, centromere midpoint ~200kb → p=200kb, q=400kb, ratio=2.0

        // Act
        var result = ChromosomeAnalyzer.AnalyzeCentromere(
            "chr2", sequence, windowSize: 10000, minAlphaSatelliteContent: 0.05);

        // Assert
        Assert.That(result.Start.HasValue && result.End.HasValue, Is.True,
            "40kb repetitive region must be detected as centromere");

        int midpoint = (result.Start!.Value + result.End!.Value) / 2;
        int pArm = Math.Min(midpoint, sequence.Length - midpoint);
        int qArm = Math.Max(midpoint, sequence.Length - midpoint);
        double armRatio = (double)qArm / pArm;

        Assert.That(armRatio, Is.GreaterThan(1.7).And.LessThanOrEqualTo(3.0),
            $"Arm ratio should be 1.7-3.0 for submetacentric, got {armRatio:F2}");
        Assert.That(result.CentromereType, Is.EqualTo("Submetacentric"),
            $"Arm ratio {armRatio:F2} in (1.7, 3.0] → Submetacentric per Levan (1964)");
    }

    [Test]
    public void AnalyzeCentromere_CentromereNearEnd_ClassifiesAsSubtelocentric()
    {
        // Arrange - Subtelocentric: arm ratio 3.0-7.0
        // Target ratio ~5.0 → p:q = 1:5 → short flank : long flank = 1:5
        string repeatUnit = "CACACACACACACACAC";
        string repetitiveRegion = string.Concat(Enumerable.Repeat(repeatUnit, 2500)); // 40kb
        string shortFlank = GenerateRandomSequence(seed: 20, length: 80000);
        string longFlank = GenerateRandomSequence(seed: 21, length: 480000);

        string sequence = shortFlank + repetitiveRegion + longFlank;
        // Total ~600kb, centromere midpoint ~100kb → p=100kb, q=500kb, ratio=5.0

        // Act
        var result = ChromosomeAnalyzer.AnalyzeCentromere(
            "chr4", sequence, windowSize: 10000, minAlphaSatelliteContent: 0.05);

        // Assert
        Assert.That(result.Start.HasValue && result.End.HasValue, Is.True,
            "40kb repetitive region must be detected as centromere");

        int midpoint = (result.Start!.Value + result.End!.Value) / 2;
        int pArm = Math.Min(midpoint, sequence.Length - midpoint);
        int qArm = Math.Max(midpoint, sequence.Length - midpoint);
        double armRatio = (double)qArm / pArm;

        Assert.That(armRatio, Is.GreaterThan(3.0).And.LessThanOrEqualTo(7.0),
            $"Arm ratio should be 3.0-7.0 for subtelocentric, got {armRatio:F2}");
        Assert.That(result.CentromereType, Is.EqualTo("Subtelocentric"),
            $"Arm ratio {armRatio:F2} in (3.0, 7.0] → Subtelocentric per Levan (1964)");
    }

    [Test]
    public void AnalyzeCentromere_CentromereNearStart_ClassifiesAsAcrocentric()
    {
        // Arrange - Acrocentric: arm ratio > 7.0
        // 80kb repetitive region at start, followed by 800kb random sequence
        // Centromere midpoint ~40kb → p=40kb, q=840kb, ratio=21
        string repeatUnit = "TATATATATATATATA";
        string repetitiveRegion = string.Concat(Enumerable.Repeat(repeatUnit, 5000)); // 80kb
        string longFlank = GenerateRandomSequence(seed: 2, length: 800000);

        string sequence = repetitiveRegion + longFlank;

        // Act
        var result = ChromosomeAnalyzer.AnalyzeCentromere(
            "chr13", sequence, windowSize: 10000, minAlphaSatelliteContent: 0.05);

        // Assert
        Assert.That(result.Start.HasValue && result.End.HasValue, Is.True,
            "80kb repetitive region at chromosome start must be detected");

        int midpoint = (result.Start!.Value + result.End!.Value) / 2;
        int pArm = Math.Min(midpoint, sequence.Length - midpoint);
        int qArm = Math.Max(midpoint, sequence.Length - midpoint);
        double armRatio = (double)qArm / pArm;

        Assert.That(armRatio, Is.GreaterThan(7.0),
            $"Arm ratio should be > 7.0 for acrocentric, got {armRatio:F2}");
        Assert.That(result.CentromereType, Is.EqualTo("Acrocentric"),
            $"Arm ratio {armRatio:F2} > 7.0 → Acrocentric per Levan (1964)");
    }

    [Test]
    public void AnalyzeCentromere_AcrocentricResult_IsAcrocentricFlagIsTrue()
    {
        // Arrange - M8: IsAcrocentric must be true when Type is Acrocentric
        // Re-use acrocentric geometry: 80kb repetitive at start + 800kb random
        string repeatUnit = "TATATATATATATATA";
        string repetitiveRegion = string.Concat(Enumerable.Repeat(repeatUnit, 5000));
        string longFlank = GenerateRandomSequence(seed: 2, length: 800000);
        string sequence = repetitiveRegion + longFlank;

        // Act
        var result = ChromosomeAnalyzer.AnalyzeCentromere(
            "chr13", sequence, windowSize: 10000, minAlphaSatelliteContent: 0.05);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(result.CentromereType, Is.EqualTo("Acrocentric"));
            Assert.That(result.IsAcrocentric, Is.True,
                "IsAcrocentric must be true when Type is Acrocentric");
        });
    }

    [Test]
    public void AnalyzeCentromere_MetacentricResult_IsAcrocentricFlagIsFalse()
    {
        // Arrange - M8: IsAcrocentric must be false for non-Acrocentric types
        string flank = GenerateRandomSequence(seed: 1, length: 100000);
        string repeatUnit = "GCGCGCGCGCGCGCGC";
        string repetitiveRegion = string.Concat(Enumerable.Repeat(repeatUnit, 6000));
        string sequence = flank + repetitiveRegion + flank;

        // Act
        var result = ChromosomeAnalyzer.AnalyzeCentromere(
            "chr1", sequence, windowSize: 10000, minAlphaSatelliteContent: 0.05);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(result.CentromereType, Is.EqualTo("Metacentric"));
            Assert.That(result.IsAcrocentric, Is.False,
                "IsAcrocentric must be false when Type is Metacentric");
        });
    }

    #endregion

    #region Should Tests - Parameter Behavior

    [Test]
    [TestCase(10000)]
    [TestCase(50000)]
    [TestCase(100000)]
    public void AnalyzeCentromere_DifferentWindowSizes_AllDetectCentromere(int windowSize)
    {
        // Arrange - S1: 160kb repetitive region should be detected regardless of window size
        string repeatUnit = "ACGTACGTACGTACGT";
        string repetitiveRegion = string.Concat(Enumerable.Repeat(repeatUnit, 10000)); // 160kb
        string flank = GenerateRandomSequence(seed: 55, length: windowSize * 2);
        string sequence = flank + repetitiveRegion + flank;

        // Act
        var result = ChromosomeAnalyzer.AnalyzeCentromere(
            "chr1", sequence, windowSize: windowSize, minAlphaSatelliteContent: 0.05);

        // Assert - All window sizes must detect the 160kb repetitive region
        Assert.Multiple(() =>
        {
            Assert.That(result.Chromosome, Is.EqualTo("chr1"));
            Assert.That(result.Start, Is.Not.Null,
                $"Window size {windowSize} must detect 160kb repetitive region");
            Assert.That(result.End, Is.Not.Null);
            Assert.That(result.CentromereType, Is.Not.EqualTo("Unknown"));
        });
    }

    [Test]
    public void AnalyzeCentromere_LowThreshold_DetectsCentromere()
    {
        // Arrange - S2: Low threshold should detect centromere
        string repeatUnit = "CGCGCGCGCGCGCGCG";
        string repetitiveRegion = string.Concat(Enumerable.Repeat(repeatUnit, 8000)); // 128kb
        string flank = GenerateRandomSequence(seed: 66, length: 100000);
        string sequence = flank + repetitiveRegion + flank;

        // Act
        var result = ChromosomeAnalyzer.AnalyzeCentromere(
            "chr1", sequence, windowSize: 10000, minAlphaSatelliteContent: 0.1);

        // Assert - Low threshold detects the centromere
        Assert.Multiple(() =>
        {
            Assert.That(result.Start, Is.Not.Null,
                "128kb repetitive region with 10% threshold must be detected");
            Assert.That(result.CentromereType, Is.Not.EqualTo("Unknown"));
        });
    }

    [Test]
    public void AnalyzeCentromere_HighThreshold_ReducesSensitivity()
    {
        // Arrange - S2: High threshold should be harder to satisfy
        // Use moderate-repeat sequence that passes low threshold but not high
        string sequence = GenerateRandomSequence(seed: 66, length: 300000);

        // Act
        var lowThresholdResult = ChromosomeAnalyzer.AnalyzeCentromere(
            "chr1", sequence, windowSize: 10000, minAlphaSatelliteContent: 0.01);
        var highThresholdResult = ChromosomeAnalyzer.AnalyzeCentromere(
            "chr1", sequence, windowSize: 10000, minAlphaSatelliteContent: 0.9);

        // Assert - Higher threshold produces equal or less detection
        Assert.That(
            highThresholdResult.AlphaSatelliteContent,
            Is.LessThanOrEqualTo(lowThresholdResult.AlphaSatelliteContent + 0.001),
            "High threshold must not produce higher score than low threshold on same input");
    }

    [Test]
    public void AnalyzeCentromere_MixedCaseInput_ProducesSameResultAsUppercase()
    {
        // Arrange - S3: Case insensitivity
        string repeatUnit = "ATCGATCGATCGATCG";
        string repetitiveRegion = string.Concat(Enumerable.Repeat(repeatUnit, 8000)); // 128kb
        string flank = GenerateRandomSequence(seed: 33, length: 100000);
        string upperSequence = flank + repetitiveRegion + flank;
        string lowerSequence = upperSequence.ToLowerInvariant();

        // Act
        var upperResult = ChromosomeAnalyzer.AnalyzeCentromere(
            "chr1", upperSequence, windowSize: 10000, minAlphaSatelliteContent: 0.05);
        var lowerResult = ChromosomeAnalyzer.AnalyzeCentromere(
            "chr1", lowerSequence, windowSize: 10000, minAlphaSatelliteContent: 0.05);

        // Assert - Lowercase must produce identical results to uppercase
        Assert.Multiple(() =>
        {
            Assert.That(lowerResult.CentromereType, Is.EqualTo(upperResult.CentromereType),
                "Case must not affect classification");
            Assert.That(lowerResult.Start, Is.EqualTo(upperResult.Start),
                "Case must not affect Start position");
            Assert.That(lowerResult.End, Is.EqualTo(upperResult.End),
                "Case must not affect End position");
            Assert.That(lowerResult.AlphaSatelliteContent,
                Is.EqualTo(upperResult.AlphaSatelliteContent).Within(0.001),
                "Case must not affect score");
        });
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
        // Arrange - All types per Levan (1964) + Unknown for undetected
        string[] validTypes =
        {
            "Metacentric", "Submetacentric", "Subtelocentric",
            "Acrocentric", "Telocentric", "Unknown"
        };
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
    public void AlphaSatelliteConsensus_IsValidDnaSequence()
    {
        // Assert - Alpha-satellite consensus per Wikipedia (human centromeric repeat)
        string consensus = ChromosomeAnalyzer.AlphaSatelliteConsensus;

        Assert.Multiple(() =>
        {
            Assert.That(consensus, Is.Not.Null.And.Not.Empty);
            Assert.That(consensus.Length, Is.GreaterThan(50),
                "Alpha-satellite monomer is ~171bp; consensus should be substantial");
            Assert.That(consensus.ToUpperInvariant().All(c => "ACGT".Contains(c)), Is.True,
                "Consensus must contain only valid DNA bases");
        });
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
