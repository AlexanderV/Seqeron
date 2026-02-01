using NUnit.Framework;
using Seqeron.Genomics;
using System.Linq;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Tests for GcSkewCalculator.
/// 
/// Test Unit: SEQ-GCSKEW-001
/// 
/// Evidence Sources:
/// - Wikipedia: GC skew, formula (G-C)/(G+C), range [-1,+1]
/// - Lobry (1996) Mol. Biol. Evol. 13:660-665 - original GC skew observations
/// - Grigoriev (1998) Nucleic Acids Res. 26:2286-2290 - cumulative GC skew method
/// 
/// Key invariant: -1 ≤ GC skew ≤ +1
/// </summary>
[TestFixture]
public class GcSkewCalculatorTests
{
    #region Formula Verification Tests (Evidence: Wikipedia, Lobry 1996)

    /// <summary>
    /// Evidence: GC skew = (G - C) / (G + C)
    /// Test: 4G, 1C → (4-1)/(4+1) = 3/5 = 0.6
    /// </summary>
    [Test]
    public void CalculateGcSkew_MoreG_ReturnsPositive()
    {
        // Arrange: GGGGC has 4G and 1C
        var sequence = new DnaSequence("GGGGC");
        const double expected = (4.0 - 1.0) / (4.0 + 1.0); // 0.6

        // Act
        double skew = GcSkewCalculator.CalculateGcSkew(sequence);

        // Assert
        Assert.That(skew, Is.EqualTo(expected).Within(0.0001));
    }

    /// <summary>
    /// Evidence: GC skew = (G - C) / (G + C)
    /// Test: 0G, 5C → (0-5)/(0+5) = -5/5 = -1.0 (boundary minimum)
    /// </summary>
    [Test]
    public void CalculateGcSkew_MoreC_ReturnsNegative()
    {
        // Arrange: CCCCC has 0G and 5C
        var sequence = new DnaSequence("CCCCC");
        const double expected = (0.0 - 5.0) / (0.0 + 5.0); // -1.0

        // Act
        double skew = GcSkewCalculator.CalculateGcSkew(sequence);

        // Assert
        Assert.That(skew, Is.EqualTo(expected).Within(0.0001));
    }

    /// <summary>
    /// Evidence: GC skew = (G - C) / (G + C)
    /// Test: 5G, 0C → (5-0)/(5+0) = 5/5 = +1.0 (boundary maximum)
    /// </summary>
    [Test]
    public void CalculateGcSkew_AllG_ReturnsPlusOne()
    {
        // Arrange: GGGGG has 5G and 0C
        var sequence = new DnaSequence("GGGGG");
        const double expected = 1.0;

        // Act
        double skew = GcSkewCalculator.CalculateGcSkew(sequence);

        // Assert
        Assert.That(skew, Is.EqualTo(expected).Within(0.0001));
    }

    /// <summary>
    /// Evidence: GC skew = (G - C) / (G + C)
    /// Test: 2G, 2C → (2-2)/(2+2) = 0/4 = 0.0
    /// </summary>
    [Test]
    public void CalculateGcSkew_EqualGC_ReturnsZero()
    {
        // Arrange: GCGC has 2G and 2C
        var sequence = new DnaSequence("GCGC");
        const double expected = 0.0;

        // Act
        double skew = GcSkewCalculator.CalculateGcSkew(sequence);

        // Assert
        Assert.That(skew, Is.EqualTo(expected).Within(0.0001));
    }

    /// <summary>
    /// Evidence: Division by zero protection when G+C=0
    /// Test: No G or C bases → return 0 (protected)
    /// </summary>
    [Test]
    public void CalculateGcSkew_NoGC_ReturnsZero()
    {
        // Arrange: AAATTT has 0G and 0C
        var sequence = new DnaSequence("AAATTT");

        // Act
        double skew = GcSkewCalculator.CalculateGcSkew(sequence);

        // Assert: Protected from division by zero
        Assert.That(skew, Is.EqualTo(0));
    }

    /// <summary>
    /// Test string overload produces same result as DnaSequence overload.
    /// </summary>
    [Test]
    public void CalculateGcSkew_StringOverload_Works()
    {
        // Arrange
        const string seq = "GGGGC";
        const double expected = (4.0 - 1.0) / (4.0 + 1.0); // 0.6

        // Act
        double skew = GcSkewCalculator.CalculateGcSkew(seq);

        // Assert
        Assert.That(skew, Is.EqualTo(expected).Within(0.0001));
    }

    /// <summary>
    /// Evidence: Wikipedia - GC skew range is [-1, +1]
    /// Invariant test: All results must be within valid range.
    /// </summary>
    [Test]
    public void CalculateGcSkew_RangeIsMinusOneToOne()
    {
        // Arrange
        var sequence = new DnaSequence("ATGCATGCATGC");

        // Act
        double skew = GcSkewCalculator.CalculateGcSkew(sequence);

        // Assert: Invariant check
        Assert.Multiple(() =>
        {
            Assert.That(skew, Is.GreaterThanOrEqualTo(-1.0), "GC skew must be >= -1");
            Assert.That(skew, Is.LessThanOrEqualTo(1.0), "GC skew must be <= +1");
        });
    }

    /// <summary>
    /// Evidence: Wikipedia - GC skew range is [-1, +1]
    /// Property-based invariant test with multiple sequences.
    /// </summary>
    [Test]
    public void CalculateGcSkew_AllSequences_RespectRangeInvariant()
    {
        // Arrange: Various sequence compositions
        var testSequences = new[]
        {
            "GGGGG",      // All G → +1
            "CCCCC",      // All C → -1
            "GCGCGC",     // Equal → 0
            "AAATTT",     // No G/C → 0
            "ATGCATGC",   // Mixed
            "GGGC",       // Positive skew
            "GCCC",       // Negative skew
            "ATATATATA",  // No G/C
            "GCATGCAT",   // Balanced
        };

        foreach (string seq in testSequences)
        {
            // Act
            double skew = GcSkewCalculator.CalculateGcSkew(seq);

            // Assert
            Assert.Multiple(() =>
            {
                Assert.That(skew, Is.GreaterThanOrEqualTo(-1.0),
                    $"Sequence '{seq}': GC skew {skew} must be >= -1");
                Assert.That(skew, Is.LessThanOrEqualTo(1.0),
                    $"Sequence '{seq}': GC skew {skew} must be <= +1");
            });
        }
    }

    /// <summary>
    /// Test case insensitivity - lowercase should work identically.
    /// </summary>
    [Test]
    public void CalculateGcSkew_LowercaseInput_HandledCorrectly()
    {
        // Arrange
        const string uppercase = "GGGGC";
        const string lowercase = "ggggc";

        // Act
        double skewUpper = GcSkewCalculator.CalculateGcSkew(uppercase);
        double skewLower = GcSkewCalculator.CalculateGcSkew(lowercase);

        // Assert
        Assert.That(skewLower, Is.EqualTo(skewUpper).Within(0.0001));
    }

    #endregion

    #region Windowed GC Skew Tests (Evidence: Grigoriev 1998)

    /// <summary>
    /// Evidence: Grigoriev 1998 - sliding window analysis for GC skew.
    /// </summary>
    [Test]
    public void CalculateWindowedGcSkew_ReturnsMultiplePoints()
    {
        // Arrange: Two windows with opposite skew
        var sequence = new DnaSequence("GGGGCCCC" + "CCCCGGGG");

        // Act
        var points = GcSkewCalculator.CalculateWindowedGcSkew(sequence, windowSize: 8, stepSize: 8).ToList();

        // Assert
        Assert.That(points.Count, Is.EqualTo(2));
    }

    /// <summary>
    /// Evidence: Grigoriev 1998 - positions reported at window centers.
    /// </summary>
    [Test]
    public void CalculateWindowedGcSkew_CorrectPositions()
    {
        // Arrange
        var sequence = new DnaSequence("GGGGCCCCGGGGCCCC");

        // Act
        var points = GcSkewCalculator.CalculateWindowedGcSkew(sequence, windowSize: 4, stepSize: 4).ToList();

        // Assert: Positions are at the center of each window (WindowStart + WindowSize/2)
        Assert.Multiple(() =>
        {
            Assert.That(points[0].Position, Is.EqualTo(2), "Window 0-3, center at 2");
            Assert.That(points[1].Position, Is.EqualTo(6), "Window 4-7, center at 6");
            Assert.That(points[2].Position, Is.EqualTo(10), "Window 8-11, center at 10");
        });
    }

    /// <summary>
    /// Evidence: Formula (G-C)/(G+C) applied per window.
    /// Test: GGGG = +1.0, CCCC = -1.0
    /// </summary>
    [Test]
    public void CalculateWindowedGcSkew_CorrectSkewValues()
    {
        // Arrange
        var sequence = new DnaSequence("GGGGCCCC");

        // Act
        var points = GcSkewCalculator.CalculateWindowedGcSkew(sequence, windowSize: 4, stepSize: 4).ToList();

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(points[0].GcSkew, Is.EqualTo(1.0).Within(0.0001), "GGGG = (4-0)/(4+0) = 1.0");
            Assert.That(points[1].GcSkew, Is.EqualTo(-1.0).Within(0.0001), "CCCC = (0-4)/(0+4) = -1.0");
        });
    }

    /// <summary>
    /// Test overlapping windows produce more data points.
    /// </summary>
    [Test]
    public void CalculateWindowedGcSkew_OverlappingWindows()
    {
        // Arrange
        var sequence = new DnaSequence("GGGGCCCC");

        // Act
        var points = GcSkewCalculator.CalculateWindowedGcSkew(sequence, windowSize: 4, stepSize: 2).ToList();

        // Assert: More points due to overlap
        Assert.That(points.Count, Is.GreaterThan(2));
    }

    /// <summary>
    /// Edge case: Empty sequence returns empty collection.
    /// </summary>
    [Test]
    public void CalculateWindowedGcSkew_EmptySequence_ReturnsEmpty()
    {
        // Act
        var points = GcSkewCalculator.CalculateWindowedGcSkew("", windowSize: 10).ToList();

        // Assert
        Assert.That(points, Is.Empty);
    }

    /// <summary>
    /// Edge case: Sequence shorter than window returns empty collection.
    /// </summary>
    [Test]
    public void CalculateWindowedGcSkew_SequenceShorterThanWindow_ReturnsEmpty()
    {
        // Act
        var points = GcSkewCalculator.CalculateWindowedGcSkew("ATGC", windowSize: 10).ToList();

        // Assert
        Assert.That(points, Is.Empty);
    }

    #endregion

    #region Cumulative GC Skew Tests (Evidence: Grigoriev 1998)

    /// <summary>
    /// Evidence: Grigoriev 1998 - cumulative GC skew is sum of window skews.
    /// Test: GGGG(+1) + CCCC(-1) + GGGG(+1) + CCCC(-1) = oscillating pattern
    /// </summary>
    [Test]
    public void CalculateCumulativeGcSkew_AccumulatesCorrectly()
    {
        // Arrange: Alternating high G and high C windows
        var sequence = new DnaSequence("GGGGCCCCGGGGCCCC");

        // Act
        var points = GcSkewCalculator.CalculateCumulativeGcSkew(sequence, windowSize: 4).ToList();

        // Assert: Verify cumulative accumulation
        Assert.Multiple(() =>
        {
            Assert.That(points.Count, Is.EqualTo(4));
            Assert.That(points[0].CumulativeGcSkew, Is.EqualTo(1.0).Within(0.0001), "+1");
            Assert.That(points[1].CumulativeGcSkew, Is.EqualTo(0.0).Within(0.0001), "+1 + (-1) = 0");
            Assert.That(points[2].CumulativeGcSkew, Is.EqualTo(1.0).Within(0.0001), "+1 + (-1) + 1 = 1");
            Assert.That(points[3].CumulativeGcSkew, Is.EqualTo(0.0).Within(0.0001), "+1 + (-1) + 1 + (-1) = 0");
        });
    }

    /// <summary>
    /// Evidence: Grigoriev 1998 - positions at window centers.
    /// </summary>
    [Test]
    public void CalculateCumulativeGcSkew_PositionsAreCorrect()
    {
        // Arrange
        var sequence = new DnaSequence("GGGGCCCCGGGG");

        // Act
        var points = GcSkewCalculator.CalculateCumulativeGcSkew(sequence, windowSize: 4).ToList();

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(points[0].Position, Is.EqualTo(2), "Window 0-3, center at 2");
            Assert.That(points[1].Position, Is.EqualTo(6), "Window 4-7, center at 6");
            Assert.That(points[2].Position, Is.EqualTo(10), "Window 8-11, center at 10");
        });
    }

    /// <summary>
    /// Test that cumulative points also include the per-window skew value.
    /// </summary>
    [Test]
    public void CalculateCumulativeGcSkew_IncludesWindowSkew()
    {
        // Arrange
        var sequence = new DnaSequence("GGGGCCCC");

        // Act
        var points = GcSkewCalculator.CalculateCumulativeGcSkew(sequence, windowSize: 4).ToList();

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(points[0].GcSkew, Is.EqualTo(1.0).Within(0.0001), "GGGG");
            Assert.That(points[1].GcSkew, Is.EqualTo(-1.0).Within(0.0001), "CCCC");
        });
    }

    /// <summary>
    /// Edge case: Empty sequence returns empty collection.
    /// </summary>
    [Test]
    public void CalculateCumulativeGcSkew_EmptySequence_ReturnsEmpty()
    {
        // Act
        var points = GcSkewCalculator.CalculateCumulativeGcSkew("").ToList();

        // Assert
        Assert.That(points, Is.Empty);
    }

    #endregion

    #region AT Skew Tests (Related metric)

    /// <summary>
    /// AT skew = (A - T) / (A + T)
    /// Test: 4A, 1T → (4-1)/(4+1) = 0.6
    /// </summary>
    [Test]
    public void CalculateAtSkew_MoreA_ReturnsPositive()
    {
        // Arrange
        var sequence = new DnaSequence("AAAAT");
        const double expected = (4.0 - 1.0) / (4.0 + 1.0); // 0.6

        // Act
        double skew = GcSkewCalculator.CalculateAtSkew(sequence);

        // Assert
        Assert.That(skew, Is.EqualTo(expected).Within(0.0001));
    }

    /// <summary>
    /// AT skew = (A - T) / (A + T)
    /// Test: 0A, 5T → (0-5)/(0+5) = -1.0
    /// </summary>
    [Test]
    public void CalculateAtSkew_MoreT_ReturnsNegative()
    {
        // Arrange
        var sequence = new DnaSequence("TTTTT");
        const double expected = -1.0;

        // Act
        double skew = GcSkewCalculator.CalculateAtSkew(sequence);

        // Assert
        Assert.That(skew, Is.EqualTo(expected).Within(0.0001));
    }

    /// <summary>
    /// AT skew with equal A and T → 0
    /// </summary>
    [Test]
    public void CalculateAtSkew_EqualAT_ReturnsZero()
    {
        // Arrange
        var sequence = new DnaSequence("ATAT");

        // Act
        double skew = GcSkewCalculator.CalculateAtSkew(sequence);

        // Assert
        Assert.That(skew, Is.EqualTo(0).Within(0.0001));
    }

    /// <summary>
    /// Division by zero protection for AT skew.
    /// </summary>
    [Test]
    public void CalculateAtSkew_NoAT_ReturnsZero()
    {
        // Arrange
        var sequence = new DnaSequence("GCGC");

        // Act
        double skew = GcSkewCalculator.CalculateAtSkew(sequence);

        // Assert
        Assert.That(skew, Is.EqualTo(0));
    }

    #endregion

    #region Replication Origin Prediction Tests (Evidence: Grigoriev 1998)

    /// <summary>
    /// Evidence: Grigoriev 1998 - minimum of cumulative GC skew = origin of replication.
    /// </summary>
    [Test]
    public void PredictReplicationOrigin_FindsMinimum()
    {
        // Arrange: Simulate bacterial genome pattern
        // GC skew decreases to origin (in C-rich region), then increases
        var sequence = new DnaSequence(
            new string('G', 50) + new string('C', 100) + new string('G', 50));

        // Act
        var prediction = GcSkewCalculator.PredictReplicationOrigin(sequence, windowSize: 10);

        // Assert: Origin should be in the C-rich region (around position 100)
        Assert.That(prediction.PredictedOrigin, Is.InRange(40, 160));
    }

    /// <summary>
    /// Evidence: Grigoriev 1998 - maximum of cumulative GC skew = terminus.
    /// </summary>
    [Test]
    public void PredictReplicationOrigin_FindsMaximum()
    {
        // Arrange: Terminus should be at maximum cumulative skew
        var sequence = new DnaSequence(
            new string('C', 50) + new string('G', 100) + new string('C', 50));

        // Act
        var prediction = GcSkewCalculator.PredictReplicationOrigin(sequence, windowSize: 10);

        // Assert: Terminus in the G-rich region
        Assert.That(prediction.PredictedTerminus, Is.InRange(40, 160));
    }

    /// <summary>
    /// Predictions must be valid positions within sequence bounds.
    /// </summary>
    [Test]
    public void PredictReplicationOrigin_ReturnsValidPositions()
    {
        // Arrange
        var sequence = new DnaSequence("GGGGGGCCCCCCGGGGGGCCCCCC");

        // Act
        var prediction = GcSkewCalculator.PredictReplicationOrigin(sequence, windowSize: 4);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(prediction.PredictedOrigin, Is.GreaterThanOrEqualTo(0));
            Assert.That(prediction.PredictedOrigin, Is.LessThan(sequence.Length));
            Assert.That(prediction.PredictedTerminus, Is.GreaterThanOrEqualTo(0));
            Assert.That(prediction.PredictedTerminus, Is.LessThan(sequence.Length));
        });
    }

    #endregion

    #region GC Content Analysis Tests

    /// <summary>
    /// Comprehensive analysis returns windowed data points.
    /// </summary>
    [Test]
    public void AnalyzeGcContent_ReturnsMultipleWindowedPoints()
    {
        // Arrange
        var sequence = new DnaSequence("GCGCATATGCGC");

        // Act
        var result = GcSkewCalculator.AnalyzeGcContent(sequence, windowSize: 4, stepSize: 4);

        // Assert
        Assert.That(result.WindowedGcContent.Count, Is.EqualTo(3));
    }

    /// <summary>
    /// GC content correctly calculated per window.
    /// Window 1: GCGC = 100% GC, Window 2: ATAT = 0% GC
    /// </summary>
    [Test]
    public void AnalyzeGcContent_CalculatesGcPercent()
    {
        // Arrange
        var sequence = new DnaSequence("GCGCATAT");

        // Act
        var result = GcSkewCalculator.AnalyzeGcContent(sequence, windowSize: 4, stepSize: 4);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(result.WindowedGcContent[0].GcContent, Is.EqualTo(100).Within(0.01));
            Assert.That(result.WindowedGcContent[1].GcContent, Is.EqualTo(0).Within(0.01));
        });
    }

    /// <summary>
    /// ATGC sequence has 50% GC content (2 of 4 bases are G or C).
    /// </summary>
    [Test]
    public void AnalyzeGcContent_50PercentGc()
    {
        // Arrange
        var sequence = new DnaSequence("ATGC");

        // Act
        var result = GcSkewCalculator.AnalyzeGcContent(sequence, windowSize: 4, stepSize: 4);

        // Assert
        Assert.That(result.OverallGcContent, Is.EqualTo(50).Within(0.01));
    }

    /// <summary>
    /// Analysis returns all expected overall metrics.
    /// </summary>
    [Test]
    public void AnalyzeGcContent_ReturnsOverallMetrics()
    {
        // Arrange
        var sequence = new DnaSequence("GGGGCCCCATAT");

        // Act
        var result = GcSkewCalculator.AnalyzeGcContent(sequence, windowSize: 4, stepSize: 4);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(result.OverallGcContent, Is.GreaterThan(0));
            Assert.That(result.SequenceLength, Is.EqualTo(12));
        });
    }

    #endregion

    #region Edge Cases and Exception Handling

    /// <summary>
    /// Empty string input returns 0.
    /// </summary>
    [Test]
    public void CalculateGcSkew_EmptySequence_ReturnsZero()
    {
        // Act
        double skew = GcSkewCalculator.CalculateGcSkew("");

        // Assert
        Assert.That(skew, Is.EqualTo(0));
    }

    /// <summary>
    /// Null DnaSequence throws ArgumentNullException.
    /// </summary>
    [Test]
    public void CalculateGcSkew_NullSequence_ThrowsException()
    {
        Assert.Throws<ArgumentNullException>(() =>
            GcSkewCalculator.CalculateGcSkew((DnaSequence)null!));
    }

    /// <summary>
    /// Null DnaSequence for windowed analysis throws ArgumentNullException.
    /// </summary>
    [Test]
    public void CalculateWindowedGcSkew_NullSequence_ThrowsException()
    {
        Assert.Throws<ArgumentNullException>(() =>
            GcSkewCalculator.CalculateWindowedGcSkew((DnaSequence)null!, 10).ToList());
    }

    /// <summary>
    /// Window size ≤ 0 throws ArgumentOutOfRangeException.
    /// </summary>
    [Test]
    public void CalculateWindowedGcSkew_ZeroWindowSize_ThrowsException()
    {
        // Arrange
        var sequence = new DnaSequence("ATGC");

        // Act & Assert
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            GcSkewCalculator.CalculateWindowedGcSkew(sequence, windowSize: 0).ToList());
    }

    /// <summary>
    /// Negative step size throws ArgumentOutOfRangeException.
    /// </summary>
    [Test]
    public void CalculateWindowedGcSkew_NegativeStepSize_ThrowsException()
    {
        // Arrange
        var sequence = new DnaSequence("ATGC");

        // Act & Assert
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            GcSkewCalculator.CalculateWindowedGcSkew(sequence, windowSize: 2, stepSize: -1).ToList());
    }

    /// <summary>
    /// Null DnaSequence for cumulative analysis throws ArgumentNullException.
    /// </summary>
    [Test]
    public void CalculateCumulativeGcSkew_NullSequence_ThrowsException()
    {
        Assert.Throws<ArgumentNullException>(() =>
            GcSkewCalculator.CalculateCumulativeGcSkew((DnaSequence)null!).ToList());
    }

    /// <summary>
    /// Null DnaSequence for origin prediction throws ArgumentNullException.
    /// </summary>
    [Test]
    public void PredictReplicationOrigin_NullSequence_ThrowsException()
    {
        Assert.Throws<ArgumentNullException>(() =>
            GcSkewCalculator.PredictReplicationOrigin((DnaSequence)null!, 10));
    }

    /// <summary>
    /// Null DnaSequence for analysis throws ArgumentNullException.
    /// </summary>
    [Test]
    public void AnalyzeGcContent_NullSequence_ThrowsException()
    {
        Assert.Throws<ArgumentNullException>(() =>
            GcSkewCalculator.AnalyzeGcContent((DnaSequence)null!, 10));
    }

    #endregion
}
