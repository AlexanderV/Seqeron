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
    /// Case insensitivity: lowercase produces same result.
    /// Source: Biopython counts both cases: s.count("G") + s.count("g").
    /// GGGGC / ggggc: G=4, C=1 → (4−1)/(4+1) = 0.6
    /// </summary>
    [Test]
    public void CalculateGcSkew_LowercaseInput_HandledCorrectly()
    {
        const double expected = 0.6;

        Assert.Multiple(() =>
        {
            Assert.That(GcSkewCalculator.CalculateGcSkew("GGGGC"), Is.EqualTo(expected).Within(0.0001));
            Assert.That(GcSkewCalculator.CalculateGcSkew("ggggc"), Is.EqualTo(expected).Within(0.0001));
        });
    }

    #endregion

    #region Windowed GC Skew Tests (Evidence: Grigoriev 1998)

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
    /// Overlapping windows produce more data points with correct values.
    /// GGGGCCCC, window=4, step=2 → 3 windows:
    ///   [0..3] GGGG → (4−0)/(4+0) = +1.0
    ///   [2..5] GGCC → (2−2)/(2+2) = 0.0
    ///   [4..7] CCCC → (0−4)/(0+4) = −1.0
    /// </summary>
    [Test]
    public void CalculateWindowedGcSkew_OverlappingWindows()
    {
        var sequence = new DnaSequence("GGGGCCCC");

        var points = GcSkewCalculator.CalculateWindowedGcSkew(sequence, windowSize: 4, stepSize: 2).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(points, Has.Count.EqualTo(3));
            Assert.That(points[0].GcSkew, Is.EqualTo(1.0).Within(0.0001), "GGGG → +1.0");
            Assert.That(points[1].GcSkew, Is.EqualTo(0.0).Within(0.0001), "GGCC → 0.0");
            Assert.That(points[2].GcSkew, Is.EqualTo(-1.0).Within(0.0001), "CCCC → −1.0");
        });
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
    /// Evidence: Grigoriev 1998 — cumulative GC skew is sum of window skews.
    /// GGGGCCCCGGGGCCCC, window=4 → 4 windows:
    ///   GGGG → skew=+1, cumulative=+1
    ///   CCCC → skew=−1, cumulative=0
    ///   GGGG → skew=+1, cumulative=+1
    ///   CCCC → skew=−1, cumulative=0
    /// </summary>
    [Test]
    public void CalculateCumulativeGcSkew_AccumulatesCorrectly()
    {
        var sequence = new DnaSequence("GGGGCCCCGGGGCCCC");

        var points = GcSkewCalculator.CalculateCumulativeGcSkew(sequence, windowSize: 4).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(points, Has.Count.EqualTo(4));
            // Per-window GcSkew
            Assert.That(points[0].GcSkew, Is.EqualTo(1.0).Within(0.0001), "GGGG → +1");
            Assert.That(points[1].GcSkew, Is.EqualTo(-1.0).Within(0.0001), "CCCC → −1");
            Assert.That(points[2].GcSkew, Is.EqualTo(1.0).Within(0.0001), "GGGG → +1");
            Assert.That(points[3].GcSkew, Is.EqualTo(-1.0).Within(0.0001), "CCCC → −1");
            // Cumulative
            Assert.That(points[0].CumulativeGcSkew, Is.EqualTo(1.0).Within(0.0001), "cum: +1");
            Assert.That(points[1].CumulativeGcSkew, Is.EqualTo(0.0).Within(0.0001), "cum: +1−1=0");
            Assert.That(points[2].CumulativeGcSkew, Is.EqualTo(1.0).Within(0.0001), "cum: 0+1=+1");
            Assert.That(points[3].CumulativeGcSkew, Is.EqualTo(0.0).Within(0.0001), "cum: +1−1=0");
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
    /// Evidence: Grigoriev 1998 — minimum of cumulative GC skew = origin.
    /// G×50 + C×100 + G×50, window=10:
    ///   5 G-windows (+1 each) → cumulative peaks at +5
    ///   10 C-windows (−1 each) → cumulative drops to −5 at window [140..149], center=145
    ///   5 G-windows (+1 each) → cumulative returns to 0
    /// Global minimum −5 at position 145 → PredictedOrigin = 145.
    /// </summary>
    [Test]
    public void PredictReplicationOrigin_FindsMinimum()
    {
        var sequence = new DnaSequence(
            new string('G', 50) + new string('C', 100) + new string('G', 50));

        var prediction = GcSkewCalculator.PredictReplicationOrigin(sequence, windowSize: 10);

        Assert.Multiple(() =>
        {
            Assert.That(prediction.PredictedOrigin, Is.EqualTo(145), "Min cumulative at center of last C-window");
            Assert.That(prediction.OriginSkew, Is.EqualTo(-5.0).Within(0.0001), "Cumulative minimum = −5");
            Assert.That(prediction.IsSignificant, Is.True, "Amplitude 10 >> threshold 0.2");
        });
    }

    /// <summary>
    /// Evidence: Grigoriev 1998 — maximum of cumulative GC skew = terminus.
    /// C×50 + G×100 + C×50, window=10:
    ///   5 C-windows (−1 each) → cumulative drops to −5
    ///   10 G-windows (+1 each) → cumulative rises to +5 at window [140..149], center=145
    ///   5 C-windows (−1 each) → cumulative returns to 0
    /// Global maximum +5 at position 145 → PredictedTerminus = 145.
    /// </summary>
    [Test]
    public void PredictReplicationOrigin_FindsMaximum()
    {
        var sequence = new DnaSequence(
            new string('C', 50) + new string('G', 100) + new string('C', 50));

        var prediction = GcSkewCalculator.PredictReplicationOrigin(sequence, windowSize: 10);

        Assert.Multiple(() =>
        {
            Assert.That(prediction.PredictedTerminus, Is.EqualTo(145), "Max cumulative at center of last G-window");
            Assert.That(prediction.TerminusSkew, Is.EqualTo(5.0).Within(0.0001), "Cumulative maximum = +5");
            Assert.That(prediction.IsSignificant, Is.True, "Amplitude 10 >> threshold 0.2");
        });
    }

    /// <summary>
    /// Exact positions and skew values for a deterministic small genome.
    /// GGGGGGCCCCCCGGGGGGCCCCCC (24 nt), window=4, step=4:
    ///   GGGG(+1), GGCC(0), CCCC(−1), GGGG(+1), GGCC(0), CCCC(−1)
    ///   Cumulative: +1, +1, 0, +1, +1, 0
    ///   Min=0 at pos 10 (first occurrence), Max=+1 at pos 2 (first occurrence)
    /// </summary>
    [Test]
    public void PredictReplicationOrigin_ExactPositionsAndSkew()
    {
        var sequence = new DnaSequence("GGGGGGCCCCCCGGGGGGCCCCCC");

        var prediction = GcSkewCalculator.PredictReplicationOrigin(sequence, windowSize: 4);

        Assert.Multiple(() =>
        {
            Assert.That(prediction.PredictedOrigin, Is.EqualTo(10), "Min cumulative at CCCC window center");
            Assert.That(prediction.PredictedTerminus, Is.EqualTo(2), "Max cumulative at first GGGG window center");
            Assert.That(prediction.OriginSkew, Is.EqualTo(0.0).Within(0.0001));
            Assert.That(prediction.TerminusSkew, Is.EqualTo(1.0).Within(0.0001));
        });
    }

    #endregion

    #region GC Content Analysis Tests

    /// <summary>
    /// GC content correctly calculated per window.
    /// GCGCATATGCGC, window=4, step=4 → 3 windows:
    ///   GCGC → 4/4 = 100%, ATAT → 0/4 = 0%, GCGC → 4/4 = 100%
    /// </summary>
    [Test]
    public void AnalyzeGcContent_CalculatesGcPercent()
    {
        var sequence = new DnaSequence("GCGCATATGCGC");

        var result = GcSkewCalculator.AnalyzeGcContent(sequence, windowSize: 4, stepSize: 4);

        Assert.Multiple(() =>
        {
            Assert.That(result.WindowedGcContent, Has.Count.EqualTo(3));
            Assert.That(result.WindowedGcContent[0].GcContent, Is.EqualTo(100).Within(0.01), "GCGC → 100%");
            Assert.That(result.WindowedGcContent[1].GcContent, Is.EqualTo(0).Within(0.01), "ATAT → 0%");
            Assert.That(result.WindowedGcContent[2].GcContent, Is.EqualTo(100).Within(0.01), "GCGC → 100%");
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
    /// Overall metrics with exact values.
    /// GGGGCCCCATAT (12 nt): GC=8 → 8/12×100 = 66.667%
    /// GC skew: G=4, C=4 → (4−4)/(4+4) = 0.0
    /// AT skew: A=2, T=2 → (2−2)/(2+2) = 0.0
    /// </summary>
    [Test]
    public void AnalyzeGcContent_ReturnsOverallMetrics()
    {
        var sequence = new DnaSequence("GGGGCCCCATAT");

        var result = GcSkewCalculator.AnalyzeGcContent(sequence, windowSize: 4, stepSize: 4);

        Assert.Multiple(() =>
        {
            Assert.That(result.OverallGcContent, Is.EqualTo(100.0 * 8 / 12).Within(0.01), "8/12 GC");
            Assert.That(result.OverallGcSkew, Is.EqualTo(0.0).Within(0.0001), "G=C → skew=0");
            Assert.That(result.OverallAtSkew, Is.EqualTo(0.0).Within(0.0001), "A=T → skew=0");
            Assert.That(result.SequenceLength, Is.EqualTo(12));
            Assert.That(result.WindowedGcSkew, Has.Count.EqualTo(3));
            Assert.That(result.WindowedGcContent, Has.Count.EqualTo(3));
        });
    }

    #endregion

    #region Biopython Cross-Verification (Bio.SeqUtils.GC_skew)

    /// <summary>
    /// Cross-verification: Biopython GC_skew("ATGCATGC", window=4) → [0.0, 0.0]
    /// Each window ATGC has G=1, C=1, (1-1)/(1+1) = 0.
    /// </summary>
    [Test]
    public void CalculateWindowedGcSkew_BiopythonCrossVerification_ATGCATGC()
    {
        var points = GcSkewCalculator.CalculateWindowedGcSkew("ATGCATGC", windowSize: 4, stepSize: 4).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(points, Has.Count.EqualTo(2));
            Assert.That(points[0].GcSkew, Is.EqualTo(0.0).Within(0.0001), "Biopython: ATGC → 0.0");
            Assert.That(points[1].GcSkew, Is.EqualTo(0.0).Within(0.0001), "Biopython: ATGC → 0.0");
        });
    }

    /// <summary>
    /// Cross-verification: Biopython GC_skew("AAAAAAAA", window=4) → [0.0, 0.0]
    /// Source: Biopython returns 0.0 on ZeroDivisionError (G+C=0).
    /// Equivalent to Biopython: test_GC_skew("A"*50) → first element 0.
    /// </summary>
    [Test]
    public void CalculateWindowedGcSkew_BiopythonCrossVerification_AllA()
    {
        var points = GcSkewCalculator.CalculateWindowedGcSkew("AAAAAAAA", windowSize: 4, stepSize: 4).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(points, Has.Count.EqualTo(2));
            Assert.That(points[0].GcSkew, Is.EqualTo(0.0).Within(0.0001), "Biopython: AAAA → 0.0 (ZeroDivisionError)");
            Assert.That(points[1].GcSkew, Is.EqualTo(0.0).Within(0.0001), "Biopython: AAAA → 0.0 (ZeroDivisionError)");
        });
    }

    /// <summary>
    /// Cross-verification: Biopython calc_gc_skew (GenomeDiagram test helper).
    /// Formula: (g-c)/(g+c), returns 0 when g+c==0.
    /// GCCC: G=1, C=3 → (1-3)/(1+3) = -0.5
    /// </summary>
    [Test]
    public void CalculateGcSkew_BiopythonCrossVerification_GCCC()
    {
        double skew = GcSkewCalculator.CalculateGcSkew("GCCC");
        Assert.That(skew, Is.EqualTo(-0.5).Within(0.0001), "Biopython: (1-3)/(1+3) = -0.5");
    }

    #endregion

    #region Edge Cases and Exception Handling

    /// <summary>
    /// Empty string: G=0, C=0, G+C=0 → zero-division protection → 0.
    /// Source: Biopython returns 0.0 on ZeroDivisionError.
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

    /// <summary>
    /// Window size ≤ 0 for cumulative analysis throws ArgumentOutOfRangeException.
    /// </summary>
    [Test]
    public void CalculateCumulativeGcSkew_ZeroWindowSize_ThrowsException()
    {
        var sequence = new DnaSequence("ATGCATGC");

        Assert.Throws<ArgumentOutOfRangeException>(() =>
            GcSkewCalculator.CalculateCumulativeGcSkew(sequence, windowSize: 0).ToList());
    }

    /// <summary>
    /// Null DnaSequence for AT skew throws ArgumentNullException.
    /// </summary>
    [Test]
    public void CalculateAtSkew_NullSequence_ThrowsException()
    {
        Assert.Throws<ArgumentNullException>(() =>
            GcSkewCalculator.CalculateAtSkew((DnaSequence)null!));
    }

    /// <summary>
    /// Empty string: A=0, T=0, A+T=0 → zero-division protection → 0.
    /// Mirrors GC skew empty-sequence behavior.
    /// </summary>
    [Test]
    public void CalculateAtSkew_EmptySequence_ReturnsZero()
    {
        double skew = GcSkewCalculator.CalculateAtSkew("");
        Assert.That(skew, Is.EqualTo(0.0));
    }

    /// <summary>
    /// Sequence shorter than window size produces no cumulative points,
    /// so PredictReplicationOrigin returns default (zeros, not significant).
    /// </summary>
    [Test]
    public void PredictReplicationOrigin_TooShortSequence_ReturnsDefault()
    {
        var sequence = new DnaSequence("ATGC");

        var result = GcSkewCalculator.PredictReplicationOrigin(sequence, windowSize: 100);

        Assert.Multiple(() =>
        {
            Assert.That(result.PredictedOrigin, Is.EqualTo(0));
            Assert.That(result.PredictedTerminus, Is.EqualTo(0));
            Assert.That(result.IsSignificant, Is.False);
        });
    }

    #endregion
}
