using NUnit.Framework;
using SuffixTree.Genomics;
using System.Linq;

namespace SuffixTree.Genomics.Tests;

[TestFixture]
public class GcSkewCalculatorTests
{
    #region Basic GC Skew Tests

    [Test]
    public void CalculateGcSkew_MoreG_ReturnsPositive()
    {
        var sequence = new DnaSequence("GGGGC"); // 4G, 1C = (4-1)/(4+1) = 0.6
        double skew = GcSkewCalculator.CalculateGcSkew(sequence);

        Assert.That(skew, Is.EqualTo(0.6).Within(0.01));
    }

    [Test]
    public void CalculateGcSkew_MoreC_ReturnsNegative()
    {
        var sequence = new DnaSequence("CCCCC"); // 0G, 5C = (0-5)/(0+5) = -1.0
        double skew = GcSkewCalculator.CalculateGcSkew(sequence);

        Assert.That(skew, Is.EqualTo(-1.0).Within(0.01));
    }

    [Test]
    public void CalculateGcSkew_EqualGC_ReturnsZero()
    {
        var sequence = new DnaSequence("GCGC"); // 2G, 2C = (2-2)/(2+2) = 0
        double skew = GcSkewCalculator.CalculateGcSkew(sequence);

        Assert.That(skew, Is.EqualTo(0).Within(0.01));
    }

    [Test]
    public void CalculateGcSkew_NoGC_ReturnsZero()
    {
        var sequence = new DnaSequence("AAATTT"); // No G or C
        double skew = GcSkewCalculator.CalculateGcSkew(sequence);

        Assert.That(skew, Is.EqualTo(0));
    }

    [Test]
    public void CalculateGcSkew_StringOverload_Works()
    {
        double skew = GcSkewCalculator.CalculateGcSkew("GGGGC");
        Assert.That(skew, Is.EqualTo(0.6).Within(0.01));
    }

    [Test]
    public void CalculateGcSkew_RangeIsMinusOneToOne()
    {
        var sequence = new DnaSequence("ATGCATGCATGC");
        double skew = GcSkewCalculator.CalculateGcSkew(sequence);

        Assert.That(skew, Is.GreaterThanOrEqualTo(-1));
        Assert.That(skew, Is.LessThanOrEqualTo(1));
    }

    #endregion

    #region Windowed GC Skew Tests

    [Test]
    public void CalculateWindowedGcSkew_ReturnsMultiplePoints()
    {
        var sequence = new DnaSequence("GGGGCCCC" + "CCCCGGGG"); // Two windows with opposite skew
        var points = GcSkewCalculator.CalculateWindowedGcSkew(sequence, windowSize: 8, stepSize: 8).ToList();

        Assert.That(points.Count, Is.EqualTo(2));
    }

    [Test]
    public void CalculateWindowedGcSkew_CorrectPositions()
    {
        var sequence = new DnaSequence("GGGGCCCCGGGGCCCC");
        var points = GcSkewCalculator.CalculateWindowedGcSkew(sequence, windowSize: 4, stepSize: 4).ToList();

        // Positions are at the center of each window
        Assert.That(points[0].Position, Is.EqualTo(2));  // Window 0-3, center at 2
        Assert.That(points[1].Position, Is.EqualTo(6));  // Window 4-7, center at 6
        Assert.That(points[2].Position, Is.EqualTo(10)); // Window 8-11, center at 10
    }

    [Test]
    public void CalculateWindowedGcSkew_CorrectSkewValues()
    {
        var sequence = new DnaSequence("GGGGCCCC"); // First: all G, Second: all C
        var points = GcSkewCalculator.CalculateWindowedGcSkew(sequence, windowSize: 4, stepSize: 4).ToList();

        Assert.That(points[0].GcSkew, Is.EqualTo(1.0).Within(0.01)); // GGGG = 1.0
        Assert.That(points[1].GcSkew, Is.EqualTo(-1.0).Within(0.01)); // CCCC = -1.0
    }

    [Test]
    public void CalculateWindowedGcSkew_OverlappingWindows()
    {
        var sequence = new DnaSequence("GGGGCCCC");
        var points = GcSkewCalculator.CalculateWindowedGcSkew(sequence, windowSize: 4, stepSize: 2).ToList();

        Assert.That(points.Count, Is.GreaterThan(2)); // More points due to overlap
    }

    [Test]
    public void CalculateWindowedGcSkew_EmptySequence_ReturnsEmpty()
    {
        var points = GcSkewCalculator.CalculateWindowedGcSkew("", windowSize: 10).ToList();
        Assert.That(points, Is.Empty);
    }

    [Test]
    public void CalculateWindowedGcSkew_SequenceShorterThanWindow_ReturnsEmpty()
    {
        var points = GcSkewCalculator.CalculateWindowedGcSkew("ATGC", windowSize: 10).ToList();
        Assert.That(points, Is.Empty);
    }

    #endregion

    #region Cumulative GC Skew Tests

    [Test]
    public void CalculateCumulativeGcSkew_AccumulatesCorrectly()
    {
        // Create sequence with alternating high G and high C windows
        var sequence = new DnaSequence("GGGGCCCCGGGGCCCC"); // Windows: GGGG(+1), CCCC(-1), GGGG(+1), CCCC(-1)
        var points = GcSkewCalculator.CalculateCumulativeGcSkew(sequence, windowSize: 4).ToList();

        Assert.That(points.Count, Is.EqualTo(4));
        Assert.That(points[0].CumulativeGcSkew, Is.EqualTo(1.0).Within(0.01));  // +1
        Assert.That(points[1].CumulativeGcSkew, Is.EqualTo(0.0).Within(0.01));  // +1 + (-1) = 0
        Assert.That(points[2].CumulativeGcSkew, Is.EqualTo(1.0).Within(0.01));  // +1 + (-1) + 1 = 1
        Assert.That(points[3].CumulativeGcSkew, Is.EqualTo(0.0).Within(0.01));  // +1 + (-1) + 1 + (-1) = 0
    }

    [Test]
    public void CalculateCumulativeGcSkew_PositionsAreCorrect()
    {
        var sequence = new DnaSequence("GGGGCCCCGGGG");
        var points = GcSkewCalculator.CalculateCumulativeGcSkew(sequence, windowSize: 4).ToList();

        // Positions should be at center of each window
        Assert.That(points[0].Position, Is.EqualTo(2)); // Window 0-3, center at 2
        Assert.That(points[1].Position, Is.EqualTo(6)); // Window 4-7, center at 6
        Assert.That(points[2].Position, Is.EqualTo(10)); // Window 8-11, center at 10
    }

    [Test]
    public void CalculateCumulativeGcSkew_IncludesWindowSkew()
    {
        var sequence = new DnaSequence("GGGGCCCC");
        var points = GcSkewCalculator.CalculateCumulativeGcSkew(sequence, windowSize: 4).ToList();

        Assert.That(points[0].GcSkew, Is.EqualTo(1.0).Within(0.01));  // GGGG
        Assert.That(points[1].GcSkew, Is.EqualTo(-1.0).Within(0.01)); // CCCC
    }

    [Test]
    public void CalculateCumulativeGcSkew_EmptySequence_ReturnsEmpty()
    {
        var points = GcSkewCalculator.CalculateCumulativeGcSkew("").ToList();
        Assert.That(points, Is.Empty);
    }

    #endregion

    #region AT Skew Tests

    [Test]
    public void CalculateAtSkew_MoreA_ReturnsPositive()
    {
        var sequence = new DnaSequence("AAAAT"); // 4A, 1T = (4-1)/(4+1) = 0.6
        double skew = GcSkewCalculator.CalculateAtSkew(sequence);

        Assert.That(skew, Is.EqualTo(0.6).Within(0.01));
    }

    [Test]
    public void CalculateAtSkew_MoreT_ReturnsNegative()
    {
        var sequence = new DnaSequence("TTTTT"); // 0A, 5T = (0-5)/(0+5) = -1.0
        double skew = GcSkewCalculator.CalculateAtSkew(sequence);

        Assert.That(skew, Is.EqualTo(-1.0).Within(0.01));
    }

    [Test]
    public void CalculateAtSkew_EqualAT_ReturnsZero()
    {
        var sequence = new DnaSequence("ATAT");
        double skew = GcSkewCalculator.CalculateAtSkew(sequence);

        Assert.That(skew, Is.EqualTo(0).Within(0.01));
    }

    [Test]
    public void CalculateAtSkew_NoAT_ReturnsZero()
    {
        var sequence = new DnaSequence("GCGC");
        double skew = GcSkewCalculator.CalculateAtSkew(sequence);

        Assert.That(skew, Is.EqualTo(0));
    }

    #endregion

    #region Replication Origin Prediction Tests

    [Test]
    public void PredictReplicationOrigin_FindsMinimum()
    {
        // Simulate bacterial genome pattern: GC skew decreases to origin, then increases
        var sequence = new DnaSequence(
            new string('G', 50) + new string('C', 100) + new string('G', 50)); // Min at position ~100

        var prediction = GcSkewCalculator.PredictReplicationOrigin(sequence, windowSize: 10);

        Assert.That(prediction.PredictedOrigin, Is.InRange(40, 160));
    }

    [Test]
    public void PredictReplicationOrigin_FindsMaximum()
    {
        // Terminus should be at maximum cumulative skew
        var sequence = new DnaSequence(
            new string('C', 50) + new string('G', 100) + new string('C', 50));

        var prediction = GcSkewCalculator.PredictReplicationOrigin(sequence, windowSize: 10);

        Assert.That(prediction.PredictedTerminus, Is.InRange(40, 160));
    }

    [Test]
    public void PredictReplicationOrigin_ReturnsValidPositions()
    {
        var sequence = new DnaSequence("GGGGGGCCCCCCGGGGGGCCCCCC");
        var prediction = GcSkewCalculator.PredictReplicationOrigin(sequence, windowSize: 4);

        Assert.That(prediction.PredictedOrigin, Is.GreaterThanOrEqualTo(0));
        Assert.That(prediction.PredictedOrigin, Is.LessThan(sequence.Length));
        Assert.That(prediction.PredictedTerminus, Is.GreaterThanOrEqualTo(0));
        Assert.That(prediction.PredictedTerminus, Is.LessThan(sequence.Length));
    }

    #endregion

    #region GC Content Analysis Tests

    [Test]
    public void AnalyzeGcContent_ReturnsMultipleWindowedPoints()
    {
        var sequence = new DnaSequence("GCGCATATGCGC");
        var result = GcSkewCalculator.AnalyzeGcContent(sequence, windowSize: 4, stepSize: 4);

        Assert.That(result.WindowedGcContent.Count, Is.EqualTo(3));
    }

    [Test]
    public void AnalyzeGcContent_CalculatesGcPercent()
    {
        var sequence = new DnaSequence("GCGCATAT"); // First window: 100% GC, Second: 0% GC
        var result = GcSkewCalculator.AnalyzeGcContent(sequence, windowSize: 4, stepSize: 4);

        Assert.That(result.WindowedGcContent[0].GcContent, Is.EqualTo(100).Within(0.01));
        Assert.That(result.WindowedGcContent[1].GcContent, Is.EqualTo(0).Within(0.01));
    }

    [Test]
    public void AnalyzeGcContent_50PercentGc()
    {
        var sequence = new DnaSequence("ATGC");
        var result = GcSkewCalculator.AnalyzeGcContent(sequence, windowSize: 4, stepSize: 4);

        Assert.That(result.OverallGcContent, Is.EqualTo(50).Within(0.01));
    }

    [Test]
    public void AnalyzeGcContent_ReturnsOverallMetrics()
    {
        var sequence = new DnaSequence("GGGGCCCCATAT");
        var result = GcSkewCalculator.AnalyzeGcContent(sequence, windowSize: 4, stepSize: 4);

        Assert.That(result.OverallGcContent, Is.GreaterThan(0));
        Assert.That(result.SequenceLength, Is.EqualTo(12));
    }

    #endregion

    #region Edge Cases

    [Test]
    public void CalculateGcSkew_EmptySequence_ReturnsZero()
    {
        double skew = GcSkewCalculator.CalculateGcSkew("");
        Assert.That(skew, Is.EqualTo(0));
    }

    [Test]
    public void CalculateGcSkew_NullSequence_ThrowsException()
    {
        Assert.Throws<ArgumentNullException>(() =>
            GcSkewCalculator.CalculateGcSkew((DnaSequence)null!));
    }

    [Test]
    public void CalculateWindowedGcSkew_NullSequence_ThrowsException()
    {
        Assert.Throws<ArgumentNullException>(() =>
            GcSkewCalculator.CalculateWindowedGcSkew((DnaSequence)null!, 10).ToList());
    }

    [Test]
    public void CalculateWindowedGcSkew_ZeroWindowSize_ThrowsException()
    {
        var sequence = new DnaSequence("ATGC");
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            GcSkewCalculator.CalculateWindowedGcSkew(sequence, windowSize: 0).ToList());
    }

    [Test]
    public void CalculateWindowedGcSkew_NegativeStepSize_ThrowsException()
    {
        var sequence = new DnaSequence("ATGC");
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            GcSkewCalculator.CalculateWindowedGcSkew(sequence, windowSize: 2, stepSize: -1).ToList());
    }

    [Test]
    public void CalculateCumulativeGcSkew_NullSequence_ThrowsException()
    {
        Assert.Throws<ArgumentNullException>(() =>
            GcSkewCalculator.CalculateCumulativeGcSkew((DnaSequence)null!).ToList());
    }

    [Test]
    public void PredictReplicationOrigin_NullSequence_ThrowsException()
    {
        Assert.Throws<ArgumentNullException>(() =>
            GcSkewCalculator.PredictReplicationOrigin((DnaSequence)null!, 10));
    }

    [Test]
    public void AnalyzeGcContent_NullSequence_ThrowsException()
    {
        Assert.Throws<ArgumentNullException>(() =>
            GcSkewCalculator.AnalyzeGcContent((DnaSequence)null!, 10));
    }

    #endregion
}
