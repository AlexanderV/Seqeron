using System;
using System.Collections.Generic;
using System.Linq;

namespace Seqeron.Genomics;

/// <summary>
/// Calculates GC skew and related metrics for identifying replication origins and termini.
/// GC skew = (G - C) / (G + C), useful for finding origin of replication in bacterial genomes.
/// </summary>
public static class GcSkewCalculator
{
    #region GC Skew Calculation

    /// <summary>
    /// Calculates GC skew for a single sequence or window.
    /// GC skew = (G - C) / (G + C).
    /// </summary>
    /// <param name="sequence">DNA sequence.</param>
    /// <returns>GC skew value (-1 to 1).</returns>
    public static double CalculateGcSkew(DnaSequence sequence)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        return CalculateGcSkewCore(sequence.Sequence);
    }

    /// <summary>
    /// Calculates GC skew from a raw sequence string.
    /// </summary>
    public static double CalculateGcSkew(string sequence)
    {
        if (string.IsNullOrEmpty(sequence))
            return 0;

        return CalculateGcSkewCore(sequence.ToUpperInvariant());
    }

    private static double CalculateGcSkewCore(string seq)
    {
        int gCount = seq.Count(c => c == 'G');
        int cCount = seq.Count(c => c == 'C');
        int total = gCount + cCount;

        return total > 0 ? (double)(gCount - cCount) / total : 0;
    }

    #endregion

    #region Sliding Window GC Skew

    /// <summary>
    /// Calculates GC skew using a sliding window across the sequence.
    /// </summary>
    /// <param name="sequence">DNA sequence.</param>
    /// <param name="windowSize">Size of the sliding window (default: 1000).</param>
    /// <param name="stepSize">Step size for window movement (default: 100).</param>
    /// <returns>Collection of GC skew values with positions.</returns>
    public static IEnumerable<GcSkewPoint> CalculateWindowedGcSkew(
        DnaSequence sequence,
        int windowSize = 1000,
        int stepSize = 100)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        if (windowSize < 1) throw new ArgumentOutOfRangeException(nameof(windowSize));
        if (stepSize < 1) throw new ArgumentOutOfRangeException(nameof(stepSize));

        return CalculateWindowedGcSkewCore(sequence.Sequence, windowSize, stepSize);
    }

    /// <summary>
    /// Calculates windowed GC skew from a raw sequence string.
    /// </summary>
    public static IEnumerable<GcSkewPoint> CalculateWindowedGcSkew(
        string sequence,
        int windowSize = 1000,
        int stepSize = 100)
    {
        if (string.IsNullOrEmpty(sequence))
            yield break;

        foreach (var point in CalculateWindowedGcSkewCore(sequence.ToUpperInvariant(), windowSize, stepSize))
            yield return point;
    }

    private static IEnumerable<GcSkewPoint> CalculateWindowedGcSkewCore(
        string seq,
        int windowSize,
        int stepSize)
    {
        for (int i = 0; i + windowSize <= seq.Length; i += stepSize)
        {
            string window = seq.Substring(i, windowSize);
            double skew = CalculateGcSkewCore(window);

            yield return new GcSkewPoint(
                Position: i + windowSize / 2,
                GcSkew: skew,
                WindowStart: i,
                WindowEnd: i + windowSize - 1);
        }
    }

    #endregion

    #region Cumulative GC Skew

    /// <summary>
    /// Calculates cumulative GC skew across the sequence.
    /// Useful for identifying origin and terminus of replication.
    /// Minimum = origin of replication, Maximum = terminus.
    /// </summary>
    /// <param name="sequence">DNA sequence.</param>
    /// <param name="windowSize">Size of the window for cumulative calculation (default: 1000).</param>
    /// <returns>Collection of cumulative GC skew values.</returns>
    public static IEnumerable<CumulativeGcSkewPoint> CalculateCumulativeGcSkew(
        DnaSequence sequence,
        int windowSize = 1000)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        if (windowSize < 1) throw new ArgumentOutOfRangeException(nameof(windowSize));

        return CalculateCumulativeGcSkewCore(sequence.Sequence, windowSize);
    }

    /// <summary>
    /// Calculates cumulative GC skew from a raw sequence string.
    /// </summary>
    public static IEnumerable<CumulativeGcSkewPoint> CalculateCumulativeGcSkew(
        string sequence,
        int windowSize = 1000)
    {
        if (string.IsNullOrEmpty(sequence))
            yield break;

        foreach (var point in CalculateCumulativeGcSkewCore(sequence.ToUpperInvariant(), windowSize))
            yield return point;
    }

    private static IEnumerable<CumulativeGcSkewPoint> CalculateCumulativeGcSkewCore(
        string seq,
        int windowSize)
    {
        double cumulative = 0;
        int stepSize = windowSize;

        for (int i = 0; i + windowSize <= seq.Length; i += stepSize)
        {
            string window = seq.Substring(i, windowSize);
            double skew = CalculateGcSkewCore(window);
            cumulative += skew;

            yield return new CumulativeGcSkewPoint(
                Position: i + windowSize / 2,
                GcSkew: skew,
                CumulativeGcSkew: cumulative);
        }
    }

    #endregion

    #region AT Skew Calculation

    /// <summary>
    /// Calculates AT skew for a sequence.
    /// AT skew = (A - T) / (A + T).
    /// </summary>
    public static double CalculateAtSkew(DnaSequence sequence)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        return CalculateAtSkewCore(sequence.Sequence);
    }

    /// <summary>
    /// Calculates AT skew from a raw sequence string.
    /// </summary>
    public static double CalculateAtSkew(string sequence)
    {
        if (string.IsNullOrEmpty(sequence))
            return 0;

        return CalculateAtSkewCore(sequence.ToUpperInvariant());
    }

    private static double CalculateAtSkewCore(string seq)
    {
        int aCount = seq.Count(c => c == 'A');
        int tCount = seq.Count(c => c == 'T');
        int total = aCount + tCount;

        return total > 0 ? (double)(aCount - tCount) / total : 0;
    }

    #endregion

    #region Origin/Terminus Prediction

    /// <summary>
    /// Predicts the origin and terminus of replication based on cumulative GC skew.
    /// </summary>
    /// <param name="sequence">DNA sequence (should be complete circular genome).</param>
    /// <param name="windowSize">Window size for analysis (default: 1000).</param>
    /// <returns>Predicted origin and terminus positions.</returns>
    public static ReplicationOriginPrediction PredictReplicationOrigin(
        DnaSequence sequence,
        int windowSize = 1000)
    {
        ArgumentNullException.ThrowIfNull(sequence);

        var cumulativePoints = CalculateCumulativeGcSkewCore(sequence.Sequence, windowSize).ToList();

        if (cumulativePoints.Count == 0)
        {
            return new ReplicationOriginPrediction(0, 0, 0, 0, false);
        }

        // Find minimum (origin) and maximum (terminus)
        var minPoint = cumulativePoints.MinBy(p => p.CumulativeGcSkew);
        var maxPoint = cumulativePoints.MaxBy(p => p.CumulativeGcSkew);

        double skewAmplitude = maxPoint!.CumulativeGcSkew - minPoint!.CumulativeGcSkew;
        bool isSignificant = Math.Abs(skewAmplitude) > cumulativePoints.Count * 0.01;

        return new ReplicationOriginPrediction(
            PredictedOrigin: minPoint.Position,
            PredictedTerminus: maxPoint.Position,
            OriginSkew: minPoint.CumulativeGcSkew,
            TerminusSkew: maxPoint.CumulativeGcSkew,
            IsSignificant: isSignificant);
    }

    #endregion

    #region Comprehensive GC Analysis

    /// <summary>
    /// Gets comprehensive GC analysis including skew, content, and variability.
    /// </summary>
    public static GcAnalysisResult AnalyzeGcContent(
        DnaSequence sequence,
        int windowSize = 1000,
        int stepSize = 100)
    {
        ArgumentNullException.ThrowIfNull(sequence);

        var windowedSkew = CalculateWindowedGcSkewCore(sequence.Sequence, windowSize, stepSize).ToList();
        var windowedContent = CalculateWindowedGcContentCore(sequence.Sequence, windowSize, stepSize).ToList();

        double overallGcContent = CalculateGcContent(sequence.Sequence);
        double overallGcSkew = CalculateGcSkewCore(sequence.Sequence);
        double overallAtSkew = CalculateAtSkewCore(sequence.Sequence);

        double gcContentVariance = windowedContent.Count > 0
            ? CalculateVariance(windowedContent.Select(w => w.GcContent).ToList())
            : 0;

        double gcSkewVariance = windowedSkew.Count > 0
            ? CalculateVariance(windowedSkew.Select(w => w.GcSkew).ToList())
            : 0;

        return new GcAnalysisResult(
            OverallGcContent: overallGcContent,
            OverallGcSkew: overallGcSkew,
            OverallAtSkew: overallAtSkew,
            GcContentVariance: gcContentVariance,
            GcSkewVariance: gcSkewVariance,
            WindowedGcSkew: windowedSkew,
            WindowedGcContent: windowedContent,
            SequenceLength: sequence.Length);
    }

    private static IEnumerable<GcContentPoint> CalculateWindowedGcContentCore(
        string seq,
        int windowSize,
        int stepSize)
    {
        for (int i = 0; i + windowSize <= seq.Length; i += stepSize)
        {
            string window = seq.Substring(i, windowSize);
            double gcContent = CalculateGcContent(window);

            yield return new GcContentPoint(
                Position: i + windowSize / 2,
                GcContent: gcContent,
                WindowStart: i,
                WindowEnd: i + windowSize - 1);
        }
    }

    private static double CalculateGcContent(string seq)
    {
        if (string.IsNullOrEmpty(seq)) return 0;
        int gcCount = seq.Count(c => c is 'G' or 'C');
        return (double)gcCount / seq.Length * 100;
    }

    private static double CalculateVariance(IList<double> values)
    {
        if (values.Count == 0) return 0;
        double mean = values.Average();
        return values.Sum(v => (v - mean) * (v - mean)) / values.Count;
    }

    #endregion
}

/// <summary>
/// A point in GC skew analysis with position and value.
/// </summary>
public readonly record struct GcSkewPoint(
    int Position,
    double GcSkew,
    int WindowStart,
    int WindowEnd);

/// <summary>
/// A point in cumulative GC skew analysis.
/// </summary>
public readonly record struct CumulativeGcSkewPoint(
    int Position,
    double GcSkew,
    double CumulativeGcSkew);

/// <summary>
/// A point in GC content analysis.
/// </summary>
public readonly record struct GcContentPoint(
    int Position,
    double GcContent,
    int WindowStart,
    int WindowEnd);

/// <summary>
/// Predicted origin and terminus of replication.
/// </summary>
public readonly record struct ReplicationOriginPrediction(
    int PredictedOrigin,
    int PredictedTerminus,
    double OriginSkew,
    double TerminusSkew,
    bool IsSignificant);

/// <summary>
/// Comprehensive GC analysis results.
/// </summary>
public sealed record GcAnalysisResult(
    double OverallGcContent,
    double OverallGcSkew,
    double OverallAtSkew,
    double GcContentVariance,
    double GcSkewVariance,
    IReadOnlyList<GcSkewPoint> WindowedGcSkew,
    IReadOnlyList<GcContentPoint> WindowedGcContent,
    int SequenceLength);
