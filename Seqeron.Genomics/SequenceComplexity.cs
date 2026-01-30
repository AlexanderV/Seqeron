using System;
using System.Collections.Generic;
using System.Linq;

namespace Seqeron.Genomics;

/// <summary>
/// Calculates various sequence complexity metrics for detecting low-complexity regions,
/// repetitive sequences, and information content.
/// </summary>
public static class SequenceComplexity
{
    #region Linguistic Complexity

    /// <summary>
    /// Calculates linguistic complexity (LC) as the ratio of observed to possible subwords.
    /// LC = 1.0 for maximum complexity, lower values indicate repeats/low complexity.
    /// </summary>
    /// <param name="sequence">DNA sequence.</param>
    /// <param name="maxWordLength">Maximum word length to consider.</param>
    /// <returns>Linguistic complexity (0 to 1).</returns>
    public static double CalculateLinguisticComplexity(DnaSequence sequence, int maxWordLength = 10)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        if (maxWordLength < 1) throw new ArgumentOutOfRangeException(nameof(maxWordLength));

        return CalculateLinguisticComplexityCore(sequence.Sequence, maxWordLength);
    }

    /// <summary>
    /// Calculates linguistic complexity from a raw sequence string.
    /// </summary>
    public static double CalculateLinguisticComplexity(string sequence, int maxWordLength = 10)
    {
        if (string.IsNullOrEmpty(sequence)) return 0;
        return CalculateLinguisticComplexityCore(sequence.ToUpperInvariant(), maxWordLength);
    }

    private static double CalculateLinguisticComplexityCore(string seq, int maxWordLength)
    {
        if (seq.Length == 0) return 0;

        long observedTotal = 0;
        long possibleTotal = 0;

        for (int wordLen = 1; wordLen <= Math.Min(maxWordLength, seq.Length); wordLen++)
        {
            var observedWords = new HashSet<string>();

            for (int i = 0; i <= seq.Length - wordLen; i++)
            {
                observedWords.Add(seq.Substring(i, wordLen));
            }

            observedTotal += observedWords.Count;

            // Maximum possible words of length wordLen
            long maxPossible = Math.Min(
                (long)Math.Pow(4, wordLen),       // 4^wordLen possible DNA words
                seq.Length - wordLen + 1);         // Can't observe more than available positions

            possibleTotal += maxPossible;
        }

        return possibleTotal > 0 ? (double)observedTotal / possibleTotal : 0;
    }

    #endregion

    #region Shannon Entropy

    /// <summary>
    /// Calculates Shannon entropy for the sequence (bits per base).
    /// Maximum entropy for DNA is 2 bits (log2(4)).
    /// </summary>
    /// <param name="sequence">DNA sequence.</param>
    /// <returns>Shannon entropy (0 to 2 for DNA).</returns>
    public static double CalculateShannonEntropy(DnaSequence sequence)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        return CalculateShannonEntropyCore(sequence.Sequence);
    }

    /// <summary>
    /// Calculates Shannon entropy from a raw sequence string.
    /// </summary>
    public static double CalculateShannonEntropy(string sequence)
    {
        if (string.IsNullOrEmpty(sequence)) return 0;
        return CalculateShannonEntropyCore(sequence.ToUpperInvariant());
    }

    private static double CalculateShannonEntropyCore(string seq)
    {
        if (seq.Length == 0) return 0;

        var frequencies = new Dictionary<char, int> { ['A'] = 0, ['T'] = 0, ['G'] = 0, ['C'] = 0 };

        foreach (char c in seq)
        {
            if (frequencies.ContainsKey(c))
                frequencies[c]++;
        }

        double entropy = 0;
        int total = frequencies.Values.Sum();

        if (total == 0) return 0;

        foreach (int count in frequencies.Values)
        {
            if (count > 0)
            {
                double p = (double)count / total;
                entropy -= p * Math.Log2(p);
            }
        }

        return entropy;
    }

    /// <summary>
    /// Calculates Shannon entropy using k-mers.
    /// </summary>
    /// <param name="sequence">DNA sequence.</param>
    /// <param name="k">K-mer size (default: 2 for dinucleotides).</param>
    /// <returns>Shannon entropy based on k-mer frequencies.</returns>
    public static double CalculateKmerEntropy(DnaSequence sequence, int k = 2)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        if (k < 1) throw new ArgumentOutOfRangeException(nameof(k));

        return CalculateKmerEntropyCore(sequence.Sequence, k);
    }

    private static double CalculateKmerEntropyCore(string seq, int k)
    {
        if (seq.Length < k) return 0;

        var kmerCounts = new Dictionary<string, int>();
        int total = 0;

        for (int i = 0; i <= seq.Length - k; i++)
        {
            string kmer = seq.Substring(i, k);
            if (kmerCounts.ContainsKey(kmer))
                kmerCounts[kmer]++;
            else
                kmerCounts[kmer] = 1;
            total++;
        }

        double entropy = 0;
        foreach (int count in kmerCounts.Values)
        {
            double p = (double)count / total;
            entropy -= p * Math.Log2(p);
        }

        return entropy;
    }

    #endregion

    #region Sliding Window Complexity

    /// <summary>
    /// Calculates complexity across the sequence using a sliding window.
    /// </summary>
    /// <param name="sequence">DNA sequence.</param>
    /// <param name="windowSize">Size of the sliding window (default: 64).</param>
    /// <param name="stepSize">Step size for window movement (default: 10).</param>
    /// <returns>Complexity values with positions.</returns>
    public static IEnumerable<ComplexityPoint> CalculateWindowedComplexity(
        DnaSequence sequence,
        int windowSize = 64,
        int stepSize = 10)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        if (windowSize < 1) throw new ArgumentOutOfRangeException(nameof(windowSize));
        if (stepSize < 1) throw new ArgumentOutOfRangeException(nameof(stepSize));

        return CalculateWindowedComplexityCore(sequence.Sequence, windowSize, stepSize);
    }

    private static IEnumerable<ComplexityPoint> CalculateWindowedComplexityCore(
        string seq,
        int windowSize,
        int stepSize)
    {
        for (int i = 0; i + windowSize <= seq.Length; i += stepSize)
        {
            string window = seq.Substring(i, windowSize);
            double entropy = CalculateShannonEntropyCore(window);
            double lc = CalculateLinguisticComplexityCore(window, Math.Min(6, windowSize));

            yield return new ComplexityPoint(
                Position: i + windowSize / 2,
                ShannonEntropy: entropy,
                LinguisticComplexity: lc,
                WindowStart: i,
                WindowEnd: i + windowSize - 1);
        }
    }

    #endregion

    #region Low Complexity Regions

    /// <summary>
    /// Finds low-complexity regions in the sequence.
    /// Uses a combination of entropy and linguistic complexity.
    /// </summary>
    /// <param name="sequence">DNA sequence.</param>
    /// <param name="windowSize">Window size for analysis (default: 64).</param>
    /// <param name="entropyThreshold">Entropy threshold below which regions are considered low complexity (default: 1.0).</param>
    /// <returns>Low-complexity regions.</returns>
    public static IEnumerable<LowComplexityRegion> FindLowComplexityRegions(
        DnaSequence sequence,
        int windowSize = 64,
        double entropyThreshold = 1.0)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        if (windowSize < 1) throw new ArgumentOutOfRangeException(nameof(windowSize));

        return FindLowComplexityRegionsCore(sequence.Sequence, windowSize, entropyThreshold);
    }

    private static IEnumerable<LowComplexityRegion> FindLowComplexityRegionsCore(
        string seq,
        int windowSize,
        double entropyThreshold)
    {
        if (seq.Length < windowSize) yield break;

        int? regionStart = null;
        double minEntropy = double.MaxValue;

        for (int i = 0; i + windowSize <= seq.Length; i++)
        {
            string window = seq.Substring(i, windowSize);
            double entropy = CalculateShannonEntropyCore(window);

            if (entropy < entropyThreshold)
            {
                if (regionStart == null)
                {
                    regionStart = i;
                    minEntropy = entropy;
                }
                else
                {
                    minEntropy = Math.Min(minEntropy, entropy);
                }
            }
            else if (regionStart != null)
            {
                // End of low-complexity region
                int end = i + windowSize - 1;
                yield return new LowComplexityRegion(
                    Start: regionStart.Value,
                    End: end,
                    Length: end - regionStart.Value + 1,
                    MinEntropy: minEntropy,
                    Sequence: seq.Substring(regionStart.Value, end - regionStart.Value + 1));

                regionStart = null;
                minEntropy = double.MaxValue;
            }
        }

        // Handle region at end of sequence
        if (regionStart != null)
        {
            int end = seq.Length - 1;
            yield return new LowComplexityRegion(
                Start: regionStart.Value,
                End: end,
                Length: end - regionStart.Value + 1,
                MinEntropy: minEntropy,
                Sequence: seq.Substring(regionStart.Value));
        }
    }

    #endregion

    #region Dust Score

    /// <summary>
    /// Calculates DUST score for low-complexity filtering (as used in BLAST).
    /// Higher scores indicate lower complexity.
    /// </summary>
    /// <param name="sequence">DNA sequence.</param>
    /// <param name="wordSize">Triplet word size (default: 3).</param>
    /// <returns>DUST score.</returns>
    public static double CalculateDustScore(DnaSequence sequence, int wordSize = 3)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        return CalculateDustScoreCore(sequence.Sequence, wordSize);
    }

    /// <summary>
    /// Calculates DUST score from a raw sequence string.
    /// </summary>
    public static double CalculateDustScore(string sequence, int wordSize = 3)
    {
        if (string.IsNullOrEmpty(sequence)) return 0;
        return CalculateDustScoreCore(sequence.ToUpperInvariant(), wordSize);
    }

    private static double CalculateDustScoreCore(string seq, int wordSize)
    {
        if (seq.Length < wordSize) return 0;

        var tripletCounts = new Dictionary<string, int>();
        int total = seq.Length - wordSize + 1;

        for (int i = 0; i <= seq.Length - wordSize; i++)
        {
            string triplet = seq.Substring(i, wordSize);
            if (tripletCounts.ContainsKey(triplet))
                tripletCounts[triplet]++;
            else
                tripletCounts[triplet] = 1;
        }

        // DUST score is sum of (count * (count - 1) / 2) for each triplet
        double score = 0;
        foreach (int count in tripletCounts.Values)
        {
            score += count * (count - 1) / 2.0;
        }

        // Normalize by window length
        return total > 1 ? score / (total - 1) : 0;
    }

    /// <summary>
    /// Masks low-complexity regions using DUST algorithm.
    /// </summary>
    /// <param name="sequence">DNA sequence.</param>
    /// <param name="windowSize">Window size for masking (default: 64).</param>
    /// <param name="threshold">DUST threshold above which to mask (default: 2.0).</param>
    /// <param name="maskChar">Character to use for masking (default: 'N').</param>
    /// <returns>Masked sequence.</returns>
    public static string MaskLowComplexity(
        DnaSequence sequence,
        int windowSize = 64,
        double threshold = 2.0,
        char maskChar = 'N')
    {
        ArgumentNullException.ThrowIfNull(sequence);

        return MaskLowComplexityCore(sequence.Sequence, windowSize, threshold, maskChar);
    }

    private static string MaskLowComplexityCore(string seq, int windowSize, double threshold, char maskChar)
    {
        if (seq.Length < windowSize) return seq;

        var masked = new char[seq.Length];
        seq.CopyTo(0, masked, 0, seq.Length);

        for (int i = 0; i + windowSize <= seq.Length; i++)
        {
            string window = seq.Substring(i, windowSize);
            double dustScore = CalculateDustScoreCore(window, 3);

            if (dustScore > threshold)
            {
                for (int j = i; j < i + windowSize; j++)
                {
                    masked[j] = maskChar;
                }
            }
        }

        return new string(masked);
    }

    #endregion

    #region Compression Ratio

    /// <summary>
    /// Estimates sequence complexity using compression ratio.
    /// Lower ratios indicate more repetitive/less complex sequences.
    /// </summary>
    /// <param name="sequence">DNA sequence.</param>
    /// <returns>Estimated compression ratio (0 to 1).</returns>
    public static double EstimateCompressionRatio(DnaSequence sequence)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        return EstimateCompressionRatioCore(sequence.Sequence);
    }

    /// <summary>
    /// Estimates compression ratio from a raw sequence string.
    /// </summary>
    public static double EstimateCompressionRatio(string sequence)
    {
        if (string.IsNullOrEmpty(sequence)) return 0;
        return EstimateCompressionRatioCore(sequence.ToUpperInvariant());
    }

    private static double EstimateCompressionRatioCore(string seq)
    {
        if (seq.Length == 0) return 0;

        // Use LZ77-like approach: count unique substrings
        var seen = new HashSet<string>();
        int uniqueCount = 0;

        for (int len = 1; len <= Math.Min(10, seq.Length); len++)
        {
            for (int i = 0; i <= seq.Length - len; i++)
            {
                string sub = seq.Substring(i, len);
                if (!seen.Contains(sub))
                {
                    seen.Add(sub);
                    uniqueCount++;
                }
            }
        }

        // Calculate expected unique substrings for random sequence
        double expected = 0;
        for (int len = 1; len <= Math.Min(10, seq.Length); len++)
        {
            expected += Math.Min(Math.Pow(4, len), seq.Length - len + 1);
        }

        return expected > 0 ? (double)uniqueCount / expected : 0;
    }

    #endregion
}

/// <summary>
/// A point in complexity analysis.
/// </summary>
public readonly record struct ComplexityPoint(
    int Position,
    double ShannonEntropy,
    double LinguisticComplexity,
    int WindowStart,
    int WindowEnd);

/// <summary>
/// A low-complexity region detected in the sequence.
/// </summary>
public readonly record struct LowComplexityRegion(
    int Start,
    int End,
    int Length,
    double MinEntropy,
    string Sequence);
