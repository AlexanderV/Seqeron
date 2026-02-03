using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading;
using System.Threading.Tasks;

namespace Seqeron.Genomics.Analysis;

/// <summary>
/// Analyzes k-mers (substrings of length k) in DNA/RNA sequences.
/// </summary>
public static class KmerAnalyzer
{
    /// <summary>
    /// Counts all k-mers in a sequence.
    /// </summary>
    /// <param name="sequence">The sequence to analyze.</param>
    /// <param name="k">The k-mer length.</param>
    /// <returns>Dictionary mapping k-mers to their counts.</returns>
    public static Dictionary<string, int> CountKmers(string sequence, int k)
    {
        if (string.IsNullOrEmpty(sequence))
            return new Dictionary<string, int>();

        if (k <= 0)
            throw new ArgumentOutOfRangeException(nameof(k), "K must be positive.");

        if (k > sequence.Length)
            return new Dictionary<string, int>();

        var seq = sequence.ToUpperInvariant();
        var counts = new Dictionary<string, int>();

        for (int i = 0; i <= seq.Length - k; i++)
        {
            string kmer = seq.Substring(i, k);
            if (!counts.TryAdd(kmer, 1))
                counts[kmer]++;
        }

        return counts;
    }

    /// <summary>
    /// Counts all k-mers in a sequence with cancellation support.
    /// </summary>
    /// <param name="sequence">The sequence to analyze.</param>
    /// <param name="k">The k-mer length.</param>
    /// <param name="cancellationToken">Cancellation token for long-running operations.</param>
    /// <param name="progress">Optional progress reporter (0.0 to 1.0).</param>
    /// <returns>Dictionary mapping k-mers to their counts.</returns>
    public static Dictionary<string, int> CountKmers(
        string sequence,
        int k,
        CancellationToken cancellationToken,
        IProgress<double>? progress = null)
    {
        if (string.IsNullOrEmpty(sequence))
            return new Dictionary<string, int>();

        if (k <= 0)
            throw new ArgumentOutOfRangeException(nameof(k), "K must be positive.");

        if (k > sequence.Length)
            return new Dictionary<string, int>();

        var seq = sequence.AsSpan();
        var counts = new Dictionary<string, int>();
        int total = sequence.Length - k + 1;
        const int checkInterval = 1000;

        for (int i = 0; i <= sequence.Length - k; i++)
        {
            if (i % checkInterval == 0)
            {
                cancellationToken.ThrowIfCancellationRequested();
                progress?.Report((double)i / total);
            }

            var kmer = new string(seq.Slice(i, k));
            if (!counts.TryAdd(kmer, 1))
                counts[kmer]++;
        }

        progress?.Report(1.0);
        return counts;
    }

    /// <summary>
    /// Counts all k-mers in a sequence asynchronously.
    /// </summary>
    public static Task<Dictionary<string, int>> CountKmersAsync(
        string sequence,
        int k,
        CancellationToken cancellationToken = default,
        IProgress<double>? progress = null)
    {
        return Task.Run(() => CountKmers(sequence, k, cancellationToken, progress), cancellationToken);
    }

    /// <summary>
    /// Counts k-mers in a DNA sequence.
    /// </summary>
    public static Dictionary<string, int> CountKmers(DnaSequence dna, int k)
    {
        return CountKmers(dna.Sequence, k);
    }

    /// <summary>
    /// Counts k-mers in a DNA sequence with cancellation support.
    /// </summary>
    public static Dictionary<string, int> CountKmers(
        DnaSequence dna,
        int k,
        CancellationToken cancellationToken,
        IProgress<double>? progress = null)
    {
        return CountKmers(dna.Sequence, k, cancellationToken, progress);
    }

    /// <summary>
    /// Counts k-mers using Span-based optimization (more memory efficient).
    /// </summary>
    public static Dictionary<string, int> CountKmersSpan(ReadOnlySpan<char> sequence, int k)
    {
        return sequence.CountKmersSpan(k);
    }

    /// <summary>
    /// Gets the k-mer spectrum (frequency distribution) of a sequence.
    /// </summary>
    /// <param name="sequence">The sequence to analyze.</param>
    /// <param name="k">The k-mer length.</param>
    /// <returns>Dictionary mapping frequency to count of k-mers with that frequency.</returns>
    public static Dictionary<int, int> GetKmerSpectrum(string sequence, int k)
    {
        var counts = CountKmers(sequence, k);
        var spectrum = new Dictionary<int, int>();

        foreach (var count in counts.Values)
        {
            if (!spectrum.TryAdd(count, 1))
                spectrum[count]++;
        }

        return spectrum;
    }

    /// <summary>
    /// Finds the most frequent k-mers in a sequence.
    /// </summary>
    /// <param name="sequence">The sequence to analyze.</param>
    /// <param name="k">The k-mer length.</param>
    /// <returns>List of most frequent k-mers.</returns>
    public static IEnumerable<string> FindMostFrequentKmers(string sequence, int k)
    {
        var counts = CountKmers(sequence, k);

        if (counts.Count == 0)
            yield break;

        int maxCount = counts.Values.Max();

        foreach (var kvp in counts.Where(kvp => kvp.Value == maxCount))
        {
            yield return kvp.Key;
        }
    }

    /// <summary>
    /// Gets the k-mer frequency (normalized count).
    /// </summary>
    /// <param name="sequence">The sequence to analyze.</param>
    /// <param name="k">The k-mer length.</param>
    /// <returns>Dictionary mapping k-mers to their frequencies (0.0 to 1.0).</returns>
    public static Dictionary<string, double> GetKmerFrequencies(string sequence, int k)
    {
        var counts = CountKmers(sequence, k);
        var total = counts.Values.Sum();

        if (total == 0)
            return new Dictionary<string, double>();

        return counts.ToDictionary(
            kvp => kvp.Key,
            kvp => (double)kvp.Value / total
        );
    }

    /// <summary>
    /// Computes the k-mer distance between two sequences using Euclidean distance.
    /// </summary>
    public static double KmerDistance(string seq1, string seq2, int k)
    {
        var freq1 = GetKmerFrequencies(seq1, k);
        var freq2 = GetKmerFrequencies(seq2, k);

        var allKmers = new HashSet<string>(freq1.Keys);
        allKmers.UnionWith(freq2.Keys);

        double sumSquares = 0;
        foreach (var kmer in allKmers)
        {
            double f1 = freq1.GetValueOrDefault(kmer, 0);
            double f2 = freq2.GetValueOrDefault(kmer, 0);
            sumSquares += (f1 - f2) * (f1 - f2);
        }

        return Math.Sqrt(sumSquares);
    }

    /// <summary>
    /// Finds unique k-mers that appear only once in the sequence.
    /// </summary>
    public static IEnumerable<string> FindUniqueKmers(string sequence, int k)
    {
        var counts = CountKmers(sequence, k);
        return counts.Where(kvp => kvp.Value == 1).Select(kvp => kvp.Key);
    }

    /// <summary>
    /// Finds k-mers that appear at least a minimum number of times.
    /// </summary>
    public static IEnumerable<(string Kmer, int Count)> FindKmersWithMinCount(
        string sequence, int k, int minCount)
    {
        var counts = CountKmers(sequence, k);
        return counts
            .Where(kvp => kvp.Value >= minCount)
            .Select(kvp => (kvp.Key, kvp.Value))
            .OrderByDescending(x => x.Value);
    }

    /// <summary>
    /// Generates all possible k-mers for a given alphabet.
    /// </summary>
    /// <param name="k">The k-mer length.</param>
    /// <param name="alphabet">The alphabet to use (default: DNA = "ACGT").</param>
    public static IEnumerable<string> GenerateAllKmers(int k, string alphabet = "ACGT")
    {
        if (k <= 0)
            throw new ArgumentOutOfRangeException(nameof(k), "K must be positive.");

        if (string.IsNullOrEmpty(alphabet))
            throw new ArgumentException("Alphabet cannot be empty.", nameof(alphabet));

        return GenerateKmersRecursive("", k, alphabet);
    }

    private static IEnumerable<string> GenerateKmersRecursive(string prefix, int k, string alphabet)
    {
        if (prefix.Length == k)
        {
            yield return prefix;
            yield break;
        }

        foreach (char c in alphabet)
        {
            foreach (var kmer in GenerateKmersRecursive(prefix + c, k, alphabet))
            {
                yield return kmer;
            }
        }
    }

    /// <summary>
    /// Calculates k-mer entropy (Shannon entropy) of the sequence.
    /// Higher entropy = more diverse k-mer composition.
    /// </summary>
    public static double CalculateKmerEntropy(string sequence, int k)
    {
        var frequencies = GetKmerFrequencies(sequence, k);

        if (frequencies.Count == 0)
            return 0;

        double entropy = 0;
        foreach (var freq in frequencies.Values)
        {
            if (freq > 0)
                entropy -= freq * Math.Log2(freq);
        }

        return Math.Round(entropy, 4);
    }

    /// <summary>
    /// Finds clumps of k-mers: regions where a k-mer appears at least t times within a window of size L.
    /// </summary>
    /// <param name="sequence">The sequence to analyze.</param>
    /// <param name="k">K-mer length.</param>
    /// <param name="windowSize">Size of the sliding window (L).</param>
    /// <param name="minOccurrences">Minimum occurrences within window (t).</param>
    /// <returns>Set of k-mers that form clumps.</returns>
    public static IEnumerable<string> FindClumps(string sequence, int k, int windowSize, int minOccurrences)
    {
        if (string.IsNullOrEmpty(sequence) || k <= 0 || windowSize < k || minOccurrences <= 0)
            yield break;

        var seq = sequence.ToUpperInvariant();
        var clumps = new HashSet<string>();

        if (windowSize > seq.Length)
            yield break;

        // Initialize window
        var windowCounts = new Dictionary<string, int>();
        for (int i = 0; i < windowSize - k + 1; i++)
        {
            string kmer = seq.Substring(i, k);
            if (!windowCounts.TryAdd(kmer, 1))
                windowCounts[kmer]++;
        }

        // Check initial window
        foreach (var kvp in windowCounts)
        {
            if (kvp.Value >= minOccurrences)
                clumps.Add(kvp.Key);
        }

        // Slide window
        for (int i = 1; i <= seq.Length - windowSize; i++)
        {
            // Remove k-mer leaving window
            string leaving = seq.Substring(i - 1, k);
            windowCounts[leaving]--;
            if (windowCounts[leaving] == 0)
                windowCounts.Remove(leaving);

            // Add k-mer entering window
            string entering = seq.Substring(i + windowSize - k, k);
            if (!windowCounts.TryAdd(entering, 1))
                windowCounts[entering]++;

            // Check for clumps
            foreach (var kvp in windowCounts)
            {
                if (kvp.Value >= minOccurrences)
                    clumps.Add(kvp.Key);
            }
        }

        foreach (var kmer in clumps)
            yield return kmer;
    }

    /// <summary>
    /// Finds the positions of all occurrences of a k-mer in a sequence.
    /// </summary>
    public static IEnumerable<int> FindKmerPositions(string sequence, string kmer)
    {
        if (string.IsNullOrEmpty(sequence) || string.IsNullOrEmpty(kmer))
            yield break;

        var seq = sequence.ToUpperInvariant();
        var km = kmer.ToUpperInvariant();

        for (int i = 0; i <= seq.Length - km.Length; i++)
        {
            if (seq.Substring(i, km.Length) == km)
                yield return i;
        }
    }

    /// <summary>
    /// Counts k-mers on both strands (forward and reverse complement).
    /// </summary>
    public static Dictionary<string, int> CountKmersBothStrands(DnaSequence dna, int k)
    {
        var forwardCounts = CountKmers(dna.Sequence, k);
        var revCompCounts = CountKmers(dna.ReverseComplement().Sequence, k);

        var combined = new Dictionary<string, int>(forwardCounts);

        foreach (var kvp in revCompCounts)
        {
            if (combined.ContainsKey(kvp.Key))
                combined[kvp.Key] += kvp.Value;
            else
                combined[kvp.Key] = kvp.Value;
        }

        return combined;
    }

    /// <summary>
    /// Analyzes k-mer composition and returns statistics.
    /// </summary>
    public static KmerStatistics AnalyzeKmers(string sequence, int k)
    {
        var counts = CountKmers(sequence, k);

        if (counts.Count == 0)
        {
            return new KmerStatistics(
                TotalKmers: 0,
                UniqueKmers: 0,
                MaxCount: 0,
                MinCount: 0,
                AverageCount: 0,
                Entropy: 0
            );
        }

        var values = counts.Values.ToList();
        int totalKmers = values.Sum();
        int uniqueKmers = counts.Count;
        int maxCount = values.Max();
        int minCount = values.Min();
        double averageCount = values.Average();
        double entropy = CalculateKmerEntropy(sequence, k);

        return new KmerStatistics(
            TotalKmers: totalKmers,
            UniqueKmers: uniqueKmers,
            MaxCount: maxCount,
            MinCount: minCount,
            AverageCount: Math.Round(averageCount, 2),
            Entropy: entropy
        );
    }
}

/// <summary>
/// Statistics about k-mers in a sequence.
/// </summary>
public readonly record struct KmerStatistics(
    int TotalKmers,
    int UniqueKmers,
    int MaxCount,
    int MinCount,
    double AverageCount,
    double Entropy);
