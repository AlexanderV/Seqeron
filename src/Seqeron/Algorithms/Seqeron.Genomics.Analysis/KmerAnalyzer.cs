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

        sequence = sequence.ToUpperInvariant();
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
    /// Computes the alignment-free k-mer distance between two sequences as the Euclidean
    /// distance between their normalized k-mer frequency vectors.
    /// </summary>
    /// <remarks>
    /// Each sequence is mapped to a vector of k-mer frequencies, where a frequency is the
    /// k-mer count divided by the total number of k-mer windows (sequence length − k + 1).
    /// The distance is √(Σ (f1[w] − f2[w])²) taken over the union of k-mers occurring in
    /// either sequence; a k-mer absent from a sequence contributes a 0 component.
    /// Identical sequences yield 0. This is the frequency (relative-count) variant of the
    /// word-composition Euclidean distance: counts are normalized per Lau et al. (2022) and
    /// the Euclidean metric is applied to the relative-frequency vectors per Boden et al.
    /// (2014); the word-vector model follows Zielezinski et al. (2017) Fig. 1 and
    /// Vinga &amp; Almeida (2003).
    /// </remarks>
    /// <param name="seq1">First sequence. Null/empty or sequences shorter than <paramref name="k"/>
    /// produce an empty frequency vector (treated as the zero vector).</param>
    /// <param name="seq2">Second sequence, same conventions as <paramref name="seq1"/>.</param>
    /// <param name="k">K-mer length; must be positive.</param>
    /// <returns>Non-negative Euclidean distance between the two frequency vectors; 0 when both are empty or equal.</returns>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="k"/> is not positive.</exception>
    public static double KmerDistance(string seq1, string seq2, int k)
    {
        if (k <= 0)
            throw new ArgumentOutOfRangeException(nameof(k), "K must be positive.");

        var freq1 = GetKmerFrequencies(seq1, k);
        var freq2 = GetKmerFrequencies(seq2, k);

        // Distance spans the union of k-mers in either sequence; absent k-mers are 0.
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

    // A k-mer is "unique" when it appears exactly once in the sequence
    // (frequency = 1), as opposed to "distinct" (each different k-mer counted
    // once). See BioInfoLogics — k-mer counting, part I (2018):
    // "Unique k-mers are those that appear only once."
    private const int UniqueKmerCount = 1;

    /// <summary>
    /// Finds the unique k-mers of length <paramref name="k"/> — those that occur
    /// exactly once (overlapping occurrence count = 1) in <paramref name="sequence"/>.
    /// </summary>
    /// <param name="sequence">The sequence to analyze (case-insensitive; upper-cased internally).</param>
    /// <param name="k">The k-mer length. Must be positive.</param>
    /// <returns>
    /// The k-mers whose occurrence count equals 1. Empty when the sequence is
    /// null/empty or when <paramref name="k"/> exceeds the sequence length
    /// (L − k + 1 ≤ 0). Order is unspecified.
    /// </returns>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="k"/> ≤ 0.</exception>
    public static IEnumerable<string> FindUniqueKmers(string sequence, int k)
    {
        var counts = CountKmers(sequence, k);
        return counts.Where(kvp => kvp.Value == UniqueKmerCount).Select(kvp => kvp.Key);
    }

    /// <summary>
    /// Finds k-mers of length <paramref name="k"/> whose overlapping occurrence
    /// count is at least <paramref name="minCount"/> (recurrent k-mers,
    /// Count(Text, Pattern) ≥ t per Compeau &amp; Pevzner), ordered by count descending.
    /// </summary>
    /// <param name="sequence">The sequence to analyze (case-insensitive; upper-cased internally).</param>
    /// <param name="k">The k-mer length. Must be positive.</param>
    /// <param name="minCount">Inclusive minimum occurrence count threshold.</param>
    /// <returns>
    /// (k-mer, Count) pairs with Count ≥ <paramref name="minCount"/>, ordered by
    /// Count descending. Empty when the sequence is null/empty or k exceeds the
    /// sequence length. With <paramref name="minCount"/> ≤ 1 every distinct k-mer
    /// qualifies.
    /// </returns>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="k"/> ≤ 0.</exception>
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
    /// Generates all possible k-mers over a given alphabet: the k-fold Cartesian
    /// product of the alphabet. The number of k-mers produced is
    /// <c>alphabet.Length^k</c> (4^k for the default DNA alphabet), per the k-mer
    /// universe size n^k (Wikipedia — K-mer; Clavijo 2018, BioInfoLogics).
    /// When the alphabet is supplied in sorted order the k-mers are emitted in
    /// lexicographic order, with the rightmost position advancing fastest
    /// (odometer ordering, cf. Python itertools.product). The default alphabet
    /// "ACGT" is already sorted, so DNA k-mers are produced AAA, AAC, ..., TTT.
    /// </summary>
    /// <param name="k">The k-mer length. Must be positive.</param>
    /// <param name="alphabet">The alphabet to use (default: DNA = "ACGT"). Must be non-empty.</param>
    /// <returns>All <c>alphabet.Length^k</c> distinct k-mers.</returns>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="k"/> is not positive.</exception>
    /// <exception cref="ArgumentException">Thrown when <paramref name="alphabet"/> is null or empty.</exception>
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

        return entropy;
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
    /// Counts every k-mer over BOTH strands of double-stranded DNA: the reported
    /// count of a k-mer <c>w</c> is its overlapping occurrences on the forward
    /// strand plus its overlapping occurrences on the reverse-complement strand.
    /// </summary>
    /// <remarks>
    /// Because double-stranded DNA carries information on both strands, a strand-aware
    /// profile sums each k-mer's count with that of its complementary reading. Counting
    /// <c>w</c> on the reverse-complement strand (read 5'→3') equals counting its reverse
    /// complement <c>RC(w)</c> on the forward strand, so this method yields
    /// <c>count[w] = forward[w] + forward[RC(w)]</c>. This is the "balance" operation of
    /// kPAL — "adding the values of each k-mer to its reverse complement" (Anvar et al.,
    /// 2014) — and reflects the generalized second Chargaff rule / inversion symmetry,
    /// whereby a k-mer's count on one strand equals its reverse-complement's count on the
    /// other (Shporer et al., 2016). Unlike canonical k-mer counting (Marçais &amp;
    /// Kingsford, 2011), which collapses the pair onto a single lexicographically-smaller
    /// key, this method retains a key for every observed k-mer; w and RC(w) carry equal
    /// counts. The total over all k-mers is 2·(L − k + 1).
    /// </remarks>
    /// <param name="sequence">The DNA sequence (case-insensitive; upper-cased internally).</param>
    /// <param name="k">The k-mer length. Must be positive.</param>
    /// <returns>
    /// Dictionary mapping each observed k-mer to its summed forward + reverse-complement
    /// strand occurrence count. Empty when the sequence is null/empty or k exceeds the
    /// sequence length (L − k + 1 ≤ 0).
    /// </returns>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="k"/> ≤ 0.</exception>
    public static Dictionary<string, int> CountKmersBothStrands(string sequence, int k)
    {
        var forwardCounts = CountKmers(sequence, k);
        var revCompCounts = CountKmers(DnaSequence.GetReverseComplementString(sequence ?? string.Empty), k);

        var combined = new Dictionary<string, int>(forwardCounts);

        foreach (var kvp in revCompCounts)
        {
            if (!combined.TryAdd(kvp.Key, kvp.Value))
                combined[kvp.Key] += kvp.Value;
        }

        return combined;
    }

    /// <summary>
    /// Counts k-mers on both strands (forward and reverse complement) of a
    /// <see cref="DnaSequence"/>. Convenience overload of
    /// <see cref="CountKmersBothStrands(string, int)"/>.
    /// </summary>
    /// <param name="dna">The DNA sequence.</param>
    /// <param name="k">The k-mer length. Must be positive.</param>
    /// <returns>Dictionary mapping each observed k-mer to its summed forward + reverse-complement count.</returns>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="k"/> ≤ 0.</exception>
    public static Dictionary<string, int> CountKmersBothStrands(DnaSequence dna, int k)
        => CountKmersBothStrands(dna.Sequence, k);

    /// <summary>
    /// Computes comprehensive k-mer composition statistics for a sequence: the
    /// total number of (overlapping) k-mers, the number of distinct k-mers, the
    /// maximum/minimum/average multiplicity, and the Shannon entropy of the k-mer
    /// frequency distribution.
    /// </summary>
    /// <remarks>
    /// The total k-mer count is the number of overlapping length-k windows,
    /// L − k + 1 (Wikipedia — K-mer; BioInfoLogics, k-mer counting part I, 2018).
    /// <see cref="KmerStatistics.UniqueKmers"/> reports the number of <i>distinct</i>
    /// k-mers (each different k-mer counted once) — not the count-1 "unique" set of
    /// <see cref="FindUniqueKmers"/>. <see cref="KmerStatistics.AverageCount"/> is the
    /// mean multiplicity total/distinct, rounded to two decimals for display.
    /// <see cref="KmerStatistics.Entropy"/> is the k-mer Shannon entropy
    /// E_k = −Σ p(α) log₂ p(α) with p(α) = mult(α) / (L − k + 1), the relative
    /// frequency of k-mer α over the L − k + 1 windows (Manca et al., 2021,
    /// "Spectral concepts in genome informational analysis", arXiv:2106.15351;
    /// log base 2 ⇒ bits). An empty sequence or k exceeding the length yields an
    /// all-zero result (L − k + 1 ≤ 0 ⇒ no k-mers).
    /// </remarks>
    /// <param name="sequence">The sequence to analyze (case-insensitive; upper-cased internally).</param>
    /// <param name="k">The k-mer length. Must be positive.</param>
    /// <returns>A <see cref="KmerStatistics"/> record with the composition statistics; all-zero when no k-mers exist.</returns>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="k"/> ≤ 0.</exception>
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
/// Comprehensive k-mer composition statistics for a sequence.
/// </summary>
/// <param name="TotalKmers">Total number of overlapping k-mers, L − k + 1.</param>
/// <param name="UniqueKmers">Number of <i>distinct</i> k-mers (each different k-mer counted once).</param>
/// <param name="MaxCount">Maximum k-mer multiplicity observed.</param>
/// <param name="MinCount">Minimum k-mer multiplicity observed.</param>
/// <param name="AverageCount">Mean multiplicity (TotalKmers / UniqueKmers), rounded to two decimals.</param>
/// <param name="Entropy">Shannon entropy of the k-mer frequency distribution, −Σ p log₂ p, in bits.</param>
public readonly record struct KmerStatistics(
    int TotalKmers,
    int UniqueKmers,
    int MaxCount,
    int MinCount,
    double AverageCount,
    double Entropy);
