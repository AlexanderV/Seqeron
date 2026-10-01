namespace Seqeron.Genomics.Analysis;

/// <summary>
/// Analyzes k-mers (substrings of length k) in DNA/RNA sequences.
/// </summary>
public static class KmerAnalyzer
{
    /// <summary>
    /// Counts all overlapping k-mers in a sequence (sliding window of length k, step 1).
    /// </summary>
    /// <remarks>
    /// Count(w) = number of start positions i in [0, L − k] with sequence[i..i+k−1] = w, so the
    /// counts sum to L − k + 1 (Wikipedia — K-mer; Rosalind KMER k-mer composition). Input is
    /// upper-cased (case-insensitive); every symbol, including IUPAC ambiguity codes such as N,
    /// is counted literally (no alphabet filtering — unlike Jellyfish, which drops windows
    /// containing a non-ACGT base). Only observed k-mers are keys (zero counts are absent).
    /// Single canonical loop: delegates to the cancellation-aware overload.
    /// </remarks>
    /// <param name="sequence">The sequence to analyze. Null/empty returns an empty dictionary.</param>
    /// <param name="k">The k-mer length. Must be positive for non-empty input.</param>
    /// <returns>Dictionary mapping k-mers to their counts; empty when k exceeds the length.</returns>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="k"/> ≤ 0 and the sequence is non-empty.</exception>
    public static Dictionary<string, int> CountKmers(string sequence, int k)
        => CountKmers(sequence, k, CancellationToken.None);

    /// <summary>
    /// Counts all k-mers in a sequence with cancellation support.
    /// </summary>
    /// <remarks>
    /// Cooperative cancellation (Microsoft Learn — "Cancellation in managed threads"): the token is
    /// polled with <see cref="CancellationToken.ThrowIfCancellationRequested"/> every 1000 windows,
    /// so the thrown <see cref="OperationCanceledException"/> carries this token. Progress is
    /// reported synchronously to <paramref name="progress"/> (TAP guidance) as the fraction of
    /// windows scanned, i/(L − k + 1) at each checkpoint (strictly increasing, in [0, 1)), and
    /// exactly one final 1.0 on every successful completion — including the trivial empty results.
    /// </remarks>
    /// <param name="sequence">The sequence to analyze.</param>
    /// <param name="k">The k-mer length.</param>
    /// <param name="cancellationToken">Cancellation token for long-running operations.</param>
    /// <param name="progress">Optional progress reporter (0.0 to 1.0); null reports nothing.</param>
    /// <returns>Dictionary mapping k-mers to their counts.</returns>
    public static Dictionary<string, int> CountKmers(
        string sequence,
        int k,
        CancellationToken cancellationToken,
        IProgress<double>? progress = null)
    {
        ValidateKmerLength(sequence, k);

        if (string.IsNullOrEmpty(sequence) || k > sequence.Length)
        {
            progress?.Report(1.0);
            return new Dictionary<string, int>();
        }

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
    /// Counts all k-mers in a sequence asynchronously (thread-pool offload of
    /// <see cref="CountKmers(string, int, CancellationToken, IProgress{double}?)"/>).
    /// </summary>
    /// <remarks>
    /// Task-based Asynchronous Pattern (Microsoft Learn — "Task-based asynchronous pattern (TAP)"
    /// and "Implementing the TAP"): the usage error k ≤ 0 (non-empty input) is validated
    /// synchronously and thrown directly from this call; every other outcome is carried by the
    /// returned task. A token already signaled at call time yields a <see cref="TaskStatus.Canceled"/>
    /// task (Task.Run does not start the delegate); a token signaled during the scan is observed at
    /// the next checkpoint and the task ends Canceled with an <see cref="OperationCanceledException"/>
    /// carrying the same token. Note: TAP advises exposing purely compute-bound work only
    /// synchronously; this wrapper is kept for API compatibility and adds no parallelism.
    /// </remarks>
    /// <exception cref="ArgumentOutOfRangeException">Thrown synchronously when <paramref name="k"/> ≤ 0 and the sequence is non-empty.</exception>
    public static Task<Dictionary<string, int>> CountKmersAsync(
        string sequence,
        int k,
        CancellationToken cancellationToken = default,
        IProgress<double>? progress = null)
    {
        ValidateKmerLength(sequence, k);
        return Task.Run(() => CountKmers(sequence, k, cancellationToken, progress), cancellationToken);
    }

    /// <summary>k must be positive for non-empty input (null/empty input yields an empty count for any k).</summary>
    private static void ValidateKmerLength(string sequence, int k)
    {
        if (!string.IsNullOrEmpty(sequence) && k <= 0)
            throw new ArgumentOutOfRangeException(nameof(k), "K must be positive.");
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
    /// <remarks>
    /// The spectrum is the count-of-counts histogram: for each multiplicity m it gives the number
    /// of distinct k-mers occurring exactly m times (Chor et al., 2009, Genome Biol. 10:R108;
    /// same semantics as Jellyfish <c>histo</c>, which increments <c>histo[count]</c> once per
    /// distinct k-mer). Only non-zero bins are returned (Jellyfish default, no <c>--full</c>), and
    /// there is no upper cap bin. Σ m·spectrum[m] = L − k + 1. Counts come from
    /// <see cref="CountKmers(string, int)"/> (single strand, case-insensitive, no alphabet filtering).
    /// </remarks>
    /// <param name="sequence">The sequence to analyze. Null/empty returns an empty dictionary.</param>
    /// <param name="k">The k-mer length.</param>
    /// <returns>Dictionary mapping multiplicity to the number of distinct k-mers with that multiplicity.</returns>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="k"/> ≤ 0 and the sequence is non-empty.</exception>
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
    /// Finds all most frequent k-mers in a sequence (Frequent Words Problem).
    /// </summary>
    /// <remarks>
    /// Returns every k-mer Pattern maximizing Count(Text, Pattern) over all k-mers of Text, where
    /// Count uses overlapping occurrences (Rosalind BA1B; Compeau &amp; Pevzner, <i>Bioinformatics
    /// Algorithms</i>, ch. 1 — the "better" frequent-words algorithm built on a single frequency map,
    /// here the canonical <see cref="CountKmers(string,int)"/>). O(|Text|·k) time.
    /// Example (BA1B sample): ACGTTGCATGTCGCATGATGCATGAGAGCT, k=4 → {CATG, GCAT} (3 occurrences each).
    /// </remarks>
    /// <param name="sequence">The sequence to analyze (case-insensitive; upper-cased internally).</param>
    /// <param name="k">The k-mer length. Must be positive for non-empty input.</param>
    /// <returns>
    /// All k-mers tied at the maximum count, in order of first occurrence (order is not part of the
    /// contract). Empty when the sequence is null/empty or k exceeds its length.
    /// </returns>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="k"/> ≤ 0 and the sequence is non-empty.</exception>
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
    /// <remarks>
    /// f(w) = count(w) / (L − k + 1): the denominator is the number of overlapping windows, which
    /// equals the sum of all counts (scikit-bio <c>Sequence.kmer_frequencies(k, overlap=True,
    /// relative=True)</c>). Only observed k-mers are keys; values sum to 1.
    /// </remarks>
    /// <param name="sequence">The sequence to analyze. Null/empty returns an empty dictionary.</param>
    /// <param name="k">The k-mer length.</param>
    /// <returns>Dictionary mapping k-mers to their frequencies (0.0 to 1.0).</returns>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="k"/> ≤ 0 and the sequence is non-empty.</exception>
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
    /// <remarks>
    /// "Unique" follows Jellyfish <c>stats</c> ("Unique" = number of k-mers with count 1,
    /// <c>uniq += val == 1</c>) and BioInfoLogics (2018); it equals <c>jellyfish dump -L 1 -U 1</c>
    /// (KMC <c>-ci1 -cx1</c>) and is the set counted by <see cref="KmerStatistics.SingletonKmers"/>.
    /// It is <b>not</b> the number of distinct k-mers (<see cref="KmerStatistics.DistinctKmers"/>).
    /// Implemented as <see cref="FindKmersWithMinCount(string, int, int, int)"/> with the range [1, 1]
    /// over one canonical <see cref="CountKmers(string, int)"/> table.
    /// </remarks>
    /// <param name="sequence">The sequence to analyze (case-insensitive; upper-cased internally).</param>
    /// <param name="k">The k-mer length. Must be positive.</param>
    /// <returns>
    /// The k-mers whose occurrence count equals 1, in ascending ordinal (lexicographic) order.
    /// Empty when the sequence is null/empty or when <paramref name="k"/> exceeds the sequence
    /// length (L − k + 1 ≤ 0).
    /// </returns>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="k"/> ≤ 0 and the sequence is non-empty.</exception>
    public static IEnumerable<string> FindUniqueKmers(string sequence, int k)
        => FindKmersWithMinCount(sequence, k, UniqueKmerCount, UniqueKmerCount).Select(p => p.Kmer);

    /// <summary>
    /// Finds k-mers of length <paramref name="k"/> whose overlapping occurrence
    /// count is at least <paramref name="minCount"/> (recurrent k-mers,
    /// Count(Text, Pattern) ≥ t per Compeau &amp; Pevzner), ordered by count descending.
    /// </summary>
    /// <remarks>
    /// Equivalent to <see cref="FindKmersWithMinCount(string, int, int, int)"/> with no upper bound
    /// (<c>jellyfish dump -L minCount</c>).
    /// </remarks>
    /// <param name="sequence">The sequence to analyze (case-insensitive; upper-cased internally).</param>
    /// <param name="k">The k-mer length. Must be positive.</param>
    /// <param name="minCount">Inclusive minimum occurrence count threshold.</param>
    /// <returns>
    /// (k-mer, Count) pairs with Count ≥ <paramref name="minCount"/>, ordered by
    /// Count descending, ties in ascending ordinal k-mer order. Empty when the sequence is
    /// null/empty or k exceeds the sequence length. With <paramref name="minCount"/> ≤ 1
    /// (including zero or negative values) every distinct k-mer qualifies.
    /// </returns>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="k"/> ≤ 0 and the sequence is non-empty.</exception>
    public static IEnumerable<(string Kmer, int Count)> FindKmersWithMinCount(
        string sequence, int k, int minCount)
        => FindKmersWithMinCount(sequence, k, minCount, int.MaxValue);

    /// <summary>
    /// Finds k-mers of length <paramref name="k"/> whose overlapping occurrence count lies in the
    /// inclusive range [<paramref name="minCount"/>, <paramref name="maxCount"/>] — the
    /// <c>-L/--lower-count</c> and <c>-U/--upper-count</c> filters of <c>jellyfish dump</c>
    /// (KMC <c>-ci</c>/<c>-cx</c>).
    /// </summary>
    /// <remarks>
    /// Mirrors Jellyfish <c>dump</c> (<c>sub_commands/dump_main.cc</c>:
    /// <c>if(it.val() &lt; lower_count || it.val() &gt; upper_count) continue;</c>), applied to the one
    /// canonical <see cref="CountKmers(string, int)"/> table through the same range filter that
    /// <see cref="AnalyzeKmers(string, int, int, int)"/> uses. Jellyfish emits its hash order; this
    /// method fixes a deterministic order instead: count descending, ties by ascending ordinal k-mer
    /// (the order of <c>jellyfish dump -c | sort -k2,2nr -k1,1</c>). <paramref name="maxCount"/> &lt;
    /// <paramref name="minCount"/> selects nothing, as in Jellyfish. O(L·k + d log d) time, O(d·k) space.
    /// </remarks>
    /// <param name="sequence">The sequence to analyze (case-insensitive; upper-cased internally).</param>
    /// <param name="k">The k-mer length. Must be positive.</param>
    /// <param name="minCount">Inclusive lower count bound (≤ 1 imposes no lower bound).</param>
    /// <param name="maxCount">Inclusive upper count bound; <see cref="int.MaxValue"/> = unbounded.</param>
    /// <returns>(k-mer, Count) pairs in the range, count descending then ordinal k-mer ascending.</returns>
    /// <exception cref="ArgumentOutOfRangeException">
    /// <paramref name="maxCount"/> is negative (Jellyfish takes unsigned values), or <paramref name="k"/> ≤ 0
    /// and the sequence is non-empty.
    /// </exception>
    public static IEnumerable<(string Kmer, int Count)> FindKmersWithMinCount(
        string sequence, int k, int minCount, int maxCount)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(maxCount);

        var counts = CountKmers(sequence, k);
        return SelectByCountRange(counts, minCount, maxCount)
            .OrderByDescending(kvp => kvp.Value)
            .ThenBy(kvp => kvp.Key, StringComparer.Ordinal)
            .Select(kvp => (kvp.Key, kvp.Value));
    }

    /// <summary>
    /// Canonical count-range filter shared by <see cref="FindKmersWithMinCount(string, int, int, int)"/>,
    /// <see cref="FindUniqueKmers"/> and <see cref="AnalyzeKmers(string, int, int, int)"/>: keeps the entries
    /// whose count is in [<paramref name="lowerCount"/>, <paramref name="upperCount"/>], i.e. skips
    /// <c>count &lt; lower || count &gt; upper</c> exactly as Jellyfish <c>dump</c>/<c>stats</c> do.
    /// </summary>
    private static IEnumerable<KeyValuePair<string, int>> SelectByCountRange(
        Dictionary<string, int> counts, int lowerCount, int upperCount)
    {
        foreach (var kvp in counts)
        {
            if (kvp.Value < lowerCount || kvp.Value > upperCount)
                continue;
            yield return kvp;
        }
    }

    /// <summary>
    /// Generates all possible k-mers over a given alphabet: the k-fold Cartesian
    /// product of the alphabet. The number of k-mers produced is
    /// <c>alphabet.Length^k</c> (4^k for the default DNA alphabet), per the k-mer
    /// universe size n^k (Wikipedia — K-mer; Clavijo 2018, BioInfoLogics).
    /// The k-mers are emitted in lexicographic order <em>with respect to the order in
    /// which the symbols appear in <paramref name="alphabet"/></em> (Rosalind LEXF,
    /// "Enumerating k-mers Lexicographically": the alphabet is an ordered permutation
    /// a1 &lt; a2 &lt; … and s &lt;Lex t iff the first mismatching symbol of s precedes
    /// that of t in the alphabet), with the rightmost position advancing fastest
    /// (odometer ordering, identical to Python <c>itertools.product(alphabet, repeat=k)</c>).
    /// The default alphabet "ACGT" is already sorted, so DNA k-mers are produced
    /// AAA, AAC, ..., TTT; e.g. alphabet "TAGC", k=2 gives TT, TA, TG, TC, AT, ..., CC.
    /// </summary>
    /// <remarks>
    /// Enumeration is lazy and uses an index odometer over a single k-character buffer,
    /// so the working space is O(k) and each k-mer costs O(k) (one string allocation).
    /// The alphabet is used verbatim (case-sensitive, not sorted, not de-duplicated): as with
    /// <c>itertools.product</c>, a repeated symbol yields repeated k-mers, so the output is
    /// duplicate-free exactly when the alphabet symbols are distinct.
    /// </remarks>
    /// <param name="k">The k-mer length. Must be positive.</param>
    /// <param name="alphabet">The ordered alphabet to use (default: DNA = "ACGT"). Must be non-empty.</param>
    /// <returns>All <c>alphabet.Length^k</c> k-mers (distinct when the alphabet symbols are distinct).</returns>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="k"/> is not positive.</exception>
    /// <exception cref="ArgumentException">Thrown when <paramref name="alphabet"/> is null or empty.</exception>
    public static IEnumerable<string> GenerateAllKmers(int k, string alphabet = "ACGT")
    {
        if (k <= 0)
            throw new ArgumentOutOfRangeException(nameof(k), "K must be positive.");

        if (string.IsNullOrEmpty(alphabet))
            throw new ArgumentException("Alphabet cannot be empty.", nameof(alphabet));

        return EnumerateCartesianProduct(k, alphabet);
    }

    /// <summary>
    /// Odometer enumeration of alphabet^k (the algorithm of CPython itertools.product):
    /// indices[i] is the alphabet position at k-mer position i; the rightmost index is
    /// advanced on every step and carries leftwards on wrap-around.
    /// </summary>
    private static IEnumerable<string> EnumerateCartesianProduct(int k, string alphabet)
    {
        int n = alphabet.Length;
        var indices = new int[k];
        var buffer = new char[k];
        Array.Fill(buffer, alphabet[0]);

        while (true)
        {
            yield return new string(buffer);

            int pos = k - 1;
            while (pos >= 0 && indices[pos] == n - 1)
            {
                indices[pos] = 0;
                buffer[pos] = alphabet[0];
                pos--;
            }

            if (pos < 0)
                yield break;

            indices[pos]++;
            buffer[pos] = alphabet[indices[pos]];
        }
    }

    /// <summary>
    /// Calculates k-mer entropy (Shannon entropy) of the sequence.
    /// Higher entropy = more diverse k-mer composition.
    /// </summary>
    /// <remarks>
    /// H = −Σ p(w) log₂ p(w) with p(w) = count(w) / (L − k + 1) (Shannon, 1948), in bits;
    /// 0 ≤ H ≤ log₂(number of distinct k-mers). Not normalized. Returns 0 when no k-mer exists
    /// (null/empty sequence or k &gt; L). Delegates to the canonical plug-in k-mer entropy
    /// <see cref="SequenceComplexity.CalculateKmerEntropy(string, int)"/> (canonical
    /// <see cref="StatisticsHelper.ShannonIndex"/> in nats ÷ ln 2 — the computation of
    /// <c>scipy.stats.entropy(counts, base=2)</c>), keeping this method's own contract: null/empty
    /// input returns 0 for any k, k ≤ 0 throws only for non-empty input.
    /// </remarks>
    /// <param name="sequence">The sequence to analyze (case-insensitive).</param>
    /// <param name="k">The k-mer length.</param>
    /// <returns>Shannon entropy in bits of the overlapping k-mer frequency distribution.</returns>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="k"/> ≤ 0 and the sequence is non-empty.</exception>
    public static double CalculateKmerEntropy(string sequence, int k)
    {
        if (string.IsNullOrEmpty(sequence))
            return 0;

        ValidateKmerLength(sequence, k);
        return SequenceComplexity.CalculateKmerEntropy(sequence, k);
    }

    /// <summary>
    /// Shannon entropy in bits of a multiplicity table: the canonical
    /// <see cref="StatisticsHelper.ShannonIndex"/> (nats) divided by ln 2 — the same operation, in the
    /// same order, as <see cref="SequenceComplexity.CalculateKmerEntropy(string, int)"/>, so a table taken
    /// from <see cref="CountKmers(string, int)"/> gives a bit-identical value. Requires a positive total.
    /// </summary>
    private static double ShannonEntropyBits(IReadOnlyList<int> counts)
        => StatisticsHelper.ShannonIndex(counts) / Ln2;

    private static readonly double Ln2 = Math.Log(2.0);

    /// <summary>
    /// Finds all distinct k-mers forming (L, t)-clumps: a k-mer forms an (L, t)-clump if some
    /// window (substring) of length L of the sequence contains at least t occurrences of it.
    /// </summary>
    /// <remarks>
    /// Clump Finding Problem (Rosalind BA1E; Compeau &amp; Pevzner, <i>Bioinformatics Algorithms</i>,
    /// ch. 1): windows are the substrings <c>Genome[i..i+L−1]</c> for i in [0, |Genome| − L], and an
    /// occurrence counts only if it lies entirely inside the window (start p with i ≤ p ≤ i + L − k),
    /// so each window holds L − k + 1 k-mer starts; overlapping occurrences are counted.
    /// Implements the textbook's <c>BetterClumpFinding</c>: the first window is counted with the
    /// canonical <see cref="CountKmers(string,int)"/>, then each slide decrements the k-mer leaving
    /// the window and increments the entering one. Only the entering k-mer's count can grow, so only
    /// it needs to be tested against t — O(|Genome|·k) time instead of rescanning the whole window.
    /// Matching is case-insensitive (upper-cased), mirroring <see cref="CountKmers(string,int)"/>.
    /// Cross-check: E. coli genome (textbook dataset), k=9, L=500, t=3 → 1904 distinct 9-mers.
    /// </remarks>
    /// <param name="sequence">The sequence to analyze (Genome).</param>
    /// <param name="k">K-mer length.</param>
    /// <param name="windowSize">Window length L.</param>
    /// <param name="minOccurrences">Minimum occurrences t within one window (inclusive).</param>
    /// <returns>
    /// Each clump-forming k-mer exactly once, in order of first detection (order is not part of
    /// the contract). Empty when the sequence is null/empty, k ≤ 0, L &lt; k, L &gt; |sequence|
    /// (no window of length L exists) or t ≤ 0.
    /// </returns>
    public static IEnumerable<string> FindClumps(string sequence, int k, int windowSize, int minOccurrences)
    {
        if (string.IsNullOrEmpty(sequence) || k <= 0 || windowSize < k || minOccurrences <= 0
            || windowSize > sequence.Length)
            yield break;

        var seq = sequence.ToUpperInvariant();
        var clumps = new HashSet<string>();

        // First window Genome[0..L−1]: canonical k-mer counting.
        var windowCounts = CountKmers(seq.Substring(0, windowSize), k);
        foreach (var kvp in windowCounts)
        {
            if (kvp.Value >= minOccurrences && clumps.Add(kvp.Key))
                yield return kvp.Key;
        }

        // Slide: window i covers k-mer starts i..i+L−k.
        for (int i = 1; i <= seq.Length - windowSize; i++)
        {
            string leaving = seq.Substring(i - 1, k);
            if (--windowCounts[leaving] == 0)
                windowCounts.Remove(leaving);

            string entering = seq.Substring(i + windowSize - k, k);
            int count = windowCounts.TryGetValue(entering, out int c) ? c + 1 : 1;
            windowCounts[entering] = count;

            if (count >= minOccurrences && clumps.Add(entering))
                yield return entering;
        }
    }

    /// <summary>
    /// Finds all starting positions where a k-mer occurs in a sequence.
    /// </summary>
    /// <remarks>
    /// Solves the Pattern Matching Problem: "find all occurrences of a pattern in a string"
    /// and report "all starting positions in Genome where Pattern appears as a substring"
    /// using 0-based indexing (Rosalind BA1D; Compeau &amp; Pevzner, <i>Bioinformatics
    /// Algorithms</i>). Occurrences may overlap and every overlapping start is reported,
    /// e.g. <c>AA</c> in <c>AAAA</c> yields 0, 1, 2. There are at most
    /// <c>L − k + 1</c> candidate positions for a length-L sequence (Wikipedia, "k-mer").
    /// Matching is case-insensitive (both arguments are upper-cased), mirroring the sibling
    /// <see cref="CountKmers(string,int)"/> methods.
    /// <para>
    /// The scan is the Knuth–Morris–Pratt matcher (Knuth, Morris &amp; Pratt 1977,
    /// <i>SIAM J. Comput.</i> 6:323–350; CLRS §32.4 KMP-MATCHER / COMPUTE-PREFIX-FUNCTION):
    /// the prefix function of the k-mer is built in O(k), then the text is read once, left
    /// to right, in O(L), so the worst case is O(L + k) instead of the Θ((L − k + 1)·k) of a
    /// window-by-window comparison (e.g. <c>A</c><sup>10⁶</sup> searched for
    /// <c>A</c><sup>99 999</sup><c>C</c>). After a full match the state falls back to
    /// π[k − 1], so overlapping occurrences are all reported, in ascending order.
    /// The repository <c>SuffixTree.FindAllOccurrences</c> (used by
    /// <c>MotifFinder.FindExactMotif</c> on <c>DnaSequence</c>) is an index for many queries
    /// against one text; for a single query it costs an O(L) tree build with a large constant
    /// and memory, and returns positions unordered (see algorithm doc §5.2).
    /// </para>
    /// </remarks>
    /// <param name="sequence">The sequence to search (case-insensitive).</param>
    /// <param name="kmer">The k-mer / pattern to locate (case-insensitive).</param>
    /// <returns>
    /// Ascending 0-based start positions of every (possibly overlapping) occurrence.
    /// Empty when <paramref name="sequence"/> or <paramref name="kmer"/> is null/empty,
    /// when the k-mer is longer than the sequence, or when the k-mer does not occur.
    /// </returns>
    public static IEnumerable<int> FindKmerPositions(string sequence, string kmer)
    {
        if (string.IsNullOrEmpty(sequence) || string.IsNullOrEmpty(kmer))
            yield break;

        if (kmer.Length > sequence.Length)
            yield break;

        var text = sequence.ToUpperInvariant();
        var pattern = kmer.ToUpperInvariant();
        int m = pattern.Length;
        int[] prefix = ComputeKmpPrefixFunction(pattern);

        // KMP-MATCHER (CLRS §32.4): q = number of pattern characters currently matched.
        int q = 0;
        for (int i = 0; i < text.Length; i++)
        {
            char c = text[i];
            while (q > 0 && pattern[q] != c)
                q = prefix[q - 1];
            if (pattern[q] == c)
                q++;
            if (q == m)
            {
                yield return i - m + 1;
                // Fall back to the longest proper border so overlapping occurrences
                // are reported too (Rosalind BA1D / SUBS).
                q = prefix[q - 1];
            }
        }
    }

    /// <summary>
    /// KMP prefix (failure) function: <c>prefix[q]</c> is the length of the longest proper
    /// prefix of <c>pattern[0..q]</c> that is also its suffix (CLRS §32.4
    /// COMPUTE-PREFIX-FUNCTION, 0-based). O(|pattern|).
    /// </summary>
    private static int[] ComputeKmpPrefixFunction(string pattern)
    {
        var prefix = new int[pattern.Length];
        int k = 0;
        for (int q = 1; q < pattern.Length; q++)
        {
            while (k > 0 && pattern[k] != pattern[q])
                k = prefix[k - 1];
            if (pattern[k] == pattern[q])
                k++;
            prefix[q] = k;
        }
        return prefix;
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
    /// number of singleton (count-1) k-mers, the maximum/minimum/mean multiplicity,
    /// and the Shannon entropy of the k-mer frequency distribution.
    /// </summary>
    /// <remarks>
    /// Equivalent to <see cref="AnalyzeKmers(string, int, int, int)"/> with no count filter
    /// (Jellyfish <c>stats</c> defaults: lower-count 0, upper-count unbounded). The total k-mer count is
    /// the number of overlapping length-k windows, L − k + 1 (Wikipedia — K-mer; BioInfoLogics,
    /// k-mer counting part I, 2018). Field definitions follow Jellyfish <c>stats</c>
    /// (<c>sub_commands/stats_main.cc</c>, <c>compute_stats</c>): Total = Σ count, Distinct = number of
    /// k-mers, Unique = number of k-mers with count 1, Max_count = max count. Note the naming:
    /// <see cref="KmerStatistics.UniqueKmers"/> is kept for backward compatibility and reports the
    /// number of <i>distinct</i> k-mers (Jellyfish "Distinct", also exposed as
    /// <see cref="KmerStatistics.DistinctKmers"/>); Jellyfish's "Unique" (count == 1, the set returned
    /// by <see cref="FindUniqueKmers"/>) is <see cref="KmerStatistics.SingletonKmers"/>.
    /// <see cref="KmerStatistics.AverageCount"/> is the exact mean multiplicity total/distinct (not rounded).
    /// <see cref="KmerStatistics.Entropy"/> is the k-mer Shannon entropy
    /// E_k = −Σ p(α) log₂ p(α) with p(α) = mult(α) / (L − k + 1), the relative
    /// frequency of k-mer α over the L − k + 1 windows (Manca et al., 2021,
    /// "Spectral concepts in genome informational analysis", arXiv:2106.15351;
    /// log base 2 ⇒ bits), bit-identical to <see cref="CalculateKmerEntropy"/>. An empty sequence or
    /// k exceeding the length yields an all-zero result (L − k + 1 ≤ 0 ⇒ no k-mers).
    /// The count table is built once with <see cref="CountKmers(string, int)"/> and every statistic is
    /// derived from it in one pass: O(L·k) time, O(D·k) space.
    /// </remarks>
    /// <param name="sequence">The sequence to analyze (case-insensitive; upper-cased internally).</param>
    /// <param name="k">The k-mer length. Must be positive.</param>
    /// <returns>A <see cref="KmerStatistics"/> record with the composition statistics; all-zero when no k-mers exist.</returns>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="k"/> ≤ 0 and the sequence is non-empty.</exception>
    public static KmerStatistics AnalyzeKmers(string sequence, int k)
        => AnalyzeKmers(sequence, k, 0, int.MaxValue);

    /// <summary>
    /// Computes k-mer composition statistics over the k-mers whose count lies in
    /// [<paramref name="lowerCount"/>, <paramref name="upperCount"/>] — the <c>-L/--lower-count</c> and
    /// <c>-U/--upper-count</c> filters of Jellyfish <c>stats</c>.
    /// </summary>
    /// <remarks>
    /// Mirrors Jellyfish <c>compute_stats</c> (<c>sub_commands/stats_main.cc</c>): a k-mer with
    /// count &lt; lower or count &gt; upper is skipped; for every retained k-mer, Unique (here
    /// <see cref="KmerStatistics.SingletonKmers"/>) += (count == 1), Total += count, Max_count = max,
    /// Distinct += 1. The remaining fields apply the same definitions to the retained multiset:
    /// MinCount = min retained count, AverageCount = Total / Distinct, Entropy = −Σ p log₂ p with
    /// p = count / Total (Total is the retained sum, which is L − k + 1 only when nothing is filtered).
    /// No retained k-mer (including upper &lt; lower) ⇒ all-zero result, as Jellyfish prints zeros.
    /// </remarks>
    /// <param name="sequence">The sequence to analyze (case-insensitive; upper-cased internally).</param>
    /// <param name="k">The k-mer length. Must be positive.</param>
    /// <param name="lowerCount">Ignore k-mers with count below this value (Jellyfish default 0).</param>
    /// <param name="upperCount">Ignore k-mers with count above this value (Jellyfish default unbounded).</param>
    /// <returns>Statistics over the retained k-mers; all-zero when none is retained.</returns>
    /// <exception cref="ArgumentOutOfRangeException">
    /// <paramref name="lowerCount"/> or <paramref name="upperCount"/> is negative (Jellyfish takes unsigned
    /// values), or <paramref name="k"/> ≤ 0 and the sequence is non-empty.
    /// </exception>
    public static KmerStatistics AnalyzeKmers(string sequence, int k, int lowerCount, int upperCount = int.MaxValue)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(lowerCount);
        ArgumentOutOfRangeException.ThrowIfNegative(upperCount);

        var counts = CountKmers(sequence, k);
        var retained = new List<int>(counts.Count);
        int total = 0, singletons = 0, maxCount = 0, minCount = int.MaxValue;

        foreach (var (_, count) in SelectByCountRange(counts, lowerCount, upperCount))
        {
            retained.Add(count);
            total += count;
            if (count == UniqueKmerCount) singletons++;
            if (count > maxCount) maxCount = count;
            if (count < minCount) minCount = count;
        }

        if (retained.Count == 0)
            return new KmerStatistics(0, 0, 0, 0, 0, 0);

        return new KmerStatistics(
            TotalKmers: total,
            UniqueKmers: retained.Count,
            MaxCount: maxCount,
            MinCount: minCount,
            AverageCount: (double)total / retained.Count,
            Entropy: ShannonEntropyBits(retained))
        {
            SingletonKmers = singletons,
        };
    }
}

/// <summary>
/// Comprehensive k-mer composition statistics for a sequence (field definitions per Jellyfish <c>stats</c>:
/// Total, Distinct, Unique = count-1, Max_count; plus min/mean multiplicity and k-mer Shannon entropy).
/// </summary>
/// <param name="TotalKmers">Total number of k-mers including multiplicity (Jellyfish "Total"), L − k + 1 when unfiltered.</param>
/// <param name="UniqueKmers">
/// Number of <i>distinct</i> k-mers (Jellyfish "Distinct"; each different k-mer counted once). The name is kept
/// for backward compatibility; it is <b>not</b> Jellyfish's "Unique" — see <see cref="SingletonKmers"/>.
/// </param>
/// <param name="MaxCount">Maximum k-mer multiplicity observed (Jellyfish "Max_count").</param>
/// <param name="MinCount">Minimum k-mer multiplicity observed.</param>
/// <param name="AverageCount">Exact mean multiplicity, TotalKmers / UniqueKmers.</param>
/// <param name="Entropy">Shannon entropy of the k-mer frequency distribution, −Σ p log₂ p, in bits.</param>
public readonly record struct KmerStatistics(
    int TotalKmers,
    int UniqueKmers,
    int MaxCount,
    int MinCount,
    double AverageCount,
    double Entropy)
{
    /// <summary>Number of distinct k-mers (Jellyfish "Distinct"); same value as <see cref="UniqueKmers"/>.</summary>
    public int DistinctKmers => UniqueKmers;

    /// <summary>
    /// Number of k-mers occurring exactly once — Jellyfish <c>stats</c> "Unique" (<c>uniq += val == 1</c>);
    /// equals the size of <see cref="KmerAnalyzer.FindUniqueKmers"/>.
    /// </summary>
    public int SingletonKmers { get; init; }
}
