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
        => CountKmersCore(sequence, k, acgtOnly: false, cancellationToken, progress);

    /// <summary>
    /// Counts k-mers under explicit <see cref="KmerCountingOptions"/>: literal (default), ACGT-only
    /// (Jellyfish window convention) or canonical (Jellyfish <c>count -C</c>).
    /// </summary>
    /// <remarks>
    /// <para><b>ACGT-only</b> (<see cref="KmerCountingOptions.AcgtOnly"/>): after case folding, a window is
    /// counted only when all k symbols are A, C, G or T. This is Jellyfish's convention: in
    /// <c>include/jellyfish/mer_iterator.hpp</c> a base whose 2-bit code (<c>mer_dna.hpp</c> <c>codes[256]</c>:
    /// A/a=0, C/c=1, G/g=2, T/t=3, every other byte negative, IUPAC codes and U included) is negative
    /// resets <c>filled_</c> to 0, so no k-mer overlapping it is emitted. Total = number of
    /// all-ACGT windows (≤ L − k + 1).</para>
    /// <para><b>Canonical</b> (<see cref="KmerCountingOptions.Canonical"/>): each k-mer w is keyed by
    /// min(w, RC(w)) — Jellyfish <c>-C</c> "count both strand, canonical representation"
    /// (<c>mer_iterator</c> returns <c>m_ &lt; rcm_ ? m_ : rcm_</c>; <c>mer_dna::get_canonical</c>;
    /// Marçais &amp; Kingsford 2011, Bioinformatics 27:764). Jellyfish's 2-bit comparison A&lt;C&lt;G&lt;T is
    /// ordinal string order on upper-case ACGT. Canonical counting is defined only over ACGT (Jellyfish
    /// cannot encode another base), so <c>Canonical = true</c> always applies the ACGT-only window rule.
    /// RC is the canonical <see cref="DnaSequence.GetReverseComplementString"/>, applied once per distinct
    /// forward k-mer. Keys are upper-case; the counts sum to the number of all-ACGT windows.</para>
    /// <para>Default options are exactly <see cref="CountKmers(string, int, CancellationToken, IProgress{double}?)"/>.
    /// Cancellation and progress semantics are those of that overload (the canonical fold runs after the final
    /// 1.0 report and is O(D·k)).</para>
    /// <para>Cross-checked against Jellyfish 2.3.1 (<c>count [-C]</c> + <c>dump -c</c>) on 20 inputs incl.
    /// N/IUPAC/lower-case/U (docs/algorithms/K-mer/K-mer_Counting.md §7.3).</para>
    /// </remarks>
    /// <param name="sequence">The sequence (case-insensitive). Null/empty returns an empty dictionary.</param>
    /// <param name="k">The k-mer length. Must be positive for non-empty input.</param>
    /// <param name="options">Counting mode.</param>
    /// <param name="cancellationToken">Cancellation token.</param>
    /// <param name="progress">Optional progress reporter (0.0 to 1.0).</param>
    /// <returns>Dictionary mapping (canonical) k-mers to their counts.</returns>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="k"/> ≤ 0 and the sequence is non-empty.</exception>
    public static Dictionary<string, int> CountKmers(
        string sequence,
        int k,
        KmerCountingOptions options,
        CancellationToken cancellationToken = default,
        IProgress<double>? progress = null)
    {
        var counts = CountKmersCore(sequence, k, options.SkipsNonAcgt, cancellationToken, progress);
        return options.Canonical ? FoldToCanonical(counts) : counts;
    }

    /// <summary>
    /// Returns the set of distinct k-mers of a sequence (literal counting, as <see cref="CountKmers(string, int)"/>).
    /// </summary>
    /// <param name="sequence">The sequence (case-insensitive). Null/empty returns an empty set.</param>
    /// <param name="k">The k-mer length. Must be positive for non-empty input.</param>
    /// <returns>A new caller-owned ordinal set of the observed k-mers (upper-case).</returns>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="k"/> ≤ 0 and the sequence is non-empty.</exception>
    public static HashSet<string> DistinctKmers(string sequence, int k)
        => DistinctKmers(sequence, k, default);

    /// <summary>
    /// Returns the set of distinct k-mers under <paramref name="options"/> — the key set of
    /// <see cref="CountKmers(string, int, KmerCountingOptions, CancellationToken, IProgress{double}?)"/>
    /// (with <c>Canonical = true</c>: the distinct canonical k-mers, Jellyfish <c>stats -C</c> "Distinct").
    /// </summary>
    /// <param name="sequence">The sequence (case-insensitive). Null/empty returns an empty set.</param>
    /// <param name="k">The k-mer length. Must be positive for non-empty input.</param>
    /// <param name="options">Counting mode.</param>
    /// <returns>A new caller-owned ordinal set of the observed (canonical) k-mers.</returns>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="k"/> ≤ 0 and the sequence is non-empty.</exception>
    public static HashSet<string> DistinctKmers(string sequence, int k, KmerCountingOptions options)
        => new(CountKmers(sequence, k, options).Keys, StringComparer.Ordinal);

    /// <summary>Folds a forward count table onto canonical keys min(w, RC(w)) (ordinal = Jellyfish 2-bit order).</summary>
    private static Dictionary<string, int> FoldToCanonical(Dictionary<string, int> forward)
    {
        var canonical = new Dictionary<string, int>(forward.Count);
        foreach (var (kmer, count) in forward)
        {
            var rc = DnaSequence.GetReverseComplementString(kmer);
            var key = string.CompareOrdinal(rc, kmer) < 0 ? rc : kmer;
            if (!canonical.TryAdd(key, count))
                canonical[key] += count;
        }
        return canonical;
    }

    private static bool IsAcgt(char c) => c is 'A' or 'C' or 'G' or 'T';

    /// <summary>
    /// The single counting loop. With <paramref name="acgtOnly"/>, windows containing a non-ACGT symbol
    /// (after upper-casing) are skipped — tracked in O(1) per window via the last non-ACGT index.
    /// </summary>
    private static Dictionary<string, int> CountKmersCore(
        string sequence,
        int k,
        bool acgtOnly,
        CancellationToken cancellationToken,
        IProgress<double>? progress)
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
        int scanned = 0, lastInvalid = -1;

        for (int i = 0; i <= sequence.Length - k; i++)
        {
            if (i % checkInterval == 0)
            {
                cancellationToken.ThrowIfCancellationRequested();
                progress?.Report((double)i / total);
            }

            if (acgtOnly)
            {
                for (; scanned < i + k; scanned++)
                {
                    if (!IsAcgt(seq[scanned]))
                        lastInvalid = scanned;
                }

                if (lastInvalid >= i)
                    continue;
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
        => GetKmerSpectrum(sequence, k, default);

    /// <summary>
    /// Gets the k-mer spectrum over counts produced under <paramref name="options"/> — with
    /// <c>Canonical = true</c> this is <c>jellyfish count -C</c> followed by <c>jellyfish histo</c>, the usual
    /// genome-profiling input (e.g. GenomeScope).
    /// </summary>
    /// <param name="sequence">The sequence to analyze. Null/empty returns an empty dictionary.</param>
    /// <param name="k">The k-mer length.</param>
    /// <param name="options">Counting mode (see <see cref="KmerCountingOptions"/>).</param>
    /// <returns>Dictionary mapping multiplicity to the number of distinct k-mers with that multiplicity.</returns>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="k"/> ≤ 0 and the sequence is non-empty.</exception>
    public static Dictionary<int, int> GetKmerSpectrum(string sequence, int k, KmerCountingOptions options)
    {
        var counts = CountKmers(sequence, k, options);
        var spectrum = new Dictionary<int, int>();

        foreach (var count in counts.Values)
        {
            if (!spectrum.TryAdd(count, 1))
                spectrum[count]++;
        }

        return spectrum;
    }

    /// <summary>Jellyfish <c>histo</c> default <c>--low</c> (1).</summary>
    public const long JellyfishHistoDefaultLow = 1;

    /// <summary>Jellyfish <c>histo</c> default <c>--high</c> (10000).</summary>
    public const long JellyfishHistoDefaultHigh = 10000;

    /// <summary>
    /// Jellyfish-<c>histo</c>-compatible k-mer histogram over counts produced under <paramref name="options"/>:
    /// <c>jellyfish count [-C] -m k</c> followed by <c>jellyfish histo -l low -h high -i increment [-f]</c>.
    /// </summary>
    /// <remarks>
    /// Delegates to <see cref="GetKmerHistogram(IEnumerable{int}, long, long, long, bool)"/> on the count table of
    /// <see cref="CountKmers(string, int, KmerCountingOptions, CancellationToken, IProgress{double}?)"/>; see that
    /// overload for the exact bucket rule. With the defaults (low 1, high 10000, increment 1, not full) the rows
    /// are the non-zero bins of <see cref="GetKmerSpectrum(string, int, KmerCountingOptions)"/> in ascending order,
    /// except that multiplicities above 10000 are pooled in the cap bin 10001.
    /// </remarks>
    /// <param name="sequence">The sequence to analyze. Null/empty yields no k-mers.</param>
    /// <param name="k">The k-mer length.</param>
    /// <param name="options">Counting mode (literal, ACGT-only or canonical <c>-C</c>).</param>
    /// <param name="low">Jellyfish <c>-l/--low</c> (default 1).</param>
    /// <param name="high">Jellyfish <c>-h/--high</c> (default 10000).</param>
    /// <param name="increment">Jellyfish <c>-i/--increment</c> (default 1).</param>
    /// <param name="full">Jellyfish <c>-f/--full</c>: also emit empty bins.</param>
    /// <returns>(Bin, Frequency) rows in ascending bin order, exactly the lines <c>jellyfish histo</c> prints.</returns>
    /// <exception cref="ArgumentOutOfRangeException">k ≤ 0 for non-empty input, or invalid histo parameters.</exception>
    public static IReadOnlyList<KmerHistogramBin> GetKmerHistogram(
        string sequence,
        int k,
        KmerCountingOptions options = default,
        long low = JellyfishHistoDefaultLow,
        long high = JellyfishHistoDefaultHigh,
        long increment = 1,
        bool full = false)
    {
        var layout = HistoLayout.Create(low, high, increment, full);
        return BuildHistogram(CountKmers(sequence, k, options).Values, layout);
    }

    /// <summary>
    /// Jellyfish-<c>histo</c>-compatible histogram of a k-mer count table (one entry per distinct k-mer),
    /// e.g. the values of a <see cref="CountKmers(string, int, KmerCountingOptions, CancellationToken, IProgress{double}?)"/> table.
    /// </summary>
    /// <remarks>
    /// Mirrors <c>sub_commands/histo_main.cc</c> (gmarcais/Jellyfish) exactly:
    /// <code>
    /// base = increment &gt;= low ? 0 : low − increment
    /// ceil = high + increment
    /// nb_buckets = (ceil + increment − base) / increment        (integer division)
    /// count &lt; base → bucket 0;  count &gt; ceil → bucket nb_buckets − 1;  else bucket (count − base) / increment
    /// bucket i is labelled base + i·increment; rows with frequency 0 are printed only with --full
    /// </code>
    /// So the last bucket (label ≥ high) is the catch-all cap bin for every count above <paramref name="high"/>,
    /// and counts below <paramref name="low"/> are pooled in the first bucket (when increment ≥ low the base is 0
    /// and the first label is 0). Cross-checked against the Jellyfish 2.3.1 binary on 72 runs
    /// (docs/algorithms/K-mer/K-mer_Frequency_Analysis.md §7.3). O(D) time for D counts; O(number of non-empty
    /// buckets) memory, O(nb_buckets) output with <paramref name="full"/>.
    /// </remarks>
    /// <param name="kmerCounts">Multiplicity of each distinct k-mer (each ≥ 0).</param>
    /// <param name="low">Jellyfish <c>-l/--low</c> (≥ 0, default 1).</param>
    /// <param name="high">Jellyfish <c>-h/--high</c> (≥ <paramref name="low"/>, default 10000).</param>
    /// <param name="increment">Jellyfish <c>-i/--increment</c> (≥ 1, default 1).</param>
    /// <param name="full">Jellyfish <c>-f/--full</c>: also emit empty bins.</param>
    /// <returns>(Bin, Frequency) rows in ascending bin order.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="kmerCounts"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException">
    /// <paramref name="low"/> &lt; 0, <paramref name="high"/> &lt; <paramref name="low"/> (Jellyfish: "High count value
    /// must be &gt;= to low count value"), <paramref name="increment"/> &lt; 1 (Jellyfish divides by it),
    /// high + 2·increment overflows, <paramref name="full"/> with more buckets than an array can hold,
    /// or a negative count.
    /// </exception>
    public static IReadOnlyList<KmerHistogramBin> GetKmerHistogram(
        IEnumerable<int> kmerCounts,
        long low = JellyfishHistoDefaultLow,
        long high = JellyfishHistoDefaultHigh,
        long increment = 1,
        bool full = false)
    {
        ArgumentNullException.ThrowIfNull(kmerCounts);
        return BuildHistogram(kmerCounts, HistoLayout.Create(low, high, increment, full));
    }

    /// <summary>Bucket layout of Jellyfish <c>histo_main.cc</c> (base, ceil, nb_buckets, inc).</summary>
    private readonly record struct HistoLayout(long Base, long Ceil, long Increment, long BucketCount, bool Full)
    {
        public static HistoLayout Create(long low, long high, long increment, bool full)
        {
            ArgumentOutOfRangeException.ThrowIfNegative(low);
            ArgumentOutOfRangeException.ThrowIfLessThan(increment, 1L);
            if (high < low)
                throw new ArgumentOutOfRangeException(nameof(high), high, "High count value must be >= low count value.");
            if (increment > (long.MaxValue - high) / 2)
                throw new ArgumentOutOfRangeException(nameof(high), high, "high + 2·increment overflows.");

            long baseCount = increment >= low ? 0 : low - increment;
            long ceil = high + increment;
            long buckets = (ceil + increment - baseCount) / increment;
            if (full && buckets > Array.MaxLength)
                throw new ArgumentOutOfRangeException(nameof(full), full, "A full histogram would have more buckets than an array can hold.");
            return new HistoLayout(baseCount, ceil, increment, buckets, full);
        }

        public long BucketOf(long count)
        {
            if (count < Base)
                return 0;
            if (count > Ceil)
                return BucketCount - 1;
            return (count - Base) / Increment;
        }
    }

    private static List<KmerHistogramBin> BuildHistogram(IEnumerable<int> kmerCounts, HistoLayout layout)
    {
        var tally = new Dictionary<long, long>();
        foreach (int count in kmerCounts)
        {
            if (count < 0)
                throw new ArgumentOutOfRangeException(nameof(kmerCounts), count, "k-mer counts must be non-negative.");
            long bucket = layout.BucketOf(count);
            tally[bucket] = tally.TryGetValue(bucket, out long n) ? n + 1 : 1;
        }

        if (!layout.Full)
        {
            return tally.OrderBy(kvp => kvp.Key)
                .Select(kvp => new KmerHistogramBin(layout.Base + kvp.Key * layout.Increment, kvp.Value))
                .ToList();
        }

        var rows = new List<KmerHistogramBin>((int)layout.BucketCount);
        for (long i = 0; i < layout.BucketCount; i++)
            rows.Add(new KmerHistogramBin(layout.Base + i * layout.Increment, tally.GetValueOrDefault(i)));
        return rows;
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

        // Single word-vector loop: the metric overload (frequency Euclidean over the union of k-mers).
        return KmerDistance(seq1, seq2, k, KmerDistanceMetric.Euclidean);
    }

    /// <summary>
    /// Alignment-free word-vector dissimilarity between two sequences under an explicit
    /// <see cref="KmerDistanceMetric"/> (literal k-mer counting, as <see cref="CountKmers(string, int)"/>).
    /// </summary>
    /// <remarks>
    /// Counts both sequences with the canonical counter and delegates to
    /// <see cref="KmerDistance(IReadOnlyDictionary{string, int}, IReadOnlyDictionary{string, int}, KmerDistanceMetric)"/>;
    /// the per-metric vector (raw counts or relative frequencies) is documented on <see cref="KmerDistanceMetric"/>.
    /// <see cref="KmerDistanceMetric.Euclidean"/> is bit-for-bit <see cref="KmerDistance(string, string, int)"/>.
    /// </remarks>
    /// <param name="seq1">First sequence (case-insensitive). Null/empty or shorter than <paramref name="k"/> gives the zero vector.</param>
    /// <param name="seq2">Second sequence, same conventions.</param>
    /// <param name="k">K-mer length; must be positive.</param>
    /// <param name="metric">The word-vector metric.</param>
    /// <returns>The metric value (a dissimilarity for every member except <see cref="KmerDistanceMetric.D2"/>, a similarity).</returns>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="k"/> is not positive or <paramref name="metric"/> is undefined.</exception>
    public static double KmerDistance(string seq1, string seq2, int k, KmerDistanceMetric metric)
    {
        if (k <= 0)
            throw new ArgumentOutOfRangeException(nameof(k), "K must be positive.");

        return KmerDistance(CountKmers(seq1, k), CountKmers(seq2, k), metric);
    }

    /// <summary>
    /// Word-vector dissimilarity between two count tables (k-mer counts, Jellyfish <c>-C</c> canonical counts
    /// from <see cref="CountKmers(string, int, KmerCountingOptions, CancellationToken, IProgress{double}?)"/>,
    /// or spaced-word counts from <see cref="CountSpacedWords(string, string)"/>).
    /// </summary>
    /// <remarks>
    /// <para>The vectors span the union of keys of the two tables; a word absent from a table is a 0 component,
    /// which equals the full |Σ|^k-dimensional vector because words absent from both contribute 0 to every
    /// metric. Frequencies are f(w) = c(w) / Σc (Σc = L − k + 1 windows for contiguous k-mers, L − ℓ + 1 for a
    /// spaced pattern of length ℓ; scikit-bio <c>kmer_frequencies(relative=True)</c>, alfpy <c>word_vector.Freqs</c>);
    /// an empty table is the zero vector.</para>
    /// <para>Definitions (see <see cref="KmerDistanceMetric"/>): Euclidean √Σ(f₁−f₂)² (Vinga &amp; Almeida 2003;
    /// Zielezinski et al. 2017); squared Euclidean on counts Σ(c₁−c₂)² (Blaisdell 1986, d_E); Manhattan Σ|f₁−f₂|;
    /// Chebyshev max|f₁−f₂|; Canberra Σ|f₁−f₂|/(f₁+f₂) (0/0 terms omitted, as scipy and alfpy); cosine distance
    /// 1 − c₁·c₂/(‖c₁‖‖c₂‖) clipped to [0, 2] as scipy (scale-invariant, so counts and frequencies agree; a zero
    /// vector has cosine similarity 0, so its distance is 1); D2 = Σ c₁(w)·c₂(w) (Torney et al. 1990; Lippert
    /// et al. 2005; Reinert et al. 2009). Cross-checked against scipy.spatial.distance and alfpy 1.0.6
    /// (docs/algorithms/K-mer/K-mer_Euclidean_Distance.md §7.2).</para>
    /// </remarks>
    /// <param name="counts1">First count table (non-negative counts).</param>
    /// <param name="counts2">Second count table (non-negative counts).</param>
    /// <param name="metric">The word-vector metric.</param>
    /// <returns>The metric value.</returns>
    /// <exception cref="ArgumentNullException">A table is null.</exception>
    /// <exception cref="ArgumentException">A table contains a negative count.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="metric"/> is undefined.</exception>
    public static double KmerDistance(
        IReadOnlyDictionary<string, int> counts1,
        IReadOnlyDictionary<string, int> counts2,
        KmerDistanceMetric metric)
    {
        ArgumentNullException.ThrowIfNull(counts1);
        ArgumentNullException.ThrowIfNull(counts2);
        if (!Enum.IsDefined(metric))
            throw new ArgumentOutOfRangeException(nameof(metric), metric, "Unknown k-mer distance metric.");

        double total1 = SumNonNegative(counts1, nameof(counts1));
        double total2 = SumNonNegative(counts2, nameof(counts2));

        double acc = 0, dot = 0, norm1 = 0, norm2 = 0;

        // Union of keys: counts1's keys in table order, then counts2's keys absent from counts1
        // (the iteration order of the original frequency-Euclidean implementation, so its sum is bit-identical).
        foreach (var (word, c1) in counts1)
            Accumulate(c1, counts2.TryGetValue(word, out var c2) ? c2 : 0);
        foreach (var (word, c2) in counts2)
        {
            if (!counts1.ContainsKey(word))
                Accumulate(0, c2);
        }

        return metric switch
        {
            KmerDistanceMetric.Euclidean => Math.Sqrt(acc),
            KmerDistanceMetric.Cosine => norm1 == 0 || norm2 == 0
                ? 1.0
                : Math.Clamp(1.0 - dot / (Math.Sqrt(norm1) * Math.Sqrt(norm2)), 0.0, 2.0),
            KmerDistanceMetric.D2 => dot,
            _ => acc,
        };

        void Accumulate(int c1, int c2)
        {
            double f1 = total1 == 0 ? 0 : c1 / total1;
            double f2 = total2 == 0 ? 0 : c2 / total2;
            switch (metric)
            {
                case KmerDistanceMetric.Euclidean:
                    acc += (f1 - f2) * (f1 - f2);
                    break;
                case KmerDistanceMetric.SquaredEuclideanCounts:
                    acc += ((double)c1 - c2) * ((double)c1 - c2);
                    break;
                case KmerDistanceMetric.Manhattan:
                    acc += Math.Abs(f1 - f2);
                    break;
                case KmerDistanceMetric.Chebyshev:
                    acc = Math.Max(acc, Math.Abs(f1 - f2));
                    break;
                case KmerDistanceMetric.Canberra:
                    if (f1 + f2 > 0)
                        acc += Math.Abs(f1 - f2) / (f1 + f2);
                    break;
                default: // Cosine, D2: inner product and norms of the count vectors
                    dot += (double)c1 * c2;
                    norm1 += (double)c1 * c1;
                    norm2 += (double)c2 * c2;
                    break;
            }
        }
    }

    private static double SumNonNegative(IReadOnlyDictionary<string, int> counts, string paramName)
    {
        long total = 0;
        foreach (var count in counts.Values)
        {
            if (count < 0)
                throw new ArgumentException("Counts must be non-negative.", paramName);
            total += count;
        }
        return total;
    }

    /// <summary>
    /// Exact k-mer Jaccard similarity J(A, B) = |K(A) ∩ K(B)| / |K(A) ∪ K(B)| of the two distinct k-mer sets
    /// (literal counting, as <see cref="DistinctKmers(string, int)"/>), as a fraction in [0, 1].
    /// </summary>
    /// <remarks>
    /// Jaccard (1901, 1912) coefficient of community applied to k-mer sets, the quantity Mash estimates by
    /// MinHash (Ondov et al. 2016, Genome Biol 17:132, eq. 1). Canonical site of the k-mer Jaccard index
    /// (DUP_MAP §16). Case-insensitive. Both sets empty (both sequences null/empty or shorter than
    /// <paramref name="k"/>): the index is undefined (0/0) and 0 is returned — the convention of
    /// <c>GenomicAnalyzer.CalculateSimilarity</c> and <c>ComparativeGenomics</c>' k-mer identity.
    /// </remarks>
    /// <param name="a">First sequence; null is treated as empty.</param>
    /// <param name="b">Second sequence; null is treated as empty.</param>
    /// <param name="k">K-mer length; must be positive.</param>
    /// <returns>The Jaccard index in [0, 1]; 1 for identical non-empty k-mer sets.</returns>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="k"/> is not positive.</exception>
    public static double JaccardSimilarity(string a, string b, int k)
        => JaccardSimilarity(a, b, k, default);

    /// <summary>
    /// Exact k-mer Jaccard similarity over the distinct k-mer sets taken under <paramref name="options"/>.
    /// With <c>Canonical = true</c> this is the exact (unsketched) Jaccard index that <c>mash dist</c>
    /// and sourmash (DNA MinHash, <c>scaled=1</c>) estimate: k-mers upper-cased, windows with a non-ACGT
    /// base skipped, each k-mer keyed by min(w, RC(w)) (Mash <c>Sketch.cpp</c> <c>addMinHashes</c>).
    /// </summary>
    /// <param name="a">First sequence; null is treated as empty.</param>
    /// <param name="b">Second sequence; null is treated as empty.</param>
    /// <param name="k">K-mer length; must be positive.</param>
    /// <param name="options">Counting mode for both sequences.</param>
    /// <returns>The Jaccard index in [0, 1]; 0 when both sets are empty.</returns>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="k"/> is not positive.</exception>
    public static double JaccardSimilarity(string a, string b, int k, KmerCountingOptions options)
    {
        var (shared, union) = SharedAndUnionKmers(a, b, k, options);
        return union == 0 ? 0.0 : (double)shared / union;
    }

    /// <summary>
    /// Mash distance D = −(1/k)·ln(2J/(1 + J)) (Ondov et al. 2016, eq. 4) computed from the exact k-mer
    /// Jaccard index under <paramref name="options"/>. Pass <c>new KmerCountingOptions(Canonical: true)</c>
    /// for <c>mash dist</c>'s default, <c>AcgtOnly: true</c> for <c>mash sketch -n</c> (non-canonical).
    /// </summary>
    /// <remarks>
    /// Boundary rules of Mash 2.x <c>CommandDistance.cpp</c>: shared = union (identical sets, including both empty)
    /// → 0; no shared k-mer → 1; otherwise the formula, capped at 1. Mash estimates J from a bottom-s MinHash
    /// sketch; this method uses the exact sets, so it equals <c>mash dist -s s</c> whenever s ≥ |K(A) ∪ K(B)|.
    /// </remarks>
    /// <param name="a">First sequence; null is treated as empty.</param>
    /// <param name="b">Second sequence; null is treated as empty.</param>
    /// <param name="k">K-mer length; must be positive.</param>
    /// <param name="options">Counting mode for both sequences.</param>
    /// <returns>The Mash distance in [0, 1].</returns>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="k"/> is not positive.</exception>
    public static double MashDistance(string a, string b, int k, KmerCountingOptions options)
    {
        var (shared, union) = SharedAndUnionKmers(a, b, k, options);
        if (shared == union)
            return 0.0;
        return MashDistanceFromJaccard((double)shared / union, k);
    }

    /// <summary>
    /// Converts a k-mer Jaccard index to the Mash distance D = −(1/k)·ln(2J/(1 + J)) (Ondov et al. 2016, eq. 4),
    /// with Mash's boundary rules: J = 1 → 0, J = 0 → 1, results above 1 capped at 1.
    /// </summary>
    /// <param name="jaccard">Jaccard index in [0, 1].</param>
    /// <param name="k">K-mer length; must be positive.</param>
    /// <returns>The Mash distance in [0, 1].</returns>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="jaccard"/> is outside [0, 1] (or NaN), or <paramref name="k"/> is not positive.</exception>
    public static double MashDistanceFromJaccard(double jaccard, int k)
    {
        if (k <= 0)
            throw new ArgumentOutOfRangeException(nameof(k), "K must be positive.");
        if (!(jaccard >= 0.0 && jaccard <= 1.0))
            throw new ArgumentOutOfRangeException(nameof(jaccard), jaccard, "Jaccard index must be in [0, 1].");

        if (jaccard == 1.0)
            return 0.0;
        if (jaccard == 0.0)
            return 1.0;

        double distance = -Math.Log(2 * jaccard / (1.0 + jaccard)) / k;
        return Math.Min(distance, 1.0);
    }

    private static (int Shared, int Union) SharedAndUnionKmers(string a, string b, int k, KmerCountingOptions options)
    {
        if (k <= 0)
            throw new ArgumentOutOfRangeException(nameof(k), "K must be positive.");

        var setA = DistinctKmers(a ?? string.Empty, k, options);
        var setB = DistinctKmers(b ?? string.Empty, k, options);
        var (small, large) = setA.Count <= setB.Count ? (setA, setB) : (setB, setA);

        int shared = 0;
        foreach (var kmer in small)
        {
            if (large.Contains(kmer))
                shared++;
        }
        return (shared, setA.Count + setB.Count - shared);
    }

    /// <summary>
    /// Counts the spaced words of <paramref name="sequence"/> with respect to a binary match pattern
    /// (Leimeister, Boden, Horwege, Lindner &amp; Morgenstern 2014, Bioinformatics 30:1991).
    /// </summary>
    /// <remarks>
    /// A pattern P ∈ {0,1}^ℓ with P[1] = P[ℓ] = 1 has weight k = number of '1' (match) positions; '0' positions are
    /// "don't care". The spaced word at window i (0 ≤ i ≤ L − ℓ) is the k-symbol string of
    /// sequence[i + j] over the match positions j, in order; each window contributes one occurrence, so the
    /// counts sum to L − ℓ + 1. The all-'1' pattern of length k gives exactly <see cref="CountKmers(string, int)"/>.
    /// Counting is literal and case-insensitive (upper-cased keys), like <see cref="CountKmers(string, int)"/>.
    /// Compare tables with <see cref="KmerDistance(IReadOnlyDictionary{string, int}, IReadOnlyDictionary{string, int}, KmerDistanceMetric)"/>
    /// (Leimeister et al. use the Euclidean distance of spaced-word frequency vectors).
    /// </remarks>
    /// <param name="sequence">The sequence; null/empty or shorter than the pattern gives an empty table.</param>
    /// <param name="pattern">Binary pattern over {'0','1'} that starts and ends with '1', e.g. "1101".</param>
    /// <returns>Dictionary mapping spaced words (match-position symbols) to their counts.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="pattern"/> is null.</exception>
    /// <exception cref="ArgumentException"><paramref name="pattern"/> is empty, contains a symbol other than '0'/'1',
    /// or does not start and end with '1'.</exception>
    public static Dictionary<string, int> CountSpacedWords(string sequence, string pattern)
    {
        ArgumentNullException.ThrowIfNull(pattern);
        if (pattern.Length == 0 || pattern[0] != '1' || pattern[^1] != '1' || pattern.Any(c => c is not ('0' or '1')))
            throw new ArgumentException("Pattern must be a non-empty string over {0,1} that starts and ends with '1'.", nameof(pattern));

        var counts = new Dictionary<string, int>();
        if (string.IsNullOrEmpty(sequence) || pattern.Length > sequence.Length)
            return counts;

        int[] matchPositions = Enumerable.Range(0, pattern.Length).Where(j => pattern[j] == '1').ToArray();
        var seq = sequence.ToUpperInvariant();
        var buffer = new char[matchPositions.Length];

        for (int i = 0; i <= seq.Length - pattern.Length; i++)
        {
            for (int m = 0; m < matchPositions.Length; m++)
                buffer[m] = seq[i + matchPositions[m]];

            var word = new string(buffer);
            if (!counts.TryAdd(word, 1))
                counts[word]++;
        }

        return counts;
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
        => FindKmersWithMinCount(sequence, k, minCount, maxCount, KmerCountingOptions.Default);

    /// <summary>
    /// Unique k-mers (count exactly 1) under <paramref name="options"/> — with <c>Canonical = true</c> this is
    /// <c>jellyfish count -C</c> followed by <c>jellyfish dump -L 1 -U 1</c>.
    /// </summary>
    /// <param name="sequence">The sequence to analyze (case-insensitive).</param>
    /// <param name="k">The k-mer length. Must be positive.</param>
    /// <param name="options">Counting mode (literal, ACGT-only, canonical).</param>
    /// <returns>The (canonical) k-mers whose count equals 1, in ascending ordinal order.</returns>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="k"/> ≤ 0 and the sequence is non-empty.</exception>
    public static IEnumerable<string> FindUniqueKmers(string sequence, int k, KmerCountingOptions options)
        => FindKmersWithMinCount(sequence, k, UniqueKmerCount, UniqueKmerCount, options).Select(p => p.Kmer);

    /// <summary>
    /// k-mers with count ≥ <paramref name="minCount"/> under <paramref name="options"/>
    /// (<c>jellyfish count [-C]</c> + <c>dump -L minCount</c>).
    /// </summary>
    /// <param name="sequence">The sequence to analyze (case-insensitive).</param>
    /// <param name="k">The k-mer length. Must be positive.</param>
    /// <param name="minCount">Inclusive lower count bound (≤ 1 imposes no lower bound).</param>
    /// <param name="options">Counting mode (literal, ACGT-only, canonical).</param>
    /// <returns>(k-mer, Count) pairs, count descending then ordinal k-mer ascending.</returns>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="k"/> ≤ 0 and the sequence is non-empty.</exception>
    public static IEnumerable<(string Kmer, int Count)> FindKmersWithMinCount(
        string sequence, int k, int minCount, KmerCountingOptions options)
        => FindKmersWithMinCount(sequence, k, minCount, int.MaxValue, options);

    /// <summary>
    /// k-mers whose count lies in [<paramref name="minCount"/>, <paramref name="maxCount"/>] under
    /// <paramref name="options"/> — <c>jellyfish count [-C]</c> followed by <c>jellyfish dump -L -U</c>
    /// (Jellyfish dump filters are normally applied to a <c>-C</c> database).
    /// </summary>
    /// <remarks>
    /// The single implementation behind every unique/min-count overload: the option-aware
    /// <see cref="CountKmers(string, int, KmerCountingOptions, CancellationToken, IProgress{double}?)"/> table
    /// filtered by the shared <c>SelectByCountRange</c> predicate. Default options give exactly
    /// <see cref="FindKmersWithMinCount(string, int, int, int)"/>. Cross-checked against Jellyfish 2.3.1
    /// <c>count [-C]</c> + <c>dump -c -L -U</c> (Unique_And_MinCount_Kmers.md §7.3).
    /// </remarks>
    /// <param name="sequence">The sequence to analyze (case-insensitive).</param>
    /// <param name="k">The k-mer length. Must be positive.</param>
    /// <param name="minCount">Inclusive lower count bound (≤ 1 imposes no lower bound).</param>
    /// <param name="maxCount">Inclusive upper count bound; <see cref="int.MaxValue"/> = unbounded.</param>
    /// <param name="options">Counting mode (literal, ACGT-only, canonical).</param>
    /// <returns>(k-mer, Count) pairs in the range, count descending then ordinal k-mer ascending.</returns>
    /// <exception cref="ArgumentOutOfRangeException">
    /// <paramref name="maxCount"/> is negative, or <paramref name="k"/> ≤ 0 and the sequence is non-empty.
    /// </exception>
    public static IEnumerable<(string Kmer, int Count)> FindKmersWithMinCount(
        string sequence, int k, int minCount, int maxCount, KmerCountingOptions options)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(maxCount);

        var counts = CountKmers(sequence, k, options);
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
    /// The pass is the one shared with <see cref="FindClumpWindows"/> (which also reports the qualifying
    /// windows); here only its run openings are used and results are streamed.
    /// Matching is case-insensitive (upper-cased), mirroring <see cref="CountKmers(string,int)"/>.
    /// Cross-check: E. coli genome (textbook dataset), k=9, L=500, t=3 → 1904 distinct 9-mers.
    /// </remarks>
    /// <param name="sequence">The sequence to analyze (Genome).</param>
    /// <param name="k">K-mer length.</param>
    /// <param name="windowSize">Window length L.</param>
    /// <param name="minOccurrences">Minimum occurrences t within one window (inclusive).</param>
    /// <returns>
    /// Each clump-forming k-mer exactly once, in order of first detection (first window: ordinal order;
    /// the order is not part of the contract). Empty when the sequence is null/empty, k ≤ 0, L &lt; k, L &gt; |sequence|
    /// (no window of length L exists) or t ≤ 0.
    /// </returns>
    public static IEnumerable<string> FindClumps(string sequence, int k, int windowSize, int minOccurrences)
    {
        if (!IsClumpScanPossible(sequence, k, windowSize, minOccurrences))
            yield break;

        var clumps = new HashSet<string>(StringComparer.Ordinal);
        foreach (var transition in ScanClumpTransitions(sequence.ToUpperInvariant(), k, windowSize, minOccurrences))
        {
            if (transition.Qualifies && clumps.Add(transition.Kmer))
                yield return transition.Kmer;
        }
    }

    /// <summary>
    /// Finds every (L, t)-clump k-mer together with the windows in which it forms a clump: for each k-mer, the
    /// maximal runs of consecutive window starts i such that <c>Genome[i..i+L−1]</c> contains at least t of its
    /// occurrences.
    /// </summary>
    /// <remarks>
    /// Same clump definition and window convention as <see cref="FindClumps"/> (Compeau &amp; Pevzner,
    /// <i>Bioinformatics Algorithms</i>, ch. 1; Rosalind BA1E): windows start at i ∈ [0, |Genome| − L] and an
    /// occurrence at p counts when i ≤ p ≤ i + L − k. The qualifying windows of a k-mer form a set of window starts;
    /// it is returned as its maximal runs [<see cref="ClumpWindowRun.FirstWindowStart"/>,
    /// <see cref="ClumpWindowRun.LastWindowStart"/>] (inclusive). A run covers the genomic interval
    /// [FirstWindowStart, LastWindowStart + L − 1], which holds the clump. <see cref="KmerClump.FirstWindowStart"/>
    /// is the leftmost window in which the k-mer qualifies, i.e. where the textbook algorithm first detects it.
    /// <para>One streaming pass shared with <see cref="FindClumps"/> (<c>BetterClumpFinding</c>): after each 1-bp
    /// slide only the leaving k-mer's count can fall and only the entering k-mer's count can rise, so a run
    /// opens when the entering k-mer reaches t and closes when the leaving k-mer drops to t − 1 (a slide whose
    /// leaving and entering k-mers are equal changes nothing). O(|Genome|·k) time, O(L) window state plus the
    /// output. Cross-checked against a Python brute force over every window (Find_Tests / K-mer_Search.md §7).</para>
    /// </remarks>
    /// <param name="sequence">The sequence to analyze (Genome, case-insensitive).</param>
    /// <param name="k">K-mer length.</param>
    /// <param name="windowSize">Window length L.</param>
    /// <param name="minOccurrences">Minimum occurrences t within one window (inclusive).</param>
    /// <returns>
    /// One <see cref="KmerClump"/> per clump k-mer (the same set as <see cref="FindClumps"/>), ordered by
    /// <see cref="KmerClump.FirstWindowStart"/> then ordinal k-mer; runs in ascending order. Empty under the same
    /// conditions as <see cref="FindClumps"/>.
    /// </returns>
    public static IReadOnlyList<KmerClump> FindClumpWindows(string sequence, int k, int windowSize, int minOccurrences)
    {
        var clumps = new List<KmerClump>();
        if (!IsClumpScanPossible(sequence, k, windowSize, minOccurrences))
            return clumps;

        var runsByKmer = new Dictionary<string, List<ClumpWindowRun>>(StringComparer.Ordinal);
        var openRunStart = new Dictionary<string, int>(StringComparer.Ordinal);
        foreach (var (kmer, window, qualifies) in ScanClumpTransitions(sequence.ToUpperInvariant(), k, windowSize, minOccurrences))
        {
            if (qualifies)
            {
                openRunStart[kmer] = window;
                if (!runsByKmer.ContainsKey(kmer))
                {
                    var runs = new List<ClumpWindowRun>();
                    runsByKmer[kmer] = runs;
                    clumps.Add(new KmerClump(kmer, runs));
                }
            }
            else
            {
                runsByKmer[kmer].Add(new ClumpWindowRun(openRunStart[kmer], window - 1));
                openRunStart.Remove(kmer);
            }
        }

        return clumps;
    }

    private static bool IsClumpScanPossible(string sequence, int k, int windowSize, int minOccurrences)
        => !string.IsNullOrEmpty(sequence) && k > 0 && windowSize >= k && minOccurrences > 0
           && windowSize <= sequence.Length;

    /// <summary>
    /// A change of a k-mer's clump status at a window start: <c>Qualifies = true</c> — the window
    /// <c>Window</c> is the first of a run in which the k-mer has ≥ t occurrences; <c>false</c> — window
    /// <c>Window</c> is the first after the run (the run ended at <c>Window − 1</c>).
    /// </summary>
    private readonly record struct ClumpTransition(string Kmer, int Window, bool Qualifies);

    /// <summary>
    /// The single (L, t)-clump sliding-window pass (Compeau &amp; Pevzner <c>BetterClumpFinding</c>) behind
    /// <see cref="FindClumps"/> and <see cref="FindClumpWindows"/>. Emits run openings in window order (the first
    /// window's in ordinal k-mer order), run closings when a count falls from t to t − 1, and after the last window
    /// a closing at |Genome| − L + 1 for every still-open run (ordinal order). <paramref name="seq"/> is upper case.
    /// </summary>
    private static IEnumerable<ClumpTransition> ScanClumpTransitions(string seq, int k, int windowSize, int minOccurrences)
    {
        // First window Genome[0..L−1]: canonical k-mer counting.
        var windowCounts = CountKmers(seq.Substring(0, windowSize), k);
        foreach (var kmer in windowCounts.Where(kvp => kvp.Value >= minOccurrences)
                     .Select(kvp => kvp.Key).Order(StringComparer.Ordinal).ToList())
            yield return new ClumpTransition(kmer, 0, true);

        // Slide: window i covers k-mer starts i..i+L−k.
        int lastWindow = seq.Length - windowSize;
        for (int i = 1; i <= lastWindow; i++)
        {
            string leaving = seq.Substring(i - 1, k);
            string entering = seq.Substring(i + windowSize - k, k);
            if (string.Equals(leaving, entering, StringComparison.Ordinal))
                continue;

            int left = windowCounts[leaving] - 1;
            if (left == 0)
                windowCounts.Remove(leaving);
            else
                windowCounts[leaving] = left;
            if (left == minOccurrences - 1)
                yield return new ClumpTransition(leaving, i, false);

            int count = windowCounts.TryGetValue(entering, out int c) ? c + 1 : 1;
            windowCounts[entering] = count;
            if (count == minOccurrences)
                yield return new ClumpTransition(entering, i, true);
        }

        foreach (var kmer in windowCounts.Where(kvp => kvp.Value >= minOccurrences)
                     .Select(kvp => kvp.Key).Order(StringComparer.Ordinal).ToList())
            yield return new ClumpTransition(kmer, lastWindow + 1, false);
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
        => AnalyzeKmers(sequence, k, default, lowerCount, upperCount);

    /// <summary>
    /// Computes k-mer statistics over counts produced under <paramref name="options"/>, with the optional
    /// Jellyfish <c>-L/-U</c> filters — with <c>Canonical = true</c> this is <c>jellyfish count -C</c> followed
    /// by <c>jellyfish stats</c>. Field semantics are those of <see cref="AnalyzeKmers(string, int, int, int)"/>;
    /// Total is the number of counted (all-ACGT, when filtering) windows.
    /// </summary>
    /// <param name="sequence">The sequence to analyze (case-insensitive).</param>
    /// <param name="k">The k-mer length. Must be positive.</param>
    /// <param name="options">Counting mode (see <see cref="KmerCountingOptions"/>).</param>
    /// <param name="lowerCount">Ignore k-mers with count below this value (Jellyfish default 0).</param>
    /// <param name="upperCount">Ignore k-mers with count above this value (Jellyfish default unbounded).</param>
    /// <returns>Statistics over the retained k-mers; all-zero when none is retained.</returns>
    /// <exception cref="ArgumentOutOfRangeException">
    /// A bound is negative, or <paramref name="k"/> ≤ 0 and the sequence is non-empty.
    /// </exception>
    public static KmerStatistics AnalyzeKmers(
        string sequence, int k, KmerCountingOptions options, int lowerCount = 0, int upperCount = int.MaxValue)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(lowerCount);
        ArgumentOutOfRangeException.ThrowIfNegative(upperCount);

        var counts = CountKmers(sequence, k, options);
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

/// <summary>
/// A maximal run of consecutive window starts in which a k-mer forms an (L, t)-clump
/// (<see cref="KmerAnalyzer.FindClumpWindows"/>); both ends inclusive, 0-based.
/// </summary>
/// <param name="FirstWindowStart">First window start i of the run.</param>
/// <param name="LastWindowStart">Last window start of the run; the run covers Genome[FirstWindowStart..LastWindowStart + L − 1].</param>
public readonly record struct ClumpWindowRun(int FirstWindowStart, int LastWindowStart);

/// <summary>
/// An (L, t)-clump k-mer and the windows in which it forms a clump (<see cref="KmerAnalyzer.FindClumpWindows"/>).
/// </summary>
/// <param name="Kmer">The clump-forming k-mer (upper case).</param>
/// <param name="WindowRuns">Maximal runs of qualifying window starts, ascending and disjoint (never empty).</param>
public sealed record KmerClump(string Kmer, IReadOnlyList<ClumpWindowRun> WindowRuns)
{
    /// <summary>Leftmost window start in which the k-mer has at least t occurrences.</summary>
    public int FirstWindowStart => WindowRuns[0].FirstWindowStart;
}

/// <summary>
/// One row of a Jellyfish-<c>histo</c>-compatible k-mer histogram
/// (<see cref="KmerAnalyzer.GetKmerHistogram(string, int, KmerCountingOptions, long, long, long, bool)"/>).
/// </summary>
/// <param name="Bin">Bucket label: its low end point, base + i·increment (Jellyfish first output column).</param>
/// <param name="Frequency">Number of distinct k-mers tallied in the bucket (Jellyfish second column).</param>
public readonly record struct KmerHistogramBin(long Bin, long Frequency);

/// <summary>
/// K-mer counting mode for the option-aware <see cref="KmerAnalyzer"/> overloads. The default value
/// (both flags false) is the library's literal counting: every symbol, including N/IUPAC, forms k-mers.
/// </summary>
/// <param name="Canonical">
/// Key each k-mer by min(w, RC(w)) in ordinal (= A&lt;C&lt;G&lt;T) order — Jellyfish <c>count -C</c>
/// (Marçais &amp; Kingsford 2011). Implies the ACGT-only window rule, because the canonical form is defined
/// only over ACGT.
/// </param>
/// <param name="AcgtOnly">
/// Skip every window that contains a symbol other than A/C/G/T after case folding — Jellyfish
/// <c>mer_iterator</c> resets its window on such a base (N, IUPAC codes, U, gaps).
/// </param>
public readonly record struct KmerCountingOptions(bool Canonical = false, bool AcgtOnly = false)
{
    /// <summary>Literal counting (both flags false); identical to the option-less overloads.</summary>
    public static KmerCountingOptions Default => default;

    /// <summary>True when non-ACGT windows are skipped: <see cref="AcgtOnly"/> or <see cref="Canonical"/>.</summary>
    public bool SkipsNonAcgt => AcgtOnly || Canonical;
}

/// <summary>
/// Word-vector metric for <see cref="KmerAnalyzer.KmerDistance(string, string, int, KmerDistanceMetric)"/>.
/// Each member states the vector it is applied to: raw counts c(w) or relative frequencies f(w) = c(w)/Σc.
/// </summary>
public enum KmerDistanceMetric
{
    /// <summary>√Σ(f₁−f₂)² on relative frequencies — the original <c>KmerDistance</c> (Vinga &amp; Almeida 2003; Zielezinski et al. 2017 Fig. 1).</summary>
    Euclidean = 0,

    /// <summary>Σ(c₁−c₂)² on raw counts — Blaisdell (1986) d_E as reviewed by Vinga &amp; Almeida (2003); alfpy <c>euclid_squared</c> on <c>Counts</c>.</summary>
    SquaredEuclideanCounts,

    /// <summary>Σ|f₁−f₂| (L1, city block) on relative frequencies — scipy <c>cityblock</c>, alfpy <c>manhattan</c> on <c>Freqs</c>.</summary>
    Manhattan,

    /// <summary>max|f₁−f₂| (L∞) on relative frequencies — scipy <c>chebyshev</c>, alfpy <c>chebyshev</c> on <c>Freqs</c>.</summary>
    Chebyshev,

    /// <summary>Σ|f₁−f₂|/(f₁+f₂) on relative frequencies, 0/0 terms omitted — scipy <c>canberra</c>, alfpy <c>canberra</c> on <c>Freqs</c>.</summary>
    Canberra,

    /// <summary>1 − c₁·c₂/(‖c₁‖‖c₂‖) (scale-invariant: counts = frequencies), clipped to [0, 2] — scipy <c>cosine</c>; a zero vector has similarity 0 (distance 1).</summary>
    Cosine,

    /// <summary>D2 = Σ c₁(w)·c₂(w), the inner product of the count vectors — a similarity statistic, not a distance (Torney et al. 1990; Lippert et al. 2005; Reinert et al. 2009).</summary>
    D2,
}
