using System.Runtime.InteropServices;
using System.Text;

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
    /// Serial driver of the single counting loop <see cref="CountWindowRange"/>: validates, upper-cases once and
    /// counts every window [0, L − k] (cancellation + progress at every 1000th window, final 1.0).
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

        var upper = sequence.ToUpperInvariant();
        int total = upper.Length - k + 1;
        var counts = new Dictionary<string, int>();
        CountWindowRange(upper, k, 0, total, acgtOnly, counts, cancellationToken,
            progress is null ? null : i => progress.Report((double)i / total));

        progress?.Report(1.0);
        return counts;
    }

    /// <summary>Cancellation/progress checkpoint interval (windows), shared by the serial and parallel drivers.</summary>
    private const int CountCheckInterval = 1000;

    /// <summary>
    /// The single k-mer counting loop: adds the windows starting at <paramref name="firstWindow"/> ..
    /// <paramref name="endWindow"/> − 1 of the upper-cased <paramref name="upper"/> to <paramref name="counts"/>.
    /// </summary>
    /// <remarks>
    /// A window reads <c>upper[i .. i + k − 1]</c>, so a range needs the k − 1 symbols after its last start —
    /// adjacent ranges of a partition overlap by k − 1 symbols and every window is counted exactly once.
    /// Lookup is by <see cref="ReadOnlySpan{T}"/> through the dictionary's alternate lookup, so a string is
    /// allocated only when a k-mer is seen for the first time (one allocation per distinct k-mer instead of
    /// one per window). At every window index divisible by <see cref="CountCheckInterval"/> the token is polled and
    /// <paramref name="checkpoint"/> receives the index. With <paramref name="acgtOnly"/>, windows containing a
    /// non-ACGT symbol are skipped, tracked in O(1) per window via the last non-ACGT index.
    /// </remarks>
    private static void CountWindowRange(
        string upper,
        int k,
        int firstWindow,
        int endWindow,
        bool acgtOnly,
        Dictionary<string, int> counts,
        CancellationToken cancellationToken,
        Action<int>? checkpoint)
    {
        var seq = upper.AsSpan();
        var lookup = counts.GetAlternateLookup<ReadOnlySpan<char>>();
        int scanned = firstWindow, lastInvalid = firstWindow - 1;

        for (int i = firstWindow; i < endWindow; i++)
        {
            if (i % CountCheckInterval == 0)
            {
                cancellationToken.ThrowIfCancellationRequested();
                checkpoint?.Invoke(i);
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

            CollectionsMarshal.GetValueRefOrAddDefault(lookup, seq.Slice(i, k), out _)++;
        }
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
    /// Counts k-mers in a DNA sequence (<see cref="CountKmers(string, int)"/> on <see cref="DnaSequence.Sequence"/>).
    /// </summary>
    /// <exception cref="ArgumentNullException"><paramref name="dna"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="k"/> ≤ 0 and the sequence is non-empty.</exception>
    public static Dictionary<string, int> CountKmers(DnaSequence dna, int k)
    {
        ArgumentNullException.ThrowIfNull(dna);
        return CountKmers(dna.Sequence, k);
    }

    /// <summary>
    /// Counts k-mers in a DNA sequence with cancellation support.
    /// </summary>
    /// <exception cref="ArgumentNullException"><paramref name="dna"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="k"/> ≤ 0 and the sequence is non-empty.</exception>
    public static Dictionary<string, int> CountKmers(
        DnaSequence dna,
        int k,
        CancellationToken cancellationToken,
        IProgress<double>? progress = null)
    {
        ArgumentNullException.ThrowIfNull(dna);
        return CountKmers(dna.Sequence, k, cancellationToken, progress);
    }

    /// <summary>
    /// Counts the k-mers of a character span (same result and contract as <see cref="CountKmers(string, int)"/>).
    /// </summary>
    /// <remarks>
    /// Runs the class's single counting loop on one upper-cased copy of the span: one O(L) copy plus one string per
    /// distinct k-mer (the loop looks windows up by span). The contract is <see cref="CountKmers(string, int)"/>'s:
    /// empty input yields an empty dictionary for any k, and k ≤ 0 throws only for non-empty input. (Core's
    /// <c>SequenceExtensions.CountKmersSpan</c> throws for k ≤ 0 even on an empty span and allocates a string per
    /// window; this method no longer delegates to it.)
    /// </remarks>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="k"/> ≤ 0 and the span is non-empty.</exception>
    public static Dictionary<string, int> CountKmersSpan(ReadOnlySpan<char> sequence, int k)
        => CountKmersCore(sequence.IsEmpty ? string.Empty : new string(sequence), k, acgtOnly: false, CancellationToken.None, null);

    /// <summary>
    /// Counts k-mers on several cores: the windows are partitioned into contiguous ranges (adjacent ranges share
    /// k − 1 symbols), each range is counted by the class's single counting loop into its own table, and the tables
    /// are merged. The result equals <see cref="CountKmers(string, int, KmerCountingOptions, CancellationToken, IProgress{double}?)"/>
    /// exactly (same keys and counts) for every input and option.
    /// </summary>
    /// <remarks>
    /// <para>Data-parallel partition of an embarrassingly parallel count, as in Jellyfish's multi-threaded counting
    /// (<c>count -t</c>; Marçais &amp; Kingsford 2011) — each window belongs to exactly one range, and counting is a
    /// sum, so the merged table is independent of the partition. Ranges are at least <see cref="ParallelMinWindowsPerRange"/>
    /// windows long; inputs with fewer windows than two ranges are counted serially.</para>
    /// <para>Cancellation: every range polls <paramref name="cancellationToken"/> at each window index divisible by 1000
    /// (as the serial loop) and <see cref="Parallel.For(int, int, ParallelOptions, Action{int})"/> observes it, so a
    /// cancelled call throws <see cref="OperationCanceledException"/> carrying the token. Progress: the fraction of
    /// windows scanned, reported at the same checkpoints, non-decreasing in [0, 1), then exactly one final 1.0; with
    /// several ranges the reports come from worker threads under a lock. The canonical fold (if requested) runs once,
    /// after the merge.</para>
    /// <para>Measured speed-up (Release, 4 cores, random 10 Mbp, 4 ranges): see docs/algorithms/K-mer/Asynchronous_K-mer_Counting.md §5.</para>
    /// </remarks>
    /// <param name="sequence">The sequence (case-insensitive). Null/empty returns an empty dictionary.</param>
    /// <param name="k">The k-mer length. Must be positive for non-empty input.</param>
    /// <param name="options">Counting mode (literal by default).</param>
    /// <param name="maxDegreeOfParallelism">Maximum number of ranges counted concurrently; −1 (default) =
    /// <see cref="Environment.ProcessorCount"/>. 1 counts serially.</param>
    /// <param name="cancellationToken">Cancellation token.</param>
    /// <param name="progress">Optional progress reporter (0.0 to 1.0).</param>
    /// <returns>Dictionary mapping (canonical) k-mers to their counts.</returns>
    /// <exception cref="ArgumentOutOfRangeException">k ≤ 0 with non-empty input, or <paramref name="maxDegreeOfParallelism"/>
    /// is 0 or less than −1.</exception>
    public static Dictionary<string, int> CountKmersParallel(
        string sequence,
        int k,
        KmerCountingOptions options = default,
        int maxDegreeOfParallelism = -1,
        CancellationToken cancellationToken = default,
        IProgress<double>? progress = null)
    {
        if (maxDegreeOfParallelism == 0 || maxDegreeOfParallelism < -1)
            throw new ArgumentOutOfRangeException(nameof(maxDegreeOfParallelism), maxDegreeOfParallelism,
                "Degree of parallelism must be positive or -1 (all processors).");
        ValidateKmerLength(sequence, k);

        int degree = maxDegreeOfParallelism == -1 ? Environment.ProcessorCount : maxDegreeOfParallelism;
        int total = string.IsNullOrEmpty(sequence) ? 0 : Math.Max(0, sequence.Length - k + 1);
        int ranges = Math.Min(degree, total / ParallelMinWindowsPerRange);
        if (ranges < 2)
            return CountKmers(sequence, k, options, cancellationToken, progress);

        cancellationToken.ThrowIfCancellationRequested();
        var upper = sequence.ToUpperInvariant();
        var tables = new Dictionary<string, int>[ranges];
        var gate = new object();
        long scanned = 0;
        double lastReported = -1;

        void Checkpoint(int window)
        {
            // Each range reports once per 1000 windows; the reported fraction is the windows covered so far.
            lock (gate)
            {
                double fraction = (double)scanned / total;
                scanned += CountCheckInterval;
                if (fraction > lastReported && fraction < 1.0)
                {
                    lastReported = fraction;
                    progress!.Report(fraction);
                }
            }
        }

        try
        {
            Parallel.For(0, ranges,
                new ParallelOptions { CancellationToken = cancellationToken, MaxDegreeOfParallelism = degree },
                r =>
                {
                    int first = (int)((long)total * r / ranges);
                    int end = (int)((long)total * (r + 1) / ranges);
                    var table = new Dictionary<string, int>();
                    CountWindowRange(upper, k, first, end, options.SkipsNonAcgt, table, cancellationToken,
                        progress is null ? null : Checkpoint);
                    tables[r] = table;
                });
        }
        catch (AggregateException ex) when (cancellationToken.IsCancellationRequested
                                            && ex.Flatten().InnerExceptions.All(e => e is OperationCanceledException))
        {
            throw new OperationCanceledException(cancellationToken);
        }

        var merged = tables.MaxBy(t => t.Count)!;
        merged.EnsureCapacity(tables.Sum(t => t.Count)); // upper bound of the union: no rehash while merging
        foreach (var table in tables)
        {
            if (ReferenceEquals(table, merged))
                continue;
            foreach (var (kmer, count) in table)
                CollectionsMarshal.GetValueRefOrAddDefault(merged, kmer, out _) += count;
        }

        progress?.Report(1.0);
        return options.Canonical ? FoldToCanonical(merged) : merged;
    }

    /// <summary>Minimum number of windows per range of <see cref="CountKmersParallel"/> (smaller inputs are counted serially).</summary>
    public const int ParallelMinWindowsPerRange = 65536;

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
        => FindMostFrequentKmers(sequence, k, KmerCountingOptions.Default);

    /// <summary>
    /// Finds the most frequent k-mers under explicit <see cref="KmerCountingOptions"/>: literal (default, Rosalind BA1B),
    /// ACGT-only windows, or canonical k-mers (Jellyfish <c>count -C</c>).
    /// </summary>
    /// <remarks>
    /// The arg-max of the table of
    /// <see cref="CountKmers(string, int, KmerCountingOptions, CancellationToken, IProgress{double}?)"/>: every key tied at
    /// the maximum count. With <c>Canonical = true</c> this is the top of <c>jellyfish count -C -m k</c> + <c>jellyfish dump -c</c>
    /// (keys are canonical representatives min(w, RC(w)); a k-mer and its reverse complement are one entry, so their
    /// occurrences on both strands are pooled — the exact-match, reverse-complement-aware frequent-words question; Rosalind
    /// BA1J adds mismatches, which this method does not do). Cross-checked against Jellyfish 2.3.1 (BA1B sample k = 4:
    /// canonical → {ATGC} with 4, literal → {CATG, GCAT} with 3).
    /// </remarks>
    /// <param name="sequence">The sequence (case-insensitive). Null/empty gives an empty result.</param>
    /// <param name="k">The k-mer length. Must be positive for non-empty input.</param>
    /// <param name="options">Counting mode.</param>
    /// <returns>All (canonical) k-mers tied at the maximum count, in table order (order is not part of the contract).</returns>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="k"/> ≤ 0 and the sequence is non-empty.</exception>
    public static IEnumerable<string> FindMostFrequentKmers(string sequence, int k, KmerCountingOptions options)
    {
        var counts = CountKmers(sequence, k, options);

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
        => GetKmerFrequencies(sequence, k, KmerCountingOptions.Default);

    /// <summary>
    /// Gets the k-mer frequency profile under explicit <see cref="KmerCountingOptions"/>: literal (default), ACGT-only
    /// (kPAL / Jellyfish window convention) or canonical (Jellyfish <c>count -C</c>).
    /// </summary>
    /// <remarks>
    /// f(w) = c(w) / Σc over the table of
    /// <see cref="CountKmers(string, int, KmerCountingOptions, CancellationToken, IProgress{double}?)"/>; Σc is the number of
    /// counted windows (L − k + 1 literal; the all-ACGT windows with <c>AcgtOnly</c> or <c>Canonical</c>), so the values
    /// always sum to 1. ACGT-only = a kPAL profile (<c>kpal/klib.py</c> <c>Profile.from_sequences</c> splits the sequence on
    /// every non-ACGT symbol and counts the k-mers of each part) divided by its total; canonical = Jellyfish
    /// <c>count -C</c> + <c>dump -c</c> counts divided by their total. Cross-checked against kPAL (run from source) and
    /// Jellyfish 2.3.1 (docs/algorithms/K-mer/K-mer_Frequency_Analysis.md §7.4).
    /// </remarks>
    /// <param name="sequence">The sequence (case-insensitive). Null/empty returns an empty dictionary.</param>
    /// <param name="k">The k-mer length. Must be positive for non-empty input.</param>
    /// <param name="options">Counting mode.</param>
    /// <returns>Dictionary mapping (canonical) k-mers to their relative frequencies; empty when nothing is counted.</returns>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="k"/> ≤ 0 and the sequence is non-empty.</exception>
    public static Dictionary<string, double> GetKmerFrequencies(string sequence, int k, KmerCountingOptions options)
    {
        var counts = CountKmers(sequence, k, options);
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

        return KmerDistance(seq1, seq2, k, metric, markovOrder: 0);
    }

    /// <summary>
    /// <see cref="KmerDistance(string, string, int, KmerDistanceMetric)"/> with an explicit background Markov order
    /// for <see cref="KmerDistanceMetric.D2Star"/> / <see cref="KmerDistanceMetric.D2Shepherd"/>
    /// (<see cref="BackgroundAdjustedD2(string, string, int, int)"/>); the other metrics have no background model
    /// and require <paramref name="markovOrder"/> = 0.
    /// </summary>
    /// <param name="seq1">First sequence (case-insensitive).</param>
    /// <param name="seq2">Second sequence.</param>
    /// <param name="k">K-mer length; must be positive (≤ <see cref="MaxBackgroundAdjustedK"/> for D2*/D2S).</param>
    /// <param name="metric">The metric.</param>
    /// <param name="markovOrder">Background Markov order r (0 ≤ r &lt; k, or −1 = BIC per sequence) for D2*/D2S; 0 for every other metric.</param>
    /// <returns>The metric value.</returns>
    /// <exception cref="ArgumentOutOfRangeException">k, <paramref name="metric"/> or <paramref name="markovOrder"/> is out of range.</exception>
    /// <exception cref="ArgumentException"><paramref name="markovOrder"/> ≠ 0 for a metric without background model, or (D2*/D2S)
    /// a sequence has no ACGT k-mer (a null sequence counts as empty, as for every other metric).</exception>
    public static double KmerDistance(string seq1, string seq2, int k, KmerDistanceMetric metric, int markovOrder)
    {
        if (k <= 0)
            throw new ArgumentOutOfRangeException(nameof(k), "K must be positive.");

        return KmerDistance(seq1, seq2, k, metric, markovOrder, bothStrands: false);
    }

    /// <summary>
    /// <see cref="KmerDistance(string, string, int, KmerDistanceMetric, int)"/> with the CAFE both-strand mode
    /// (<c>-R</c>) for <see cref="KmerDistanceMetric.D2Star"/> / <see cref="KmerDistanceMetric.D2Shepherd"/>
    /// (<see cref="BackgroundAdjustedD2(string, string, int, int, bool)"/>).
    /// </summary>
    /// <remarks>
    /// A null sequence is treated as the empty sequence for every metric (the word-vector metrics then see the zero
    /// vector; D2*/D2S reject it with <see cref="ArgumentException"/> because an empty sequence has no background).
    /// <paramref name="bothStrands"/> applies only to D2*/D2S and to <see cref="KmerDistanceMetric.SpacedEvolutionary"/>;
    /// for the plain word-vector metrics count canonical k-mers
    /// with <see cref="CountKmers(string, int, KmerCountingOptions, CancellationToken, IProgress{double}?)"/> and use
    /// the count-table overload. <see cref="KmerDistanceMetric.SpacedEvolutionary"/> is the <c>spaced -d EV</c> distance
    /// with the contiguous pattern 1^k: <see cref="SpacedWordDistance(string, string, IReadOnlyList{string}, KmerDistanceMetric, KmerCountingOptions, bool)"/>
    /// with that one pattern (<paramref name="bothStrands"/> = <c>spaced</c>'s reverse-complement mode, seq1 on both strands).
    /// </remarks>
    /// <param name="seq1">First sequence (case-insensitive).</param>
    /// <param name="seq2">Second sequence.</param>
    /// <param name="k">K-mer length; must be positive (≤ <see cref="MaxBackgroundAdjustedK"/> for D2*/D2S).</param>
    /// <param name="metric">The metric.</param>
    /// <param name="markovOrder">Background Markov order (D2*/D2S only; 0 for every other metric).</param>
    /// <param name="bothStrands">CAFE <c>-R</c> both-strand counts and background (D2*/D2S), or <c>spaced</c>'s reverse-complement
    /// mode (<see cref="KmerDistanceMetric.SpacedEvolutionary"/>); false for every other metric.</param>
    /// <returns>The metric value.</returns>
    /// <exception cref="ArgumentOutOfRangeException">k, <paramref name="metric"/> or <paramref name="markovOrder"/> is out of range.</exception>
    /// <exception cref="ArgumentException"><paramref name="markovOrder"/> ≠ 0 or <paramref name="bothStrands"/> set for a metric
    /// without background model, (D2*/D2S) a sequence has no ACGT k-mer, or (SpacedEvolutionary) a sequence has fewer than
    /// k letters.</exception>
    public static double KmerDistance(string seq1, string seq2, int k, KmerDistanceMetric metric, int markovOrder, bool bothStrands)
    {
        if (k <= 0)
            throw new ArgumentOutOfRangeException(nameof(k), "K must be positive.");

        seq1 ??= string.Empty;
        seq2 ??= string.Empty;
        return metric switch
        {
            KmerDistanceMetric.D2Star => BackgroundAdjustedD2(seq1, seq2, k, markovOrder, bothStrands).D2StarDistance,
            KmerDistanceMetric.D2Shepherd => BackgroundAdjustedD2(seq1, seq2, k, markovOrder, bothStrands).D2ShepherdDistance,
            _ when markovOrder != 0 => throw new ArgumentException(
                "markovOrder applies only to the background-adjusted metrics D2Star and D2Shepherd.", nameof(markovOrder)),
            KmerDistanceMetric.SpacedEvolutionary => SpacedWordDistance(
                seq1, seq2, [new string('1', k)], metric, new KmerCountingOptions(AcgtOnly: true), bothStrands),
            _ when bothStrands => throw new ArgumentException(
                "bothStrands applies only to the background-adjusted metrics D2Star and D2Shepherd; count canonical k-mers (KmerCountingOptions.Canonical) for the other metrics.",
                nameof(bothStrands)),
            _ => KmerDistance(CountKmers(seq1, k), CountKmers(seq2, k), metric),
        };
    }

    /// <summary>
    /// Parses a metric name as used by the MCP tools (case-insensitive, surrounding blanks ignored): <c>euclidean</c>
    /// (also null/empty), <c>squared_euclidean_counts</c>, <c>manhattan</c>, <c>chebyshev</c>, <c>canberra</c>,
    /// <c>cosine</c>, <c>d2</c>, <c>d2star</c>, <c>d2shepherd</c> (alias <c>d2s</c>), <c>jensen_shannon</c> (alias <c>js</c>),
    /// <c>euclidean_counts</c>, <c>ev</c> (alias <c>evolutionary</c>; <see cref="KmerDistanceMetric.SpacedEvolutionary"/>).
    /// </summary>
    /// <exception cref="ArgumentException">The name is not one of the above.</exception>
    public static KmerDistanceMetric ParseDistanceMetric(string? name) =>
        (name ?? string.Empty).Trim().ToLowerInvariant() switch
        {
            "" or "euclidean" => KmerDistanceMetric.Euclidean,
            "squared_euclidean_counts" => KmerDistanceMetric.SquaredEuclideanCounts,
            "manhattan" => KmerDistanceMetric.Manhattan,
            "chebyshev" => KmerDistanceMetric.Chebyshev,
            "canberra" => KmerDistanceMetric.Canberra,
            "cosine" => KmerDistanceMetric.Cosine,
            "d2" => KmerDistanceMetric.D2,
            "d2star" => KmerDistanceMetric.D2Star,
            "d2shepherd" or "d2s" => KmerDistanceMetric.D2Shepherd,
            "jensen_shannon" or "js" => KmerDistanceMetric.JensenShannon,
            "euclidean_counts" => KmerDistanceMetric.EuclideanCounts,
            "ev" or "evolutionary" => KmerDistanceMetric.SpacedEvolutionary,
            _ => throw new ArgumentException(
                "metric must be one of: euclidean, squared_euclidean_counts, manhattan, chebyshev, canberra, cosine, d2, d2star, d2shepherd, jensen_shannon, euclidean_counts, ev",
                nameof(name)),
        };

    /// <summary>Largest k for <see cref="BackgroundAdjustedD2"/>: the statistics sum over all 4^k DNA words (4^12 ≈ 1.7·10⁷).</summary>
    public const int MaxBackgroundAdjustedK = 12;

    /// <summary>
    /// Background-adjusted word-match statistics D2* and D2S (Reinert, Chew, Sun &amp; Waterman 2009, J Comput Biol
    /// 16:1615; Wan, Reinert, Sun &amp; Waterman 2010) and their dissimilarities d2* and d2S (Song et al. 2014,
    /// Brief Bioinform 15:343), with a Markov background of order <paramref name="markovOrder"/> estimated from each
    /// sequence.
    /// </summary>
    /// <remarks>
    /// <para>Counts X_w, Y_w are the single-strand k-mer counts over the ACGT windows (upper-cased; a window with
    /// another symbol is skipped, the Jellyfish convention used by CAFE), n̄ = Σ X_w, m̄ = Σ Y_w. The expected count is
    /// E_X(w) = n̄·p̂_X(w) under the order-r Markov chain fitted to that sequence by maximum likelihood:
    /// p̂(w) = N(w₁..w_r)/Σ N(r-mers) · Π_{i&gt;r} N(w_{i−r}..w_i)/Σ_a N(w_{i−r}..w_{i−1}a), with the r-mer and
    /// (r+1)-mer counts N taken over the sequence's ACGT windows (r = 0: p̂(w) = Π p̂(w_i), the letter frequencies).
    /// Centred counts X̃ = X − E_X, Ỹ = Y − E_Y; the sums run over all 4^k words w ∈ {A,C,G,T}^k (words absent from
    /// both sequences contribute too):</para>
    /// <para>D2* = Σ X̃Ỹ/√(E_X E_Y) (words with E_X·E_Y = 0 omitted); d2* = ½(1 − D2*/√(Σ X̃²/E_X · Σ Ỹ²/E_Y)).</para>
    /// <para>D2S = Σ X̃Ỹ/√(X̃² + Ỹ²) (words with X̃ = Ỹ = 0 omitted); d2S = ½(1 − D2S/√(Σ X̃²/√(X̃²+Ỹ²) · Σ Ỹ²/√(X̃²+Ỹ²))).</para>
    /// <para>These are the formulas of CAFE (Lu et al. 2017, <c>dist_model.cpp</c> D2starStrategy/D2sheppStrategy,
    /// single-strand mode). CAFE derives the r- and (r+1)-mer counts by marginalising the k-mer table on the prefix;
    /// this method counts them on the sequence (the textbook maximum-likelihood estimator). A Python replica of the
    /// formulas reproduces the CAFE binary to 6 digits when given CAFE's estimator, and this method to 1e-12 with the
    /// sequence estimator (docs/algorithms/K-mer/K-mer_Euclidean_Distance.md §7.4). A distance is NaN when a
    /// normaliser is 0 (e.g. a sequence whose k-mer counts equal their expectation exactly). The two
    /// dissimilarities are clamped to [0, 1]: the Cauchy–Schwarz inequality bounds them there, and rounding could otherwise
    /// return −1.1e-16 for identical sequences (the raw D2*/D2S are not clamped).</para>
    /// </remarks>
    /// <param name="seq1">First sequence (case-insensitive).</param>
    /// <param name="seq2">Second sequence (case-insensitive).</param>
    /// <param name="k">Word length, 1 ≤ k ≤ <see cref="MaxBackgroundAdjustedK"/>.</param>
    /// <param name="markovOrder">Background Markov order r, 0 ≤ r &lt; k (0 = i.i.d. letters, the default), or
    /// <see cref="AutoMarkovOrder"/> (−1) to choose each sequence's order in [0, min(k − 1, 10)] by
    /// <see cref="SelectMarkovOrder"/> (BIC; CAFE <c>-M -1</c>).</param>
    /// <returns>The two statistics, the two dissimilarities and the orders used.</returns>
    /// <exception cref="ArgumentNullException">A sequence is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException">k or <paramref name="markovOrder"/> is out of range.</exception>
    /// <exception cref="ArgumentException">A sequence has no ACGT k-mer window.</exception>
    public static D2StarStatistics BackgroundAdjustedD2(string seq1, string seq2, int k, int markovOrder = 0)
        => BackgroundAdjustedD2(seq1, seq2, k, markovOrder, bothStrands: false);

    /// <summary>
    /// <see cref="BackgroundAdjustedD2(string, string, int, int)"/> with an optional both-strand mode, the semantics of
    /// CAFE's <c>-R</c> option (Lu et al. 2017, <c>kmer.cpp</c> <c>KmerModel::load</c> and
    /// <c>KmerProbEnsembDelegate::getKmerlogProb</c>).
    /// </summary>
    /// <remarks>
    /// <para>With <paramref name="bothStrands"/> = false this is exactly the single-strand statistic.</para>
    /// <para>With <paramref name="bothStrands"/> = true each word's count is X^R(w) = X(w) + X(RC(w)) (so the total is
    /// 2n̄, palindromes w = RC(w) count twice), and the background probability is symmetrised,
    /// p^R(w) = ½(p̂(w) + p̂(RC(w))), where p̂ is the order-r Markov chain fitted to the given strand. The expected count is
    /// therefore E^R(w) = 2n̄·p^R(w) = E_X(w) + E_X(RC(w)), the expectation of X(w) + X(RC(w)). The sums still run over
    /// all 4^k words (each non-palindromic pair {w, RC(w)} contributes twice, as in CAFE). The Markov order chosen by
    /// BIC (−1) is computed on the given strand, as CAFE's <c>getEstMarkovOrder</c>.</para>
    /// <para>CAFE takes p̂ from prefix marginals of the k-mer table (and prunes probability-1 factors). This method uses
    /// the maximum-likelihood chain on the sequence, as in the single-strand mode. A Python replica reproduces the CAFE
    /// <c>-R</c> binary to its 6 printed digits on 20 runs when given CAFE's estimator, and this method to 1e-12 with
    /// the sequence estimator (docs/algorithms/K-mer/K-mer_Euclidean_Distance.md §7.5).</para>
    /// <para>Memory: the counts and the Markov tables are sparse (one entry per k-mer / r-mer / (r+1)-mer present in
    /// the sequence; dense arrays only for r + 1 ≤ 8), so high orders such as k = 12, r = 11 need no 4^12 array. The
    /// time is Θ(4^k·k) for every order.</para>
    /// </remarks>
    /// <param name="seq1">First sequence (case-insensitive).</param>
    /// <param name="seq2">Second sequence (case-insensitive).</param>
    /// <param name="k">Word length, 1 ≤ k ≤ <see cref="MaxBackgroundAdjustedK"/>.</param>
    /// <param name="markovOrder">Background Markov order r, 0 ≤ r &lt; k, or <see cref="AutoMarkovOrder"/> (−1, BIC).</param>
    /// <param name="bothStrands">Combine each word with its reverse complement (CAFE <c>-R</c>).</param>
    /// <returns>The two statistics, the two dissimilarities and the orders used.</returns>
    /// <exception cref="ArgumentNullException">A sequence is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException">k or <paramref name="markovOrder"/> is out of range.</exception>
    /// <exception cref="ArgumentException">A sequence has no ACGT k-mer window.</exception>
    public static D2StarStatistics BackgroundAdjustedD2(string seq1, string seq2, int k, int markovOrder, bool bothStrands)
    {
        ArgumentNullException.ThrowIfNull(seq1);
        ArgumentNullException.ThrowIfNull(seq2);
        if (k <= 0 || k > MaxBackgroundAdjustedK)
            throw new ArgumentOutOfRangeException(nameof(k), k, $"K must be in [1, {MaxBackgroundAdjustedK}].");
        if (markovOrder < AutoMarkovOrder || markovOrder >= k)
            throw new ArgumentOutOfRangeException(nameof(markovOrder), markovOrder, "Markov order must be in [0, k) or -1 (BIC).");

        int maxAutoOrder = Math.Min(k - 1, MaxAutoMarkovOrder);
        int order1 = markovOrder == AutoMarkovOrder ? SelectMarkovOrder(seq1, maxAutoOrder) : markovOrder;
        int order2 = markovOrder == AutoMarkovOrder ? SelectMarkovOrder(seq2, maxAutoOrder) : markovOrder;
        var x = WordBackground.Fit(seq1, k, order1, nameof(seq1));
        var y = WordBackground.Fit(seq2, k, order2, nameof(seq2));

        double starNum = 0, starX = 0, starY = 0, shepNum = 0, shepX = 0, shepY = 0;
        var word = new int[k];
        var probX = new double[k + 1];
        var probY = new double[k + 1];
        var probXRc = new double[k + 1];
        var probYRc = new double[k + 1];
        probX[0] = probY[0] = probXRc[0] = probYRc[0] = 1.0;
        int depth = 0;
        word[0] = -1;

        // Odometer over {A,C,G,T}^k in lexicographic order; prob[d] = background probability of the d-symbol prefix
        // (probRc[d]: the factors of RC(word) fixed by the first d symbols, complete at d = k).
        while (depth >= 0)
        {
            if (++word[depth] == 4)
            {
                depth--;
                continue;
            }

            probX[depth + 1] = probX[depth] * x.Factor(word, depth, order1);
            probY[depth + 1] = probY[depth] * y.Factor(word, depth, order2);
            if (bothStrands)
            {
                probXRc[depth + 1] = probXRc[depth] * x.ReverseComplementFactor(word, depth, order1);
                probYRc[depth + 1] = probYRc[depth] * y.ReverseComplementFactor(word, depth, order2);
            }

            if (depth < k - 1)
            {
                word[++depth] = -1;
                continue;
            }

            long code = WordBackground.Encode(word, 0, k);
            double cx = x.CountOf(code), cy = y.CountOf(code);
            double px = probX[k], py = probY[k];
            if (bothStrands)
            {
                long rc = WordBackground.EncodeReverseComplement(word);
                cx += x.CountOf(rc);
                cy += y.CountOf(rc);
                px += probXRc[k];
                py += probYRc[k];
            }

            // Single strand: E = n̄·p. Both strands: E = 2n̄·½(p + p_RC) = n̄·(p + p_RC).
            double ex = x.Windows * px, ey = y.Windows * py;
            double xt = cx - ex, yt = cy - ey;
            if (ex > 0 && ey > 0)
            {
                starNum += xt * yt / Math.Sqrt(ex * ey);
                starX += xt * xt / ex;
                starY += yt * yt / ey;
            }

            double norm = Math.Sqrt(xt * xt + yt * yt);
            if (norm > 0)
            {
                shepNum += xt * yt / norm;
                shepX += xt * xt / norm;
                shepY += yt * yt / norm;
            }
        }

        return new D2StarStatistics(
            starNum,
            shepNum,
            ClampUnit(0.5 * (1.0 - starNum / (Math.Sqrt(starX) * Math.Sqrt(starY)))),
            ClampUnit(0.5 * (1.0 - shepNum / (Math.Sqrt(shepX) * Math.Sqrt(shepY)))))
        {
            MarkovOrder1 = order1,
            MarkovOrder2 = order2,
        };
    }

    /// <summary>
    /// Clamps a dissimilarity ½(1 − ratio) to [0, 1]. By the Cauchy–Schwarz inequality |ratio| ≤ 1, so values outside [0, 1]
    /// are only floating-point rounding (e.g. −1.1e-16 for identical sequences). NaN (a zero normaliser) passes through.
    /// </summary>
    private static double ClampUnit(double distance) => Math.Clamp(distance, 0.0, 1.0);

    /// <summary><c>markovOrder</c> value of <see cref="BackgroundAdjustedD2"/> that selects each sequence's order by BIC.</summary>
    public const int AutoMarkovOrder = -1;

    /// <summary>Highest order tried by the BIC selection (CAFE <c>MAX_ORDER</c> = 10).</summary>
    public const int MaxAutoMarkovOrder = 10;

    /// <summary>
    /// Bayesian information criterion of the order-r Markov chain fitted to a sequence:
    /// BIC(r) = −2·ln L̂_r + (|A| − 1)·|A|^r·ln N_r, with |A| = 4, ln L̂_r = Σ N(u a)·ln(N(u a)/Σ_b N(u b)) over the
    /// sequence's ACGT (r+1)-mers u a, and N_r their number (Schwarz 1978; Katz 1981 Markov-order estimation; the
    /// criterion CAFE <c>-M -1</c> minimises).
    /// </summary>
    /// <param name="sequence">The sequence (case-insensitive; windows with a non-ACGT symbol are skipped).</param>
    /// <param name="order">Markov order r ≥ 0.</param>
    /// <returns>The BIC value (smaller is better).</returns>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="order"/> is negative.</exception>
    /// <exception cref="ArgumentException">The sequence has no ACGT (r+1)-mer.</exception>
    public static double MarkovOrderBic(string sequence, int order)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentOutOfRangeException.ThrowIfNegative(order);

        var words = CountKmers(sequence, order + 1, new KmerCountingOptions(AcgtOnly: true));
        if (words.Count == 0)
            throw new ArgumentException("Sequence has no (order+1)-mer window over A/C/G/T.", nameof(sequence));

        var contextTotals = new Dictionary<string, double>(StringComparer.Ordinal);
        double observations = 0;
        foreach (var (word, count) in words)
        {
            var context = word[..order];
            contextTotals[context] = contextTotals.GetValueOrDefault(context) + count;
            observations += count;
        }

        double logLikelihood = 0;
        foreach (var (word, count) in words)
            logLikelihood += count * Math.Log(count / contextTotals[word[..order]]);

        return -2.0 * logLikelihood + 3.0 * Math.Pow(4, order) * Math.Log(observations);
    }

    /// <summary>
    /// The Markov order r ∈ [0, <paramref name="maxOrder"/>] with the smallest <see cref="MarkovOrderBic"/>
    /// (ties → the smaller order).
    /// </summary>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="maxOrder"/> is negative.</exception>
    /// <exception cref="ArgumentException">The sequence has no ACGT (maxOrder+1)-mer.</exception>
    public static int SelectMarkovOrder(string sequence, int maxOrder)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(maxOrder);
        int best = 0;
        double bestBic = double.PositiveInfinity;
        for (int r = 0; r <= maxOrder; r++)
        {
            double bic = MarkovOrderBic(sequence, r);
            if (bic < bestBic)
            {
                bestBic = bic;
                best = r;
            }
        }
        return best;
    }

    /// <summary>Per-sequence k-mer counts (2-bit codes) and fitted order-r Markov background for <see cref="BackgroundAdjustedD2"/>.</summary>
    private sealed class WordBackground
    {
        // Largest (r+1)-mer table kept as a dense array (4^8 doubles = 512 KB); higher orders use sparse tables
        // keyed by the 2-bit code, holding only the r-mers / (r+1)-mers that occur (absent = probability 0).
        private const int MaxDenseTableWordLength = 8;

        private readonly Dictionary<long, int> _counts = new();
        private ProbabilityTable _initial = null!;
        private ProbabilityTable _transition = null!;

        public double Windows { get; private set; }

        public static WordBackground Fit(string sequence, int k, int order, string paramName)
        {
            var acgt = new KmerCountingOptions(AcgtOnly: true);
            var kmers = CountKmers(sequence, k, acgt);
            if (kmers.Count == 0)
                throw new ArgumentException("Sequence has no k-mer window over A/C/G/T.", paramName);

            var model = new WordBackground();
            foreach (var (kmer, count) in kmers)
            {
                model._counts[Encode(kmer)] = count;
                model.Windows += count;
            }

            // Initial distribution over r-mers and transition probabilities P(a | r-mer context), both maximum
            // likelihood from the sequence's ACGT r-mer and (r+1)-mer counts.
            bool dense = order + 1 <= MaxDenseTableWordLength;
            model._initial = new ProbabilityTable(dense ? 1 << (2 * order) : 0);
            if (order == 0)
            {
                model._initial.Set(0, 1.0);
            }
            else
            {
                var rmers = CountKmers(sequence, order, acgt);
                double totalR = rmers.Values.Sum();
                foreach (var (rmer, count) in rmers)
                    model._initial.Set(Encode(rmer), count / totalR);
            }

            var words = CountKmers(sequence, order + 1, acgt);
            var rowTotals = new Dictionary<long, double>();
            foreach (var (w, count) in words)
            {
                long context = Encode(w) >> 2;
                rowTotals[context] = rowTotals.GetValueOrDefault(context) + count;
            }

            model._transition = new ProbabilityTable(dense ? 1 << (2 * (order + 1)) : 0);
            foreach (var (w, count) in words)
            {
                long code = Encode(w);
                model._transition.Set(code, count / rowTotals[code >> 2]);
            }

            return model;
        }

        public int CountOf(long code) => _counts.TryGetValue(code, out var c) ? c : 0;

        /// <summary>Probability factor contributed by symbol <paramref name="word"/>[<paramref name="depth"/>].</summary>
        public double Factor(int[] word, int depth, int order)
        {
            if (depth < order - 1)
                return 1.0;
            if (depth == order - 1)
                return _initial.Get(Encode(word, 0, order));
            return _transition.Get(Encode(word, depth - order, order + 1));
        }

        /// <summary>
        /// Factor of p̂(RC(word)) fixed once <paramref name="word"/>[0..<paramref name="depth"/>] is known: the window
        /// word[depth−r..depth] is, reverse-complemented, the context c(w_depth)..c(w_{depth−r+1}) followed by
        /// c(w_{depth−r}) (CAFE's reverse-complement delegate); at depth k − 1 the initial r-mer
        /// c(w_{k−1})..c(w_{k−r}) of RC(word) is multiplied in.
        /// </summary>
        public double ReverseComplementFactor(int[] word, int depth, int order)
        {
            double factor = 1.0;
            if (depth >= order)
                factor = _transition.Get(EncodeComplementDescending(word, depth, order + 1));
            if (depth == word.Length - 1 && order > 0)
                factor *= _initial.Get(EncodeComplementDescending(word, depth, order));
            return factor;
        }

        public static long Encode(int[] word, int start, int length)
        {
            long code = 0;
            for (int i = start; i < start + length; i++)
                code = (code << 2) | (uint)word[i];
            return code;
        }

        /// <summary>Code of the reverse complement of the whole word.</summary>
        public static long EncodeReverseComplement(int[] word) => EncodeComplementDescending(word, word.Length - 1, word.Length);

        // Code of c(word[from]), c(word[from − 1]), …, c(word[from − length + 1]) with c(x) = 3 − x (A↔T, C↔G).
        private static long EncodeComplementDescending(int[] word, int from, int length)
        {
            long code = 0;
            for (int i = from; i > from - length; i--)
                code = (code << 2) | (uint)(3 - word[i]);
            return code;
        }

        private static long Encode(string kmer)
        {
            long code = 0;
            foreach (char c in kmer)
                code = (code << 2) | (uint)(c switch { 'A' => 0, 'C' => 1, 'G' => 2, _ => 3 });
            return code;
        }

        /// <summary>Probabilities indexed by 2-bit word code: a dense array, or (size 0) a sparse dictionary.</summary>
        private sealed class ProbabilityTable(int denseSize)
        {
            private readonly double[]? _dense = denseSize > 0 ? new double[denseSize] : null;
            private readonly Dictionary<long, double>? _sparse = denseSize > 0 ? null : new Dictionary<long, double>();

            public void Set(long code, double value)
            {
                if (_dense is not null)
                    _dense[code] = value;
                else
                    _sparse![code] = value;
            }

            public double Get(long code)
            {
                if (_dense is not null)
                    return _dense[code];
                return _sparse!.TryGetValue(code, out var v) ? v : 0.0;
            }
        }
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
    /// et al. 2005; Reinert et al. 2009); Jensen–Shannon divergence ½Σ f₁ log₂(f₁/m) + ½Σ f₂ log₂(f₂/m), m = ½(f₁+f₂)
    /// (Lin 1991; the JS measure of Leimeister et al. 2014 and of the <c>spaced</c> program, = scipy
    /// <c>jensenshannon(p, q, base=2)²</c>); Euclidean on counts √Σ(c₁−c₂)² (the <c>spaced</c> 1.2 <c>-d EU</c> value).
    /// Cross-checked against scipy.spatial.distance and alfpy 1.0.6
    /// (docs/algorithms/K-mer/K-mer_Euclidean_Distance.md §7.2, §7.5).</para>
    /// </remarks>
    /// <param name="counts1">First count table (non-negative counts).</param>
    /// <param name="counts2">Second count table (non-negative counts).</param>
    /// <param name="metric">The word-vector metric.</param>
    /// <returns>The metric value.</returns>
    /// <exception cref="ArgumentNullException">A table is null.</exception>
    /// <exception cref="ArgumentException">A table contains a negative count, or <paramref name="metric"/> is D2*/D2S or
    /// <see cref="KmerDistanceMetric.SpacedEvolutionary"/> (they need the sequences).</exception>
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
        if (metric is KmerDistanceMetric.D2Star or KmerDistanceMetric.D2Shepherd)
            throw new ArgumentException(
                "D2*/D2S need each sequence's background model; use KmerDistance(string, string, int, metric) or BackgroundAdjustedD2.",
                nameof(metric));
        if (metric is KmerDistanceMetric.SpacedEvolutionary)
            throw new ArgumentException(
                "The spaced EV distance needs the sequence lengths and base composition; use SpacedWordDistance or KmerDistance(string, string, int, metric).",
                nameof(metric));

        double total1 = SumNonNegative(counts1, nameof(counts1));
        double total2 = SumNonNegative(counts2, nameof(counts2));
        return WordVectorDistance(counts1, total1, counts2, total2, metric);
    }

    /// <summary>
    /// The single word-vector metric loop: frequencies are c / <paramref name="total1"/> and c / <paramref name="total2"/>
    /// (0 when a total is 0). The count-table overload passes Σc; <see cref="SpacedWordDistance(string, string, IReadOnlyList{string}, KmerDistanceMetric, KmerCountingOptions, bool)"/>
    /// passes the <c>spaced</c> window totals.
    /// </summary>
    private static double WordVectorDistance(
        IReadOnlyDictionary<string, int> counts1,
        double total1,
        IReadOnlyDictionary<string, int> counts2,
        double total2,
        KmerDistanceMetric metric)
    {
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
            KmerDistanceMetric.Euclidean or KmerDistanceMetric.EuclideanCounts => Math.Sqrt(acc),
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
                case KmerDistanceMetric.EuclideanCounts:
                    acc += ((double)c1 - c2) * ((double)c1 - c2);
                    break;
                case KmerDistanceMetric.JensenShannon:
                    double m = 0.5 * (f1 + f2);
                    if (f1 > 0)
                        acc += 0.5 * f1 * Math.Log2(f1 / m);
                    if (f2 > 0)
                        acc += 0.5 * f2 * Math.Log2(f2 / m);
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
        var (shared, countA, countB) = SharedAndSetSizes(a, b, k, options);
        return (shared, countA + countB - shared);
    }

    /// <summary>|K(A) ∩ K(B)|, |K(A)| and |K(B)| of the distinct k-mer sets under <paramref name="options"/> (null = empty).</summary>
    private static (int Shared, int CountA, int CountB) SharedAndSetSizes(string a, string b, int k, KmerCountingOptions options)
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
        return (shared, setA.Count, setB.Count);
    }

    /// <summary>
    /// Exact containment index C(A, B) = |K(A) ∩ K(B)| / |K(A)|: the fraction of the distinct k-mers of
    /// <paramref name="a"/> that also occur in <paramref name="b"/> (Koslicki &amp; Zabeti 2019, Appl Math Comput
    /// 354:206, "Improving MinHash via the containment index"; sourmash <c>compare --containment</c>, which reports
    /// C(row, column); exact for <c>scaled=1</c> sketches). Asymmetric: C(A, B)·|K(A)| = C(B, A)·|K(B)| = |K(A) ∩ K(B)|.
    /// </summary>
    /// <remarks>
    /// Pass <c>new KmerCountingOptions(Canonical: true)</c> for sourmash's DNA k-mers (strand-collapsed, non-ACGT windows
    /// skipped). K(A) empty (null/empty or shorter than <paramref name="k"/>): the ratio is 0/0 and 0 is returned,
    /// as sourmash's <c>contained_by</c> returns 0 for an empty sketch and as <see cref="JaccardSimilarity(string, string, int)"/> does.
    /// </remarks>
    /// <param name="a">The sequence whose k-mers are tested for containment; null is treated as empty.</param>
    /// <param name="b">The containing sequence; null is treated as empty.</param>
    /// <param name="k">K-mer length; must be positive.</param>
    /// <param name="options">Counting mode for both sequences.</param>
    /// <returns>C(A, B) in [0, 1]; 1 when every k-mer of <paramref name="a"/> occurs in <paramref name="b"/>.</returns>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="k"/> is not positive.</exception>
    public static double ContainmentIndex(string a, string b, int k, KmerCountingOptions options)
    {
        var (shared, countA, _) = SharedAndSetSizes(a, b, k, options);
        return countA == 0 ? 0.0 : (double)shared / countA;
    }

    /// <summary>Exact containment index C(A, B) with literal counting; see <see cref="ContainmentIndex(string, string, int, KmerCountingOptions)"/>.</summary>
    /// <param name="a">The sequence whose k-mers are tested for containment; null is treated as empty.</param>
    /// <param name="b">The containing sequence; null is treated as empty.</param>
    /// <param name="k">K-mer length; must be positive.</param>
    /// <returns>C(A, B) in [0, 1].</returns>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="k"/> is not positive.</exception>
    public static double ContainmentIndex(string a, string b, int k)
        => ContainmentIndex(a, b, k, default);

    #region MinHash (Mash) sketches

    /// <summary>Default Mash sketch size (<c>mash sketch -s</c>, 1000).</summary>
    public const int DefaultMashSketchSize = 1000;

    /// <summary>Default Mash hash seed (<c>mash sketch -S</c>, 42).</summary>
    public const uint DefaultMashSeed = 42;

    /// <summary>Largest k-mer size Mash accepts (<c>mash sketch -k</c>, 1..32).</summary>
    public const int MaxMashKmerSize = 32;

    /// <summary>
    /// Builds a bottom-s MinHash sketch of one sequence exactly as <c>mash sketch</c> 2.3 does (Ondov et al. 2016,
    /// Genome Biol 17:132): see <see cref="CreateMinHashSketch(IEnumerable{string}, int, int, bool, uint)"/>.
    /// </summary>
    /// <param name="sequence">The sequence; null is treated as empty.</param>
    /// <param name="k">K-mer size, 1..32 (Mash default 21).</param>
    /// <param name="sketchSize">Sketch size s ≥ 1 (Mash default 1000).</param>
    /// <param name="canonical">True (default): canonical k-mers min(w, RC(w)); false: forward k-mers (<c>mash sketch -n</c>).</param>
    /// <param name="seed">MurmurHash3 seed (Mash default 42).</param>
    /// <returns>The sketch.</returns>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="k"/> outside 1..32 or <paramref name="sketchSize"/> &lt; 1.</exception>
    public static MinHashSketch CreateMinHashSketch(
        string sequence, int k, int sketchSize = DefaultMashSketchSize, bool canonical = true, uint seed = DefaultMashSeed)
        => CreateMinHashSketch([sequence ?? string.Empty], k, sketchSize, canonical, seed);

    /// <summary>
    /// Builds a bottom-s MinHash sketch of a set of records (e.g. the contigs of one FASTA file) exactly as
    /// <c>mash sketch</c> 2.3 does (Ondov et al. 2016; Mash <c>Sketch.cpp</c> <c>sketchFile</c>/<c>addMinHashes</c>,
    /// <c>hash.cpp</c>, <c>MinHashHeap.cpp</c>).
    /// </summary>
    /// <remarks>
    /// <para>K-mers: upper-cased, windows containing a non-ACGT symbol skipped (Mash nucleotide alphabet), each k-mer
    /// replaced by the lexicographically smaller of itself and its reverse complement (<c>memcmp(fwd, rev) &lt;= 0</c>)
    /// unless <paramref name="canonical"/> is false — the k-mer set of
    /// <see cref="DistinctKmers(string, int, KmerCountingOptions)"/> with <c>Canonical</c> (or <c>AcgtOnly</c>) set. No
    /// k-mer spans two records.</para>
    /// <para>Hash: MurmurHash3_x64_128 of the k ASCII bytes with <paramref name="seed"/>; Mash keeps the first 64-bit
    /// word h1 when 4^k &gt; 2^32 (k ≥ 17, <c>use64</c>) and otherwise its low 32 bits (<c>hash32</c>), compared unsigned.
    /// The sketch is the <paramref name="sketchSize"/> smallest distinct hash values, ascending (<c>MinHashHeap</c> keeps
    /// non-redundant hashes).</para>
    /// <para>Length: Σ record lengths over records of length ≥ k (all symbols, N included); shorter records are skipped
    /// like Mash's <c>l &lt; kmerSize</c> rule. The length enters the p-value of <see cref="CompareMinHashSketches"/>.</para>
    /// </remarks>
    /// <param name="records">The records; null entries are treated as empty.</param>
    /// <param name="k">K-mer size, 1..32 (Mash default 21).</param>
    /// <param name="sketchSize">Sketch size s ≥ 1 (Mash default 1000).</param>
    /// <param name="canonical">True (default): canonical k-mers; false: forward k-mers (<c>mash sketch -n</c>).</param>
    /// <param name="seed">MurmurHash3 seed (Mash default 42).</param>
    /// <returns>The sketch.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="records"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="k"/> outside 1..32 or <paramref name="sketchSize"/> &lt; 1.</exception>
    public static MinHashSketch CreateMinHashSketch(
        IEnumerable<string> records, int k, int sketchSize = DefaultMashSketchSize, bool canonical = true, uint seed = DefaultMashSeed)
    {
        ArgumentNullException.ThrowIfNull(records);
        if (k < 1 || k > MaxMashKmerSize)
            throw new ArgumentOutOfRangeException(nameof(k), k, $"K must be in 1..{MaxMashKmerSize} (mash sketch -k).");
        ArgumentOutOfRangeException.ThrowIfLessThan(sketchSize, 1);

        bool use64 = k > 16; // Mash: use64 = 4^k > 2^32
        var options = canonical ? new KmerCountingOptions(Canonical: true) : new KmerCountingOptions(AcgtOnly: true);
        var kmers = new HashSet<string>(StringComparer.Ordinal);
        long length = 0;
        foreach (var record in records)
        {
            var seq = record ?? string.Empty;
            if (seq.Length < k)
                continue;
            length += seq.Length;
            kmers.UnionWith(CountKmers(seq, k, options).Keys);
        }

        var hashes = new ulong[kmers.Count];
        Span<byte> buffer = stackalloc byte[MaxMashKmerSize];
        int n = 0;
        foreach (var kmer in kmers)
        {
            var bytes = buffer[..k];
            for (int i = 0; i < k; i++)
                bytes[i] = (byte)kmer[i];
            ulong h1 = MurmurHash3X64_128(bytes, seed).H1;
            hashes[n++] = use64 ? h1 : (uint)h1;
        }

        Array.Sort(hashes);
        var bottom = new List<ulong>(Math.Min(sketchSize, hashes.Length));
        foreach (var h in hashes)
        {
            if (bottom.Count == sketchSize)
                break;
            if (bottom.Count == 0 || bottom[^1] != h)
                bottom.Add(h);
        }

        return new MinHashSketch(k, sketchSize, canonical, seed, use64, length, bottom.ToArray());
    }

    /// <summary>
    /// Compares two MinHash sketches as <c>mash dist</c> 2.3 does (Mash <c>CommandDistance.cpp</c>
    /// <c>compareSketches</c> + <c>pValue</c>; Ondov et al. 2016 eqs. 1, 4 and the p-value of the Methods).
    /// </summary>
    /// <remarks>
    /// <para>The sorted sketches are merged until s = min(s_A, s_B) union hashes have been visited (or a sketch is exhausted,
    /// in which case the rest of the other is added and the total capped at s); x = shared hashes among them, and
    /// J = x / denominator — Mash's "x/s" column. Distance: x = denominator → 0; x = 0 → 1; otherwise
    /// −ln(2J/(1 + J))/k capped at 1 (<see cref="MashDistanceFromJaccard"/>).</para>
    /// <para>p-value: P(X ≥ x) for X ~ Binomial(denominator, r) with r = p_A·p_B/(p_A + p_B − p_A·p_B) and
    /// p = 1/(1 + 4^k/length) — Mash's <c>gsl_cdf_binomial_Q(x − 1, r, denominator)</c>, evaluated by
    /// <see cref="StatisticsHelper.BinomialUpperTail"/>; x = 0 → 1.</para>
    /// <para>Both sketches empty: denominator 0, Jaccard 0 (the convention of <see cref="JaccardSimilarity(string, string, int)"/>;
    /// Mash would print nan), distance 0, p-value 1.</para>
    /// </remarks>
    /// <param name="reference">Reference sketch (Mash's first argument).</param>
    /// <param name="query">Query sketch.</param>
    /// <returns>Shared hashes, denominator, Jaccard estimate, Mash distance and p-value.</returns>
    /// <exception cref="ArgumentNullException">A sketch is null.</exception>
    /// <exception cref="ArgumentException">The sketches differ in k, seed or canonical mode (Mash refuses to compare them), or a
    /// sketch is malformed: <c>K</c> outside 1..32 (Mash <c>-k</c> range), <c>Use64</c> ≠ (<c>K</c> &gt; 16) (Mash
    /// <c>use64 = 4^k &gt; 2^32</c>), a hash &gt; 2^32 − 1 in a 32-bit sketch, <c>Hashes</c> null, not strictly ascending
    /// (sorted, distinct), longer than its <c>SketchSize</c>, <c>SketchSize</c> &lt; 1 or <c>Length</c> &lt; 0 (build sketches with <see cref="CreateMinHashSketch(IEnumerable{string}, int, int, bool, uint)"/>
    /// or <see cref="MinHashSketch.FromHashes"/>).</exception>
    public static MashComparison CompareMinHashSketches(MinHashSketch reference, MinHashSketch query)
    {
        ArgumentNullException.ThrowIfNull(reference);
        ArgumentNullException.ThrowIfNull(query);
        if (reference.K != query.K || reference.Seed != query.Seed || reference.Canonical != query.Canonical)
            throw new ArgumentException("Sketches must have the same k-mer size, seed and canonical mode.", nameof(query));
        ValidateMinHashSketch(reference, nameof(reference));
        ValidateMinHashSketch(query, nameof(query));

        int sketchSize = Math.Min(reference.SketchSize, query.SketchSize);
        var r = reference.Hashes;
        var q = query.Hashes;
        int i = 0, j = 0, common = 0, denom = 0;
        while (denom < sketchSize && i < r.Count && j < q.Count)
        {
            if (r[i] < q[j])
                i++;
            else if (q[j] < r[i])
                j++;
            else
            {
                i++;
                j++;
                common++;
            }
            denom++;
        }
        if (denom < sketchSize)
        {
            denom += (r.Count - i) + (q.Count - j);
            if (denom > sketchSize)
                denom = sketchSize;
        }

        double jaccard = denom == 0 ? 0.0 : (double)common / denom;
        double distance = common == denom ? 0.0 : MashDistanceFromJaccard(jaccard, reference.K);
        double pValue = MashPValue(common, reference.Length, query.Length, reference.K, denom);
        return new MashComparison(common, denom, jaccard, distance, pValue);
    }

    private static void ValidateMinHashSketch(MinHashSketch sketch, string paramName)
    {
        // Mash Command.cpp: -k is an integer option with range 1..32; Sketch.cpp: use64 = alphabetSize^k > 2^32,
        // i.e. k > 16 for DNA, and a 32-bit sketch stores only 32-bit hash values.
        if (sketch.K < 1 || sketch.K > MaxMashKmerSize)
            throw new ArgumentException($"K must be in 1..{MaxMashKmerSize} (mash sketch -k).", paramName);
        if (sketch.Use64 != sketch.K > 16)
            throw new ArgumentException("Use64 must equal K > 16 (Mash: use64 = 4^k > 2^32).", paramName);
        if (sketch.SketchSize < 1 || sketch.Length < 0)
            throw new ArgumentException("SketchSize must be >= 1 and Length >= 0.", paramName);
        if (sketch.Hashes is null)
            throw new ArgumentException("Hashes must not be null.", paramName);
        if (sketch.Hashes.Count > sketch.SketchSize)
            throw new ArgumentException("A bottom-s sketch holds at most SketchSize hashes.", paramName);
        ThrowIfNotStrictlyAscending(sketch.Hashes, paramName);
        if (!sketch.Use64 && sketch.Hashes.Count > 0 && sketch.Hashes[^1] > uint.MaxValue)
            throw new ArgumentException("A 32-bit sketch (Use64 = false) holds only hash values <= 2^32 - 1.", paramName);
    }

    private static void ThrowIfNotStrictlyAscending(IReadOnlyList<ulong> hashes, string paramName)
    {
        for (int i = 1; i < hashes.Count; i++)
        {
            if (hashes[i] <= hashes[i - 1])
                throw new ArgumentException("Hashes must be sorted ascending and distinct.", paramName);
        }
    }

    /// <summary>
    /// Mash p-value of observing ≥ <paramref name="sharedHashes"/> shared hashes by chance (Mash <c>CommandDistance.cpp</c>
    /// <c>pValue</c>): P(X ≥ x), X ~ Binomial(<paramref name="sketchSize"/>, r), r = p₁p₂/(p₁ + p₂ − p₁p₂),
    /// p_i = 1/(1 + 4^k/length_i); x = 0 → 1.
    /// </summary>
    /// <param name="sharedHashes">x, the shared hashes.</param>
    /// <param name="length1">Total length of the first sequence set (Mash reference length).</param>
    /// <param name="length2">Total length of the second sequence set.</param>
    /// <param name="k">K-mer size, 1..32 (k-mer space 4^k; Mash <c>Command.cpp</c> declares <c>-k</c> as an integer option with
    /// range 1..32).</param>
    /// <param name="sketchSize">The comparison denominator (number of binomial trials).</param>
    /// <remarks>
    /// Inputs that no comparison can produce are rejected instead of being passed through Mash's formula: x &gt; s (the
    /// shared hashes are a subset of the s visited union hashes; Mash's <c>gsl_cdf_binomial_Q(x − 1, r, s)</c> silently
    /// returns 0 for x − 1 ≥ s) and x ≥ 1 with a length of 0 (an empty sequence set has no hashes; Mash would compute
    /// p = 1/(1 + 4^k/0) = 0, r = 0/0 = NaN, and GSL's beta CDF of NaN is NaN). k is limited to Mash's 1..32: for
    /// k ≥ 512, 4^k overflows a double to ∞, r = 0/0 = NaN and the binomial tail is undefined (k = 32 gives 4^k ≈ 1.8·10^19,
    /// still a positive r for any length ≥ 1).
    /// </remarks>
    /// <returns>The p-value in [0, 1].</returns>
    /// <exception cref="ArgumentOutOfRangeException">A count or length is negative, or <paramref name="k"/> is outside 1..32 (Mash's
    /// k-mer size range).</exception>
    /// <exception cref="ArgumentException"><paramref name="sharedHashes"/> exceeds <paramref name="sketchSize"/>, or is
    /// positive while a length is 0.</exception>
    public static double MashPValue(long sharedHashes, long length1, long length2, int k, long sketchSize)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(sharedHashes);
        ArgumentOutOfRangeException.ThrowIfNegative(length1);
        ArgumentOutOfRangeException.ThrowIfNegative(length2);
        ArgumentOutOfRangeException.ThrowIfNegative(sketchSize);
        if (k < 1 || k > MaxMashKmerSize)
            throw new ArgumentOutOfRangeException(nameof(k), k, $"K must be in 1..{MaxMashKmerSize} (mash -k).");
        if (sharedHashes > sketchSize)
            throw new ArgumentException("Shared hashes cannot exceed the sketch size (number of compared hashes).", nameof(sharedHashes));
        if (sharedHashes == 0)
            return 1.0;
        if (length1 == 0 || length2 == 0)
            throw new ArgumentException(
                "Shared hashes require two non-empty sequence sets (length >= 1); an empty set has no hashes.",
                length1 == 0 ? nameof(length1) : nameof(length2));

        double kmerSpace = Math.Pow(4.0, k);
        double pX = 1.0 / (1.0 + kmerSpace / length1);
        double pY = 1.0 / (1.0 + kmerSpace / length2);
        double r = pX * pY / (pX + pY - pX * pY);
        return StatisticsHelper.BinomialUpperTail(sharedHashes, sketchSize, r);
    }

    /// <summary>
    /// MurmurHash3_x64_128 (Austin Appleby, public domain; the <c>MurmurHash3.cpp</c> bundled with Mash and used by
    /// sourmash / <c>mmh3.hash64</c>), little-endian block reads.
    /// </summary>
    /// <param name="data">Bytes to hash.</param>
    /// <param name="seed">32-bit seed (zero-extended to both 64-bit lanes).</param>
    /// <returns>The two 64-bit output words (h1 first, as written to <c>out[0]</c>).</returns>
    public static (ulong H1, ulong H2) MurmurHash3X64_128(ReadOnlySpan<byte> data, uint seed)
    {
        const ulong c1 = 0x87c37b91114253d5UL;
        const ulong c2 = 0x4cf5ad432745937fUL;
        int len = data.Length;
        int nblocks = len / 16;
        ulong h1 = seed, h2 = seed;

        for (int i = 0; i < nblocks; i++)
        {
            ulong k1 = System.Buffers.Binary.BinaryPrimitives.ReadUInt64LittleEndian(data.Slice(i * 16, 8));
            ulong k2 = System.Buffers.Binary.BinaryPrimitives.ReadUInt64LittleEndian(data.Slice(i * 16 + 8, 8));

            k1 *= c1; k1 = ulong.RotateLeft(k1, 31); k1 *= c2; h1 ^= k1;
            h1 = ulong.RotateLeft(h1, 27); h1 += h2; h1 = h1 * 5 + 0x52dce729;
            k2 *= c2; k2 = ulong.RotateLeft(k2, 33); k2 *= c1; h2 ^= k2;
            h2 = ulong.RotateLeft(h2, 31); h2 += h1; h2 = h2 * 5 + 0x38495ab5;
        }

        var tail = data[(nblocks * 16)..];
        int rem = len & 15;
        if (rem > 8)
        {
            ulong k2 = 0;
            for (int t = rem - 1; t >= 8; t--)
                k2 ^= (ulong)tail[t] << ((t - 8) * 8);
            k2 *= c2; k2 = ulong.RotateLeft(k2, 33); k2 *= c1; h2 ^= k2;
        }
        if (rem > 0)
        {
            ulong k1 = 0;
            for (int t = Math.Min(rem, 8) - 1; t >= 0; t--)
                k1 ^= (ulong)tail[t] << (t * 8);
            k1 *= c1; k1 = ulong.RotateLeft(k1, 31); k1 *= c2; h1 ^= k1;
        }

        h1 ^= (ulong)len; h2 ^= (ulong)len;
        h1 += h2; h2 += h1;
        h1 = FMix64(h1); h2 = FMix64(h2);
        h1 += h2; h2 += h1;
        return (h1, h2);
    }

    private static ulong FMix64(ulong k)
    {
        k ^= k >> 33;
        k *= 0xff51afd7ed558ccdUL;
        k ^= k >> 33;
        k *= 0xc4ceb9fe1a85ec53UL;
        k ^= k >> 33;
        return k;
    }

    #endregion

    #region FracMinHash (sourmash scaled sketches)

    /// <summary>Default sourmash hash seed (<c>MinHash(seed=42)</c>).</summary>
    public const uint DefaultSourmashSeed = 42;

    /// <summary>Largest sourmash scaled factor (<c>ScaledType</c> = u32: 4294967295).</summary>
    public const long MaxSourmashScaled = uint.MaxValue;

    /// <summary>
    /// The FracMinHash threshold for a <paramref name="scaled"/> value, as sourmash 4.9.4 computes it (Rust core
    /// <c>sketch/minhash.rs</c> <c>max_hash_for_scaled</c>): scaled = 1 → 2^64 − 1; otherwise
    /// <c>(u64::MAX as f64 / scaled as f64) as u64</c>, i.e. ⌊2^64 / scaled⌋ evaluated in double precision and truncated.
    /// </summary>
    /// <remarks>
    /// <c>u64::MAX as f64</c> rounds to 2^64. The Python helper <c>_get_max_hash_for_scaled</c> rounds instead of truncating;
    /// the two agree while 2^64/scaled ≥ 2^53 (scaled ≤ 2048) and differ by one above (e.g. scaled 7919: Rust
    /// 2329428472497733, Python 2329428472497734). <c>MinHash(n=0, scaled=S)</c> passes S to the Rust constructor, and
    /// <c>MinHash.downsample(scaled=S)</c> converts its rounded max_hash back to S and does the same, so the Rust value is the
    /// one applied (sourmash <c>MinHash(0, 21, scaled=7919)._max_hash</c> and <c>.downsample(scaled=7919)._max_hash</c> =
    /// 2329428472497733; scaled 4294967295 → 4294967297).
    /// </remarks>
    /// <param name="scaled">The scaled factor S, 1 ≤ S ≤ 4294967295 (sourmash <c>ScaledType</c> is u32).</param>
    /// <returns>max_hash; a sketch keeps the hashes h ≤ max_hash (about 1/S of all hashes).</returns>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="scaled"/> is outside 1..4294967295.</exception>
    public static ulong FracMinHashMaxHash(long scaled)
    {
        ThrowIfScaledOutOfRange(scaled, nameof(scaled));
        return scaled == 1 ? ulong.MaxValue : (ulong)(18446744073709551616.0 / scaled);
    }

    private static void ThrowIfScaledOutOfRange(long scaled, string paramName)
    {
        if (scaled < 1 || scaled > MaxSourmashScaled)
            throw new ArgumentOutOfRangeException(paramName, scaled, $"Scaled must be in 1..{MaxSourmashScaled} (sourmash u32 scaled).");
    }

    /// <summary>
    /// Builds a FracMinHash ("scaled") sketch of one sequence exactly as sourmash 4.9.4
    /// <c>MinHash(n=0, ksize=k, scaled=S, seed=42, track_abundance=…).add_sequence(seq, force=True)</c> does (Irber, Brooks,
    /// Reiter et al. 2022, "Lightweight compositional analysis of metagenomes with FracMinHash and minimum metagenome covers",
    /// bioRxiv 2022.01.11.475838; Hera, Pierce-Ward &amp; Koslicki 2023, Genome Res 33:1061).
    /// </summary>
    /// <remarks>
    /// <para>K-mers: upper-cased; windows with a symbol other than A/C/G/T are skipped (sourmash <c>force=True</c>; without it
    /// sourmash raises on such a window); each k-mer replaced by the lexicographically smaller of itself and its reverse
    /// complement (Rust <c>signature.rs</c> <c>SeqToHashes</c>: <c>std::cmp::min(kmer, krc)</c> on the byte strings) — the
    /// counts of <see cref="CountKmers(string, int, KmerCountingOptions, CancellationToken, IProgress{double}?)"/> with
    /// <c>Canonical</c>. With <paramref name="canonical"/> = false the forward k-mers (<c>AcgtOnly</c>) are hashed instead;
    /// sourmash has no such DNA mode (it is Mash's <c>-n</c>).</para>
    /// <para>Hash: the first 64-bit word of <see cref="MurmurHash3X64_128"/> of the k ASCII bytes with
    /// <paramref name="seed"/> (sourmash <c>_hash_murmur</c>). The sketch keeps every distinct hash h ≤
    /// <see cref="FracMinHashMaxHash"/>(<paramref name="scaled"/>) (Rust <c>add_hash</c>: <c>hash &gt; max_hash</c> is dropped),
    /// ascending. <paramref name="scaled"/> = 1 keeps all hashes, so the sketch is the exact canonical k-mer set (up to
    /// 64-bit hash collisions).</para>
    /// <para><paramref name="trackAbundance"/> (sourmash <c>track_abundance=True</c>): every k-mer occurrence adds 1 to the
    /// abundance of its hash (Rust <c>add_hash_with_abundance</c>), so <see cref="FracMinHashSketch.Abundances"/>[i] is the
    /// (canonical) count of the k-mer(s) hashing to <c>Hashes[i]</c>; the hash set is unchanged.</para>
    /// </remarks>
    /// <param name="sequence">The sequence; null is treated as empty.</param>
    /// <param name="k">K-mer size ≥ 1 (sourmash default 31).</param>
    /// <param name="scaled">Scaled factor S in 1..4294967295 (sourmash default 1000).</param>
    /// <param name="canonical">True (default, sourmash): canonical k-mers; false: forward k-mers.</param>
    /// <param name="seed">MurmurHash3 seed (sourmash default 42).</param>
    /// <param name="trackAbundance">True: also record the abundance of each kept hash (sourmash <c>track_abundance</c>).</param>
    /// <returns>The sketch.</returns>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="k"/> is less than 1, or <paramref name="scaled"/> is outside
    /// 1..4294967295.</exception>
    public static FracMinHashSketch CreateFracMinHashSketch(
        string sequence, int k, long scaled, bool canonical = true, uint seed = DefaultSourmashSeed, bool trackAbundance = false)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(k, 1);
        ulong maxHash = FracMinHashMaxHash(scaled);
        var options = canonical ? new KmerCountingOptions(Canonical: true) : new KmerCountingOptions(AcgtOnly: true);
        var counts = CountKmers(sequence ?? string.Empty, k, options);

        var kept = new Dictionary<ulong, long>();
        var bytes = new byte[k];
        foreach (var (kmer, count) in counts)
        {
            for (int i = 0; i < k; i++)
                bytes[i] = (byte)kmer[i];
            ulong h = MurmurHash3X64_128(bytes, seed).H1;
            if (h <= maxHash && !kept.TryAdd(h, count))
                kept[h] += count;
        }

        var hashes = kept.Keys.ToArray();
        Array.Sort(hashes);
        long[]? abundances = null;
        if (trackAbundance)
        {
            abundances = new long[hashes.Length];
            for (int i = 0; i < hashes.Length; i++)
                abundances[i] = kept[hashes[i]];
        }
        return new FracMinHashSketch(k, scaled, maxHash, canonical, seed, hashes, abundances);
    }

    /// <summary>
    /// Downsamples a FracMinHash sketch to a larger scaled factor, as sourmash 4.9.4 <c>MinHash.downsample(scaled=S′)</c>
    /// (Rust <c>downsample_scaled</c>): a new sketch at S′ holding the hashes h ≤ <see cref="FracMinHashMaxHash"/>(S′), with
    /// their abundances when tracked — the sketch that sketching the sequence at S′ would give.
    /// </summary>
    /// <param name="sketch">The sketch.</param>
    /// <param name="newScaled">The target scaled factor S′ ≥ <c>sketch.Scaled</c> (S′ = S returns an equal sketch).</param>
    /// <returns>The downsampled sketch.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="sketch"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="newScaled"/> is outside 1..4294967295.</exception>
    /// <exception cref="ArgumentException"><paramref name="newScaled"/> &lt; <c>sketch.Scaled</c> (sourmash: "new scaled … is lower
    /// than current sample scaled", Rust <c>CannotUpsampleScaled</c>), or the sketch is malformed (see
    /// <see cref="CompareFracMinHashSketches"/>).</exception>
    public static FracMinHashSketch DownsampleFracMinHash(FracMinHashSketch sketch, long newScaled)
    {
        ArgumentNullException.ThrowIfNull(sketch);
        ThrowIfScaledOutOfRange(newScaled, nameof(newScaled));
        ValidateFracMinHashSketch(sketch, nameof(sketch));
        if (newScaled < sketch.Scaled)
            throw new ArgumentException(
                $"New scaled {newScaled} is lower than the sketch's scaled {sketch.Scaled} (a FracMinHash sketch cannot be upsampled).",
                nameof(newScaled));

        ulong maxHash = FracMinHashMaxHash(newScaled);
        var hashes = sketch.Hashes;
        int n = 0;
        while (n < hashes.Count && hashes[n] <= maxHash)
            n++;
        var kept = new ulong[n];
        for (int i = 0; i < n; i++)
            kept[i] = hashes[i];
        long[]? abundances = null;
        if (sketch.Abundances is { } ab)
        {
            abundances = new long[n];
            for (int i = 0; i < n; i++)
                abundances[i] = ab[i];
        }
        return sketch with { Scaled = newScaled, MaxHash = maxHash, Hashes = kept, Abundances = abundances };
    }

    /// <summary>
    /// Compares two FracMinHash sketches as sourmash 4.9.4 does: Jaccard (<c>MinHash.jaccard</c>), containment of each in the
    /// other (<c>contained_by</c>), maximum containment (<c>max_containment</c>) and, when both track abundance, the angular
    /// similarity (<c>angular_similarity</c> = <c>similarity(ignore_abundance=False)</c>) and the abundance-weighted
    /// containments (<c>contained_by_weighted</c>).
    /// </summary>
    /// <remarks>
    /// <para>x = |A ∩ B| (shared hashes), u = |A ∪ B|. Jaccard = x / max(1, u) (Rust <c>KmerMinHash::jaccard</c>; 0 for two
    /// empty sketches).</para>
    /// <para>Containment of A in B (sourmash <c>A.contained_by(B)</c>, Python <c>minhash.py</c>): |A| = 0 → 0; otherwise
    /// x / (|A|·b) with the bias factor b = 1 − (1 − 1/S)^(|A|·S) (Hera et al. 2023), clamped to [0, 1]. Maximum
    /// containment uses min(|A|, |B|) in place of |A|. For S = 1, b = 1 and these are the exact ratios of
    /// <see cref="ContainmentIndex(string, string, int, KmerCountingOptions)"/>.</para>
    /// <para>Angular similarity (Rust <c>KmerMinHash::angular_similarity</c>): with abundance vectors a, b over each sketch's
    /// own hashes, cos = min(1, Σ_{h∈A∩B} a_h·b_h / (‖a‖·‖b‖)) (‖·‖ over all of each sketch's abundances; 0 when a norm
    /// is 0) and similarity = 1 − 2·acos(cos)/π. Weighted containment of A in B (<c>contained_by_weighted</c>):
    /// Σ_{h∈A∩B} a_h / Σ_{h∈A} a_h (not bias-corrected; 0 for an empty A, where sourmash divides by zero). These are null
    /// unless the sketch(es) carry <see cref="FracMinHashSketch.Abundances"/>: angular needs both (sourmash raises
    /// otherwise), the weighted containment of A needs A's.</para>
    /// <para>Different <c>scaled</c>: refused unless <paramref name="downsample"/> is true (sourmash <c>downsample=True</c>);
    /// then both sketches are downsampled to S = max(S_A, S_B) with <see cref="DownsampleFracMinHash"/> and every value is
    /// computed at S. Jaccard and angular similarity equal sourmash's <c>jaccard</c> / <c>angular_similarity(…,
    /// downsample=True)</c> (Rust downsamples the smaller-scaled sketch). The containments equal
    /// <c>sourmash compare --containment / --max-containment</c>, which downsamples every signature to the common maximum
    /// scaled first (<c>commands.py</c>), i.e. <c>A.downsample(scaled=S).contained_by(B.downsample(scaled=S))</c>. The Python
    /// methods <c>contained_by(…, downsample=True)</c> / <c>max_containment(…, downsample=True)</c> differ when <c>self</c> has
    /// the smaller scaled: only the shared count is downsampled while the denominator keeps <c>len(self)</c> and
    /// <c>self.scaled</c> of the undownsampled sketch (e.g. sourmash 4.9.4: 0.0291 instead of 0.3155), so they are not
    /// reproduced.</para>
    /// </remarks>
    /// <param name="a">First sketch (sourmash <c>self</c>).</param>
    /// <param name="b">Second sketch (<c>other</c>).</param>
    /// <param name="downsample">True: compare sketches of different scaled at the larger one (sourmash <c>downsample=True</c>).</param>
    /// <returns>Shared and union hash counts, Jaccard, both containments, the maximum containment and the abundance metrics.</returns>
    /// <exception cref="ArgumentNullException">A sketch is null.</exception>
    /// <exception cref="ArgumentException">The sketches differ in k, seed or canonical mode, or in scaled with
    /// <paramref name="downsample"/> false, or a sketch is malformed (<c>Scaled</c> outside 1..4294967295, <c>MaxHash</c> ≠
    /// max_hash(<c>Scaled</c>), <c>Hashes</c> null or not strictly ascending, a hash above <c>MaxHash</c>, or <c>Abundances</c>
    /// not one positive count per hash).</exception>
    public static FracMinHashComparison CompareFracMinHashSketches(FracMinHashSketch a, FracMinHashSketch b, bool downsample = false)
    {
        ArgumentNullException.ThrowIfNull(a);
        ArgumentNullException.ThrowIfNull(b);
        ValidateFracMinHashSketch(a, nameof(a));
        ValidateFracMinHashSketch(b, nameof(b));
        if (a.K != b.K || a.Seed != b.Seed || a.Canonical != b.Canonical)
            throw new ArgumentException("Sketches must have the same k-mer size, seed and canonical mode.", nameof(b));
        if (a.Scaled != b.Scaled)
        {
            if (!downsample)
                throw new ArgumentException(
                    "Sketches must have the same scaled (pass downsample: true to compare at the larger scaled, sourmash downsample=True).",
                    nameof(b));
            long common = Math.Max(a.Scaled, b.Scaled);
            a = DownsampleFracMinHash(a, common);
            b = DownsampleFracMinHash(b, common);
        }

        var x = a.Hashes;
        var y = b.Hashes;
        var xa = a.Abundances;
        var ya = b.Abundances;
        int i = 0, j = 0, shared = 0;
        ulong dot = 0;
        long sharedA = 0, sharedB = 0;
        while (i < x.Count && j < y.Count)
        {
            if (x[i] < y[j])
                i++;
            else if (y[j] < x[i])
                j++;
            else
            {
                shared++;
                if (xa is not null)
                    sharedA += xa[i];
                if (ya is not null)
                    sharedB += ya[j];
                if (xa is not null && ya is not null)
                    dot = unchecked(dot + (ulong)xa[i] * (ulong)ya[j]);
                i++;
                j++;
            }
        }

        int union = x.Count + y.Count - shared;
        double jaccard = (double)shared / Math.Max(1, union);
        double? angular = xa is not null && ya is not null ? AngularSimilarity(dot, xa, ya) : null;
        return new FracMinHashComparison(
            shared,
            union,
            jaccard,
            DebiasedContainment(shared, x.Count, a.Scaled),
            DebiasedContainment(shared, y.Count, a.Scaled),
            DebiasedContainment(shared, Math.Min(x.Count, y.Count), a.Scaled),
            angular,
            xa is null ? null : WeightedContainment(sharedA, xa),
            ya is null ? null : WeightedContainment(sharedB, ya));
    }

    // Rust KmerMinHash::angular_similarity: u64 sums of squares over each sketch's abundances, prod over shared hashes.
    private static double AngularSimilarity(ulong dot, IReadOnlyList<long> a, IReadOnlyList<long> b)
    {
        ulong aSq = 0, bSq = 0;
        foreach (var v in a)
            aSq = unchecked(aSq + (ulong)v * (ulong)v);
        foreach (var v in b)
            bSq = unchecked(bSq + (ulong)v * (ulong)v);
        double normA = Math.Sqrt(aSq);
        double normB = Math.Sqrt(bSq);
        if (normA == 0 || normB == 0)
            return 0.0;
        double cos = Math.Min(dot / (normA * normB), 1.0);
        return 1.0 - 2.0 * Math.Acos(cos) / Math.PI;
    }

    // sourmash minhash.py contained_by_weighted: inflate(other ∩ self).sum_abundances / self.sum_abundances.
    private static double WeightedContainment(long sharedAbundance, IReadOnlyList<long> abundances)
    {
        long total = 0;
        foreach (var v in abundances)
            total += v;
        return total == 0 ? 0.0 : (double)sharedAbundance / total;
    }

    // sourmash minhash.py contained_by / max_containment: x / (denom · (1 − (1 − 1/S)^(denom·S))), clamped to [0, 1].
    private static double DebiasedContainment(int shared, int denom, long scaled)
    {
        if (denom == 0)
            return 0.0;
        double totalDenom = (double)denom * scaled;
        double biasFactor = 1.0 - Math.Pow(1.0 - 1.0 / scaled, totalDenom);
        double containment = shared / (denom * biasFactor);
        if (containment >= 1)
            return 1.0;
        return containment <= 0 ? 0.0 : containment;
    }

    private static void ValidateFracMinHashSketch(FracMinHashSketch sketch, string paramName)
    {
        if (sketch.Scaled < 1 || sketch.Scaled > MaxSourmashScaled || sketch.MaxHash != FracMinHashMaxHash(sketch.Scaled))
            throw new ArgumentException(
                $"Scaled must be in 1..{MaxSourmashScaled} and MaxHash must equal FracMinHashMaxHash(Scaled).", paramName);
        if (sketch.Hashes is null)
            throw new ArgumentException("Hashes must not be null.", paramName);
        ThrowIfNotStrictlyAscending(sketch.Hashes, paramName);
        if (sketch.Hashes.Count > 0 && sketch.Hashes[^1] > sketch.MaxHash)
            throw new ArgumentException("A FracMinHash sketch holds only hashes <= MaxHash.", paramName);
        if (sketch.Abundances is { } ab)
        {
            if (ab.Count != sketch.Hashes.Count)
                throw new ArgumentException("Abundances must hold one count per hash.", paramName);
            foreach (var v in ab)
            {
                if (v < 1)
                    throw new ArgumentException("Abundances must be positive (sourmash removes a hash whose abundance is 0).", paramName);
            }
        }
    }

    #endregion

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
        => CountSpacedWords(sequence, pattern, KmerCountingOptions.Default);

    /// <summary>
    /// Counts the spaced words of <paramref name="sequence"/> under <see cref="KmerCountingOptions"/>: literal (default) or
    /// ACGT-only — the word rule of the <c>spaced</c> program.
    /// </summary>
    /// <remarks>
    /// <para><b>AcgtOnly</b>: a window whose symbol at some <b>match</b> position is not A/C/G/T (after case folding) gives no
    /// word; symbols at don't-care positions are ignored. This is <c>spaced</c> 1.2.0 (<c>src/sort.h</c> <c>spacedDNA</c>):
    /// every letter other than A/C/G/T is stored as 'N', and a word is kept only while <c>correctWord</c> holds, i.e. no
    /// match position reads 'N'. For the all-'1' pattern this is the k-mer ACGT-only window rule of
    /// <see cref="CountKmers(string, int, KmerCountingOptions, CancellationToken, IProgress{double}?)"/>. The counts then
    /// sum to the number of kept words (≤ L − ℓ + 1). This method windows the string as given; <c>spaced</c>'s reader also
    /// deletes every non-letter character (gaps '-', '*', digits, blanks) before windowing, which
    /// <see cref="SpacedWordDistance(string, string, IReadOnlyList{string}, KmerDistanceMetric, KmerCountingOptions, bool)"/>
    /// does in its <c>AcgtOnly</c> (spaced-faithful) mode.</para>
    /// <para><b>Canonical</b> is rejected (<see cref="ArgumentException"/>): the reverse strand's spaced word at a window is
    /// read with the mirrored pattern, so min(word, RC(word)) is not a strand-independent key unless the pattern is a
    /// palindrome, and neither the paper nor <c>spaced</c> defines it. <c>spaced</c>'s reverse-complement handling is the
    /// both-strand comparison of
    /// <see cref="SpacedWordDistance(string, string, IReadOnlyList{string}, KmerDistanceMetric, KmerCountingOptions, bool)"/>.</para>
    /// </remarks>
    /// <param name="sequence">The sequence; null/empty or shorter than the pattern gives an empty table.</param>
    /// <param name="pattern">Binary pattern over {'0','1'} that starts and ends with '1', e.g. "1101".</param>
    /// <param name="options">Counting mode; <c>Canonical</c> must be false.</param>
    /// <returns>Dictionary mapping spaced words (upper-case match-position symbols) to their counts.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="pattern"/> is null.</exception>
    /// <exception cref="ArgumentException"><paramref name="pattern"/> is malformed, or <paramref name="options"/> has
    /// <c>Canonical = true</c>.</exception>
    public static Dictionary<string, int> CountSpacedWords(string sequence, string pattern, KmerCountingOptions options)
    {
        ValidateSpacedPattern(pattern);
        ThrowIfCanonicalSpaced(options);

        var counts = new Dictionary<string, int>();
        if (string.IsNullOrEmpty(sequence) || pattern.Length > sequence.Length)
            return counts;

        int[] matchPositions = Enumerable.Range(0, pattern.Length).Where(j => pattern[j] == '1').ToArray();
        var seq = sequence.ToUpperInvariant();
        var buffer = new char[matchPositions.Length];
        var lookup = counts.GetAlternateLookup<ReadOnlySpan<char>>();

        for (int i = 0; i <= seq.Length - pattern.Length; i++)
        {
            bool keep = true;
            for (int m = 0; m < matchPositions.Length; m++)
            {
                char c = seq[i + matchPositions[m]];
                if (options.AcgtOnly && !IsAcgt(c))
                {
                    keep = false;
                    break;
                }
                buffer[m] = c;
            }

            if (keep)
                CollectionsMarshal.GetValueRefOrAddDefault(lookup, buffer, out _)++;
        }

        return counts;
    }

    private static void ValidateSpacedPattern(string pattern)
    {
        ArgumentNullException.ThrowIfNull(pattern);
        if (pattern.Length == 0 || pattern[0] != '1' || pattern[^1] != '1' || pattern.Any(c => c is not ('0' or '1')))
            throw new ArgumentException("Pattern must be a non-empty string over {0,1} that starts and ends with '1'.", nameof(pattern));
    }

    private static void ThrowIfCanonicalSpaced(KmerCountingOptions options)
    {
        if (options.Canonical)
            throw new ArgumentException(
                "Canonical spaced words are not defined (the reverse strand is read with the mirrored pattern); use bothStrands = true for the spaced program's reverse-complement mode.",
                nameof(options));
    }

    /// <summary>
    /// Multiple-pattern spaced-word distance (Leimeister, Boden, Horwege, Lindner &amp; Morgenstern 2014,
    /// Bioinformatics 30:1991): the average, over a set of patterns of equal weight, of the per-pattern word-vector
    /// distance between the two spaced-word count tables.
    /// </summary>
    /// <remarks>
    /// <para>d_P(S₁, S₂) = (1/m) Σ_{i=1..m} d(N_{P_i}(S₁), N_{P_i}(S₂)), where N_P(S) is
    /// <see cref="CountSpacedWords(string, string)"/> and d is
    /// <see cref="KmerDistance(IReadOnlyDictionary{string, int}, IReadOnlyDictionary{string, int}, KmerDistanceMetric)"/>
    /// (the paper: "the distance between two sequences is the average of the distances based on the individual
    /// patterns"). The paper applies the Euclidean distance (<see cref="KmerDistanceMetric.Euclidean"/>, relative
    /// frequencies, the default) and the Jensen–Shannon distance (<see cref="KmerDistanceMetric.JensenShannon"/>) to
    /// relative spaced-word frequencies.</para>
    /// <para>Reference program: <c>spaced</c> 1.2.0 (Debian/Ubuntu package source, <c>sort.h</c> <c>spacedDNA</c>), run
    /// with <c>-r</c> (single strand) and <c>-f</c> (fixed patterns): its <c>-d JS</c> output equals
    /// <see cref="KmerDistanceMetric.JensenShannon"/>, and its <c>-d EU</c> output equals
    /// <see cref="KmerDistanceMetric.EuclideanCounts"/> (the program takes the Euclidean distance of raw counts, not of
    /// frequencies). This overload counts literal words (<see cref="CountSpacedWords(string, string)"/>); on sequences
    /// with symbols other than A/C/G/T, <c>spaced</c> drops a word with such a symbol at a match position, which is
    /// <c>KmerCountingOptions.AcgtOnly</c> of the 6-argument overload. <c>spaced</c>'s default mode (no <c>-r</c>), which
    /// compares one sequence's forward words with both strands of the other, is that overload's <c>bothStrands</c>
    /// (docs/algorithms/K-mer/K-mer_Euclidean_Distance.md §7.5, §7.7).</para>
    /// </remarks>
    /// <param name="seq1">First sequence (case-insensitive; null/empty or shorter than a pattern gives an empty table).</param>
    /// <param name="seq2">Second sequence, same conventions.</param>
    /// <param name="patterns">One or more binary patterns (see <see cref="CountSpacedWords(string, string)"/>), all of the
    /// same weight (number of '1'). The all-'1' pattern of length k gives the contiguous k-mer distance.</param>
    /// <param name="metric">The per-pattern word-vector metric; any metric of the count-table overload (not D2*/D2S), or
    /// <see cref="KmerDistanceMetric.SpacedEvolutionary"/> (<c>spaced -d EV</c>, always read in the spaced-faithful mode; see the
    /// 6-argument overload).</param>
    /// <returns>The mean of the per-pattern values.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="patterns"/> or one of its patterns is null.</exception>
    /// <exception cref="ArgumentException"><paramref name="patterns"/> is empty, a pattern is malformed, the weights
    /// differ, or <paramref name="metric"/> is D2*/D2S.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="metric"/> is undefined.</exception>
    public static double SpacedWordDistance(
        string seq1,
        string seq2,
        IReadOnlyList<string> patterns,
        KmerDistanceMetric metric = KmerDistanceMetric.Euclidean)
        => SpacedWordDistance(seq1, seq2, patterns, metric, KmerCountingOptions.Default, bothStrands: false);

    /// <summary>
    /// Multiple-pattern spaced-word distance with the word rule of <see cref="CountSpacedWords(string, string, KmerCountingOptions)"/>
    /// and, optionally, the reverse-complement mode of the <c>spaced</c> program; also the <c>spaced -d EV</c> evolutionary
    /// distance (<see cref="KmerDistanceMetric.SpacedEvolutionary"/>).
    /// </summary>
    /// <remarks>
    /// <para>For each pattern P (length ℓ) the two word tables are compared with the word-vector metric, and the values are
    /// averaged over the patterns, as in
    /// <see cref="SpacedWordDistance(string, string, IReadOnlyList{string}, KmerDistanceMetric)"/>. Frequencies are
    /// count ÷ W with W = L − ℓ + 1, the number of pattern windows (0 when L &lt; ℓ). That is <c>spaced</c>'s denominator
    /// (<c>src/sort.h</c>: <c>row[i] /= seqWordEnd − seqStart</c>). With the default literal words every window gives one
    /// word, so W = Σc and the result equals the 4-argument overload. With <c>AcgtOnly</c> the dropped words still count in W,
    /// as in <c>spaced</c>, so the frequencies sum to less than 1.</para>
    /// <para><b>Spaced-faithful mode = <c>AcgtOnly</c></b> (with or without <paramref name="bothStrands"/>). Both sequences are
    /// first read as <c>spacedDNA</c> reads a FASTA record: every character that is not an ASCII letter (gap '-', '*',
    /// digits, blanks) is deleted, letters are upper-cased, and every letter other than A/C/G/T becomes N. L and W are taken
    /// on the read sequence, so "ACG-TACGT" is windowed as "ACGTACGT" (8 letters). The literal default mode windows the string
    /// as given (a gap is a symbol) and is not <c>spaced</c>'s rule.</para>
    /// <para><b><paramref name="bothStrands"/></b> reproduces <c>spaced</c>'s default mode (run without <c>-r</c>), in which the
    /// sequence that comes first in the input is compared on both strands and the other on its forward strand only.
    /// In <c>spacedDNA</c> the matrix entry d[i][j] (i &gt; j) uses row_i = forward counts of sequence i and, for sequence j,
    /// counts on the forward strand plus counts on the reverse-complement strand (the reverse-complement string read with
    /// the same pattern). For <c>-d EU</c> it takes |row_i − (row_j + row_j′)|; for <c>-d JS</c> it takes
    /// (row_j + row_j′) ÷ (2·W_j) against row_i ÷ W_i. Here <paramref name="seq1"/> plays sequence j (both strands; the first
    /// FASTA record) and <paramref name="seq2"/> sequence i (forward strand; the second record). For every metric, the
    /// seq1 vector becomes F₁ + R₁ with total 2·W₁, and the seq2 vector stays F₂ with total W₂. The value therefore
    /// depends on the argument order. This is a convention of the tool: the paper (Leimeister et al. 2014) does not
    /// define it. Exact <c>spaced</c> output needs <c>AcgtOnly = true</c>. Metric correspondence:
    /// <see cref="KmerDistanceMetric.JensenShannon"/> = <c>-d JS</c>, <see cref="KmerDistanceMetric.EuclideanCounts"/> =
    /// <c>-d EU</c>, <see cref="KmerDistanceMetric.SpacedEvolutionary"/> = <c>-d EV</c>.</para>
    /// <para><b><see cref="KmerDistanceMetric.SpacedEvolutionary"/></b> (Morgenstern, Zhu, Horwege &amp; Leimeister 2015,
    /// Algorithms Mol Biol 10:5), exactly as <c>spaced</c> 1.2.0 <c>sort.h</c> (EV branch) computes it; always read in the
    /// spaced-faithful mode (<c>AcgtOnly</c> implied). All patterns must have the same length ℓ (the program keeps one
    /// weight w and one don't-care count). N = Σ_P Σ_w min(c₂(w), c₁(w)) is the number of spaced-word matches summed over the
    /// patterns, with c₂ the forward counts of <paramref name="seq2"/> and c₁ those of <paramref name="seq1"/> (forward +
    /// reverse strand with <paramref name="bothStrands"/>). With L₁, L₂ the read lengths (N letters included),
    /// m = min(L₁, L₂) − ℓ + 1, M = max(L₁, L₂) − ℓ + 1, base frequencies f(a) = count(a) ÷ L and background match probability
    /// q = Σ_a f₁(a)·f₂(a) (with <paramref name="bothStrands"/>, each f(a) is averaged with f(complement a)), the value under
    /// the root is V = N ÷ (|P|·m) − s·M·q^w, s = 2 with <paramref name="bothStrands"/> and 1 without. If V ≥ 0 the
    /// match probability per site is p = V^(1/w) and the distance d = −¾·ln(4p/3 − 1/3) (Jukes–Cantor); otherwise
    /// <c>spaced</c> prints the saturation value <see cref="SpacedEvolutionarySaturationDistance"/> (1.2). As in the
    /// program, p &lt; ¼ gives NaN and p = ¼ gives +∞ (e.g. a sequence of N only).</para>
    /// <para>Cross-checked against the <c>spaced</c> 1.2.0 binary (Ubuntu archive) with and without <c>-r</c>: 248 runs on 11
    /// sequence pairs (with N, IUPAC, lower case, gaps, '*' and digits) × up to 7 pattern sets × EV/JS/EU, all equal to the
    /// 12 printed digits (docs/algorithms/K-mer/K-mer_Euclidean_Distance.md §7.7, §7.8).</para>
    /// </remarks>
    /// <param name="seq1">First sequence (both strands when <paramref name="bothStrands"/>).</param>
    /// <param name="seq2">Second sequence (forward strand).</param>
    /// <param name="patterns">One or more binary patterns of equal weight (and, for EV, equal length).</param>
    /// <param name="metric">The per-pattern word-vector metric (not D2*/D2S), or <see cref="KmerDistanceMetric.SpacedEvolutionary"/>.</param>
    /// <param name="options">Word rule; <c>AcgtOnly</c> = the spaced-faithful reader and N-word rule; <c>Canonical</c> must be false.</param>
    /// <param name="bothStrands">The <c>spaced</c> reverse-complement mode (seq1 on both strands vs seq2 forward).</param>
    /// <returns>The mean of the per-pattern values, or the EV distance.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="patterns"/> or one of its patterns is null.</exception>
    /// <exception cref="ArgumentException"><paramref name="patterns"/> is empty, a pattern is malformed, the weights differ,
    /// <paramref name="metric"/> is D2*/D2S, <paramref name="options"/> has <c>Canonical = true</c>, or (EV) the pattern
    /// lengths differ or a read sequence is shorter than the pattern.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="metric"/> is undefined.</exception>
    public static double SpacedWordDistance(
        string seq1,
        string seq2,
        IReadOnlyList<string> patterns,
        KmerDistanceMetric metric,
        KmerCountingOptions options,
        bool bothStrands = false)
    {
        ArgumentNullException.ThrowIfNull(patterns);
        if (patterns.Count == 0)
            throw new ArgumentException("At least one pattern is required.", nameof(patterns));

        int weight = -1;
        foreach (var pattern in patterns)
        {
            if (pattern is null)
                throw new ArgumentNullException(nameof(patterns), "Patterns must not contain null.");
            ValidateSpacedPattern(pattern);
            int w = pattern.Count(c => c == '1');
            if (weight >= 0 && w != weight)
                throw new ArgumentException("All patterns must have the same weight (number of '1' positions).", nameof(patterns));
            weight = w;
        }

        ThrowIfCanonicalSpaced(options);
        if (!Enum.IsDefined(metric))
            throw new ArgumentOutOfRangeException(nameof(metric), metric, "Unknown k-mer distance metric.");
        if (metric is KmerDistanceMetric.D2Star or KmerDistanceMetric.D2Shepherd)
            throw new ArgumentException(
                "D2*/D2S need each sequence's background model; use KmerDistance(string, string, int, metric) or BackgroundAdjustedD2.",
                nameof(metric));

        seq1 ??= string.Empty; // null = empty sequence (zero vector), as before
        seq2 ??= string.Empty;
        if (metric is KmerDistanceMetric.SpacedEvolutionary)
            return SpacedEvolutionaryDistance(SpacedRead(seq1), SpacedRead(seq2), patterns, weight, bothStrands);

        if (options.AcgtOnly)
        {
            seq1 = SpacedRead(seq1);
            seq2 = SpacedRead(seq2);
        }

        string reverse1 = bothStrands ? DnaSequence.GetReverseComplementString(seq1) : string.Empty;
        double sum = 0;
        foreach (var pattern in patterns)
        {
            var counts1 = CountSpacedWords(seq1, pattern, options);
            var counts2 = CountSpacedWords(seq2, pattern, options);
            double windows1 = Math.Max(0, seq1.Length - pattern.Length + 1);
            double windows2 = Math.Max(0, seq2.Length - pattern.Length + 1);
            if (bothStrands)
            {
                foreach (var (word, count) in CountSpacedWords(reverse1, pattern, options))
                    CollectionsMarshal.GetValueRefOrAddDefault(counts1, word, out _) += count;
                windows1 *= 2;
            }

            sum += WordVectorDistance(counts1, windows1, counts2, windows2, metric);
        }
        return sum / patterns.Count;
    }

    /// <summary>
    /// The value <c>spaced</c> 1.2.0 prints for <c>-d EV</c> when the estimated match rate is undefined (negative value under
    /// the root: fewer spaced-word matches than expected by chance), <c>sort.h</c> <c>dmat[0][i][j] = 1.2</c>.
    /// </summary>
    public const double SpacedEvolutionarySaturationDistance = 1.2;

    /// <summary>
    /// <c>spacedDNA</c>'s FASTA reader applied to one sequence: non-letters (ASCII <c>isalpha</c> false) are deleted, letters
    /// upper-cased, every letter other than A/C/G/T stored as N.
    /// </summary>
    private static string SpacedRead(string sequence)
    {
        var read = new StringBuilder(sequence.Length);
        foreach (char ch in sequence)
        {
            if (!char.IsAsciiLetter(ch))
                continue;
            char c = char.ToUpperInvariant(ch);
            read.Append(IsAcgt(c) ? c : 'N');
        }
        return read.ToString();
    }

    /// <summary><c>spaced</c> 1.2.0 <c>sort.h</c> EV branch on two read sequences (see <see cref="KmerDistanceMetric.SpacedEvolutionary"/>).</summary>
    private static double SpacedEvolutionaryDistance(
        string read1, string read2, IReadOnlyList<string> patterns, int weight, bool bothStrands)
    {
        int patternLength = patterns[0].Length;
        foreach (var pattern in patterns)
        {
            if (pattern.Length != patternLength)
                throw new ArgumentException(
                    "The spaced EV distance needs all patterns of the same length (spaced keeps one weight and one don't-care count).",
                    nameof(patterns));
        }
        if (read1.Length < patternLength || read2.Length < patternLength)
            throw new ArgumentException(
                "The spaced EV distance needs both sequences to have at least as many letters as the pattern length.",
                read1.Length < patternLength ? "seq1" : "seq2");

        var acgtOnly = new KmerCountingOptions(AcgtOnly: true);
        string reverse1 = bothStrands ? DnaSequence.GetReverseComplementString(read1) : string.Empty;
        double matches = 0;
        foreach (var pattern in patterns)
        {
            var counts1 = CountSpacedWords(read1, pattern, acgtOnly);
            if (bothStrands)
            {
                foreach (var (word, count) in CountSpacedWords(reverse1, pattern, acgtOnly))
                    CollectionsMarshal.GetValueRefOrAddDefault(counts1, word, out _) += count;
            }
            foreach (var (word, count2) in CountSpacedWords(read2, pattern, acgtOnly))
            {
                if (counts1.TryGetValue(word, out int count1))
                    matches += Math.Min(count1, count2);
            }
        }

        int ell = patternLength - 1;
        double min = Math.Min(read1.Length, read2.Length) - ell;
        double max = Math.Max(read1.Length, read2.Length) - ell;
        double[] f1 = BaseFrequencies(read1);
        double[] f2 = BaseFrequencies(read2);
        double q = 0;
        for (int a = 0; a < 4; a++)
        {
            if (bothStrands)
                q += ((f2[a] + f2[3 - a]) * 0.5) * ((f1[a] + f1[3 - a]) * 0.5);
            else
                q += f2[a] * f1[a];
        }

        double underRoot = bothStrands
            ? matches / (patterns.Count * min) - 2 * max * Math.Pow(q, weight)
            : matches / (patterns.Count * min) - max * Math.Pow(q, weight);
        if (underRoot >= 0)
        {
            double p = Math.Pow(underRoot, 1.0 / weight);
            return -0.75 * Math.Log((4.0 / 3.0) * p - (1.0 / 3.0));
        }
        return SpacedEvolutionarySaturationDistance;

        // A, C, G, T counts ÷ read length (N letters count in the length only), spaced's frequencies[k].
        static double[] BaseFrequencies(string read)
        {
            var f = new double[4];
            foreach (char c in read)
            {
                switch (c)
                {
                    case 'A': f[0]++; break;
                    case 'C': f[1]++; break;
                    case 'G': f[2]++; break;
                    case 'T': f[3]++; break;
                }
            }
            for (int a = 0; a < 4; a++)
                f[a] /= read.Length;
            return f;
        }
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
        => CountKmersBothStrands(sequence, k, KmerCountingOptions.Default);

    /// <summary>
    /// Counts every k-mer over both strands (count[w] = forward[w] + forward[RC(w)], palindromes doubled) under
    /// <see cref="KmerCountingOptions"/>: literal (default) or ACGT-only — the kPAL convention.
    /// </summary>
    /// <remarks>
    /// <para>The forward table is
    /// <see cref="CountKmers(string, int, KmerCountingOptions, CancellationToken, IProgress{double}?)"/>; each entry is then
    /// also added to its reverse complement's entry, which is kPAL's <c>Profile.balance()</c> (<c>kpal/klib.py</c>:
    /// <c>counts[i] += counts[i_rc]; counts[i_rc] += temp</c>, and <c>counts[i] += counts[i]</c> for a palindrome).
    /// Counting the forward strand's windows of RC(w) equals counting w on the reverse strand, so the result is the same
    /// as counting the reverse-complement string.</para>
    /// <para><b>AcgtOnly</b>: windows containing a non-ACGT symbol are skipped on both strands. kPAL does this:
    /// <c>Profile.from_sequences</c> splits every sequence on the regular expression <c>[^AaCcGgTt]</c> and counts only
    /// the k-mers inside the parts. Total = 2 × (number of all-ACGT windows). Cross-checked against kPAL run from source
    /// (e.g. ACGTNACGTAAcgtRTT, k = 3 → ACG = CGT = 6, AAC = GTA = GTT = TAA = TAC = TTA = 1, Σ 18).</para>
    /// <para><b>Canonical</b> is rejected (<see cref="ArgumentException"/>). Canonical counting already is the both-strand
    /// count: Jellyfish <c>count -C</c> (<see cref="CountKmers(string, int, KmerCountingOptions, CancellationToken, IProgress{double}?)"/>
    /// with <c>Canonical = true</c>) keys forward[w] + forward[RC(w)] by min(w, RC(w)) and counts a palindrome once.
    /// Restricting this balanced table to canonical keys would instead double the palindromes (e.g. GAATTCNNACGTTGCAGGATCCATGCRYacgtgcaNTTGCA,
    /// k = 4: AATT 2 here vs 1 in <c>jellyfish -C</c>), a combination that neither kPAL (which has no canonical profiles)
    /// nor Jellyfish defines.</para>
    /// </remarks>
    /// <param name="sequence">The DNA sequence (case-insensitive). Null/empty gives an empty dictionary.</param>
    /// <param name="k">The k-mer length. Must be positive for non-empty input.</param>
    /// <param name="options">Counting mode; <c>Canonical</c> must be false.</param>
    /// <returns>Dictionary mapping each k-mer seen on either strand to its summed forward + reverse-strand count.</returns>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="k"/> ≤ 0 and the sequence is non-empty.</exception>
    /// <exception cref="ArgumentException"><paramref name="options"/> has <c>Canonical = true</c>.</exception>
    public static Dictionary<string, int> CountKmersBothStrands(string sequence, int k, KmerCountingOptions options)
    {
        if (options.Canonical)
            throw new ArgumentException(
                "Canonical counting is already the both-strand count (Jellyfish count -C); use CountKmers(sequence, k, new KmerCountingOptions(Canonical: true)).",
                nameof(options));

        var forward = CountKmers(sequence, k, options);
        var combined = new Dictionary<string, int>(forward);
        foreach (var (kmer, count) in forward)
        {
            var rc = DnaSequence.GetReverseComplementString(kmer);
            CollectionsMarshal.GetValueRefOrAddDefault(combined, rc, out _) += count;
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
    /// <exception cref="ArgumentNullException"><paramref name="dna"/> is null.</exception>
    public static Dictionary<string, int> CountKmersBothStrands(DnaSequence dna, int k)
    {
        ArgumentNullException.ThrowIfNull(dna);
        return CountKmersBothStrands(dna.Sequence, k);
    }

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
/// Background-adjusted word-match statistics of <see cref="KmerAnalyzer.BackgroundAdjustedD2(string, string, int, int)"/>:
/// the similarity statistics D2* and D2S and their dissimilarities d2* and d2S (in [0, 1]; Song et al. 2014).
/// </summary>
public readonly record struct D2StarStatistics(double D2Star, double D2Shepherd, double D2StarDistance, double D2ShepherdDistance)
{
    /// <summary>Background Markov order used for the first sequence (the requested order, or the BIC choice for −1).</summary>
    public int MarkovOrder1 { get; init; }

    /// <summary>Background Markov order used for the second sequence.</summary>
    public int MarkovOrder2 { get; init; }
}

/// <summary>
/// A bottom-s MinHash sketch as built by <c>mash sketch</c> (<see cref="KmerAnalyzer.CreateMinHashSketch(IEnumerable{string}, int, int, bool, uint)"/>).
/// </summary>
/// <param name="K">K-mer size.</param>
/// <param name="SketchSize">Target sketch size s (the sketch holds min(s, distinct hashes) values).</param>
/// <param name="Canonical">True for canonical k-mers (Mash default), false for <c>-n</c>.</param>
/// <param name="Seed">MurmurHash3 seed.</param>
/// <param name="Use64">True when the 64-bit hash h1 is used (k ≥ 17); false for its low 32 bits.</param>
/// <param name="Length">Total length of the sketched records of length ≥ k (Mash reference length).</param>
/// <param name="Hashes">The smallest distinct hash values, ascending.</param>
public sealed record MinHashSketch(int K, int SketchSize, bool Canonical, uint Seed, bool Use64, long Length, IReadOnlyList<ulong> Hashes)
{
    /// <summary>
    /// Builds a well-formed sketch from arbitrary hash values (e.g. read from a Mash <c>.msh</c> file): the values are sorted,
    /// duplicates removed and the <paramref name="sketchSize"/> smallest kept, which is the form
    /// <see cref="KmerAnalyzer.CompareMinHashSketches"/> requires.
    /// </summary>
    /// <param name="k">K-mer size, 1..32.</param>
    /// <param name="sketchSize">Sketch size s ≥ 1.</param>
    /// <param name="canonical">Canonical mode of the hashed k-mers.</param>
    /// <param name="seed">Hash seed.</param>
    /// <param name="length">Total sequence length (Mash reference length), ≥ 0.</param>
    /// <param name="hashes">Hash values in any order.</param>
    /// <returns>The sketch, with <see cref="Use64"/> = k ≥ 17 as Mash sets it.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="hashes"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="k"/> outside 1..32, <paramref name="sketchSize"/> &lt; 1,
    /// or <paramref name="length"/> &lt; 0.</exception>
    /// <exception cref="ArgumentException">k ≤ 16 (32-bit sketch, Mash <c>use64 = 4^k &gt; 2^32</c>) and a kept hash exceeds
    /// 2^32 − 1.</exception>
    public static MinHashSketch FromHashes(int k, int sketchSize, bool canonical, uint seed, long length, IEnumerable<ulong> hashes)
    {
        ArgumentNullException.ThrowIfNull(hashes);
        if (k < 1 || k > KmerAnalyzer.MaxMashKmerSize)
            throw new ArgumentOutOfRangeException(nameof(k), k, $"K must be in 1..{KmerAnalyzer.MaxMashKmerSize} (mash sketch -k).");
        ArgumentOutOfRangeException.ThrowIfLessThan(sketchSize, 1);
        ArgumentOutOfRangeException.ThrowIfNegative(length);
        bool use64 = k > 16;
        var bottom = hashes.Distinct().Order().Take(sketchSize).ToArray();
        if (!use64 && bottom.Length > 0 && bottom[^1] > uint.MaxValue)
            throw new ArgumentException("For k <= 16 Mash stores 32-bit hashes; every value must be <= 2^32 - 1.", nameof(hashes));
        return new MinHashSketch(k, sketchSize, canonical, seed, use64, length, bottom);
    }
}

/// <summary>
/// A FracMinHash ("scaled") sketch as built by sourmash (<see cref="KmerAnalyzer.CreateFracMinHashSketch"/>).
/// </summary>
/// <param name="K">K-mer size.</param>
/// <param name="Scaled">Scaled factor S (1..4294967295, sourmash u32).</param>
/// <param name="MaxHash">Threshold <see cref="KmerAnalyzer.FracMinHashMaxHash"/>(S): the sketch holds every hash ≤ MaxHash.</param>
/// <param name="Canonical">True for canonical k-mers (sourmash), false for forward k-mers.</param>
/// <param name="Seed">MurmurHash3 seed.</param>
/// <param name="Hashes">The kept distinct hash values, ascending.</param>
/// <param name="Abundances">With <c>trackAbundance</c> (sourmash <c>track_abundance=True</c>): the count of each hash, parallel to
/// <paramref name="Hashes"/>; null otherwise.</param>
public sealed record FracMinHashSketch(
    int K, long Scaled, ulong MaxHash, bool Canonical, uint Seed, IReadOnlyList<ulong> Hashes, IReadOnlyList<long>? Abundances = null);

/// <summary>
/// Result of <see cref="KmerAnalyzer.CompareFracMinHashSketches"/> — sourmash <c>jaccard</c>, <c>contained_by</c> (both ways),
/// <c>max_containment</c> and, with abundances, <c>angular_similarity</c> / <c>contained_by_weighted</c>.
/// </summary>
/// <param name="SharedHashes">|A ∩ B| (sourmash <c>count_common</c>).</param>
/// <param name="UnionHashes">|A ∪ B|.</param>
/// <param name="Jaccard">Shared / max(1, union).</param>
/// <param name="ContainmentAInB">sourmash <c>A.contained_by(B)</c> (bias-corrected, clamped to [0, 1]).</param>
/// <param name="ContainmentBInA">sourmash <c>B.contained_by(A)</c>.</param>
/// <param name="MaxContainment">sourmash <c>A.max_containment(B)</c>.</param>
/// <param name="AngularSimilarity">sourmash <c>A.angular_similarity(B)</c> (1 − 2·acos(cos)/π of the abundance vectors) when
/// both sketches track abundance; null otherwise.</param>
/// <param name="WeightedContainmentAInB">sourmash <c>A.contained_by_weighted(B)</c> when A tracks abundance; null otherwise.</param>
/// <param name="WeightedContainmentBInA">sourmash <c>B.contained_by_weighted(A)</c> when B tracks abundance; null otherwise.</param>
public readonly record struct FracMinHashComparison(
    int SharedHashes, int UnionHashes, double Jaccard, double ContainmentAInB, double ContainmentBInA, double MaxContainment,
    double? AngularSimilarity = null, double? WeightedContainmentAInB = null, double? WeightedContainmentBInA = null);

/// <summary>
/// Result of <see cref="KmerAnalyzer.CompareMinHashSketches"/> — one <c>mash dist</c> output line.
/// </summary>
/// <param name="SharedHashes">x, the shared hashes (numerator of Mash's "x/s" column).</param>
/// <param name="Denominator">Union hashes visited (denominator of "x/s"; s unless both sketches are smaller).</param>
/// <param name="Jaccard">Jaccard estimate x / denominator.</param>
/// <param name="Distance">Mash distance −ln(2J/(1+J))/k (0 when x = denominator, 1 when x = 0, capped at 1).</param>
/// <param name="PValue">Probability of ≥ x shared hashes by chance (Mash binomial p-value).</param>
public readonly record struct MashComparison(int SharedHashes, int Denominator, double Jaccard, double Distance, double PValue);

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

    /// <summary>
    /// d2* = ½(1 − D2*/√(Σ X̃²/E_X · Σ Ỹ²/E_Y)) on background-centred counts X̃ = X − E_X, with D2* = Σ X̃Ỹ/√(E_X E_Y)
    /// (Reinert et al. 2009; Song et al. 2014; CAFE <c>D2star</c>), order-0 (i.i.d.) background estimated from each
    /// sequence; a dissimilarity in [0, 1]. Needs the sequences (string overload or
    /// <see cref="KmerAnalyzer.BackgroundAdjustedD2(string, string, int, int)"/> for a Markov order r).
    /// </summary>
    D2Star,

    /// <summary>
    /// d2S = ½(1 − D2S/√(Σ X̃²/√(X̃²+Ỹ²) · Σ Ỹ²/√(X̃²+Ỹ²))), D2S = Σ X̃Ỹ/√(X̃²+Ỹ²) ("D2 shepherd"; Reinert et al. 2009;
    /// Wan et al. 2010; Song et al. 2014; CAFE <c>D2shepp</c>), order-0 background estimated from each sequence; a
    /// dissimilarity in [0, 1]. Needs the sequences (see <see cref="D2Star"/>).
    /// </summary>
    D2Shepherd,

    /// <summary>
    /// Jensen–Shannon divergence ½Σ f₁ log₂(f₁/m) + ½Σ f₂ log₂(f₂/m), m = ½(f₁+f₂), on relative frequencies (Lin 1991;
    /// Leimeister et al. 2014 "JS"; <c>spaced -d JS</c>); in [0, 1], = scipy <c>jensenshannon(p, q, base=2)</c> squared.
    /// One empty table against a non-empty one gives ½ (the zero vector convention of this method).
    /// </summary>
    JensenShannon,

    /// <summary>√Σ(c₁−c₂)² on raw counts — the per-pattern Euclidean value of the <c>spaced</c> 1.2 program (<c>-d EU</c>).</summary>
    EuclideanCounts,

    /// <summary>
    /// Evolutionary distance (substitutions per site, Jukes–Cantor corrected) estimated from the number of spaced-word
    /// matches (Morgenstern, Zhu, Horwege &amp; Leimeister 2015, Algorithms Mol Biol 10:5; <c>spaced -d EV</c>). Not a
    /// word-vector metric: it needs the sequence lengths and base composition, so only the string overloads accept it
    /// (<see cref="KmerAnalyzer.SpacedWordDistance(string, string, IReadOnlyList{string}, KmerDistanceMetric, KmerCountingOptions, bool)"/>,
    /// and <see cref="KmerAnalyzer.KmerDistance(string, string, int, KmerDistanceMetric, int, bool)"/> with the pattern 1^k).
    /// </summary>
    SpacedEvolutionary,
}
