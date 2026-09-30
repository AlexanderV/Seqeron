namespace Seqeron.Genomics.Analysis;

/// <summary>
/// Finds various types of repeats in DNA sequences including microsatellites (STRs),
/// minisatellites (VNTRs), inverted repeats, and direct repeats.
/// </summary>
public static class RepeatFinder
{
    #region Microsatellite (STR) Detection

    /// <summary>
    /// Finds microsatellites (Short Tandem Repeats) in a DNA sequence.
    /// STRs are 1-6 bp motifs repeated consecutively.
    /// </summary>
    /// <remarks>
    /// <para>
    /// For every unit length <c>p</c> in <c>[minUnitLength, maxUnitLength]</c> the detector reports each
    /// <b>maximal perfect run</b> of period <c>p</c> exactly once: a run is the longest interval
    /// <c>S[a..e)</c> with <c>S[x] = S[x+p]</c> for all <c>a ≤ x &lt; e−p</c> that cannot be extended to the
    /// left or right (a maximal repetition, Kolpakov &amp; Kucherov 1999). The run is reported at its
    /// left end <c>a</c> with <c>RepeatUnit = S[a..a+p)</c> (the motif phase at the run start) and
    /// <c>RepeatCount = ⌊(e−a)/p⌋</c> complete copies; a trailing partial copy (&lt; p bases) is not
    /// counted, and rotations of the same run (e.g. <c>TA</c> inside <c>ATATATA</c>) are not reported again.
    /// This is the per-motif-size convention of MISA (Thiel et al. 2003; leftmost regex match
    /// <c>([acgt]{p})\1{k-1,}</c>) and pytrf/Krait (Du et al. 2018: seed at the run start, <c>repeat = length / p</c>).
    /// </para>
    /// <para>
    /// A unit that is itself a repetition of a shorter unit (e.g. <c>ATAT</c>, <c>AA</c>) is not reported —
    /// its run is reported at the primitive unit length (MISA "reject false type motifs"). Only units made of
    /// the unambiguous bases A/C/G/T are reported (MISA <c>[acgt]</c>; pytrf skips <c>N</c>), so runs of
    /// <c>N</c> (assembly gaps) or other symbols are never microsatellites. Runs of different unit lengths may
    /// overlap and are reported independently. Input is case-insensitive; positions are 0-based.
    /// </para>
    /// </remarks>
    /// <param name="sequence">DNA sequence to search.</param>
    /// <param name="minUnitLength">Minimum repeat unit length (default: 1).</param>
    /// <param name="maxUnitLength">Maximum repeat unit length (default: 6).</param>
    /// <param name="minRepeats">Minimum number of complete consecutive copies to report (default: 3).</param>
    /// <returns>Collection of microsatellite repeats found, ordered by unit length then position.</returns>
    public static IEnumerable<MicrosatelliteResult> FindMicrosatellites(
        DnaSequence sequence,
        int minUnitLength = 1,
        int maxUnitLength = 6,
        int minRepeats = 3)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ValidateMicrosatelliteParameters(minUnitLength, maxUnitLength, minRepeats);

        return FindMicrosatellitesCore(
            sequence.Sequence, minUnitLength, maxUnitLength, minRepeats, CancellationToken.None, null);
    }

    /// <summary>
    /// Finds microsatellites with cancellation support. Same semantics as
    /// <see cref="FindMicrosatellites(DnaSequence,int,int,int)"/>.
    /// </summary>
    /// <param name="sequence">DNA sequence to search.</param>
    /// <param name="minUnitLength">Minimum repeat unit length.</param>
    /// <param name="maxUnitLength">Maximum repeat unit length.</param>
    /// <param name="minRepeats">Minimum number of repeats.</param>
    /// <param name="cancellationToken">Cancellation token for long-running operations.</param>
    /// <param name="progress">Optional progress reporter (0.0 to 1.0).</param>
    /// <returns>Collection of microsatellite repeats found.</returns>
    public static IEnumerable<MicrosatelliteResult> FindMicrosatellites(
        DnaSequence sequence,
        int minUnitLength,
        int maxUnitLength,
        int minRepeats,
        CancellationToken cancellationToken,
        IProgress<double>? progress = null)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ValidateMicrosatelliteParameters(minUnitLength, maxUnitLength, minRepeats);

        return FindMicrosatellitesCore(
            sequence.Sequence, minUnitLength, maxUnitLength, minRepeats, cancellationToken, progress);
    }

    /// <summary>
    /// Finds microsatellites in a raw sequence string (case-insensitive). Same semantics as
    /// <see cref="FindMicrosatellites(DnaSequence,int,int,int)"/>; <c>null</c>/empty input yields no results.
    /// </summary>
    public static IEnumerable<MicrosatelliteResult> FindMicrosatellites(
        string sequence,
        int minUnitLength = 1,
        int maxUnitLength = 6,
        int minRepeats = 3)
    {
        ValidateMicrosatelliteParameters(minUnitLength, maxUnitLength, minRepeats);

        if (string.IsNullOrEmpty(sequence))
            return [];

        return FindMicrosatellitesCore(
            sequence.ToUpperInvariant(), minUnitLength, maxUnitLength, minRepeats, CancellationToken.None, null);
    }

    /// <summary>
    /// Finds microsatellites in a raw sequence string with cancellation support.
    /// </summary>
    public static IEnumerable<MicrosatelliteResult> FindMicrosatellites(
        string sequence,
        int minUnitLength,
        int maxUnitLength,
        int minRepeats,
        CancellationToken cancellationToken,
        IProgress<double>? progress = null)
    {
        ValidateMicrosatelliteParameters(minUnitLength, maxUnitLength, minRepeats);

        if (string.IsNullOrEmpty(sequence))
            return [];

        return FindMicrosatellitesCore(
            sequence.ToUpperInvariant(), minUnitLength, maxUnitLength, minRepeats, cancellationToken, progress);
    }

    private static void ValidateMicrosatelliteParameters(int minUnitLength, int maxUnitLength, int minRepeats)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(minUnitLength, 1);
        ArgumentOutOfRangeException.ThrowIfLessThan(maxUnitLength, minUnitLength);
        ArgumentOutOfRangeException.ThrowIfLessThan(minRepeats, 2);
    }

    /// <summary>
    /// Single scan core (upper-case input): for each unit length, visits only run starts
    /// (left-maximal positions), extends the run to its right-maximal end and reports it once.
    /// O(n) character comparisons per unit length.
    /// </summary>
    private static IEnumerable<MicrosatelliteResult> FindMicrosatellitesCore(
        string seq,
        int minUnitLength,
        int maxUnitLength,
        int minRepeats,
        CancellationToken cancellationToken,
        IProgress<double>? progress)
    {
        int n = seq.Length;
        long unitLengths = (long)maxUnitLength - minUnitLength + 1;
        double totalPositions = Math.Max(1.0, (double)n * unitLengths);
        int sinceCheck = 0;
        const int checkInterval = 1000;

        for (int unitLen = minUnitLength; unitLen <= maxUnitLength; unitLen++)
        {
            if ((long)unitLen * minRepeats > n)
                break; // no longer unit length can fit minRepeats copies either

            int i = 0;
            while (i + unitLen * minRepeats <= n)
            {
                if (++sinceCheck >= checkInterval)
                {
                    sinceCheck = 0;
                    cancellationToken.ThrowIfCancellationRequested();
                    progress?.Report(((double)(unitLen - minUnitLength) * n + i) / totalPositions);
                }

                // Run start (left-maximal): the period-p run cannot be extended one base to the left.
                // Positions inside a run (rotations / suffixes of it) are skipped.
                if (i > 0 && seq[i - 1] == seq[i - 1 + unitLen])
                {
                    i++;
                    continue;
                }

                // Right-maximal end e (exclusive): extend while S[x] == S[x - p].
                int e = i + unitLen;
                while (e < n && seq[e] == seq[e - unitLen])
                    e++;

                int repeats = (e - i) / unitLen;
                if (repeats >= minRepeats && IsReportableUnit(seq, i, unitLen))
                {
                    string unit = seq.Substring(i, unitLen);
                    yield return new MicrosatelliteResult(
                        Position: i,
                        RepeatUnit: unit,
                        RepeatCount: repeats,
                        TotalLength: repeats * unitLen,
                        RepeatType: ClassifyRepeatType(unit));
                }

                // Every position in (i, e - p] lies inside this run (not left-maximal); the next possible
                // run start of the same period is e - p + 1.
                i = e - unitLen + 1;
            }
        }

        cancellationToken.ThrowIfCancellationRequested();
        progress?.Report(1.0);
    }

    /// <summary>
    /// A unit is reportable when it consists only of A/C/G/T and is primitive (not a power of a shorter word).
    /// </summary>
    private static bool IsReportableUnit(string seq, int start, int unitLen)
    {
        for (int k = start; k < start + unitLen; k++)
        {
            if (AcgtCode(seq[k]) < 0)
                return false;
        }

        return !IsRedundantUnit(seq.Substring(start, unitLen));
    }

    private static bool IsRedundantUnit(string unit)
    {
        if (unit.Length <= 1) return false;

        // Check if unit is made of smaller repeating pattern
        for (int subLen = 1; subLen < unit.Length; subLen++)
        {
            if (unit.Length % subLen != 0) continue;

            string subUnit = unit.Substring(0, subLen);
            bool isRedundant = true;

            for (int i = subLen; i < unit.Length; i += subLen)
            {
                if (unit.Substring(i, subLen) != subUnit)
                {
                    isRedundant = false;
                    break;
                }
            }

            if (isRedundant) return true;
        }

        return false;
    }

    private static RepeatType ClassifyRepeatType(string unit)
    {
        return unit.Length switch
        {
            1 => RepeatType.Mononucleotide,
            2 => RepeatType.Dinucleotide,
            3 => RepeatType.Trinucleotide,
            4 => RepeatType.Tetranucleotide,
            5 => RepeatType.Pentanucleotide,
            6 => RepeatType.Hexanucleotide,
            _ => RepeatType.Complex
        };
    }

    #endregion

    #region Approximate (Imperfect/Interrupted) Tandem Repeat Detection — TRF model

    // --- Tandem Repeats Finder (TRF; Benson 1999) ----------------------------------------------------
    // Benson G (1999) "Tandem repeats finder: a program to analyze DNA sequences", Nucleic Acids Res
    // 27(2):573-580, https://doi.org/10.1093/nar/27.2.573; TRF 4.10.0 README (github.com/
    // Benson-Genomics-Lab/TRF: parameters, "TRF Definitions", "How does Tandem Repeats Finder work?").
    //
    // TRF has a DETECTION component (k-tuple matches at a common distance d, tested against statistical
    // criteria) and an ANALYSIS component (wraparound dynamic programming (WDP) of the sequence against
    // tandem copies of a candidate pattern, majority-rule consensus, realignment against the consensus,
    // statistics "between adjacent copies"). This region implements:
    //   * the ANALYSIS component in full (reported score / indices / period / copy number / consensus /
    //     %matches / %indels / composition identical to compiled TRF 4.10.0 on 1524/1524 analysed
    //     candidates with pattern <= 20 bp, the range where TRF itself runs the full WDP (SMALLDISTANCE);
    //     for larger patterns TRF restricts WDP to a narrow diagonal band — a speed heuristic — whereas
    //     this code keeps the full (optimal) WDP: 827/959 identical);
    //   * the DETECTION component's k-tuple trigger with Benson's tuple sizes (Table 1 / README: k = 4 for
    //     d <= 29, 5 for 30..159, 7 for >= 160 at PM = .80) and the sum-of-heads criterion R(d,k,PM),
    //     derived here from the exact mean/variance of R (normal approximation, 95% one-sided, floor
    //     k+1) — this reproduces TRF's own sumdata80 table for all d = 1..2000;
    //   * the three-best-periods ("multiples") test, the minimum copy-number rule, per-distance
    //     "already aligned" suppression, and TRF's redundancy elimination / MaxPeriod filter.
    // NOT implemented (declared residual, see FindApproximateTandemRepeats remarks): the apparent-size
    // (waiting-time) criterion, whose cut-offs TRF estimates by simulation; the random-walk distance
    // range d +/- floor(2.3*sqrt(PI*d)) summation; the narrow-band WDP for patterns > 20; and TRF's
    // best-period list for d > 250. A line-by-line port of TRF is also excluded by licence: TRF is
    // AGPL-3.0, this library is MIT.

    /// <summary>TRF match weight: "Match ... The recomended values for Match Mismatch and Delta are 2, 7, and 7" (TRF README).</summary>
    private const int TrfMatchWeight = 2;

    /// <summary>TRF mismatch weight −7 (TRF README; weights are "interpreted as negative numbers").</summary>
    private const int TrfMismatchWeight = -7;

    /// <summary>TRF indel weight (Delta) −7 per gap column (TRF README).</summary>
    private const int TrfIndelWeight = -7;

    /// <summary>Marks a dead WDP cell (TRF: a zero cell beyond the candidate is set to −1000 so that the
    /// local alignment cannot restart there).</summary>
    private const int TrfDeadCell = -1000;

    /// <summary>
    /// TRF Min_Distance_Window = 20: the smallest tandem-repeat span the detector wants to see; the k-tuple
    /// distance window and the backward WDP scan both extend at least this far (TRF 4.10.0 source, tr30dat.h).
    /// </summary>
    private const int TrfMinDistanceWindow = 20;

    /// <summary>
    /// Default minimum alignment score to report a tandem repeat. Benson (1999): "Only those repeats
    /// scoring at least 50 ... are reported"; the recommended TRF parameter set uses Minscore = 50.
    /// </summary>
    public const int DefaultApproximateMinScore = 50;

    /// <summary>
    /// Largest supported period. TRF README (4.10.0): "TRF will throw an error if a value of over 2000 is
    /// given for MaxPeriod ... very large TRs are outside the scope of the TRF statistical models."
    /// </summary>
    public const int MaxApproximatePeriod = 2000;

    /// <summary>
    /// Finds approximate (imperfect / interrupted) tandem repeats with the Tandem Repeats Finder model
    /// (Benson 1999; recommended parameters Match 2, Mismatch 7, Delta 7, PM 80, PI 10).
    /// </summary>
    /// <remarks>
    /// <para><b>Analysis (TRF-exact).</b> For every candidate (position i, distance d) the sequence is aligned
    /// by wraparound dynamic programming (local alignment against unlimited tandem copies of the pattern
    /// S[i−d+1..i]); a consensus is taken by majority rule from that alignment and the sequence is realigned
    /// against the consensus. Reported values follow the TRF table: indices, period = most common distance
    /// between matching characters of adjacent copies, copy number = aligned consensus columns / consensus
    /// size, consensus size, % matches and % indels "between adjacent copies overall" (not between the
    /// sequence and the consensus), alignment score, nucleotide composition and entropy. A repeat needs at
    /// least 1.9 copies (1.8 for large patterns) and score ≥ <paramref name="minScore"/>. Scoring: +2 for
    /// an identical A/C/G/T pair, −7 for any other pair (N and other symbols never match), −7 per gap.</para>
    /// <para><b>Detection (TRF criteria, partial).</b> A candidate is examined when a k-tuple match at distance d
    /// ends at i and the heads counted in k-runs over the last max(d, 20) positions reach the sum-of-heads
    /// cut-off; d must be among the three best periods of the aligned region (period 1 needs ≥ 80 % of one
    /// base); overlapping reports are reduced with TRF's redundancy rule (≥ 90 % overlap, same period or a
    /// multiple scoring ≤ 1.1×). TRF's simulated apparent-size criterion, random-walk distance ranges,
    /// narrow-band alignment for patterns &gt; 20 and best-period list are not reproduced, so the set of
    /// reported loci can differ from TRF (measured on random sequences with embedded repeats: 92.6 % of TRF
    /// rows identical and 96 % found at region level for periods ≤ 20; 80.5 % / 93 % for periods ≤ 100 —
    /// docs/Evidence/REP-APPROX-001-Evidence.md).</para>
    /// <para>Complexity: O(n · maxPeriod) for the k-tuple scan plus one wraparound DP over the aligned region per
    /// examined candidate (O(region · pattern) time, O(pattern) memory unless the candidate passes the copy and
    /// best-period tests, which store the traceback matrix).</para>
    /// <para>Coordinates are 0-based (<c>Start</c>; TRF prints 1-based indices). Percentages are exact
    /// (TRF truncates them to integers). Output is ordered by start, then end, then period.</para>
    /// </remarks>
    /// <param name="sequence">DNA sequence to search.</param>
    /// <param name="minPeriod">Minimum reported period (≥ 1; default 1).</param>
    /// <param name="maxPeriod">Maximum candidate distance and reported period (minPeriod..2000; default 6).</param>
    /// <param name="minScore">Minimum TRF alignment score to report (≥ 1; default 50, Benson 1999).</param>
    /// <returns>Approximate tandem repeats ordered by start position.</returns>
    public static IEnumerable<ApproximateTandemRepeatResult> FindApproximateTandemRepeats(
        DnaSequence sequence,
        int minPeriod = 1,
        int maxPeriod = 6,
        int minScore = DefaultApproximateMinScore)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ValidateApproximateParameters(minPeriod, maxPeriod, minScore);
        return FindApproximateTandemRepeatsCore(sequence.Sequence, minPeriod, maxPeriod, minScore);
    }

    /// <summary>
    /// Finds approximate (imperfect / interrupted) tandem repeats in a raw sequence string with the Tandem
    /// Repeats Finder model; case-insensitive, any non-A/C/G/T symbol never matches. See
    /// <see cref="FindApproximateTandemRepeats(DnaSequence,int,int,int)"/>. Null or empty input yields no repeats.
    /// </summary>
    public static IEnumerable<ApproximateTandemRepeatResult> FindApproximateTandemRepeats(
        string sequence,
        int minPeriod = 1,
        int maxPeriod = 6,
        int minScore = DefaultApproximateMinScore)
    {
        ValidateApproximateParameters(minPeriod, maxPeriod, minScore);
        if (string.IsNullOrEmpty(sequence))
            return Array.Empty<ApproximateTandemRepeatResult>();

        return FindApproximateTandemRepeatsCore(sequence.ToUpperInvariant(), minPeriod, maxPeriod, minScore);
    }

    private static void ValidateApproximateParameters(int minPeriod, int maxPeriod, int minScore)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(minPeriod, 1);
        ArgumentOutOfRangeException.ThrowIfLessThan(maxPeriod, minPeriod);
        ArgumentOutOfRangeException.ThrowIfGreaterThan(maxPeriod, MaxApproximatePeriod);
        ArgumentOutOfRangeException.ThrowIfLessThan(minScore, 1);
    }

    /// <summary>One column of a TRF alignment: sequence symbol (or '-'), pattern symbol (or '-'), 1-based
    /// sequence index, 0-based pattern index (for an inserted sequence symbol: the next pattern position).</summary>
    private readonly record struct TrfColumn(char Seq, char Pat, int SeqIndex, int PatIndex);

    /// <summary>A traced WDP alignment; <see cref="Columns"/> run from the RIGHTMOST column to the leftmost.</summary>
    private sealed record TrfAlignment(int Score, double CopyNumber, TrfColumn[] Columns)
    {
        public int First => Columns[^1].SeqIndex;
        public int Last => Columns[0].SeqIndex;
    }

    /// <summary>Extent of a WDP alignment without traceback (same optimum and same path as the traceback).</summary>
    private readonly record struct TrfExtent(int Score, int First, int Last, double CopyNumber);

    private static int TrfWeight(char a, char b) =>
        a == b && AcgtCode(a) >= 0 ? TrfMatchWeight : TrfMismatchWeight;

    /// <summary>TRF k-tuple size for distance d at PM = .80 (Benson 1999 Table 1; TRF 4.10.0: 4 / 5 / 7).</summary>
    private static int TrfTupleSize(int d)
    {
        if (d <= 29)
            return 4;
        return d <= 159 ? 5 : 7;
    }

    private static readonly int[] SumOfHeadsCache = new int[MaxApproximatePeriod + 1];

    /// <summary>
    /// Sum-of-heads criterion for distance d (Benson 1999): R(d,k,PM) = total heads in head runs of length ≥ k
    /// in an iid Bernoulli(PM) sequence of length d; "the distribution of R is well approximated by the normal
    /// distribution and its exact mean and variance can be calculated"; the criterion is the largest x such
    /// that R ≥ x 95 % of the time. Mean and variance are computed exactly by a run-length Markov chain; the
    /// cut-off is ⌊μ − 1.65σ⌋, never below k + 1 ("the smallest pattern for tuple size k [has] a sum-of-heads
    /// criterion of at least k+1"). Reproduces TRF 4.10.0's PM = 80 table for every d = 1..2000.
    /// </summary>
    internal static int TrfSumOfHeadsCriterion(int d)
    {
        int cached = Volatile.Read(ref SumOfHeadsCache[d]);
        if (cached != 0)
            return cached;

        int k = TrfTupleSize(d);
        const double pm = TrfDefaultMatchProbability;
        // State r = current head-run length (0..k-1) or k = inside a run already counted.
        var p = new double[k + 1];
        var m1 = new double[k + 1];
        var m2 = new double[k + 1];
        var np = new double[k + 1];
        var n1 = new double[k + 1];
        var n2 = new double[k + 1];
        p[0] = 1.0;
        for (int step = 0; step < d; step++)
        {
            Array.Clear(np);
            Array.Clear(n1);
            Array.Clear(n2);
            for (int r = 0; r <= k; r++)
            {
                if (p[r] == 0.0)
                    continue;
                // tails: run resets, nothing added
                np[0] += p[r] * (1 - pm);
                n1[0] += m1[r] * (1 - pm);
                n2[0] += m2[r] * (1 - pm);
                // heads: the k-th head of a run adds k, every later head adds 1
                int next = r < k ? r + 1 : k;
                int add = 1;
                if (r < k - 1)
                    add = 0;
                else if (r == k - 1)
                    add = k;
                np[next] += p[r] * pm;
                n1[next] += (m1[r] + add * p[r]) * pm;
                n2[next] += (m2[r] + 2.0 * add * m1[r] + (double)add * add * p[r]) * pm;
            }
            (p, np) = (np, p);
            (m1, n1) = (n1, m1);
            (m2, n2) = (n2, m2);
        }

        double mean = m1.Sum();
        double sd = Math.Sqrt(Math.Max(0.0, m2.Sum() - mean * mean));
        int criterion = Math.Max(k + 1, (int)(mean - 1.65 * sd));
        Volatile.Write(ref SumOfHeadsCache[d], criterion);
        return criterion;
    }

    private static IReadOnlyList<ApproximateTandemRepeatResult> FindApproximateTandemRepeatsCore(
        string sequence,
        int minPeriod,
        int maxPeriod,
        int minScore)
    {
        int n = sequence.Length;
        var s = new char[n + 1]; // 1-based, as in TRF
        sequence.CopyTo(0, s, 1, n);

        var found = new List<ApproximateTandemRepeatResult>();
        var runLength = new int[maxPeriod + 1];
        var seenEnd = new int[maxPeriod + 1];
        var windows = new TupleMatchWindow?[maxPeriod + 1];
        var bestPeriods = new Dictionary<(int First, int Last), int[]>();

        for (int i = 2; i <= n; i++)
        {
            bool acgt = AcgtCode(s[i]) >= 0;
            int dMax = Math.Min(maxPeriod, i - 1);
            for (int d = 1; d <= dMax; d++)
            {
                if (!acgt || s[i] != s[i - d])
                {
                    runLength[d] = 0;
                    continue;
                }

                int k = TrfTupleSize(d);
                if (++runLength[d] < k)
                    continue;

                // A k-tuple match at distance d ends at i; record it in the distance window of the last
                // max(d, 20) positions (adjacent tuple matches form one entry of growing size).
                var window = windows[d] ??= new TupleMatchWindow(Math.Max(d, TrfMinDistanceWindow) + 1);
                window.Add(i, k, i - Math.Max(d, TrfMinDistanceWindow) + 1);

                if (seenEnd[d] >= i || window.Heads < TrfSumOfHeadsCriterion(d))
                    continue;

                var repeat = AnalyzeTrfCandidate(s, n, i, d, maxPeriod, minScore, seenEnd, bestPeriods);
                if (repeat is not null)
                    found.Add(repeat.Value);
            }
        }

        // TRF: drop periods above MaxPeriod, sort by start (stable), eliminate redundancy. The minimum
        // period is a library option applied afterwards, so it never resurrects a redundant multiple.
        var reported = RemoveTrfRedundancy(found.Where(r => r.Period <= maxPeriod).OrderBy(r => r.Start).ToList());

        return reported
            .Where(r => r.Period >= minPeriod)
            .OrderBy(r => r.Start)
            .ThenBy(r => r.Start + r.SpanLength)
            .ThenBy(r => r.Period)
            .ToList();
    }

    /// <summary>
    /// Sliding window of k-tuple matches at one distance d (TRF distance list): runs of adjacent tuple
    /// matches, each stored as (end position, number of heads); runs whose end falls before the window's
    /// left end are dropped. <see cref="Heads"/> is the sum of heads in k-runs inside the window.
    /// </summary>
    private sealed class TupleMatchWindow(int capacity)
    {
        private readonly int[] _end = new int[capacity];
        private readonly int[] _size = new int[capacity];
        private int _head;
        private int _count;

        public int Heads { get; private set; }

        public void Add(int position, int tupleSize, int windowLeft)
        {
            while (_count > 0 && _end[_head] < windowLeft)
            {
                Heads -= _size[_head];
                _head = (_head + 1) % _end.Length;
                _count--;
            }

            int tail = (_head + _count - 1 + _end.Length) % _end.Length;
            if (_count > 0 && _end[tail] == position - 1)
            {
                _end[tail] = position;
                _size[tail]++;
                Heads++;
                return;
            }

            tail = (_head + _count) % _end.Length;
            _end[tail] = position;
            _size[tail] = tupleSize;
            _count++;
            Heads += tupleSize;
        }
    }

    /// <summary>
    /// TRF analysis of one candidate (position <paramref name="i"/>, distance <paramref name="d"/>): align the
    /// candidate pattern S[i−d+1..i], apply the copy-number and three-best-periods tests, build the consensus,
    /// realign against it and report the statistics. Updates <paramref name="seenEnd"/>[d] with the end of each
    /// alignment so that later matches at d inside the aligned region are not re-analysed (TRF).
    /// </summary>
    private static ApproximateTandemRepeatResult? AnalyzeTrfCandidate(
        char[] s, int n, int i, int d, int maxPeriod, int minScore, int[] seenEnd,
        Dictionary<(int First, int Last), int[]> bestPeriods)
    {
        var pattern = new char[d];
        Array.Copy(s, i - d + 1, pattern, 0, d);

        // Extent-only pass first: the multiples test needs only the aligned region, so the traceback
        // matrix is built only for candidates that pass both tests (identical path either way).
        var extent = TrfWraparoundExtent(s, n, i, pattern);
        if (extent is null)
            return null;
        MarkAligned(seenEnd, d, extent.Value.Last);
        if (!MeetsTrfCopyNumber(extent.Value.CopyNumber, d, d) ||
            !IsAmongTrfBestPeriods(s, extent.Value.First, extent.Value.Last, d, maxPeriod, bestPeriods))
            return null;

        var first = TrfWraparoundAlign(s, n, i, pattern);
        if (first is null)
            return null;
        char[] consensus = TrfConsensus(first.Columns, d);
        if (consensus.Length == 0)
            return null;

        var final = TrfWraparoundAlign(s, n, i, consensus);
        if (final is null)
            return null;
        MarkAligned(seenEnd, d, final.Last);
        // TRF quirk kept: the 50 < size <= 100 ramp uses the candidate distance d.
        if (!MeetsTrfCopyNumber(final.CopyNumber, consensus.Length, d) || final.Score < minScore)
            return null;

        return TrfStatistics(s, final, consensus);
    }

    /// <summary>Records that sequence positions up to <paramref name="last"/> were aligned at distance d.</summary>
    private static void MarkAligned(int[] seenEnd, int d, int last) => seenEnd[d] = last;

    /// <summary>
    /// TRF minimum copy number: ≥ 1.9 copies for patterns ≤ 50, ramping from 1.9 down to 1.8 for 50..100,
    /// ≥ 1.8 above 100 (TRF 4.10.0; Benson 1999: "If at least two copies of the pattern are aligned with the
    /// sequence, the tandem repeat is reported").
    /// </summary>
    private static bool MeetsTrfCopyNumber(double copies, int size, int d)
    {
        if (size <= 50)
            return copies >= 1.9;
        if (size <= 100)
            return copies >= 1.9 - 0.002 * (d - 50);
        return copies >= 1.8;
    }

    /// <summary>Result of a WDP fill: best score, its cell, and (extent mode) the traced path's start and
    /// consumed pattern columns.</summary>
    private readonly record struct TrfFill(int Score, int RealRow, int Row, int Col, int First, int Consumed);

    /// <summary>
    /// Runs the TRF wraparound DP for <paramref name="pattern"/> around candidate end <paramref name="start"/>:
    /// a backward local scan locates the leftmost row reaching the best score, then a forward local alignment
    /// starting one pattern length before it yields the optimum (first strictly greatest cell, row-major).
    /// Each row is computed in two passes because the horizontal dependency wraps from the last pattern column
    /// to the first; zero cells beyond the candidate are killed so the alignment cannot restart there.
    /// With <paramref name="rows"/> the final forward rows are stored for traceback; without, every cell carries
    /// the start row and consumed pattern columns of the path the traceback would follow (same predecessor
    /// preference — diagonal, then vertical, then horizontal — on the final values), so the extent is known
    /// without the O(rows × pattern) matrix.
    /// </summary>
    private static TrfFill TrfWraparoundFill(char[] s, int n, int start, char[] pattern, List<int[]>? rows)
    {
        int size = pattern.Length;
        var weights = TrfWeightRows(pattern);
        var up = new int[size];
        var diag = new int[size];
        var cur = new int[size];

        // Backward scan (sequence read right-to-left from the candidate end, pattern read in reverse).
        Array.Fill(up, TrfIndelWeight);
        int maxScore = 0;
        int minRow = start;
        int killBelow = start - Math.Max(size, TrfMinDistanceWindow);
        int realRow = start + 1;
        bool endOfTrace = false;
        while (!endOfTrace && realRow > 1)
        {
            realRow--;
            int[] w = weights[TrfSymbolClass(s[realRow])];
            int left = TrfIndelWeight;
            for (int c = size - 1; c >= 0; c--)
            {
                diag[c] += w[c];
                left = Math.Max(Math.Max(0, diag[c]), Math.Max(up[c], left)) + TrfIndelWeight;
            }

            endOfTrace = true;
            for (int c = size - 1; c >= 0; c--)
            {
                int v = Math.Max(Math.Max(0, diag[c]), Math.Max(up[c], left));
                left = up[c] = v + TrfIndelWeight;
                if (realRow <= killBelow && v == 0)
                {
                    v = TrfDeadCell;
                    left = up[c] = TrfDeadCell;
                }
                else
                {
                    endOfTrace = false;
                }

                cur[c] = v;
                if (v >= maxScore)
                {
                    maxScore = v;
                    minRow = realRow;
                }
            }

            for (int c = 0; c < size - 1; c++)
                diag[c] = cur[c + 1];
            diag[size - 1] = cur[0];
        }

        // Forward local alignment from one pattern length before the leftmost best row.
        Array.Fill(up, TrfIndelWeight);
        Array.Clear(diag);
        bool trackExtent = rows is null;
        var prev = new int[size];
        var prevFirst = trackExtent ? new int[size] : [];
        var prevCols = trackExtent ? new int[size] : [];
        var curFirst = trackExtent ? new int[size] : [];
        var curCols = trackExtent ? new int[size] : [];
        var pending = trackExtent ? new bool[size] : [];
        rows?.Add(new int[size]);

        realRow = Math.Max(minRow - size - 1, 0);
        int row = 0;
        maxScore = 0;
        int bestReal = -1, bestRow = -1, bestCol = -1, bestFirst = 0, bestConsumed = 0;
        endOfTrace = false;
        while (!endOfTrace && realRow < n)
        {
            row++;
            realRow++;
            int[] w = weights[TrfSymbolClass(s[realRow])];
            int[] values = rows is null ? cur : new int[size];
            int left = TrfIndelWeight;
            for (int c = 0; c < size; c++)
            {
                diag[c] += w[c];
                left = Math.Max(Math.Max(0, diag[c]), Math.Max(up[c], left)) + TrfIndelWeight;
            }

            endOfTrace = true;
            int rowMax = 0, rowMaxCol = -1;
            for (int c = 0; c < size; c++)
            {
                int v = Math.Max(Math.Max(0, diag[c]), Math.Max(up[c], left));
                left = up[c] = v + TrfIndelWeight;
                if (realRow > start && v == 0)
                {
                    v = TrfDeadCell;
                    left = up[c] = TrfDeadCell;
                }
                else
                {
                    endOfTrace = false;
                }

                values[c] = v;
                if (v > maxScore && v > rowMax)
                {
                    rowMax = v;
                    rowMaxCol = c;
                }
            }

            if (trackExtent)
                TrackExtent(values, prev, w, realRow, prevFirst, prevCols, curFirst, curCols, pending);

            if (rowMaxCol >= 0)
            {
                maxScore = rowMax;
                bestReal = realRow;
                bestRow = row;
                bestCol = rowMaxCol;
                if (trackExtent)
                {
                    bestFirst = curFirst[rowMaxCol];
                    bestConsumed = curCols[rowMaxCol];
                }
            }

            for (int c = size - 1; c > 0; c--)
                diag[c] = values[c - 1];
            diag[0] = values[size - 1];

            if (rows is null)
            {
                Array.Copy(values, prev, size);
                Array.Copy(curFirst, prevFirst, size);
                Array.Copy(curCols, prevCols, size);
            }
            else
            {
                rows.Add(values);
            }
        }

        return new TrfFill(maxScore, bestReal, bestRow, bestCol, bestFirst, bestConsumed);
    }

    /// <summary>Per-row extent bookkeeping for <see cref="TrfWraparoundFill"/> (start row and consumed pattern
    /// columns inherited along the traceback's preferred predecessor).</summary>
    private static void TrackExtent(
        int[] values, int[] prev, int[] w, int realRow,
        int[] prevFirst, int[] prevCols, int[] curFirst, int[] curCols, bool[] pending)
    {
        int size = values.Length;
        for (int c = 0; c < size; c++)
        {
            pending[c] = false;
            int v = values[c];
            if (v <= 0)
                continue;
            int jp = c == 0 ? size - 1 : c - 1;
            if (v == prev[jp] + w[c])
            {
                bool starts = prev[jp] <= 0;
                curFirst[c] = starts ? realRow : prevFirst[jp];
                curCols[c] = (starts ? 0 : prevCols[jp]) + 1;
            }
            else if (v == prev[c] + TrfIndelWeight)
            {
                curFirst[c] = prevFirst[c];
                curCols[c] = prevCols[c];
            }
            else
            {
                pending[c] = true; // horizontal move: inherits from its left neighbour in this row
            }
        }

        // Horizontal chains may wrap from the last pattern column to the first: two sweeps resolve them.
        for (int sweep = 0; sweep < 2; sweep++)
        {
            for (int c = 0; c < size; c++)
            {
                int jp = c == 0 ? size - 1 : c - 1;
                if (pending[c] && !pending[jp] && values[jp] > 0)
                {
                    curFirst[c] = curFirst[jp];
                    curCols[c] = curCols[jp] + 1;
                    pending[c] = false;
                }
            }
        }
    }

    /// <summary>Weight of each pattern column against A, C, G, T and any other symbol (TRF: +2 for an
    /// identical A/C/G/T pair, −7 otherwise; N never matches).</summary>
    private static int[][] TrfWeightRows(char[] pattern)
    {
        var weights = new int[5][];
        for (int symbol = 0; symbol < 5; symbol++)
        {
            weights[symbol] = new int[pattern.Length];
            for (int c = 0; c < pattern.Length; c++)
                weights[symbol][c] = symbol < 4 && TrfSymbolClass(pattern[c]) == symbol ? TrfMatchWeight : TrfMismatchWeight;
        }
        return weights;
    }

    /// <summary>A/C/G/T → 0..3, any other symbol → 4 (never matches).</summary>
    private static int TrfSymbolClass(char c)
    {
        int code = AcgtCode(c);
        return code < 0 ? 4 : code;
    }

    /// <summary>WDP optimum and its extent (first/last sequence index, copy number) without traceback.</summary>
    private static TrfExtent? TrfWraparoundExtent(char[] s, int n, int start, char[] pattern)
    {
        var fill = TrfWraparoundFill(s, n, start, pattern, rows: null);
        if (fill.Score <= 0)
            return null;
        return new TrfExtent(fill.Score, fill.First, fill.RealRow, TrfCopyNumber(fill.Consumed, pattern.Length));
    }

    /// <summary>TRF copy number: full passes through the pattern plus the final partial pass, i.e. aligned
    /// pattern columns / pattern size (summed as TRF does, integer part then fraction).</summary>
    private static double TrfCopyNumber(int consumedColumns, int size) =>
        consumedColumns / size + (double)(consumedColumns % size) / size;

    /// <summary>
    /// Runs the TRF wraparound DP and traces the optimal local alignment back from its best cell
    /// (predecessor preference: match/mismatch, then a sequence symbol against a gap, then a pattern symbol
    /// against a gap). Columns are returned rightmost first.
    /// </summary>
    private static TrfAlignment? TrfWraparoundAlign(char[] s, int n, int start, char[] pattern)
    {
        int size = pattern.Length;
        var rows = new List<int[]>();
        var fill = TrfWraparoundFill(s, n, start, pattern, rows);
        if (fill.Score <= 0)
            return null;

        var columns = new List<TrfColumn>();
        int i = fill.RealRow, r = fill.Row, j = fill.Col, consumed = 0;
        while (rows[r][j] > 0)
        {
            int v = rows[r][j];
            int jp = j == 0 ? size - 1 : j - 1;
            if (v == rows[r - 1][jp] + TrfWeight(s[i], pattern[j]))
            {
                columns.Add(new TrfColumn(s[i], pattern[j], i, j));
                consumed++;
                i--;
                r--;
                j = jp;
            }
            else if (v == rows[r - 1][j] + TrfIndelWeight)
            {
                columns.Add(new TrfColumn(s[i], '-', i, (j + 1) % size));
                i--;
                r--;
            }
            else if (v == rows[r][jp] + TrfIndelWeight)
            {
                columns.Add(new TrfColumn('-', pattern[j], i + 1, j));
                consumed++;
                j = jp;
            }
            else
            {
                throw new InvalidOperationException("Wraparound DP traceback is inconsistent.");
            }
        }

        return new TrfAlignment(fill.Score, TrfCopyNumber(consumed, size), columns.ToArray());
    }

    /// <summary>
    /// Majority-rule consensus from an alignment against a pattern of <paramref name="patternLength"/>
    /// (Benson 1999: "we determine a consensus pattern by majority rule from the alignment of the copies with
    /// P"). Slot 2c+1 is pattern position c, slot 2c the insertion point before it. A position takes the most
    /// frequent aligned symbol, or is deleted when gaps are at least as frequent (ties: gap, A, C, G, T); an
    /// insertion point receives its most frequent inserted base when insertions occur there in at least half
    /// of the passes. Follows the TRF 4.10.0 (non-weighted) consensus rule.
    /// </summary>
    private static char[] TrfConsensus(TrfColumn[] columns, int patternLength)
    {
        int slots = 2 * patternLength + 1;
        var counts = new int[5, slots]; // A C G T gap
        var inserts = new int[slots];
        var passes = new int[slots];

        int last = -1;
        int k = 0;
        while (k < columns.Length)
        {
            int index = columns[k].PatIndex;
            if (index != last)
            {
                int symbol = ConsensusSymbol(columns[k].Seq);
                if (symbol >= 0)
                    counts[symbol, 2 * index + 1]++;
                if (last != -1)
                    passes[index == patternLength - 1 ? 0 : 2 * index + 2]++;
                last = index;
                k++;
            }
            else
            {
                inserts[2 * index]++;
                while (k < columns.Length && columns[k].PatIndex == last)
                {
                    int symbol = ConsensusSymbol(columns[k].Seq);
                    if (symbol >= 0)
                        counts[symbol, 2 * index]++;
                    k++;
                }
            }
        }

        const string Bases = "ACGT";
        var consensus = new List<char>(patternLength + 4);
        for (int slot = 0; slot < slots; slot++)
        {
            if (slot % 2 == 1)
            {
                int best = counts[4, slot];
                char chosen = '-';
                for (int b = 0; b < 4; b++)
                {
                    if (counts[b, slot] > best)
                    {
                        best = counts[b, slot];
                        chosen = Bases[b];
                    }
                }
                if (chosen != '-')
                    consensus.Add(chosen);
            }
            else if (passes[slot] != 0 && (float)inserts[slot] / passes[slot] >= 0.5f)
            {
                int best = counts[0, slot];
                char chosen = 'A';
                for (int b = 1; b < 4; b++)
                {
                    if (counts[b, slot] > best)
                    {
                        best = counts[b, slot];
                        chosen = Bases[b];
                    }
                }
                consensus.Add(chosen);
            }
        }

        return consensus.ToArray();
    }

    /// <summary>Consensus count slot: A/C/G/T → 0..3, gap → 4, other symbols are not counted (−1).</summary>
    private static int ConsensusSymbol(char c) => c == '-' ? 4 : AcgtCode(c);

    /// <summary>Match / mismatch / indel counts between adjacent copies and the TRF period (most common
    /// distance between matching characters).</summary>
    private readonly record struct TrfCopyComparison(int Matches, int Mismatches, int Indels, int Period);

    /// <summary>
    /// Compares ADJACENT copies through the consensus alignment (TRF: statistics refer to "the matches,
    /// mismatches and indels overall between adjacent copies in the sequence, not between the sequence and the
    /// consensus pattern"): two cursors one consensus period apart walk the alignment; aligned symbols of the
    /// two copies are a match or a mismatch, a symbol against a gap is an indel. The period is "the most
    /// common matching distance between corresponding characters in the alignment" (ties: smallest).
    /// </summary>
    private static TrfCopyComparison CompareAdjacentCopies(TrfColumn[] columns)
    {
        int length = columns.Length;
        int lp = 0;
        while (lp < length && columns[lp].Pat == '-')
            lp++;
        int rp = lp + 1;
        while (rp < length && columns[rp].PatIndex != columns[lp].PatIndex)
            rp++;
        while (rp < length && columns[rp].Pat == '-')
            rp++;
        if (rp >= length)
            return new TrfCopyComparison(0, 0, 0, 0);

        int matches = 0, mismatches = 0, indels = 0;
        var distances = new Dictionary<int, int>();
        void Match(int a, int b)
        {
            matches++;
            int distance = Math.Abs(columns[b].SeqIndex - columns[a].SeqIndex);
            distances[distance] = distances.GetValueOrDefault(distance) + 1;
        }

        while (rp < length && lp < rp)
        {
            bool leftGapPat = columns[lp].Pat == '-';
            bool rightGapPat = columns[rp].Pat == '-';
            if (!leftGapPat && !rightGapPat)
            {
                bool leftGapSeq = columns[lp].Seq == '-';
                bool rightGapSeq = columns[rp].Seq == '-';
                if (!leftGapSeq && !rightGapSeq)
                {
                    if (columns[lp].Seq == columns[rp].Seq) Match(lp, rp);
                    else mismatches++;
                }
                else if (leftGapSeq != rightGapSeq)
                {
                    indels++;
                }
                lp++;
                rp++;
            }
            else if (leftGapPat && rightGapPat)
            {
                if (columns[lp].Seq == columns[rp].Seq) Match(lp, rp);
                else mismatches++;
                lp++;
                rp++;
            }
            else if (leftGapPat)
            {
                indels++;
                lp++;
            }
            else
            {
                indels++;
                rp++;
            }
        }

        int period = 0, bestCount = 0;
        foreach (var (distance, count) in distances.OrderBy(kv => kv.Key))
        {
            if (count > bestCount)
            {
                bestCount = count;
                period = distance;
            }
        }

        return new TrfCopyComparison(matches, mismatches, indels, period);
    }

    /// <summary>Builds the reported TRF statistics from the final (consensus) alignment.</summary>
    private static ApproximateTandemRepeatResult TrfStatistics(char[] s, TrfAlignment alignment, char[] consensus)
    {
        var copies = CompareAdjacentCopies(alignment.Columns);
        int trials = copies.Matches + copies.Mismatches + copies.Indels;

        int first = alignment.First, last = alignment.Last, span = last - first + 1;
        var region = new string(s, first, span);
        int a = 0, c = 0, g = 0, t = 0;
        foreach (char ch in region)
        {
            switch (ch)
            {
                case 'A': a++; break;
                case 'C': c++; break;
                case 'G': g++; break;
                case 'T': t++; break;
            }
        }

        // TRF prints the consensus starting at the pattern position aligned with the first repeat base.
        int phase = alignment.Columns[^1].PatIndex;
        string rotated = new string(consensus, phase, consensus.Length - phase) + new string(consensus, 0, phase);

        return new ApproximateTandemRepeatResult(
            Start: first - 1,
            SpanLength: span,
            Period: copies.Period,
            ConsensusSize: consensus.Length,
            Consensus: rotated,
            CopyNumber: alignment.CopyNumber,
            PercentMatches: trials > 0 ? 100.0 * copies.Matches / trials : 0.0,
            PercentIndels: trials > 0 ? 100.0 * copies.Indels / trials : 0.0,
            AlignmentScore: alignment.Score)
        {
            PercentA = 100.0 * a / span,
            PercentC = 100.0 * c / span,
            PercentG = 100.0 * g / span,
            PercentT = 100.0 * t / span,
            Entropy = SequenceComplexity.CalculateShannonEntropy(region),
        };
    }

    /// <summary>
    /// TRF multiples test: the candidate distance must be one of the three best periods of the aligned region
    /// (Benson 1999 / TRF README: redundant reporting at multiples of the pattern size is limited "to, at most,
    /// three pattern sizes"). Best periods = the largest counts of distances between identical dinucleotides in
    /// the region after removing their least-squares linear trend (TRF 4.10.0 method; non-ACGT symbols count as
    /// A as in TRF). Period 1 instead requires ≥ 80 % of the region to be one base.
    /// </summary>
    private static bool IsAmongTrfBestPeriods(
        char[] s, int first, int last, int d, int maxPeriod, Dictionary<(int First, int Last), int[]> cache)
    {
        int length = last - first + 1;
        if (d == 1)
        {
            var composition = new int[4];
            for (int p = first; p <= last; p++)
                composition[AcgtIndexOrA(s[p])]++;
            return composition.Max() * 100.0f / length >= 80.0f;
        }

        // The best periods depend only on the region; different distances often align the same region.
        if (!cache.TryGetValue((first, last), out int[]? best))
        {
            best = TrfBestPeriods(s, first, length, maxPeriod);
            cache[(first, last)] = best;
        }

        return Array.IndexOf(best, d) >= 0;
    }

    /// <summary>The three best periods of s[first..first+length−1] (see <see cref="IsAmongTrfBestPeriods"/>).</summary>
    private static int[] TrfBestPeriods(char[] s, int first, int length, int maxPeriod)
    {
        var best = new int[3];
        int end = length - 2;
        if (end < 1)
            return best;

        var counts = new double[length];
        var history = new int[length];
        var heads = new int[16];
        Array.Fill(heads, -1);
        const int MaxCountedDistance = 3 * MaxApproximatePeriod;
        for (int p = 0; p <= end; p++)
        {
            int tuple = AcgtIndexOrA(s[first + p]) * 4 + AcgtIndexOrA(s[first + p + 1]);
            history[p] = heads[tuple];
            heads[tuple] = p;
            int distance = 0;
            for (int cur = p; history[cur] != -1 && distance < MaxCountedDistance; cur = history[cur])
            {
                distance = p - history[cur];
                counts[distance] += 1.0;
            }
        }

        double xy = 0, x = 0, y = 0, x2 = 0;
        for (int q = 1; q <= end; q++)
        {
            xy += q * counts[q];
            x += q;
            y += counts[q];
            x2 += (double)q * q;
        }
        double slope = (end * xy - x * y) / (end * x2 - x * x);
        for (int q = 1; q <= end; q++)
            counts[q] -= q * slope;

        int top = Math.Min(end, maxPeriod);
        for (int pick = 0; pick < best.Length; pick++)
        {
            int bestIndex = 0;
            double bestValue = 0.0;
            for (int q = 1; q <= top; q++)
            {
                if (counts[q] > bestValue)
                {
                    bestIndex = q;
                    bestValue = counts[q];
                }
            }
            best[pick] = bestIndex;
            counts[bestIndex] = 0.0;
        }

        return best;
    }

    /// <summary>A/C/G/T → 0..3 with any other symbol counted as A (TRF's zero-initialised index table).</summary>
    private static int AcgtIndexOrA(char c) => Math.Max(0, AcgtCode(c));

    /// <summary>
    /// TRF redundancy elimination over repeats sorted by start: of two repeats overlapping by ≥ 90 % of one of
    /// them, that one is dropped when it has the same period and no higher score, or a period that is a
    /// multiple of the other's and a score ≤ 1.1× (TRF README "Redundancy": the same repeat detected at
    /// several period sizes / "the same period size may be detected more than once").
    /// </summary>
    private static List<ApproximateTandemRepeatResult> RemoveTrfRedundancy(List<ApproximateTandemRepeatResult> repeats)
    {
        int i = 0;
        while (i < repeats.Count)
        {
            bool removedI = false;
            int j = i + 1;
            while (j < repeats.Count)
            {
                var a = repeats[i];
                var b = repeats[j];
                int overlap = Math.Min(a.Start + a.SpanLength, b.Start + b.SpanLength) - Math.Max(a.Start, b.Start);
                if (overlap <= 0)
                    break;
                if (!(overlap / (double)a.SpanLength < 0.9) && IsTrfRedundant(a, b))
                {
                    repeats.RemoveAt(i);
                    removedI = true;
                    break;
                }
                if (!(overlap / (double)b.SpanLength < 0.9) && IsTrfRedundant(b, a))
                {
                    repeats.RemoveAt(j);
                    continue;
                }
                j++;
            }

            if (!removedI)
                i++;
        }

        return repeats;
    }

    private static bool IsTrfRedundant(ApproximateTandemRepeatResult x, ApproximateTandemRepeatResult y) =>
        (x.Period > y.Period && y.Period > 0 && x.Period % y.Period == 0 && x.AlignmentScore <= 1.1 * y.AlignmentScore) ||
        (x.Period == y.Period && x.AlignmentScore <= y.AlignmentScore);

    #endregion

    #region TRF Bernoulli statistics (Benson 1999)

    // Benson (1999) / TRF README "Probabilistic Model of Tandem Repeats": "We model alignment of two tandem
    // copies of a pattern of length n by a sequence of n independent Bernoulli trials ... P(Heads), which we
    // also call PM or matching probability, represents the average percent identity between the copies ...
    // PI or indel probability specifies the average percentage of insertions and deletions between the
    // copies." Default "(PM = .80, PI = .10)"; "Probabilistic data is available for PM values of 80 and 75 and
    // PI values of 10 and 20."

    /// <summary>Benson (1999) default Bernoulli matching probability PM = 0.80.</summary>
    public const double TrfDefaultMatchProbability = 0.80;

    /// <summary>Benson (1999) default Bernoulli indel probability PI = 0.10.</summary>
    public const double TrfDefaultIndelProbability = 0.10;

    /// <summary>
    /// Estimates the TRF Bernoulli-model parameters of a tandem-repeat tract: PM (match probability) and PI
    /// (indel probability) between ADJACENT copies. The tract is analysed exactly as TRF analyses a detected
    /// repeat: wraparound-DP alignment against tandem copies of the candidate pattern (the last
    /// <paramref name="period"/> bases of the tract), majority-rule consensus, realignment against the
    /// consensus, then comparison of each copy with the next one through that alignment (TRF: statistics refer
    /// to "the matches, mismatches and indels overall between adjacent copies in the sequence, not between the
    /// sequence and the consensus pattern"). Each compared column is one Bernoulli trial: heads = match,
    /// tails = mismatch or indel. PM and PI therefore equal TRF's reported % matches / % indels (/100) for
    /// the same region.
    /// </summary>
    /// <remarks>
    /// The statistics cover the locally aligned part of the tract (flanks that do not align are ignored, as in
    /// TRF). A tract without two aligned copies yields zero trials and PM = PI = 0.
    /// <see cref="TandemRepeatBernoulliStatistics.ExpectedMatches"/> is PM × trials (= matches).
    /// </remarks>
    /// <param name="repeatTract">The tandem-repeat tract (≥ 2 × period symbols; case-insensitive).</param>
    /// <param name="period">Candidate period of the tract (1..2000).</param>
    /// <param name="expectedMatchProbability">PM the tract is compared with (default 0.80, Benson 1999).</param>
    public static TandemRepeatBernoulliStatistics ComputeBernoulliStatistics(
        string repeatTract,
        int period,
        double expectedMatchProbability = TrfDefaultMatchProbability)
    {
        ArgumentNullException.ThrowIfNull(repeatTract);
        ArgumentOutOfRangeException.ThrowIfLessThan(period, 1);
        ArgumentOutOfRangeException.ThrowIfGreaterThan(period, MaxApproximatePeriod);
        if (!(expectedMatchProbability >= 0.0 && expectedMatchProbability <= 1.0))
            throw new ArgumentOutOfRangeException(nameof(expectedMatchProbability));
        if (repeatTract.Length < period * 2)
            throw new ArgumentException(
                "A tandem repeat needs at least two contiguous copies of the period.", nameof(repeatTract));

        int n = repeatTract.Length;
        var s = new char[n + 1];
        repeatTract.ToUpperInvariant().CopyTo(0, s, 1, n);
        var pattern = new char[period];
        Array.Copy(s, n - period + 1, pattern, 0, period);

        TrfCopyComparison copies = default;
        double copyNumber = 0.0;
        var initial = TrfWraparoundAlign(s, n, n, pattern);
        if (initial is not null)
        {
            char[] consensus = TrfConsensus(initial.Columns, period);
            var final = consensus.Length > 0 ? TrfWraparoundAlign(s, n, n, consensus) : null;
            if (final is not null)
            {
                copies = CompareAdjacentCopies(final.Columns);
                copyNumber = final.CopyNumber;
            }
        }

        int trials = copies.Matches + copies.Mismatches + copies.Indels;
        double pm = trials > 0 ? (double)copies.Matches / trials : 0.0;
        double pi = trials > 0 ? (double)copies.Indels / trials : 0.0;

        return new TandemRepeatBernoulliStatistics(
            Period: period,
            AdjacentCopyPairs: Math.Max(0, (int)Math.Ceiling(copyNumber) - 1),
            BernoulliTrials: trials,
            Matches: copies.Matches,
            Mismatches: copies.Mismatches,
            Indels: copies.Indels,
            MatchProbability: pm,
            IndelProbability: pi,
            PercentMatches: pm * 100.0,
            PercentIndels: pi * 100.0,
            ExpectedMatches: pm * trials,
            MeetsExpectedMatchProbability: trials > 0 && pm >= expectedMatchProbability);
    }

    #endregion

    #region Inverted Repeat Detection

    /// <summary>
    /// Finds inverted repeats: a left arm followed, after a loop, by a right arm equal to the reverse complement of
    /// the left arm — exact (perfect) stems by default, optionally with up to <c>maxMismatches</c> mismatched pairs,
    /// a maximum arm length and G·U wobble pairs. Such structures can form hairpin/stem-loop (single strand) or
    /// cruciform (duplex) structures.
    /// </summary>
    /// <remarks>
    /// <para>
    /// <b>Model (default: exact, maximal stems).</b> A stem is a triple (<c>LeftArmStart = i</c>,
    /// <c>RightArmStart = j</c>, <c>ArmLength = A</c>) such that base <c>s[i+k]</c> is the Watson–Crick
    /// complement of <c>s[j+A−1−k]</c> for every <c>0 ≤ k &lt; A</c>, with
    /// <c>A ≥ minArmLength</c> and <c>minLoopLength ≤ j − (i + A) ≤ maxLoopLength</c>.
    /// Only stems that are <b>not contained, in both arms, in another such stem</b> are reported. This is
    /// the reporting rule of EMBOSS <c>palindrome</c> (M. Faller; <c>palindrome.c</c>: every stem is extended
    /// inward from its outer pair and a stem that is a subset of an already found stem in both halves —
    /// <c>palindrome_AInB</c> — is dropped, default <c>-overlap Y</c>), run with <c>-nummismatches 0</c>.
    /// Consequently every reported stem is maximal: it cannot be extended outward, and it can be extended
    /// inward only by making the loop shorter than <paramref name="minLoopLength"/>; sub-stems and
    /// "slipped" re-pairings lying inside a longer stem are not reported. With
    /// <paramref name="minLoopLength"/> = 0 the result set equals EMBOSS <c>palindrome</c>
    /// (<c>-nummismatches 0 -overlap Y</c>, <c>-gaplimit = maxLoopLength</c>, <c>-minpallen = minArmLength</c>,
    /// unbounded <c>-maxpallen</c>), cross-checked on random sequences.
    /// </para>
    /// <para>
    /// <b>Pairing.</b> By default only the unambiguous bases A, C, G, T pair, via the canonical
    /// <see cref="SequenceExtensions.GetComplementBase(char)"/> (A↔T, C↔G). As in EMBOSS <c>einverted</c>
    /// (which scores a match only for a/c/g/t), N and other IUPAC ambiguity codes never form a pair, so a run
    /// of N is not reported as a stem.
    /// </para>
    /// <para>
    /// <b>Options.</b> <c>maxMismatches &gt; 0</c> and a finite <c>maxArmLength</c> follow EMBOSS <c>palindrome</c>
    /// <c>-nummismatches</c> / <c>-maxpallen</c> exactly (see the parameters; with <c>minLoopLength = 0</c> the result
    /// set equals the EMBOSS 6.6.0 binary, cross-checked on thousands of random cases); <c>allowWobble</c> adds G·U
    /// pairs under the same maximal-stem rule. No gaps (bulges): for scored, gap-tolerant inverted repeats use
    /// <see cref="FindInvertedRepeatsScored(string, int, int, int, int, int)"/> (EMBOSS <c>einverted</c>).
    /// </para>
    /// <para>Coordinates are 0-based; results are ordered by <c>LeftArmStart</c>, then <c>RightArmStart</c>.
    /// Parameters are validated eagerly.</para>
    /// </remarks>
    /// <param name="sequence">DNA sequence to search.</param>
    /// <param name="minArmLength">Minimum length of each arm (default: 4, must be ≥ 2).</param>
    /// <param name="maxLoopLength">Maximum loop length between arms (default: 50, must be ≥ <paramref name="minLoopLength"/>).</param>
    /// <param name="minLoopLength">Minimum loop length (default: 3, must be ≥ 0).</param>
    /// <param name="maxMismatches">
    /// Maximum number of mismatched pairs inside a stem (default 0 = exact stems; must be ≥ 0). EMBOSS
    /// <c>palindrome -nummismatches</c> semantics: every stem starts at a pairing outer pair and is extended inward
    /// pair by pair until the (<paramref name="maxMismatches"/>+1)-th mismatch (or until the loop would become shorter
    /// than <paramref name="minLoopLength"/>); mismatches at the inner end of the stem are trimmed, interior mismatches
    /// are kept and counted in <see cref="InvertedRepeatResult.Mismatches"/>; a stem lying inside another stem in both
    /// arms is dropped (<c>-overlap Y</c>).
    /// </param>
    /// <param name="maxArmLength">
    /// Maximum arm length (default <see cref="int.MaxValue"/> = unbounded; must be ≥ <paramref name="minArmLength"/>).
    /// EMBOSS <c>palindrome -maxpallen</c> semantics: stems are only started from outer pairs spanning at most
    /// <c>2·maxArmLength + maxLoopLength + 1</c> bases, and stems longer than <paramref name="maxArmLength"/> are not
    /// reported but still suppress the stems they contain (so a longer stem is <b>not</b> split into shorter pieces;
    /// only parts of it that do not lie inside a longer candidate can surface).
    /// </param>
    /// <param name="allowWobble">
    /// When true, G·T (G·U) wobble pairs count as pairs in addition to Watson–Crick pairs, using the canonical
    /// <see cref="RnaSecondaryStructure.CanPair(char, char)"/> (A·U, G·C, G·U — T is read as U; Crick 1966,
    /// Varani &amp; McClain 2000). Default false (Watson–Crick A·T/C·G only, U never pairs).
    /// </param>
    /// <returns>The maximal (non-nested) stems, ordered by <c>LeftArmStart</c>, then <c>RightArmStart</c>.</returns>
    public static IEnumerable<InvertedRepeatResult> FindInvertedRepeats(
        DnaSequence sequence,
        int minArmLength = 4,
        int maxLoopLength = 50,
        int minLoopLength = 3,
        int maxMismatches = 0,
        int maxArmLength = int.MaxValue,
        bool allowWobble = false)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ValidateInvertedRepeatParameters(minArmLength, maxLoopLength, minLoopLength, maxMismatches, maxArmLength);

        return FindInvertedRepeatsDispatch(
            sequence.Sequence, minArmLength, maxLoopLength, minLoopLength, maxMismatches, maxArmLength, allowWobble);
    }

    /// <summary>
    /// Finds maximal inverted repeats in a raw sequence string (case-insensitive).
    /// Same model, options and validation as
    /// <see cref="FindInvertedRepeats(DnaSequence, int, int, int, int, int, bool)"/>;
    /// <c>null</c> or empty input yields no results. Characters other than A/C/G/T never pair
    /// (with <paramref name="allowWobble"/>, U pairs as T).
    /// </summary>
    public static IEnumerable<InvertedRepeatResult> FindInvertedRepeats(
        string sequence,
        int minArmLength = 4,
        int maxLoopLength = 50,
        int minLoopLength = 3,
        int maxMismatches = 0,
        int maxArmLength = int.MaxValue,
        bool allowWobble = false)
    {
        ValidateInvertedRepeatParameters(minArmLength, maxLoopLength, minLoopLength, maxMismatches, maxArmLength);

        if (string.IsNullOrEmpty(sequence))
            return Array.Empty<InvertedRepeatResult>();

        return FindInvertedRepeatsDispatch(
            sequence.ToUpperInvariant(), minArmLength, maxLoopLength, minLoopLength, maxMismatches, maxArmLength, allowWobble);
    }

    private static void ValidateInvertedRepeatParameters(
        int minArmLength, int maxLoopLength, int minLoopLength, int maxMismatches, int maxArmLength)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(minArmLength, 2);
        ArgumentOutOfRangeException.ThrowIfNegative(minLoopLength);
        ArgumentOutOfRangeException.ThrowIfLessThan(maxLoopLength, minLoopLength);
        ArgumentOutOfRangeException.ThrowIfNegative(maxMismatches);
        ArgumentOutOfRangeException.ThrowIfLessThan(maxArmLength, minArmLength);
    }

    private static List<InvertedRepeatResult> FindInvertedRepeatsDispatch(
        string seq, int minArmLength, int maxLoopLength, int minLoopLength,
        int maxMismatches, int maxArmLength, bool allowWobble)
    {
        // An arm can never exceed n/2, so maxArmLength ≥ n/2 removes neither a start pair nor a stem.
        bool armBounded = maxArmLength < seq.Length / 2;
        if (maxMismatches == 0 && !armBounded)
            return FindInvertedRepeatsCore(seq, minArmLength, maxLoopLength, minLoopLength, allowWobble);

        int armLimit = armBounded ? maxArmLength : int.MaxValue;
        return FindInvertedRepeatsMismatchCore(
            seq, minArmLength, maxLoopLength, minLoopLength, maxMismatches, armLimit, allowWobble);
    }

    /// <summary>True when <paramref name="a"/> and <paramref name="b"/> form a Watson–Crick pair (ACGT only).</summary>
    private static bool IsWatsonCrickPair(char a, char b) =>
        AcgtCode(a) >= 0 && SequenceExtensions.GetComplementBase(a) == b;

    /// <summary>Stem pairing rule: Watson–Crick (ACGT), or the canonical RNA pair set incl. G·U when <paramref name="wobble"/>.</summary>
    private static bool IsStemPair(char a, char b, bool wobble) =>
        wobble ? RnaSecondaryStructure.CanPair(a, b) : IsWatsonCrickPair(a, b);

    /// <summary>
    /// Maximal exact stems. A stem is identified by its innermost pair (<c>iIn</c>, <c>rIn</c>) with loop
    /// <c>L = rIn − iIn − 1</c>; it is inward-maximal when the next inner pair does not pair or would leave a loop
    /// shorter than <c>minLoop</c>. Each inward-maximal stem is extended outward as far as the bases pair
    /// (EMBOSS palindrome's candidate set with 0 mismatches); then stems contained in both arms of another stem
    /// are removed (<c>palindrome_AInB</c>, see <see cref="IsContainedInShiftedStem"/>).
    /// </summary>
    private static List<InvertedRepeatResult> FindInvertedRepeatsCore(
        string seq,
        int minArmLength,
        int maxLoopLength,
        int minLoopLength,
        bool wobble)
    {
        int n = seq.Length;
        var results = new List<InvertedRepeatResult>();
        if (n < 2 * minArmLength + minLoopLength)
            return results;

        // Largest loop that can still fit two minimum arms; caps a huge maxLoopLength (no overflow / no hang).
        int maxLoop = (int)Math.Min(maxLoopLength, (long)n - 2 * minArmLength);

        for (int iIn = minArmLength - 1; iIn < n; iIn++)
        {
            for (int loop = minLoopLength; loop <= maxLoop; loop++)
            {
                int rIn = iIn + loop + 1;
                if (rIn + minArmLength - 1 >= n)
                    break;
                if (!IsStemPair(seq[iIn], seq[rIn], wobble))
                    continue;
                // Inward-maximal: the next inner pair must fail, unless it would make the loop < minLoop.
                if (loop - 2 >= minLoopLength && IsStemPair(seq[iIn + 1], seq[rIn - 1], wobble))
                    continue;

                int arm = 1;
                while (iIn - arm >= 0 && rIn + arm < n && IsStemPair(seq[iIn - arm], seq[rIn + arm], wobble))
                    arm++;
                if (arm < minArmLength)
                    continue;

                int iOut = iIn - arm + 1;
                if (IsContainedInShiftedStem(seq, iOut, rIn, arm, loop - minLoopLength, wobble))
                    continue;

                results.Add(new InvertedRepeatResult(
                    LeftArmStart: iOut,
                    RightArmStart: rIn,
                    ArmLength: arm,
                    LoopLength: loop,
                    LeftArm: seq.Substring(iOut, arm),
                    RightArm: seq.Substring(rIn, arm),
                    Loop: seq.Substring(iIn + 1, loop),
                    CanFormHairpin: loop >= 3));
            }
        }

        results.Sort(static (x, y) =>
        {
            int c = x.LeftArmStart.CompareTo(y.LeftArmStart);
            return c != 0 ? c : x.RightArmStart.CompareTo(y.RightArmStart);
        });
        return results;
    }

    /// <summary>
    /// True when the stem (<paramref name="left"/>, <paramref name="right"/>, <paramref name="arm"/>) lies, in both
    /// arms, inside a stem on another diagonal. Covering both arms of the stem from a diagonal shifted by
    /// <c>m ≥ 1</c> needs a stem with arm <c>arm + m</c> and loop <c>loop − m</c> that either keeps the left start and
    /// the right start (<c>(left, right, arm + m)</c>) or keeps both arm ends (<c>(left − m, right − m, arm + m)</c>);
    /// every containing stem contains one of these, so testing <c>1 ≤ m ≤ loop − minLoop</c> is exact.
    /// </summary>
    private static bool IsContainedInShiftedStem(string seq, int left, int right, int arm, int maxShift, bool wobble)
    {
        int n = seq.Length;
        for (int m = 1; m <= maxShift; m++)
        {
            int len = arm + m;
            if (right + len <= n && IsExactStem(seq, left, right, len, wobble))
                return true;
            if (left - m >= 0 && IsExactStem(seq, left - m, right - m, len, wobble))
                return true;
        }
        return false;
    }

    /// <summary>True when s[i+k] pairs with s[j+len−1−k] for all k (checked from the innermost pair out).</summary>
    private static bool IsExactStem(string seq, int i, int j, int len, bool wobble)
    {
        for (int k = len - 1; k >= 0; k--)
        {
            if (!IsStemPair(seq[i + k], seq[j + len - 1 - k], wobble))
                return false;
        }
        return true;
    }

    /// <summary>
    /// Mismatch-tolerant and/or arm-bounded stems: the candidate set of EMBOSS <c>palindrome</c>
    /// (<c>-nummismatches k</c>, <c>-maxpallen</c>) followed by the <c>palindrome_AInB</c> nesting filter.
    /// All pairs of a stem lie on one anti-diagonal <c>D = left + right</c>; along a diagonal the pairs are indexed
    /// from the innermost admissible pair (loop ≥ <paramref name="minLoop"/>, index 0) outward. A candidate starts at a
    /// pairing outer index <c>u</c> and walks inward to just before the (k+1)-th mismatch (or to index 0); trailing
    /// mismatches are trimmed, so its inner end is the first pairing index above that stop. Kept when
    /// arm ≥ <paramref name="minArm"/>, loop ≤ <paramref name="maxLoopLength"/>, and (bounded arm) the outer span is at
    /// most <c>2·maxArm + maxLoopLength + 1</c>; a candidate lying inside another candidate in both arms is dropped
    /// (any container differs by at most <c>loop − minLoop</c> diagonals); candidates longer than
    /// <paramref name="maxArm"/> are then not reported (<c>palindrome_Print</c> filter).
    /// </summary>
    private static List<InvertedRepeatResult> FindInvertedRepeatsMismatchCore(
        string seq, int minArm, int maxLoopLength, int minLoop, int maxMismatches, int maxArm, bool wobble)
    {
        int n = seq.Length;
        var results = new List<InvertedRepeatResult>();
        if (n < 2 * minArm + minLoop)
            return results;

        long maxLoop = Math.Min(maxLoopLength, (long)n);
        long maxSpan = maxArm == int.MaxValue ? long.MaxValue : 2L * maxArm + maxLoopLength + 1;
        int diagonals = 2 * n - 1;

        // Candidates grouped by diagonal; within a diagonal stored with increasing outer-left position.
        var outerLeft = new List<int>();
        var innerLeft = new List<int>();
        var mismatchCount = new List<int>();
        var diagStart = new int[diagonals + 1];
        var mismatchAt = new List<int>();
        var nextPairAfter = new List<int>();
        var mismatchesBelow = new List<int>(); // per index u: mismatches at indices < u
        var diagOuter = new List<int>();
        var diagInner = new List<int>();
        var diagMism = new List<int>();

        for (int d = 0; d < diagonals; d++)
        {
            diagStart[d] = outerLeft.Count;
            long a0Long = (long)d - 1 - minLoop;
            if (a0Long < 0)
                continue;
            int a0 = (int)(a0Long / 2);
            int loop0 = d - 2 * a0 - 1;
            if (loop0 > maxLoop)
                continue;
            int innerLimit = (int)((maxLoop - loop0) / 2);
            long uMaxLong = Math.Min(a0, (long)n - 1 - d + a0);
            if (maxSpan != long.MaxValue)
            {
                long spanRoom = maxSpan - loop0 - 2;
                uMaxLong = spanRoom < 0 ? -1 : Math.Min(uMaxLong, spanRoom / 2);
            }
            int uMax = (int)uMaxLong;

            mismatchAt.Clear();
            nextPairAfter.Clear();
            mismatchesBelow.Clear();
            diagOuter.Clear();
            diagInner.Clear();
            diagMism.Clear();
            int pending = 0;
            int firstPair = -1;
            for (int u = 0; u <= uMax; u++)
            {
                mismatchesBelow.Add(mismatchAt.Count);
                if (!IsStemPair(seq[a0 - u], seq[d - a0 + u], wobble))
                {
                    mismatchAt.Add(u);
                    nextPairAfter.Add(-1);
                    // Every later start stops at or above this (k+1)-th mismatch → its inner end exceeds innerLimit.
                    if (mismatchAt.Count > maxMismatches && mismatchAt[^(maxMismatches + 1)] >= innerLimit)
                        break;
                    continue;
                }

                if (firstPair < 0)
                    firstPair = u;
                for (; pending < mismatchAt.Count; pending++)
                    nextPairAfter[pending] = u;

                // Inner end: first pairing index above the (k+1)-th mismatch met walking inward (trailing mismatches trimmed).
                int inner = mismatchAt.Count > maxMismatches
                    ? nextPairAfter[mismatchAt.Count - 1 - maxMismatches]
                    : firstPair;
                int mism = mismatchAt.Count - mismatchesBelow[inner];

                if (u - inner + 1 >= minArm && inner <= innerLimit)
                {
                    diagOuter.Add(a0 - u);
                    diagInner.Add(a0 - inner);
                    diagMism.Add(mism);
                }
            }

            for (int c = diagOuter.Count - 1; c >= 0; c--)
            {
                outerLeft.Add(diagOuter[c]);
                innerLeft.Add(diagInner[c]);
                mismatchCount.Add(diagMism[c]);
            }
        }
        diagStart[diagonals] = outerLeft.Count;

        // Prefix maximum of the inner-left end per diagonal (outer-left ascending) for containment queries, and the
        // longest candidate arm per diagonal / overall (a container shifted by m diagonals is ≥ m bases longer).
        var prefixMaxInner = new int[outerLeft.Count];
        var longestOnDiagonal = new int[diagonals];
        int longestArm = 0;
        for (int d = 0; d < diagonals; d++)
        {
            int best = int.MinValue;
            for (int c = diagStart[d]; c < diagStart[d + 1]; c++)
            {
                prefixMaxInner[c] = best = Math.Max(best, innerLeft[c]);
                longestOnDiagonal[d] = Math.Max(longestOnDiagonal[d], innerLeft[c] - outerLeft[c] + 1);
            }
            longestArm = Math.Max(longestArm, longestOnDiagonal[d]);
        }

        // True when diagonal dd holds a candidate with outer-left ≤ x and inner-left ≥ y.
        bool Covered(int dd, int x, int y)
        {
            if (longestOnDiagonal[dd] < y - x + 1)
                return false;
            int lo = diagStart[dd], hi = diagStart[dd + 1] - 1, found = -1;
            while (lo <= hi)
            {
                int mid = (lo + hi) >>> 1;
                if (outerLeft[mid] <= x) { found = mid; lo = mid + 1; }
                else hi = mid - 1;
            }
            return found >= 0 && prefixMaxInner[found] >= y;
        }

        for (int d = 0; d < diagonals; d++)
        {
            for (int c = diagStart[d]; c < diagStart[d + 1]; c++)
            {
                int fs = outerLeft[c], fe = innerLeft[c];
                int arm = fe - fs + 1;
                if (arm > maxArm)
                    continue;
                int loop = d - 2 * fe - 1;
                int rightOuter = d - fs, rightInner = d - fe;

                bool contained = Covered(d, fs - 1, fe);
                int maxShift = Math.Min(loop - minLoop, longestArm - arm);
                for (int m = 1; !contained && m <= maxShift; m++)
                {
                    int up = d + m, down = d - m;
                    if (up < diagonals && Covered(up, Math.Min(fs, up - rightOuter), Math.Max(fe, up - rightInner)))
                        contained = true;
                    else if (down >= 0 && Covered(down, Math.Min(fs, down - rightOuter), Math.Max(fe, down - rightInner)))
                        contained = true;
                }
                if (contained)
                    continue;

                results.Add(new InvertedRepeatResult(
                    LeftArmStart: fs,
                    RightArmStart: rightInner,
                    ArmLength: arm,
                    LoopLength: loop,
                    LeftArm: seq.Substring(fs, arm),
                    RightArm: seq.Substring(rightInner, arm),
                    Loop: seq.Substring(fe + 1, loop),
                    CanFormHairpin: loop >= 3)
                { Mismatches = mismatchCount[c] });
            }
        }

        results.Sort(static (x, y) =>
        {
            int c = x.LeftArmStart.CompareTo(y.LeftArmStart);
            return c != 0 ? c : x.RightArmStart.CompareTo(y.RightArmStart);
        });
        return results;
    }

    #endregion

    #region Scored (einverted) Inverted Repeat Detection

    /// <summary>EMBOSS <c>einverted</c> default gap penalty (<c>-gap 12</c>).</summary>
    public const int EinvertedDefaultGapPenalty = 12;

    /// <summary>EMBOSS <c>einverted</c> default minimum score (<c>-threshold 50</c>).</summary>
    public const int EinvertedDefaultThreshold = 50;

    /// <summary>EMBOSS <c>einverted</c> default match score (<c>-match 3</c>).</summary>
    public const int EinvertedDefaultMatchScore = 3;

    /// <summary>EMBOSS <c>einverted</c> default mismatch score (<c>-mismatch -4</c>).</summary>
    public const int EinvertedDefaultMismatchScore = -4;

    /// <summary>
    /// EMBOSS 6.6 <c>einverted</c> default maximum repeat extent (<c>-maxrepeat 2000</c>, einverted.acd). Durbin's
    /// original compile-time value was 4000 (comment in <c>einverted.c</c>); pass it explicitly to reproduce that.
    /// </summary>
    public const int EinvertedDefaultMaxRepeatLength = 2000;

    /// <summary>Sentinel score of <c>einverted</c> (<c>rogue</c>) marking cells outside the sequence/window.</summary>
    private const int EinvertedRogue = 1_000_000;

    /// <summary>
    /// Finds imperfect (mismatch- and gap-tolerant) inverted repeats by score, reproducing EMBOSS <c>einverted</c>
    /// (Durbin &amp; Thierry-Mieg 1993, "Inverted repeats by dynamic programming").
    /// </summary>
    /// <remarks>
    /// <para>
    /// <b>Scoring.</b> A local alignment of the sequence against its own reverse complement, grown outward from the
    /// loop: cell (<c>i</c>, <c>k</c>) pairs right-arm base <c>i</c> with left-arm base <c>i − 1 − k</c> and scores
    /// <c>H(i,k) = max(s(i, i−1−k) + max(0, H(i−1, k−2)), max(H(i, k−1), H(i−1, k−1)) − gap)</c>, where <c>s</c> is
    /// <paramref name="matchScore"/> for a Watson–Crick pair of A/C/G/T and <paramref name="mismatchScore"/>
    /// otherwise (N and other symbols never pair). Only cells with <c>k &lt; maxRepeatLength − 2</c> (outer span at
    /// most <paramref name="maxRepeatLength"/>) are scored.
    /// </para>
    /// <para>
    /// <b>Reporting.</b> Exactly the <c>einverted</c> procedure: for each right end <c>i</c> the first best cell with
    /// score ≥ <paramref name="threshold"/> is remembered against its left start (a later end with the same left start
    /// replaces it only with a higher score; left starts at or before the last reported right end are ignored); once
    /// the scan is <paramref name="maxRepeatLength"/> past a remembered left start, the best-scoring end in the window
    /// is traced back (same-row gap, previous-row gap, then diagonal) and reported, and the window is cleared.
    /// Results are returned in einverted's report order (not necessarily sorted by position). Degenerate parameters
    /// (threshold ≤ match score, zero gap penalty) can make einverted 6.6.0 abort (SIGFPE: its percentage divides by
    /// matches + mismatches): a trace-back that records no column is consumed without a result, and a repeat made
    /// only of gap columns is returned with <see cref="ScoredInvertedRepeatResult.PercentMatches"/> = 0 (the scan
    /// otherwise continues exactly as einverted's code).
    /// </para>
    /// <para>
    /// Cross-checked against the EMBOSS 6.6.0 <c>einverted</c> binary (coordinates, score, matches, mismatches,
    /// gaps and both alignment rows) on random sequences with planted imperfect inverted repeats. Output differences
    /// by design: coordinates are 0-based and inclusive; alignment rows use the upper-cased input characters (einverted
    /// prints lower case and shows non-ACGT symbols as '-').
    /// </para>
    /// </remarks>
    /// <param name="sequence">DNA sequence (case-insensitive; symbols other than A/C/G/T never pair).</param>
    /// <param name="gapPenalty">Gap penalty (≥ 0, default 12).</param>
    /// <param name="threshold">Minimum reported score (≥ 0, default 50).</param>
    /// <param name="matchScore">Score of a Watson–Crick pair (≥ 0, default 3).</param>
    /// <param name="mismatchScore">Score of any other pair (≤ 0, default −4).</param>
    /// <param name="maxRepeatLength">
    /// Maximum extent from the start of the repeat to the end of its inverted copy (einverted <c>-maxrepeat</c>,
    /// ≥ 2, default 2000). Memory is O(min(maxRepeatLength, n)²).
    /// </param>
    /// <returns>Scored inverted repeats in einverted report order.</returns>
    public static IEnumerable<ScoredInvertedRepeatResult> FindInvertedRepeatsScored(
        string sequence,
        int gapPenalty = EinvertedDefaultGapPenalty,
        int threshold = EinvertedDefaultThreshold,
        int matchScore = EinvertedDefaultMatchScore,
        int mismatchScore = EinvertedDefaultMismatchScore,
        int maxRepeatLength = EinvertedDefaultMaxRepeatLength)
    {
        ValidateScoredInvertedRepeatParameters(gapPenalty, threshold, matchScore, mismatchScore, maxRepeatLength);
        if (string.IsNullOrEmpty(sequence))
            return Array.Empty<ScoredInvertedRepeatResult>();
        return new EinvertedScanner(sequence.ToUpperInvariant(), gapPenalty, threshold, matchScore, mismatchScore, maxRepeatLength).Run();
    }

    /// <summary>
    /// Scored inverted repeats of a <see cref="DnaSequence"/>; see
    /// <see cref="FindInvertedRepeatsScored(string, int, int, int, int, int)"/>.
    /// </summary>
    public static IEnumerable<ScoredInvertedRepeatResult> FindInvertedRepeatsScored(
        DnaSequence sequence,
        int gapPenalty = EinvertedDefaultGapPenalty,
        int threshold = EinvertedDefaultThreshold,
        int matchScore = EinvertedDefaultMatchScore,
        int mismatchScore = EinvertedDefaultMismatchScore,
        int maxRepeatLength = EinvertedDefaultMaxRepeatLength)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        return FindInvertedRepeatsScored(sequence.Sequence, gapPenalty, threshold, matchScore, mismatchScore, maxRepeatLength);
    }

    private static void ValidateScoredInvertedRepeatParameters(
        int gapPenalty, int threshold, int matchScore, int mismatchScore, int maxRepeatLength)
    {
        // Ranges of einverted.acd; scores must stay below the rogue sentinel (einverted's arrays assume it).
        ArgumentOutOfRangeException.ThrowIfNegative(gapPenalty);
        ArgumentOutOfRangeException.ThrowIfGreaterThanOrEqual(gapPenalty, EinvertedRogue);
        ArgumentOutOfRangeException.ThrowIfNegative(threshold);
        ArgumentOutOfRangeException.ThrowIfGreaterThanOrEqual(threshold, EinvertedRogue);
        ArgumentOutOfRangeException.ThrowIfNegative(matchScore);
        ArgumentOutOfRangeException.ThrowIfGreaterThanOrEqual(matchScore, EinvertedRogue);
        ArgumentOutOfRangeException.ThrowIfGreaterThan(mismatchScore, 0);
        ArgumentOutOfRangeException.ThrowIfLessThanOrEqual(mismatchScore, -EinvertedRogue);
        ArgumentOutOfRangeException.ThrowIfLessThan(maxRepeatLength, 2);
    }

    /// <summary>
    /// The einverted scan: a ring of <c>W</c> DP rows (row <c>i</c> in slot <c>i mod W</c>), per-row best scores and a
    /// left-start → right-end table, both ring-indexed. With <c>W ≥ n + 2</c> no slot is ever reused and every row is
    /// scored in full, so <c>W</c> is capped at <c>n + 2</c> without changing the output.
    /// </summary>
    private sealed class EinvertedScanner
    {
        private readonly string _seq;
        private readonly int[] _code;
        private readonly int _n, _w, _gap, _threshold, _match, _mismatch;
        private readonly int[][] _rows;
        private readonly int[] _bestEnd;   // einverted "back": slot of a left start → right end (row) of its best repeat
        private readonly int[] _rowBest;   // einverted "localMax": best score of a row (0 = none)
        private readonly List<ScoredInvertedRepeatResult> _results = new();

        public EinvertedScanner(string seq, int gap, int threshold, int match, int mismatch, int maxRepeat)
        {
            _seq = seq;
            _n = seq.Length;
            _w = (int)Math.Min(maxRepeat, (long)_n + 2);
            _gap = gap;
            _threshold = threshold;
            _match = match;
            _mismatch = mismatch;
            _code = new int[_n];
            for (int p = 0; p < _n; p++)
            {
                int c = AcgtCode(seq[p]);
                _code[p] = c < 0 ? 4 : c;
            }
            _rows = new int[_w][];
            for (int r = 0; r < _w; r++)
                _rows[r] = new int[_w];
            _bestEnd = new int[_w];
            _rowBest = new int[_w];
        }

        /// <summary>Pair score of right base <paramref name="i"/> with left base <c>i − 1 − k</c>; rogue off the sequence.</summary>
        private int PairScore(int i, int k)
        {
            int partner = i - 1 - k;
            if (partner < 0)
                return EinvertedRogue;
            int ci = _code[i];
            return ci < 4 && _code[partner] == 3 - ci ? _match : _mismatch;
        }

        public List<ScoredInvertedRepeatResult> Run()
        {
            int w = _w;
            int lastReported = -1;
            for (int i = 0; i < _n + w; i++)
            {
                int slot = i % w;
                if (_bestEnd[slot] != 0)
                    lastReported = ReportWindow(i, slot);

                if (i >= _n)
                    continue;

                int[] prev = _rows[(i + w - 1) % w];
                int[] cur = _rows[slot];
                for (int k = 0; k < w - 2; k++)
                    cur[k] = PairScore(i, k);
                cur[w - 2] = cur[w - 1] = EinvertedRogue;
                if (i == 0)
                    Array.Clear(prev); // row "−1" of einverted is the zero-initialised matrix

                int best = _threshold - 1, bestK = 0;
                int c = cur[0], diag = -EinvertedRogue;
                // Always ends at the rogue column (k ≤ W − 2 or the first column past the sequence start).
                for (int k = 1; k < w; k++)
                {
                    int d = cur[k];
                    if (diag > 0)
                        d += diag;
                    diag = prev[k - 1];
                    if (diag > c)
                        c = diag;
                    c -= _gap;
                    if (d > c)
                        c = d;
                    cur[k] = c;
                    if (c >= EinvertedRogue)
                        break;
                    if (c > best)
                    {
                        best = c;
                        bestK = k;
                    }
                }

                if (bestK != 0)
                {
                    int leftStart = i - bestK - 1;
                    int startSlot = leftStart % w;
                    int previousEnd = _bestEnd[startSlot];
                    if (leftStart > lastReported && (previousEnd == 0 || _rowBest[previousEnd % w] < best))
                    {
                        _bestEnd[startSlot] = i;
                        _rowBest[slot] = best;
                    }
                }
                else
                {
                    _rowBest[slot] = 0;
                }
            }
            return _results;
        }

        /// <summary>Reports the best repeat of the window (i − W, end] and clears it; returns the new last-reported end.</summary>
        private int ReportWindow(int i, int slot)
        {
            int w = _w;
            int windowEnd = _bestEnd[slot];
            int bestRow = windowEnd;
            bool done = false;
            while (!done)
            {
                int bestScore = 0;
                bestRow = windowEnd;
                for (int j = windowEnd; j > i - w; j--)
                {
                    int score = _rowBest[j % w];
                    if (score > bestScore)
                    {
                        bestRow = j;
                        bestScore = score;
                    }
                }
                if (bestScore == 0)
                    break;
                done = TraceBack(bestScore, bestRow, i - w);
                if (!done)
                    _rowBest[bestRow % w] = 0;
            }

            for (int j = bestRow; j >= i - w; j--)
            {
                _bestEnd[j % w] = 0;
                _rowBest[j % w] = 0;
            }
            return bestRow;
        }

        private bool TraceBack(int score, int endRow, int minRow)
        {
            int w = _w;
            int[] row = _rows[endRow % w];
            int k = 0;
            while (k < w && row[k] != score)
                k++;

            var right = new List<int>();
            var left = new List<int>();
            int remaining = score, i = endRow, matches = 0, mismatches = 0, gaps = 0;
            while (remaining > 0 && k >= 1)
            {
                if (i < minRow)
                    return false;
                right.Add(i);
                left.Add(i - 1 - k);
                if (row[k - 1] == remaining + _gap)
                {
                    remaining += _gap;
                    gaps++;
                    k--;
                    continue;
                }
                row = _rows[(i - 1 + w) % w];
                if (row[k - 1] == remaining + _gap)
                {
                    remaining += _gap;
                    gaps++;
                    i--;
                    k--;
                    continue;
                }
                int s = PairScore(i, k);
                remaining -= s;
                if (s == _match)
                    matches++;
                else
                    mismatches++;
                i--;
                k -= 2;
            }

            // einverted divides by (matches + mismatches) and prints coordinates here: with no aligned column, or a
            // column before the sequence start (possible only after a ring slot was reused), it aborts or prints
            // garbage. Such a trace-back still counts as reported (scan state as in einverted) but yields no result.
            if (left.Count == 0 || left[0] < 0)
                return true;

            int len = left.Count;
            var leftRow = new char[len];
            var midRow = new char[len];
            var rightRow = new char[len];
            for (int t = 0; t < len; t++)
            {
                bool leftGap = t + 1 < len && left[t] == left[t + 1];
                bool rightGap = t + 1 < len && right[t] == right[t + 1];
                leftRow[t] = leftGap ? '-' : _seq[left[t]];
                rightRow[t] = rightGap ? '-' : _seq[right[t]];
                midRow[t] = !leftGap && !rightGap && _code[right[t]] + _code[left[t]] == 3 ? '|' : ' ';
            }

            _results.Add(new ScoredInvertedRepeatResult(
                LeftArmStart: left[0],
                LeftArmEnd: left[len - 1],
                RightArmStart: right[len - 1],
                RightArmEnd: right[0],
                Score: score,
                Matches: matches,
                Mismatches: mismatches,
                Gaps: gaps,
                LeftArmAlignment: new string(leftRow),
                MatchLine: new string(midRow),
                RightArmAlignment: new string(rightRow)));
            return true;
        }
    }

    #endregion


    #region Direct Repeat Detection

    /// <summary>
    /// Finds exact direct repeats reported as <b>maximal repeated pairs</b> (Gusfield 1997, §7.12;
    /// Kurtz &amp; Schleiermacher 1999 REPuter forward repeats; MUMmer <c>repeat-match -f</c>,
    /// Kurtz et al. 2004).
    /// </summary>
    /// <remarks>
    /// <para>
    /// A pair of 0-based positions <c>i &lt; j</c> with length <c>L</c> is reported when
    /// <c>S[i..i+L) = S[j..j+L)</c>, the pair is <b>right-maximal</b> (<c>L</c> is the full length of the
    /// common prefix of the suffixes at <c>i</c> and <c>j</c>) and <b>left-maximal</b> (<c>i = 0</c> or
    /// <c>S[i−1] ≠ S[j−1]</c>). Every exact direct repeat is therefore reported exactly once, at its full
    /// extent — its shorter sub-copies (nested windows of the same pair of copies) are not repeated.
    /// This is exactly the forward-strand output of MUMmer <c>repeat-match -f -n minLength</c>
    /// (converted to 0-based positions), then filtered by <paramref name="maxLength"/> and
    /// <paramref name="minSpacing"/>.
    /// </para>
    /// <para>
    /// Only A/C/G/T match (case-insensitive); any other symbol (N, IUPAC codes, gaps) never matches and
    /// terminates a repeat — the MUMmer <c>mummer -n</c> convention ("match only the characters a, c, g,
    /// or t"), consistent with the ACGT-only rule of the other <see cref="RepeatFinder"/> methods.
    /// </para>
    /// <para>
    /// Filters: <c>minLength ≤ L ≤ maxLength</c> — a maximal repeat longer than <paramref name="maxLength"/>
    /// is <b>not</b> reported (it is not truncated into sub-windows; pass a larger <paramref name="maxLength"/>
    /// to see it); <c>Spacing = j − i − L ≥ minSpacing</c>. <paramref name="minSpacing"/> may be negative to
    /// admit overlapping copies (e.g. <c>int.MinValue</c> returns every maximal pair, as repeat-match does);
    /// <c>0</c> admits abutting (tandem) copies.
    /// </para>
    /// <para>
    /// Algorithm: suffix array + Kasai LCP array (shared with <see cref="SequenceComplexity"/>), bottom-up
    /// traversal of the lcp-interval tree with per-left-character position lists (Gusfield 1997 §7.12.3;
    /// Abouelhoda, Kurtz &amp; Ohlebusch 2004, maximal repeated pairs on enhanced suffix arrays):
    /// O(n log² n + z) time for z maximal pairs of length within [minLength, maxLength], O(n) extra space.
    /// Results are ordered by (FirstPosition, SecondPosition); each position pair occurs at most once.
    /// </para>
    /// </remarks>
    /// <param name="sequence">DNA sequence to search.</param>
    /// <param name="minLength">Minimum repeat length (default: 5, must be ≥ 2).</param>
    /// <param name="maxLength">Maximum repeat length (default: 50, must be ≥ <paramref name="minLength"/>).</param>
    /// <param name="minSpacing">Minimum number of bases between the copies (default: 1; negative admits overlap).</param>
    /// <returns>Maximal direct-repeat pairs, sorted by (FirstPosition, SecondPosition).</returns>
    public static IEnumerable<DirectRepeatResult> FindDirectRepeats(
        DnaSequence sequence,
        int minLength = 5,
        int maxLength = 50,
        int minSpacing = 1)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentOutOfRangeException.ThrowIfLessThan(minLength, 2);
        ArgumentOutOfRangeException.ThrowIfLessThan(maxLength, minLength);

        return FindDirectRepeatsCore(sequence.Sequence, minLength, maxLength, minSpacing);
    }

    /// <summary>
    /// Finds maximal exact direct-repeat pairs in a raw sequence string (case-insensitive; non-ACGT symbols
    /// never match). <c>null</c> or empty input yields no results. See
    /// <see cref="FindDirectRepeats(DnaSequence,int,int,int)"/> for the reporting convention.
    /// </summary>
    public static IEnumerable<DirectRepeatResult> FindDirectRepeats(
        string sequence,
        int minLength = 5,
        int maxLength = 50,
        int minSpacing = 1)
    {
        // Same numeric validation as the DnaSequence overload, eager (at the call, not on enumeration).
        ArgumentOutOfRangeException.ThrowIfLessThan(minLength, 2);
        ArgumentOutOfRangeException.ThrowIfLessThan(maxLength, minLength);

        if (string.IsNullOrEmpty(sequence))
            return Array.Empty<DirectRepeatResult>();

        return FindDirectRepeatsCore(sequence.ToUpperInvariant(), minLength, maxLength, minSpacing);
    }

    /// <summary>Left-character class for a suffix start with no matchable left neighbour (p = 0 or non-ACGT).</summary>
    private const int UniqueLeftClass = 4;
    private const int LeftClassCount = 5;

    /// <summary>A/C/G/T → 0..3, other symbols → −1 (shared <see cref="SequenceComplexity.AcgtCode"/>; inputs here are upper-cased).</summary>
    private static int AcgtCode(char c) => SequenceComplexity.AcgtCode(c);

    private static List<DirectRepeatResult> FindDirectRepeatsCore(
        string seq,
        int minLength,
        int maxLength,
        int minSpacing)
    {
        var results = new List<DirectRepeatResult>();
        int n = seq.Length;
        if (n <= minLength)
            return results;

        // Symbols: A/C/G/T → 0..3; every other symbol gets a unique code 4 + p, so it never matches.
        var symbols = new int[n];
        var leftClass = new int[n];
        for (int p = 0; p < n; p++)
        {
            int code = AcgtCode(seq[p]);
            symbols[p] = code >= 0 ? code : 4 + p;
            leftClass[p] = UniqueLeftClass;
        }
        for (int p = 1; p < n; p++)
        {
            if (symbols[p - 1] < 4)
                leftClass[p] = symbols[p - 1];
        }

        int[] sa = SequenceComplexity.BuildSuffixArray(symbols);
        int[] lcp = SequenceComplexity.BuildLcpArray(symbols, sa);

        // Bottom-up lcp-interval traversal. Each interval keeps its suffix positions in one linked list
        // per left-character class; merging a child into its parent at string depth ℓ emits every pair
        // (p from the child, q already in the parent) whose left characters differ (or are undefined):
        // those pairs have LCP exactly ℓ (right-maximal) and are left-maximal.
        var next = new int[n];
        var stackLcp = new int[n + 1];
        var stackHead = new int[(n + 1) * LeftClassCount];
        var stackTail = new int[(n + 1) * LeftClassCount];
        var childHead = new int[LeftClassCount];
        var childTail = new int[LeftClassCount];
        int top = 0;
        stackLcp[0] = 0;
        Array.Fill(stackHead, -1, 0, LeftClassCount);

        for (int r = 1; r <= n; r++)
        {
            // Pending child: the leaf for suffix sa[r − 1].
            Array.Fill(childHead, -1);
            int leaf = sa[r - 1];
            next[leaf] = -1;
            childHead[leftClass[leaf]] = leaf;
            childTail[leftClass[leaf]] = leaf;

            int h = r < n ? lcp[r] : 0;
            while (stackLcp[top] > h)
            {
                MergeDirectRepeatLists(seq, top, childHead, childTail, stackLcp, stackHead, stackTail, next,
                    minLength, maxLength, minSpacing, results);
                for (int c = 0; c < LeftClassCount; c++)
                {
                    childHead[c] = stackHead[top * LeftClassCount + c];
                    childTail[c] = stackTail[top * LeftClassCount + c];
                }
                top--;
            }

            if (stackLcp[top] < h)
            {
                top++;
                stackLcp[top] = h;
                Array.Fill(stackHead, -1, top * LeftClassCount, LeftClassCount);
            }

            MergeDirectRepeatLists(seq, top, childHead, childTail, stackLcp, stackHead, stackTail, next,
                minLength, maxLength, minSpacing, results);
        }

        results.Sort(static (a, b) =>
        {
            int c = a.FirstPosition.CompareTo(b.FirstPosition);
            return c != 0 ? c : a.SecondPosition.CompareTo(b.SecondPosition);
        });
        return results;
    }

    /// <summary>
    /// Emits the maximal pairs between a child interval's lists and the lists already accumulated in the
    /// interval at stack slot <paramref name="node"/>, then concatenates the child lists into the node.
    /// </summary>
    private static void MergeDirectRepeatLists(
        string seq,
        int node,
        int[] childHead,
        int[] childTail,
        int[] stackLcp,
        int[] stackHead,
        int[] stackTail,
        int[] next,
        int minLength,
        int maxLength,
        int minSpacing,
        List<DirectRepeatResult> results)
    {
        int length = stackLcp[node];
        int baseIdx = node * LeftClassCount;

        if (length >= minLength && length <= maxLength)
        {
            for (int a = 0; a < LeftClassCount; a++)
            {
                if (childHead[a] < 0) continue;
                for (int b = 0; b < LeftClassCount; b++)
                {
                    if (a == b && a != UniqueLeftClass) continue; // same left character: not left-maximal
                    int nodeListHead = stackHead[baseIdx + b];
                    if (nodeListHead < 0) continue;

                    for (int p = childHead[a]; p >= 0; p = next[p])
                    {
                        for (int q = nodeListHead; q >= 0; q = next[q])
                        {
                            int i = Math.Min(p, q);
                            int j = Math.Max(p, q);
                            long spacing = (long)j - i - length;
                            if (spacing < minSpacing) continue;
                            results.Add(new DirectRepeatResult(
                                FirstPosition: i,
                                SecondPosition: j,
                                RepeatSequence: seq.Substring(i, length),
                                Length: length,
                                Spacing: (int)spacing));
                        }
                    }
                }
            }
        }

        for (int c = 0; c < LeftClassCount; c++)
        {
            if (childHead[c] < 0) continue;
            int idx = baseIdx + c;
            if (stackHead[idx] < 0)
            {
                stackHead[idx] = childHead[c];
            }
            else
            {
                next[stackTail[idx]] = childHead[c];
            }
            stackTail[idx] = childTail[c];
        }
    }

    #endregion

    #region Tandem Repeat Summary

    /// <summary>
    /// Gets a summary of all perfect microsatellites (1–6 bp units) in a sequence.
    /// </summary>
    /// <remarks>
    /// <para>
    /// The summary aggregates exactly the list returned by
    /// <see cref="FindMicrosatellites(DnaSequence,int,int,int)"/> with unit lengths 1–6 and the given
    /// <paramref name="minRepeats"/> (one maximal primitive run per unit length; no second scan), following the
    /// statistics conventions of MISA (Thiel et al. 2003, <c>misa.pl</c> <c>.statistics</c> output) and Krait
    /// (Du et al. 2018, <c>statistics.py</c>):
    /// </para>
    /// <list type="bullet">
    /// <item><description><c>TotalRepeats</c> — number of reported microsatellites (MISA "Total number of identified SSRs").</description></item>
    /// <item><description>Per-class counts for all six unit sizes, mono … hexa (MISA "Distribution to different repeat
    /// type classes": one row per unit size; Krait <c>motifTypeStatis</c>: Mono, Di, Tri, Tetra, Penta, Hexa). The six
    /// counts always sum to <c>TotalRepeats</c>.</description></item>
    /// <item><description><c>TotalRepeatBases</c> — sum of the repeat lengths (Krait "Length (bp)" = <c>SUM(length)</c>);
    /// runs of different unit lengths that overlap are each counted in full.</description></item>
    /// <item><description><c>PercentageOfSequence</c> — percent of the sequence's bases covered by at least one reported
    /// microsatellite: |∪ [Position, Position+TotalLength)| / sequence length × 100, so always in [0, 100]. When no two
    /// reported runs overlap this equals <c>TotalRepeatBases</c> / length × 100 (Krait's relative density expressed as a
    /// percentage; DnaSequence contains only A/C/G/T, so Krait's valid-base denominator equals the length).</description></item>
    /// <item><description><c>LongestRepeat</c> — the run with the largest <c>TotalLength</c>; ties go to the shorter unit,
    /// then the leftmost position (the <see cref="FindMicrosatellites(DnaSequence,int,int,int)"/> order).</description></item>
    /// <item><description><c>MostFrequentUnit</c> — the repeat unit string (as reported, i.e. the motif phase at the run
    /// start, not rotation- or strand-canonicalized; MISA table "Frequency of identified SSR motifs") that occurs in the
    /// most runs; ties go to the unit whose first run appears first in that order.</description></item>
    /// </list>
    /// </remarks>
    /// <param name="sequence">DNA sequence to summarize.</param>
    /// <param name="minRepeats">Minimum number of complete copies (≥ 2; default 3), applied to every unit length.</param>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> is <c>null</c>.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="minRepeats"/> is less than 2.</exception>
    public static TandemRepeatSummary GetTandemRepeatSummary(
        DnaSequence sequence,
        int minRepeats = 3)
    {
        ArgumentNullException.ThrowIfNull(sequence);

        var microsatellites = FindMicrosatellites(sequence, 1, 6, minRepeats).ToList();

        var byType = microsatellites
            .GroupBy(m => m.RepeatType)
            .ToDictionary(g => g.Key, g => g.ToList());

        int totalBases = microsatellites.Sum(m => m.TotalLength);

        // PercentageOfSequence is the fraction of the sequence that is tandem-repeat, so it must use the
        // DISTINCT bases covered by any microsatellite (the union of their spans), not the raw sum of repeat
        // lengths: overlapping repeats (e.g. a homopolymer run also matched as a dinucleotide repeat) would
        // otherwise double-count bases and push the percentage above 100. Covered bases ≤ sequence length, so
        // the percentage is always in [0, 100]. TotalRepeatBases keeps the (possibly overlapping) repeat content.
        long coveredBases = CountCoveredBases(microsatellites, sequence.Length);
        double percentageOfSequence = sequence.Length > 0
            ? (double)coveredBases / sequence.Length * 100
            : 0;

        return new TandemRepeatSummary(
            TotalRepeats: microsatellites.Count,
            TotalRepeatBases: totalBases,
            PercentageOfSequence: percentageOfSequence,
            MononucleotideRepeats: byType.GetValueOrDefault(RepeatType.Mononucleotide)?.Count ?? 0,
            DinucleotideRepeats: byType.GetValueOrDefault(RepeatType.Dinucleotide)?.Count ?? 0,
            TrinucleotideRepeats: byType.GetValueOrDefault(RepeatType.Trinucleotide)?.Count ?? 0,
            TetranucleotideRepeats: byType.GetValueOrDefault(RepeatType.Tetranucleotide)?.Count ?? 0,
            PentanucleotideRepeats: byType.GetValueOrDefault(RepeatType.Pentanucleotide)?.Count ?? 0,
            HexanucleotideRepeats: byType.GetValueOrDefault(RepeatType.Hexanucleotide)?.Count ?? 0,
            // MicrosatelliteResult is a struct: FirstOrDefault() on the struct sequence would return default(...)
            // (Position 0, RepeatUnit null) instead of null when nothing is found, so lift to the nullable type first.
            LongestRepeat: microsatellites
                .OrderByDescending(m => m.TotalLength)
                .Select(m => (MicrosatelliteResult?)m)
                .FirstOrDefault(),
            MostFrequentUnit: microsatellites
                .GroupBy(m => m.RepeatUnit)
                .OrderByDescending(g => g.Count())
                .FirstOrDefault()?.Key);
    }

    /// <summary>
    /// Counts the number of distinct sequence bases covered by at least one microsatellite, i.e. the length
    /// of the union of the half-open spans <c>[Position, Position + TotalLength)</c>. Spans are clamped to the
    /// sequence length and merged in start order, so the result never exceeds <paramref name="sequenceLength"/>.
    /// </summary>
    private static long CountCoveredBases(
        IReadOnlyList<MicrosatelliteResult> microsatellites, int sequenceLength)
    {
        if (microsatellites.Count == 0)
        {
            return 0;
        }

        var intervals = microsatellites
            .Select(m => (Start: m.Position, End: Math.Min(sequenceLength, m.Position + m.TotalLength)))
            .Where(iv => iv.End > iv.Start)
            .OrderBy(iv => iv.Start)
            .ToList();

        long covered = 0;
        int currentStart = -1;
        int currentEnd = -1;
        foreach (var iv in intervals)
        {
            if (iv.Start > currentEnd)
            {
                // Disjoint from the current run: close it out and start a new one.
                covered += currentEnd - currentStart;
                currentStart = iv.Start;
                currentEnd = iv.End;
            }
            else
            {
                // Overlapping or adjacent: extend the current run.
                currentEnd = Math.Max(currentEnd, iv.End);
            }
        }

        covered += currentEnd - currentStart;
        return covered;
    }

    #endregion

    #region Palindrome Detection

    /// <summary>
    /// Finds reverse-complement palindromes (sequences that read the same 5'→3' on both strands, i.e. equal
    /// their own reverse complement — the recognition-site form of most Type II restriction enzymes).
    /// </summary>
    /// <remarks>
    /// <para>
    /// Reports <b>every</b> window <c>S[i..i+len)</c> with <c>minLength ≤ len ≤ maxLength</c> that equals its
    /// reverse complement, including windows nested inside longer palindromes and overlapping ones — the
    /// Rosalind REVP ("Locating Restriction Sites") convention: "the position and length of every reverse
    /// palindrome in the string having length between 4 and 12". This differs from
    /// <see cref="FindInvertedRepeats(DnaSequence,int,int,int,int,int,bool)"/>, which reports only maximal stems (EMBOSS
    /// <c>palindrome</c>).
    /// </para>
    /// <para>
    /// A reverse-complement palindrome always has even length (the middle base of an odd window would have to be
    /// its own complement), so only even lengths are scanned; an odd <paramref name="maxLength"/> is effectively
    /// rounded down. Only the unambiguous bases A, C, G, T pair (canonical
    /// <see cref="SequenceExtensions.GetComplementBase(char)"/>): a window containing N, another IUPAC ambiguity
    /// code, U, a gap or any other symbol is never reported (as in the REVP alphabet and EMBOSS
    /// <c>einverted</c>; EMBOSS <c>palindrome</c> likewise rejects all-N stems) — e.g. <c>NNNN</c> or <c>SSSS</c>
    /// equal their symbolic reverse complement but need not be palindromic once resolved (<c>SS</c> may be CC).
    /// </para>
    /// <para>
    /// Positions are 0-based (REVP positions are 1-based: add 1). Results are ordered by <c>Position</c>,
    /// then <c>Length</c> (the order of the REVP sample output). Parameters are validated eagerly on both
    /// overloads. Cost O(n·maxLength) comparisons (one capped centre expansion per inter-base gap) plus output.
    /// </para>
    /// </remarks>
    /// <param name="sequence">DNA sequence to search.</param>
    /// <param name="minLength">Minimum palindrome length (default: 4; must be even and ≥ 4).</param>
    /// <param name="maxLength">Maximum palindrome length (default: 12; must be ≥ <paramref name="minLength"/>).</param>
    /// <returns>Every palindromic window, ordered by position then length.</returns>
    public static IEnumerable<PalindromeResult> FindPalindromes(
        DnaSequence sequence,
        int minLength = 4,
        int maxLength = 12)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ValidatePalindromeParameters(minLength, maxLength);

        return FindPalindromesCore(sequence.Sequence, minLength, maxLength);
    }

    /// <summary>
    /// Finds reverse-complement palindromes in a raw sequence string (case-insensitive).
    /// Same model and validation as <see cref="FindPalindromes(DnaSequence, int, int)"/>; <c>null</c> or empty
    /// input yields no results. Characters other than A/C/G/T never pair.
    /// </summary>
    public static IEnumerable<PalindromeResult> FindPalindromes(
        string sequence,
        int minLength = 4,
        int maxLength = 12)
    {
        ValidatePalindromeParameters(minLength, maxLength);

        return string.IsNullOrEmpty(sequence)
            ? Array.Empty<PalindromeResult>()
            : FindPalindromesCore(sequence.ToUpperInvariant(), minLength, maxLength);
    }

    private static void ValidatePalindromeParameters(int minLength, int maxLength)
    {
        if (minLength < 4 || minLength % 2 != 0)
            throw new ArgumentOutOfRangeException(nameof(minLength), minLength, "Must be even and >= 4");
        ArgumentOutOfRangeException.ThrowIfLessThan(maxLength, minLength);
    }

    private static IEnumerable<PalindromeResult> FindPalindromesCore(
        string seq,
        int minLength,
        int maxLength)
    {
        int n = seq.Length;
        int minHalf = minLength / 2;
        int maxHalf = Math.Min(maxLength, n) / 2; // capping by n also avoids int overflow for huge maxLength
        if (maxHalf < minHalf)
            yield break;

        // radius[c] = number of consecutive Watson–Crick pairs (seq[c−1−k], seq[c+k]) around the gap before c,
        // capped at maxHalf. The window (i, 2h) is a palindrome iff radius[i + h] ≥ h.
        var radius = new int[n];
        for (int c = 1; c < n; c++)
        {
            int r = 0;
            while (r < maxHalf && c - r - 1 >= 0 && c + r < n && IsWatsonCrickPair(seq[c - r - 1], seq[c + r]))
                r++;
            radius[c] = r;
        }

        for (int i = 0; i + 2 * minHalf <= n; i++)
        {
            for (int h = minHalf; h <= maxHalf && i + 2 * h <= n; h++)
            {
                if (radius[i + h] >= h)
                    yield return new PalindromeResult(Position: i, Sequence: seq.Substring(i, 2 * h), Length: 2 * h);
            }
        }
    }

    #endregion
}

/// <summary>
/// Type of repeat unit.
/// </summary>
public enum RepeatType
{
    Mononucleotide,  // A, T, G, C
    Dinucleotide,    // AT, GC, CA, etc.
    Trinucleotide,   // CAG, CGG, etc.
    Tetranucleotide, // GATA, AAAT, etc.
    Pentanucleotide, // AAAAT, etc.
    Hexanucleotide,  // AAAAAG, etc.
    Complex          // Longer than 6bp
}

/// <summary>
/// Result of microsatellite (STR) detection.
/// </summary>
public readonly record struct MicrosatelliteResult(
    int Position,
    string RepeatUnit,
    int RepeatCount,
    int TotalLength,
    RepeatType RepeatType)
{
    /// <summary>
    /// Gets the full repeat sequence.
    /// </summary>
    public string FullSequence => string.Concat(Enumerable.Repeat(RepeatUnit, RepeatCount));
}

/// <summary>
/// Result of approximate (imperfect / interrupted) tandem-repeat detection with the Tandem Repeats Finder
/// model (Benson 1999). Fields mirror the TRF table: <see cref="Start"/> (0-based; TRF prints 1-based) and
/// <see cref="SpanLength"/> give the indices; <see cref="Period"/> is the most common distance between matching
/// characters of adjacent copies (may differ from <see cref="ConsensusSize"/>); <see cref="CopyNumber"/> is the
/// number of copies aligned with the consensus; <see cref="PercentMatches"/> / <see cref="PercentIndels"/> are
/// between ADJACENT copies (exact values; TRF truncates to integers); <see cref="AlignmentScore"/> is the WDP
/// score; <see cref="Consensus"/> starts at the phase of the first repeat base.
/// </summary>
public readonly record struct ApproximateTandemRepeatResult(
    int Start,
    int SpanLength,
    int Period,
    int ConsensusSize,
    string Consensus,
    double CopyNumber,
    double PercentMatches,
    double PercentIndels,
    int AlignmentScore)
{
    /// <summary>Percentage of A in the repeat region (TRF "A" column, exact; denominator = region length).</summary>
    public double PercentA { get; init; }

    /// <summary>Percentage of C in the repeat region (TRF "C" column, exact).</summary>
    public double PercentC { get; init; }

    /// <summary>Percentage of G in the repeat region (TRF "G" column, exact).</summary>
    public double PercentG { get; init; }

    /// <summary>Percentage of T in the repeat region (TRF "T" column, exact).</summary>
    public double PercentT { get; init; }

    /// <summary>Shannon entropy of the region's A/C/G/T composition in bits, 0–2 (TRF "Entropy (0-2)").</summary>
    public double Entropy { get; init; }
}

/// <summary>
/// Tandem Repeats Finder Bernoulli-model statistical measures (Benson 1999) for a detected repeat tract.
/// Benson models the alignment of two adjacent copies as a sequence of independent Bernoulli trials;
/// <see cref="MatchProbability"/> = P(Heads) = PM (matching probability) is the average percent identity
/// between adjacent copies, and <see cref="IndelProbability"/> = PI (indel probability) is the average
/// percentage of insertions and deletions between them. The statistics are computed between ADJACENT
/// copies, not between the tract and the consensus pattern.
/// </summary>
public readonly record struct TandemRepeatBernoulliStatistics(
    int Period,
    int AdjacentCopyPairs,
    int BernoulliTrials,
    int Matches,
    int Mismatches,
    int Indels,
    double MatchProbability,
    double IndelProbability,
    double PercentMatches,
    double PercentIndels,
    double ExpectedMatches,
    bool MeetsExpectedMatchProbability);

/// <summary>
/// Result of inverted repeat detection.
/// </summary>
public readonly record struct InvertedRepeatResult(
    int LeftArmStart,
    int RightArmStart,
    int ArmLength,
    int LoopLength,
    string LeftArm,
    string RightArm,
    string Loop,
    bool CanFormHairpin)
{
    /// <summary>
    /// Total length of the inverted repeat structure.
    /// </summary>
    public int TotalLength => 2 * ArmLength + LoopLength;

    /// <summary>
    /// Number of mismatched (non-pairing) positions inside the stem; always 0 for exact stems
    /// (<c>maxMismatches = 0</c>). Mismatches never occur at either end of a stem.
    /// </summary>
    public int Mismatches { get; init; }
}

/// <summary>
/// A scored imperfect inverted repeat (EMBOSS <c>einverted</c> semantics). Coordinates are 0-based and inclusive:
/// <c>LeftArmStart ≤ LeftArmEnd &lt; RightArmStart ≤ RightArmEnd</c>. The three alignment rows are column-aligned:
/// <see cref="LeftArmAlignment"/> reads the left arm 5'→3' (outer → inner), <see cref="RightArmAlignment"/> reads the
/// right arm from its outer end inward (3'→5'), '-' marks a gap, and <see cref="MatchLine"/> has '|' for each
/// Watson–Crick pair. <see cref="Score"/>, <see cref="Matches"/>, <see cref="Mismatches"/> and <see cref="Gaps"/> are
/// einverted's figures: a final trace-back gap step is counted in <see cref="Gaps"/> but not drawn, and an innermost
/// loop-0 pair reached by the last diagonal step contributes to <see cref="Score"/> without being counted or drawn.
/// </summary>
public readonly record struct ScoredInvertedRepeatResult(
    int LeftArmStart,
    int LeftArmEnd,
    int RightArmStart,
    int RightArmEnd,
    int Score,
    int Matches,
    int Mismatches,
    int Gaps,
    string LeftArmAlignment,
    string MatchLine,
    string RightArmAlignment)
{
    /// <summary>Length of the left arm (bases).</summary>
    public int LeftArmLength => LeftArmEnd - LeftArmStart + 1;

    /// <summary>Length of the right arm (bases).</summary>
    public int RightArmLength => RightArmEnd - RightArmStart + 1;

    /// <summary>Bases between the arms.</summary>
    public int LoopLength => RightArmStart - LeftArmEnd - 1;

    /// <summary>Percentage of aligned (non-gap) columns that pair, 100·matches / (matches + mismatches); 0 when none.</summary>
    public double PercentMatches => Matches + Mismatches == 0 ? 0 : 100.0 * Matches / (Matches + Mismatches);
}

/// <summary>
/// Result of direct repeat detection.
/// </summary>
public readonly record struct DirectRepeatResult(
    int FirstPosition,
    int SecondPosition,
    string RepeatSequence,
    int Length,
    int Spacing);

/// <summary>
/// Result of palindrome detection.
/// </summary>
public readonly record struct PalindromeResult(
    int Position,
    string Sequence,
    int Length);

/// <summary>
/// Summary of perfect microsatellites (1–6 bp units) in a sequence; see
/// <see cref="RepeatFinder.GetTandemRepeatSummary(DnaSequence,int)"/> for each field's definition.
/// The six per-class counts (mono … hexa) sum to <see cref="TotalRepeats"/>.
/// </summary>
public readonly record struct TandemRepeatSummary(
    int TotalRepeats,
    int TotalRepeatBases,
    double PercentageOfSequence,
    int MononucleotideRepeats,
    int DinucleotideRepeats,
    int TrinucleotideRepeats,
    int TetranucleotideRepeats,
    int PentanucleotideRepeats,
    int HexanucleotideRepeats,
    MicrosatelliteResult? LongestRepeat,
    string? MostFrequentUnit);
