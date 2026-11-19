namespace Seqeron.Genomics.Analysis;

/// <summary>
/// All <b>maximal repetitions</b> ("runs") of a string in O(n log² n) time (suffix-array construction dominates; the
/// scan itself is O(n) longest-common-extension queries), independent of the period range — the engine behind
/// <see cref="RepeatFinder.FindMicrosatellites(DnaSequence,int,int,int)"/> (and so of <see cref="GenomicAnalyzer.FindTandemRepeats(DnaSequence,int,int)"/>)
/// when the unit-length range is too wide for one linear scan per unit length (more than <c>RepeatFinder.RunsEngineMinUnitLengths</c> admissible unit lengths).
/// </summary>
/// <remarks>
/// <para>
/// <b>Definition</b> (Kolpakov &amp; Kucherov 1999): a run is a triple <c>(a, e, p)</c> — an interval <c>S[a..e)</c>
/// whose <i>smallest</i> period is <c>p</c>, with <c>e − a ≥ 2p</c>, that cannot be extended to the left or right
/// without changing its smallest period (equivalently: maximal with period p). A string has fewer than n runs
/// (the "runs theorem", Bannai, I, Inenaga, Nakashima, Takeda &amp; Tsuruta 2017, SIAM J. Comput. 46:1501).
/// </para>
/// <para>
/// <b>Algorithm</b> (the constructive proof of the runs theorem, arXiv:1406.0263 sections 3 and 4): for each of the two
/// alphabet orders ≺0 (natural) and ≺1 (reversed), the string is read as w$ with the end marker $ smaller than every
/// letter under ≺0 and therefore larger than every letter under ≺1 (a proper prefix is smaller under ≺0, larger under
/// ≺1); let
/// <c>λ(i)</c> = the longest Lyndon word starting at <c>i</c> = <c>S[i..j)</c> where <c>j</c> is the next position with
/// <c>S[j..] ≺ S[i..]</c> (next-smaller-suffix, one monotonic-stack pass). Every run <c>(a, e, p)</c> has, for one of
/// the two orders, a Lyndon root <c>S[i..i+p)</c> equal to <c>λ(i)</c> with <c>i &gt; a</c>. So for every (order, i) take
/// <c>p = |λ(i)|</c>, extend the period <c>p</c> to the right with the longest common extension
/// <c>LCE(i, i+p)</c> and to the left with the longest common suffix of the prefixes ending at <c>i−1</c> and
/// <c>i+p−1</c>, and keep the candidate when it spans at least 2p. A candidate is always a genuine run (a Lyndon word is
/// primitive, so the smallest period of an interval of length ≥ 2p with period p is p), and the set union over both
/// orders contains every run; duplicates are removed. LCE queries are answered in O(1) from the suffix array of the text
/// and of its reverse (prefix doubling, <see cref="SequenceComplexity.BuildSuffixArray"/>), the Kasai LCP array
/// (<see cref="SequenceComplexity.BuildLcpArray"/>) and a block-decomposed range-minimum structure.
/// </para>
/// </remarks>
internal static class MaximalRepetitions
{
    /// <summary>A maximal repetition: <c>S[Start..End)</c> has smallest period <see cref="Period"/> and length ≥ 2·Period.</summary>
    internal readonly record struct Run(int Start, int End, int Period);

    /// <summary>All runs of <paramref name="s"/>, ordered by period, then start (each run once).</summary>
    internal static List<Run> FindRuns(string s, CancellationToken cancellationToken = default)
    {
        int n = s.Length;
        var runs = new List<Run>();
        if (n < 2)
            return runs;

        cancellationToken.ThrowIfCancellationRequested();

        var forward = new LceIndex(Encode(s, reverse: false));
        cancellationToken.ThrowIfCancellationRequested();
        var backward = new LceIndex(Encode(s, reverse: true));
        cancellationToken.ThrowIfCancellationRequested();
        var seen = new HashSet<(int, int, int)>();

        for (int order = 0; order < 2; order++)
        {
            // next-smaller-suffix (w.r.t. the order) by a monotonic stack: nsv[i] = j, the end of lambda(i).
            int noSmaller = order == 0 ? n : n + 1; // the end marker is the smaller suffix under order 0 only
            var nsv = new int[n];
            var stack = new int[n];
            int top = -1;
            for (int i = n - 1; i >= 0; i--)
            {
                if ((i & 4095) == 0)
                    cancellationToken.ThrowIfCancellationRequested();
                while (top >= 0 && !SuffixSmaller(s, forward, stack[top], i, order))
                    top--;
                nsv[i] = top >= 0 ? stack[top] : noSmaller;
                stack[++top] = i;
            }

            for (int i = 0; i < n; i++)
            {
                if ((i & 4095) == 0)
                    cancellationToken.ThrowIfCancellationRequested();
                int j = nsv[i];
                int p = j - i;
                // j == n (order 0: lambda(i) reaches the end of the text) is a valid last root: the run then only extends
                // left. j == n + 1 (order 1) means lambda(i) contains the end marker, which no run root can.
                if (j > n)
                    continue;

                int right = j + forward.Lce(i, j);                            // first mismatch of period p, exclusive end
                int leftExt = i > 0 ? backward.Lce(n - i, n - j) : 0;         // common suffix of prefixes ending at i-1 / j-1
                int left = i - leftExt;
                if (right - left >= 2 * p && seen.Add((left, right, p)))
                    runs.Add(new Run(left, right, p));
            }
        }

        runs.Sort((x, y) => x.Period != y.Period ? x.Period.CompareTo(y.Period) : x.Start.CompareTo(y.Start));
        return runs;
    }

    /// <summary>True when suffix S[a..] ≺ S[b..] (a ≠ b) in the given order; end of string = $: a proper prefix is smaller under ≺0 and larger under ≺1.</summary>
    private static bool SuffixSmaller(string s, LceIndex index, int a, int b, int order)
    {
        int n = s.Length;
        int l = index.Lce(a, b);
        if (a + l == n) return order == 0; // S[a..] is a proper prefix of S[b..]
        if (b + l == n) return order == 1; // S[b..] is a proper prefix of S[a..]
        return order == 0 ? s[a + l] < s[b + l] : s[a + l] > s[b + l];
    }

    private static int[] Encode(string s, bool reverse)
    {
        int n = s.Length;
        var symbols = new int[n];
        for (int i = 0; i < n; i++)
            symbols[i] = reverse ? s[n - 1 - i] : s[i];
        return symbols;
    }

    /// <summary>Longest common extension of two suffixes in O(1): suffix array + Kasai LCP + block sparse-table RMQ.</summary>
    private sealed class LceIndex
    {
        private const int Block = 16;
        private readonly int _n;
        private readonly int[] _rank;
        private readonly int[] _lcp;
        private readonly int[][] _sparse; // sparse[k][b] = min of block minima b .. b + 2^k - 1

        internal LceIndex(int[] symbols)
        {
            _n = symbols.Length;
            int[] sa = SequenceComplexity.BuildSuffixArray(symbols);
            _lcp = SequenceComplexity.BuildLcpArray(symbols, sa);
            _rank = new int[_n];
            for (int r = 0; r < _n; r++) _rank[sa[r]] = r;

            int blocks = (_n + Block - 1) / Block;
            var blockMin = new int[blocks];
            for (int b = 0; b < blocks; b++)
            {
                int min = int.MaxValue;
                for (int x = b * Block, end = Math.Min(_n, x + Block); x < end; x++)
                    if (_lcp[x] < min) min = _lcp[x];
                blockMin[b] = min;
            }

            int levels = 1;
            while ((1 << levels) <= blocks) levels++;
            _sparse = new int[levels][];
            _sparse[0] = blockMin;
            for (int k = 1; k < levels; k++)
            {
                int len = blocks - (1 << k) + 1;
                if (len <= 0) { _sparse[k] = []; continue; }
                var prev = _sparse[k - 1];
                var cur = new int[len];
                for (int b = 0; b < len; b++)
                    cur[b] = Math.Min(prev[b], prev[b + (1 << (k - 1))]);
                _sparse[k] = cur;
            }
        }

        /// <summary>Length of the longest common prefix of the suffixes at <paramref name="i"/> and <paramref name="j"/> (i ≠ j).</summary>
        internal int Lce(int i, int j)
        {
            if (i >= _n || j >= _n) return 0;
            int a = _rank[i], b = _rank[j];
            if (a > b) (a, b) = (b, a);
            return RangeMin(a + 1, b);
        }

        private int RangeMin(int l, int r) // inclusive, l <= r
        {
            int bl = l / Block, br = r / Block;
            int min = int.MaxValue;
            if (bl == br)
            {
                for (int x = l; x <= r; x++)
                    if (_lcp[x] < min) min = _lcp[x];
                return min;
            }

            for (int x = l, end = (bl + 1) * Block; x < end; x++)
                if (_lcp[x] < min) min = _lcp[x];
            for (int x = br * Block; x <= r; x++)
                if (_lcp[x] < min) min = _lcp[x];
            if (bl + 1 <= br - 1)
            {
                int from = bl + 1, to = br - 1;
                int k = 31 - System.Numerics.BitOperations.LeadingZeroCount((uint)(to - from + 1));
                int m = Math.Min(_sparse[k][from], _sparse[k][to - (1 << k) + 1]);
                if (m < min) min = m;
            }

            return min;
        }
    }
}
