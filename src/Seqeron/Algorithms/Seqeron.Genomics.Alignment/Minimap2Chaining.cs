namespace Seqeron.Genomics.Alignment;

/// <summary>A seed hit (anchor) of minimap2: a minimizer of the query that occurs in a target sequence.</summary>
/// <param name="TargetId">Index of the target (reference) sequence.</param>
/// <param name="Reverse">The hit is on the reverse strand of the query.</param>
/// <param name="TargetEnd">0-based position of the LAST base of the seed on the target.</param>
/// <param name="QueryEnd">0-based position of the last base of the seed on the query (reverse strand: on the reverse-complemented query).</param>
/// <param name="Span">Length of the seed (the k-mer span).</param>
public readonly record struct MinimapAnchor(int TargetId, bool Reverse, int TargetEnd, int QueryEnd, int Span)
{
    /// <summary><c>mm128_t.x</c>: <c>rev &lt;&lt; 63 | tid &lt;&lt; 32 | tpos</c>.</summary>
    internal ulong X => (Reverse ? 1UL << 63 : 0UL) | ((ulong)(uint)TargetId << 32) | (uint)TargetEnd;

    /// <summary><c>mm128_t.y</c>: <c>span &lt;&lt; 32 | qpos</c>.</summary>
    internal ulong Y => ((ulong)(uint)Span << 32) | (uint)QueryEnd;
}

/// <summary>
/// Chaining parameters of minimap2 (defaults of <c>mm_mapopt_init</c> 2.26 for DNA-vs-DNA mapping; <c>-k -g -r -z -n -m</c> and the chain scales).
/// </summary>
/// <remarks>
/// Validation follows the program (A4-4, measured on the installed 2.26): <c>mm_check_opt</c> refuses only a non-positive k ("-k and -w must be positive") - and <c>mm_sketch</c> asserts k &lt;= 28 -, so
/// <see cref="KmerSize"/> outside 1 - 28 is an <see cref="ArgumentOutOfRangeException"/>. Every other value is ACCEPTED by the program and by this port with the program's arithmetic: a negative or zero bandwidth
/// (<c>-r -5</c>: no anchor pair is within the bandwidth, no chain), a negative gap limit (<c>-g -1</c>: the bandwidth is the limit), <c>-n 0</c> / <c>-n -1</c>, <c>-m -5</c>, a zero or negative number of predecessors
/// (<c>--max-chain-iter 0</c>: every anchor starts its own chain), negative scales, and a NaN / infinite / huge scale, where the C float-to-int conversion of the penalty gives <c>INT_MIN</c> (x86-64
/// <c>cvttss2si</c>), reproduced by the port (the chain scores equal the installed program's, see <c>minimap2_option_boundary_oracle.json</c>).
/// </remarks>
public sealed record MinimapChainOptions
{
    /// <summary><c>is_cdna</c> of <c>mg_lchain_dp</c> (<c>--splice</c> / <c>-x splice</c>): a long gap on the target (an intron) costs <c>min(linear, log)</c> penalty, the query gap limit is not raised to the bandwidth and the chain is never cut by the score drop.</summary>
    public bool IsCdna { get; init; }

    /// <summary>k-mer size of the index (<c>-k</c>, 15): scales the gap penalties.</summary>
    public int KmerSize { get; init; } = 15;

    /// <summary>Largest distance on the target between two anchors of a chain (<c>-g</c>, 5000; at least the bandwidth).</summary>
    public int MaxTargetGap { get; init; } = 5000;

    /// <summary>Largest distance on the query (<c>-g</c>, 5000; at least the bandwidth).</summary>
    public int MaxQueryGap { get; init; } = 5000;

    /// <summary>Chaining bandwidth (<c>-r</c>, 500): largest difference between the target and query distance of two anchors.</summary>
    public int Bandwidth { get; init; } = 500;

    /// <summary>Stop the search for a predecessor after this many consecutive anchors that do not improve it (<c>max_chain_skip</c>, 25).</summary>
    public int MaxChainSkip { get; init; } = 25;

    /// <summary>Largest number of predecessors examined per anchor (<c>max_chain_iter</c>, 5000).</summary>
    public int MaxChainIterations { get; init; } = 5000;

    /// <summary>Fewest anchors on a chain (<c>-n</c>, 3).</summary>
    public int MinCount { get; init; } = 3;

    /// <summary>Smallest chain score (<c>-m</c>, 40).</summary>
    public int MinScore { get; init; } = 40;

    /// <summary>Gap penalty scale (<c>chain_gap_scale</c>, 0.8): the linear penalty per base of diagonal difference is <c>scale * 0.01 * k</c>.</summary>
    public float ChainGapScale { get; init; } = 0.8f;

    /// <summary>Skip penalty scale (<c>chain_skip_scale</c>, 0).</summary>
    public float ChainSkipScale { get; init; }
}

/// <summary>A chain of <see cref="Minimap2Chaining.MapChains"/>: its score, target, strand and anchors.</summary>
/// <param name="Score">Chain score.</param>
/// <param name="TargetId">Index of the target sequence.</param>
/// <param name="Reverse">The chain is on the reverse strand of the query.</param>
/// <param name="Anchors">The anchors of the chain in increasing order along the target.</param>
public sealed record MinimapChainHit(int Score, int TargetId, bool Reverse, IReadOnlyList<MinimapAnchor> Anchors);

/// <summary>A chain of anchors found by <see cref="Minimap2Chaining.Chain"/>.</summary>
/// <param name="Score">Chain score (minimap2's <c>s1</c> tag of the corresponding hit).</param>
/// <param name="Anchors">Indices (into the input anchor list) of the anchors of the chain, in increasing order along the chain.</param>
public sealed record MinimapChain(int Score, IReadOnlyList<int> Anchors);

/// <summary>
/// The chaining step of minimap2 (Li 2018, Bioinformatics 34:3094; Li 2021, Bioinformatics 37:4572): <c>mg_lchain_dp</c> of <c>lchain.c</c> (2.26), the dynamic program that gives every anchor the best score of a chain ending in it
/// — <c>min(min(dq, dr), span)</c> of the step minus a gap penalty <c>0.8 * 0.01 * k * |dr - dq| + 0.5 * log2(|dr - dq| + 1)</c> — with minimap2's heuristics (at most <c>max_iter</c> predecessors, a
/// stop after <c>max_skip</c> non-improving predecessors, the bandwidth and distance limits, the "max-f" shortcut), then the backtracking <c>mg_chain_backtrack</c> with the Z-drop-like peak rule
/// (<c>max_drop = bandwidth</c>) and the compaction <c>compact_a</c> that orders the chains by the target position of their first anchor. Port line by line; the float arithmetic of the program is kept.
/// </summary>
public static partial class Minimap2Chaining
{
    private const int MinInt = int.MinValue;

    // mm_check_opt: "-k and -w must be positive"; mm_sketch: assert(k > 0 && k <= 28). Everything else is accepted by the program (see MinimapChainOptions).
    private static void ValidateOptions(MinimapChainOptions o, string paramName)
    {
        if (o.KmerSize < 1 || o.KmerSize > 28)
            throw new ArgumentOutOfRangeException(paramName, o.KmerSize, "KmerSize must be 1 - 28: minimap2 stops with \"-k and -w must be positive\" (mm_check_opt) for k <= 0 and fails the assertion k <= 28 of mm_sketch above 28.");
    }

    // (int)f in the C program on x86-64 (cvttss2si): NaN and every value outside the int range give INT_MIN (0x80000000); .NET gives 0 / saturates
    private static int CInt(float f) => float.IsNaN(f) || f >= 2147483648f || f < -2147483648f ? int.MinValue : (int)f;

    /// <summary>
    /// Chains the anchors. The anchors must be sorted as minimap2 sorts them (by <c>rev | target id | target position</c>, the order of <c>--print-seeds</c>; <see cref="SortAnchors"/>).
    /// </summary>
    /// <param name="anchors">The anchors in minimap2's order.</param>
    /// <param name="options">Chaining options (null = minimap2 defaults).</param>
    /// <param name="cancellationToken">Cancellation token, polled every 64 anchors (A4-1).</param>
    /// <returns>The chains in the order minimap2 stores them (by the target position of the first anchor); empty when no chain reaches the minimum score and anchor count.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="anchors"/> is null.</exception>
    /// <exception cref="ArgumentException">The anchors are not sorted by (reverse, target id, target position).</exception>
    public static IReadOnlyList<MinimapChain> Chain(IReadOnlyList<MinimapAnchor> anchors, MinimapChainOptions? options = null, CancellationToken cancellationToken = default)
    {
        ArgumentNullException.ThrowIfNull(anchors);
        using var scope = AlignmentCancellation.Enter(cancellationToken);
        var o = options ?? new MinimapChainOptions();
        ValidateOptions(o, nameof(options));
        int n = anchors.Count;
        if (n == 0)
            return [];
        if (n > SequenceAligner.MinimapMaxAnchors)
            throw SequenceAligner.ResourceLimit("minimap2 chaining", "keeps four arrays (score, predecessor, peak score, mark: 36 bytes) per anchor besides the anchors", $"{n} anchors", SequenceAligner.MinimapMaxAnchors,
                nameof(SequenceAligner) + "." + nameof(SequenceAligner.MinimapMaxAnchors), "The program has no such limit and dies at the allocation", "chain the anchors of the query in pieces", nameof(anchors));
        var ax = new ulong[n];
        var ay = new ulong[n];
        for (int i = 0; i < n; i++)
        {
            ax[i] = anchors[i].X;
            ay[i] = anchors[i].Y;
            if (i > 0 && ax[i] < ax[i - 1])
                throw new ArgumentException("The anchors must be sorted by (reverse, target id, target position).", nameof(anchors));
        }

        return ChainCore(ax, ay, o);
    }

    /// <summary><c>mg_lchain_dp</c> + <c>mg_chain_backtrack</c> + <c>compact_a</c> over the anchor arrays <c>x</c> / <c>y</c> (sorted by <c>x</c>); the chain members are indices into the arrays.</summary>
    internal static IReadOnlyList<MinimapChain> ChainCore(ulong[] ax, ulong[] ay, MinimapChainOptions o)
    {
        int n = ax.Length;
        float penGap = (float)(o.ChainGapScale * 0.01 * o.KmerSize);
        float penSkip = (float)(o.ChainSkipScale * 0.01 * o.KmerSize);
        int bw = o.Bandwidth, maxDistX = Math.Max(o.MaxTargetGap, bw), maxDistY = o.IsCdna ? o.MaxQueryGap : Math.Max(o.MaxQueryGap, bw), maxDrop = o.IsCdna ? int.MaxValue : bw;

        var f = new int[n];
        var p = new long[n];
        var v = new int[n];
        var t = new int[n];
        long st = 0, maxIi = -1;
        for (long i = 0; i < n; i++)
        {
            if ((i & 63) == 0)
                AlignmentCancellation.ThrowIfRequested();
            long maxJ = -1;
            int maxF = (int)(ay[i] >> 32 & 0xff), nSkip = 0;
            while (st < i && (ax[i] >> 32 != ax[st] >> 32 || ax[i] > ax[st] + (ulong)maxDistX))
                st++;
            if (i - st > o.MaxChainIterations)
                st = i - o.MaxChainIterations;
            long j;
            for (j = i - 1; j >= st; j--)
            {
                int sc = Score(ax, ay, i, j, maxDistX, maxDistY, bw, penGap, penSkip, o.IsCdna);
                if (sc == MinInt)
                    continue;
                sc += f[j];
                if (sc > maxF)
                {
                    maxF = sc;
                    maxJ = j;
                    if (nSkip > 0)
                        nSkip--;
                }
                else if (t[j] == (int)i && ++nSkip > o.MaxChainSkip)
                {
                    break;
                }

                if (p[j] >= 0)
                    t[p[j]] = (int)i;
            }

            long endJ = j;
            if (maxIi < 0 || ax[i] - ax[maxIi] > (ulong)maxDistX)
            {
                int max = MinInt;
                maxIi = -1;
                for (j = i - 1; j >= st; j--)
                    if (max < f[j])
                    {
                        max = f[j];
                        maxIi = j;
                    }
            }

            if (maxIi >= 0 && maxIi < endJ)
            {
                int tmp = Score(ax, ay, i, maxIi, maxDistX, maxDistY, bw, penGap, penSkip, o.IsCdna);
                if (tmp != MinInt && maxF < tmp + f[maxIi])
                {
                    maxF = tmp + f[maxIi];
                    maxJ = maxIi;
                }
            }

            f[i] = maxF;
            p[i] = maxJ;
            v[i] = maxJ >= 0 && v[maxJ] > maxF ? v[maxJ] : maxF;
            if (maxIi < 0 || (ax[i] - ax[maxIi] <= (ulong)maxDistX && f[maxIi] < f[i]))
                maxIi = i;
        }

        return Backtrack(n, f, p, t, o.MinCount, o.MinScore, maxDrop, ax);
    }

    /// <summary>
    /// Sorts anchors as minimap2 does before chaining: <c>radix_sort_128x</c> of klib on <c>reverse | target id | target position</c>. That sort is an in-place MSD radix sort and is <b>not stable</b> above 64
    /// elements (anchors with the same key come out in the order the buckets leave them): it is ported exactly, so equal keys end up in minimap2's order.
    /// </summary>
    public static IReadOnlyList<MinimapAnchor> SortAnchors(IEnumerable<MinimapAnchor> anchors)
    {
        ArgumentNullException.ThrowIfNull(anchors);
        var list = anchors.ToList();
        var a = list.Select(x => new Mm128(x.X, x.Y)).ToArray();
        RadixSort(a);
        return a.Select(m => new MinimapAnchor((int)(m.X >> 32 & 0x7fffffff), (m.X >> 63) != 0, (int)(uint)m.X, (int)(uint)m.Y, (int)(m.Y >> 32 & 0xff))).ToList();
    }

    internal struct Mm128(ulong x, ulong y)
    {
        public ulong X = x;
        public ulong Y = y;
    }

    // klib krsort / ksort.h KRADIX_SORT_INIT(128x, mm128_t, key = x, 8): insertion sort up to 64 elements, otherwise the in-place MSD radix sort (8 bits per digit)
    private static void RadixSort(Mm128[] a) => RadixSort(a, 0, a.Length);

    private static void RadixSort(Mm128[] a, int beg, int end)
    {
        if (end - beg <= 64)
            InsertSort(a, beg, end);
        else
            RadixPass(a, beg, end, 8, 56);
    }

    private static void InsertSort(Mm128[] a, int beg, int end)
    {
        for (int i = beg + 1; i < end; i++)
        {
            if (a[i].X >= a[i - 1].X)
                continue;
            var tmp = a[i];
            int j = i;
            for (; j > beg && tmp.X < a[j - 1].X; j--)
                a[j] = a[j - 1];
            a[j] = tmp;
        }
    }

    private static void RadixPass(Mm128[] a, int beg, int end, int nBits, int shift)
    {
        int size = 1 << nBits, mask = size - 1;
        var bb = new int[size];
        var be = new int[size];
        Array.Fill(bb, beg);
        Array.Fill(be, beg);
        for (int i = beg; i != end; i++)
            be[(int)(a[i].X >> shift & (ulong)mask)]++;
        for (int k = 1; k != size; k++)
        {
            be[k] += be[k - 1] - beg;
            bb[k] = be[k - 1];
        }

        for (int k = 0; k != size;)
        {
            if (bb[k] != be[k])
            {
                int l = (int)(a[bb[k]].X >> shift & (ulong)mask);
                if (l != k)
                {
                    var tmp = a[bb[k]];
                    do
                    {
                        var swap = tmp;
                        tmp = a[bb[l]];
                        a[bb[l]++] = swap;
                        l = (int)(tmp.X >> shift & (ulong)mask);
                    }
                    while (l != k);
                    a[bb[k]++] = tmp;
                }
                else
                {
                    bb[k]++;
                }
            }
            else
            {
                k++;
            }
        }

        bb[0] = beg;
        for (int k = 1; k != size; k++)
            bb[k] = be[k - 1];
        if (shift != 0)
        {
            int next = shift > nBits ? shift - nBits : 0;
            for (int k = 0; k != size; k++)
            {
                if (be[k] - bb[k] > 64)
                    RadixPass(a, bb[k], be[k], nBits, next);
                else if (be[k] - bb[k] > 1)
                    InsertSort(a, bb[k], be[k]);
            }
        }
    }

    // comput_sc (one segment)
    private static int Score(ulong[] ax, ulong[] ay, long i, long j, int maxDistX, int maxDistY, int bw, float penGap, float penSkip, bool isCdna)
    {
        int dq = (int)ay[i] - (int)ay[j];
        if (dq <= 0 || dq > maxDistX)
            return MinInt;
        int dr = (int)(ax[i] - ax[j]);
        if (dr == 0 || dq > maxDistY)
            return MinInt;
        int dd = dr > dq ? dr - dq : dq - dr;
        if (dd > bw)
            return MinInt;
        int dg = dr < dq ? dr : dq;
        int qSpan = (int)(ay[j] >> 32 & 0xff);
        int sc = qSpan < dg ? qSpan : dg;
        if (dd != 0 || dg > qSpan)
        {
            float linPen = penGap * dd + penSkip * dg;
            float logPen = dd >= 1 ? Log2(dd + 1) : 0.0f;
            if (isCdna && dr > dq)
                sc = unchecked(sc - CInt(linPen < logPen ? linPen : logPen));         // a deletion (an intron on the target): the cheaper of the two penalties
            else
                sc = unchecked(sc - CInt(linPen + .5f * logPen));
        }

        return sc;
    }

    // mg_log2: a quadratic fit of log2 on the mantissa (only valid for x >= 2)
    private static float Log2(float x)
    {
        uint bits = BitConverter.SingleToUInt32Bits(x);
        float log2 = (int)(bits >> 23 & 255) - 128;
        bits &= ~(255u << 23);
        bits += 127u << 23;
        float m = BitConverter.UInt32BitsToSingle(bits);
        log2 += (-0.34484843f * m + 2.02466578f) * m - 0.67487759f;
        return log2;
    }

    // mg_chain_bk_end
    private static long BackEnd(int maxDrop, Mm128[] z, int[] f, long[] p, int[] t, int k)
    {
        long i = (long)z[k].Y, endI = -1, maxI = i;
        int maxS = 0;
        if (i < 0 || t[i] != 0)
            return i;
        do
        {
            t[i] = 2;
            endI = i = p[i];
            int s = i < 0 ? (int)z[k].X : (int)z[k].X - f[i];
            if (s > maxS)
            {
                maxS = s;
                maxI = i;
            }
            else if (maxS - s > maxDrop)
            {
                break;
            }
        }
        while (i >= 0 && t[i] == 0);

        for (i = (long)z[k].Y; i >= 0 && i != endI; i = p[i])
            t[i] = 0;
        return maxI;
    }

    // mg_chain_backtrack + compact_a
    private static List<MinimapChain> Backtrack(int n, int[] f, long[] p, int[] t, int minCount, int minScore, int maxDrop, ulong[] ax)
    {
        var zList = new List<Mm128>();
        for (int i = 0; i < n; i++)
            if (f[i] >= minScore)
                zList.Add(new Mm128((ulong)f[i], (ulong)i));
        if (zList.Count == 0)
            return [];
        var z = zList.ToArray();
        RadixSort(z);

        Array.Clear(t);
        var chains = new List<(int Score, List<int> Anchors)>();
        for (int k = z.Length - 1; k >= 0; k--)
        {
            if ((k & 255) == 0)
                AlignmentCancellation.ThrowIfRequested();
            if (t[(int)z[k].Y] != 0)
                continue;
            long endI = BackEnd(maxDrop, z, f, p, t, k);
            var members = new List<int>();
            long i;
            for (i = (long)z[k].Y; i != endI; i = p[i])
            {
                members.Add((int)i);
                t[i] = 1;
            }

            int sc = i < 0 ? (int)z[k].X : (int)z[k].X - f[i];
            if (sc >= minScore && members.Count > 0 && members.Count >= minCount)
                chains.Add((sc, members));
        }

        if (chains.Count == 0)
            return [];

        // compact_a: every chain is written in ascending order along the chain; the chains are sorted (radix_sort_128x on the target position of the first anchor) so that adjacent chains may be joined
        var asc = chains.Select(c => (c.Score, Anchors: Enumerable.Reverse(c.Anchors).ToList())).ToList();
        var w = new Mm128[asc.Count];
        for (int i = 0; i < w.Length; i++)
            w[i] = new Mm128(ax[asc[i].Anchors[0]], (ulong)i);
        RadixSort(w);
        return w.Select(m => new MinimapChain(asc[(int)m.Y].Score, asc[(int)m.Y].Anchors)).ToList();
    }

    // ------------------------------------------------------------------------------------------------ minimizer seeds (sketch.c, index.c, seed.c, map.c)

    /// <summary>
    /// The seed hits minimap2 collects for a query against a set of targets with a minimizer index of k-mer size <paramref name="k"/> and window <paramref name="w"/>
    /// (<c>mm_sketch</c> of <c>sketch.c</c> on every target and on the query, the hash index, <c>mm_seed_mz_flt</c>, <c>mm_collect_matches</c> with <c>mm_seed_select</c> and <c>collect_seed_hits</c>), sorted as minimap2 sorts them.
    /// </summary>
    /// <remarks>
    /// Frequent minimizers (A4-7): a minimizer that occurs more than <c>mid_occ</c> times (<c>max(10, the (1 - 2e-4) quantile of the occurrence counts + 1)</c>) is dropped except for the least frequent ones of a streak of
    /// frequent minimizers (<c>mm_seed_select</c>, one per 500 query bases, never above 4095 occurrences). The options of <see cref="CollectSeeds"/> (<c>-H</c>, <c>-T</c>, <c>-f</c>, <c>-U</c>, <c>-X</c>, ...) are available there. The query is searched on both strands.
    /// </remarks>
    /// <param name="targets">The target (reference) sequences, A/C/G/T/U (other symbols break the k-mer), any case.</param>
    /// <param name="query">The query sequence.</param>
    /// <param name="k">k-mer size (1-28, <c>-k</c>).</param>
    /// <param name="w">Window size (1-255, <c>-w</c>).</param>
    /// <param name="queryOccurrenceFraction">Query minimizers that occur more than <c>mid_occ</c> times and in more than this fraction of the query's minimizers are dropped (<c>q_occ_frac</c>, 0.01). As in the program every
    /// value is accepted: a zero, negative or NaN fraction switches the filter off (<c>q_occ_frac &gt; 0</c> is false), an infinite one drops nothing (measured: <c>--q-occ-frac nan / -1 / inf</c> give the chains of the default).</param>
    /// <param name="cancellationToken">Cancellation token, polled every 64 k-mer positions of a sequence and per query minimizer (A4-1).</param>
    /// <returns>The anchors, sorted by (reverse, target id, target position).</returns>
    /// <exception cref="ArgumentNullException">An argument is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="k"/> or <paramref name="w"/> is out of range.</exception>
    public static IReadOnlyList<MinimapAnchor> CollectAnchors(IReadOnlyList<string> targets, string query, int k = 15, int w = 10, float queryOccurrenceFraction = 0.01f, CancellationToken cancellationToken = default)
    {
        ArgumentNullException.ThrowIfNull(targets);
        ArgumentNullException.ThrowIfNull(query);
        return CollectSeeds(targets, query, new MinimapSeedOptions { KmerSize = k, WindowSize = w, QueryOccurrenceFraction = queryOccurrenceFraction }, cancellationToken: cancellationToken).Anchors;
    }

    /// <summary>
    /// The chains minimap2 builds for a query: <see cref="CollectAnchors"/> then <see cref="Chain"/> with the options of the index (<c>k</c> scales the gap penalties). Every chain is returned with the
    /// target, strand and anchors it consists of.
    /// </summary>
    /// <param name="targets">The target sequences.</param>
    /// <param name="query">The query sequence.</param>
    /// <param name="k">k-mer size.</param>
    /// <param name="w">Window size.</param>
    /// <param name="options">Chaining options; <see cref="MinimapChainOptions.KmerSize"/> is replaced by <paramref name="k"/> (null = defaults).</param>
    /// <param name="cancellationToken">Cancellation token (see <see cref="CollectAnchors"/> and <see cref="Chain"/>).</param>
    /// <returns>The chains, ordered by the position of their first anchor.</returns>
        public static IReadOnlyList<MinimapChainHit> MapChains(IReadOnlyList<string> targets, string query, int k = 15, int w = 10, MinimapChainOptions? options = null, CancellationToken cancellationToken = default)
    {
        using var scope = AlignmentCancellation.Enter(cancellationToken);
        var anchors = CollectAnchors(targets, query, k, w);
        var chains = Chain(anchors, (options ?? new MinimapChainOptions()) with { KmerSize = k });
        return chains.Select(c =>
        {
            var members = c.Anchors.Select(i => anchors[i]).ToList();
            return new MinimapChainHit(c.Score, members[0].TargetId, members[0].Reverse, members);
        }).ToList();
    }

    private static List<(ulong X, ulong Y)> Sketch(string str, int w, int k, uint rid, bool hpc = false)
    {
        var result = new List<(ulong X, ulong Y)>();
        int len = str.Length;
        if (len == 0)
            return result;
        ulong shift1 = (ulong)(2 * (k - 1)), mask = (1UL << (2 * k)) - 1;
        ulong[] kmer = [0, 0];
        int l = 0, bufPos = 0, minPos = 0, kmerSpanHpc = 0;
        var tq = new int[32]; // tiny_queue_t: the run lengths of the last k homopolymer runs
        int tqFront = 0, tqCount = 0;
        var bufX = new ulong[w];
        var bufY = new ulong[w];
        Array.Fill(bufX, ulong.MaxValue);
        Array.Fill(bufY, ulong.MaxValue);
        ulong minX = ulong.MaxValue, minY = ulong.MaxValue;

#pragma warning disable S127 // -H: the sketch moves to the end of a homopolymer run (sketch.c is_hpc)
        for (int i = 0; i < len; i++)
        {
            if ((i & 0xFFFF) == 0)
                AlignmentCancellation.ThrowIfRequested();
            int c = Nt4(str[i]);
            ulong infoX = ulong.MaxValue, infoY = ulong.MaxValue;
            if (c < 4)
            {
                int kmerSpan;
                if (hpc)
                {
                    // -H: homopolymer-compressed k-mers (sketch.c is_hpc): $i is put at the end of the current run, the span is the sum of the last k run lengths
                    int skipLen = 1;
                    if (i + 1 < len && Nt4(str[i + 1]) == c)
                    {
                        for (skipLen = 2; i + skipLen < len; ++skipLen)
                            if (Nt4(str[i + skipLen]) != c)
                                break;
                        i += skipLen - 1;
                    }

                    tq[(tqCount++ + tqFront) & 0x1f] = skipLen;
                    kmerSpanHpc += skipLen;
                    if (tqCount > k)
                    {
                        kmerSpanHpc -= tq[tqFront++];
                        tqFront &= 0x1f;
                        --tqCount;
                    }

                    kmerSpan = kmerSpanHpc;
                }
                else
                {
                    kmerSpan = l + 1 < k ? l + 1 : k;
                }

                kmer[0] = (kmer[0] << 2 | (uint)c) & mask;
                kmer[1] = (kmer[1] >> 2) | (3UL ^ (uint)c) << (int)shift1;
                if (kmer[0] == kmer[1])
                    continue; // symmetric k-mer: the strand is unknown
                int z = kmer[0] < kmer[1] ? 0 : 1;
                l++;
                if (l >= k && kmerSpan < 256)
                {
                    infoX = Hash64(kmer[z], mask) << 8 | (uint)kmerSpan;
                    infoY = (ulong)rid << 32 | (uint)i << 1 | (uint)z;
                }
            }
            else
            {
                l = 0;
                tqCount = tqFront = 0;
                kmerSpanHpc = 0;
            }

            bufX[bufPos] = infoX;
            bufY[bufPos] = infoY;
            if (l == w + k - 1 && minX != ulong.MaxValue)
            {
                for (int j = bufPos + 1; j < w; j++)
                    if (minX == bufX[j] && bufY[j] != minY)
                        result.Add((bufX[j], bufY[j]));
                for (int j = 0; j < bufPos; j++)
                    if (minX == bufX[j] && bufY[j] != minY)
                        result.Add((bufX[j], bufY[j]));
            }

            if (infoX <= minX)
            {
                if (l >= w + k && minX != ulong.MaxValue)
                    result.Add((minX, minY));
                minX = infoX;
                minY = infoY;
                minPos = bufPos;
            }
            else if (bufPos == minPos)
            {
                if (l >= w + k - 1 && minX != ulong.MaxValue)
                    result.Add((minX, minY));
                minX = ulong.MaxValue;
                for (int j = bufPos + 1; j < w; j++)
                    if (minX >= bufX[j])
                    {
                        minX = bufX[j];
                        minY = bufY[j];
                        minPos = j;
                    }

                for (int j = 0; j <= bufPos; j++)
                    if (minX >= bufX[j])
                    {
                        minX = bufX[j];
                        minY = bufY[j];
                        minPos = j;
                    }

                if (l >= w + k - 1 && minX != ulong.MaxValue)
                {
                    for (int j = bufPos + 1; j < w; j++)
                        if (minX == bufX[j] && minY != bufY[j])
                            result.Add((bufX[j], bufY[j]));
                    for (int j = 0; j <= bufPos; j++)
                        if (minX == bufX[j] && minY != bufY[j])
                            result.Add((bufX[j], bufY[j]));
                }
            }

            if (++bufPos == w)
                bufPos = 0;
        }

#pragma warning restore S127
        if (minX != ulong.MaxValue)
            result.Add((minX, minY));
        return result;
    }

    // seq_nt4_table: A C G T (U) in either case = 0 1 2 3, everything else 4
    private static int Nt4(char ch) => ch switch
    {
        'A' or 'a' => 0,
        'C' or 'c' => 1,
        'G' or 'g' => 2,
        'T' or 't' or 'U' or 'u' => 3,
        _ => 4,
    };

    private static ulong Hash64(ulong key, ulong mask)
    {
        key = (~key + (key << 21)) & mask;
        key ^= key >> 24;
        key = (key + (key << 3) + (key << 8)) & mask;
        key ^= key >> 14;
        key = (key + (key << 2) + (key << 4)) & mask;
        key ^= key >> 28;
        key = (key + (key << 31)) & mask;
        return key;
    }
}
