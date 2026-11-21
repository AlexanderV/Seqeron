namespace Seqeron.Genomics.Alignment;

/// <summary>
/// Seeding options of minimap2 2.26 (A4-7): the minimizer sketch (<c>-k -w -H</c>, the <c>sdust</c> masking <c>-T</c>), the frequent-minimizer rules (<c>-f -U</c>, <c>--q-occ-frac</c>, <c>mid_occ</c> / <c>max_occ</c> /
/// <c>max_max_occ</c> / <c>occ_dist</c> of <c>mm_mapopt_t</c>) and the all-vs-all / strand switches (<c>-X</c>, <c>--for-only</c>, <c>--rev-only</c>). The defaults are those of <c>mm_mapopt_init</c>.
/// </summary>
public sealed record MinimapSeedOptions
{
    /// <summary>k-mer size of the index (<c>-k</c>, 1-28, default 15).</summary>
    public int KmerSize { get; init; } = 15;

    /// <summary>Minimizer window (<c>-w</c>, 1-255, default 10).</summary>
    public int WindowSize { get; init; } = 10;

    /// <summary><c>-H</c>: homopolymer-compressed k-mers (the index and the query; positions and spans stay in the original coordinates).</summary>
    public bool HomopolymerCompressed { get; init; }

    /// <summary><c>-T</c> (hidden option of 2.26, <c>sdust_thres</c>, default 0 = off): query minimizers of which more than half falls into a low-complexity region found by <c>sdust</c> (window 64, this threshold) are dropped.</summary>
    public int SdustThreshold { get; init; }

    /// <summary><c>-f</c> (<c>mid_occ_frac</c>, default 2e-4): the occurrence count above which a minimizer is frequent is the <c>(1 - f)</c> quantile of the occurrence counts of the index plus one (<c>mm_idx_cal_max_occ</c>); <c>f &lt;= 0</c> = no limit.</summary>
    public float MidOccFraction { get; init; } = 2e-4f;

    /// <summary><c>-U</c> lower bound (<c>min_mid_occ</c>, default 10) of <c>mid_occ</c>.</summary>
    public int MinMidOcc { get; init; } = 10;

    /// <summary><c>-U</c> upper bound (<c>max_mid_occ</c>, default 1 000 000) of <c>mid_occ</c> (only when larger than <see cref="MinMidOcc"/>).</summary>
    public int MaxMidOcc { get; init; } = 1_000_000;

    /// <summary>An explicit <c>mid_occ</c> (&gt; 0): replaces the value computed from <see cref="MidOccFraction"/> (<c>mm_mapopt_update</c> computes it only when it is not positive).</summary>
    public int MidOcc { get; init; }

    /// <summary><c>max_occ</c> (default 0 = unused; a preset value above <c>mid_occ</c> re-chains with this limit when the first chaining found no chain on all segments of a short-read pair).</summary>
    public int MaxOcc { get; init; }

    /// <summary><c>max_max_occ</c> (default 4095): minimizers with more occurrences are never used, even when <c>mm_seed_select</c> would keep one of a high-occurrence streak.</summary>
    public int MaxMaxOcc { get; init; } = 4095;

    /// <summary><c>occ_dist</c> (default 500; 0 = off, as in the <c>ava-*</c> presets): in a streak of frequent minimizers up to one per <c>occ_dist</c> query bases (the least frequent ones) is kept (<c>mm_seed_select</c>).</summary>
    public int OccDist { get; init; } = 500;

    /// <summary><c>--q-occ-frac</c> (default 0.01): query minimizers that occur more than <c>mid_occ</c> times in the query and in more than this fraction of its minimizers are dropped (<c>mm_seed_mz_flt</c>).</summary>
    public float QueryOccurrenceFraction { get; init; } = 0.01f;

    /// <summary><c>-X</c> (<c>MM_F_NO_DIAG</c>): hits of a target of the same name and length on the diagonal are skipped (and the same-strand hits of such a target are marked self).</summary>
    public bool SkipSelf { get; init; }

    /// <summary><c>-X</c> (<c>MM_F_NO_DUAL</c>): a target whose name sorts BEFORE the query name is skipped, so that a pair of sequences is mapped once (needs the names).</summary>
    public bool SkipDual { get; init; }

    /// <summary><c>--for-only</c>: only hits on the forward strand of the target.</summary>
    public bool ForwardOnly { get; init; }

    /// <summary><c>--rev-only</c>: only hits on the reverse strand.</summary>
    public bool ReverseOnly { get; init; }

    /// <summary><c>--heap-sort</c> (<c>MM_F_HEAP_SORT</c>, set by the <c>sr</c> preset): the anchors are collected by a k-way merge of the occurrence lists (<c>collect_seed_hits_heap</c>) instead of a radix sort; anchors with the same reference position
    /// come out in the heap order, not in the radix-sort order.</summary>
    public bool HeapSort { get; init; }

    /// <summary><c>--qstrand</c> (<c>MM_F_QSTRAND</c>): the query-strand mode: anchors of the reverse strand keep the query coordinate and take the reverse-complemented reference coordinate.</summary>
    public bool QueryStrand { get; init; }
}

/// <summary>The seed stage of minimap2 for one query: the sorted anchor array and the values the later stages use.</summary>
/// <param name="Anchors">The anchors (<c>--print-seeds</c>), sorted as minimap2 sorts them.</param>
/// <param name="RepeatLength">Length of the query regions covered by frequent minimizers (<c>rep_len</c>, the PAF <c>rl</c> tag).</param>
/// <param name="MidOcc">The occurrence limit <c>mid_occ</c> that was applied (the line <c>mid_occ = N</c> of the program's log).</param>
public sealed record MinimapSeeds(IReadOnlyList<MinimapAnchor> Anchors, int RepeatLength, int MidOcc);

public static partial class Minimap2Chaining
{
    // MM_SEED_TANDEM / MM_SEED_SELF: flag bits of mm128_t.y of an anchor (they enter the hash of a chain, mm_gen_regs)
    private const ulong SeedTandem = 1UL << 42;
    private const ulong SeedSelf = 1UL << 43;

    /// <summary>The query minimizers after sdust masking and <c>mm_seed_mz_flt</c>, the index and <c>mid_occ</c>; reused by the re-chaining with another occurrence limit.</summary>
    internal sealed class SeedContext
    {
        internal required Dictionary<ulong, List<ulong>> Index { get; init; }
        internal required List<(ulong X, ulong Y)> Minimizers { get; init; }
        internal required int MidOcc { get; init; }
        internal required int QueryLength { get; init; }
    }

    /// <summary>The anchor array <c>a</c> of <c>collect_seed_hits</c> with its repeat length and the minimizer list of <c>mm_est_err</c>.</summary>
    internal sealed class SeedArray
    {
        internal required Mm128[] A { get; init; }
        internal required int RepeatLength { get; init; }
        internal required ulong[] MiniPos { get; init; }
    }

    private static void ValidateSeedOptions(MinimapSeedOptions o)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(o.KmerSize, 1);
        ArgumentOutOfRangeException.ThrowIfGreaterThan(o.KmerSize, 28);
        ArgumentOutOfRangeException.ThrowIfLessThan(o.WindowSize, 1);
        ArgumentOutOfRangeException.ThrowIfGreaterThan(o.WindowSize, 255);
    }

    /// <summary>The index of <c>mm_idx_gen</c>: minimizer hash -> the positions <c>rid &lt;&lt; 32 | pos &lt;&lt; 1 | strand</c> in ascending order (<c>worker_post</c>).</summary>
    private static Dictionary<ulong, List<ulong>> BuildIndex(IReadOnlyList<string> targets, int k, int w, bool hpc)
    {
        var index = new Dictionary<ulong, List<ulong>>();
        for (int rid = 0; rid < targets.Count; rid++)
        {
            ArgumentNullException.ThrowIfNull(targets[rid]);
            int indexed = 0;
            foreach (var (x, y) in Sketch(targets[rid], w, k, (uint)rid, hpc))
            {
                if ((++indexed & 0xFFFF) == 0)
                    AlignmentCancellation.ThrowIfRequested();
                if (!index.TryGetValue(x >> 8, out var list))
                    index[x >> 8] = list = [];
                list.Add(y);
            }
        }

        foreach (var list in index.Values)
        {
            AlignmentCancellation.ThrowIfRequested();
            list.Sort(); // worker_post: radix_sort_64 of the positions of a minimizer (ascending: rid, position, strand)
        }

        return index;
    }

    /// <summary><c>mm_idx_cal_max_occ</c>: the <c>(1 - f)</c> quantile of the occurrence counts (a minimizer seen once counts 1) plus one; <c>int.MaxValue</c> for <c>f &lt;= 0</c>.</summary>
    private static int CalcMaxOcc(Dictionary<ulong, List<ulong>> index, float f)
    {
        if (f <= 0f)
            return int.MaxValue;
        if (index.Count == 0)
            return 1; // the program sorts an empty array and reads element 0: undefined; 1 is below every min_mid_occ
        var counts = index.Values.Select(l => l.Count).OrderBy(c => c).ToArray();
        return counts[(int)((1.0 - f) * counts.Length)] + 1;
    }

    private static int ResolveMidOcc(Dictionary<ulong, List<ulong>> index, MinimapSeedOptions o)
    {
        if (o.MidOcc > 0)
            return o.MidOcc;
        int midOcc = CalcMaxOcc(index, o.MidOccFraction);
        if (midOcc < o.MinMidOcc)
            midOcc = o.MinMidOcc;
        if (o.MaxMidOcc > o.MinMidOcc && midOcc > o.MaxMidOcc)
            midOcc = o.MaxMidOcc;
        return midOcc;
    }

    /// <summary>The minimizer index of the targets and <c>mid_occ</c> (<c>mm_idx_gen</c> + <c>mm_mapopt_update</c>), shared by all queries of a run.</summary>
    internal sealed class SeedIndex
    {
        internal required Dictionary<ulong, List<ulong>> Index { get; init; }
        internal required int MidOcc { get; init; }
    }

    internal static SeedIndex BuildSeedIndex(IReadOnlyList<string> targets, MinimapSeedOptions o)
    {
        ArgumentNullException.ThrowIfNull(targets);
        ValidateSeedOptions(o);
        var index = BuildIndex(targets, o.KmerSize, o.WindowSize, o.HomopolymerCompressed);
        return new SeedIndex { Index = index, MidOcc = ResolveMidOcc(index, o) };
    }

    /// <summary>The query minimizers of <c>collect_minimizers</c> (sketch of every segment, the positions of the later segments shifted by the length of the earlier ones, sdust) and <c>mm_seed_mz_flt</c>.</summary>
    internal static SeedContext PrepareQuerySeeds(SeedIndex ix, IReadOnlyList<string> segments, MinimapSeedOptions o)
    {
        var mv = new List<(ulong X, ulong Y)>();
        int sum = 0;
        for (int s = 0; s < segments.Count; s++)
        {
            var part = Sketch(segments[s], o.WindowSize, o.KmerSize, (uint)s, o.HomopolymerCompressed);
            for (int j = 0; j < part.Count; j++)
                part[j] = (part[j].X, part[j].Y + ((ulong)sum << 1));
            if (o.SdustThreshold > 0)
                part = DustMinimizers(part, segments[s], o.SdustThreshold);
            mv.AddRange(part);
            sum += segments[s].Length;
        }

        int midOcc = ix.MidOcc;
        if (mv.Count > midOcc && o.QueryOccurrenceFraction > 0f && midOcc > 0)
        {
            // mm_seed_mz_flt
            var drop = new HashSet<ulong>();
            foreach (var grp in mv.GroupBy(m => m.X))
            {
                int cnt = grp.Count();
                if (cnt > midOcc && cnt > mv.Count * o.QueryOccurrenceFraction)
                    drop.Add(grp.Key);
            }

            mv = mv.Where(m => !drop.Contains(m.X)).ToList();
        }

        return new SeedContext { Index = ix.Index, Minimizers = mv, MidOcc = midOcc, QueryLength = sum };
    }

    /// <summary>The index, <c>mid_occ</c> and the query minimizers of <c>collect_minimizers</c> (sketch, sdust) and <c>mm_seed_mz_flt</c>.</summary>
    internal static SeedContext PrepareSeeds(IReadOnlyList<string> targets, string query, MinimapSeedOptions o)
    {
        ArgumentNullException.ThrowIfNull(targets);
        ArgumentNullException.ThrowIfNull(query);
        return PrepareQuerySeeds(BuildSeedIndex(targets, o), [query], o);
    }

    // ---------------------------------------------------------------------------------------------------------- sdust (sdust.c)

    private const int SdWlen = 3;
    private const int SdWtot = 1 << (SdWlen << 1);
    private const int SdWmsk = SdWtot - 1;

    private struct PerfectInterval
    {
        public int Start, Finish, R, L;
    }

    /// <summary><c>mm_dust_minier</c>: squeeze out the minimizers of which more than half overlaps a low-complexity region (<c>sdust_core</c>, window 64).</summary>
    private static List<(ulong X, ulong Y)> DustMinimizers(List<(ulong X, ulong Y)> a, string seq, int threshold)
    {
        var dreg = SdustCore(seq, threshold, 64);
        var kept = new List<(ulong X, ulong Y)>(a.Count);
        int u = 0, nDreg = dreg.Count;
        foreach (var m in a)
        {
            int qpos = (int)((uint)m.Y >> 1), span = (int)(m.X & 0xff);
            int s = qpos - (span - 1), e = s + span;
            while (u < nDreg && (int)dreg[u] <= s)
                ++u;
            if (u < nDreg && (int)(dreg[u] >> 32) < e)
            {
                int l = 0;
                for (int v = u; v < nDreg && (int)(dreg[v] >> 32) < e; ++v)
                {
                    int ss = s > (int)(dreg[v] >> 32) ? s : (int)(dreg[v] >> 32);
                    int ee = e < (int)dreg[v] ? e : (int)(uint)dreg[v];
                    l += ee - ss;
                }

                if (l <= span >> 1)
                    kept.Add(m);
            }
            else
            {
                kept.Add(m);
            }
        }

        return kept;
    }

    /// <summary><c>sdust_core</c> (Morgulis et al. 2006, the symmetric DUST of minimap2's <c>sdust.c</c>): the masked intervals <c>start &lt;&lt; 32 | finish</c> (0-based, half-open) of a sequence.</summary>
    internal static List<ulong> SdustCore(string seq, int threshold, int window)
    {
        var res = new List<ulong>();
        var perfect = new List<PerfectInterval>();
        var w = new List<int>(); // kdq: the words of the current window
        var cv = new int[SdWtot];
        var cw = new int[SdWtot];
        int rv = 0, rw = 0, bigL = 0, l = 0, start;
        uint t = 0;
        int lSeq = seq.Length;
        for (int i = 0; i <= lSeq; ++i)
        {
            int b = i < lSeq ? Nt4(seq[i]) : 4;
            if (b < 4)
            {
                ++l;
                t = (t << 2 | (uint)b) & SdWmsk;
                if (l >= SdWlen)
                {
                    start = (l - window > 0 ? l - window : 0) + (i + 1 - l);
                    SaveMaskedRegions(res, perfect, start);
                    ShiftWindow((int)t, w, threshold, window, ref bigL, ref rw, ref rv, cw, cv);
                    if (rw * 10 > bigL * threshold)
                        FindPerfect(perfect, w, threshold, start, bigL, rv, cv);
                }
            }
            else
            {
                start = (l - window + 1 > 0 ? l - window + 1 : 0) + (i + 1 - l);
                while (perfect.Count > 0)
                    SaveMaskedRegions(res, perfect, start++);
                l = 0;
                t = 0;
            }
        }

        return res;
    }

    private static void ShiftWindow(int t, List<int> w, int tThreshold, int window, ref int bigL, ref int rw, ref int rv, int[] cw, int[] cv)
    {
        int s;
        if (w.Count >= window - SdWlen + 1)
        {
            s = w[0];
            w.RemoveAt(0);
            rw -= --cw[s];
            if (bigL > w.Count)
            {
                --bigL;
                rv -= --cv[s];
            }
        }

        w.Add(t);
        ++bigL;
        rw += cw[t]++;
        rv += cv[t]++;
        if (cv[t] * 10 > tThreshold << 1)
        {
            do
            {
                s = w[w.Count - bigL];
                rv -= --cv[s];
                --bigL;
            }
            while (s != t);
        }
    }

    private static void SaveMaskedRegions(List<ulong> res, List<PerfectInterval> p, int start)
    {
        if (p.Count == 0 || p[^1].Start >= start)
            return;
        var last = p[^1];
        bool saved = false;
        if (res.Count > 0)
        {
            int s = (int)(res[^1] >> 32), f = (int)(uint)res[^1];
            if (last.Start <= f)
            {
                saved = true;
                res[^1] = (ulong)(uint)s << 32 | (uint)(f > last.Finish ? f : last.Finish);
            }
        }

        if (!saved)
            res.Add((ulong)(uint)last.Start << 32 | (uint)last.Finish);
        int keep = p.Count;
        while (keep > 0 && p[keep - 1].Start < start) // remove the perfect intervals that have fallen out of the window
            --keep;
        p.RemoveRange(keep, p.Count - keep);
    }

    private static void FindPerfect(List<PerfectInterval> pList, List<int> w, int tThreshold, int start, int bigL, int rv, int[] cv)
    {
        var c = (int[])cv.Clone();
        int r = rv, maxR = 0, maxL = 0;
        for (int i = w.Count - bigL - 1; i >= 0; --i)
        {
            int t = w[i];
            r += c[t]++;
            int newR = r, newL = w.Count - i - 1;
            if (newR * 10 > tThreshold * newL)
            {
                int j;
                for (j = 0; j < pList.Count && pList[j].Start >= i + start; ++j)
                {
                    var p = pList[j];
                    if (maxR == 0 || p.R * maxL > maxR * p.L)
                    {
                        maxR = p.R;
                        maxL = p.L;
                    }
                }

                if (maxR == 0 || newR * maxL >= maxR * newL)
                {
                    maxR = newR;
                    maxL = newL;
                    pList.Insert(j, new PerfectInterval { Start = i + start, Finish = w.Count + (SdWlen - 1) + start, R = newR, L = newL });
                }
            }
        }
    }

    // ---------------------------------------------------------------------------------------------------------- mm_seed_select / mm_collect_matches / collect_seed_hits (seed.c, map.c)

    private sealed class Seed
    {
        public uint QPos;
        public int QSpan;
        public List<ulong> Cr = [];
        public int N;
        public bool Tandem, Flt;
        public int SegId;
    }

    private const int MaxMaxHighOcc = 128;

    // ks_heapadjust of klib ksort.h for uint64_t: a MAX-heap (a parent is never smaller than a child)
    private static void HeapAdjust(int i, int n, ulong[] l)
    {
        int k = i;
        ulong tmp = l[i];
        while ((k = (k << 1) + 1) < n)
        {
            if (k != n - 1 && l[k] < l[k + 1])
                ++k;
            if (l[k] < tmp)
                break;
            l[i] = l[k];
            i = k;
        }

        l[i] = tmp;
    }

    private static void HeapMake(int n, ulong[] l)
    {
        for (int i = (n >> 1) - 1; i >= 0; --i)
            HeapAdjust(i, n, l);
    }

    /// <summary><c>mm_seed_select</c>: for the high-occurrence minimizers choose up to <c>(pe - ps) / dist</c> of each streak (the ones with the fewest occurrences).</summary>
    private static void SeedSelect(List<Seed> a, int len, int maxOcc, int maxMaxOcc, int dist)
    {
        int n = a.Count;
        var b = new ulong[MaxMaxHighOcc];
        if (n is 0 or 1)
            return;
        int m = 0;
        for (int i = 0; i < n; ++i)
            if (a[i].N > maxOcc)
                ++m;
        if (m == 0)
            return;
        int last0 = -1;
        for (int i = 0; i <= n; ++i)
        {
            if (i == n || a[i].N <= maxOcc)
            {
                if (i - last0 > 1)
                {
                    int ps = last0 < 0 ? 0 : (int)(a[last0].QPos >> 1);
                    int pe = i == n ? len : (int)(a[i].QPos >> 1);
                    int st = last0 + 1, en = i;
                    int maxHighOcc = (int)((double)(pe - ps) / dist + .499);
                    if (maxHighOcc > 0)
                    {
                        if (maxHighOcc > MaxMaxHighOcc)
                            maxHighOcc = MaxMaxHighOcc;
                        int j, k;
                        for (j = st, k = 0; j < en && k < maxHighOcc; ++j, ++k)
                            b[k] = (ulong)a[j].N << 32 | (uint)j;
                        HeapMake(k, b);
                        for (; j < en; ++j)
                        {
                            if (a[j].N < (int)(b[0] >> 32))
                            {
                                b[0] = (ulong)a[j].N << 32 | (uint)j;
                                HeapAdjust(0, k, b);
                            }
                        }

                        for (j = 0; j < k; ++j)
                            a[(int)(uint)b[j]].Flt = true;
                    }

                    for (int j = st; j < en; ++j)
                        a[j].Flt = !a[j].Flt;
                    for (int j = st; j < en; ++j)
                        if (a[j].N > maxMaxOcc)
                            a[j].Flt = true;
                }

                last0 = i;
            }
        }
    }

    /// <summary><c>mm_collect_matches</c> + <c>collect_seed_hits</c> with the occurrence limit <paramref name="maxOcc"/>: the sorted anchor array (with the tandem / self flags), the repeat length and the minimizer positions.</summary>
    internal static SeedArray CollectSeedArray(SeedContext ctx, MinimapSeedOptions o, int maxOcc, IReadOnlyList<string>? targetNames, IReadOnlyList<int>? targetLengths, string? queryName)
    {
        var mv = ctx.Minimizers;
        int qlen = ctx.QueryLength;
        // mm_seed_collect_all
        var m = new List<Seed>();
        for (int i = 0; i < mv.Count; i++)
        {
            var p = mv[i];
            if (!ctx.Index.TryGetValue(p.X >> 8, out var cr))
                continue;
            var q = new Seed { QPos = (uint)p.Y, QSpan = (int)(p.X & 0xff), Cr = cr, N = cr.Count, SegId = (int)(p.Y >> 32) };
            if (i > 0 && p.X >> 8 == mv[i - 1].X >> 8)
                q.Tandem = true;
            if (i < mv.Count - 1 && p.X >> 8 == mv[i + 1].X >> 8)
                q.Tandem = true;
            m.Add(q);
        }

        if (o.OccDist > 0 && o.MaxMaxOcc > maxOcc)
        {
            SeedSelect(m, qlen, maxOcc, o.MaxMaxOcc, o.OccDist);
        }
        else
        {
            foreach (var q in m)
                if (q.N > maxOcc)
                    q.Flt = true;
        }

        int repStart = 0, repEnd = 0, repLen = 0;
        long nA = 0;
        var kept = new List<Seed>();
        var miniPos = new List<ulong>();
        foreach (var q in m)
        {
            if (q.Flt)
            {
                int en = (int)(q.QPos >> 1) + 1, st = en - q.QSpan;
                if (st > repEnd)
                {
                    repLen += repEnd - repStart;
                    repStart = st;
                    repEnd = en;
                }
                else
                {
                    repEnd = en;
                }
            }
            else
            {
                nA += q.N;
                miniPos.Add((ulong)q.QSpan << 32 | q.QPos >> 1);
                kept.Add(q);
            }
        }

        repLen += repEnd - repStart;
        if (nA > SequenceAligner.MinimapMaxAnchors)
            throw SequenceAligner.ResourceLimit("minimap2 seeding", "stores one anchor (20 bytes) for every occurrence of every query minimizer and four arrays of them for the chaining",
                string.Create(System.Globalization.CultureInfo.InvariantCulture, $"{nA} anchors (the occurrences of {kept.Count} query minimizers in the targets)"), SequenceAligner.MinimapMaxAnchors, nameof(SequenceAligner) + "." + nameof(SequenceAligner.MinimapMaxAnchors),
                "The program has no such limit: it dies at the allocation, and drops the minimizers that occur more than mid_occ times (mm_seed_select)", "use a longer k-mer or window, mask the repeat, or map a shorter query", "query");

        bool filterNames = (o.SkipSelf || o.SkipDual) && queryName is not null && targetNames is not null;
        // skip_seed: false = keep; isSelf is the MM_SEED_SELF flag
        bool SkipSeed(ulong r, Seed q, out bool isSelf)
        {
            isSelf = false;
            if (filterNames)
            {
                int tid = (int)(r >> 32);
                int cmp = string.CompareOrdinal(queryName, targetNames![tid]);
                if (o.SkipSelf && cmp == 0 && targetLengths![tid] == qlen)
                {
                    if ((uint)r >> 1 == q.QPos >> 1)
                        return true;
                    if ((r & 1) == (q.QPos & 1))
                        isSelf = true;
                }

                if (o.SkipDual && cmp > 0)
                    return true;
            }

            if (o.ForwardOnly || o.ReverseOnly)
            {
                bool forward = (r & 1) == (q.QPos & 1);
                if (forward && o.ReverseOnly)
                    return true;
                if (!forward && o.ForwardOnly)
                    return true;
            }

            return false;
        }

        if (o.HeapSort)
        {
            // collect_seed_hits_heap: a k-way merge of the (sorted) occurrence lists; forward anchors are stored from the front, reverse anchors from the back (then reversed)
            long total = nA;
            var heapA = new Mm128[total];
            int nFor = 0, nRev = 0;
            var heap = new Mm128[kept.Count];
            int heapSize = 0;
            for (int i = 0; i < kept.Count; ++i)
            {
                if (kept[i].N > 0)
                    heap[heapSize++] = new Mm128(kept[i].Cr[0], (ulong)i << 32);
            }

            HeapMakeMm(heapSize, heap);
            int polledH = 0;
            while (heapSize > 0)
            {
                if ((++polledH & 1023) == 0)
                    AlignmentCancellation.ThrowIfRequested();
                var q = kept[(int)(heap[0].Y >> 32)];
                ulong r = heap[0].X;
                int rpos = (int)((uint)r >> 1);
                if (!SkipSeed(r, q, out bool isSelf))
                {
                    Mm128 pa;
                    if ((r & 1) == (q.QPos & 1))
                    {
                        pa = new Mm128((r & 0xffffffff00000000UL) | (uint)rpos, (ulong)q.QSpan << 32 | q.QPos >> 1);
                        pa.Y |= (ulong)q.SegId << 48;
                        if (q.Tandem)
                            pa.Y |= SeedTandem;
                        if (isSelf)
                            pa.Y |= SeedSelf;
                        heapA[nFor++] = pa;
                    }
                    else
                    {
                        pa = new Mm128(1UL << 63 | (r & 0xffffffff00000000UL) | (uint)rpos, (ulong)q.QSpan << 32 | (uint)(qlen - ((int)(q.QPos >> 1) + 1 - q.QSpan) - 1));
                        pa.Y |= (ulong)q.SegId << 48;
                        if (q.Tandem)
                            pa.Y |= SeedTandem;
                        if (isSelf)
                            pa.Y |= SeedSelf;
                        heapA[total - (++nRev)] = pa;
                    }
                }

                if ((uint)heap[0].Y < (uint)(q.N - 1))
                {
                    heap[0] = new Mm128(kept[(int)(heap[0].Y >> 32)].Cr[(int)(uint)heap[0].Y + 1], heap[0].Y + 1);
                }
                else
                {
                    heap[0] = heap[heapSize - 1];
                    --heapSize;
                }

                HeapDownMm(0, heapSize, heap);
            }

            for (long jj = 0; jj < nRev >> 1; ++jj)
            {
                var t = heapA[total - 1 - jj];
                heapA[total - 1 - jj] = heapA[total - (nRev - jj)];
                heapA[total - (nRev - jj)] = t;
            }

            var res = new Mm128[nFor + nRev];
            Array.Copy(heapA, 0, res, 0, nFor);
            Array.Copy(heapA, total - nRev, res, nFor, nRev);
            return new SeedArray { A = res, RepeatLength = repLen, MiniPos = miniPos.ToArray() };
        }

        var a = new List<Mm128>((int)Math.Min(nA, int.MaxValue));
        int polled = 0;
        foreach (var q in kept)
        {
            if ((++polled & 255) == 0)
                AlignmentCancellation.ThrowIfRequested();
            foreach (ulong r in q.Cr)
            {
                int rpos = (int)((uint)r >> 1);
                if (SkipSeed(r, q, out bool isSelf))
                    continue;
                ulong x, y;
                if ((r & 1) == (q.QPos & 1))
                {
                    x = (r & 0xffffffff00000000UL) | (uint)rpos;
                    y = (ulong)q.QSpan << 32 | q.QPos >> 1;
                }
                else if (!o.QueryStrand)
                {
                    x = 1UL << 63 | (r & 0xffffffff00000000UL) | (uint)rpos;
                    y = (ulong)q.QSpan << 32 | (uint)(qlen - ((int)(q.QPos >> 1) + 1 - q.QSpan) - 1);
                }
                else
                {
                    int len = targetLengths![(int)(r >> 32)];
                    x = 1UL << 63 | (r & 0xffffffff00000000UL) | (uint)(len - (rpos + 1 - q.QSpan) - 1);
                    y = (ulong)q.QSpan << 32 | q.QPos >> 1;
                }

                y |= (ulong)q.SegId << 48;
                if (q.Tandem)
                    y |= SeedTandem;
                if (isSelf)
                    y |= SeedSelf;
                a.Add(new Mm128(x, y));
            }
        }

        var arr = a.ToArray();
        RadixSort(arr);
        return new SeedArray { A = arr, RepeatLength = repLen, MiniPos = miniPos.ToArray() };
    }

    // ks_heapdown / ks_heapmake of ksort.h for KSORT_INIT(heap, mm128_t, heap_lt) with heap_lt(a, b) = a.x > b.x (a min-heap on x)
    private static void HeapDownMm(int i, int n, Mm128[] l)
    {
        int k = i;
        var tmp = l[i];
        while ((k = (k << 1) + 1) < n)
        {
            if (k != n - 1 && l[k].X > l[k + 1].X)
                ++k;
            if (l[k].X > tmp.X)
                break;
            l[i] = l[k];
            i = k;
        }

        l[i] = tmp;
    }

    private static void HeapMakeMm(int n, Mm128[] l)
    {
        for (int i = (n >> 1) - 1; i >= 0; --i)
            HeapDownMm(i, n, l);
    }

    private static IReadOnlyList<MinimapAnchor> ToAnchors(Mm128[] a)
    {
        var list = new List<MinimapAnchor>(a.Length);
        foreach (var m in a)
            list.Add(new MinimapAnchor((int)(m.X >> 32 & 0x7fffffff), (m.X >> 63) != 0, (int)(uint)m.X, (int)(uint)m.Y, (int)(m.Y >> 32 & 0xff)));
        return list;
    }

    /// <summary>
    /// The seed stage of minimap2 for a query with the options of <see cref="MinimapSeedOptions"/>: the sorted anchors (<c>--print-seeds</c>), the repeat length (<c>rl</c>) and <c>mid_occ</c>. Frequent minimizers
    /// (more than <c>mid_occ</c> occurrences) are dropped except for the least frequent ones of a streak (<c>mm_seed_select</c>); the dropped query bases are the repeat length.
    /// </summary>
    /// <param name="targets">The target (reference) sequences, A/C/G/T/U (other symbols break the k-mer), any case.</param>
    /// <param name="query">The query sequence.</param>
    /// <param name="options">Seeding options (null = the defaults of <c>minimap2 2.26</c>).</param>
    /// <param name="targetNames">Names of the targets (needed by <see cref="MinimapSeedOptions.SkipSelf"/> / <see cref="MinimapSeedOptions.SkipDual"/>), null = none.</param>
    /// <param name="queryName">Name of the query (the first word of the FASTA header), null = none.</param>
    /// <param name="cancellationToken">Cancellation token.</param>
    /// <exception cref="ArgumentNullException">An argument is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException">k or w is out of range.</exception>
    public static MinimapSeeds CollectSeeds(IReadOnlyList<string> targets, string query, MinimapSeedOptions? options = null, IReadOnlyList<string>? targetNames = null, string? queryName = null, CancellationToken cancellationToken = default)
    {
        using var scope = AlignmentCancellation.Enter(cancellationToken);
        var o = options ?? new MinimapSeedOptions();
        var ctx = PrepareSeeds(targets, query, o);
        var arr = CollectSeedArray(ctx, o, ctx.MidOcc, targetNames, targets.Select(t => t.Length).ToArray(), queryName);
        return new MinimapSeeds(ToAnchors(arr.A), arr.RepeatLength, ctx.MidOcc);
    }
}
