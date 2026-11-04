using System.Numerics;
using System.Runtime.CompilerServices;
using System.Runtime.Intrinsics;

namespace SuffixTree.Esa;

/// <summary>Which resolution produced the result of a direct LCP attempt (<see cref="DirectReport"/>).</summary>
internal enum DirectOutcome
{
    /// <summary>The attempt gave up (too many long ranks) and produced nothing: the three-pass kernel is the better algorithm.</summary>
    GaveUp,

    /// <summary>Every long rank was resolved on its own, vector by vector, within the step budget.</summary>
    ResolvedPerRank,

    /// <summary>The step budget was exhausted (a repeat is long) and the long ranks were resolved in text order with the Kasai carry.</summary>
    ResolvedInTextOrder,
}

/// <summary>What the direct LCP kernel did, for tests and diagnostics (production callers pass none); <see cref="ForceTextOrder"/> is the only input.</summary>
internal sealed class DirectReport
{
    /// <summary>Input (test hook): skip the per-rank resolution and resolve in text order whatever the text.</summary>
    internal bool ForceTextOrder { get; init; }

    /// <summary>How the attempt ended.</summary>
    internal DirectOutcome Outcome { get; set; }

    /// <summary>Number of long ranks the scan found (0 when it gave up).</summary>
    internal int LongRanks { get; set; }

    /// <summary>Number of distinct text positions among the long ranks (text-order resolution only; equal to <see cref="LongRanks"/> for a valid suffix array).</summary>
    internal int LongPositions { get; set; }

    /// <summary>Vector steps the per-rank resolution used (all of the budget plus one when it was exhausted).</summary>
    internal long PerRankSteps { get; set; }
}

/// <summary>
/// Iteration 1o.3: the LCP array computed in suffix-array order, in one sequential pass, instead of the three passes of Phi / PLCP / gather (whose
/// random accesses into two 4-byte-per-character arrays cost about 0.25 s of the 0.33 s at 10 Mbp: the phase is bound by cache misses, not by instructions).
/// <para>
/// <b>Scan.</b> For r = 1..n the suffixes <c>SA[r-1]</c> and <c>SA[r]</c> are compared over one 16-byte vector (8 symbols for 16-bit texts): one random access into
/// the text per rank (1 byte per character, a quarter of the footprint of the old passes), prefetched from the suffix array a fixed distance ahead, the result is the
/// number of leading equal symbols from a compare mask, with no branch on the data. A value below the vector width is exact and is stored. A full-width match is
/// recorded as "long" (its rank goes to a list, the byte holds the vector width as a placeholder). On random DNA about 0.1 % of the ranks are long.
/// </para>
/// <para>
/// <b>Per-rank resolution (the fast path).</b> Each long rank is resolved by continuing the comparison from the vector width, vector by vector. This is not amortised (Kasai's
/// argument needs text order): a repeat of length L costs about L^2 / 32 vector steps. It is the cheapest way when the long matches are short, which is the case on texts
/// without long repeats (random DNA: no step at all for 99.9 % of the long ranks), so it runs first, with an allowance of <see cref="ResolveStepsPerLong"/> steps per long rank
/// (plus <see cref="ResolveBaseSteps"/>, cumulative: it is checked against the ranks processed so far). Iteration 1o.3b gave the whole attempt up when this budget (then 2 n steps) ran
/// out, which one repeat longer than about 8 sqrt(n) symbols did (REVIEW-S1a-opt MAJ-1). With <see cref="TextOrderMinLongs"/> or more long ranks the per-rank loop is not tried at all.
/// </para>
/// <para>
/// <b>Resolution in text order (iteration 1o.3c, when the budget runs out).</b> The long ranks are sorted by text position with a bitmap (set the bit of every position, prefix
/// counts per 64-bit word, scatter the records to their index: a counting sort in one pass, no comparisons, the bitmap is 1/8 byte per character) and walked in increasing text
/// position with the Kasai invariant: the LCP of position <c>j</c> with its suffix-array predecessor is at least <c>PLCP[j-1] - 1</c> (<see cref="LcpReference.ComputePlcp{TReader}"/>), so when the
/// previous long position is <c>j - 1</c> the comparison starts at <c>max(PLCP[j-1] - 1, w)</c>, otherwise at <c>w</c> (the scan has seen the first <c>w</c> symbols, and a position whose
/// predecessor in text order is not long has a PLCP below <c>w</c>). This is the Kasai run restricted to the long positions: the carried length drops by at most one per position, and a
/// run of long positions ends only where the length is back at <c>w</c>, so the symbols compared in all runs together number at most the long positions plus the lengths of the
/// runs' first extensions: <b>O(n / w + long ranks) vector steps for every text</b>; a repeat of length L costs L / w steps once and 1 step per position. The values go to an array
/// indexed by the order of the long ranks (rank ascending), then the kind is chosen from the exact counts (<see cref="CompactLcp.SelectKind"/>) and the long entries are written in
/// rank order into the byte array (kind <c>U8</c> with the escape list) or into a widened <c>ushort</c> / <c>int</c> array (kinds <c>U16</c> / <c>U32</c>), so the kind never forces a fallback.
/// Its memory is 20 bytes per long rank (rank list, text position, predecessor, sorted record), which is why the number of long ranks is capped.
/// </para>
/// <para>
/// <b>When it gives up.</b> The scan gives up when more than a quarter of the ranks are long, or after <c>max(n / 32, EarlyAbortMinRanks)</c> ranks when more than half of them
/// are long (a repeated text, a tiny alphabet): the three-pass kernel (<see cref="Build{TSym}"/>, linear in the worst case) is then the better algorithm and runs, so the result is always
/// identical to <see cref="LcpReference"/> for a valid suffix array (same kind, same escape list, every entry) and the time is O(n).
/// </para>
/// <para>
/// Why two resolutions: the text-order machinery is about 1,000 bytes of IL more to compile than the per-rank loop, and a cold process pays for it (about 4 ms of JIT at 1 Mbp, 6 % of the build;
/// measured in attempt 1 of 1o.3c), while on the common input (no repeat of more than a few hundred symbols) the per-rank loop is cheaper. The text-order code is compiled only when it runs.
/// </para>
/// <para>
/// Memory safety: every suffix-array entry is range-checked when it is read, the vector loads happen only when both suffixes have at least a vector of symbols left
/// (<c>n - max(a, b) &gt;= width</c>), the scalar compare and the continuation never pass <c>n - max(a, b)</c>, the start length of a text-order continuation is clamped to that limit.
/// For an input that is not a permutation of the suffix positions the result is unspecified (it can differ from the reference, which is garbage in garbage out as well; a repeated
/// text position keeps the minimum <c>w</c> for the records it displaced), the access pattern stays inside the buffers.
/// </para>
/// </summary>
internal static unsafe partial class LcpKernel
{
    /// <summary>The scan gives up (and the three-pass algorithm runs) when more than this fraction of the ranks, 1 / <see cref="LongFractionDivisor"/>, is long.</summary>
    internal const int LongFractionDivisor = 4;

    /// <summary>
    /// Early exit of the scan on repetitive texts: once at least <c>max(n / 32, EarlyAbortMinRanks)</c> ranks are done and more than half of them are long, the scan gives up at once
    /// (the quarter limit would only be reached after about 30 % of the pass; on a text of 100 copies of a block that is 27 ms of wasted work at 10 Mbp).
    /// </summary>
    internal const int EarlyAbortMinRanks = 1024;

    /// <summary>
    /// Allowance of the per-rank resolution in vector steps for every long rank it has processed (random DNA needs 0.0 to 0.1; a repeat of length L needs L^2 / 32 for L ranks). The
    /// allowance is cumulative (<see cref="ResolveBaseSteps"/> plus this per processed rank), so the attempt is abandoned within a few thousand ranks when the long matches are
    /// long, wherever they are in the rank order: every processed rank costs about 100 ns (two random cache lines), a budget per rank of the whole list wasted 0.11 s at 10 Mbp.
    /// </summary>
    internal const int ResolveStepsPerLong = 4;

    /// <summary>Allowance of the per-rank resolution in vector steps before the first rank.</summary>
    internal const int ResolveBaseSteps = 4096;

    /// <summary>
    /// With at least this many long ranks the text-order resolution is used at once: it costs about 30 ns per long rank, the per-rank loop about 100 ns (two random cache
    /// lines each, however short the match), and the text-order code costs about 4 ms of JIT in a cold process, which is what 40,000 per-rank resolutions cost. Random DNA has
    /// 11,000 long ranks at 10 Mbp (per-rank) and 2.3 million at 100 Mbp (text order).
    /// </summary>
    internal const int TextOrderMinLongs = 1 << 16;

    /// <summary>Builds the compact LCP in suffix-array order with a three-pass fallback; <paramref name="trace"/> gets the scan in <c>PlcpTicks</c>, the per-rank resolution in <c>CountTicks</c> and, for the text-order resolution, its output in <c>OutputTicks</c>.</summary>
    internal static CompactLcp BuildDirect<TSym>(Buf<TSym> text, int n, Buf<int> sa, NativeArena arena, LcpTrace? trace = null, DirectReport? report = null)
        where TSym : unmanaged, IEquatable<TSym>
    {
        ArgumentNullException.ThrowIfNull(arena);
        if (sa[0] != n) throw new ArgumentException("sa[0] must be n, the empty suffix", nameof(sa));
        CompactLcp? direct = TryDirect(text, n, sa, arena, trace, report);
        return direct ?? Build(text, n, sa, arena, trace);
    }

    /// <summary>The direct algorithm; null when it gave up (see the type remarks), in which case nothing but arena memory was used (and it has been released again).</summary>
    internal static CompactLcp? TryDirect<TSym>(Buf<TSym> text, int n, Buf<int> sa, NativeArena arena, LcpTrace? trace = null, DirectReport? report = null)
        where TSym : unmanaged, IEquatable<TSym>
    {
        int longCap = n / LongFractionDivisor;
        int escapeCap = n / 16; // 16 * escapes <= n: the U8 rule of CompactLcp.SelectKind
        int mark = arena.Mark();
        Buf<int> longRanks = arena.Allocate<int>(longCap + 1, zero: false, "longRanks");
        int rankMark = arena.Mark();
        Buf<int> escapeRanks = arena.Allocate<int>(escapeCap + 1, zero: false, "escapeRanks");
        Buf<int> escapeValues = arena.Allocate<int>(escapeCap + 1, zero: false, "escapeValues");
        byte[] lcp = GC.AllocateUninitializedArray<byte>(n + 1);

        long tick = trace is null ? 0 : System.Diagnostics.Stopwatch.GetTimestamp();
        CompactLcp result;
        fixed (byte* p = lcp)
        {
            var dst = new Buf<byte>(p, n + 1, "lcp8");
            int longs = ScanDirect(text, n, sa, dst, longRanks, longCap);
            if (trace is not null)
            {
                long now = System.Diagnostics.Stopwatch.GetTimestamp();
                trace.PlcpTicks += now - tick;
                tick = now;
            }

            if (longs < 0)
            {
                arena.ReleaseTo(mark);
                if (report is not null) report.Outcome = DirectOutcome.GaveUp;
                return null;
            }

            int escapes = -1;
            long steps = 0;
            if (longs < TextOrderMinLongs && (report is null || !report.ForceTextOrder))
                escapes = ResolveLong(text, n, sa, dst, longRanks, longs, escapeRanks, escapeValues, escapeCap, out steps);
            if (report is not null)
            {
                report.LongRanks = longs;
                report.PerRankSteps = steps;
            }

            if (trace is not null)
            {
                long now = System.Diagnostics.Stopwatch.GetTimestamp();
                trace.CountTicks += now - tick;
                tick = now;
            }

            if (escapes >= 0)
            {
                var overflowRank = new int[escapes];
                var overflowValue = new int[escapes];
                escapeRanks.AsSpan(escapes).CopyTo(overflowRank);
                escapeValues.AsSpan(escapes).CopyTo(overflowValue);
                result = CompactLcp.FromParts(LcpKind.U8, n + 1, lcp, null, null, overflowRank, overflowValue);
                if (report is not null) report.Outcome = DirectOutcome.ResolvedPerRank;
                if (trace is not null) trace.OutputTicks += System.Diagnostics.Stopwatch.GetTimestamp() - tick;
            }
            else
            {
                arena.ReleaseTo(rankMark); // the escape lists of the failed attempt are dead
                result = ResolveInTextOrder(text, n, sa, arena, dst, lcp, longs, longRanks, report);
                if (trace is not null) trace.OutputTicks += System.Diagnostics.Stopwatch.GetTimestamp() - tick;
            }
        }

        arena.ReleaseTo(mark);
        return result;
    }

    /// <summary>
    /// The scan: <c>dst[r]</c> = the common prefix of <c>SA[r-1]</c> and <c>SA[r]</c> when it is shorter than the vector width, else the width (a placeholder) and the
    /// rank is appended to <paramref name="longRanks"/>. Returns the number of long ranks, or -1 when there would be more than <paramref name="longCap"/> or the early exit (<see cref="EarlyAbortMinRanks"/>) fires.
    /// </summary>
    /// <exception cref="ArgumentException">An entry <c>sa[r]</c>, <c>r &gt;= 1</c>, is outside <c>[0, n)</c>.</exception>
    private static int ScanDirect<TSym>(Buf<TSym> text, int n, Buf<int> sa, Buf<byte> dst, Buf<int> longRanks, int longCap)
        where TSym : unmanaged, IEquatable<TSym>
    {
        int width = Vector128<TSym>.Count;
        int distance = LcpPrefetch.Direct;
        int previous = n;
        int longs = 0;
        int probe = Math.Max(n >> 5, EarlyAbortMinRanks);
        dst[0] = 0;
        for (int r = 1; r <= n; r++)
        {
            int s = sa[r];
            int ahead = sa[Math.Min(r + distance, n)];
            text.Prefetch(ahead);
            text.Prefetch(ahead + (width - 1)); // the vector at `ahead` straddles two cache lines in one case of four (a miss there is not covered by the first prefetch)
            if ((uint)s >= (uint)n) throw BadEntry(r, s, n, nameof(sa));
            int top = previous ^ ((previous ^ s) & ((previous - s) >> 31)); // max(previous, s) without a branch: which one is larger is a coin flip
            int limit = n - top;
            int v;
            if (limit >= width)
            {
                v = EqualPrefix(text, previous, s);
                if (v == width)
                {
                    if (longs == longCap || (r >= probe && longs > (r >> 1))) return -1;
                    longRanks[longs++] = r;
                }
            }
            else
            {
                v = CommonPrefix(text, previous, s, limit);
            }

            dst[r] = (byte)v;
            previous = s;
        }

        return longs;
    }

    /// <summary>
    /// Continues the comparison of every long rank from the vector width, writes the exact value (or the escape 255 plus a list entry). Returns the number of escapes,
    /// or -1 when the cumulative step allowance (<see cref="ResolveBaseSteps"/> + <see cref="ResolveStepsPerLong"/> per processed rank) is exhausted or more than
    /// <paramref name="escapeCap"/> escapes would be needed; <paramref name="steps"/> is the number of full-vector steps taken (the first one over the allowance included).
    /// </summary>
    internal static int ResolveLong<TSym>(Buf<TSym> text, int n, Buf<int> sa, Buf<byte> dst, Buf<int> longRanks, int longs, Buf<int> escapeRanks, Buf<int> escapeValues, int escapeCap, out long steps)
        where TSym : unmanaged, IEquatable<TSym>
    {
        int width = Vector128<TSym>.Count;
        int distance = LcpPrefetch.Resolve;
        int escapes = 0;
        long work = 0;
        long budget = ResolveBaseSteps;
        steps = 0;
        for (int k = 0; k < longs; k++)
        {
            budget += ResolveStepsPerLong;
            if (k + distance < longs)
            {
                int ahead = longRanks[k + distance];
                text.Prefetch(sa[ahead - 1] + width);
                text.Prefetch(sa[ahead] + width);
            }
            int r = longRanks[k];
            int a = sa[r - 1], b = sa[r];
            int limit = n - Math.Max(a, b);
            int i = width; // the first `width` symbols are equal: the scan saw them
            while (true)
            {
                if (i + width > limit)
                {
                    while (i < limit && text[a + i].Equals(text[b + i])) i++;
                    break;
                }
                int c = EqualPrefix(text, a + i, b + i);
                i += c;
                if (c < width) break;
                if (++work > budget)
                {
                    steps = work;
                    return -1;
                }
            }
            if (i < byte.MaxValue)
            {
                dst[r] = (byte)i;
                continue;
            }
            if (escapes == escapeCap)
            {
                steps = work;
                return -1;
            }
            dst[r] = byte.MaxValue;
            escapeRanks[escapes] = r;
            escapeValues[escapes] = i;
            escapes++;
        }
        steps = work;
        return escapes;
    }

    /// <summary>
    /// The resolution in text order (see the type remarks): gathers the text position and the predecessor of every long rank, sorts the positions with a bitmap, walks them in
    /// increasing order with the Kasai carry, chooses the kind from the exact counts and writes the long entries in rank order. The workspace comes from <paramref name="arena"/>
    /// and is released again.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    private static CompactLcp ResolveInTextOrder<TSym>(Buf<TSym> text, int n, Buf<int> sa, NativeArena arena, Buf<byte> dst, byte[] lcp8, int longs, Buf<int> longRanks, DirectReport? report)
        where TSym : unmanaged, IEquatable<TSym>
    {
        int width = Vector128<TSym>.Count;
        int words = (n + 63) >> 6;
        int mark = arena.Mark();
        Buf<int> values = arena.Allocate<int>(longs + 1, zero: false, "longValues"); // first the text positions, then (same memory) the LCP values
        int workMark = arena.Mark();
        Buf<int> predecessors = arena.Allocate<int>(longs + 1, zero: false, "longPredecessors");
        Buf<ulong> bits = arena.Allocate<ulong>(words + 1, zero: true, "positionBits");
        Buf<int> before = arena.Allocate<int>(words + 1, zero: false, "bitsBefore");
        Buf<ulong> sorted = arena.Allocate<ulong>(longs + 1, zero: false, "sortedLong");

        // 1. the text position and the predecessor of every long rank (the scan checked that they are inside the text and the suffix array is not modified during a build; the
        //    check is repeated because a wrong value here would be an out-of-bounds write into the bitmap), and the bitmap of the positions
        int lookAhead = LcpPrefetch.Resolve;
        for (int k = 0; k < longs; k++)
        {
            if (k + lookAhead < longs) sa.Prefetch(longRanks[k + lookAhead] - 1);
            int r = longRanks[k];
            int j = sa[r], p = sa[r - 1];
            if ((uint)j >= (uint)n) throw BadEntry(r, j, n, nameof(sa));
            if ((uint)p >= (uint)n) throw BadEntry(r - 1, p, n, nameof(sa));
            values[k] = j;
            predecessors[k] = p;
            bits[j >> 6] |= 1UL << (j & 63);
        }

        // 2. the number of set bits before every word
        int running = 0;
        for (int w = 0; w < words; w++)
        {
            before[w] = running;
            running += BitOperations.PopCount(bits[w]);
        }

        int distinct = running;

        // 3. counting sort: the record (predecessor, k) goes to the index of its position among the set bits; the result array starts at the minimum every long rank has (the width)
        for (int k = 0; k < longs; k++)
        {
            int j = values[k];
            int w = j >> 6;
            int index = before[w] + BitOperations.PopCount(bits[w] & ((1UL << (j & 63)) - 1));
            sorted[index] = (uint)predecessors[k] | ((ulong)(uint)k << 32);
            values[k] = width;
        }

        // 4. the long positions in increasing order, each continuing from the carried length
        int t = 0;
        int previousPosition = -2, previousLength = 0;
        int ahead = LcpPrefetch.Resolve;
        for (int w = 0; w < words; w++)
        {
            ulong word = bits[w];
            while (word != 0)
            {
                int j = (w << 6) + BitOperations.TrailingZeroCount(word);
                word &= word - 1;
                ulong record = sorted[t];
                if (t + ahead < distinct) text.Prefetch((int)(uint)sorted[t + ahead] + width);
                t++;
                int p = (int)(uint)record;
                int k = (int)(record >> 32);
                int top = j ^ ((j ^ p) & ((j - p) >> 31)); // max(j, p) without a branch
                int limit = n - top;
                int l = j == previousPosition + 1 ? Math.Max(previousLength - 1, width) : width;
                if (l > limit) l = limit; // never needed for a valid suffix array (the true value is at most the limit); keeps a corrupt one inside the text
                while (true)
                {
                    if (l + width > limit)
                    {
                        while (l < limit && text[j + l].Equals(text[p + l])) l++;
                        break;
                    }

                    int c = EqualPrefix(text, j + l, p + l);
                    l += c;
                    if (c < width) break;
                }

                values[k] = l;
                previousPosition = j;
                previousLength = l;
            }
        }

        arena.ReleaseTo(workMark); // predecessors, bitmap, counts and sorted records are dead
        if (report is not null)
        {
            report.Outcome = DirectOutcome.ResolvedInTextOrder;
            report.LongPositions = distinct;
        }

        CompactLcp result = Materialize(n, dst, lcp8, longs, longRanks, values);
        arena.ReleaseTo(mark);
        return result;
    }

    /// <summary>
    /// Chooses the kind from the exact counts of values at or above 255 and 65,535 among the long ranks (the others are below the vector width), then writes the long entries in
    /// rank order: into the byte array (<c>U8</c>, with the escape list), or into a widened <c>ushort</c> / <c>int</c> copy of it.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    private static CompactLcp Materialize(int n, Buf<byte> dst, byte[] lcp8, int longs, Buf<int> longRanks, Buf<int> values)
    {
        long escapes8 = 0, escapes16 = 0;
        for (int k = 0; k < longs; k++)
        {
            int v = values[k];
            if (v >= byte.MaxValue)
            {
                escapes8++;
                if (v >= ushort.MaxValue) escapes16++;
            }
        }

        LcpKind kind = CompactLcp.SelectKind(n, escapes8, escapes16);
        int escapes = kind switch { LcpKind.U8 => (int)escapes8, LcpKind.U16 => (int)escapes16, _ => 0 };
        var overflowRank = new int[escapes];
        var overflowValue = new int[escapes];
        int e = 0;
        switch (kind)
        {
            case LcpKind.U8:
                for (int k = 0; k < longs; k++)
                {
                    int v = values[k];
                    int r = longRanks[k];
                    if (v < byte.MaxValue)
                    {
                        dst[r] = (byte)v;
                        continue;
                    }

                    dst[r] = byte.MaxValue;
                    overflowRank[e] = r;
                    overflowValue[e] = v;
                    e++;
                }

                return CompactLcp.FromParts(kind, n + 1, lcp8, null, null, overflowRank, overflowValue);
            case LcpKind.U16:
                {
                    ushort[] lcp16 = GC.AllocateUninitializedArray<ushort>(n + 1);
                    Widen(lcp8, lcp16);
                    for (int k = 0; k < longs; k++)
                    {
                        int v = values[k];
                        int r = longRanks[k];
                        if (v < ushort.MaxValue)
                        {
                            lcp16[r] = (ushort)v;
                            continue;
                        }

                        lcp16[r] = ushort.MaxValue;
                        overflowRank[e] = r;
                        overflowValue[e] = v;
                        e++;
                    }

                    return CompactLcp.FromParts(kind, n + 1, null, lcp16, null, overflowRank, overflowValue);
                }

            default:
                {
                    int[] lcp32 = GC.AllocateUninitializedArray<int>(n + 1);
                    Widen(lcp8, lcp32);
                    for (int k = 0; k < longs; k++) lcp32[longRanks[k]] = values[k];
                    return CompactLcp.FromParts(kind, n + 1, null, null, lcp32, overflowRank, overflowValue);
                }
        }
    }

    /// <summary><c>dst[i] = src[i]</c> with 16 bytes per step (safe code; the lengths are checked by the slices).</summary>
    private static void Widen(ReadOnlySpan<byte> src, Span<ushort> dst)
    {
        int i = 0;
        for (; i + 16 <= src.Length; i += 16)
        {
            (Vector128<ushort> lower, Vector128<ushort> upper) = Vector128.Widen(Vector128.Create(src.Slice(i, 16)));
            lower.CopyTo(dst.Slice(i, 8));
            upper.CopyTo(dst.Slice(i + 8, 8));
        }

        for (; i < src.Length; i++) dst[i] = src[i];
    }

    /// <summary><c>dst[i] = src[i]</c> with 16 bytes per step (safe code; the lengths are checked by the slices).</summary>
    private static void Widen(ReadOnlySpan<byte> src, Span<int> dst)
    {
        int i = 0;
        for (; i + 16 <= src.Length; i += 16)
        {
            (Vector128<ushort> lower, Vector128<ushort> upper) = Vector128.Widen(Vector128.Create(src.Slice(i, 16)));
            (Vector128<uint> a, Vector128<uint> b) = Vector128.Widen(lower);
            (Vector128<uint> c, Vector128<uint> d) = Vector128.Widen(upper);
            a.AsInt32().CopyTo(dst.Slice(i, 4));
            b.AsInt32().CopyTo(dst.Slice(i + 4, 4));
            c.AsInt32().CopyTo(dst.Slice(i + 8, 4));
            d.AsInt32().CopyTo(dst.Slice(i + 12, 4));
        }

        for (; i < src.Length; i++) dst[i] = src[i];
    }

    /// <summary>The number of leading equal symbols of the two vectors at <paramref name="a"/> and <paramref name="b"/>, the vector width when all are equal. Both vectors must lie inside the text.</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static int EqualPrefix<TSym>(Buf<TSym> t, int a, int b)
        where TSym : unmanaged
    {
        uint equal = Vector128.Equals(t.LoadVector128(a), t.LoadVector128(b)).ExtractMostSignificantBits();
        return BitOperations.TrailingZeroCount(~equal); // the bits above the vector width are zero in `equal`, so an all-equal vector gives exactly the width
    }

    [MethodImpl(MethodImplOptions.NoInlining)]
    private static ArgumentException BadEntry(int rank, int value, int n, string parameter)
        => new($"sa[{rank}] = {value} is outside [0, {n})", parameter);
}
