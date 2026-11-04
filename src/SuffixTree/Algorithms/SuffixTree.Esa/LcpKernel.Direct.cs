using System.Numerics;
using System.Runtime.CompilerServices;
using System.Runtime.Intrinsics;

namespace SuffixTree.Esa;

/// <summary>How an attempt of the direct LCP kernel ended (<see cref="DirectReport"/>).</summary>
internal enum DirectOutcome
{
    /// <summary>The compact LCP was produced.</summary>
    Succeeded,

    /// <summary>A sample of the ranks showed that clearly more than a quarter of them are long: nothing was allocated or scanned, the three-pass kernel is the better algorithm.</summary>
    GaveUpBySample,

    /// <summary>The scan met more than a quarter of long ranks (the sample had been too optimistic): the scan was abandoned.</summary>
    GaveUpByCap,
}

/// <summary>What the direct LCP kernel did, for tests and diagnostics (not used by production callers).</summary>
internal sealed class DirectReport
{
    /// <summary>How the attempt ended.</summary>
    internal DirectOutcome Outcome { get; set; }

    /// <summary>Number of long ranks the scan found (0 when it did not run to the end).</summary>
    internal int LongRanks { get; set; }

    /// <summary>Number of distinct text positions among the long ranks (equal to <see cref="LongRanks"/> for a valid suffix array).</summary>
    internal int LongPositions { get; set; }
}

/// <summary>
/// Iteration 1o.3: the LCP array computed in suffix-array order, in one sequential pass, instead of the three passes of Phi / PLCP / gather (whose
/// random accesses into two 4-byte-per-character arrays cost about 0.25 s of the 0.33 s at 10 Mbp: the phase is bound by cache misses, not by instructions).
/// Iteration 1o.3c (REVIEW-S1a-opt MAJ-1, MIN-4) replaced the bounded vector-by-vector continuation by the Kasai argument restricted to the long ranks, so a repeat of any
/// length costs only its own length and the kernel no longer gives up on texts with long repeats.
/// <para>
/// <b>Scan.</b> For r = 1..n the suffixes <c>SA[r-1]</c> and <c>SA[r]</c> are compared over one 16-byte vector (8 symbols for 16-bit texts): one random access into
/// the text per rank (1 byte per character, a quarter of the footprint of the old passes), prefetched from the suffix array a fixed distance ahead, the result is the
/// number of leading equal symbols from a compare mask, with no branch on the data. A value below the vector width <c>w</c> is exact and is stored. A full-width match is
/// a "long" rank: the byte holds <c>w</c> as a placeholder and the triple (rank, text position, position of the predecessor) goes to three lists. On random DNA about 0.1 % of
/// the ranks are long, on a text with a repeat of 1 % of its length 1 %, on a text with many copied blocks up to a quarter and more.
/// </para>
/// <para>
/// <b>Resolution in text order.</b> The long ranks are sorted by text position with a bitmap (set the bit of every position, prefix counts per 64-bit word, scatter the
/// records to their index: a counting sort in one pass, no comparisons, the bitmap is 1/8 byte per character) and then walked in increasing text position with the Kasai
/// invariant: the LCP of position <c>j</c> with its suffix-array predecessor is at least <c>PLCP[j-1] - 1</c> (<see cref="LcpReference.ComputePlcp{TReader}"/>), so when the previous long position is
/// <c>j - 1</c> the comparison starts at <c>max(PLCP[j-1] - 1, w)</c>, otherwise at <c>w</c> (the scan has seen the first <c>w</c> symbols, and a position whose predecessor in text
/// order is not long has a PLCP below <c>w</c>). This is the Kasai run restricted to the long positions: the carried length drops by at most one per position, and a run of long
/// positions ends only where the length is back at <c>w</c>, so the symbols compared in all runs together number at most the long positions plus the lengths of the runs' first
/// extensions: <b>O(n / w + long ranks) vector steps for every text</b>, a repeat of length L costs L / w steps once and 1 step per position, not L^2 / (2 w). There is no budget and
/// no fallback for long repeats. The values go to an array indexed by the order of the long ranks (rank ascending), then the kind is chosen from the exact counts (<see cref="CompactLcp.SelectKind"/>)
/// and the long entries are written in rank order into the byte array (kind <c>U8</c> with the escape list) or into a widened <c>ushort</c> / <c>int</c> array (kinds <c>U16</c> / <c>U32</c>).
/// </para>
/// <para>
/// <b>When it does not run.</b> The memory of the resolution is 20 bytes per long rank (the three lists, the sorted records), so the number of long ranks is capped at a quarter of n
/// (5 bytes per character at the cap). More long ranks mean a highly repetitive text (a repeat of a quarter of the text, many copies of a block, a tiny alphabet) on which the
/// three-pass kernel is the better algorithm anyway. The decision is made before anything is allocated, from a stratified sample of <see cref="SampleRanks"/> ranks (exact for n up to that
/// number): the attempt is abandoned when the sample estimates more than 9/8 of the cap. A text whose sample was too optimistic is stopped by the cap in the scan itself (the cost is one scan,
/// about 6 ns per character, bounded). In every such case the three-pass kernel runs (<see cref="Build{TSym}"/>, linear in the worst case), so the result is always identical to <see cref="LcpReference"/> for a
/// valid suffix array (same kind, same escape list, every entry) and the time is O(n).
/// </para>
/// <para>
/// Memory safety: every suffix-array entry is range-checked when it is read, the vector loads happen only when both suffixes have at least a vector of symbols left
/// (<c>n - max(a, b) &gt;= width</c>), the scalar compare and the continuation never pass <c>n - max(a, b)</c>, the start length of a continuation is clamped to that limit. For an input that
/// is not a permutation of the suffix positions the result is unspecified (it can differ from the reference, which is garbage in garbage out as well; a repeated text position keeps
/// the minimum <c>w</c> for the records it displaced), the access pattern stays inside the buffers.
/// </para>
/// </summary>
internal static unsafe partial class LcpKernel
{
    /// <summary>The attempt is abandoned when more than this fraction of the ranks, 1 / <see cref="LongFractionDivisor"/>, is long (memory: 20 bytes per long rank).</summary>
    internal const int LongFractionDivisor = 4;

    /// <summary>Ranks examined by the up-front sample (stratified: one pseudo-random rank per block of n / SampleRanks); a text with fewer ranks is examined completely.</summary>
    internal const int SampleRanks = 4096;

    /// <summary>The sample gives up when its estimate exceeds the cap by more than 1 / this (the statistical error of 4,096 ranks is under 0.8 % points).</summary>
    internal const int SampleSlackDivisor = 8;

    /// <summary>Builds the compact LCP in suffix-array order with a three-pass fallback; <paramref name="trace"/> gets the scan in <c>PlcpTicks</c>, the resolution in <c>CountTicks</c> and the output in <c>OutputTicks</c>.</summary>
    internal static CompactLcp BuildDirect<TSym>(Buf<TSym> text, int n, Buf<int> sa, NativeArena arena, LcpTrace? trace = null)
        where TSym : unmanaged, IEquatable<TSym>
    {
        ArgumentNullException.ThrowIfNull(arena);
        if (sa[0] != n) throw new ArgumentException("sa[0] must be n, the empty suffix", nameof(sa));
        CompactLcp? direct = TryDirect(text, n, sa, arena, trace);
        return direct ?? Build(text, n, sa, arena, trace);
    }

    /// <summary>The direct algorithm; null when it gave up (see the type remarks), in which case nothing but arena memory was used (and it has been released again).</summary>
    internal static CompactLcp? TryDirect<TSym>(Buf<TSym> text, int n, Buf<int> sa, NativeArena arena, LcpTrace? trace = null, DirectReport? report = null)
        where TSym : unmanaged, IEquatable<TSym>
    {
        int longCap = n / LongFractionDivisor;
        long tick = trace is null ? 0 : System.Diagnostics.Stopwatch.GetTimestamp();
        if (SampleSaysTooManyLongRanks(text, n, sa, longCap))
        {
            if (trace is not null) trace.PlcpTicks += System.Diagnostics.Stopwatch.GetTimestamp() - tick;
            if (report is not null) report.Outcome = DirectOutcome.GaveUpBySample;
            return null;
        }

        int mark = arena.Mark();
        Buf<int> longRanks = arena.Allocate<int>(longCap + 1, zero: false, "longRanks");
        Buf<int> longPositions = arena.Allocate<int>(longCap + 1, zero: false, "longPositions");
        int resolutionMark = arena.Mark(); // everything allocated from here on is dead when the values are known: it is released before the output arrays are allocated
        Buf<int> longPredecessors = arena.Allocate<int>(longCap + 1, zero: false, "longPredecessors");
        byte[] lcp = GC.AllocateUninitializedArray<byte>(n + 1);
        CompactLcp result;
        fixed (byte* p = lcp)
        {
            var dst = new Buf<byte>(p, n + 1, "lcp8");
            int longs = ScanDirect(text, n, sa, dst, longRanks, longPositions, longPredecessors, longCap);
            if (trace is not null)
            {
                long now = System.Diagnostics.Stopwatch.GetTimestamp();
                trace.PlcpTicks += now - tick;
                tick = now;
            }
            if (longs < 0)
            {
                arena.ReleaseTo(mark);
                if (report is not null) report.Outcome = DirectOutcome.GaveUpByCap;
                return null;
            }

            Buf<int> values = ResolveInTextOrder(text, n, arena, longs, longPositions, longPredecessors, out int distinct);
            arena.ReleaseTo(resolutionMark);
            if (trace is not null)
            {
                long now = System.Diagnostics.Stopwatch.GetTimestamp();
                trace.CountTicks += now - tick;
                tick = now;
            }

            result = Materialize(n, dst, lcp, longs, longRanks, values);
            if (report is not null)
            {
                report.Outcome = DirectOutcome.Succeeded;
                report.LongRanks = longs;
                report.LongPositions = distinct;
            }
        }

        arena.ReleaseTo(mark);
        if (trace is not null) trace.OutputTicks += System.Diagnostics.Stopwatch.GetTimestamp() - tick;
        return result;
    }

    /// <summary>
    /// True when a sample of ranks estimates more long ranks than the cap allows (plus the slack of <see cref="SampleSlackDivisor"/>); exact (no slack) when the text has at most
    /// <see cref="SampleRanks"/> ranks. The sample takes one pseudo-random rank out of each of <see cref="SampleRanks"/> equal blocks of the ranks, so a rank region that is
    /// entirely long is seen in proportion to its size whatever its position.
    /// </summary>
    /// <exception cref="ArgumentException">A sampled entry <c>sa[r]</c>, <c>r &gt;= 1</c>, is outside <c>[0, n)</c>.</exception>
    private static bool SampleSaysTooManyLongRanks<TSym>(Buf<TSym> text, int n, Buf<int> sa, int longCap)
        where TSym : unmanaged, IEquatable<TSym>
    {
        if (n < 1) return false;
        int width = Vector128<TSym>.Count;
        int samples = Math.Min(n, SampleRanks);
        int hits = 0;
        for (int i = 0; i < samples; i++)
        {
            long low = (long)i * n / samples, high = (long)(i + 1) * n / samples; // block i of the ranks 1..n is [1 + low, 1 + high), at least one rank
            int r = 1 + (int)low + (int)(unchecked((uint)i * 2654435761u) % (uint)(high - low));
            int s = sa[r], previous = sa[r - 1];
            if ((uint)s >= (uint)n) throw BadEntry(r, s, n, nameof(sa));
            if ((uint)previous > (uint)n) throw BadEntry(r - 1, previous, n, nameof(sa));
            int top = Math.Max(previous, s);
            if (n - top >= width && EqualPrefix(text, previous, s) == width) hits++;
        }

        if (samples == n) return hits > longCap;
        return (long)hits * n > (long)samples * (longCap + (longCap / SampleSlackDivisor));
    }

    /// <summary>
    /// The scan: <c>dst[r]</c> = the common prefix of <c>SA[r-1]</c> and <c>SA[r]</c> when it is shorter than the vector width, else the width (a placeholder) and the
    /// rank, the text position <c>SA[r]</c> and the predecessor <c>SA[r-1]</c> are appended to the three lists. Returns the number of long ranks, or -1 when there would be more than <paramref name="longCap"/>.
    /// </summary>
    /// <exception cref="ArgumentException">An entry <c>sa[r]</c>, <c>r &gt;= 1</c>, is outside <c>[0, n)</c>.</exception>
    private static int ScanDirect<TSym>(Buf<TSym> text, int n, Buf<int> sa, Buf<byte> dst, Buf<int> longRanks, Buf<int> longPositions, Buf<int> longPredecessors, int longCap)
        where TSym : unmanaged, IEquatable<TSym>
    {
        int width = Vector128<TSym>.Count;
        int distance = LcpPrefetch.Direct;
        int previous = n;
        int longs = 0;
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
                    if (longs == longCap) return -1;
                    longRanks[longs] = r;
                    longPositions[longs] = s;
                    longPredecessors[longs] = previous;
                    longs++;
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
    /// Computes the exact LCP of every long rank (the array <c>values[k]</c>, k = index in rank order) by walking the long text positions in increasing order with the Kasai
    /// carry (see the type remarks). <paramref name="distinct"/> is the number of distinct positions (= <paramref name="longs"/> for a valid suffix array). The returned buffer is
    /// the memory of <paramref name="longPositions"/>, which is dead after the sorting pass; the blocks allocated here are the caller's to release (they are dead on return).
    /// </summary>
    private static Buf<int> ResolveInTextOrder<TSym>(Buf<TSym> text, int n, NativeArena arena, int longs, Buf<int> longPositions, Buf<int> longPredecessors, out int distinct)
        where TSym : unmanaged, IEquatable<TSym>
    {
        int width = Vector128<TSym>.Count;
        int words = (n + 63) >> 6;
        Buf<ulong> bits = arena.Allocate<ulong>(words + 1, zero: true, "positionBits");
        Buf<int> before = arena.Allocate<int>(words + 1, zero: false, "bitsBefore");
        Buf<ulong> sorted = arena.Allocate<ulong>(longs + 1, zero: false, "sortedLong");

        // 1. the bitmap of the long positions and the number of set bits before every word
        for (int k = 0; k < longs; k++)
        {
            int j = longPositions[k];
            bits[j >> 6] |= 1UL << (j & 63);
        }

        int running = 0;
        for (int w = 0; w < words; w++)
        {
            before[w] = running;
            running += BitOperations.PopCount(bits[w]);
        }

        distinct = running;

        // 2. counting sort: the record (predecessor, k) goes to the index of its position among the set bits; the result array starts at the minimum every long rank has (the width)
        Buf<int> values = longPositions;
        for (int k = 0; k < longs; k++)
        {
            int j = longPositions[k];
            int w = j >> 6;
            int index = before[w] + BitOperations.PopCount(bits[w] & ((1UL << (j & 63)) - 1));
            sorted[index] = (uint)longPredecessors[k] | ((ulong)(uint)k << 32);
            values[k] = width;
        }

        // 3. the long positions in increasing order, each continuing from the carried length
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

        return values;
    }

    /// <summary>
    /// Chooses the kind from the exact counts of values at or above 255 and 65,535 among the long ranks (the others are below the vector width), then writes the long entries in
    /// rank order: into the byte array (<c>U8</c>, with the escape list), or into a widened <c>ushort</c> / <c>int</c> copy of it.
    /// </summary>
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
