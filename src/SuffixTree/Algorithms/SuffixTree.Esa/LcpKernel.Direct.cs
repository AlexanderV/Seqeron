using System.Numerics;
using System.Runtime.CompilerServices;
using System.Runtime.Intrinsics;

namespace SuffixTree.Esa;

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
/// <b>Resolution.</b> Each long rank is resolved by continuing the comparison from the vector width, vector by vector. This is not amortised (Kasai's argument needs
/// text order), so it is bounded: the scan gives up when more than a quarter of the ranks are long, the resolution gives up after <c>2 n</c> vector steps or when the
/// escapes no longer fit the U8 representation (<c>16 * escapes &gt; n</c>), and in every such case the kernel falls back to the three-pass algorithm
/// (<see cref="Build{TSym}"/>), which is linear in the worst case. The result is therefore always identical to <see cref="LcpReference"/> for a valid suffix array
/// (same kind, same escape list, every entry), and the running time is O(n) in all cases.
/// </para>
/// <para>
/// Memory safety: every suffix-array entry is range-checked when it is read, the vector loads happen only when both suffixes have at least a vector of symbols left
/// (<c>n - max(a, b) &gt;= width</c>), the scalar compare and the continuation never pass <c>n - max(a, b)</c>. For an input that is not a permutation of the
/// suffix positions the result is unspecified (it can differ from the reference, which is garbage in garbage out as well), the access pattern stays inside the buffers.
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

    /// <summary>Budget of the resolution in vector steps, per text character.</summary>
    internal const int ResolveStepsPerChar = 2;

    /// <summary>Builds the compact LCP in suffix-array order with a three-pass fallback; <paramref name="trace"/> gets the scan in <c>PlcpTicks</c> and the resolution in <c>CountTicks</c>.</summary>
    internal static CompactLcp BuildDirect<TSym>(Buf<TSym> text, int n, Buf<int> sa, NativeArena arena, LcpTrace? trace = null)
        where TSym : unmanaged, IEquatable<TSym>
    {
        ArgumentNullException.ThrowIfNull(arena);
        if (sa[0] != n) throw new ArgumentException("sa[0] must be n, the empty suffix", nameof(sa));
        CompactLcp? direct = TryDirect(text, n, sa, arena, trace);
        return direct ?? Build(text, n, sa, arena, trace);
    }

    /// <summary>The direct algorithm; null when it gave up (see the type remarks), in which case nothing but arena memory was used.</summary>
    internal static CompactLcp? TryDirect<TSym>(Buf<TSym> text, int n, Buf<int> sa, NativeArena arena, LcpTrace? trace = null)
        where TSym : unmanaged, IEquatable<TSym>
    {
        int longCap = n / LongFractionDivisor;
        int escapeCap = n / 16; // 16 * escapes <= n: the U8 rule of CompactLcp.SelectKind
        Buf<int> longRanks = arena.Allocate<int>(longCap + 1, zero: false, "longRanks");
        Buf<int> escapeRanks = arena.Allocate<int>(escapeCap + 1, zero: false, "escapeRanks");
        Buf<int> escapeValues = arena.Allocate<int>(escapeCap + 1, zero: false, "escapeValues");
        byte[] lcp = GC.AllocateUninitializedArray<byte>(n + 1);

        long tick = trace is null ? 0 : System.Diagnostics.Stopwatch.GetTimestamp();
        int escapes;
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
            if (longs < 0) return null;
            escapes = ResolveLong(text, n, sa, dst, longRanks, longs, escapeRanks, escapeValues, escapeCap, (long)ResolveStepsPerChar * n);
            if (trace is not null)
            {
                long now = System.Diagnostics.Stopwatch.GetTimestamp();
                trace.CountTicks += now - tick;
                tick = now;
            }
            if (escapes < 0) return null;
        }

        var overflowRank = new int[escapes];
        var overflowValue = new int[escapes];
        escapeRanks.AsSpan(escapes).CopyTo(overflowRank);
        escapeValues.AsSpan(escapes).CopyTo(overflowValue);
        CompactLcp result = CompactLcp.FromParts(LcpKind.U8, n + 1, lcp, null, null, overflowRank, overflowValue);
        if (trace is not null) trace.OutputTicks += System.Diagnostics.Stopwatch.GetTimestamp() - tick;
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
    /// or -1 when the step budget is exhausted or more than <paramref name="escapeCap"/> escapes would be needed.
    /// </summary>
    private static int ResolveLong<TSym>(Buf<TSym> text, int n, Buf<int> sa, Buf<byte> dst, Buf<int> longRanks, int longs, Buf<int> escapeRanks, Buf<int> escapeValues, int escapeCap, long budget)
        where TSym : unmanaged, IEquatable<TSym>
    {
        int width = Vector128<TSym>.Count;
        int distance = LcpPrefetch.Resolve;
        int escapes = 0;
        long work = 0;
        for (int k = 0; k < longs; k++)
        {
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
                if (++work > budget) return -1;
            }
            if (i < byte.MaxValue)
            {
                dst[r] = (byte)i;
                continue;
            }
            if (escapes == escapeCap) return -1;
            dst[r] = byte.MaxValue;
            escapeRanks[escapes] = r;
            escapeValues[escapes] = i;
            escapes++;
        }
        return escapes;
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
