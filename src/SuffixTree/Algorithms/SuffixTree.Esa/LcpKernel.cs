using System.Numerics;
using System.Runtime.CompilerServices;

namespace SuffixTree.Esa;

/// <summary>Which LCP implementation a harness runs: the safe permanent oracle or the pointer kernel (iteration 1o.1a).</summary>
internal enum LcpKernelKind
{
    /// <summary><see cref="LcpReference"/>: safe code, generic over the index type, never optimized.</summary>
    Reference,

    /// <summary><see cref="LcpKernel"/>: <c>int</c>-only pointer kernel behind <see cref="Buf{T}"/>.</summary>
    Pointer,
}

/// <summary>
/// Iteration 1o.1a (REVIEW-S1a section 6): the LCP construction of <see cref="LcpReference"/> as an unsafe <c>int</c>-only kernel. Same three passes and the
/// same result bit for bit (kind, escape list, every entry): Phi over a zeroed native workspace, PLCP in place with the Kasai invariant and the escape counting
/// fused into the loop, and pass 3 written straight into the compact arrays. Every access goes through <see cref="Buf{T}"/> (a raw pointer in Release, bounds-checked
/// under <c>ESA_BOUNDS</c>); the work buffer comes from a <see cref="NativeArena"/>. Differences to the reference that matter for safety:
/// <list type="bullet">
/// <item>the symbol comparison runs over 8-byte words and never reads at or beyond <c>text[n]</c> (the compared length is capped by <c>n - max(j, p)</c> before any access);</item>
/// <item>every suffix-array entry is range-checked once in pass 1 (<c>0 &lt;= sa[i] &lt; n</c> for <c>i &gt;= 1</c>, <c>sa[0] = n</c>), because pass 2 and pass 3 use the entries as indices and the
///   array may come from an untrusted container (S1b); a wrong value throws <see cref="ArgumentException"/>; a non-permutation that passes the range check gives garbage values but is memory safe
///   (Phi is zero-initialised, as the reference's array is);</item>
/// <item>pass 3 cannot write more escapes than were counted: the announced number is checked in the (rare) escape branch, and the final count is compared (<see cref="InvalidOperationException"/> like the reference's builder).</item>
/// </list>
/// </summary>
internal static unsafe class LcpKernel
{
    /// <summary>Builds the compact LCP with the selected implementation.</summary>
    internal static CompactLcp Build(LcpKernelKind kind, PreparedText prepared, int[] sa, LcpTrace? trace = null)
    {
        ArgumentNullException.ThrowIfNull(prepared);
        if (kind == LcpKernelKind.Pointer) return Build(prepared, sa, trace);
        return prepared.Width == SymbolWidth.Byte
            ? LcpReference.Build(prepared.ByteReader(), sa, trace)
            : LcpReference.Build(prepared.UShortReader(), sa, trace);
    }

    /// <summary>Builds the compact LCP of a prepared text from its suffix array (n + 1 layout) with its own workspace.</summary>
    internal static CompactLcp Build(PreparedText prepared, ReadOnlySpan<int> sa, LcpTrace? trace = null)
    {
        using var arena = new NativeArena();
        return Build(prepared, sa, arena, trace);
    }

    /// <summary>Builds the compact LCP with the workspace taken from <paramref name="arena"/> (the tests pass a guarded arena and verify its canaries afterwards).</summary>
    internal static CompactLcp Build(PreparedText prepared, ReadOnlySpan<int> sa, NativeArena arena, LcpTrace? trace = null)
    {
        ArgumentNullException.ThrowIfNull(prepared);
        ArgumentNullException.ThrowIfNull(arena);
        int n = prepared.Length;
        if (sa.Length != n + 1)
            throw new ArgumentException($"the suffix array needs n + 1 = {n + 1} entries, got {sa.Length}", nameof(sa));
        fixed (int* pSa = sa)
        {
            var saBuf = new Buf<int>(pSa, n + 1, "sa");
            if (prepared.Width == SymbolWidth.Byte)
            {
                fixed (byte* pText = prepared.Bytes) return Build(new Buf<byte>(pText, n, "text"), n, saBuf, arena, trace);
            }
            fixed (ushort* pText = prepared.UShorts) return Build(new Buf<ushort>(pText, n, "text"), n, saBuf, arena, trace);
        }
    }

    /// <summary>The kernel on accessor-wrapped input: <paramref name="text"/> has n symbols, <paramref name="sa"/> n + 1 entries.</summary>
    internal static CompactLcp Build<TSym>(Buf<TSym> text, int n, Buf<int> sa, NativeArena arena, LcpTrace? trace = null)
        where TSym : unmanaged, IEquatable<TSym>
    {
        Buf<int> plcp = ComputePlcp(text, n, sa, arena, out long escapes8, out long escapes16, trace);
        long tick = trace is null ? 0 : System.Diagnostics.Stopwatch.GetTimestamp();
        LcpKind kind = CompactLcp.SelectKind(n, escapes8, escapes16);
        int escapes = kind switch { LcpKind.U8 => (int)escapes8, LcpKind.U16 => (int)escapes16, _ => 0 };
        if (trace is not null)
        {
            long now = System.Diagnostics.Stopwatch.GetTimestamp();
            trace.CountTicks += now - tick;
            tick = now;
        }
        CompactLcp result = WriteCompact(plcp, sa, n, kind, escapes);
        if (trace is not null) trace.OutputTicks += System.Diagnostics.Stopwatch.GetTimestamp() - tick;
        return result;
    }

    /// <summary>
    /// Passes 1 and 2: the permuted LCP array (length n) in a zeroed arena block, and the numbers of values at or above 255 and 65,535 (counted in the loop).
    /// </summary>
    /// <exception cref="ArgumentException"><c>sa[0] != n</c> or an entry <c>sa[i]</c>, <c>i &gt;= 1</c>, outside <c>[0, n)</c>.</exception>
    internal static Buf<int> ComputePlcp<TSym>(Buf<TSym> text, int n, Buf<int> sa, NativeArena arena, out long escapes8, out long escapes16, LcpTrace? trace = null)
        where TSym : unmanaged, IEquatable<TSym>
    {
        ArgumentNullException.ThrowIfNull(arena);
        if (sa[0] != n) throw new ArgumentException("sa[0] must be n, the empty suffix", nameof(sa));
        Buf<int> phi = arena.Allocate<int>(n, zero: true, "phi");

        long tick = trace is null ? 0 : System.Diagnostics.Stopwatch.GetTimestamp();
        int bad = FillPhi(sa, phi, n);
        if (bad >= 0) throw new ArgumentException($"sa[{bad}] = {sa[bad]} is outside [0, {n})", nameof(sa));
        if (trace is not null)
        {
            long now = System.Diagnostics.Stopwatch.GetTimestamp();
            trace.PhiTicks += now - tick;
            tick = now;
        }

        (escapes8, escapes16) = PlcpLoop(text, n, phi);
        if (trace is not null) trace.PlcpTicks += System.Diagnostics.Stopwatch.GetTimestamp() - tick;
        return phi;
    }

    // The hot loops are separate methods marked AggressiveOptimization (compiled once to full optimization, no tiering): a loop-heavy method that runs once
    // would otherwise be left in Tier0 code and enter its loop through on-stack replacement, which keeps address-exposed struct locals such as a Buf<T> in the
    // Tier0 frame (a reload of the base pointer on every access) and was measured at 2.3 times the time of the safe reference loop (1o.1a, first attempt).

    /// <summary><c>Phi[SA[i]] = SA[i-1]</c> for <c>i = 1..n</c>; returns the first rank whose entry is outside <c>[0, n)</c>, or -1.</summary>
    private static int FillPhi(Buf<int> sa, Buf<int> phi, int n)
    {
        int previous = n;
        for (int i = 1; i <= n; i++)
        {
            int s = sa[i];
            if ((uint)s >= (uint)n) return i;
            phi[s] = previous;
            previous = s;
        }
        return -1;
    }

    /// <summary>PLCP in place over Phi (Kasai invariant, at most 2n symbol comparisons) with the escape counting fused into the loop.</summary>
    private static (long Escapes8, long Escapes16) PlcpLoop<TSym>(Buf<TSym> text, int n, Buf<int> phi)
        where TSym : unmanaged, IEquatable<TSym>
    {
        long e8 = 0, e16 = 0;
        int l = 0;
        for (int j = 0; j < n; j++)
        {
            int p = phi[j];
            if (p == n)
            {
                l = 0;
                phi[j] = 0;
                continue;
            }
            // symbols that both suffixes still have: n - max(j, p). The maximum is computed without a branch: j < p is a coin flip on random text, and
            // the JIT turned Math.Max into a conditional jump whose mispredictions cost about as much as the whole gain of the kernel (1o.1a).
            int limit = n - (j ^ ((j ^ p) & ((j - p) >> 31)));
            if (l < limit) l += CommonPrefix(text, j + l, p + l, limit - l);
            phi[j] = l;
            if (l >= 255)
            {
                e8++;
                if (l >= 65535) e16++;
            }
            if (l > 0) l--;
        }
        return (e8, e16);
    }

    /// <summary>
    /// Pass 3: <c>LCP[r] = PLCP[SA[r]]</c> for <c>r = 1..n</c> straight into the arrays of <paramref name="kind"/>; escaped entries go to the two parallel lists.
    /// </summary>
    private static CompactLcp WriteCompact(Buf<int> plcp, Buf<int> sa, int n, LcpKind kind, int escapes)
    {
        var overflowRank = new int[escapes];
        var overflowValue = new int[escapes];
        fixed (int* pRank = overflowRank)
        fixed (int* pValue = overflowValue)
        {
            var ranks = new Buf<int>(pRank, escapes, "overflowRank");
            var values = new Buf<int>(pValue, escapes, "overflowValue");
            switch (kind)
            {
                case LcpKind.U8:
                    {
                        byte[] lcp = GC.AllocateUninitializedArray<byte>(n + 1);
                        fixed (byte* p = lcp) WriteU8(plcp, sa, n, new Buf<byte>(p, n + 1, "lcp8"), ranks, values, escapes);
                        return CompactLcp.FromParts(kind, n + 1, lcp, null, null, overflowRank, overflowValue);
                    }
                case LcpKind.U16:
                    {
                        ushort[] lcp = GC.AllocateUninitializedArray<ushort>(n + 1);
                        fixed (ushort* p = lcp) WriteU16(plcp, sa, n, new Buf<ushort>(p, n + 1, "lcp16"), ranks, values, escapes);
                        return CompactLcp.FromParts(kind, n + 1, null, lcp, null, overflowRank, overflowValue);
                    }
                default:
                    {
                        int[] lcp = GC.AllocateUninitializedArray<int>(n + 1);
                        fixed (int* p = lcp) WriteU32(plcp, sa, n, new Buf<int>(p, n + 1, "lcp32"));
                        return CompactLcp.FromParts(kind, n + 1, null, null, lcp, overflowRank, overflowValue);
                    }
            }
        }
    }

    private static void WriteU32(Buf<int> plcp, Buf<int> sa, int n, Buf<int> dst)
    {
        dst[0] = 0;
        for (int r = 1; r <= n; r++) dst[r] = plcp[sa[r]];
    }

    private static void WriteU8(Buf<int> plcp, Buf<int> sa, int n, Buf<byte> dst, Buf<int> ranks, Buf<int> values, int escapes)
    {
        dst[0] = 0;
        int k = 0;
        for (int r = 1; r <= n; r++)
        {
            int v = plcp[sa[r]];
            if (v < byte.MaxValue)
            {
                dst[r] = (byte)v;
                continue;
            }
            dst[r] = byte.MaxValue;
            if ((uint)k >= (uint)escapes) ThrowEscapeMismatch(k, escapes);
            ranks[k] = r;
            values[k] = v;
            k++;
        }
        if (k != escapes) ThrowEscapeMismatch(k, escapes);
    }

    private static void WriteU16(Buf<int> plcp, Buf<int> sa, int n, Buf<ushort> dst, Buf<int> ranks, Buf<int> values, int escapes)
    {
        dst[0] = 0;
        int k = 0;
        for (int r = 1; r <= n; r++)
        {
            int v = plcp[sa[r]];
            if (v < ushort.MaxValue)
            {
                dst[r] = (ushort)v;
                continue;
            }
            dst[r] = ushort.MaxValue;
            if ((uint)k >= (uint)escapes) ThrowEscapeMismatch(k, escapes);
            ranks[k] = r;
            values[k] = v;
            k++;
        }
        if (k != escapes) ThrowEscapeMismatch(k, escapes);
    }

    private static void ThrowEscapeMismatch(int written, int announced)
        => throw new InvalidOperationException($"{written} escapes met, {announced} announced from PLCP: the suffix array is not a permutation");

    /// <summary>
    /// The length of the common prefix of the suffixes at <paramref name="a"/> and <paramref name="b"/>, at most <paramref name="max"/> (the caller guarantees that
    /// both have that many symbols left). The first <see cref="HeadSymbols"/> symbols are compared one by one with a branch on each: on random text almost every
    /// match ends there, and a branch lets the CPU speculate the carried length <c>l</c> and overlap the cache misses of the following iterations (a word compare with
    /// a trailing-zero count makes <c>l</c> depend on the loaded data and serialises them: measured 2.3 times slower on 10 Mbp DNA). Long matches (repeats) continue
    /// with 8-byte words; the tail is compared symbol by symbol and nothing at or beyond <c>a + max</c> / <c>b + max</c> is read.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static int CommonPrefix<TSym>(Buf<TSym> t, int a, int b, int max)
        where TSym : unmanaged, IEquatable<TSym>
    {
        int i = 0;
        int head = max < HeadSymbols ? max : HeadSymbols;
        while (i < head)
        {
            if (!t[a + i].Equals(t[b + i])) return i;
            i++;
        }
        int perWord = sizeof(ulong) / sizeof(TSym);
        while (i + perWord <= max)
        {
            ulong x = t.ReadUnaligned<ulong>(a + i) ^ t.ReadUnaligned<ulong>(b + i);
            if (x != 0)
            {
                int bits = BitConverter.IsLittleEndian ? BitOperations.TrailingZeroCount(x) : BitOperations.LeadingZeroCount(x);
                return i + (bits >> (sizeof(TSym) == 1 ? 3 : 4)); // bits / (8 * sizeof(TSym)) for the two widths
            }
            i += perWord;
        }
        while (i < max && t[a + i].Equals(t[b + i])) i++;
        return i;
    }

    /// <summary>Symbols compared one by one before the word loop takes over.</summary>
    private const int HeadSymbols = 8;
}
