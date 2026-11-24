using System.Numerics;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;

// Stored-position layouts of ADR-005 section 5: how one array of n entries (values below V) is stored and read back.
//   w4   fixed 32-bit unsigned entries (the SESA POS_WIDTH 4 of today, n + 1 <= 2^32 - 1)
//   w5   fixed 40-bit entries, byte-packed (POS_WIDTH 5: unaligned 8-byte load and a mask; 6.25 % of the entries straddle a 64-byte line)
//   w8   fixed 64-bit entries (POS_WIDTH 8)
//   hyb  hybrid with escapes: 32-bit slots, the value 0xFFFFFFFF = "long", the long values (40-bit) in a second array reached through a
//        bitmap + cumulative count per 64 slots (the nearest equivalent for a rank-indexed array of the zoned layout of the v6 tree)
//   p34  34-bit entries packed 15 to a 64-byte line (510 of 512 bits; one cache line per read, never a straddle; covers values below 2^34)
// All readers are readonly structs over raw pointers, so that a generic loop is specialised by the JIT per layout (as MappedSearch is).

internal interface IReader
{
    ulong Get(long i);

    /// <summary>The sum of the entries [i, i + m): a locate reads the slice of the suffix array sequentially.</summary>
    ulong SumRange(long i, long m)
    {
        ulong s = 0;
        for (long k = 0; k < m; k++) s += Get(i + k);
        return s;
    }
}

internal static unsafe class Mem
{
    internal static byte* Alloc(long bytes)
    {
        byte* p = (byte*)NativeMemory.AlignedAlloc((nuint)(bytes + 64), 4096);
        if (p == null) throw new OutOfMemoryException();
        return p;
    }

    internal static void Free(void* p) => NativeMemory.AlignedFree(p);
}

internal static class Rng
{
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    internal static ulong Mix(ulong x)
    {
        x += 0x9E3779B97F4A7C15UL;
        x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9UL;
        x = (x ^ (x >> 27)) * 0x94D049BB133111EBUL;
        return x ^ (x >> 31);
    }

    /// <summary>The canonical value of entry i: uniform in [0, v).</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    internal static ulong Value(ulong seed, long i, ulong v) => Math.BigMul(Mix(seed ^ (ulong)i), v, out _);

    /// <summary>The next index of a chain: the low 32 bits of the value, scaled to [0, n).</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    internal static long Next(ulong value, long n) => (long)(((value & 0xFFFFFFFFUL) * (ulong)n) >> 32);
}

internal readonly unsafe struct R4(byte* p) : IReader
{
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public ulong Get(long i) => ((uint*)p)[i];

    public ulong SumRange(long i, long m)
    {
        ulong s = 0;
        uint* q = (uint*)p + i;
        for (long k = 0; k < m; k++) s += q[k];
        return s;
    }
}

internal readonly unsafe struct R5(byte* p) : IReader
{
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public ulong Get(long i) => Unsafe.ReadUnaligned<ulong>(p + (5 * i)) & 0xFFFFFFFFFFUL;

    public ulong SumRange(long i, long m)
    {
        ulong s = 0;
        byte* q = p + (5 * i);
        for (long k = 0; k < m; k++, q += 5) s += Unsafe.ReadUnaligned<ulong>(q) & 0xFFFFFFFFFFUL;
        return s;
    }
}

internal readonly unsafe struct R8(byte* p) : IReader
{
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public ulong Get(long i) => ((ulong*)p)[i];

    public ulong SumRange(long i, long m)
    {
        ulong s = 0;
        ulong* q = (ulong*)p + i;
        for (long k = 0; k < m; k++) s += q[k];
        return s;
    }
}

/// <summary>34 bits per entry, 15 entries per 64-byte line.</summary>
internal readonly unsafe struct R34(byte* p) : IReader
{
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public ulong Get(long i)
    {
        long line = i / 15;
        int k = (int)(i - (line * 15));
        int bit = k * 34;
        int b = Math.Min(bit >> 3, 56);
        ulong x = Unsafe.ReadUnaligned<ulong>(p + (line * 64) + b);
        return (x >> (bit - (b << 3))) & 0x3FFFFFFFFUL;
    }

    public ulong SumRange(long i, long m)
    {
        ulong s = 0;
        for (long k = 0; k < m; k++) s += Get(i + k);
        return s;
    }
}

/// <summary>The hybrid: 32-bit slots, an escape value, a 16-byte record (bitmap, cumulative count) per 64 slots and the long values at 40 bits.</summary>
internal readonly unsafe struct RHyb(byte* slots, byte* blocks, byte* longs) : IReader
{
    internal const uint Escape = 0xFFFFFFFFu;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public ulong Get(long i)
    {
        uint s = ((uint*)slots)[i];
        if (s != Escape) return s;
        return GetLong(i);
    }

    [MethodImpl(MethodImplOptions.NoInlining)]
    private ulong GetLong(long i)
    {
        ulong* rec = (ulong*)(blocks + ((i >> 6) * 16));
        long j = (long)rec[1] + BitOperations.PopCount(rec[0] & ((1UL << (int)(i & 63)) - 1));
        return Unsafe.ReadUnaligned<ulong>(longs + (5 * j)) & 0xFFFFFFFFFFUL;
    }

    /// <summary>The sum of the entries [i, i + m) by a sequential decode (one rank lookup, then the long values in order).</summary>
    public ulong SumRange(long i, long m)
    {
        ulong sum = 0;
        long j = -1;
        for (long x = i; x < i + m; x++)
        {
            uint s = ((uint*)slots)[x];
            if (s != Escape) { sum += s; continue; }
            if (j < 0)
            {
                ulong* rec = (ulong*)(blocks + ((x >> 6) * 16));
                j = (long)rec[1] + BitOperations.PopCount(rec[0] & ((1UL << (int)(x & 63)) - 1));
            }
            sum += Unsafe.ReadUnaligned<ulong>(longs + (5 * j)) & 0xFFFFFFFFFFUL;
            j++;
        }
        return sum;
    }
}

/// <summary>An allocated array in one of the layouts, with the native blocks it owns.</summary>
internal sealed unsafe class Packed : IDisposable
{
    internal string Layout { get; }
    internal long N { get; }
    internal long Bytes { get; private set; }
    internal long Escapes { get; private set; }
    internal byte* Main;
    internal byte* Blocks;
    internal byte* Longs;

    private Packed(string layout, long n) { Layout = layout; N = n; }

    internal static long MainBytes(string layout, long n) => layout switch
    {
        "w4" or "hyb" => 4 * n,
        "w5" => 5 * n,
        "w8" => 8 * n,
        "p34" => ((n + 14) / 15) * 64,
        _ => throw new ArgumentException(layout),
    };

    /// <summary>The array whose entry i holds <c>Rng.Value(seed, i, v)</c> (w4: the low 32 bits; w5: the low 40 bits; the others the whole value; hyb: escapes above 2^32 - 2).</summary>
    internal static Packed Create(string layout, long n, ulong v, ulong seed, int threads = 4)
        => CreateFrom(layout, n, i => Rng.Value(seed, i, v), RHyb.Escape, threads);

    /// <summary>The array whose entry i holds <c>value(i)</c>; for the hybrid the entries at or above <paramref name="escapeAt"/> (at most 2^32 - 1) are the long ones.</summary>
    internal static Packed CreateFrom(string layout, long n, Func<long, ulong> value, ulong escapeAt, int threads = 4)
    {
        var a = new Packed(layout, n);
        a.Main = Mem.Alloc(MainBytes(layout, n));
        a.Bytes = MainBytes(layout, n);
        Unsafe.InitBlock(a.Main + a.Bytes, 0, 64);
        if (layout == "hyb") a.BuildHybrid(value, escapeAt, threads);
        else ParallelFor(n, threads, (lo, hi) => a.FillSimple(lo, hi, value));
        return a;
    }

    private static void ParallelFor(long n, int threads, Action<long, long> body)
    {
        long chunk = (n + threads - 1) / threads;
        chunk = (chunk + 959) / 960 * 960; // a multiple of 64 and of 15 (every chunk starts at a line boundary of every layout)
        var ts = new List<Thread>();
        for (long lo = 0; lo < n; lo += chunk)
        {
            long l = lo, h = Math.Min(n, lo + chunk);
            var t = new Thread(() => body(l, h));
            ts.Add(t);
            t.Start();
        }
        foreach (Thread t in ts) t.Join();
    }

    private void FillSimple(long lo, long hi, Func<long, ulong> value)
    {
        switch (Layout)
        {
            case "w4": for (long i = lo; i < hi; i++) ((uint*)Main)[i] = (uint)value(i); break;
            case "w5": for (long i = lo; i < hi; i++) { ulong x = value(i); byte* q = Main + (5 * i); *(uint*)q = (uint)x; q[4] = (byte)(x >> 32); } break;
            case "w8": for (long i = lo; i < hi; i++) ((ulong*)Main)[i] = value(i); break;
            case "p34":
                for (long line = lo / 15; line * 15 < hi; line++)
                {
                    ulong* w = (ulong*)(Main + (line * 64));
                    for (int q = 0; q < 8; q++) w[q] = 0;
                    for (int k = 0; k < 15 && (line * 15) + k < N; k++)
                    {
                        ulong x = value((line * 15) + k) & 0x3FFFFFFFFUL;
                        int bit = k * 34, wi = bit >> 6, sh = bit & 63;
                        w[wi] |= x << sh;
                        if (sh + 34 > 64) w[wi + 1] |= x >> (64 - sh);
                    }
                }
                break;
            default: throw new InvalidOperationException();
        }
    }

    private void BuildHybrid(Func<long, ulong> value, ulong escapeAt, int threads)
    {
        long blocks = (N + 63) / 64;
        Blocks = Mem.Alloc(blocks * 16);
        var counts = new long[blocks];
        ParallelFor(N, threads, (lo, hi) =>
        {
            for (long b = lo / 64; b * 64 < hi; b++)
            {
                ulong bitmap = 0;
                long end = Math.Min(N, (b + 1) * 64);
                for (long i = b * 64; i < end; i++)
                {
                    ulong x = value(i);
                    if (x >= escapeAt) { bitmap |= 1UL << (int)(i & 63); ((uint*)Main)[i] = RHyb.Escape; }
                    else ((uint*)Main)[i] = (uint)x;
                }
                ((ulong*)Blocks)[2 * b] = bitmap;
                counts[b] = BitOperations.PopCount(bitmap);
            }
        });
        long total = 0;
        for (long b = 0; b < blocks; b++) { ((ulong*)Blocks)[(2 * b) + 1] = (ulong)total; total += counts[b]; }
        Escapes = total;
        Longs = Mem.Alloc((total * 5) + 8);
        ParallelFor(N, threads, (lo, hi) =>
        {
            for (long b = lo / 64; b * 64 < hi; b++)
            {
                long j = (long)((ulong*)Blocks)[(2 * b) + 1];
                long end = Math.Min(N, (b + 1) * 64);
                for (long i = b * 64; i < end; i++)
                {
                    if (((uint*)Main)[i] != RHyb.Escape) continue;
                    ulong x = value(i);
                    byte* q = Longs + (5 * j++);
                    *(uint*)q = (uint)x;
                    q[4] = (byte)(x >> 32);
                }
            }
        });
        Bytes = (4 * N) + (blocks * 16) + (total * 5);
    }

    /// <summary>Touches every page once (a sequential pass), so that the first random reads do not pay the page fault of anonymous memory.</summary>
    internal ulong Touch()
    {
        ulong s = 0;
        s += TouchRange(Main, MainBytes(Layout, N));
        if (Blocks != null) s += TouchRange(Blocks, ((N + 63) / 64) * 16);
        if (Longs != null) s += TouchRange(Longs, (Escapes * 5) + 8);
        return s;
    }

    private static ulong TouchRange(byte* p, long len)
    {
        ulong s = 0;
        for (long o = 0; o < len; o += 4096) s += p[o];
        return s;
    }

    public void Dispose()
    {
        if (Main != null) Mem.Free(Main);
        if (Blocks != null) Mem.Free(Blocks);
        if (Longs != null) Mem.Free(Longs);
        Main = Blocks = Longs = null;
    }
}
