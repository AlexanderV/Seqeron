using System.Runtime.CompilerServices;

// The delta forms of the child table (ADR-005 section 5, the owner's candidate "a delta-from-index form with escapes for CHILD"): a cell stores the value of the field as a small signed distance from its own index
// and the few fields that are too far (the top of the tree) go to a sorted list of (index, value) pairs found by binary search, as the LCP escapes do (CompactLcp).
//   d16  two bytes: the kind of the field in bits 14-15 (as CHILD codec 2), the distance in bits 0-13 (-8191..8191; -8192 = "see the list")        2 B per cell, the kind needs no LCP read
//   d8   one byte: the distance (-126..127; -128 = empty cell, -127 = "see the list"), the kind told by comparing LCP values as in CHILD codec 1   1 B per cell
// The readers return what the 4-byte layouts return (the value, the kind above bit 40 for d16), so the descent code is the one of the other variants.

internal sealed unsafe class DeltaCells : IDisposable
{
    internal byte* Cells;
    internal uint* Pairs;
    internal long PairCount;
    internal long Bytes;
    internal readonly int Width;
    internal readonly bool Hashed;
    internal ulong* Table;     // open addressing, linear probing: (key + 1) << 32 | value, 0 = free
    internal int TableShift;

    internal DeltaCells(int width, uint* tagged4, long n, bool packedRaw, bool hashed = false)
    {
        Width = width;
        Hashed = hashed;
        long cells = n + 1;
        Cells = Mem.Alloc(cells * width);
        var pairs = new List<(uint, uint)>();
        for (long k = 0; k < cells; k++)
        {
            uint cell = tagged4[k];
            uint kind = cell >> 30;
            long value = cell & ((1u << 30) - 1);
            if (width == 2)
            {
                ushort x;
                if (kind == 0) x = 0;
                else
                {
                    long d = value - k;
                    if (d >= -8191 && d <= 8191) x = (ushort)((kind << 14) | ((uint)d & 0x3FFF));
                    else { x = (ushort)((kind << 14) | 0x2000); pairs.Add(((uint)k, (uint)value)); }
                }
                ((ushort*)Cells)[k] = x;
            }
            else
            {
                // packed: a next link to the following rank is not stored; empty = -128
                long raw = kind == 0 || (kind == 1 && value == k + 1) ? 0 : value;
                sbyte b;
                if (raw == 0) b = -128;
                else
                {
                    long d = raw - k;
                    if (d >= -126 && d <= 127) b = (sbyte)d;
                    else { b = -127; pairs.Add(((uint)k, (uint)raw)); }
                }
                ((sbyte*)Cells)[k] = b;
            }
        }
        PairCount = pairs.Count;
        Pairs = (uint*)Mem.Alloc((PairCount * 8) + 8);
        for (int i = 0; i < pairs.Count; i++) { Pairs[2 * i] = pairs[i].Item1; Pairs[(2 * i) + 1] = pairs[i].Item2; }
        Bytes = (cells * width) + (PairCount * 8);
        if (hashed)
        {
            int bits = 4;
            while ((1L << bits) < PairCount * 2) bits++;
            TableShift = 64 - bits;
            long slots = 1L << bits;
            Table = (ulong*)Mem.Alloc(slots * 8);
            new Span<ulong>(Table, (int)slots).Clear();
            for (int i = 0; i < PairCount; i++)
            {
                ulong h = ((ulong)Pairs[2 * i] * 0x9E3779B97F4A7C15UL) >> TableShift;
                while (Table[h] != 0) h = (h + 1) & (ulong)(slots - 1);
                Table[h] = ((ulong)(Pairs[2 * i] + 1UL) << 32) | Pairs[(2 * i) + 1];
            }
            Bytes = (cells * width) + (slots * 8);
        }
        for (long o = 0; o < cells * width; o += 4096) _ = Cells[o];
    }

    internal long Escaped(long k)
    {
        if (Hashed)
        {
            ulong mask = (1UL << (64 - TableShift)) - 1;
            ulong h = ((ulong)k * 0x9E3779B97F4A7C15UL) >> TableShift;
            while (true)
            {
                ulong e = Table[h];
                if ((e >> 32) == (ulong)k + 1) return (long)(e & 0xFFFFFFFFUL);
                if (e == 0) throw new InvalidOperationException("no escape for cell " + k);
                h = (h + 1) & mask;
            }
        }
        long lo = 0, hi = PairCount - 1;
        while (lo <= hi)
        {
            long mid = (lo + hi) >> 1;
            uint r = Pairs[2 * mid];
            if (r == k) return Pairs[(2 * mid) + 1];
            if (r < k) lo = mid + 1; else hi = mid - 1;
        }
        throw new InvalidOperationException("no escape for cell " + k);
    }

    public void Dispose()
    {
        Mem.Free(Cells);
        Mem.Free(Pairs);
        if (Table != null) Mem.Free(Table);
    }
}

/// <summary>d16: Get(i) = (kind &lt;&lt; 40) | value, so that <see cref="TaggedCellsR{TR}"/> with 40 value bits reads it.</summary>
internal readonly unsafe struct RD16(DeltaCells d) : IReader
{
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public ulong Get(long i)
    {
        ushort x = ((ushort*)d.Cells)[i];
        int delta = ((short)(x << 2)) >> 2;
        ulong kind = (ulong)(x >> 14) << 40;
        if (delta == -8192) return kind | (ulong)d.Escaped(i);
        return kind | (ulong)(i + delta);
    }
}

/// <summary>d8: Get(i) = the packed raw value (0 for an empty cell), so that <see cref="PackedCellsR{TR, TL}"/> reads it.</summary>
internal readonly unsafe struct RD8(DeltaCells d) : IReader
{
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public ulong Get(long i)
    {
        sbyte b = ((sbyte*)d.Cells)[i];
        if (b == -128) return 0;
        if (b == -127) return (ulong)d.Escaped(i);
        return (ulong)(i + b);
    }
}
