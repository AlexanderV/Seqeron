using System.Runtime.CompilerServices;

namespace SuffixTree.Esa;

/// <summary>Storage kind of a <see cref="CompactLcp"/> (ADR-002 section 4; the container kinds of ADR-003).</summary>
internal enum LcpKind : byte
{
    /// <summary>One byte per entry; the value 255 is an escape, the true value (>= 255) is in the overflow list.</summary>
    U8,

    /// <summary>Two bytes per entry; the value 65,535 is an escape resolved through the overflow list.</summary>
    U16,

    /// <summary>Four bytes per entry, no overflow list.</summary>
    U32,
}

/// <summary>
/// The LCP array (<c>n + 1</c> entries, <c>LCP[0] = LCP[1] = 0</c>, <c>LCP[r] = lcp(SA[r-1], SA[r])</c>) in the smallest of three
/// representations: <c>U8</c> + overflow, <c>U16</c> + overflow, <c>U32</c>. An entry whose value reaches the escape (255 resp.
/// 65,535, so the escape value itself is an overflow entry too) is stored as the escape and its true value as a pair
/// <c>(rank, value)</c> in two parallel ascending arrays. Sequential consumers use a <see cref="Cursor"/> (amortised O(1)
/// per entry), random access binary-searches the list, which is only needed for the rare large values.
/// </summary>
internal sealed unsafe class CompactLcp
{
    private const int U8Escape = byte.MaxValue;
    private const int U16Escape = ushort.MaxValue;

    private readonly byte[]? _u8;
    private readonly ushort[]? _u16;
    private readonly int[]? _u32;
    private readonly int[] _overflowRank;
    private readonly int[] _overflowValue;

    // the mapped form (work package 1b.7f): the base array and the interleaved (rank, value) pairs of the overflow list in a SESA container; the managed arrays above are then null / empty
    private readonly byte* _mappedBase;
    private readonly uint* _mappedPairs;
    private readonly int _mappedPairCount;
    private readonly bool _mappedChecked;
    private int[]? _mappedRanks;
    private int[]? _mappedValues;

    private CompactLcp(LcpKind kind, int length, byte[]? u8, ushort[]? u16, int[]? u32, int[] overflowRank, int[] overflowValue)
    {
        Kind = kind;
        Length = length;
        _u8 = u8;
        _u16 = u16;
        _u32 = u32;
        _overflowRank = overflowRank;
        _overflowValue = overflowValue;
    }

    private CompactLcp(LcpKind kind, int length, byte* baseArray, uint* pairs, int pairCount, bool isChecked)
    {
        Kind = kind;
        Length = length;
        _overflowRank = [];
        _overflowValue = [];
        _mappedBase = baseArray;
        _mappedPairs = pairs;
        _mappedPairCount = pairCount;
        _mappedChecked = isChecked;
        IsMapped = true;
    }

    /// <summary>
    /// The LCP array of a SESA container, read in place (work package 1b.7f): <paramref name="baseArray"/> holds <paramref name="length"/> entries of the width of <paramref name="kind"/>, <paramref name="pairs"/> the
    /// <paramref name="pairCount"/> (rank, value) pairs of the overflow list as consecutive 32-bit words. The memory belongs to the container and must outlive this object. The entries are untrusted: a read that the
    /// arrays cannot serve is an <see cref="InvalidDataException"/>, and with <paramref name="isChecked"/> every value is also checked against <c>0..n+1</c>.
    /// </summary>
    internal static CompactLcp FromMapped(LcpKind kind, int length, byte* baseArray, uint* pairs, int pairCount, bool isChecked)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(length, 1);
        ArgumentOutOfRangeException.ThrowIfNegative(pairCount);
        if (kind == LcpKind.U32 && pairCount != 0) throw new InvalidDataException("an LCP array of 32 bits per entry has no overflow list");
        return new CompactLcp(kind, length, baseArray, pairs, pairCount, isChecked);
    }

    /// <summary>True when the array is read in place from a container.</summary>
    internal bool IsMapped { get; }

    /// <summary>The mapped base array (null for a managed one).</summary>
    internal byte* MappedBase => _mappedBase;

    /// <summary>The representation.</summary>
    internal LcpKind Kind { get; }

    /// <summary>Number of entries, <c>n + 1</c>.</summary>
    internal int Length { get; }

    /// <summary>Number of escaped entries (always 0 for U32).</summary>
    internal int OverflowCount => IsMapped ? _mappedPairCount : _overflowRank.Length;

    /// <summary>Bytes of the base array plus 8 bytes per overflow entry.</summary>
    internal long PayloadBytes => ((long)Length * Kind switch { LcpKind.U8 => 1, LcpKind.U16 => 2, _ => 4 }) + (8L * OverflowCount);

    /// <summary>The corruption of a mapped file is an <see cref="InvalidDataException"/>, that of a managed table the <see cref="InvalidOperationException"/> it always was.</summary>
    private Exception Corrupt(string message) => IsMapped ? new InvalidDataException(message + ": the file is damaged") : new InvalidOperationException(message);

    /// <summary>The escape value of the base array, or -1 when there is none (U32).</summary>
    private int EscapeValue => Kind switch { LcpKind.U8 => U8Escape, LcpKind.U16 => U16Escape, _ => -1 };

    /// <summary>The base array of kind <see cref="LcpKind.U8"/> (null for the other kinds). Read access for the pointer readers of the navigation pack (work package 1b.4); never written.</summary>
    internal byte[]? BaseU8 => _u8;

    /// <summary>The base array of kind <see cref="LcpKind.U16"/> (null for the other kinds); see <see cref="BaseU8"/>.</summary>
    internal ushort[]? BaseU16 => _u16;

    /// <summary>The base array of kind <see cref="LcpKind.U32"/> (null for the other kinds); see <see cref="BaseU8"/>.</summary>
    internal int[]? BaseU32 => _u32;

    /// <summary>The ranks of the escaped entries, strictly ascending (empty for U32); see <see cref="BaseU8"/>.</summary>
    internal int[] OverflowRanks => IsMapped ? MappedArrays().Ranks : _overflowRank;

    /// <summary>The true values of the escaped entries, parallel to <see cref="OverflowRanks"/>; see <see cref="BaseU8"/>.</summary>
    internal int[] OverflowValues => IsMapped ? MappedArrays().Values : _overflowValue;

    /// <summary>The LCP value at <paramref name="rank"/> (random access; binary search for escaped entries).</summary>
    internal int this[int rank]
    {
        get
        {
            int b = BaseValue(rank);
            return b != EscapeValue ? b : EscapedValue(rank);
        }
    }

    /// <summary>The value of an escaped entry (the slow arm of the indexer, out of line so that the arm that most reads take stays small enough to be inlined into the pattern descent).</summary>
    [MethodImpl(MethodImplOptions.NoInlining)]
    private int EscapedValue(int rank)
    {
        int k = FindOverflow(rank);
        if (k < 0) throw Corrupt($"overflow entry for rank {rank} is missing");
        return OverflowValueAt(k);
    }

    /// <summary>The index in the overflow list of <paramref name="rank"/>, or a negative number.</summary>
    private int FindOverflow(int rank)
    {
        if (!IsMapped) return Array.BinarySearch(_overflowRank, rank);
        int lo = 0, hi = _mappedPairCount - 1;
        while (lo <= hi)
        {
            int mid = (int)(((uint)lo + (uint)hi) >> 1);
            uint v = _mappedPairs[2 * (long)mid];
            if (v == (uint)rank) return mid;
            if (v < (uint)rank) lo = mid + 1;
            else hi = mid - 1;
        }
        return -1;
    }

    private int OverflowRankAt(int index) => IsMapped ? MappedWord(2 * (long)index) : _overflowRank[index];

    private int OverflowValueAt(int index) => IsMapped ? MappedWord((2 * (long)index) + 1) : _overflowValue[index];

    /// <summary>A word of the mapped overflow list as a non-negative <see cref="int"/> (a word above <c>int.MaxValue</c> cannot be a rank or an LCP value).</summary>
    private int MappedWord(long word)
    {
        uint v = _mappedPairs[word];
        if (v > int.MaxValue) throw Corrupt($"the overflow list holds {v}");
        return (int)v;
    }

    /// <summary>The overflow list as two arrays (the navigation pack's readers want arrays; it is built once, on the first suffix link).</summary>
    private (int[] Ranks, int[] Values) MappedArrays()
    {
        if (_mappedRanks is null || _mappedValues is null)
        {
            var ranks = new int[_mappedPairCount];
            var values = new int[_mappedPairCount];
            for (int i = 0; i < ranks.Length; i++)
            {
                ranks[i] = OverflowRankAt(i);
                values[i] = OverflowValueAt(i);
            }
            _mappedValues = values;
            _mappedRanks = ranks;
        }
        return (_mappedRanks, _mappedValues);
    }

    /// <summary>A forward cursor starting at the first overflow entry.</summary>
    internal Cursor CreateCursor() => new(this);

    /// <summary>All entries as a dense array (tests, diagnostics).</summary>
    internal int[] ToArray()
    {
        var result = new int[Length];
        var cursor = CreateCursor();
        for (int i = 0; i < result.Length; i++) result[i] = cursor.Get(i);
        return result;
    }

    /// <summary>
    /// The selection rule of ADR-002 section 4, from the text length <paramref name="n"/> and the numbers of values
    /// <c>&gt;= 255</c> (<paramref name="escapes8"/>) and <c>&gt;= 65,535</c> (<paramref name="escapes16"/>):
    /// <c>U8</c> + overflow iff <c>8 * escapes8 &lt;= n / 2</c> (the overflow list is at most half of the base array), written
    /// <c>16 * escapes8 &lt;= n</c>; else <c>U16</c> + overflow iff <c>8 * escapes16 &lt;= n</c> (the same condition one level up:
    /// overflow list at most half of the 2n-byte base array; the ADR states this tier only as "else U16(+overflow)"); else <c>U32</c>.
    /// </summary>
    internal static LcpKind SelectKind(int n, long escapes8, long escapes16)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(n);
        ArgumentOutOfRangeException.ThrowIfNegative(escapes8);
        ArgumentOutOfRangeException.ThrowIfNegative(escapes16);
        if (escapes16 > escapes8)
            throw new ArgumentException("every value >= 65,535 is also >= 255, so escapes16 cannot exceed escapes8");
        if (16 * escapes8 <= n) return LcpKind.U8;
        if (8 * escapes16 <= n) return LcpKind.U16;
        return LcpKind.U32;
    }

    /// <summary>Stores a dense LCP array (n + 1 entries) in the representation chosen by <see cref="SelectKind"/>.</summary>
    internal static CompactLcp FromDense(ReadOnlySpan<int> lcp)
    {
        if (lcp.IsEmpty) throw new ArgumentException("an LCP array has n + 1 >= 1 entries", nameof(lcp));
        long e8 = 0, e16 = 0;
        foreach (int v in lcp)
        {
            if (v >= U8Escape) e8++;
            if (v >= U16Escape) e16++;
        }
        return FromDense(lcp, SelectKind(lcp.Length - 1, e8, e16));
    }

    /// <summary>Stores a dense LCP array (n + 1 entries) in a given representation.</summary>
    internal static CompactLcp FromDense(ReadOnlySpan<int> lcp, LcpKind kind)
    {
        if (lcp.IsEmpty) throw new ArgumentException("an LCP array has n + 1 >= 1 entries", nameof(lcp));
        int escape = kind switch { LcpKind.U8 => U8Escape, LcpKind.U16 => U16Escape, _ => int.MaxValue };
        int escapes = 0;
        foreach (int v in lcp)
        {
            ArgumentOutOfRangeException.ThrowIfNegative(v);
            if (kind != LcpKind.U32 && v >= escape) escapes++;
        }
        var builder = new Builder(kind, lcp.Length, escapes);
        for (int r = 0; r < lcp.Length; r++) builder.Set(r, lcp[r]);
        return builder.Complete();
    }

    /// <summary>
    /// Wraps the arrays that an optimized kernel has written (no copy). Checks the shape the readers rely on: exactly the base array of <paramref name="kind"/>
    /// with <paramref name="length"/> entries, an overflow list of equal lengths with strictly ascending ranks inside the array and values at or above the
    /// escape, no overflow list for U32, and that every escape slot of the base array is listed (the base array is scanned once for the escape value).
    /// </summary>
    internal static CompactLcp FromParts(LcpKind kind, int length, byte[]? u8, ushort[]? u16, int[]? u32, int[] overflowRank, int[] overflowValue)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(length, 1);
        ArgumentNullException.ThrowIfNull(overflowRank);
        ArgumentNullException.ThrowIfNull(overflowValue);
        bool shape = kind switch
        {
            LcpKind.U8 => u8 is { } a && a.Length == length && u16 is null && u32 is null,
            LcpKind.U16 => u16 is { } b && b.Length == length && u8 is null && u32 is null,
            _ => u32 is { } c && c.Length == length && u8 is null && u16 is null,
        };
        if (!shape) throw new ArgumentException("exactly one base array of the stated kind and length is required");
        if (overflowRank.Length != overflowValue.Length) throw new ArgumentException("overflow ranks and values differ in length");
        if (kind == LcpKind.U32 && overflowRank.Length != 0) throw new ArgumentException("U32 has no overflow list");
        int escape = kind switch { LcpKind.U8 => U8Escape, LcpKind.U16 => U16Escape, _ => int.MaxValue };
        for (int i = 0; i < overflowRank.Length; i++)
        {
            if ((uint)overflowRank[i] >= (uint)length || (i > 0 && overflowRank[i - 1] >= overflowRank[i]))
                throw new ArgumentException("overflow ranks must be strictly ascending and inside the array");
            if (overflowValue[i] < escape) throw new ArgumentException("an overflow value must be at least the escape value");
        }
        return new CompactLcp(kind, length, u8, u16, u32, overflowRank, overflowValue);
    }

    [MethodImpl(MethodImplOptions.NoInlining)]
    private void ThrowRankOutside(int rank)
    {
        if (IsMapped) throw new InvalidDataException($"the LCP array is read at rank {rank}, but it has {Length} entries: the file is damaged");
        throw new ArgumentOutOfRangeException(nameof(rank));
    }

    private int BaseValue(int rank)
    {
        if ((uint)rank >= (uint)Length) ThrowRankOutside(rank);
        if (IsMapped) return MappedValue(rank);
        return Kind switch
        {
            LcpKind.U8 => _u8![rank],
            LcpKind.U16 => _u16![rank],
            _ => _u32![rank],
        };
    }

    private int MappedValue(int rank)
    {
        switch (Kind)
        {
            case LcpKind.U8:
                return _mappedBase[rank];
            case LcpKind.U16:
                return ((ushort*)_mappedBase)[rank];
            default:
                int v = ((int*)_mappedBase)[rank];
                if (_mappedChecked && (uint)v > (uint)Length) throw Corrupt($"the LCP value {(uint)v} at rank {rank} is larger than the text");
                return v;
        }
    }

    /// <summary>Sequential reader that resolves escaped entries with a moving index into the overflow list.</summary>
    internal struct Cursor
    {
        private readonly CompactLcp _owner;
        private int _next;

        internal Cursor(CompactLcp owner)
        {
            _owner = owner;
            _next = 0;
        }

        /// <summary>
        /// The LCP value at <paramref name="rank"/>. Amortised O(1) when the ranks are non-decreasing; a smaller rank than the
        /// previous escaped one repositions by binary search.
        /// </summary>
        internal int Get(int rank)
        {
            int b = _owner.BaseValue(rank);
            if (b != _owner.EscapeValue) return b;

            int count = _owner.OverflowCount;
            if (_next > 0 && _owner.OverflowRankAt(_next - 1) >= rank)
            {
                _next = _owner.FindOverflow(rank);
                if (_next < 0) throw _owner.Corrupt($"overflow entry for rank {rank} is missing");
            }
            while (_next < count && _owner.OverflowRankAt(_next) < rank) _next++;
            if (_next == count || _owner.OverflowRankAt(_next) != rank)
                throw _owner.Corrupt($"overflow entry for rank {rank} is missing");
            return _owner.OverflowValueAt(_next);
        }
    }

    /// <summary>
    /// Writes entries in ascending rank order straight into the chosen representation, so that no dense 4-byte-per-entry array is
    /// needed (ADR-002 section 4, pass 3). The number of escapes must be known in advance.
    /// </summary>
    internal sealed class Builder
    {
        private readonly LcpKind _kind;
        private readonly int _length;
        private readonly byte[]? _u8;
        private readonly ushort[]? _u16;
        private readonly int[]? _u32;
        private readonly int[] _overflowRank;
        private readonly int[] _overflowValue;
        private int _overflowFilled;

        internal Builder(LcpKind kind, int length, int overflowCount)
        {
            ArgumentOutOfRangeException.ThrowIfLessThan(length, 1);
            ArgumentOutOfRangeException.ThrowIfNegative(overflowCount);
            if (kind == LcpKind.U32 && overflowCount != 0)
                throw new ArgumentException("U32 has no overflow list", nameof(overflowCount));
            _kind = kind;
            _length = length;
            _u8 = kind == LcpKind.U8 ? new byte[length] : null;
            _u16 = kind == LcpKind.U16 ? new ushort[length] : null;
            _u32 = kind == LcpKind.U32 ? new int[length] : null;
            _overflowRank = new int[overflowCount];
            _overflowValue = new int[overflowCount];
        }

        /// <summary>Stores <paramref name="value"/> at <paramref name="rank"/>; escaped entries must arrive in ascending rank order.</summary>
        internal void Set(int rank, int value)
        {
            ArgumentOutOfRangeException.ThrowIfNegative(value);
            if ((uint)rank >= (uint)_length) throw new ArgumentOutOfRangeException(nameof(rank));
            switch (_kind)
            {
                case LcpKind.U8:
                    if (value >= U8Escape)
                    {
                        _u8![rank] = U8Escape;
                        AddOverflow(rank, value);
                    }
                    else
                    {
                        _u8![rank] = (byte)value;
                    }
                    break;
                case LcpKind.U16:
                    if (value >= U16Escape)
                    {
                        _u16![rank] = U16Escape;
                        AddOverflow(rank, value);
                    }
                    else
                    {
                        _u16![rank] = (ushort)value;
                    }
                    break;
                default:
                    _u32![rank] = value;
                    break;
            }
        }

        /// <summary>Finishes the build; fails if the announced number of escapes was not written.</summary>
        internal CompactLcp Complete()
        {
            if (_overflowFilled != _overflowRank.Length)
                throw new InvalidOperationException($"{_overflowFilled} escapes written, {_overflowRank.Length} announced");
            return new CompactLcp(_kind, _length, _u8, _u16, _u32, _overflowRank, _overflowValue);
        }

        private void AddOverflow(int rank, int value)
        {
            if (_overflowFilled == _overflowRank.Length)
                throw new InvalidOperationException("more escapes than announced");
            if (_overflowFilled > 0 && _overflowRank[_overflowFilled - 1] >= rank)
                throw new InvalidOperationException("escaped entries must be written in ascending rank order");
            _overflowRank[_overflowFilled] = rank;
            _overflowValue[_overflowFilled] = value;
            _overflowFilled++;
        }
    }
}
