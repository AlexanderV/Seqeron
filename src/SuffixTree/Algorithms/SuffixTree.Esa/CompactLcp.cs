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
internal sealed class CompactLcp
{
    private const int U8Escape = byte.MaxValue;
    private const int U16Escape = ushort.MaxValue;

    private readonly byte[]? _u8;
    private readonly ushort[]? _u16;
    private readonly int[]? _u32;
    private readonly int[] _overflowRank;
    private readonly int[] _overflowValue;

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

    /// <summary>The representation.</summary>
    internal LcpKind Kind { get; }

    /// <summary>Number of entries, <c>n + 1</c>.</summary>
    internal int Length { get; }

    /// <summary>Number of escaped entries (always 0 for U32).</summary>
    internal int OverflowCount => _overflowRank.Length;

    /// <summary>Bytes of the base array plus 8 bytes per overflow entry.</summary>
    internal long PayloadBytes => ((long)Length * Kind switch { LcpKind.U8 => 1, LcpKind.U16 => 2, _ => 4 }) + (8L * OverflowCount);

    /// <summary>The escape value of the base array, or -1 when there is none (U32).</summary>
    private int EscapeValue => Kind switch { LcpKind.U8 => U8Escape, LcpKind.U16 => U16Escape, _ => -1 };

    /// <summary>The base array of kind <see cref="LcpKind.U8"/> (null for the other kinds). Read access for the pointer readers of the navigation pack (work package 1b.4); never written.</summary>
    internal byte[]? BaseU8 => _u8;

    /// <summary>The base array of kind <see cref="LcpKind.U16"/> (null for the other kinds); see <see cref="BaseU8"/>.</summary>
    internal ushort[]? BaseU16 => _u16;

    /// <summary>The base array of kind <see cref="LcpKind.U32"/> (null for the other kinds); see <see cref="BaseU8"/>.</summary>
    internal int[]? BaseU32 => _u32;

    /// <summary>The ranks of the escaped entries, strictly ascending (empty for U32); see <see cref="BaseU8"/>.</summary>
    internal int[] OverflowRanks => _overflowRank;

    /// <summary>The true values of the escaped entries, parallel to <see cref="OverflowRanks"/>; see <see cref="BaseU8"/>.</summary>
    internal int[] OverflowValues => _overflowValue;

    /// <summary>The LCP value at <paramref name="rank"/> (random access; binary search for escaped entries).</summary>
    internal int this[int rank]
    {
        get
        {
            int b = BaseValue(rank);
            if (b != EscapeValue) return b;
            int k = Array.BinarySearch(_overflowRank, rank);
            if (k < 0) throw new InvalidOperationException($"overflow entry for rank {rank} is missing");
            return _overflowValue[k];
        }
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

    private int BaseValue(int rank)
    {
        if ((uint)rank >= (uint)Length) throw new ArgumentOutOfRangeException(nameof(rank));
        return Kind switch
        {
            LcpKind.U8 => _u8![rank],
            LcpKind.U16 => _u16![rank],
            _ => _u32![rank],
        };
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

            int[] ranks = _owner._overflowRank;
            if (_next > 0 && ranks[_next - 1] >= rank)
            {
                _next = Array.BinarySearch(ranks, rank);
                if (_next < 0) throw new InvalidOperationException($"overflow entry for rank {rank} is missing");
            }
            while (_next < ranks.Length && ranks[_next] < rank) _next++;
            if (_next == ranks.Length || ranks[_next] != rank)
                throw new InvalidOperationException($"overflow entry for rank {rank} is missing");
            return _owner._overflowValue[_next];
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
