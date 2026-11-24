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

    /// <summary>Eight bytes per entry, no overflow list (LCP codec 4: a text above 2^32 with a repeat longer than 2^32 - 2); read by the wide forms only.</summary>
    U64,
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

    // the wide form (work package S2.3): any stored width, the U64 kind, any length; every read goes through the store, the fields of the forms above are unused
    private readonly LcpStore? _store;

    private CompactLcp(LcpKind kind, int length, byte[]? u8, ushort[]? u16, int[]? u32, int[] overflowRank, int[] overflowValue)
    {
        if (kind == LcpKind.U64) throw new ArgumentOutOfRangeException(nameof(kind), kind, "the managed form has no U64 kind (LCP codec 4 is read from a container)");
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
        if (kind == LcpKind.U64) throw new ArgumentOutOfRangeException(nameof(kind), kind, "the 32-bit form reads no U64 array: use FromMappedWide");
        if (kind == LcpKind.U32 && pairCount != 0) throw new InvalidDataException("an LCP array of 32 bits per entry has no overflow list");
        return new CompactLcp(kind, length, baseArray, pairs, pairCount, isChecked);
    }

    private CompactLcp(LcpStore store)
    {
        Kind = store.Kind;
        Length = store.Length;
        _overflowRank = [];
        _overflowValue = [];
        _store = store;
        IsMapped = true;
    }

    /// <summary>
    /// The LCP array of a container (or of any memory) at any stored width and size (work package S2.3): <paramref name="baseArray"/> holds <paramref name="length"/> entries of the width of <paramref name="kind"/> (1, 2, 4 or 8 bytes),
    /// <paramref name="pairs"/> the <paramref name="pairCount"/> pairs <c>(rank, value)</c> of the overflow list, each field <paramref name="posWidth"/> (4 or 5) bytes, so a pair is <c>2 x POS_WIDTH</c> bytes (8 or 10) and a list of 5-byte fields
    /// ends with the 3 bytes of padding that belong to its section. The memory belongs to the container and must outlive this object. The entries are untrusted, as in <see cref="FromMapped"/>.
    /// </summary>
    /// <exception cref="ArgumentOutOfRangeException">The length is below 1, the pair count is negative or the width is not 4 or 5.</exception>
    /// <exception cref="InvalidDataException">An array of 32 or 64 bits per entry has no overflow list.</exception>
    internal static CompactLcp FromMappedWide(LcpKind kind, long length, byte* baseArray, byte* pairs, long pairCount, int posWidth, bool isChecked)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(length, 1);
        ArgumentOutOfRangeException.ThrowIfNegative(pairCount);
        if (posWidth is not (4 or 5)) throw new ArgumentOutOfRangeException(nameof(posWidth), posWidth, "the stored width of a position is 4 or 5 bytes");
        if (kind is LcpKind.U32 or LcpKind.U64 && pairCount != 0) throw new InvalidDataException("an LCP array of 32 or 64 bits per entry has no overflow list");
        long fields = 2 * pairCount;
        return kind switch
        {
            LcpKind.U8 => WideOf(kind, new Pos8<TrustedAccess>(baseArray, length), posWidth, pairs, fields, isChecked),
            LcpKind.U16 => WideOf(kind, new Pos16<TrustedAccess>(baseArray, length), posWidth, pairs, fields, isChecked),
            LcpKind.U32 => WideOf(kind, new Pos32<TrustedAccess>(baseArray, length), posWidth, pairs, fields, isChecked),
            _ => WideOf(kind, new Pos64<TrustedAccess>(baseArray, length), posWidth, pairs, fields, isChecked),
        };
    }

    private static CompactLcp WideOf<TBase>(LcpKind kind, TBase baseEntries, int posWidth, byte* pairs, long fields, bool isChecked)
        where TBase : struct, IPosReader =>
        posWidth == 5
            ? FromReaders(kind, baseEntries, new Pos40<TrustedAccess>(pairs, fields), isChecked)
            : FromReaders(kind, baseEntries, new Pos32<TrustedAccess>(pairs, fields), isChecked);

    /// <summary>
    /// The LCP array over readers (work package S2.3; the generator-backed tests run it at n above 2^31): <paramref name="baseEntries"/> has <c>n + 1</c> elements of the width of <paramref name="kind"/>, <paramref name="pairs"/> has
    /// <c>2 x pairCount</c> elements (rank, value, rank, value, ...). A rank outside the array is an <see cref="InvalidDataException"/> in both modes; with <paramref name="isChecked"/> every 32- and 64-bit value is also checked against the text.
    /// </summary>
    /// <exception cref="ArgumentException">The reader's element width is not the width of the kind, or the list has an odd number of elements.</exception>
    /// <exception cref="InvalidDataException">An array of 32 or 64 bits per entry has no overflow list.</exception>
    internal static CompactLcp FromReaders<TBase, TPairs>(LcpKind kind, TBase baseEntries, TPairs pairs, bool isChecked)
        where TBase : struct, IPosReader
        where TPairs : struct, IPosReader
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(baseEntries.Count, 1);
        if (TBase.ElementBytes != LcpStore.WidthOf(kind)) throw new ArgumentException($"the base array of {kind} has {LcpStore.WidthOf(kind)} bytes per entry, the reader has {TBase.ElementBytes}", nameof(baseEntries));
        if ((pairs.Count & 1) != 0) throw new ArgumentException("the overflow list has pairs: an even number of elements", nameof(pairs));
        if (kind is LcpKind.U32 or LcpKind.U64 && pairs.Count != 0) throw new InvalidDataException("an LCP array of 32 or 64 bits per entry has no overflow list");
        return new CompactLcp(new LcpStore<TBase, TPairs>(kind, baseEntries, pairs, isChecked));
    }

    /// <summary>True when the array is read in place from a container (or through readers).</summary>
    internal bool IsMapped { get; }

    /// <summary>True when the array is the wide form: read through a store, at any width, size and kind.</summary>
    internal bool IsWide => _store is not null;

    /// <summary>The mapped base array (null for a managed one).</summary>
    internal byte* MappedBase => _mappedBase;

    /// <summary>The representation.</summary>
    internal LcpKind Kind { get; }

    /// <summary>Number of entries, <c>n + 1</c>.</summary>
    internal long Length { get; }

    /// <summary>Number of escaped entries (always 0 for U32 and U64).</summary>
    internal long OverflowCount
    {
        get
        {
            if (_store is not null) return _store.PairCount;
            return IsMapped ? _mappedPairCount : _overflowRank.Length;
        }
    }

    /// <summary>The overflow count as the 32-bit engine reads it (the wide form can hold more than <see cref="int.MaxValue"/> pairs).</summary>
    private int OverflowCountInt
    {
        get
        {
            if (_store is not null) return IntBridge.Length(_store.PairCount, "the overflow list of the LCP array");
            return IsMapped ? _mappedPairCount : _overflowRank.Length;
        }
    }

    /// <summary>Bytes of the base array plus one pair per overflow entry (8 bytes, or <c>2 x POS_WIDTH</c> of the wide form).</summary>
    internal long PayloadBytes => (Length * LcpStore.WidthOf(Kind)) + ((_store is null ? 8L : 2L * _store.PosWidth) * OverflowCount);

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

    /// <summary>The LCP value at a 64-bit rank (any kind, any stored width; the wide form reads values above <see cref="int.MaxValue"/>).</summary>
    /// <exception cref="InvalidDataException">The rank is outside the array (a mapped form), an escape has no pair, or a value is damaged.</exception>
    /// <exception cref="ArgumentOutOfRangeException">The rank is outside the array (a managed form).</exception>
    internal long this[long rank]
    {
        get
        {
            if (_store is not null) return _store.Value(rank);
            if ((ulong)rank >= (ulong)Length) ThrowRankOutside(rank);
            return this[(int)rank];
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
        if (_store is not null)
        {
            long found = _store.Find(rank);
            return found < 0 ? -1 : IntBridge.Length(found, "an index of the overflow list of the LCP array");
        }
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

    private int OverflowRankAt(int index)
    {
        if (_store is not null) return IntBridge.Value(_store.PairRank(index), "a rank of the overflow list of the LCP array");
        return IsMapped ? MappedWord(2 * (long)index) : _overflowRank[index];
    }

    private int OverflowValueAt(int index)
    {
        if (_store is not null) return IntBridge.Value(_store.PairValue(index), "a value of the overflow list of the LCP array");
        return IsMapped ? MappedWord((2 * (long)index) + 1) : _overflowValue[index];
    }

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
            var ranks = new int[OverflowCountInt];
            var values = new int[ranks.Length];
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
        var result = new int[IntBridge.Length(Length, "the LCP array (ToArray)")];
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
        if (kind == LcpKind.U64) throw new ArgumentOutOfRangeException(nameof(kind), kind, "the managed form has no U64 kind (LCP codec 4 is read from a container)");
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
    private void ThrowRankOutside(long rank)
    {
        if (IsMapped) throw new InvalidDataException($"the LCP array is read at rank {rank}, but it has {Length} entries: the file is damaged");
        throw new ArgumentOutOfRangeException(nameof(rank));
    }

    private int BaseValue(int rank)
    {
        if ((ulong)rank >= (ulong)Length) ThrowRankOutside(rank);
        if (IsMapped) return _store is null ? MappedValue(rank) : IntBridge.Value(_store.BaseEntry(rank), "an LCP value");
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
                if (_mappedChecked && (uint)v > (ulong)Length) throw Corrupt($"the LCP value {(uint)v} at rank {rank} is larger than the text");
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

            int count = _owner.OverflowCountInt;
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
            if (kind == LcpKind.U64) throw new ArgumentOutOfRangeException(nameof(kind), kind, "the managed form has no U64 kind (LCP codec 4 is read from a container)");
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
