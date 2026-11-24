namespace SuffixTree.Esa;

/// <summary>
/// The wide form of a <see cref="CompactLcp"/> (work package S2.3): an LCP array read through position readers, at any stored width, in any of the four kinds (U8, U16, U32, U64), with 64-bit ranks and values and an overflow list of pairs
/// <c>(rank, value)</c> whose fields are <c>POS_WIDTH</c> wide (a pair is 8 bytes at width 4, 10 at width 5). The class is abstract so that <see cref="CompactLcp"/> (a non-generic class) can hold one: the cold paths of the container reader and
/// the generator-backed tests; the pattern descent of <see cref="MappedSearch"/> reads the base array through its own pointer readers and calls <see cref="Value"/> only for an escaped entry.
/// </summary>
internal abstract class LcpStore
{
    /// <summary>The kind.</summary>
    internal abstract LcpKind Kind { get; }

    /// <summary>Number of entries, <c>n + 1</c>.</summary>
    internal abstract long Length { get; }

    /// <summary>Number of pairs of the overflow list.</summary>
    internal abstract long PairCount { get; }

    /// <summary>Bytes of one field of a pair (4 or 5).</summary>
    internal abstract int PosWidth { get; }

    /// <summary>The base entry at <paramref name="rank"/> as it is stored (an escape is not resolved); a 32- or 64-bit value is checked against the text in Checked mode.</summary>
    /// <exception cref="InvalidDataException">The rank is outside the array, or a checked value is larger than the text.</exception>
    internal abstract long BaseEntry(long rank);

    /// <summary>The LCP value at <paramref name="rank"/> (an escape resolved through the overflow list).</summary>
    /// <exception cref="InvalidDataException">The rank is outside the array, the escape has no pair, or a value is damaged.</exception>
    internal abstract long Value(long rank);

    /// <summary>The rank field of pair <paramref name="index"/>.</summary>
    internal abstract long PairRank(long index);

    /// <summary>The value field of pair <paramref name="index"/>.</summary>
    internal abstract long PairValue(long index);

    /// <summary>The index of the pair of <paramref name="rank"/> (binary search over the ranks of the pairs, which are ascending), or -1.</summary>
    internal abstract long Find(long rank);

    /// <summary>Bytes of one entry of the base array of a kind.</summary>
    internal static int WidthOf(LcpKind kind) => kind switch
    {
        LcpKind.U8 => 1,
        LcpKind.U16 => 2,
        LcpKind.U32 => 4,
        _ => 8,
    };
}

/// <inheritdoc cref="LcpStore"/>
internal sealed class LcpStore<TBase, TPairs> : LcpStore
    where TBase : struct, IPosReader
    where TPairs : struct, IPosReader
{
    private readonly LcpKind _kind;
    private readonly TBase _base;
    private readonly TPairs _pairs;
    private readonly long _length;
    private readonly long _pairCount;
    private readonly bool _checked;
    private readonly ulong _escape;
    private readonly bool _hasEscape;

    internal LcpStore(LcpKind kind, TBase baseEntries, TPairs pairs, bool isChecked)
    {
        _kind = kind;
        _base = baseEntries;
        _pairs = pairs;
        _length = baseEntries.Count;
        _pairCount = pairs.Count >> 1;
        _checked = isChecked;
        _hasEscape = kind is LcpKind.U8 or LcpKind.U16;
        _escape = kind == LcpKind.U8 ? byte.MaxValue : ushort.MaxValue;
    }

    internal override LcpKind Kind => _kind;

    internal override long Length => _length;

    internal override long PairCount => _pairCount;

    internal override int PosWidth => TPairs.ElementBytes;

    internal override long BaseEntry(long rank)
    {
        if ((ulong)rank >= (ulong)_length) throw new InvalidDataException($"the LCP array is read at rank {rank}, but it has {_length} entries: the file is damaged");
        ulong v = _base.Get(rank);
        if (_hasEscape) return (long)v;
        if (v > long.MaxValue || (_checked && v > (ulong)_length)) throw new InvalidDataException($"the LCP value {v} at rank {rank} is larger than the text: the file is damaged");
        return (long)v;
    }

    internal override long Value(long rank)
    {
        long b = BaseEntry(rank);
        if (!_hasEscape || (ulong)b != _escape) return b;
        long k = Find(rank);
        if (k < 0) throw new InvalidDataException($"overflow entry for rank {rank} is missing: the file is damaged");
        return PairValue(k);
    }

    internal override long PairRank(long index) => Field(2 * index, "the overflow list");

    internal override long PairValue(long index) => Field((2 * index) + 1, "the overflow list");

    internal override long Find(long rank)
    {
        long lo = 0, hi = _pairCount - 1;
        while (lo <= hi)
        {
            long mid = lo + ((hi - lo) >> 1);
            ulong v = _pairs.Get(2 * mid);
            if (v == (ulong)rank) return mid;
            if (v < (ulong)rank) lo = mid + 1;
            else hi = mid - 1;
        }
        return -1;
    }

    /// <summary>A field of the overflow list as a non-negative <see cref="long"/> (a field above <c>long.MaxValue</c> cannot be a rank or an LCP value).</summary>
    private long Field(long element, string what)
    {
        ulong v = _pairs.Get(element);
        if (v > long.MaxValue) throw new InvalidDataException($"{what} holds {v}: the file is damaged");
        return (long)v;
    }
}
