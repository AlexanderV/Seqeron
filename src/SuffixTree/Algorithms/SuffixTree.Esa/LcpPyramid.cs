using System.Runtime.CompilerServices;

namespace SuffixTree.Esa;

/// <summary>
/// The range-minimum / previous-smaller / next-smaller helper over the LCP array (ADR-001 section 5, work package 1b.4; the production counterpart of the prototype
/// <c>LcpRangeHelper</c> of 1b.4a): a pyramid of block minima (fanout F, a power of two; level k holds the minimum of F^k consecutive entries) answers <see cref="Min"/>,
/// <see cref="PrevLess"/> and <see cref="NextLess"/> in O(F x levels) sequential reads, with about 4 / (F - 1) bytes per entry (0.13 at F = 32). The levels are one pinned
/// <c>int</c> array; the LCP entries themselves are NOT copied: they are read through a pointer to the base array of the <see cref="CompactLcp"/> (U8, U16 or U32), pinned for the life of
/// this object. An escaped entry (U8: 255, U16: 65,535, whose true value is in the overflow list) is resolved only when a comparison needs it: <c>value &lt; bound</c> is decided by the
/// escape alone while <c>bound &lt;= escape</c>, which is every comparison of a suffix-link search at depth up to 255 (the common case); deeper nodes pay one binary search per escaped entry
/// they look at. Every loop runs inside a block or a level group, so a search reads at most F x (levels + 1) entries; an inconsistent pyramid ends in an
/// <see cref="InvalidOperationException"/>, never in a read outside the arrays. Immutable after <see cref="Build"/>; safe for concurrent readers.
/// </summary>
internal sealed unsafe class LcpPyramid : IDisposable
{
    /// <summary>The fanout used by the navigation pack (0.13 bytes per entry; measured speed does not depend on it in 8-64, Log 1b.4a).</summary>
    internal const int DefaultFanout = 32;

    private const int U8Escape = byte.MaxValue;
    private const int U16Escape = ushort.MaxValue;
    private const int MaxLevels = 32;

    private readonly CompactLcp _lcp;
    private readonly LcpKind _kind;
    private readonly int _length;
    private readonly int _fanout;
    private readonly int _logFanout;
    private readonly int _levelCount;
    private readonly PinnedArray<byte>? _pinned8;
    private readonly PinnedArray<ushort>? _pinned16;
    private readonly PinnedArray<int>? _pinned32;
    private readonly PinnedArray<int> _levelData;
    private readonly Buf<int> _levels;
#pragma warning disable S3459, S1144 // the fields of an inline array are only reachable through its indexer
    private LevelTable _offset;
    private LevelTable _count;

    [InlineArray(MaxLevels)]
    private struct LevelTable
    {
        private int _element;
    }
#pragma warning restore S3459, S1144

    private LcpPyramid(CompactLcp lcp, int fanout)
    {
        _lcp = lcp;
        _kind = lcp.Kind;
        _length = lcp.Length;
        _fanout = fanout;
        _logFanout = int.TrailingZeroCount(fanout);

        // level sizes: ceil(length / F), ceil(that / F), ... down to one block; none for a single entry
        long total = 0;
        int levels = 0;
        if (_length > 1)
        {
            int size = _length;
            do
            {
                size = (int)(((long)size + fanout - 1) / fanout);
                if (levels == MaxLevels) throw new InvalidOperationException("more pyramid levels than the table holds");
                _offset[levels] = (int)total;
                _count[levels] = size;
                total += size;
                levels++;
            }
            while (size > 1);
        }
        _levelCount = levels;
        _levelData = PinnedArray<int>.Allocate(checked((int)total), zero: false, "lcp-pyramid");
        _levels = _levelData.Buf;

        switch (_kind)
        {
            case LcpKind.U8:
                _pinned8 = PinnedArray<byte>.Pin(lcp.BaseU8!, "lcp-u8");
                Fill(new Reader8(_pinned8.Buf, lcp.OverflowRanks, lcp.OverflowValues));
                break;
            case LcpKind.U16:
                _pinned16 = PinnedArray<ushort>.Pin(lcp.BaseU16!, "lcp-u16");
                Fill(new Reader16(_pinned16.Buf, lcp.OverflowRanks, lcp.OverflowValues));
                break;
            default:
                _pinned32 = PinnedArray<int>.Pin(lcp.BaseU32!, "lcp-u32");
                Fill(new Reader32(_pinned32.Buf));
                break;
        }
    }

    /// <summary>Builds the pyramid over <paramref name="lcp"/> (n + 1 entries); the base array is pinned until <see cref="Dispose"/>.</summary>
    /// <exception cref="ArgumentNullException"><paramref name="lcp"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException">The fanout is not a power of two of at least 2.</exception>
    internal static LcpPyramid Build(CompactLcp lcp, int fanout = DefaultFanout)
    {
        ArgumentNullException.ThrowIfNull(lcp);
        if (fanout < 2 || (fanout & (fanout - 1)) != 0)
            throw new ArgumentOutOfRangeException(nameof(fanout), fanout, "the fanout must be a power of two of at least 2");
        return new LcpPyramid(lcp, fanout);
    }

    /// <summary>Number of LCP entries (n + 1).</summary>
    internal int Length => _length;

    /// <summary>The fanout F.</summary>
    internal int Fanout => _fanout;

    /// <summary>Number of pyramid levels above the LCP array.</summary>
    internal int Levels => _levelCount;

    /// <summary>Bytes of the level arrays (the LCP array is not counted).</summary>
    internal long LevelBytes => 4L * _levelData.Length;

    /// <summary>The true LCP value at <paramref name="rank"/> (an escaped entry is resolved).</summary>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="rank"/> is outside <c>0..Length-1</c>.</exception>
    internal int ValueAt(int rank)
    {
        if ((uint)rank >= (uint)_length) throw new ArgumentOutOfRangeException(nameof(rank));
        return _kind switch
        {
            LcpKind.U8 => new Reader8(_pinned8!.Buf, _lcp.OverflowRanks, _lcp.OverflowValues).Value(rank),
            LcpKind.U16 => new Reader16(_pinned16!.Buf, _lcp.OverflowRanks, _lcp.OverflowValues).Value(rank),
            _ => new Reader32(_pinned32!.Buf).Value(rank),
        };
    }

    /// <summary>Minimum of <c>LCP[lo..hi]</c> (inclusive).</summary>
    /// <exception cref="ArgumentOutOfRangeException">The range is not inside <c>0..Length-1</c>.</exception>
    internal int Min(int lo, int hi)
    {
        if ((uint)lo > (uint)hi || hi >= _length) throw new ArgumentOutOfRangeException(nameof(lo), $"[{lo}, {hi}] is not a range of {_length} entries");
        return _kind switch
        {
            LcpKind.U8 => MinCore(new Reader8(_pinned8!.Buf, _lcp.OverflowRanks, _lcp.OverflowValues), lo, hi),
            LcpKind.U16 => MinCore(new Reader16(_pinned16!.Buf, _lcp.OverflowRanks, _lcp.OverflowValues), lo, hi),
            _ => MinCore(new Reader32(_pinned32!.Buf), lo, hi),
        };
    }

    /// <summary>The largest <c>j &lt;= pos</c> with <c>LCP[j] &lt; bound</c>, or -1 (also for a negative <paramref name="pos"/>).</summary>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="pos"/> is at or beyond <see cref="Length"/>.</exception>
    internal int PrevLess(int pos, int bound)
    {
        if (pos < 0) return -1;
        ArgumentOutOfRangeException.ThrowIfGreaterThanOrEqual(pos, _length);
        return _kind switch
        {
            LcpKind.U8 => PrevLessCore(new Reader8(_pinned8!.Buf, _lcp.OverflowRanks, _lcp.OverflowValues), pos, bound),
            LcpKind.U16 => PrevLessCore(new Reader16(_pinned16!.Buf, _lcp.OverflowRanks, _lcp.OverflowValues), pos, bound),
            _ => PrevLessCore(new Reader32(_pinned32!.Buf), pos, bound),
        };
    }

    /// <summary>The smallest <c>j &gt;= pos</c> with <c>LCP[j] &lt; bound</c>, or <see cref="Length"/> (also for a <paramref name="pos"/> at or beyond it).</summary>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="pos"/> is negative.</exception>
    internal int NextLess(int pos, int bound)
    {
        if (pos >= _length) return _length;
        ArgumentOutOfRangeException.ThrowIfNegative(pos);
        return _kind switch
        {
            LcpKind.U8 => NextLessCore(new Reader8(_pinned8!.Buf, _lcp.OverflowRanks, _lcp.OverflowValues), pos, bound),
            LcpKind.U16 => NextLessCore(new Reader16(_pinned16!.Buf, _lcp.OverflowRanks, _lcp.OverflowValues), pos, bound),
            _ => NextLessCore(new Reader32(_pinned32!.Buf), pos, bound),
        };
    }

    /// <summary>Previous smaller value (strict): the largest <c>j &lt; p</c> with <c>LCP[j] &lt; LCP[p]</c>, or -1.</summary>
    internal int PreviousSmaller(int p) => PrevLess(p - 1, ValueAt(p));

    /// <summary>Next smaller value (strict): the smallest <c>j &gt; p</c> with <c>LCP[j] &lt; LCP[p]</c>, or <see cref="Length"/>.</summary>
    internal int NextSmaller(int p) => NextLess(p + 1, ValueAt(p));

    /// <summary>The suffix-link kernel (c2) in one call: <c>PrevLess(a, t)</c> and <c>NextLess(b + 1, t)</c> with the readers resolved once, so that the kind switch is paid once per link.</summary>
    /// <param name="a">Left endpoint, <c>Psi(l)</c>.</param>
    /// <param name="b">Right endpoint, <c>Psi(r)</c> (at least <paramref name="a"/>).</param>
    /// <param name="t">The value of the target interval, <c>depth - 1</c>.</param>
    /// <param name="left">The largest <c>j &lt;= a</c> with <c>LCP[j] &lt; t</c>.</param>
    /// <param name="right">The smallest <c>j &gt; b</c> with <c>LCP[j] &lt; t</c>, or <see cref="Length"/>.</param>
    /// <param name="leftValue">The true <c>LCP[left]</c>.</param>
    /// <param name="rightValue">The true <c>LCP[right]</c>, or -1 when <paramref name="right"/> is <see cref="Length"/>.</param>
    internal void Borders(int a, int b, int t, out int left, out int right, out int leftValue, out int rightValue)
    {
        if ((uint)a >= (uint)_length || b < a || b >= _length) throw new ArgumentOutOfRangeException(nameof(a));
        switch (_kind)
        {
            case LcpKind.U8:
                BordersCore(new Reader8(_pinned8!.Buf, _lcp.OverflowRanks, _lcp.OverflowValues), a, b, t, out left, out right, out leftValue, out rightValue);
                break;
            case LcpKind.U16:
                BordersCore(new Reader16(_pinned16!.Buf, _lcp.OverflowRanks, _lcp.OverflowValues), a, b, t, out left, out right, out leftValue, out rightValue);
                break;
            default:
                BordersCore(new Reader32(_pinned32!.Buf), a, b, t, out left, out right, out leftValue, out rightValue);
                break;
        }
    }

    /// <summary>Releases the pins; idempotent. The pyramid must not be used afterwards.</summary>
    public void Dispose()
    {
        _pinned8?.Dispose();
        _pinned16?.Dispose();
        _pinned32?.Dispose();
        _levelData.Dispose();
    }

    private void Fill<TR>(TR reader)
        where TR : struct, IBaseReader
    {
        if (_levelCount == 0) return;
        int f = _fanout;
        for (int b = 0, blocks = _count[0]; b < blocks; b++)
        {
            int start = b * f;
            int end = Math.Min(_length, start + f);
            int min = int.MaxValue;
            for (int j = start; j < end; j++) min = Math.Min(min, reader.Value(j));
            _levels[b] = min;
        }
        for (int k = 1; k < _levelCount; k++)
        {
            int below = _offset[k - 1], belowCount = _count[k - 1], here = _offset[k];
            for (int b = 0, blocks = _count[k]; b < blocks; b++)
            {
                int start = b * f;
                int end = Math.Min(belowCount, start + f);
                int min = int.MaxValue;
                for (int j = start; j < end; j++) min = Math.Min(min, _levels[below + j]);
                _levels[here + b] = min;
            }
        }
    }

    private void BordersCore<TR>(TR reader, int a, int b, int t, out int left, out int right, out int leftValue, out int rightValue)
        where TR : struct, IBaseReader
    {
        left = PrevLessCore(reader, a, t);
        right = NextLessCore(reader, b + 1, t);
        leftValue = left >= 0 ? reader.Value(left) : -1;
        rightValue = right < _length ? reader.Value(right) : -1;
    }

    private int MinCore<TR>(TR reader, int lo, int hi)
        where TR : struct, IBaseReader
    {
        int f = _fanout;
        int min = int.MaxValue;
        int loBlock = lo >> _logFanout, hiBlock = hi >> _logFanout;
        if (loBlock == hiBlock)
        {
            for (int j = lo; j <= hi; j++) min = Math.Min(min, reader.Value(j));
            return min;
        }

        // the partial blocks at both ends are read from the LCP array itself
        for (int j = lo, end = lo | (f - 1); j <= end; j++) min = Math.Min(min, reader.Value(j));
        for (int j = hi & ~(f - 1); j <= hi; j++) min = Math.Min(min, reader.Value(j));

        // the whole blocks in between, level by level: both ends contribute their partial group, then the range moves one level up
        int l = loBlock + 1, r = hiBlock - 1;
        int k = 1;
        while (l <= r)
        {
            int off = _offset[k - 1];
            if ((l >> _logFanout) == (r >> _logFanout))
            {
                for (int b = l; b <= r; b++) min = Math.Min(min, _levels[off + b]);
                break;
            }
            for (int b = l, end = l | (f - 1); b <= end; b++) min = Math.Min(min, _levels[off + b]);
            for (int b = r & ~(f - 1); b <= r; b++) min = Math.Min(min, _levels[off + b]);
            l = (l >> _logFanout) + 1;
            r = (r >> _logFanout) - 1;
            k++;
        }
        return min;
    }

    private int PrevLessCore<TR>(TR reader, int pos, int bound)
        where TR : struct, IBaseReader
    {
        int f = _fanout;
        for (int j = pos, start = pos & ~(f - 1); j >= start; j--)
            if (reader.Less(j, bound)) return j;

        int idx = pos >> _logFanout;
        for (int k = 1; k <= _levelCount; k++)
        {
            int off = _offset[k - 1];
            for (int b = idx - 1, start = idx & ~(f - 1); b >= start; b--)
                if (_levels[off + b] < bound) return RightmostBelow(reader, k, b, bound);
            idx >>= _logFanout;
        }
        return -1;
    }

    private int NextLessCore<TR>(TR reader, int pos, int bound)
        where TR : struct, IBaseReader
    {
        if (pos >= _length) return _length;
        int f = _fanout;
        for (int j = pos, end = Math.Min(pos | (f - 1), _length - 1); j <= end; j++)
            if (reader.Less(j, bound)) return j;

        int idx = pos >> _logFanout;
        for (int k = 1; k <= _levelCount; k++)
        {
            int off = _offset[k - 1];
            for (int b = idx + 1, end = Math.Min(idx | (f - 1), _count[k - 1] - 1); b <= end; b++)
                if (_levels[off + b] < bound) return LeftmostBelow(reader, k, b, bound);
            idx >>= _logFanout;
        }
        return _length;
    }

    /// <summary>The rightmost position below <paramref name="bound"/> inside block <paramref name="b"/> of level <paramref name="k"/>, whose minimum is known to be below it.</summary>
    private int RightmostBelow<TR>(TR reader, int k, int b, int bound)
        where TR : struct, IBaseReader
    {
        int f = _fanout;
        while (k > 1)
        {
            k--;
            int off = _offset[k - 1];
            int lowest = b * f;
            int c = Math.Min(lowest + f - 1, _count[k - 1] - 1);
            while (_levels[off + c] >= bound)
            {
                if (--c < lowest) throw new InvalidOperationException("the LCP pyramid is inconsistent: no entry below the bound in a block whose minimum is below it");
            }
            b = c;
        }
        int first = b * f;
        for (int position = Math.Min(first + f - 1, _length - 1); position >= first; position--)
            if (reader.Less(position, bound)) return position;
        throw new InvalidOperationException("the LCP pyramid is inconsistent: no entry below the bound in a block whose minimum is below it");
    }

    /// <summary>The leftmost position below <paramref name="bound"/> inside block <paramref name="b"/> of level <paramref name="k"/>, whose minimum is known to be below it.</summary>
    private int LeftmostBelow<TR>(TR reader, int k, int b, int bound)
        where TR : struct, IBaseReader
    {
        int f = _fanout;
        while (k > 1)
        {
            k--;
            int off = _offset[k - 1];
            int c = b * f;
            int last = Math.Min(c + f - 1, _count[k - 1] - 1);
            while (_levels[off + c] >= bound)
            {
                if (++c > last) throw new InvalidOperationException("the LCP pyramid is inconsistent: no entry below the bound in a block whose minimum is below it");
            }
            b = c;
        }
        int end = Math.Min((b * f) + f - 1, _length - 1);
        for (int position = b * f; position <= end; position++)
            if (reader.Less(position, bound)) return position;
        throw new InvalidOperationException("the LCP pyramid is inconsistent: no entry below the bound in a block whose minimum is below it");
    }

    /// <summary>Read access to the base entries of one <see cref="LcpKind"/>; a struct, so that each search is specialised per kind.</summary>
    private interface IBaseReader
    {
        /// <summary>The true value at <paramref name="rank"/>.</summary>
        int Value(int rank);

        /// <summary>Whether the true value at <paramref name="rank"/> is below <paramref name="bound"/>.</summary>
        bool Less(int rank, int bound);
    }

    private readonly struct Reader8(Buf<byte> entries, int[] overflowRanks, int[] overflowValues) : IBaseReader
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public int Value(int rank)
        {
            int v = entries[rank];
            return v != U8Escape ? v : Resolve(overflowRanks, overflowValues, rank);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public bool Less(int rank, int bound)
        {
            int v = entries[rank];
            if (v != U8Escape) return v < bound;
            return bound > U8Escape && Resolve(overflowRanks, overflowValues, rank) < bound;
        }
    }

    private readonly struct Reader16(Buf<ushort> entries, int[] overflowRanks, int[] overflowValues) : IBaseReader
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public int Value(int rank)
        {
            int v = entries[rank];
            return v != U16Escape ? v : Resolve(overflowRanks, overflowValues, rank);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public bool Less(int rank, int bound)
        {
            int v = entries[rank];
            if (v != U16Escape) return v < bound;
            return bound > U16Escape && Resolve(overflowRanks, overflowValues, rank) < bound;
        }
    }

    private readonly struct Reader32(Buf<int> entries) : IBaseReader
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public int Value(int rank) => entries[rank];

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public bool Less(int rank, int bound) => entries[rank] < bound;
    }

    private static int Resolve(int[] overflowRanks, int[] overflowValues, int rank)
    {
        int k = Array.BinarySearch(overflowRanks, rank);
        if (k < 0) throw new InvalidOperationException($"overflow entry for rank {rank} is missing");
        return overflowValues[k];
    }
}
