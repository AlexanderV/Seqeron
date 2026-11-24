using System.Runtime.CompilerServices;

namespace SuffixTree.Esa;

/// <summary>
/// The range-minimum / previous-smaller / next-smaller helper over the LCP array (ADR-001 section 5, work package 1b.4; the production counterpart of the prototype
/// <c>LcpRangeHelper</c> of 1b.4a): a pyramid of block minima (fanout F, a power of two; level k holds the minimum of F^k consecutive entries) answers <see cref="Min(int, int)"/>,
/// <see cref="PrevLess(int, int)"/> and <see cref="NextLess(int, int)"/> in O(F x levels) sequential reads, with about 4 / (F - 1) bytes per entry (0.13 at F = 32). The levels are one pinned
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

    private readonly LcpKind _kind;
    private readonly bool _mapped;
    private readonly int _length;
    private readonly long _longLength;
    private readonly CompactLcp? _wideLcp;
    private readonly int _fanout;
    private readonly int _logFanout;
    private readonly int _levelCount;
    private readonly PinnedArray<byte>? _pinned8;
    private readonly PinnedArray<ushort>? _pinned16;
    private readonly PinnedArray<int>? _pinned32;
    private readonly PinnedArray<int> _levelData;
    private readonly Buf<int> _levels;
    private readonly Reader8 _reader8;
    private readonly Reader16 _reader16;
    private readonly Reader32 _reader32;
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
        _kind = lcp.Kind;
        _mapped = lcp.IsMapped;
        _longLength = lcp.Length;
        _length = (int)Math.Min(lcp.Length, int.MaxValue); // what the 32-bit cores compare with; a wide pyramid of more entries never runs them
        _fanout = fanout;
        _logFanout = int.TrailingZeroCount(fanout);

        // level sizes: ceil(length / F), ceil(that / F), ... down to one block; none for a single entry. The levels are one array of 32-bit minima: a text whose levels do not fit it is refused (a bridge of S2.3, removed by S2.5)
        long total = 0;
        int levels = 0;
        if (_longLength > 1)
        {
            long size = _longLength;
            do
            {
                size = (size + fanout - 1) / fanout;
                if (levels == MaxLevels) throw new InvalidOperationException("more pyramid levels than the table holds");
                if (total + size > int.MaxValue) throw new NotSupportedException($"the pyramid of an LCP array of {_longLength} entries has more than {int.MaxValue} 32-bit block minima: the 64-bit navigation of stage S2 (S2.5) holds it");
                _offset[levels] = (int)total;
                _count[levels] = (int)size;
                total += size;
                levels++;
            }
            while (size > 1);
        }
        _levelCount = levels;
        _levelData = PinnedArray<int>.Allocate((int)total, zero: false, "lcp-pyramid");
        _levels = _levelData.Buf;

        if (lcp.IsWide)
        {
            // an LCP array of any width, kind and size is read through its store; a block minimum above int.MaxValue cannot be a 32-bit level
            _wideLcp = lcp;
            try
            {
                FillWide(lcp);
            }
            catch
            {
                _levelData.Dispose();
                throw;
            }
            return;
        }

        if (lcp.IsMapped)
        {
            // the LCP array of a container is read where it lies; its overflow list becomes two arrays (built once, with the pack), and a list that does not match the escapes of the base array is a damaged file
            try
            {
                switch (_kind)
                {
                    case LcpKind.U8:
                        _reader8 = new Reader8(new Buf<byte>(lcp.MappedBase, _length, "lcp-u8"), lcp.OverflowRanks, lcp.OverflowValues);
                        Fill(_reader8);
                        break;
                    case LcpKind.U16:
                        _reader16 = new Reader16(new Buf<ushort>((ushort*)lcp.MappedBase, _length, "lcp-u16"), lcp.OverflowRanks, lcp.OverflowValues);
                        Fill(_reader16);
                        break;
                    default:
                        _reader32 = new Reader32(new Buf<int>((int*)lcp.MappedBase, _length, "lcp-u32"));
                        Fill(_reader32);
                        break;
                }
            }
            catch (InvalidOperationException e)
            {
                _levelData.Dispose();
                throw new InvalidDataException(e.Message + ": the file is damaged", e);
            }
            return;
        }

        switch (_kind)
        {
            case LcpKind.U8:
                _pinned8 = PinnedArray<byte>.Pin(lcp.BaseU8!, "lcp-u8");
                _reader8 = new Reader8(_pinned8.Buf, lcp.OverflowRanks, lcp.OverflowValues);
                Fill(_reader8);
                break;
            case LcpKind.U16:
                _pinned16 = PinnedArray<ushort>.Pin(lcp.BaseU16!, "lcp-u16");
                _reader16 = new Reader16(_pinned16.Buf, lcp.OverflowRanks, lcp.OverflowValues);
                Fill(_reader16);
                break;
            default:
                _pinned32 = PinnedArray<int>.Pin(lcp.BaseU32!, "lcp-u32");
                _reader32 = new Reader32(_pinned32.Buf);
                Fill(_reader32);
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
    internal long Length => _longLength;

    /// <summary>The fanout F.</summary>
    internal int Fanout => _fanout;

    /// <summary>Number of pyramid levels above the LCP array.</summary>
    internal int Levels => _levelCount;

    /// <summary>Bytes of the level arrays (the LCP array is not counted).</summary>
    internal long LevelBytes => 4L * _levelData.Length;

    /// <summary>The true LCP value at <paramref name="rank"/> (an escaped entry is resolved).</summary>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="rank"/> is outside <c>0..Length-1</c>.</exception>
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    internal int ValueAt(int rank)
    {
        if (_wideLcp is not null) return IntBridge.Value(ValueAtWide(rank), "an LCP value");
        if ((uint)rank >= (uint)_length) throw new ArgumentOutOfRangeException(nameof(rank));
        return _kind switch
        {
            LcpKind.U8 => _reader8.Value(rank),
            LcpKind.U16 => _reader16.Value(rank),
            _ => _reader32.Value(rank),
        };
    }

    /// <summary>Minimum of <c>LCP[lo..hi]</c> (inclusive).</summary>
    /// <exception cref="ArgumentOutOfRangeException">The range is not inside <c>0..Length-1</c>.</exception>
    internal int Min(int lo, int hi)
    {
        if (_wideLcp is not null) return IntBridge.Value(Min((long)lo, hi), "an LCP value");
        if ((uint)lo > (uint)hi || hi >= _length) throw new ArgumentOutOfRangeException(nameof(lo), $"[{lo}, {hi}] is not a range of {_length} entries");
        return _kind switch
        {
            LcpKind.U8 => MinCore(_reader8, lo, hi),
            LcpKind.U16 => MinCore(_reader16, lo, hi),
            _ => MinCore(_reader32, lo, hi),
        };
    }

    /// <summary>The largest <c>j &lt;= pos</c> with <c>LCP[j] &lt; bound</c>, or -1 (also for a negative <paramref name="pos"/>).</summary>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="pos"/> is at or beyond <see cref="Length"/>.</exception>
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    internal int PrevLess(int pos, int bound)
    {
        if (pos < 0) return -1;
        if (_wideLcp is not null) return IntBridge.Value(PrevLess((long)pos, bound), "a rank");
        ArgumentOutOfRangeException.ThrowIfGreaterThanOrEqual(pos, _length);
        return _kind switch
        {
            LcpKind.U8 => PrevLessCore(_reader8, pos, bound),
            LcpKind.U16 => PrevLessCore(_reader16, pos, bound),
            _ => PrevLessCore(_reader32, pos, bound),
        };
    }

    /// <summary>The smallest <c>j &gt;= pos</c> with <c>LCP[j] &lt; bound</c>, or <see cref="Length"/> (also for a <paramref name="pos"/> at or beyond it).</summary>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="pos"/> is negative.</exception>
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    internal int NextLess(int pos, int bound)
    {
        if (_wideLcp is not null) return IntBridge.Value(NextLess((long)pos, bound), "a rank");
        if (pos >= _length) return _length;
        ArgumentOutOfRangeException.ThrowIfNegative(pos);
        return _kind switch
        {
            LcpKind.U8 => NextLessCore(_reader8, pos, bound),
            LcpKind.U16 => NextLessCore(_reader16, pos, bound),
            _ => NextLessCore(_reader32, pos, bound),
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
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    internal void Borders(int a, int b, int t, out int left, out int right, out int leftValue, out int rightValue)
    {
        if (_wideLcp is not null)
        {
            BordersWide(a, b, t, out long l, out long r, out long lv, out long rv);
            left = IntBridge.Value(l, "a rank");
            right = IntBridge.Value(r, "a rank");
            leftValue = lv < 0 ? -1 : IntBridge.Value(lv, "an LCP value");
            rightValue = rv < 0 ? -1 : IntBridge.Value(rv, "an LCP value");
            return;
        }
        if ((uint)a >= (uint)_length || b < a || b >= _length) throw new ArgumentOutOfRangeException(nameof(a));
        switch (_kind)
        {
            case LcpKind.U8:
                BordersCore(_reader8, a, b, t, out left, out right, out leftValue, out rightValue);
                break;
            case LcpKind.U16:
                BordersCore(_reader16, a, b, t, out left, out right, out leftValue, out rightValue);
                break;
            default:
                BordersCore(_reader32, a, b, t, out left, out right, out leftValue, out rightValue);
                break;
        }
    }

    // ---- the 64-bit members (work package S2.3): the same queries at 64-bit ranks and bounds. A pyramid of an LCP array that the 32-bit cores can address answers them through those cores (a rank above int.MaxValue is out of range);
    // a pyramid of a wide LCP array (any width, kind and size: CompactLcp.IsWide) runs the cores below, which read the base entries through the store of the array and the block minima from the 32-bit levels.

    /// <summary>The true LCP value at a 64-bit rank.</summary>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="rank"/> is outside <c>0..Length-1</c>.</exception>
    internal long ValueAt(long rank) => _wideLcp is not null ? ValueAtWide(rank) : ValueAt(Rank(rank, nameof(rank)));

    /// <summary>Minimum of <c>LCP[lo..hi]</c> (inclusive) at 64-bit ranks.</summary>
    /// <exception cref="ArgumentOutOfRangeException">The range is not inside <c>0..Length-1</c>.</exception>
    internal long Min(long lo, long hi)
    {
        if (_wideLcp is null)
        {
            if ((ulong)hi >= (ulong)_longLength) throw new ArgumentOutOfRangeException(nameof(hi), $"[{lo}, {hi}] is not a range of {_longLength} entries");
            return Min(Rank(lo, nameof(lo)), (int)hi);
        }
        if ((ulong)lo > (ulong)hi || hi >= _longLength) throw new ArgumentOutOfRangeException(nameof(lo), $"[{lo}, {hi}] is not a range of {_longLength} entries");
        return MinWide(lo, hi);
    }

    /// <summary>The largest <c>j &lt;= pos</c> with <c>LCP[j] &lt; bound</c>, or -1 (also for a negative <paramref name="pos"/>), at 64-bit ranks and bounds.</summary>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="pos"/> is at or beyond <see cref="Length"/>.</exception>
    internal long PrevLess(long pos, long bound)
    {
        if (pos < 0) return -1;
        ArgumentOutOfRangeException.ThrowIfGreaterThanOrEqual(pos, _longLength);
        return _wideLcp is not null ? PrevLessWide(pos, bound) : PrevLess((int)pos, Bound(bound));
    }

    /// <summary>The smallest <c>j &gt;= pos</c> with <c>LCP[j] &lt; bound</c>, or <see cref="Length"/> (also for a <paramref name="pos"/> at or beyond it), at 64-bit ranks and bounds.</summary>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="pos"/> is negative.</exception>
    internal long NextLess(long pos, long bound)
    {
        if (pos >= _longLength) return _longLength;
        ArgumentOutOfRangeException.ThrowIfNegative(pos);
        return _wideLcp is not null ? NextLessWide(pos, bound) : NextLess((int)pos, Bound(bound));
    }

    /// <summary>Previous smaller value (strict) at 64-bit ranks.</summary>
    internal long PreviousSmaller(long p) => PrevLess(p - 1, ValueAt(p));

    /// <summary>Next smaller value (strict) at 64-bit ranks.</summary>
    internal long NextSmaller(long p) => NextLess(p + 1, ValueAt(p));

    /// <summary>The suffix-link kernel at 64-bit ranks and values (see the 32-bit overload).</summary>
    internal void Borders(long a, long b, long t, out long left, out long right, out long leftValue, out long rightValue)
    {
        if (_wideLcp is not null)
        {
            BordersWide(a, b, t, out left, out right, out leftValue, out rightValue);
            return;
        }
        if ((ulong)a >= (ulong)_longLength || b < a || b >= _longLength) throw new ArgumentOutOfRangeException(nameof(a));
        Borders((int)a, (int)b, Bound(t), out int l, out int r, out int lv, out int rv);
        left = l;
        right = r;
        leftValue = lv;
        rightValue = rv;
    }

    /// <summary>A rank as an <see cref="int"/> for the 32-bit cores: a pyramid that they can address has fewer than 2^31 entries, so a rank above is out of range.</summary>
    private static int Rank(long rank, string name) => (ulong)rank <= int.MaxValue ? (int)rank : throw new ArgumentOutOfRangeException(name, rank, "a rank outside the LCP array");

    /// <summary>A bound as an <see cref="int"/>: the values of an array that the 32-bit cores address are below 2^31, so a bound outside the integers decides every comparison the way its clamped value does.</summary>
    private static int Bound(long bound) => (int)Math.Clamp(bound, int.MinValue, int.MaxValue);

    private long ValueAtWide(long rank)
    {
        if ((ulong)rank >= (ulong)_longLength) throw new ArgumentOutOfRangeException(nameof(rank));
        return _wideLcp![rank];
    }

    /// <summary>The levels of a wide pyramid: the minimum of every block, read through the store; a minimum above <c>int.MaxValue</c> does not fit a 32-bit level.</summary>
    private void FillWide(CompactLcp lcp)
    {
        if (_levelCount == 0) return;
        long f = _fanout;
        for (int b = 0, blocks = _count[0]; b < blocks; b++)
        {
            long start = b * f;
            long end = Math.Min(_longLength, start + f);
            long min = long.MaxValue;
            for (long j = start; j < end; j++) min = Math.Min(min, lcp[j]);
            if (min > int.MaxValue) throw new NotSupportedException($"the LCP array has a block whose minimum is {min}, more than a 32-bit level of the pyramid holds: the 64-bit navigation of stage S2 (S2.5) holds it");
            _levels[b] = (int)min;
        }
        FillUpperLevels();
    }

    private void FillUpperLevels()
    {
        int f = _fanout;
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

    private void BordersWide(long a, long b, long t, out long left, out long right, out long leftValue, out long rightValue)
    {
        if ((ulong)a >= (ulong)_longLength || b < a || b >= _longLength) throw new ArgumentOutOfRangeException(nameof(a));
        CompactLcp lcp = _wideLcp!;
        left = PrevLessWide(a, t);
        right = NextLessWide(b + 1, t);
        leftValue = left >= 0 ? lcp[left] : -1;
        rightValue = right < _longLength ? lcp[right] : -1;
    }

    private long MinWide(long lo, long hi)
    {
        CompactLcp lcp = _wideLcp!;
        int f = _fanout;
        long min = long.MaxValue;
        long loBlock = lo >> _logFanout, hiBlock = hi >> _logFanout;
        if (loBlock == hiBlock)
        {
            for (long j = lo; j <= hi; j++) min = Math.Min(min, lcp[j]);
            return min;
        }

        // the partial blocks at both ends are read from the LCP array itself
        for (long j = lo, end = lo | (long)(f - 1); j <= end; j++) min = Math.Min(min, lcp[j]);
        for (long j = hi & ~(long)(f - 1); j <= hi; j++) min = Math.Min(min, lcp[j]);

        // the whole blocks in between, level by level: both ends contribute their partial group, then the range moves one level up (block indices are below the level sizes, which are 32-bit numbers)
        int l = (int)(loBlock + 1), r = (int)(hiBlock - 1);
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

    private long PrevLessWide(long pos, long bound)
    {
        if (pos < 0) return -1;
        CompactLcp lcp = _wideLcp!;
        int f = _fanout;
        for (long j = pos, start = pos & ~(long)(f - 1); j >= start; j--)
            if (lcp[j] < bound) return j;

        int idx = (int)(pos >> _logFanout);
        for (int k = 1; k <= _levelCount; k++)
        {
            int off = _offset[k - 1];
            for (int b = idx - 1, start = idx & ~(f - 1); b >= start; b--)
                if (_levels[off + b] < bound) return RightmostBelowWide(k, b, bound);
            idx >>= _logFanout;
        }
        return -1;
    }

    private long NextLessWide(long pos, long bound)
    {
        if (pos >= _longLength) return _longLength;
        CompactLcp lcp = _wideLcp!;
        int f = _fanout;
        for (long j = pos, end = Math.Min(pos | (long)(f - 1), _longLength - 1); j <= end; j++)
            if (lcp[j] < bound) return j;

        int idx = (int)(pos >> _logFanout);
        for (int k = 1; k <= _levelCount; k++)
        {
            int off = _offset[k - 1];
            for (int b = idx + 1, end = Math.Min(idx | (f - 1), _count[k - 1] - 1); b <= end; b++)
                if (_levels[off + b] < bound) return LeftmostBelowWide(k, b, bound);
            idx >>= _logFanout;
        }
        return _longLength;
    }

    private long RightmostBelowWide(int k, int b, long bound)
    {
        CompactLcp lcp = _wideLcp!;
        int f = _fanout;
        while (k > 1)
        {
            k--;
            int off = _offset[k - 1];
            int lowest = b * f;
            int c = Math.Min(lowest + f - 1, _count[k - 1] - 1);
            while (_levels[off + c] >= bound)
            {
                if (--c < lowest) throw Inconsistent();
            }
            b = c;
        }
        long first = (long)b * f;
        for (long position = Math.Min(first + f - 1, _longLength - 1); position >= first; position--)
            if (lcp[position] < bound) return position;
        throw Inconsistent();
    }

    private long LeftmostBelowWide(int k, int b, long bound)
    {
        CompactLcp lcp = _wideLcp!;
        int f = _fanout;
        while (k > 1)
        {
            k--;
            int off = _offset[k - 1];
            int c = b * f;
            int last = Math.Min(c + f - 1, _count[k - 1] - 1);
            while (_levels[off + c] >= bound)
            {
                if (++c > last) throw Inconsistent();
            }
            b = c;
        }
        long end = Math.Min(((long)b * f) + f - 1, _longLength - 1);
        for (long position = (long)b * f; position <= end; position++)
            if (lcp[position] < bound) return position;
        throw Inconsistent();
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
        FillUpperLevels();
    }

    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    private void BordersCore<TR>(TR reader, int a, int b, int t, out int left, out int right, out int leftValue, out int rightValue)
        where TR : struct, IBaseReader
    {
        left = PrevLessCore(reader, a, t);
        right = NextLessCore(reader, b + 1, t);
        leftValue = left >= 0 ? reader.Value(left) : -1;
        rightValue = right < _length ? reader.Value(right) : -1;
    }

    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
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

    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
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

    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
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
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
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
                if (--c < lowest) throw Inconsistent();
            }
            b = c;
        }
        int first = b * f;
        for (int position = Math.Min(first + f - 1, _length - 1); position >= first; position--)
            if (reader.Less(position, bound)) return position;
        throw Inconsistent();
    }

    /// <summary>The leftmost position below <paramref name="bound"/> inside block <paramref name="b"/> of level <paramref name="k"/>, whose minimum is known to be below it.</summary>
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
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
                if (++c > last) throw Inconsistent();
            }
            b = c;
        }
        int end = Math.Min((b * f) + f - 1, _length - 1);
        for (int position = b * f; position <= end; position++)
            if (reader.Less(position, bound)) return position;
        throw Inconsistent();
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

    /// <summary>The pyramid and the entries it was built from disagree (a block whose minimum is below the bound holds no entry below it): never for an LCP array in memory, a damaged overflow list for one over a file.</summary>
    private Exception Inconsistent()
    {
        const string message = "the LCP pyramid is inconsistent: no entry below the bound in a block whose minimum is below it";
        return _mapped ? new InvalidDataException(message + ": the file is damaged") : new InvalidOperationException(message);
    }

    private static int Resolve(int[] overflowRanks, int[] overflowValues, int rank)
    {
        int k = Array.BinarySearch(overflowRanks, rank);
        if (k < 0) throw new InvalidOperationException($"overflow entry for rank {rank} is missing");
        return overflowValues[k];
    }
}
