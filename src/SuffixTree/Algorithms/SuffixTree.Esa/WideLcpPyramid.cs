using System.Runtime.CompilerServices;

namespace SuffixTree.Esa;

/// <summary>
/// The range-minimum / previous-smaller / next-smaller queries of the navigation pack at 64-bit ranks and values (work package S2.5): the members of <see cref="LcpPyramid"/> that the long contract needs, with no 32-bit bridge. A text
/// above 2^31 symbols has ranks, block indices and minima above <see cref="int.MaxValue"/> (<c>a^n</c>: the LCP value at rank r is r - 1).
/// </summary>
internal interface ILcpRanges : IDisposable
{
    /// <summary>Number of LCP entries (<c>n + 1</c>).</summary>
    long Length { get; }

    /// <summary>The fanout of the pyramid.</summary>
    int Fanout { get; }

    /// <summary>Bytes of the levels that the pyramid holds in memory (0 when it reads the persisted RMQ section where it lies).</summary>
    long LevelBytes { get; }

    /// <summary>The true LCP value at <paramref name="rank"/>.</summary>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="rank"/> is outside <c>0..Length-1</c>.</exception>
    long ValueAt(long rank);

    /// <summary>Minimum of <c>LCP[lo..hi]</c> (inclusive).</summary>
    /// <exception cref="ArgumentOutOfRangeException">The range is not inside <c>0..Length-1</c>.</exception>
    long Min(long lo, long hi);

    /// <summary>The largest <c>j &lt;= pos</c> with <c>LCP[j] &lt; bound</c>, or -1 (also for a negative <paramref name="pos"/>).</summary>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="pos"/> is at or beyond <see cref="Length"/>.</exception>
    long PrevLess(long pos, long bound);

    /// <summary>The smallest <c>j &gt;= pos</c> with <c>LCP[j] &lt; bound</c>, or <see cref="Length"/> (also for a <paramref name="pos"/> at or beyond it).</summary>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="pos"/> is negative.</exception>
    long NextLess(long pos, long bound);

    /// <summary>The suffix-link kernel in one call: <c>PrevLess(a, t)</c> and <c>NextLess(b + 1, t)</c> with the true LCP values at the two borders (-1 when there is none).</summary>
    void Borders(long a, long b, long t, out long left, out long right, out long leftValue, out long rightValue);
}

/// <summary>
/// The pyramid of block minima (fanout F, a power of two; level k holds the minimum of F^(k+1) consecutive entries) over an LCP array of any width, kind and length, with 64-bit levels (work package S2.5). The levels are either read where
/// they lie (the persisted RMQ section of a container, <see cref="PosSource"/> over its entries: 4 or 5 bytes each) or built in memory (<c>long[]</c>, 8 bytes per block: 0.25 bytes per entry; a text whose first level does not fit an array
/// must have the RMQ section). The LCP entries are read through <see cref="CompactLcp"/> (its escape is resolved there). Every loop runs inside a block or a level group, so a query reads at most F x (levels + 1) entries, and a pyramid that
/// disagrees with the array it describes ends in <see cref="InvalidDataException"/> (a container) or <see cref="InvalidOperationException"/> (memory), never in a read outside the arrays. Immutable after construction; safe for concurrent readers.
/// This is the form for every text that the 32-bit <see cref="LcpPyramid"/> of S1 cannot hold (and the form that a test builds over generated sections).
/// </summary>
internal sealed class WideLcpPyramid : ILcpRanges
{
    private const int MaxLevels = 64;

    private readonly CompactLcp _lcp;
    private readonly bool _mapped;
    private readonly long _length;
    private readonly int _fanout;
    private readonly int _logFanout;
    private readonly int _levelCount;
    private readonly long[] _offset = new long[MaxLevels];
    private readonly long[] _count = new long[MaxLevels];
    private readonly long[]? _built;
    private readonly PosSource? _persisted;

    private WideLcpPyramid(CompactLcp lcp, int fanout, PosSource? persisted)
    {
        _lcp = lcp;
        _mapped = lcp.IsMapped;
        _length = lcp.Length;
        _fanout = fanout;
        _logFanout = int.TrailingZeroCount(fanout);
        long total = 0;
        int levels = 0;
        if (_length > 1)
        {
            long size = _length;
            do
            {
                size = (size + fanout - 1) / fanout;
                if (levels == MaxLevels) throw new InvalidOperationException("more pyramid levels than the table holds");
                _offset[levels] = total;
                _count[levels] = size;
                total += size;
                levels++;
            }
            while (size > 1);
        }
        _levelCount = levels;
        if (persisted is not null)
        {
            if (persisted.Count != total) throw new InvalidDataException($"the RMQ section has {persisted.Count} entries, the pyramid of {_length} LCP entries at fanout {fanout} has {total}: the file is damaged");
            _persisted = persisted;
            return;
        }
        if (total > Array.MaxLength) throw new NotSupportedException($"the pyramid of an LCP array of {_length} entries has {total} block minima, more than an array of memory holds ({Array.MaxLength}): the container needs its RMQ section (write it with the builder option)");
        _built = new long[total];
        Fill();
    }

    /// <summary>The pyramid of <paramref name="lcp"/> (n + 1 entries) with its levels built in memory.</summary>
    /// <exception cref="ArgumentOutOfRangeException">The fanout is not a power of two of at least 2.</exception>
    /// <exception cref="NotSupportedException">The levels do not fit an array.</exception>
    internal static WideLcpPyramid Build(CompactLcp lcp, int fanout = LcpPyramid.DefaultFanout)
    {
        ArgumentNullException.ThrowIfNull(lcp);
        LcpPyramid.CheckFanout(fanout);
        return new WideLcpPyramid(lcp, fanout, null);
    }

    /// <summary>The pyramid of <paramref name="lcp"/> whose levels are the entries of <paramref name="levels"/> (the RMQ section, level 0 first): nothing is built or copied.</summary>
    /// <exception cref="InvalidDataException">The section does not hold the levels of this array at this fanout.</exception>
    internal static WideLcpPyramid FromLevels(CompactLcp lcp, int fanout, PosSource levels)
    {
        ArgumentNullException.ThrowIfNull(lcp);
        ArgumentNullException.ThrowIfNull(levels);
        LcpPyramid.CheckFanout(fanout);
        return new WideLcpPyramid(lcp, fanout, levels);
    }

    /// <inheritdoc />
    public long Length => _length;

    /// <inheritdoc />
    public int Fanout => _fanout;

    /// <summary>Number of pyramid levels above the LCP array.</summary>
    internal int Levels => _levelCount;

    /// <inheritdoc />
    public long LevelBytes => _built is null ? 0 : 8L * _built.Length;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private long Level(long index) => _built is not null ? _built[index] : _persisted![index];

    private void Fill()
    {
        if (_levelCount == 0) return;
        long f = _fanout;
        for (long b = 0, blocks = _count[0]; b < blocks; b++)
        {
            long start = b * f;
            long end = Math.Min(_length, start + f);
            long min = long.MaxValue;
            for (long j = start; j < end; j++) min = Math.Min(min, _lcp[j]);
            _built![b] = min;
        }
        for (int k = 1; k < _levelCount; k++)
        {
            long below = _offset[k - 1], belowCount = _count[k - 1], here = _offset[k];
            for (long b = 0, blocks = _count[k]; b < blocks; b++)
            {
                long start = b * f;
                long end = Math.Min(belowCount, start + f);
                long min = long.MaxValue;
                for (long j = start; j < end; j++) min = Math.Min(min, _built![below + j]);
                _built![here + b] = min;
            }
        }
    }

    /// <inheritdoc />
    public long ValueAt(long rank)
    {
        if ((ulong)rank >= (ulong)_length) throw new ArgumentOutOfRangeException(nameof(rank));
        return _lcp[rank];
    }

    /// <inheritdoc />
    public long Min(long lo, long hi)
    {
        if ((ulong)lo > (ulong)hi || hi >= _length) throw new ArgumentOutOfRangeException(nameof(lo), $"[{lo}, {hi}] is not a range of {_length} entries");
        long f = _fanout;
        long min = long.MaxValue;
        long loBlock = lo >> _logFanout, hiBlock = hi >> _logFanout;
        if (loBlock == hiBlock)
        {
            for (long j = lo; j <= hi; j++) min = Math.Min(min, _lcp[j]);
            return min;
        }

        // the partial blocks at both ends are read from the LCP array itself
        for (long j = lo, end = lo | (f - 1); j <= end; j++) min = Math.Min(min, _lcp[j]);
        for (long j = hi & ~(f - 1); j <= hi; j++) min = Math.Min(min, _lcp[j]);

        // the whole blocks in between, level by level: both ends contribute their partial group, then the range moves one level up
        long l = loBlock + 1, r = hiBlock - 1;
        int k = 1;
        while (l <= r)
        {
            long off = _offset[k - 1];
            if ((l >> _logFanout) == (r >> _logFanout))
            {
                for (long b = l; b <= r; b++) min = Math.Min(min, Level(off + b));
                break;
            }
            for (long b = l, end = l | (f - 1); b <= end; b++) min = Math.Min(min, Level(off + b));
            for (long b = r & ~(f - 1); b <= r; b++) min = Math.Min(min, Level(off + b));
            l = (l >> _logFanout) + 1;
            r = (r >> _logFanout) - 1;
            k++;
        }
        return min;
    }

    /// <inheritdoc />
    public long PrevLess(long pos, long bound)
    {
        if (pos < 0) return -1;
        ArgumentOutOfRangeException.ThrowIfGreaterThanOrEqual(pos, _length);
        long f = _fanout;
        for (long j = pos, start = pos & ~(f - 1); j >= start; j--)
            if (_lcp[j] < bound) return j;

        long idx = pos >> _logFanout;
        for (int k = 1; k <= _levelCount; k++)
        {
            long off = _offset[k - 1];
            for (long b = idx - 1, start = idx & ~(f - 1); b >= start; b--)
                if (Level(off + b) < bound) return RightmostBelow(k, b, bound);
            idx >>= _logFanout;
        }
        return -1;
    }

    /// <inheritdoc />
    public long NextLess(long pos, long bound)
    {
        if (pos >= _length) return _length;
        ArgumentOutOfRangeException.ThrowIfNegative(pos);
        long f = _fanout;
        for (long j = pos, end = Math.Min(pos | (f - 1), _length - 1); j <= end; j++)
            if (_lcp[j] < bound) return j;

        long idx = pos >> _logFanout;
        for (int k = 1; k <= _levelCount; k++)
        {
            long off = _offset[k - 1];
            for (long b = idx + 1, end = Math.Min(idx | (f - 1), _count[k - 1] - 1); b <= end; b++)
                if (Level(off + b) < bound) return LeftmostBelow(k, b, bound);
            idx >>= _logFanout;
        }
        return _length;
    }

    /// <inheritdoc />
    public void Borders(long a, long b, long t, out long left, out long right, out long leftValue, out long rightValue)
    {
        if ((ulong)a >= (ulong)_length || b < a || b >= _length) throw new ArgumentOutOfRangeException(nameof(a));
        left = PrevLess(a, t);
        right = NextLess(b + 1, t);
        leftValue = left >= 0 ? _lcp[left] : -1;
        rightValue = right < _length ? _lcp[right] : -1;
    }

    /// <summary>Previous smaller value (strict): the largest <c>j &lt; p</c> with <c>LCP[j] &lt; LCP[p]</c>, or -1.</summary>
    internal long PreviousSmaller(long p) => PrevLess(p - 1, ValueAt(p));

    /// <summary>Next smaller value (strict): the smallest <c>j &gt; p</c> with <c>LCP[j] &lt; LCP[p]</c>, or <see cref="Length"/>.</summary>
    internal long NextSmaller(long p) => NextLess(p + 1, ValueAt(p));

    private long RightmostBelow(int k, long b, long bound)
    {
        long f = _fanout;
        while (k > 1)
        {
            k--;
            long off = _offset[k - 1];
            long lowest = b * f;
            long c = Math.Min(lowest + f - 1, _count[k - 1] - 1);
            while (Level(off + c) >= bound)
            {
                if (--c < lowest) throw Inconsistent();
            }
            b = c;
        }
        long first = b * f;
        for (long position = Math.Min(first + f - 1, _length - 1); position >= first; position--)
            if (_lcp[position] < bound) return position;
        throw Inconsistent();
    }

    private long LeftmostBelow(int k, long b, long bound)
    {
        long f = _fanout;
        while (k > 1)
        {
            k--;
            long off = _offset[k - 1];
            long c = b * f;
            long last = Math.Min(c + f - 1, _count[k - 1] - 1);
            while (Level(off + c) >= bound)
            {
                if (++c > last) throw Inconsistent();
            }
            b = c;
        }
        long end = Math.Min((b * f) + f - 1, _length - 1);
        for (long position = b * f; position <= end; position++)
            if (_lcp[position] < bound) return position;
        throw Inconsistent();
    }

    /// <summary>The pyramid and the entries it was built from disagree (a block whose minimum is below the bound holds no entry below it): never for a pyramid built in memory over an array that does not change, a damaged section for one read in place.</summary>
    private Exception Inconsistent()
    {
        const string message = "the LCP pyramid is inconsistent: no entry below the bound in a block whose minimum is below it";
        return _mapped ? new InvalidDataException(message + ": the file is damaged") : new InvalidOperationException(message);
    }

    /// <summary>Releases nothing: the levels are an array of memory or a section of the container, which another owner holds.</summary>
    public void Dispose()
    {
    }
}
