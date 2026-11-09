using System.Runtime.CompilerServices;

namespace SuffixTree.Esa.Bench.Links;

/// <summary>
/// The RMQ + PSV/NSV helper of link kernel (c) (ADR-001 section 5, prototype P1): a pyramid of block minima over the LCP array (fanout F, a power of two; level k holds the
/// minimum of F^k consecutive entries) answers range-minimum queries with the leftmost minimum position and the previous / next strictly smaller value in O(F * levels)
/// sequential reads, with about 4 / (F - 1) bytes per entry. The LCP entries themselves are read through the <see cref="CompactLcp"/> of the index (no copy).
/// </summary>
internal sealed class LcpRangeHelper
{
    private readonly CompactLcp _lcp;
    private readonly int _fanout;
    private readonly int _logFanout;
    private readonly int _length;

    /// <summary><c>_levels[k - 1]</c> is level k: entry <c>b</c> is the minimum of the F entries <c>b*F .. b*F + F - 1</c> of the level below (the LCP array for k = 1).</summary>
    private readonly int[][] _levels;

    /// <param name="lcp">The LCP array (n + 1 entries).</param>
    /// <param name="fanout">F, a power of two of at least 2.</param>
    internal LcpRangeHelper(CompactLcp lcp, int fanout = 16)
    {
        ArgumentNullException.ThrowIfNull(lcp);
        if (fanout < 2 || (fanout & (fanout - 1)) != 0) throw new ArgumentOutOfRangeException(nameof(fanout), fanout, "the fanout must be a power of two of at least 2");
        _lcp = lcp;
        _fanout = fanout;
        _logFanout = int.TrailingZeroCount(fanout);
        _length = lcp.Length;

        var levels = new List<int[]>();
        var cursor = lcp.CreateCursor();
        int size = (_length + fanout - 1) / fanout;
        if (_length > 1)
        {
            var first = new int[size];
            for (int b = 0; b < size; b++)
            {
                int end = Math.Min(_length, (b + 1) * fanout);
                int min = int.MaxValue;
                for (int j = b * fanout; j < end; j++) min = Math.Min(min, cursor.Get(j));
                first[b] = min;
            }
            levels.Add(first);
            while (levels[^1].Length > 1)
            {
                int[] below = levels[^1];
                var next = new int[(below.Length + fanout - 1) / fanout];
                for (int b = 0; b < next.Length; b++)
                {
                    int end = Math.Min(below.Length, (b + 1) * fanout);
                    int min = int.MaxValue;
                    for (int j = b * fanout; j < end; j++) min = Math.Min(min, below[j]);
                    next[b] = min;
                }
                levels.Add(next);
            }
        }
        _levels = [.. levels];
    }

    /// <summary>Number of LCP entries (n + 1).</summary>
    internal int Length => _length;

    /// <summary>The fanout F.</summary>
    internal int Fanout => _fanout;

    /// <summary>Number of pyramid levels above the LCP array.</summary>
    internal int Levels => _levels.Length;

    /// <summary>Bytes of the pyramid (the LCP array is not counted).</summary>
    internal long PayloadBytes
    {
        get
        {
            long bytes = 0;
            foreach (int[] level in _levels) bytes += 4L * level.Length;
            return bytes;
        }
    }

    /// <summary>Minimum of <c>LCP[lo..hi]</c> (inclusive).</summary>
    internal int Min(int lo, int hi)
    {
        var none = new NoCount();
        return Min(lo, hi, ref none);
    }

    /// <summary>Leftmost position of the minimum of <c>LCP[lo..hi]</c> (inclusive).</summary>
    internal int ArgMin(int lo, int hi)
    {
        var none = new NoCount();
        return ArgMin(lo, hi, ref none);
    }

    /// <summary>The largest <c>j &lt;= pos</c> with <c>LCP[j] &lt; bound</c>, or -1.</summary>
    internal int PrevLess(int pos, int bound)
    {
        var none = new NoCount();
        return PrevLess(pos, bound, ref none);
    }

    /// <summary>The smallest <c>j &gt;= pos</c> with <c>LCP[j] &lt; bound</c>, or <see cref="Length"/>.</summary>
    internal int NextLess(int pos, int bound)
    {
        var none = new NoCount();
        return NextLess(pos, bound, ref none);
    }

    /// <summary>Previous smaller value (strict): the largest <c>j &lt; p</c> with <c>LCP[j] &lt; LCP[p]</c>, or -1.</summary>
    internal int PreviousSmaller(int p) => PrevLess(p - 1, _lcp[p]);

    /// <summary>Next smaller value (strict): the smallest <c>j &gt; p</c> with <c>LCP[j] &lt; LCP[p]</c>, or <see cref="Length"/>.</summary>
    internal int NextSmaller(int p) => NextLess(p + 1, _lcp[p]);

    /// <summary>The minimum with step accounting (one step per entry read).</summary>
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    internal int Min<TCnt>(int lo, int hi, ref TCnt cnt)
        where TCnt : struct, IStepCounter
    {
        if ((uint)lo > (uint)hi || hi >= _length) throw new ArgumentOutOfRangeException(nameof(lo), $"[{lo}, {hi}] is not a range of {_length} entries");
        int f = _fanout;
        int min = int.MaxValue;
        int loBlock = lo >> _logFanout, hiBlock = hi >> _logFanout;
        if (loBlock == hiBlock)
        {
            for (int j = lo; j <= hi; j++)
            {
                cnt.Add(1);
                min = Math.Min(min, _lcp[j]);
            }
            return min;
        }

        // the partial blocks at both ends are read from the LCP array itself
        for (int j = lo, end = lo | (f - 1); j <= end; j++)
        {
            cnt.Add(1);
            min = Math.Min(min, _lcp[j]);
        }
        for (int j = hi & ~(f - 1); j <= hi; j++)
        {
            cnt.Add(1);
            min = Math.Min(min, _lcp[j]);
        }

        // the whole blocks in between, level by level: both ends contribute their partial group, then the range moves one level up
        int l = loBlock + 1, r = hiBlock - 1;
        int k = 1;
        while (l <= r)
        {
            int[] level = _levels[k - 1];
            if ((l >> _logFanout) == (r >> _logFanout))
            {
                for (int b = l; b <= r; b++)
                {
                    cnt.Add(1);
                    min = Math.Min(min, level[b]);
                }
                break;
            }
            for (int b = l, end = l | (f - 1); b <= end; b++)
            {
                cnt.Add(1);
                min = Math.Min(min, level[b]);
            }
            for (int b = r & ~(f - 1); b <= r; b++)
            {
                cnt.Add(1);
                min = Math.Min(min, level[b]);
            }
            l = (l >> _logFanout) + 1;
            r = (r >> _logFanout) - 1;
            k++;
        }
        return min;
    }

    /// <summary><see cref="ArgMin(int, int)"/> with step accounting: the leftmost position holding the minimum is the first position at or after <c>lo</c> that is not above it.</summary>
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    internal int ArgMin<TCnt>(int lo, int hi, ref TCnt cnt)
        where TCnt : struct, IStepCounter
    {
        int min = Min(lo, hi, ref cnt);
        return NextLess(lo, min + 1, ref cnt);
    }

    /// <summary><see cref="PrevLess(int, int)"/> with step accounting.</summary>
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    internal int PrevLess<TCnt>(int pos, int bound, ref TCnt cnt)
        where TCnt : struct, IStepCounter
    {
        if (pos < 0) return -1;
        int f = _fanout;
        for (int j = pos, start = pos & ~(f - 1); j >= start; j--)
        {
            cnt.Add(1);
            if (_lcp[j] < bound) return j;
        }

        int idx = pos >> _logFanout;
        for (int k = 1; k <= _levels.Length; k++)
        {
            int[] level = _levels[k - 1];
            for (int b = idx - 1, start = idx & ~(f - 1); b >= start; b--)
            {
                cnt.Add(1);
                if (level[b] < bound) return RightmostBelow(k, b, bound, ref cnt);
            }
            idx >>= _logFanout;
        }
        return -1;
    }

    /// <summary><see cref="NextLess(int, int)"/> with step accounting.</summary>
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    internal int NextLess<TCnt>(int pos, int bound, ref TCnt cnt)
        where TCnt : struct, IStepCounter
    {
        if (pos >= _length) return _length;
        int f = _fanout;
        for (int j = pos, end = Math.Min(pos | (f - 1), _length - 1); j <= end; j++)
        {
            cnt.Add(1);
            if (_lcp[j] < bound) return j;
        }

        int idx = pos >> _logFanout;
        for (int k = 1; k <= _levels.Length; k++)
        {
            int[] level = _levels[k - 1];
            for (int b = idx + 1, end = Math.Min(idx | (f - 1), level.Length - 1); b <= end; b++)
            {
                cnt.Add(1);
                if (level[b] < bound) return LeftmostBelow(k, b, bound, ref cnt);
            }
            idx >>= _logFanout;
        }
        return _length;
    }

    /// <summary>The rightmost position below <paramref name="bound"/> inside block <paramref name="b"/> of level <paramref name="k"/>, whose minimum is known to be below it.</summary>
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    private int RightmostBelow<TCnt>(int k, int b, int bound, ref TCnt cnt)
        where TCnt : struct, IStepCounter
    {
        int f = _fanout;
        while (k > 1)
        {
            k--;
            int[] level = _levels[k - 1];
            int c = Math.Min((b * f) + f - 1, level.Length - 1);
            while (true)
            {
                cnt.Add(1);
                if (level[c] < bound) break;
                c--;
            }
            b = c;
        }
        int position = Math.Min((b * f) + f - 1, _length - 1);
        while (true)
        {
            cnt.Add(1);
            if (_lcp[position] < bound) return position;
            position--;
        }
    }

    /// <summary>The leftmost position below <paramref name="bound"/> inside block <paramref name="b"/> of level <paramref name="k"/>, whose minimum is known to be below it.</summary>
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    private int LeftmostBelow<TCnt>(int k, int b, int bound, ref TCnt cnt)
        where TCnt : struct, IStepCounter
    {
        int f = _fanout;
        while (k > 1)
        {
            k--;
            int[] level = _levels[k - 1];
            int c = b * f;
            while (true)
            {
                cnt.Add(1);
                if (level[c] < bound) break;
                c++;
            }
            b = c;
        }
        int position = b * f;
        while (true)
        {
            cnt.Add(1);
            if (_lcp[position] < bound) return position;
            position++;
        }
    }
}
