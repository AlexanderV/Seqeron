namespace SuffixTree.Esa.Tests.Large;

/// <summary>
/// The closed form of the prefix of length n of <c>(a^k b)^infinity</c> for any n and any k >= 1 (work package S2.1, ADR-005 section 3.7; <c>(ab)^m</c> is k = 1, n = 2m). With the period <c>p = k + 1</c>, <c>n = q p + t</c>
/// (0 &lt;= t &lt;= k) the text is q full periods and a last partial run <c>a^t</c>. The suffixes sort as follows (an end of text is smaller than every symbol; a suffix is a prefix of every earlier suffix of its phase):
/// <list type="bullet">
/// <item>rank 0 the empty suffix; ranks 1..t the pure runs <c>a^s</c> in the last period (shortest first);</item>
/// <item>then k + 1 blocks of q ranks, block <c>phase</c> = 0..k holding the suffixes that start at <c>phase</c> in a full period (a leading run of k - phase letters a, then b; block k starts with b), the later position first:
/// <c>Rank(phase, j) = 1 + t + phase q + j</c> is the suffix at <c>(q - 1 - j) p + phase</c>.</item>
/// </list>
/// So <c>LCP(phase, 0) = k - phase</c> (the first of a block), <c>LCP(phase, j) = t + j p - phase</c> (j >= 1: the shorter suffix is a prefix), <c>LCP(0, 0) = t</c>, and the pure runs have <c>LCP = s - 1</c>.
/// Internal nodes: the chains <c>(phase, j)</c> for j >= 1 (k + 1 of length q - 1, each interval from the rank before to the end of its block, nested) and, for each value v from 1 to max(t, k - 1), one interval I_v of value v whose
/// l-indices are the pure rank v + 1 (when v &lt;= t) and the first rank of block k - v. The child table follows from the nesting (see <see cref="Down"/>, <see cref="Up"/>, <see cref="Next"/>); when q = 0 the text is <c>a^t</c>.
/// </summary>
internal sealed class PeriodicAkB : ClosedFormText
{
    private readonly long _k;
    private readonly long _p;
    private readonly long _q;
    private readonly long _t;
    private readonly PowerOfA? _pure;

    internal PeriodicAkB(long n, long k)
        : base(n)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(k, 1);
        _k = k;
        _p = k + 1;
        _q = n / _p;
        _t = n % _p;
        if (_q == 0) _pure = new PowerOfA(n);
    }

    /// <summary><c>(ab)^m</c>: the case k = 1 with n = 2m.</summary>
    internal static PeriodicAkB AbRepeat(long m) => new(2 * m, 1);

    internal override string Name => _k == 1 && Length % 2 == 0 ? $"(ab)^{Length / 2}" : $"(a^{_k} b)* n={Length}";

    internal override long Period => _p;

    internal override int Sigma => Length == 0 ? 0 : _q == 0 ? 1 : 2;

    internal override long SymbolCount(int code) => code switch { 0 => Length - _q, 1 => _q, _ => 0 };

    internal override int InfiniteCodeAt(long position) => position % _p == _k ? 1 : 0;

    /// <summary>The rank of the suffix of block <paramref name="phase"/> that is the j-th (later positions first).</summary>
    private long Rank(long phase, long j) => 1 + _t + (phase * _q) + j;

    internal override long SuffixAt(long rank)
    {
        CheckRank(rank);
        if (rank == 0) return Length;
        if (rank <= _t) return Length - rank; // the pure runs (all ranks when q = 0)
        long x = rank - 1 - _t;
        long phase = x / _q, j = x % _q;
        return ((_q - 1 - j) * _p) + phase;
    }

    internal override long RankOf(long position)
    {
        CheckRank(position, nameof(position));
        if (position == Length) return 0;
        if (position >= _q * _p) return Length - position; // the last partial period: the run a^s with s = n - position
        long phase = position % _p, c = position / _p;
        return Rank(phase, _q - 1 - c);
    }

    internal override long Lcp(long rank)
    {
        CheckRank(rank);
        if (rank <= 1) return 0;
        if (rank <= _t) return rank - 1;
        if (rank == _t + 1) return _t;
        long x = rank - 1 - _t;
        long phase = x / _q, j = x % _q;
        return j == 0 ? _k - phase : _t + (j * _p) - phase;
    }

    internal override Int128 LcpSum
    {
        get
        {
            if (_q == 0) return _pure!.LcpSum;
            Int128 t = _t, q = _q, k = _k, p = _p;
            return (t * (t - 1) / 2) + (q * t) + (p * p * q * (q - 1) / 2) + (k * (k - 1) / 2) + ((q - 1) * k * t) - ((q - 1) * k * (k + 1) / 2);
        }
    }

    internal override long InternalCount => _q == 0 ? _pure!.InternalCount : (_p * (_q - 1)) + Math.Max(_t, _k - 1);

    internal override long LrsDepth
    {
        get
        {
            if (_q == 0) return _pure!.LrsDepth;
            if (_q >= 2) return Length - _p;
            long v = Math.Max(_t, _k - 1);
            return v >= 1 ? v : 0;
        }
    }

    internal override long LrsRank
    {
        get
        {
            if (_q == 0) return _pure!.LrsRank;
            if (_q >= 2) return _t + _q - 1; // the rank before the last element of block 0
            long v = LrsDepth;
            if (v == 0) return -1;
            return v <= _t ? v : _t + 1;
        }
    }

    internal override RankRange RangeOf(ReadOnlySpan<char> pattern)
    {
        if (pattern.IsEmpty) return RangeOfRun(0);
        int run = 0;
        while (run < pattern.Length && pattern[run] == 'a') run++;
        if (run == pattern.Length) return RangeOfRun(run);
        if (_pure is not null || run > _k) return RankRange.Empty;
        long phase = _k - run; // the pattern starts with a^run b: its occurrences are the positions of the phase k - run
        for (int x = 0; x < pattern.Length; x++)
        {
            char expected = (phase + x) % _p == _k ? 'b' : 'a';
            if (pattern[x] != expected) return RankRange.Empty;
        }
        return Occurrences(phase, pattern.Length);
    }

    internal override RankRange RangeOfRun(long length)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(length);
        if (length == 0) return Length == 0 ? RankRange.Empty : new RankRange(1, Length);
        if (_pure is not null) return _pure.RangeOfRun(length);
        if (length > _k) return RankRange.Empty;
        // a^length starts the pure runs a^s (s >= length) and the blocks of the phases 0..k - length
        long first = length <= _t ? length : _t + 1;
        long last = _t + ((_k - length + 1) * _q);
        return new RankRange(first, last);
    }

    internal override RankRange RangeOfFactor(long start, long length)
    {
        if (start < 0 || length < 0 || start > Length - length) throw new ArgumentOutOfRangeException(nameof(start), start, $"[{start}, {start} + {length}) is not a factor of a text of {Length} symbols");
        if (length == 0) return RangeOfRun(0);
        if (_pure is not null) return _pure.RangeOfRun(length);
        long phase = start % _p;
        long leadingRun = phase == _k ? 0 : _k - phase;
        return length <= leadingRun ? RangeOfRun(length) : Occurrences(phase, length);
    }

    /// <summary>The suffixes of block <paramref name="phase"/> that leave room for a pattern of <paramref name="length"/> symbols that starts at that phase and contains a b: the earliest positions, the last elements of the block.</summary>
    private RankRange Occurrences(long phase, long length)
    {
        long room = Length - length - phase;
        if (room < 0) return RankRange.Empty;
        long cMax = room / _p;
        return new RankRange(Rank(phase, _q - 1 - cMax), Rank(phase, _q - 1));
    }

    /// <summary>The first l-index of the interval I_v of value v.</summary>
    private long FirstOfValue(long v) => v <= _t ? v + 1 : 1 + _t + ((_k - v) * _q);

    /// <summary>
    /// <c>Down[i]</c>: the root (<c>Down[0] = 1</c>); the chain intervals (phase, j) with j >= 2 are the last child of the chain (phase, j - 1), so <c>Down[Rank(phase, j - 1)] = Rank(phase, j)</c>; the chain (phase, 1) is the last child
    /// of I_(k - phase) (or of the root for phase k), and of I_k when phase 0 and t = k, so <c>Down[Rank(phase, 0)] = Rank(phase, 1)</c> there.
    /// </summary>
    internal override long Down(long i)
    {
        CheckRank(i, nameof(i));
        if (_pure is not null) return _pure.Down(i);
        if (i == 0) return 1;
        if (i <= _t) return 0;
        long x = i - 1 - _t;
        long phase = x / _q, j = x % _q;
        if (j == 0) return _q >= 2 && (phase >= 1 || _t == _k) ? i + 1 : 0;
        return j <= _q - 2 ? i + 1 : 0;
    }

    /// <summary>
    /// <c>Up[j]</c>: I_v ends at <c>Rank(k - v, q - 1) = t + (k - v + 1) q</c> and is never the last child of its parent I_(v-1) (or of the root), so <c>Up[t + (k - v + 1) q] = FirstOfValue(v)</c>; the chain (0, 1) ends at
    /// <c>t + q</c> and is not the last child when t &lt; k (its parent I_t ends later; the root for t = 0).
    /// </summary>
    internal override long Up(long j)
    {
        CheckRank(j, nameof(j));
        if (_pure is not null) return 0;
        if (j > _t && (j - _t) % _q == 0)
        {
            long v = _k + 1 - ((j - _t) / _q);
            if (v >= 1 && v <= Math.Max(_t, _k - 1)) return FirstOfValue(v);
        }
        return _t < _k && _q >= 2 && j == _t + _q ? _t + 2 : 0;
    }

    /// <summary><c>Next[k]</c>: the root has the l-indices 1 and <c>Rank(k, 0)</c>; I_v with v &lt;= min(t, k - 1) has the l-indices v + 1 and <c>Rank(k - v, 0)</c>.</summary>
    internal override long Next(long k)
    {
        CheckRank(k, nameof(k));
        if (_pure is not null) return 0;
        if (k == 1) return Rank(_k, 0);
        return k >= 2 && k - 1 <= Math.Min(_t, _k - 1) ? Rank(_k - (k - 1), 0) : 0;
    }
}
