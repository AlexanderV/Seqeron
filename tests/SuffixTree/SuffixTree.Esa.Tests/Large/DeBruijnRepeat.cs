using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.Large;

/// <summary>
/// The closed form of <c>u^m</c> (m >= 2) for the de Bruijn word <c>u</c> of order d over sigma letters (work package S2.1, ADR-005 end of section 3.7): fan-out sigma at every level down to depth d, then the periodic tail.
/// With <c>P = sigma^d</c> every word of length d is the cyclic window <c>V(r)</c> of exactly one rotation r of u, so the suffixes of length at least d are ordered by the window value (a digit string of base sigma) and, within a
/// rotation, by length (the shorter is a prefix of the longer): the suffix array is a table lookup in <c>V</c> and its inverse (one period, <c>P</c> entries, computed once from the word of
/// <see cref="StringFamilies.DeBruijn"/>) and a division. The <c>d - 1</c> suffixes shorter than d (lengths 1 .. d - 1, the last symbols of u) are ordered by the end of the text, not by the window: each is placed by counting
/// (the suffixes of the windows below its prefix, plus the shorter ones that are smaller) and the rank maps skip and re-insert them. Two suffixes of different rotations have the LCP of their windows, bounded by their lengths.
/// There is no closed form of the child table (it is the scan of the engine over the generated LCP): <see cref="HasChildTable"/> is false.
/// </summary>
internal sealed class DeBruijnRepeat : ClosedFormText
{
    private readonly int _sigma;
    private readonly int _order;
    private readonly long _period;
    private readonly long _m;
    private readonly byte[] _u;
    private readonly int[] _windowOfRotation;
    private readonly int[] _rotationOfWindow;
    private readonly long[] _specialRank;   // [l - 1]: the rank of the suffix of length l (1 <= l < d)
    private readonly long[] _specialPure;   // [l - 1]: the index of that suffix in the order of the blocks (m V(P - l))
    private readonly long[] _sortedRanks;
    private readonly long[] _sortedPure;

    internal DeBruijnRepeat(int sigma, int order, long m)
        : base(checked(Pow(sigma, order) * m))
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(sigma, 2);
        ArgumentOutOfRangeException.ThrowIfGreaterThan(sigma, 26);
        ArgumentOutOfRangeException.ThrowIfLessThan(order, 1);
        ArgumentOutOfRangeException.ThrowIfLessThan(m, 2);
        _sigma = sigma;
        _order = order;
        _m = m;
        _period = Pow(sigma, order);
        (_u, _windowOfRotation, _rotationOfWindow) = TablesOf(sigma, order);

        int specials = order - 1;
        _specialRank = new long[specials];
        _specialPure = new long[specials];
        for (int l = 1; l <= specials; l++)
        {
            long value = PrefixValue(_period - l, l); // the suffix of length l is the first l symbols of the rotation P - l
            long below = value * Pow(sigma, order - l);
            long smallerLong = (_m * below) - Enumerable.Range(1, specials).Count(other => _windowOfRotation[_period - other] < below);
            long smallerShort = Enumerable.Range(1, specials).Count(other => other != l && CompareSuffixesOfU(other, l) < 0);
            _specialRank[l - 1] = 1 + smallerLong + smallerShort;
            _specialPure[l - 1] = _m * _windowOfRotation[_period - l];
        }
        _sortedRanks = [.. _specialRank.Order()];
        _sortedPure = [.. _specialPure.Order()];
    }

    /// <summary>The order d of the word (the length of the windows that tell the rotations apart).</summary>
    internal int Order => _order;

    internal override bool HasChildTable => false;

    internal override string Name => $"u^{_m} (de Bruijn {_sigma}^{_order})";

    internal override long Period => _period;

    internal override int Sigma => _sigma;

    internal override long SymbolCount(int code) => code >= 0 && code < _sigma ? _m * (_period / _sigma) : 0;

    internal override int InfiniteCodeAt(long position) => _u[position % _period];

    private static readonly System.Collections.Concurrent.ConcurrentDictionary<(int Sigma, int Order), (byte[] U, int[] Window, int[] Rotation)> Cache = new();

    /// <summary>The word u of <see cref="StringFamilies.DeBruijn"/>, the value of the cyclic window of each rotation and its inverse (one period, shared by every copy count).</summary>
    private static (byte[] U, int[] Window, int[] Rotation) TablesOf(int sigma, int order) => Cache.GetOrAdd((sigma, order), key =>
    {
        long period = Pow(key.Sigma, key.Order);
        string word = StringFamilies.DeBruijn(key.Sigma, key.Order);
        if (word.Length != period) throw new InvalidOperationException("the de Bruijn word has the wrong length");
        var u = new byte[period];
        for (int i = 0; i < period; i++) u[i] = (byte)(word[i] - 'a');
        var windowOf = new int[period];
        var rotationOf = new int[period];
        long window = 0;
        for (int j = 0; j < key.Order; j++) window = (window * key.Sigma) + u[j];
        long high = period / key.Sigma;
        for (int r = 0; r < period; r++)
        {
            windowOf[r] = (int)window;
            rotationOf[window] = r;
            window = ((window % high) * key.Sigma) + u[(r + key.Order) % period];
        }
        return (u, windowOf, rotationOf);
    });

    private static long Pow(int b, int e)
    {
        long r = 1;
        for (int i = 0; i < e; i++) r = checked(r * b);
        return r;
    }

    /// <summary>The value (base sigma) of the first <paramref name="length"/> symbols of the rotation <paramref name="rotation"/> of u (cyclic).</summary>
    private long PrefixValue(long rotation, int length)
    {
        long value = 0;
        for (int j = 0; j < length; j++) value = (value * _sigma) + _u[(rotation + j) % _period];
        return value;
    }

    /// <summary>Compares the suffixes of u of lengths a and b (a proper prefix is smaller).</summary>
    private int CompareSuffixesOfU(int a, int b) => new ReadOnlySpan<byte>(_u, (int)_period - a, a).SequenceCompareTo(new ReadOnlySpan<byte>(_u, (int)_period - b, b));

    /// <summary>The number of equal leading digits of two windows (as strings of d digits of base sigma).</summary>
    private int WindowLcp(long rotationA, long rotationB)
    {
        int a = _windowOfRotation[rotationA], b = _windowOfRotation[rotationB];
        long divisor = _period / _sigma;
        int common = 0;
        while (common < _order && (a / divisor) % _sigma == (b / divisor) % _sigma)
        {
            common++;
            divisor /= _sigma;
        }
        return common;
    }

    private bool IsSpecialPosition(long position) => position < Length && Length - position < _order;

    internal override long SuffixAt(long rank)
    {
        CheckRank(rank);
        if (rank == 0) return Length;
        int special = Array.IndexOf(_specialRank, rank);
        if (special >= 0) return Length - (special + 1);
        long k = rank - 1;
        foreach (long s in _sortedRanks)
            if (s < rank) k--;
        long pure = k;
        foreach (long q in _sortedPure)
            if (q <= pure) pure++;
        long window = pure / _m, j = pure % _m;
        return ((_m - 1 - j) * _period) + _rotationOfWindow[window];
    }

    internal override long RankOf(long position)
    {
        CheckRank(position, nameof(position));
        if (position == Length) return 0;
        if (IsSpecialPosition(position)) return _specialRank[(int)(Length - position) - 1];
        long c = position / _period, r = position % _period;
        long pure = (_m * _windowOfRotation[r]) + (_m - 1 - c);
        long k = pure;
        foreach (long q in _sortedPure)
            if (q < pure) k--;
        long rank = 1 + k;
        foreach (long s in _sortedRanks)
            if (s <= rank) rank++;
        return rank;
    }

    /// <summary>The LCP of two different suffixes: within a rotation the shorter one; across rotations the windows' LCP, bounded by both lengths.</summary>
    private long PairLcp(long a, long b)
    {
        long la = Length - a, lb = Length - b;
        if (a == b) return la;
        if ((a - b) % _period == 0) return Math.Min(la, lb);
        return Math.Min(Math.Min(la, lb), WindowLcp(a % _period, b % _period));
    }

    internal override long Lcp(long rank)
    {
        CheckRank(rank);
        return rank <= 1 ? 0 : PairLcp(SuffixAt(rank - 1), SuffixAt(rank));
    }

    internal override Int128 LcpSum
    {
        get
        {
            if (Length > 1L << 26) throw new NotSupportedException("the sum of the LCP of u^m is a pass over the ranks: refused above 2^26 symbols");
            Int128 sum = 0;
            for (long r = 2; r <= Length; r++) sum += Lcp(r);
            return sum;
        }
    }

    /// <summary>The non-empty right-special factors: every word shorter than d (sigma^1 + ... + sigma^(d-1)) and the suffixes of the text of length d .. n - P (they occur again one period earlier).</summary>
    internal override long InternalCount
    {
        get
        {
            long words = 0;
            for (int l = 1; l < _order; l++) words += Pow(_sigma, l);
            return words + (Length - _period - _order + 1);
        }
    }

    internal override long LrsDepth => Length - _period;

    internal override long LrsRank => RankOf(_period);

    internal override RankRange RangeOf(ReadOnlySpan<char> pattern)
    {
        if (pattern.IsEmpty) return new RankRange(1, Length);
        var digits = new byte[pattern.Length];
        for (int i = 0; i < pattern.Length; i++)
        {
            int d = pattern[i] - 'a';
            if (d < 0 || d >= _sigma) return RankRange.Empty;
            digits[i] = (byte)d;
        }
        int len = digits.Length;
        if (len >= _order)
        {
            long window = 0;
            for (int j = 0; j < _order; j++) window = (window * _sigma) + digits[j];
            long r = _rotationOfWindow[window];
            for (int x = _order; x < len; x++)
                if (digits[x] != _u[(r + x) % _period]) return RankRange.Empty;
            return InRotation(r, len);
        }
        long value = 0;
        for (int j = 0; j < len; j++) value = (value * _sigma) + digits[j];
        long span = Pow(_sigma, _order - len);
        long wrapped = 0;
        for (long r = _period - len + 1; r < _period; r++)
            if (PrefixValue(r, len) == value) wrapped++;
        long count = (_m * span) - wrapped;
        long top = (value * span) + span - 1;
        long last = RankOf(_rotationOfWindow[top]);
        return new RankRange(last - count + 1, last);
    }

    /// <summary>The occurrences of a pattern of <paramref name="length"/> >= d symbols that starts at the rotation <paramref name="rotation"/>: its positions that leave room, the shortest suffixes of the rotation's block.</summary>
    private RankRange InRotation(long rotation, long length)
    {
        long room = Length - length - rotation;
        if (room < 0) return RankRange.Empty;
        long cMax = room / _period;
        return new RankRange(RankOf((cMax * _period) + rotation), RankOf(rotation));
    }

    internal override RankRange RangeOfRun(long length)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(length);
        if (length == 0) return new RankRange(1, Length);
        return length > _order ? RankRange.Empty : RangeOf(new string('a', (int)length));
    }

    internal override RankRange RangeOfFactor(long start, long length)
    {
        if (start < 0 || length < 0 || start > Length - length) throw new ArgumentOutOfRangeException(nameof(start), start, $"[{start}, {start} + {length}) is not a factor of a text of {Length} symbols");
        if (length == 0) return new RankRange(1, Length);
        if (length < _order) return RangeOf(new string(Enumerable.Range(0, (int)length).Select(i => CharAt(start + i)).ToArray()));
        return InRotation(start % _period, length);
    }

    internal override long Down(long i) => throw NoChildTable();

    internal override long Up(long j) => throw NoChildTable();

    internal override long Next(long k) => throw NoChildTable();

    private static NotSupportedException NoChildTable() => new("the child table of u^m has no closed form: it is the scan of the engine over the generated LCP (ADR-005 section 3.7)");
}
