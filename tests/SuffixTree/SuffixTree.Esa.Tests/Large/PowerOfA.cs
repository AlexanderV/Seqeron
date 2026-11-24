namespace SuffixTree.Esa.Tests.Large;

/// <summary>
/// The closed form of <c>a^n</c> (work package S2.1, ADR-005 section 3.7). <c>SA[r] = n - r</c>, <c>LCP[r] = r - 1</c> (0 for the first two ranks): the suffix of length r is a prefix of the one of length r + 1, so the tree is a chain
/// of n - 1 internal nodes <c>a^k</c>, each with two children (a leaf and the next node). The lcp-interval of value k is <c>[k..n]</c> with the single l-index k + 1; the root has the l-index 1. Every interval is the last child of its
/// parent, so the child table is <c>Down[i] = i + 1</c> and nothing else.
/// </summary>
internal sealed class PowerOfA : ClosedFormText
{
    internal PowerOfA(long n) : base(n)
    {
    }

    internal override string Name => $"a^{Length}";

    internal override long Period => 1;

    internal override int Sigma => Length == 0 ? 0 : 1;

    internal override long SymbolCount(int code) => code == 0 ? Length : 0;

    internal override int InfiniteCodeAt(long position) => 0;

    internal override long SuffixAt(long rank)
    {
        CheckRank(rank);
        return Length - rank;
    }

    internal override long RankOf(long position)
    {
        CheckRank(position, nameof(position));
        return Length - position;
    }

    internal override long Lcp(long rank)
    {
        CheckRank(rank);
        return rank <= 1 ? 0 : rank - 1;
    }

    /// <summary>The sum of 0, 0, 1, ..., n - 1.</summary>
    internal override Int128 LcpSum => Length == 0 ? 0 : (Int128)Length * (Length - 1) / 2;

    internal override long InternalCount => Length == 0 ? 0 : Length - 1;

    internal override long LrsDepth => Length >= 2 ? Length - 1 : 0;

    internal override long LrsRank => Length >= 2 ? Length - 1 : -1;

    internal override RankRange RangeOf(ReadOnlySpan<char> pattern)
    {
        foreach (char c in pattern)
            if (c != 'a') return RankRange.Empty;
        return RangeOfRun(pattern.Length);
    }

    /// <summary>The suffixes of length at least L are the ranks L..n.</summary>
    internal override RankRange RangeOfRun(long length)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(length);
        if (length == 0) return Length == 0 ? RankRange.Empty : new RankRange(1, Length);
        return length > Length ? RankRange.Empty : new RankRange(length, Length);
    }

    internal override RankRange RangeOfFactor(long start, long length)
    {
        if (start < 0 || length < 0 || start > Length - length) throw new ArgumentOutOfRangeException(nameof(start), start, $"[{start}, {start} + {length}) is not a factor of a text of {Length} symbols");
        return RangeOfRun(length);
    }

    internal override long Down(long i)
    {
        CheckRank(i, nameof(i));
        return i < Length ? i + 1 : 0;
    }

    internal override long Up(long j)
    {
        CheckRank(j, nameof(j));
        return 0;
    }

    internal override long Next(long k)
    {
        CheckRank(k, nameof(k));
        return 0;
    }
}
