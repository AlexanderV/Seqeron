using SuffixTree.Indexing;

namespace SuffixTree.Esa;

internal sealed partial class EsaSuffixIndex
{
    /// <inheritdoc />
    public long SuffixAt(long rank)
    {
        CheckRank(rank, nameof(rank));
        return _index.SuffixArray[(int)rank];
    }

    /// <summary><c>ISA[position]</c> from the navigation pack (built on first use).</summary>
    public long RankOf(long position)
    {
        CheckRank(position, nameof(position));
        return _index.Navigation.Isa((int)position);
    }

    /// <summary><c>LCP[rank]</c> from the compact LCP array (no pack needed).</summary>
    public long LcpAt(long rank)
    {
        CheckRank(rank, nameof(rank));
        return _index.Lcp[(int)rank];
    }

    /// <summary>The range minimum <c>min LCP[min + 1 .. max]</c> on the LCP pyramid of the navigation pack; equal ranks give the length of the suffix.</summary>
    public long Lcp(long rankA, long rankB)
    {
        CheckRank(rankA, nameof(rankA));
        CheckRank(rankB, nameof(rankB));
        if (rankA == rankB) return _n - _index.SuffixArray[(int)rankA];
        return _index.Navigation.Pyramid.Min((int)Math.Min(rankA, rankB) + 1, (int)Math.Max(rankA, rankB));
    }

    /// <summary>
    /// The lowest common ancestor of the leaves of two ranks: the lcp-interval of value <c>m = min LCP[lo + 1 .. hi]</c> around them, whose borders are the previous and the next LCP value below <c>m</c> (the root
    /// when <c>m = 0</c>); the leaf itself for equal ranks, the root for the empty text.
    /// </summary>
    public IndexLocus Lca(long rankA, long rankB)
    {
        CheckRank(rankA, nameof(rankA));
        CheckRank(rankB, nameof(rankB));
        if (_n == 0) return Root;
        if (rankA == rankB)
        {
            int rank = (int)rankA;
            long depth = _n - _index.SuffixArray[rank] + 1;
            return new IndexLocus(rank, rank, depth, depth);
        }
        int lo = (int)Math.Min(rankA, rankB);
        int hi = (int)Math.Max(rankA, rankB);
        LcpPyramid pyramid = _index.Navigation.Pyramid;
        int value = pyramid.Min(lo + 1, hi);
        if (value == 0) return Root;
        return new IndexLocus(pyramid.PrevLess(lo, value), pyramid.NextLess(hi + 1, value) - 1, value, value);
    }
}
