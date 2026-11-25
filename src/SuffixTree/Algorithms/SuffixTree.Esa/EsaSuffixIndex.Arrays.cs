using SuffixTree.Indexing;

namespace SuffixTree.Esa;

internal sealed partial class EsaSuffixIndex
{
    /// <inheritdoc />
    public long SuffixAt(long rank)
    {
        CheckRank(rank, nameof(rank));
        return _index.Sa[rank];
    }

    /// <summary><c>ISA[position]</c> from the navigation pack (built on first use).</summary>
    /// <exception cref="NotSupportedException">A container of more than <see cref="EsaNavigationPack.PackHeapLimit"/> symbols: it has no ISA section and the heap does not hold ISA at that size (the message names the limit).</exception>
    public long RankOf(long position)
    {
        CheckRank(position, nameof(position));
        return _nav.Pack.Isa(position, "EsaSuffixIndex.RankOf");
    }

    /// <summary><c>LCP[rank]</c> from the compact LCP array (no pack needed).</summary>
    public long LcpAt(long rank)
    {
        CheckRank(rank, nameof(rank));
        return _index.Lcp[rank];
    }

    /// <summary>The range minimum <c>min LCP[min + 1 .. max]</c> on the LCP pyramid of the navigation pack; equal ranks give the length of the suffix.</summary>
    public long Lcp(long rankA, long rankB)
    {
        CheckRank(rankA, nameof(rankA));
        CheckRank(rankB, nameof(rankB));
        return _nav.Lcp(rankA, rankB);
    }

    /// <summary>
    /// The lowest common ancestor of the leaves of two ranks: the lcp-interval of value <c>m = min LCP[lo + 1 .. hi]</c> around them, whose borders are the previous and the next LCP value below <c>m</c> (the root
    /// when <c>m = 0</c>); the leaf itself for equal ranks, the root for the empty text.
    /// </summary>
    public IndexLocus Lca(long rankA, long rankB)
    {
        CheckRank(rankA, nameof(rankA));
        CheckRank(rankB, nameof(rankB));
        (long lb, long rb, long depth) = _nav.Lca(rankA, rankB);
        return new IndexLocus(lb, rb, depth, depth);
    }
}
