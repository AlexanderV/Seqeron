using SuffixTree.Indexing;

namespace SuffixTree.Esa;

public sealed partial class EsaSuffixTree
{
    /// <inheritdoc />
    public string Engine
    {
        get
        {
            using ReadScope scope = Enter();
            return _contract.Engine;
        }
    }

    /// <inheritdoc />
    public long TextLength
    {
        get
        {
            using ReadScope scope = Enter();
            return _contract.TextLength;
        }
    }

    /// <remarks>The 64-bit view of the text; <see cref="Text"/> is the <see cref="ITextSource"/> view of the same string.</remarks>
    IIndexText ISuffixIndex.Text
    {
        get
        {
            using ReadScope scope = Enter();
            return _guardedIndexText ?? _contract.Text;
        }
    }

    /// <inheritdoc />
    public IndexLocus Root
    {
        get
        {
            using ReadScope scope = Enter();
            return _contract.Root;
        }
    }

    /// <inheritdoc />
    public IndexSupport Support(IndexOperation operation)
    {
        using ReadScope scope = Enter();
        return _contract.Support(operation);
    }

    /// <inheritdoc />
    public string? SupportNote(IndexOperation operation)
    {
        using ReadScope scope = Enter();
        return _contract.SupportNote(operation);
    }

    /// <inheritdoc />
    public long Count(in IndexLocus locus)
    {
        using ReadScope scope = Enter();
        return _contract.Count(in locus);
    }

    /// <inheritdoc />
    public bool TryExtendRight(in IndexLocus locus, int symbol, out IndexLocus result)
    {
        using ReadScope scope = Enter();
        return _contract.TryExtendRight(in locus, symbol, out result);
    }

    /// <inheritdoc />
    public bool TryExtendRight(in IndexLocus locus, ReadOnlySpan<char> pattern, out IndexLocus result)
    {
        using ReadScope scope = Enter();
        return _contract.TryExtendRight(in locus, pattern, out result);
    }

    /// <inheritdoc />
    public bool TryExtendLeft(in IndexLocus locus, int symbol, out IndexLocus result)
    {
        using ReadScope scope = Enter();
        return _contract.TryExtendLeft(in locus, symbol, out result);
    }

    /// <inheritdoc />
    public int Locate(in IndexLocus locus, long skip, Span<long> positions)
    {
        using ReadScope scope = Enter();
        return _contract.Locate(in locus, skip, positions);
    }

    /// <inheritdoc />
    public bool IsExplicitNode(in IndexLocus locus)
    {
        using ReadScope scope = Enter();
        return _contract.IsExplicitNode(in locus);
    }

    /// <inheritdoc />
    public bool TryParent(in IndexLocus locus, out IndexLocus parent)
    {
        using ReadScope scope = Enter();
        return _contract.TryParent(in locus, out parent);
    }

    /// <inheritdoc />
    public bool TryChild(in IndexLocus node, int firstSymbol, out IndexLocus child)
    {
        using ReadScope scope = Enter();
        return _contract.TryChild(in node, firstSymbol, out child);
    }

    /// <inheritdoc />
    public int GetChildren(in IndexLocus node, Span<IndexLocus> children)
    {
        using ReadScope scope = Enter();
        return _contract.GetChildren(in node, children);
    }

    /// <inheritdoc />
    public bool TrySuffixLink(in IndexLocus node, out IndexLocus link)
    {
        using ReadScope scope = Enter();
        return _contract.TrySuffixLink(in node, out link);
    }

    /// <inheritdoc />
    public int EdgeSymbol(in IndexLocus node, long offset)
    {
        using ReadScope scope = Enter();
        return _contract.EdgeSymbol(in node, offset);
    }

    /// <inheritdoc />
    public long SuffixAt(long rank)
    {
        using ReadScope scope = Enter();
        return _contract.SuffixAt(rank);
    }

    /// <inheritdoc />
    public long RankOf(long position)
    {
        using ReadScope scope = Enter();
        return _contract.RankOf(position);
    }

    /// <inheritdoc />
    public long LcpAt(long rank)
    {
        using ReadScope scope = Enter();
        return _contract.LcpAt(rank);
    }

    /// <inheritdoc />
    public long Lcp(long rankA, long rankB)
    {
        using ReadScope scope = Enter();
        return _contract.Lcp(rankA, rankB);
    }

    /// <inheritdoc />
    public IndexLocus Lca(long rankA, long rankB)
    {
        using ReadScope scope = Enter();
        return _contract.Lca(rankA, rankB);
    }
}
