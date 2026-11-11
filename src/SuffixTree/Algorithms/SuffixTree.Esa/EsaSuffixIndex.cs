using SuffixTree.Indexing;

namespace SuffixTree.Esa;

/// <summary>
/// STUB of work package 1b.5 (red-first state): every member throws <see cref="NotImplementedException"/>. Replaced by the implementation in the next commit.
/// </summary>
internal sealed class EsaSuffixIndex : ISuffixTreeNavigation, ISuffixArrayAccess
{
    internal EsaSuffixIndex(EsaIndex index) => throw new NotImplementedException();

    internal static EsaSuffixIndex Build(string text, ChildLayoutKind? layout = null) => throw new NotImplementedException();

    public string Engine => throw new NotImplementedException();

    public long TextLength => throw new NotImplementedException();

    public IIndexText Text => throw new NotImplementedException();

    public IndexLocus Root => throw new NotImplementedException();

    public IndexSupport Support(IndexOperation operation) => throw new NotImplementedException();

    public string? SupportNote(IndexOperation operation) => throw new NotImplementedException();

    public long Count(in IndexLocus locus) => throw new NotImplementedException();

    public bool TryExtendRight(in IndexLocus locus, int symbol, out IndexLocus result) => throw new NotImplementedException();

    public bool TryExtendRight(in IndexLocus locus, ReadOnlySpan<char> pattern, out IndexLocus result) => throw new NotImplementedException();

    public bool TryExtendLeft(in IndexLocus locus, int symbol, out IndexLocus result) => throw new NotImplementedException();

    public int Locate(in IndexLocus locus, long skip, Span<long> positions) => throw new NotImplementedException();

    public bool IsExplicitNode(in IndexLocus locus) => throw new NotImplementedException();

    public bool TryParent(in IndexLocus locus, out IndexLocus parent) => throw new NotImplementedException();

    public bool TryChild(in IndexLocus node, int firstSymbol, out IndexLocus child) => throw new NotImplementedException();

    public int GetChildren(in IndexLocus node, Span<IndexLocus> children) => throw new NotImplementedException();

    public bool TrySuffixLink(in IndexLocus node, out IndexLocus link) => throw new NotImplementedException();

    public int EdgeSymbol(in IndexLocus node, long offset) => throw new NotImplementedException();

    public long SuffixAt(long rank) => throw new NotImplementedException();

    public long RankOf(long position) => throw new NotImplementedException();

    public long LcpAt(long rank) => throw new NotImplementedException();

    public long Lcp(long rankA, long rankB) => throw new NotImplementedException();

    public IndexLocus Lca(long rankA, long rankB) => throw new NotImplementedException();
}
