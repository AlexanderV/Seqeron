#pragma warning disable S2325, S4144, S1144, S3453, CA1822 // stub of the red run (work package 1b.4): every member throws NotImplementedException; the implementation replaces this file
namespace SuffixTree.Esa;

/// <summary>
/// <see cref="ISuffixTreeNavigator{TNode}"/> over an in-memory <see cref="EsaIndex"/> (work package 1b.4): the generic algorithms of <c>SuffixTree.Core</c> (LCS family, anchors, MEM, MUM) run on
/// it unchanged. Nodes are lcp-intervals (<see cref="EsaNode"/>), node for node the terminated Ukkonen tree. A value type like the other navigators, so that the JIT specialises the generic
/// algorithms. Stub of the red run.
/// </summary>
internal readonly struct EsaNavigator : ISuffixTreeNavigator<EsaNode>
{
    /// <summary>Creates the navigator of <paramref name="index"/> (builds the navigation pack on first use).</summary>
    internal EsaNavigator(EsaIndex index)
    {
        _ = index;
        throw new NotImplementedException();
    }

    public ITextSource Text => throw new NotImplementedException();

    public EsaNode Root => throw new NotImplementedException();

    public EsaNode NullNode => throw new NotImplementedException();

    public bool IsNull(EsaNode node) => throw new NotImplementedException();

    public bool IsRoot(EsaNode node) => throw new NotImplementedException();

    public int GetEdgeSymbol(EsaNode node, int offset) => throw new NotImplementedException();

    public int LengthOf(EsaNode node) => throw new NotImplementedException();

    public EsaNode GetSuffixLink(EsaNode node) => throw new NotImplementedException();

    public bool TryGetChild(EsaNode node, int key, out EsaNode child) => throw new NotImplementedException();

    public void CollectLeaves(EsaNode node, int depth, List<int> results) => throw new NotImplementedException();

    public int FindAnyLeafPosition(EsaNode node, int depthFromRoot) => throw new NotImplementedException();

    public void GetChildren(EsaNode node, List<EsaNode> children) => throw new NotImplementedException();
}
