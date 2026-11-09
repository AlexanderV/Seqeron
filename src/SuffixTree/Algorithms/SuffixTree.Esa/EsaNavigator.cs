using System.Runtime.CompilerServices;

namespace SuffixTree.Esa;

/// <summary>
/// <see cref="ISuffixTreeNavigator{TNode}"/> over an in-memory <see cref="EsaIndex"/> (work package 1b.4): the generic algorithms of <c>SuffixTree.Core</c> (LCS family, anchors, MEM, MUM) run on it
/// unchanged. Nodes are lcp-intervals (<see cref="EsaNode"/>), node for node the terminated Ukkonen tree (ADR-001 section 4). A value type like the other navigators, so that the JIT specialises
/// the generic algorithms. What each member costs: <see cref="GetEdgeSymbol"/> and <see cref="LengthOf"/> read the text at <c>Pos</c> (one random read) or nothing; <see cref="TryGetChild"/> walks the
/// l-indices of the node (the children are in ascending key order, the terminator first) and compares the symbol at the node's depth in the first suffix of each candidate (two dependent reads per
/// candidate); <see cref="GetSuffixLink"/> is the kernel of <see cref="EsaNavigationPack"/>; <see cref="CollectLeaves"/> is a contiguous copy of a suffix-array slice. The node handle needs no
/// parent pointer: the algorithms never ask for one (ADR-001 section 3).
/// </summary>
internal readonly struct EsaNavigator : ISuffixTreeNavigator<EsaNode>
{
    private readonly EsaNavigationPack _pack;
    private readonly string _text;
    private readonly ITextSource _source;
    private readonly int _n;
    private readonly Buf<int> _sa;
    private readonly LcpPyramid _lcp;
    private readonly ChildLayoutKind _kind;
    private readonly TaggedChild _tagged;
    private readonly PackedChild _packed;
    private readonly Reference3Child _reference3;

    /// <summary>Creates the navigator of <paramref name="index"/> (builds the navigation pack on first use).</summary>
    internal EsaNavigator(EsaIndex index)
    {
        ArgumentNullException.ThrowIfNull(index);
        _pack = index.Navigation;
        _text = index.Text;
        _source = new StringTextSource(index.Text);
        _n = index.Length;
        _sa = _pack.Sa;
        _lcp = _pack.Pyramid;
        _kind = index.Child.Kind;
        _tagged = index.Child.Tagged;
        _packed = index.Child.Packed;
        _reference3 = index.Child.Reference3;
    }

    /// <inheritdoc />
    public ITextSource Text => _source;

    /// <inheritdoc />
    public EsaNode Root
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get => new(0, _n, 0, 0, _n, true);
    }

    /// <inheritdoc />
    public EsaNode NullNode
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get => EsaNode.Null;
    }

    /// <inheritdoc />
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public bool IsNull(EsaNode node) => node.Lb < 0;

    /// <inheritdoc />
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public bool IsRoot(EsaNode node) => node.Depth == 0;

    /// <inheritdoc />
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public int GetEdgeSymbol(EsaNode node, int offset)
    {
        // Pos + ParentDepth + offset <= 2n: unsigned, so that it cannot wrap for a text of more than 2^30 symbols
        uint index = (uint)node.Pos + (uint)node.ParentDepth + (uint)offset;
        return index < (uint)_n ? _text[(int)index] : -1;
    }

    /// <inheritdoc />
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public int LengthOf(EsaNode node) => node.Depth - node.ParentDepth;

    /// <inheritdoc />
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    public EsaNode GetSuffixLink(EsaNode node) => _pack.SuffixLink(node);

    /// <inheritdoc />
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    public bool TryGetChild(EsaNode node, int key, out EsaNode child)
    {
        if (node.Lb == node.Rb)
        {
            child = EsaNode.Null;
            return false;
        }
        return _kind switch
        {
            ChildLayoutKind.Tagged => TryGetChildCore(_tagged, node, key, out child),
            ChildLayoutKind.Packed => TryGetChildCore(_packed, node, key, out child),
            _ => TryGetChildCore(_reference3, node, key, out child),
        };
    }

    /// <inheritdoc />
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    public void CollectLeaves(EsaNode node, int depth, List<int> results)
    {
        // rank 0 is the empty suffix: it belongs to the root (and is the root's terminator leaf), and it is no position
        int from = Math.Max(node.Lb, 1);
        int count = node.Rb - from + 1;
        if (count > 0) results.AddRange(new ReadOnlySpan<int>(_pack.Index.SuffixArray, from, count));
    }

    /// <inheritdoc />
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    public int FindAnyLeafPosition(EsaNode node, int depthFromRoot)
    {
        if (node.Lb == 0) return node.Rb >= 1 ? _sa[1] : -1; // the root, or its terminator leaf
        return node.Pos;
    }

    /// <inheritdoc />
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    public void GetChildren(EsaNode node, List<EsaNode> children)
    {
        children.Clear();
        if (node.Lb == node.Rb) return;
        switch (_kind)
        {
            case ChildLayoutKind.Tagged:
                GetChildrenCore(_tagged, node, children);
                break;
            case ChildLayoutKind.Packed:
                GetChildrenCore(_packed, node, children);
                break;
            default:
                GetChildrenCore(_reference3, node, children);
                break;
        }
    }

    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    private bool TryGetChildCore<TChild>(TChild table, in EsaNode node, int key, out EsaNode child)
        where TChild : struct, IChildTable
    {
        int rb = node.Rb;
        int a = node.Lb;
        int boundary = table.FirstLIndex(node.Lb, rb, node.LastChild);
        if (boundary >= 0 && (boundary <= a || boundary > rb)) ThrowCorruptChildTable();
        while (true)
        {
            int saA = _sa[a];
            int p = saA + node.Depth;
            int symbol = p < _n ? _text[p] : -1;
            if (symbol == key)
            {
                int b = boundary >= 0 ? boundary - 1 : rb;
                child = MakeChild(table, a, b, saA, node.Depth, boundary < 0);
                return true;
            }
            if (symbol > key || boundary < 0)
            {
                child = EsaNode.Null;
                return false;
            }
            a = boundary;
            boundary = table.NextLIndex(boundary, rb);
            if (boundary >= 0 && (boundary <= a || boundary > rb)) ThrowCorruptChildTable();
        }
    }

    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    private void GetChildrenCore<TChild>(TChild table, in EsaNode node, List<EsaNode> children)
        where TChild : struct, IChildTable
    {
        int rb = node.Rb;
        int a = node.Lb;
        int boundary = table.FirstLIndex(node.Lb, rb, node.LastChild);
        if (boundary >= 0 && (boundary <= a || boundary > rb)) ThrowCorruptChildTable();
        while (true)
        {
            int b = boundary >= 0 ? boundary - 1 : rb;
            children.Add(MakeChild(table, a, b, _sa[a], node.Depth, boundary < 0));
            if (boundary < 0) return;
            a = boundary;
            boundary = table.NextLIndex(boundary, rb);
            if (boundary >= 0 && (boundary <= a || boundary > rb)) ThrowCorruptChildTable();
        }
    }

    /// <summary>The child <c>[a..b]</c> whose first suffix starts at <paramref name="saA"/>: a leaf (depth = rest of the suffix + the terminator) or an interval whose value is the LCP at its first l-index.</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private EsaNode MakeChild<TChild>(TChild table, int a, int b, int saA, int parentDepth, bool childIsLast)
        where TChild : struct, IChildTable
        => a == b
            ? new EsaNode(a, a, _n - saA + 1, parentDepth, saA, childIsLast)
            : new EsaNode(a, b, _lcp.ValueAt(table.FirstLIndex(a, b, childIsLast)), parentDepth, saA, childIsLast);

    private static void ThrowCorruptChildTable()
        => throw new InvalidOperationException("the child table is inconsistent: the l-indices of an interval must increase and stay inside it");
}
