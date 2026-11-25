namespace SuffixTree.Esa;

/// <summary>
/// <see cref="ISuffixTreeNavigator{TNode}"/> over an index that the 32-bit navigator of S1 does not serve (work package S2.5): a container with 5-byte entries or a 64-bit LCP array, a view of a test, a text of any size and width. The generic
/// algorithms of <c>SuffixTree.Core</c> (the LCS family, anchors, MEM, MUM) are 32-bit (positions, lengths and counts are <see cref="int"/> in their contracts; S3 adds their 64-bit copies, ADR-001 section 8), so this navigator hands them
/// <see cref="EsaNode"/> and answers from the long navigator (<see cref="EsaNavigator64"/>), at any stored width: every number is narrowed at this one place, with <see cref="IntBridge"/>, and a text whose numbers do not fit is
/// <see cref="NotSupportedException"/> that names the operations and the 64-bit navigator, at the creation of the navigator (the length of the text) or at the first node that does not fit. The navigator of S1 (<see cref="EsaNavigator"/>)
/// stays the code for the indexes that the 32-bit forms hold; the selection is made once per index (<see cref="EsaIndex.IsNarrow"/>).
/// </summary>
internal readonly struct EsaNavigatorWide : ISuffixTreeNavigator<EsaNode>
{
    private const string Operations = "MEM, MUM, anchors and the longest common substrings (the generic algorithms of SuffixTree.Core are 32-bit; the long copies are stage S3's)";

    private readonly EsaNavigator64 _nav;
    private readonly EsaText _text;
    private readonly SaView _sa;
    private readonly ITextSource _source;
    private readonly int _n;
    private readonly bool _mapped;

    /// <summary>Creates the navigator of <paramref name="index"/>.</summary>
    /// <exception cref="NotSupportedException">The text has more than <see cref="int.MaxValue"/> symbols (the generic algorithms are 32-bit), or the container cannot link.</exception>
    internal EsaNavigatorWide(EsaIndex index)
    {
        ArgumentNullException.ThrowIfNull(index);
        _n = IntBridge.Length(index.TextLength, "the text of " + Operations);
        _nav = new EsaNavigator64(index);
        _nav.Pack.EnsureLinkable();
        _text = index.TextView;
        _sa = index.Sa;
        _source = index.TextView.AsTextSource();
        _mapped = index.IsMapped;
    }

    /// <inheritdoc />
    public ITextSource Text => _source;

    /// <inheritdoc />
    public EsaNode Root => new(0, _n, 0, 0, _n, true);

    /// <inheritdoc />
    public EsaNode NullNode => EsaNode.Null;

    /// <inheritdoc />
    public bool IsNull(EsaNode node) => node.Lb < 0;

    /// <inheritdoc />
    public bool IsRoot(EsaNode node) => node.Depth == 0;

    /// <inheritdoc />
    public int GetEdgeSymbol(EsaNode node, int offset)
    {
        long index = (long)node.Pos + node.ParentDepth + offset;
        return (ulong)index < (ulong)_n ? _text[index] : -1;
    }

    /// <inheritdoc />
    public int LengthOf(EsaNode node)
    {
        if (_mapped && node.Lb < 0) throw new InvalidDataException("the tree of the container has no child where its text and tables say there must be one: the file is damaged");
        return node.Depth - node.ParentDepth;
    }

    /// <inheritdoc />
    public EsaNode GetSuffixLink(EsaNode node) => Narrow(_nav.SuffixLink(EsaNode64.Of(node)));

    /// <inheritdoc />
    public bool TryGetChild(EsaNode node, int key, out EsaNode child)
    {
        if (_nav.TryChild(EsaNode64.Of(node), key, out EsaNode64 found))
        {
            child = Narrow(found);
            return true;
        }
        child = EsaNode.Null;
        return false;
    }

    /// <inheritdoc />
    public void GetChildren(EsaNode node, List<EsaNode> children)
    {
        children.Clear();
        var sink = new NarrowSink(children);
        _ = _nav.Children(EsaNode64.Of(node), ref sink);
    }

    /// <inheritdoc />
    public void CollectLeaves(EsaNode node, int depth, List<int> results)
    {
        // rank 0 is the empty suffix: it belongs to the root, and it is no position
        for (long rank = Math.Max(node.Lb, 1); rank <= node.Rb; rank++) results.Add(Position(_sa[rank]));
    }

    /// <inheritdoc />
    public int FindAnyLeafPosition(EsaNode node, int depthFromRoot)
    {
        if (node.Lb == 0) return node.Rb >= 1 ? Position(_sa[1]) : -1; // the root, or its terminator leaf
        return node.Pos;
    }

    /// <summary>A suffix array entry as a position: a number that is no position of the text is damage (the entry of a container is untrusted).</summary>
    private int Position(long entry)
    {
        if ((ulong)entry >= (ulong)_n) throw new InvalidDataException($"the suffix array holds {entry}, which is no position of a text of {_n} symbols: the file is damaged");
        return checked((int)entry);
    }

    private static EsaNode Narrow(in EsaNode64 node)
    {
        if (node.IsNull) return EsaNode.Null;
        return new EsaNode(IntBridge.Length(node.Lb, Operations), IntBridge.Length(node.Rb, Operations), IntBridge.Length(node.Depth, Operations), IntBridge.Length(node.ParentDepth, Operations), IntBridge.Length(node.Pos, Operations), node.LastChild);
    }

    private readonly struct NarrowSink(List<EsaNode> list) : IChildSink
    {
        public void Add(in EsaNode64 child) => list.Add(Narrow(child));
    }
}
