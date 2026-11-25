namespace SuffixTree.Esa;

/// <summary>The text of the long navigator: the symbol at a position in the alphabet of the keys that <see cref="EsaNavigator64.TryChild"/> is asked for (code units for an index, whatever the test says for generated readers).</summary>
internal interface ILetters
{
    /// <summary>The symbol at <paramref name="position"/> (<c>0 &lt;= position &lt; n</c>).</summary>
    int At(long position);
}

/// <summary>The text of an index as <see cref="ILetters"/> (a string or the packed section through <see cref="EsaText"/>).</summary>
internal sealed class EsaLetters(EsaText text) : ILetters
{
    /// <inheritdoc />
    public int At(long position) => text[position];
}

/// <summary>The child table and the LCP array of an index as the cells that the long navigator reads (<see cref="IWideCells"/>).</summary>
internal sealed class ChildTableCells(ChildTable table, CompactLcp lcp) : IWideCells
{
    /// <inheritdoc />
    public long First(long lb, long rb, bool lastChild) => table.FirstLIndex(lb, rb, lastChild);

    /// <inheritdoc />
    public long Next(long k, long rb) => table.NextLIndex(k, rb);

    /// <inheritdoc />
    public long Depth(long first) => lcp[first];
}

/// <summary>Receives the children of a node in key order (<see cref="EsaNavigator64.Children{TSink}"/>); a struct, so that the walk is specialised and allocates nothing.</summary>
internal interface IChildSink
{
    /// <summary>One child.</summary>
    void Add(in EsaNode64 child);
}

/// <summary>Counts the children.</summary>
internal struct CountSink : IChildSink
{
    /// <summary>Number of children seen.</summary>
    internal int Count;

    /// <inheritdoc />
    public void Add(in EsaNode64 child) => Count++;
}

/// <summary>Collects the children into a list.</summary>
internal readonly struct ListSink(List<EsaNode64> list) : IChildSink
{
    /// <inheritdoc />
    public void Add(in EsaNode64 child) => list.Add(child);
}

/// <summary>
/// The navigator of the tree at 64 bits (work package S2.5): the nodes are <see cref="EsaNode64"/>, every rank, depth and position is <see cref="long"/>, and nothing is narrowed. It is what <see cref="EsaSuffixIndex"/> asks for the
/// children, the parent, the suffix link and the range queries of a node (the contract members), and what the navigator of the generic algorithms is built on for an index that the 32-bit navigator cannot serve
/// (<see cref="EsaNavigatorWide"/>). The 32-bit <see cref="EsaNode"/> stays the handle of the generic algorithms of <c>SuffixTree.Core</c> (they are 32-bit); this is the one navigator surface at 64 bits: <see cref="EsaNode64"/> converges
/// with it, the 32-bit navigator of S1 (<see cref="EsaNavigator"/>) is its specialisation for the indexes that the 32-bit forms hold. The sources are abstract (<see cref="PosSource"/>, <see cref="IWideCells"/>, <see cref="ILetters"/>), so the
/// same code runs over the views of an index at any stored width and over the generated readers of a test at 2^32 symbols and more. The pack (<see cref="EsaNavigationPack"/>) is asked for only by the operations that need it (a node
/// is never built from it): the children of a node cost no pack.
/// </summary>
internal sealed class EsaNavigator64
{
    private readonly long _n;
    private readonly PosSource _sa;
    private readonly IWideCells _cells;
    private readonly ILetters _text;
    private readonly Func<EsaNavigationPack> _pack;
    private readonly bool _mapped;

    /// <summary>The navigator of <paramref name="index"/> (the pack is built when an operation first needs it).</summary>
    internal EsaNavigator64(EsaIndex index)
        : this(index.TextLength, PosSource.Of(index.Sa), new ChildTableCells(index.Child, index.Lcp), new EsaLetters(index.TextView), () => index.Navigation, index.IsMapped)
    {
    }

    /// <summary>The navigator over sources.</summary>
    /// <param name="n">The length of the text.</param>
    /// <param name="sa">The suffix array, <c>n + 1</c> entries.</param>
    /// <param name="cells">The first and next l-indices of the intervals and the LCP value at an l-index.</param>
    /// <param name="text">The text.</param>
    /// <param name="pack">The pack, asked for by the suffix link, the parent and the range queries.</param>
    /// <param name="mapped">True when the sources are untrusted (a container, views of a test): an inconsistent table is <see cref="InvalidDataException"/>, not <see cref="InvalidOperationException"/>.</param>
    internal EsaNavigator64(long n, PosSource sa, IWideCells cells, ILetters text, Func<EsaNavigationPack> pack, bool mapped)
    {
        _n = n;
        _sa = sa;
        _cells = cells;
        _text = text;
        _pack = pack;
        _mapped = mapped;
    }

    /// <summary>The pack (built on first use).</summary>
    internal EsaNavigationPack Pack => _pack();

    /// <summary>The root node <c>(0, n, 0, 0, n, true)</c>.</summary>
    internal EsaNode64 Root => new(0, _n, 0, 0, _n, true);

    /// <summary>The child of <paramref name="node"/> whose edge starts with <paramref name="key"/> (-1: the terminator); false when there is none (a leaf has no child).</summary>
    internal bool TryChild(in EsaNode64 node, int key, out EsaNode64 child)
    {
        if (node.IsLeaf)
        {
            child = EsaNode64.Null;
            return false;
        }
        long rb = node.Rb;
        long a = node.Lb;
        long boundary = _cells.First(node.Lb, rb, node.LastChild);
        CheckBoundary(boundary, a, rb);
        while (true)
        {
            long saA = _sa[a];
            long at = saA + node.Depth;
            int symbol = (ulong)at < (ulong)_n ? _text.At(at) : -1;
            if (symbol == key)
            {
                long b = boundary >= 0 ? boundary - 1 : rb;
                child = MakeChild(a, b, saA, node.Depth, boundary < 0);
                return true;
            }
            if (symbol > key || boundary < 0)
            {
                child = EsaNode64.Null;
                return false;
            }
            a = boundary;
            boundary = _cells.Next(boundary, rb);
            CheckBoundary(boundary, a, rb);
        }
    }

    /// <summary>Walks the children of <paramref name="node"/> in key order, hands each to <paramref name="sink"/> and returns their number (a leaf has none).</summary>
    internal int Children<TSink>(in EsaNode64 node, ref TSink sink)
        where TSink : struct, IChildSink, allows ref struct
    {
        if (node.IsLeaf) return 0;
        long rb = node.Rb;
        long a = node.Lb;
        int count = 0;
        long boundary = _cells.First(node.Lb, rb, node.LastChild);
        CheckBoundary(boundary, a, rb);
        while (true)
        {
            long b = boundary >= 0 ? boundary - 1 : rb;
            EsaNode64 child = MakeChild(a, b, _sa[a], node.Depth, boundary < 0);
            sink.Add(in child);
            count++;
            if (boundary < 0) return count;
            a = boundary;
            boundary = _cells.Next(boundary, rb);
            CheckBoundary(boundary, a, rb);
        }
    }

    /// <summary>The suffix link (see <see cref="EsaNavigationPack.SuffixLink(in EsaNode64)"/>).</summary>
    internal EsaNode64 SuffixLink(in EsaNode64 node) => _pack().SuffixLink(node);

    /// <summary>The parent of an explicit node (see <see cref="EsaNavigationPack.Parent"/>).</summary>
    internal EsaNode64 Parent(in EsaNode64 node) => _pack().Parent(node);

    /// <summary>The Lcp of two ranks (see <see cref="EsaNavigationPack.Lcp"/>).</summary>
    internal long Lcp(long rankA, long rankB) => _pack().Lcp(rankA, rankB);

    /// <summary>The lowest common ancestor of two ranks (see <see cref="EsaNavigationPack.Lca"/>).</summary>
    internal (long Lb, long Rb, long Depth) Lca(long rankA, long rankB) => _pack().Lca(rankA, rankB);

    /// <summary>The child <c>[a..b]</c> whose first suffix starts at <paramref name="saA"/>: a leaf (depth = rest of the suffix + the terminator) or an interval whose value is the LCP at its first l-index.</summary>
    private EsaNode64 MakeChild(long a, long b, long saA, long parentDepth, bool childIsLast)
    {
        if (a == b)
        {
            if (_mapped && _n - saA + 1 <= parentDepth) ThrowCorruptChildTable();
            return new EsaNode64(a, a, _n - saA + 1, parentDepth, saA, childIsLast);
        }
        long first = _cells.First(a, b, childIsLast);
        if (first <= a || first > b) ThrowCorruptChildTable();
        long depth = _cells.Depth(first);
        if (_mapped && depth <= parentDepth) ThrowCorruptChildTable();
        return new EsaNode64(a, b, depth, parentDepth, saA, childIsLast);
    }

    private void CheckBoundary(long boundary, long previous, long rb)
    {
        if (boundary >= 0 && (boundary <= previous || boundary > rb)) ThrowCorruptChildTable();
    }

    private void ThrowCorruptChildTable()
        => throw (_mapped
            ? new InvalidDataException("the child table is inconsistent: the l-indices of an interval must increase and stay inside it: the file is damaged")
            : new InvalidOperationException("the child table is inconsistent: the l-indices of an interval must increase and stay inside it"));
}
