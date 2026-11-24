using SuffixTree.Indexing;

namespace SuffixTree.Esa;

internal sealed partial class EsaSuffixIndex
{
    private readonly bool _isMapped;

    /// <inheritdoc />
    public bool IsExplicitNode(in IndexLocus locus)
    {
        Resolve(locus, nameof(locus));
        return locus.Length == locus.Aux;
    }

    /// <summary>
    /// The nearest explicit ancestor above the locus (above the node below it for a point inside an edge). The parent depth is <c>max(LCP[lb], LCP[rb + 1])</c>; its borders are the previous and the
    /// next LCP value below that depth, found on the LCP pyramid of the navigation pack.
    /// </summary>
    public bool TryParent(in IndexLocus locus, out IndexLocus parent)
    {
        EsaNode64 node = Resolve(locus, nameof(locus)); // the locus is copied into the handle before parent is written: parent may be the same variable as locus
        parent = default;
        if (node.Depth == 0) return false;
        long parentDepth = node.ParentDepth;
        if (parentDepth == 0)
        {
            parent = Root;
            return true;
        }
        LcpPyramid pyramid = _index.Navigation.Pyramid;
        long left = pyramid.PrevLess(node.Lb, parentDepth);
        long right = pyramid.NextLess(node.Rb + 1, parentDepth);
        parent = new IndexLocus(left, right - 1, parentDepth, parentDepth);
        return true;
    }

    /// <inheritdoc />
    public bool TryChild(in IndexLocus node, int firstSymbol, out IndexLocus child)
    {
        EsaNode64 resolved = ResolveNode(node, nameof(node)); // copied before child is written: child may be the same variable as node
        child = default;
        if (!FindChild(resolved, firstSymbol, out EsaNode64 found)) return false;
        child = ToLocus(found);
        return true;
    }

    /// <inheritdoc />
    public int GetChildren(in IndexLocus node, Span<IndexLocus> children)
    {
        EsaNode64 resolved = ResolveNode(node, nameof(node)); // copied before the first child is written: children may hold the variable that holds node
        if (resolved.IsLeaf) return 0;
        int count = WalkChildren(resolved, default, false);
        if (count > children.Length) return count;
        _ = WalkChildren(resolved, children, true);
        return count;
    }

    /// <summary>
    /// The suffix link of an internal node: for <c>aY</c> the node <c>Y</c> (kernel (c2) of the navigation pack); false for the root and for leaves.
    /// </summary>
    public bool TrySuffixLink(in IndexLocus node, out IndexLocus link)
    {
        EsaNode64 resolved = ResolveNode(node, nameof(node)); // copied before link is written: link may be the same variable as node
        link = default;
        if (resolved.Depth == 0 || resolved.IsLeaf) return false;
        EsaNavigationPack pack = _index.Navigation; // the pack is 32-bit until S2.5: it refuses a text above 2^31 - 2 symbols by name
        link = ToLocus(EsaNode64.Of(pack.SuffixLink(resolved.ToNode("EsaSuffixIndex.TrySuffixLink (the suffix-link kernel of the navigation pack is 32-bit until S2.5)"))));
        return true;
    }

    /// <summary>
    /// The text symbol at <c>SA[lb] + parentDepth + offset</c>, the terminator -1 from the end of the text on. Inside the incoming edge every occurrence agrees; behind it (offsets at or beyond the edge
    /// length) the first suffix of the interval is read, as the navigator of the generic algorithms does.
    /// </summary>
    public int EdgeSymbol(in IndexLocus node, long offset)
    {
        EsaNode64 resolved = ResolveNode(node, nameof(node));
        ArgumentOutOfRangeException.ThrowIfNegative(offset);
        if (offset >= _n) return -1;
        long index = resolved.Pos + resolved.ParentDepth + offset;
        return index < _n ? _index.TextView[index] : -1;
    }

    /// <summary>The child of <paramref name="node"/> whose edge starts with <paramref name="key"/> (-1: the terminator).</summary>
    private bool FindChild(in EsaNode64 node, int key, out EsaNode64 child)
    {
        if (node.IsLeaf)
        {
            child = EsaNode64.Null;
            return false;
        }
        SaView sa = _index.Sa;
        EsaText text = _index.TextView;
        ChildTable table = _index.Child;
        long rb = node.Rb;
        long a = node.Lb;
        long boundary = table.FirstLIndex(node.Lb, rb, node.LastChild);
        CheckBoundary(boundary, a, rb);
        while (true)
        {
            long saA = sa[a];
            long at = saA + node.Depth;
            int symbol = (ulong)at < (ulong)_n ? text[at] : -1;
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
            boundary = table.NextLIndex(boundary, rb);
            CheckBoundary(boundary, a, rb);
        }
    }

    /// <summary>Walks the children of <paramref name="node"/> in key order; counts them, and writes them when <paramref name="write"/> is set (the caller has checked that they fit).</summary>
    private int WalkChildren(in EsaNode64 node, Span<IndexLocus> destination, bool write)
    {
        SaView sa = _index.Sa;
        ChildTable table = _index.Child;
        long rb = node.Rb;
        long a = node.Lb;
        int count = 0;
        long boundary = table.FirstLIndex(node.Lb, rb, node.LastChild);
        CheckBoundary(boundary, a, rb);
        while (true)
        {
            if (write)
            {
                long b = boundary >= 0 ? boundary - 1 : rb;
                destination[count] = ToLocus(MakeChild(a, b, sa[a], node.Depth, boundary < 0));
            }
            count++;
            if (boundary < 0) return count;
            a = boundary;
            boundary = table.NextLIndex(boundary, rb);
            CheckBoundary(boundary, a, rb);
        }
    }

    /// <summary>The child <c>[a..b]</c> whose first suffix starts at <paramref name="saA"/>: a leaf (depth = rest of the suffix + the terminator) or an interval whose value is the LCP at its first l-index.</summary>
    private EsaNode64 MakeChild(long a, long b, long saA, long parentDepth, bool childIsLast)
    {
        if (a == b)
        {
            if (_isMapped && _n - saA + 1 <= parentDepth) ThrowCorruptChildTable();
            return new EsaNode64(a, a, _n - saA + 1, parentDepth, saA, childIsLast);
        }
        long first = _index.Child.FirstLIndex(a, b, childIsLast);
        if (first <= a || first > b) ThrowCorruptChildTable();
        long depth = _index.Lcp[first];
        if (_isMapped && depth <= parentDepth) ThrowCorruptChildTable();
        return new EsaNode64(a, b, depth, parentDepth, saA, childIsLast);
    }

    private void CheckBoundary(long boundary, long previous, long rb)
    {
        if (boundary >= 0 && (boundary <= previous || boundary > rb)) ThrowCorruptChildTable();
    }

    private void ThrowCorruptChildTable()
        => throw (_isMapped
            ? new InvalidDataException("the child table is inconsistent: the l-indices of an interval must increase and stay inside it: the file is damaged")
            : new InvalidOperationException("the child table is inconsistent: the l-indices of an interval must increase and stay inside it"));
}
