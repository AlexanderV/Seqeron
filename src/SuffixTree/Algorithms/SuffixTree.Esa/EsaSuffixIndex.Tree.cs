using SuffixTree.Indexing;

namespace SuffixTree.Esa;

internal sealed partial class EsaSuffixIndex
{
    private readonly ChildLayoutKind _kind;
    private readonly TaggedChild _tagged;
    private readonly PackedChild _packed;
    private readonly Reference3Child _reference3;

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
        EsaNode node = Resolve(locus, nameof(locus));
        parent = default;
        if (node.Depth == 0) return false;
        int parentDepth = node.ParentDepth;
        if (parentDepth == 0)
        {
            parent = Root;
            return true;
        }
        LcpPyramid pyramid = _index.Navigation.Pyramid;
        int left = pyramid.PrevLess(node.Lb, parentDepth);
        int right = pyramid.NextLess(node.Rb + 1, parentDepth);
        parent = new IndexLocus(left, right - 1, parentDepth, parentDepth);
        return true;
    }

    /// <inheritdoc />
    public bool TryChild(in IndexLocus node, int firstSymbol, out IndexLocus child)
    {
        EsaNode resolved = ResolveNode(node, nameof(node));
        child = default;
        if (!FindChild(resolved, firstSymbol, out EsaNode found)) return false;
        child = ToLocus(found);
        return true;
    }

    /// <inheritdoc />
    public int GetChildren(in IndexLocus node, Span<IndexLocus> children)
    {
        EsaNode resolved = ResolveNode(node, nameof(node));
        if (resolved.IsLeaf) return 0;
        int count = _kind switch
        {
            ChildLayoutKind.Tagged => WalkChildren(_tagged, resolved, default, false),
            ChildLayoutKind.Packed => WalkChildren(_packed, resolved, default, false),
            _ => WalkChildren(_reference3, resolved, default, false),
        };
        if (count > children.Length) return count;
        _ = _kind switch
        {
            ChildLayoutKind.Tagged => WalkChildren(_tagged, resolved, children, true),
            ChildLayoutKind.Packed => WalkChildren(_packed, resolved, children, true),
            _ => WalkChildren(_reference3, resolved, children, true),
        };
        return count;
    }

    /// <summary>
    /// The suffix link of an internal node: for <c>aY</c> the node <c>Y</c> (kernel (c2) of the navigation pack); false for the root and for leaves.
    /// </summary>
    public bool TrySuffixLink(in IndexLocus node, out IndexLocus link)
    {
        EsaNode resolved = ResolveNode(node, nameof(node));
        link = default;
        if (resolved.Depth == 0 || resolved.IsLeaf) return false;
        link = ToLocus(_index.Navigation.SuffixLink(resolved));
        return true;
    }

    /// <summary>
    /// The text symbol at <c>SA[lb] + parentDepth + offset</c>, the terminator -1 from the end of the text on. Inside the incoming edge every occurrence agrees; behind it (offsets at or beyond the edge
    /// length) the first suffix of the interval is read, as the navigator of the generic algorithms does.
    /// </summary>
    public int EdgeSymbol(in IndexLocus node, long offset)
    {
        EsaNode resolved = ResolveNode(node, nameof(node));
        ArgumentOutOfRangeException.ThrowIfNegative(offset);
        if (offset >= _n) return -1;
        long index = (long)resolved.Pos + resolved.ParentDepth + offset;
        return index < _n ? _index.Text[(int)index] : -1;
    }

    /// <summary>The child of <paramref name="node"/> whose edge starts with <paramref name="key"/> (-1: the terminator).</summary>
    private bool FindChild(in EsaNode node, int key, out EsaNode child)
    {
        if (node.IsLeaf)
        {
            child = EsaNode.Null;
            return false;
        }
        return _kind switch
        {
            ChildLayoutKind.Tagged => FindChildCore(_tagged, node, key, out child),
            ChildLayoutKind.Packed => FindChildCore(_packed, node, key, out child),
            _ => FindChildCore(_reference3, node, key, out child),
        };
    }

    private bool FindChildCore<TChild>(TChild table, in EsaNode node, int key, out EsaNode child)
        where TChild : struct, IChildTable
    {
        int[] sa = _index.SuffixArray;
        string text = _index.Text;
        int rb = node.Rb;
        int a = node.Lb;
        int boundary = table.FirstLIndex(node.Lb, rb, node.LastChild);
        CheckBoundary(boundary, a, rb);
        while (true)
        {
            int saA = sa[a];
            uint at = (uint)saA + (uint)node.Depth;
            int symbol = at < (uint)_n ? text[(int)at] : -1;
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
            CheckBoundary(boundary, a, rb);
        }
    }

    /// <summary>Walks the children of <paramref name="node"/> in key order; counts them, and writes them when <paramref name="write"/> is set (the caller has checked that they fit).</summary>
    private int WalkChildren<TChild>(TChild table, in EsaNode node, Span<IndexLocus> destination, bool write)
        where TChild : struct, IChildTable
    {
        int[] sa = _index.SuffixArray;
        int rb = node.Rb;
        int a = node.Lb;
        int count = 0;
        int boundary = table.FirstLIndex(node.Lb, rb, node.LastChild);
        CheckBoundary(boundary, a, rb);
        while (true)
        {
            if (write)
            {
                int b = boundary >= 0 ? boundary - 1 : rb;
                destination[count] = ToLocus(MakeChild(table, a, b, sa[a], node.Depth, boundary < 0));
            }
            count++;
            if (boundary < 0) return count;
            a = boundary;
            boundary = table.NextLIndex(boundary, rb);
            CheckBoundary(boundary, a, rb);
        }
    }

    /// <summary>The child <c>[a..b]</c> whose first suffix starts at <paramref name="saA"/>: a leaf (depth = rest of the suffix + the terminator) or an interval whose value is the LCP at its first l-index.</summary>
    private EsaNode MakeChild<TChild>(TChild table, int a, int b, int saA, int parentDepth, bool childIsLast)
        where TChild : struct, IChildTable
    {
        if (a == b) return new EsaNode(a, a, _n - saA + 1, parentDepth, saA, childIsLast);
        int first = table.FirstLIndex(a, b, childIsLast);
        if (first <= a || first > b) ThrowCorruptChildTable();
        return new EsaNode(a, b, _index.Lcp[first], parentDepth, saA, childIsLast);
    }

    private static void CheckBoundary(int boundary, int previous, int rb)
    {
        if (boundary >= 0 && (boundary <= previous || boundary > rb)) ThrowCorruptChildTable();
    }

    private static void ThrowCorruptChildTable()
        => throw new InvalidOperationException("the child table is inconsistent: the l-indices of an interval must increase and stay inside it");
}
