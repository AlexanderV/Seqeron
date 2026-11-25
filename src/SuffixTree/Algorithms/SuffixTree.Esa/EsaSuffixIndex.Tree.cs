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
    /// The nearest explicit ancestor above the locus (above the node below it for a point inside an edge). The parent depth is <c>max(LCP[lb], LCP[rb + 1])</c>; its borders are the previous and the next LCP value below that depth, found on
    /// the LCP pyramid of the navigation pack (<see cref="EsaNavigator64.Parent"/>).
    /// </summary>
    public bool TryParent(in IndexLocus locus, out IndexLocus parent)
    {
        EsaNode64 node = Resolve(locus, nameof(locus)); // the locus is copied into the handle before parent is written: parent may be the same variable as locus
        parent = default;
        if (node.Depth == 0) return false;
        parent = ToLocus(_nav.Parent(node));
        return true;
    }

    /// <inheritdoc />
    public bool TryChild(in IndexLocus node, int firstSymbol, out IndexLocus child)
    {
        EsaNode64 resolved = ResolveNode(node, nameof(node)); // copied before child is written: child may be the same variable as node
        child = default;
        if (!_nav.TryChild(resolved, firstSymbol, out EsaNode64 found)) return false;
        child = ToLocus(found);
        return true;
    }

    /// <inheritdoc />
    public int GetChildren(in IndexLocus node, Span<IndexLocus> children)
    {
        EsaNode64 resolved = ResolveNode(node, nameof(node)); // copied before the first child is written: children may hold the variable that holds node
        if (resolved.IsLeaf) return 0;
        var counter = default(CountSink);
        int count = _nav.Children(resolved, ref counter);
        if (count > children.Length) return count;
        var writer = new LocusSink(children);
        _ = _nav.Children(resolved, ref writer);
        return count;
    }

    /// <summary>
    /// The suffix link of an internal node: for <c>aY</c> the node <c>Y</c> (kernel (c2) of the navigation pack, at 64 bits); false for the root and for leaves.
    /// </summary>
    public bool TrySuffixLink(in IndexLocus node, out IndexLocus link)
    {
        EsaNode64 resolved = ResolveNode(node, nameof(node)); // copied before link is written: link may be the same variable as node
        link = default;
        if (resolved.Depth == 0 || resolved.IsLeaf) return false;
        link = ToLocus(_nav.SuffixLink(resolved));
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

    /// <summary>Writes the children into the caller's span (the caller has checked that they fit).</summary>
    private ref struct LocusSink : IChildSink
    {
        private readonly Span<IndexLocus> _destination;
        private int _next;

        internal LocusSink(Span<IndexLocus> destination)
        {
            _destination = destination;
            _next = 0;
        }

        public void Add(in EsaNode64 child) => _destination[_next++] = ToLocus(child);
    }
}
