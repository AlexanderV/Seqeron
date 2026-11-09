using System;

namespace SuffixTree.Indexing;

/// <summary>Navigation in the virtual suffix tree of an index (ADR-001 section 6). The tree is the terminated tree: key -1 is the terminator child.</summary>
public interface ISuffixTreeNavigation : ISuffixIndex
{
    /// <summary>True when <paramref name="locus"/> is an explicit node (the root, an internal node or a leaf) and not a point inside an edge.</summary>
    bool IsExplicitNode(in IndexLocus locus);

    /// <summary>The nearest explicit ancestor above <paramref name="locus"/>; false for the root.</summary>
    bool TryParent(in IndexLocus locus, out IndexLocus parent);

    /// <summary>The child of <paramref name="node"/> whose edge starts with <paramref name="firstSymbol"/> (-1 selects the terminator child); false if there is none.</summary>
    bool TryChild(in IndexLocus node, int firstSymbol, out IndexLocus child);

    /// <summary>
    /// Writes the children of <paramref name="node"/> in ascending key order (terminator first) when they fit and returns their number; when the
    /// number exceeds <paramref name="children"/>.Length nothing is written.
    /// </summary>
    int GetChildren(in IndexLocus node, Span<IndexLocus> children);

    /// <summary>The suffix link of an internal node; false for the root and for leaves.</summary>
    bool TrySuffixLink(in IndexLocus node, out IndexLocus link);

    /// <summary>The symbol at <paramref name="offset"/> on the incoming edge of <paramref name="node"/>; -1 for an offset at or beyond the text end.</summary>
    int EdgeSymbol(in IndexLocus node, long offset);
}
