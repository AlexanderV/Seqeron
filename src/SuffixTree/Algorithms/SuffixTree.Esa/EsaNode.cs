namespace SuffixTree.Esa;

/// <summary>
/// A node of the terminated suffix tree as the ESA navigator hands it to the generic algorithms (work package 1b.4, ADR-001 section 4): the lcp-interval
/// <c>[Lb..Rb]</c> of its path label (a singleton for a leaf), the string depth <see cref="Depth"/> of the label (a leaf counts the terminator, so its depth is
/// <c>n - SA[Lb] + 1</c>), the string depth of the parent, <c>Pos = SA[Lb]</c> (an occurrence of the label, cached so that reading an edge symbol needs no suffix-array
/// access) and whether the node is the last child of its parent (the child-table field of an interval sits at <c>Rb</c> or at <c>Lb</c> by this flag). The root is
/// <c>(0, n, 0, 0, n, true)</c>; the null node has a negative <see cref="Lb"/>.
/// </summary>
internal readonly record struct EsaNode(int Lb, int Rb, int Depth, int ParentDepth, int Pos, bool LastChild)
{
    /// <summary>The null node (no node).</summary>
    internal static EsaNode Null => new(-1, -1, -1, -1, -1, false);

    /// <summary>True for a leaf (a singleton interval), and for the root of the empty text.</summary>
    internal bool IsLeaf => Lb == Rb;

    /// <summary>True for the null node.</summary>
    internal bool IsNull => Lb < 0;

    /// <summary>Length of the incoming edge, <c>Depth - ParentDepth</c> (0 for the root).</summary>
    internal int EdgeLength => Depth - ParentDepth;

    /// <summary>Number of ranks in the interval (the leaf count; the root counts the empty suffix too).</summary>
    internal int RankCount => Rb - Lb + 1;
}
