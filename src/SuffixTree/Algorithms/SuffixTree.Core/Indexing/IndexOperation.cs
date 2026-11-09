namespace SuffixTree.Indexing;

/// <summary>
/// The operations of the index contract whose availability differs between engines. <see cref="ISuffixIndex.Support"/> reports, per
/// operation, whether an engine answers it natively, answers it correctly but slowly (or after building an auxiliary structure), or not at all
/// (ADR-001 section 6).
/// </summary>
public enum IndexOperation
{
    /// <summary>Number of start positions of a locus string (<see cref="ISuffixIndex.Count"/>).</summary>
    Count,

    /// <summary>Appending a symbol to the locus string (<see cref="ISuffixIndex.TryExtendRight(in IndexLocus, int, out IndexLocus)"/>).</summary>
    ExtendRight,

    /// <summary>Prepending a symbol to the locus string (<see cref="ISuffixIndex.TryExtendLeft"/>).</summary>
    ExtendLeft,

    /// <summary>Paged enumeration of the start positions (<see cref="ISuffixIndex.Locate"/>).</summary>
    LazyLocate,

    /// <summary>Random access to the indexed text (<see cref="ISuffixIndex.Text"/>).</summary>
    TextAccess,

    /// <summary>Nearest explicit ancestor (<see cref="ISuffixTreeNavigation.TryParent"/>).</summary>
    Parent,

    /// <summary>Child by first edge symbol, enumeration of children (<see cref="ISuffixTreeNavigation.TryChild"/>).</summary>
    Child,

    /// <summary>Suffix link of an internal node (<see cref="ISuffixTreeNavigation.TrySuffixLink"/>).</summary>
    SuffixLink,

    /// <summary>String depth of a locus (<see cref="IndexLocus.Length"/>) and edge access (<see cref="ISuffixTreeNavigation.EdgeSymbol"/>).</summary>
    StringDepth,

    /// <summary><c>SA[rank]</c> (<see cref="ISuffixArrayAccess.SuffixAt"/>).</summary>
    SuffixArrayAccess,

    /// <summary><c>ISA[position]</c> (<see cref="ISuffixArrayAccess.RankOf"/>).</summary>
    InverseSuffixArrayAccess,

    /// <summary><c>LCP[rank]</c> (<see cref="ISuffixArrayAccess.LcpAt"/>).</summary>
    LcpAccess,

    /// <summary>Longest common prefix of two ranks (<see cref="ISuffixArrayAccess.Lcp"/>).</summary>
    Lcp,

    /// <summary>Lowest common ancestor of two leaves (<see cref="ISuffixArrayAccess.Lca"/>).</summary>
    Lca,
}
