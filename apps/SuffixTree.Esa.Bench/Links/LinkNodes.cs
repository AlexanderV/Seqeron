#pragma warning disable S2325 // stubs of the red run: the bodies are written after the tests
namespace SuffixTree.Esa.Bench.Links;

/// <summary>All internal nodes (lcp-intervals of value &gt;= 1; the root is not a link source) of an index, with string depth and node depth (number of edges from the root).</summary>
internal sealed class LinkNodeSet
{
    /// <summary>Number of internal nodes.</summary>
    internal int Count => throw new NotImplementedException();

    /// <summary>Left boundaries.</summary>
    internal int[] Lb => throw new NotImplementedException();

    /// <summary>Right boundaries.</summary>
    internal int[] Rb => throw new NotImplementedException();

    /// <summary>String depths (the lcp-interval values).</summary>
    internal int[] Depth => throw new NotImplementedException();

    /// <summary>Node depths: 1 for a child of the root.</summary>
    internal int[] NodeDepth => throw new NotImplementedException();

    /// <summary>Enumerates the internal nodes by one stack scan over the LCP array.</summary>
    internal static LinkNodeSet Enumerate(CompactLcp lcp) => throw new NotImplementedException();

    /// <summary>The handle of node <paramref name="i"/>.</summary>
    internal LinkNode Node(int i, int[] sa) => throw new NotImplementedException();

    /// <summary>A deterministic random sample of <paramref name="count"/> node indices without repetition (all of them, shuffled, when fewer exist).</summary>
    internal int[] Sample(int count, ulong seed) => throw new NotImplementedException();

    /// <summary>The indices of the <paramref name="count"/> largest keys, largest first (ties: smaller index first).</summary>
    internal static int[] Largest(ReadOnlySpan<int> keys, int count) => throw new NotImplementedException();
}
