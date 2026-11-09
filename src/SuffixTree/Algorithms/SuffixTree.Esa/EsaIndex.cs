namespace SuffixTree.Esa;

/// <summary>
/// The in-memory enhanced suffix array of a text (work package 1b.1): the suffix array and the compact LCP array of stage S1a, the child table, the header
/// statistics, and the lcp-interval tree on top of them. The index is node-for-node the terminated suffix tree of the text (ADR-001 section 4). It shares
/// the caller's string (no copy); the dense symbols of the preparation step are kept for the pattern search. Immutable after construction and safe for
/// concurrent readers.
/// </summary>
internal sealed partial class EsaIndex
{
    /// <summary>
    /// The child-table layout that <see cref="Build"/> uses by default for a text of <paramref name="n"/> symbols (decided by prototype P3, Log of PROGRESS-S1.md): the tagged
    /// single array while its 30 value bits suffice (it is the fastest: no LCP comparison on a read), else the LCP-disambiguated single array (the same 4 bytes per rank, text up to
    /// <see cref="Array.MaxLength"/>).
    /// </summary>
    internal static ChildLayoutKind DefaultLayoutFor(int n) => n <= TaggedChild.MaxLength ? ChildLayoutKind.Tagged : ChildLayoutKind.Packed;

    private EsaIndex(string text, PreparedText prepared, int[] sa, CompactLcp lcp, ChildTable child, EsaStatistics statistics)
    {
        Text = text;
        Prepared = prepared;
        SuffixArray = sa;
        Lcp = lcp;
        Child = child;
        Statistics = statistics;
    }

    /// <summary>The indexed text (the string given to <see cref="Build"/>).</summary>
    internal string Text { get; }

    /// <summary>Number of symbols n.</summary>
    internal int Length => Text.Length;

    /// <summary>The text after the order-preserving dense remap.</summary>
    internal PreparedText Prepared { get; }

    /// <summary><c>n + 1</c> entries, <c>SuffixArray[0] = n</c> (the empty suffix), then the suffixes in ascending order.</summary>
    internal int[] SuffixArray { get; }

    /// <summary><c>n + 1</c> entries, <c>LCP[0] = LCP[1] = 0</c>, <c>LCP[r] = lcp(SA[r-1], SA[r])</c>.</summary>
    internal CompactLcp Lcp { get; }

    /// <summary>The child table.</summary>
    internal ChildTable Child { get; }

    /// <summary>NodeCount, LeafCount, MaxDepth, the deepest internal node.</summary>
    internal EsaStatistics Statistics { get; }

    /// <summary>
    /// Builds the index: text preparation, SA-IS (pointer kernel), the direct LCP kernel, then one scan for the child table and the statistics.
    /// </summary>
    /// <exception cref="ArgumentNullException"><paramref name="text"/> is null.</exception>
    /// <exception cref="NotSupportedException">The text does not fit the 32-bit path (or the layout).</exception>
    internal static EsaIndex Build(string text, ChildLayoutKind? layout = null)
    {
        ArgumentNullException.ThrowIfNull(text);
        PreparedText prepared = TextPreparation.Prepare(text);
        int[] sa = SaisKernel.Build(prepared);
        CompactLcp lcp = LcpKernel.Build(LcpKernelKind.Direct, prepared, sa);
        return FromArrays(text, prepared, sa, lcp, layout ?? DefaultLayoutFor(text.Length));
    }

    /// <summary>
    /// The index of already built arrays (<paramref name="sa"/> and <paramref name="lcp"/> of <paramref name="prepared"/>): runs the child-table scan only. The arrays are shared, not
    /// copied. The measurement harness uses it to time the scan on its own and to compare the layouts on the same SA and LCP.
    /// </summary>
    internal static EsaIndex FromArrays(string text, PreparedText prepared, int[] sa, CompactLcp lcp, ChildLayoutKind layout)
    {
        ArgumentNullException.ThrowIfNull(text);
        ArgumentNullException.ThrowIfNull(prepared);
        ArgumentNullException.ThrowIfNull(sa);
        ArgumentNullException.ThrowIfNull(lcp);
        ChildTable child = ChildTable.Build(layout, text.Length, lcp, out EsaStatistics statistics);
        return new EsaIndex(text, prepared, sa, lcp, child, statistics);
    }

    /// <summary>The same index with a child table of another layout (the suffix array and the LCP array are shared).</summary>
    internal EsaIndex WithLayout(ChildLayoutKind layout) => FromArrays(Text, Prepared, SuffixArray, Lcp, layout);
}
