using SuffixTree.Esa.External;

namespace SuffixTree.Esa;

/// <summary>
/// The in-memory enhanced suffix array of a text (work package 1b.1): the suffix array and the compact LCP array of stage S1a, the child table, the header
/// statistics, and the lcp-interval tree on top of them. The index is node-for-node the terminated suffix tree of the text (ADR-001 section 4). It shares
/// the caller's string (no copy); the dense symbols of the preparation step are kept for the pattern search. Immutable after construction and safe for
/// concurrent readers. Since work package 1b.7f an index can also be the view of a mapped SESA container (<see cref="FromFile"/>): the suffix array, the LCP array, the child table and the packed text are then read in place
/// through <see cref="Sa"/>, <see cref="Lcp"/>, <see cref="Child"/> and <see cref="TextView"/>, which serve both kinds of index; the members that hand out the managed arrays (<see cref="SuffixArray"/>, <see cref="Text"/>,
/// <see cref="Prepared"/>) exist for the in-memory index only.
/// </summary>
internal sealed partial class EsaIndex
{
    private readonly string? _text;
    private readonly PreparedText? _prepared;
    private readonly SaView _sa;
    private readonly SesaFile? _file;
    private readonly MappedSearch? _mappedSearch;
    private readonly int _length;

    /// <summary>
    /// The child-table layout that <see cref="Build"/> uses by default for a text of <paramref name="n"/> symbols (decided by prototype P3, Log of PROGRESS-S1.md): the tagged
    /// single array while its 30 value bits suffice (it is the fastest: no LCP comparison on a read), else the LCP-disambiguated single array (the same 4 bytes per rank, text up to
    /// <see cref="Array.MaxLength"/>).
    /// </summary>
    internal static ChildLayoutKind DefaultLayoutFor(int n) => n <= TaggedChild.MaxLength ? ChildLayoutKind.Tagged : ChildLayoutKind.Packed;

    private EsaIndex(string text, PreparedText prepared, int[] sa, CompactLcp lcp, ChildTable child, EsaStatistics statistics)
    {
        _text = text;
        _prepared = prepared;
        _sa = new SaView(sa);
        _length = text.Length;
        TextView = EsaText.FromString(text);
        Lcp = lcp;
        Child = child;
        Statistics = statistics;
    }

    private EsaIndex(SesaFile file, EsaText text, SaView sa, CompactLcp lcp, ChildTable child, EsaStatistics statistics)
    {
        _file = file;
        _sa = sa;
        _length = text.Length;
        TextView = text;
        Lcp = lcp;
        Child = child;
        Statistics = statistics;
        // a container of 2^31 - 1 symbols or more can be opened (its header read) but has no tree here: no child table, nothing to search
        _mappedSearch = file.Header.TextLength < int.MaxValue ? new MappedSearch(file) : null;
    }

    /// <summary>
    /// The index that reads the sections of an opened container in place (work package 1b.7f): nothing is copied, the arrays stay in the mapping, which <paramref name="file"/> owns (disposing the index does not dispose it; the
    /// facade does).
    /// </summary>
    internal static EsaIndex FromFile(SesaFile file)
    {
        ArgumentNullException.ThrowIfNull(file);
        return new EsaIndex(file, file.TextView, file.Sa, file.Lcp, file.Child, file.Statistics);
    }

    /// <summary>The container this index reads in place, or null for an index that was built in memory.</summary>
    internal SesaFile? File => _file;

    /// <summary>True when the index reads a mapped container (the managed arrays do not exist then).</summary>
    internal bool IsMapped => _file is not null;

    /// <summary>The indexed text (the string given to <see cref="Build"/>); an index that reads a container has no string.</summary>
    /// <exception cref="NotSupportedException">The index reads a container: use <see cref="TextView"/>.</exception>
    internal string Text => _text ?? throw new NotSupportedException("an index that reads a container has no managed text: use TextView");

    /// <summary>The text as code units, for both kinds of index (a string or the packed section).</summary>
    internal EsaText TextView { get; }

    /// <summary>Number of symbols n.</summary>
    internal int Length => _length;

    /// <summary>The text after the order-preserving dense remap.</summary>
    /// <exception cref="NotSupportedException">The index reads a container: the packed section is the dense text.</exception>
    internal PreparedText Prepared => _prepared ?? throw new NotSupportedException("an index that reads a container has no prepared text: the packed section is the dense text");

    /// <summary><c>n + 1</c> entries, <c>SuffixArray[0] = n</c> (the empty suffix), then the suffixes in ascending order.</summary>
    /// <exception cref="NotSupportedException">The index reads a container: use <see cref="Sa"/>.</exception>
    internal int[] SuffixArray => _sa.Array ?? throw new NotSupportedException("an index that reads a container has no managed suffix array: use Sa");

    /// <summary>The suffix array as a view, for both kinds of index.</summary>
    internal SaView Sa => _sa;

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
