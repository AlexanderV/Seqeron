namespace SuffixTree.Esa;

/// <summary>
/// Header statistics of an ESA, produced by the child-table scan (formulas of ADR-001 section 4): <c>NodeCount = 1 + InternalCount + (Length + 1)</c> for
/// <c>Length &gt; 0</c> (root, internal nodes with value &gt; 0, one leaf per suffix including the empty one) and 1 for the empty text, exactly as the Ukkonen tree
/// reports; <c>LeafCount = MaxDepth = Length</c>. <c>LrsDepth</c> and <c>LrsRank</c> describe the deepest internal node (the longest repeated substring is
/// <c>text[SA[LrsRank] .. +LrsDepth)</c>; on ties the interval that closes first). <c>StackPeak</c> is the peak height of the scan stack.
/// </summary>
internal readonly record struct EsaStatistics(int Length, long NodeCount, int InternalCount, int LeafCount, int MaxDepth, int LrsDepth, int LrsRank, int StackPeak);

/// <summary>Counts the intervals and finds the deepest one while the scan runs.</summary>
internal struct StatisticsSink : IIntervalSink
{
    internal int Internal;
    internal int LrsDepth;
    internal int LrsRank;

    public readonly void LIndex(int previous, int k)
    {
    }

    public void Close(int ell, int lb, int rb, int first, bool lastChild)
    {
        if (ell <= 0) return;
        Internal++;
        if (ell > LrsDepth)
        {
            LrsDepth = ell;
            LrsRank = lb;
        }
    }
}

/// <summary>Forwards every scan event to two sinks.</summary>
internal struct TeeSink<TA, TB>(TA a, TB b) : IIntervalSink
    where TA : struct, IIntervalSink
    where TB : struct, IIntervalSink
{
    private TA _a = a;
    private TB _b = b;

    internal readonly TA First => _a;

    internal readonly TB Second => _b;

    public void LIndex(int previous, int k)
    {
        _a.LIndex(previous, k);
        _b.LIndex(previous, k);
    }

    public void Close(int ell, int lb, int rb, int first, bool lastChild)
    {
        _a.Close(ell, lb, rb, first, lastChild);
        _b.Close(ell, lb, rb, first, lastChild);
    }
}

/// <summary>
/// A built child table of one of the three layouts, with the slow, layout-switching accessors used by tests and the diagnostic paths
/// (<c>Traverse</c>). The hot paths use the layout structs directly through generics so that the JIT specialises them.
/// </summary>
internal sealed class ChildTable
{
    private ChildTable(ChildLayoutKind kind, Reference3Child reference, PackedChild packed, TaggedChild tagged)
    {
        Kind = kind;
        Reference3 = reference;
        Packed = packed;
        Tagged = tagged;
    }

    /// <summary>The layout.</summary>
    internal ChildLayoutKind Kind { get; }

    /// <summary>The three-array layout (default value unless <see cref="Kind"/> is <see cref="ChildLayoutKind.Reference3"/>).</summary>
    internal Reference3Child Reference3 { get; }

    /// <summary>The LCP-disambiguated single-array layout (default value unless <see cref="Kind"/> is <see cref="ChildLayoutKind.Packed"/>).</summary>
    internal PackedChild Packed { get; }

    /// <summary>The tagged single-array layout (default value unless <see cref="Kind"/> is <see cref="ChildLayoutKind.Tagged"/>).</summary>
    internal TaggedChild Tagged { get; }

    /// <summary>Bytes of the table arrays.</summary>
    internal long PayloadBytes => Kind switch
    {
        ChildLayoutKind.Reference3 => Reference3.PayloadBytes,
        ChildLayoutKind.Packed => Packed.PayloadBytes,
        _ => Tagged.PayloadBytes,
    };

    /// <summary>See <see cref="IChildTable.FirstLIndex"/>.</summary>
    internal int FirstLIndex(int lb, int rb, bool lastChild) => Kind switch
    {
        ChildLayoutKind.Reference3 => Reference3.FirstLIndex(lb, rb, lastChild),
        ChildLayoutKind.Packed => Packed.FirstLIndex(lb, rb, lastChild),
        _ => Tagged.FirstLIndex(lb, rb, lastChild),
    };

    /// <summary>See <see cref="IChildTable.NextLIndex"/>.</summary>
    internal int NextLIndex(int k, int rb) => Kind switch
    {
        ChildLayoutKind.Reference3 => Reference3.NextLIndex(k, rb),
        ChildLayoutKind.Packed => Packed.NextLIndex(k, rb),
        _ => Tagged.NextLIndex(k, rb),
    };

    /// <summary>
    /// Builds the child table of the given layout by one scan over the LCP array and returns the header statistics of the same scan.
    /// </summary>
    /// <exception cref="NotSupportedException">The layout cannot hold a text of this length (<see cref="TaggedChild.MaxLength"/>).</exception>
    internal static ChildTable Build(ChildLayoutKind kind, int n, CompactLcp lcp, out EsaStatistics statistics)
    {
        ArgumentNullException.ThrowIfNull(lcp);
        ArgumentOutOfRangeException.ThrowIfNegative(n);
        if (lcp.Length != n + 1) throw new ArgumentException($"the LCP array needs n + 1 = {n + 1} entries, got {lcp.Length}", nameof(lcp));
        switch (kind)
        {
            case ChildLayoutKind.Reference3:
                {
                    var layout = Reference3Child.Create(n);
                    statistics = Scan(n, lcp, layout);
                    return new ChildTable(kind, layout, default, default);
                }
            case ChildLayoutKind.Packed:
                {
                    var layout = PackedChild.Create(n, lcp);
                    statistics = Scan(n, lcp, layout);
                    return new ChildTable(kind, default, layout, default);
                }
            case ChildLayoutKind.Tagged:
                {
                    var layout = TaggedChild.Create(n);
                    statistics = Scan(n, lcp, layout);
                    return new ChildTable(kind, default, default, layout);
                }
            default:
                throw new ArgumentOutOfRangeException(nameof(kind), kind, "not a child-table layout");
        }
    }

    private static EsaStatistics Scan<TLayout>(int n, CompactLcp lcp, TLayout layout)
        where TLayout : struct, IIntervalSink
    {
        var stream = new CompactLcpStream(lcp);
        var sink = new TeeSink<TLayout, StatisticsSink>(layout, default);
        int peak = LcpIntervalScan.Run(n, ref stream, ref sink);
        StatisticsSink stats = sink.Second;
        return n == 0
            ? new EsaStatistics(0, 1, 0, 0, 0, 0, 0, 0)
            : new EsaStatistics(n, 1L + stats.Internal + (n + 1L), stats.Internal, n, n, stats.LrsDepth, stats.LrsRank, peak);
    }
}
