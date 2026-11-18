#pragma warning disable S2325 // stub: the members become instance members of the implementation
using System.Runtime.InteropServices;

namespace SuffixTree.Esa.External;

/// <summary>How the single-array child table (layout A) encodes a field: <see cref="Tagged"/> = kind in bits 30-31 (the in-RAM <c>TaggedChild</c>), <see cref="Packed"/> = the raw value (the in-RAM <c>PackedChild</c>).</summary>
internal enum ChildCodec
{
    Tagged,
    Packed,
}

/// <summary>Parameters of PH5 (ADR-004 section 3.2). The output never depends on any of them.</summary>
internal sealed record ChildBuildOptions
{
    /// <summary>Directory of the spill file of the scan stack (must exist).</summary>
    internal required string ScratchDir { get; init; }

    /// <summary>Entries per stack block (ADR: 256 K); the two top blocks are in RAM.</summary>
    internal int BlockEntries { get; init; } = 256 * 1024;

    /// <summary>Size of one mapped window of the CHILD writer (a power of two).</summary>
    internal int WindowBytes { get; init; } = 1 << 20;

    /// <summary>Windows the CHILD writer keeps mapped at a time.</summary>
    internal int WindowSlots { get; init; } = 8;

    /// <summary>msync a window when it is released (bounds the dirty pages).</summary>
    internal bool FlushWindows { get; init; } = true;

    /// <summary>fsync the spill file after this many unflushed bytes (0 = never).</summary>
    internal long SpillFlushBytes { get; init; } = 16L << 20;

    internal ChildCodec Codec { get; init; } = ChildCodec.Packed;

    /// <summary>Shared cap of the native buffers of the scan (null = unlimited).</summary>
    internal ArenaBudget? Budget { get; init; }

    /// <summary>Test hook: options of the arena (guard bands, poison pattern); null = production arena.</summary>
    internal ArenaOptions? ArenaGuards { get; init; }

    /// <summary>Receives the spill's I/O and scratch and the CHILD windows (null = not counted).</summary>
    internal IoCounters? Counters { get; init; }
}

/// <summary>What the forward scan found and did.</summary>
internal sealed class ChildScanResult
{
    internal long N { get; set; }

    /// <summary>Internal nodes with a value above 0 (the root is not counted).</summary>
    internal long InternalCount { get; set; }

    /// <summary>NODE_COUNT = <c>1 + I + n + 1</c> (root, internal nodes, one leaf per suffix including the empty one); 1 for the empty text.</summary>
    internal long NodeCount { get; set; }

    internal uint LrsDepth { get; set; }

    /// <summary>Left boundary of the first-closing deepest interval: the longest repeated substring starts at <c>SA[LrsRank]</c>.</summary>
    internal uint LrsRank { get; set; }

    internal long MaxLcp { get; set; }

    /// <summary>Peak height of the scan stack (at most <c>1 + max LCP</c>).</summary>
    internal long StackPeak { get; set; }

    internal SpillStats Spill { get; set; }

    internal WindowStats ChildWindows { get; set; }

    /// <summary>Seconds spent inside <see cref="ChildScan.Consume"/> and <see cref="ChildScan.Complete"/>.</summary>
    internal double Seconds { get; set; }
}

/// <summary>An open interval of the scan: value, left boundary, first and last l-index (<see cref="ChildScan.None"/> = none yet). 16 bytes (ADR-004 section 3.2: "12-16 B").</summary>
[StructLayout(LayoutKind.Sequential)]
internal struct IntervalEntry
{
    internal uint Ell;
    internal uint Lb;
    internal uint First;
    internal uint Last;
}

/// <summary>
/// The forward scan of PH5 (ADR-004 sections 3.2 and 4.2, row PH5): one pass over the LCP in rank order with the external stack, writing layout A (the child table in the single-array form, writes behind the cursor through
/// bounded mapped windows) and counting the internal nodes and the deepest interval (NODE_COUNT, LRS_DEPTH, LRS_RANK). It is an <see cref="ILcpChunkConsumer"/>: the PH4 hook
/// (<c>Consume(firstRank, values)</c>, dense resolved values in ascending rank order) feeds it fused with PH4-B3.
/// </summary>
internal sealed class ChildScan : ILcpChunkConsumer, IDisposable
{
    /// <summary>"No index": the stored value of a missing l-index.</summary>
    internal const uint None = uint.MaxValue;

    /// <summary>Creates the scan of a text of <paramref name="n"/> symbols writing into <paramref name="file"/> at <paramref name="start"/> (<c>4 x (n + 1)</c> bytes: the range is checked, the caller creates the file at its final length).</summary>
    /// <exception cref="NotSupportedException">Tagged codec for a text above 2^30 - 1 symbols, or a text above <see cref="PackedTextFormat.MaxSymbols"/>.</exception>
    internal static ChildScan Create(long n, MappedOutputFile file, long start, ChildBuildOptions options) => throw new NotImplementedException();

    /// <summary>The next block of LCP values: ranks <c>firstRank ..</c>, which must continue where the previous call ended (the first call starts at rank 0).</summary>
    /// <exception cref="InvalidOperationException">The ranks are not in order or exceed n + 1.</exception>
    /// <exception cref="InvalidDataException">LCP[0] or LCP[1] is not 0, or a value exceeds the text length.</exception>
    public void Consume(long firstRank, ReadOnlySpan<uint> lcp) => throw new NotImplementedException();

    /// <summary>Closes the intervals still open, closes the output writer and returns what the scan found.</summary>
    /// <exception cref="InvalidOperationException">Not all <c>n + 1</c> values were delivered.</exception>
    internal ChildScanResult Complete() => throw new NotImplementedException();

    public void Dispose() => throw new NotImplementedException();
}
