#pragma warning disable S2325, S1144, S1172, S3459, CS9113 // stub of the tests-first step: every member throws
using SuffixTree.Esa;
using SuffixTree.EmProto.Harness;
using SuffixTree.EmProto.Ph4;

namespace SuffixTree.EmProto.Ph5;

/// <summary>How the single-array child table (layout A) encodes a field: <see cref="Tagged"/> = kind in bits 30-31 (the in-RAM <c>TaggedChild</c>), <see cref="Packed"/> = the raw value (the in-RAM <c>PackedChild</c>).</summary>
internal enum ChildCodec
{
    Tagged,
    Packed,
}

/// <summary>Parameters of PH5 (ADR-004 section 3.2). The output never depends on any of them.</summary>
internal sealed record ChildBuildOptions
{
    internal required string ScratchDir { get; init; }

    /// <summary>Entries per stack block (ADR: 256 K); the two top blocks are in RAM.</summary>
    internal int BlockEntries { get; init; } = 256 * 1024;

    /// <summary>Size of one mapped window of the CHILD, PSV and NSV writers.</summary>
    internal int WindowBytes { get; init; } = 1 << 20;

    /// <summary>Windows the layout A writer keeps mapped at a time.</summary>
    internal int WindowSlots { get; init; } = 8;

    /// <summary>msync a window when it is released (bounds the dirty pages).</summary>
    internal bool FlushWindows { get; init; } = true;

    /// <summary>fsync the spill file after this many unflushed bytes (0 = never).</summary>
    internal long SpillFlushBytes { get; init; } = 16L << 20;

    internal int ReadBufferBytes { get; init; } = 1 << 20;

    internal ChildCodec Codec { get; init; } = ChildCodec.Tagged;

    internal ArenaBudget? Budget { get; init; }

    internal ArenaOptions? ArenaGuards { get; init; }

    internal bool TouchBuffers { get; init; }
}

/// <summary>The files the forward scan writes (null = that part is off): the single-array child table (layout A) and the PSV array of layout F.</summary>
internal sealed record ChildScanOutputs
{
    internal string? ChildPath { get; init; }

    internal string? PsvPath { get; init; }
}

/// <summary>What the forward scan found and did.</summary>
internal sealed class ChildScanResult
{
    internal long N { get; set; }

    internal long InternalCount { get; set; }

    internal long NodeCount { get; set; }

    internal uint LrsDepth { get; set; }

    internal uint LrsRank { get; set; }

    internal long MaxLcp { get; set; }

    internal long StackPeak { get; set; }

    internal SpillStats Spill { get; set; }

    internal WindowStats ChildWindows { get; set; }

    internal WindowStats PsvWindows { get; set; }

    internal double Seconds { get; set; }
}

/// <summary>STUB (tests first): the forward scan of PH5 as a consumer of the PH4 hook (dense LCP values in ascending rank order).</summary>
internal sealed class ChildScan : ILcpChunkConsumer, IDisposable
{
    internal static ChildScan Create(long n, ChildScanOutputs outputs, ChildBuildOptions options) => throw new NotImplementedException();

    public void Consume(long firstRank, ReadOnlySpan<uint> lcp) => throw new NotImplementedException();

    internal ChildScanResult Complete() => throw new NotImplementedException();

    public void Dispose() => throw new NotImplementedException();
}
