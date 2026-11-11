using SuffixTree.EmProto.Harness;
using SuffixTree.EmProto.Text;

namespace SuffixTree.EmProto.Dcs;

/// <summary>
/// Parameters of the external-memory PH1 (ADR-004 section 4.2, prototype "E3b"): everything that sets a buffer size, a bucket capacity or a recursion threshold. The output never depends on
/// any of them (ADR-004 A8, owner answer Q5); the fuzz tests scale them down to K = 3..50, N_base = 4..200, X = 7..64 so that every multi-bucket and multi-level path runs on small texts.
/// </summary>
internal sealed record ExtDcsOptions
{
    /// <summary>K_s: records (sample positions) per sample bucket; a bucket of at most K_s records is sorted in RAM.</summary>
    internal long BucketCapacity { get; init; } = 1 << 20;

    /// <summary>N_base: a reduced string of at most this many symbols is sorted by SA-IS in RAM; longer ones go through the external levels.</summary>
    internal long BaseCase { get; init; } = 4 << 20;

    /// <summary>F_max: largest number of sample buckets one split creates.</summary>
    internal int FanOut { get; init; } = 1024;

    /// <summary>Records kept by the reservoir that chooses the splitters of one split.</summary>
    internal int ReservoirSize { get; init; } = 16384;

    /// <summary>Worker threads of the bucket tasks and of the suffix sort of the external levels.</summary>
    internal int Threads { get; init; } = 1;

    /// <summary>Seed of the reservoir sampling; the output does not depend on it, the shape of the bucket tree does.</summary>
    internal ulong Seed { get; init; } = 1;

    /// <summary>R: slots per merge window and per chunk of the suffix-array inversion (the staging array holds R four-byte entries).</summary>
    internal long ChunkSlots { get; init; } = 4 << 20;

    /// <summary>A bucket run of at most this many bytes stays in memory, larger runs go to a scratch file.</summary>
    internal int RunSpillBytes { get; init; } = 64 * 1024;

    /// <summary>Size of one read or write buffer of the merge, of the tuple files and of the class cursors.</summary>
    internal int IoBufferBytes { get; init; } = 64 * 1024;

    /// <summary>Size of one sliding window of the mapped rank file writer.</summary>
    internal int WindowBytes { get; init; } = 1 << 20;

    /// <summary>K of the suffix sort of a reduced string at the external levels (records per leaf bucket of the position sort).</summary>
    internal long SortLeafCapacity { get; init; } = 1 << 20;

    /// <summary>Directory of the scratch files (created when missing; every file of a build is removed when it ends).</summary>
    internal required string ScratchDir { get; init; }

    /// <summary>Shared cap of the native buffers (null = unlimited).</summary>
    internal ArenaBudget? Budget { get; init; }
}

/// <summary>What one level of the recursion did.</summary>
internal sealed class ExtDcsLevelStats
{
    /// <summary>0 = the input text, 1 = the reduced string of level 0, ...</summary>
    internal int Level { get; set; }

    /// <summary>Symbols of the text of this level.</summary>
    internal long TextLength { get; set; }

    /// <summary>Sample positions of this level.</summary>
    internal long SampleSize { get; set; }

    /// <summary>Distinct names (distinct X-prefixes).</summary>
    internal long DistinctNames { get; set; }

    /// <summary>Sample positions whose name is shared with at least one other sample position.</summary>
    internal long NonUniqueRecords { get; set; }

    /// <summary>All names distinct: the names are the ranks and there is no level below.</summary>
    internal bool NamesUnique => DistinctNames == SampleSize;

    /// <summary>The reduced string of this level was sorted by SA-IS in RAM (the recursion ended here).</summary>
    internal bool BaseCase { get; set; }

    /// <summary>Leaf buckets that were sorted.</summary>
    internal long Leaves { get; set; }

    /// <summary>Scans of the sample positions by tasks that split a bucket instead of sorting it (the root and the buckets that overflowed).</summary>
    internal long SplitScans { get; set; }

    /// <summary>Scans of the sample positions by all tasks (one per leaf, one per split, one per overflow).</summary>
    internal long Scans { get; set; }

    /// <summary>Leaves whose run went to a scratch file.</summary>
    internal long FileRuns { get; set; }

    /// <summary>Seconds of the bucket tasks, of the merge into the sink, and of everything below this level.</summary>
    internal double NamingSeconds { get; set; }

    internal double MergeSeconds { get; set; }

    internal double BelowSeconds { get; set; }
}

/// <summary>What <see cref="ExternalDcs.BuildRankFile"/> did, level by level.</summary>
internal sealed class ExtDcsStats
{
    internal List<ExtDcsLevelStats> Levels { get; } = [];

    /// <summary>Index of the deepest level reached.</summary>
    internal int MaxLevel => Levels.Count - 1;
}

/// <summary>
/// PH1 of ADR-004 built in external memory (prototype "E3b"): the rank file <c>rank.L0[slot]</c> of the difference-cover sample, byte for byte the table of <see cref="DcsBuilder"/>.
/// </summary>
internal static class ExternalDcs
{
    /// <summary>Builds the rank file of <paramref name="text"/> at <paramref name="rankPath"/> (length <c>Slots * 4</c>), written through bounded mapped windows.</summary>
    internal static ExtDcsStats BuildRankFile(PackedText text, DifferenceCover cover, string rankPath, ExtDcsOptions options)
        => throw new NotImplementedException("E3b stub: the external-memory PH1 is written after its tests");
}

/// <summary>Counters of the writes through mapped windows (the explicit pread / pwrite calls are counted by <see cref="IoCounters"/>).</summary>
internal static class MappedIo
{
    private static long s_writeBytes;
    private static long s_windows;
    private static long s_maxWindow;

    internal static long WriteBytes => Interlocked.Read(ref s_writeBytes);

    /// <summary>Windows (views) created.</summary>
    internal static long Windows => Interlocked.Read(ref s_windows);

    /// <summary>Largest window created, in bytes.</summary>
    internal static long MaxWindowBytes => Interlocked.Read(ref s_maxWindow);

    internal static void Reset()
    {
        Interlocked.Exchange(ref s_writeBytes, 0);
        Interlocked.Exchange(ref s_windows, 0);
        Interlocked.Exchange(ref s_maxWindow, 0);
    }

    internal static void AddWrite(long bytes) => Interlocked.Add(ref s_writeBytes, bytes);

    internal static void AddWindow(long bytes)
    {
        Interlocked.Increment(ref s_windows);
        Atomic.RaiseTo(ref s_maxWindow, bytes);
    }
}
