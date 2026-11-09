using SuffixTree.EmProto.Harness;

namespace SuffixTree.EmProto.Build;

/// <summary>How a run of positions with equal 64-bit keys is ordered by the comparator.</summary>
internal enum TieSort
{
    /// <summary>Stable merge sort with binary insertion for tiny runs (default; fewest comparator calls, n / 2 per merge on reverse-ordered runs).</summary>
    Merge,

    /// <summary>The introsort of <see cref="MemoryExtensions.Sort{T, TComparer}(Span{T}, TComparer)"/> (the version of the first E1 measurements).</summary>
    Intro,
}

/// <summary>Parameters of the bucket partition and the bucket sort (ADR-004 sections 4.2 and 5; the fuzz tests scale them down to K = 3..50, F = 2..4, X = 7..64).</summary>
internal sealed record BuildOptions
{
    /// <summary>K: records per leaf bucket (a bucket of at most K positions is sorted in RAM).</summary>
    internal long LeafCapacity { get; init; } = 1 << 20;

    /// <summary>F: largest number of children of one partition step (the fan-out of ADR-004: 64).</summary>
    internal int FanOut { get; init; } = 64;

    /// <summary>Random sample positions per child when splitters are chosen.</summary>
    internal int Oversample { get; init; } = 16;

    /// <summary>Worker threads (1 = the caller's thread).</summary>
    internal int Threads { get; init; } = 1;

    /// <summary>Seed of the splitter sampling; the output does not depend on it, the shape of the partition does.</summary>
    internal ulong Seed { get; init; } = 1;

    /// <summary>Directory for the position files of the buckets; null keeps every bucket in memory (E1, E2).</summary>
    internal string? ScratchDir { get; init; }

    /// <summary>Write buffer per child bucket in file mode (a bucket that fits stays in memory).</summary>
    internal int SpillBytes { get; init; } = 64 * 1024;

    /// <summary>Ordering of key-tie runs in the leaf buckets.</summary>
    internal TieSort TieSort { get; init; } = TieSort.Merge;

    /// <summary>Shared cap of the native working buffers (null = unlimited).</summary>
    internal ArenaBudget? Budget { get; init; }

    /// <summary>Reduced strings of the DCS recursion with at most this many symbols are sorted by SA-IS (N_base of ADR-004: 16 M).</summary>
    internal long DcsBaseCase { get; init; } = 16_000_000;
}

/// <summary>What a build did, for the report.</summary>
internal sealed class BuildStats
{
    internal WorkCounters Work { get; } = new();

    internal long Buckets { get; set; }

    internal long Leaves { get; set; }

    internal long Partitions { get; set; }

    internal int MaxDepth { get; set; }

    internal long MaxLeaf { get; set; }

    internal long Retries { get; set; }

    internal double PartitionSeconds { get; set; }

    internal double LeafSeconds { get; set; }
}
