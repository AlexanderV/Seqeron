#pragma warning disable S2325, S2933, S1144, CA1823, IDE0052, CA1822, S1172, IDE0060, CA1801, S3251, S1186 // STUB: removed when the implementation lands
namespace SuffixTree.Esa.External;

/// <summary>STUB (work package 1b.7d, tests first). Parameters of PH2 / PH3.</summary>
internal sealed record SortOptions
{
    internal long LeafCapacity { get; init; } = 1 << 20;

    internal int FanOut { get; init; } = 64;

    internal int Oversample { get; init; } = 16;

    internal int Threads { get; init; } = 1;

    internal ulong Seed { get; init; } = 1;

    internal required string ScratchDir { get; init; }

    internal int SpillBytes { get; init; } = 64 * 1024;

    internal ArenaBudget? Budget { get; init; }

    internal ArenaOptions? ArenaGuards { get; init; }

    internal IoCounters? Counters { get; init; }

    internal Action<string>? Trace { get; init; }
}

/// <summary>STUB (work package 1b.7d, tests first). What a sort did.</summary>
internal sealed class SortStats
{
    internal SortWork Work { get; } = new();

    internal long Buckets { get; set; }

    internal long Leaves { get; set; }

    internal long Partitions { get; set; }

    internal int MaxDepth { get; set; }

    internal long MaxLeaf { get; set; }

    internal long Retries { get; set; }

    internal double PartitionSeconds { get; set; }

    internal double LeafSeconds { get; set; }
}

/// <summary>STUB (work package 1b.7d, tests first). PH2 and PH3: the partition into lexicographic buckets and the bucket sort.</summary>
internal static class BucketSorter
{
    internal static SortStats Run(MappedTextReader text, SortOptions options, ISaSliceSink sink, RankTable ranks) => throw new NotImplementedException();
}
