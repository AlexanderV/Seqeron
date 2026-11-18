#pragma warning disable S2325, S2933, S1144, CA1823, IDE0052, CA1822, S1172, IDE0060, CA1801, S3251, S1186 // STUB: removed when the implementation lands
namespace SuffixTree.Esa.External;

/// <summary>STUB (work package 1b.7d, tests first). Parameters of PH1.</summary>
internal sealed record ExtDcsOptions
{
    internal long BucketCapacity { get; init; } = 1 << 20;

    internal long BaseCase { get; init; } = 16 << 20;

    internal int FanOut { get; init; } = 1024;

    internal int ReservoirSize { get; init; } = 16384;

    internal int Threads { get; init; } = 1;

    internal ulong Seed { get; init; } = 1;

    internal long ChunkSlots { get; init; } = 4 << 20;

    internal int RunSpillBytes { get; init; } = 64 * 1024;

    internal int IoBufferBytes { get; init; } = 64 * 1024;

    internal long MergeBufferBudget { get; init; } = 16 << 20;

    internal int WindowBytes { get; init; } = MappedWindowWriter.DefaultWindowBytes;

    internal long SortLeafCapacity { get; init; } = 1 << 20;

    internal required string ScratchDir { get; init; }

    internal ArenaBudget? Budget { get; init; }

    internal ArenaOptions? ArenaGuards { get; init; }

    internal IoCounters? Counters { get; init; }

    internal MappedFileOptions? Files { get; init; }

    internal Action<string>? Trace { get; init; }
}

/// <summary>STUB (work package 1b.7d, tests first). What one level of the recursion did.</summary>
internal sealed class ExtDcsLevelStats
{
    internal int Level { get; set; }

    internal long TextLength { get; set; }

    internal long SampleSize { get; set; }

    internal long DistinctNames { get; set; }

    internal long NonUniqueRecords { get; set; }

    internal bool NamesUnique => DistinctNames == SampleSize;

    internal bool BaseCase { get; set; }

    internal long Leaves { get; set; }

    internal long SplitScans { get; set; }

    internal long Scans { get; set; }

    internal long FileRuns { get; set; }

    internal double NamingSeconds { get; set; }

    internal double MergeSeconds { get; set; }

    internal double BelowSeconds { get; set; }
}

/// <summary>STUB (work package 1b.7d, tests first). What <see cref="ExternalDcs.BuildRankFile"/> did, level by level.</summary>
internal sealed class ExtDcsStats
{
    internal List<ExtDcsLevelStats> Levels { get; } = [];

    internal int MaxLevel => Levels.Count - 1;

    /// <summary>What the creation of the rank file did about the disk space.</summary>
    internal ReservationOutcome Reservation { get; set; }
}

/// <summary>STUB (work package 1b.7d, tests first). PH1: the rank file of the difference-cover sample.</summary>
internal static class ExternalDcs
{
    internal static ExtDcsStats BuildRankFile(MappedTextReader text, DifferenceCover cover, string rankPath, ExtDcsOptions options) => throw new NotImplementedException();
}
