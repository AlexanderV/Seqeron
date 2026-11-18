#pragma warning disable S2325, S2933, S1144, CA1823, IDE0052, CA1822, S1172, IDE0060, CA1801, S3251, S1186 // STUB: removed when the implementation lands
namespace SuffixTree.Esa.External;

/// <summary>STUB (work package 1b.7d, tests first). Explicit values for the sizes that the rule would otherwise derive (null = derive).</summary>
internal sealed record SaBuildTuning
{
    internal long? DcsBucketCapacity { get; init; }

    internal long? DcsBaseCase { get; init; }

    internal long? DcsChunkSlots { get; init; }

    internal long? DcsSortCapacity { get; init; }

    internal int? DcsRunSpillBytes { get; init; }

    internal int? DcsFanOut { get; init; }

    internal int? DcsReservoir { get; init; }

    internal int? DcsIoBytes { get; init; }

    internal long? DcsMergeBudget { get; init; }

    internal long? BucketCapacity { get; init; }

    internal int? FanOut { get; init; }

    internal int? Oversample { get; init; }

    internal int? SpillBytes { get; init; }

    internal int? WindowBytes { get; init; }
}

/// <summary>STUB (work package 1b.7d, tests first). What <see cref="SaBuilder.Build"/> needs besides the files.</summary>
internal sealed record SaBuildOptions
{
    internal required string ScratchDir { get; init; }

    internal int Threads { get; init; } = 4;

    internal int Period { get; init; } = 1024;

    internal bool ExactPeriod { get; init; }

    internal long ScaleN { get; init; }

    internal long RuleN { get; init; }

    internal long ArenaCapBytes { get; init; } = 160L << 20;

    internal SaBuildTuning Tuning { get; init; } = new();

    internal ulong Seed { get; init; } = 1;

    internal IoCounters? Counters { get; init; }

    internal MappedFileOptions? Files { get; init; }

    internal ArenaOptions? ArenaGuards { get; init; }

    internal bool SyncOutput { get; init; }

    internal bool VerifyRanks { get; init; }

    internal string? RankPath { get; init; }

    internal bool MeasureRss { get; init; } = true;

    internal IProcSource? ProcSource { get; init; }

    internal double RssAnonCapMb { get; init; } = 320;

    internal Action<string>? Trace { get; init; }
}

/// <summary>STUB (work package 1b.7d, tests first). Every size of one run.</summary>
internal sealed record SaBuildParameters
{
    internal long N { get; init; }

    internal double Factor { get; init; }

    internal bool Scaled { get; init; }

    internal long DcsBucketCapacity { get; init; }

    internal long DcsBaseCase { get; init; }

    internal long DcsChunkSlots { get; init; }

    internal long DcsSortCapacity { get; init; }

    internal int DcsRunSpillBytes { get; init; }

    internal int DcsFanOut { get; init; }

    internal int DcsReservoir { get; init; }

    internal int DcsIoBytes { get; init; }

    internal long DcsMergeBudget { get; init; }

    internal long BucketCapacity { get; init; }

    internal int FanOut { get; init; }

    internal int Oversample { get; init; }

    internal int SpillBytes { get; init; }

    internal int WindowBytes { get; init; }
}

/// <summary>STUB (work package 1b.7d, tests first). The scale-down rule of ADR-004 (E3, E3b, E6) for PH1 and PH2 / PH3.</summary>
internal static class SaBuildRules
{
    internal const long ReferenceN = 3_000_000_000;

    internal const long BucketReference = 1_000_000;

    internal const long BucketUnscaled = 1 << 20;

    internal const long DcsBucketReference = 1 << 20;

    internal const long DcsBaseReference = 16 << 20;

    internal const long DcsChunkReference = 4 << 20;

    internal const long DcsSortReference = 1 << 20;

    internal const long RunSpillReference = 64 * 1024;

    internal static SaBuildParameters Resolve(long n, SaBuildOptions options) => throw new NotImplementedException();
}

/// <summary>STUB (work package 1b.7d, tests first). What a finished build knows.</summary>
internal sealed class SaBuildResult
{
    internal required SaBuildParameters Parameters { get; init; }

    internal long N { get; init; }

    internal int Period { get; init; }

    internal int CoverResidues { get; init; }

    internal CoverKind CoverKind { get; init; }

    internal required ExtDcsStats Dcs { get; init; }

    internal required SortStats Sort { get; init; }

    internal IoSnapshot Io { get; init; }

    internal IoSnapshot Ph1Io { get; init; }

    internal IoSnapshot Ph23Io { get; init; }

    internal double Ph1Seconds { get; init; }

    internal double Ph23Seconds { get; init; }

    internal double Seconds => Ph1Seconds + Ph23Seconds;

    internal long Ph1ArenaPeakBytes { get; init; }

    internal long Ph23ArenaPeakBytes { get; init; }

    internal long RankFileBytes { get; init; }

    internal long Ph1ScratchPeakBytes { get; init; }

    internal long Ph23ScratchPeakBytes { get; init; }

    internal double? RssAnonPeakMb { get; init; }

    internal double? Ph1RssAnonPeakMb { get; init; }

    internal double? Ph23RssAnonPeakMb { get; init; }

    internal RssVerdict RssVerdict { get; init; }

    internal ReservationOutcome Reservation { get; init; }

    /// <summary>Levels of the recursion below the text (0 = the names were unique or the first reduced string fitted the base case).</summary>
    internal int RecursionLevels => Dcs.MaxLevel;
}

/// <summary>STUB (work package 1b.7d, tests first). PH1 and PH2 / PH3 of ADR-004 on a packed text: writes the suffix array file.</summary>
internal static class SaBuilder
{
    internal static SaBuildResult Build(MappedTextReader text, string saPath, SaBuildOptions options) => throw new NotImplementedException();
}
