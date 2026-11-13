using SuffixTree.Esa;
using SuffixTree.EmProto.Experiments;
using SuffixTree.EmProto.Harness;
using SuffixTree.EmProto.Ph4;
using SuffixTree.EmProto.Ph5;

namespace SuffixTree.EmProto.Pipeline;

/// <summary>
/// Explicit values for the sizes that the scale-down rule would otherwise derive (null = derive). The tests set every one of them to tiny values (K = 3..50, F = 2..4, X = 7..64, windows of a few bytes) so that
/// every multi-bucket, multi-level, multi-chunk and spill path runs on texts of 10^3-10^5 symbols; the output never depends on any of them (ADR-004 section 5, "Oracles per stage").
/// </summary>
internal sealed record EndToEndTuning
{
    /// <summary>PH2 / PH3: records per leaf bucket (K).</summary>
    internal long? BucketCapacity { get; init; }

    /// <summary>PH2: fan-out of one partition step (F).</summary>
    internal int? FanOut { get; init; }

    /// <summary>PH2: random sample positions per child when splitters are chosen.</summary>
    internal int? Oversample { get; init; }

    /// <summary>PH2: write buffer per child bucket.</summary>
    internal int? SpillBytes { get; init; }

    /// <summary>PH1: records per sample bucket (K_s).</summary>
    internal long? DcsBucketCapacity { get; init; }

    /// <summary>PH1: reduced strings of at most this many symbols are sorted by SA-IS in RAM (N_base).</summary>
    internal long? DcsBaseCase { get; init; }

    /// <summary>PH1: slots per merge window and per inversion chunk (R of PH1).</summary>
    internal long? DcsChunkSlots { get; init; }

    /// <summary>PH1: K of the suffix sort of a reduced string.</summary>
    internal long? DcsSortCapacity { get; init; }

    /// <summary>PH1: a bucket run of at most this many bytes stays in memory.</summary>
    internal int? DcsRunSpillBytes { get; init; }

    /// <summary>PH1: fan-out of a split of the sample buckets.</summary>
    internal int? DcsFanOut { get; init; }

    /// <summary>PH1: reservoir of a split.</summary>
    internal int? DcsReservoir { get; init; }

    /// <summary>PH1: size of a read or write buffer of the merge and of the tuple files.</summary>
    internal int? DcsIoBytes { get; init; }

    /// <summary>PH1: total of the read buffers of the merge.</summary>
    internal long? DcsMergeBudget { get; init; }

    /// <summary>PH4: R, positions per Phi chunk and ranks per LCP chunk.</summary>
    internal long? LcpChunk { get; init; }

    /// <summary>PH4: parallel B2 tasks (null = what the arena cap allows at the reference geometry).</summary>
    internal int? LcpWorkers { get; init; }

    /// <summary>PH4: sequential read buffer.</summary>
    internal int? LcpReadBytes { get; init; }

    /// <summary>PH4: write buffer of the tuple distributions.</summary>
    internal int? LcpTupleBytes { get; init; }

    /// <summary>PH5: entries per block of the stack.</summary>
    internal int? BlockEntries { get; init; }

    /// <summary>PH5: window of the CHILD writer (a power of two).</summary>
    internal int? ChildWindowBytes { get; init; }

    /// <summary>PH5: windows of the CHILD writer.</summary>
    internal int? ChildSlots { get; init; }

    /// <summary>PH5: fsync of the spill file after this many bytes.</summary>
    internal long? SpillFlushBytes { get; init; }

    /// <summary>PH0 / PH4 / PH1: size of a sliding window of the mapped writers.</summary>
    internal int? WindowBytes { get; init; }

    /// <summary>PH0: symbols per block of the streaming pack (a multiple of 8).</summary>
    internal int? PrepareBlockSymbols { get; init; }
}

/// <summary>What <see cref="EndToEnd.Run"/> needs: the input text (one byte per symbol), the prefix of the output files, a scratch directory and the sizes.</summary>
internal sealed record EndToEndOptions
{
    /// <summary>The input text file: one byte per symbol, no terminator (the benchmark inputs); the distinct byte values become the symbols <c>0..sigma-1</c> in increasing byte order.</summary>
    internal required string InputPath { get; init; }

    /// <summary>The outputs are <c>&lt;prefix&gt;.pk</c> (+ <c>.pk.meta</c>), <c>.sa</c>, <c>.lcp</c>, <c>.ovf</c>, <c>.cld</c> and <c>.stats</c>.</summary>
    internal required string OutPrefix { get; init; }

    /// <summary>The scratch directory (created; deleted with everything in it when the run ends). Must not lie under <see cref="OutPrefix"/> as a string prefix, or the dirty-page view of the outputs would count it.</summary>
    internal required string ScratchDir { get; init; }

    /// <summary>Worker threads of the bucket sorts and of the B2 tasks.</summary>
    internal int Threads { get; init; } = 4;

    /// <summary>
    /// The scale-down emulation (ADR-004 section 5, E3 / E3b / E4 / E5): when set (3,000,000,000 for the 100 Mbp runs), every size that decides the number of buckets, levels or chunks is its reference value
    /// for 3x10^9 symbols multiplied by n / ScaleN; 0 = unscaled (the reference values themselves).
    /// </summary>
    internal long ScaleN { get; init; }

    /// <summary>X of the difference cover (nominal; a Singer period nearby is used unless <see cref="ExactPeriod"/>).</summary>
    internal int Period { get; init; } = 1024;

    internal bool ExactPeriod { get; init; }

    internal ChildCodec Codec { get; init; } = ChildCodec.Packed;

    /// <summary>Cap of the native buffers of each of PH1, PH3 and PH4 (the stages run one after the other).</summary>
    internal long ArenaCapBytes { get; init; } = 160L << 20;

    /// <summary>Cap of the native buffers of PH5 (it lives next to PH4 in the fused phase).</summary>
    internal long ChildArenaCapBytes { get; init; } = 64L << 20;

    /// <summary>fsync every output at the end (reported as its own phase, not counted in the stage times).</summary>
    internal bool SyncOutputs { get; init; } = true;

    /// <summary>Makes every native buffer resident when it is allocated (the worst case of the working set, for the RssAnon figures of the scale-down runs).</summary>
    internal bool TouchBuffers { get; init; }

    internal EndToEndTuning Tuning { get; init; } = new();

    /// <summary>Test hook: called as <c>("begin" | "end", phase)</c> around every phase; an exception thrown here fails the run like a failing stage.</summary>
    internal Action<string, string>? Trace { get; init; }

    /// <summary>Test hook: options of the arenas of PH4 and PH5 (guard bands, poison pattern); null = production arenas.</summary>
    internal ArenaOptions? ArenaGuards { get; init; }
}

/// <summary>Every size of one run, resolved from <see cref="EndToEndOptions"/> and the text length (what the report prints, and what the tests pin).</summary>
internal sealed record EndToEndParameters
{
    internal long N { get; init; }

    /// <summary>n / ScaleN, or 1 when unscaled.</summary>
    internal double Factor { get; init; }

    internal bool Scaled { get; init; }

    // PH1 (E3b)
    internal long DcsBucketCapacity { get; init; }

    internal long DcsBaseCase { get; init; }

    internal long DcsChunkSlots { get; init; }

    internal long DcsSortCapacity { get; init; }

    internal int DcsRunSpillBytes { get; init; }

    internal int DcsFanOut { get; init; }

    internal int DcsReservoir { get; init; }

    internal int DcsIoBytes { get; init; }

    internal long DcsMergeBudget { get; init; }

    // PH2 / PH3 (E3)
    internal long BucketCapacity { get; init; }

    internal int FanOut { get; init; }

    internal int Oversample { get; init; }

    internal int SpillBytes { get; init; }

    // PH4 (E4)
    internal long LcpChunk { get; init; }

    internal int LcpWorkers { get; init; }

    internal int LcpReadBytes { get; init; }

    internal int LcpTupleBytes { get; init; }

    // PH5 (E5)
    internal int BlockEntries { get; init; }

    internal int ChildWindowBytes { get; init; }

    internal int ChildSlots { get; init; }

    internal long SpillFlushBytes { get; init; }

    // shared
    internal int WindowBytes { get; init; }

    internal int PrepareBlockSymbols { get; init; }

    public override string ToString() => string.Create(System.Globalization.CultureInfo.InvariantCulture,
        $"n={N} f={Factor:E3} scaled={Scaled}: PH1 K_s={DcsBucketCapacity} N_base={DcsBaseCase} R={DcsChunkSlots} K_sort={DcsSortCapacity} run-spill={DcsRunSpillBytes}; PH3 K={BucketCapacity} F={FanOut}; PH4 R={LcpChunk} P={(N + LcpChunk - 1) / LcpChunk} B2 workers={LcpWorkers}; PH5 block={BlockEntries} window={ChildWindowBytes} slots={ChildSlots}");
}

/// <summary>
/// The scale-down rule of E3, E3b and E4 applied to every stage (read from their commands, not invented): with <c>f = n / ScaleN</c> PH1 takes K_s, N_base, R and K_sort as their 3x10^9 reference values times f
/// (<c>Math.Round</c>, floors 3 / 4 / 16 / 3) and the run-spill bytes as 64 KiB x f (E3b); PH3 takes K = n / 3000 for ScaleN = 3x10^9, i.e. <c>n x 10^6 / ScaleN</c>, because the reference K is 10^6 records (E3 <c>--k-div 3000</c>);
/// PH4 takes R = <c>ceil(8,000,000 x f)</c> with a floor of 3 and the number of parallel B2 tasks that the arena cap allows at the reference geometry (E4); PH5 has no size that depends on n (E5). Unscaled, every
/// stage uses its reference value (PH3: 2^20, the default of <see cref="Build.BuildOptions"/>). An explicit <see cref="EndToEndTuning"/> value wins.
/// </summary>
internal static class ScaleRules
{
    /// <summary>The reference K of PH3 at 3x10^9 symbols: 10^6 records (3x10^9 / 3000).</summary>
    internal const long BucketReference = 1_000_000;

    /// <summary>The default K of an unscaled run (<see cref="Build.BuildOptions.LeafCapacity"/>).</summary>
    internal const long BucketUnscaled = 1 << 20;

    /// <summary>The reference run-spill bytes of PH1 (E3b <c>--run-spill-kb 64</c>).</summary>
    internal const long RunSpillReference = 64 * 1024;

    internal static EndToEndParameters Resolve(long n, EndToEndOptions o)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(n);
        ArgumentOutOfRangeException.ThrowIfNegative(o.ScaleN);
        EndToEndTuning t = o.Tuning;
        bool scaled = o.ScaleN > 0;
        double f = scaled ? n / (double)o.ScaleN : 1.0;
        long Scale(long reference, long floor) => scaled ? Math.Max(floor, (long)Math.Round(reference * f)) : reference;

        // PH4: workers from the reference geometry (R = 8,000,000 at 3x10^9), as E4 and the fused E5 do
        long chunk = t.LcpChunk ?? (scaled ? Math.Max(3, (long)Math.Ceiling(E4.ChunkRef * f)) : E4.ChunkRef);
        var reference = new ExtLcpOptions { ScratchDir = string.Empty, ChunkPositions = E4.ChunkRef, Threads = o.Threads, Mode = LcpMode.Distributed };
        long perTask = ExternalLcp.WorkerBytes(E4.ReferenceN, reference);
        int workers = t.LcpWorkers ?? (int)Math.Clamp(o.ArenaCapBytes / Math.Max(1, perTask), 1, o.Threads);

        return new EndToEndParameters
        {
            N = n,
            Factor = f,
            Scaled = scaled,
            DcsBucketCapacity = Math.Max(3, t.DcsBucketCapacity ?? Scale(E3b.KsRef, 3)),
            DcsBaseCase = Math.Max(4, t.DcsBaseCase ?? Scale(E3b.BaseRef, 4)),
            DcsChunkSlots = Math.Max(16, t.DcsChunkSlots ?? Scale(E3b.ChunkRef, 16)),
            DcsSortCapacity = Math.Max(3, t.DcsSortCapacity ?? Scale(E3b.SortKRef, 3)),
            DcsRunSpillBytes = t.DcsRunSpillBytes ?? (int)Scale(RunSpillReference, 0),
            DcsFanOut = t.DcsFanOut ?? 1024,
            DcsReservoir = t.DcsReservoir ?? 16384,
            DcsIoBytes = t.DcsIoBytes ?? 64 * 1024,
            DcsMergeBudget = t.DcsMergeBudget ?? 16L << 20,
            BucketCapacity = Math.Max(2, t.BucketCapacity ?? (scaled ? (long)((decimal)n * BucketReference / o.ScaleN) : BucketUnscaled)),
            FanOut = t.FanOut ?? 64,
            Oversample = t.Oversample ?? 16,
            SpillBytes = t.SpillBytes ?? 64 * 1024,
            LcpChunk = chunk,
            LcpWorkers = workers,
            LcpReadBytes = t.LcpReadBytes ?? 1 << 20,
            LcpTupleBytes = t.LcpTupleBytes ?? 32 * 1024,
            BlockEntries = t.BlockEntries ?? 256 * 1024,
            ChildWindowBytes = t.ChildWindowBytes ?? 1 << 20,
            ChildSlots = t.ChildSlots ?? 8,
            SpillFlushBytes = t.SpillFlushBytes ?? 16L << 20,
            WindowBytes = t.WindowBytes ?? 1 << 20,
            PrepareBlockSymbols = t.PrepareBlockSymbols ?? 1 << 20,
        };
    }
}
