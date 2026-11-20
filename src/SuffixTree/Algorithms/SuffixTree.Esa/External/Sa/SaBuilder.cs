using System.Diagnostics;

namespace SuffixTree.Esa.External;

/// <summary>
/// Explicit values for the sizes that the rule would otherwise derive (null = derive). The tests set every one of them to tiny values (K = 3..50, F = 2..4, X = 7..64, windows of a few bytes) so that every multi-bucket,
/// multi-level, spill and window path runs on texts of 10^1-10^4 symbols; the output never depends on any of them (ADR-004 section 5, "Oracles per stage").
/// </summary>
internal sealed record SaBuildTuning
{
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

    /// <summary>PH2 / PH3: records per leaf bucket (K).</summary>
    internal long? BucketCapacity { get; init; }

    /// <summary>PH2: fan-out of one partition step (F).</summary>
    internal int? FanOut { get; init; }

    /// <summary>PH2: random sample positions per child when splitters are chosen.</summary>
    internal int? Oversample { get; init; }

    /// <summary>PH2: write buffer per child bucket.</summary>
    internal int? SpillBytes { get; init; }

    /// <summary>PH1 / PH3: size of a sliding window of the mapped writers (the suffix array, the rank file).</summary>
    internal int? WindowBytes { get; init; }
}

/// <summary>What <see cref="SaBuilder.Build(MappedTextReader, string, SaBuildOptions)"/> needs besides the files.</summary>
internal sealed record SaBuildOptions
{
    /// <summary>The scratch directory (created when missing; <c>ph1</c> and <c>ph3</c> inside it and the rank file are created and removed by the build, the directory itself is removed when this build created it and it is empty).</summary>
    internal required string ScratchDir { get; init; }

    /// <summary>Worker threads of the bucket tasks of PH1 and of PH3 (the calling thread is one of them).</summary>
    internal int Threads { get; init; } = 4;

    /// <summary>X of the difference cover (nominal; a Singer period nearby is used unless <see cref="ExactPeriod"/>: 993 for the default 1024).</summary>
    internal int Period { get; init; } = 1024;

    /// <summary>Use <see cref="Period"/> itself (the pruned cover where it is not a Singer period).</summary>
    internal bool ExactPeriod { get; init; }

    /// <summary>
    /// The scale-down emulation of ADR-004 section 5 (E3, E3b, E6): when set (3,000,000,000 for the 100 Mbp runs), every size that decides the number of buckets, levels or chunks is its reference value for 3x10^9 symbols multiplied by
    /// n / ScaleN; 0 = the reference values themselves.
    /// </summary>
    internal long ScaleN { get; init; }

    /// <summary>E7 series A: the n at which the scale-down rule is evaluated (0 = the text length). Ignored when <see cref="ScaleN"/> is 0.</summary>
    internal long RuleN { get; init; }

    /// <summary>Cap of the native buffers of each of PH1 and PH2 / PH3 (the two phases run one after the other; owner decision: 160 MiB).</summary>
    internal long ArenaCapBytes { get; init; } = 160L << 20;

    internal SaBuildTuning Tuning { get; init; } = new();

    /// <summary>Seed of the sampling; the output does not depend on it.</summary>
    internal ulong Seed { get; init; } = 1;

    /// <summary>The counters of the build (a private instance when null; the result carries a snapshot either way).</summary>
    internal IoCounters? Counters { get; init; }

    /// <summary>How the output file, the rank file and the scratch files are created (reservation of the disk space); null = <see cref="ReservationMode.Auto"/>.</summary>
    internal MappedFileOptions? Files { get; init; }

    /// <summary>Test hook: options of the arenas (guard bands, poison pattern); null = production arenas.</summary>
    internal ArenaOptions? ArenaGuards { get; init; }

    /// <summary>fsync the suffix array when the build is done (not part of the stage time of the prototype).</summary>
    internal bool SyncOutput { get; init; }

    /// <summary>Read every slot of the rank file once between PH1 and PH2 / PH3 and check it (distinct ranks, the unused slots zero): an O(S) pass that the prototype does not make.</summary>
    internal bool VerifyRanks { get; init; }

    /// <summary>The rank file (scratch: created by PH1, consumed by PH3, deleted when PH3 is done); null = <c>ranks.bin</c> in the scratch directory.</summary>
    internal string? RankPath { get; init; }

    /// <summary>Sample RssAnon while the build runs (the figure of the gate); unavailable where there is no <c>/proc</c>.</summary>
    internal bool MeasureRss { get; init; } = true;

    /// <summary>Where the kernel's text files come from (null = the real <c>/proc</c>); a test passes a source that fails like a platform without one.</summary>
    internal IProcSource? ProcSource { get; init; }

    /// <summary>The gate of the owner decision of 2026-11-18: RssAnon of the process at most 320 MB (10^6 bytes) at every n with the default parameters and four threads.</summary>
    internal double RssAnonCapMb { get; init; } = 320;

    /// <summary>
    /// Test hook, called with <c>"ph1"</c> and <c>"ph23"</c> when a phase begins, <c>"verify"</c> before the rank file check, and <c>"naming"</c>, <c>"merge"</c>, <c>"rank-written"</c>, <c>"partition"</c>, <c>"scan"</c>, <c>"handover"</c> and <c>"leaf"</c> from the workers; an exception thrown
    /// here fails the build like a failing step.
    /// </summary>
    internal Action<string>? Trace { get; init; }
}

/// <summary>Every size of one run, resolved from <see cref="SaBuildOptions"/> and the text length (what the report prints, and what the tests pin).</summary>
internal sealed record SaBuildParameters
{
    internal long N { get; init; }

    /// <summary>n / ScaleN (or RuleN / ScaleN), or 1 when unscaled.</summary>
    internal double Factor { get; init; }

    internal bool Scaled { get; init; }

    // PH1
    internal long DcsBucketCapacity { get; init; }

    internal long DcsBaseCase { get; init; }

    internal long DcsChunkSlots { get; init; }

    internal long DcsSortCapacity { get; init; }

    internal int DcsRunSpillBytes { get; init; }

    internal int DcsFanOut { get; init; }

    internal int DcsReservoir { get; init; }

    internal int DcsIoBytes { get; init; }

    internal long DcsMergeBudget { get; init; }

    // PH2 / PH3
    internal long BucketCapacity { get; init; }

    internal int FanOut { get; init; }

    internal int Oversample { get; init; }

    internal int SpillBytes { get; init; }

    // shared
    internal int WindowBytes { get; init; }

    public override string ToString() => string.Create(System.Globalization.CultureInfo.InvariantCulture,
        $"n={N} f={Factor:E3} scaled={Scaled}: PH1 K_s={DcsBucketCapacity} N_base={DcsBaseCase} R={DcsChunkSlots} K_sort={DcsSortCapacity} run-spill={DcsRunSpillBytes}; PH3 K={BucketCapacity} F={FanOut} window={WindowBytes}");
}

/// <summary>
/// The scale-down rule of E3, E3b and E6 for PH1 and PH2 / PH3 (read from the prototype's <c>ScaleRules</c>, not invented): with <c>f = n / ScaleN</c> PH1 takes K_s, N_base, R and K_sort as their 3x10^9 reference values times f
/// (<c>Math.Round</c>, floors 3 / 4 / 16 / 3) and the run-spill bytes as 64 KiB x f (E3b); PH3 takes K = n x 10^6 / ScaleN, because the reference K is 10^6 records (E3: n / 3000 for ScaleN = 3x10^9). Unscaled, every stage uses
/// its reference value (PH3: 2^20; N_base: 16 Mi symbols, the owner's decision of 2026-11-18). An explicit value of <see cref="SaBuildOptions.Tuning"/> wins.
/// </summary>
internal static class SaBuildRules
{
    /// <summary>The reference text length of the sizes: 3x10^9 symbols.</summary>
    internal const long ReferenceN = 3_000_000_000;

    /// <summary>The reference K of PH3 at 3x10^9 symbols: 10^6 records (3x10^9 / 3000).</summary>
    internal const long BucketReference = 1_000_000;

    /// <summary>The K of an unscaled run.</summary>
    internal const long BucketUnscaled = 1 << 20;

    /// <summary>The reference K_s of PH1.</summary>
    internal const long DcsBucketReference = 1 << 20;

    /// <summary>The reference N_base of PH1: 16 Mi symbols.</summary>
    internal const long DcsBaseReference = 16 << 20;

    /// <summary>The reference R of PH1.</summary>
    internal const long DcsChunkReference = 4 << 20;

    /// <summary>The reference K of the suffix sort of a reduced string.</summary>
    internal const long DcsSortReference = 1 << 20;

    /// <summary>The reference run-spill bytes of PH1 (E3b <c>--run-spill-kb 64</c>).</summary>
    internal const long RunSpillReference = 64 * 1024;

    /// <summary>The sizes of a run.</summary>
    /// <exception cref="ArgumentOutOfRangeException">A negative n, scale or rule, fewer than one thread, or an explicit value that its stage cannot use.</exception>
    internal static SaBuildParameters Resolve(long n, SaBuildOptions options)
    {
        ArgumentNullException.ThrowIfNull(options);
        ArgumentOutOfRangeException.ThrowIfNegative(n);
        ArgumentOutOfRangeException.ThrowIfNegative(options.ScaleN);
        ArgumentOutOfRangeException.ThrowIfNegative(options.RuleN);
        ArgumentOutOfRangeException.ThrowIfLessThan(options.Threads, 1);
        SaBuildTuning t = options.Tuning;
        bool scaled = options.ScaleN > 0;
        long ruleN = options.RuleN > 0 ? options.RuleN : n;
        double f = scaled ? ruleN / (double)options.ScaleN : 1.0;
        long Scale(long reference, long floor) => scaled ? Math.Max(floor, (long)Math.Round(reference * f)) : reference;
        var p = new SaBuildParameters
        {
            N = n,
            Factor = f,
            Scaled = scaled,
            DcsBucketCapacity = t.DcsBucketCapacity ?? Math.Max(3, Scale(DcsBucketReference, 3)),
            DcsBaseCase = t.DcsBaseCase ?? Math.Max(4, Scale(DcsBaseReference, 4)),
            DcsChunkSlots = t.DcsChunkSlots ?? Math.Max(16, Scale(DcsChunkReference, 16)),
            DcsSortCapacity = t.DcsSortCapacity ?? Math.Max(3, Scale(DcsSortReference, 3)),
            DcsRunSpillBytes = t.DcsRunSpillBytes ?? (int)Scale(RunSpillReference, 0),
            DcsFanOut = t.DcsFanOut ?? 1024,
            DcsReservoir = t.DcsReservoir ?? 16384,
            DcsIoBytes = t.DcsIoBytes ?? 64 * 1024,
            DcsMergeBudget = t.DcsMergeBudget ?? 16L << 20,
            BucketCapacity = t.BucketCapacity ?? Math.Max(2, scaled ? (long)((decimal)ruleN * BucketReference / options.ScaleN) : BucketUnscaled),
            FanOut = t.FanOut ?? 64,
            Oversample = t.Oversample ?? 16,
            SpillBytes = t.SpillBytes ?? 64 * 1024,
            WindowBytes = t.WindowBytes ?? MappedWindowWriter.DefaultWindowBytes,
        };
        ArgumentOutOfRangeException.ThrowIfLessThan(p.DcsBucketCapacity, 2);
        ArgumentOutOfRangeException.ThrowIfGreaterThan(p.DcsBucketCapacity, int.MaxValue / 2);
        ArgumentOutOfRangeException.ThrowIfLessThan(p.DcsBaseCase, 1);
        ArgumentOutOfRangeException.ThrowIfLessThan(p.DcsChunkSlots, 1);
        ArgumentOutOfRangeException.ThrowIfLessThan(p.DcsSortCapacity, 2);
        ArgumentOutOfRangeException.ThrowIfNegative(p.DcsRunSpillBytes);
        ArgumentOutOfRangeException.ThrowIfLessThan(p.DcsFanOut, 2);
        ArgumentOutOfRangeException.ThrowIfLessThan(p.DcsReservoir, 3);
        ArgumentOutOfRangeException.ThrowIfLessThan(p.DcsIoBytes, 16);
        ArgumentOutOfRangeException.ThrowIfLessThan(p.DcsMergeBudget, 8);
        ArgumentOutOfRangeException.ThrowIfLessThan(p.BucketCapacity, 2);
        ArgumentOutOfRangeException.ThrowIfGreaterThan(p.BucketCapacity, int.MaxValue / 2);
        ArgumentOutOfRangeException.ThrowIfLessThan(p.FanOut, 2);
        ArgumentOutOfRangeException.ThrowIfLessThan(p.Oversample, 1);
        ArgumentOutOfRangeException.ThrowIfLessThan(p.SpillBytes, 1);
        ArgumentOutOfRangeException.ThrowIfLessThan(p.WindowBytes, 1);
        return p;
    }
}

/// <summary>What a finished build knows.</summary>
internal sealed class SaBuildResult
{
    internal required SaBuildParameters Parameters { get; init; }

    internal long N { get; init; }

    /// <summary>The period of the cover that was used, its number of residues and how it was constructed.</summary>
    internal int Period { get; init; }

    internal int CoverResidues { get; init; }

    internal CoverKind CoverKind { get; init; }

    /// <summary>PH1 level by level.</summary>
    internal required ExtDcsStats Dcs { get; init; }

    /// <summary>PH2 / PH3.</summary>
    internal required SortStats Sort { get; init; }

    /// <summary>The counters of the build at its end (the shared counters, whatever they held before included).</summary>
    internal IoSnapshot Io { get; init; }

    /// <summary>The I/O of PH1 and of PH2 / PH3 each (differences of the counters at the phase boundaries).</summary>
    internal IoSnapshot Ph1Io { get; init; }

    internal IoSnapshot Ph23Io { get; init; }

    /// <summary>Wall seconds of PH1 (up to and including the optional check of the rank file) and of PH2 / PH3 (without the optional fsync).</summary>
    internal double Ph1Seconds { get; init; }

    internal double Ph23Seconds { get; init; }

    internal double Seconds => Ph1Seconds + Ph23Seconds;

    /// <summary>Peak of the native buffers of PH1 and of PH2 / PH3.</summary>
    internal long Ph1ArenaPeakBytes { get; init; }

    internal long Ph23ArenaPeakBytes { get; init; }

    /// <summary>Bytes of the rank file (scratch while PH3 reads it).</summary>
    internal long RankFileBytes { get; init; }

    /// <summary>Peak scratch bytes at the end of PH1 (the rank file not yet counted, as in the prototype) and of the whole build (the rank file and the position files of PH3); peaks of the counters, so with a counter that served other stages they include those.</summary>
    internal long Ph1ScratchPeakBytes { get; init; }

    internal long Ph23ScratchPeakBytes { get; init; }

    /// <summary>Peak RssAnon of the whole process in MB (10^6 bytes) while the build ran, in PH1 and in PH2 / PH3, or null where it is unavailable.</summary>
    internal double? RssAnonPeakMb { get; init; }

    internal double? Ph1RssAnonPeakMb { get; init; }

    internal double? Ph23RssAnonPeakMb { get; init; }

    internal RssVerdict RssVerdict { get; init; }

    /// <summary>What the creation of the suffix array file did about the disk space.</summary>
    internal ReservationOutcome Reservation { get; init; }

    /// <summary>Levels of the recursion below the text (0 = the names were unique or the first reduced string fitted the base case).</summary>
    internal int RecursionLevels => Dcs.MaxLevel;
}

/// <summary>
/// PH1 (the difference-cover rank file), PH2 (the partition) and PH3 (the bucket sort) of ADR-004 on a packed text, in one process, as E6 ran them: the suffix array (<c>n + 1</c> little-endian <c>uint</c>, <c>SA[0] = n</c>) is
/// created at its final length through <see cref="MappedOutputFile.Create"/> (a full disk is an <see cref="IOException"/> at creation, found before PH1 starts) and written through bounded mapped windows at the final offsets of
/// the buckets; nothing is built in RAM and written afterwards. The rank file is scratch: PH1 creates it, PH3 reads it through a read-only mapping, and it is deleted as soon as PH3 is done. Every failure ends the build with its own
/// exception type, removes the suffix array, the rank file and the scratch subdirectories (the scratch root only when the build created it and it is empty) and leaves no descriptor or mapping behind. The RssAnon gate is a
/// verdict in the result, never an exception.
/// </summary>
internal static class SaBuilder
{
    /// <summary>Builds the suffix array of <paramref name="text"/> into the file <paramref name="saPath"/>.</summary>
    /// <exception cref="NotSupportedException">The text is longer than <see cref="PackedTextFormat.MaxSymbols"/>.</exception>
    /// <exception cref="IOException">A file cannot be created, written or spilled to (a full disk included); the suffix array and the rank file are removed.</exception>
    /// <exception cref="InvalidDataException">The rank file or a level of it is inconsistent.</exception>
    /// <exception cref="ArenaBudgetExceededException">The native buffers of a phase would exceed <see cref="SaBuildOptions.ArenaCapBytes"/>.</exception>
    internal static SaBuildResult Build(MappedTextReader text, string saPath, SaBuildOptions options)
    {
        ArgumentNullException.ThrowIfNull(text);
        ArgumentNullException.ThrowIfNull(saPath);
        ArgumentNullException.ThrowIfNull(options);
        return BuildCore(text, saPath, null, 0, options);
    }

    private static SaBuildResult BuildCore(MappedTextReader text, string? saPath, MappedOutputFile? section, long saOffset, SaBuildOptions options)
    {
        ArgumentNullException.ThrowIfNull(options.ScratchDir);
        ArgumentOutOfRangeException.ThrowIfNegative(options.ArenaCapBytes);
        ArgumentOutOfRangeException.ThrowIfNegative(options.RssAnonCapMb);
        long n = text.Length;
        if (n > PackedTextFormat.MaxSymbols) throw new NotSupportedException($"{n} symbols: positions are 32-bit, the limit is {PackedTextFormat.MaxSymbols}");
        if (text.Bits > PackedTextFormat.MaxBits) throw new ArgumentException($"a text of {text.Bits} bits per symbol is not a packed text of the file codec", nameof(text));
        if (section is not null)
        {
            ArgumentOutOfRangeException.ThrowIfNegative(saOffset);
            if (saOffset > section.Length - (4 * (n + 1))) throw new ArgumentOutOfRangeException(nameof(saOffset), saOffset, $"the suffix array of {4 * (n + 1)} bytes does not fit the file of {section.Length} bytes at this offset");
        }
        SaBuildParameters p = SaBuildRules.Resolve(n, options);
        DifferenceCover cover = DifferenceCover.Create(options.ExactPeriod ? options.Period : DifferenceCover.NearestSingerPeriod(options.Period));

        string scratch = Path.GetFullPath(options.ScratchDir);
        bool createdRoot = !Directory.Exists(scratch);
        string ph1 = Path.Combine(scratch, "ph1"), ph3 = Path.Combine(scratch, "ph3");
        string rankPath = options.RankPath ?? Path.Combine(scratch, "ranks.bin");
        IoCounters counters = options.Counters ?? new IoCounters();
        var ph1Budget = new ArenaBudget(options.ArenaCapBytes);
        var ph3Budget = new ArenaBudget(options.ArenaCapBytes);
        var ledger = new ScratchLedger(counters);
        using var whole = options.MeasureRss ? new RssSampler(new RssReader(options.ProcSource), 10) : null;
        bool saCreated = false;
        bool rankTouched = false;
        try
        {
            Directory.CreateDirectory(ph1);
            Directory.CreateDirectory(ph3);
            whole?.Start();
            using MappedOutputFile? owned = section is null ? MappedOutputFile.Create(saPath!, 4 * (n + 1), options.Files) : null;
            saCreated = owned is not null;
            MappedOutputFile saFile = section ?? owned!;

            // ---- PH1
            options.Trace?.Invoke("ph1");
            IoSnapshot io0 = counters.Snapshot();
            using var rss1 = options.MeasureRss ? new RssSampler(new RssReader(options.ProcSource), 10) : null;
            rss1?.Start();
            long t0 = Stopwatch.GetTimestamp();
            rankTouched = true;
            ExtDcsStats dcs = ExternalDcs.BuildRankFile(text, cover, rankPath, new ExtDcsOptions
            {
                ScratchDir = ph1,
                BucketCapacity = p.DcsBucketCapacity,
                BaseCase = p.DcsBaseCase,
                ChunkSlots = p.DcsChunkSlots,
                SortLeafCapacity = p.DcsSortCapacity,
                Threads = options.Threads,
                Seed = options.Seed,
                FanOut = p.DcsFanOut,
                ReservoirSize = p.DcsReservoir,
                RunSpillBytes = p.DcsRunSpillBytes,
                IoBufferBytes = p.DcsIoBytes,
                MergeBufferBudget = p.DcsMergeBudget,
                WindowBytes = p.WindowBytes,
                Budget = ph1Budget,
                ArenaGuards = options.ArenaGuards,
                Counters = counters,
                Files = options.Files,
                Trace = options.Trace,
            });
            long rankBytes = new FileInfo(rankPath).Length;
            long ph1ScratchPeak = counters.ScratchPeakBytes;
            ledger.Add(rankBytes); // the rank file is scratch of the pipeline from here on
            if (options.VerifyRanks)
            {
                options.Trace?.Invoke("verify");
                using RankTable table = RankTable.Open(rankPath, cover, n);
                table.Verify();
            }
            double ph1Seconds = Seconds(t0);
            RssStats rssPh1 = rss1?.Stop() ?? default;
            IoSnapshot io1 = counters.Snapshot();

            // ---- PH2 / PH3
            options.Trace?.Invoke("ph23");
            using var rss3 = options.MeasureRss ? new RssSampler(new RssReader(options.ProcSource), 10) : null;
            rss3?.Start();
            t0 = Stopwatch.GetTimestamp();
            SortStats sort;
            using (RankTable ranks = RankTable.Open(rankPath, cover, n))
            {
                sort = BucketSorter.Run(text, new SortOptions
                {
                    LeafCapacity = p.BucketCapacity,
                    FanOut = p.FanOut,
                    Oversample = p.Oversample,
                    Threads = options.Threads,
                    Seed = options.Seed,
                    ScratchDir = ph3,
                    SpillBytes = p.SpillBytes,
                    Ledger = ledger,
                    Budget = ph3Budget,
                    ArenaGuards = options.ArenaGuards,
                    Counters = counters,
                    Trace = options.Trace,
                }, section is null ? new MappedSaSliceSink(saFile, p.WindowBytes, counters) : new MappedSaSectionSink(saFile, saOffset, 4 * (n + 1), p.WindowBytes, counters), ranks);
            }
            File.Delete(rankPath); // the rank file has served its only consumer
            ledger.Remove(rankBytes);
            double ph23Seconds = Seconds(t0);
            RssStats rssPh3 = rss3?.Stop() ?? default;
            IoSnapshot io2 = counters.Snapshot();
            if (section is null)
            {
                saFile.Complete(options.SyncOutput);
                CheckLength(saPath!, n);
            }
            RssStats all = whole?.Stop() ?? default;
            double? anonMb = all.Available ? all.MaxRssAnonKb * 1024.0 / 1e6 : null;
            RssVerdict verdict = RssVerdict.Unavailable;
            if (anonMb is { } mb) verdict = mb <= options.RssAnonCapMb ? RssVerdict.Within : RssVerdict.Exceeded;
            return new SaBuildResult
            {
                Parameters = p,
                N = n,
                Period = cover.Period,
                CoverResidues = cover.Residues.Length,
                CoverKind = cover.Kind,
                Dcs = dcs,
                Sort = sort,
                Io = io2,
                Ph1Io = io1.Since(io0),
                Ph23Io = io2.Since(io1),
                Ph1Seconds = ph1Seconds,
                Ph23Seconds = ph23Seconds,
                Ph1ArenaPeakBytes = ph1Budget.PeakBytes,
                Ph23ArenaPeakBytes = ph3Budget.PeakBytes,
                RankFileBytes = rankBytes,
                Ph1ScratchPeakBytes = ph1ScratchPeak,
                Ph23ScratchPeakBytes = io2.ScratchPeakBytes,
                RssAnonPeakMb = anonMb,
                Ph1RssAnonPeakMb = rssPh1.Available ? rssPh1.MaxRssAnonKb * 1024.0 / 1e6 : null,
                Ph23RssAnonPeakMb = rssPh3.Available ? rssPh3.MaxRssAnonKb * 1024.0 / 1e6 : null,
                RssVerdict = verdict,
                Reservation = saFile.Outcome,
            };
        }
        catch
        {
            ledger.ReleaseAll();
            if (saCreated) DeleteQuietly(saPath!);
            if (rankTouched) DeleteQuietly(rankPath);
            throw;
        }
        finally
        {
            DeleteDirectoryQuietly(ph1);
            DeleteDirectoryQuietly(ph3);
            if (createdRoot) DeleteEmptyDirectoryQuietly(scratch);
        }
    }

    /// <summary>
    /// Builds the suffix array of <paramref name="text"/> into the range <c>[saOffset, saOffset + 4 (n + 1))</c> of <paramref name="file"/> (a section of the container): the file exists at that length already, the stage
    /// writes through windows of it and neither completes nor deletes it (the owner of the file does, also on a failure); everything else is as in <see cref="Build(MappedTextReader, string, SaBuildOptions)"/>.
    /// </summary>
    /// <exception cref="ArgumentOutOfRangeException">The range does not lie inside the file.</exception>
    internal static SaBuildResult Build(MappedTextReader text, MappedOutputFile file, long saOffset, SaBuildOptions options)
    {
        ArgumentNullException.ThrowIfNull(text);
        ArgumentNullException.ThrowIfNull(file);
        ArgumentNullException.ThrowIfNull(options);
        return BuildCore(text, null, file, saOffset, options);
    }

    /// <summary>The suffix array has its final length: <c>4 (n + 1)</c> bytes.</summary>
    internal static void CheckLength(string saPath, long n)
    {
        long got = new FileInfo(saPath).Length;
        if (got != 4 * (n + 1)) throw new InvalidOperationException($"the suffix array has {got} bytes, expected {4 * (n + 1)}");
    }

    private static double Seconds(long t0) => (Stopwatch.GetTimestamp() - t0) / (double)Stopwatch.Frequency;

    private static void DeleteQuietly(string path)
    {
        try
        {
            File.Delete(path);
        }
        catch (Exception e) when (e is IOException or UnauthorizedAccessException)
        {
            // not ours to clean up if it cannot be deleted; the failure of the build is what the caller needs
        }
    }

    private static void DeleteDirectoryQuietly(string path)
    {
        try
        {
            if (Directory.Exists(path)) Directory.Delete(path, true);
        }
        catch (Exception e) when (e is IOException or UnauthorizedAccessException)
        {
            // reported by the files that are left
        }
    }

    private static void DeleteEmptyDirectoryQuietly(string path)
    {
        try
        {
            if (Directory.Exists(path) && !Directory.EnumerateFileSystemEntries(path).Any()) Directory.Delete(path);
        }
        catch (Exception e) when (e is IOException or UnauthorizedAccessException)
        {
            // someone else put a file there between the check and the delete: it is theirs
        }
    }
}
