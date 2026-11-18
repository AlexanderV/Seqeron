using System.Diagnostics;
using System.Globalization;
using System.Numerics;

namespace SuffixTree.Esa.External;

/// <summary>The files PH4 and PH5 write: the base array of the LCP, its overflow list, the child table (layout A) and, when a path is given, the statistics text.</summary>
/// <param name="Lcp">Base array: <c>n + 1</c> entries of the width the rule of ADR-002 section 4 chooses (little endian).</param>
/// <param name="Overflow">Overflow list: pairs <c>(rank, value)</c> of two 4-byte fields (empty for U32 and when no value escapes).</param>
/// <param name="Child">Child table: <c>n + 1</c> entries of 4 bytes.</param>
/// <param name="Stats">The statistics (<see cref="LcpChildResult.StatsText"/>), or null for none.</param>
internal sealed record LcpChildOutputs(string Lcp, string Overflow, string Child, string? Stats = null);

/// <summary>
/// Explicit values for the sizes that the rule would otherwise derive (null = derive). The tests set every one of them to tiny values (R = 3..50, blocks of 2..5 entries, windows of a few bytes) so that every
/// multi-chunk, spill and window path runs on texts of 10^3-10^5 symbols; the output never depends on any of them.
/// </summary>
internal sealed record LcpChildTuning
{
    /// <summary>PH4: R, positions per Phi chunk and ranks per LCP chunk.</summary>
    internal long? LcpChunk { get; init; }

    /// <summary>PH4: parallel B2 tasks (null = what the arena cap allows at the reference geometry).</summary>
    internal int? LcpWorkers { get; init; }

    /// <summary>PH4: sequential read buffer.</summary>
    internal int? LcpReadBytes { get; init; }

    /// <summary>PH4: write buffer of the tuple distributions.</summary>
    internal int? LcpTupleBytes { get; init; }

    /// <summary>PH4: sliding window of the LCP and overflow writers.</summary>
    internal int? WindowBytes { get; init; }

    /// <summary>PH5: entries per block of the stack.</summary>
    internal int? BlockEntries { get; init; }

    /// <summary>PH5: window of the CHILD writer (a power of two).</summary>
    internal int? ChildWindowBytes { get; init; }

    /// <summary>PH5: windows of the CHILD writer.</summary>
    internal int? ChildSlots { get; init; }

    /// <summary>PH5: fsync of the spill file after this many bytes.</summary>
    internal long? SpillFlushBytes { get; init; }
}

/// <summary>What <see cref="LcpChildBuilder.Build"/> needs besides the files.</summary>
internal sealed record LcpChildOptions
{
    /// <summary>The scratch directory (created when missing; <c>ph4</c> and <c>ph5</c> inside it are created and removed by the build, the directory itself is removed when this build created it and it is empty).</summary>
    internal required string ScratchDir { get; init; }

    /// <summary>Worker threads of B2 (the B1, B3 and PH5 work is sequential).</summary>
    internal int Threads { get; init; } = 4;

    /// <summary>The codec of the child table: the owner's decision of 2026-11-13 is <see cref="ChildCodec.Packed"/>; <see cref="ChildCodec.Tagged"/> holds ranks in 30 bits (texts up to 2^30 - 1 symbols).</summary>
    internal ChildCodec Codec { get; init; } = ChildCodec.Packed;

    /// <summary>
    /// The scale-down emulation of ADR-004 section 5 (E3 to E7): when set (3,000,000,000 for the 100 Mbp runs), R is its reference value for 3x10^9 symbols multiplied by n / ScaleN; 0 = the reference values themselves.
    /// </summary>
    internal long ScaleN { get; init; }

    /// <summary>E7 series A: the n at which the scale-down rule is evaluated (0 = the text length). Ignored when <see cref="ScaleN"/> is 0.</summary>
    internal long RuleN { get; init; }

    /// <summary>Cap of the native buffers of PH4 (owner decision: 160 MiB).</summary>
    internal long ArenaCapBytes { get; init; } = 160L << 20;

    /// <summary>Cap of the native buffers of PH5 (it lives next to PH4 in the fused phase; owner decision: 64 MiB).</summary>
    internal long ChildArenaCapBytes { get; init; } = 64L << 20;

    internal LcpChildTuning Tuning { get; init; } = new();

    /// <summary>The counters of the build (a private instance when null; the result carries a snapshot either way).</summary>
    internal IoCounters? Counters { get; init; }

    /// <summary>How the output files are created (reservation of the disk space); null = <see cref="ReservationMode.Auto"/>.</summary>
    internal MappedFileOptions? Files { get; init; }

    /// <summary>Test hook: options of the arenas of PH4 and PH5 (guard bands, poison pattern); null = production arenas.</summary>
    internal ArenaOptions? ArenaGuards { get; init; }

    /// <summary>The alphabet size, which the statistics text records; required when <see cref="LcpChildOutputs.Stats"/> is set.</summary>
    internal int? Sigma { get; init; }

    /// <summary>fsync every output file when the build is done (not part of the stage time of the prototype).</summary>
    internal bool SyncOutputs { get; init; }

    /// <summary>Sample RssAnon while the build runs (the figure of the gate); unavailable where there is no <c>/proc</c>.</summary>
    internal bool MeasureRss { get; init; } = true;

    /// <summary>Where the kernel's text files come from (null = the real <c>/proc</c>); a test passes a source that fails like a platform without one.</summary>
    internal IProcSource? ProcSource { get; init; }

    /// <summary>The gate of the owner decision of 2026-11-18: RssAnon of the process at most 320 MB (10^6 bytes) at every n with the default parameters and four threads.</summary>
    internal double RssAnonCapMb { get; init; } = 320;

    /// <summary>Test hook, called with <c>"b1"</c>, <c>"b2"</c>, <c>"b3"</c> when a step of PH4 begins; an exception thrown here fails the build like a failing step.</summary>
    internal Action<string>? Trace { get; init; }
}

/// <summary>Every size of one run, resolved from <see cref="LcpChildOptions"/> and the text length (what the report prints, and what the tests pin).</summary>
internal sealed record LcpChildParameters
{
    internal long N { get; init; }

    /// <summary>n / ScaleN (or RuleN / ScaleN), or 1 when unscaled.</summary>
    internal double Factor { get; init; }

    internal bool Scaled { get; init; }

    // PH4
    internal long LcpChunk { get; init; }

    internal int LcpWorkers { get; init; }

    internal int LcpReadBytes { get; init; }

    internal int LcpTupleBytes { get; init; }

    internal int WindowBytes { get; init; }

    // PH5
    internal int BlockEntries { get; init; }

    internal int ChildWindowBytes { get; init; }

    internal int ChildSlots { get; init; }

    internal long SpillFlushBytes { get; init; }
}

/// <summary>The scale-down rule of ADR-004 (E4, E5, E6) for PH4 and PH5 only.</summary>
internal static class LcpChildRules
{
    /// <summary>The reference text length of the sizes: 3x10^9 symbols.</summary>
    internal const long ReferenceN = 3_000_000_000;

    /// <summary>The reference R: 8,000,000 positions per chunk (P = 375 chunks at <see cref="ReferenceN"/>).</summary>
    internal const long ChunkReference = 8_000_000;

    /// <summary>
    /// The sizes of a run: with <c>f = RuleN / ScaleN</c> (RuleN = the text length when it is 0) R is <c>ceil(8,000,000 x f)</c> with a floor of 3 when <see cref="LcpChildOptions.ScaleN"/> is set, else 8,000,000; the number of parallel B2
    /// tasks is what the PH4 arena cap allows at the reference geometry (R = 8,000,000 at 3x10^9 symbols: about 77 MB per task, two of four under 160 MiB); the buffers and windows have their reference values. An explicit value of
    /// <see cref="LcpChildOptions.Tuning"/> wins.
    /// </summary>
    /// <exception cref="ArgumentOutOfRangeException">A negative n, scale or cap, fewer than one thread, or an explicit value that its stage cannot use.</exception>
    internal static LcpChildParameters Resolve(long n, LcpChildOptions options)
    {
        ArgumentNullException.ThrowIfNull(options);
        ArgumentOutOfRangeException.ThrowIfNegative(n);
        ArgumentOutOfRangeException.ThrowIfNegative(options.ScaleN);
        ArgumentOutOfRangeException.ThrowIfNegative(options.RuleN);
        ArgumentOutOfRangeException.ThrowIfLessThan(options.Threads, 1);
        LcpChildTuning t = options.Tuning;
        bool scaled = options.ScaleN > 0;
        long ruleN = options.RuleN > 0 ? options.RuleN : n;
        double f = scaled ? ruleN / (double)options.ScaleN : 1.0;
        long chunk = t.LcpChunk ?? (scaled ? Math.Max(3, (long)Math.Ceiling(ChunkReference * f)) : ChunkReference);
        var reference = new ExtLcpOptions { ScratchDir = string.Empty, ChunkPositions = ChunkReference, Threads = options.Threads };
        long perTask = ExternalLcp.WorkerBytes(ReferenceN, reference);
        int workers = t.LcpWorkers ?? (int)Math.Clamp(options.ArenaCapBytes / Math.Max(1, perTask), 1, options.Threads);
        var p = new LcpChildParameters
        {
            N = n,
            Factor = f,
            Scaled = scaled,
            LcpChunk = chunk,
            LcpWorkers = workers,
            LcpReadBytes = t.LcpReadBytes ?? 1 << 20,
            LcpTupleBytes = t.LcpTupleBytes ?? 32 * 1024,
            WindowBytes = t.WindowBytes ?? MappedWindowWriter.DefaultWindowBytes,
            BlockEntries = t.BlockEntries ?? 256 * 1024,
            ChildWindowBytes = t.ChildWindowBytes ?? 1 << 20,
            ChildSlots = t.ChildSlots ?? 8,
            SpillFlushBytes = t.SpillFlushBytes ?? 16L << 20,
        };
        ArgumentOutOfRangeException.ThrowIfLessThan(p.LcpChunk, 1);
        ArgumentOutOfRangeException.ThrowIfLessThan(p.LcpWorkers, 1);
        ArgumentOutOfRangeException.ThrowIfLessThan(p.LcpReadBytes, 1);
        ArgumentOutOfRangeException.ThrowIfLessThan(p.LcpTupleBytes, 1);
        ArgumentOutOfRangeException.ThrowIfLessThan(p.WindowBytes, 1);
        ArgumentOutOfRangeException.ThrowIfLessThan(p.BlockEntries, 2);
        if (p.ChildWindowBytes < 4 || !BitOperations.IsPow2(p.ChildWindowBytes)) throw new ArgumentOutOfRangeException(nameof(options), p.ChildWindowBytes, "a window of the child table is a power of two of at least 4 bytes");
        ArgumentOutOfRangeException.ThrowIfLessThan(p.ChildSlots, 1);
        ArgumentOutOfRangeException.ThrowIfNegative(p.SpillFlushBytes);
        return p;
    }
}

/// <summary>The verdict of the RssAnon gate of the builder.</summary>
internal enum RssVerdict
{
    /// <summary>The platform has no <c>/proc</c> (or sampling was off): nothing was measured, the gate was neither passed nor failed.</summary>
    Unavailable,

    /// <summary>The peak RssAnon of the process was at most the cap.</summary>
    Within,

    /// <summary>The peak RssAnon of the process exceeded the cap.</summary>
    Exceeded,
}

/// <summary>What a finished build knows.</summary>
internal sealed class LcpChildResult
{
    internal required LcpChildParameters Parameters { get; init; }

    internal required ExtLcpStats Lcp { get; init; }

    internal required ChildScanResult Scan { get; init; }

    /// <summary>The counters of the build at its end.</summary>
    internal IoSnapshot Io { get; init; }

    /// <summary>Wall seconds of PH4 + PH5 (without the optional fsync).</summary>
    internal double Seconds { get; init; }

    /// <summary>Peak of the native buffers of PH4 and of PH5.</summary>
    internal long ArenaPeakBytes { get; init; }

    internal long ChildArenaPeakBytes { get; init; }

    /// <summary>Peak RssAnon of the whole process in MB (10^6 bytes) while the build ran, or null where it is unavailable.</summary>
    internal double? RssAnonPeakMb { get; init; }

    internal RssVerdict RssVerdict { get; init; }

    /// <summary>What the creation of the LCP base file did about the disk space (the other outputs use the same mode).</summary>
    internal ReservationOutcome Reservation { get; init; }

    /// <summary>The line-oriented statistics text (<c>.stats</c> of the prototype): NODE_COUNT, the internal nodes, LRS_DEPTH and LRS_RANK and the shape of the output.</summary>
    internal string StatsText(int sigma) => string.Create(CultureInfo.InvariantCulture,
        $"nodecount={Scan.NodeCount}\ninternal={Scan.InternalCount}\nlrsdepth={Scan.LrsDepth}\nlrsrank={Scan.LrsRank}\nn={Lcp.N}\nsigma={sigma}\nlcpwidth={LcpWidth.Bytes(Lcp.Kind)}\noverflow={Lcp.OverflowCount}\nmaxlcp={Scan.MaxLcp}\n");
}

/// <summary>
/// PH4 (base mode) and PH5 (layout A) of ADR-004 on a given suffix-array file, in one process, fused through the B3 hook exactly as in E6: the compact LCP, its overflow list and the child table are written through
/// bounded mapped windows (every big output is created with <see cref="MappedOutputFile.Create"/>, so a full disk is an <see cref="IOException"/> at creation), nothing is built in RAM and written afterwards.
/// </summary>
internal static class LcpChildBuilder
{
    /// <summary>Builds the outputs of <paramref name="text"/> from the suffix array file <paramref name="saPath"/> (<c>n + 1</c> little-endian <c>uint</c> entries, <c>SA[0] = n</c>).</summary>
    /// <exception cref="InvalidDataException">The suffix array file is not the file of a permutation of the text positions (wrong length, <c>SA[0] != n</c>, an entry out of range or repeated).</exception>
    /// <exception cref="NotSupportedException">The text is longer than <see cref="PackedTextFormat.MaxSymbols"/>, or longer than the tagged codec holds.</exception>
    /// <exception cref="IOException">A file cannot be created, written or spilled to (a full disk included); every output of this build is removed.</exception>
    internal static LcpChildResult Build(string saPath, MappedTextReader text, LcpChildOutputs outputs, LcpChildOptions options)
    {
        ArgumentNullException.ThrowIfNull(saPath);
        ArgumentNullException.ThrowIfNull(text);
        ArgumentNullException.ThrowIfNull(outputs);
        ArgumentNullException.ThrowIfNull(options);
        ArgumentNullException.ThrowIfNull(outputs.Lcp);
        ArgumentNullException.ThrowIfNull(outputs.Overflow);
        ArgumentNullException.ThrowIfNull(outputs.Child);
        ArgumentOutOfRangeException.ThrowIfNegative(options.ArenaCapBytes);
        ArgumentOutOfRangeException.ThrowIfNegative(options.ChildArenaCapBytes);
        if (options.Sigma is < 0) throw new ArgumentOutOfRangeException(nameof(options), options.Sigma, "sigma is not negative");
        if (outputs.Stats is not null && options.Sigma is null) throw new ArgumentException("the statistics file records the alphabet size: set LcpChildOptions.Sigma", nameof(options));
        long n = text.Length;
        if (n > PackedTextFormat.MaxSymbols) throw new NotSupportedException($"{n} symbols: positions are 32-bit, the limit is {PackedTextFormat.MaxSymbols}");
        if (options.Codec == ChildCodec.Tagged && n > TaggedChild.MaxLength) throw new NotSupportedException($"the tagged codec holds ranks in 30 bits; a text of {n} symbols needs the packed codec");
        LcpChildParameters p = LcpChildRules.Resolve(n, options);

        string scratch = Path.GetFullPath(options.ScratchDir);
        bool createdRoot = !Directory.Exists(scratch);
        string ph4 = Path.Combine(scratch, "ph4"), ph5 = Path.Combine(scratch, "ph5");
        IoCounters counters = options.Counters ?? new IoCounters();
        var lcpBudget = new ArenaBudget(options.ArenaCapBytes);
        var childBudget = new ArenaBudget(options.ChildArenaCapBytes);
        using var sampler = options.MeasureRss ? new RssSampler(new RssReader(options.ProcSource), 10) : null;
        var childPath = outputs.Child;
        bool childCreated = false;
        bool lcpCreated = false;
        bool statsStarted = false;
        try
        {
            Directory.CreateDirectory(ph4);
            Directory.CreateDirectory(ph5);
            sampler?.Start();
            var lcpOptions = new ExtLcpOptions
            {
                ScratchDir = ph4,
                ChunkPositions = p.LcpChunk,
                Threads = options.Threads,
                Workers = p.LcpWorkers,
                ReadBufferBytes = p.LcpReadBytes,
                TupleBufferBytes = p.LcpTupleBytes,
                WindowBytes = p.WindowBytes,
                Budget = lcpBudget,
                ArenaGuards = options.ArenaGuards,
                Counters = counters,
                Files = options.Files,
                Trace = options.Trace,
            };
            var childOptions = new ChildBuildOptions
            {
                ScratchDir = ph5,
                BlockEntries = p.BlockEntries,
                WindowBytes = p.ChildWindowBytes,
                WindowSlots = p.ChildSlots,
                FlushWindows = true,
                SpillFlushBytes = p.SpillFlushBytes,
                Codec = options.Codec,
                Budget = childBudget,
                ArenaGuards = options.ArenaGuards,
                Counters = counters,
            };
            ExtLcpStats lcp;
            ChildScanResult scan;
            long t0 = Stopwatch.GetTimestamp();
            using (var lazy = new LazyChildScan(n, childPath, childOptions, options.Files))
            {
                try
                {
                    lcp = ExternalLcp.Build(saPath, text, outputs.Lcp, outputs.Overflow, lcpOptions, lazy);
                    scan = lazy.Complete();
                }
                finally
                {
                    childCreated = lazy.Created;
                }
            }
            lcpCreated = true; // PH4 removes its own files when it fails; from here on they are this build's to remove
            double seconds = (Stopwatch.GetTimestamp() - t0) / (double)Stopwatch.Frequency;
            if (options.SyncOutputs)
            {
                foreach (string path in new[] { outputs.Lcp, outputs.Overflow, childPath }) Sync(path);
            }
            CheckLengths(outputs, n, lcp, childPath);
            RssStats rss = sampler?.Stop() ?? default;
            double? anonMb = rss.Available ? rss.MaxRssAnonKb * 1024.0 / 1e6 : null;
            RssVerdict verdict = RssVerdict.Unavailable;
            if (anonMb is { } mb) verdict = mb <= options.RssAnonCapMb ? RssVerdict.Within : RssVerdict.Exceeded;
            var result = new LcpChildResult
            {
                Parameters = p,
                Lcp = lcp,
                Scan = scan,
                Io = counters.Snapshot(),
                Seconds = seconds,
                ArenaPeakBytes = lcpBudget.PeakBytes,
                ChildArenaPeakBytes = childBudget.PeakBytes,
                RssAnonPeakMb = anonMb,
                RssVerdict = verdict,
                Reservation = lcp.Reservation,
            };
            if (outputs.Stats is not null)
            {
                statsStarted = true;
                File.WriteAllText(outputs.Stats, result.StatsText(options.Sigma!.Value));
            }
            return result;
        }
        catch
        {
            if (lcpCreated)
            {
                DeleteQuietly(outputs.Lcp);
                DeleteQuietly(outputs.Overflow);
            }
            if (childCreated) DeleteQuietly(childPath);
            if (statsStarted && !Directory.Exists(outputs.Stats)) DeleteQuietly(outputs.Stats!);
            throw;
        }
        finally
        {
            DeleteDirectoryQuietly(ph4);
            DeleteDirectoryQuietly(ph5);
            if (createdRoot) DeleteEmptyDirectoryQuietly(scratch);
        }
    }

    /// <summary>The child table of PH5 is created when PH4 delivers its first block, i.e. after the width of the LCP is known and the LCP files exist (the creation order of E6, which keeps the peak of the scratch on the disk where the ADR measured it).</summary>
    private sealed class LazyChildScan(long n, string childPath, ChildBuildOptions options, MappedFileOptions? files) : ILcpChunkConsumer, IDisposable
    {
        private MappedOutputFile? _file;
        private ChildScan? _scan;

        internal bool Created { get; private set; }

        public void Consume(long firstRank, ReadOnlySpan<uint> lcp)
        {
            if (_scan is null)
            {
                _file = MappedOutputFile.Create(childPath, 4 * (n + 1), files);
                Created = true;
                _scan = ChildScan.Create(n, _file, 0, options);
            }
            _scan.Consume(firstRank, lcp);
        }

        internal ChildScanResult Complete()
        {
            ChildScan scan = _scan ?? throw new InvalidOperationException("PH4 delivered no LCP block");
            ChildScanResult result = scan.Complete();
            _file!.Complete(fsync: false);
            return result;
        }

        public void Dispose()
        {
            _scan?.Dispose();
            _file?.Dispose();
        }
    }

    /// <summary>Every output has its final length: the base LCP of <c>(n + 1)</c> entries of the chosen width, the overflow pairs, the child table of <c>n + 1</c> entries.</summary>
    private static void CheckLengths(LcpChildOutputs outputs, long n, ExtLcpStats lcp, string childPath)
    {
        static void Expect(string path, long bytes, string what)
        {
            long got = new FileInfo(path).Length;
            if (got != bytes) throw new InvalidOperationException($"{what} has {got} bytes, expected {bytes}");
        }
        Expect(outputs.Lcp, LcpWidth.Bytes(lcp.Kind) * (n + 1), "the base LCP");
        Expect(outputs.Overflow, 8 * lcp.OverflowCount, "the overflow list");
        Expect(childPath, 4 * (n + 1), "the child table");
    }

    private static void Sync(string path)
    {
        using var fs = new FileStream(path, FileMode.Open, FileAccess.ReadWrite);
        fs.Flush(true);
    }

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
