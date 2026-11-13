using System.Diagnostics;
using System.Globalization;
using SuffixTree.EmProto.Build;
using SuffixTree.EmProto.Dcs;
using SuffixTree.EmProto.Experiments;
using SuffixTree.EmProto.Harness;
using SuffixTree.EmProto.Ph4;
using SuffixTree.EmProto.Ph5;
using SuffixTree.EmProto.Text;

namespace SuffixTree.EmProto.Pipeline;

/// <summary>What one phase of the pipeline cost: wall time, the program's own byte counters, the kernel's byte counters of the process, RssAnon, scratch and dirty pages.</summary>
internal sealed record PhaseRecord(
    string Name,
    double Seconds,
    long ExplicitReadBytes,
    long ExplicitWriteBytes,
    long MappedWriteBytes,
    ProcIoSnapshot Proc,
    double RssAnonPeakMb,
    double RssFilePeakMb,
    long ScratchPeakBytes,
    long ScratchLiveBytes,
    long DirtyPeakKb,
    long ProcessDirtyPeakKb);

/// <summary>The output files of a run.</summary>
internal sealed record EndToEndFiles(string Packed, string PackedMeta, string Sa, string Lcp, string Overflow, string Child, string Stats)
{
    internal static EndToEndFiles Of(string prefix) => new(prefix + ".pk", prefix + ".pk.meta", prefix + ".sa", prefix + ".lcp", prefix + ".ovf", prefix + ".cld", prefix + ".stats");

    internal IEnumerable<string> All => [Packed, PackedMeta, Sa, Lcp, Overflow, Child, Stats];

    /// <summary>The files that are big and written through mappings (everything but the two small text files).</summary>
    internal IEnumerable<string> Big => [Packed, Sa, Lcp, Overflow, Child];
}

/// <summary>What a finished run knows.</summary>
internal sealed class EndToEndResult
{
    internal required EndToEndParameters Parameters { get; init; }

    internal required PrepareResult Prepare { get; init; }

    internal required EndToEndFiles Files { get; init; }

    internal List<PhaseRecord> Phases { get; } = [];

    internal ExtDcsStats Dcs { get; set; } = new();

    internal BuildStats Build { get; set; } = new();

    internal ExtLcpStats Lcp { get; set; } = new();

    internal ChildScanResult Scan { get; set; } = new();

    /// <summary>The options each stage was given (what the tests compare with <see cref="ScaleRules"/>).</summary>
    internal ExtDcsOptions DcsOptions { get; set; } = null!;

    internal BuildOptions BuildOptions { get; set; } = null!;

    internal ExtLcpOptions LcpOptions { get; set; } = null!;

    internal ChildBuildOptions ChildOptions { get; set; } = null!;

    internal int CoverPeriod { get; set; }

    internal int CoverResidues { get; set; }

    /// <summary>Peak RssAnon over the whole run, MB (0 where /proc is not available).</summary>
    internal double RssAnonPeakMb { get; set; }

    /// <summary>Peak of the scratch bytes that the stages track (files of PH1, PH2 / PH3, PH4, PH5, and the rank file as long as it lives).</summary>
    internal long ScratchPeakBytes { get; set; }

    internal long ScratchLiveBytes { get; set; }

    /// <summary>Entries (files or directories) still in the scratch directory after the clean-up; the directory itself is removed.</summary>
    internal int ScratchEntriesLeft { get; set; }

    internal bool ScratchDirectoryExists { get; set; }

    /// <summary>Bytes of the native arenas of the four stages still reserved at the end (0 when every arena was released).</summary>
    internal long ArenaLiveBytes { get; set; }

    internal long ArenaPeakBytes { get; set; }

    internal Dictionary<string, long> FileBytes { get; } = [];

    internal PhaseRecord Phase(string name) => Phases.First(p => p.Name == name);

    /// <summary>Wall seconds of PH0 .. PH5 (without the final sync).</summary>
    internal double StageSeconds => Phases.Where(p => p.Name != EndToEnd.PhaseSync).Sum(p => p.Seconds);

    /// <summary>Wall seconds of PH0 .. PH3 (the suffix array stage of the gate).</summary>
    internal double SaStageSeconds => Phases.Where(p => p.Name is EndToEnd.PhasePrepare or EndToEnd.PhaseDcs or EndToEnd.PhaseSa).Sum(p => p.Seconds);

    /// <summary>The run left nothing behind: no scratch entry, no tracked scratch byte, no live arena byte.</summary>
    internal bool Clean => ScratchEntriesLeft == 0 && !ScratchDirectoryExists && ScratchLiveBytes == 0 && ArenaLiveBytes == 0;
}

/// <summary>
/// E6 of ADR-004 section 5: the whole construction in one process. Text in, then PH0 (alphabet scan, counts, SHA-256, the packed text streamed into its file), PH1 (the difference-cover rank file in external memory, E3b),
/// PH2 / PH3 (partition and bucket sort, the suffix array streamed to its final offsets through mapped windows, E3, the rank file is deleted as soon as it has been used) and PH4 with PH5 as the consumer of its B3 hook
/// (the compact LCP in base mode, E4, and the child table of layout A with the codec of <see cref="EndToEndOptions.Codec"/>, NODE_COUNT, LRS, E5), all with the scale-down rule of <see cref="ScaleRules"/>. Every big output is
/// created at its final length and filled through bounded mapped windows; nothing is built in RAM and written afterwards. The scratch directory is deleted at the end (also when a stage fails).
/// </summary>
internal static unsafe class EndToEnd
{
    internal const string PhasePrepare = "ph0_prepare";
    internal const string PhaseDcs = "ph1_dcs";
    internal const string PhaseSa = "ph23_sa";
    internal const string PhaseLcpChild = "ph45_lcp_child";
    internal const string PhaseSync = "sync";

    /// <summary>The stage phases in the order they run (the sync comes after them, when asked for).</summary>
    internal static readonly string[] StagePhases = [PhasePrepare, PhaseDcs, PhaseSa, PhaseLcpChild];

    internal static EndToEndResult Run(EndToEndOptions o, DirtySampler? dirty = null)
    {
        ArgumentNullException.ThrowIfNull(o);
        ArgumentOutOfRangeException.ThrowIfLessThan(o.Threads, 1);
        if (!File.Exists(o.InputPath)) throw new FileNotFoundException("the input text does not exist", o.InputPath);
        long n = new FileInfo(o.InputPath).Length;
        if (n == 0) throw new ArgumentException("the input text is empty", nameof(o));
        if (n >= uint.MaxValue - 2) throw new NotSupportedException("positions are 32-bit in the prototype");
        string full = Path.GetFullPath(o.ScratchDir);
        if (Directory.Exists(full) && Directory.EnumerateFileSystemEntries(full).Any()) throw new ArgumentException($"the scratch directory '{o.ScratchDir}' is not empty", nameof(o));
        EndToEndParameters p = ScaleRules.Resolve(n, o);
        EndToEndFiles files = EndToEndFiles.Of(o.OutPrefix);
        string dirPh1 = Path.Combine(full, "ph1"), dirPh3 = Path.Combine(full, "ph3"), dirPh4 = Path.Combine(full, "ph4"), dirPh5 = Path.Combine(full, "ph5");
        string ranksPath = Path.Combine(full, "ranks.bin");
        Directory.CreateDirectory(full);
        IoCounters.Reset();
        MappedIo.Reset();
        ScratchTracker.Reset();
        var budgets = new List<ArenaBudget>();
        ArenaBudget Budget(long cap)
        {
            var b = new ArenaBudget(cap);
            budgets.Add(b);
            return b;
        }

        bool measure = RssSampler.IsSupported;
        using var whole = new RssSampler(10);
        if (measure) whole.Start();
        var phases = new List<PhaseRecord>();
        EndToEndResult? result = null;
        try
        {
            // ---- PH0
            PrepareResult prep = default!;
            Phase(o, phases, dirty, measure, PhasePrepare, () =>
            {
                prep = Prepare.Run(o.InputPath, files.Packed, p.PrepareBlockSymbols, p.WindowBytes);
                File.WriteAllText(files.PackedMeta, prep.MetaText(Path.GetFileName(o.InputPath)));
            });
            result = new EndToEndResult { Parameters = p, Prepare = prep, Files = files };

            using var mapped = new MappedFile(files.Packed, writable: false);
            using var packed = new PackedText(mapped.Pointer, n, prep.Bits, mapped, ownsOwner: false);
            DifferenceCover cover = o.ExactPeriod ? DifferenceCover.Create(o.Period) : DifferenceCover.Create(DifferenceCover.NearestSingerPeriod(o.Period));
            result.CoverPeriod = cover.Period;
            result.CoverResidues = cover.Residues.Length;

            // ---- PH1
            ArenaBudget dcsBudget = Budget(o.ArenaCapBytes);
            Phase(o, phases, dirty, measure, PhaseDcs, () =>
            {
                var dcs = new ExtDcsOptions
                {
                    ScratchDir = dirPh1,
                    BucketCapacity = p.DcsBucketCapacity,
                    BaseCase = p.DcsBaseCase,
                    ChunkSlots = p.DcsChunkSlots,
                    SortLeafCapacity = p.DcsSortCapacity,
                    Threads = o.Threads,
                    Seed = 1,
                    FanOut = p.DcsFanOut,
                    ReservoirSize = p.DcsReservoir,
                    RunSpillBytes = p.DcsRunSpillBytes,
                    IoBufferBytes = p.DcsIoBytes,
                    MergeBufferBudget = p.DcsMergeBudget,
                    WindowBytes = p.WindowBytes,
                    Budget = dcsBudget,
                    TouchBuffers = o.TouchBuffers,
                };
                result.DcsOptions = dcs;
                result.Dcs = ExternalDcs.BuildRankFile(packed, cover, ranksPath, dcs);
                ScratchTracker.Add(new FileInfo(ranksPath).Length); // the rank file is scratch of the pipeline from here on
            });

            // ---- PH2 / PH3
            ArenaBudget saBudget = Budget(o.ArenaCapBytes);
            Phase(o, phases, dirty, measure, PhaseSa, () =>
            {
                var build = new BuildOptions
                {
                    LeafCapacity = p.BucketCapacity,
                    FanOut = p.FanOut,
                    Threads = o.Threads,
                    Seed = 1,
                    Oversample = p.Oversample,
                    ScratchDir = dirPh3,
                    SpillBytes = p.SpillBytes,
                    Budget = saBudget,
                };
                result.BuildOptions = build;
                Directory.CreateDirectory(dirPh3);
                long rankBytes = new FileInfo(ranksPath).Length;
                using (var ranks = new DcsRankTable(cover, new MappedFile(ranksPath, writable: false), n))
                using (var sink = new MappedSaSink(files.Sa, n + 1))
                {
                    result.Build = EmBuilder.Run(packed, build, sink, ranks);
                }
                File.Delete(ranksPath); // the rank file has served its only consumer
                ScratchTracker.Remove(rankBytes);
            });

            // ---- PH4 + PH5 (fused through the B3 hook)
            ArenaBudget lcpBudget = Budget(o.ArenaCapBytes);
            ArenaBudget childBudget = Budget(o.ChildArenaCapBytes);
            Phase(o, phases, dirty, measure, PhaseLcpChild, () =>
            {
                var lcp = new ExtLcpOptions
                {
                    ScratchDir = dirPh4,
                    ChunkPositions = p.LcpChunk,
                    Threads = o.Threads,
                    Workers = p.LcpWorkers,
                    ReadBufferBytes = p.LcpReadBytes,
                    TupleBufferBytes = p.LcpTupleBytes,
                    WindowBytes = p.WindowBytes,
                    Mode = LcpMode.Distributed,
                    TouchBuffers = o.TouchBuffers,
                    Budget = lcpBudget,
                    ArenaGuards = o.ArenaGuards,
                };
                var child = new ChildBuildOptions
                {
                    ScratchDir = dirPh5,
                    BlockEntries = p.BlockEntries,
                    WindowBytes = p.ChildWindowBytes,
                    WindowSlots = p.ChildSlots,
                    FlushWindows = true,
                    SpillFlushBytes = p.SpillFlushBytes,
                    ReadBufferBytes = p.LcpReadBytes,
                    Codec = o.Codec,
                    TouchBuffers = o.TouchBuffers,
                    Budget = childBudget,
                    ArenaGuards = o.ArenaGuards,
                };
                result.LcpOptions = lcp;
                result.ChildOptions = child;
                using var scan = new E5.LazyScan(n, new ChildScanOutputs { ChildPath = files.Child }, child, null);
                result.Lcp = ExternalLcp.Build(files.Sa, packed, files.Lcp, files.Overflow, lcp, scan);
                result.Scan = scan.Complete();
            });
            File.WriteAllText(files.Stats, StatsText(result));
            if (o.SyncOutputs) Phase(o, phases, dirty, measure, PhaseSync, () => SyncAll(files));
        }
        catch
        {
            DeleteScratch(full);
            throw;
        }
        finally
        {
            if (measure) whole.Stop();
        }

        EndToEndResult r = result!;
        r.Phases.AddRange(phases);
        RssStats all = measure ? whole.Stop() : default;
        r.RssAnonPeakMb = all.MaxRssAnonMb;
        r.ScratchPeakBytes = ScratchTracker.Peak;
        r.ScratchLiveBytes = ScratchTracker.Live;
        r.ArenaLiveBytes = budgets.Sum(b => b.LiveBytes);
        r.ArenaPeakBytes = budgets.Max(b => b.PeakBytes);
        DeleteScratch(full);
        r.ScratchDirectoryExists = Directory.Exists(full);
        r.ScratchEntriesLeft = r.ScratchDirectoryExists ? Directory.GetFileSystemEntries(full).Length : 0;
        foreach (string f in files.All) r.FileBytes[f] = new FileInfo(f).Length;
        CheckLengths(r, n);
        return r;
    }

    /// <summary>The line-oriented statistics file: NODE_COUNT, the internal nodes, LRS_DEPTH and LRS_RANK (the keys of <c>e5-oracle</c>'s <c>.stats</c>) and the shape of the output.</summary>
    internal static string StatsText(EndToEndResult r) => string.Create(CultureInfo.InvariantCulture,
        $"nodecount={r.Scan.NodeCount}\ninternal={r.Scan.InternalCount}\nlrsdepth={r.Scan.LrsDepth}\nlrsrank={r.Scan.LrsRank}\nn={r.Prepare.N}\nsigma={r.Prepare.Sigma}\nlcpwidth={LcpWidth.Bytes(r.Lcp.Kind)}\noverflow={r.Lcp.OverflowCount}\nmaxlcp={r.Scan.MaxLcp}\n");

    private static void Phase(EndToEndOptions o, List<PhaseRecord> phases, DirtySampler? dirty, bool measure, string name, Action body)
    {
        o.Trace?.Invoke("begin", name);
        long r0 = IoCounters.ReadBytes, w0 = IoCounters.WriteBytes, m0 = MappedIo.WriteBytes;
        ProcIoSnapshot io0 = ProcIoSnapshot.Read();
        using var rss = new RssSampler(10);
        if (measure) rss.Start();
        dirty?.TakePhaseMax();
        long t0 = Stopwatch.GetTimestamp();
        body();
        double seconds = (Stopwatch.GetTimestamp() - t0) / (double)Stopwatch.Frequency;
        RssStats rs = measure ? rss.Stop() : default;
        (long dirtyKb, long processKb) = dirty?.TakePhaseMax() ?? default;
        ProcIoSnapshot io = ProcIoSnapshot.Read().Since(io0);
        phases.Add(new PhaseRecord(name, seconds, IoCounters.ReadBytes - r0, IoCounters.WriteBytes - w0, MappedIo.WriteBytes - m0, io, rs.MaxRssAnonMb, rs.MaxRssFileMb, ScratchTracker.Peak, ScratchTracker.Live, dirtyKb, processKb));
        o.Trace?.Invoke("end", name);
    }

    private static void SyncAll(EndToEndFiles files)
    {
        foreach (string f in files.All)
        {
            using var fs = new FileStream(f, FileMode.Open, FileAccess.ReadWrite);
            fs.Flush(true);
        }
    }

    private static void DeleteScratch(string dir)
    {
        try
        {
            if (Directory.Exists(dir)) Directory.Delete(dir, true);
        }
        catch (IOException)
        {
            // reported through ScratchEntriesLeft
        }
    }

    /// <summary>Every output has its final length: the packed text with its padding, SA and CHILD of <c>4 (n + 1)</c> bytes, the base LCP of <c>(n + 1)</c> entries of the chosen width, the overflow pairs.</summary>
    private static void CheckLengths(EndToEndResult r, long n)
    {
        void Expect(string path, long bytes, string what)
        {
            long got = r.FileBytes[path];
            if (got != bytes) throw new InvalidOperationException($"{what} has {got} bytes, expected {bytes}");
        }
        Expect(r.Files.Packed, PackedText.TotalBytes(n, r.Prepare.Bits), "the packed text");
        Expect(r.Files.Sa, 4 * (n + 1), "the suffix array");
        Expect(r.Files.Lcp, LcpWidth.Bytes(r.Lcp.Kind) * (n + 1), "the base LCP");
        Expect(r.Files.Overflow, 8 * r.Lcp.OverflowCount, "the overflow list");
        Expect(r.Files.Child, 4 * (n + 1), "the child table");
    }
}
