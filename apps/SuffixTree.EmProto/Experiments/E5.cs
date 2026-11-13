using System.Diagnostics;
using System.Globalization;
using System.Reflection;
using System.Runtime.InteropServices;
using SuffixTree.Esa;
using SuffixTree.EmProto.Build;
using SuffixTree.EmProto.Dcs;
using SuffixTree.EmProto.Harness;
using SuffixTree.EmProto.Ph4;
using SuffixTree.EmProto.Ph5;
using SuffixTree.EmProto.Text;

namespace SuffixTree.EmProto.Experiments;

/// <summary>
/// E5 of ADR-004 section 5 (owner decision 2026-11-11, section 6c; the yes for E5 of 2026-11-13): PH5, the lcp-interval scan with the external stack, producing the child table in layout A (Abouelhoda up / down / nextlIndex
/// in the single-array form, writes behind the cursor through bounded mapped windows) and in layout F (PSV forward, NSV by a backward scan of the LCP file, every write sequential), NODE_COUNT and LRS_DEPTH / LRS_RANK.
/// <c>e5</c> reads the LCP from a given compact LCP file (<c>--lcp &lt;prefix&gt;</c>, the E4 output: <c>.lcp</c> base array, <c>.ovf</c> overflow pairs) or, with <c>--fused --sa</c>, runs PH4 (base mode) with the scan as the
/// consumer of its B3 hook; <c>e5-oracle</c> writes the reference files from the in-RAM structures of the product (the oracle of gate G1 above the size that the tests reach, <c>--verify</c> compares byte for byte).
/// Measured: time per phase, RssAnon, spill volume and I/O, stack height, windows, and the dirty pages (system-wide <c>/proc/meminfo</c> Dirty every 20 ms with the baseline subtracted, and the per-process <c>smaps_rollup</c> view).
/// </summary>
internal static unsafe class E5
{
    /// <summary>The names of the output files of a run, <c>&lt;prefix&gt;.cld</c> (layout A), <c>.psv</c> and <c>.nsv</c> (layout F).</summary>
    internal static string ChildFile(string prefix) => prefix + ".cld";

    internal static string PsvFile(string prefix) => prefix + ".psv";

    internal static string NsvFile(string prefix) => prefix + ".nsv";

    internal static int Build(Options o, Report r)
    {
        string textPath = o.Require("text");
        string outPrefix = o.Require("out");
        string scratch = o.Require("scratch");
        string caseName = o.Str("case", "e5");
        string layout = o.Str("layout", "both");
        if (layout is not ("a" or "f" or "both" or "stats")) throw new ArgumentException($"unknown --layout '{layout}' (a, f, both or stats)");
        ChildCodec codec = o.Str("codec", "tagged") switch
        {
            "tagged" => ChildCodec.Tagged,
            "packed" => ChildCodec.Packed,
            var other => throw new ArgumentException($"unknown --codec '{other}' (tagged or packed)"),
        };
        bool fused = o.Flag("fused");
        (long n, int bits) = E3.ReadMeta(textPath);
        long ph5Budget = o.Long("ph5-budget-mb", 64) << 20;
        string spillDir = o.Get("spill-dir") ?? scratch;
        var options = new ChildBuildOptions
        {
            ScratchDir = spillDir,
            BlockEntries = o.Int("block-entries", 256 * 1024),
            WindowBytes = o.Int("window-kb", 1024) * 1024,
            WindowSlots = o.Int("slots", 8),
            FlushWindows = !o.Flag("no-flush"),
            SpillFlushBytes = o.Long("spill-flush-mb", 16) << 20,
            ReadBufferBytes = o.Int("io-kb", 1024) * 1024,
            Codec = codec,
            TouchBuffers = o.Flag("touch-buffers"),
            Budget = new ArenaBudget(ph5Budget),
            ParallelBackward = o.Flag("parallel-nsv"),
        };
        bool wantA = layout is "a" or "both";
        bool wantF = layout is "f" or "both";
        var outputs = new ChildBuildOutputs
        {
            Child = wantA ? ChildFile(outPrefix) : null,
            Psv = wantF ? PsvFile(outPrefix) : null,
            Nsv = wantF ? NsvFile(outPrefix) : null,
        };
        Directory.CreateDirectory(scratch);
        Directory.CreateDirectory(spillDir);
        foreach (string entry in Directory.GetFileSystemEntries(scratch))
        {
            if (Directory.Exists(entry)) Directory.Delete(entry, true);
            else File.Delete(entry);
        }
        string lcpPrefix = fused ? outPrefix + ".ph4" : o.Require("lcp");

        r.Row("E5", caseName, "n", n, "symbols", $"b={bits} bits, layout={layout}, codec={codec}, fused with PH4={fused}");
        r.Row("E5", caseName, "parameters", string.Create(CultureInfo.InvariantCulture, $"parallel NSV={options.ParallelBackward} block={options.BlockEntries} entries ({options.BlockEntries * 16L >> 10} KiB) window={options.WindowBytes} slots={options.WindowSlots} flush windows={options.FlushWindows} spill flush={options.SpillFlushBytes} read={options.ReadBufferBytes} spill dir={spillDir}"), "", "entry 16 B (l, lb, first, last); two blocks in RAM");

        IoCounters.Reset();
        MappedIo.Reset();
        ScratchTracker.Reset();
        ProcIoSnapshot io0 = ProcIoSnapshot.Read();
        using var rss = new RssSampler(10, keepSeries: false);
        using var dirty = new DirtySampler(20, filePrefix: Path.GetFullPath(outPrefix), keepSeries: o.Get("series") is not null);
        dirty.Start();
        rss.Start();
        long t0 = Stopwatch.GetTimestamp();
        ChildBuildStats stats;
        ExtLcpStats? ph4 = null;
        double ph4Seconds = 0;
        if (fused)
        {
            string saPath = o.Require("sa");
            ExtLcpOptions lo = Ph4Options(o, n, scratch, out string description);
            r.Row("E5", caseName, "ph4_parameters", description, "", "PH4 base mode; the scan is the consumer of its B3 hook");
            using var mapped = new MappedFile(textPath, writable: false);
            using var packed = new PackedText(mapped.Pointer, n, bits, mapped, ownsOwner: false);
            using var lazy = new LazyScan(n, new ChildScanOutputs { ChildPath = outputs.Child, PsvPath = outputs.Psv }, options, dirty);
            long p0 = Stopwatch.GetTimestamp();
            ph4 = ExternalLcp.Build(saPath, packed, lcpPrefix + ".lcp", lcpPrefix + ".ovf", lo, lazy);
            ChildScanResult forward = lazy.Complete();
            ph4Seconds = (Stopwatch.GetTimestamp() - p0) / (double)Stopwatch.Frequency;
            NsvResult? backward = null;
            if (outputs.Nsv is not null)
            {
                using var file = new CompactLcpFile(lcpPrefix + ".lcp", lcpPrefix + ".ovf", n, options.ReadBufferBytes);
                backward = NsvScan.Run(n, file, outputs.Nsv, options);
            }
            stats = new ChildBuildStats { Forward = forward, Backward = backward, ForwardWallSeconds = ph4Seconds, LcpBytes = new FileInfo(lcpPrefix + ".lcp").Length + new FileInfo(lcpPrefix + ".ovf").Length };
        }
        else
        {
            dirty.Mark();
            stats = ExternalChild.BuildFromFiles(n, lcpPrefix + ".lcp", lcpPrefix + ".ovf", outputs, options);
        }
        double total = (Stopwatch.GetTimestamp() - t0) / (double)Stopwatch.Frequency;
        RssStats rs = rss.Stop();
        ProcIoSnapshot io = ProcIoSnapshot.Read().Since(io0);
        long f0 = Stopwatch.GetTimestamp();
        foreach (string? p in new[] { outputs.Child, outputs.Psv, outputs.Nsv })
        {
            if (p is null) continue;
            using var fs = new FileStream(p, FileMode.Open, FileAccess.ReadWrite);
            fs.Flush(true);
        }
        double fsync = (Stopwatch.GetTimestamp() - f0) / (double)Stopwatch.Frequency;
        DirtyStats ds = dirty.Stop();
        if (o.Get("series") is { } seriesPath)
            File.WriteAllLines(seriesPath, dirty.Series.Select(s => string.Create(CultureInfo.InvariantCulture, $"{s.Seconds:F3}\t{s.DirtyKb}\t{s.ProcessKb}")));
        Report(r, caseName, o, n, options, stats, outputs, total, fsync, rs, io, ds, ph4, fused);
        string ph4Dir = Path.Combine(scratch, "ph4");
        if (Directory.Exists(ph4Dir) && !Directory.EnumerateFileSystemEntries(ph4Dir).Any()) Directory.Delete(ph4Dir);
        string[] left = Directory.GetFileSystemEntries(scratch).Concat(Directory.GetFileSystemEntries(spillDir)).Distinct().ToArray();
        bool ok = Healthy(rs.MaxRssAnonMb, o.Double("rss-cap-mb", 256), left.Length, ScratchTracker.Live, options.Budget!.LiveBytes);
        if (o.Get("verify") is { } oracle) ok &= Verify(r, caseName, oracle, outputs, codec, stats.Forward, n);
        return ok ? 0 : 1;
    }

    /// <summary>The exit code of a run before the comparison with the oracle: the RssAnon cap held and no scratch entry, tracked scratch byte or arena byte is left behind.</summary>
    internal static bool Healthy(double rssAnonMb, double capMb, int scratchEntriesLeft, long scratchLiveBytes, long arenaLiveBytes)
        => rssAnonMb <= capMb && scratchEntriesLeft == 0 && scratchLiveBytes == 0 && arenaLiveBytes == 0;

    private static void Report(Report r, string c, Options o, long n, ChildBuildOptions options, ChildBuildStats st, ChildBuildOutputs outputs, double total, double fsync, RssStats rs, ProcIoSnapshot io, DirtyStats ds, ExtLcpStats? ph4, bool fused)
    {
        ChildScanResult f = st.Forward;
        NsvResult? b = st.Backward;
        double entries = n + 1.0;
        r.Row("E5", c, "ph5_s", total, "s", fused ? "wall of PH4 + PH5 fused (ph5_scan_s is the share inside the hook), plus the NSV pass" : "wall of the forward scan including the LCP read, plus the NSV pass when asked for (a second pass, or concurrent with --parallel-nsv; without fsync)");
        r.Row("E5", c, "forward_s", st.ForwardWallSeconds, "s", fused ? "PH4 with the hook (PH4 + forward scan)" : "forward scan including reading the LCP file");
        r.Row("E5", c, "scan_s", f.Seconds, "s", "time inside Consume + Complete of the forward scan (without reading the LCP)");
        r.Row("E5", c, "scan_ns_per_entry", f.Seconds * 1e9 / entries, "ns", "scan_s / (n + 1)");
        r.Row("E5", c, "forward_ns_per_entry", st.ForwardWallSeconds * 1e9 / entries, "ns", "forward_s / (n + 1)");
        if (b is not null)
        {
            r.Row("E5", c, "nsv_s", b.Seconds, "s", "backward scan of the LCP file: NSV written in descending rank order");
            r.Row("E5", c, "nsv_ns_per_entry", b.Seconds * 1e9 / entries, "ns", "nsv_s / (n + 1)");
        }
        r.Row("E5", c, "fsync_s", fsync, "s", "flush of the output files after the run");
        if (ph4 is not null) r.Row("E5", c, "ph4_width", ph4.Kind.ToString(), "", $"PH4 in the same process: B1 {ph4.B1Seconds:F2} s, B2 {ph4.B2Seconds:F2} s, B3 {ph4.B3Seconds:F2} s (B3 contains the scan)");
        r.Row("E5", c, "rssanon_peak", rs.MaxRssAnonMb, "MB", $"max RssAnon of {rs.Samples} samples at 10 ms; gate G2 <= 256 MB incl. the .NET runtime and thread stacks (max threads {rs.MaxThreads})");
        r.Row("E5", c, "rssfile_peak", rs.MaxRssFileMb, "MB", "mapped file pages resident (windows, text): reported, not counted");
        r.Row("E5", c, "vmhwm", rs.MaxVmHwmKb * 1024.0 / 1e6, "MB", "VmHWM (anonymous + file-backed resident), reported only");
        GCMemoryInfo gc = GC.GetGCMemoryInfo();
        r.Row("E5", c, "gc_heap", gc.HeapSizeBytes / 1e6, "MB", $"managed heap at the end (DOTNET_GCHeapHardLimit = {Environment.GetEnvironmentVariable("DOTNET_GCHeapHardLimit") ?? "not set"}); collections gen0/1/2 = {GC.CollectionCount(0)}/{GC.CollectionCount(1)}/{GC.CollectionCount(2)}");
        r.Row("E5", c, "arena_peak", options.Budget!.PeakBytes / 1e6, "MB", $"peak native buffers charged to the ArenaBudget (cap {options.Budget.CapBytes >> 20} MiB), live at the end {options.Budget.LiveBytes}");
        r.Row("E5", c, "max_lcp", f.MaxLcp, "symbols", "largest LCP value");
        r.Row("E5", c, "stack_peak", f.StackPeak, "entries", $"peak height of the forward stack (<= 1 + max LCP = {f.MaxLcp + 1}); {f.StackPeak * 16.0 / 1e6:F2} MB at 16 B per entry");
        if (b is not null) r.Row("E5", c, "nsv_stack_peak", b.StackPeak, "entries", $"peak height of the NSV stack (8 B per entry): {b.StackPeak * 8.0 / 1e6:F2} MB");
        SpillStats s = f.Spill;
        long expected = 16L * (f.MaxLcp + 1);
        r.Row("E5", c, "spill_written", s.BytesWritten, "bytes", $"forward stack: {s.BlocksWritten} blocks written, {s.BlocksRead} read back ({s.BytesRead} bytes), {s.CleanEvictions} clean evictions, {s.Flushes} fsyncs; peak file {s.PeakFileBytes} bytes");
        r.Row("E5", c, "spill_expected", expected, "bytes", $"16 B x (max LCP + 1); spill_written / expected = {(expected == 0 ? 0 : s.BytesWritten / (double)expected):F6}; difference {(expected - s.BytesWritten) / 16.0:F0} entries = {(expected - s.BytesWritten) / (16.0 * options.BlockEntries):F2} blocks of {options.BlockEntries}");
        if (b is not null) r.Row("E5", c, "nsv_spill_written", b.Spill.BytesWritten, "bytes", $"NSV stack: {b.Spill.BlocksWritten} blocks written, {b.Spill.BlocksRead} read back; expected 8 B x (max LCP + 1) = {8L * (f.MaxLcp + 1)}");
        // the spill file counts in write_bytes (its pages are dirtied) unless it lives on a memory file system (the isolation runs put it on tmpfs)
        bool spillOnDisk = !string.Equals(new DriveInfo(options.ScratchDir).DriveFormat, "tmpfs", StringComparison.Ordinal);
        long spillTotal = s.BytesWritten + (b?.Spill.BytesWritten ?? 0);
        long spillCounted = spillOnDisk ? spillTotal : 0;
        WindowStats cw = f.ChildWindows;
        if (outputs.Child is not null) r.Row("E5", c, "child_windows", cw.Opened, "n", $"CHILD: {cw.Opened} windows mapped ({cw.DistinctWindows} distinct, {cw.Reopened} mapped again after release), {cw.Evictions} evictions, {cw.Flushes} flushes covering {cw.FlushedBytes} bytes, {cw.Writes} writes ({cw.Writes / entries:F3} per rank)");
        if (outputs.Psv is not null) r.Row("E5", c, "psv_windows", f.PsvWindows.Opened, "n", $"PSV: {f.PsvWindows.Opened} windows ({f.PsvWindows.Reopened} again), {f.PsvWindows.Flushes} flushes");
        if (b is not null) r.Row("E5", c, "nsv_windows", b.Windows.Opened, "n", $"NSV: {b.Windows.Opened} windows ({b.Windows.Reopened} again), {b.Windows.Flushes} flushes");
        r.Row("E5", c, "io_explicit_read", IoCounters.ReadBytes, "bytes", $"{IoCounters.ReadBytes / entries:F3} B/char: the LCP file(s) (forward and backward) and the spill read-back ({IoCounters.ReadCalls} calls)");
        r.Row("E5", c, "io_explicit_write", IoCounters.WriteBytes, "bytes", $"{IoCounters.WriteBytes / entries:F3} B/char: the spill file(s) ({IoCounters.WriteCalls} calls)" + (fused ? " plus the PH4 tuple files" : string.Empty));
        r.Row("E5", c, "blocklayer_read", io.ReadBytes, "bytes", "/proc/self/io read_bytes (page-cache misses)");
        r.Row("E5", c, "blocklayer_write", io.WriteBytes, "bytes", "/proc/self/io write_bytes: bytes of pages that went from clean to dirty (every re-dirtying of a flushed page counts again)");
        long final = 0;
        if (outputs.Child is not null) final += new FileInfo(outputs.Child).Length;
        if (outputs.Psv is not null) final += new FileInfo(outputs.Psv).Length;
        if (outputs.Nsv is not null) final += new FileInfo(outputs.Nsv).Length;
        r.Row("E5", c, "output_bytes", final, "bytes", $"{final / entries:F3} B/char: CHILD (4 B/char) for layout A, PSV + NSV (8 B/char) for layout F");
        long dirtied = io.WriteBytes - spillCounted - (fused && ph4 is not null ? IoCounters.WriteBytes - spillTotal + MappedIo.WriteBytes : 0);
        r.Row("E5", c, "wa_dirtied", final == 0 ? 0 : dirtied / (double)final, "x", $"write amplification, definition 1: (write_bytes of /proc/self/io - the spill bytes written {spillCounted}{(spillOnDisk ? string.Empty : " (the spill is on tmpfs: not counted)")}{(fused ? " - the PH4 explicit writes and its mapped LCP output" : string.Empty)}) / final size of the output files; 1.0 = every page was dirtied once; the device-level definition 2 is taken by the run script (sectors written by the block device after sync)");
        r.Row("E5", c, "dirty_peak", ds.MaxDeltaKb * 1024.0 / 1e6, "MB", $"max /proc/meminfo Dirty minus the baseline ({ds.BaselineKb * 1024.0 / 1e6:F2} MB, median of 5 samples before the start) over {ds.Samples} samples at 20 ms; absolute max {ds.MaxDirtyKb * 1024.0 / 1e6:F2} MB; gate G3 <= 64 MB");
        r.Row("E5", c, "dirty_mean", ds.MeanDeltaKb * 1024.0 / 1e6, "MB", "mean of Dirty minus the baseline over the run");
        r.Row("E5", c, "dirty_peak_after_mark", ds.MaxDeltaAfterMarkKb * 1024.0 / 1e6, "MB", $"max Dirty after the PH5 start minus Dirty at that moment ({ds.MarkBaselineKb * 1024.0 / 1e6:F2} MB)");
        r.Row("E5", c, "writeback_peak", ds.MaxWritebackKb * 1024.0 / 1e6, "MB", "max /proc/meminfo Writeback");
        r.Row("E5", c, "proc_dirty_peak", ds.MaxProcessDirtyKb * 1024.0 / 1e6, "MB", "max Private_Dirty + Shared_Dirty of the mappings of /proc/self/smaps whose path is an output file of this run: the dirty pages this process has mapped in them (the per-process view; anonymous memory is not in it)");
        r.Row("E5", c, "node_count", f.NodeCount, "n", $"NODE_COUNT = 1 + I + n + 1 with I = {f.InternalCount} internal nodes");
        r.Row("E5", c, "lrs", f.LrsDepth, "symbols", $"LRS_DEPTH {f.LrsDepth}, LRS_RANK {f.LrsRank}");
    }

    private static bool Verify(Report r, string c, string oracle, ChildBuildOutputs outputs, ChildCodec codec, ChildScanResult f, long n)
    {
        bool ok = true;
        void Compare(string what, string? produced, string reference)
        {
            if (produced is null) return;
            (bool equal, long at) = E3b.CompareFiles(produced, reference);
            r.Row("E5", c, what + "_equals_oracle", equal ? "yes" : "NO", "", equal ? $"byte-identical to {Path.GetFileName(reference)} ({new FileInfo(produced).Length} bytes)" : $"first difference at byte {at}");
            ok &= equal;
        }
        Compare("child", outputs.Child, oracle + (codec == ChildCodec.Tagged ? ".cldT" : ".cldP"));
        Compare("psv", outputs.Psv, oracle + ".psv");
        Compare("nsv", outputs.Nsv, oracle + ".nsv");
        Dictionary<string, string> expected = File.ReadAllLines(oracle + ".stats").Select(l => l.Split('=')).ToDictionary(p => p[0], p => p[1]);
        bool stats = expected["nodecount"] == f.NodeCount.ToString(CultureInfo.InvariantCulture)
            && expected["internal"] == f.InternalCount.ToString(CultureInfo.InvariantCulture)
            && expected["lrsdepth"] == f.LrsDepth.ToString(CultureInfo.InvariantCulture)
            && expected["lrsrank"] == f.LrsRank.ToString(CultureInfo.InvariantCulture);
        r.Row("E5", c, "stats_equal_oracle", stats ? "yes" : "NO", "", $"NODE_COUNT {f.NodeCount} (oracle {expected["nodecount"]}), internal {f.InternalCount} ({expected["internal"]}), LRS_DEPTH {f.LrsDepth} ({expected["lrsdepth"]}), LRS_RANK {f.LrsRank} ({expected["lrsrank"]}); in-RAM stack peak {expected["stackpeak"]}");
        return ok && stats;
    }

    /// <summary>PH4 options for the fused mode, as <see cref="E4.Build"/> makes them (base mode, scale-down by <c>--scale-n</c>).</summary>
    private static ExtLcpOptions Ph4Options(Options o, long n, string scratch, out string description)
    {
        int threads = o.Int("threads", 4);
        long budgetCap = o.Long("budget-mb", 160) << 20;
        bool scaled = o.Get("scale-n") is not null;
        double f = scaled ? n / (double)o.Long("scale-n", E4.ReferenceN) : 1.0;
        long chunk = scaled ? Math.Max(3, (long)Math.Ceiling(E4.ChunkRef * f)) : E4.ChunkRef;
        if (o.Get("chunk") is { } explicitChunk) chunk = long.Parse(explicitChunk, CultureInfo.InvariantCulture);
        ExtLcpOptions options = new()
        {
            ScratchDir = Path.Combine(scratch, "ph4"),
            ChunkPositions = chunk,
            Threads = threads,
            Mode = LcpMode.Distributed,
            Budget = new ArenaBudget(budgetCap),
        };
        ExtLcpOptions reference = options with { ChunkPositions = E4.ChunkRef };
        long perTask = ExternalLcp.WorkerBytes(E4.ReferenceN, reference);
        int workers = o.Int("b2-workers", (int)Math.Clamp(budgetCap / Math.Max(1, perTask), 1, threads));
        description = string.Create(CultureInfo.InvariantCulture, $"R={chunk} P={(n + chunk - 1) / chunk} threads={threads} B2 workers={workers} arena cap {budgetCap >> 20} MiB");
        return options with { Workers = workers };
    }

    /// <summary>Creates the scan at the first chunk of the hook, so that its stack blocks and windows are not allocated while PH4 is in its B1 / B2 peak.</summary>
    private sealed class LazyScan(long n, ChildScanOutputs outputs, ChildBuildOptions options, DirtySampler? sampler) : ILcpChunkConsumer, IDisposable
    {
        private ChildScan? _scan;

        public void Consume(long firstRank, ReadOnlySpan<uint> lcp)
        {
            if (_scan is null)
            {
                sampler?.Mark();
                _scan = ChildScan.Create(n, outputs, options);
            }
            _scan.Consume(firstRank, lcp);
        }

        internal ChildScanResult Complete() => (_scan ?? throw new InvalidOperationException("PH4 delivered no LCP chunk")).Complete();

        public void Dispose() => _scan?.Dispose();
    }

    /// <summary>
    /// <c>e5-oracle</c>: the reference files of a text from the in-RAM structures of the product (needs n below 2^31): the compact LCP of the E4 files, the child table of the tagged and of the packed layout of
    /// <see cref="ChildTable"/> (raw field arrays, <c>.cldT</c> and <c>.cldP</c>), the previous and next smaller values of the <see cref="LcpPyramid"/> (<c>.psv</c>, <c>.nsv</c>: 32-bit, 0xFFFFFFFF for none on the left, n + 1 for none on the right) and the header
    /// statistics (<c>.stats</c>).
    /// </summary>
    internal static int Oracle(Options o, Report r)
    {
        string textPath = o.Require("text");
        string lcpPrefix = o.Require("lcp");
        string outPrefix = o.Require("out");
        (long n, _) = E3.ReadMeta(textPath);
        if (n >= int.MaxValue) throw new ArgumentException("the in-RAM oracle needs n < 2^31");
        int length = checked((int)n + 1);
        long baseLength = new FileInfo(lcpPrefix + ".lcp").Length;
        int width = (int)(baseLength / length);
        byte[] baseBytes = File.ReadAllBytes(lcpPrefix + ".lcp");
        ReadOnlySpan<uint> pairs = MemoryMarshal.Cast<byte, uint>(File.ReadAllBytes(lcpPrefix + ".ovf"));
        var ranks = new int[pairs.Length / 2];
        var values = new int[pairs.Length / 2];
        for (int i = 0; i < ranks.Length; i++)
        {
            ranks[i] = (int)pairs[2 * i];
            values[i] = (int)pairs[(2 * i) + 1];
        }
        CompactLcp lcp = width switch
        {
            1 => CompactLcp.FromParts(LcpKind.U8, length, baseBytes, null, null, ranks, values),
            2 => CompactLcp.FromParts(LcpKind.U16, length, null, MemoryMarshal.Cast<byte, ushort>(baseBytes).ToArray(), null, ranks, values),
            _ => CompactLcp.FromParts(LcpKind.U32, length, null, null, MemoryMarshal.Cast<byte, int>(baseBytes).ToArray(), ranks, values),
        };
#pragma warning disable S3011 // the oracle reads the private field arrays of the in-RAM layouts
        long t0 = Stopwatch.GetTimestamp();
        EsaStatistics stats;
        {
            ChildTable tagged = ChildTable.Build(ChildLayoutKind.Tagged, (int)n, lcp, out stats);
            WriteRaw(outPrefix + ".cldT", (uint[])typeof(TaggedChild).GetField("_cld", BindingFlags.NonPublic | BindingFlags.Instance)!.GetValue(tagged.Tagged)!);
        }
        {
            ChildTable packed = ChildTable.Build(ChildLayoutKind.Packed, (int)n, lcp, out _);
            WriteRaw(outPrefix + ".cldP", (uint[])typeof(PackedChild).GetField("_cld", BindingFlags.NonPublic | BindingFlags.Instance)!.GetValue(packed.Packed)!);
        }
#pragma warning restore S3011
        double child = (Stopwatch.GetTimestamp() - t0) / (double)Stopwatch.Frequency;
        t0 = Stopwatch.GetTimestamp();
        using (LcpPyramid pyramid = LcpPyramid.Build(lcp, 32))
        {
            using var psv = new FileStream(outPrefix + ".psv", FileMode.Create, FileAccess.Write, FileShare.None, 1 << 20);
            using var nsv = new FileStream(outPrefix + ".nsv", FileMode.Create, FileAccess.Write, FileShare.None, 1 << 20);
            var pb = new uint[1 << 18];
            var nb = new uint[1 << 18];
            for (int start = 0; start < length; start += pb.Length)
            {
                int count = Math.Min(pb.Length, length - start);
                for (int i = 0; i < count; i++)
                {
                    pb[i] = (uint)pyramid.PreviousSmaller(start + i); // -1 -> 0xFFFFFFFF
                    nb[i] = (uint)pyramid.NextSmaller(start + i);
                }
                psv.Write(MemoryMarshal.AsBytes(pb.AsSpan(0, count)));
                nsv.Write(MemoryMarshal.AsBytes(nb.AsSpan(0, count)));
            }
        }
        double smaller = (Stopwatch.GetTimestamp() - t0) / (double)Stopwatch.Frequency;
        File.WriteAllLines(outPrefix + ".stats", new[]
        {
            $"nodecount={stats.NodeCount}",
            $"internal={stats.InternalCount}",
            $"lrsdepth={stats.LrsDepth}",
            $"lrsrank={stats.LrsRank}",
            $"stackpeak={stats.StackPeak}",
        });
        string c = o.Str("case", "e5-oracle");
        r.Row("E5", c, "oracle_child_s", child, "s", $"in-RAM ChildTable.Build, tagged and packed, n={n}, LCP kind {lcp.Kind}, overflow entries {lcp.OverflowCount}");
        r.Row("E5", c, "oracle_smaller_s", smaller, "s", "PSV and NSV of every rank from the LcpPyramid (fanout 32)");
        r.Row("E5", c, "oracle_stats", stats.NodeCount, "n", $"NODE_COUNT {stats.NodeCount}, internal {stats.InternalCount}, LRS_DEPTH {stats.LrsDepth}, LRS_RANK {stats.LrsRank}, in-RAM stack peak {stats.StackPeak}");
        return 0;
    }

    private static void WriteRaw(string path, uint[] data)
    {
        using var fs = new FileStream(path, FileMode.Create, FileAccess.Write, FileShare.None, 1 << 20);
        fs.Write(MemoryMarshal.AsBytes(data.AsSpan()));
    }
}
