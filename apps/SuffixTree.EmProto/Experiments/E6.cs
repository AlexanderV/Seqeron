using System.Diagnostics;
using System.Globalization;
using SuffixTree.EmProto.Build;
using SuffixTree.EmProto.Dcs;
using SuffixTree.EmProto.Harness;
using SuffixTree.EmProto.Ph4;
using SuffixTree.EmProto.Ph5;
using SuffixTree.EmProto.Pipeline;

namespace SuffixTree.EmProto.Experiments;

/// <summary>
/// E6 of ADR-004 section 5: the first end-to-end run, PH0 to PH5, text in, packed text + suffix array + compact LCP + CHILD (layout A, packed codec) + statistics out, in one process with bounded RAM, every big output streamed
/// through mapped windows, the scratch deleted at the end. <c>e6</c> runs <see cref="EndToEnd"/> with the scale-down rule of <see cref="ScaleRules"/> (<c>--scale-n 3000000000</c> for the 100 Mbp runs) and reports every phase;
/// <c>--verify-*</c> compare the outputs byte for byte with reference files (the packed file of <c>pack</c>, the kernel suffix array of <c>oracle --out</c>, the files of <c>lcp-oracle</c> and <c>e5-oracle</c>);
/// <c>e6-mktext</c> writes the input text of a family (setup, outside the pipeline).
/// </summary>
internal static class E6
{
    /// <summary><c>e6-mktext --family &lt;spec&gt; --out &lt;file&gt; [--data dir]</c>: the ASCII text of a family, one byte per symbol, as the benchmark inputs.</summary>
    internal static int MakeText(Options o, Report r)
    {
        string spec = o.Require("family");
        string outPath = o.Require("out");
        Gen.SymbolText text = FamilyLoader.Load(spec, o.Str("data", string.Empty));
        File.WriteAllBytes(outPath, Gen.Families.ToAscii(text));
        r.Row("E6", spec, "text_file", text.Length, "bytes", $"{text.Length} symbols, sigma={text.Sigma}, letters '{text.Letters}', written to {Path.GetFileName(outPath)}");
        return 0;
    }

    /// <summary>The command line as <see cref="EndToEndOptions"/>.</summary>
    internal static EndToEndOptions ParseOptions(Options o)
    {
        ChildCodec codec = o.Str("codec", "packed") switch
        {
            "tagged" => ChildCodec.Tagged,
            "packed" => ChildCodec.Packed,
            var other => throw new ArgumentException($"unknown --codec '{other}' (tagged or packed)"),
        };
        long? L(string key) => o.Get(key) is { } v ? long.Parse(v, CultureInfo.InvariantCulture) : null;
        int? I(string key) => o.Get(key) is { } v ? int.Parse(v, CultureInfo.InvariantCulture) : null;
        int? Kb(string key) => I(key) is { } v ? v * 1024 : null;
        return new EndToEndOptions
        {
            InputPath = o.Require("input"),
            OutPrefix = o.Require("out"),
            ScratchDir = o.Require("scratch"),
            Threads = o.Int("threads", 4),
            ScaleN = o.Long("scale-n", 0),
            Period = o.Int("period", 1024),
            ExactPeriod = o.Flag("exact"),
            Codec = codec,
            ArenaCapBytes = o.Long("budget-mb", 160) << 20,
            ChildArenaCapBytes = o.Long("ph5-budget-mb", 64) << 20,
            SyncOutputs = !o.Flag("no-sync"),
            TouchBuffers = o.Flag("touch-buffers"),
            Tuning = new EndToEndTuning
            {
                BucketCapacity = L("k"),
                FanOut = I("fan"),
                Oversample = I("oversample"),
                SpillBytes = Kb("spill-kb"),
                DcsBucketCapacity = L("k-s"),
                DcsBaseCase = L("nbase"),
                DcsChunkSlots = L("dcs-chunk"),
                DcsSortCapacity = L("sort-k"),
                DcsRunSpillBytes = Kb("run-spill-kb"),
                DcsFanOut = I("dcs-fan"),
                DcsReservoir = I("reservoir"),
                DcsIoBytes = Kb("dcs-io-kb"),
                DcsMergeBudget = L("merge-budget-mb") is { } mb ? mb << 20 : null,
                LcpChunk = L("chunk"),
                LcpWorkers = I("b2-workers"),
                LcpReadBytes = Kb("io-kb"),
                LcpTupleBytes = Kb("tuple-kb"),
                BlockEntries = I("block-entries"),
                ChildWindowBytes = Kb("child-window-kb"),
                ChildSlots = I("slots"),
                SpillFlushBytes = L("spill-flush-mb") is { } sf ? sf << 20 : null,
                WindowBytes = Kb("window-kb"),
                PrepareBlockSymbols = I("prepare-block"),
            },
        };
    }

    internal static int Run(Options o, Report r)
    {
        EndToEndOptions eo = ParseOptions(o);
        string c = o.Str("case", "e6");
        Directory.CreateDirectory(Path.GetDirectoryName(Path.GetFullPath(eo.OutPrefix))!);
        using var dirty = new DirtySampler(20, filePrefix: Path.GetFullPath(eo.OutPrefix), keepSeries: o.Get("series") is not null);
        dirty.Start();
        ProcIoSnapshot io0 = ProcIoSnapshot.Read();
        long t0 = Stopwatch.GetTimestamp();
        EndToEndResult res = EndToEnd.Run(eo, dirty);
        double wall = (Stopwatch.GetTimestamp() - t0) / (double)Stopwatch.Frequency;
        ProcIoSnapshot io = ProcIoSnapshot.Read().Since(io0);
        DirtyStats ds = dirty.Stop();
        if (o.Get("series") is { } seriesPath)
            File.WriteAllLines(seriesPath, dirty.Series.Select(s => string.Create(CultureInfo.InvariantCulture, $"{s.Seconds:F3}\t{s.DirtyKb}\t{s.ProcessKb}")));
        bool ok = Report(r, c, o, eo, res, wall, io, ds);
        ok &= Verify(r, c, o, res);
        return ok ? 0 : 1;
    }

    private static bool Report(Report r, string c, Options o, EndToEndOptions eo, EndToEndResult res, double wall, ProcIoSnapshot io, DirtyStats ds)
    {
        long n = res.Prepare.N;
        double perChar = 1.0 / n;
        EndToEndParameters p = res.Parameters;
        r.Row("E6", c, "n", n, "symbols", $"b={res.Prepare.Bits} bits, sigma={res.Prepare.Sigma}, threads={eo.Threads}, codec={eo.Codec}, X={res.CoverPeriod} (|D|={res.CoverResidues}), touch buffers={eo.TouchBuffers}");
        r.Row("E6", c, "parameters", p.ToString(), "", eo.ScaleN > 0 ? $"scale-down by n / {eo.ScaleN} = {p.Factor:E3}: every stage by the rule of E3 / E3b / E4 (ScaleRules)" : "unscaled (reference values and explicit values)");
        r.Row("E6", c, "stage_parameters", string.Create(CultureInfo.InvariantCulture, $"PH1 fan-out {p.DcsFanOut}, reservoir {p.DcsReservoir}, io {p.DcsIoBytes}, merge budget {p.DcsMergeBudget}; PH3 fan-out {p.FanOut}, oversample {p.Oversample}, spill {p.SpillBytes}; PH4 read {p.LcpReadBytes}, tuple {p.LcpTupleBytes}; PH5 block {p.BlockEntries} entries, window {p.ChildWindowBytes} x {p.ChildSlots}, spill flush {p.SpillFlushBytes}; windows {p.WindowBytes}; arena caps {eo.ArenaCapBytes >> 20} + {eo.ChildArenaCapBytes >> 20} MiB"), "");
        foreach (PhaseRecord ph in res.Phases)
        {
            string k = ph.Name;
            r.Row("E6", c, k + "_s", ph.Seconds, "s", ph.Name == EndToEnd.PhaseSync ? "fsync of all outputs (not part of the stage times)" : "wall time of the phase");
            r.Row("E6", c, k + "_io_explicit_read", ph.ExplicitReadBytes, "bytes", $"{ph.ExplicitReadBytes * perChar:F3} B/char (program counter, explicit reads)");
            r.Row("E6", c, k + "_io_explicit_write", ph.ExplicitWriteBytes, "bytes", $"{ph.ExplicitWriteBytes * perChar:F3} B/char (program counter, explicit writes)");
            r.Row("E6", c, k + "_io_mapped_write", ph.MappedWriteBytes, "bytes", $"{ph.MappedWriteBytes * perChar:F3} B/char (program counter, writes through mapped windows; CHILD windows are counted at the run level)");
            r.Row("E6", c, k + "_proc_rchar", ph.Proc.RChar, "bytes", $"/proc/self/io rchar during the phase, {ph.Proc.RChar * perChar:F3} B/char");
            r.Row("E6", c, k + "_proc_wchar", ph.Proc.WChar, "bytes", $"/proc/self/io wchar during the phase, {ph.Proc.WChar * perChar:F3} B/char");
            r.Row("E6", c, k + "_proc_read_bytes", ph.Proc.ReadBytes, "bytes", "/proc/self/io read_bytes (block layer, page-cache misses)");
            r.Row("E6", c, k + "_proc_write_bytes", ph.Proc.WriteBytes, "bytes", "/proc/self/io write_bytes (pages that went from clean to dirty)");
            r.Row("E6", c, k + "_rssanon_peak", ph.RssAnonPeakMb, "MB", "max RssAnon of the phase at 10 ms");
            r.Row("E6", c, k + "_rssfile_peak", ph.RssFilePeakMb, "MB", "mapped file pages resident (reported, not counted)");
            r.Row("E6", c, k + "_scratch_peak_cum", ph.ScratchPeakBytes, "bytes", $"{ph.ScratchPeakBytes * perChar:F3} B/char: peak of the tracked scratch bytes up to the end of this phase");
            r.Row("E6", c, k + "_scratch_live_end", ph.ScratchLiveBytes, "bytes", "tracked scratch bytes alive at the end of this phase");
            r.Row("E6", c, k + "_dirty_peak", Math.Max(0, ph.DirtyPeakKb - ds.BaselineKb) * 1024.0 / 1e6, "MB", $"max /proc/meminfo Dirty in the phase minus the baseline {ds.BaselineKb * 1024.0 / 1e6:F2} MB (system-wide: scratch included)");
            r.Row("E6", c, k + "_proc_dirty_peak", ph.ProcessDirtyPeakKb * 1024.0 / 1e6, "MB", "max Private_Dirty + Shared_Dirty of the mappings of the output files in the phase (per-process view)");
        }
        double total = res.StageSeconds;
        r.Row("E6", c, "total_s", total, "s", "PH0 + PH1 + PH2/PH3 + PH4/PH5: wall time of the stages (gate G3: <= 90 s at 100 Mbp, 4 threads, cold cache)");
        r.Row("E6", c, "sa_stage_s", res.SaStageSeconds, "s", "PH0 + PH1 + PH2/PH3 (gate G3: <= 30 s)");
        r.Row("E6", c, "wall_s", wall, "s", "EndToEnd.Run from start to end, including the clean-up and the final sync");
        r.Row("E6", c, "ns_per_symbol", total * 1e9 / n, "ns", "total_s / n");
        r.Row("E6", c, "throughput", n / total / 1e6, "Msym/s", "n / total_s");
        foreach (string ph in EndToEnd.StagePhases) r.Row("E6", c, ph + "_share", res.Phase(ph).Seconds / total, "", "share of total_s");
        r.Row("E6", c, "rssanon_peak", res.RssAnonPeakMb, "MB", "max RssAnon over the whole run at 10 ms; gate G2 <= 256 MB incl. the .NET runtime and thread stacks");
        GCMemoryInfo gc = GC.GetGCMemoryInfo();
        r.Row("E6", c, "gc_heap", gc.HeapSizeBytes / 1e6, "MB", $"managed heap at the end (DOTNET_GCHeapHardLimit = {Environment.GetEnvironmentVariable("DOTNET_GCHeapHardLimit") ?? "not set"}); collections gen0/1/2 = {GC.CollectionCount(0)}/{GC.CollectionCount(1)}/{GC.CollectionCount(2)}");
        r.Row("E6", c, "arena_peak", res.ArenaPeakBytes / 1e6, "MB", $"largest peak of the native buffers of one stage charged to its ArenaBudget (caps {eo.ArenaCapBytes >> 20} / {eo.ChildArenaCapBytes >> 20} MiB), live at the end {res.ArenaLiveBytes}");
        long child = 4 * (n + 1);
        long explicitRead = res.Phases.Where(x => x.Name != EndToEnd.PhaseSync).Sum(x => x.ExplicitReadBytes);
        long explicitWrite = res.Phases.Where(x => x.Name != EndToEnd.PhaseSync).Sum(x => x.ExplicitWriteBytes);
        long mapped = res.Phases.Where(x => x.Name != EndToEnd.PhaseSync).Sum(x => x.MappedWriteBytes) + child;
        long counted = explicitRead + explicitWrite + mapped;
        r.Row("E6", c, "io_explicit_read", explicitRead, "bytes", $"{explicitRead * perChar:F3} B/char (program counters, all phases)");
        r.Row("E6", c, "io_explicit_write", explicitWrite, "bytes", $"{explicitWrite * perChar:F3} B/char");
        r.Row("E6", c, "io_mapped_write", mapped, "bytes", $"{mapped * perChar:F3} B/char: packed text, rank file, SA, LCP through the counted windows plus CHILD {child} bytes (its windows are not counted: the file is written once, 4 (n + 1) bytes)");
        r.Row("E6", c, "io_total_counted", counted, "bytes", $"{counted * perChar:F3} B/char: explicit I/O plus mapped writes (gate G3: <= 100 B/char); reads through the mappings (text, ranks) are not in it");
        r.Row("E6", c, "proc_rchar", io.RChar, "bytes", $"{io.RChar * perChar:F3} B/char: /proc/self/io rchar over EndToEnd.Run (read syscalls)");
        r.Row("E6", c, "proc_wchar", io.WChar, "bytes", $"{io.WChar * perChar:F3} B/char: /proc/self/io wchar over EndToEnd.Run (write syscalls)");
        r.Row("E6", c, "proc_rw_total", io.RChar + io.WChar, "bytes", $"{(io.RChar + io.WChar) * perChar:F3} B/char: rchar + wchar (syscalls only; the mapped windows are not in it)");
        r.Row("E6", c, "blocklayer_read", io.ReadBytes, "bytes", $"{io.ReadBytes * perChar:F3} B/char: /proc/self/io read_bytes over the run (page-cache misses)");
        r.Row("E6", c, "blocklayer_write", io.WriteBytes, "bytes", $"{io.WriteBytes * perChar:F3} B/char: /proc/self/io write_bytes over the run (pages dirtied)");
        r.Row("E6", c, "scratch_peak", res.ScratchPeakBytes, "bytes", $"{res.ScratchPeakBytes * perChar:F3} B/char: peak of the tracked scratch files over the whole pipeline (apparent size); live at the end {res.ScratchLiveBytes}; entries left {res.ScratchEntriesLeft}; directory left {res.ScratchDirectoryExists}");
        long footprint = 0;
        foreach (string f in res.Files.All)
        {
            long len = res.FileBytes[f];
            footprint += len;
            r.Row("E6", c, "file_" + Path.GetExtension(f).TrimStart('.'), len, "bytes", $"{len * perChar:F3} B/char: {Path.GetFileName(f)}");
        }
        r.Row("E6", c, "output_bytes", footprint, "bytes", $"{footprint * perChar:F3} B/char: disk footprint of the outputs (text, SA, LCP + overflow, CHILD, metadata)");
        r.Row("E6", c, "dirty_peak", ds.MaxDeltaKb * 1024.0 / 1e6, "MB", $"max /proc/meminfo Dirty minus the baseline ({ds.BaselineKb * 1024.0 / 1e6:F2} MB) over {ds.Samples} samples at 20 ms; system-wide, scratch of PH2/PH3 and PH4 included; absolute max {ds.MaxDirtyKb * 1024.0 / 1e6:F2} MB");
        r.Row("E6", c, "dirty_mean", ds.MeanDeltaKb * 1024.0 / 1e6, "MB", "mean of Dirty minus the baseline over the run");
        r.Row("E6", c, "writeback_peak", ds.MaxWritebackKb * 1024.0 / 1e6, "MB", "max /proc/meminfo Writeback");
        r.Row("E6", c, "proc_dirty_peak", ds.MaxProcessDirtyKb * 1024.0 / 1e6, "MB", "max Private_Dirty + Shared_Dirty of the mappings of the output files of this run (per-process view of the mapped dirty pages of the outputs; gate G3 <= 64 MB)");
        ExtDcsLevelStats l0 = res.Dcs.Levels[0];
        r.Row("E6", c, "ph1_levels", res.Dcs.MaxLevel, "n", $"deepest level of the DCS recursion; level 0: sample {l0.SampleSize}, distinct names {l0.DistinctNames}, shared {100.0 * l0.NonUniqueRecords / Math.Max(1, l0.SampleSize):F2} %, leaf buckets {l0.Leaves}, split scans {l0.SplitScans}");
        BuildStats b = res.Build;
        r.Row("E6", c, "ph3_buckets", b.Buckets, "n", $"leaves={b.Leaves}, partitions={b.Partitions}, max depth={b.MaxDepth}, largest leaf={b.MaxLeaf}, retries={b.Retries}; tie share {100.0 * b.Work.TiePositions / Math.Max(1, b.Work.LeafPositions):F4} %, rank lookups {b.Work.RankLookups}");
        ExtLcpStats l = res.Lcp;
        r.Row("E6", c, "ph4_width", l.Kind.ToString(), "", $"chosen before any entry was written; escapes >= 255: {l.Escapes8} ({100.0 * l.Escapes8 / n:F4} %), >= 65535: {l.Escapes16}, overflow entries {l.OverflowCount}; B1 {l.B1Seconds:F2} s, B2 {l.B2Seconds:F2} s, B3 {l.B3Seconds:F2} s (B3 contains the PH5 scan); chunks {l.TextChunks}/{l.RankChunks}, B2 tasks {l.Workers}, rounds {l.Rounds}");
        ChildScanResult s = res.Scan;
        r.Row("E6", c, "ph5_scan_s", s.Seconds, "s", "time inside the PH5 scan (Consume + Complete): the share of the fused phase that belongs to PH5");
        r.Row("E6", c, "node_count", s.NodeCount, "n", $"NODE_COUNT = 1 + I + n + 1 with I = {s.InternalCount} internal nodes");
        r.Row("E6", c, "lrs", s.LrsDepth, "symbols", $"LRS_DEPTH {s.LrsDepth}, LRS_RANK {s.LrsRank}; max LCP {s.MaxLcp}, stack peak {s.StackPeak}, spill written {s.Spill.BytesWritten} bytes");
        r.Row("E6", c, "sha256", res.Prepare.Sha256, "", "SHA-256 of the input text (PH0)");
        bool ok = res.RssAnonPeakMb <= o.Double("rss-cap-mb", 256) && res.Clean;
        r.Row("E6", c, "clean", res.Clean ? "yes" : "NO", "", $"scratch directory removed, no tracked scratch byte, no live arena byte; RssAnon within {o.Double("rss-cap-mb", 256)} MB: {(res.RssAnonPeakMb <= o.Double("rss-cap-mb", 256) ? "yes" : "NO")}");
        return ok;
    }

    /// <summary>The comparisons of <c>--verify-pk</c>, <c>--verify-sa</c>, <c>--verify-sa-hash</c>, <c>--verify-lcp</c>, <c>--verify-child</c> and <c>--verify-stats</c>.</summary>
    private static bool Verify(Report r, string c, Options o, EndToEndResult res)
    {
        bool ok = true;
        void Compare(string what, string produced, string? reference)
        {
            if (reference is null) return;
            (bool equal, long at) = E3b.CompareFiles(produced, reference);
            r.Row("E6", c, what + "_equals_ref", equal ? "yes" : "NO", "", equal ? $"byte-identical to {Path.GetFileName(reference)} ({new FileInfo(produced).Length} bytes)" : $"first difference at byte {at}");
            ok &= equal;
        }
        Compare("pk", res.Files.Packed, o.Get("verify-pk"));
        Compare("sa", res.Files.Sa, o.Get("verify-sa"));
        if (o.Get("verify-sa-hash") is { } expectedHash)
        {
            string hash = SaHash.ComputeFile(res.Files.Sa);
            bool equal = string.Equals(hash, expectedHash, StringComparison.OrdinalIgnoreCase);
            r.Row("E6", c, "sa_hash_equals_ref", equal ? "yes" : "NO", "", $"SHA-256 of SA[1..n] {hash}");
            ok &= equal;
        }
        if (o.Get("verify-lcp") is { } lcp)
        {
            Compare("lcp", res.Files.Lcp, lcp + ".lcp");
            Compare("ovf", res.Files.Overflow, lcp + ".ovf");
        }
        Compare("child", res.Files.Child, o.Get("verify-child"));
        if (o.Get("verify-stats") is { } statsPath)
        {
            Dictionary<string, string> expected = File.ReadAllLines(statsPath).Select(x => x.Split('=')).ToDictionary(x => x[0], x => x[1]);
            Dictionary<string, string> got = File.ReadAllLines(res.Files.Stats).Select(x => x.Split('=')).ToDictionary(x => x[0], x => x[1]);
            bool stats = new[] { "nodecount", "internal", "lrsdepth", "lrsrank" }.All(k => expected[k] == got[k]);
            r.Row("E6", c, "stats_equal_ref", stats ? "yes" : "NO", "", $"NODE_COUNT {got["nodecount"]} (ref {expected["nodecount"]}), internal {got["internal"]} ({expected["internal"]}), LRS_DEPTH {got["lrsdepth"]} ({expected["lrsdepth"]}), LRS_RANK {got["lrsrank"]} ({expected["lrsrank"]})");
            ok &= stats;
        }
        return ok;
    }
}
