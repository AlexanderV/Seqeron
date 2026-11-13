using System.Diagnostics;
using System.Globalization;
using SuffixTree.Esa;
using SuffixTree.EmProto.Build;
using SuffixTree.EmProto.Dcs;
using SuffixTree.EmProto.Harness;
using SuffixTree.EmProto.Text;

namespace SuffixTree.EmProto.Experiments;

/// <summary>
/// E3 of ADR-004 section 5: scale-down emulation of the memory ratio. The packed text and the sample-rank file are read through read-only memory mappings, the position buckets go to scratch files,
/// the suffix array is written at its final offsets through bounded mapped windows, all with the bucket capacity K = n / 3000 (the 3x10^9 / 1 M ratio of the design is 3000 leaves). The commands here are
/// the pieces: <c>pack</c> (text to the packed file), <c>mkranks</c> (the rank file), <c>e3</c> (the build under measurement), <c>hash-sa</c> and <c>oracle</c> (the correctness comparison).
/// </summary>
internal static unsafe class E3
{
    internal static int Pack(Options o, Report r)
    {
        string spec = o.Require("family");
        string outPath = o.Require("out");
        Gen.SymbolText text = FamilyLoader.Load(spec, o.Str("data", string.Empty));
        int bits = Math.Max(1, PackedText.BitsFor(text.Sigma));
        long total = PackedText.TotalBytes(text.Length, bits);
        using (var fs = new FileStream(outPath, FileMode.Create, FileAccess.Write))
        {
            var bytes = new byte[total];
            fixed (byte* p = bytes) PackedText.PackInto(text.Symbols, bits, p);
            fs.Write(bytes);
        }
        File.WriteAllText(outPath + ".meta", string.Create(CultureInfo.InvariantCulture, $"n={text.Length}\nbits={bits}\nsigma={text.Sigma}\nletters={text.Letters}\nname={text.Name}\n"));
        r.Row("E3", spec, "packed_file", total, "bytes", $"n={text.Length}, sigma={text.Sigma}, b={bits}, {total / (double)text.Length:F3} B/char");
        return 0;
    }

    internal static int MakeRanks(Options o, Report r)
    {
        string textPath = o.Require("text");
        string outPath = o.Require("out");
        (long n, int bits) = ReadMeta(textPath);
        int requested = o.Int("period", 1024);
        DifferenceCover cover = o.Flag("exact") ? DifferenceCover.Create(requested) : DifferenceCover.Create(DifferenceCover.NearestSingerPeriod(requested));
        using var mapped = new MappedFile(textPath, writable: false);
        using var packed = new PackedText(mapped.Pointer, n, bits, mapped, ownsOwner: false);
        long t0 = Stopwatch.GetTimestamp();
        using DcsRankTable table = DcsBuilder.Build(packed, cover, new WorkCounters(), o.Long("dcs-base", 16_000_000), out DcsStats st);
        double seconds = (Stopwatch.GetTimestamp() - t0) / (double)Stopwatch.Frequency;
        long bytes = table.Slots * 4;
        using (var fs = new FileStream(outPath, FileMode.Create, FileAccess.Write))
        {
            const long chunk = 1 << 24;
            for (long off = 0; off < bytes; off += chunk) fs.Write(new ReadOnlySpan<byte>(table.Owner.Pointer + off, (int)Math.Min(chunk, bytes - off)));
        }
        File.WriteAllText(outPath + ".meta", string.Create(CultureInfo.InvariantCulture, $"period={cover.Period}\nn={n}\nkind={cover.Kind}\nresidues={cover.Residues.Length}\n"));
        r.Row("E3", "ranks", "mkranks_s", seconds, "s", $"period {cover.Period} ({cover.Kind}), |D|={cover.Residues.Length}, sample={st.SampleSize}, unique names={st.NamesUnique}, levels={st.Levels}");
        r.Row("E3", "ranks", "rank_file", bytes, "bytes", $"{bytes / (double)n:F3} B/char");
        r.Row("E3", "ranks", "vmhwm", E1.ReadStatusKb("VmHWM:") * 1024.0 / 1e6, "MB", "peak resident set of the rank construction (in RAM here, external in the design)");
        return 0;
    }

    internal static int HashSa(Options o, Report r)
    {
        string path = o.Require("file");
        string hash = SaHash.ComputeFile(path);
        r.Row("E3", o.Str("case", "sa"), "sa_hash", hash, "", $"SHA-256 of entries 1..n of {Path.GetFileName(path)} (libsais convention)");
        return 0;
    }

    internal static int OracleHash(Options o, Report r)
    {
        string spec = o.Require("family");
        Gen.SymbolText text = FamilyLoader.Load(spec, o.Str("data", string.Empty));
        int[] sa = Oracle.SaisKernelSa(text, out double seconds);
        string hash;
        fixed (int* p = sa) hash = SaHash.Compute((uint*)p + 1, sa.Length - 1L);
        if (o.Get("out") is { } saOut)
        {
            using var fs = new FileStream(saOut, FileMode.Create, FileAccess.Write, FileShare.None, 1 << 20);
            fs.Write(System.Runtime.InteropServices.MemoryMarshal.AsBytes(sa.AsSpan()));
        }
        r.Row("E3", spec, "oracle_sa_hash", hash, "", $"SaisKernel in RAM on the same text, {seconds:F2} s");
        return 0;
    }

    internal static int Build(Options o, Report r)
    {
        string textPath = o.Require("text");
        string saPath = o.Require("sa");
        string scratch = o.Require("scratch");
        string caseName = o.Str("case", "e3");
        (long n, int bits) = ReadMeta(textPath);
        long k = o.Long("k", Math.Max(2, n / o.Long("k-div", 3000)));
        int threads = o.Int("threads", 4);
        long budgetCap = o.Long("budget-mb", 160) << 20;
        Directory.CreateDirectory(scratch);
        foreach (string f in Directory.GetFiles(scratch)) File.Delete(f);

        using var mapped = new MappedFile(textPath, writable: false);
        using var packed = new PackedText(mapped.Pointer, n, bits, mapped, ownsOwner: false);
        DcsRankTable? ranks = null;
        MappedFile? rankMap = null;
        string rankKind = "none";
        if (o.Get("ranks") is { } rankPath)
        {
            string[] meta = File.ReadAllLines(rankPath + ".meta");
            int period = int.Parse(meta.First(l => l.StartsWith("period=", StringComparison.Ordinal))[7..], CultureInfo.InvariantCulture);
            DifferenceCover cover = DifferenceCover.Create(period);
            rankMap = new MappedFile(rankPath, writable: false);
            ranks = new DcsRankTable(cover, rankMap, n);
            rankKind = $"X={period} {cover.Kind} |D|={cover.Residues.Length}";
        }

        var budget = new ArenaBudget(budgetCap);
        var options = new BuildOptions { LeafCapacity = k, FanOut = o.Int("fan", 64), Threads = threads, Seed = 1, Oversample = o.Int("oversample", 16), ScratchDir = scratch, SpillBytes = o.Int("spill-kb", 64) * 1024, Budget = budget };
        using var sink = new MappedSaSink(saPath, n + 1);

        IoCounters.Reset();
        ScratchTracker.Reset();
        ProcIoSnapshot io0 = ProcIoSnapshot.Read();
        (long minflt0, long majflt0) = Faults();
        using var rss = new RssSampler(10, keepSeries: o.Flag("series"));
        using var mem = new MemSampler(o.Get("cgroup-usage"));
        mem.Start();
        rss.Start();
        long t0 = Stopwatch.GetTimestamp();
        BuildStats st = EmBuilder.Run(packed, options, sink, ranks);
        double build = (Stopwatch.GetTimestamp() - t0) / (double)Stopwatch.Frequency;
        RssStats rs = rss.Stop();
        MemStats ms = mem.Stop();
        ProcIoSnapshot io = ProcIoSnapshot.Read().Since(io0);
        (long minflt1, long majflt1) = Faults();
        long scratchPeak = ScratchTracker.Peak;
        long scratchLeft = Directory.GetFiles(scratch).Length;
        sink.Dispose();
        long f0 = Stopwatch.GetTimestamp();
        using (var fs = new FileStream(saPath, FileMode.Open, FileAccess.ReadWrite)) fs.Flush(true);
        double fsync = (Stopwatch.GetTimestamp() - f0) / (double)Stopwatch.Frequency;

        r.Row("E3", caseName, "n", n, "symbols", $"b={bits} bits, K={k} (n/{n / k}), fan-out {options.FanOut}, threads={threads}, ranks: {rankKind}");
        r.Row("E3", caseName, "sa_s", build, "s", "partition + sort, SA written through mapped windows (page cache), without fsync");
        r.Row("E3", caseName, "fsync_s", fsync, "s", "flush of the SA file after the build");
        r.Row("E3", caseName, "ns_per_symbol", build * 1e9 / n, "ns", "sa_s / n");
        r.Row("E3", caseName, "rssanon_peak", rs.MaxRssAnonMb, "MB", $"max RssAnon of {rs.Samples} samples at 10 ms; gate <= 256 MB incl. the .NET runtime and thread stacks (max threads {rs.MaxThreads})");
        r.Row("E3", caseName, "rssfile_peak", rs.MaxRssFileMb, "MB", "mapped file pages resident (text, ranks, SA windows): reported, not counted (ADR-004 section 6a)");
        r.Row("E3", caseName, "vmhwm", rs.MaxVmHwmKb * 1024.0 / 1e6, "MB", "VmHWM (anonymous + file-backed resident), reported only");
        GCMemoryInfo gc = GC.GetGCMemoryInfo();
        r.Row("E3", caseName, "gc_heap", gc.HeapSizeBytes / 1e6, "MB", $"managed heap at the end; GC total available memory {gc.TotalAvailableMemoryBytes / 1e6:F0} MB (DOTNET_GCHeapHardLimit = {Environment.GetEnvironmentVariable("DOTNET_GCHeapHardLimit") ?? "not set"}); collections gen0/1/2 = {GC.CollectionCount(0)}/{GC.CollectionCount(1)}/{GC.CollectionCount(2)}");
        r.Row("E3", caseName, "arena_peak", budget.PeakBytes / 1e6, "MB", $"peak native buffers charged to the ArenaBudget (cap {budgetCap >> 20} MiB)");
        r.Row("E3", caseName, "io_explicit_read", IoCounters.ReadBytes, "bytes", $"{IoCounters.ReadBytes / (double)n:F2} B/char; scratch position files and samples ({IoCounters.ReadCalls} calls)");
        r.Row("E3", caseName, "io_explicit_write", IoCounters.WriteBytes, "bytes", $"{IoCounters.WriteBytes / (double)n:F2} B/char ({IoCounters.WriteCalls} calls)");
        r.Row("E3", caseName, "io_mapped_write", (n + 1) * 4, "bytes", $"{(n + 1) * 4 / (double)n:F2} B/char: the SA, once, through the mapped windows (largest window {sink.MaxWindowBytes} bytes)");
        r.Row("E3", caseName, "io_total_counted", IoCounters.ReadBytes + IoCounters.WriteBytes + ((n + 1) * 4), "bytes", $"{(IoCounters.ReadBytes + IoCounters.WriteBytes + ((n + 1) * 4)) / (double)n:F2} B/char explicit I/O plus SA (text and rank reads through mappings come on top: see block-layer read_bytes)");
        r.Row("E3", caseName, "blocklayer_read", io.ReadBytes, "bytes", $"/proc/self/io read_bytes during the build (page-cache misses: cold text/rank pages and re-reads under a limit); {io.ReadBytes / (double)n:F2} B/char");
        r.Row("E3", caseName, "blocklayer_write", io.WriteBytes, "bytes", "/proc/self/io write_bytes during the build (writeback caused by the process)");
        r.Row("E3", caseName, "page_faults", majflt1 - majflt0, "major", $"minor {minflt1 - minflt0}");
        r.Row("E3", caseName, "scratch_peak", scratchPeak, "bytes", $"{scratchPeak / (double)n:F2} B/char; files left in the scratch directory: {scratchLeft}");
        r.Row("E3", caseName, "dirty_peak", ms.MaxDirtyKb * 1024.0 / 1e6, "MB", $"max /proc/meminfo Dirty (system-wide); max Cached {ms.MaxCachedKb * 1024.0 / 1e6:F0} MB; cgroup usage max {ms.MaxCgroupUsageBytes / 1e6:F0} MB");
        if (ms.CgSamples > 0)
        {
            r.Row("E3", caseName, "cgroup_rss", ms.MaxCgRssBytes / 1e6, "MB", $"memory.stat rss of the sub-group: max, mean {ms.AvgCgRssBytes / 1e6:F1} MB over {ms.CgSamples} samples at 50 ms");
            r.Row("E3", caseName, "cgroup_mapped_file", ms.AvgCgMappedFileBytes / 1e6, "MB", "mean mapped file pages charged to the sub-group (text + rank + SA windows)");
            r.Row("E3", caseName, "cgroup_cache", ms.AvgCgCacheBytes / 1e6, "MB", "mean page cache charged to the sub-group (incl. scratch files and the SA file)");
        }
        r.Row("E3", caseName, "buckets", st.Buckets, "n", $"leaves={st.Leaves}, partitions={st.Partitions}, max depth={st.MaxDepth}, largest leaf={st.MaxLeaf}, retries={st.Retries}; tie share {100.0 * st.Work.TiePositions / Math.Max(1, st.Work.LeafPositions):F4} %, rank lookups {st.Work.RankLookups}");
        r.Row("E3", caseName, "phase_cpu_s", st.PartitionSeconds + st.LeafSeconds, "s", $"summed over workers: partition {st.PartitionSeconds:F2}, leaf sort {st.LeafSeconds:F2}");
        rankMap?.Dispose();
        return rs.MaxRssAnonMb <= o.Double("rss-cap-mb", 256) && scratchLeft == 0 ? 0 : 1;
    }

    private static (long Minor, long Major) Faults()
    {
        string[] f = File.ReadAllText("/proc/self/stat").Split(')')[^1].Split(' ', StringSplitOptions.RemoveEmptyEntries);
        // after "(comm)": state is f[0], minflt is field 10 (f[7]), majflt is field 12 (f[9])
        return (long.Parse(f[7], CultureInfo.InvariantCulture), long.Parse(f[9], CultureInfo.InvariantCulture));
    }

    internal static (long N, int Bits) ReadMeta(string textPath)
    {
        string[] meta = File.ReadAllLines(textPath + ".meta");
        long n = long.Parse(meta.First(l => l.StartsWith("n=", StringComparison.Ordinal))[2..], CultureInfo.InvariantCulture);
        int bits = int.Parse(meta.First(l => l.StartsWith("bits=", StringComparison.Ordinal))[5..], CultureInfo.InvariantCulture);
        return (n, bits);
    }
}
