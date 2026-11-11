using System.Diagnostics;
using System.Globalization;
using SuffixTree.EmProto.Build;
using SuffixTree.EmProto.Dcs;
using SuffixTree.EmProto.Harness;
using SuffixTree.EmProto.Text;

namespace SuffixTree.EmProto.Experiments;

/// <summary>
/// E3b of ADR-004 (owner decision 2026-11-11, section 6c): the external-memory PH1, i.e. the difference-cover rank file built with bounded RAM where E3 built it in RAM and mapped it. The command
/// <c>e3b</c> builds the rank file of a packed text (<c>pack</c> makes one) and reports wall time, RssAnon, I/O counters, scratch, the recursion and the share of repeated names; <c>--verify</c> compares the
/// result byte for byte with a rank file of the in-RAM builder (<c>mkranks</c>). The scale-down emulation of gate G3 is <c>--scale-n 3000000000</c>: every size that decides the number of buckets,
/// levels or windows is its 3x10^9 reference value (the <c>--*-ref</c> options) multiplied by n / 3x10^9, so that the build has the structure that 3x10^9 symbols would have under the full-size parameters.
/// </summary>
internal static unsafe class E3b
{
    /// <summary>Reference values for 3x10^9 symbols (ADR-004 section 4.2 and the choices of the E3b results): K_s, N_base, R, K of the external sort.</summary>
    internal const long ReferenceN = 3_000_000_000;
    internal const long KsRef = 1 << 20;
    internal const long BaseRef = 16 << 20;
    internal const long ChunkRef = 4 << 20;
    internal const long SortKRef = 1 << 20;

    internal static int Build(Options o, Report r)
    {
        string textPath = o.Require("text");
        string outPath = o.Require("out");
        string scratch = o.Require("scratch");
        string caseName = o.Str("case", "e3b");
        (long n, int bits) = E3.ReadMeta(textPath);
        int requested = o.Int("period", 1024);
        DifferenceCover cover = o.Flag("exact") ? DifferenceCover.Create(requested) : DifferenceCover.Create(DifferenceCover.NearestSingerPeriod(requested));
        int threads = o.Int("threads", 4);
        long budgetCap = o.Long("budget-mb", 160) << 20;
        bool scaled = o.Get("scale-n") is not null;
        double f = scaled ? n / (double)o.Long("scale-n", ReferenceN) : 1.0;
        long Scale(string name, long reference, long floor = 3)
        {
            if (o.Get(name) is { } v) return long.Parse(v, CultureInfo.InvariantCulture);
            if (!scaled) return reference;
            return Math.Max(floor, (long)Math.Round(reference * f));
        }
        ExtDcsOptions options = new()
        {
            ScratchDir = scratch,
            BucketCapacity = Math.Max(3, Scale("k-s", KsRef)),
            BaseCase = Math.Max(4, Scale("nbase", BaseRef)),
            ChunkSlots = Math.Max(16, Scale("chunk", ChunkRef)),
            SortLeafCapacity = Math.Max(3, Scale("sort-k", SortKRef)),
            Threads = threads,
            Seed = (ulong)o.Long("seed", 1),
            FanOut = o.Int("fan", 1024),
            ReservoirSize = o.Int("reservoir", 16384),
            RunSpillBytes = o.Get("run-spill-kb") is { } spill ? int.Parse(spill, CultureInfo.InvariantCulture) * 1024 : (int)Scale("run-spill", 64 * 1024, 0),
            IoBufferBytes = o.Int("io-kb", 64) * 1024,
            MergeBufferBudget = o.Long("merge-budget-mb", 16) << 20,
            WindowBytes = o.Int("window-kb", 1024) * 1024,
            Budget = new ArenaBudget(budgetCap),
        };
        Directory.CreateDirectory(scratch);
        foreach (string file in Directory.GetFileSystemEntries(scratch)) File.Delete(file);

        long sRef = cover.SampleCount(ReferenceN);
        long sHere = cover.SampleCount(n);
        r.Row("E3b", caseName, "n", n, "symbols", $"b={bits} bits, X={cover.Period} {cover.Kind} |D|={cover.Residues.Length}, sample {sHere} ({sHere * 100.0 / n:F2} %), threads={threads}");
        r.Row("E3b", caseName, "parameters", string.Create(CultureInfo.InvariantCulture, $"K_s={options.BucketCapacity} N_base={options.BaseCase} R={options.ChunkSlots} K_sort={options.SortLeafCapacity}"), "", scaled
            ? $"scaled by n/3e9 = {f:E3} from K_s={KsRef} N_base={BaseRef} R={ChunkRef} K_sort={SortKRef} run-spill=64 KiB (run spill {options.RunSpillBytes} bytes: every bucket run of the 3e9 build is a scratch file); the 3e9 build has sample {sRef}, about {2 * sRef / KsRef} sample buckets (2 x sample / K_s), this one about {2 * sHere / options.BucketCapacity}"
            : "unscaled (explicit values or the reference values)");

        IoCounters.Reset();
        MappedIo.Reset();
        ScratchTracker.Reset();
        ProcIoSnapshot io0 = ProcIoSnapshot.Read();
        using var rss = new RssSampler(10, keepSeries: o.Flag("series"));
        using var mem = new MemSampler(o.Get("cgroup-usage"));
        using var mapped = new MappedFile(textPath, writable: false);
        using var packed = new PackedText(mapped.Pointer, n, bits, mapped, ownsOwner: false);
        mem.Start();
        rss.Start();
        long t0 = Stopwatch.GetTimestamp();
        ExtDcsStats st = ExternalDcs.BuildRankFile(packed, cover, outPath, options);
        double build = (Stopwatch.GetTimestamp() - t0) / (double)Stopwatch.Frequency;
        RssStats rs = rss.Stop();
        MemStats ms = mem.Stop();
        ProcIoSnapshot io = ProcIoSnapshot.Read().Since(io0);
        long scratchPeak = ScratchTracker.Peak;
        long scratchLive = ScratchTracker.Live;
        string[] left = Directory.GetFileSystemEntries(scratch);
        long f0 = Stopwatch.GetTimestamp();
        using (var fs = new FileStream(outPath, FileMode.Open, FileAccess.ReadWrite)) fs.Flush(true);
        double fsync = (Stopwatch.GetTimestamp() - f0) / (double)Stopwatch.Frequency;
        long rankBytes = new FileInfo(outPath).Length;

        r.Row("E3b", caseName, "ph1_s", build, "s", "external-memory PH1: sample buckets, merge, levels below, rank file through mapped windows (without fsync)");
        r.Row("E3b", caseName, "fsync_s", fsync, "s", "flush of the rank file after the build");
        r.Row("E3b", caseName, "ns_per_symbol", build * 1e9 / n, "ns", "ph1_s / n");
        r.Row("E3b", caseName, "rssanon_peak", rs.MaxRssAnonMb, "MB", $"max RssAnon of {rs.Samples} samples at 10 ms; gate G2 <= 256 MB incl. the .NET runtime and thread stacks (max threads {rs.MaxThreads})");
        r.Row("E3b", caseName, "rssfile_peak", rs.MaxRssFileMb, "MB", "mapped file pages resident (text, level files, rank windows): reported, not counted (ADR-004 section 6a)");
        r.Row("E3b", caseName, "vmhwm", rs.MaxVmHwmKb * 1024.0 / 1e6, "MB", "VmHWM (anonymous + file-backed resident), reported only");
        GCMemoryInfo gc = GC.GetGCMemoryInfo();
        r.Row("E3b", caseName, "gc_heap", gc.HeapSizeBytes / 1e6, "MB", $"managed heap at the end; GC total available memory {gc.TotalAvailableMemoryBytes / 1e6:F0} MB (DOTNET_GCHeapHardLimit = {Environment.GetEnvironmentVariable("DOTNET_GCHeapHardLimit") ?? "not set"}); collections gen0/1/2 = {GC.CollectionCount(0)}/{GC.CollectionCount(1)}/{GC.CollectionCount(2)}");
        r.Row("E3b", caseName, "arena_peak", options.Budget!.PeakBytes / 1e6, "MB", $"peak native buffers charged to the ArenaBudget (cap {budgetCap >> 20} MiB), live at the end {options.Budget.LiveBytes}");
        r.Row("E3b", caseName, "io_explicit_read", IoCounters.ReadBytes, "bytes", $"{IoCounters.ReadBytes / (double)n:F3} B/char; run, level, suffix-array and tuple files ({IoCounters.ReadCalls} calls)");
        r.Row("E3b", caseName, "io_explicit_write", IoCounters.WriteBytes, "bytes", $"{IoCounters.WriteBytes / (double)n:F3} B/char ({IoCounters.WriteCalls} calls)");
        r.Row("E3b", caseName, "io_mapped_write", MappedIo.WriteBytes, "bytes", $"{MappedIo.WriteBytes / (double)n:F3} B/char through mapped windows of at most {options.WindowBytes} bytes ({MappedIo.Windows} windows, largest {MappedIo.MaxWindowBytes}); the suffix arrays of the external levels (EmBuilder) are not in this number");
        long counted = IoCounters.ReadBytes + IoCounters.WriteBytes + MappedIo.WriteBytes;
        r.Row("E3b", caseName, "io_total_counted", counted, "bytes", $"{counted / (double)n:F3} B/char explicit I/O plus mapped writes (text reads and the SA writes of the external levels come on top: see block-layer read_bytes)");
        r.Row("E3b", caseName, "blocklayer_read", io.ReadBytes, "bytes", $"/proc/self/io read_bytes during the build (page-cache misses); {io.ReadBytes / (double)n:F3} B/char");
        r.Row("E3b", caseName, "blocklayer_write", io.WriteBytes, "bytes", "/proc/self/io write_bytes during the build (writeback caused by the process)");
        r.Row("E3b", caseName, "scratch_peak", scratchPeak, "bytes", $"{scratchPeak / (double)n:F3} B/char peak of all scratch files (runs, reduced strings, level rank files, suffix arrays, tuple files, position buckets), excluding the output; live at the end {scratchLive}; entries left in the scratch directory: {left.Length}");
        r.Row("E3b", caseName, "rank_file", rankBytes, "bytes", $"{rankBytes / (double)n:F3} B/char (the output)");
        r.Row("E3b", caseName, "dirty_peak", ms.MaxDirtyKb * 1024.0 / 1e6, "MB", $"max /proc/meminfo Dirty (system-wide); max Cached {ms.MaxCachedKb * 1024.0 / 1e6:F0} MB");
        r.Row("E3b", caseName, "levels", st.MaxLevel, "n", $"deepest level reached (0 = the names of the text are all distinct); levels {st.Levels.Count}");
        ExtDcsLevelStats l0 = st.Levels[0];
        r.Row("E3b", caseName, "nonunique_share", 100.0 * l0.NonUniqueRecords / Math.Max(1, l0.SampleSize), "%", $"sample positions at level 0 whose name is shared: {l0.NonUniqueRecords} of {l0.SampleSize}; distinct names {l0.DistinctNames}");
        foreach (ExtDcsLevelStats l in st.Levels)
        {
            r.Row("E3b", caseName, $"level{l.Level}", l.TextLength, "symbols", string.Create(CultureInfo.InvariantCulture, $"sample {l.SampleSize}, distinct names {l.DistinctNames}, shared {l.NonUniqueRecords} ({100.0 * l.NonUniqueRecords / Math.Max(1, l.SampleSize):F2} %), unique={l.NamesUnique}, base case={l.BaseCase}, leaf buckets {l.Leaves}, split scans {l.SplitScans}, scans {l.Scans}, file runs {l.FileRuns}, naming {l.NamingSeconds:F2} s, merge {l.MergeSeconds:F2} s, below {l.BelowSeconds:F2} s"));
        }
        bool ok = rs.MaxRssAnonMb <= o.Double("rss-cap-mb", 256) && left.Length == 0 && scratchLive == 0 && options.Budget.LiveBytes == 0;
        if (o.Get("verify") is { } oracle)
        {
            (bool equal, long diffAt) = CompareFiles(outPath, oracle);
            r.Row("E3b", caseName, "rank_equals_oracle", equal ? "yes" : "NO", "", equal ? $"byte-identical to {Path.GetFileName(oracle)} ({rankBytes} bytes)" : $"first difference at byte {diffAt} of {rankBytes} (oracle {new FileInfo(oracle).Length} bytes)");
            ok &= equal;
        }
        return ok ? 0 : 1;
    }

    /// <summary>Byte-for-byte comparison of two files in 1 MB pieces: equality and the offset of the first difference (or of the shorter end).</summary>
    internal static (bool Equal, long FirstDifference) CompareFiles(string a, string b)
    {
        using var fa = new FileStream(a, FileMode.Open, FileAccess.Read, FileShare.Read, 1 << 20, FileOptions.SequentialScan);
        using var fb = new FileStream(b, FileMode.Open, FileAccess.Read, FileShare.Read, 1 << 20, FileOptions.SequentialScan);
        var x = new byte[1 << 20];
        var y = new byte[1 << 20];
        long at = 0;
        while (true)
        {
            int gx = ReadFull(fa, x);
            int gy = ReadFull(fb, y);
            int common = Math.Min(gx, gy);
            int diff = x.AsSpan(0, common).CommonPrefixLength(y.AsSpan(0, common));
            if (diff < common) return (false, at + diff);
            if (gx != gy) return (false, at + common);
            if (gx == 0) return (true, -1);
            at += gx;
        }
    }

    private static int ReadFull(Stream s, byte[] buffer)
    {
        int total = 0;
        while (total < buffer.Length)
        {
            int got = s.Read(buffer, total, buffer.Length - total);
            if (got <= 0) break;
            total += got;
        }
        return total;
    }
}
