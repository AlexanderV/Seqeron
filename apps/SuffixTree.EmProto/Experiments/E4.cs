using System.Diagnostics;
using System.Globalization;
using System.Runtime.InteropServices;
using SuffixTree.Esa;
using SuffixTree.EmProto.Build;
using SuffixTree.EmProto.Dcs;
using SuffixTree.EmProto.Harness;
using SuffixTree.EmProto.Ph4;
using SuffixTree.EmProto.Text;

namespace SuffixTree.EmProto.Experiments;

/// <summary>
/// E4 of ADR-004 section 5 (owner decision 2026-11-11, section 6c): PH4, the LCP array from a given suffix-array file with bounded RAM. The command <c>e4</c> builds the compact LCP (a base file and an overflow
/// file, the layout of <see cref="ExternalLcp"/>) from a packed text and the SA file of <c>e3</c> and reports wall time per phase B1 / B2 / B3, RssAnon, I/O counters, scratch, the width and the escape
/// share; <c>lcp-oracle</c> writes the same two files from the in-RAM <see cref="LcpKernel"/> (the oracle of gate G1 up to 100 Mbp, <c>--verify</c> compares byte for byte), <c>check-lcp</c> checks sampled ranks
/// against a direct comparison of the text and the structure of the overflow list (the check above the size that fits in RAM). The scale-down emulation of gate G2 is <c>--scale-n 3000000000</c>: R is its
/// reference value 8,000,000 (ADR-004 section 4.2: P = n / R = 375 chunks at 3x10^9 symbols) multiplied by n / 3x10^9 (rounded up), so that the chunk count, and with <c>--rounds</c> the round count, are those of 3x10^9.
/// </summary>
internal static unsafe class E4
{
    /// <summary>The text length that the scale-down emulation reproduces.</summary>
    internal const long ReferenceN = 3_000_000_000;

    /// <summary>R at <see cref="ReferenceN"/>: positions per chunk (ADR-004 section 4.2), P = 375.</summary>
    internal const long ChunkRef = 8_000_000;

    internal static int Build(Options o, Report r)
    {
        string textPath = o.Require("text");
        string saPath = o.Require("sa");
        string outPrefix = o.Require("out");
        string scratch = o.Require("scratch");
        string caseName = o.Str("case", "e4");
        LcpMode mode = o.Str("mode", "distributed") switch
        {
            "distributed" => LcpMode.Distributed,
            "bounded" => LcpMode.Bounded,
            var other => throw new ArgumentException($"unknown --mode '{other}' (distributed or bounded)"),
        };
        (long n, int bits) = E3.ReadMeta(textPath);
        int threads = o.Int("threads", 4);
        long budgetCap = o.Long("budget-mb", 160) << 20;
        bool scaled = o.Get("scale-n") is not null;
        double f = scaled ? n / (double)o.Long("scale-n", ReferenceN) : 1.0;
        long chunk = scaled ? Math.Max(3, (long)Math.Ceiling(ChunkRef * f)) : ChunkRef;
        if (o.Get("chunk") is { } explicitChunk) chunk = long.Parse(explicitChunk, CultureInfo.InvariantCulture);
        ExtLcpOptions options = new()
        {
            ScratchDir = scratch,
            ChunkPositions = chunk,
            Threads = threads,
            ReadBufferBytes = o.Int("io-kb", 1024) * 1024,
            TupleBufferBytes = o.Int("tuple-kb", 32) * 1024,
            WindowBytes = o.Int("window-kb", 1024) * 1024,
            Mode = mode,
            TouchBuffers = o.Flag("touch-buffers"),
            Budget = new ArenaBudget(budgetCap),
        };
        // rounds of the bounded mode: --rounds r gives ceil(P / r) chunks per round, so that the round count is r whatever n is
        long chunks = (n + chunk - 1) / chunk;
        if (mode == LcpMode.Bounded) options = options with { RoundChunks = o.Get("rounds") is { } rounds ? (chunks + long.Parse(rounds, CultureInfo.InvariantCulture) - 1) / long.Parse(rounds, CultureInfo.InvariantCulture) : o.Long("round-chunks", 0) };
        // the number of parallel B2 tasks is what the arena cap allows at the REFERENCE geometry (R = 8,000,000 at 3x10^9), so that a scaled run has the task count of the real one
        ExtLcpOptions reference = options with { ChunkPositions = ChunkRef };
        long perTask = ExternalLcp.WorkerBytes(ReferenceN, reference);
        int workers = o.Int("b2-workers", (int)Math.Clamp(budgetCap / Math.Max(1, perTask), 1, threads));
        options = options with { Workers = workers };

        Directory.CreateDirectory(scratch);
        foreach (string file in Directory.GetFileSystemEntries(scratch)) File.Delete(file);
        string basePath = outPrefix + ".lcp";
        string overflowPath = outPrefix + ".ovf";

        r.Row("E4", caseName, "n", n, "symbols", $"b={bits} bits, mode={mode}, threads={threads}, B2 workers={workers}, touch buffers={options.TouchBuffers}");
        long roundCount = 1;
        if (mode == LcpMode.Bounded && options.RoundChunks > 0) roundCount = (chunks + options.RoundChunks - 1) / options.RoundChunks;
        r.Row("E4", caseName, "parameters", string.Create(CultureInfo.InvariantCulture, $"R={chunk} P={chunks} rounds={roundCount} tuple={options.TupleBufferBytes} read={options.ReadBufferBytes} window={options.WindowBytes}"), "", scaled
            ? $"R scaled by n/3e9 = {f:E3} from R={ChunkRef} (P = 375 at 3e9); B2 workers from the reference geometry: one task holds {perTask / 1048576.0:F1} MiB at 3e9, arena cap {budgetCap >> 20} MiB"
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
        ExtLcpStats st = ExternalLcp.Build(saPath, packed, basePath, overflowPath, options);
        double build = (Stopwatch.GetTimestamp() - t0) / (double)Stopwatch.Frequency;
        RssStats rs = rss.Stop();
        MemStats ms = mem.Stop();
        ProcIoSnapshot io = ProcIoSnapshot.Read().Since(io0);
        long scratchPeak = ScratchTracker.Peak;
        long scratchLive = ScratchTracker.Live;
        string[] left = Directory.GetFileSystemEntries(scratch);
        long f0 = Stopwatch.GetTimestamp();
        foreach (string p in new[] { basePath, overflowPath })
        {
            using var fs = new FileStream(p, FileMode.Open, FileAccess.ReadWrite);
            fs.Flush(true);
        }
        double fsync = (Stopwatch.GetTimestamp() - f0) / (double)Stopwatch.Frequency;
        long baseBytes = new FileInfo(basePath).Length;
        long overflowBytes = new FileInfo(overflowPath).Length;

        r.Row("E4", caseName, "ph4_s", build, "s", "PH4 total: B1 + B2 + B3 (without fsync)");
        r.Row("E4", caseName, "b1_s", st.B1Seconds, "s", mode == LcpMode.Bounded ? "B1 of all rounds: suffix array scanned once per round, Phi pairs of the round distributed" : "B1: suffix array scanned once, Phi tuples distributed by position / R");
        r.Row("E4", caseName, "b2_s", st.B2Seconds, "s", mode == LcpMode.Bounded ? "B2 of all rounds: scatter, PLCP loop, PLCP chunk to the text-order file" : "B2: scatter, PLCP loop, rank tuples distributed by rank / R");
        r.Row("E4", caseName, "b3_s", st.B3Seconds, "s", mode == LcpMode.Bounded ? "B3: gather PLCP[SA[i]] in rank order, output through mapped windows" : "B3: rank chunks scattered and written through mapped windows");
        r.Row("E4", caseName, "fsync_s", fsync, "s", "flush of the two output files after the build");
        r.Row("E4", caseName, "ns_per_symbol", build * 1e9 / n, "ns", "ph4_s / n");
        r.Row("E4", caseName, "throughput", n / build / 1e6, "Msym/s", "n / ph4_s; the gate of ADR-004 row E4 is 5 MB/s = 5 M symbols/s on 4 threads");
        r.Row("E4", caseName, "rssanon_peak", rs.MaxRssAnonMb, "MB", $"max RssAnon of {rs.Samples} samples at 10 ms; gate G2 <= 256 MB incl. the .NET runtime and thread stacks (max threads {rs.MaxThreads})");
        r.Row("E4", caseName, "rssfile_peak", rs.MaxRssFileMb, "MB", "mapped file pages resident (text, output windows, PLCP file in the bounded mode): reported, not counted (ADR-004 section 6a)");
        r.Row("E4", caseName, "vmhwm", rs.MaxVmHwmKb * 1024.0 / 1e6, "MB", "VmHWM (anonymous + file-backed resident), reported only");
        GCMemoryInfo gc = GC.GetGCMemoryInfo();
        r.Row("E4", caseName, "gc_heap", gc.HeapSizeBytes / 1e6, "MB", $"managed heap at the end; GC total available memory {gc.TotalAvailableMemoryBytes / 1e6:F0} MB (DOTNET_GCHeapHardLimit = {Environment.GetEnvironmentVariable("DOTNET_GCHeapHardLimit") ?? "not set"}); collections gen0/1/2 = {GC.CollectionCount(0)}/{GC.CollectionCount(1)}/{GC.CollectionCount(2)}");
        r.Row("E4", caseName, "arena_peak", options.Budget!.PeakBytes / 1e6, "MB", $"peak native buffers charged to the ArenaBudget (cap {budgetCap >> 20} MiB), live at the end {options.Budget.LiveBytes}");
        r.Row("E4", caseName, "io_explicit_read", IoCounters.ReadBytes, "bytes", $"{IoCounters.ReadBytes / (double)n:F3} B/char; suffix array and tuple files ({IoCounters.ReadCalls} calls)");
        r.Row("E4", caseName, "io_explicit_write", IoCounters.WriteBytes, "bytes", $"{IoCounters.WriteBytes / (double)n:F3} B/char; tuple files and the PLCP file ({IoCounters.WriteCalls} calls)");
        r.Row("E4", caseName, "io_mapped_write", MappedIo.WriteBytes, "bytes", $"{MappedIo.WriteBytes / (double)n:F3} B/char: the LCP section through mapped windows of at most {options.WindowBytes} bytes ({MappedIo.Windows} windows, largest {MappedIo.MaxWindowBytes})");
        long counted = IoCounters.ReadBytes + IoCounters.WriteBytes + MappedIo.WriteBytes;
        long gather = mode == LcpMode.Bounded ? 4 * n : 0;
        r.Row("E4", caseName, "io_total_counted", counted, "bytes", $"{counted / (double)n:F3} B/char explicit I/O plus mapped writes; the bounded mode's gather reads {gather / (double)n:F1} B/char more through a read-only mapping ({(counted + gather) / (double)n:F3} B/char with it); the text reads come on top (see block-layer read_bytes)");
        r.Row("E4", caseName, "blocklayer_read", io.ReadBytes, "bytes", $"/proc/self/io read_bytes during the build (page-cache misses); {io.ReadBytes / (double)n:F3} B/char");
        r.Row("E4", caseName, "blocklayer_write", io.WriteBytes, "bytes", "/proc/self/io write_bytes during the build (writeback caused by the process)");
        r.Row("E4", caseName, "scratch_peak", scratchPeak, "bytes", $"{scratchPeak / (double)n:F3} B/char peak of all scratch files (tuple files, PLCP file), excluding the output; live at the end {scratchLive}; entries left in the scratch directory: {left.Length}");
        r.Row("E4", caseName, "lcp_file", baseBytes + overflowBytes, "bytes", $"{(baseBytes + overflowBytes) / (double)n:F3} B/char (the output): base {baseBytes}, overflow {overflowBytes}");
        r.Row("E4", caseName, "dirty_peak", ms.MaxDirtyKb * 1024.0 / 1e6, "MB", $"max /proc/meminfo Dirty (system-wide); max Cached {ms.MaxCachedKb * 1024.0 / 1e6:F0} MB");
        r.Row("E4", caseName, "width", st.Kind.ToString(), "", $"chosen before any entry was written; escapes >= 255: {st.Escapes8} ({100.0 * st.Escapes8 / n:F4} %), >= 65535: {st.Escapes16} ({100.0 * st.Escapes16 / n:F4} %), overflow entries {st.OverflowCount}");
        r.Row("E4", caseName, "max_plcp", st.MaxPlcp, "symbols", "largest LCP value");
        r.Row("E4", caseName, "compare_steps", st.CompareSteps, "steps", $"vector (or word) compare steps of all PLCP loops, {st.CompareSteps / (double)n:F3} per position; restart cost (symbols compared at the first position of each of the {st.TextChunks} chunks, no carry): {st.RestartSymbols}");
        r.Row("E4", caseName, "chunks", st.TextChunks, "n", $"text chunks {st.TextChunks}, rank chunks {st.RankChunks}, rounds {st.Rounds}, B2 tasks in parallel {st.Workers}");
        bool ok = rs.MaxRssAnonMb <= o.Double("rss-cap-mb", 256) && left.Length == 0 && scratchLive == 0 && options.Budget.LiveBytes == 0;
        if (o.Get("verify") is { } oracle)
        {
            (bool eqBase, long diffBase) = E3b.CompareFiles(basePath, oracle + ".lcp");
            (bool eqOverflow, long diffOverflow) = E3b.CompareFiles(overflowPath, oracle + ".ovf");
            bool equal = eqBase && eqOverflow;
            r.Row("E4", caseName, "lcp_equals_oracle", equal ? "yes" : "NO", "", equal ? $"base ({baseBytes} bytes) and overflow ({overflowBytes} bytes) byte-identical to the LcpKernel files {Path.GetFileName(oracle)}" : $"first difference: base at byte {diffBase}, overflow at byte {diffOverflow}");
            ok &= equal;
        }
        return ok ? 0 : 1;
    }

    /// <summary><c>lcp-oracle</c>: the compact LCP of the in-RAM <see cref="LcpKernel"/> written in the layout of <see cref="ExternalLcp"/> (base file, overflow pairs); needs n below 2^31.</summary>
    internal static int Oracle(Options o, Report r)
    {
        string textPath = o.Require("text");
        string saPath = o.Require("sa");
        string outPrefix = o.Require("out");
        (long n, int bits) = E3.ReadMeta(textPath);
        if (n >= int.MaxValue || bits > 8) throw new ArgumentException("the kernel oracle needs n < 2^31 and at most 8 bits per symbol");
        using var mapped = new MappedFile(textPath, writable: false);
        using var packed = new PackedText(mapped.Pointer, n, bits, mapped, ownsOwner: false);
        var symbols = new byte[n];
        for (long i = 0; i < n; i++) symbols[i] = (byte)packed.Symbol(i);
        var sa = new int[n + 1];
        using (var fs = new FileStream(saPath, FileMode.Open, FileAccess.Read, FileShare.Read, 1 << 20, FileOptions.SequentialScan)) fs.ReadExactly(MemoryMarshal.AsBytes(sa.AsSpan()));
        var prepared = new PreparedText(new char[1 << bits], symbols, null);
        long t0 = Stopwatch.GetTimestamp();
        CompactLcp lcp = LcpKernel.Build(prepared, sa);
        double seconds = (Stopwatch.GetTimestamp() - t0) / (double)Stopwatch.Frequency;
        using (var fs = new FileStream(outPrefix + ".lcp", FileMode.Create, FileAccess.Write, FileShare.None, 1 << 20))
        {
            switch (lcp.Kind)
            {
                case LcpKind.U8: fs.Write(lcp.BaseU8!); break;
                case LcpKind.U16: fs.Write(MemoryMarshal.AsBytes(lcp.BaseU16.AsSpan())); break;
                default: fs.Write(MemoryMarshal.AsBytes(lcp.BaseU32.AsSpan())); break;
            }
        }
        var pairs = new uint[lcp.OverflowCount * 2];
        for (int i = 0; i < lcp.OverflowCount; i++)
        {
            pairs[2 * i] = (uint)lcp.OverflowRanks[i];
            pairs[(2 * i) + 1] = (uint)lcp.OverflowValues[i];
        }
        File.WriteAllBytes(outPrefix + ".ovf", MemoryMarshal.AsBytes(pairs.AsSpan()).ToArray());
        r.Row("E4", o.Str("case", "lcp-oracle"), "kernel_s", seconds, "s", $"LcpKernel (three passes in RAM, one thread) on n={n}; kind {lcp.Kind}, overflow entries {lcp.OverflowCount}");
        return 0;
    }

    /// <summary>
    /// <c>check-lcp</c>: the output of <c>e4</c> against the text itself where no in-RAM oracle fits: <c>--samples</c> random ranks (and the first and last) compared with the longest common prefix of
    /// <c>SA[r-1]</c> and <c>SA[r]</c> computed directly on the packed text, and a scan of the whole structure (every entry at the escape value has exactly one overflow pair, ranks ascending, values at or above the escape).
    /// </summary>
    internal static int CheckSample(Options o, Report r)
    {
        string textPath = o.Require("text");
        string saPath = o.Require("sa");
        string prefix = o.Require("lcp");
        string caseName = o.Str("case", "check-lcp");
        (long n, int bits) = E3.ReadMeta(textPath);
        long baseLength = new FileInfo(prefix + ".lcp").Length;
        int width = (int)(baseLength / (n + 1));
        if (width is not (1 or 2 or 4) || baseLength != (n + 1) * width) throw new ArgumentException($"the base file has {baseLength} bytes, not (n + 1) x 1, 2 or 4");
        uint escape = width switch { 1 => byte.MaxValue, 2 => ushort.MaxValue, _ => uint.MaxValue };
        long overflowCount = new FileInfo(prefix + ".ovf").Length / 8;
        using var mapped = new MappedFile(textPath, writable: false);
        using var packed = new PackedText(mapped.Pointer, n, bits, mapped, ownsOwner: false);
        using var baseFile = new CountedFile(prefix + ".lcp", FileMode.Open, FileAccess.Read);
        using var overflowFile = new CountedFile(prefix + ".ovf", FileMode.Open, FileAccess.Read);
        using var sa = new CountedFile(saPath, FileMode.Open, FileAccess.Read);

        uint BaseAt(long rank)
        {
            Span<byte> b = stackalloc byte[4];
            baseFile.ReadExactAt(b[..width], rank * width);
            return width switch { 1 => b[0], 2 => BitConverter.ToUInt16(b), _ => BitConverter.ToUInt32(b) };
        }

        (uint Rank, uint Value) PairAt(long i)
        {
            Span<byte> b = stackalloc byte[8];
            overflowFile.ReadExactAt(b, i * 8);
            return (BitConverter.ToUInt32(b), BitConverter.ToUInt32(b[4..]));
        }

        // structure: one pass over the base file
        long escapesInBase = 0;
        long problems = 0;
        {
            var block = new byte[1 << 20];
            long at = 0;
            long total = baseLength;
            long k = 0; // next overflow pair
            long previousRank = -1;
            var pairBuffer = new byte[8 * 4096];
            long pairsBuffered = 0, pairsStart = 0;
            while (at < total)
            {
                int take = (int)Math.Min(block.Length, total - at);
                baseFile.ReadExactAt(block.AsSpan(0, take), at);
                for (int i = 0; i < take; i += width)
                {
                    long rank = (at + i) / width;
                    uint v = width switch { 1 => block[i], 2 => BitConverter.ToUInt16(block, i), _ => BitConverter.ToUInt32(block, i) };
                    if (width < 4 && v == escape)
                    {
                        escapesInBase++;
                        if (k >= overflowCount)
                        {
                            problems++; // an escape slot without a pair
                            continue;
                        }
                        if (k >= pairsStart + pairsBuffered)
                        {
                            pairsStart = k;
                            pairsBuffered = Math.Min(4096, overflowCount - k);
                            overflowFile.ReadExactAt(pairBuffer.AsSpan(0, (int)pairsBuffered * 8), k * 8);
                        }
                        uint pr = BitConverter.ToUInt32(pairBuffer, (int)(k - pairsStart) * 8);
                        uint pv = BitConverter.ToUInt32(pairBuffer, ((int)(k - pairsStart) * 8) + 4);
                        if (pr != rank || pr <= previousRank || pv < escape) problems++;
                        previousRank = pr;
                        k++;
                    }
                }
                at += take;
            }
            if (width == 4 ? overflowCount != 0 : escapesInBase != overflowCount) problems++;
        }

        // samples
        long samples = o.Long("samples", 1000);
        var rng = new SplitMix64((ulong)o.Long("seed", 20240011));
        long mismatches = 0, checkedRanks = 0;
        long steps = 0;
        Span<uint> two = stackalloc uint[2];
        for (long s = 0; s < samples + 2; s++)
        {
            if (n == 0) break;
            long rank = 1 + (long)(rng.Next() % (ulong)n);
            if (s == 0) rank = 1;
            if (s == 1) rank = n;
            sa.ReadExactAt(MemoryMarshal.AsBytes(two), (rank - 1) * 4);
            long expected = 0;
            if (rank > 1)
            {
                long a = two[0], b = two[1];
                expected = packed.Lce(a, b, n - Math.Max(a, b), ref steps);
            }
            long got = BaseAt(rank);
            if (width < 4 && got == escape)
            {
                long lo = 0, hi = overflowCount - 1;
                got = -1;
                while (lo <= hi)
                {
                    long mid = (lo + hi) / 2;
                    (uint pr, uint pv) = PairAt(mid);
                    if (pr == rank) { got = pv; break; }
                    if (pr < rank) lo = mid + 1;
                    else hi = mid - 1;
                }
            }
            checkedRanks++;
            if (got != expected) mismatches++;
        }
        r.Row("E4", caseName, "structure_problems", problems, "n", $"width {width}, escapes in the base array {escapesInBase}, overflow pairs {overflowCount}: every escape slot has its pair, ranks ascending, values at or above the escape");
        r.Row("E4", caseName, "sample_mismatches", mismatches, "n", $"{checkedRanks} ranks (first, last and {samples} random, seed {o.Long("seed", 20240011)}) compared with the LCE of SA[r-1], SA[r] on the packed text");
        return problems == 0 && mismatches == 0 ? 0 : 1;
    }
}
