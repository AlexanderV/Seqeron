using System.Diagnostics;
using System.Globalization;
using System.Text;
using SuffixTree.EmProto.Dcs;
using SuffixTree.EmProto.Harness;
using SuffixTree.EmProto.Build;
using SuffixTree.EmProto.Ph6;
using SuffixTree.EmProto.Text;
using SuffixTree.Esa;
using SuffixTree.Esa.Bench.Links;

namespace SuffixTree.EmProto.Experiments;

/// <summary>
/// E9 of ADR-004 section 5: PH6, the navigation pack input Psi, built from the BWT by sigma-way distribution (<c>e9-psi</c>) against the route through the inverse permutation (<c>e9-isa</c> = ISA by DIST, then <c>e9-psi-isa</c>),
/// the BWT pass that feeds it (<c>e9-bwt</c>, route (ii): a separate pass over the SA file and the packed text), the in-RAM tables of the suffix-link prototype P1 as oracle (<c>e9-oracle</c>) and the smoke test of the P1 link kernels
/// on the streamed tables (<c>e9-smoke</c>). Every command is one measured phase: wall time, the program's own byte counters (explicit reads and writes, writes through mapped windows), <c>/proc/self/io</c>
/// (rchar + wchar and the block-layer counters), RssAnon and the dirty pages of the output.
/// </summary>
internal static class E9
{
    private const string Exp = "E9";

    /// <summary>One measured phase: counters are zeroed at the start, the rows are written at the end.</summary>
    private sealed class Probe : IDisposable
    {
        private readonly RssSampler _rss = new(10);
        private readonly DirtySampler? _dirty;
        private readonly bool _rssOn;
        private readonly ProcIoSnapshot _io0;
        private readonly long _t0;

        /// <param name="outputPrefix">Path of the output file whose mappings the dirty-page sampler follows.</param>
        /// <param name="samplers">False (<c>--no-sampler</c>): neither the RssAnon sampler nor the dirty-page sampler runs inside the process, so that rchar + wchar of <c>/proc/self/io</c> contain no read of <c>/proc/self/status</c> or <c>smaps</c>.</param>
        internal Probe(string outputPrefix, bool samplers = true)
        {
            IoCounters.Reset();
            MappedIo.Reset();
            ScratchTracker.Reset();
            if (samplers)
            {
                _dirty = new DirtySampler(20, filePrefix: Path.GetFullPath(outputPrefix));
                _dirty.Start();
                _rssOn = RssSampler.IsSupported;
                if (_rssOn) _rss.Start();
            }
            _io0 = ProcIoSnapshot.Read();
            _t0 = Stopwatch.GetTimestamp();
        }

        internal void Finish(Report r, string c, string phase, long n, string notes, ArenaBudget? budget = null)
        {
            double seconds = (Stopwatch.GetTimestamp() - _t0) / (double)Stopwatch.Frequency;
            ProcIoSnapshot io = ProcIoSnapshot.Read().Since(_io0);
            RssStats rs = _rssOn ? _rss.Stop() : default;
            DirtyStats ds = _dirty?.Stop() ?? default;
            double perChar = 1.0 / n;
            long counted = IoCounters.ReadBytes + IoCounters.WriteBytes + MappedIo.WriteBytes;
            r.Row(Exp, c, phase + "_s", seconds, "s", notes);
            r.Row(Exp, c, phase + "_ns_per_char", seconds * 1e9 / n, "ns", "seconds / n");
            r.Row(Exp, c, phase + "_io_explicit_read", IoCounters.ReadBytes, "bytes", $"{IoCounters.ReadBytes * perChar:F4} B/char (program counter, explicit reads, {IoCounters.ReadCalls} calls)");
            r.Row(Exp, c, phase + "_io_explicit_write", IoCounters.WriteBytes, "bytes", $"{IoCounters.WriteBytes * perChar:F4} B/char (program counter, explicit writes, {IoCounters.WriteCalls} calls)");
            r.Row(Exp, c, phase + "_io_mapped_write", MappedIo.WriteBytes, "bytes", $"{MappedIo.WriteBytes * perChar:F4} B/char (program counter, writes through mapped windows: {MappedIo.Windows} windows, largest {MappedIo.MaxWindowBytes} bytes)");
            r.Row(Exp, c, phase + "_io_counted_total", counted, "bytes", $"{counted * perChar:F4} B/char: explicit reads + explicit writes + mapped writes (reads through mappings are not counted)");
            r.Row(Exp, c, phase + "_proc_rchar", io.RChar, "bytes", $"/proc/self/io rchar during the phase, {io.RChar * perChar:F4} B/char");
            r.Row(Exp, c, phase + "_proc_wchar", io.WChar, "bytes", $"/proc/self/io wchar during the phase, {io.WChar * perChar:F4} B/char (mapped writes are not syscalls)");
            r.Row(Exp, c, phase + "_proc_read_bytes", io.ReadBytes, "bytes", $"/proc/self/io read_bytes (block layer, page-cache misses), {io.ReadBytes * perChar:F4} B/char");
            r.Row(Exp, c, phase + "_proc_write_bytes", io.WriteBytes, "bytes", $"/proc/self/io write_bytes (pages that went from clean to dirty), {io.WriteBytes * perChar:F4} B/char");
            if (_dirty is not null)
            {
                r.Row(Exp, c, phase + "_rssanon_peak", rs.MaxRssAnonMb, "MB", "max RssAnon of the phase at 10 ms (program sampler)");
                r.Row(Exp, c, phase + "_rssfile_peak", rs.MaxRssFileMb, "MB", "mapped file pages resident (reported, not counted)");
                r.Row(Exp, c, phase + "_dirty_peak", ds.MaxDeltaKb * 1024.0 / 1e6, "MB", $"max /proc/meminfo Dirty minus the baseline ({ds.BaselineKb * 1024.0 / 1e6:F2} MB) at 20 ms, system-wide");
                r.Row(Exp, c, phase + "_proc_dirty_peak", ds.MaxProcessDirtyKb * 1024.0 / 1e6, "MB", "max Private_Dirty + Shared_Dirty of the mappings of the output file (per-process view)");
            }
            if (budget is not null) r.Row(Exp, c, phase + "_arena_peak", budget.PeakBytes, "bytes", $"peak charged to the ArenaBudget (cap {budget.CapBytes >> 20} MiB), live at the end {budget.LiveBytes}");
            r.Row(Exp, c, phase + "_scratch_peak", ScratchTracker.Peak, "bytes", $"{ScratchTracker.Peak * perChar:F4} B/char: peak of the tracked scratch files");
        }

        public void Dispose()
        {
            _rss.Dispose();
            _dirty?.Dispose();
        }
    }

    private static string KernelLabel(LinkKernel k) => k switch
    {
        LinkKernel.Scan => "a_scan",
        LinkKernel.Descent => "b_descent",
        LinkKernel.RmqPsvNsv => "c_rmq",
        _ => "c2_psv",
    };

    private static (long N, int Bits, int Sigma) ReadMeta(string textPath)
    {
        (long n, int bits) = E3.ReadMeta(textPath);
        string[] meta = File.ReadAllLines(textPath + ".meta");
        string? sigma = meta.FirstOrDefault(l => l.StartsWith("sigma=", StringComparison.Ordinal));
        if (sigma is null) throw new ArgumentException($"{textPath}.meta has no sigma=");
        return (n, bits, int.Parse(sigma[6..], CultureInfo.InvariantCulture));
    }

    private static long FileRows(string saPath)
    {
        long bytes = new FileInfo(saPath).Length;
        if (bytes % 4 != 0 || bytes < 4) throw new ArgumentException($"{saPath} is not a suffix array file (n + 1 entries of 4 bytes)");
        return (bytes / 4) - 1;
    }

    private static bool Verify(Report r, string c, string what, string produced, string? reference)
    {
        if (reference is null) return true;
        (bool equal, long at) = E3b.CompareFiles(produced, reference);
        r.Row(Exp, c, what + "_equals_ref", equal ? "yes" : "NO", "", equal ? $"byte-identical to {Path.GetFileName(reference)} ({new FileInfo(produced).Length} bytes)" : $"first difference at byte {at}");
        return equal;
    }

    /// <summary><c>e9-bwt --text &lt;packed&gt; --sa &lt;sa file&gt; --out &lt;bwt file&gt;</c>: the BWT and its meta file from the suffix array file.</summary>
    internal static unsafe int Bwt(Options o, Report r)
    {
        string textPath = o.Require("text"), saPath = o.Require("sa"), outPath = o.Require("out");
        string c = o.Str("case", "e9");
        (long n, int bits, int sigma) = ReadMeta(textPath);
        var budget = new ArenaBudget(o.Long("budget-mb", 160) << 20);
        using var mapped = new MappedFile(textPath, writable: false);
        using var text = new PackedText(mapped.Pointer, n, bits, mapped, ownsOwner: false);
        BwtMeta meta;
        using (var probe = new Probe(outPath, !o.Flag("no-sampler")))
        {
            meta = BwtPass.Build(saPath, text, sigma, outPath, new BwtOptions { ReadBytes = o.Int("read-kb", 1024) * 1024, WindowBytes = o.Int("window-kb", 1024) * 1024, Budget = budget });
            probe.Finish(r, c, "bwt", n, "wall time of the BWT pass: BWT[j] = T[SA[j] - 1], SA read once, text read at random through its mapping", budget);
        }
        r.Row(Exp, c, "bwt_meta", meta.Text().Replace('\n', ' ').Trim(), "", $"n={n} sigma={sigma} bits={bits}; primary row = the row of the suffix 0");
        return 0;
    }

    /// <summary><c>e9-psi --bwt &lt;bwt&gt; --out &lt;psi&gt;</c>: PH6, Psi by sigma-way distribution of the BWT.</summary>
    internal static int Psi(Options o, Report r)
    {
        string bwtPath = o.Require("bwt"), outPath = o.Require("out");
        string c = o.Str("case", "e9");
        BwtMeta meta = BwtMeta.Read(bwtPath);
        var budget = new ArenaBudget(o.Long("budget-mb", 160) << 20);
        var options = new PsiOptions
        {
            WindowBytes = o.Int("window-kb", 64) * 1024,
            SymbolsPerPass = o.Int("symbols-per-pass", 64),
            ReadBytes = o.Int("read-kb", 1024) * 1024,
            CyclicEntryZero = o.Flag("cyclic"),
            WindowCapBytes = o.Long("budget-mb", 160) << 20,
            Budget = budget,
        };
        PsiStats st;
        using (var probe = new Probe(outPath, !o.Flag("no-sampler")))
        {
            st = PsiStream.Build(bwtPath, meta, outPath, options);
            probe.Finish(r, c, "psi", meta.N, "wall time of the Psi pass: one BWT scan per group of symbols, one output cursor per symbol", budget);
        }
        r.Row(Exp, c, "psi_passes", st.Passes, "n", $"sigma={meta.Sigma}, symbols per pass {st.SymbolsPerPass} (limit {options.SymbolsPerPass}), window {Math.Max(1, options.WindowBytes / 4) * 4} bytes, {st.Windows} windows of the cursors, mapped working set bound {st.MappedBoundBytes} bytes");
        r.Row(Exp, c, "psi_bwt_bytes_read", st.BwtBytesRead, "bytes", $"{st.BwtBytesRead / (double)meta.N:F4} B/char: passes x (n + 1)");
        r.Row(Exp, c, "psi_bytes_written", st.PsiBytesWritten, "bytes", $"{st.PsiBytesWritten / (double)meta.N:F4} B/char: 4 (n + 1)");
        if (o.Get("pk-meta") is { } pkMeta)
        {
            string? line = File.ReadAllLines(pkMeta).FirstOrDefault(l => l.StartsWith("counts=", StringComparison.Ordinal));
            bool same = line is not null && line[7..] == string.Join(',', meta.Counts);
            string verdict = "n/a";
            if (line is not null) verdict = same ? "yes" : "NO";
            r.Row(Exp, c, "counts_equal_ph0", verdict, "", $"the symbol counts of the BWT meta against the C[] of PH0 in {Path.GetFileName(pkMeta)}");
            if (line is not null && !same) return 1;
        }
        return Verify(r, c, "psi", outPath, o.Get("verify")) ? 0 : 1;
    }

    /// <summary><c>e9-isa --sa &lt;sa file&gt; --out &lt;isa&gt; --scratch &lt;dir&gt;</c>: ISA by DIST.</summary>
    internal static int Isa(Options o, Report r)
    {
        string saPath = o.Require("sa"), outPath = o.Require("out"), scratch = o.Require("scratch");
        string c = o.Str("case", "e9");
        long n = FileRows(saPath);
        long scaleN = o.Long("scale-n", 0);
        long chunk = scaleN > 0 ? Math.Max(3, (long)Math.Ceiling(E4.ChunkRef * n / (double)scaleN)) : E4.ChunkRef;
        if (o.Get("chunk") is { } cv) chunk = long.Parse(cv, CultureInfo.InvariantCulture);
        var budget = new ArenaBudget(o.Long("budget-mb", 160) << 20);
        var options = new DistOptions { ScratchDir = scratch, ChunkPositions = chunk, IoBufferBytes = o.Int("io-kb", 64) * 1024, ReadBytes = o.Int("read-kb", 1024) * 1024, WindowBytes = o.Int("window-kb", 1024) * 1024, Budget = budget };
        DistStats st;
        using (var probe = new Probe(outPath, !o.Flag("no-sampler")))
        {
            st = IsaByDist.Build(saPath, n, outPath, options);
            probe.Finish(r, c, "isa", n, "wall time of ISA by DIST: distribute (SA[i], i) by SA[i] / R, then scatter chunk by chunk", budget);
        }
        r.Row(Exp, c, "isa_chunks", st.Chunks, "n", $"R = {chunk} positions per chunk, P = {st.Chunks}, tuple buffer {options.IoBufferBytes} bytes; tuples written {st.TupleBytesWritten} bytes ({st.TupleBytesWritten / (double)n:F4} B/char), read {st.TupleBytesRead} bytes");
        return Verify(r, c, "isa", outPath, o.Get("verify")) ? 0 : 1;
    }

    /// <summary><c>e9-psi-isa --sa &lt;sa file&gt; --isa &lt;isa file&gt; --out &lt;psi&gt;</c>: Psi = ISA[SA[i] + 1] from the ISA file.</summary>
    internal static int PsiFromIsaCommand(Options o, Report r)
    {
        string saPath = o.Require("sa"), isaPath = o.Require("isa"), outPath = o.Require("out");
        string c = o.Str("case", "e9");
        long n = FileRows(saPath);
        var budget = new ArenaBudget(o.Long("budget-mb", 160) << 20);
        using (var probe = new Probe(outPath, !o.Flag("no-sampler")))
        {
            PsiFromIsa.Build(saPath, isaPath, n, outPath, o.Int("read-kb", 1024) * 1024, o.Int("window-kb", 1024) * 1024, budget);
            probe.Finish(r, c, "psi_from_isa", n, "wall time of Psi = ISA[SA[i] + 1]: SA read sequentially, ISA read at random through its mapping", budget);
        }
        return Verify(r, c, "psi_from_isa", outPath, o.Get("verify")) ? 0 : 1;
    }

    /// <summary><c>e9-oracle --sa &lt;sa file&gt; --out &lt;prefix&gt;</c>: the in-RAM ISA and Psi of P1 (<c>LinkTables</c>) as files <c>prefix.isa</c> and <c>prefix.psi</c> (PSI[0] = 0).</summary>
    internal static int Oracle(Options o, Report r)
    {
        string saPath = o.Require("sa"), prefix = o.Require("out");
        string c = o.Str("case", "e9");
        long n = FileRows(saPath);
        long t0 = Stopwatch.GetTimestamp();
        PsiFiles.WriteOracle(saPath, prefix + ".isa", prefix + ".psi");
        r.Row(Exp, c, "oracle_s", (Stopwatch.GetTimestamp() - t0) / (double)Stopwatch.Frequency, "s", $"LinkTables.BuildIsa / BuildPsi in RAM on the SA file of n = {n}; written as {Path.GetFileName(prefix)}.isa and .psi (PSI[0] = -1 of P1 stored as 0)");
        return 0;
    }

    /// <summary><c>e9-smoke --input &lt;text file&gt; --sa &lt;sa&gt; --isa &lt;isa&gt; --psi &lt;psi&gt; [--psi2 &lt;psi from isa&gt;]</c>: every kernel x source of P1 with the in-RAM tables and with the streamed ones.</summary>
    internal static int Smoke(Options o, Report r)
    {
        string input = o.Require("input"), saPath = o.Require("sa"), isaPath = o.Require("isa"), psiPath = o.Require("psi");
        string c = o.Str("case", "e9");
        int sample = o.Int("sample", 200_000);
        ulong seed = (ulong)o.Long("seed", 1);
        long budget = o.Long("link-budget", 3_000_000);
        int fanout = o.Int("fanout", 16);
        string text = Encoding.Latin1.GetString(File.ReadAllBytes(input));
        long t0 = Stopwatch.GetTimestamp();
        EsaIndex index = EsaIndex.Build(text, ChildLayoutKind.Tagged);
        r.Row(Exp, c, "smoke_index_s", (Stopwatch.GetTimestamp() - t0) / (double)Stopwatch.Frequency, "s", $"EsaIndex.Build (tagged child table) of n = {text.Length}");
        int[] sa = index.SuffixArray;
        bool saSame = sa.AsSpan().SequenceEqual(PsiFiles.ReadInts(saPath));
        r.Row(Exp, c, "smoke_sa_file_equals_index", saSame ? "yes" : "NO", "", "the suffix array file the streamed tables were built from against the SA of the in-RAM index");
        int[] isaRam = LinkTables.BuildIsa(sa);
        int[] psiRam = LinkTables.BuildPsi(sa, isaRam);
        int[] isaDist = PsiFiles.ReadInts(isaPath);
        bool isaSame = isaDist.AsSpan().SequenceEqual(isaRam);
        r.Row(Exp, c, "smoke_isa_equals_ram", isaSame ? "yes" : "NO", "", "ISA from DIST against LinkTables.BuildIsa");
        bool ok = saSame && isaSame;
        var variants = new List<(string Name, string Path)> { ("psi_stream", psiPath) };
        if (o.Get("psi2") is { } psi2) variants.Add(("psi_from_isa", psi2));
        foreach ((string name, string path) in variants)
        {
            int[] psi = PsiFiles.LoadForP1(path);
            bool psiSame = psi.AsSpan(1).SequenceEqual(psiRam.AsSpan(1)) && psi[0] == psiRam[0];
            r.Row(Exp, c, $"smoke_{name}_equals_ram", psiSame ? "yes" : "NO", "", $"{name} (entry 0 mapped 0 -> -1) against LinkTables.BuildPsi");
            long t1 = Stopwatch.GetTimestamp();
            SmokeResult res = LinkSmoke.Run(index, isaRam, psiRam, isaDist, psi, sample, seed, budget, fanout, p1Sampling: true, null);
            r.Row(Exp, c, $"smoke_{name}_s", (Stopwatch.GetTimestamp() - t1) / (double)Stopwatch.Frequency, "s", "all kernels x sources, two runs each (in-RAM tables, streamed tables)");
            r.Row(Exp, c, $"smoke_{name}_internal_nodes", res.NodeCount, "n", $"internal nodes of the index; header statistic {index.Statistics.InternalCount}; sample {res.SampleCount} (sample A of the P1 benchmark: --sample {sample} --seed {seed})");
            foreach (SmokeCell cell in res.Cells)
                r.Row(Exp, c, $"smoke_{name}_{KernelLabel(cell.Kernel)}_{cell.Source.ToString().ToLowerInvariant()}", cell.Mismatches, "mismatches", $"{cell.Nodes} nodes, {cell.GaveUp} links over the step budget {budget} (the same in both runs)");
            r.Row(Exp, c, $"smoke_{name}_mismatches", res.Mismatches, "n", "sum over the eight cells: success flag, target and step count of every link, tables under test against the in-RAM tables");
            ok &= psiSame && res.Mismatches == 0 && res.NodeCount == index.Statistics.InternalCount;
        }
        r.Row(Exp, c, "smoke_ok", ok ? "yes" : "NO", "", "all comparisons equal");
        return ok ? 0 : 1;
    }
}
