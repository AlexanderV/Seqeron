using System.Buffers.Binary;
using System.Diagnostics;
using System.Globalization;
using System.Security.Cryptography;
using SuffixTree;
using SuffixTree.Esa;

// The in-memory build, entry and Open measurements of the S1 gate (work package 1b.9). One process = one measurement; every command prints
// `R<TAB>op<TAB>metric<TAB>value<TAB>unit<TAB>notes` lines (the protocol of lib/benchlib.py). The input is loaded first and is not part of any timer.
//
//   build  <text file> <name> [sa_vectors.tsv]     the production in-memory pipeline stage by stage, exactly the three calls of EsaIndex.Build:
//                                                  TextPreparation.Prepare, SaisKernel.Build(prepared), LcpKernel.Build(Direct), then the child-table scan (EsaIndex.FromArrays
//                                                  with the default layout of the length); build_sa_s = prepare + SA-IS, build_total_s = SA + LCP + child (the gate metric)
//   lean   <text file> <name> [sa_vectors.tsv]     the same stages, but nothing is read or counted before the build (no /proc read, no GC counter): see the note at `Lean` in this file
//   steady <text file> <name> [builds=5]           the same stages <builds> times in one process, the last three reported
//   entry  <text file> <name>                      one call of the public factory EsaSuffixTreeFactory.CreateEsa(text) (what a user calls), one process
//   open   <container> <pattern file>              the first Open of a container in this process (cold JIT; nothing before it but the arguments), the first query after it, then the warm Open at the three levels
//   openp  <container> <pattern file>              the same with the pattern file read before the first Open (the order of the first session of 1b.9)
//
// The `Checked` and `Trusted` throughput and the facade are measured by abapp (analysis/bench/abapp), the child-table-free query rows by the Bench command `query`.
return Gate.Run(args);

static class Gate
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    private static void Row(string op, string metric, double v, string unit, string notes = "") => Console.WriteLine(string.Join('\t', "R", op, metric, v.ToString("R", Inv), unit, notes));
    private static void Row(string op, string metric, long v, string unit, string notes = "") => Console.WriteLine(string.Join('\t', "R", op, metric, v.ToString(Inv), unit, notes));
    private static void Row(string op, string metric, string v, string unit, string notes = "") => Console.WriteLine(string.Join('\t', "R", op, metric, v, unit, notes));

    private static double Seconds(long t0) => (Stopwatch.GetTimestamp() - t0) / (double)Stopwatch.Frequency;

    private static long StatusBytes(string key)
    {
        foreach (string line in File.ReadLines("/proc/self/status"))
            if (line.StartsWith(key + ":", StringComparison.Ordinal))
                return long.Parse(new string(line.Where(char.IsAsciiDigit).ToArray()), Inv) * 1024;
        return 0;
    }

    private static string LoadText(string path)
    {
        using var stream = File.OpenRead(path);
        return string.Create(checked((int)stream.Length), stream, static (span, s) =>
        {
            var buffer = new byte[1 << 16];
            int pos = 0, read;
            while ((read = s.Read(buffer)) > 0)
                for (int i = 0; i < read; i++) span[pos++] = (char)buffer[i];
        });
    }

    private sealed record Stages(PreparedText Prepared, int[] Sa, CompactLcp Lcp, EsaIndex Index, double Prepare, double Sais, double Lcpt, double Child)
    {
        public double SaSeconds => Prepare + Sais;
        public double Total => Prepare + Sais + Lcpt + Child;
    }

    private static Stages Stage(string text)
    {
        long t0 = Stopwatch.GetTimestamp();
        PreparedText prepared = TextPreparation.Prepare(text);
        double tp = Seconds(t0);
        long t1 = Stopwatch.GetTimestamp();
        int[] sa = SaisKernel.Build(prepared);
        double ts = Seconds(t1);
        long t2 = Stopwatch.GetTimestamp();
        CompactLcp lcp = LcpKernel.Build(LcpKernelKind.Direct, prepared, sa);
        double tl = Seconds(t2);
        long t3 = Stopwatch.GetTimestamp();
        EsaIndex index = EsaIndex.FromArrays(text, prepared, sa, lcp, EsaIndex.DefaultLayoutFor(text.Length));
        double tc = Seconds(t3);
        return new Stages(prepared, sa, lcp, index, tp, ts, tl, tc);
    }

    private static void Append(IncrementalHash hash, ReadOnlySpan<int> values)
    {
        Span<byte> buffer = stackalloc byte[4 * 1024];
        while (!values.IsEmpty)
        {
            int k = Math.Min(values.Length, 1024);
            for (int j = 0; j < k; j++) BinaryPrimitives.WriteInt32LittleEndian(buffer.Slice(4 * j, 4), values[j]);
            hash.AppendData(buffer.Slice(0, 4 * k));
            values = values.Slice(k);
        }
    }

    private static string SaHash(int[] sa)
    {
        using var hash = IncrementalHash.CreateHash(HashAlgorithmName.SHA256);
        Append(hash, sa.AsSpan(1));
        return Convert.ToHexStringLower(hash.GetHashAndReset());
    }

    private static string LcpHash(CompactLcp lcp, int n)
    {
        using var hash = IncrementalHash.CreateHash(HashAlgorithmName.SHA256);
        var block = new int[4096];
        var cursor = lcp.CreateCursor();
        for (int r = 1; r <= n;)
        {
            int k = Math.Min(block.Length, n - r + 1);
            for (int j = 0; j < k; j++) block[j] = cursor.Get(r + j);
            Append(hash, block.AsSpan(0, k));
            r += k;
        }
        return Convert.ToHexStringLower(hash.GetHashAndReset());
    }

    private static (string Sa, string Lcp)? Expected(string vectors, string name)
    {
        foreach (string line in File.ReadLines(vectors))
        {
            if (line.StartsWith('#')) continue;
            string[] f = line.Split('\t');
            if (f.Length == 6 && f[0] == name) return (f[4], f[5]);
        }
        return null;
    }

    // `build` reads /proc/self/status and the GC counters BEFORE the first stage (as the Bench command `build` does); the code of those calls is compiled cold just before the build. `lean` moves them behind the build:
    // the first cold build of a process then starts from the load of the text, as a call of the public factory does (`entry`). The two are reported side by side (1b.9: the cold build of a fresh process is
    // sensitive to what the process compiled just before it: the tiering of the kernels).
    private static int Build(string[] a, bool lean = false)
    {
        string text = LoadText(a[0]);
        string name = a[1];
        int n = text.Length;
        long allocated0 = 0;
        int g0 = 0, g1 = 0, g2 = 0;
        if (!lean)
        {
            Row("build", "input_chars", n, "chars");
            Row("build", "peak_rss_hwm_after_load", StatusBytes("VmHWM"), "bytes", "VmHWM after loading the text as a managed string (2 B/char), before the build");
            allocated0 = GC.GetTotalAllocatedBytes(precise: true);
            g0 = GC.CollectionCount(0); g1 = GC.CollectionCount(1); g2 = GC.CollectionCount(2);
        }
        Stages s = Stage(text);
        long allocated = lean ? 0 : GC.GetTotalAllocatedBytes(precise: true) - allocated0;
        long hwm = StatusBytes("VmHWM");
        if (lean) Row("build", "input_chars", n, "chars");
        const string note = "Stopwatch around the library call only; input already in memory; production kernels (SA-IS pointer, direct LCP) and the child scan of the default layout";
        Row("build", "build_prepare_s", s.Prepare, "s", note);
        Row("build", "build_sais_s", s.Sais, "s", "SaisKernel.Build(prepared)");
        Row("build", "build_sa_s", s.SaSeconds, "s", "prepare + SA-IS");
        Row("build", "build_lcp_s", s.Lcpt, "s", "LcpKernel.Build(Direct): Phi / PLCP into the compact LCP");
        Row("build", "build_child_s", s.Child, "s", $"EsaIndex.FromArrays: one scan for the child table ({EsaIndex.DefaultLayoutFor(n)}) and the statistics");
        Row("build", "build_total_s", s.Total, "s", "SA + LCP + child table: the S1 gate metric (<= 3.3 s at 10 Mbp)");
        Row("build", "throughput", n / s.Total, "chars/s", "SA + LCP + child");
        Row("build", "peak_rss_hwm_after_build", hwm, "bytes", "VmHWM of the whole process after the build (runtime, text as a string, SA, LCP, child table, transient arrays)");
        Row("build", "peak_rss_bytes_per_char", hwm / (double)n, "bytes/char", "gate of PROGRESS-S1: <= 24 B/char");
        if (!lean)
        {
            Row("build", "gc_allocated_bytes", allocated, "bytes", "GC.GetTotalAllocatedBytes delta during the build");
            Row("build", "gc_collections_gen0", GC.CollectionCount(0) - g0, "count");
            Row("build", "gc_collections_gen1", GC.CollectionCount(1) - g1, "count");
            Row("build", "gc_collections_gen2", GC.CollectionCount(2) - g2, "count");
        }

        long saBytes = 4L * (n + 1);
        long lcpBytes = s.Lcp.PayloadBytes;
        long textBytes = ((long)n * (s.Prepared.Width == SymbolWidth.Byte ? 1 : 2)) + (2L * s.Prepared.AlphabetSize);
        long childBytes = s.Index.Child.PayloadBytes;
        long core = saBytes + lcpBytes + textBytes + childBytes;
        Row("build", "index_bytes_sa", saBytes, "bytes", "int32 x (n + 1)");
        Row("build", "index_bytes_lcp", lcpBytes, "bytes", $"compact LCP kind {s.Lcp.Kind}, {s.Lcp.OverflowCount} escapes");
        Row("build", "index_bytes_text", textBytes, "bytes", $"remapped symbols ({s.Prepared.Width}) + alphabet table of {s.Prepared.AlphabetSize} code units");
        Row("build", "index_bytes_child", childBytes, "bytes", $"child table payload, layout {EsaIndex.DefaultLayoutFor(n)}");
        Row("build", "index_bytes_core", core, "bytes", "in-memory Core: text + SA + LCP + child table (no navigation pack)");
        Row("build", "index_bytes_core_per_char", core / (double)n, "bytes/char", "gate: <= 16 B/char incl. text (Core)");
        Row("build", "internal_nodes", s.Index.Statistics.InternalCount, "count");
        Row("build", "alphabet_size", s.Prepared.AlphabetSize, "count");
        string saHash = SaHash(s.Sa), lcpHash = LcpHash(s.Lcp, n);
        Row("build", "sa_sha256", saHash, "hex", "SHA-256 of SA[1..n] as little-endian int32 (= libsais SA)");
        Row("build", "lcp_sha256", lcpHash, "hex", "SHA-256 of LCP[1..n] as little-endian int32 (= libsais LCP)");
        if (a.Length > 2)
        {
            var e = Expected(a[2], name);
            if (e is null) Row("build", "vector_match", "NA", "bool", $"{name} is not listed in {Path.GetFileName(a[2])}");
            else
            {
                bool ok = e.Value.Sa == saHash && e.Value.Lcp == lcpHash;
                Row("build", "vector_match", ok ? 1 : 0, "bool", "SA and LCP hashes equal the libsais vectors of analysis/sa_vectors.tsv");
                if (!ok) Environment.ExitCode = 3;
            }
        }
        GC.KeepAlive(s);
        return Environment.ExitCode;
    }

    private static int Steady(string[] a)
    {
        string text = LoadText(a[0]);
        int builds = a.Length > 2 ? int.Parse(a[2], Inv) : 5;
        if (builds < 4) throw new ArgumentException("builds must be at least 4");
        for (int k = 1; k <= builds; k++)
        {
            Stages s = Stage(text);
            string op = k > builds - 3 ? "steady" : "warmup";
            string note = $"build {k} of {builds} in one process" + (k <= 2 ? " (JIT warm-up, not reported as steady state)" : "");
            Row(op, "build_sa_s", s.SaSeconds, "s", note);
            Row(op, "build_lcp_s", s.Lcpt, "s", note);
            Row(op, "build_child_s", s.Child, "s", note);
            Row(op, "build_total_s", s.Total, "s", note);
            GC.KeepAlive(s);
        }
        Row("steady", "peak_rss_hwm_after_builds", StatusBytes("VmHWM"), "bytes", $"VmHWM after {builds} builds");
        return 0;
    }

    private static int Entry(string[] a)
    {
        string text = LoadText(a[0]);
        int n = text.Length;
        Row("entry", "input_chars", n, "chars");
        long t0 = Stopwatch.GetTimestamp();
        EsaSuffixTree tree = EsaSuffixTreeFactory.CreateEsa(text);
        double t = Seconds(t0);
        Row("entry", "entry_build_s", t, "s", "EsaSuffixTreeFactory.CreateEsa(text): the public entry (prepare + SA-IS + LCP + child + facade), cold JIT");
        Row("entry", "peak_rss_hwm_after_entry", StatusBytes("VmHWM"), "bytes");
        long t1 = Stopwatch.GetTimestamp();
        int c = tree.CountOccurrences(text.Substring(n / 2, 20));
        Row("entry", "first_count_s", Seconds(t1), "s", $"first CountOccurrences of a 20-mer of the text, cold: result {c}");
        GC.KeepAlive(tree);
        return 0;
    }

    private static string[] FirstLines(string path, int k) => File.ReadLines(path).Take(k).ToArray();

    // `open` makes the very first Open of the process before anything else is read or compiled (the pattern file is read after it); `openp` reads the 1000 patterns first (the order of the first session of 1b.9).
    private static int Open(string[] a, bool patternsFirst)
    {
        string container = a[0];
        string[] patterns = patternsFirst ? FirstLines(a[1], 1000) : [];
        long t0 = Stopwatch.GetTimestamp();
        ISuffixTree first = SuffixIndexFactory.Open(container);
        double firstOpen = Seconds(t0);
        if (!patternsFirst) patterns = FirstLines(a[1], 1000);
        long size = new FileInfo(container).Length;
        Row("open", "container_bytes", size, "bytes");
        Row("open", "first_open_s", firstOpen, "s", "SuffixIndexFactory.Open at the Header level, Checked, the FIRST call of this process (cold JIT; the page cache is as the runner left it)" + (patternsFirst ? "; the pattern file was read before" : "; nothing but the arguments before it"));
        long t1 = Stopwatch.GetTimestamp();
        int c = first.CountOccurrences(patterns[0]);
        Row("open", "first_query_s", Seconds(t1), "s", $"the first CountOccurrences of a 20-mer after the open (cold): result {c}");
        long t2 = Stopwatch.GetTimestamp();
        long sum = 0;
        for (int i = 0; i < patterns.Length; i++) sum += first.CountOccurrences(patterns[i]);
        Row("open", "first_1000_queries_s", Seconds(t2), "s", $"1000 counts of present 20-mers after the open (warming the tiers): sum {sum}");
        (first as IDisposable)?.Dispose();

        const int reps = 200;
        var samples = new double[reps];
        for (int i = 0; i < reps; i++)
        {
            long t = Stopwatch.GetTimestamp();
            ISuffixTree tree = SuffixIndexFactory.Open(container);
            samples[i] = Seconds(t);
            (tree as IDisposable)?.Dispose();
        }
        Array.Sort(samples);
        Row("open", "warm_open_median_s", samples[reps / 2], "s", $"SuffixIndexFactory.Open, Header, Checked; median of {reps} opens in this process after the first (JIT warm, file cached)");
        Row("open", "warm_open_min_s", samples[0], "s");
        Row("open", "warm_open_max_s", samples[^1], "s");

        var facade = new double[reps];
        for (int i = 0; i < reps; i++)
        {
            long t = Stopwatch.GetTimestamp();
            EsaSuffixTree tree = EsaSuffixTreeFactory.Open(container);
            facade[i] = Seconds(t);
            tree.Dispose();
        }
        Array.Sort(facade);
        Row("open", "warm_open_facade_median_s", facade[reps / 2], "s", $"EsaSuffixTreeFactory.Open (the facade), median of {reps}");

        long t3 = Stopwatch.GetTimestamp();
        (SuffixIndexFactory.Open(container, new EsaOpenOptions { Validation = EsaValidationLevel.Sections }) as IDisposable)?.Dispose();
        Row("open", "open_sections_s", Seconds(t3), "s", "Open at the Sections level: one pass over the file with the CRC-32C of every section (warm file cache; the JIT of this path is cold)");
        long t4 = Stopwatch.GetTimestamp();
        (SuffixIndexFactory.Open(container, new EsaOpenOptions { Validation = EsaValidationLevel.Structure }) as IDisposable)?.Dispose();
        Row("open", "open_structure_s", Seconds(t4), "s", "Open at the Structure level: the structure re-derived in one pass (warm file cache; the JIT of this path is cold)");
        Row("open", "peak_rss_hwm", StatusBytes("VmHWM"), "bytes");
        return 0;
    }

    public static int Run(string[] args)
    {
        try
        {
            return args.Length == 0 ? 2 : args[0] switch
            {
                "build" => Build(args[1..]),
                "lean" => Build(args[1..], lean: true),
                "steady" => Steady(args[1..]),
                "entry" => Entry(args[1..]),
                "open" => Open(args[1..], patternsFirst: false),
                "openp" => Open(args[1..], patternsFirst: true),
                _ => 2,
            };
        }
        catch (Exception e)
        {
            Console.Error.WriteLine(e);
            return 1;
        }
    }
}
