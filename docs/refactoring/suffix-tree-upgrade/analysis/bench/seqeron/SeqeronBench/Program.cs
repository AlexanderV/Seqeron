using System.Diagnostics;
using System.Globalization;
using SuffixTree;
using SuffixTree.Persistent;
using ST = SuffixTree.SuffixTree;

// Seqeron suffix-tree benchmark harness.  One process = one measurement run; the Python driver
// (../../run_bench.py) repeats runs, takes medians and writes bench.tsv.
// Output protocol: lines "R\t<op>\t<metric>\t<value>\t<unit>\t<notes>" on stdout.
internal static class Program
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    private static void Emit(string op, string metric, double value, string unit, string notes = "")
        => Console.WriteLine($"R\t{op}\t{metric}\t{value.ToString("R", Inv)}\t{unit}\t{notes}");

    private static void EmitI(string op, string metric, long value, string unit, string notes = "")
        => Console.WriteLine($"R\t{op}\t{metric}\t{value.ToString(Inv)}\t{unit}\t{notes}");

    // ---------- /proc based memory sampling ----------
    private static long _peakAnon, _peakFile, _peakRss;
    private static volatile bool _stopSampler;

    private static void ReadStatus(out long rss, out long anon, out long file, out long hwm)
    {
        rss = anon = file = hwm = 0;
        foreach (var line in File.ReadLines("/proc/self/status"))
        {
            if (line.StartsWith("VmRSS:", StringComparison.Ordinal)) rss = Kb(line);
            else if (line.StartsWith("RssAnon:", StringComparison.Ordinal)) anon = Kb(line);
            else if (line.StartsWith("RssFile:", StringComparison.Ordinal)) file = Kb(line);
            else if (line.StartsWith("VmHWM:", StringComparison.Ordinal)) hwm = Kb(line);
        }
    }

    private static long Kb(string l)
    {
        var p = l.Split((char[]?)null, StringSplitOptions.RemoveEmptyEntries);
        return long.Parse(p[1], Inv) * 1024;
    }

    private static Thread StartSampler()
    {
        var t = new Thread(() =>
        {
            while (!_stopSampler)
            {
                try
                {
                    ReadStatus(out var rss, out var anon, out var file, out _);
                    if (rss > _peakRss) _peakRss = rss;
                    if (anon > _peakAnon) _peakAnon = anon;
                    if (file > _peakFile) _peakFile = file;
                }
                catch (IOException) { }
                Thread.Sleep(20);
            }
        })
        { IsBackground = true, Priority = ThreadPriority.BelowNormal };
        t.Start();
        return t;
    }

    private static void EmitProcessMemory(string op)
    {
        _stopSampler = true;
        ReadStatus(out var rss, out var anon, out var file, out var hwm);
        EmitI(op, "peak_rss_hwm", hwm, "bytes", "VmHWM at end of process (includes file-backed mmap pages)");
        EmitI(op, "peak_rss_anon_sampled", Math.Max(_peakAnon, anon), "bytes", "max RssAnon sampled every 20ms");
        EmitI(op, "peak_rss_file_sampled", Math.Max(_peakFile, file), "bytes", "max RssFile sampled every 20ms (mmap'd file pages)");
        using var p = Process.GetCurrentProcess();
        EmitI(op, "peak_working_set64", p.PeakWorkingSet64, "bytes", "Process.PeakWorkingSet64");
        EmitI(op, "gc_total_committed_end", GC.GetGCMemoryInfo().TotalCommittedBytes, "bytes", "managed heap committed at end");
    }

    // ---------- helpers ----------
    private static string[] Lines(string path)
        => File.ReadAllLines(path).Where(s => s.Length > 0).ToArray();

    private static string ReadQ(string path) => File.ReadAllText(path).Trim();

    /// Runs body (one full pass over the queries) repeatedly until >= minSeconds elapsed (>=1 pass).
    private static (double secPerPass, int passes, long check) Passes(Func<long> body, double minSeconds = 1.0)
    {
        long check = body(); // untimed full warm-up pass (tier-1 JIT, caches)
        var sw = Stopwatch.StartNew();
        int passes = 0;
        do { check = body(); passes++; } while (sw.Elapsed.TotalSeconds < minSeconds);
        sw.Stop();
        return (sw.Elapsed.TotalSeconds / passes, passes, check);
    }

    private static long CountPass(ISuffixTree t, string[] qs)
    {
        long s = 0;
        for (int i = 0; i < qs.Length; i++) s += t.CountOccurrences(qs[i]);
        return s;
    }

    private static long ContainsPass(ISuffixTree t, string[] qs)
    {
        long s = 0;
        for (int i = 0; i < qs.Length; i++) if (t.Contains(qs[i])) s++;
        return s;
    }

    private static long LocatePass(ISuffixTree t, string[] qs, out long hits)
    {
        long sum = 0; hits = 0;
        for (int i = 0; i < qs.Length; i++)
        {
            var r = t.FindAllOccurrences(qs[i]);
            hits += r.Count;
            for (int j = 0; j < r.Count; j++) sum += r[j];
        }
        return sum;
    }

    private static void Warm(ISuffixTree t, string[] present, string[] loc)
    {
        var w = present.Length > 20000 ? present[..20000] : present;
        CountPass(t, w); ContainsPass(t, w); CountPass(t, w);
        LocatePass(t, loc.Length > 200 ? loc[..200] : loc, out _);
    }

    private static void QuerySuite(ISuffixTree t, string qp, string tag)
    {
        var present = Lines(qp + ".present20.txt");
        var absent = Lines(qp + ".absent20.txt");
        var loc = Lines(qp + ".locate8.txt");
        var loc4Path = qp + ".locate4.txt";
        Warm(t, present, loc);

        void Rate(string op, int nq, Func<long> body, string unit = "queries/s", string notes = "")
        {
            var (spp, passes, chk) = Passes(body);
            Emit(op, "throughput", nq / spp, unit, $"{notes}passes={passes}".Trim());
            Emit(op, "time_per_pass", spp, "s", $"queries/pass={nq}");
            EmitI(op, "checksum", chk, "count", "sum of results (cross-tool oracle)");
        }

        Rate("count_present20", present.Length, () => CountPass(t, present));
        Rate("count_absent20", absent.Length, () => CountPass(t, absent));
        Rate("contains_present20", present.Length, () => ContainsPass(t, present));
        Rate("contains_absent20", absent.Length, () => ContainsPass(t, absent));
        {
            long hits = 0;
            var (spp, passes, chk) = Passes(() => LocatePass(t, loc, out hits));
            Emit("locate8", "throughput", loc.Length / spp, "patterns/s", $"passes={passes}");
            Emit("locate8", "positions_per_s", hits / spp, "positions/s", $"total_hits_per_pass={hits}");
            EmitI("locate8", "total_hits", hits, "count", "occurrences returned per pass (oracle)");
            EmitI("locate8", "checksum", chk, "count", "sum of positions (oracle)");
        }
        if (File.Exists(loc4Path))
        {
            var loc4 = Lines(loc4Path);
            long hits = 0;
            var (spp, passes, chk) = Passes(() => LocatePass(t, loc4, out hits));
            Emit("locate4", "throughput", loc4.Length / spp, "patterns/s", $"passes={passes}");
            Emit("locate4", "positions_per_s", hits / spp, "positions/s", $"total_hits_per_pass={hits}");
            EmitI("locate4", "total_hits", hits, "count", "occurrences returned per pass (oracle)");
            EmitI("locate4", "checksum", chk, "count", "sum of positions (oracle)");
        }

        // per-query latency of CountOccurrences (present 20-mers), single thread
        {
            var ticks = new long[present.Length];
            for (int i = 0; i < present.Length; i++)
            {
                long a = Stopwatch.GetTimestamp();
                t.CountOccurrences(present[i]);
                ticks[i] = Stopwatch.GetTimestamp() - a;
            }
            Array.Sort(ticks);
            double ns(long tk) => tk * 1e9 / Stopwatch.Frequency;
            Emit("count_present20", "latency_p50", ns(ticks[ticks.Length / 2]), "ns", "per-query Stopwatch (includes ~25ns timer overhead)");
            Emit("count_present20", "latency_p99", ns(ticks[(int)(ticks.Length * 0.99)]), "ns", "");
            Emit("count_present20", "latency_mean", ticks.Select(ns).Average(), "ns", "");
        }

        // concurrent read throughput (count_present20, locate8) on ONE shared tree
        foreach (int th in new[] { 1, 4 })
        {
            Concurrent(t, present, th, $"count_present20_{th}t", (tree, q) => CountPass(tree, q));
            Concurrent(t, loc, th, $"locate8_{th}t", (tree, q) => LocatePass(tree, q, out _), "patterns/s");
        }
        _ = tag;
    }

    private static void Concurrent(ISuffixTree t, string[] qs, int threads, string op,
        Func<ISuffixTree, string[], long> pass, string unit = "queries/s")
    {
        // every thread runs the FULL query set `reps` times (identical per-thread working set for 1 and 4 threads)
        var sw0 = Stopwatch.StartNew();
        long expect = pass(t, qs);
        double single = sw0.Elapsed.TotalSeconds;
        int reps = Math.Min(400, Math.Max(1, (int)Math.Ceiling(1.0 / Math.Max(single, 1e-6))));
        long[] sums = new long[threads];
        using var gate = new ManualResetEventSlim(false);
        using var ready = new CountdownEvent(threads);
        var ths = new Thread[threads];
        for (int i = 0; i < threads; i++)
        {
            int k = i;
            ths[i] = new Thread(() =>
            {
                ready.Signal(); gate.Wait();
                long s = 0;
                for (int r = 0; r < reps; r++) s = pass(t, qs);
                sums[k] = s;
            });
            ths[i].Start();
        }
        ready.Wait();
        var sw = Stopwatch.StartNew();
        gate.Set();
        foreach (var th in ths) th.Join();
        sw.Stop();
        bool ok = sums.All(x => x == expect);
        Emit(op, "throughput", (double)threads * qs.Length * reps / sw.Elapsed.TotalSeconds, unit,
            $"threads={threads} reps={reps} each thread runs the full query set; one shared tree instance");
        EmitI(op, "checksum_ok", ok ? 1 : 0, "bool", "every thread's result equals the single-thread result");
    }

    private static void AnalyticsSuite(ISuffixTree t, string qp)
    {
        var sw = Stopwatch.StartNew();
        var lrs = t.LongestRepeatedSubstring();
        Emit("lrs", "first_call_time", sw.Elapsed.TotalSeconds, "s", "LongestRepeatedSubstring first call");
        EmitI("lrs", "length", lrs.Length, "chars", "");
        sw.Restart();
        var lrs2 = t.LongestRepeatedSubstring();
        Emit("lrs", "second_call_time", sw.Elapsed.TotalSeconds, "s", "cached");
        _ = lrs2;

        var q = ReadQ(qp + ".memq100k.txt");
        sw.Restart();
        var mem = t.FindMaximalExactMatches(q, 20);
        Emit("mem_q100k_min20", "time", sw.Elapsed.TotalSeconds, "s", "FindMaximalExactMatches, 100 kbp query w/ 1% substitutions, minLength 20 (first call, JIT included)");
        EmitI("mem_q100k_min20", "matches", mem.Count, "count", "complete MEM set size (oracle)");
        long ck = 0; foreach (var m in mem) ck += m.PositionInText * 31L + m.PositionInQuery * 17L + m.Length;
        EmitI("mem_q100k_min20", "checksum", ck, "count", "");
        sw.Restart();
        mem = t.FindMaximalExactMatches(q, 20);
        Emit("mem_q100k_min20", "time_second_call", sw.Elapsed.TotalSeconds, "s", "second call (JIT warm)");
        sw.Restart();
        var mum = t.FindMaximalUniqueMatches(q, 20);
        Emit("mum_q100k_min20", "time", sw.Elapsed.TotalSeconds, "s", "FindMaximalUniqueMatches(Both), same query");
        EmitI("mum_q100k_min20", "matches", mum.Count, "count", "");
        sw.Restart();
        var anchors = t.FindExactMatchAnchors(q, 20);
        Emit("anchors_q100k_min20", "time", sw.Elapsed.TotalSeconds, "s", "FindExactMatchAnchors (matching-statistics anchors)");
        EmitI("anchors_q100k_min20", "matches", anchors.Count, "count", "");
    }

    private sealed class StageProgress : IProgress<(string Stage, double Percent)>
    {
        public readonly List<(string Stage, double At)> Marks = new();
        private readonly Stopwatch _sw;
        private string _last = "";
        public StageProgress(Stopwatch sw) => _sw = sw;
        public void Report((string Stage, double Percent) v)
        {
            if (v.Stage == _last) return;
            _last = v.Stage;
            Marks.Add((v.Stage, _sw.Elapsed.TotalSeconds));
        }
    }

    private static void EmitStages(StageProgress sp, double total, string op)
    {
        for (int i = 0; i < sp.Marks.Count; i++)
        {
            double end = i + 1 < sp.Marks.Count ? sp.Marks[i + 1].At : total;
            string name = sp.Marks[i].Stage.Replace(' ', '_').Replace(':', '_');
            Emit(op, "stage_" + name, end - sp.Marks[i].At, "s", "builder progress-callback stage");
        }
    }

    // ---------- modes ----------
    private static int Main(string[] a)
    {
        if (a.Length < 2) { Console.Error.WriteLine("usage: mem|pers-build|pers-load|pers-query ..."); return 2; }
        StartSampler();
        string mode = a[0];
        switch (mode)
        {
            case "mem": return Mem(a);              // mem <text> <qprefix>
            case "pers-build": return PersBuild(a); // pers-build <text> <treefile>
            case "pers-load": return PersLoad(a);   // pers-load <treefile> <qprefix>   (driver drops caches for cold)
            case "pers-query": return PersQuery(a); // pers-query <treefile> <qprefix> [analytics]
            default: Console.Error.WriteLine("bad mode"); return 2;
        }
    }

    private static int Mem(string[] a)
    {
        string text = File.ReadAllText(a[1]);
        string qp = a[2];
        GC.Collect(); GC.WaitForPendingFinalizers(); GC.Collect();
        long h0 = GC.GetTotalMemory(true);
        long alloc0 = GC.GetTotalAllocatedBytes(true);
        var sw = Stopwatch.StartNew();
        var tree = ST.Build(text);
        sw.Stop();
        EmitProcessMemoryAfterBuild();
        Emit("build", "time", sw.Elapsed.TotalSeconds, "s", "SuffixTree.Build(string) wall time (Ukkonen, in-memory)");
        Emit("build", "throughput", text.Length / sw.Elapsed.TotalSeconds, "chars/s", "");
        long alloc1 = GC.GetTotalAllocatedBytes(true);
        GC.Collect(); GC.WaitForPendingFinalizers(); GC.Collect();
        long h1 = GC.GetTotalMemory(true);
        EmitI("build", "index_bytes_heap_delta", h1 - h0, "bytes", "GC.GetTotalMemory(true) after minus before Build; text string (2 B/char) excluded");
        Emit("build", "bytes_per_char", (double)(h1 - h0) / text.Length, "bytes/char", "heap delta / n");
        EmitI("build", "total_allocated_during_build", alloc1 - alloc0, "bytes", "GC.GetTotalAllocatedBytes delta");
        EmitI("build", "text_string_bytes", 2L * text.Length, "bytes", "UTF-16 string of the input (kept alive; not in index_bytes)");
        EmitI("build", "node_count", tree.NodeCount, "count", "");
        EmitI("build", "leaf_count", tree.LeafCount, "count", "");
        QuerySuite(tree, qp, "mem");
        AnalyticsSuite(tree, qp);
        EmitProcessMemory("process");
        GC.KeepAlive(tree);
        return 0;
    }

    private static void EmitProcessMemoryAfterBuild()
    {
        ReadStatus(out var rss, out var anon, out _, out var hwm);
        EmitI("build", "peak_rss_hwm_after_build", hwm, "bytes", "VmHWM right after build, before any query");
        EmitI("build", "peak_rss_anon_after_build", Math.Max(_peakAnon, anon), "bytes", "max RssAnon sampled up to end of build");
        _ = rss;
    }

    private static int PersBuild(string[] a)
    {
        string text = File.ReadAllText(a[1]);
        string path = a[2];
        foreach (var f in new[] { path, path + ".children.tmp", path + ".depth.tmp" }) if (File.Exists(f)) File.Delete(f);
        GC.Collect();
        long h0 = GC.GetTotalMemory(true);
        var sw = Stopwatch.StartNew();
        var sp = new StageProgress(sw);
        var tree = PersistentSuffixTreeFactory.CreatePersistent(new StringTextSource(text), path, sp);
        sw.Stop();
        double total = sw.Elapsed.TotalSeconds;
        long h1 = GC.GetTotalMemory(false);
        Emit("build", "time", total, "s", "PersistentSuffixTreeFactory.CreatePersistent(StringTextSource, file) wall time incl. finalize+trim");
        Emit("build", "throughput", text.Length / total, "chars/s", "");
        EmitStages(sp, total, "build");
        EmitI("build", "managed_heap_delta", h1 - h0, "bytes", "GC heap growth during build (graph lives in mmap, not heap); text string excluded");
        EmitI("build", "node_count", tree.NodeCount, "count", "");
        EmitI("build", "leaf_count", tree.LeafCount, "count", "");
        tree.Dispose();
        long size = new FileInfo(path).Length;
        EmitI("build", "index_bytes_file", size, "bytes", "final file size (nodes+child arrays+text; text stored as 1 B/char ASCII)");
        Emit("build", "bytes_per_char", (double)size / text.Length, "bytes/char", "file bytes / n");
        EmitI("build", "temp_files_left", new[] { path + ".children.tmp", path + ".depth.tmp" }.Count(File.Exists), "count", "");
        EmitProcessMemory("process");
        return 0;
    }

    private static int PersLoad(string[] a)
    {
        string path = a[1];
        string qp = a[2];
        var present = Lines(qp + ".present20.txt");
        var sw = Stopwatch.StartNew();
        var tree = PersistentSuffixTreeFactory.LoadPersistent(path);
        double load = sw.Elapsed.TotalSeconds;
        Emit("load", "time", load, "s", "PersistentSuffixTreeFactory.LoadPersistent (mmap+header validation+madvise prefetch)");
        sw.Restart();
        int c = tree.CountOccurrences(present[0]);
        Emit("load", "first_query_time", sw.Elapsed.TotalSeconds, "s", "first CountOccurrences after load (includes JIT)");
        sw.Restart();
        long s = CountPass(tree, present[..1000]);
        Emit("load", "first_1000_queries_time", sw.Elapsed.TotalSeconds, "s", "next 1000 count queries");
        sw.Restart();
        s += CountPass(tree, present);
        Emit("load", "first_100k_queries_time", sw.Elapsed.TotalSeconds, "s", "one pass of 100k count queries right after load (cache state as set by driver)");
        EmitI("load", "checksum", s, "count", $"c0={c}");
        EmitProcessMemory("process");
        tree.Dispose();
        return 0;
    }

    private static int PersQuery(string[] a)
    {
        string path = a[1];
        string qp = a[2];
        bool analytics = a.Length > 3 && a[3] == "analytics";
        using var tree = PersistentSuffixTreeFactory.LoadPersistent(path);
        QuerySuite(tree, qp, "pers");
        if (analytics) AnalyticsSuite(tree, qp);
        EmitProcessMemory("process");
        return 0;
    }
}
