using System.Diagnostics;
using System.Globalization;
using SuffixTree.Persistent;
using ST = SuffixTree.SuffixTree;

namespace SuffixTree.PerfHarness;

/// <summary>
/// The measurements of one child process. Line protocol (stdout): <c>R\top\tmetric\tvalue\tunit</c>, identical to the
/// harness that produced BENCHMARK.md (<c>analysis/bench/seqeron/SeqeronBench</c>), so rows are comparable.
/// One process = one measurement run: memory metrics need a fresh process.
/// </summary>
internal static class Measurements
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    private static void Emit(string op, string metric, double value, string unit)
        => Console.WriteLine($"R\t{op}\t{metric}\t{value.ToString("R", Inv)}\t{unit}");

    // ───────────────────────── /proc based memory (Linux) ─────────────────────────

    private static long ReadHighWaterMarkBytes()
    {
        try
        {
            foreach (string line in File.ReadLines("/proc/self/status"))
            {
                if (line.StartsWith("VmHWM:", StringComparison.Ordinal))
                {
                    string[] parts = line.Split((char[]?)null, StringSplitOptions.RemoveEmptyEntries);
                    return long.Parse(parts[1], Inv) * 1024;
                }
            }
        }
        catch (IOException)
        {
            // /proc is not available (not Linux): the memory metric is simply not reported.
        }
        return 0;
    }

    private static void EmitPeakRss(string op, string metric)
    {
        long hwm = ReadHighWaterMarkBytes();
        if (hwm > 0)
            Emit(op, metric, hwm, "bytes");
    }

    // ───────────────────────── query helpers ─────────────────────────

    private static string[] Lines(string path) => File.ReadAllLines(path).Where(s => s.Length > 0).ToArray();

    /// <summary>Runs <paramref name="body"/> (one pass over the queries) after an untimed warm-up pass until at least one second elapsed.</summary>
    private static (double SecondsPerPass, long Check) Passes(Func<long> body)
    {
        body(); // untimed warm-up pass
        long check;
        var sw = Stopwatch.StartNew();
        int passes = 0;
        do
        {
            check = body();
            passes++;
        }
        while (sw.Elapsed.TotalSeconds < 1.0);
        return (sw.Elapsed.TotalSeconds / passes, check);
    }

    private static long CountPass(ISuffixTree tree, string[] queries)
    {
        long sum = 0;
        foreach (string q in queries)
            sum += tree.CountOccurrences(q);
        return sum;
    }

    private static long LocatePass(ISuffixTree tree, string[] queries, out long hits)
    {
        long sum = 0;
        hits = 0;
        foreach (string q in queries)
        {
            var r = tree.FindAllOccurrences(q);
            hits += r.Count;
            for (int j = 0; j < r.Count; j++)
                sum += r[j];
        }
        return sum;
    }

    private static void QuerySuite(ISuffixTree tree, string queryPrefix)
    {
        string[] present = Lines(queryPrefix + ".present20.txt");
        string[] absent = Lines(queryPrefix + ".absent20.txt");
        string[] locate = Lines(queryPrefix + ".locate8.txt");

        // warm-up of the JIT and the caches, as in the original harness
        string[] warm = present.Length > 20000 ? present[..20000] : present;
        CountPass(tree, warm);
        CountPass(tree, warm);
        LocatePass(tree, locate.Length > 200 ? locate[..200] : locate, out _);

        var (spp, check) = Passes(() => CountPass(tree, present));
        Emit("count_present20", "throughput", present.Length / spp, "queries/s");
        Emit("count_present20", "checksum", check, "count");

        (spp, check) = Passes(() => CountPass(tree, absent));
        Emit("count_absent20", "throughput", absent.Length / spp, "queries/s");
        Emit("count_absent20", "checksum", check, "count");

        long hits = 0;
        (spp, check) = Passes(() => LocatePass(tree, locate, out hits));
        Emit("locate8", "throughput", locate.Length / spp, "patterns/s");
        Emit("locate8", "total_hits", hits, "count");
        Emit("locate8", "checksum", check, "count");

        // per-query latency of CountOccurrences (present 20-mers), one thread
        var ticks = new long[present.Length];
        for (int i = 0; i < present.Length; i++)
        {
            long a = Stopwatch.GetTimestamp();
            tree.CountOccurrences(present[i]);
            ticks[i] = Stopwatch.GetTimestamp() - a;
        }
        Array.Sort(ticks);
        double ToNs(long tk) => tk * 1e9 / Stopwatch.Frequency;
        Emit("count_present20", "latency_p50", ToNs(ticks[ticks.Length / 2]), "ns");
        Emit("count_present20", "latency_p99", ToNs(ticks[(int)(ticks.Length * 0.99)]), "ns");

        Concurrent(tree, present, 4, "count_present20_4t");
    }

    /// <summary>Every thread runs the full query set <c>reps</c> times on ONE shared tree instance.</summary>
    private static void Concurrent(ISuffixTree tree, string[] queries, int threads, string op)
    {
        var sw0 = Stopwatch.StartNew();
        long expect = CountPass(tree, queries);
        double single = sw0.Elapsed.TotalSeconds;
        int reps = Math.Min(400, Math.Max(1, (int)Math.Ceiling(1.0 / Math.Max(single, 1e-6))));
        var sums = new long[threads];
        using var gate = new ManualResetEventSlim(false);
        using var ready = new CountdownEvent(threads);
        var workers = new Thread[threads];
        for (int i = 0; i < threads; i++)
        {
            int k = i;
            workers[i] = new Thread(() =>
            {
                ready.Signal();
                gate.Wait();
                long s = 0;
                for (int r = 0; r < reps; r++)
                    s = CountPass(tree, queries);
                sums[k] = s;
            });
            workers[i].Start();
        }
        ready.Wait();
        var sw = Stopwatch.StartNew();
        gate.Set();
        foreach (var w in workers)
            w.Join();
        sw.Stop();
        Emit(op, "throughput", (double)threads * queries.Length * reps / sw.Elapsed.TotalSeconds, "queries/s");
        Emit(op, "checksum_ok", sums.All(x => x == expect) ? 1 : 0, "bool");
    }

    private static void AnalyticsSuite(ISuffixTree tree, string queryPrefix)
    {
        Emit("lrs", "length", tree.LongestRepeatedSubstring().Length, "chars");
        string query = File.ReadAllText(queryPrefix + ".memq100k.txt").Trim();
        var mem = tree.FindMaximalExactMatches(query, 20);
        Emit("mem_q100k_min20", "matches", mem.Count, "count");
        long ck = 0;
        foreach (var m in mem)
            ck += m.PositionInText * 31L + m.PositionInQuery * 17L + m.Length;
        Emit("mem_q100k_min20", "checksum", ck, "count");
        Emit("mum_q100k_min20", "matches", tree.FindMaximalUniqueMatches(query, 20).Count, "count");
        Emit("anchors_q100k_min20", "matches", tree.FindExactMatchAnchors(query, 20).Count, "count");
    }

    // ───────────────────────── child modes ─────────────────────────

    /// <summary>In-memory tree: build, memory, queries. Arguments: text file, query prefix.</summary>
    public static int InMemory(string textPath, string queryPrefix)
    {
        string text = File.ReadAllText(textPath);
        GC.Collect();
        GC.WaitForPendingFinalizers();
        GC.Collect();
        long heapBefore = GC.GetTotalMemory(true);

        var sw = Stopwatch.StartNew();
        var tree = ST.Build(text);
        sw.Stop();

        EmitPeakRss("build", "peak_rss_hwm_after_build");
        Emit("build", "time", sw.Elapsed.TotalSeconds, "s");
        GC.Collect();
        GC.WaitForPendingFinalizers();
        GC.Collect();
        long heapAfter = GC.GetTotalMemory(true);
        Emit("build", "bytes_per_char", (double)(heapAfter - heapBefore) / text.Length, "B/char");
        Emit("build", "node_count", tree.NodeCount, "count");
        Emit("build", "leaf_count", tree.LeafCount, "count");

        QuerySuite(tree, queryPrefix);
        AnalyticsSuite(tree, queryPrefix);
        GC.KeepAlive(tree);
        return 0;
    }

    /// <summary>Persistent build into a file. Arguments: text file, tree file.</summary>
    public static int PersistentBuild(string textPath, string treePath)
    {
        string text = File.ReadAllText(textPath);
        foreach (string f in new[] { treePath, treePath + ".children.tmp", treePath + ".depth.tmp" })
        {
            if (File.Exists(f))
                File.Delete(f);
        }

        var sw = Stopwatch.StartNew();
        var tree = PersistentSuffixTreeFactory.CreatePersistent(new StringTextSource(text), treePath);
        sw.Stop();

        Emit("build", "time", sw.Elapsed.TotalSeconds, "s");
        Emit("build", "node_count", tree.NodeCount, "count");
        Emit("build", "leaf_count", tree.LeafCount, "count");
        tree.Dispose();
        Emit("build", "bytes_per_char", (double)new FileInfo(treePath).Length / text.Length, "B/char");
        EmitPeakRss("build_process", "peak_rss_hwm");
        return 0;
    }

    /// <summary>Persistent queries on an existing file. Arguments: tree file, query prefix.</summary>
    public static int PersistentQuery(string treePath, string queryPrefix)
    {
        using var tree = PersistentSuffixTreeFactory.LoadPersistent(treePath);
        QuerySuite(tree, queryPrefix);
        AnalyticsSuite(tree, queryPrefix);
        return 0;
    }

    /// <summary>
    /// Load time of the default (header-only) path, then of the opt-in strict validation. Arguments: tree file.
    /// </summary>
    public static int PersistentLoad(string treePath)
    {
        var sw = Stopwatch.StartNew();
        using (var tree = PersistentSuffixTreeFactory.LoadPersistent(treePath))
            Emit("warm_load", "time", sw.Elapsed.TotalSeconds, "s");

        sw.Restart();
        using (var tree = PersistentSuffixTreeFactory.LoadPersistent(treePath, validateStructure: true))
            Emit("strict_load", "time", sw.Elapsed.TotalSeconds, "s");
        return 0;
    }
}
