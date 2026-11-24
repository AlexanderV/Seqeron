using System.Diagnostics;
using System.Runtime.CompilerServices;

// The synthetic access patterns of ADR-005 section 5 on arrays of n entries in one layout (Layouts.cs):
//   chase   a chain of dependent random reads (the next index is the value just read): the latency of one random read, the pattern of an ISA / PSI hop and of one step of a descent
//   thr     independent random reads (the indices come from a generator, not from the data): the throughput that the memory-level parallelism allows, the pattern of a batch of queries
//   scan    m consecutive entries from a random start (the locate of an SA slice): ns per entry
//   null    the same loops on a cache-resident array of 4,096 entries (the decode and loop floor) and with no read at all
// Each measurement is repeated and reported as median / min / max of ns per operation; 4 threads run the same loop on one array from different seeds.

internal static unsafe class Synthetic
{
    internal sealed record Row(string Kind, string Layout, long N, int Threads, string Metric, double Median, double Min, double Max, string Unit, int Reps, string Notes);

    internal static List<Row> RunLayout(string layout, long n, ulong v, string[] kinds, int[] threadCounts, int reps, long ops, int[] scanLengths, string tag)
    {
        var rows = new List<Row>();
        using Packed a = Packed.Create(layout, n, v, 0x1234UL);
        ulong touch = a.Touch();
        double bytesPerEntry = (double)a.Bytes / n;
        string notes = $"{tag} bytes={a.Bytes} B/entry={bytesPerEntry:F3} escapes={a.Escapes} V={v} touch={touch & 0xFF}";
        foreach (int t in threadCounts)
        {
            foreach (string kind in kinds)
            {
                if (kind == "scan")
                {
                    foreach (int m in scanLengths)
                        rows.Add(Measure(layout, a, kind, n, t, reps, ops * 10, m, notes));
                }
                else rows.Add(Measure(layout, a, kind, n, t, reps, kind == "thr" ? ops * 10 : ops, 0, notes));
            }
        }
        return rows;
    }

    private static Row Measure(string layout, Packed a, string kind, long n, int threads, int reps, long ops, int m, string notes)
    {
        var samples = new List<double>();
        for (int r = 0; r < reps; r++)
            samples.Add(RunOnce(layout, a, kind, n, threads, ops, m, (ulong)(r + 1)));
        samples.Sort();
        string metric = kind == "scan" ? $"scan{m}_ns_per_entry" : kind + "_ns_per_read";
        return new Row(kind, layout, n, threads, metric, samples[samples.Count / 2], samples[0], samples[^1], kind == "scan" ? "ns/entry" : "ns/read", reps, notes);
    }

    // returns ns per operation of the whole run (wall clock / operations summed over the threads x threads, i.e. the per-thread latency when threads > 1 for chase; for thr and scan the aggregate time per op)
    private static double RunOnce(string layout, Packed a, string kind, long n, int threads, long ops, int m, ulong seed)
    {
        long perThread = ops / threads;
        double[] ns = new double[threads];
        var barrier = new Barrier(threads + 1);
        var ts = new Thread[threads];
        ulong[] sink = new ulong[threads * 8];
        for (int t = 0; t < threads; t++)
        {
            int id = t;
            ts[t] = new Thread(() =>
            {
                barrier.SignalAndWait();
                var sw = Stopwatch.StartNew();
                ulong s = Dispatch(layout, a, kind, n, perThread, m, (seed << 8) + (ulong)id);
                sw.Stop();
                ns[id] = sw.Elapsed.TotalMilliseconds * 1e6;
                sink[id * 8] = s;
            });
            ts[t].Start();
        }
        var wall = Stopwatch.StartNew();
        barrier.SignalAndWait();
        foreach (Thread t in ts) t.Join();
        wall.Stop();
        GC.KeepAlive(sink);
        double totalOps = (double)perThread * threads;
        // chase: report the latency seen by each thread (time of the slowest thread / its own reads); thr and scan: aggregate time per operation
        if (kind == "chase") return ns.Max() / perThread;
        return wall.Elapsed.TotalMilliseconds * 1e6 / totalOps;
    }

    private static ulong Dispatch(string layout, Packed a, string kind, long n, long ops, int m, ulong seed)
    {
        switch (layout)
        {
            case "w4": return Loop(new R4(a.Main), kind, n, ops, m, seed);
            case "w5": return Loop(new R5(a.Main), kind, n, ops, m, seed);
            case "w8": return Loop(new R8(a.Main), kind, n, ops, m, seed);
            case "p34": return Loop(new R34(a.Main), kind, n, ops, m, seed);
            case "hyb": return LoopHyb(new RHyb(a.Main, a.Blocks, a.Longs), kind, n, ops, m, seed);
            default: throw new ArgumentException(layout);
        }
    }

    private static ulong LoopHyb(RHyb r, string kind, long n, long ops, int m, ulong seed)
    {
        if (kind != "scan") return Loop(r, kind, n, ops, m, seed);
        // the scan of the hybrid is the sequential decode (one rank lookup per slice), which is what a locate would do
        ulong s = 0, x = seed | 1;
        long entries = 0;
        while (entries < ops)
        {
            x = Rng.Mix(x);
            long start = (long)(((x & 0xFFFFFFFFUL) * (ulong)(n - m)) >> 32);
            s += r.SumRange(start, m);
            entries += m;
        }
        return s;
    }

    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    private static ulong Loop<TR>(TR r, string kind, long n, long ops, int m, ulong seed) where TR : struct, IReader
    {
        ulong s = 0;
        switch (kind)
        {
            case "chase":
            {
                // a chain of dependent reads; it restarts at a new random index every 1,024 steps, because the trajectory of a random function falls into a cycle of about sqrt(n) entries
                // (30,000 at n = 10^9), which then stays in the cache and would measure the cache and not the memory
                ulong x = seed | 1;
                long idx = (long)(Rng.Mix(x) % (ulong)n);
                for (long k = 0; k < ops; k++)
                {
                    ulong v = r.Get(idx);
                    s += v;
                    idx = Rng.Next(v, n);
                    if ((k & 1023) == 1023)
                    {
                        x = Rng.Mix(x);
                        idx = (long)(((x & 0xFFFFFFFFUL) * (ulong)n) >> 32);
                    }
                }
                return s;
            }
            case "thr":
            {
                ulong x = seed | 1;
                for (long k = 0; k < ops; k++)
                {
                    x ^= x << 13; x ^= x >> 7; x ^= x << 17;
                    long idx = (long)(((x & 0xFFFFFFFFUL) * (ulong)n) >> 32);
                    s += r.Get(idx);
                }
                return s;
            }
            case "scan":
            {
                ulong x = seed | 1;
                long entries = 0;
                while (entries < ops)
                {
                    x ^= x << 13; x ^= x >> 7; x ^= x << 17;
                    long start = (long)(((x & 0xFFFFFFFFUL) * (ulong)(n - m)) >> 32);
                    for (int i = 0; i < m; i++) s += r.Get(start + i);
                    entries += m;
                }
                return s;
            }
            default: throw new ArgumentException(kind);
        }
    }

    /// <summary>The null experiment: the loops with no read (the index generation alone) on the layout's cheapest read, a cache-resident array.</summary>
    internal static List<Row> Null(string layout, int reps, long ops, string tag)
    {
        var rows = new List<Row>();
        using Packed a = Packed.Create(layout, 4096, 6_200_000_000UL, 0x77UL, 1);
        foreach (string kind in new[] { "chase", "thr", "scan" })
        {
            var samples = new List<double>();
            int m = kind == "scan" ? 1500 : 0;
            for (int r = 0; r < reps; r++)
            {
                var sw = Stopwatch.StartNew();
                ulong s = Dispatch(layout, a, kind, 4096, ops, m, (ulong)(r + 1));
                sw.Stop();
                GC.KeepAlive(s);
                samples.Add(sw.Elapsed.TotalMilliseconds * 1e6 / ops);
            }
            samples.Sort();
            rows.Add(new Row("null", layout, 4096, 1, kind == "scan" ? "scan1500_ns_per_entry" : kind + "_ns_per_read", samples[samples.Count / 2], samples[0], samples[^1], "ns/op", reps, tag + " cache-resident 4096 entries"));
        }
        return rows;
    }

    /// <summary>The loop with no memory read at all (the generator and the arithmetic of the three loops).</summary>
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    internal static List<Row> Empty(int reps, long ops, string tag)
    {
        var rows = new List<Row>();
        var samples = new List<double>();
        for (int r = 0; r < reps; r++)
        {
            var sw = Stopwatch.StartNew();
            ulong x = (ulong)(r + 1) | 1, s = 0;
            long n = 1_000_000_000;
            for (long k = 0; k < ops; k++)
            {
                x ^= x << 13; x ^= x >> 7; x ^= x << 17;
                s += ((x & 0xFFFFFFFFUL) * (ulong)n) >> 32;
            }
            sw.Stop();
            GC.KeepAlive(s);
            samples.Add(sw.Elapsed.TotalMilliseconds * 1e6 / ops);
        }
        samples.Sort();
        rows.Add(new Row("null", "none", 0, 1, "generator_ns_per_op", samples[samples.Count / 2], samples[0], samples[^1], "ns/op", reps, tag + " no read"));
        return rows;
    }
}
