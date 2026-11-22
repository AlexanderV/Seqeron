using System.Diagnostics;
using SuffixTree.Esa;

// The cost of the registration of a call on a tree over a container, in nanoseconds per call (work package 1b.9g).
//   FacadeCost facade <container> <text file>   IsEmpty (the registration, one comparison, the exit) on a tree over the container and on an in-memory tree of the first 100,000 symbols of the text,
//                                               1 and 4 threads on the one tree; the difference between the two trees is the registration
//   FacadeCost gate                             (compile with -p:Gate=true against a build of 1b.9g or later) the pieces of the new protocol on their own: the cell of the thread, the store and the clearing of an entry, a whole
//                                               Enter + Leave of a gate, and the striped count of 1b.8 emulated (Interlocked.Increment, read of the flag, Interlocked.Decrement, read of the flag)
// Each figure is the median of 7 (facade) or 5 (gate) runs of 10^8 (5 x 10^7 for the count) calls per thread; the range is the minimum and the maximum.
internal static class M
{
    private static long s_sink;

    private static double Median(List<double> values)
    {
        values.Sort();
        return values[values.Count / 2];
    }

    private static double Run(int threads, Func<long, long> body, long iterations, int repeats)
    {
        var results = new List<double>();
        for (int rep = 0; rep < repeats; rep++)
        {
            using var go = new ManualResetEventSlim();
            var nanos = new double[threads];
            List<Thread> workers = Enumerable.Range(0, threads).Select(k => new Thread(() =>
            {
                go.Wait();
                long t0 = Stopwatch.GetTimestamp();
                long s = body(iterations);
                nanos[k] = Stopwatch.GetElapsedTime(t0).TotalSeconds * 1e9 / iterations;
                Interlocked.Add(ref s_sink, s);
            })).ToList();
            workers.ForEach(t => t.Start());
            Thread.Sleep(30);
            go.Set();
            workers.ForEach(t => t.Join());
            results.Add(nanos.Average());
        }

        double median = Median(results);
        Console.Write($"{median,7:F2} ns/call (min {results.Min():F2} max {results.Max():F2})");
        return median;
    }

    private static int Main(string[] args)
    {
        if (args.Length >= 3 && args[0] == "facade")
        {
            EsaSuffixTree file = EsaSuffixTreeFactory.Open(args[1]);
            EsaSuffixTree memory = EsaSuffixTreeFactory.CreateEsa(File.ReadAllText(args[2])[..100_000]);
            foreach ((string name, EsaSuffixTree tree) in new[] { ("file-backed", file), ("in-memory", memory) })
                foreach (int threads in new[] { 1, 4 })
                {
                    Console.Write($"facade {name,-12} threads {threads}  IsEmpty ");
                    Run(threads, n =>
                    {
                        long s = 0;
                        for (long i = 0; i < n; i++)
                            if (tree.IsEmpty) s++;
                        return s;
                    }, 100_000_000, 7);
                    Console.WriteLine();
                }

            return 0;
        }

#if GATE
        if (args.Length >= 1 && args[0] == "gate")
        {
            var gate = new ReaderGate<ReaderGateSeams>(new object(), () => { });
            int flag = 0;
            var counts = new int[8 * 16];
            foreach (int threads in new[] { 1, 4 })
            {
                Console.Write($"gate threads {threads}  volatile read of a static      ");
                Run(threads, n => { long s = 0; for (long i = 0; i < n; i++) s += Volatile.Read(ref flag); return s; }, 200_000_000, 5);
                Console.WriteLine();
                Console.Write($"gate threads {threads}  ReaderCell.Current             ");
                Run(threads, n => { long s = 0; for (long i = 0; i < n; i++) s += ReaderCell.Current.Depth; return s; }, 200_000_000, 5);
                Console.WriteLine();
                Console.Write($"gate threads {threads}  cell.Publish + Retract         ");
                Run(threads, n => { ReaderCell c = ReaderCell.Current; for (long i = 0; i < n; i++) { c.Publish(7); c.Retract(); } return 0; }, 200_000_000, 5);
                Console.WriteLine();
                Console.Write($"gate threads {threads}  gate.Enter + Leave (1b.9g)     ");
                Run(threads, n => { for (long i = 0; i < n; i++) gate.Leave(gate.Enter()); return 0; }, 200_000_000, 5);
                Console.WriteLine();
                Console.Write($"gate threads {threads}  striped count of 1b.8 (emulated) ");
                Run(threads, n =>
                {
                    int stripe = Environment.CurrentManagedThreadId & 7;
                    for (long i = 0; i < n; i++)
                    {
                        Interlocked.Increment(ref counts[stripe * 16]);
                        if (Volatile.Read(ref flag) != 0) return 1;
                        Interlocked.Decrement(ref counts[stripe * 16]);
                        if (Volatile.Read(ref flag) != 0) return 1;
                    }

                    return 0;
                }, 50_000_000, 5);
                Console.WriteLine();
            }

            return 0;
        }
#endif

        Console.Error.WriteLine("usage: FacadeCost facade <container> <text file> | FacadeCost gate");
        return 2;
    }
}
