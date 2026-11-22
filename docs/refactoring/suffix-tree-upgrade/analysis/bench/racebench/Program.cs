using System.Collections.Concurrent;
using System.Diagnostics;
using System.Globalization;
using SuffixTree;
using SuffixTree.Esa;
using SuffixTree.Indexing;

// The Dispose-race client of the suffix-tree upgrade (S1b, work package 1b.9g): readers that run calls on a tree while the tree is disposed after a random delay, one process = many rounds.
//   racebench mem  <rounds> <seed>                              an in-memory tree: four readers walk suffix links (longest common substring, maximal exact matches, ranks), a thread forces compacting collections and refills
//                                                                 the freed memory, the main thread disposes after a random 0-2 ms (the race of 1b.9b: the pins of the navigation pack)
//   racebench file <rounds> <seed> <container> <text file>      a tree over a container: four readers run every kind of member (count, locate, contains, the longest common substring, the maximal exact matches, a traversal,
//                                                                 the suffixes, the suffix array, the text), the main thread disposes after a random 0-40 ms (the race of 1b.8: nothing may read the unmapped view)
// Every answer is compared with the answer of an in-memory tree made here; a call must complete with the right answer or throw ObjectDisposedException, anything else is a defect (a wrong answer, another exception, a crash of the
// process, a hang: the driver reports an exit code that is not 0 or a timeout). Output: one line, "RESULT mode rounds threads completed disposed wrong other"; exit code 0 when clean, 1 for a wrong answer or another exception.
internal static class Program
{
    private const int Readers = 4;

    private static int Main(string[] args)
    {
        if (args.Length < 3) return Usage();
        string mode = args[0];
        int rounds = int.Parse(args[1], CultureInfo.InvariantCulture);
        ulong seed = ulong.Parse(args[2], CultureInfo.InvariantCulture);
        var rng = new Random(unchecked((int)seed));
        long completed = 0, disposed = 0, wrong = 0;
        var others = new ConcurrentBag<string>();
        if (mode == "mem") RunMemory(rounds, rng, ref completed, ref disposed, ref wrong, others);
        else if (mode == "file" && args.Length >= 5) RunFile(rounds, rng, args[3], File.ReadAllText(args[4], System.Text.Encoding.Latin1), ref completed, ref disposed, ref wrong, others);
        else return Usage();
        Console.WriteLine($"RESULT mode={mode} rounds={rounds} threads={Readers} completed={completed} disposed={disposed} wrong={wrong} other={others.Count}");
        foreach (string o in others.Take(3)) Console.WriteLine("OTHER " + o);
        return wrong == 0 && others.IsEmpty && disposed == (long)rounds * Readers ? 0 : 1;
    }

    private static int Usage()
    {
        Console.Error.WriteLine("usage: racebench mem <rounds> <seed> | racebench file <rounds> <seed> <container> <text file>");
        return 2;
    }

    private static string RandomDna(int n, Random rng)
    {
        const string letters = "ACGT";
        return string.Create(n, rng, (span, r) =>
        {
            for (int i = 0; i < span.Length; i++) span[i] = letters[r.Next(4)];
        });
    }

    private static long TraverseCount(ISuffixTree tree)
    {
        var visitor = new Counting();
        tree.Traverse(visitor);
        return visitor.Nodes;
    }

    private sealed class Counting : ISuffixTreeVisitor
    {
        internal long Nodes;

        public void VisitNode(int startIndex, int endIndex, int leafCount, int childCount, int depth) => Nodes++;

        public void EnterBranch(int key)
        {
        }

        public void ExitBranch()
        {
        }
    }

    private static void Junk(int arrays)
    {
        // allocations that reuse the memory that a compacting collection freed, filled with a value that is no suffix array entry of any tree here
        var keep = new int[arrays][];
        for (int i = 0; i < keep.Length; i++)
        {
            keep[i] = new int[1500 + (i % 7 * 100)];
            Array.Fill(keep[i], int.MaxValue - 7);
        }

        GC.KeepAlive(keep);
    }

    private static void SpinFor(int microseconds)
    {
        long end = Stopwatch.GetTimestamp() + (microseconds * Stopwatch.Frequency / 1_000_000);
        while (Stopwatch.GetTimestamp() < end) Thread.SpinWait(20);
    }

    private static void RunMemory(int rounds, Random rng, ref long completedTotal, ref long disposedTotal, ref long wrongTotal, ConcurrentBag<string> others)
    {
        long completed = 0, disposed = 0, wrong = 0;
        for (int round = 0; round < rounds; round++)
        {
            string text = RandomDna(1_500 + (round * 13), rng);
            string[] patterns = Enumerable.Range(0, 32).Select(i => text.Substring((i * 97) % (text.Length - 30), 4 + (i % 9))).ToArray();
            string query = string.Concat(text.AsSpan(200, 300), "TTTTT", text.AsSpan(900, 200));
            using EsaSuffixTree reference = EsaSuffixTree.Build(text);
            int[] counts = patterns.Select(reference.CountOccurrences).ToArray();
            string lcs = reference.LongestCommonSubstring(query);
            var mems = reference.FindMaximalExactMatches(query, 8);
            long nodes = TraverseCount(reference);
            Junk(100);
            EsaSuffixTree tree = EsaSuffixTree.Build(text);
            _ = tree.LongestCommonSubstring(query); // the pack exists: the readers walk suffix links
            using var stop = new ManualResetEventSlim();
            var collector = new Thread(() =>
            {
                while (!stop.IsSet)
                {
                    Junk(40);
                    GC.Collect(GC.MaxGeneration, GCCollectionMode.Forced, blocking: true, compacting: true);
                }
            });
            using var start = new Barrier(Readers + 1);
            var threads = Enumerable.Range(0, Readers).Select(t => new Thread(() =>
            {
                start.SignalAndWait();
                long k = t;
                while (true)
                {
                    int i = (int)(k % patterns.Length);
                    try
                    {
                        bool ok = (k % 6) switch
                        {
                            0 => tree.CountOccurrences(patterns[i]) == counts[i],
                            1 => tree.LongestCommonSubstring(query) == lcs,
                            2 => tree.FindMaximalExactMatches(query, 8).SequenceEqual(mems),
                            3 => TraverseCount(tree) == nodes,
                            4 => ((ISuffixArrayAccess)tree).RankOf(1 + (k % text.Length)) >= 0,
                            _ => tree.Contains(patterns[i]),
                        };
                        if (ok) Interlocked.Increment(ref completed);
                        else Interlocked.Increment(ref wrong);
                    }
                    catch (ObjectDisposedException)
                    {
                        Interlocked.Increment(ref disposed);
                        return;
                    }
                    catch (Exception e)
                    {
                        others.Add($"round {round}: {e.GetType().Name}: {e.Message}");
                        return;
                    }

                    k += Readers;
                }
            })).ToList();
            collector.Start();
            foreach (Thread thread in threads) thread.Start();
            start.SignalAndWait();
            SpinFor(rng.Next(0, 7) * 300);
            tree.Dispose();
            SpinFor(2_000);
            foreach (Thread thread in threads)
                if (!thread.Join(TimeSpan.FromSeconds(120))) Environment.FailFast("a reader did not finish after Dispose (a hang)");
            stop.Set();
            collector.Join();
        }

        completedTotal += completed;
        disposedTotal += disposed;
        wrongTotal += wrong;
    }

    private static void RunFile(int rounds, Random rng, string container, string text, ref long completedTotal, ref long disposedTotal, ref long wrongTotal, ConcurrentBag<string> others)
    {
        long completed = 0, disposed = 0, wrong = 0;
        string[] patterns = Enumerable.Range(0, 64).Select(i => text.Substring((i * 613) % (text.Length - 30), 6 + (i % 17))).ToArray();
        string query = string.Concat(text.AsSpan(5_000, 600), "TTTTT", text.AsSpan(9_000, 400));
        using EsaSuffixTree memory = EsaSuffixTree.Build(text);
        int[] counts = patterns.Select(memory.CountOccurrences).ToArray();
        int[][] positions = patterns.Select(p => memory.FindAllOccurrences(p).Order().ToArray()).ToArray();
        var mems = memory.FindMaximalExactMatches(query, 12);
        string lcs = memory.LongestCommonSubstring(query);
        long nodes = TraverseCount(memory);
        for (int round = 0; round < rounds; round++)
        {
            EsaSuffixTree tree = EsaSuffixTreeFactory.Open(container, new EsaOpenOptions
            {
                Validation = round % 3 == 2 ? EsaValidationLevel.Structure : EsaValidationLevel.Header,
                Access = round % 2 == 0 ? EsaAccessMode.Checked : EsaAccessMode.Trusted,
            });
            using var start = new Barrier(Readers + 1);
            var threads = Enumerable.Range(0, Readers).Select(t => new Thread(() =>
            {
                start.SignalAndWait();
                long k = t;
                while (true)
                {
                    int i = (int)(k % patterns.Length);
                    try
                    {
                        bool ok = (k % 9) switch
                        {
                            0 => tree.CountOccurrences(patterns[i]) == counts[i],
                            1 => tree.FindAllOccurrences(patterns[i]).Order().SequenceEqual(positions[i]),
                            2 => tree.Contains(patterns[i]),
                            3 => tree.LongestCommonSubstring(query) == lcs,
                            4 => tree.FindMaximalExactMatches(query, 12).SequenceEqual(mems),
                            5 => TraverseCount(tree) == nodes,
                            6 => tree.EnumerateSuffixes().Take(700).Count() == 700,
                            7 => ((ISuffixArrayAccess)tree).SuffixAt(1 + (k % text.Length)) >= 0,
                            _ => tree.Text[(int)(k % text.Length)] == text[(int)(k % text.Length)],
                        };
                        if (ok) Interlocked.Increment(ref completed);
                        else Interlocked.Increment(ref wrong);
                    }
                    catch (ObjectDisposedException)
                    {
                        Interlocked.Increment(ref disposed);
                        return;
                    }
                    catch (Exception e)
                    {
                        others.Add($"round {round}: {e.GetType().Name}: {e.Message}");
                        return;
                    }

                    k += Readers;
                }
            })).ToList();
            foreach (Thread thread in threads) thread.Start();
            start.SignalAndWait();
            SpinFor(rng.Next(0, 40_000));
            tree.Dispose();
            foreach (Thread thread in threads)
                if (!thread.Join(TimeSpan.FromSeconds(120))) Environment.FailFast("a reader did not finish after Dispose (a hang)");
        }

        completedTotal += completed;
        disposedTotal += disposed;
        wrongTotal += wrong;
    }
}
