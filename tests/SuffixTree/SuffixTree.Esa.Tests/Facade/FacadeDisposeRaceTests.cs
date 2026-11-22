using System.Collections.Concurrent;
using System.Diagnostics;
using System.Globalization;
using System.Text.RegularExpressions;
using SuffixTree.Esa.Tests.External;
using SuffixTree.Esa.Tests.TestSupport;
using SuffixTree.Indexing;

namespace SuffixTree.Esa.Tests.Facade;

/// <summary>
/// A read that races <c>Dispose</c> of a file-backed tree (work package 1b.8, the open point of 1b.7f: "Dispose during a read is not safe by design and not tested"). The contract that is tested: a call that is in
/// flight when <c>Dispose</c> runs either completes with the right answer or throws <see cref="ObjectDisposedException"/>, and no call touches the mapping after it was unmapped. The mapping is the one thing in this
/// process that can kill it (an access violation is not catchable), so the race runs in a child process: this test starts <c>dotnet test</c> on the explicit class <see cref="FacadeDisposeRaceChildTests"/> of the same
/// assembly and reads its verdict and its result line; a crash of the child is a failure of this test with the evidence attached.
/// </summary>
[TestFixture]
public class FacadeDisposeRaceTests
{
    internal const string ChildCategory = "DisposeRaceChild";
    internal const string ResultFileVariable = "ESA_DISPOSE_RACE_OUT";
    private const int ChildTests = 7;

    [Test]
    public void ReadsRacingDispose_OfAFileBackedTree_CompleteOrThrowObjectDisposed_AndNeverTouchTheUnmappedView()
    {
        string dll = typeof(FacadeDisposeRaceTests).Assembly.Location;
        using var dir = new TempDir();
        string result = dir.File("result.txt");
        var start = new ProcessStartInfo(Environment.ProcessPath ?? "dotnet")
        {
            RedirectStandardOutput = true,
            RedirectStandardError = true,
            UseShellExecute = false,
        };
        foreach (string arg in new[] { "test", dll, "--nologo", "--filter", "Category=" + ChildCategory }) start.ArgumentList.Add(arg);
        start.Environment[ResultFileVariable] = result;

        using Process child = Process.Start(start) ?? throw new InvalidOperationException("the child process did not start");
        Task<string> stdout = child.StandardOutput.ReadToEndAsync();
        Task<string> stderr = child.StandardError.ReadToEndAsync();
        if (!child.WaitForExit(TimeSpan.FromMinutes(6)))
        {
            child.Kill(entireProcessTree: true);
            Assert.Fail("the child process of the Dispose race did not finish within 6 minutes (a hang inside a read or inside Dispose)");
        }
        string output = stdout.Result + stderr.Result;
        string lines = File.Exists(result) ? File.ReadAllText(result) : "(no result file: the child died before it wrote one)";
        TestContext.Out.WriteLine(lines);
        Assert.That(output, Does.Not.Contain("crashed"), $"the test host of the race crashed (a read touched the unmapped view?)\n{lines}\n{Tail(output)}");
        Match passed = Regex.Match(output, @"Passed:\s+(\d+)");
        Match failed = Regex.Match(output, @"Failed:\s+(\d+)");
        Assert.That(passed.Success && failed.Success, Is.True, $"no verdict from the child\n{Tail(output)}");
        Assert.That(int.Parse(failed.Groups[1].Value, CultureInfo.InvariantCulture), Is.Zero, $"a test of the race failed\n{lines}\n{Tail(output)}");
        Assert.That(int.Parse(passed.Groups[1].Value, CultureInfo.InvariantCulture), Is.EqualTo(ChildTests), $"the child ran {passed.Groups[1].Value} of {ChildTests} tests\n{Tail(output)}");
    }

    private static string Tail(string text) => text.Length <= 3000 ? text : text[^3000..];
}

/// <summary>
/// The tests that run in the child process of <see cref="FacadeDisposeRaceTests"/> (explicit: the default run never starts this class itself, a failure or a crash here is the finding). Every test but one opens a file-backed tree,
/// lets threads read it, and disposes it while they do; the other one (work package 1b.9b) does the same with an in-memory tree while a collector thread moves memory.
/// </summary>
[TestFixture]
[Category(FacadeDisposeRaceTests.ChildCategory)]
[Explicit("child process of FacadeDisposeRaceTests; run with --filter Category=DisposeRaceChild (a crash of the test host is the finding it detects)")]
public class FacadeDisposeRaceChildTests
{
    private const int ReaderThreads = 4;
    private const int Rounds = 80;

    private sealed class CountingVisitor : ISuffixTreeVisitor
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

    private sealed class DisposingVisitor(IDisposable tree) : ISuffixTreeVisitor
    {
        internal long Nodes;

        public void VisitNode(int startIndex, int endIndex, int leafCount, int childCount, int depth)
        {
            if (Nodes++ == 3) tree.Dispose();
        }

        public void EnterBranch(int key)
        {
        }

        public void ExitBranch()
        {
        }
    }

    private sealed class BlockingVisitor(ManualResetEventSlim inside, ManualResetEventSlim go) : ISuffixTreeVisitor
    {
        internal long Nodes;

        public void VisitNode(int startIndex, int endIndex, int leafCount, int childCount, int depth)
        {
            if (Nodes++ == 3)
            {
                inside.Set();
                Assert.That(go.Wait(TimeSpan.FromSeconds(60)), Is.True, "the reader was never released");
            }
        }

        public void EnterBranch(int key)
        {
        }

        public void ExitBranch()
        {
        }
    }

    private static void Report(string line)
    {
        TestContext.Out.WriteLine(line);
        string? file = Environment.GetEnvironmentVariable(FacadeDisposeRaceTests.ResultFileVariable);
        if (file is not null) File.AppendAllText(file, line + "\n");
    }

    private static EsaSuffixTree OpenFile(string path, int round) => SesaTestFile.Open(
        path,
        round % 3 == 2 ? EsaValidationLevel.Structure : EsaValidationLevel.Header,
        round % 2 == 0 ? EsaAccessMode.Checked : EsaAccessMode.Trusted);

    [Test]
    public void ReaderThreads_RacingDispose_EitherCompleteTheCall_OrThrowObjectDisposed()
    {
        string text = StringFamilies.RandomDna(40_000, 3) + StringFamilies.RandomDna(2_000, 9);
        string[] patterns = Enumerable.Range(0, 64).Select(i => text.Substring((i * 613) % (text.Length - 30), 6 + (i % 17))).ToArray();
        string query = string.Concat(text.AsSpan(5_000, 600), "TTTTT", text.AsSpan(9_000, 400));
        using EsaSuffixTree memory = EsaSuffixTree.Build(text);
        int[] counts = patterns.Select(memory.CountOccurrences).ToArray();
        int[][] positions = patterns.Select(p => memory.FindAllOccurrences(p).Order().ToArray()).ToArray();
        var mems = memory.FindMaximalExactMatches(query, 12);
        string lcs = memory.LongestCommonSubstring(query);
        long nodes = 0;
        var counting = new CountingVisitor();
        memory.Traverse(counting);
        nodes = counting.Nodes;
        using var dir = new TempDir();
        string path = SesaTestFile.Write(dir, text);

        long completed = 0, disposed = 0, wrong = 0;
        var others = new ConcurrentBag<string>();
        for (int round = 0; round < Rounds; round++)
        {
            EsaSuffixTree tree = OpenFile(path, round);
            var start = new Barrier(ReaderThreads + 1);
            var threads = Enumerable.Range(0, ReaderThreads).Select(t => new Thread(() =>
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
                    k += ReaderThreads;
                }
            })).ToList();
            foreach (Thread thread in threads) thread.Start();
            start.SignalAndWait();
            SpinFor(round % 7 * 400);
            tree.Dispose();
            foreach (Thread thread in threads) Assert.That(thread.Join(TimeSpan.FromSeconds(120)), Is.True, "a reader did not finish after Dispose (a hang)");
        }

        Report($"RESULT race rounds={Rounds} threads={ReaderThreads} completed={completed} disposed={disposed} wrong={wrong} other={others.Count}");
        Assert.Multiple(() =>
        {
            Assert.That(wrong, Is.Zero, "a call that completed returned a wrong answer (a read of freed memory?)");
            Assert.That(others, Is.Empty, "a call failed with something other than ObjectDisposedException: " + string.Join(" | ", others.Take(3)));
            Assert.That(disposed, Is.EqualTo(Rounds * ReaderThreads), "every reader ends at the disposal");
        });
    }

    private static long TraverseCount(ISuffixTree tree)
    {
        var visitor = new CountingVisitor();
        tree.Traverse(visitor);
        return visitor.Nodes;
    }

    private static void SpinFor(int microseconds)
    {
        long end = Stopwatch.GetTimestamp() + (microseconds * Stopwatch.Frequency / 1_000_000);
        while (Stopwatch.GetTimestamp() < end) Thread.SpinWait(20);
    }

    [Test]
    public void ADisposeThatArrivesInsideAReadOfAnotherThread_DoesNotUnmapUnderIt()
    {
        string text = StringFamilies.RandomDna(20_000, 5);
        using var dir = new TempDir();
        string path = SesaTestFile.Write(dir, text);
        EsaSuffixTree tree = OpenFile(path, 0);
        using EsaSuffixTree memory = EsaSuffixTree.Build(text);
        long expected = TraverseCount(memory);
        using var inside = new ManualResetEventSlim();
        using var go = new ManualResetEventSlim();
        var visitor = new BlockingVisitor(inside, go);
        Task reader = Task.Run(() => tree.Traverse(visitor));
        Assert.That(inside.Wait(TimeSpan.FromSeconds(60)), Is.True);
        Task disposer = Task.Run(tree.Dispose);
        Assert.That(reader.Wait(TimeSpan.FromMilliseconds(300)), Is.False, "the read is blocked inside the visitor");
        go.Set();
        Assert.That(reader.Wait(TimeSpan.FromSeconds(60)), Is.True, "the read that was in flight did not finish");
        Assert.That(disposer.Wait(TimeSpan.FromSeconds(60)), Is.True, "Dispose did not return once the read had finished");
        Report($"RESULT in-flight nodes={visitor.Nodes} expected={expected}");
        Assert.That(visitor.Nodes, Is.EqualTo(expected), "the read in flight finished with every node");
        Assert.That(() => tree.CountOccurrences("A"), Throws.TypeOf<ObjectDisposedException>());
    }

    [Test]
    public void ADisposeFromInsideAVisitor_DoesNotDeadlock_AndTheTraversalFinishes()
    {
        string text = StringFamilies.RandomDna(8_000, 11);
        using var dir = new TempDir();
        string path = SesaTestFile.Write(dir, text);
        EsaSuffixTree tree = OpenFile(path, 1);
        using EsaSuffixTree memory = EsaSuffixTree.Build(text);
        long expected = TraverseCount(memory);
        var visitor = new DisposingVisitor(tree);
        Task traversal = Task.Run(() => tree.Traverse(visitor));
        Assert.That(traversal.Wait(TimeSpan.FromSeconds(60)), Is.True, "a Dispose from inside the traversal never returned");
        Report($"RESULT dispose-in-visitor nodes={visitor.Nodes} expected={expected}");
        Assert.That(visitor.Nodes, Is.EqualTo(expected));
        Assert.That(() => tree.Contains("A"), Throws.TypeOf<ObjectDisposedException>());
        tree.Dispose();
    }

    [Test]
    public void TheLazyResults_OfADisposedTree_ThrowObjectDisposed_NeverReadTheUnmappedView()
    {
        string text = StringFamilies.RandomDna(5_000, 13);
        using var dir = new TempDir();
        string path = SesaTestFile.Write(dir, text);
        EsaSuffixTree tree = OpenFile(path, 0);
        ITextSource source = tree.Text;
        IIndexText indexText = ((ISuffixIndex)tree).Text;
        IEnumerator<string> suffixes = tree.EnumerateSuffixes().GetEnumerator();
        Assert.That(suffixes.MoveNext(), Is.True);
        Assert.That(source[7], Is.EqualTo(text[7]));
        Assert.That(indexText.SymbolAt(7), Is.EqualTo(text[7]));
        tree.Dispose();
        Assert.Multiple(() =>
        {
            Assert.That(() => suffixes.MoveNext(), Throws.TypeOf<ObjectDisposedException>(), "the enumerator");
            Assert.That(() => source[7], Throws.TypeOf<ObjectDisposedException>(), "ITextSource indexer");
            Assert.That(() => source.Substring(3, 5), Throws.TypeOf<ObjectDisposedException>(), "ITextSource.Substring");
            Assert.That(() => source.Slice(3, 5).ToString(), Throws.TypeOf<ObjectDisposedException>(), "ITextSource.Slice");
            Assert.That(() => source.Length, Throws.Nothing, "the length is a value");
            Assert.That(() => indexText.SymbolAt(7), Throws.TypeOf<ObjectDisposedException>(), "IIndexText.SymbolAt");
            Assert.That(() => indexText.CopyTo(0, new char[4]), Throws.TypeOf<ObjectDisposedException>(), "IIndexText.CopyTo");
        });
        Report("RESULT lazy results after Dispose: ObjectDisposedException");
    }

    // ---------------------------------------------------------------- in-memory trees (work package 1b.9b)

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

    /// <summary>
    /// The question of 1b.9b: does <c>Dispose</c> of an in-memory tree release anything that a reader which is still inside a call uses? Yes, before this package: the navigation pack pins the suffix array and the LCP base array
    /// with GC handles and reads them through pointers, and <c>Dispose</c> freed the handles, so a compacting collection could move an array while a reader held its old address; the reader then took garbage for a position and used it
    /// as an index into the other pinned arrays, which ended the process (3 of 3 runs of this test against the behaviour of <c>e4a0b2da</c>: the test host crashed). The collector thread forces compacting collections and fills the
    /// freed memory while the readers walk suffix links (the longest common substring, the maximal exact matches, the ranks) and the main thread disposes the tree. A call must complete with the right answer or throw
    /// <see cref="ObjectDisposedException"/>.
    /// </summary>
    [Test]
    public void ReaderThreads_RacingDispose_OfAnInMemoryTree_CompleteOrThrowObjectDisposed_WhileTheCollectorMovesMemory()
    {
        const int rounds = 120;
        long completed = 0, disposed = 0, wrong = 0;
        var others = new ConcurrentBag<string>();
        for (int round = 0; round < rounds; round++)
        {
            string text = StringFamilies.RandomDna(1_500 + (round * 13), (ulong)round + 100);
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
            var start = new Barrier(ReaderThreads + 1);
            var threads = Enumerable.Range(0, ReaderThreads).Select(t => new Thread(() =>
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
                    k += ReaderThreads;
                }
            })).ToList();
            collector.Start();
            foreach (Thread thread in threads) thread.Start();
            start.SignalAndWait();
            SpinFor(round % 7 * 300);
            tree.Dispose();
            SpinFor(2_000);
            foreach (Thread thread in threads) Assert.That(thread.Join(TimeSpan.FromSeconds(120)), Is.True, "a reader did not finish after Dispose (a hang)");
            stop.Set();
            collector.Join();
        }

        Report($"RESULT in-memory race rounds={rounds} threads={ReaderThreads} completed={completed} disposed={disposed} wrong={wrong} other={others.Count}");
        Assert.Multiple(() =>
        {
            Assert.That(wrong, Is.Zero, "a call that completed returned a wrong answer (a read of memory that the collector moved?)");
            Assert.That(others, Is.Empty, "a call failed with something other than ObjectDisposedException: " + string.Join(" | ", others.Take(3)));
            Assert.That(disposed, Is.EqualTo(rounds * ReaderThreads), "every reader ends at the disposal");
        });
    }

    /// <summary>
    /// The registry of the reader cells (work package 1b.9g) under a race: threads that start, register a cell (or take over the cell of a thread that has just ended), call the tree a few times and end, while the tree is disposed. The
    /// releaser scans the list while it grows and while cells change owner; a call must complete with the right answer or throw <see cref="ObjectDisposedException"/>, and nothing may read the unmapped view.
    /// </summary>
    [Test]
    public void ThreadsThatStartAndEnd_WhileTheTreeIsDisposed_RegisterAndTakeOverCellsRacingTheScan_AndNoCallReadsTheUnmappedView()
    {
        string text = StringFamilies.RandomDna(30_000, 21);
        string[] patterns = Enumerable.Range(0, 32).Select(i => text.Substring((i * 311) % (text.Length - 30), 6 + (i % 11))).ToArray();
        using EsaSuffixTree memory = EsaSuffixTree.Build(text);
        int[] counts = patterns.Select(memory.CountOccurrences).ToArray();
        using var dir = new TempDir();
        string path = SesaTestFile.Write(dir, text);
        long completed = 0, disposed = 0, wrong = 0, threadsStarted = 0;
        var others = new ConcurrentBag<string>();
        const int rounds = 60;
        for (int round = 0; round < rounds; round++)
        {
            EsaSuffixTree tree = OpenFile(path, round);
            using var stop = new ManualResetEventSlim();
            void Reader(int seed)
            {
                try
                {
                    for (int k = 0; k < 3; k++)
                    {
                        int i = (seed + k) % patterns.Length;
                        if (tree.CountOccurrences(patterns[i]) == counts[i]) Interlocked.Increment(ref completed);
                        else Interlocked.Increment(ref wrong);
                    }
                }
                catch (ObjectDisposedException)
                {
                    Interlocked.Increment(ref disposed);
                }
                catch (Exception e)
                {
                    others.Add($"round {round}: {e.GetType().Name}: {e.Message}");
                }
            }

            var spawners = Enumerable.Range(0, 3).Select(s => new Thread(() =>
            {
                int seed = s * 7;
                while (!stop.IsSet)
                {
                    var t = new Thread(() => Reader(seed++));
                    t.Start();
                    Interlocked.Increment(ref threadsStarted);
                    t.Join();
                }
            })).ToList();
            foreach (Thread spawner in spawners) spawner.Start();
            SpinFor(round % 9 * 500);
            tree.Dispose();
            SpinFor(1_500);
            stop.Set();
            foreach (Thread spawner in spawners) Assert.That(spawner.Join(TimeSpan.FromSeconds(120)), Is.True, "a spawner did not finish (a hang)");
            Assert.That(() => tree.CountOccurrences(patterns[0]), Throws.TypeOf<ObjectDisposedException>(), $"round {round}");
        }

        Report($"RESULT short-lived readers rounds={rounds} threads={threadsStarted} completed={completed} disposed={disposed} wrong={wrong} other={others.Count}");
        Assert.Multiple(() =>
        {
            Assert.That(wrong, Is.Zero, "a call that completed returned a wrong answer");
            Assert.That(others, Is.Empty, "a call failed with something other than ObjectDisposedException: " + string.Join(" | ", others.Take(3)));
            Assert.That(completed, Is.GreaterThan(0));
        });
    }

    [Test]
    public void ManyTreesOfOneFile_DisposedWhileTheirReadersRun_LeaveTheOthersIntact()
    {
        string text = StringFamilies.RandomDna(10_000, 17);
        using var dir = new TempDir();
        string path = SesaTestFile.Write(dir, text);
        string pattern = text.Substring(100, 12);
        using EsaSuffixTree memory = EsaSuffixTree.Build(text);
        int expected = memory.CountOccurrences(pattern);
        var trees = Enumerable.Range(0, 6).Select(i => OpenFile(path, i)).ToArray();
        var failures = new ConcurrentBag<string>();
        var threads = trees.Select((tree, i) => new Thread(() =>
        {
            try
            {
                for (int k = 0; k < 3000; k++)
                    if (tree.CountOccurrences(pattern) != expected) failures.Add($"tree {i} gave a wrong count");
            }
            catch (ObjectDisposedException)
            {
                // the tree of this thread was disposed: the end of its reads
            }
            catch (Exception e)
            {
                failures.Add($"tree {i}: {e.GetType().Name}: {e.Message}");
            }
        })).ToList();
        foreach (Thread thread in threads) thread.Start();
        for (int i = 0; i < trees.Length; i += 2) trees[i].Dispose();
        foreach (Thread thread in threads) Assert.That(thread.Join(TimeSpan.FromSeconds(120)), Is.True);
        foreach (EsaSuffixTree tree in trees) tree.Dispose();
        Report($"RESULT many trees failures={failures.Count}");
        Assert.That(failures, Is.Empty);
    }
}
