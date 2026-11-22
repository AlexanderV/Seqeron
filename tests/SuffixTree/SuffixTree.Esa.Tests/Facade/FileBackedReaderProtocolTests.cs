using SuffixTree.Esa.External;
using SuffixTree.Esa.Tests.External;
using SuffixTree.Esa.Tests.TestSupport;
using SuffixTree.Indexing;

namespace SuffixTree.Esa.Tests.Facade;

/// <summary>
/// The reader protocol of 1b.9g at the level of the tree: the calls of a tree over a container are entries in the cell of the calling thread (the in-memory tree takes none), they nest, they leave with the call, a call after
/// <c>Dispose</c> leaves nothing behind, and trees that are many and threads that are short-lived add nothing to the registry. The lifetime contract itself (a call in flight completes, a later call throws, nothing reads the unmapped view) is
/// tested by <see cref="FacadeDisposeRaceTests"/>, <see cref="FacadeReaderCountTests"/> and <see cref="FacadeLifetimeTests"/>.
/// </summary>
[TestFixture]
public class FileBackedReaderProtocolTests
{
    private static readonly TimeSpan Wait = TimeSpan.FromSeconds(30);

    private sealed class DepthProbe : ISuffixTreeVisitor
    {
        internal readonly List<int> Depths = [];
        internal Action? OnNode;

        public void VisitNode(int startIndex, int endIndex, int leafCount, int childCount, int depth)
        {
            Depths.Add(ReaderCell.Current.Depth);
            OnNode?.Invoke();
        }

        public void EnterBranch(int key)
        {
        }

        public void ExitBranch()
        {
        }
    }

    private static string Container(TempDir dir, string text = "mississippi$banana$abracadabra") => SesaTestFile.Write(dir, text);

    [Test]
    public void ACallOfATreeOverAContainer_IsAnEntryInTheCellOfItsThread_AndAnInMemoryTreeTakesNone()
    {
        using var dir = new TempDir();
        string text = StringFamilies.RandomDna(500, 5);
        using EsaSuffixTree file = SesaTestFile.Open(Container(dir, text));
        using EsaSuffixTree memory = EsaSuffixTree.Build(text);
        var inFile = new DepthProbe();
        var inMemory = new DepthProbe();
        file.Traverse(inFile);
        memory.Traverse(inMemory);
        Assert.Multiple(() =>
        {
            Assert.That(inFile.Depths, Is.Not.Empty);
            Assert.That(inFile.Depths.Distinct(), Is.EqualTo(new[] { 1 }), "inside the traversal of a tree over a container this thread is inside one call");
            Assert.That(inMemory.Depths.Distinct(), Is.EqualTo(new[] { 0 }), "an in-memory tree takes no entry");
            Assert.That(file.ReaderId, Is.Not.Zero);
            Assert.That(memory.ReaderId, Is.Zero);
            Assert.That(ReaderCell.Current.Depth, Is.Zero, "nothing is left after the calls");
            Assert.That(ReaderCell.AnyInside(file.ReaderId), Is.False);
        });
    }

    [Test]
    public void TheCallsOfAGuardedResult_NestInsideTheCall_AndEveryMemberOfTheTreeCountsItself()
    {
        using var dir = new TempDir();
        string text = StringFamilies.RandomDna(300, 6);
        using EsaSuffixTree file = SesaTestFile.Open(Container(dir, text));
        ITextSource source = file.Text;
        var probe = new DepthProbe();
        // from inside a visitor: the traversal is call 1; the text read, the count and a second traversal nest
        probe.OnNode = () =>
        {
            if (probe.Depths.Count != 3) return;
            Assert.That(source[3], Is.EqualTo(text[3]));
            Assert.That(file.CountOccurrences("A"), Is.GreaterThan(0));
            var inner = new DepthProbe();
            file.Traverse(inner);
            Assert.That(inner.Depths.Distinct(), Is.EqualTo(new[] { 2 }), "a call inside a call is the second entry");
            Assert.That(ReaderCell.Current.Depth, Is.EqualTo(1));
        };
        file.Traverse(probe);
        Assert.That(ReaderCell.Current.Depth, Is.Zero);
    }

    [Test]
    public void TwoTreesOneInsideTheOther_AreSeparateEntries_AndLeaveInTheOrderOfTheCalls()
    {
        using var dir = new TempDir();
        using EsaSuffixTree a = SesaTestFile.Open(Container(dir, StringFamilies.RandomDna(200, 7)));
        using EsaSuffixTree b = SesaTestFile.Open(Container(dir, StringFamilies.RandomDna(200, 8)));
        var seen = new List<(bool InA, bool InB, int Depth)>();
        var outer = new DepthProbe();
        outer.OnNode = () =>
        {
            if (outer.Depths.Count != 2) return;
            var inner = new DepthProbe();
            inner.OnNode = () =>
            {
                if (inner.Depths.Count == 2) seen.Add((ReaderCell.AnyInside(a.ReaderId), ReaderCell.AnyInside(b.ReaderId), ReaderCell.Current.Depth));
            };
            b.Traverse(inner);
            seen.Add((ReaderCell.AnyInside(a.ReaderId), ReaderCell.AnyInside(b.ReaderId), ReaderCell.Current.Depth));
        };
        a.Traverse(outer);
        Assert.Multiple(() =>
        {
            Assert.That(seen[0], Is.EqualTo((true, true, 2)), "inside both");
            Assert.That(seen[1], Is.EqualTo((true, false, 1)), "the inner call left, the outer one is still inside");
            Assert.That(ReaderCell.AnyInside(a.ReaderId) || ReaderCell.AnyInside(b.ReaderId), Is.False);
        });
    }

    [Test]
    public void ACallAfterDispose_ThrowsObjectDisposed_AndLeavesNothingInTheCell_ForEveryKindOfMember()
    {
        using var dir = new TempDir();
        EsaSuffixTree tree = SesaTestFile.Open(Container(dir));
        ITextSource source = tree.Text;
        long id = tree.ReaderId;
        tree.Dispose();
        ReaderCell cell = ReaderCell.Current;
        Assert.Multiple(() =>
        {
            Assert.That(() => tree.CountOccurrences("ss"), Throws.TypeOf<ObjectDisposedException>());
            Assert.That(() => tree.Contains("ss"), Throws.TypeOf<ObjectDisposedException>());
            Assert.That(() => tree.FindAllOccurrences("ss"), Throws.TypeOf<ObjectDisposedException>());
            Assert.That(() => tree.NodeCount, Throws.TypeOf<ObjectDisposedException>());
            Assert.That(() => tree.Traverse(new DepthProbe()), Throws.TypeOf<ObjectDisposedException>());
            Assert.That(() => source[1], Throws.TypeOf<ObjectDisposedException>());
            Assert.That(cell.Depth, Is.Zero);
            Assert.That(ReaderCell.AnyInside(id), Is.False);
        });
    }

    [Test]
    public void ManyTreesOverContainers_OneThread_AddNothingToTheRegistry_AndEveryOneIsReleasedByItsDispose()
    {
        using var dir = new TempDir();
        string path = Container(dir);
        _ = ReaderCell.Current;
        int cells = ReaderCell.CellCount;
        var ids = new HashSet<long>();
        for (int i = 0; i < 150; i++)
        {
            EsaSuffixTree tree = SesaTestFile.Open(path, access: i % 2 == 0 ? EsaAccessMode.Checked : EsaAccessMode.Trusted);
            SesaFile file = tree.File!;
            Assert.That(ids.Add(tree.ReaderId), Is.True, "a gate id is never used twice");
            Assert.That(tree.CountOccurrences("ss"), Is.EqualTo(2));
            tree.Dispose();
            Assert.That(() => file.Payload(file.Parsed.Text), Throws.TypeOf<ObjectDisposedException>(), $"tree {i}: the mapping was released by the disposal");
            Assert.That(ReaderCell.AnyInside(tree.ReaderId), Is.False);
        }

        Assert.That(ReaderCell.CellCount, Is.EqualTo(cells), "150 trees read from one thread register nothing but the thread");
    }

    [Test]
    public void ShortLivedThreads_ThatReadSeveralTrees_DoNotGrowTheRegistry()
    {
        using var dir = new TempDir();
        string path = Container(dir);
        EsaSuffixTree[] trees = Enumerable.Range(0, 3).Select(_ => SesaTestFile.Open(path)).ToArray();
        _ = ReaderCell.Current;
        int before = ReaderCell.CellCount;
        int wrong = 0;
        for (int i = 0; i < 200; i++)
        {
            var thread = new Thread(() =>
            {
                try
                {
                    foreach (EsaSuffixTree tree in trees)
                        if (tree.CountOccurrences("ss") != 2) Interlocked.Increment(ref wrong);
                }
                catch (Exception)
                {
                    Interlocked.Increment(ref wrong);
                }
            }) { IsBackground = true };
            thread.Start();
            Assert.That(thread.Join(Wait), Is.True);
        }

        Assert.That(wrong, Is.Zero, "every call of every short-lived thread answered");

        int after = ReaderCell.CellCount;
        foreach (EsaSuffixTree tree in trees) tree.Dispose();
        Assert.That(after - before, Is.LessThanOrEqualTo(1), "200 short-lived threads one after the other are one cell");
    }

    [Test]
    public void SeveralReadersInsideAtOnce_HoldTheMappingUntilTheLastOneLeaves()
    {
        using var dir = new TempDir();
        EsaSuffixTree tree = SesaTestFile.Open(Container(dir, StringFamilies.RandomDna(400, 9)));
        SesaFile file = tree.File!;
        const int readers = 3;
        using var inside = new CountdownEvent(readers);
        var leave = Enumerable.Range(0, readers).Select(_ => new ManualResetEventSlim()).ToArray();
        var tasks = Enumerable.Range(0, readers).Select(i => Task.Factory.StartNew(() =>
        {
            tree.Traverse(new Blocking(inside, leave[i]));
        }, TaskCreationOptions.LongRunning)).ToArray();
        Assert.That(inside.Wait(Wait), Is.True);
        tree.Dispose();
        for (int i = 0; i < readers; i++)
        {
            Assert.That(() => file.Payload(file.Parsed.Text), Throws.Nothing, $"{readers - i} reader(s) still inside: the mapping stays");
            leave[i].Set();
            Assert.That(tasks[i].Wait(Wait), Is.True);
        }

        Assert.That(() => file.Payload(file.Parsed.Text), Throws.TypeOf<ObjectDisposedException>(), "the last reader released the mapping");
        foreach (ManualResetEventSlim e in leave) e.Dispose();
    }

    [Test]
    public void TheProcessWideBarrier_IsAvailableAndReturns_OnThisPlatform()
    {
        // the protocol needs it on every Dispose of a tree that other threads read; it must not throw, whatever the platform does underneath
        Assert.DoesNotThrow(Interlocked.MemoryBarrierProcessWide);
    }

    private sealed class Blocking(CountdownEvent inside, ManualResetEventSlim go) : ISuffixTreeVisitor
    {
        private int _nodes;

        public void VisitNode(int startIndex, int endIndex, int leafCount, int childCount, int depth)
        {
            if (_nodes++ != 2) return;
            inside.Signal();
            go.Wait(Wait);
        }

        public void EnterBranch(int key)
        {
        }

        public void ExitBranch()
        {
        }
    }
}
