using System.Reflection;
using System.Runtime.InteropServices;

namespace SuffixTree.Esa.Tests.Facade;

/// <summary>
/// The registry of the thread-exclusive reader slots (work package 1b.9g): one cell per thread, found without a lock or a fence, whatever the number of trees; the cell of a thread that has ended is taken over, so that many short-lived threads
/// do not grow the registry; an entry is seen by a scan from another thread while it is stored and not after it is cleared, for the outermost call and for calls nested to any depth; two cells never share a cache line.
/// </summary>
[TestFixture]
public class ReaderCellTests
{
    private static readonly TimeSpan Wait = TimeSpan.FromSeconds(30);

    private static void OnThread(Action body)
    {
        Exception? failure = null;
        var thread = new Thread(() =>
        {
            try
            {
                body();
            }
            catch (Exception e)
            {
                failure = e;
            }
        }) { IsBackground = true };
        thread.Start();
        Assert.That(thread.Join(Wait), Is.True);
        if (failure is not null) throw new AssertionException("on the thread: " + failure);
    }

    private static bool ScanFromAnotherThread(long id)
    {
        bool seen = false;
        OnThread(() => seen = ReaderCell.AnyInside(id));
        return seen;
    }

    [Test]
    public void Current_IsOneCellPerThread_StableOnThatThread_AndAnotherCellOnAnotherThread()
    {
        ReaderCell mine = ReaderCell.Current;
        Assert.That(ReaderCell.Current, Is.SameAs(mine));
        ReaderCell? other = null;
        ReaderCell? again = null;
        OnThread(() =>
        {
            other = ReaderCell.Current;
            again = ReaderCell.Current;
        });
        Assert.Multiple(() =>
        {
            Assert.That(other, Is.Not.Null);
            Assert.That(other, Is.Not.SameAs(mine));
            Assert.That(again, Is.SameAs(other));
        });
    }

    [Test]
    public void AnEntry_IsSeenByAScanFromAnotherThread_WhileItIsStored_AndNotAfterItIsCleared()
    {
        long id = ReaderCell.NextId();
        ReaderCell cell = ReaderCell.Current;
        Assert.That(ScanFromAnotherThread(id), Is.False);
        cell.Publish(id);
        Assert.That(ScanFromAnotherThread(id), Is.True);
        Assert.That(ScanFromAnotherThread(ReaderCell.NextId()), Is.False, "another tree's id is not held");
        cell.Retract();
        Assert.That(ScanFromAnotherThread(id), Is.False);
    }

    [Test]
    public void NestedEntries_OfManyLevelsAndTrees_AreAllSeen_AndLeaveOneByOneInTheReverseOrder()
    {
        const int levels = 37; // 1 outermost + 36 nested: four full segments of 8 and a part
        long[] ids = Enumerable.Range(0, levels).Select(_ => ReaderCell.NextId()).ToArray();
        ReaderCell cell = ReaderCell.Current;
        Assert.That(cell.Depth, Is.Zero);
        for (int i = 0; i < levels; i++)
        {
            cell.Publish(ids[i]);
            Assert.That(cell.Depth, Is.EqualTo(i + 1));
            for (int k = 0; k <= i; k++) Assert.That(ReaderCell.AnyInside(ids[k]), Is.True, $"level {k} while {i} are inside");
            if (i + 1 < levels) Assert.That(ReaderCell.AnyInside(ids[i + 1]), Is.False);
        }

        Assert.That(ScanFromAnotherThread(ids[levels - 1]), Is.True, "the deepest entry is seen from another thread");
        for (int i = levels - 1; i >= 0; i--)
        {
            cell.Retract();
            Assert.That(cell.Depth, Is.EqualTo(i));
            Assert.That(ReaderCell.AnyInside(ids[i]), Is.False, $"level {i} left");
            for (int k = 0; k < i; k++) Assert.That(ReaderCell.AnyInside(ids[k]), Is.True, $"level {k} stays while level {i} leaves");
        }

        Assert.That(ids.Any(ReaderCell.AnyInside), Is.False);
    }

    [Test]
    public void TheSameTreeEntered_AtSeveralLevels_IsInsideUntilTheOutermostCallLeaves()
    {
        long id = ReaderCell.NextId();
        ReaderCell cell = ReaderCell.Current;
        for (int i = 0; i < 12; i++) cell.Publish(id);
        for (int i = 12; i > 1; i--)
        {
            cell.Retract();
            Assert.That(ReaderCell.AnyInside(id), Is.True, $"{i - 1} calls are still inside");
        }

        cell.Retract();
        Assert.That(ReaderCell.AnyInside(id), Is.False);
    }

    [Test]
    public void ALevelThatIsReusedAfterItLeft_IsStoredAgain()
    {
        long a = ReaderCell.NextId(), b = ReaderCell.NextId(), c = ReaderCell.NextId();
        ReaderCell cell = ReaderCell.Current;
        cell.Publish(a);
        cell.Publish(b);
        cell.Retract();
        cell.Publish(c);
        Assert.Multiple(() =>
        {
            Assert.That(ReaderCell.AnyInside(a), Is.True);
            Assert.That(ReaderCell.AnyInside(b), Is.False);
            Assert.That(ReaderCell.AnyInside(c), Is.True);
        });
        cell.Retract();
        cell.Retract();
        Assert.That(new[] { a, b, c }.Any(ReaderCell.AnyInside), Is.False);
    }

    [Test]
    public void ThreadsThatEnd_HandTheirCellOn_SoTheRegistryStaysAsSmallAsTheLargestNumberOfLiveThreads()
    {
        _ = ReaderCell.Current; // this thread has a cell
        int before = ReaderCell.CellCount;
        long id = ReaderCell.NextId();
        for (int i = 0; i < 300; i++)
        {
            OnThread(() =>
            {
                ReaderCell cell = ReaderCell.Current;
                cell.Publish(id);
                cell.Retract();
            });
        }

        int after = ReaderCell.CellCount;
        TestContext.Out.WriteLine($"cells before {before}, after 300 sequential threads {after}");
        Assert.That(after - before, Is.LessThanOrEqualTo(1), "300 threads one after the other need one cell, not 300");
    }

    [Test]
    public void ThreadsAliveAtTheSameTime_NeedACellEach_AndTheirCellsAreTakenOverWhenTheyEnd()
    {
        const int n = 12;
        _ = ReaderCell.Current;
        using var go = new ManualResetEventSlim();
        using var holding = new CountdownEvent(n);
        var cells = new ReaderCell[n];
        long[] ids = Enumerable.Range(0, n).Select(_ => ReaderCell.NextId()).ToArray();
        var failures = new System.Collections.Concurrent.ConcurrentBag<Exception>();
        var threads = Enumerable.Range(0, n).Select(i =>
        {
            var t = new Thread(() =>
            {
                try
                {
                    cells[i] = ReaderCell.Current;
                    cells[i].Publish(ids[i]);
                }
                catch (Exception e)
                {
                    failures.Add(e);
                }

                holding.Signal();
                go.Wait(Wait);
                try
                {
                    cells[i].Retract();
                }
                catch (Exception e)
                {
                    failures.Add(e);
                }
            }) { IsBackground = true };
            t.Start();
            return t;
        }).ToList();
        Assert.That(holding.Wait(Wait), Is.True);
        Assert.Multiple(() =>
        {
            Assert.That(ReaderCell.CellCount, Is.GreaterThanOrEqualTo(n + 1), "twelve threads and this one are alive: the registry has at least thirteen cells");
            Assert.That(cells.Distinct().Count(), Is.EqualTo(n), "no two live threads share a cell");
            foreach (long id in ids) Assert.That(ReaderCell.AnyInside(id), Is.True);
        });
        go.Set();
        threads.ForEach(t => Assert.That(t.Join(Wait), Is.True));
        Assert.That(failures, Is.Empty);
        Assert.That(ids.Any(ReaderCell.AnyInside), Is.False, "every thread has left");
        int peak = ReaderCell.CellCount;

        // twelve more threads, one at a time: they take over the cells of the twelve that ended
        for (int i = 0; i < n * 3; i++) OnThread(() => ReaderCell.Current.Publish(ReaderCell.NextId()));
        Assert.That(ReaderCell.CellCount, Is.EqualTo(peak), "the cells of the ended threads were reused");
    }

    [Test]
    public void ACellThatIsTakenOver_StartsEmpty_EvenWhenItsThreadEndedInsideACall()
    {
        long stale = ReaderCell.NextId();
        ReaderCell? first = null;
        OnThread(() =>
        {
            first = ReaderCell.Current;
            first.Publish(stale);
            for (int i = 0; i < 10; i++) first.Publish(ReaderCell.NextId()); // the thread ends with eleven calls "inside"
        });
        ReaderCell? second = null;
        long fresh = ReaderCell.NextId();
        OnThread(() =>
        {
            second = ReaderCell.Current;
            Assert.That(second.Depth, Is.Zero, "a taken-over cell is empty");
            second.Publish(fresh);
            Assert.That(ReaderCell.AnyInside(stale), Is.False, "the entries of the ended thread are gone");
            second.Retract();
        });
        Assert.That(second, Is.SameAs(first), "the cell of the ended thread was the one that was free");
    }

    [Test]
    public void ManyIds_OneThread_AddNoCell()
    {
        ReaderCell cell = ReaderCell.Current;
        int before = ReaderCell.CellCount;
        for (int i = 0; i < 5_000; i++)
        {
            long id = ReaderCell.NextId();
            cell.Publish(id);
            Assert.That(ReaderCell.AnyInside(id), Is.True);
            cell.Retract();
        }

        Assert.That(ReaderCell.CellCount, Is.EqualTo(before));
    }

    [Test]
    public void TheIdsOfTheGates_AreStrictlyIncreasing_AndNeverZero()
    {
        long previous = 0;
        for (int i = 0; i < 100; i++)
        {
            long id = ReaderCell.NextId();
            Assert.That(id, Is.GreaterThan(previous));
            previous = id;
        }

        Assert.That(previous, Is.GreaterThan(0));
    }

    [Test]
    public void TheWrittenFields_HaveRoomForAFullCacheLineOnBothSides_SoNoTwoCellsShareALine()
    {
        var layout = typeof(ReaderCell).StructLayoutAttribute!;
        Assert.That(layout.Value, Is.EqualTo(LayoutKind.Explicit));
        int size = layout.Size;
        FieldInfo[] instance = typeof(ReaderCell).GetFields(BindingFlags.Instance | BindingFlags.NonPublic | BindingFlags.Public);
        var offsets = instance.ToDictionary(f => f.Name, f => f.GetCustomAttribute<FieldOffsetAttribute>()!.Value);
        // the fields that the owner writes on every call
        int[] hot = [offsets["_first"], offsets["_nestedTop"], offsets["_segments"]];
        Assert.Multiple(() =>
        {
            Assert.That(size, Is.GreaterThanOrEqualTo(256), "at least two lines of room around the hot fields");
            Assert.That(hot.Min(), Is.GreaterThanOrEqualTo(64 + 64), "the hot fields start 128 bytes into the object");
            Assert.That(size - (hot.Max() + 8), Is.GreaterThanOrEqualTo(64 + 64), "and end 128 bytes before its end");
            foreach (KeyValuePair<string, int> cold in offsets.Where(o => o.Key is "_next" or "_owner")) Assert.That(cold.Value, Is.LessThan(64), $"{cold.Key} is read by a scan, away from the hot fields");
        });
    }
}
