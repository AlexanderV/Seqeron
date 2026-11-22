using System.Collections.Concurrent;
using System.Diagnostics;

namespace SuffixTree.Esa.Tests.Facade;

/// <summary>Records the points of the reader protocol that the gate reports (work package 1b.9g) and lets a test stop a named thread at one of them.</summary>
internal sealed class GateProbe
{
    private readonly List<string> _events = [];

    /// <summary>Called at every point, on the thread that reaches it, with <c>"&lt;thread name&gt;:&lt;point&gt;"</c>.</summary>
    internal Action<string>? Hook { get; set; }

    internal int Count(string point) => Events.Count(e => e.EndsWith(":" + point, StringComparison.Ordinal));

    internal string[] Events
    {
        get
        {
            lock (_events) return [.. _events];
        }
    }

    internal string[] EventsOf(string thread) => [.. Events.Where(e => e.StartsWith(thread + ":", StringComparison.Ordinal)).Select(e => e[(thread.Length + 1)..])];

    internal void Log(string point)
    {
        string entry = (Thread.CurrentThread.Name ?? "?") + ":" + point;
        lock (_events) _events.Add(entry);
        Hook?.Invoke(entry);
    }
}

internal readonly struct ProbeSeams(GateProbe probe) : IReaderGateSeams
{
    public void ProcessWideBarrier()
    {
        probe.Log("barrier");
        Interlocked.MemoryBarrierProcessWide();
    }

    public void AfterPublish() => probe.Log("published");

    public void AfterRetract() => probe.Log("retracted");

    public void BeforeScan() => probe.Log("before-scan");

    public void AfterScan(bool inside) => probe.Log(inside ? "scan-inside" : "scan-clear");
}

/// <summary>
/// The protocol of the readers of a file-backed tree (work package 1b.9g), step by step: a gate over the thread cells with the seams of <see cref="GateProbe"/>, so that a thread can be stopped between the store of its entry and the read of the
/// flag, or between the clearing of its entry and the read of the flag, and the releaser can run in that gap. These are the interleavings that the hardware of a weak memory model makes real and that no race test reaches reliably; the tests
/// run them in a sequentially consistent order, which checks the logic of the protocol (what is read, what is stored, in which order, who releases); the memory-model argument is in ADR-004 section 5a ("1b.9g").
/// </summary>
[TestFixture]
public class ReaderGateTests
{
    private static readonly TimeSpan Wait = TimeSpan.FromSeconds(30);

    private sealed class Setup
    {
        internal readonly GateProbe Probe = new();
        internal int Releases;
        internal string? Releaser;
        internal readonly ReaderGate<ProbeSeams> Gate;

        internal Setup()
        {
            Gate = new ReaderGate<ProbeSeams>(this, () =>
            {
                Releaser = Thread.CurrentThread.Name;
                Interlocked.Increment(ref Releases);
            }, new ProbeSeams(Probe));
        }
    }

    private static readonly ConcurrentDictionary<Thread, Exception> Failures = new();

    // a failure inside a thread of a test is the failure of the test, not of the host: the thread records it and Join reports it
    private static Thread Start(string name, Action body)
    {
        Thread? thread = null;
        thread = new Thread(() =>
        {
            try
            {
                body();
            }
            catch (Exception e)
            {
                Failures[thread!] = e;
            }
        }) { Name = name, IsBackground = true };
        thread.Start();
        return thread;
    }

    private static void Join(Thread thread)
    {
        Assert.That(thread.Join(Wait), Is.True, $"thread {thread.Name} did not finish");
        if (Failures.TryRemove(thread, out Exception? failure)) throw new AssertionException($"thread {thread.Name} failed: {failure}");
    }

    /// <summary>A reader that enters, signals that it is inside, stays until told to leave, and leaves.</summary>
    private static Thread Inside(Setup s, string name, ManualResetEventSlim inside, ManualResetEventSlim leave)
    {
        return Start(name, () =>
        {
            ReaderCell cell = s.Gate.Enter();
            inside.Set();
            leave.Wait(Wait);
            s.Gate.Leave(cell);
        });
    }

    [Test]
    public void ACallOfALiveGate_TakesNoBarrier_NoScan_AndNoRelease()
    {
        var s = new Setup();
        for (int i = 0; i < 500; i++) s.Gate.Leave(s.Gate.Enter());
        Assert.Multiple(() =>
        {
            Assert.That(s.Probe.Count("barrier"), Is.Zero);
            Assert.That(s.Probe.Count("before-scan"), Is.Zero);
            Assert.That(s.Probe.Count("published"), Is.EqualTo(500), "one entry per call");
            Assert.That(s.Probe.Count("retracted"), Is.EqualTo(500), "one clearing per call");
            Assert.That(s.Releases, Is.Zero);
            Assert.That(s.Gate.IsDisposed, Is.False);
            Assert.That(s.Gate.IsReleased, Is.False);
            Assert.That(ReaderCell.AnyInside(s.Gate.Id), Is.False);
        });
    }

    [Test]
    public void TheEntryOfACall_IsInTheCell_BetweenEnterAndLeave_AndOnlyThen()
    {
        var s = new Setup();
        Assert.That(ReaderCell.AnyInside(s.Gate.Id), Is.False);
        ReaderCell cell = s.Gate.Enter();
        Assert.That(ReaderCell.AnyInside(s.Gate.Id), Is.True);
        Assert.That(cell, Is.SameAs(ReaderCell.Current), "the cell of the calling thread");
        s.Gate.Leave(cell);
        Assert.That(ReaderCell.AnyInside(s.Gate.Id), Is.False);
    }

    [Test]
    public void Dispose_OfAnIdleGate_RaisesTheFlag_RunsTheBarrierThenTheScan_AndReleasesOnce()
    {
        var s = new Setup();
        Thread t = Start("disposer", s.Gate.Dispose);
        Join(t);
        Assert.Multiple(() =>
        {
            Assert.That(s.Probe.EventsOf("disposer"), Is.EqualTo(new[] { "barrier", "before-scan", "scan-clear" }), "the flag, then the barrier, then the scan, then the release");
            Assert.That(s.Releases, Is.EqualTo(1));
            Assert.That(s.Releaser, Is.EqualTo("disposer"));
            Assert.That(s.Gate.IsDisposed, Is.True);
            Assert.That(s.Gate.IsReleased, Is.True);
        });
    }

    [Test]
    public void Dispose_IsIdempotent_AndReleasesOnce_WhateverTheNumberOfDisposers()
    {
        var s = new Setup();
        using var go = new ManualResetEventSlim();
        var threads = Enumerable.Range(0, 8).Select(i => Start("d" + i, () =>
        {
            go.Wait(Wait);
            s.Gate.Dispose();
            s.Gate.Dispose();
        })).ToList();
        go.Set();
        threads.ForEach(Join);
        Assert.Multiple(() =>
        {
            Assert.That(s.Releases, Is.EqualTo(1));
            Assert.That(s.Probe.Count("barrier"), Is.EqualTo(1), "only the disposer that raised the flag runs the barrier: the others find it raised and return");
            Assert.That(s.Probe.Count("before-scan"), Is.EqualTo(1));
        });
    }

    [Test]
    public void Dispose_WithAReaderInside_DoesNotRelease_AndTheReaderThatLeavesReleases_AfterItsOwnBarrierAndScan()
    {
        var s = new Setup();
        using var inside = new ManualResetEventSlim();
        using var leave = new ManualResetEventSlim();
        Thread reader = Inside(s, "reader", inside, leave);
        Assert.That(inside.Wait(Wait), Is.True);
        Thread d = Start("disposer", s.Gate.Dispose);
        Join(d);
        Assert.Multiple(() =>
        {
            Assert.That(s.Probe.EventsOf("disposer"), Is.EqualTo(new[] { "barrier", "before-scan", "scan-inside" }));
            Assert.That(s.Releases, Is.Zero, "a call is inside: nothing is released");
            Assert.That(s.Gate.IsDisposed, Is.True);
            Assert.That(s.Gate.IsReleased, Is.False);
        });
        leave.Set();
        Join(reader);
        Assert.Multiple(() =>
        {
            Assert.That(s.Probe.EventsOf("reader"), Is.EqualTo(new[] { "published", "retracted", "barrier", "before-scan", "scan-clear" }), "the reader clears its entry, sees the flag, and runs the barrier and the scan itself");
            Assert.That(s.Releases, Is.EqualTo(1));
            Assert.That(s.Releaser, Is.EqualTo("reader"));
            Assert.That(s.Gate.IsReleased, Is.True);
        });
    }

    [Test]
    public void TwoReadersInside_TheLastOneToLeaveReleases_NotTheFirst()
    {
        var s = new Setup();
        using var insideA = new ManualResetEventSlim();
        using var insideB = new ManualResetEventSlim();
        using var leaveA = new ManualResetEventSlim();
        using var leaveB = new ManualResetEventSlim();
        Thread a = Inside(s, "a", insideA, leaveA);
        Thread b = Inside(s, "b", insideB, leaveB);
        Assert.That(insideA.Wait(Wait) && insideB.Wait(Wait), Is.True);
        s.Gate.Dispose();
        Assert.That(s.Releases, Is.Zero);
        leaveA.Set();
        Join(a);
        Assert.Multiple(() =>
        {
            Assert.That(s.Releases, Is.Zero, "the other reader is still inside");
            Assert.That(s.Probe.EventsOf("a").Last(), Is.EqualTo("scan-inside"));
        });
        leaveB.Set();
        Join(b);
        Assert.Multiple(() =>
        {
            Assert.That(s.Releases, Is.EqualTo(1));
            Assert.That(s.Releaser, Is.EqualTo("b"));
        });
    }

    [Test]
    public void ADisposeBetweenTheStoreAndTheCheckOfAnEntry_IsSeenByTheScan_AndTheReaderThatThenSeesTheFlagGoesAwayAndReleases()
    {
        var s = new Setup();
        using var published = new ManualResetEventSlim();
        using var resume = new ManualResetEventSlim();
        s.Probe.Hook = e =>
        {
            if (e != "reader:published") return;
            published.Set();
            resume.Wait(Wait);
        };
        Exception? thrown = null;
        Thread reader = Start("reader", () =>
        {
            try
            {
                s.Gate.Enter();
            }
            catch (Exception e)
            {
                thrown = e;
            }
        });
        Assert.That(published.Wait(Wait), Is.True);
        Thread d = Start("disposer", s.Gate.Dispose);
        Join(d);
        Assert.Multiple(() =>
        {
            Assert.That(s.Probe.EventsOf("disposer"), Is.EqualTo(new[] { "barrier", "before-scan", "scan-inside" }), "the scan finds the entry that was stored before the flag was read");
            Assert.That(s.Releases, Is.Zero);
        });
        resume.Set();
        Join(reader);
        Assert.Multiple(() =>
        {
            Assert.That(thrown, Is.TypeOf<ObjectDisposedException>(), "the second read sees the flag");
            Assert.That(s.Probe.EventsOf("reader"), Is.EqualTo(new[] { "published", "retracted", "barrier", "before-scan", "scan-clear" }), "it takes its entry back and releases, because the disposer left the release to it");
            Assert.That(s.Releases, Is.EqualTo(1));
            Assert.That(s.Releaser, Is.EqualTo("reader"));
            Assert.That(ReaderCell.AnyInside(s.Gate.Id), Is.False, "no entry is left behind");
        });
    }

    [Test]
    public void ADisposeBetweenTheClearingAndTheCheckOfAnExit_ReleasesInTheDisposer_AndTheReaderDoesNotReleaseAgain()
    {
        var s = new Setup();
        using var retracted = new ManualResetEventSlim();
        using var resume = new ManualResetEventSlim();
        s.Probe.Hook = e =>
        {
            if (e != "reader:retracted") return;
            retracted.Set();
            resume.Wait(Wait);
        };
        using var inside = new ManualResetEventSlim();
        using var leave = new ManualResetEventSlim();
        Thread reader = Inside(s, "reader", inside, leave);
        Assert.That(inside.Wait(Wait), Is.True);
        leave.Set();
        Assert.That(retracted.Wait(Wait), Is.True, "the reader has cleared its entry and is stopped before it reads the flag");
        Thread d = Start("disposer", s.Gate.Dispose);
        Join(d);
        Assert.Multiple(() =>
        {
            Assert.That(s.Probe.EventsOf("disposer"), Is.EqualTo(new[] { "barrier", "before-scan", "scan-clear" }));
            Assert.That(s.Releases, Is.EqualTo(1), "the entry was gone: the disposer releases");
        });
        resume.Set();
        Join(reader);
        Assert.Multiple(() =>
        {
            Assert.That(s.Releases, Is.EqualTo(1), "the reader sees the flag and finds the release claimed");
            Assert.That(s.Probe.Count("barrier"), Is.EqualTo(1), "and does not run another barrier");
            Assert.That(s.Releaser, Is.EqualTo("disposer"));
        });
    }

    [Test]
    public void ACallThatStartsAfterDispose_ThrowsObjectDisposed_WithoutPublishingAnEntry_AndWithoutABarrier()
    {
        var s = new Setup();
        using var inside = new ManualResetEventSlim();
        using var leave = new ManualResetEventSlim();
        Thread reader = Inside(s, "reader", inside, leave);
        Assert.That(inside.Wait(Wait), Is.True);
        s.Gate.Dispose();
        int barriers = s.Probe.Count("barrier");
        Exception? thrown = null;
        Join(Start("late", () =>
        {
            for (int i = 0; i < 50; i++)
            {
                try
                {
                    s.Gate.Enter();
                }
                catch (Exception e)
                {
                    thrown = e;
                }
            }
        }));
        Assert.Multiple(() =>
        {
            Assert.That(thrown, Is.TypeOf<ObjectDisposedException>());
            Assert.That(s.Probe.EventsOf("late"), Is.Empty, "the late caller never published an entry, took none back and ran no barrier");
            Assert.That(s.Probe.Count("barrier"), Is.EqualTo(barriers));
            Assert.That(s.Releases, Is.Zero, "the reader that was inside still holds the release back");
        });
        leave.Set();
        Join(reader);
        Assert.That(s.Releases, Is.EqualTo(1));
    }

    [Test]
    public void AfterTheRelease_NoCallRunsTheBarrierAgain()
    {
        var s = new Setup();
        s.Gate.Dispose();
        Assert.That(s.Probe.Count("barrier"), Is.EqualTo(1));
        for (int i = 0; i < 20; i++) Assert.Throws<ObjectDisposedException>(() => s.Gate.Enter());
        s.Gate.Dispose();
        Assert.Multiple(() =>
        {
            Assert.That(s.Probe.Count("barrier"), Is.EqualTo(1));
            Assert.That(s.Releases, Is.EqualTo(1));
        });
    }

    [Test]
    public void ASecondDispose_WhileAReaderIsStillInside_RunsNoBarrierAndNoScan_OnlyTheDisposeThatRaisedTheFlagDoes()
    {
        var s = new Setup();
        using var inside = new ManualResetEventSlim();
        using var leave = new ManualResetEventSlim();
        Thread reader = Inside(s, "reader", inside, leave);
        Assert.That(inside.Wait(Wait), Is.True);
        s.Gate.Dispose();
        Assert.That(s.Probe.Count("barrier"), Is.EqualTo(1));
        s.Gate.Dispose();
        s.Gate.Dispose();
        Assert.Multiple(() =>
        {
            Assert.That(s.Probe.Count("barrier"), Is.EqualTo(1), "a Dispose that finds the flag raised returns at once");
            Assert.That(s.Probe.Count("before-scan"), Is.EqualTo(1));
            Assert.That(s.Releases, Is.Zero);
        });
        leave.Set();
        Join(reader);
        Assert.That(s.Releases, Is.EqualTo(1));
    }

    [Test]
    public void TheScan_FollowsTheBarrier_AndSeesAnEntryThatIsStoredWhileTheBarrierRuns()
    {
        // an entry that is still in a store buffer becomes visible at the barrier: the hook stores one in the cell of the disposing thread when the barrier runs; a scan that came before the barrier would miss it and release
        var s = new Setup();
        s.Probe.Hook = e =>
        {
            if (e.EndsWith(":barrier", StringComparison.Ordinal)) ReaderCell.Current.Publish(s.Gate.Id);
        };
        s.Gate.Dispose();
        try
        {
            Assert.Multiple(() =>
            {
                Assert.That(s.Releases, Is.Zero, "the scan saw the entry: it came after the barrier");
                Assert.That(s.Probe.Count("scan-inside"), Is.EqualTo(1));
            });
        }
        finally
        {
            ReaderCell.Current.Retract();
        }
    }

    [Test]
    public void ADisposeFromInsideACall_ReleasesWhenTheOutermostCallLeaves()
    {
        var s = new Setup();
        ReaderCell outer = s.Gate.Enter();
        ReaderCell inner = s.Gate.Enter();
        Assert.That(inner, Is.SameAs(outer));
        Assert.That(outer.Depth, Is.EqualTo(2));
        s.Gate.Dispose();
        Assert.That(s.Releases, Is.Zero, "this very thread is inside");
        s.Gate.Leave(inner);
        Assert.Multiple(() =>
        {
            Assert.That(s.Releases, Is.Zero, "the outer call is still inside");
            Assert.That(ReaderCell.AnyInside(s.Gate.Id), Is.True);
        });
        s.Gate.Leave(outer);
        Assert.Multiple(() =>
        {
            Assert.That(s.Releases, Is.EqualTo(1));
            Assert.That(ReaderCell.AnyInside(s.Gate.Id), Is.False);
            Assert.That(outer.Depth, Is.Zero);
        });
    }

    [Test]
    public void ManyReadersThatLeaveTogether_ReleaseOnce()
    {
        for (int round = 0; round < 30; round++)
        {
            var s = new Setup();
            const int readers = 6;
            using var go = new ManualResetEventSlim();
            using var inside = new CountdownEvent(readers);
            var threads = Enumerable.Range(0, readers).Select(i => Start("r" + i, () =>
            {
                ReaderCell cell = s.Gate.Enter();
                inside.Signal();
                go.Wait(Wait);
                s.Gate.Leave(cell);
            })).ToList();
            Assert.That(inside.Wait(Wait), Is.True);
            s.Gate.Dispose();
            go.Set();
            threads.ForEach(Join);
            Assert.That(s.Releases, Is.EqualTo(1), $"round {round}");
        }
    }

    [Test]
    public void TheIdsOfTwoGates_DifferAndAreNeverZero()
    {
        var a = new Setup();
        var b = new Setup();
        Assert.Multiple(() =>
        {
            Assert.That(a.Gate.Id, Is.Not.Zero);
            Assert.That(b.Gate.Id, Is.Not.Zero);
            Assert.That(a.Gate.Id, Is.Not.EqualTo(b.Gate.Id));
        });
    }

    [Test]
    public void TheEntriesOfTwoGates_AreSeparate_WhenAThreadIsInsideBoth()
    {
        var a = new Setup();
        var b = new Setup();
        ReaderCell outer = a.Gate.Enter();
        ReaderCell inner = b.Gate.Enter();
        Assert.Multiple(() =>
        {
            Assert.That(ReaderCell.AnyInside(a.Gate.Id), Is.True);
            Assert.That(ReaderCell.AnyInside(b.Gate.Id), Is.True);
        });
        b.Gate.Leave(inner);
        Assert.Multiple(() =>
        {
            Assert.That(ReaderCell.AnyInside(a.Gate.Id), Is.True);
            Assert.That(ReaderCell.AnyInside(b.Gate.Id), Is.False);
        });
        b.Gate.Dispose();
        Assert.That(b.Releases, Is.EqualTo(1), "a thread that is inside another gate does not hold this one");
        a.Gate.Leave(outer);
        a.Gate.Dispose();
        Assert.That(a.Releases, Is.EqualTo(1));
    }

    [Test]
    public void TheArgumentsOfTheConstructor_AreChecked()
    {
        Assert.Throws<ArgumentNullException>(() => _ = new ReaderGate<ReaderGateSeams>(null!, () => { }));
        Assert.Throws<ArgumentNullException>(() => _ = new ReaderGate<ReaderGateSeams>(new object(), null!));
    }

    [Test]
    public void TheExceptionOfAFailedEnter_NamesTheOwner()
    {
        var owner = new Setup();
        var gate = new ReaderGate<ReaderGateSeams>(owner, () => { });
        gate.Dispose();
        var e = Assert.Throws<ObjectDisposedException>(() => gate.Enter());
        Assert.That(e!.ObjectName, Is.EqualTo(typeof(Setup).FullName));
    }

    /// <summary>
    /// The invariant of the whole protocol at the level of the gate: readers loop on a resource that the release invalidates, a disposer raises the flag after a random delay, and no reader may find the resource released while it is inside a call
    /// (at the start of the call or at its end). Under a sequentially consistent order this fails when the second read of the flag, the scan, the release claim or the barrier-then-scan order is wrong; the weak-order cases are the
    /// subject of the proof in ADR-004.
    /// </summary>
    [Test]
    public void ReadersNeverFindTheResourceReleased_WhileTheyAreInsideACall_WhateverTheMomentOfTheDispose()
    {
        const int rounds = 400;
        const int readers = 4;
        long calls = 0, refused = 0;
        var violations = new ConcurrentBag<string>();
        var rng = new Random(20261122);
        for (int round = 0; round < rounds; round++)
        {
            int released = 0;
            var gate = new ReaderGate<ReaderGateSeams>(this, () => Volatile.Write(ref released, 1));
            using var start = new Barrier(readers + 1);
            var threads = Enumerable.Range(0, readers).Select(_ => new Thread(() =>
            {
                start.SignalAndWait();
                while (true)
                {
                    ReaderCell cell;
                    try
                    {
                        cell = gate.Enter();
                    }
                    catch (ObjectDisposedException)
                    {
                        Interlocked.Increment(ref refused);
                        return;
                    }
                    catch (Exception e)
                    {
                        violations.Add($"round {round}: {e.GetType().Name}");
                        Interlocked.Increment(ref refused);
                        return;
                    }

                    if (Volatile.Read(ref released) != 0) violations.Add($"round {round}: released at the start of a call");
                    Thread.SpinWait(20);
                    if (Volatile.Read(ref released) != 0) violations.Add($"round {round}: released inside a call");
                    gate.Leave(cell);
                    Interlocked.Increment(ref calls);
                }
            })).ToList();
            threads.ForEach(t => t.Start());
            start.SignalAndWait();
            long end = Stopwatch.GetTimestamp() + (rng.Next(0, 300) * Stopwatch.Frequency / 1_000_000);
            while (Stopwatch.GetTimestamp() < end) Thread.SpinWait(10);
            gate.Dispose();
            threads.ForEach(t => Assert.That(t.Join(Wait), Is.True));
            Assert.That(Volatile.Read(ref released), Is.EqualTo(1), $"round {round}: every reader left, so the release must have run");
        }

        TestContext.Out.WriteLine($"calls={calls} refused={refused}");
        Assert.That(violations, Is.Empty, string.Join(" | ", violations.Take(3)));
        Assert.That(refused, Is.EqualTo((long)rounds * readers), "every reader ends at the disposal");
    }
}
