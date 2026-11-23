using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.Scaffolding;

/// <summary>
/// The allocation meter of the test project (work package 1b.9h, the flaky 1b.3 test): what it reads for code of known allocation, that it reports a collection in its window, and the three verdicts of
/// <see cref="AllocationReading.AssertLessThan"/>.
/// </summary>
[TestFixture]
public class AllocationMeterTests
{
    private static object? s_sink;          // a static field: the allocation cannot be elided or put on the stack

    private const int ArrayHeaderBytes = 24;
    private const long Slack = 64;

    [Test]
    public void CodeThatAllocatesNothing_ReadsZero_InAWindowWithoutCollections()
    {
        long x = 0;
        AllocationReading reading = AllocationMeter.Measure(() => x += 3);
        Assert.That(reading.Bytes, Is.Zero);
        Assert.That(reading.GcFree, Is.True);
        Assert.That(x, Is.EqualTo(3), "the code ran once");
    }

    [Test]
    public void ASmallArray_ReadsItsSize()
    {
        AllocationReading reading = AllocationMeter.Measure(() => s_sink = new long[11]);
        Assert.That(reading.Bytes, Is.InRange(8 * 11 + ArrayHeaderBytes - 8, 8 * 11 + ArrayHeaderBytes + Slack), "a long[11]: 88 bytes of elements and the header");
        Assert.That(reading.GcFree, Is.True);
    }

    [Test]
    public void ALargeArray_ReadsItsSize_WhateverTheStateOfTheThreadsAllocationContext()
    {
        for (int round = 0; round < 8; round++)
        {
            for (int k = 0; k < round * 37; k++) s_sink = new byte[64 + k];     // leave the rest of the allocation context at a different size each round
            AllocationReading reading = AllocationMeter.Measure(() => s_sink = new long[200_001]);
            Assert.That(reading.Bytes, Is.InRange(8 * 200_001L + ArrayHeaderBytes - 8, 8 * 200_001L + ArrayHeaderBytes + Slack), $"round {round}");
            Assert.That(reading.GcFree, Is.True, $"round {round}");
        }
    }

    [Test]
    public void ACollectionInTheWindow_IsReported()
    {
        AllocationReading reading = AllocationMeter.Measure(() => GC.Collect());
        Assert.That(reading.GcFree, Is.False, "an induced collection leaves the no-GC region, and the window is not free of collections");
    }

    [Test]
    public void TheRegionIsLeft_AfterTheWindow_AlsoWhenTheCodeThrows()
    {
        Assert.Throws<InvalidOperationException>(() => AllocationMeter.Measure(() => throw new InvalidOperationException("the code under test failed")));
        Assert.That(System.Runtime.GCSettings.LatencyMode, Is.Not.EqualTo(System.Runtime.GCLatencyMode.NoGCRegion), "a later window can open its own region");
        Assert.That(AllocationMeter.Measure(() => s_sink = new byte[10]).GcFree, Is.True);
    }

    [Test]
    public void ARegionThatIsAlreadyOpen_IsUsedAsItIs_AndLeftOpen()
    {
        Assert.That(GC.TryStartNoGCRegion(1 << 20), Is.True, "the outer region");
        try
        {
            AllocationReading reading = AllocationMeter.Measure(() => s_sink = new byte[10]);
            Assert.That(reading.GcFree, Is.True, "the window is in the region of the caller, which holds");
            Assert.That(System.Runtime.GCSettings.LatencyMode, Is.EqualTo(System.Runtime.GCLatencyMode.NoGCRegion), "the region of the caller is left as it was found (a start inside it would have ended it)");
            AllocationReading lost = AllocationMeter.Measure(() => GC.Collect());
            Assert.That(lost.GcFree, Is.False, "a collection that ends the region of the caller is reported");
        }
        finally
        {
            if (System.Runtime.GCSettings.LatencyMode == System.Runtime.GCLatencyMode.NoGCRegion) GC.EndNoGCRegion();
        }
    }

    [Test]
    public void NoCode_IsAnArgumentError()
    {
        Assert.Throws<ArgumentNullException>(() => AllocationMeter.Measure(null!));
    }

    [Test]
    public void TheCodeRunsOnceOnTheCallingThread()
    {
        int runs = 0;
        int thread = -1;
        AllocationMeter.Measure(() => { runs++; thread = Environment.CurrentManagedThreadId; });
        Assert.That(runs, Is.EqualTo(1));
        Assert.That(thread, Is.EqualTo(Environment.CurrentManagedThreadId));
    }

    [Test]
    public void TheMeterReadsTheCallingThreadOnly()
    {
        using var stop = new CancellationTokenSource();
        var helper = new Thread(() =>
        {
            var keep = new byte[8][];
            int i = 0;
            while (!stop.IsCancellationRequested)
            {
                keep[i++ & 7] = new byte[20_000];          // another thread allocates 20 KB at a time, all the time
                Thread.Yield();
            }
        })
        { IsBackground = true };
        helper.Start();
        try
        {
            for (int k = 0; k < 20; k++)
            {
                AllocationReading reading = AllocationMeter.Measure(() => Thread.SpinWait(200_000));
                Assert.That(reading.Bytes, Is.LessThan(1024), $"window {k}: the allocations of the other thread are not in the reading");
            }
        }
        finally
        {
            stop.Cancel();
            helper.Join();
        }
    }

    [Test]
    public void TheReading_IsExact_WhileAnotherThreadAllocatesLargeObjects_AndStartsBackgroundCollections()
    {
        using var stop = new CancellationTokenSource();
        int gen2Before = GC.CollectionCount(2);
        var helper = new Thread(() =>
        {
            var keep = new byte[8][];
            int i = 0;
            while (!stop.IsCancellationRequested)
            {
                keep[i++ & 7] = new byte[100_000];         // large objects (the threshold is 85,000 bytes): what starts the background collections of the concurrent collector
                stop.Token.WaitHandle.WaitOne(1);
            }
        })
        { IsBackground = true };
        helper.Start();
        try
        {
            for (int k = 0; k < 120; k++)
            {
                AllocationReading nothing = AllocationMeter.Measure(() => Thread.Yield());
                AllocationReading small = AllocationMeter.Measure(() => s_sink = new long[11]);
                AllocationReading large = AllocationMeter.Measure(() => s_sink = new long[200_001]);
                Assert.That(nothing.Bytes, Is.Zero, $"window {k}: nothing");
                Assert.That(small.Bytes, Is.InRange(8 * 11 + ArrayHeaderBytes - 8, 8 * 11 + ArrayHeaderBytes + Slack), $"window {k}: small");
                Assert.That(large.Bytes, Is.InRange(8 * 200_001L + ArrayHeaderBytes - 8, 8 * 200_001L + ArrayHeaderBytes + Slack), $"window {k}: large");
                Assert.That(nothing.GcFree && small.GcFree && large.GcFree, Is.True, $"window {k}: no collection ran in a window");
            }
        }
        finally
        {
            stop.Cancel();
            helper.Join();
        }

        TestContext.Out.WriteLine($"gen 2 collections during the run (the blocking ones that open each window included): {GC.CollectionCount(2) - gen2Before}");
    }

    [Test]
    public void AReadingBelowTheLimit_Passes_WhetherOrNotACollectionRan()
    {
        Assert.DoesNotThrow(() => new AllocationReading(100, true).AssertLessThan(1024, "a"));
        Assert.DoesNotThrow(() => new AllocationReading(100, false).AssertLessThan(1024, "b"));
    }

    [Test]
    public void AReadingAboveTheLimit_WithoutACollection_IsAFailureOfTheCode()
    {
        var e = Assert.Throws<AssertionException>(() => new AllocationReading(200_000, true).AssertLessThan(1024, "what the code allocates"));
        Assert.That(e!.Message, Does.Contain("what the code allocates").And.Contain("200,000").And.Contain("no garbage collection ran"));
    }

    [Test]
    public void AReadingAboveTheLimit_WithACollection_IsInconclusive_NotAFailure()
    {
        var e = Assert.Throws<InconclusiveException>(() => new AllocationReading(200_000, false).AssertLessThan(1024, "what the code allocates"));
        Assert.That(e!.Message, Does.Contain("what the code allocates").And.Contain("a garbage collection ran"));
    }

    [Test]
    public void TheLimitIsExclusive()
    {
        Assert.DoesNotThrow(() => new AllocationReading(1023, true).AssertLessThan(1024, "edge"));
        Assert.Throws<AssertionException>(() => new AllocationReading(1024, true).AssertLessThan(1024, "edge"));
    }
}
