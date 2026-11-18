using SuffixTree.Esa.External;

namespace SuffixTree.Esa.Tests.External;

/// <summary>The byte and call counters of the explicit I/O, of the mapped windows and of the scratch files.</summary>
[TestFixture]
public class IoCountersTests
{
    [Test]
    public void EveryCounter_StartsAtZero_AndCountsWhatItIsGiven()
    {
        var c = new IoCounters();
        Assert.That(c.Snapshot(), Is.EqualTo(default(IoSnapshot)));
        c.AddRead(100);
        c.AddRead(28);
        c.AddWrite(7);
        c.AddMappedWrite(1000);
        c.AddMappedWrite(24);
        c.AddWindow(4096);
        c.AddWindow(512);
        Assert.That(c.ReadBytes, Is.EqualTo(128));
        Assert.That(c.ReadCalls, Is.EqualTo(2));
        Assert.That(c.WriteBytes, Is.EqualTo(7));
        Assert.That(c.WriteCalls, Is.EqualTo(1));
        Assert.That(c.MappedWriteBytes, Is.EqualTo(1024));
        Assert.That(c.MappedWindows, Is.EqualTo(2));
        Assert.That(c.MappedMaxWindowBytes, Is.EqualTo(4096), "the largest window, not the last");
        IoSnapshot s = c.Snapshot();
        Assert.That((s.ReadBytes, s.ReadCalls, s.WriteBytes, s.WriteCalls, s.MappedWriteBytes, s.MappedWindows, s.MappedMaxWindowBytes), Is.EqualTo((128L, 2L, 7L, 1L, 1024L, 2L, 4096L)));
        Assert.That(s.TotalBytes, Is.EqualTo(128 + 7 + 1024));
    }

    [Test]
    public void Scratch_TracksLiveAndPeak()
    {
        var c = new IoCounters();
        c.AddScratch(10);
        c.AddScratch(5);
        Assert.That((c.ScratchLiveBytes, c.ScratchPeakBytes), Is.EqualTo((15L, 15L)));
        c.RemoveScratch(12);
        Assert.That((c.ScratchLiveBytes, c.ScratchPeakBytes), Is.EqualTo((3L, 15L)));
        c.AddScratch(20);
        Assert.That((c.ScratchLiveBytes, c.ScratchPeakBytes), Is.EqualTo((23L, 23L)));
        c.RemoveScratch(23);
        Assert.That((c.ScratchLiveBytes, c.ScratchPeakBytes), Is.EqualTo((0L, 23L)));
    }

    [Test]
    public void Reset_ClearsEverything()
    {
        var c = new IoCounters();
        c.AddRead(1);
        c.AddWrite(2);
        c.AddMappedWrite(3);
        c.AddWindow(4);
        c.AddScratch(5);
        c.Reset();
        Assert.That(c.Snapshot(), Is.EqualTo(default(IoSnapshot)));
        c.AddWindow(2);
        Assert.That(c.MappedMaxWindowBytes, Is.EqualTo(2), "the maximum starts again after a reset");
    }

    [Test]
    public void Since_SubtractsTheCumulativeFields_AndKeepsTheMaximaAndTheLiveScratch()
    {
        var c = new IoCounters();
        c.AddRead(10);
        c.AddWrite(20);
        c.AddMappedWrite(30);
        c.AddWindow(64);
        c.AddScratch(100);
        IoSnapshot first = c.Snapshot();
        c.AddRead(1);
        c.AddRead(2);
        c.AddWrite(4);
        c.AddMappedWrite(8);
        c.AddWindow(32);
        c.AddScratch(50);
        c.RemoveScratch(120);
        IoSnapshot delta = c.Snapshot().Since(first);
        Assert.That((delta.ReadBytes, delta.ReadCalls, delta.WriteBytes, delta.WriteCalls, delta.MappedWriteBytes, delta.MappedWindows), Is.EqualTo((3L, 2L, 4L, 1L, 8L, 1L)));
        Assert.That(delta.MappedMaxWindowBytes, Is.EqualTo(64), "a maximum is not a difference");
        Assert.That(delta.ScratchLiveBytes, Is.EqualTo(30));
        Assert.That(delta.ScratchPeakBytes, Is.EqualTo(150));
    }

    [Test]
    public void Threads_AddUpExactly()
    {
        var c = new IoCounters();
        Parallel.For(0, 8, new ParallelOptions { MaxDegreeOfParallelism = 8 }, _ =>
        {
            for (int i = 0; i < 50_000; i++)
            {
                c.AddRead(3);
                c.AddWrite(5);
                c.AddMappedWrite(7);
                c.AddScratch(11);
                c.RemoveScratch(11);
                c.AddWindow(i);
            }
        });
        Assert.That(c.ReadBytes, Is.EqualTo(8L * 50_000 * 3));
        Assert.That(c.ReadCalls, Is.EqualTo(8L * 50_000));
        Assert.That(c.WriteBytes, Is.EqualTo(8L * 50_000 * 5));
        Assert.That(c.WriteCalls, Is.EqualTo(8L * 50_000));
        Assert.That(c.MappedWriteBytes, Is.EqualTo(8L * 50_000 * 7));
        Assert.That(c.MappedWindows, Is.EqualTo(8L * 50_000));
        Assert.That(c.MappedMaxWindowBytes, Is.EqualTo(49_999));
        Assert.That(c.ScratchLiveBytes, Is.Zero);
        Assert.That(c.ScratchPeakBytes, Is.InRange(11, 8 * 11));
    }
}

/// <summary>The running maximum used by the peaks.</summary>
[TestFixture]
public class AtomicTests
{
    [Test]
    public void RaiseTo_OnlyRaises()
    {
        long x = 5;
        Atomic.RaiseTo(ref x, 3);
        Assert.That(x, Is.EqualTo(5));
        Atomic.RaiseTo(ref x, 5);
        Assert.That(x, Is.EqualTo(5));
        Atomic.RaiseTo(ref x, 9);
        Assert.That(x, Is.EqualTo(9));
        Atomic.RaiseTo(ref x, long.MinValue);
        Assert.That(x, Is.EqualTo(9));
    }

    [Test]
    public void RaiseTo_UnderContention_EndsAtTheMaximum()
    {
        long x = 0;
        Parallel.For(0, 1000, i => Atomic.RaiseTo(ref x, i * 7L));
        Assert.That(x, Is.EqualTo(999 * 7L));
    }
}

/// <summary>The kernel text parser.</summary>
[TestFixture]
public class ProcTextTests
{
    private static bool Get(string text, string key, out long value) => ProcText.TryGetLong(System.Text.Encoding.ASCII.GetBytes(text), System.Text.Encoding.ASCII.GetBytes(key), out value);

    [Test]
    public void ParsesAKeyedLine_WithTabsAndSpaces_AndATrailingUnit()
    {
        Assert.That(Get("Name:\tdotnet\nRssAnon:\t   12345 kB\nRssFile:\t 77 kB\n", "RssAnon:", out long v), Is.True);
        Assert.That(v, Is.EqualTo(12345));
        Assert.That(Get("rchar: 99\nwchar: 7\n", "wchar:", out v), Is.True);
        Assert.That(v, Is.EqualTo(7));
        Assert.That(Get("rchar: 0", "rchar:", out v), Is.True, "last line without a newline");
        Assert.That(v, Is.Zero);
        Assert.That(Get("A:5\n", "A:", out v), Is.True, "no blank at all");
        Assert.That(v, Is.EqualTo(5));
    }

    [Test]
    public void MatchesWholeKeysAtTheStartOfALine_Only()
    {
        // the trap of /proc/self/io: cancelled_write_bytes ends in the key write_bytes
        const string io = "rchar: 1\nwchar: 2\nsyscr: 3\nsyscw: 4\nread_bytes: 5\nwrite_bytes: 6\ncancelled_write_bytes: 700\n";
        Assert.That(Get(io, "write_bytes:", out long v), Is.True);
        Assert.That(v, Is.EqualTo(6));
        Assert.That(Get("cancelled_write_bytes: 700\n", "write_bytes:", out _), Is.False, "a key inside another key is not a match");
        Assert.That(Get("xRssAnon: 5\n", "RssAnon:", out _), Is.False);
        Assert.That(Get("Name: RssAnon: 9\n", "RssAnon:", out _), Is.False, "a key inside a value is not a match");
        Assert.That(Get("RssFile: 1\nRssAnon: 8\n", "RssAnon:", out v), Is.True, "second line");
        Assert.That(v, Is.EqualTo(8));
        Assert.That(Get("RssAnon: 8\nRssAnon: 9\n", "RssAnon:", out v), Is.True);
        Assert.That(v, Is.EqualTo(8), "the first match wins");
    }

    [Test]
    public void AnEmptyKey_NeverMatches()
    {
        Assert.That(ProcText.TryGetLong("5 apples\n"u8, ""u8, out long v), Is.False);
        Assert.That(v, Is.Zero);
    }

    [Test]
    public void ReportsFalse_ForAMissingKey_ANonNumber_AnEmptyText_AndAnOverflow()
    {
        Assert.That(Get("", "A:", out long v), Is.False);
        Assert.That(v, Is.Zero);
        Assert.That(Get("B: 1\n", "A:", out _), Is.False);
        Assert.That(Get("A: kB\n", "A:", out _), Is.False, "no digits");
        Assert.That(Get("A:\n7\n", "A:", out _), Is.False, "the number must be on the line");
        Assert.That(Get("A:", "A:", out _), Is.False, "key at the very end");
        Assert.That(Get("A: 9223372036854775807\n", "A:", out v), Is.True);
        Assert.That(v, Is.EqualTo(long.MaxValue));
        Assert.That(Get("A: 9223372036854775808\n", "A:", out _), Is.False, "above long.MaxValue");
        Assert.That(Get("A: 99999999999999999999999\n", "A:", out _), Is.False);
        Assert.That(Get("A: -5\n", "A:", out _), Is.False, "the kernel prints unsigned numbers");
    }
}

/// <summary>/proc/self/io through an injectable source (G4: the counters degrade where /proc does not exist).</summary>
[TestFixture]
public class ProcIoTests
{
    private const string Io = "rchar: 1000\nwchar: 2000\nsyscr: 30\nsyscw: 40\nread_bytes: 5000\nwrite_bytes: 6000\ncancelled_write_bytes: 77\n";

    [Test]
    public void ParsesAllSixFields()
    {
        ProcIoSnapshot? s = ProcIoSnapshot.TryRead(new FakeProcSource().With("/proc/self/io", Io));
        Assert.That(s, Is.EqualTo(new ProcIoSnapshot(1000, 2000, 30, 40, 5000, 6000)));
    }

    [Test]
    public void Since_IsTheFieldWiseDifference()
    {
        var a = new ProcIoSnapshot(10, 20, 30, 40, 50, 60);
        var b = new ProcIoSnapshot(11, 23, 35, 44, 55, 66);
        Assert.That(b.Since(a), Is.EqualTo(new ProcIoSnapshot(1, 3, 5, 4, 5, 6)));
    }

    [Test]
    public void NoProcFile_IsUnavailable_NotAnException()
    {
        Assert.That(ProcIoSnapshot.TryRead(new FakeProcSource()), Is.Null, "a macOS or Windows process: no such file");
    }

    [TestCase("rchar: 1\nwchar: 2\nsyscr: 3\nsyscw: 4\nread_bytes: 5\n", TestName = "a key is missing")]
    [TestCase("", TestName = "empty file")]
    [TestCase("rchar: x\nwchar: 2\nsyscr: 3\nsyscw: 4\nread_bytes: 5\nwrite_bytes: 6\n", TestName = "a value is not a number")]
    [TestCase("garbage garbage garbage", TestName = "garbage")]
    public void APartialOrGarbledFile_IsUnavailable(string text)
    {
        Assert.That(ProcIoSnapshot.TryRead(new FakeProcSource().With("/proc/self/io", text)), Is.Null);
    }

    private sealed class DishonestSource(string text) : IProcSource
    {
        public int TryRead(string path, Span<byte> buffer)
        {
            System.Text.Encoding.ASCII.GetBytes(text).AsSpan().CopyTo(buffer);
            return int.MaxValue; // claims far more than the buffer holds
        }
    }

    [Test]
    public void ASourceThatClaimsMoreBytesThanTheBufferHolds_IsReadWithinTheBuffer()
    {
        Assert.That(ProcIoSnapshot.TryRead(new DishonestSource(Io)), Is.EqualTo(new ProcIoSnapshot(1000, 2000, 30, 40, 5000, 6000)));
        Assert.That(new RssReader(new DishonestSource("RssAnon: 77 kB\n")).Read(), Is.EqualTo(new RssSnapshot(true, 77, 0, 0, 0)));
    }

    [Test]
    public void AHostileSource_IsUnavailable_NotAnException()
    {
        foreach (Exception e in new Exception[] { new IOException("x"), new UnauthorizedAccessException("x"), new NotSupportedException("x"), new ArgumentException("x") })
        {
            var source = new FakeProcSource().With("/proc/self/io", Io);
            source.Throws = e;
            Assert.That(ProcIoSnapshot.TryRead(source), Is.Null, e.GetType().Name);
        }
    }

    [Test]
    public void TheRealFileSystem_ReadsWhereProcExists_AndNeverThrowsWhereItDoesNot()
    {
        ProcIoSnapshot? s = ProcIoSnapshot.TryRead();
        if (File.Exists("/proc/self/io"))
        {
            Assert.That(s, Is.Not.Null);
            Assert.That(s!.Value.RChar, Is.GreaterThan(0));
            Assert.That(s.Value.SysCr, Is.GreaterThan(0));
        }
        var buffer = new byte[64];
        Assert.That(FileProcSource.Instance.TryRead("/no/such/dir/never", buffer), Is.EqualTo(-1));
        Assert.That(FileProcSource.Instance.TryRead(Path.GetTempPath(), buffer), Is.EqualTo(-1), "a directory is not a readable file");
        Assert.That(FileProcSource.Instance.TryRead("", buffer), Is.EqualTo(-1), "an empty path");
    }

    [Test]
    public void TheRealSource_ReadsTheStartOfALongFile_IntoASmallBuffer()
    {
        using var dir = new TempDir();
        string path = dir.File("p.txt");
        File.WriteAllText(path, new string('x', 100));
        var small = new byte[10];
        Assert.That(FileProcSource.Instance.TryRead(path, small), Is.EqualTo(10));
        var big = new byte[500];
        Assert.That(FileProcSource.Instance.TryRead(path, big), Is.EqualTo(100));
        File.WriteAllBytes(path, []);
        Assert.That(FileProcSource.Instance.TryRead(path, big), Is.Zero, "an empty file reads zero bytes, it does not fail");
    }

    [Test]
    public void IoMeasurement_ReportsTheCountersAndTheProcDelta_OfTheInterval()
    {
        var counters = new IoCounters();
        counters.AddRead(100);
        var source = new FakeProcSource().WithSeries("/proc/self/io", call => $"rchar: {1000 + (call * 500)}\nwchar: {call * 10}\nsyscr: 1\nsyscw: 2\nread_bytes: {call * 4096}\nwrite_bytes: 0\n");
        var m = new IoMeasurement(counters, source);
        counters.AddRead(28);
        counters.AddWrite(64);
        counters.AddMappedWrite(1024);
        IoReport r = m.Finish();
        Assert.That(r.Counted.ReadBytes, Is.EqualTo(28), "the counters of the interval, not since the start of the process");
        Assert.That(r.Counted.WriteBytes, Is.EqualTo(64));
        Assert.That(r.Counted.MappedWriteBytes, Is.EqualTo(1024));
        Assert.That(r.Proc, Is.EqualTo(new ProcIoSnapshot(500, 10, 0, 0, 4096, 0)));
    }

    [Test]
    public void IoMeasurement_WithoutProc_ReportsTheCountersAndNoProcDelta()
    {
        var counters = new IoCounters();
        var m = new IoMeasurement(counters, new FakeProcSource());
        counters.AddWrite(5);
        IoReport r = m.Finish();
        Assert.That(r.Counted.WriteBytes, Is.EqualTo(5));
        Assert.That(r.Proc, Is.Null);
        // a /proc that appears or vanishes between the two readings is unavailable as well
        var flaky = new FakeProcSource().WithSeries("/proc/self/io", call => call == 0 ? "rchar: 1\nwchar: 1\nsyscr: 1\nsyscw: 1\nread_bytes: 1\nwrite_bytes: 1\n" : "garbage");
        Assert.That(new IoMeasurement(counters, flaky).Finish().Proc, Is.Null);
        Assert.Throws<ArgumentNullException>(() => _ = new IoMeasurement(null!));
    }
}
