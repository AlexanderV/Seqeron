using SuffixTree.EmProto.Ph5;

namespace SuffixTree.EmProto.Tests.Ph5;

/// <summary>
/// The writer of the CHILD, PSV and NSV files: a file of known length through a few aligned mapped windows (the owner's rule of 2026-11-04: stream into memory-mapped files), released least recently used first,
/// flushed when released so that the dirty pages stay bounded.
/// </summary>
[TestFixture]
public class WindowWriterTests
{
    private string _dir = "";

    [SetUp]
    public void SetUp()
    {
        _dir = Path.Combine(Path.GetTempPath(), "win-" + Guid.NewGuid().ToString("N"));
        Directory.CreateDirectory(_dir);
    }

    [TearDown]
    public void TearDown() => Directory.Delete(_dir, true);

    [Test]
    public void RandomWrites_MatchAModel_ForEveryWindowAndSlotCount()
    {
        var rng = new Random(77);
        foreach (int window in new[] { 8, 16, 64, 4096 })
        {
            foreach (int slots in new[] { 1, 2, 3, 8 })
            {
                long entries = 1 + rng.Next(3000);
                var model = new uint[entries];
                string path = Path.Combine(_dir, $"r-{window}-{slots}.bin");
                using (var w = new AlignedWindowWriter(path, entries * 4, window, slots, flushOnRelease: slots % 2 == 0))
                {
                    Assert.That(new FileInfo(path).Length, Is.EqualTo(entries * 4), "the file has its final length at once");
                    for (int k = 0; k < 5000; k++)
                    {
                        long i = rng.Next((int)entries);
                        uint v = (uint)rng.Next();
                        w.Write32(i, v);
                        model[i] = v;
                    }
                }
                Assert.That(ChildFiles.ReadUInts(path), Is.EqualTo(model), $"window {window}, slots {slots}");
            }
        }
    }

    [Test]
    public void SequentialAndBackwardWrites_OpenEveryWindowOnce()
    {
        const int entries = 1000;
        const int window = 64; // 16 entries
        long windows = ((entries * 4) + window - 1) / window;
        foreach (bool backward in new[] { false, true })
        {
            string path = Path.Combine(_dir, backward ? "b.bin" : "f.bin");
            WindowStats stats;
            using (var w = new AlignedWindowWriter(path, entries * 4L, window, 2, true))
            {
                for (int k = 0; k < entries; k++) w.Write32(backward ? entries - 1 - k : k, (uint)k);
                stats = w.Stats;
            }
            Assert.That(stats.Opened, Is.EqualTo(windows), "every window mapped once");
            Assert.That(stats.Reopened, Is.Zero);
            Assert.That(stats.DistinctWindows, Is.EqualTo(windows));
            Assert.That(stats.Writes, Is.EqualTo(entries));
            uint[] got = ChildFiles.ReadUInts(path);
            for (int k = 0; k < entries; k++) Assert.That(got[backward ? entries - 1 - k : k], Is.EqualTo((uint)k));
        }
    }

    [Test]
    public void LeastRecentlyUsedWindow_IsReleasedFirst_AndReopeningIsCounted()
    {
        const int window = 16; // 4 entries
        using var w = new AlignedWindowWriter(Path.Combine(_dir, "l.bin"), 4 * window * 4L, window, 2, true);
        w.Write32(0, 1);   // window 0
        w.Write32(4, 2);   // window 1
        w.Write32(0, 3);   // window 0 is the most recent now
        w.Write32(8, 4);   // window 2: window 1 goes
        Assert.That(w.Stats.Evictions, Is.EqualTo(1));
        w.Write32(0, 5);   // still open: no new mapping
        Assert.That(w.Stats.Opened, Is.EqualTo(3));
        w.Write32(4, 6);   // window 1 again: mapped a second time, window 2 goes (0 was used more recently)
        Assert.That(w.Stats.Opened, Is.EqualTo(4));
        Assert.That(w.Stats.Reopened, Is.EqualTo(1));
        Assert.That(w.Stats.DistinctWindows, Is.EqualTo(3));
    }

    [Test]
    public void Flush_OnRelease_CoversTheWrittenWindows_AndIsOffWhenNotAsked()
    {
        const int window = 4096;
        foreach (bool flush in new[] { true, false })
        {
            WindowStats stats;
            using (var w = new AlignedWindowWriter(Path.Combine(_dir, $"fl-{flush}.bin"), 10 * window, window, 2, flush))
            {
                for (int k = 0; k < 10; k++) w.Write32(k * (window / 4), 1);
                stats = w.Stats;
            }
            if (flush)
            {
                Assert.That(stats.Flushes, Is.EqualTo(8), "eight windows were released while writing");
                Assert.That(stats.FlushedBytes, Is.EqualTo(8L * window));
            }
            else
            {
                Assert.That(stats.Flushes, Is.Zero);
            }
        }
    }

    [Test]
    public void ByteWrites_CrossWindowBoundaries()
    {
        const int window = 16;
        string path = Path.Combine(_dir, "bw.bin");
        var model = new byte[100];
        using (var w = new AlignedWindowWriter(path, model.Length, window, 2, true))
        {
            var rng = new Random(5);
            for (int k = 0; k < 200; k++)
            {
                int offset = rng.Next(model.Length - 1);
                int length = 1 + rng.Next(Math.Min(40, model.Length - offset));
                var data = new byte[length];
                rng.NextBytes(data);
                w.Write(offset, data);
                data.CopyTo(model, offset);
            }
        }
        Assert.That(File.ReadAllBytes(path), Is.EqualTo(model));
    }

    [Test]
    public void ByteWrites_CrossWindowsOfAWholePage_WithoutLeavingTheMapping()
    {
        // windows of one page each: a piece that runs past the end of its window must be cut there, the rest goes to the next window (a window smaller than a page hides the error: the page is mapped around it)
        const int window = 4096;
        string path = Path.Combine(_dir, "bp.bin");
        var model = new byte[3 * window];
        using (var w = new AlignedWindowWriter(path, model.Length, window, 2, true))
        {
            var rng = new Random(9);
            foreach (int offset in new[] { window - 6, (2 * window) - 1, 100, window - 1, (2 * window) - 20 })
            {
                var data = new byte[Math.Min(30, model.Length - offset)];
                rng.NextBytes(data);
                w.Write(offset, data);
                data.CopyTo(model, offset);
            }
        }
        Assert.That(File.ReadAllBytes(path), Is.EqualTo(model));
    }

    [Test]
    public void WritesOutsideTheFile_AndBadParameters_AreRefused()
    {
        using var w = new AlignedWindowWriter(Path.Combine(_dir, "o.bin"), 40, 16, 2, true);
        w.Write32(9, 1);
        Assert.Throws<ArgumentOutOfRangeException>(() => w.Write32(10, 1));
        Assert.Throws<ArgumentOutOfRangeException>(() => w.Write32(-1, 1));
        Assert.Throws<ArgumentOutOfRangeException>(() => w.Write(38, new byte[3]));
        Assert.Throws<ArgumentOutOfRangeException>(() => new AlignedWindowWriter(Path.Combine(_dir, "p.bin"), 40, 0, 2, true));
        Assert.Throws<ArgumentOutOfRangeException>(() => new AlignedWindowWriter(Path.Combine(_dir, "q.bin"), 40, 16, 0, true));
        foreach (int window in new[] { 3, 6, 24, 4097, -16 }) Assert.Throws<ArgumentOutOfRangeException>(() => new AlignedWindowWriter(Path.Combine(_dir, "r.bin"), 40, window, 2, true), $"window {window}");
        Assert.Throws<ArgumentOutOfRangeException>(() => new AlignedWindowWriter(Path.Combine(_dir, "s.bin"), -1, 16, 2, true));
    }

    [Test]
    public void ExplicitFlush_FlushesEveryOpenWindow_AndDisposeFlushesTheRest()
    {
        const int window = 64;
        var w = new AlignedWindowWriter(Path.Combine(_dir, "ef.bin"), 4 * window, window, 3, flushOnRelease: true);
        w.Write32(0, 1);
        w.Write32(16, 2);
        Assert.That(w.Stats.Flushes, Is.Zero);
        w.Flush();
        Assert.That(w.Stats.Flushes, Is.EqualTo(2), "two windows are open");
        Assert.That(w.Stats.FlushedBytes, Is.EqualTo(2L * window));
        w.Write32(32, 3);
        w.Dispose();
        Assert.That(w.Stats.Flushes, Is.EqualTo(2 + 3), "dispose flushes the three open windows");
        w.Dispose();
        Assert.That(w.Stats.Flushes, Is.EqualTo(5), "a second dispose does nothing");
        Assert.Throws<ObjectDisposedException>(() => w.Write32(0, 1));
    }

    [Test]
    public void EmptyFile_IsCreated()
    {
        string path = Path.Combine(_dir, "z.bin");
        using (new AlignedWindowWriter(path, 0, 16, 2, true))
        {
        }
        Assert.That(new FileInfo(path).Length, Is.Zero);
    }
}
