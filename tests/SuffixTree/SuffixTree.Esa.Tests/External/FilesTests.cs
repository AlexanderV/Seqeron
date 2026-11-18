using SuffixTree.Esa.External;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.External;

/// <summary>The counting positional file.</summary>
[TestFixture]
public class CountedFileTests
{
    [Test]
    public void WriteThenRead_RoundTrips_AndCountsBytesAndCalls()
    {
        using var dir = new TempDir();
        var c = new IoCounters();
        byte[] data = Enumerable.Range(0, 1000).Select(i => (byte)(i * 7)).ToArray();
        using (var f = new CountedFile(dir.File("a.bin"), FileMode.Create, FileAccess.ReadWrite, c))
        {
            Assert.That(f.Path, Is.EqualTo(dir.File("a.bin")));
            f.WriteAt(data.AsSpan(0, 600), 0);
            f.WriteAt(data.AsSpan(600), 600);
            Assert.That(f.Length, Is.EqualTo(1000));
            var back = new byte[400];
            f.ReadExactAt(back, 300);
            Assert.That(back, Is.EqualTo(data.AsSpan(300, 400).ToArray()));
            f.Flush();
        }
        Assert.That(c.WriteBytes, Is.EqualTo(1000));
        Assert.That(c.WriteCalls, Is.EqualTo(2));
        Assert.That(c.ReadBytes, Is.EqualTo(400));
        Assert.That(c.ReadCalls, Is.GreaterThanOrEqualTo(1));
        Assert.That(File.ReadAllBytes(dir.File("a.bin")), Is.EqualTo(data));
    }

    [Test]
    public void ReadAt_ReturnsFewerBytesAtTheEnd_AndReadExactAtThrows()
    {
        using var dir = new TempDir();
        File.WriteAllBytes(dir.File("b.bin"), [1, 2, 3, 4, 5]);
        using var f = new CountedFile(dir.File("b.bin"), FileMode.Open, FileAccess.Read);
        var buf = new byte[10];
        Assert.That(f.ReadAt(buf, 2), Is.EqualTo(3));
        Assert.That(buf.AsSpan(0, 3).ToArray(), Is.EqualTo(new byte[] { 3, 4, 5 }));
        Assert.That(f.ReadAt(buf, 5), Is.Zero, "at the end");
        Assert.That(f.ReadAt(buf, 100), Is.Zero, "beyond the end");
        Assert.That(f.ReadAt(buf.AsSpan(0, 0), 0), Is.Zero, "an empty buffer");
        Assert.Throws<EndOfStreamException>(() => f.ReadExactAt(buf, 0), "10 bytes of a 5-byte file");
        Assert.Throws<EndOfStreamException>(() => f.ReadExactAt(buf.AsSpan(0, 4), 3), "a truncated tail");
        f.ReadExactAt(buf.AsSpan(0, 5), 0);
        f.ReadExactAt(buf.AsSpan(0, 0), 99); // nothing to read, nothing to fail
    }

    [Test]
    public void SetLength_ExtendsAndTruncates_AndAWriteBeyondTheEndExtends()
    {
        using var dir = new TempDir();
        using var f = new CountedFile(dir.File("c.bin"), FileMode.Create, FileAccess.ReadWrite);
        f.SetLength(100);
        Assert.That(f.Length, Is.EqualTo(100));
        f.WriteAt([9], 150);
        Assert.That(f.Length, Is.EqualTo(151));
        f.SetLength(10);
        Assert.That(f.Length, Is.EqualTo(10));
    }

    [Test]
    public void WorksWithoutCounters_AndFailsLoudly_OnMissingFilesAndBadOffsets()
    {
        using var dir = new TempDir();
        Assert.Throws<FileNotFoundException>(() => _ = new CountedFile(dir.File("missing.bin"), FileMode.Open, FileAccess.Read));
        Assert.Throws<DirectoryNotFoundException>(() => _ = new CountedFile(dir.File("no/such/x.bin"), FileMode.Create, FileAccess.ReadWrite));
        using var f = new CountedFile(dir.File("d.bin"), FileMode.Create, FileAccess.ReadWrite);
        f.WriteAt([1, 2, 3], 0);
        Assert.Throws<ArgumentOutOfRangeException>(() => f.WriteAt([1], -1));
        Assert.Throws<ArgumentOutOfRangeException>(() => f.ReadAt(new byte[1], -1));
    }

    [Test]
    public void ADiskWriteError_IsSurfaced_AsAnIoException()
    {
        if (!File.Exists("/dev/full")) Assert.Ignore("no /dev/full on this platform");
        using var f = new CountedFile("/dev/full", FileMode.Open, FileAccess.Write);
        Assert.Throws<IOException>(() => f.WriteAt(new byte[4096], 0), "ENOSPC");
    }

    [Test]
    public void SharesTheFile_WithAnotherHandle()
    {
        using var dir = new TempDir();
        using var w = new CountedFile(dir.File("e.bin"), FileMode.Create, FileAccess.ReadWrite);
        using var r = new CountedFile(dir.File("e.bin"), FileMode.Open, FileAccess.Read);
        w.WriteAt([5, 6, 7], 0);
        var buf = new byte[3];
        r.ReadExactAt(buf, 0);
        Assert.That(buf, Is.EqualTo(new byte[] { 5, 6, 7 }));
    }
}

/// <summary>The mapped sliding-window writer (default window 1 MiB: E9 found that the .NET view accessor flushes every window on disposal, so 64 KB windows doubled the device writes).</summary>
[TestFixture]
public class MappedWindowWriterTests
{
    [Test]
    public void TheDefaultWindow_IsAtLeastOneMiB_AndIsUsed()
    {
        Assert.That(MappedWindowWriter.DefaultWindowBytes, Is.GreaterThanOrEqualTo(1 << 20));
        using var dir = new TempDir();
        var c = new IoCounters();
        const int length = 3 * (1 << 20) + 17;
        using (var w = new MappedWindowWriter(dir.File("w.bin"), length, counters: c))
        {
            w.Append(new byte[length]);
            Assert.That(w.Position, Is.EqualTo(length));
        }
        Assert.That(c.MappedMaxWindowBytes, Is.EqualTo(MappedWindowWriter.DefaultWindowBytes));
        Assert.That(c.MappedWindows, Is.EqualTo(4), "3 full windows and the 17-byte tail");
        Assert.That(c.MappedWriteBytes, Is.EqualTo(length));
        Assert.That(new FileInfo(dir.File("w.bin")).Length, Is.EqualTo(length));
    }

    [Test]
    public void Append_AcrossWindowsOfEverySize_WritesTheBytes()
    {
        using var dir = new TempDir();
        var rng = new SplitMix64(5);
        foreach (int window in new[] { 1, 2, 3, 7, 16, 100, 1 << 20 })
        {
            foreach (int length in new[] { 0, 1, 2, 15, 16, 17, 100, 333 })
            {
                var data = new byte[length];
                for (int i = 0; i < length; i++) data[i] = (byte)rng.NextInt(256);
                string path = dir.File($"w{window}_{length}.bin");
                var c = new IoCounters();
                using (var w = new MappedWindowWriter(path, length, window, c))
                {
                    int at = 0;
                    while (at < length)
                    {
                        int take = Math.Min(length - at, 1 + rng.NextInt(40));
                        w.Append(data.AsSpan(at, take));
                        at += take;
                    }
                    Assert.That(w.Length, Is.EqualTo(length));
                    Assert.That(w.Path, Is.EqualTo(path));
                }
                Assert.That(File.ReadAllBytes(path), Is.EqualTo(data), $"window {window} length {length}");
                Assert.That(c.MappedWriteBytes, Is.EqualTo(length));
                Assert.That(c.MappedWindows, Is.GreaterThanOrEqualTo(length == 0 ? 0 : (length + window - 1) / window));
                if (length > 0) Assert.That(c.MappedMaxWindowBytes, Is.EqualTo(Math.Min(window, length)));
            }
        }
    }

    [Test]
    public void WriteAt_MayJump_AndOverwrites()
    {
        using var dir = new TempDir();
        string path = dir.File("j.bin");
        using (var w = new MappedWindowWriter(path, 64, 8))
        {
            w.WriteAt(40, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]);
            w.WriteAt(0, [9, 9]);
            w.WriteAt(41, [77]);
            Assert.That(w.Position, Is.Zero, "WriteAt does not move the append position");
            w.Append([5, 5, 5]);
            Assert.That(w.Position, Is.EqualTo(3));
        }
        byte[] file = File.ReadAllBytes(path);
        Assert.That(file.Length, Is.EqualTo(64));
        Assert.That(file[40..52], Is.EqualTo(new byte[] { 1, 77, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 }));
        Assert.That(file[0..3], Is.EqualTo(new byte[] { 5, 5, 5 }));
        Assert.That(file[3..40], Is.All.Zero);
    }

    [Test]
    public void RangesOutsideTheFile_AreRefused_AndBadArgumentsToo()
    {
        using var dir = new TempDir();
        using var w = new MappedWindowWriter(dir.File("r.bin"), 10, 4);
        Assert.Throws<ArgumentOutOfRangeException>(() => w.WriteAt(-1, [1]));
        Assert.Throws<ArgumentOutOfRangeException>(() => w.WriteAt(10, [1]));
        Assert.Throws<ArgumentOutOfRangeException>(() => w.WriteAt(9, [1, 2]));
        w.WriteAt(9, [1]);
        w.WriteAt(10, []); // an empty write at the end is nothing
        w.Append(new byte[10]);
        Assert.Throws<ArgumentOutOfRangeException>(() => w.Append([1]), "appending past the length");
        Assert.Throws<ArgumentOutOfRangeException>(() => _ = new MappedWindowWriter(dir.File("x.bin"), -1, 4));
        Assert.Throws<ArgumentOutOfRangeException>(() => _ = new MappedWindowWriter(dir.File("y.bin"), 10, 0));
    }

    [Test]
    public void ARecreatedFile_IsTruncatedToTheNewLength()
    {
        using var dir = new TempDir();
        string path = dir.File("t.bin");
        File.WriteAllBytes(path, new byte[1000]);
        using (var w = new MappedWindowWriter(path, 10, 4)) w.Append(new byte[10]);
        Assert.That(new FileInfo(path).Length, Is.EqualTo(10));
    }

    [Test]
    public void ATargetInAMissingDirectory_IsSurfaced()
    {
        using var dir = new TempDir();
        Assert.Throws<DirectoryNotFoundException>(() => _ = new MappedWindowWriter(dir.File("no/such/w.bin"), 10, 4));
    }

    [Test]
    public void ADeviceThatCannotBeSized_IsSurfaced_NotHung()
    {
        if (!File.Exists("/dev/full")) Assert.Ignore("no /dev/full on this platform");
        Assert.Catch<Exception>(() =>
        {
            using var w = new MappedWindowWriter("/dev/full", 4096, 4096);
            w.Append(new byte[4096]);
        });
    }

    [Test]
    public void Dispose_IsIdempotent_AndAnEmptyFileHasNoMapping()
    {
        using var dir = new TempDir();
        var w = new MappedWindowWriter(dir.File("z.bin"), 0, 4);
        w.Dispose();
        w.Dispose();
        Assert.That(new FileInfo(dir.File("z.bin")).Length, Is.Zero);
    }
}
