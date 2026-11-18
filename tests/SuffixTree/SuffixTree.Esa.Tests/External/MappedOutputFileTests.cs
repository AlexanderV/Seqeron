using System.Diagnostics;
using SuffixTree.Esa.External;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.External;

/// <summary>
/// The shared output file and its range writers (work package 1b.7b): disjoint ranges written in any order by several writers, offsets above 2^32 bytes on sparse files, completion and fsync, and the reservation of the
/// disk space (preallocation where the platform has it, an explicit free-space check where it has not), through an injectable source.
/// </summary>
[TestFixture]
public class MappedOutputFileTests
{
    private static byte[] Bytes(int n, ulong seed)
    {
        var rng = new SplitMix64(seed);
        var b = new byte[n];
        for (int i = 0; i < n; i++) b[i] = (byte)rng.NextInt(256);
        return b;
    }

    // ---- ranges

    [Test]
    public void RangeWriters_UseOffsetsRelativeToTheirRange_AndTheFileHoldsEveryRange()
    {
        using var dir = new TempDir();
        string path = dir.File("r.bin");
        byte[] a = Bytes(300, 1), b = Bytes(500, 2), c = Bytes(200, 3);
        using (MappedOutputFile f = MappedOutputFile.Create(path, 1000))
        {
            Assert.That(f.Length, Is.EqualTo(1000));
            Assert.That(f.Path, Is.EqualTo(path));
            using MappedWindowWriter wa = f.CreateWriter(0, 300, 64);
            using MappedWindowWriter wb = f.CreateWriter(300, 500, 100);
            using MappedWindowWriter wc = f.CreateWriter(800, 200, 7);
            Assert.That((wa.FileOffset, wa.Length, wa.Position), Is.EqualTo((0L, 300L, 0L)));
            Assert.That((wb.FileOffset, wb.Length, wb.Position), Is.EqualTo((300L, 500L, 0L)));
            Assert.That((wc.FileOffset, wc.Length), Is.EqualTo((800L, 200L)));
            wc.WriteAt(100, c.AsSpan(100)); // random order: the second half first
            wb.Append(b.AsSpan(0, 250));
            wa.Append(a);
            wc.WriteAt(0, c.AsSpan(0, 100));
            wb.Append(b.AsSpan(250));
            Assert.That(wb.Position, Is.EqualTo(500));
            Assert.That(wc.Position, Is.Zero, "WriteAt does not move the append position");
            Assert.That(wa.Path, Is.EqualTo(path));
        }
        Assert.That(File.ReadAllBytes(path), Is.EqualTo(a.Concat(b).Concat(c).ToArray()));
    }

    [Test]
    public void ARangeWriter_RefusesWritesOutsideItsRange_EvenInsideTheFile()
    {
        using var dir = new TempDir();
        using MappedOutputFile f = MappedOutputFile.Create(dir.File("o.bin"), 100);
        using MappedWindowWriter w = f.CreateWriter(40, 20, 8);
        Assert.Throws<ArgumentOutOfRangeException>(() => w.WriteAt(-1, [1]));
        Assert.Throws<ArgumentOutOfRangeException>(() => w.WriteAt(20, [1]));
        Assert.Throws<ArgumentOutOfRangeException>(() => w.WriteAt(19, [1, 2]));
        Assert.Throws<ArgumentOutOfRangeException>(() => w.WriteAt(long.MaxValue, [1]));
        w.WriteAt(19, [1]);
        w.WriteAt(20, []);
        w.Append(new byte[20]);
        Assert.Throws<ArgumentOutOfRangeException>(() => w.Append([1]));
    }

    [Test]
    public void CreateWriter_RefusesARangeOutsideTheFile_AndBadArguments()
    {
        using var dir = new TempDir();
        using MappedOutputFile f = MappedOutputFile.Create(dir.File("o.bin"), 100);
        Assert.Throws<ArgumentOutOfRangeException>(() => f.CreateWriter(-1, 10));
        Assert.Throws<ArgumentOutOfRangeException>(() => f.CreateWriter(0, -1));
        Assert.Throws<ArgumentOutOfRangeException>(() => f.CreateWriter(95, 10));
        Assert.Throws<ArgumentOutOfRangeException>(() => f.CreateWriter(101, 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => f.CreateWriter(0, 10, 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => f.CreateWriter(long.MaxValue, 10));
        f.CreateWriter(100, 0).Dispose();
        f.CreateWriter(90, 10).Dispose();
        f.CreateWriter().Dispose();
        Assert.Throws<ArgumentOutOfRangeException>(() => MappedOutputFile.Create(dir.File("neg.bin"), -1));
        Assert.Throws<ArgumentNullException>(() => MappedOutputFile.Create(null!, 1));
    }

    [Test]
    public void TheWholeFileWriter_HasTheDefaultWindowOfAtLeastOneMiB_AndCountsIt()
    {
        using var dir = new TempDir();
        var c = new IoCounters();
        const int length = (2 * (1 << 20)) + 5;
        using (MappedOutputFile f = MappedOutputFile.Create(dir.File("w.bin"), length))
        {
            using MappedWindowWriter w = f.CreateWriter(counters: c);
            Assert.That((w.FileOffset, w.Length), Is.EqualTo((0L, (long)length)));
            w.Append(new byte[length]);
        }
        Assert.That(MappedWindowWriter.DefaultWindowBytes, Is.GreaterThanOrEqualTo(1 << 20));
        Assert.That(c.MappedMaxWindowBytes, Is.EqualTo(MappedWindowWriter.DefaultWindowBytes));
        Assert.That(c.MappedWindows, Is.EqualTo(3));
        Assert.That(c.MappedWriteBytes, Is.EqualTo(length));
    }

    [Test]
    public void FourThreads_WriteTheirOwnRanges_ThroughOneMapping()
    {
        using var dir = new TempDir();
        string path = dir.File("t.bin");
        const int each = 3_000_001;
        var expected = new byte[4 * each];
        using (MappedOutputFile f = MappedOutputFile.Create(path, expected.Length))
        {
            Parallel.For(0, 4, new ParallelOptions { MaxDegreeOfParallelism = 4 }, t =>
            {
                byte[] data = Bytes(each, (ulong)t);
                data.CopyTo(expected, t * each);
                using MappedWindowWriter w = f.CreateWriter((long)t * each, each, 1 << 12);
                var rng = new SplitMix64((ulong)t + 100);
                int at = 0;
                while (at < each)
                {
                    int take = Math.Min(each - at, 1 + rng.NextInt(20_000));
                    w.Append(data.AsSpan(at, take));
                    at += take;
                }
            });
            f.Complete();
        }
        Assert.That(File.ReadAllBytes(path).AsSpan().SequenceEqual(expected), Is.True);
    }

    [Test]
    public void ADisposedWriter_HasUnmappedItsWindow_WhileTheFileStaysOpen()
    {
        using var dir = new TempDir();
        string path = dir.File("m.bin");
        using MappedOutputFile f = MappedOutputFile.Create(path, 1 << 20);
        MappedWindowWriter w = f.CreateWriter(0, 1 << 20, 4096);
        w.Append([1, 2, 3]);
        bool? open = ProcMaps.Contains(path);
        if (open is null) Assert.Ignore("no /proc/self/maps here");
        Assert.That(open, Is.True, "the window is mapped while the writer is open");
        w.Flush();
        Assert.That(ProcMaps.Contains(path), Is.False, "Flush releases the window");
        w.Append([4]);
        Assert.That(ProcMaps.Contains(path), Is.True);
        w.Dispose();
        Assert.That(ProcMaps.Contains(path), Is.False, "Dispose releases the window");
    }

    [Test]
    public void AnOwningWriter_UnmapsEverything_WhenDisposed()
    {
        using var dir = new TempDir();
        string path = dir.File("o.bin");
        var w = new MappedWindowWriter(path, 1 << 20, 4096);
        w.Append([1, 2, 3]);
        if (ProcMaps.Contains(path) is null) Assert.Ignore("no /proc/self/maps here");
        Assert.That(ProcMaps.Contains(path), Is.True);
        w.Dispose();
        Assert.That(ProcMaps.Contains(path), Is.False);
    }

    [Test]
    public void Flush_ReleasesTheWindow_SoTheNextWriteOpensAnother()
    {
        using var dir = new TempDir();
        var c = new IoCounters();
        using MappedOutputFile f = MappedOutputFile.Create(dir.File("f.bin"), 100);
        using MappedWindowWriter w = f.CreateWriter(0, 100, 100, c);
        w.Append([1, 2, 3]);
        Assert.That(c.MappedWindows, Is.EqualTo(1));
        w.Flush();
        w.Flush();
        w.Append([4]);
        Assert.That(c.MappedWindows, Is.EqualTo(2));
        w.Dispose();
        Assert.That(File.ReadAllBytes(dir.File("f.bin"))[..4], Is.EqualTo(new byte[] { 1, 2, 3, 4 }));
    }

    [Test]
    public void Complete_NeedsEveryWriterToBeDisposed_AndMayBeRepeated_AndSyncsOrNot()
    {
        using var dir = new TempDir();
        MappedOutputFile f = MappedOutputFile.Create(dir.File("c.bin"), 10);
        MappedWindowWriter w = f.CreateWriter();
        w.Append([1, 2, 3]);
        InvalidOperationException? e = Assert.Throws<InvalidOperationException>(() => f.Complete());
        Assert.That(e!.Message, Does.Contain("writer"));
        w.Dispose();
        w.Dispose();
        f.Complete();
        f.Complete(false);
        f.Complete(true);
        MappedWindowWriter w2 = f.CreateWriter(5, 5);
        Assert.Throws<InvalidOperationException>(() => f.Complete(false), "a writer created after the first completion counts as well");
        w2.Dispose();
        f.Dispose();
        f.Dispose();
        Assert.Throws<ObjectDisposedException>(() => f.Complete());
        Assert.Throws<ObjectDisposedException>(() => f.CreateWriter());
        Assert.That(File.ReadAllBytes(dir.File("c.bin"))[..3], Is.EqualTo(new byte[] { 1, 2, 3 }));
        using var exclusive = new FileStream(dir.File("c.bin"), FileMode.Open, FileAccess.ReadWrite, FileShare.None);
    }

    [Test]
    public void AnEmptyFile_HasNoMapping_AndAnEmptyRangeIsFine()
    {
        using var dir = new TempDir();
        using (MappedOutputFile f = MappedOutputFile.Create(dir.File("z.bin"), 0))
        {
            using (MappedWindowWriter w = f.CreateWriter())
            {
                w.Append([]);
                w.WriteAt(0, []);
                Assert.Throws<ArgumentOutOfRangeException>(() => w.Append([1]));
            }
            f.Complete();
        }
        Assert.That(new FileInfo(dir.File("z.bin")).Length, Is.Zero);
    }

    [Test]
    public void ATargetInAMissingDirectory_IsSurfaced_AndNothingIsLeft()
    {
        using var dir = new TempDir();
        Assert.Throws<DirectoryNotFoundException>(() => MappedOutputFile.Create(dir.File("no/such/w.bin"), 10));
        Assert.That(Directory.GetFileSystemEntries(dir.Path), Is.Empty);
    }

    [Test]
    public void ADeviceThatCannotBeSized_IsSurfaced_AndTheLinkIsNotLeftBehind()
    {
        if (!File.Exists("/dev/full") || !OperatingSystem.IsLinux()) Assert.Ignore("no /dev/full on this platform");
        foreach (ReservationMode mode in new[] { ReservationMode.Sparse, ReservationMode.Preallocate, ReservationMode.Auto })
        {
            using var dir = new TempDir();
            string link = dir.File("full.bin");
            File.CreateSymbolicLink(link, "/dev/full");
            Assert.Throws<IOException>(() => MappedOutputFile.Create(link, 4096, new MappedFileOptions { Reservation = mode }), mode.ToString());
            Assert.That(File.Exists("/dev/full"), Is.True);
            Assert.That(Directory.GetFileSystemEntries(dir.Path), Is.Empty, $"{mode}: a failed creation leaves no file");
        }
    }

    // ---- offsets above 2^32 on sparse files

    [Test]
    public void RangeWriters_AtFileOffsets_AboveTwoToThe32_OnASparseFile()
    {
        using var dir = new TempDir();
        string path = dir.File("big.bin");
        const long length = 5_000_000_000L;
        using (MappedOutputFile f = MappedOutputFile.Create(path, length, new MappedFileOptions { Reservation = ReservationMode.Sparse }))
        {
            using (MappedWindowWriter hi = f.CreateWriter(4_900_000_000L, 100_000_000L, 1 << 20))
            {
                hi.WriteAt(0, [1, 2, 3, 4]);
                hi.WriteAt(99_999_998L, [7, 7]);
            }
            using (MappedWindowWriter straddle = f.CreateWriter(4_294_967_290L, 100))
            {
                straddle.WriteAt(0, [9, 8, 7, 6, 5, 4, 3, 2, 1, 0, 11, 12]); // crosses 2^32 - 1 and 2^32
            }
        }
        Assert.That(new FileInfo(path).Length, Is.EqualTo(length));
        using var fs = new FileStream(path, FileMode.Open, FileAccess.Read);
        var buf = new byte[12];
        fs.Position = 4_900_000_000L;
        fs.ReadExactly(buf.AsSpan(0, 4));
        Assert.That(buf[..4], Is.EqualTo(new byte[] { 1, 2, 3, 4 }));
        fs.Position = 4_999_999_998L;
        fs.ReadExactly(buf.AsSpan(0, 2));
        Assert.That(buf[..2], Is.EqualTo(new byte[] { 7, 7 }));
        fs.Position = 4_294_967_290L;
        fs.ReadExactly(buf);
        Assert.That(buf, Is.EqualTo(new byte[] { 9, 8, 7, 6, 5, 4, 3, 2, 1, 0, 11, 12 }));
        File.Delete(path);
    }

    // ---- the reservation

    [Test]
    public void TheDefaultMode_IsAuto_AndAutoResolvesToTheMeasuredChoice()
    {
        var o = new MappedFileOptions();
        Assert.That(o.Reservation, Is.EqualTo(ReservationMode.Auto));
        Assert.That(o.Source, Is.Null);
        Assert.That(MappedOutputFile.AutoMode, Is.EqualTo(ReservationMode.Preallocate), "see ADR-004 section 5a, 1b.7b, for the measurement");
        using var dir = new TempDir();
        using MappedOutputFile f = MappedOutputFile.Create(dir.File("a.bin"), 1000, new MappedFileOptions { Source = new FakeReservationSource { CanPreallocate = true } });
        Assert.That(f.Outcome, Is.EqualTo(ReservationOutcome.Preallocated));
        Assert.That(new FileInfo(dir.File("a.bin")).Length, Is.EqualTo(1000));
    }

    [Test]
    public void Sparse_CreatesTheFileAtItsLength_AndTouchesNeitherFallocateNorTheFreeSpace()
    {
        using var dir = new TempDir();
        var source = new FakeReservationSource { CanPreallocate = true, Free = 0 };
        using MappedOutputFile f = MappedOutputFile.Create(dir.File("s.bin"), 5000, new MappedFileOptions { Reservation = ReservationMode.Sparse, Source = source });
        Assert.That(f.Outcome, Is.EqualTo(ReservationOutcome.Sparse));
        Assert.That((source.PreallocateCalls, source.FreeCalls), Is.EqualTo((0, 0)));
        Assert.That(new FileInfo(dir.File("s.bin")).Length, Is.EqualTo(5000));
    }

    [Test]
    public void Preallocate_WherePossible_CallsTheSourceOnce_AndTheFileHasItsLength()
    {
        using var dir = new TempDir();
        var source = new FakeReservationSource { CanPreallocate = true };
        using MappedOutputFile f = MappedOutputFile.Create(dir.File("p.bin"), 7777, new MappedFileOptions { Reservation = ReservationMode.Preallocate, Source = source });
        Assert.That(f.Outcome, Is.EqualTo(ReservationOutcome.Preallocated));
        Assert.That(source.PreallocateCalls, Is.EqualTo(1));
        Assert.That(source.FreeCalls, Is.Zero);
        Assert.That(new FileInfo(dir.File("p.bin")).Length, Is.EqualTo(7777));
        using MappedWindowWriter w = f.CreateWriter();
        w.Append(new byte[7777]);
    }

    [Test]
    public void Preallocate_WhereThePlatformHasNone_ChecksTheFreeSpace_ExplicitlyAndThenCreatesSparse()
    {
        using var dir = new TempDir();
        string path = dir.File("n.bin");
        var enough = new FakeReservationSource { CanPreallocate = false, Free = 5000 };
        using (MappedOutputFile f = MappedOutputFile.Create(path, 5000, new MappedFileOptions { Reservation = ReservationMode.Preallocate, Source = enough }))
        {
            Assert.That(f.Outcome, Is.EqualTo(ReservationOutcome.FreeSpaceChecked));
        }
        Assert.That(enough.LastDirectory, Is.EqualTo(Path.GetDirectoryName(Path.GetFullPath(path))), "the directory of the file is asked, the file does not exist yet");
        var tight = new FakeReservationSource { CanPreallocate = false, Free = 4999 };
        IOException? e = Assert.Throws<IOException>(() => MappedOutputFile.Create(path, 5000, new MappedFileOptions { Reservation = ReservationMode.Preallocate, Source = tight }));
        Assert.That(e!.Message, Does.Contain("5000").And.Contain("4999"));
        Assert.That(File.Exists(path), Is.False, "a refused creation removes the file it would have replaced");
        var unknown = new FakeReservationSource { CanPreallocate = false, Free = -1 };
        using (MappedOutputFile f = MappedOutputFile.Create(path, 5000, new MappedFileOptions { Reservation = ReservationMode.Preallocate, Source = unknown }))
        {
            Assert.That(f.Outcome, Is.EqualTo(ReservationOutcome.Unchecked), "free space unknown: sparse, said so");
        }
        var empty = new FakeReservationSource { CanPreallocate = false, Free = 0 };
        using (MappedOutputFile f = MappedOutputFile.Create(path, 0, new MappedFileOptions { Reservation = ReservationMode.Preallocate, Source = empty }))
        {
            Assert.That(f.Outcome, Is.EqualTo(ReservationOutcome.FreeSpaceChecked), "a file of length 0 needs no space");
        }
    }

    [Test]
    public void ANegativeLength_IsRefusedBeforeAnythingIsTouched()
    {
        using var dir = new TempDir();
        string path = dir.File("keep.bin");
        File.WriteAllBytes(path, [1, 2, 3]);
        Assert.Throws<ArgumentOutOfRangeException>(() => MappedOutputFile.Create(path, -1));
        Assert.That(File.ReadAllBytes(path), Is.EqualTo(new byte[] { 1, 2, 3 }), "the file that was there is still there");
        Assert.Throws<ArgumentNullException>(() => MappedOutputFile.Create(null!, 10));
    }

    [Test]
    public void ARecreatedFile_IsEmptyOfTheOldContent_InEveryMode()
    {
        using var dir = new TempDir();
        string path = dir.File("old.bin");
        foreach (ReservationMode mode in new[] { ReservationMode.Sparse, ReservationMode.Preallocate })
        {
            foreach (bool system in new[] { false, true })
            {
                File.WriteAllBytes(path, Enumerable.Repeat((byte)0xAB, 3000).ToArray());
                var options = new MappedFileOptions { Reservation = mode, Source = system ? null : new FakeReservationSource { CanPreallocate = true } };
                using (MappedOutputFile.Create(path, 2000, options)) { }
                Assert.That(File.ReadAllBytes(path), Is.EqualTo(new byte[2000]), $"{mode} system={system}");
            }
            if (mode == ReservationMode.Preallocate)
            {
                File.WriteAllBytes(path, Enumerable.Repeat((byte)0xAB, 3000).ToArray());
                using (MappedOutputFile.Create(path, 2000, new MappedFileOptions { Reservation = mode, Source = new FakeReservationSource { CanPreallocate = false, Free = 1 << 20 } })) { }
                Assert.That(File.ReadAllBytes(path), Is.EqualTo(new byte[2000]), "the explicit fallback creates the file anew as well");
            }
        }
    }

    [Test]
    public void APreallocationThatFails_IsSurfaced_AndNoFileIsLeft()
    {
        using var dir = new TempDir();
        var source = new FakeReservationSource { CanPreallocate = true, PreallocateThrows = new IOException("Failed to create because the disk was full.") };
        Assert.Throws<IOException>(() => MappedOutputFile.Create(dir.File("x.bin"), 100, new MappedFileOptions { Reservation = ReservationMode.Preallocate, Source = source }));
        Assert.That(Directory.GetFileSystemEntries(dir.Path), Is.Empty);
        File.WriteAllBytes(dir.File("old.bin"), [1, 2, 3]);
        Assert.Throws<IOException>(() => MappedOutputFile.Create(dir.File("old.bin"), 100, new MappedFileOptions { Reservation = ReservationMode.Preallocate, Source = source }));
        Assert.That(File.Exists(dir.File("old.bin")), Is.False, "the file that was to be replaced is gone: it was truncated by the creation");
    }

    [Test]
    public void TheSystemSource_OnLinux_ReservesTheBlocks_AndTheSparseFileDoesNot()
    {
        if (!OperatingSystem.IsLinux()) Assert.Ignore("the system source preallocates on Linux only");
        using var dir = new TempDir();
        const long length = 256L << 20;
        long free0 = SystemReservationSource.Instance.AvailableFreeBytes(dir.Path);
        if (free0 < 0) Assert.Ignore("free space unknown here");
        Assert.That(SystemReservationSource.Instance.CanPreallocate, Is.True);
        using (MappedOutputFile f = MappedOutputFile.Create(dir.File("pre.bin"), length, new MappedFileOptions { Reservation = ReservationMode.Preallocate }))
        {
            Assert.That(f.Outcome, Is.EqualTo(ReservationOutcome.Preallocated));
            long free1 = SystemReservationSource.Instance.AvailableFreeBytes(dir.Path);
            if (free0 - free1 < 16L << 20) Assert.Ignore($"this file system does not reserve blocks (free space fell by {free0 - free1} bytes)");
            Assert.That(free0 - free1, Is.GreaterThan(length * 9 / 10), "the blocks are taken at creation");
            Assert.That(new FileInfo(dir.File("pre.bin")).Length, Is.EqualTo(length));
        }
        File.Delete(dir.File("pre.bin"));
        long free2 = SystemReservationSource.Instance.AvailableFreeBytes(dir.Path);
        using (var legacy = new MappedWindowWriter(dir.File("legacy.bin"), length))
        {
            long free4 = SystemReservationSource.Instance.AvailableFreeBytes(dir.Path);
            Assert.That(free2 - free4, Is.LessThan(length / 10), "the constructor with a path keeps creating sparse files (Dist is unchanged)");
            Assert.That(legacy.Length, Is.EqualTo(length));
        }
        File.Delete(dir.File("legacy.bin"));
        free2 = SystemReservationSource.Instance.AvailableFreeBytes(dir.Path);
        using (MappedOutputFile f = MappedOutputFile.Create(dir.File("sp.bin"), length, new MappedFileOptions { Reservation = ReservationMode.Sparse }))
        {
            long free3 = SystemReservationSource.Instance.AvailableFreeBytes(dir.Path);
            Assert.That(free2 - free3, Is.LessThan(length / 10), "a sparse file takes no blocks until they are written");
            Assert.That(f.Outcome, Is.EqualTo(ReservationOutcome.Sparse));
        }
    }

    [Test]
    public void TheSystemSource_ReturnsUnknown_ForADirectoryThatDoesNotExist_AndNeverThrows()
    {
        Assert.That(SystemReservationSource.Instance.AvailableFreeBytes("/no/such/directory/at/all"), Is.EqualTo(-1));
        Assert.That(SystemReservationSource.Instance.AvailableFreeBytes(string.Empty), Is.EqualTo(-1));
        Assert.That(SystemReservationSource.Instance.AvailableFreeBytes(Path.GetTempPath()), Is.GreaterThan(0).Or.EqualTo(-1));
    }

    [Test]
    public void ARealFullFileSystem_RefusesAPreallocatedFileUpFront_WithAnIoException()
    {
        if (!OperatingSystem.IsLinux()) Assert.Ignore("needs Linux to mount a small tmpfs");
        using var dir = new TempDir();
        string mnt = dir.Sub("mnt");
        Directory.CreateDirectory(mnt);
        if (!TryMountTmpfs(mnt, 4)) Assert.Ignore("could not mount a tmpfs here");
        try
        {
            string path = Path.Combine(mnt, "big.bin");
            IOException? e = Assert.Throws<IOException>(() => MappedOutputFile.Create(path, 32L << 20, new MappedFileOptions { Reservation = ReservationMode.Preallocate }));
            Assert.That(e!.Message, Does.Contain("full").IgnoreCase.Or.Contain("space").IgnoreCase);
            Assert.That(File.Exists(path), Is.False);
            string path2 = Path.Combine(mnt, "fits.bin");
            using (MappedOutputFile f = MappedOutputFile.Create(path2, 1L << 20, new MappedFileOptions { Reservation = ReservationMode.Preallocate }))
            {
                Assert.That(f.Outcome, Is.EqualTo(ReservationOutcome.Preallocated));
            }
            // the free-space check of the fallback sees the same: a file that does not fit is refused before it is created
            var noFallocate = new FakeReservationSource { CanPreallocate = false, Free = SystemReservationSource.Instance.AvailableFreeBytes(mnt) };
            Assert.Throws<IOException>(() => MappedOutputFile.Create(Path.Combine(mnt, "big2.bin"), 32L << 20, new MappedFileOptions { Reservation = ReservationMode.Preallocate, Source = noFallocate }));
        }
        finally
        {
            Unmount(mnt);
        }
    }

    private static bool TryMountTmpfs(string dir, int megabytes)
    {
        try
        {
            using Process? p = Process.Start(new ProcessStartInfo("mount", ["-t", "tmpfs", "-o", $"size={megabytes}m", "tmpfs", dir]) { RedirectStandardError = true, RedirectStandardOutput = true });
            p!.WaitForExit(10_000);
            return p.ExitCode == 0;
        }
        catch (System.ComponentModel.Win32Exception)
        {
            return false;
        }
    }

    private static void Unmount(string dir)
    {
        using Process? p = Process.Start(new ProcessStartInfo("umount", [dir]) { RedirectStandardError = true, RedirectStandardOutput = true });
        p!.WaitForExit(10_000);
    }
}
