using Microsoft.Win32.SafeHandles;
using SuffixTree.Esa.External;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.External;

/// <summary>
/// A source that records the reservation of a range that was added to a file and can fail inside it like a disk that fills up: the range is written as zeros (so that the partial work of a failed call is real) and
/// the call fails after <see cref="FailAfterBytes"/> bytes of it.
/// </summary>
internal sealed class TailSource : IReservationSource
{
    /// <summary>What <see cref="CanPreallocate"/> says.</summary>
    internal bool Can { get; init; }

    /// <summary>Free bytes reported (-1: unknown).</summary>
    internal long Free { get; set; } = -1;

    /// <summary>A call whose range is longer than this writes that many bytes and then throws an <see cref="IOException"/>.</summary>
    internal long FailAfterBytes { get; set; } = long.MaxValue;

    /// <summary>The source records the range and returns without writing (it reserves by another means).</summary>
    internal bool WritesNothing { get; init; }

    /// <summary>Every range that was asked for, in order.</summary>
    internal List<(long Offset, long Length)> Reserved { get; } = [];

    public bool CanPreallocate => Can;

    public FileStream CreatePreallocated(string path, long length) => new(path, FileMode.Create, FileAccess.ReadWrite, FileShare.ReadWrite);

    public long AvailableFreeBytes(string directory) => Free;

    public void ReserveRange(SafeFileHandle file, long offset, long length, IoCounters? counters)
    {
        lock (Reserved) Reserved.Add((offset, length));
        if (WritesNothing) return;
        if (FailAfterBytes < length)
        {
            ZeroFill.Write(file, offset, FailAfterBytes, counters);
            throw new IOException($"disk full (injected {FailAfterBytes} bytes into the reservation of {length})");
        }
        ZeroFill.Write(file, offset, length, counters);
    }
}

/// <summary>A source whose reservation fails with an exception that is not an <see cref="IOException"/>.</summary>
internal sealed class ThrowingSource : IReservationSource
{
    public bool CanPreallocate => true;

    public FileStream CreatePreallocated(string path, long length) => new(path, FileMode.Create, FileAccess.ReadWrite, FileShare.ReadWrite);

    public long AvailableFreeBytes(string directory) => -1;

    public void ReserveRange(SafeFileHandle file, long offset, long length, IoCounters? counters) => throw new InvalidOperationException("the reservation failed in a way that is not an IOException");
}

/// <summary>
/// A source that does what the system does, except that it cannot tell the free space: the pre-flights of the builder (the early one and the late one) pass whatever the volume holds, so that only the volume itself
/// can refuse a file or a write (a disk that fills after the pre-flight, which is the window of work package 1b.9d).
/// </summary>
internal sealed class BlindSource : IReservationSource
{
    public bool CanPreallocate => SystemReservationSource.Instance.CanPreallocate;

    public FileStream CreatePreallocated(string path, long length) => SystemReservationSource.Instance.CreatePreallocated(path, length);

    public long AvailableFreeBytes(string directory) => -1;
}

/// <summary>
/// The reservation of what is added to a file (work package 1b.9d): <see cref="MappedOutputFile.Grow"/> reserves the blocks of the new range before it maps it, a failure of the reservation leaves the file as it
/// was, the mode <see cref="ReservationMode.Sparse"/> still gives the sparse growth, and <see cref="Dist"/> no longer creates its target sparse.
/// </summary>
[TestFixture]
public sealed class Reservation1b9dTests
{
    private static MappedFileOptions Opts(TailSource source, ReservationMode mode = ReservationMode.Auto) => new() { Source = source, Reservation = mode };

    // ------------------------------------------------------------------ the range goes through the source

    [Test]
    public void Grow_ReservesTheNewRange_ThroughTheSource_AndSaysFilled()
    {
        using var dir = new TempDir();
        string path = dir.File("g.bin");
        foreach (bool can in new[] { false, true })
        {
            File.Delete(path);
            var source = new TailSource { Can = can, Free = can ? 1L << 40 : -1 };
            var counters = new IoCounters();
            using MappedOutputFile file = MappedOutputFile.Create(path, 1000, Opts(source));
            using (MappedWindowWriter w = file.CreateWriter(0, 1000, 64)) w.Append(Enumerable.Repeat((byte)0x5A, 1000).ToArray());
            Assert.That(file.Grow(5000, Opts(source), counters), Is.EqualTo(ReservationOutcome.Filled), $"can {can}");
            Assert.That(source.Reserved, Is.EqualTo(new[] { (1000L, 4000L) }), "exactly the range that was added");
            Assert.That(counters.ReservedWriteBytes, Is.EqualTo(4000));
            Assert.That(counters.WriteBytes, Is.Zero, "the zeros are reported on their own");
            Assert.That((file.Length, new FileInfo(path).Length), Is.EqualTo((5000L, 5000L)));
            Assert.That(file.Grow(9000, Opts(source), counters), Is.EqualTo(ReservationOutcome.Filled));
            Assert.That(source.Reserved, Is.EqualTo(new[] { (1000L, 4000L), (5000L, 4000L) }), "and again for the next growth");
            Assert.That(counters.ReservedWriteBytes, Is.EqualTo(8000));
            file.Complete(false);
            byte[] bytes = File.ReadAllBytes(path);
            Assert.That(bytes.AsSpan(0, 1000).ToArray(), Is.All.EqualTo((byte)0x5A), "what was there stays");
            Assert.That(bytes.AsSpan(1000).ToArray(), Is.All.Zero, "what is added is zero");
        }
    }

    [Test]
    public void Grow_ExtendsTheFile_EvenWhenTheSourceWritesNothing()
    {
        // a source that reserves by other means (a platform call that a later version of the BCL may offer) is not obliged to write: the length is the caller's business, not the source's
        using var dir = new TempDir();
        var source = new TailSource { Can = true, Free = -1, WritesNothing = true };
        using MappedOutputFile file = MappedOutputFile.Create(dir.File("g.bin"), 1000, Opts(source));
        Assert.That(file.Grow(5000, Opts(source)), Is.EqualTo(ReservationOutcome.Filled));
        Assert.That(source.Reserved, Is.EqualTo(new[] { (1000L, 4000L) }));
        Assert.That((file.Length, new FileInfo(file.Path).Length), Is.EqualTo((5000L, 5000L)));
        using (MappedWindowWriter w = file.CreateWriter(1000, 4000, 4096)) w.Append(new byte[4000]);
    }

    [Test]
    public void Grow_InTheExplicitPreallocateMode_Reserves_AsTheAutomaticModeDoes()
    {
        using var dir = new TempDir();
        var source = new TailSource { Can = false, Free = -1 };
        using MappedOutputFile file = MappedOutputFile.Create(dir.File("g.bin"), 1000, Opts(source, ReservationMode.Preallocate));
        Assert.That(file.Grow(3000, Opts(source, ReservationMode.Preallocate)), Is.EqualTo(ReservationOutcome.Filled));
        Assert.That(source.Reserved, Is.EqualTo(new[] { (1000L, 2000L) }));
    }

    [Test]
    public void AReservationThatFailsInAnyWay_LeavesTheFileAsItWas()
    {
        using var dir = new TempDir();
        var source = new ThrowingSource();
        using MappedOutputFile file = MappedOutputFile.Create(dir.File("g.bin"), 1000, new MappedFileOptions { Source = source });
        Assert.Throws<InvalidOperationException>(() => file.Grow(9000, new MappedFileOptions { Source = source }));
        Assert.That((file.Length, new FileInfo(file.Path).Length), Is.EqualTo((1000L, 1000L)), "whatever the exception, the length is given back");
    }

    [Test]
    public void ZeroFill_OnAFullVolume_CountsOnlyTheBytesThatWereWritten()
    {
        using var dir = new TempDir();
        using SmallVolume volume = SmallVolume.MountOrIgnore(dir, "mnt", 3L << 20);
        var counters = new IoCounters();
        using SafeFileHandle h = File.OpenHandle(volume.File("z.bin"), FileMode.Create, FileAccess.ReadWrite, FileShare.ReadWrite);
        Assert.Throws<IOException>(() => ZeroFill.Write(h, 0, 8L << 20, counters));
        Assert.That(counters.ReservedWriteBytes, Is.EqualTo(3L << 20), "the three chunks that fit were written; the fourth met the end of the volume and is not counted");
    }

    [Test]
    public void Grow_WithNoGrowth_ReservesNothing()
    {
        using var dir = new TempDir();
        var source = new TailSource { Can = true, Free = 1L << 40 };
        using MappedOutputFile file = MappedOutputFile.Create(dir.File("g.bin"), 1000, Opts(source));
        Assert.That(file.Grow(1000, Opts(source)), Is.EqualTo(ReservationOutcome.Sparse), "no growth, no check, no reservation");
        Assert.That(source.Reserved, Is.Empty);
    }

    [Test]
    public void Grow_ThroughTheDefaultOfTheInterface_WritesTheZerosInChunksOfOneMiB_AndCountsEveryByte()
    {
        using var dir = new TempDir();
        string path = dir.File("g.bin");
        var source = new FakeReservationSource { CanPreallocate = false, Free = -1 }; // no override of ReserveRange: the default of the interface runs
        var counters = new IoCounters();
        const long growth = (2 * ZeroFill.ChunkBytes) + 12_345;
        using MappedOutputFile file = MappedOutputFile.Create(path, 4096, new MappedFileOptions { Source = source });
        using (MappedWindowWriter w = file.CreateWriter(0, 4096, 4096)) w.Append(Enumerable.Repeat((byte)0xC3, 4096).ToArray());
        Assert.That(file.Grow(4096 + growth, new MappedFileOptions { Source = source }, counters), Is.EqualTo(ReservationOutcome.Filled));
        Assert.That(ZeroFill.ChunkBytes, Is.EqualTo(1 << 20));
        Assert.That(counters.ReservedWriteBytes, Is.EqualTo(growth));
        file.Complete(false);
        byte[] bytes = File.ReadAllBytes(path);
        Assert.That(bytes.Length, Is.EqualTo(4096 + growth));
        Assert.That(bytes.AsSpan(0, 4096).ToArray(), Is.All.EqualTo((byte)0xC3));
        Assert.That(bytes.AsSpan(4096).ToArray(), Is.All.Zero);
    }

    [Test]
    public void ZeroFill_WritesExactlyTheRange_NothingBeforeIt_AndNothingBehindIt()
    {
        using var dir = new TempDir();
        string path = dir.File("z.bin");
        byte[] before = Enumerable.Range(0, 3_000_000).Select(i => (byte)(i % 251 + 1)).ToArray();
        File.WriteAllBytes(path, before);
        foreach ((long offset, long length) in new[] { (0L, 0L), (7L, 1L), (1000L, (long)ZeroFill.ChunkBytes), (1000L, ZeroFill.ChunkBytes + 1L), (500_000L, 2_000_001L) })
        {
            File.WriteAllBytes(path, before);
            var counters = new IoCounters();
            using (SafeFileHandle h = File.OpenHandle(path, FileMode.Open, FileAccess.ReadWrite, FileShare.ReadWrite)) ZeroFill.Write(h, offset, length, counters);
            byte[] after = File.ReadAllBytes(path);
            Assert.That(after.Length, Is.EqualTo(before.Length));
            for (int i = 0; i < after.Length; i++)
            {
                byte expected = i >= offset && i < offset + length ? (byte)0 : before[i];
                if (after[i] != expected) Assert.Fail($"offset {offset} length {length}: byte {i} is {after[i]}, expected {expected}");
            }
            Assert.That(counters.ReservedWriteBytes, Is.EqualTo(length));
        }
        using (SafeFileHandle h = File.OpenHandle(path, FileMode.Open, FileAccess.ReadWrite, FileShare.ReadWrite))
        {
            Assert.That(Assert.Throws<ArgumentOutOfRangeException>(() => ZeroFill.Write(h, -1, 1, null))!.ParamName, Is.EqualTo("offset"), "refused by ZeroFill, not by the call below it");
            Assert.That(Assert.Throws<ArgumentOutOfRangeException>(() => ZeroFill.Write(h, 0, -1, null))!.ParamName, Is.EqualTo("length"));
            Assert.That(Assert.Throws<ArgumentNullException>(() => ZeroFill.Write(null!, 0, 1, null))!.ParamName, Is.EqualTo("file"));
            ZeroFill.Write(h, 4_000_000, 0, null);
            Assert.That(RandomAccess.GetLength(h), Is.EqualTo(3_000_000), "no bytes, no write, nothing behind the end");
            ZeroFill.Write(h, 3_000_000, 10, null); // a range behind the end extends the file with zeros
        }
        Assert.That(new FileInfo(path).Length, Is.EqualTo(3_000_010));
    }

    [Test]
    public void Grow_InTheSparseMode_StaysSparse_AsBefore()
    {
        using var dir = new TempDir();
        var source = new TailSource { Can = true, Free = 1L << 40 };
        var counters = new IoCounters();
        using MappedOutputFile file = MappedOutputFile.Create(dir.File("g.bin"), 1000, Opts(source, ReservationMode.Sparse));
        Assert.That(file.Grow(5_000_000, Opts(source, ReservationMode.Sparse), counters), Is.EqualTo(ReservationOutcome.FreeSpaceChecked), "the check is made, nothing is reserved");
        Assert.That(source.Reserved, Is.Empty);
        Assert.That(counters.ReservedWriteBytes, Is.Zero);
        Assert.That(new FileInfo(file.Path).Length, Is.EqualTo(5_000_000));
        source.Free = -1;
        Assert.That(file.Grow(6_000_000, Opts(source, ReservationMode.Sparse), counters), Is.EqualTo(ReservationOutcome.Unchecked));
        Assert.That(source.Reserved, Is.Empty);
    }

    [Test]
    public void TheReservedBytes_AreCountedOnTheirOwn_AndReset()
    {
        var c = new IoCounters();
        c.AddReserved(5);
        c.AddReserved(7);
        Assert.That(c.ReservedWriteBytes, Is.EqualTo(12));
        Assert.That((c.WriteBytes, c.WriteCalls, c.MappedWriteBytes, c.Snapshot().TotalBytes), Is.EqualTo((0L, 0L, 0L, 0L)), "they are not explicit writes, not mapped writes and not part of the snapshot");
        c.AddWrite(3);
        Assert.That((c.ReservedWriteBytes, c.WriteBytes), Is.EqualTo((12L, 3L)));
        c.Reset();
        Assert.That((c.ReservedWriteBytes, c.WriteBytes), Is.EqualTo((0L, 0L)));
    }

    // ------------------------------------------------------------------ a failure inside the reservation

    [Test]
    public void AReservationThatFails_AtAnyByte_IsAnIoException_AndTheFileIsAsItWas()
    {
        using var dir = new TempDir();
        string path = dir.File("g.bin");
        const long growth = (2 * ZeroFill.ChunkBytes) + 500;
        foreach (long failAfter in new long[] { 0, 1, ZeroFill.ChunkBytes - 1, ZeroFill.ChunkBytes, ZeroFill.ChunkBytes + 1, growth - 1 })
        {
            var source = new TailSource { Can = true, Free = 1L << 40, FailAfterBytes = failAfter };
            var counters = new IoCounters();
            using MappedOutputFile file = MappedOutputFile.Create(path, 1000, Opts(source));
            using (MappedWindowWriter w = file.CreateWriter(0, 1000, 64)) w.Append(Enumerable.Repeat((byte)0x77, 1000).ToArray());
            IOException? e = Assert.Throws<IOException>(() => file.Grow(1000 + growth, Opts(source), counters), $"failure after {failAfter}");
            Assert.That(e!.Message, Does.Contain("disk full"));
            Assert.That((file.Length, new FileInfo(path).Length), Is.EqualTo((1000L, 1000L)), $"after {failAfter}: the length is given back");
            Assert.That(counters.ReservedWriteBytes, Is.EqualTo(failAfter), "what was written before the failure is counted");
            // the file is usable and the volume's blocks of the failed range are gone with the length: a good source reserves the same range now
            source.FailAfterBytes = long.MaxValue;
            Assert.That(file.Grow(1000 + growth, Opts(source)), Is.EqualTo(ReservationOutcome.Filled));
            Assert.That(source.Reserved, Is.EqualTo(new[] { (1000L, growth), (1000L, growth) }));
            using (MappedWindowWriter w = file.CreateWriter(1000, growth, 4096)) w.Append(new byte[growth]);
            file.Complete(false);
            Assert.That(File.ReadAllBytes(path).AsSpan(0, 1000).ToArray(), Is.All.EqualTo((byte)0x77));
            File.Delete(path);
        }
    }

    [Test]
    public void AReservationThatFails_LeavesNoMappingOfTheFileThatIsNotTheOldOne()
    {
        using var dir = new TempDir();
        var source = new TailSource { Can = true, Free = 1L << 40, FailAfterBytes = 0 };
        MappedOutputFile file = MappedOutputFile.Create(dir.File("g.bin"), 4096, Opts(source));
        Assert.Throws<IOException>(() => file.Grow(8192, Opts(source)));
        using (MappedWindowWriter w = file.CreateWriter(0, 4096, 4096)) w.Append(new byte[4096]); // the old mapping still serves
        file.Dispose();
        Assert.That(ProcMaps.Contains(dir.File("g.bin")) ?? false, Is.False, "nothing stays mapped");
    }

    [Test]
    public void TheFreeSpaceIsCheckedBefore_AnyZeroIsWritten()
    {
        using var dir = new TempDir();
        var source = new TailSource { Can = true, Free = 1999 };
        var counters = new IoCounters();
        using MappedOutputFile file = MappedOutputFile.Create(dir.File("g.bin"), 1000, Opts(source));
        var e = Assert.Throws<IOException>(() => file.Grow(3000, Opts(source), counters));
        Assert.That(e!.Message, Does.Contain("free space"));
        Assert.That(source.Reserved, Is.Empty, "the check refuses first, the zeros are not written");
        Assert.That(counters.ReservedWriteBytes, Is.Zero);
    }

    // ------------------------------------------------------------------ a real volume

    [Test]
    public void OnARealVolume_TheGrowthTakesItsBlocksAtOnce()
    {
        using var dir = new TempDir();
        using SmallVolume volume = SmallVolume.MountOrIgnore(dir, "mnt", 16L << 20);
        using MappedOutputFile file = MappedOutputFile.Create(volume.File("g.bin"), 1L << 20);
        long before = volume.Free;
        Assert.That(file.Grow(9L << 20), Is.EqualTo(ReservationOutcome.Filled));
        Assert.That(before - volume.Free, Is.GreaterThanOrEqualTo(8L << 20), "the 8 MiB that were added are on the volume's books now, not when somebody writes them");
        // and it is the page of the mapping that the disk could not have given later: fill the rest of the volume, then write the whole range
        volume.FillToTheBrim();
        using (MappedWindowWriter w = file.CreateWriter(0, file.Length, 1 << 20)) w.Append(new byte[file.Length]);
        file.Complete(false);
    }

    [Test]
    public void OnARealVolume_AGrowthThatDoesNotFit_IsAnIoException_EvenWhenNobodyKnowsTheFreeSpace_AndGivesTheBlocksBack()
    {
        using var dir = new TempDir();
        using SmallVolume volume = SmallVolume.MountOrIgnore(dir, "mnt", 8L << 20);
        var blind = new BlindSource();
        using MappedOutputFile file = MappedOutputFile.Create(volume.File("g.bin"), 1L << 20, new MappedFileOptions { Source = blind });
        long before = volume.Free;
        IOException? e = Assert.Throws<IOException>(() => file.Grow(40L << 20, new MappedFileOptions { Source = blind }));
        Assert.That(e!.Message, Does.Contain("space").IgnoreCase.Or.Contain("full").IgnoreCase);
        Assert.That((file.Length, new FileInfo(volume.File("g.bin")).Length), Is.EqualTo((1L << 20, 1L << 20)));
        Assert.That(volume.Free, Is.GreaterThanOrEqualTo(before - 4096), "the blocks that were written before the disk filled are given back with the length");
    }

    // ------------------------------------------------------------------ Dist

    private static (string Sa, int N) WritePermutation(TempDir dir, int n)
    {
        string sa = dir.File("sa.bin");
        var rng = new SplitMix64(9);
        FileBits.WriteUInts(sa, [.. Enumerable.Range(0, n + 1).Select(i => (uint)i).OrderBy(_ => rng.Next())]);
        return (sa, n);
    }

    [Test]
    public void Dist_WithNoFileOptions_CreatesItsTargetWithAReservation_AndTheSparseModeStillGivesASparseOne()
    {
        using var dir = new TempDir();
        (string sa, int n) = WritePermutation(dir, 5000);
        DistOptions Options(MappedFileOptions? files) => new() { ScratchDir = dir.Sub("scratch"), ChunkEntries = 700, Files = files };
        DistStats stats = Dist.BuildIsa(sa, n, dir.File("isa.bin"), Options(null));
        Assert.That(stats.Reservation, Is.Not.EqualTo(ReservationOutcome.Sparse), "the default is Auto");
        Assert.That(stats.Reservation, Is.EqualTo(OperatingSystem.IsLinux() ? ReservationOutcome.Preallocated : stats.Reservation));
        Assert.That(new FileInfo(dir.File("isa.bin")).Length, Is.EqualTo(4L * (n + 1)));
        var fake = new FakeReservationSource { CanPreallocate = false, Free = 1L << 40 };
        Assert.That(Dist.BuildIsa(sa, n, dir.File("isa2.bin"), Options(new MappedFileOptions { Source = fake })).Reservation, Is.EqualTo(ReservationOutcome.FreeSpaceChecked));
        Assert.That(fake.FreeCalls, Is.EqualTo(1));
        Assert.That(Dist.BuildIsa(sa, n, dir.File("isa3.bin"), Options(new MappedFileOptions { Reservation = ReservationMode.Sparse })).Reservation, Is.EqualTo(ReservationOutcome.Sparse));
        Assert.That(File.ReadAllBytes(dir.File("isa2.bin")), Is.EqualTo(File.ReadAllBytes(dir.File("isa.bin"))));
        Assert.That(File.ReadAllBytes(dir.File("isa3.bin")), Is.EqualTo(File.ReadAllBytes(dir.File("isa.bin"))));
    }

    [Test]
    public void Dist_WithNoFileOptions_OnAVolumeThatCannotHoldTheTarget_IsAnIoException_AtCreation_AndLeavesNoTarget()
    {
        using var dir = new TempDir();
        (string sa, int n) = WritePermutation(dir, 1_000_000);
        using SmallVolume volume = SmallVolume.MountOrIgnore(dir, "mnt", 2L << 20);
        string target = volume.File("isa.bin");
        var options = new DistOptions { ScratchDir = dir.Sub("scratch"), ChunkEntries = 200_000 };
        IOException? e = Assert.Throws<IOException>(() => Dist.BuildIsa(sa, n, target, options));
        Assert.That(e!.Message, Does.Contain("space").IgnoreCase.Or.Contain("full").IgnoreCase);
        Assert.That(Directory.GetFileSystemEntries(volume.Path), Is.Empty);
        Assert.That(Directory.GetFileSystemEntries(dir.Sub("scratch")), Is.Empty);
    }
}
