using System.Runtime.InteropServices;
using SuffixTree.Esa.External;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.External;

/// <summary>The scratch side of DIST: chunk files of fixed-size tuples, a write buffer per chunk, sequential readers, deletion as soon as a chunk is consumed.</summary>
[TestFixture]
[HangLimit(30)]
public class DistScratchTests
{
    [StructLayout(LayoutKind.Sequential)]
    private struct Tup(uint a, uint b, uint c)
    {
        internal uint A = a;
        internal uint B = b;
        internal uint C = c;
    }

    private static BudgetedArena Arena(ArenaBudget? b = null) => new(b ?? new ArenaBudget(long.MaxValue), ArenaOptions.Guarded(0xCD));

    [Test]
    public void AppendAndRead_RoundTrip_PerChunk_InOrder()
    {
        using var dir = new TempDir();
        var c = new IoCounters();
        using var scratch = new DistScratch<Tup>(dir.Path, "t", 3, c);
        Assert.That(scratch.Chunks, Is.EqualTo(3));
        scratch.Append(1, [new Tup(1, 2, 3), new Tup(4, 5, 6)]);
        scratch.Append(1, [new Tup(7, 8, 9)]);
        scratch.Append(0, [new Tup(10, 11, 12)]);
        scratch.Append(2, ReadOnlySpan<Tup>.Empty); // nothing, and no file
        Assert.That((scratch.Entries(0), scratch.Entries(1), scratch.Entries(2)), Is.EqualTo((1L, 3L, 0L)));
        Assert.That(File.Exists(scratch.PathOf(2)), Is.False);
        Assert.That(File.Exists(scratch.PathOf(1)), Is.True);
        Assert.That(Path.GetDirectoryName(scratch.PathOf(1)), Is.EqualTo(dir.Path));
        Assert.That(Path.GetFileName(scratch.PathOf(1)), Does.StartWith("t-"));
        Assert.That(scratch.PathOf(0), Is.Not.EqualTo(scratch.PathOf(1)));
        Assert.That(c.WriteBytes, Is.EqualTo(4 * 12));
        Assert.That(c.ScratchLiveBytes, Is.EqualTo(4 * 12));

        using (DistReader<Tup> r = scratch.OpenReader(1))
        {
            Assert.That(r.Count, Is.EqualTo(3));
            var buf = new Tup[2];
            Assert.That(r.Read(buf), Is.EqualTo(2));
            Assert.That((buf[0].A, buf[1].C), Is.EqualTo((1u, 6u)));
            Assert.That(r.Read(buf), Is.EqualTo(1));
            Assert.That(buf[0].B, Is.EqualTo(8u));
            Assert.That(r.Read(buf), Is.Zero, "the end of the chunk");
            Assert.That(r.Read(Span<Tup>.Empty), Is.Zero);
        }
        using (DistReader<Tup> empty = scratch.OpenReader(2))
        {
            Assert.That(empty.Count, Is.Zero);
            Assert.That(empty.Read(new Tup[4]), Is.Zero, "a chunk that never got a tuple has no file and reads nothing");
        }
        Assert.That(c.ReadBytes, Is.EqualTo(3 * 12));
    }

    [Test]
    public void Delete_RemovesTheFile_AndGivesTheBytesBack_Once()
    {
        using var dir = new TempDir();
        var c = new IoCounters();
        using var scratch = new DistScratch<uint>(dir.Path, "d", 2, c);
        scratch.Append(0, [1, 2, 3]);
        scratch.Append(1, [4]);
        Assert.That(c.ScratchLiveBytes, Is.EqualTo(16));
        scratch.Delete(0);
        Assert.That(File.Exists(scratch.PathOf(0)), Is.False);
        Assert.That(scratch.Entries(0), Is.Zero);
        Assert.That(c.ScratchLiveBytes, Is.EqualTo(4));
        scratch.Delete(0); // nothing left: no second give-back
        scratch.Delete(1);
        Assert.That(c.ScratchLiveBytes, Is.Zero);
        Assert.That(c.ScratchPeakBytes, Is.EqualTo(16));
    }

    [Test]
    public void Dispose_RemovesTheFilesThatAreLeft()
    {
        using var dir = new TempDir();
        var c = new IoCounters();
        var scratch = new DistScratch<uint>(dir.Path, "z", 4, c);
        scratch.Append(0, [1]);
        scratch.Append(3, [1, 2]);
        scratch.Dispose();
        Assert.That(Directory.GetFileSystemEntries(dir.Path), Is.Empty);
        Assert.That(c.ScratchLiveBytes, Is.Zero);
    }

    [Test]
    public void ATruncatedChunkFile_IsSurfaced_NotAHangOrGarbage()
    {
        using var dir = new TempDir();
        using var scratch = new DistScratch<uint>(dir.Path, "tr", 1);
        scratch.Append(0, [1, 2, 3, 4, 5, 6, 7, 8]);
        using (var f = File.OpenWrite(scratch.PathOf(0))) f.SetLength(13);
        using DistReader<uint> r = scratch.OpenReader(0);
        Assert.That(r.Count, Is.EqualTo(8), "the count is the writer's, the file is what is damaged");
        Assert.Throws<EndOfStreamException>(() => r.Read(new uint[8]));
    }

    [Test]
    public void Appends_FromManyThreads_NeverOverlap()
    {
        using var dir = new TempDir();
        using var scratch = new DistScratch<uint>(dir.Path, "mt", 3);
        const int threads = 8;
        const int batches = 200;
        Parallel.For(0, threads, new ParallelOptions { MaxDegreeOfParallelism = threads }, t =>
        {
            for (int b = 0; b < batches; b++)
            {
                uint[] items = Enumerable.Range(0, 1 + (b % 5)).Select(i => (uint)((t << 24) | (b << 8) | i)).ToArray();
                scratch.Append((t + b) % 3, items);
            }
        });
        long total = 0;
        var seen = new HashSet<uint>();
        for (long chunk = 0; chunk < 3; chunk++)
        {
            using DistReader<uint> r = scratch.OpenReader(chunk);
            total += r.Count;
            var buf = new uint[7];
            int got;
            while ((got = r.Read(buf)) > 0) foreach (uint v in buf.AsSpan(0, got)) Assert.That(seen.Add(v), Is.True, "no tuple twice, none lost to an overwrite");
        }
        long expected = 0;
        for (int b = 0; b < batches; b++) expected += 1 + (b % 5);
        Assert.That(total, Is.EqualTo(threads * expected));
        Assert.That(seen, Has.Count.EqualTo(threads * expected));
    }

    [Test]
    public void Arguments_AreValidated()
    {
        using var dir = new TempDir();
        Assert.Throws<ArgumentOutOfRangeException>(() => _ = new DistScratch<uint>(dir.Path, "x", -1));
        using var s = new DistScratch<uint>(dir.Path, "x", 0);
        Assert.That(s.Chunks, Is.Zero);
    }

    // ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // the writer

    [Test]
    public void Writer_FlushesAChunkWhenItsBufferIsFull_AndKeepsTheRemainder()
    {
        using var dir = new TempDir();
        var budget = new ArenaBudget(1 << 20);
        using BudgetedArena arena = Arena(budget);
        using var scratch = new DistScratch<Tup>(dir.Path, "w", 3);
        var w = new DistWriter<Tup>(scratch, arena, 0, 3, 48); // 48 bytes = 4 tuples of 12
        Assert.That(w.Capacity, Is.EqualTo(4));
        Assert.That(budget.LiveBytes, Is.EqualTo(3 * 4 * 12), "one block for all the buffers, charged to the budget");
        for (uint i = 0; i < 3; i++) w.Add(1, new Tup(i, 0, 0));
        Assert.That((w.Pending(1), scratch.Entries(1)), Is.EqualTo((3, 0L)));
        w.Add(1, new Tup(3, 0, 0)); // the fourth fills the buffer: flushed
        Assert.That((w.Pending(1), scratch.Entries(1)), Is.EqualTo((0, 4L)));
        w.Add(1, new Tup(4, 0, 0));
        w.Add(2, new Tup(100, 0, 0));
        Assert.That(w.PendingSpan(1).ToArray().Select(t => t.A), Is.EqualTo(new uint[] { 4 }));
        Assert.That(w.PendingSpan(2).Length, Is.EqualTo(1));
        Assert.That(w.PendingSpan(0).Length, Is.Zero);
        w.FlushAll();
        Assert.That((scratch.Entries(0), scratch.Entries(1), scratch.Entries(2)), Is.EqualTo((0L, 5L, 1L)));
        Assert.That((w.Pending(0), w.Pending(1), w.Pending(2)), Is.EqualTo((0, 0, 0)));
        using DistReader<Tup> r = scratch.OpenReader(1);
        var back = new Tup[5];
        Assert.That(r.Read(back), Is.EqualTo(5));
        Assert.That(back.Select(t => t.A), Is.EqualTo(new uint[] { 0, 1, 2, 3, 4 }), "in the order added");
        arena.VerifyCanaries();
    }

    [Test]
    public void Writer_ChunkRange_IsAWindowOfTheChunks()
    {
        using var dir = new TempDir();
        using BudgetedArena arena = Arena();
        using var scratch = new DistScratch<uint>(dir.Path, "r", 10);
        var w = new DistWriter<uint>(scratch, arena, 4, 3, 16); // chunks 4, 5, 6
        w.Add(4, 1);
        w.Add(6, 2);
        w.FlushAll();
        Assert.That((scratch.Entries(4), scratch.Entries(5), scratch.Entries(6)), Is.EqualTo((1L, 0L, 1L)));
        Assert.Throws<ArgumentOutOfRangeException>(() => w.Add(3, 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => w.Add(7, 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => w.Pending(7));
        Assert.Throws<ArgumentOutOfRangeException>(() => w.Add(4L + (1L << 32), 0), "a distance that would wrap an int");
    }

    [Test]
    public void Writer_BufferSmallerThanOneTuple_HoldsOneTuple()
    {
        using var dir = new TempDir();
        using BudgetedArena arena = Arena();
        using var scratch = new DistScratch<ulong>(dir.Path, "s", 1);
        var w = new DistWriter<ulong>(scratch, arena, 0, 1, 1);
        Assert.That(w.Capacity, Is.EqualTo(1));
        w.Add(0, 5);
        Assert.That((w.Pending(0), scratch.Entries(0)), Is.EqualTo((0, 1L)), "flushed at once");
    }

    [Test]
    public void Writer_OverTheBudget_Throws_AndTheCountsOfThePropertyHold()
    {
        using var dir = new TempDir();
        using BudgetedArena arena = Arena(new ArenaBudget(100));
        using var scratch = new DistScratch<uint>(dir.Path, "b", 1);
        Assert.Throws<ArenaBudgetExceededException>(() => _ = new DistWriter<uint>(scratch, arena, 0, 3, 40), "3 buffers of 40 bytes = 120 > 100");
        Assert.DoesNotThrow(() => _ = new DistWriter<uint>(scratch, arena, 0, 2, 40));
    }

    [Test]
    public void Writer_RandomTuples_ArriveExactlyOnce_AndInOrderPerChunk()
    {
        var rng = new SplitMix64(4242);
        for (int round = 0; round < 40; round++)
        {
            using var dir = new TempDir();
            int chunks = 1 + rng.NextInt(7);
            int bufBytes = 4 * (1 + rng.NextInt(9));
            using BudgetedArena arena = Arena();
            using var scratch = new DistScratch<uint>(dir.Path, "p", chunks);
            var w = new DistWriter<uint>(scratch, arena, 0, chunks, bufBytes);
            var expected = new List<uint>[chunks];
            for (int c = 0; c < chunks; c++) expected[c] = [];
            int count = rng.NextInt(300);
            for (int i = 0; i < count; i++)
            {
                int c = rng.NextInt(chunks);
                w.Add(c, (uint)i);
                expected[c].Add((uint)i);
            }
            for (int c = 0; c < chunks; c++)
            {
                // flushed + pending = everything, in order
                var got = new List<uint>();
                using (DistReader<uint> r = scratch.OpenReader(c))
                {
                    var buf = new uint[5];
                    int n;
                    while ((n = r.Read(buf)) > 0) got.AddRange(buf.AsSpan(0, n).ToArray());
                }
                got.AddRange(w.PendingSpan(c).ToArray());
                Assert.That(got, Is.EqualTo(expected[c]), $"round {round} chunk {c}");
                Assert.That(w.Pending(c), Is.LessThan(w.Capacity));
            }
        }
    }
}
