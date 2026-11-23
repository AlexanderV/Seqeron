using System.Buffers;
using SuffixTree.Esa.External;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.External;

/// <summary>
/// The gaps that the audit of work package 1b.9c showed in the position store of the sort (<c>External/Sort/PositionStore.cs</c>) and in the pooled buffers of the partition: every pooled array that a stage rents goes back to the
/// pool, a bucket in a file closes and truncates its file. The pool of <c>ArrayPool&lt;T&gt;.Shared</c> keeps the last array that a thread returned for each size in a slot of that thread, so a test that rents, returns and
/// rents again on its own thread sees the very same array exactly when the code between the two calls returned its own. Each test names the mutants it kills (ids of <c>analysis/mutation/1b7d_*_listed.json</c>).
/// </summary>
[TestFixture]
[HangLimit(30)]
public class SaPositionGap1b9cTests
{
    private TempDir _dir = null!;

    [SetUp]
    public void SetUp() => _dir = new TempDir();

    [TearDown]
    public void TearDown() => _dir.Dispose();

    private PositionFiles NewFiles(string name, out string dir, out ScratchLedger ledger, out IoCounters counters)
    {
        dir = _dir.Sub(name);
        Directory.CreateDirectory(dir);
        counters = new IoCounters();
        ledger = new ScratchLedger(counters);
        return new PositionFiles(dir, "t", ledger, counters);
    }

    /// <summary>An array of the smallest size of the pool that this thread has just returned, so that the next Rent of that size on this thread gives it back.</summary>
    private static uint[] Probe(int size)
    {
        uint[] probe = ArrayPool<uint>.Shared.Rent(size);
        ArrayPool<uint>.Shared.Return(probe);
        return probe;
    }

    // ---------------------------------------------------------------- T04: the array of a memory source goes back to the pool, once

    [Test]
    public void TheMemorySource_GivesItsArrayBackToThePool_OnceWhateverTheNumberOfDisposals()
    {
        // T04 (MemoryPositions: the array is not returned)
        uint[] array = ArrayPool<uint>.Shared.Rent(100);
        var source = new MemoryPositions(array, 5);
        source.Dispose();
        uint[] again = ArrayPool<uint>.Shared.Rent(100);
        Assert.That(again, Is.SameAs(array), "the array that the source held is back in the pool when it is disposed");
        source.Dispose();
        uint[] third = ArrayPool<uint>.Shared.Rent(100);
        Assert.That(third, Is.Not.SameAs(array), "a second Dispose returns nothing: the array is in the hands of the caller of the first Rent now");
        ArrayPool<uint>.Shared.Return(again);
        ArrayPool<uint>.Shared.Return(third);
    }

    // ---------------------------------------------------------------- T18, T22, T28: a bucket in a file

    [Test]
    public void ThePosWriter_GivesItsBufferBackToThePool_WhenItHandsOverABucketInAFile_AndWhenItIsDisposedUncompleted()
    {
        // T18 (PosWriter: the file bucket buffer is not returned): the buffer of a bucket that went to a file is not needed after the last flush
        PositionFiles files = NewFiles("pool", out _, out _, out _);
        uint[] probe = Probe(16);
        var writer = new PosWriter(files, 64, 16); // takes the probe: capacity 16 positions
        for (uint p = 0; p < 40; p++) writer.Append(p);
        using (IPosSource source = writer.Complete())
        {
            Assert.That(source, Is.TypeOf<FilePositions>(), "40 positions do not fit 16");
            uint[] after = ArrayPool<uint>.Shared.Rent(16);
            Assert.That(after, Is.SameAs(probe), "the buffer of the writer was returned by Complete");
            ArrayPool<uint>.Shared.Return(after);
        }

        var abandoned = new PosWriter(files, 64, 16);
        for (uint p = 0; p < 40; p++) abandoned.Append(p);
        abandoned.Dispose();
        uint[] last = ArrayPool<uint>.Shared.Rent(16);
        Assert.That(last, Is.SameAs(probe), "the buffer of a writer that is disposed uncompleted is returned as well");
        ArrayPool<uint>.Shared.Return(last);
    }

    [Test]
    public void AWriterThatIsDisposedUncompleted_ClosesItsFile()
    {
        // T22 (PosWriter: Dispose does not close the file): on Linux the deletion of an open file succeeds, so the directory is empty either way; the descriptor is what shows the difference
        PositionFiles files = NewFiles("fd", out string dir, out ScratchLedger ledger, out _);
        var writer = new PosWriter(files, 64, 16);
        for (uint p = 0; p < 40; p++) writer.Append(p);
        Assert.That(Directory.GetFiles(dir), Has.Length.EqualTo(1), "the bucket went to a file");
        ProcFdsOpenOn(dir, expected: true);
        writer.Dispose();
        Assert.That(Directory.GetFiles(dir), Is.Empty);
        Assert.That(ledger.Live, Is.Zero);
        ProcFds.AssertClosed("after the disposal of the writer", [], [dir]);
    }

    private static void ProcFdsOpenOn(string dir, bool expected)
    {
        List<string>? open = ProcFds.Open();
        if (open is null) return;
        string prefix = Path.GetFullPath(dir) + "/";
        Assert.That(open.Any(f => f.StartsWith(prefix, StringComparison.Ordinal)), Is.EqualTo(expected), "the file of the bucket is open while the writer lives");
    }

    [Test]
    public void ThePosWriter_CreatesItsFileAnew_AndKeepsNothingThatWasInItsPlace()
    {
        // T28 (PosWriter: the file is opened for append): a file that the name of the bucket already points to is truncated, not written into; the bytes that were in it are not part of the bucket
        PositionFiles files = NewFiles("stale", out string dir, out _, out _);
        string path = Path.Combine(dir, "t-00000001.bin"); // the first name that the files hand out
        File.WriteAllBytes(path, [.. Enumerable.Repeat((byte)0xEE, 4000)]);
        var writer = new PosWriter(files, 64, 16);
        for (uint p = 0; p < 20; p++) writer.Append(p + 1000);
        using IPosSource source = writer.Complete();
        Assert.That(source, Is.TypeOf<FilePositions>());
        Assert.That(new FileInfo(path).Length, Is.EqualTo(80), "20 positions of 4 bytes: nothing of the 4,000 bytes that were there");
        uint[] got = new uint[20];
        source.Read(0, got);
        Assert.That(got, Is.EqualTo(Enumerable.Range(1000, 20).Select(i => (uint)i).ToArray()));
    }

    // ---------------------------------------------------------------- T14: the buffer grows up to the capacity before the bucket goes to a file

    [Test]
    public void ABucketThatFitsTheCapacity_StaysInMemory_WhateverTheSizeOfItsFirstBuffer()
    {
        // T14 (PosWriter: the buffer is flushed instead of grown): the hint of 16 gives a first buffer of 16 positions; a bucket of 1,000 positions in a capacity of 262,144 grows it, it does not spill
        PositionFiles files = NewFiles("grow", out string dir, out ScratchLedger ledger, out IoCounters counters);
        using var writer = new PosWriter(files, 1 << 20, 16);
        for (uint p = 0; p < 1000; p++) writer.Append((p * 7) + 1);
        Assert.That(Directory.GetFiles(dir), Is.Empty, "no file");
        using IPosSource source = writer.Complete();
        Assert.That(source, Is.TypeOf<MemoryPositions>());
        Assert.That(source.Count, Is.EqualTo(1000));
        for (int i = 0; i < 1000; i += 37) Assert.That(source.At(i), Is.EqualTo((uint)((i * 7) + 1)), $"position {i}: the grown buffer kept what the small one held");
        Assert.That(source.At(0), Is.EqualTo(1u));
        Assert.That(source.At(999), Is.EqualTo(6994u));
        Assert.That(ledger.Live, Is.Zero);
        Assert.That(counters.WriteBytes, Is.Zero);
    }

    [Test]
    public void ThePosWriter_SpillsAtTheSameCount_WhateverTheHintOfTheSizeAndWhateverBelowSixteenPositionsTheCapacityIs()
    {
        // not a kill of its own: the documentation of the equivalence of H71 (the writers of a partition expect the whole bucket) and T12 (the capacity has no floor of 16). The first buffer is the pool's array for
        // min(capacity, hint) and the buffer doubles up to the capacity, so the count at which a bucket spills is the smallest size of the pool at or above the capacity (32 for 20; 16 for anything up to 16) and no hint changes it
        foreach (int hint in new[] { 1, 16, 20, 1000, 1 << 20 })
        {
            foreach ((int spillBytes, int fits) in new[] { (80, 32), (8, 16), (0, 16), (64, 16) })
            {
                PositionFiles files = NewFiles($"hint-{hint}-{spillBytes}", out string dir, out _, out _);
                var writer = new PosWriter(files, spillBytes, hint);
                for (uint p = 0; p < fits; p++) writer.Append(p);
                Assert.That(Directory.GetFiles(dir), Is.Empty, $"hint {hint}, {spillBytes} bytes: {fits} positions fit");
                writer.Append(99);
                Assert.That(Directory.GetFiles(dir), Has.Length.EqualTo(1), $"hint {hint}, {spillBytes} bytes: one more goes to the file");
                writer.Dispose();
            }
        }
    }

    // ---------------------------------------------------------------- H65: the chunk of the scan of a partition goes back to the pool

    [Test]
    public void ThePartition_GivesItsReadChunkBackToThePool_WhenItEnds()
    {
        // H65 (Partition: the chunk is not returned): one worker, so the partition runs on this thread; a root of 2,000 positions in leaves of 1,500 is split once and its children are leaves. The hook "partition" returns an array of the size of the chunk
        // (16,384 positions) to this thread's slot, the partition rents it as its chunk, and the first leaf finds it in the slot again exactly when the partition returned it
        SaCase c = SaCase.Of(StringFamilies.RandomDna(2000, 3));
        DifferenceCover cover = DifferenceCover.Create(7);
        string textPath = SaRun.WriteText(_dir, c);
        using MappedTextReader text = MappedTextReader.Open(textPath, c.N, c.Bits);
        using RankTable ranks = RankTable.Open(SaInputs.WriteRankFile(_dir, c, cover), cover, c.N);
        uint[]? probe = null, atLeaf = null;
        int partitions = 0, leaves = 0;
        SortOptions options = SaTiny.Sort(_dir.Sub("scratch"), 1, 1) with
        {
            LeafCapacity = 1500,
            Trace = step =>
            {
                if (step == "partition" && ++partitions == 1) probe = Probe(1 << 14);
                else if (step == "leaf" && ++leaves == 1)
                {
                    atLeaf = ArrayPool<uint>.Shared.Rent(1 << 14);
                    ArrayPool<uint>.Shared.Return(atLeaf);
                }
            },
        };
        var sink = new MemorySaSink(c.N + 1);
        SortStats stats = BucketSorter.Run(text, options, sink, ranks);
        Assert.That(sink.Data, Is.EqualTo(c.Sa));
        Assert.That((partitions, stats.Partitions), Is.EqualTo((1, 1L)), "the root is split once and its children are leaves");
        Assert.That(leaves, Is.GreaterThan(1));
        Assert.That(probe, Is.Not.Null);
        Assert.That(atLeaf, Is.SameAs(probe), "the partition took the array of the probe as its chunk and gave it back");
    }
}
