using SuffixTree.Esa.External;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.External;

/// <summary>The window writer of the child table: 32-bit writes at arbitrary entries of a range of a shared output file through a small set of aligned windows, least recently used released first.</summary>
[TestFixture]
public class AlignedWindowWriterTests
{
    private TempDir _dir = null!;

    [SetUp]
    public void SetUp() => _dir = new TempDir();

    [TearDown]
    public void TearDown() => _dir.Dispose();

    /// <summary>The counts of an LRU set of <paramref name="slots"/> windows over the access sequence (an independent model of the stats).</summary>
    private static WindowStats Model(long[] windowOfWrite, long windowCount, int slots, bool flush, long windowBytes, long length)
    {
        var resident = new Dictionary<long, long>(); // window -> last use
        var seen = new HashSet<long>();
        long clock = 0, opened = 0, reopened = 0, flushes = 0, flushed = 0, evictions = 0, distinct = 0;
        long SizeOf(long w) => Math.Min(windowBytes, length - (w * windowBytes));
        foreach (long w in windowOfWrite)
        {
            clock++;
            if (resident.ContainsKey(w))
            {
                resident[w] = clock;
                continue;
            }
            if (resident.Count == slots)
            {
                long victim = resident.OrderBy(kv => kv.Value).First().Key;
                if (flush)
                {
                    flushes++;
                    flushed += SizeOf(victim);
                }
                resident.Remove(victim);
                evictions++;
            }
            resident[w] = clock;
            opened++;
            if (!seen.Add(w)) reopened++;
            else distinct++;
        }
        if (flush)
        {
            foreach (long w in resident.Keys)
            {
                flushes++;
                flushed += SizeOf(w);
            }
        }
        _ = windowCount;
        return new WindowStats(opened, reopened, flushes, flushed, windowOfWrite.Length, evictions, distinct);
    }

    [Test]
    public void RandomWrites_LandAtTheirEntries_ForEveryWindowSlotsAndFlushChoice_AndTheStatsFollowTheLruModel()
    {
        var rng = new SplitMix64(3);
        foreach (int window in new[] { 4, 8, 16, 64 })
        {
            foreach (int slots in new[] { 1, 2, 3, 8 })
            {
                foreach (bool flush in new[] { true, false })
                {
                    long entries = 70 + rng.NextInt(40); // not a multiple of the window: the last window is short
                    long length = entries * 4;
                    string path = _dir.File($"w{window}-{slots}-{flush}");
                    var model = new uint[entries];
                    var windows = new List<long>();
                    var counters = new IoCounters();
                    WindowStats stats;
                    using (MappedOutputFile f = MappedOutputFile.Create(path, length))
                    {
                        var w = new AlignedWindowWriter(f, 0, length, window, slots, flush, counters);
                        Assert.That(w.Length, Is.EqualTo(length));
                        for (int i = 0; i < 400; i++)
                        {
                            long idx = rng.NextInt((int)entries);
                            if (i % 3 == 0) idx = Math.Clamp((i * 7L) % entries, 0, entries - 1);
                            uint v = (uint)rng.Next();
                            w.Write32(idx, v);
                            model[idx] = v;
                            windows.Add((idx * 4) / window);
                        }
                        w.Dispose();
                        stats = w.Stats;
                        f.Complete(fsync: false);
                    }
                    Assert.That(FileBits.ReadUInts(path), Is.EqualTo(model), $"window {window} slots {slots} flush {flush}");
                    WindowStats expected = Model([.. windows], (length + window - 1) / window, slots, flush, window, length);
                    Assert.That(stats, Is.EqualTo(expected), $"window {window} slots {slots} flush {flush}");
                    Assert.That(counters.MappedWriteBytes, Is.EqualTo(4 * 400L), "the bytes of the 32-bit writes");
                    Assert.That(counters.MappedWindows, Is.EqualTo(stats.Opened));
                    Assert.That(counters.MappedMaxWindowBytes, Is.EqualTo(window));
                }
            }
        }
    }

    [Test]
    public void SequentialWrites_OpenEachWindowOnce_AndNeverReopen()
    {
        string path = _dir.File("seq");
        long length = 4 * 100;
        using MappedOutputFile f = MappedOutputFile.Create(path, length);
        var w = new AlignedWindowWriter(f, 0, length, 16, 2, true);
        for (long i = 0; i < 100; i++) w.Write32(i, (uint)(i * 3));
        w.Dispose();
        Assert.That(w.Stats, Is.EqualTo(new WindowStats(25, 0, 25, 400, 100, 23, 25)));
        Assert.That(FileBits.ReadUInts(path), Is.EqualTo(Enumerable.Range(0, 100).Select(i => (uint)(i * 3))));
    }

    [Test]
    public void OneSlot_AlternatingWindows_ReopensEveryTime()
    {
        string path = _dir.File("alt");
        using MappedOutputFile f = MappedOutputFile.Create(path, 64);
        var w = new AlignedWindowWriter(f, 0, 64, 16, 1, true);
        for (int round = 0; round < 5; round++)
        {
            w.Write32(0, 1);
            w.Write32(8, 2);
        }
        w.Dispose();
        Assert.That(w.Stats.Opened, Is.EqualTo(10));
        Assert.That(w.Stats.Reopened, Is.EqualTo(8));
        Assert.That(w.Stats.DistinctWindows, Is.EqualTo(2));
        Assert.That(w.Stats.Evictions, Is.EqualTo(9));
    }

    [Test]
    public void ARangeInsideALargerFile_IsWrittenAtItsOffset_AndNothingElse()
    {
        string path = _dir.File("range");
        using MappedOutputFile f = MappedOutputFile.Create(path, 200);
        using (MappedWindowWriter fill = f.CreateWriter(16)) fill.Append(Enumerable.Repeat((byte)0xEE, 200).ToArray());
        var w = new AlignedWindowWriter(f, 40, 80, 16, 2, true);
        for (long i = 0; i < 20; i++) w.Write32(i, (uint)i + 1);
        Assert.Throws<ArgumentOutOfRangeException>(() => w.Write32(20, 1), "the range ends after 20 entries");
        w.Dispose();
        f.Complete(fsync: false);
        byte[] all = File.ReadAllBytes(path);
        Assert.That(all[..40], Is.All.EqualTo((byte)0xEE));
        Assert.That(all[120..], Is.All.EqualTo((byte)0xEE));
        Assert.That(System.Runtime.InteropServices.MemoryMarshal.Cast<byte, uint>(all.AsSpan(40, 80)).ToArray(), Is.EqualTo(Enumerable.Range(1, 20).Select(i => (uint)i)));
    }

    [Test]
    public void ABadIndex_IsRefused_NotWritten()
    {
        string path = _dir.File("bad");
        using MappedOutputFile f = MappedOutputFile.Create(path, 40);
        using var w = new AlignedWindowWriter(f, 0, 40, 8, 2, true);
        Assert.Throws<ArgumentOutOfRangeException>(() => w.Write32(-1, 1));
        Assert.Throws<ArgumentOutOfRangeException>(() => w.Write32(10, 1));
        Assert.Throws<ArgumentOutOfRangeException>(() => w.Write32(long.MaxValue / 4, 1));
        w.Write32(9, 7);
    }

    [Test]
    public void TheConstructor_ChecksItsArguments()
    {
        string path = _dir.File("args");
        using MappedOutputFile f = MappedOutputFile.Create(path, 64);
        Assert.Throws<ArgumentNullException>(() => new AlignedWindowWriter(null!, 0, 8, 8, 1, true));
        Assert.Throws<ArgumentOutOfRangeException>(() => new AlignedWindowWriter(f, 0, 8, 3, 1, true), "a window is a power of two");
        Assert.Throws<ArgumentOutOfRangeException>(() => new AlignedWindowWriter(f, 0, 8, 2, 1, true), "a window holds at least one entry");
        Assert.Throws<ArgumentOutOfRangeException>(() => new AlignedWindowWriter(f, 0, 8, 0, 1, true));
        Assert.Throws<ArgumentOutOfRangeException>(() => new AlignedWindowWriter(f, 0, 8, 8, 0, true));
        Assert.Throws<ArgumentOutOfRangeException>(() => new AlignedWindowWriter(f, 0, -4, 8, 1, true));
        Assert.Throws<ArgumentOutOfRangeException>(() => new AlignedWindowWriter(f, -4, 8, 8, 1, true));
        Assert.Throws<ArgumentOutOfRangeException>(() => new AlignedWindowWriter(f, 60, 8, 8, 1, true), "the range ends after the file");
        Assert.Throws<ArgumentOutOfRangeException>(() => new AlignedWindowWriter(f, 0, 68, 8, 1, true));
        new AlignedWindowWriter(f, 56, 8, 8, 1, true).Dispose(); // exactly the end of the file is fine
        new AlignedWindowWriter(f, 0, 0, 8, 1, true).Dispose(); // an empty range
    }

    [Test]
    public void TheWriter_CountsItselfAsAnOpenWriter_UntilItIsDisposed()
    {
        string path = _dir.File("count");
        using MappedOutputFile f = MappedOutputFile.Create(path, 16);
        var w = new AlignedWindowWriter(f, 0, 16, 8, 1, true);
        w.Write32(0, 1);
        Assert.Throws<InvalidOperationException>(() => f.Complete(fsync: false), "a writer is still open");
        w.Dispose();
        w.Dispose(); // idempotent
        f.Complete(fsync: false);
        Assert.Throws<ObjectDisposedException>(() => w.Write32(0, 2), "a write after the disposal");
        Assert.That(ProcMaps.Contains(path) ?? false, Is.False, "the windows are unmapped");
    }

    [Test]
    public void DisposeBeforeAnyWrite_IsFine_AndOpensNothing()
    {
        string path = _dir.File("none");
        using MappedOutputFile f = MappedOutputFile.Create(path, 16);
        var w = new AlignedWindowWriter(f, 0, 16, 8, 2, true);
        w.Dispose();
        Assert.That(w.Stats, Is.EqualTo(default(WindowStats)));
        f.Complete(fsync: false);
    }

    [Test]
    public void WithoutFlushOnRelease_NoFlushIsCounted_ButTheBytesAreStillWritten()
    {
        string path = _dir.File("noflush");
        using MappedOutputFile f = MappedOutputFile.Create(path, 64);
        var w = new AlignedWindowWriter(f, 0, 64, 8, 1, false);
        for (long i = 0; i < 16; i++) w.Write32(i, (uint)i + 100);
        w.Dispose();
        Assert.That(w.Stats.Flushes, Is.Zero);
        Assert.That(w.Stats.FlushedBytes, Is.Zero);
        Assert.That(FileBits.ReadUInts(path), Is.EqualTo(Enumerable.Range(100, 16).Select(i => (uint)i)));
    }
}

/// <summary>The forward scan of PH5 (layout A, both codecs) fed from the dense LCP in blocks of any size, against the in-RAM child tables and statistics.</summary>
[TestFixture]
public class ChildScanTests
{
    private TempDir _dir = null!;

    [SetUp]
    public void SetUp() => _dir = new TempDir();

    [TearDown]
    public void TearDown() => _dir.Dispose();

    private ChildBuildOptions Options(ChildCodec codec, int block, int window = 16, int slots = 2, IoCounters? counters = null) => new()
    {
        ScratchDir = _dir.Path,
        BlockEntries = block,
        WindowBytes = window,
        WindowSlots = slots,
        SpillFlushBytes = 64,
        Codec = codec,
        ArenaGuards = ArenaOptions.Guarded(0xCD),
        Counters = counters,
    };

    private (ChildScanResult Result, uint[] Child) Run(LcpChildCase c, ChildCodec codec, int block, int feed, int window = 16, int slots = 2)
    {
        string path = _dir.File($"child-{codec}-{block}-{feed}-{window}-{slots}");
        using MappedOutputFile f = MappedOutputFile.Create(path, 4L * (c.N + 1));
        ChildScanResult result;
        using (ChildScan scan = ChildScan.Create(c.N, f, 0, Options(codec, block, window, slots)))
        {
            uint[] dense = c.DenseLcp();
            for (int at = 0; at < dense.Length; at += feed) scan.Consume(at, dense.AsSpan(at, Math.Min(feed, dense.Length - at)));
            result = scan.Complete();
        }
        f.Complete(fsync: false);
        return (result, FileBits.ReadUInts(path));
    }

    private static void Check(LcpChildCase c, ChildCodec codec, ChildScanResult r, uint[] child, string label)
    {
        Assert.That(child, Is.EqualTo(c.Child(codec)), $"{label}: CHILD ({codec})");
        Assert.That((r.NodeCount, r.InternalCount, (int)r.LrsDepth, (int)r.LrsRank, r.MaxLcp, r.N), Is.EqualTo((c.Statistics.NodeCount, (long)c.Statistics.InternalCount, c.Statistics.LrsDepth, c.Statistics.LrsRank, (long)c.MaxLcp, (long)c.N)), $"{label}: statistics");
    }

    [Test]
    public void SmallTexts_BothCodecs_EveryFeedSize_EqualTheInRamTables()
    {
        foreach (string text in new[] { "a", "aa", "ab", "banana", "mississippi", "abracadabra", "aaaaaaaaaaaa", "abababababab", "abcabcabcabc", "aabaabaab", "ba", "cabbage", StringFamilies.Fibonacci(60), StringFamilies.ThueMorse(50), StringFamilies.RandomDna(80, 5) })
        {
            LcpChildCase c = LcpChildCase.Of(text);
            foreach (ChildCodec codec in new[] { ChildCodec.Tagged, ChildCodec.Packed })
            {
                int v = 0;
                foreach (int feed in new[] { 1, 2, 3, 5, 1000 })
                {
                    (ChildScanResult r, uint[] child) = Run(c, codec, 2 + (v % 4), feed, 16 << (v % 3), 1 + (v % 3));
                    Check(c, codec, r, child, $"'{text}' feed {feed}");
                    v++;
                }
            }
        }
    }

    [Test]
    public void TheEmptyText_HasOneNode_AndAOneEntryTable()
    {
        LcpChildCase c = LcpChildCase.Of(string.Empty);
        foreach (ChildCodec codec in new[] { ChildCodec.Tagged, ChildCodec.Packed })
        {
            (ChildScanResult r, uint[] child) = Run(c, codec, 2, 1);
            Assert.That(child, Is.EqualTo(new uint[] { 0 }));
            Assert.That((r.NodeCount, r.InternalCount, r.LrsDepth, r.LrsRank, r.MaxLcp, r.StackPeak), Is.EqualTo((1L, 0L, 0u, 0u, 0L, 0L)));
        }
    }

    [Test]
    public void Families_EqualTheInRamTables_AtSeveralSizes()
    {
        foreach (int n in new[] { 150, 1200 })
        {
            foreach ((string name, string text) in LcpChildRun.Families(n))
            {
                LcpChildCase c = LcpChildCase.Of(text);
                (ChildScanResult r1, uint[] t) = Run(c, ChildCodec.Tagged, 3, 7);
                Check(c, ChildCodec.Tagged, r1, t, $"{name} {n}");
                (ChildScanResult r2, uint[] p) = Run(c, ChildCodec.Packed, 5, 1000, 64, 3);
                Check(c, ChildCodec.Packed, r2, p, $"{name} {n}");
            }
        }
    }

    [Test]
    public void ACubeOfN_ChainsTheWholeStack_AndSpillsExactlyWhenItOutgrowsTwoBlocks()
    {
        foreach (int block in new[] { 2, 3, 5 })
        {
            foreach (int n in new[] { block - 1, block, block + 1, (2 * block) - 1, 2 * block, (2 * block) + 1, (3 * block) + 1, 10 * block })
            {
                if (n < 1) continue;
                LcpChildCase c = LcpChildCase.Of(new string('a', n));
                (ChildScanResult r, uint[] child) = Run(c, ChildCodec.Packed, block, 4);
                Check(c, ChildCodec.Packed, r, child, $"a^{n} block {block}");
                Assert.That(r.StackPeak, Is.EqualTo(n), $"a^{n}: the root and the chain of intervals of depth 1 .. n - 1");
                Assert.That(r.Spill.BlocksWritten > 0, Is.EqualTo(n > 2 * block), $"a^{n} block {block}: spill");
                Assert.That(r.Spill.BytesWritten, Is.EqualTo(r.Spill.BlocksWritten * block * 16L));
                Assert.That(r.Spill.PeakFileBytes % (block * 16L), Is.Zero);
            }
        }
    }

    [Test]
    public void TheHook_AcceptsAnEmptyBlock_AnywhereInTheStream()
    {
        LcpChildCase c = LcpChildCase.Of("banana");
        string path = _dir.File("hook");
        using MappedOutputFile f = MappedOutputFile.Create(path, 4L * (c.N + 1));
        using ChildScan scan = ChildScan.Create(c.N, f, 0, Options(ChildCodec.Packed, 2));
        uint[] dense = c.DenseLcp();
        scan.Consume(0, []);
        scan.Consume(0, dense.AsSpan(0, 3));
        scan.Consume(3, []);
        scan.Consume(3, dense.AsSpan(3));
        ChildScanResult r = scan.Complete();
        Assert.That(FileBits.ReadUInts(path), Is.EqualTo(c.Packed));
        Assert.That(r.NodeCount, Is.EqualTo(c.Statistics.NodeCount));
    }

    [Test]
    public void Blocks_MustArriveInRankOrder_AndNotBeyondNPlusOne()
    {
        LcpChildCase c = LcpChildCase.Of("banana");
        using MappedOutputFile f = MappedOutputFile.Create(_dir.File("order"), 4L * (c.N + 1));
        using ChildScan scan = ChildScan.Create(c.N, f, 0, Options(ChildCodec.Packed, 2));
        uint[] d = c.DenseLcp();
        Assert.Throws<InvalidOperationException>(() => scan.Consume(1, d.AsSpan(1, 2)), "the first block must start at rank 0");
        scan.Consume(0, d.AsSpan(0, 3));
        Assert.Throws<InvalidOperationException>(() => scan.Consume(4, d.AsSpan(4, 1)), "a gap");
        Assert.Throws<InvalidOperationException>(() => scan.Consume(2, d.AsSpan(2, 2)), "a repeat");
        scan.Consume(3, d.AsSpan(3, 3));
        Assert.Throws<InvalidOperationException>(() => scan.Consume(6, [1]), "rank n + 1 does not exist");
        scan.Consume(6, []);
        Assert.That(scan.Complete().NodeCount, Is.EqualTo(c.Statistics.NodeCount));
    }

    [Test]
    public void ABlockPastNPlusOne_IsRefusedWhole()
    {
        using MappedOutputFile f = MappedOutputFile.Create(_dir.File("past"), 4L * 4);
        using ChildScan scan = ChildScan.Create(3, f, 0, Options(ChildCodec.Packed, 2));
        Assert.Throws<InvalidOperationException>(() => scan.Consume(0, [0, 0, 1, 1, 0]), "five values for four ranks");
        scan.Consume(0, [0, 0, 1, 1]);
        scan.Complete();
    }

    [Test]
    public void TheValuesAtRanks0And1_MustBeZero_AndNoValueMayExceedN()
    {
        using (MappedOutputFile f = MappedOutputFile.Create(_dir.File("v0"), 16))
        {
            using ChildScan scan = ChildScan.Create(3, f, 0, Options(ChildCodec.Packed, 2));
            Assert.Throws<InvalidDataException>(() => scan.Consume(0, [1, 0, 0, 0]));
        }
        using (MappedOutputFile f = MappedOutputFile.Create(_dir.File("v1"), 16))
        {
            using ChildScan scan = ChildScan.Create(3, f, 0, Options(ChildCodec.Packed, 2));
            Assert.Throws<InvalidDataException>(() => scan.Consume(0, [0, 1, 0, 0]));
        }
        using (MappedOutputFile f = MappedOutputFile.Create(_dir.File("v2"), 16))
        {
            using ChildScan scan = ChildScan.Create(3, f, 0, Options(ChildCodec.Packed, 2));
            Assert.Throws<InvalidDataException>(() => scan.Consume(0, [0, 0, 4, 0]), "a value above n");
        }
        using (MappedOutputFile f = MappedOutputFile.Create(_dir.File("v3"), 16))
        {
            using ChildScan scan = ChildScan.Create(3, f, 0, Options(ChildCodec.Packed, 2));
            scan.Consume(0, [0, 0, 3, 0]); // a value equal to n is not refused by the scan (no LCP of a text of length n reaches it, but the scan does not judge the array)
            scan.Complete();
        }
    }

    [Test]
    public void Complete_NeedsEveryValue_AndHappensOnce()
    {
        using MappedOutputFile f = MappedOutputFile.Create(_dir.File("complete"), 16);
        using ChildScan scan = ChildScan.Create(3, f, 0, Options(ChildCodec.Packed, 2));
        scan.Consume(0, [0, 0, 1]);
        Assert.Throws<InvalidOperationException>(() => scan.Complete(), "one value short");
        scan.Consume(3, [0]);
        scan.Complete();
        Assert.Throws<InvalidOperationException>(() => scan.Complete(), "twice");
        Assert.Throws<InvalidOperationException>(() => scan.Consume(4, [0]), "after the end");
    }

    [Test]
    public void AfterDispose_EverythingIsRefused_AndDisposeIsIdempotent()
    {
        using MappedOutputFile f = MappedOutputFile.Create(_dir.File("disposed"), 16);
        ChildScan scan = ChildScan.Create(3, f, 0, Options(ChildCodec.Packed, 2));
        scan.Consume(0, [0, 0]);
        scan.Dispose();
        scan.Dispose();
        Assert.Throws<ObjectDisposedException>(() => scan.Consume(2, [1, 0]));
        Assert.Throws<ObjectDisposedException>(() => scan.Complete());
        f.Complete(fsync: false); // the writer of the scan is closed
    }

    [Test]
    public void Create_RefusesWhatItCannotHold()
    {
        using MappedOutputFile f = MappedOutputFile.Create(_dir.File("create"), 16);
        Assert.Throws<ArgumentNullException>(() => ChildScan.Create(3, null!, 0, Options(ChildCodec.Packed, 2)));
        Assert.Throws<ArgumentNullException>(() => ChildScan.Create(3, f, 0, null!));
        Assert.Throws<ArgumentOutOfRangeException>(() => ChildScan.Create(-1, f, 0, Options(ChildCodec.Packed, 2)));
        Assert.Throws<ArgumentOutOfRangeException>(() => ChildScan.Create(4, f, 0, Options(ChildCodec.Packed, 2)), "the file holds 4 entries, n = 4 needs 5");
        Assert.Throws<ArgumentOutOfRangeException>(() => ChildScan.Create(3, f, 4, Options(ChildCodec.Packed, 2)), "the range starts too late");
        Assert.Throws<NotSupportedException>(() => ChildScan.Create(1L << 30, f, 0, Options(ChildCodec.Tagged, 2)), "ranks need more than 30 bits");
        Assert.Throws<NotSupportedException>(() => ChildScan.Create(PackedTextFormat.MaxSymbols + 1, f, 0, Options(ChildCodec.Packed, 2)), "positions are 32-bit");
        Assert.Throws<ArgumentOutOfRangeException>(() => ChildScan.Create(3, f, 0, Options(ChildCodec.Packed, 1)), "a block holds at least 2 entries");
        ChildScan.Create(3, f, 0, Options(ChildCodec.Tagged, 2)).Dispose();
    }

    [Test]
    public void ARefusedCreate_LeavesNoWriterOpenAndNothingAllocated()
    {
        using MappedOutputFile f = MappedOutputFile.Create(_dir.File("refused"), 16);
        var budget = new ArenaBudget(long.MaxValue);
        Assert.Throws<NotSupportedException>(() => ChildScan.Create(1L << 30, f, 0, Options(ChildCodec.Tagged, 2) with { Budget = budget }));
        Assert.Throws<ArgumentOutOfRangeException>(() => ChildScan.Create(3, f, 0, Options(ChildCodec.Packed, 1) with { Budget = budget }));
        Assert.That(budget.LiveBytes, Is.Zero);
        f.Complete(fsync: false);
    }

    [Test]
    public void TheArenaBudget_IsChargedForTheTwoStackBlocks_AndGivenBackOnDispose()
    {
        using MappedOutputFile f = MappedOutputFile.Create(_dir.File("budget"), 16);
        var budget = new ArenaBudget(long.MaxValue);
        ChildScan scan = ChildScan.Create(3, f, 0, Options(ChildCodec.Packed, 10) with { Budget = budget });
        Assert.That(budget.LiveBytes, Is.EqualTo(2 * 10 * 16));
        scan.Dispose();
        Assert.That(budget.LiveBytes, Is.Zero);
        Assert.Throws<ArenaBudgetExceededException>(() => ChildScan.Create(3, f, 0, Options(ChildCodec.Packed, 10) with { Budget = new ArenaBudget((2 * 10 * 16) - 1) }));
    }

    [Test]
    public void ASpillThatFails_EndsTheScan_WithTheIoException_AndDisposeStillCleansUp()
    {
        string scratch = _dir.Sub("scratch");
        Directory.CreateDirectory(scratch);
        LcpChildCase c = LcpChildCase.Of(new string('a', 40));
        string path = _dir.File("spillfail");
        using MappedOutputFile f = MappedOutputFile.Create(path, 4L * (c.N + 1));
        var counters = new IoCounters();
        ChildScan scan = ChildScan.Create(c.N, f, 0, Options(ChildCodec.Packed, 2) with { ScratchDir = scratch, Counters = counters });
        uint[] dense = c.DenseLcp();
        Directory.Delete(scratch); // the spill file cannot be created any more
        Assert.Throws<DirectoryNotFoundException>(() => scan.Consume(0, dense));
        scan.Dispose();
        f.Complete(fsync: false);
        Assert.That(counters.ScratchLiveBytes, Is.Zero);
    }

    [Test]
    public void TheCounters_SeeTheSpill_AndTheWindowWrites()
    {
        LcpChildCase c = LcpChildCase.Of(new string('a', 200));
        var counters = new IoCounters();
        string path = _dir.File("counted");
        using MappedOutputFile f = MappedOutputFile.Create(path, 4L * (c.N + 1));
        ChildScanResult r;
        using (ChildScan scan = ChildScan.Create(c.N, f, 0, Options(ChildCodec.Packed, 4, counters: counters)))
        {
            scan.Consume(0, c.DenseLcp());
            r = scan.Complete();
        }
        Assert.That(counters.WriteBytes, Is.EqualTo(r.Spill.BytesWritten));
        Assert.That(counters.ReadBytes, Is.EqualTo(r.Spill.BytesRead));
        Assert.That(counters.ScratchPeakBytes, Is.EqualTo(r.Spill.PeakFileBytes));
        Assert.That(counters.ScratchLiveBytes, Is.Zero, "the spill file is removed with the scan");
        Assert.That(counters.MappedWriteBytes, Is.EqualTo(4 * r.ChildWindows.Writes));
        Assert.That(r.Spill.BytesWritten, Is.GreaterThan(0));
    }

    [Test]
    public void LrsDepth_AndRank_AreTheFirstClosingDeepestInterval()
    {
        foreach (string text in new[] { "banana", "mississippi", "abcabxabcd", "aaaa", "abcd", "xyzxyzxyz", "a" })
        {
            LcpChildCase c = LcpChildCase.Of(text);
            (ChildScanResult r, _) = Run(c, ChildCodec.Packed, 2, 3);
            Assert.That(((int)r.LrsDepth, (int)r.LrsRank), Is.EqualTo((c.Statistics.LrsDepth, c.Statistics.LrsRank)), text);
        }
    }

    [Test]
    public void TheTwoCodecs_DifferExactlyWhereTheyShould_AndTheTaggedKindsAreInBits30And31()
    {
        LcpChildCase c = LcpChildCase.Of("mississippi");
        (_, uint[] tagged) = Run(c, ChildCodec.Tagged, 3, 4);
        (_, uint[] packed) = Run(c, ChildCodec.Packed, 3, 4);
        Assert.That(tagged.Any(v => v >> 30 != 0), Is.True);
        Assert.That(packed.All(v => v < c.N + 1 || v == 0), Is.True, "the packed codec stores raw ranks");
        for (int k = 0; k < tagged.Length; k++)
        {
            if (tagged[k] >> 30 is 2 or 3) Assert.That(packed[k], Is.EqualTo(tagged[k] & ((1u << 30) - 1)), $"up / down fields hold the same rank at {k}");
        }
    }
}
