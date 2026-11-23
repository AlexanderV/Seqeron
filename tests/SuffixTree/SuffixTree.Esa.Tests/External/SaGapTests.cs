using SuffixTree.Esa.External;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.External;

/// <summary>
/// Tests that close the gaps which the mutation audit of work package 1b.7d showed (round 1): the position store read at every offset, the radix sort on keys that differ in one digit only, a leaf above the capacity,
/// the scratch root and the account of the callers, the boundary of the base case, a failure at every call of the hook and under every arena cap, the two budgets of the phases, the layout check of a reduced string in memory, and the
/// difference covers pinned.
/// </summary>
[TestFixture]
[HangLimit(30)]
public class SaGapTests
{
    private TempDir _dir = null!;

    [SetUp]
    public void SetUp() => _dir = new TempDir();

    [TearDown]
    public void TearDown() => _dir.Dispose();

    private string Scratch => _dir.Sub("scratch");

    // ---- the position store

    [Test]
    public void ThePositionSources_ReadAtEveryOffset_AndAtAnyIndex_InMemoryAndInAFile()
    {
        string dir = _dir.Sub("pos");
        Directory.CreateDirectory(dir);
        var counters = new IoCounters();
        var ledger = new ScratchLedger(counters);
        var files = new PositionFiles(dir, "t", ledger, counters);
        Assert.That(files.NextPath(), Does.EndWith("t-00000001.bin"), "the numbering starts at one");
        Assert.That(files.NextPath(), Does.EndWith("t-00000002.bin"));
        uint[] data = Enumerable.Range(0, 100).Select(i => (uint)((i * 7919) + 3)).ToArray();
        foreach (int spill in new[] { 1 << 20, 64 })
        {
            using var writer = new PosWriter(files, spill, 100);
            foreach (uint p in data) writer.Append(p);
            Assert.That(writer.Count, Is.EqualTo(100));
            using IPosSource source = writer.Complete();
            Assert.That(source.Count, Is.EqualTo(100));
            Assert.That(source is FilePositions, Is.EqualTo(spill == 64), $"spill {spill}: a bucket that fits stays in memory, a larger one goes to a file");
            foreach (int start in new[] { 0, 1, 17, 90 })
            {
                uint[] got = new uint[Math.Min(10, 100 - start)];
                source.Read(start, got);
                Assert.That(got, Is.EqualTo(data.Skip(start).Take(got.Length).ToArray()), $"spill {spill}, from {start}");
            }
            for (int i = 0; i < 100; i += 11) Assert.That(source.At(i), Is.EqualTo(data[i]), $"spill {spill}, at {i}");
            Assert.That(source.At(0), Is.EqualTo(data[0]));
            Assert.That(source.At(1), Is.EqualTo(data[1]));
            Assert.That(source.At(99), Is.EqualTo(data[99]));
        }
        Assert.That(ledger.Live, Is.Zero, "the bytes are given back when the sources are disposed");
        var again = new PosWriter(files, 64, 100);
        for (uint p = 0; p < 100; p++) again.Append(p);
        IPosSource twice = again.Complete();
        ledger.Add(5000);
        twice.Dispose();
        twice.Dispose();
        Assert.That(ledger.Live, Is.EqualTo(5000), "a second Dispose gives nothing back again");
        Assert.That(Directory.GetFiles(dir), Is.Empty, "and the files are deleted");
    }

    [Test]
    public void AWriterThatIsDisposedUncompleted_DeletesItsFile_AndGivesTheBytesBack()
    {
        string dir = _dir.Sub("pos2");
        Directory.CreateDirectory(dir);
        var ledger = new ScratchLedger(new IoCounters());
        var files = new PositionFiles(dir, "t", ledger, null);
        var writer = new PosWriter(files, 64, 100);
        for (uint p = 0; p < 100; p++) writer.Append(p);
        Assert.That(Directory.GetFiles(dir), Is.Not.Empty, "the bucket went to a file");
        Assert.That(ledger.Live, Is.GreaterThan(0));
        writer.Dispose();
        Assert.That(Directory.GetFiles(dir), Is.Empty);
        Assert.That(ledger.Live, Is.Zero);
    }

    // ---- the sorts of the leaves

    [Test]
    public void TheRadixSort_SortsKeysThatDifferInOneDigitOnly_WhateverTheDigit()
    {
        foreach (int digit in new[] { 0, 1, 3, 4, 6, 7 })
        {
            foreach (bool outlierBelow in new[] { true, false })
            {
                const int m = 30;
                using var arena = new BudgetedArena(new ArenaBudget(long.MaxValue), ArenaOptions.Guarded(0xCD));
                Buf<ulong> k1 = arena.Allocate<ulong>(m, zero: false, "k1");
                Buf<ulong> k2 = arena.Allocate<ulong>(m, zero: false, "k2");
                Buf<uint> p1 = arena.Allocate<uint>(m, zero: false, "p1");
                Buf<uint> p2 = arena.Allocate<uint>(m, zero: false, "p2");
                const ulong majority = 0x1111111111111111UL;
                ulong outlier = (majority & ~(0xFFUL << (8 * digit))) | ((outlierBelow ? 0x05UL : 0x29UL) << (8 * digit));
                var expected = new List<(ulong Key, uint Pos)>();
                for (int i = 0; i < m; i++)
                {
                    ulong key = i == 10 ? outlier : majority;
                    k1.AsSpan(m)[i] = key;
                    p1.AsSpan(m)[i] = (uint)i;
                    expected.Add((key, (uint)i));
                }
                expected = expected.OrderBy(x => x.Key).ToList(); // stable
                KeySort.RadixSort(k1, k2, p1, p2, m);
                for (int i = 0; i < m; i++)
                {
                    Assert.That((k1.AsSpan(m)[i], p1.AsSpan(m)[i]), Is.EqualTo(expected[i]), $"digit {digit}, outlier {(outlierBelow ? "below" : "above")}, entry {i}");
                }
            }
        }
    }

    [Test]
    public void ALeafAboveTheCapacity_IsRefused()
    {
        SaCase c = SaCase.Of(StringFamilies.RandomDna(60, 3));
        DifferenceCover cover = DifferenceCover.Create(7);
        string textPath = SaRun.WriteText(_dir, c);
        using MappedTextReader reader = MappedTextReader.Open(textPath, c.N, c.Bits);
        using RankTable ranks = RankTable.Open(SaInputs.WriteRankFile(_dir, c, cover), cover, c.N);
        var order = new SuffixOrder(reader, ranks, new SortWork());
        using var arena = new BudgetedArena(new ArenaBudget(long.MaxValue), ArenaOptions.Guarded(0xCD));
        var sorter = new LeafSorter(arena, 8);
        uint[] rented = System.Buffers.ArrayPool<uint>.Shared.Rent(9);
        for (uint i = 0; i < 9; i++) rented[i] = i;
        using var bucket = new MemoryPositions(rented, 9);
        var sink = new MemorySaSink(c.N + 1);
        using ISaSliceWriter writer = sink.CreateWriter();
        Assert.Throws<InvalidOperationException>(() => sorter.Sort(bucket, 0, order, writer));
    }

    // ---- PH2 / PH3 on its own: the scratch root and the account of a caller

    [Test]
    public void TheSort_KeepsAScratchRootThatExists_UsesTheCallersLedger_AndLeavesItsBytes()
    {
        SaCase c = SaCase.Of(SaRun.GenomeLike(3000, 4));
        DifferenceCover cover = DifferenceCover.Create(13);
        Directory.CreateDirectory(Scratch);
        var counters = new IoCounters();
        var ledger = new ScratchLedger(counters);
        ledger.Add(1000); // the bytes of an outer build
        long highest = 0;
        string textPath = SaRun.WriteText(_dir, c);
        using MappedTextReader reader = MappedTextReader.Open(textPath, c.N, c.Bits);
        using RankTable ranks = RankTable.Open(SaInputs.WriteRankFile(_dir, c, cover), cover, c.N);
        var sink = new MemorySaSink(c.N + 1);
        SortOptions o = SaTiny.Sort(Scratch, 3, 1) with { Ledger = ledger, Counters = counters, SpillBytes = 16, Trace = step => { if (step == "leaf") highest = Math.Max(highest, ledger.Live); } };
        SortStats stats = BucketSorter.Run(reader, o, sink, ranks);
        Assert.That(sink.Data, Is.EqualTo(c.Sa));
        Assert.That(Directory.Exists(Scratch), Is.True, "a directory that existed is not the sort's to remove, empty or not");
        Assert.That(Directory.GetFileSystemEntries(Scratch), Is.Empty);
        Assert.That(ledger.Live, Is.EqualTo(1000), "the sort gives back what it took and nothing of the outer build");
        Assert.That(highest, Is.GreaterThan(1000), "its scratch files were booked in the ledger of the caller");
        Assert.That(stats.PartitionSeconds, Is.GreaterThan(0));
        Assert.That(stats.MaxLeaf, Is.GreaterThan(0));

        // a sort that fails does not release a ledger that is not its own, and keeps the root that existed
        int leaves = 0;
        var sink2 = new MemorySaSink(c.N + 1);
        SortOptions failing = o with { Trace = step => { if (step == "leaf" && ++leaves == 3) throw new FormatException("hook"); } };
        Assert.Throws<FormatException>(() => BucketSorter.Run(reader, failing, sink2, ranks));
        Assert.That(ledger.Live, Is.GreaterThanOrEqualTo(1000), "after a failure the caller's ledger still holds the bytes of the outer build (what the failed sort left is the caller's to settle: the builder releases the whole ledger)");
        Assert.That(Directory.Exists(Scratch), Is.True);
        Assert.That(Directory.GetFileSystemEntries(Scratch), Is.Empty);
    }

    // ---- PH1 on its own: the base case, the scratch root, the tag, a failure at every point

    [Test]
    public void TheBaseCase_IsTakenAtTheExactLength_OfTheReducedString()
    {
        SaCase c = SaCase.Of(SaRun.GenomeLike(3000, 3));
        DifferenceCover cover = DifferenceCover.Create(13);
        long red = new SampleGeometry(c.N, cover).RedLength;
        using var dir1 = new TempDir();
        ExtDcsStats atBase = ExternalDcsTests.RunAndCompare(dir1, c, cover, SaTiny.Dcs(dir1.Sub("scratch"), 4) with { BaseCase = red }, "base case at the length", out _);
        Assert.That(atBase.Levels[0].BaseCase, Is.True, "a reduced string of exactly N_base symbols is sorted in RAM");
        Assert.That(atBase.MaxLevel, Is.Zero);
        using var dir2 = new TempDir();
        ExtDcsStats below = ExternalDcsTests.RunAndCompare(dir2, c, cover, SaTiny.Dcs(dir2.Sub("scratch"), 4) with { BaseCase = red - 1 }, "base case one below", out _);
        Assert.That(below.Levels[0].BaseCase, Is.False, "one symbol more than N_base goes to the external sort");
        Assert.That(below.MaxLevel, Is.GreaterThan(0));
    }

    private static Exception? BuildAndCatch(TempDir dir, SaCase c, DifferenceCover cover, ExtDcsOptions o, string rank)
    {
        string textPath = SaRun.WriteText(dir, c, "-t" + Guid.NewGuid().ToString("N")[..4]);
        using MappedTextReader reader = MappedTextReader.Open(textPath, c.N, c.Bits);
        try
        {
            ExternalDcs.BuildRankFile(reader, cover, rank, o);
            return null;
        }
        catch (Exception e)
        {
            return e;
        }
    }

    [Test]
    public void AnOptionThatIsOutOfRange_IsRefusedBeforeAnyStepRuns_AndNothingIsCreated()
    {
        SaCase c = SaCase.Of("banana");
        DifferenceCover cover = DifferenceCover.Create(7);
        foreach ((string what, Func<ExtDcsOptions, ExtDcsOptions> bad) in new (string, Func<ExtDcsOptions, ExtDcsOptions>)[]
        {
            ("window 0", o => o with { WindowBytes = 0 }),
            ("bucket capacity 1", o => o with { BucketCapacity = 1 }),
            ("base case 0", o => o with { BaseCase = 0 }),
            ("threads 0", o => o with { Threads = 0 }),
            ("chunk slots 0", o => o with { ChunkSlots = 0 }),
            ("fan-out 1", o => o with { FanOut = 1 }),
            ("reservoir 2", o => o with { ReservoirSize = 2 }),
            ("buffers 15", o => o with { IoBufferBytes = 15 }),
        })
        {
            using var dir = new TempDir();
            int steps = 0;
            string rank = dir.File("r");
            Exception? e = BuildAndCatch(dir, c, cover, bad(SaTiny.Dcs(dir.Sub("scratch"), 2) with { Trace = _ => steps++ }), rank);
            Assert.That(e, Is.TypeOf<ArgumentOutOfRangeException>(), what);
            Assert.That(steps, Is.Zero, $"{what}: refused before the first step");
            Assert.That(Directory.Exists(dir.Sub("scratch")), Is.False, what);
            Assert.That(File.Exists(rank), Is.False, what);
        }
    }

    [Test]
    public void AScratchRootThatExists_AndIsEmpty_StaysAfterASuccessAndAFailure_AndAForeignFileOfAnotherBuildToo()
    {
        SaCase c = SaCase.Of(SaRun.GenomeLike(2500, 6));
        DifferenceCover cover = DifferenceCover.Create(13);
        Directory.CreateDirectory(Scratch);
        Assert.That(BuildAndCatch(_dir, c, cover, SaTiny.Dcs(Scratch, 4), _dir.File("r1")), Is.Null);
        Assert.That(Directory.Exists(Scratch), Is.True, "an empty scratch root that existed stays after a success");
        Assert.That(Directory.GetFileSystemEntries(Scratch), Is.Empty);
        int calls = 0;
        Exception? failure = BuildAndCatch(_dir, c, cover, SaTiny.Dcs(Scratch, 4) with { Trace = _ => { if (++calls == 5) throw new FormatException("hook"); } }, _dir.File("r2"));
        Assert.That(failure, Is.TypeOf<FormatException>());
        Assert.That(Directory.Exists(Scratch), Is.True, "and after a failure");
        Assert.That(Directory.GetFileSystemEntries(Scratch), Is.Empty);
        string[] foreign = ["e3b-00000000.tmp", "e3b-00000000-L1-pos-00000001.bin", "e3b-00000000-ranks.bin"];
        foreach (string f in foreign) File.WriteAllText(Path.Combine(Scratch, f), "another build's file");
        calls = 0;
        failure = BuildAndCatch(_dir, c, cover, SaTiny.Dcs(Scratch, 4) with { Trace = _ => { if (++calls == 5) throw new FormatException("hook"); } }, _dir.File("r3"));
        Assert.That(failure, Is.TypeOf<FormatException>());
        Assert.That(Directory.GetFileSystemEntries(Scratch).Select(Path.GetFileName).OrderBy(x => x), Is.EqualTo(foreign.OrderBy(x => x)), "the files of other builds, even of a tag that looks alike, are not touched");
    }

    [Test]
    [HangLimit(40)]
    public void APh1ThatFailsAtAnyCallOfTheHook_OrUnderAnyArenaCap_LeavesNothing_AndGivesTheAccountBack()
    {
        SaCase c = SaCase.Of(SaRun.GenomeLike(2500, 6));
        DifferenceCover cover = DifferenceCover.Create(13);
        int total = 0;
        using (var dir = new TempDir())
        {
            Assert.That(BuildAndCatch(dir, c, cover, SaTiny.Dcs(dir.Sub("scratch"), 4) with { Trace = _ => total++ }, dir.File("r")), Is.Null);
        }
        Assert.That(total, Is.GreaterThan(20), "the build calls the hook at its steps");
        for (int k = 1; k <= total; k += Math.Max(1, total / 45))
        {
            using var dir = new TempDir();
            int calls = 0;
            var counters = new IoCounters();
            string rank = dir.File("r");
            Exception? e = BuildAndCatch(dir, c, cover, SaTiny.Dcs(dir.Sub("scratch"), 4) with { Counters = counters, Trace = _ => { if (++calls == k) throw new FormatException("hook"); } }, rank);
            Assert.That(e, Is.TypeOf<FormatException>(), $"hook call {k} of {total}");
            Assert.That(File.Exists(rank), Is.False, $"hook call {k}: the rank file");
            Assert.That(Directory.Exists(dir.Sub("scratch")), Is.False, $"hook call {k}: the scratch directory");
            Assert.That(counters.ScratchLiveBytes, Is.Zero, $"hook call {k}: the account of the scratch bytes");
        }
        foreach (long cap in new long[] { 0, 100, 300, 600, 1000, 1500, 2000, 3000, 4000, 6000, 8000, 12000, 20000, 50000 })
        {
            using var dir = new TempDir();
            var counters = new IoCounters();
            string rank = dir.File("r");
            Exception? e = BuildAndCatch(dir, c, cover, SaTiny.Dcs(dir.Sub("scratch"), 4) with { Counters = counters, Budget = new ArenaBudget(cap) }, rank);
            Assert.That(counters.ScratchLiveBytes, Is.Zero, $"cap {cap}: the account of the scratch bytes");
            Assert.That(Directory.Exists(dir.Sub("scratch")), Is.False, $"cap {cap}: the scratch directory");
            if (e is null)
            {
                using RankTable t = RankTable.Open(rank, cover, c.N);
                Assert.DoesNotThrow(t.Verify, $"cap {cap}");
            }
            else
            {
                Assert.That(e, Is.TypeOf<ArenaBudgetExceededException>(), $"cap {cap}: {e}");
                Assert.That(File.Exists(rank), Is.False, $"cap {cap}: the rank file");
            }
        }
    }

    // ---- the stage: one budget per phase

    [Test]
    public void TheArenaCap_IsOneBudgetPerPhase_PH1IsCappedBeforePh23IsReached()
    {
        SaCase c = SaCase.Of(SaRun.GenomeLike(3000, 4));
        SaRun.Outcome free = SaRun.Run(_dir, c, SaTiny.Build(Scratch, 1, 1));
        long p1 = free.Result.Ph1ArenaPeakBytes, p3 = free.Result.Ph23ArenaPeakBytes;
        Assert.That(Math.Min(p1, p3), Is.GreaterThan(0));
        foreach (long cap in new[] { Math.Min(p1, p3) - 1, Math.Max(p1, p3) - 1 })
        {
            using var dir = new TempDir();
            var steps = new List<string>();
            SaBuildOptions o = SaTiny.Build(dir.Sub("scratch"), 1, 1) with { ArenaCapBytes = cap, Trace = steps.Add };
            SaRun.AssertFailsCleanly(dir, c, o, typeof(ArenaBudgetExceededException), $"cap {cap} (PH1 needs {p1}, PH2 / PH3 {p3})");
            Assert.That(steps.Contains("ph23"), Is.EqualTo(cap >= p1), $"cap {cap}: PH2 / PH3 is reached exactly when PH1 fits");
        }
    }

    // ---- a reduced string in memory

    private static unsafe void OpenInMemory(byte[] data, long n)
    {
        fixed (byte* p = data)
        {
            using MappedTextReader r = MappedTextReader.FromMemoryReduced(p, data.Length, n);
            Assert.That(r.Length, Is.EqualTo(n));
        }
    }

    [Test]
    public void AReducedStringInMemory_HasTheLayoutOfTheFile_OrIsRefused()
    {
        byte[] good = new byte[(4 * 3) + 64];
        Assert.DoesNotThrow(() => OpenInMemory(good, 3));
        byte[] badPadding = (byte[])good.Clone();
        badPadding[^1] = 1;
        Assert.Throws<InvalidDataException>(() => OpenInMemory(badPadding, 3), "a padding byte that is not zero");
        Assert.Throws<InvalidDataException>(() => OpenInMemory(new byte[(4 * 3) + 63], 3), "a truncated buffer");
        Assert.Throws<InvalidDataException>(() => OpenInMemory(new byte[(4 * 3) + 65], 3), "a buffer that is too long");
    }

    // ---- the difference covers, pinned

    [Test]
    public void TheCoversOfEveryPeriodUpTo300_AreThePinnedOnes_AndSoAreTheNearestSingerPeriods()
    {
        var text = new System.Text.StringBuilder();
        for (int x = 2; x <= 300; x++)
        {
            DifferenceCover cover = DifferenceCover.Create(x);
            text.Append(x).Append(':').Append(cover.Kind).Append(':').Append(DifferenceCover.NearestSingerPeriod(x)).Append(':').AppendJoin(',', cover.Residues).Append('\n');
        }
        string digest = Convert.ToHexString(System.Security.Cryptography.SHA256.HashData(System.Text.Encoding.ASCII.GetBytes(text.ToString())));
        TestContext.Out.WriteLine(digest);
        Assert.That(digest, Is.EqualTo(PinnedCoverDigest), "the residues, the kinds and the Singer periods of the prototype (the differential test compares every period up to 2,300 with it)");
    }

    private const string PinnedCoverDigest = "DE3523DB16CC6E08467991386F86A53CC411AFC12EA00DBA6546E765419022FD";
}
