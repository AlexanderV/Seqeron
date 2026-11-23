using SuffixTree.Esa.External;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.External;

/// <summary>A suffix-array sink in RAM that records every write (the entries written, how often, by how many writers), for the tests of the partition and the bucket sort, which would otherwise need a file.</summary>
internal sealed class MemorySaSink(int entries) : ISaSliceSink
{
    private int _created;
    private int _disposed;
    private int _writes;

    internal uint[] Data { get; } = new uint[entries];

    internal int[] Writes { get; } = new int[entries];

    internal int Created => Volatile.Read(ref _created);

    internal int Disposed => Volatile.Read(ref _disposed);

    internal int WriteCalls => Volatile.Read(ref _writes);

    /// <summary>When set, called with the call number and the offset before every write; an exception it throws is the failure of the write.</summary>
    internal Action<int, long>? Hook { get; set; }

    public ISaSliceWriter CreateWriter()
    {
        Interlocked.Increment(ref _created);
        return new Writer(this);
    }

    private sealed class Writer(MemorySaSink owner) : ISaSliceWriter
    {
        public void Write(long entryOffset, ReadOnlySpan<uint> data)
        {
            int call = Interlocked.Increment(ref owner._writes);
            owner.Hook?.Invoke(call, entryOffset);
            for (int i = 0; i < data.Length; i++)
            {
                owner.Data[entryOffset + i] = data[i];
                Interlocked.Increment(ref owner.Writes[entryOffset + i]);
            }
        }

        public void Dispose() => Interlocked.Increment(ref owner._disposed);
    }
}

/// <summary>The comparator of PH3 (the difference-cover comparison) against the definition, on every pair of positions.</summary>
[TestFixture]
[HangLimit(30)]
public class SuffixOrderTests
{
    private TempDir _dir = null!;

    [SetUp]
    public void SetUp() => _dir = new TempDir();

    [TearDown]
    public void TearDown() => _dir.Dispose();

    private static IEnumerable<string> Texts()
    {
        yield return new string('a', 90);
        yield return StringFamilies.Fibonacci(130);
        yield return StringFamilies.ThueMorse(120);
        yield return StringFamilies.RandomDna(150, 3);
        yield return "abababababababababababababababababababab";
        yield return "banana";
        yield return "a";
        yield return LcpChildRun.RepeatRich(200, 4);
    }

    [Test]
    public void Compare_IsTheDefinition_ForEveryPairOfPositions_EveryPeriod_AndBothCommonPrefixPaths()
    {
        foreach (int x in new[] { 2, 3, 7, 9, 13, 31 })
        {
            DifferenceCover cover = DifferenceCover.Create(x);
            foreach (string text in Texts())
            {
                SaCase c = SaCase.Of(text);
                using var dir = new TempDir();
                string textPath = SaRun.WriteText(dir, c);
                using MappedTextReader reader = MappedTextReader.Open(textPath, c.N, c.Bits);
                using RankTable ranks = RankTable.Open(SaInputs.WriteRankFile(dir, c, cover), cover, c.N);
                foreach (bool vectors in new[] { true, false })
                {
                    reader.UseVectors = vectors;
                    var work = new SortWork();
                    var order = new SuffixOrder(reader, ranks, work);
                    for (uint i = 0; i < c.N; i++)
                    {
                        for (uint j = 0; j < c.N; j++)
                        {
                            int expected = DcsOracle.Compare(c.Symbols, (int)i, (int)j);
                            if (order.Compare(i, j) != expected) Assert.Fail($"X = {x}, '{text[..Math.Min(20, text.Length)]}': Compare({i}, {j}) = {order.Compare(i, j)}, expected {expected} (vectors {vectors})");
                        }
                    }
                    Assert.That(work.Compares, Is.EqualTo((long)c.N * c.N), "one call per pair");
                }
            }
        }
    }

    [Test]
    public void CompareTied_IsTheDefinition_ForPositionsWhoseKeysAreEqual()
    {
        foreach (int x in new[] { 3, 7, 13, 31 })
        {
            DifferenceCover cover = DifferenceCover.Create(x);
            foreach (string text in Texts())
            {
                SaCase c = SaCase.Of(text);
                using var dir = new TempDir();
                string textPath = SaRun.WriteText(dir, c);
                using MappedTextReader reader = MappedTextReader.Open(textPath, c.N, c.Bits);
                using RankTable ranks = RankTable.Open(SaInputs.WriteRankFile(dir, c, cover), cover, c.N);
                var order = new SuffixOrder(reader, ranks, new SortWork());
                int pairs = 0;
                for (uint i = 0; i < c.N; i++)
                {
                    for (uint j = 0; j < c.N; j++)
                    {
                        if (reader.Key(i) != reader.Key(j)) continue;
                        pairs++;
                        int expected = DcsOracle.Compare(c.Symbols, (int)i, (int)j);
                        if (order.CompareTied(i, j) != expected) Assert.Fail($"X = {x}, '{text[..Math.Min(20, text.Length)]}': CompareTied({i}, {j}) = {order.CompareTied(i, j)}, expected {expected}");
                    }
                }
                Assert.That(pairs, Is.GreaterThan(0));
            }
        }
    }

    [Test]
    public void ThePositionAtTheEnd_AndOnePastIt_IsTheEmptySuffix_TheSmallest()
    {
        SaCase c = SaCase.Of("banana");
        DifferenceCover cover = DifferenceCover.Create(7);
        string textPath = SaRun.WriteText(_dir, c);
        using MappedTextReader reader = MappedTextReader.Open(textPath, c.N, c.Bits);
        using RankTable ranks = RankTable.Open(SaInputs.WriteRankFile(_dir, c, cover), cover, c.N);
        var order = new SuffixOrder(reader, ranks, new SortWork());
        for (uint i = 0; i < c.N; i++)
        {
            Assert.That(order.Compare(i, (uint)c.N), Is.EqualTo(1), $"suffix {i} against the empty suffix");
            Assert.That(order.Compare((uint)c.N, i), Is.EqualTo(-1));
        }
        Assert.That(order.Compare((uint)c.N, (uint)c.N), Is.Zero);
    }

    [Test]
    public void TheCounters_CountTheCallsTheLceAndTheRankLookups_AndAddUp()
    {
        SaCase c = SaCase.Of(new string('a', 200) + StringFamilies.RandomDna(100, 2));
        DifferenceCover cover = DifferenceCover.Create(13);
        string textPath = SaRun.WriteText(_dir, c);
        using MappedTextReader reader = MappedTextReader.Open(textPath, c.N, c.Bits);
        using RankTable ranks = RankTable.Open(SaInputs.WriteRankFile(_dir, c, cover), cover, c.N);
        var work = new SortWork();
        var order = new SuffixOrder(reader, ranks, work);
        order.Compare(3, 50);
        Assert.That(work.Compares, Is.EqualTo(1));
        Assert.That(work.TiedCompares, Is.Zero);
        Assert.That(work.RankLookups, Is.EqualTo(2), "two sample ranks decide a long tie");
        order.CompareTied(3, 50);
        Assert.That((work.Compares, work.TiedCompares), Is.EqualTo((2L, 1L)));
        order.Compare(7, 7);
        Assert.That(work.Compares, Is.EqualTo(3), "a call with equal positions is a call");
        Assert.That(work.RankLookups, Is.EqualTo(4));
        var total = new SortWork { Compares = 10, RankLookups = 3, LceCalls = 4, LceSteps = 5, LceSymbols = 6, TiedCompares = 7, TiePositions = 8, LeafPositions = 9, ScanPositions = 11 };
        total.Add(work);
        Assert.That((total.Compares, total.RankLookups, total.LceCalls, total.TiedCompares), Is.EqualTo((13L, 7L, 4L + work.LceCalls, 8L)));
        Assert.That((total.LceSteps, total.LceSymbols, total.TiePositions, total.LeafPositions, total.ScanPositions), Is.EqualTo((5 + work.LceSteps, 6 + work.LceSymbols, 8L, 9L, 11L)));
    }
}

/// <summary>PH2 and PH3 of work package 1b.7d (the partition into lexicographic buckets and the bucket sort) against the definition of the suffix array.</summary>
[TestFixture]
[HangLimit(30)]
public class BucketSorterTests
{
    private TempDir _dir = null!;

    [SetUp]
    public void SetUp() => _dir = new TempDir();

    [TearDown]
    public void TearDown() => _dir.Dispose();

    private string Scratch => _dir.Sub("scratch");

    /// <summary>Sorts the case with the oracle rank file and checks the suffix array, that every entry was written exactly once, and that nothing is left.</summary>
    private static SortStats SortAndCompare(TempDir dir, SaCase c, DifferenceCover cover, SortOptions o, string label, MemorySaSink? sink = null)
    {
        string textPath = SaRun.WriteText(dir, c, "s" + Guid.NewGuid().ToString("N")[..5]);
        using MappedTextReader reader = MappedTextReader.Open(textPath, c.N, c.Bits);
        using RankTable ranks = RankTable.Open(SaInputs.WriteRankFile(dir, c, cover, "r" + Guid.NewGuid().ToString("N")[..5]), cover, c.N);
        sink ??= new MemorySaSink(c.N + 1);
        SortStats stats = BucketSorter.Run(reader, o, sink, ranks);
        uint[] expected = c.Sa;
        for (int k = 0; k <= c.N; k++)
        {
            if (sink.Data[k] != expected[k]) Assert.Fail($"{label}: SA[{k}] = {sink.Data[k]}, expected {expected[k]} (text '{(c.Text.Length > 40 ? c.Text[..40] + "..." : c.Text)}')");
        }
        Assert.That(sink.Writes.All(w => w == 1), Is.True, $"{label}: every entry is written exactly once");
        Assert.That(sink.Created, Is.EqualTo(sink.Disposed), $"{label}: every writer is disposed");
        Assert.That(o.Counters?.ScratchLiveBytes ?? 0, Is.Zero, $"{label}: no scratch byte left in the account");
        if (Directory.Exists(o.ScratchDir)) Assert.That(Directory.GetFileSystemEntries(o.ScratchDir), Is.Empty, $"{label}: scratch files are left behind");
        if (o.Budget is { } b) Assert.That(b.LiveBytes, Is.Zero, $"{label}: arena bytes are still reserved");
        return stats;
    }

    [Test]
    public void SmallTexts_AreTheSuffixArray_WithRotatingTinyParameters()
    {
        int v = 0;
        foreach (string text in new[] { "a", "ab", "ba", "aa", "banana", "mississippi", "abracadabra", new string('a', 60), "abababababababababab", "abcabcabcabcabcabcabc", StringFamilies.Fibonacci(90), StringFamilies.ThueMorse(75), StringFamilies.RandomDna(120, 4), StringFamilies.RandomProtein(100, 5) })
        {
            SaCase c = SaCase.Of(text);
            foreach (int variant in new[] { v, v + 5, v + 9 })
            {
                using var dir = new TempDir();
                DifferenceCover cover = DifferenceCover.Create(SaTiny.Periods[variant % SaTiny.Periods.Length]);
                SortAndCompare(dir, c, cover, SaTiny.Sort(dir.Sub("scratch"), variant, 1 + (variant % 4)), $"'{(text.Length > 20 ? text[..20] : text)}' variant {variant}");
            }
            v++;
        }
    }

    [Test]
    public void TheEmptyText_HasOnlyTheEmptySuffix()
    {
        SaCase c = SaCase.Of(string.Empty);
        SortStats stats = SortAndCompare(_dir, c, DifferenceCover.Create(7), SaTiny.Sort(Scratch, 1), "empty");
        Assert.That(stats.Leaves, Is.Zero);
    }

    [Test]
    [HangLimit(100)]
    public void EveryString_OverSigma1To3_OfTheDefaultTier_IsTheSuffixArray()
    {
        int count = 0;
        foreach ((int sigma, int max) in new[] { (1, 12), (2, 9), (3, 5) })
        {
            foreach (string s in ExhaustiveStrings.All(sigma, max))
            {
                if (s.Length == 0) continue;
                SaCase c = SaCase.Of(s);
                using var dir = new TempDir();
                int v = count++;
                SortAndCompare(dir, c, DifferenceCover.Create(SaTiny.Periods[v % SaTiny.Periods.Length]), SaTiny.Sort(dir.Sub("scratch"), v, 1 + (v % 3)), $"'{s}' variant {v}");
            }
        }
        Assert.That(count, Is.EqualTo(12 + 1022 + 363));
    }

    [Test]
    public void Families_AreTheSuffixArray_AtEveryTinyVariant_AndThePartitionNeedsSeveralLevels()
    {
        int v = 0;
        int deepest = 0;
        foreach (int n in new[] { 600, 3000 })
        {
            foreach ((string name, string text) in SaRun.Families(n))
            {
                SaCase c = SaCase.Of(text);
                using var dir = new TempDir();
                SortStats stats = SortAndCompare(dir, c, DifferenceCover.Create(SaTiny.Periods[v % SaTiny.Periods.Length]), SaTiny.Sort(dir.Sub("scratch"), v, 1 + (v % 4)), $"{name} {n} variant {v}");
                deepest = Math.Max(deepest, stats.MaxDepth);
                Assert.That(stats.Leaves, Is.GreaterThan(0));
                Assert.That(stats.MaxLeaf, Is.LessThanOrEqualTo(SaTiny.Sort("x", v).LeafCapacity), $"{name}: no leaf is larger than K");
                Assert.That(stats.Buckets, Is.EqualTo(stats.Leaves + stats.Partitions));
                v++;
            }
        }
        Assert.That(deepest, Is.GreaterThanOrEqualTo(2), "a tiny K and fan-out split the buckets more than once");
    }

    [Test]
    [HangLimit(40)]
    public void TheOutput_IsIdenticalForEveryThreadCount_LeafCapacity_FanOut_AndSeed()
    {
        foreach ((string name, string text) in SaRun.Families(9000).Where(f => f.Name is "a^n" or "fibonacci" or "tandem" or "random-dna" or "genome-like"))
        {
            SaCase c = SaCase.Of(text);
            DifferenceCover cover = DifferenceCover.Create(31);
            int v = 0;
            foreach (int threads in new[] { 1, 2, 4 })
            {
                foreach ((long k, int fan) in new[] { (17L, 2), (300L, 3), (4000L, 64), (1L << 20, 64) })
                {
                    using var dir = new TempDir();
                    SortOptions o = SaTiny.Sort(dir.Sub("scratch"), v++, threads) with { LeafCapacity = k, FanOut = fan, Seed = (ulong)v * 7 };
                    SortAndCompare(dir, c, cover, o, $"{name}: {threads} threads, K {k}, F {fan}");
                }
            }
        }
    }

    [Test]
    public void ThePartition_NeverNeedsTheRetry_BecauseTheSplittersAreDistinctSuffixesOfTheSample()
    {
        // the splitters are real suffixes ordered strictly, so with at least as many samples as children both sides of a splitter hold a sample record: a split is never degenerate (the retry is a safety net)
        SaCase c = SaCase.Of(StringFamilies.RandomDna(400, 8));
        DifferenceCover cover = DifferenceCover.Create(13);
        for (int seed = 1; seed <= 60; seed++)
        {
            using var dir = new TempDir();
            SortOptions o = SaTiny.Sort(dir.Sub("scratch"), 3) with { LeafCapacity = 20, FanOut = 2 + (seed % 5), Oversample = 1, Seed = (ulong)seed };
            Assert.That(SortAndCompare(dir, c, cover, o, $"seed {seed}").Retries, Is.Zero);
        }
    }

    [Test]
    public void ARepetitiveText_WithATinyK_IsSplitIntoEqualBuckets_BecauseTheSplittersAreRealSuffixes()
    {
        SaCase c = SaCase.Of(new string('a', 2000));
        SortStats stats = SortAndCompare(_dir, c, DifferenceCover.Create(31), SaTiny.Sort(Scratch, 2) with { LeafCapacity = 50, FanOut = 4, Oversample = 8 }, "a^n");
        Assert.That(stats.MaxLeaf, Is.LessThanOrEqualTo(50));
        Assert.That(stats.Leaves, Is.GreaterThanOrEqualTo(2000 / 50));
    }

    [Test]
    public void TheSlices_AreWrittenAtTheirFinalOffsets_ByWritersOfTheWorkers()
    {
        SaCase c = SaCase.Of(StringFamilies.Fibonacci(2000));
        var sink = new MemorySaSink(c.N + 1);
        SortStats stats = SortAndCompare(_dir, c, DifferenceCover.Create(13), SaTiny.Sort(Scratch, 4, 3), "slices", sink);
        Assert.That(sink.Created, Is.InRange(1, 3), "one writer per worker");
        Assert.That(sink.WriteCalls, Is.EqualTo(stats.Leaves + 1), "one write per leaf and the write of SA[0]");
        Assert.That(sink.Data[0], Is.EqualTo((uint)c.N));
    }

    [Test]
    public void ASinkThatThrows_EndsTheBuild_WithItsOwnException_AndEveryWriterIsDisposed()
    {
        SaCase c = SaCase.Of(StringFamilies.RandomDna(1500, 4) + new string('a', 400));
        DifferenceCover cover = DifferenceCover.Create(13);
        foreach (int failAt in new[] { 1, 2, 5, 20 })
        {
            foreach (int threads in new[] { 1, 3 })
            {
                using var dir = new TempDir();
                var sink = new MemorySaSink(c.N + 1) { Hook = (call, _) => { if (call == failAt) throw new IOException("the disk is gone"); } };
                SortOptions o = SaTiny.Sort(dir.Sub("scratch"), 2, threads) with { Counters = new IoCounters() };
                string textPath = SaRun.WriteText(dir, c);
                string rank = SaInputs.WriteRankFile(dir, c, cover);
                using MappedTextReader reader = MappedTextReader.Open(textPath, c.N, c.Bits);
                using RankTable ranks = RankTable.Open(rank, cover, c.N);
                IOException e = Assert.Throws<IOException>(() => BucketSorter.Run(reader, o, sink, ranks), $"write {failAt}, {threads} threads")!;
                Assert.That(e.Message, Is.EqualTo("the disk is gone"));
                Assert.That(sink.Created, Is.EqualTo(sink.Disposed), "every writer is disposed");
                Assert.That(Directory.Exists(o.ScratchDir) ? Directory.GetFileSystemEntries(o.ScratchDir) : [], Is.Empty, "no scratch file is left");
                Assert.That(o.Counters!.ScratchLiveBytes, Is.Zero);
                ProcFds.AssertClosed("sink failure", [], [o.ScratchDir]);
            }
        }
    }

    [Test]
    public void AHookThatThrows_AtALeafOrAPartition_ReachesTheCaller_NotWrapped()
    {
        SaCase c = SaCase.Of(StringFamilies.RandomDna(1500, 4));
        DifferenceCover cover = DifferenceCover.Create(13);
        foreach (string stage in new[] { "leaf", "partition" })
        {
            using var dir = new TempDir();
            SortOptions o = SaTiny.Sort(dir.Sub("scratch"), 2, 2) with { Trace = step => { if (step == stage) throw new FormatException(stage); } };
            string textPath = SaRun.WriteText(dir, c);
            string rank = SaInputs.WriteRankFile(dir, c, cover);
            using MappedTextReader reader = MappedTextReader.Open(textPath, c.N, c.Bits);
            using RankTable ranks = RankTable.Open(rank, cover, c.N);
            var sink = new MemorySaSink(c.N + 1);
            FormatException e = Assert.Throws<FormatException>(() => BucketSorter.Run(reader, o, sink, ranks), stage)!;
            Assert.That(e.Message, Is.EqualTo(stage));
            Assert.That(sink.Created, Is.EqualTo(sink.Disposed));
            Assert.That(Directory.Exists(o.ScratchDir) ? Directory.GetFileSystemEntries(o.ScratchDir) : [], Is.Empty);
        }
    }

    [Test]
    public void APartitionSpillThatFails_IsAnIoException_AndNoPositionFileIsLeft()
    {
        SaCase c = SaCase.Of(StringFamilies.RandomDna(3000, 4));
        DifferenceCover cover = DifferenceCover.Create(13);
        using var dir = new TempDir();
        string scratch = dir.Sub("scratch");
        bool done = false;
        SortOptions o = SaTiny.Sort(scratch, 2, 2) with
        {
            SpillBytes = 16,
            Trace = step =>
            {
                if (step == "partition" && !done)
                {
                    done = true;
                    Directory.Delete(scratch, true);
                }
            },
        };
        string textPath = SaRun.WriteText(dir, c);
        string rank = SaInputs.WriteRankFile(dir, c, cover);
        using MappedTextReader reader = MappedTextReader.Open(textPath, c.N, c.Bits);
        using RankTable ranks = RankTable.Open(rank, cover, c.N);
        var sink = new MemorySaSink(c.N + 1);
        Assert.Throws(Is.InstanceOf<IOException>(), () => BucketSorter.Run(reader, o, sink, ranks));
        Assert.That(Directory.Exists(scratch), Is.False);
        Assert.That(sink.Created, Is.EqualTo(sink.Disposed));
        ProcFds.AssertClosed("spill", [], [scratch]);
    }

    [Test]
    public void AnArenaBudgetThatIsTooSmall_EndsWithTheBudgetException_AndGivesEverythingBack()
    {
        SaCase c = SaCase.Of(StringFamilies.RandomDna(800, 3));
        DifferenceCover cover = DifferenceCover.Create(13);
        foreach (long cap in new long[] { 0, 50, (4 * (8 + 8 + 4 + 4)) - 1 })
        {
            using var dir = new TempDir();
            var budget = new ArenaBudget(cap);
            SortOptions o = SaTiny.Sort(dir.Sub("scratch"), 3, 2) with { Budget = budget, LeafCapacity = 4 };
            string textPath = SaRun.WriteText(dir, c);
            string rank = SaInputs.WriteRankFile(dir, c, cover);
            using MappedTextReader reader = MappedTextReader.Open(textPath, c.N, c.Bits);
            using RankTable ranks = RankTable.Open(rank, cover, c.N);
            var sink = new MemorySaSink(c.N + 1);
            Assert.Throws<ArenaBudgetExceededException>(() => BucketSorter.Run(reader, o, sink, ranks), $"cap {cap}");
            Assert.That(budget.LiveBytes, Is.Zero);
            Assert.That(sink.Created, Is.EqualTo(sink.Disposed));
        }
    }

    [Test]
    public void TheLeafBuffers_AreFourBlocksOfK_AndTheArenaPeakIsWhatTheBudgetAllows()
    {
        SaCase c = SaCase.Of(StringFamilies.RandomDna(2000, 3));
        DifferenceCover cover = DifferenceCover.Create(13);
        var budget = new ArenaBudget(long.MaxValue);
        SortOptions o = SaTiny.Sort(Scratch, 3, 1) with { Budget = budget, LeafCapacity = 64, ArenaGuards = null };
        SortAndCompare(_dir, c, cover, o, "peak");
        Assert.That(budget.PeakBytes, Is.EqualTo(64L * (8 + 8 + 4 + 4)), "keys, keys2, positions, positions2 of K entries");
    }

    [Test]
    public void TheCounters_CountThePartitionScansAndTheLeafPositions()
    {
        SaCase c = SaCase.Of(StringFamilies.RandomDna(3000, 5));
        var counters = new IoCounters();
        SortStats stats = SortAndCompare(_dir, c, DifferenceCover.Create(13), SaTiny.Sort(Scratch, 3, 2) with { Counters = counters, SpillBytes = 16 }, "counters");
        Assert.That(stats.Work.LeafPositions, Is.EqualTo(c.N), "every position is sorted in exactly one leaf");
        Assert.That(stats.Work.ScanPositions, Is.GreaterThanOrEqualTo(c.N), "every position goes through at least one partition scan");
        Assert.That(counters.WriteBytes, Is.GreaterThan(0), "the position files");
        Assert.That(counters.ReadBytes, Is.GreaterThan(0));
        Assert.That(counters.ScratchPeakBytes, Is.GreaterThan(0));
        Assert.That(counters.ScratchLiveBytes, Is.Zero);
        Assert.That(stats.PartitionSeconds + stats.LeafSeconds, Is.GreaterThan(0));
    }

    [Test]
    public void TheOptions_AreChecked()
    {
        SaCase c = SaCase.Of("banana");
        DifferenceCover cover = DifferenceCover.Create(7);
        string textPath = SaRun.WriteText(_dir, c);
        using MappedTextReader reader = MappedTextReader.Open(textPath, c.N, c.Bits);
        using RankTable ranks = RankTable.Open(SaInputs.WriteRankFile(_dir, c, cover), cover, c.N);
        // a small arena budget on every option set: a bound that is not enforced would allocate its buffers from the native heap, 8 GB for the capacity below (the mutant that dropped the bound killed the test host, work package 1b.9h, G1); with the budget the
        // allocation is refused with an exception that is not the one expected here, so the test fails by an assertion
        SortOptions ok = SaTiny.Sort(Scratch, 1) with { Budget = new ArenaBudget(8 << 20) };
        var sink = new MemorySaSink(c.N + 1);
        foreach ((string what, SortOptions bad) in new (string, SortOptions)[]
        {
            ("leaf capacity 1", ok with { LeafCapacity = 1 }),
            ("leaf capacity too large", ok with { LeafCapacity = (int.MaxValue / 2) + 1L }),
            ("fan-out 1", ok with { FanOut = 1 }),
            ("threads 0", ok with { Threads = 0 }),
            ("oversample 0", ok with { Oversample = 0 }),
            ("spill 0", ok with { SpillBytes = 0 }),
        })
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => BucketSorter.Run(reader, bad, sink, ranks), what);
        }
        Assert.Throws<ArgumentNullException>(() => BucketSorter.Run(null!, ok, sink, ranks));
        Assert.Throws<ArgumentNullException>(() => BucketSorter.Run(reader, null!, sink, ranks));
        Assert.Throws<ArgumentNullException>(() => BucketSorter.Run(reader, ok, null!, ranks));
        Assert.Throws<ArgumentNullException>(() => BucketSorter.Run(reader, ok, sink, null!));
    }

    [Test]
    public void TheMappedSink_WritesEverySliceAtItsOffset_ThroughWindowsOfItsOwn()
    {
        SaCase c = SaCase.Of(StringFamilies.RandomDna(1200, 7));
        string saPath = _dir.File("sa");
        var counters = new IoCounters();
        using (MappedOutputFile file = MappedOutputFile.Create(saPath, 4L * (c.N + 1)))
        {
            var sink = new MappedSaSliceSink(file, 256, counters);
            using (ISaSliceWriter a = sink.CreateWriter())
            using (ISaSliceWriter b = sink.CreateWriter())
            {
                uint[] sa = c.Sa;
                a.Write(0, sa.AsSpan(0, 300));
                b.Write(300, sa.AsSpan(300, 500));
                a.Write(800, sa.AsSpan(800));
                b.Write(0, sa.AsSpan(0, 0));
            }
            file.Complete();
        }
        Assert.That(FileBits.ReadUInts(saPath), Is.EqualTo(c.Sa));
        Assert.That(counters.MappedWriteBytes, Is.EqualTo(4L * (c.N + 1)));
        Assert.That(counters.MappedMaxWindowBytes, Is.EqualTo(256));
        Assert.That(ProcMaps.Contains(saPath) ?? false, Is.False);
    }

    [Test]
    public void TheMappedSink_RefusesAWriteOutsideTheFile()
    {
        string saPath = _dir.File("sa");
        using MappedOutputFile file = MappedOutputFile.Create(saPath, 40);
        var sink = new MappedSaSliceSink(file, 16, null);
        using ISaSliceWriter w = sink.CreateWriter();
        Assert.Throws<ArgumentOutOfRangeException>(() => w.Write(8, new uint[5]));
        Assert.Throws<ArgumentOutOfRangeException>(() => w.Write(-1, new uint[1]));
        Assert.DoesNotThrow(() => w.Write(5, new uint[5]));
    }
}
