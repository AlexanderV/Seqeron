using System.Globalization;
using SuffixTree.EmProto.Build;
using SuffixTree.EmProto.Dcs;
using SuffixTree.EmProto.Experiments;
using SuffixTree.EmProto.Harness;
using SuffixTree.EmProto.Tests.Build;
using SuffixTree.EmProto.Text;
using SuffixTree.Esa.Tests.TestSupport;
using SplitMix64 = SuffixTree.EmProto.Harness.SplitMix64;

namespace SuffixTree.EmProto.Tests.Ph1;

/// <summary>The oracle of the external-memory PH1 (gate G1 of the prototype E3b): the in-RAM <see cref="DcsBuilder"/> table, byte for byte, and the file builder to compare with it.</summary>
internal static class ExtDcsHarness
{
    /// <summary>The slots of the in-RAM rank table as bytes (little endian uint per slot, zero in the unused slots of the last block).</summary>
    internal static unsafe byte[] Oracle(PackedText text, DifferenceCover cover)
    {
        using DcsRankTable table = DcsBuilder.Build(text, cover, new WorkCounters(), 16_000_000, out _);
        return new ReadOnlySpan<byte>(table.Owner.Pointer, checked((int)(table.Slots * 4))).ToArray();
    }

    internal static ExtDcsOptions Options(string scratch, long k, long nbase, int threads = 1, ulong seed = 1) => new()
    {
        ScratchDir = scratch,
        BucketCapacity = k,
        BaseCase = nbase,
        Threads = threads,
        Seed = seed,
        FanOut = 4,
        ReservoirSize = 64,
        ChunkSlots = 40,
        RunSpillBytes = 64,
        IoBufferBytes = 128,
        MergeBufferBudget = 1 << 20,
        WindowBytes = 256,
        SortLeafCapacity = 7,
    };

    /// <summary>Builds the rank file and returns its bytes and the stats; the scratch directory must be empty afterwards.</summary>
    internal static (byte[] Bytes, ExtDcsStats Stats) Run(PackedText text, DifferenceCover cover, ExtDcsOptions options, string outDir)
    {
        string rank = Path.Combine(outDir, "rank-" + Guid.NewGuid().ToString("N") + ".bin");
        ExtDcsStats stats = ExternalDcs.BuildRankFile(text, cover, rank, options);
        byte[] bytes = File.ReadAllBytes(rank);
        File.Delete(rank);
        Assert.That(Directory.Exists(options.ScratchDir) ? Directory.GetFileSystemEntries(options.ScratchDir) : [], Is.Empty, "scratch files left behind");
        return (bytes, stats);
    }

    /// <summary>Names by definition (a dictionary over the first min(X, ...) symbols and the zero padded key), for the share of non-unique names.</summary>
    internal static (long Distinct, long NonUnique) NaiveNames(byte[] sym, int bits, DifferenceCover cover)
    {
        int n = sym.Length;
        int x = cover.Period;
        int ks = 64 / bits;
        var groups = new Dictionary<string, int>();
        long s = 0;
        for (int p = 0; p <= n; p++)
        {
            if (cover.IndexOf[p % x] < 0) continue;
            s++;
            string id;
            if (n - p < x) id = "u" + p.ToString(CultureInfo.InvariantCulture);
            else
            {
                var sb = new System.Text.StringBuilder();
                for (int i = 0; i < x; i++) sb.Append((char)(sym[p + i] + 1));
                sb.Append('|');
                for (int i = 0; i < ks; i++) sb.Append((char)(p + i < n ? sym[p + i] + 1 : 1));
                id = sb.ToString();
            }
            groups[id] = groups.GetValueOrDefault(id) + 1;
        }
        long nonUnique = groups.Values.Where(v => v >= 2).Sum();
        return (groups.Count, nonUnique);
    }
}

[TestFixture]
public class ExternalDcsExhaustiveTests
{
    private string _dir = string.Empty;

    [SetUp]
    public void SetUp()
    {
        _dir = Path.Combine(Path.GetTempPath(), "emproto-ph1-" + Guid.NewGuid().ToString("N"));
        Directory.CreateDirectory(_dir);
    }

    [TearDown]
    public void TearDown() => Directory.Delete(_dir, true);

    private static IEnumerable<byte[]> AllStrings(int sigma, int maxLength)
    {
        for (int length = 0; length <= maxLength; length++)
        {
            long total = 1;
            for (int i = 0; i < length; i++) total *= sigma;
            for (long code = 0; code < total; code++)
            {
                var s = new byte[length];
                long c = code;
                for (int i = length - 1; i >= 0; i--)
                {
                    s[i] = (byte)(c % sigma);
                    c /= sigma;
                }
                yield return s;
            }
        }
    }

    private void CheckAll(int sigma, int maxLength, int bits, int x, long k, long nbase, int threads, int minBaseCases, int minExternalLevels)
    {
        DifferenceCover cover = DifferenceCover.Create(x);
        string scratch = Path.Combine(_dir, "scratch");
        ExtDcsOptions o = ExtDcsHarness.Options(scratch, k, nbase, threads);
        long cases = 0, external = 0, baseCase = 0, splits = 0;
        foreach (byte[] s in AllStrings(sigma, maxLength))
        {
            using PackedText t = PackedText.Pack(s, bits);
            byte[] expected = ExtDcsHarness.Oracle(t, cover);
            (byte[] got, ExtDcsStats st) = ExtDcsHarness.Run(t, cover, o, _dir);
            if (!got.AsSpan().SequenceEqual(expected)) Assert.Fail($"sigma={sigma} X={x} K={k} base={nbase} b={bits} text=[{string.Join(',', s)}]: rank file differs from the in-RAM table");
            cases++;
            if (st.MaxLevel >= 1) external++;
            if (st.Levels.Any(l => l.BaseCase)) baseCase++;
            if (st.Levels.Sum(l => l.SplitScans) > 0) splits++;
        }
        TestContext.Out.WriteLine($"sigma={sigma} len<={maxLength} X={x} K={k} base={nbase}: {cases} strings, {baseCase} with a base case, {external} with an external level, {splits} with a split");
        Assert.That(cases, Is.GreaterThan(0));
        Assert.That(baseCase, Is.GreaterThanOrEqualTo(minBaseCases), "strings whose names are not unique end in the base case");
        Assert.That(external, Is.GreaterThanOrEqualTo(minExternalLevels), "the external levels must be exercised");
        Assert.That(splits, Is.GreaterThan(cases / 2), "most strings have a sample larger than K and split");
    }

    [Test]
    public void AllBinaryStrings_Length0To11_X3_TinyBuckets_BaseCaseInRam_EqualTheInRamTable() => CheckAll(2, 11, 1, 3, 3, 1 << 20, 1, 100, 0);

    [Test]
    public void AllBinaryStrings_Length0To12_X7_TinyBuckets_ExternalLevels_EqualTheInRamTable() => CheckAll(2, 12, 1, 7, 3, 3, 2, 0, 4);

    [Test]
    public void AllBinaryStrings_Length0To10_X9_PrunedCover_TwoBitSymbols_EqualTheInRamTable() => CheckAll(2, 10, 2, 9, 3, 5, 3, 0, 0);

    [Test]
    public void AllTernaryStrings_Length0To6_X5_EqualTheInRamTable() => CheckAll(3, 6, 2, 5, 2, 4, 4, 0, 0);
}

[TestFixture]
public class ExternalDcsDifferentialTests
{
    private string _dir = string.Empty;

    [SetUp]
    public void SetUp()
    {
        _dir = Path.Combine(Path.GetTempPath(), "emproto-ph1-" + Guid.NewGuid().ToString("N"));
        Directory.CreateDirectory(_dir);
    }

    [TearDown]
    public void TearDown() => Directory.Delete(_dir, true);

    [Test]
    public void AdversarialFamilies_ScaledDownParameters_EqualTheInRamTable_AndEveryPathRuns()
    {
        var rng = new SplitMix64(7070);
        string scratch = Path.Combine(_dir, "scratch");
        int cases = 0, multiLevel = 0, splits = 0, baseCases = 0, fileRuns = 0, deep = 0, multiBucket = 0;
        foreach (int n in new[] { 1000, 5000 })
        {
            foreach (DenseCase c in BuilderTestSupport.Families(n, 5))
            {
                for (int rep = 0; rep < 2; rep++)
                {
                    long k = 3 + rng.NextInt(48);
                    int x = 7 + rng.NextInt(58);
                    int bits = new[] { c.MinBits, 2, 3, 4, 8, 16, 32 }.Where(b => b >= c.MinBits).OrderBy(_ => rng.Next()).First();
                    int threads = new[] { 1, 2, 4 }[rng.NextInt(3)];
                    long nbase = rng.NextInt(3) == 0 ? 16_000_000 : 4 + rng.NextInt(197);
                    ExtDcsOptions o = ExtDcsHarness.Options(scratch, k, nbase, threads, rng.Next()) with
                    {
                        FanOut = 2 + rng.NextInt(5),
                        ReservoirSize = 8 + rng.NextInt(200),
                        ChunkSlots = 4 + rng.NextInt(100),
                        RunSpillBytes = 16 + rng.NextInt(300),
                        IoBufferBytes = 16 + rng.NextInt(400),
                        MergeBufferBudget = 8 + rng.NextInt(3000),
                        WindowBytes = 16 + rng.NextInt(4000),
                        SortLeafCapacity = 3 + rng.NextInt(48),
                    };
                    DifferenceCover cover = DifferenceCover.Create(x);
                    using PackedText t = PackedText.Pack(c.Symbols, bits);
                    byte[] expected = ExtDcsHarness.Oracle(t, cover);
                    (byte[] got, ExtDcsStats st) = ExtDcsHarness.Run(t, cover, o, _dir);
                    if (!got.AsSpan().SequenceEqual(expected))
                        Assert.Fail($"{c.Name} K={k} X={x} b={bits} threads={threads} base={nbase} F={o.FanOut} reservoir={o.ReservoirSize} R={o.ChunkSlots}: rank file differs from the in-RAM table");
                    Assert.That(st.Levels[0].SampleSize, Is.EqualTo(cover.SampleCount(n)));
                    (long distinct, long nonUnique) = ExtDcsHarness.NaiveNames(c.Symbols, bits, cover);
                    Assert.That(st.Levels[0].DistinctNames, Is.EqualTo(distinct), $"{c.Name} X={x} b={bits}: distinct names");
                    Assert.That(st.Levels[0].NonUniqueRecords, Is.EqualTo(nonUnique), $"{c.Name} X={x} b={bits}: non-unique records");
                    cases++;
                    if (st.MaxLevel >= 1) multiLevel++;
                    if (st.MaxLevel >= 2) deep++;
                    if (st.Levels.Sum(l => l.SplitScans) > 0) splits++;
                    if (st.Levels.Any(l => l.BaseCase)) baseCases++;
                    if (st.Levels.Sum(l => l.FileRuns) > 0) fileRuns++;
                    if (st.Levels[0].Leaves > 1) multiBucket++;
                }
            }
        }
        Assert.That(cases, Is.GreaterThan(100));
        Assert.That(multiBucket, Is.GreaterThan(cases / 2), "most cases must have several sample buckets");
        Assert.That(splits, Is.GreaterThan(cases / 2));
        Assert.That(multiLevel, Is.GreaterThan(20), "names are not unique on the periodic and run families: the recursion must run");
        Assert.That(deep, Is.GreaterThan(5), "levels 2 and below must run");
        Assert.That(baseCases, Is.GreaterThan(20));
        Assert.That(fileRuns, Is.GreaterThan(20), "runs above the spill size must go to scratch files");
    }

    [Test]
    public void LargerTexts_ManyLevels_EqualTheInRamTable()
    {
        var rng = new SplitMix64(11);
        string scratch = Path.Combine(_dir, "scratch");
        foreach (DenseCase c in new[]
        {
            BuilderTestSupport.Dense("dna", StringFamilies.RandomDna(60_000, 9)),
            BuilderTestSupport.Dense("same", StringFamilies.AllSame(40_000)),
            BuilderTestSupport.Dense("fib", StringFamilies.Fibonacci(30_000)),
            BuilderTestSupport.Dense("p64", StringFamilies.Periodic(30_000, 64, 3)),
            BuilderTestSupport.Dense("runs", StringFamilies.Runs(30_000, 300, 4)),
        })
        {
            int x = 13 + rng.NextInt(50);
            DifferenceCover cover = DifferenceCover.Create(x);
            using PackedText t = PackedText.Pack(c.Symbols, c.MinBits);
            byte[] expected = ExtDcsHarness.Oracle(t, cover);
            ExtDcsOptions o = ExtDcsHarness.Options(scratch, 40 + rng.NextInt(200), 300, 4) with { ChunkSlots = 700, RunSpillBytes = 512, IoBufferBytes = 1024, WindowBytes = 4096, SortLeafCapacity = 50 };
            (byte[] got, ExtDcsStats st) = ExtDcsHarness.Run(t, cover, o, _dir);
            Assert.That(got.AsSpan().SequenceEqual(expected), Is.True, $"{c.Name} X={x}");
            TestContext.Out.WriteLine($"{c.Name} X={x}: levels {st.Levels.Count}, leaves {st.Levels[0].Leaves}, unique at level 0: {st.Levels[0].NamesUnique}");
        }
    }

    [Test]
    public void Output_DoesNotDependOnThreadsOrBucketCapacity()
    {
        string scratch = Path.Combine(_dir, "scratch");
        foreach (DenseCase c in new[]
        {
            BuilderTestSupport.Dense("binary", StringFamilies.Random(3000, 2, 6)),
            BuilderTestSupport.Dense("fib", StringFamilies.Fibonacci(3000)),
            BuilderTestSupport.Dense("same", StringFamilies.AllSame(3000)),
            BuilderTestSupport.Dense("dna", StringFamilies.RandomDna(3000, 7)),
        })
        {
            DifferenceCover cover = DifferenceCover.Create(31);
            using PackedText t = PackedText.Pack(c.Symbols, c.MinBits);
            byte[] expected = ExtDcsHarness.Oracle(t, cover);
            foreach (int threads in new[] { 1, 2, 4 })
            {
                foreach (long k in new long[] { 3, 17, 200, 1 << 20 })
                {
                    foreach (long nbase in new long[] { 6, 150, 1 << 20 })
                    {
                        (byte[] got, _) = ExtDcsHarness.Run(t, cover, ExtDcsHarness.Options(scratch, k, nbase, threads, (ulong)(threads * 1000) + (ulong)k) with { TouchBuffers = k == 17 }, _dir);
                        Assert.That(got.AsSpan().SequenceEqual(expected), Is.True, $"{c.Name} threads={threads} K={k} base={nbase}");
                    }
                }
            }
        }
    }

    [TestCase("dna:300000", 100)]
    [TestCase("genome:400000", 80)]
    [TestCase("tandem:171x400", 30)]
    [TestCase("nrun:300000:40000", 60)]
    [TestCase("a:200000", 50)]
    public void ProductionPeriod993_ScaledBucketsAndBaseCase_EqualTheInRamTable(string family, int bucketsTarget)
    {
        SuffixTree.EmProto.Gen.SymbolText text = FamilyLoader.Load(family, string.Empty);
        DifferenceCover cover = DifferenceCover.Create(993);
        int bits = Math.Max(1, PackedText.BitsFor(text.Sigma));
        using PackedText t = PackedText.Pack(text.Symbols, bits);
        byte[] expected = ExtDcsHarness.Oracle(t, cover);
        long s = cover.SampleCount(text.Length);
        ExtDcsOptions o = new()
        {
            ScratchDir = Path.Combine(_dir, "scratch"),
            BucketCapacity = Math.Max(3, s / bucketsTarget * 2),
            BaseCase = Math.Max(10, s / 20),
            Threads = 4,
            ChunkSlots = Math.Max(10, s / 6),
            RunSpillBytes = 4096,
            IoBufferBytes = 4096,
            WindowBytes = 16384,
            SortLeafCapacity = 1000,
            Budget = new ArenaBudget(64L << 20),
        };
        (byte[] got, ExtDcsStats st) = ExtDcsHarness.Run(t, cover, o, _dir);
        Assert.That(got.AsSpan().SequenceEqual(expected), Is.True, family);
        Assert.That(st.Levels[0].Leaves, Is.GreaterThan(bucketsTarget / 2), "the scaled capacity must give many buckets");
        Assert.That(o.Budget!.LiveBytes, Is.EqualTo(0), "every native buffer returned to the budget");
        Assert.That(o.Budget.PeakBytes, Is.GreaterThan(0));
        TestContext.Out.WriteLine($"{family}: levels {st.Levels.Count}, leaves {st.Levels[0].Leaves}, unique {st.Levels[0].NamesUnique}, nonunique {st.Levels[0].NonUniqueRecords}/{st.Levels[0].SampleSize}");
    }

    [Test]
    public void RecursionDepth_OfAllSame_ShrinksByAboutPeriodOverSampleDensityPerLevel()
    {
        string scratch = Path.Combine(_dir, "scratch");
        DenseCase c = BuilderTestSupport.Dense("same", StringFamilies.AllSame(20_000));
        DifferenceCover cover = DifferenceCover.Create(31); // |D| = 6, density 6/31
        using PackedText t = PackedText.Pack(c.Symbols, 1);
        (_, ExtDcsStats st) = ExtDcsHarness.Run(t, cover, ExtDcsHarness.Options(scratch, 100, 50, 2), _dir);
        foreach (ExtDcsLevelStats l in st.Levels) TestContext.Out.WriteLine($"level {l.Level}: n={l.TextLength} S={l.SampleSize} distinct={l.DistinctNames} nonunique={l.NonUniqueRecords} base={l.BaseCase} leaves={l.Leaves}");
        Assert.That(st.Levels[0].NamesUnique, Is.False);
        Assert.That(st.Levels[0].TextLength, Is.EqualTo(20_000));
        Assert.That(st.MaxLevel, Is.EqualTo(2), "levels 0 and 1 have repeated names, level 2 (758 symbols, runs shorter than X) has unique names");
        Assert.That(st.Levels[^1].NamesUnique, Is.True);
        Assert.That(st.Levels.All(l => !l.BaseCase), Is.True, "N_base = 50 is never reached: the chain ends with unique names");
        for (int l = 1; l <= st.MaxLevel; l++)
        {
            Assert.That(st.Levels[l].Level, Is.EqualTo(l));
            Assert.That(st.Levels[l].TextLength, Is.EqualTo(st.Levels[l - 1].SampleSize + cover.Residues.Length), "the reduced string has one symbol per sample position and one separator per class");
        }
        (_, ExtDcsStats inRam) = ExtDcsHarness.Run(t, cover, ExtDcsHarness.Options(scratch, 100, 1000, 2), _dir);
        Assert.That(inRam.MaxLevel, Is.EqualTo(1), "with N_base = 1000 the 3,878 symbol reduced string is an external level whose own reduced string (758) is sorted in RAM");
        Assert.That(inRam.Levels[1].BaseCase, Is.True);
        Assert.That(inRam.Levels[0].BaseCase, Is.False);
    }

    [Test]
    public void RandomDna_WithProductionPeriod_HasUniqueNames_AndNoRecursion()
    {
        SuffixTree.EmProto.Gen.SymbolText text = FamilyLoader.Load("dna:200000", string.Empty);
        DifferenceCover cover = DifferenceCover.Create(993);
        using PackedText t = PackedText.Pack(text.Symbols, 2);
        (byte[] got, ExtDcsStats st) = ExtDcsHarness.Run(t, cover, ExtDcsHarness.Options(Path.Combine(_dir, "scratch"), 500, 100), _dir);
        Assert.That(st.MaxLevel, Is.EqualTo(0));
        Assert.That(st.Levels[0].NamesUnique, Is.True);
        Assert.That(st.Levels[0].NonUniqueRecords, Is.EqualTo(0));
        Assert.That(st.Levels[0].Leaves, Is.GreaterThan(5));
        Assert.That(got.AsSpan().SequenceEqual(ExtDcsHarness.Oracle(t, cover)), Is.True);
    }

    [Test]
    public void RankFile_IsWrittenThroughBoundedMappedWindows_AndScratchIsTracked()
    {
        DenseCase c = BuilderTestSupport.Dense("fib", StringFamilies.Fibonacci(30_000));
        DifferenceCover cover = DifferenceCover.Create(31);
        using PackedText t = PackedText.Pack(c.Symbols, 1);
        MappedIo.Reset();
        IoCounters.Reset();
        ScratchTracker.Reset();
        ExtDcsOptions o = ExtDcsHarness.Options(Path.Combine(_dir, "scratch"), 300, 200, 2) with { WindowBytes = 4096, RunSpillBytes = 256, ChunkSlots = 500 };
        (byte[] got, _) = ExtDcsHarness.Run(t, cover, o, _dir);
        Assert.That(MappedIo.WriteBytes, Is.GreaterThanOrEqualTo(cover.SampleCount(30_000) * 4), "every rank entry goes through a window");
        Assert.That(MappedIo.Windows, Is.GreaterThan(got.Length / 8192), "the file is written window by window");
        Assert.That(MappedIo.MaxWindowBytes, Is.LessThanOrEqualTo(4096));
        Assert.That(ScratchTracker.Peak, Is.GreaterThan(0));
        Assert.That(ScratchTracker.Live, Is.EqualTo(0), "every scratch file was deleted and untracked");
        Assert.That(IoCounters.ReadBytes, Is.GreaterThan(0));
        Assert.That(IoCounters.WriteBytes, Is.GreaterThan(0));
    }

    [Test]
    public void MergeBuffers_ShareOneBudget_InsteadOfGrowingWithTheNumberOfRuns()
    {
        // written after the first version of the merge, which gave every run file its own 64 KB buffer: the part of the working memory that grows with n / K_s is now capped
        DenseCase c = BuilderTestSupport.Dense("dna", StringFamilies.RandomDna(50_000, 12));
        DifferenceCover cover = DifferenceCover.Create(31);
        using PackedText t = PackedText.Pack(c.Symbols, 2);
        byte[] expected = ExtDcsHarness.Oracle(t, cover);
        long[] peak = new long[2];
        for (int i = 0; i < 2; i++)
        {
            ExtDcsOptions o = ExtDcsHarness.Options(Path.Combine(_dir, "scratch"), 100, 1 << 20) with { IoBufferBytes = 4096, RunSpillBytes = 64, MergeBufferBudget = i == 0 ? 1 << 20 : 19_000, Budget = new ArenaBudget(64 << 20) };
            (byte[] got, ExtDcsStats st) = ExtDcsHarness.Run(t, cover, o, _dir);
            Assert.That(got.AsSpan().SequenceEqual(expected), Is.True, $"budget case {i}");
            Assert.That(st.Levels[0].FileRuns, Is.GreaterThan(80), "many run files");
            peak[i] = o.Budget!.PeakBytes;
        }
        Assert.That(peak[1] + 300_000, Is.LessThan(peak[0]), $"the small budget must shrink the merge buffers ({peak[1]} against {peak[0]} bytes)");
    }

    [Test]
    public void TinyAndDegenerateTexts_EqualTheInRamTable()
    {
        string scratch = Path.Combine(_dir, "scratch");
        foreach (int x in new[] { 7, 13, 993 })
        {
            DifferenceCover cover = DifferenceCover.Create(x);
            foreach (byte[] s in new byte[][] { [], [0], [1], [0, 0], [1, 0], new byte[x - 1], new byte[x], new byte[x + 1], new byte[(2 * x) + 3] })
            {
                using PackedText t = PackedText.Pack(s, 1);
                byte[] expected = ExtDcsHarness.Oracle(t, cover);
                (byte[] got, _) = ExtDcsHarness.Run(t, cover, ExtDcsHarness.Options(scratch, 2, 3, 2), _dir);
                Assert.That(got, Is.EqualTo(expected), $"X={x} n={s.Length}");
            }
        }
    }

    [Test]
    public void InvalidOptions_AreRefused()
    {
        string scratch = Path.Combine(_dir, "scratch");
        DifferenceCover cover = DifferenceCover.Create(7);
        using PackedText t = PackedText.Pack(new byte[] { 0, 1, 0 }, 1);
        string rank = Path.Combine(_dir, "r.bin");
        Assert.Throws<ArgumentOutOfRangeException>(() => ExternalDcs.BuildRankFile(t, cover, rank, ExtDcsHarness.Options(scratch, 1, 4)));
        Assert.Throws<ArgumentOutOfRangeException>(() => ExternalDcs.BuildRankFile(t, cover, rank, ExtDcsHarness.Options(scratch, 4, 0)));
        Assert.Throws<ArgumentOutOfRangeException>(() => ExternalDcs.BuildRankFile(t, cover, rank, ExtDcsHarness.Options(scratch, 4, 4, 0)));
        Assert.Throws<ArgumentOutOfRangeException>(() => ExternalDcs.BuildRankFile(t, cover, rank, ExtDcsHarness.Options(scratch, 4, 4) with { ChunkSlots = 0 }));
        ExtDcsOptions ok = ExtDcsHarness.Options(scratch, 4, 4);
        foreach (ExtDcsOptions bad in new[]
        {
            ok with { FanOut = 1 },
            ok with { ReservoirSize = 2 },
            ok with { IoBufferBytes = 15 },
            ok with { WindowBytes = 0 },
            ok with { MergeBufferBudget = 7 },
            ok with { SortLeafCapacity = 1 },
            ok with { RunSpillBytes = -1 },
            ok with { BucketCapacity = (int.MaxValue / 2) + 1L },
        })
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => ExternalDcs.BuildRankFile(t, cover, rank, bad), $"{bad}");
        }
        // the smallest accepted value of every option works and gives the right table
        byte[] expected = ExtDcsHarness.Oracle(t, cover);
        foreach (ExtDcsOptions edge in new[]
        {
            ok with { FanOut = 2 },
            ok with { ReservoirSize = 3 },
            ok with { IoBufferBytes = 16 },
            ok with { WindowBytes = 1 },
            ok with { MergeBufferBudget = 8 },
            ok with { SortLeafCapacity = 2 },
            ok with { RunSpillBytes = 0 },
            ok with { BucketCapacity = 2, BaseCase = 1, ChunkSlots = 1, Threads = 1 },
            ok with { BucketCapacity = int.MaxValue / 2 },
        })
        {
            (byte[] got, _) = ExtDcsHarness.Run(t, cover, edge, _dir);
            Assert.That(got, Is.EqualTo(expected), $"{edge}");
        }
    }

    [Test]
    public void BaseCaseBoundary_ReducedLengthEqualToNBase_IsTheBaseCase_OneLessIsAnExternalLevel()
    {
        string scratch = Path.Combine(_dir, "scratch");
        DenseCase c = BuilderTestSupport.Dense("same", StringFamilies.AllSame(2000));
        DifferenceCover cover = DifferenceCover.Create(31);
        using PackedText t = PackedText.Pack(c.Symbols, 1);
        long reduced = new SampleGeometry(2000, cover).RedLength;
        (_, ExtDcsStats atLimit) = ExtDcsHarness.Run(t, cover, ExtDcsHarness.Options(scratch, 100, reduced, 2), _dir);
        Assert.That(atLimit.MaxLevel, Is.EqualTo(0));
        Assert.That(atLimit.Levels[0].BaseCase, Is.True, "a reduced string of exactly N_base symbols is sorted in RAM");
        (_, ExtDcsStats below) = ExtDcsHarness.Run(t, cover, ExtDcsHarness.Options(scratch, 100, reduced - 1, 2), _dir);
        Assert.That(below.Levels[0].BaseCase, Is.False);
        Assert.That(below.MaxLevel, Is.GreaterThanOrEqualTo(1), "one symbol more than N_base goes through an external level");
    }

    [Test]
    public void ReducedStringNotShorterThanTheText_IsSortedInRam_InsteadOfRecursingForEver()
    {
        string scratch = Path.Combine(_dir, "scratch");
        // a period and a length at which the reduced string (Samples + |D| symbols) is exactly as long as the text, with at least two sample suffixes of X symbols (a^n has repeated names then)
        (int x, int n) = (from period in Enumerable.Range(3, 30)
                          from len in Enumerable.Range(period + 2, 90)
                          let g = new SampleGeometry(len, DifferenceCover.Create(period))
                          where g.RedLength == len && Enumerable.Range(0, len - period + 1).Count(p => g.Cover.IndexOf[p % period] >= 0) >= 2
                          select (period, len)).First();
        DifferenceCover cover = DifferenceCover.Create(x);
        using PackedText t = PackedText.Pack(new byte[n], 1);
        (byte[] got, ExtDcsStats st) = ExtDcsHarness.Run(t, cover, ExtDcsHarness.Options(scratch, 3, 3, 1), _dir);
        Assert.That(st.Levels[0].NamesUnique, Is.False);
        Assert.That(st.MaxLevel, Is.EqualTo(0), $"a^{n} with X = {x} has a reduced string of {n} symbols: no gain, so the base case");
        Assert.That(st.Levels[0].BaseCase, Is.True);
        Assert.That(got, Is.EqualTo(ExtDcsHarness.Oracle(t, cover)));
        int longer = Enumerable.Range(n + 1, 200).First(len => new SampleGeometry(len, cover).RedLength < len);
        using PackedText down = PackedText.Pack(new byte[longer], 1);
        (byte[] got2, ExtDcsStats st2) = ExtDcsHarness.Run(down, cover, ExtDcsHarness.Options(scratch, 3, 3, 1), _dir);
        Assert.That(st2.MaxLevel, Is.GreaterThanOrEqualTo(1), $"a^{longer}: the reduced string is shorter than the text, so the external level is used");
        Assert.That(got2, Is.EqualTo(ExtDcsHarness.Oracle(down, cover)));
    }

    [Test]
    public void BucketBuffers_AreSizedByTheSampleCount_NotByTheCapacity()
    {
        DenseCase c = BuilderTestSupport.Dense("dna", StringFamilies.RandomDna(5000, 3));
        DifferenceCover cover = DifferenceCover.Create(31);
        using PackedText t = PackedText.Pack(c.Symbols, 2);
        var budget = new ArenaBudget(64L << 20);
        ExtDcsOptions o = ExtDcsHarness.Options(Path.Combine(_dir, "scratch"), 1 << 20, 1 << 20) with { Budget = budget };
        (byte[] got, _) = ExtDcsHarness.Run(t, cover, o, _dir);
        Assert.That(got.AsSpan().SequenceEqual(ExtDcsHarness.Oracle(t, cover)), Is.True);
        Assert.That(budget.PeakBytes, Is.LessThan(1 << 20), "a text with about 1,000 sample positions must not reserve 24 MB for a capacity of 1 Mi");
    }

    [Test]
    public void SplitScan_IsUsedForTheRoot_AndTheLeavesCountIsNearTwiceTheSampleOverTheCapacity()
    {
        DenseCase c = BuilderTestSupport.Dense("dna", StringFamilies.RandomDna(60_000, 8));
        DifferenceCover cover = DifferenceCover.Create(31);
        using PackedText t = PackedText.Pack(c.Symbols, 2);
        long s = cover.SampleCount(60_000);
        const long k = 200;
        (_, ExtDcsStats st) = ExtDcsHarness.Run(t, cover, ExtDcsHarness.Options(Path.Combine(_dir, "scratch"), k, 1 << 20) with { FanOut = 1024, ReservoirSize = 4096 }, _dir);
        ExtDcsLevelStats l = st.Levels[0];
        Assert.That(l.SplitScans, Is.GreaterThanOrEqualTo(1));
        Assert.That(l.Leaves, Is.InRange((2 * s / k) - 5, (2 * s / k) + 60), $"the root splits into about 2 x sample / K = {2 * s / k} buckets");
        Assert.That(l.Scans, Is.GreaterThanOrEqualTo(l.Leaves + l.SplitScans), "every leaf and every split scans the sample once; an overflowing collect scan adds one more");
    }
}

/// <summary>
/// Differential fuzz of the external-memory PH1 (ADR-004 section 5, "Oracles per stage"): scaled-down K = 3..50, N_base = 4..200, X = 7..64, F = 2..6 on texts of 10^3..10^5 symbols from the adversarial
/// families, against the in-RAM table. Category <c>Fuzz</c>, <c>[Explicit]</c>; environment <c>EMPROTO_FUZZ_SEED</c> (default 1) and <c>EMPROTO_FUZZ_SECONDS</c> (the budget, default 120).
/// </summary>
[TestFixture]
[Category("Fuzz")]
[Explicit("Differential fuzz of the external-memory PH1; run with --filter Category=Fuzz (see PROGRESS-S1.md, Log E3b)")]
public class ExternalDcsFuzzTests
{
    [Test]
    public void ScaledDownParameters_Texts1kTo100k_EqualTheInRamTable_WithinTheBudget()
    {
        ulong seed = ulong.Parse(Environment.GetEnvironmentVariable("EMPROTO_FUZZ_SEED") ?? "1", CultureInfo.InvariantCulture);
        double seconds = double.Parse(Environment.GetEnvironmentVariable("EMPROTO_FUZZ_SECONDS") ?? "120", CultureInfo.InvariantCulture);
        string dir = Path.Combine(Path.GetTempPath(), "emproto-ph1-fuzz-" + Guid.NewGuid().ToString("N"));
        Directory.CreateDirectory(dir);
        try
        {
            var clock = System.Diagnostics.Stopwatch.StartNew();
            var rng = new SplitMix64(seed * 0x9E3779B97F4A7C15UL);
            long cases = 0, multiLevel = 0;
            while (clock.Elapsed.TotalSeconds < seconds)
            {
                int n = new[] { 1000, 3000, 10_000, 30_000, 100_000 }[rng.NextInt(5)];
                var all = BuilderTestSupport.Families(n, rng.Next()).ToList();
                DenseCase c = all[rng.NextInt(all.Count)];
                int x = 7 + rng.NextInt(58);
                DifferenceCover cover = DifferenceCover.Create(x);
                long s = cover.SampleCount(n);
                long k = Math.Max(3 + rng.NextInt(48), s / 400); // at most about 800 buckets: the scans cost buckets x sample
                int bits = new[] { c.MinBits, 2, 3, 4, 8, 16, 32 }.Where(b => b >= c.MinBits).OrderBy(_ => rng.Next()).First();
                ExtDcsOptions o = ExtDcsHarness.Options(Path.Combine(dir, "scratch"), k, 4 + rng.NextInt(197), new[] { 1, 2, 4 }[rng.NextInt(3)], rng.Next()) with
                {
                    FanOut = 2 + rng.NextInt(5),
                    ReservoirSize = 8 + rng.NextInt(300),
                    ChunkSlots = 4 + rng.NextInt(400),
                    RunSpillBytes = 16 + rng.NextInt(600),
                    IoBufferBytes = 16 + rng.NextInt(600),
                    MergeBufferBudget = 8 + rng.NextInt(6000),
                    WindowBytes = 16 + rng.NextInt(8000),
                    SortLeafCapacity = 3 + rng.NextInt(48),
                };
                using PackedText t = PackedText.Pack(c.Symbols, bits);
                byte[] expected = ExtDcsHarness.Oracle(t, cover);
                (byte[] got, ExtDcsStats st) = ExtDcsHarness.Run(t, cover, o, dir);
                if (!got.AsSpan().SequenceEqual(expected)) Assert.Fail($"fuzz case {cases} (seed {seed}): {c.Name} X={x} K={k} b={bits} threads={o.Threads} base={o.BaseCase} F={o.FanOut}: rank file differs");
                cases++;
                if (st.MaxLevel >= 1) multiLevel++;
            }
            TestContext.Out.WriteLine($"fuzz: {cases} cases, {multiLevel} with at least one recursion level, {clock.Elapsed.TotalSeconds:F0} s");
            Assert.That(cases, Is.GreaterThan(0));
        }
        finally
        {
            Directory.Delete(dir, true);
        }
    }
}
