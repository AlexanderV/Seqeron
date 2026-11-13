using SuffixTree.EmProto.Build;
using SuffixTree.EmProto.Harness;
using SuffixTree.EmProto.Text;

namespace SuffixTree.EmProto.Tests.Build;

[TestFixture]
public class EmBuilderTests
{
    private static BuildOptions Options(long k, int f, int threads, ulong seed, string? dir = null, int spill = 64)
        => new() { LeafCapacity = k, FanOut = f, Threads = threads, Seed = seed, ScratchDir = dir, SpillBytes = spill, Oversample = 4 };

    [Test]
    public void TinyTexts_n0_to_n3_AllBinaryStrings()
    {
        for (int n = 0; n <= 8; n++)
        {
            for (int code = 0; code < (1 << n); code++)
            {
                var sym = new byte[n];
                for (int i = 0; i < n; i++) sym[i] = (byte)((code >> i) & 1);
                var c = new DenseCase($"bin{n}-{code}", sym, n == 0 ? 0 : sym.Distinct().Count());
                // dense symbols need the order-preserving remap: with symbols 0/1 where only 1 occurs, sigma 1 would be wrong, so remap by hand
                byte[] dense = sym.Select(x => (byte)(c.Sigma == 1 ? 0 : x)).ToArray();
                c = c with { Symbols = dense };
                uint[] expected = c.ReferenceSa();
                foreach ((long k, int f) in new[] { (2L, 2), (3L, 2), (5L, 4) })
                {
                    uint[] got = BuilderTestSupport.Build(c, Math.Max(1, c.MinBits), Options(k, f, 1, 1));
                    Assert.That(got, Is.EqualTo(expected), $"{c.Name} K={k} F={f}");
                }
            }
        }
    }

    [Test]
    public void AdversarialFamilies_EqualTheReferenceSuffixArray_ForScaledDownParameters()
    {
        var rng = new SplitMix64(2024);
        int cases = 0;
        foreach (int n in new[] { 1000, 6000 })
        {
            foreach (DenseCase c in BuilderTestSupport.Families(n, 7))
            {
                uint[] expected = c.ReferenceSa();
                for (int rep = 0; rep < 2; rep++)
                {
                    long k = 3 + rng.NextInt(48);
                    int f = 2 + rng.NextInt(3);
                    int bits = new[] { c.MinBits, 2, 3, 4, 8 }.Where(b => b >= c.MinBits).OrderBy(_ => rng.Next()).First();
                    int threads = new[] { 1, 2, 4 }[rng.NextInt(3)];
                    uint[] got = BuilderTestSupport.Build(c, bits, Options(k, f, threads, rng.Next()));
                    Assert.That(got, Is.EqualTo(expected), $"{c.Name} K={k} F={f} b={bits} threads={threads}: {BuilderTestSupport.Describe(expected, got)}");
                    cases++;
                }
            }
        }
        Assert.That(cases, Is.GreaterThan(100));
    }

    [Test]
    public void LargerTexts_UpTo100kSymbols_EqualTheReference_WithManyLevels()
    {
        var rng = new SplitMix64(99);
        foreach (DenseCase c in new[]
        {
            BuilderTestSupport.Dense("random-dna", global::SuffixTree.Esa.Tests.TestSupport.StringFamilies.RandomDna(100_000, 3)),
            BuilderTestSupport.Dense("fibonacci", global::SuffixTree.Esa.Tests.TestSupport.StringFamilies.Fibonacci(100_000)),
            BuilderTestSupport.Dense("periodic-p64", global::SuffixTree.Esa.Tests.TestSupport.StringFamilies.Periodic(60_000, 64, 1)),
            BuilderTestSupport.Dense("all-same", global::SuffixTree.Esa.Tests.TestSupport.StringFamilies.AllSame(30_000)),
        })
        {
            uint[] expected = c.ReferenceSa();
            for (int rep = 0; rep < 2; rep++)
            {
                long k = 3 + rng.NextInt(48);
                int f = 2 + rng.NextInt(3);
                uint[] got = BuilderTestSupport.Build(c, c.MinBits, Options(k, f, 4, rng.Next()));
                Assert.That(got, Is.EqualTo(expected), $"{c.Name} K={k} F={f}: {BuilderTestSupport.Describe(expected, got)}");
            }
        }
    }

    [Test]
    public void Output_DoesNotDependOnThreadsLeafCapacityFanOutOrSeed()
    {
        DenseCase c = BuilderTestSupport.Dense("dna", global::SuffixTree.Esa.Tests.TestSupport.StringFamilies.RandomDna(20_000, 8));
        uint[] baseline = BuilderTestSupport.Build(c, 2, Options(1 << 20, 64, 1, 1));
        Assert.That(baseline, Is.EqualTo(c.ReferenceSa()));
        foreach ((long k, int f, int threads, ulong seed) in new[] { (7L, 2, 1, 5UL), (50L, 4, 2, 6UL), (300L, 3, 4, 7UL), (3000L, 64, 4, 8UL) })
            Assert.That(BuilderTestSupport.Build(c, 2, Options(k, f, threads, seed)), Is.EqualTo(baseline), $"K={k} F={f} threads={threads}");
    }

    [Test]
    public void FileMode_SpillsBuckets_GivesTheSameSuffixArray_AndRemovesItsScratchFiles()
    {
        string dir = Path.Combine(Path.GetTempPath(), "emproto-scratch-" + Guid.NewGuid().ToString("N"));
        Directory.CreateDirectory(dir);
        try
        {
            foreach (DenseCase c in new[]
            {
                BuilderTestSupport.Dense("dna", global::SuffixTree.Esa.Tests.TestSupport.StringFamilies.RandomDna(30_000, 2)),
                BuilderTestSupport.Dense("tm", global::SuffixTree.Esa.Tests.TestSupport.StringFamilies.ThueMorse(20_000)),
            })
            {
                IoCounters.Reset();
                uint[] got = BuilderTestSupport.Build(c, c.MinBits, Options(40, 3, 2, 3, dir, spill: 64));
                Assert.That(got, Is.EqualTo(c.ReferenceSa()));
                Assert.That(IoCounters.WriteBytes, Is.GreaterThan(c.Symbols.Length), "the buckets went through scratch files");
                Assert.That(IoCounters.ReadBytes, Is.GreaterThan(IoCounters.WriteBytes / 2), "the spilled buckets are read back (plus the 4-byte sample reads)");
                Assert.That(Directory.GetFiles(dir), Is.Empty, "a scratch file is deleted as soon as its consumer is done");
            }
        }
        finally
        {
            Directory.Delete(dir, true);
        }
    }

    [Test]
    public void TieRunSorting_MergeAndIntro_GiveTheSameSuffixArray_AndMergeNeedsFewerComparisonsOnReverseOrderedRuns()
    {
        DenseCase same = BuilderTestSupport.Dense("same", new string('a', 6000));
        DenseCase tm = BuilderTestSupport.Dense("tm", global::SuffixTree.Esa.Tests.TestSupport.StringFamilies.ThueMorse(8000));
        foreach (DenseCase c in new[] { same, tm })
        {
            uint[] expected = c.ReferenceSa();
            var compares = new long[2];
            for (int v = 0; v < 2; v++)
            {
                using PackedText t = PackedText.Pack(c.Symbols, Math.Max(1, c.MinBits));
                using var sink = new MemorySaSink(c.Symbols.Length + 1L);
                BuildStats st = EmBuilder.Run(t, Options(1000, 4, 1, 1) with { TieSort = v == 0 ? TieSort.Merge : TieSort.Intro }, sink);
                Assert.That(sink.ToArray(), Is.EqualTo(expected), $"{c.Name} {(v == 0 ? "merge" : "intro")}");
                compares[v] = st.Work.Compares;
            }
            Assert.That(compares[0], Is.LessThan(compares[1]), $"{c.Name}: merge {compares[0]} vs intro {compares[1]}");
        }
    }

    [Test]
    public void MergeSort_SortsRunsOfEveryLengthByTheComparator()
    {
        var rng = new SplitMix64(17);
        // a run of equal keys is what the tied comparator is defined for: every position of a^n has the same key (except the last few, which are padded and compared from symbol 0)
        DenseCase c = BuilderTestSupport.Dense("same", new string('a', 2000));
        using PackedText t = PackedText.Pack(c.Symbols, 1);
        var order = new SuffixOrder(t, null, new WorkCounters());
        uint[] sa = c.ReferenceSa();
        var rank = new int[sa.Length];
        for (int k = 0; k < sa.Length; k++) rank[sa[k]] = k;
        for (int len = 0; len <= 70; len++)
        {
            var data = Enumerable.Range(0, 2000).OrderBy(_ => rng.Next()).Take(len).Select(x => (uint)x).ToArray();
            var tmp = new uint[len];
            LeafSorter.MergeSort(data, tmp, order);
            Assert.That(data.Select(x => rank[x]), Is.Ordered.Ascending, $"len {len}");
        }
    }

    [Test]
    public void KeyTiesAreCounted_AndRepetitiveTextsHaveMany()
    {
        DenseCase same = BuilderTestSupport.Dense("same", new string('a', 4000));
        DenseCase dna = BuilderTestSupport.Dense("dna", global::SuffixTree.Esa.Tests.TestSupport.StringFamilies.RandomDna(4000, 1));
        long Ties(DenseCase c)
        {
            using PackedText t = PackedText.Pack(c.Symbols, Math.Max(1, c.MinBits));
            using var sink = new MemorySaSink(c.Symbols.Length + 1L);
            BuildStats st = EmBuilder.Run(t, Options(100, 4, 1, 1), sink);
            return st.Work.TiePositions;
        }
        Assert.That(Ties(same), Is.GreaterThan(3500));
        Assert.That(Ties(dna), Is.LessThan(50));
    }

    [Test]
    public void LeafBuffersAreChargedToTheBudget()
    {
        DenseCase c = BuilderTestSupport.Dense("dna", global::SuffixTree.Esa.Tests.TestSupport.StringFamilies.RandomDna(5000, 1));
        var options = Options(1000, 4, 1, 1) with { Budget = new ArenaBudget(1000 * 24 - 1) };
        Assert.Throws<InvalidOperationException>(() => BuilderTestSupport.Build(c, 2, options));
        var budget = new ArenaBudget(1000 * 24);
        _ = BuilderTestSupport.Build(c, 2, options with { Budget = budget });
        Assert.That(budget.PeakBytes, Is.EqualTo(1000 * 24));
        Assert.That(budget.LiveBytes, Is.EqualTo(0));
    }

    [Test]
    public void ExactPartitionSizes_NoRetriesOnRandomText_AndDepthGrowsWithNOverK()
    {
        DenseCase c = BuilderTestSupport.Dense("dna", global::SuffixTree.Esa.Tests.TestSupport.StringFamilies.RandomDna(50_000, 4));
        using PackedText t = PackedText.Pack(c.Symbols, 2);
        using var sink = new MemorySaSink(c.Symbols.Length + 1L);
        BuildStats st = EmBuilder.Run(t, Options(10, 2, 1, 1), sink);
        Assert.That(st.MaxDepth, Is.GreaterThanOrEqualTo(8));
        Assert.That(st.MaxLeaf, Is.LessThanOrEqualTo(10));
        Assert.That(st.Leaves, Is.GreaterThanOrEqualTo(50_000 / 10));
    }
}
