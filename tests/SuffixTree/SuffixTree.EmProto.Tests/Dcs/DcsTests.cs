using SuffixTree.EmProto.Build;
using SuffixTree.EmProto.Dcs;
using SuffixTree.EmProto.Harness;
using SuffixTree.EmProto.Tests.Build;
using SuffixTree.EmProto.Text;
using SuffixTree.Esa.Tests.TestSupport;
using SplitMix64 = SuffixTree.EmProto.Harness.SplitMix64;

namespace SuffixTree.EmProto.Tests.Dcs;

[TestFixture]
public class DifferenceCoverTests
{
    /// <summary>Independent check of the definition: every residue d has a, b in D with b - a = d (mod X).</summary>
    private static bool CoversEveryDifference(int x, int[] d)
    {
        var covered = new bool[x];
        foreach (int a in d)
        {
            foreach (int b in d) covered[(((b - a) % x) + x) % x] = true;
        }
        return covered.All(c => c);
    }

    [Test]
    public void EveryResidueDifferenceIsCovered_ForEveryPeriodFrom2To2300()
    {
        for (int x = 2; x <= 2300; x++)
        {
            DifferenceCover c = DifferenceCover.Create(x);
            Assert.That(c.Period, Is.EqualTo(x));
            Assert.That(CoversEveryDifference(x, c.Residues), Is.True, $"X={x} ({c.Kind}, |D|={c.Residues.Length})");
            Assert.That(c.Verify(), Is.True, $"X={x}");
            Assert.That(c.Residues, Is.Ordered.Ascending.And.Unique);
            Assert.That(c.Residues.All(r => r >= 0 && r < x), Is.True);
            for (int delta = 0; delta < x; delta++)
            {
                int a = c.Witness[delta];
                Assert.That(a, Is.GreaterThanOrEqualTo(0), $"X={x} delta={delta}");
                Assert.That(c.IndexOf[a], Is.GreaterThanOrEqualTo(0));
                Assert.That(c.IndexOf[(a + delta) % x], Is.GreaterThanOrEqualTo(0), $"witness X={x} delta={delta}");
            }
        }
    }

    [TestCase(512)]
    [TestCase(1024)]
    [TestCase(2048)]
    [TestCase(4096)]
    public void PowersOfTwoUsedByAdr004_AreCoveredAndSmall(int x)
    {
        DifferenceCover c = DifferenceCover.Create(x);
        Assert.That(CoversEveryDifference(x, c.Residues), Is.True);
        Assert.That(c.Residues.Length, Is.LessThan(3 * Math.Sqrt(x)), $"|D|={c.Residues.Length} for X={x} ({c.Kind})");
        TestContext.Out.WriteLine($"X={x}: |D|={c.Residues.Length} ({c.Density:P2}), {c.Kind}; sqrt(1.5X)={Math.Sqrt(1.5 * x):F1}");
    }

    [Test]
    public void SingerSets_ArePerfectDifferenceSets_OfSizeQPlusOne()
    {
        foreach (int q in new[] { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47 })
        {
            int x = (q * q) + q + 1;
            DifferenceCover c = DifferenceCover.Create(x);
            Assert.That(c.Kind, Is.EqualTo(CoverKind.Singer), $"q={q}");
            Assert.That(c.Residues, Has.Length.EqualTo(q + 1));
            Assert.That(c.IsPerfect(), Is.True, $"q={q} X={x}: every non-zero difference exactly once");
        }
    }

    [Test]
    public void NearestSingerPeriod_ForTheNominalPeriodsOfAdr004()
    {
        Assert.That(DifferenceCover.NearestSingerPeriod(512), Is.EqualTo(553));
        Assert.That(DifferenceCover.NearestSingerPeriod(1024), Is.EqualTo(993));
        Assert.That(DifferenceCover.NearestSingerPeriod(2048), Is.EqualTo(1893));
    }

    [Test]
    public void SampleCount_MatchesTheEnumerationOfSamplePositions()
    {
        foreach (int x in new[] { 7, 13, 20, 57 })
        {
            DifferenceCover c = DifferenceCover.Create(x);
            for (int n = 0; n < 300; n++)
            {
                long expected = 0;
                for (int p = 0; p <= n; p++) expected += c.IndexOf[p % x] >= 0 ? 1 : 0;
                Assert.That(c.SampleCount(n), Is.EqualTo(expected), $"X={x} n={n}");
            }
        }
    }
}

[TestFixture]
public class DcsRankTests
{
    private static string Text(string family, int n) => family switch
    {
        "dna" => StringFamilies.RandomDna(n, 4),
        "same" => StringFamilies.AllSame(n),
        "fib" => StringFamilies.Fibonacci(n),
        "tm" => StringFamilies.ThueMorse(n),
        "p5" => StringFamilies.Periodic(n, 5, 1),
        "runs" => StringFamilies.Runs(n, 16, 3),
        "ab" => StringFamilies.AbRepeat(n),
        _ => throw new ArgumentException(family),
    };

    [Test]
    public void SampleRanks_OrderTheSampleSuffixesLikeTheSuffixArray_IncludingRecursion()
    {
        foreach (string family in new[] { "dna", "same", "fib", "tm", "p5", "runs", "ab" })
        {
            foreach (int n in new[] { 1, 2, 5, 40, 333 })
            {
                DenseCase c = BuilderTestSupport.Dense($"{family}{n}", Text(family, n));
                uint[] sa = c.ReferenceSa();
                var rankInSa = new int[n + 1];
                for (int k = 0; k <= n; k++) rankInSa[sa[k]] = k;
                foreach (int x in new[] { 7, 13, 20, 31 })
                {
                    foreach (long baseCase in new long[] { 16_000_000, 12 })
                    {
                        DifferenceCover cover = DifferenceCover.Create(x);
                        using PackedText t = PackedText.Pack(c.Symbols, Math.Max(1, c.MinBits));
                        using DcsRankTable table = DcsBuilder.Build(t, cover, new WorkCounters(), baseCase, out DcsStats st);
                        var sample = Enumerable.Range(0, n + 1).Where(p => cover.IndexOf[p % x] >= 0).ToList();
                        Assert.That(st.SampleSize, Is.EqualTo(sample.Count));
                        foreach (int p in sample)
                        {
                            foreach (int q in sample)
                            {
                                if (p == q) continue;
                                Assert.That(table.Rank(p).CompareTo(table.Rank(q)), Is.EqualTo(rankInSa[p].CompareTo(rankInSa[q])), $"{family} n={n} X={x} base={baseCase} p={p} q={q}");
                            }
                        }
                    }
                }
            }
        }
    }

    [Test]
    public void Recursion_RunsOnRepetitiveTexts_WhenTheBaseCaseIsSmall()
    {
        DenseCase c = BuilderTestSupport.Dense("same", StringFamilies.AllSame(3000));
        DifferenceCover cover = DifferenceCover.Create(13);
        using PackedText t = PackedText.Pack(c.Symbols, 1);
        using DcsRankTable table = DcsBuilder.Build(t, cover, new WorkCounters(), 40, out DcsStats st);
        Assert.That(st.NamesUnique, Is.False);
        Assert.That(st.Levels, Is.GreaterThanOrEqualTo(2), "the reduced string of a^3000 has about 970 symbols, still above the base case of 40 after one more level");
        Assert.That(table.Rank(3000 - (3000 % 13) + cover.Residues.First(r => r <= 3000 % 13)), Is.GreaterThan(0u));
    }

    [Test]
    public void Comparator_EqualsTheNaiveSuffixOrder_ForAllPairs()
    {
        foreach (string family in new[] { "dna", "same", "fib", "tm", "p5", "runs", "ab" })
        {
            DenseCase c = BuilderTestSupport.Dense(family, Text(family, 220));
            int n = c.Symbols.Length;
            uint[] sa = c.ReferenceSa();
            var rankInSa = new int[n + 1];
            for (int k = 0; k <= n; k++) rankInSa[sa[k]] = k;
            foreach (int x in new[] { 7, 13, 24, 57 })
            {
                DifferenceCover cover = DifferenceCover.Create(x);
                foreach (int bits in new[] { Math.Max(1, c.MinBits), 8 })
                {
                    using PackedText t = PackedText.Pack(c.Symbols, bits);
                    using DcsRankTable table = DcsBuilder.Build(t, cover, new WorkCounters(), 16_000_000, out _);
                    var plain = new SuffixOrder(t, null, new WorkCounters());
                    var dcs = new SuffixOrder(t, table, new WorkCounters());
                    for (int i = 0; i < n; i++)
                    {
                        for (int j = 0; j < n; j++)
                        {
                            int expected = rankInSa[i].CompareTo(rankInSa[j]);
                            Assert.That(Math.Sign(plain.Compare((uint)i, (uint)j)), Is.EqualTo(expected), $"direct {family} i={i} j={j}");
                            Assert.That(Math.Sign(dcs.Compare((uint)i, (uint)j)), Is.EqualTo(expected), $"dcs {family} X={x} b={bits} i={i} j={j}");
                        }
                    }
                }
            }
        }
    }
}

[TestFixture]
public class DcsBuilderDifferentialTests
{
    [Test]
    public void BuildWithDcs_EqualsTheReferenceSuffixArray_ForScaledDownParameters()
    {
        var rng = new SplitMix64(4242);
        int cases = 0;
        foreach (int n in new[] { 1000, 5000 })
        {
            foreach (DenseCase c in BuilderTestSupport.Families(n, 3))
            {
                uint[] expected = c.ReferenceSa();
                for (int rep = 0; rep < 2; rep++)
                {
                    long k = 3 + rng.NextInt(48);
                    int f = 2 + rng.NextInt(3);
                    int x = 7 + rng.NextInt(58);
                    int bits = new[] { c.MinBits, 2, 3, 4, 8 }.Where(b => b >= c.MinBits).OrderBy(_ => rng.Next()).First();
                    int threads = new[] { 1, 2, 4 }[rng.NextInt(3)];
                    long baseCase = rng.NextInt(2) == 0 ? 16_000_000 : 30 + rng.NextInt(200);
                    DifferenceCover cover = DifferenceCover.Create(x);
                    using PackedText t = PackedText.Pack(c.Symbols, bits);
                    using DcsRankTable table = DcsBuilder.Build(t, cover, new WorkCounters(), baseCase, out _);
                    using var sink = new MemorySaSink(c.Symbols.Length + 1L);
                    var options = new BuildOptions { LeafCapacity = k, FanOut = f, Threads = threads, Seed = rng.Next(), Oversample = 4 };
                    _ = EmBuilder.Run(t, options, sink, table);
                    uint[] got = sink.ToArray();
                    Assert.That(got, Is.EqualTo(expected), $"{c.Name} K={k} F={f} X={x} b={bits} threads={threads} base={baseCase}: {BuilderTestSupport.Describe(expected, got)}");
                    cases++;
                }
            }
        }
        Assert.That(cases, Is.GreaterThan(100));
    }

    [Test]
    public void LargerTextsWithDcs_100kSymbols_ManyLevels()
    {
        var rng = new SplitMix64(11);
        foreach (DenseCase c in new[]
        {
            BuilderTestSupport.Dense("dna", StringFamilies.RandomDna(100_000, 9)),
            BuilderTestSupport.Dense("same", StringFamilies.AllSame(40_000)),
            BuilderTestSupport.Dense("fib", StringFamilies.Fibonacci(60_000)),
            BuilderTestSupport.Dense("p64", StringFamilies.Periodic(50_000, 64, 3)),
        })
        {
            uint[] expected = c.ReferenceSa();
            int x = 13 + rng.NextInt(50);
            DifferenceCover cover = DifferenceCover.Create(x);
            using PackedText t = PackedText.Pack(c.Symbols, c.MinBits);
            using DcsRankTable table = DcsBuilder.Build(t, cover, new WorkCounters(), 2000, out _);
            using var sink = new MemorySaSink(c.Symbols.Length + 1L);
            _ = EmBuilder.Run(t, new BuildOptions { LeafCapacity = 20, FanOut = 3, Threads = 4, Seed = 5, Oversample = 4 }, sink, table);
            Assert.That(sink.ToArray(), Is.EqualTo(expected), $"{c.Name} X={x}");
        }
    }
}
