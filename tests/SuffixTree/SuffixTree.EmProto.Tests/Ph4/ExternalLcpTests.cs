using System.Globalization;
using SuffixTree.Esa;
using SuffixTree.EmProto.Harness;
using SuffixTree.EmProto.Ph4;
using SuffixTree.EmProto.Tests.Build;
using SuffixTree.EmProto.Text;
using SuffixTree.Esa.Tests.TestSupport;
using SplitMix64 = SuffixTree.EmProto.Harness.SplitMix64;

namespace SuffixTree.EmProto.Tests.Ph4;

/// <summary>
/// Gate G1 of the prototype E4 (ADR-004 section 5, row E4, "LCP equal to LcpKernel"): every string of a small alphabet up to a length, in both modes, with tiny chunks so that the chunk boundaries of the
/// text and of the ranks are crossed; the compact form (width, base values, overflow ranks and values) equal to the kernel and the decoded values equal to the safe reference.
/// </summary>
[TestFixture]
public class ExternalLcpExhaustiveTests
{
    private string _dir = string.Empty;

    [SetUp]
    public void SetUp()
    {
        _dir = Path.Combine(Path.GetTempPath(), "emproto-ph4-" + Guid.NewGuid().ToString("N"));
        Directory.CreateDirectory(_dir);
    }

    [TearDown]
    public void TearDown() => Directory.Delete(_dir, true);

    private void CheckAll(int sigma, int maxLength, LcpMode mode, long[] chunkSizes)
    {
        string scratch = Path.Combine(_dir, "scratch");
        int bits = Math.Max(1, PackedText.BitsFor(sigma));
        long cases = 0, multiChunk = 0, multiRound = 0;
        foreach (byte[] s in LcpHarness.AllStrings(sigma, 1, maxLength))
        {
            var c = new DenseCase("s", s, sigma);
            long r = chunkSizes[(int)(cases % chunkSizes.Length)];
            int threads = new[] { 1, 2, 4 }[(int)(cases % 3)];
            long rounds = mode == LcpMode.Bounded ? 1 + (cases % 2) : 0;
            LcpOutcome got = LcpHarness.CheckAgainstOracles(c, bits, LcpHarness.Options(scratch, r, threads, mode, rounds), _dir, $"{mode} sigma={sigma} text=[{string.Join(',', s)}] R={r} threads={threads} rounds={rounds}");
            cases++;
            if (got.Stats.TextChunks > 1) multiChunk++;
            if (got.Stats.Rounds > 1) multiRound++;
        }
        TestContext.Out.WriteLine($"{mode} sigma={sigma} length<={maxLength}: {cases} strings, {multiChunk} with several chunks, {multiRound} with several rounds");
        Assert.That(cases, Is.GreaterThan(0));
        Assert.That(multiChunk, Is.GreaterThan(cases / 2), "most strings must cross chunk boundaries");
        if (mode == LcpMode.Bounded) Assert.That(multiRound, Is.GreaterThan(cases / 5), "the bounded mode must run several rounds");
    }

    [Test]
    public void AllBinaryStrings_Length1To12_Distributed_TinyChunks_EqualKernelAndReference() => CheckAll(2, 12, LcpMode.Distributed, [3, 5]);

    [Test]
    public void AllBinaryStrings_Length1To11_Bounded_TinyChunksAndRounds_EqualKernelAndReference() => CheckAll(2, 11, LcpMode.Bounded, [3, 4]);

    [Test]
    public void AllTernaryStrings_Length1To6_Distributed_EqualKernelAndReference() => CheckAll(3, 6, LcpMode.Distributed, [2, 3, 7]);

    [Test]
    public void AllTernaryStrings_Length1To6_Bounded_EqualKernelAndReference() => CheckAll(3, 6, LcpMode.Bounded, [2, 5]);
}

[TestFixture]
public class ExternalLcpDifferentialTests
{
    private string _dir = string.Empty;

    [SetUp]
    public void SetUp()
    {
        _dir = Path.Combine(Path.GetTempPath(), "emproto-ph4-" + Guid.NewGuid().ToString("N"));
        Directory.CreateDirectory(_dir);
    }

    [TearDown]
    public void TearDown() => Directory.Delete(_dir, true);

    private string Scratch => Path.Combine(_dir, "scratch");

    [Test]
    public void AdversarialFamilies_ScaledDownParameters_EqualBothOracles()
    {
        var rng = new SplitMix64(4040);
        int cases = 0, distributed = 0, bounded = 0, multiChunk = 0, scalar = 0, multiRound = 0;
        var kinds = new HashSet<LcpKind>();
        foreach (int n in new[] { 1000, 5000 })
        {
            foreach (DenseCase c in BuilderTestSupport.Families(n, 5))
            {
                for (int rep = 0; rep < 2; rep++)
                {
                    long r = Math.Max(3 + rng.NextInt(48), n / 400);
                    LcpMode mode = rng.NextInt(2) == 0 ? LcpMode.Distributed : LcpMode.Bounded;
                    int bits = new[] { c.MinBits, 2, 3, 4, 8, 16, 32 }.Where(b => b >= c.MinBits).OrderBy(_ => rng.Next()).First();
                    ExtLcpOptions o = LcpHarness.Options(Scratch, r, new[] { 1, 2, 4 }[rng.NextInt(3)], mode, mode == LcpMode.Bounded ? 1 + rng.NextInt(4) : 0) with
                    {
                        ReadBufferBytes = 16 + rng.NextInt(400),
                        TupleBufferBytes = 12 + rng.NextInt(300),
                        WindowBytes = 16 + rng.NextInt(4000),
                        TouchBuffers = rng.NextInt(4) == 0,
                    };
                    bool useScalar = rng.NextInt(3) == 0;
                    LcpOutcome got = LcpHarness.CheckAgainstOracles(c, bits, o, _dir, $"{c.Name} {mode} R={r} b={bits} threads={o.Threads} rounds={o.RoundChunks} read={o.ReadBufferBytes} tuple={o.TupleBufferBytes} window={o.WindowBytes} scalar={useScalar}", scalar: useScalar);
                    cases++;
                    if (mode == LcpMode.Distributed) distributed++;
                    else bounded++;
                    if (got.Stats.TextChunks > 1) multiChunk++;
                    if (got.Stats.Rounds > 1) multiRound++;
                    if (useScalar) scalar++;
                    kinds.Add(got.Stats.Kind);
                }
            }
        }
        Assert.That(cases, Is.GreaterThan(100));
        Assert.That(distributed, Is.GreaterThan(30));
        Assert.That(bounded, Is.GreaterThan(30));
        Assert.That(multiChunk, Is.GreaterThan(cases * 9 / 10), "many chunks");
        Assert.That(multiRound, Is.GreaterThan(10), "several rounds of the bounded mode");
        Assert.That(scalar, Is.GreaterThan(10), "the scalar compare path");
        Assert.That(kinds, Does.Contain(LcpKind.U8));
        Assert.That(kinds, Does.Contain(LcpKind.U16), "periodic texts of 5000 symbols have values above 255 in plenty and none above 65,535: U16 without overflow");
    }

    /// <summary>A random 5000-symbol text with a 600-symbol block copied once: a few hundred values above 255, far fewer than n / 16.</summary>
    private static string U8WithOverflowText()
    {
        string head = StringFamilies.RandomDna(5000, 3);
        return string.Concat(head, head.AsSpan(1000, 600), StringFamilies.RandomDna(800, 4));
    }

    /// <summary>X X with a random X of 70,000 symbols: 4,466 values above 65,535 (at most n / 8), about 70,000 above 255 (more than n / 16).</summary>
    private static string U16WithOverflowText()
    {
        string x = StringFamilies.RandomDna(70_000, 5);
        return x + x;
    }

    [TestCase("u8", false, 400)]
    [TestCase("u8", true, 300)]
    [TestCase("u16", false, 700)]
    [TestCase("u16", true, 9000)]
    [TestCase("u32", false, 900)]
    [TestCase("u32", true, 20000)]
    public void Widths_U8U16U32_WithOverflowLists_EqualKernel(string width, bool boundedMode, long r)
    {
        LcpMode mode = boundedMode ? LcpMode.Bounded : LcpMode.Distributed;
        (string text, LcpKind expected, bool overflow) = width switch
        {
            "u8" => (U8WithOverflowText(), LcpKind.U8, true),
            "u16" => (U16WithOverflowText(), LcpKind.U16, true),
            _ => (StringFamilies.AllSame(150_000), LcpKind.U32, false),
        };
        DenseCase c = BuilderTestSupport.Dense(width, text);
        LcpOutcome got = LcpHarness.CheckAgainstOracles(c, c.MinBits, LcpHarness.Options(Scratch, r, 4, mode, mode == LcpMode.Bounded ? 3 : 0), _dir, $"{width} {mode} R={r}", reference: false);
        Assert.That(got.Stats.Kind, Is.EqualTo(expected));
        if (overflow) Assert.That(got.Stats.OverflowCount, Is.GreaterThan(100));
        else Assert.That(got.Stats.OverflowCount, Is.Zero, "U32 has no overflow list");
        Assert.That(got.Stats.OverflowCount, Is.EqualTo(got.Compact.OverflowCount));
        TestContext.Out.WriteLine($"{width} {mode}: kind {got.Stats.Kind}, escapes8 {got.Stats.Escapes8}, escapes16 {got.Stats.Escapes16}, overflow {got.Stats.OverflowCount}, chunks {got.Stats.TextChunks}");
    }

    [Test]
    public void LargerTexts_ManyChunks_EqualKernel()
    {
        var rng = new SplitMix64(12);
        foreach (DenseCase c in new[]
        {
            BuilderTestSupport.Dense("dna", StringFamilies.RandomDna(60_000, 9)),
            BuilderTestSupport.Dense("fib", StringFamilies.Fibonacci(30_000)),
            BuilderTestSupport.Dense("p64", StringFamilies.Periodic(30_000, 64, 3)),
            BuilderTestSupport.Dense("runs", StringFamilies.Runs(30_000, 300, 4)),
            BuilderTestSupport.Dense("tm", StringFamilies.ThueMorse(30_000)),
        })
        {
            foreach (LcpMode mode in new[] { LcpMode.Distributed, LcpMode.Bounded })
            {
                long r = 100 + rng.NextInt(900);
                LcpOutcome got = LcpHarness.CheckAgainstOracles(c, c.MinBits, LcpHarness.Options(Scratch, r, 4, mode, 4), _dir, $"{c.Name} {mode} R={r}", reference: false);
                Assert.That(got.Stats.TextChunks, Is.GreaterThan(20));
            }
        }
    }

    [Test]
    public void Output_DoesNotDependOnThreadsChunkSizeBuffersModeOrRounds()
    {
        foreach (DenseCase c in new[]
        {
            BuilderTestSupport.Dense("binary", StringFamilies.Random(3000, 2, 6)),
            BuilderTestSupport.Dense("fib", StringFamilies.Fibonacci(3000)),
            BuilderTestSupport.Dense("same", StringFamilies.AllSame(3000)),
            BuilderTestSupport.Dense("dna", StringFamilies.RandomDna(3000, 7)),
        })
        {
            uint[] sa = c.ReferenceSa();
            using PackedText t = PackedText.Pack(c.Symbols, c.MinBits);
            LcpOutcome first = LcpHarness.Run(t, sa, LcpHarness.Options(Scratch, 3000), _dir);
            foreach (LcpMode mode in new[] { LcpMode.Distributed, LcpMode.Bounded })
            {
                foreach (int threads in new[] { 1, 2, 4 })
                {
                    foreach (long r in new long[] { 3, 17, 200, 1 << 20 })
                    {
                        foreach (int buffer in new[] { 12, 100, 4096 })
                        {
                            ExtLcpOptions o = LcpHarness.Options(Scratch, r, threads, mode, threads) with { ReadBufferBytes = Math.Max(16, buffer), TupleBufferBytes = buffer, WindowBytes = buffer, TouchBuffers = buffer == 100 };
                            if (r == 3 && (buffer != 12 || threads != 4)) continue; // 1000 chunks per spool: enough once
                            LcpOutcome got = LcpHarness.Run(t, sa, o, _dir);
                            Assert.That(got.BaseBytes.AsSpan().SequenceEqual(first.BaseBytes), Is.True, $"{c.Name} {mode} threads={threads} R={r} buffer={buffer}: base file");
                            Assert.That(got.OverflowBytes.AsSpan().SequenceEqual(first.OverflowBytes), Is.True, $"{c.Name} {mode} threads={threads} R={r} buffer={buffer}: overflow file");
                            Assert.That(got.Stats.Kind, Is.EqualTo(first.Stats.Kind));
                        }
                    }
                }
            }
        }
    }

    private sealed class RecordingConsumer : ILcpChunkConsumer
    {
        internal long Next;
        internal int Calls;
        internal readonly List<uint> Values = [];
        internal readonly List<string> Violations = [];

        public void Consume(long firstRank, ReadOnlySpan<uint> lcp)
        {
            if (firstRank != Next) Violations.Add($"chunk starts at rank {firstRank}, expected {Next}");
            Next = firstRank + lcp.Length;
            Calls++;
            Values.AddRange(lcp.ToArray());
        }
    }

    [TestCase(false)]
    [TestCase(true)]
    public void ChunkConsumer_ReceivesEveryRankOnceInAscendingOrder_WithResolvedValues(bool boundedMode)
    {
        LcpMode mode = boundedMode ? LcpMode.Bounded : LcpMode.Distributed;
        DenseCase c = BuilderTestSupport.Dense("periodic", StringFamilies.Periodic(2000, 37, 8) + StringFamilies.AllSame(900) + StringFamilies.RandomDna(500, 2));
        uint[] sa = c.ReferenceSa();
        using PackedText t = PackedText.Pack(c.Symbols, c.MinBits);
        var consumer = new RecordingConsumer();
        LcpOutcome got = LcpHarness.Run(t, sa, LcpHarness.Options(Scratch, 77, 4, mode, 2), _dir, consumer);
        Assert.That(consumer.Violations, Is.Empty);
        Assert.That(consumer.Next, Is.EqualTo(c.Symbols.Length + 1L));
        Assert.That(consumer.Calls, Is.EqualTo(got.Stats.RankChunks));
        int[] expected = LcpHarness.Reference(c, sa).ToArray();
        Assert.That(consumer.Values.Select(v => (int)v).ToArray().AsSpan().SequenceEqual(expected), Is.True, "resolved LCP values in rank order");
        Assert.That(got.Stats.Escapes8, Is.GreaterThan(0), "the text has values above 255, so the stream carries resolved values, not escapes");
    }

    [Test]
    public void Stats_CountEscapesAndChunks()
    {
        DenseCase c = BuilderTestSupport.Dense("mix", StringFamilies.Periodic(3000, 41, 8) + StringFamilies.RandomDna(2000, 2));
        uint[] sa = c.ReferenceSa();
        int[] lcp = LcpHarness.Reference(c, sa).ToArray();
        using PackedText t = PackedText.Pack(c.Symbols, c.MinBits);
        foreach (LcpMode mode in new[] { LcpMode.Distributed, LcpMode.Bounded })
        {
            LcpOutcome got = LcpHarness.Run(t, sa, LcpHarness.Options(Scratch, 330, 2, mode, 3), _dir);
            Assert.That(got.Stats.N, Is.EqualTo(5000));
            Assert.That(got.Stats.Escapes8, Is.EqualTo(lcp.Count(v => v >= 255)));
            Assert.That(got.Stats.Escapes16, Is.EqualTo(lcp.Count(v => v >= 65535)));
            Assert.That(got.Stats.MaxPlcp, Is.EqualTo(lcp.Max()));
            Assert.That(got.Stats.TextChunks, Is.EqualTo((5000 + 329) / 330));
            Assert.That(got.Stats.RankChunks, Is.EqualTo((5001 + 329) / 330));
            Assert.That(got.Stats.Rounds, Is.EqualTo(mode == LcpMode.Bounded ? (16 + 2) / 3 : 1));
            Assert.That(got.Stats.CompareSteps, Is.GreaterThan(0));
        }
    }

    [Test]
    public void ArenaBudget_IsBalancedAndEnforced_AndFailureLeavesNoScratch()
    {
        DenseCase c = BuilderTestSupport.Dense("dna", StringFamilies.RandomDna(4000, 3));
        uint[] sa = c.ReferenceSa();
        using PackedText t = PackedText.Pack(c.Symbols, c.MinBits);
        foreach (LcpMode mode in new[] { LcpMode.Distributed, LcpMode.Bounded })
        {
            var budget = new ArenaBudget(8 << 20);
            LcpHarness.Run(t, sa, LcpHarness.Options(Scratch, 100, 4, mode, 5) with { Budget = budget }, _dir);
            Assert.That(budget.LiveBytes, Is.Zero, "every native buffer returned to the budget");
            Assert.That(budget.PeakBytes, Is.GreaterThan(0));

            var tiny = new ArenaBudget(64);
            string saPath = LcpHarness.WriteSa(sa, _dir);
            Assert.Throws<ArenaBudgetExceededException>(() => ExternalLcp.Build(saPath, t, Path.Combine(_dir, "b.bin"), Path.Combine(_dir, "o.bin"), LcpHarness.Options(Scratch, 100, 4, mode, 5) with { Budget = tiny }));
            Assert.That(tiny.LiveBytes, Is.Zero);
            Assert.That(Directory.Exists(Scratch) ? Directory.GetFileSystemEntries(Scratch) : [], Is.Empty);
        }
    }

    [TestCase(false)]
    [TestCase(true)]
    public void InvalidSuffixArrays_Throw_AndLeaveNoScratchAndNoOutput(bool boundedMode)
    {
        LcpMode mode = boundedMode ? LcpMode.Bounded : LcpMode.Distributed;
        DenseCase c = BuilderTestSupport.Dense("dna", StringFamilies.RandomDna(900, 3));
        uint[] good = c.ReferenceSa();
        using PackedText t = PackedText.Pack(c.Symbols, c.MinBits);
        string basePath = Path.Combine(_dir, "b.bin");
        string ovfPath = Path.Combine(_dir, "o.bin");

        uint[] noEmpty = (uint[])good.Clone();
        noEmpty[0] = 5; // sa[0] must be n
        uint[] outOfRange = (uint[])good.Clone();
        outOfRange[^1] = (uint)c.Symbols.Length + 5; // the last tuple, after most of the scratch was written
        uint[] negative = (uint[])good.Clone();
        negative[10] = uint.MaxValue;
        uint[] equalsN = (uint[])good.Clone();
        equalsN[10] = (uint)c.Symbols.Length; // n is the empty suffix, it belongs to rank 0 only
        foreach (uint[] bad in new[] { noEmpty, outOfRange, negative, equalsN })
        {
            string saPath = LcpHarness.WriteSa(bad, _dir);
            Assert.Throws<ArgumentException>(() => ExternalLcp.Build(saPath, t, basePath, ovfPath, LcpHarness.Options(Scratch, 60, 2, mode, 3)));
            Assert.That(Directory.Exists(Scratch) ? Directory.GetFileSystemEntries(Scratch) : [], Is.Empty, "scratch files removed after a failure");
            Assert.That(File.Exists(basePath) || File.Exists(ovfPath), Is.False, "no output before the width is known");
        }

        string shortSa = LcpHarness.WriteSa(good[..^1], _dir);
        Assert.Throws<ArgumentException>(() => ExternalLcp.Build(shortSa, t, basePath, ovfPath, LcpHarness.Options(Scratch, 60, 2, mode, 3)), "n + 1 entries are required");
    }

    [TestCase(0)]
    [TestCase(1)]
    [TestCase(2)]
    [TestCase(3)]
    public void EmptyAndTinyTexts_AllSame(int n)
    {
        foreach (LcpMode mode in new[] { LcpMode.Distributed, LcpMode.Bounded })
        {
            var symbols = new byte[n];
            uint[] sa = Enumerable.Range(0, n + 1).Select(i => (uint)(n - i)).ToArray(); // the suffix array of a^n
            using PackedText t = PackedText.Pack(symbols, 1);
            LcpOutcome got = LcpHarness.Run(t, sa, LcpHarness.Options(Scratch, 2, 2, mode, 1), _dir);
            Assert.That(got.Stats.Kind, Is.EqualTo(LcpKind.U8));
            int[] expected = Enumerable.Range(0, n + 1).Select(r => Math.Max(0, r - 1)).ToArray(); // LCP[r] = r - 1 for r >= 1
            Assert.That(got.Compact.ToArray(), Is.EqualTo(expected));
            Assert.That(got.Stats.OverflowCount, Is.Zero);
            Assert.That(got.Stats.N, Is.EqualTo(n));
        }
    }
}

/// <summary>
/// Differential fuzz of PH4 (ADR-004 section 5, "Oracles per stage"): scaled-down R = 3..50 on texts of 10^3..10^5 symbols from the adversarial families, a periodic text with long repeats and two blocks that
/// force the widths U16 and U32, against <see cref="LcpKernel"/> and <see cref="LcpReference"/>. Category <c>Fuzz</c>, <c>[Explicit]</c>; environment <c>EMPROTO_FUZZ_SEED</c> (default 1) and
/// <c>EMPROTO_FUZZ_SECONDS</c> (the budget, default 120).
/// </summary>
[TestFixture]
[Category("Fuzz")]
[Explicit("Differential fuzz of PH4; run with --filter Category=Fuzz (see PROGRESS-S1.md, Log E4)")]
public class ExternalLcpFuzzTests
{
    [Test]
    public void ScaledDownParameters_Texts1kTo100k_EqualBothOracles_AllWidths()
    {
        ulong seed = ulong.Parse(Environment.GetEnvironmentVariable("EMPROTO_FUZZ_SEED") ?? "1", CultureInfo.InvariantCulture);
        double seconds = double.Parse(Environment.GetEnvironmentVariable("EMPROTO_FUZZ_SECONDS") ?? "120", CultureInfo.InvariantCulture);
        string dir = Path.Combine(Path.GetTempPath(), "emproto-ph4-fuzz-" + Guid.NewGuid().ToString("N"));
        Directory.CreateDirectory(dir);
        try
        {
            var clock = System.Diagnostics.Stopwatch.StartNew();
            var rng = new SplitMix64(seed * 0x9E3779B97F4A7C15UL);
            long cases = 0, bounded = 0;
            var kinds = new Dictionary<LcpKind, long>();
            while (clock.Elapsed.TotalSeconds < seconds)
            {
                DenseCase c;
                int pick = rng.NextInt(40);
                if (pick == 0)
                {
                    string x = StringFamilies.RandomDna(66_000 + rng.NextInt(20_000), rng.Next());
                    c = BuilderTestSupport.Dense("xx", x + x);
                }
                else if (pick == 1)
                {
                    c = BuilderTestSupport.Dense("same", StringFamilies.AllSame(140_000 + rng.NextInt(30_000)));
                }
                else
                {
                    int n = new[] { 1000, 3000, 10_000, 30_000, 100_000 }[rng.NextInt(5)];
                    var all = BuilderTestSupport.Families(n, rng.Next()).ToList();
                    c = all[rng.NextInt(all.Count)];
                }
                long r = Math.Max(3 + rng.NextInt(48), c.Symbols.Length / 3000);
                LcpMode mode = rng.NextInt(2) == 0 ? LcpMode.Distributed : LcpMode.Bounded;
                int bits = new[] { c.MinBits, 2, 3, 4, 8, 16, 32 }.Where(b => b >= c.MinBits).OrderBy(_ => rng.Next()).First();
                ExtLcpOptions o = LcpHarness.Options(Path.Combine(dir, "scratch"), r, new[] { 1, 2, 4 }[rng.NextInt(3)], mode, mode == LcpMode.Bounded ? 1 + rng.NextInt(6) : 0) with
                {
                    ReadBufferBytes = 16 + rng.NextInt(600),
                    TupleBufferBytes = 12 + rng.NextInt(400),
                    WindowBytes = 16 + rng.NextInt(8000),
                };
                LcpOutcome got = LcpHarness.CheckAgainstOracles(c, bits, o, dir, $"fuzz case {cases} (seed {seed}): {c.Name} n={c.Symbols.Length} {mode} R={r} b={bits} threads={o.Threads}", reference: c.Symbols.Length <= 30_000, scalar: rng.NextInt(4) == 0);
                cases++;
                if (mode == LcpMode.Bounded) bounded++;
                kinds[got.Stats.Kind] = kinds.GetValueOrDefault(got.Stats.Kind) + 1;
            }
            TestContext.Out.WriteLine($"fuzz: {cases} cases ({bounded} bounded), kinds {string.Join(", ", kinds.Select(k => $"{k.Key} {k.Value}"))}, {clock.Elapsed.TotalSeconds:F0} s");
            Assert.That(cases, Is.GreaterThan(0));
        }
        finally
        {
            Directory.Delete(dir, true);
        }
    }
}
