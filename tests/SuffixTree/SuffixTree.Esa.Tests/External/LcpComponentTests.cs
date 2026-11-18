using SuffixTree.Esa.External;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.External;

/// <summary>The width rule, the chunk arithmetic and the parameter rule of PH4 + PH5 (work package 1b.7c).</summary>
[TestFixture]
public class LcpRuleTests
{
    // ---- LcpWidth

    [TestCase(0, 0, 0, "U8")]
    [TestCase(16, 1, 0, "U8")]
    [TestCase(15, 1, 0, "U16")] // 16 * 1 > 15, 8 * 0 <= 15
    [TestCase(16, 1, 1, "U8")] // the first test passes: 16 * 1 <= 16
    [TestCase(17, 2, 1, "U16")]
    [TestCase(1_000, 62, 0, "U8")]
    [TestCase(1_000, 63, 0, "U16")]
    [TestCase(1_000, 63, 125, "U16")]
    [TestCase(1_000, 200, 126, "U32")]
    [TestCase(1_000, 200, 125, "U16")]
    [TestCase(1_000, 1_000, 1_000, "U32")]
    public void Select_FollowsTheRuleOfAdr002_AtItsBoundaries(long n, long e8, long e16, string expected)
    {
        // 16 * escapes8 <= n selects U8; else 8 * escapes16 <= n selects U16; else U32
        Assert.That(LcpWidth.Select(n, e8, e16), Is.EqualTo(Enum.Parse<LcpKind>(expected)));
    }

    [Test]
    public void Select_EqualsCompactLcpSelectKind_OnAGrid()
    {
        var rng = new SplitMix64(7);
        for (int i = 0; i < 20_000; i++)
        {
            int n = rng.NextInt(3000);
            long e8 = rng.NextInt(n + 1);
            long e16 = rng.NextInt((int)e8 + 1);
            Assert.That(LcpWidth.Select(n, e8, e16), Is.EqualTo(CompactLcp.SelectKind(n, e8, e16)), $"n={n} e8={e8} e16={e16}");
        }
        for (int n = 0; n <= 40; n++)
        {
            for (int e8 = 0; e8 <= n + 2; e8++)
            {
                for (int e16 = 0; e16 <= e8; e16++) Assert.That(LcpWidth.Select(n, e8, e16), Is.EqualTo(CompactLcp.SelectKind(n, e8, e16)), $"n={n} e8={e8} e16={e16}");
            }
        }
    }

    [Test]
    public void Select_WorksAbove2Pow31_WhereCompactLcpTakesAnInt()
    {
        long n = 3_000_000_000;
        Assert.That(LcpWidth.Select(n, n / 16, 0), Is.EqualTo(LcpKind.U8));
        Assert.That(LcpWidth.Select(n, (n / 16) + 1, n / 8), Is.EqualTo(LcpKind.U16));
        Assert.That(LcpWidth.Select(n, n, (n / 8) + 1), Is.EqualTo(LcpKind.U32));
        Assert.That(LcpWidth.Select(4_294_967_292, 268_435_455, 0), Is.EqualTo(LcpKind.U8), "16 x (2^28 - 1) = 4,294,967,280 <= n");
        Assert.That(LcpWidth.Select(4_294_967_292, 268_435_456, 0), Is.EqualTo(LcpKind.U16), "16 x 2^28 = 2^32 > n, and no value reaches 65,535");
    }

    [Test]
    public void Select_RefusesMoreLongValuesThanShortOnes()
        => Assert.Throws<ArgumentException>(() => LcpWidth.Select(100, 3, 4));

    [Test]
    public void Bytes_Escape_OverflowCount_PerKind()
    {
        Assert.That(LcpWidth.Bytes(LcpKind.U8), Is.EqualTo(1));
        Assert.That(LcpWidth.Bytes(LcpKind.U16), Is.EqualTo(2));
        Assert.That(LcpWidth.Bytes(LcpKind.U32), Is.EqualTo(4));
        Assert.That(LcpWidth.Escape(LcpKind.U8), Is.EqualTo(255u));
        Assert.That(LcpWidth.Escape(LcpKind.U16), Is.EqualTo(65535u));
        Assert.That(LcpWidth.Escape(LcpKind.U32), Is.EqualTo(uint.MaxValue));
        Assert.That(LcpWidth.OverflowCount(LcpKind.U8, 11, 5), Is.EqualTo(11));
        Assert.That(LcpWidth.OverflowCount(LcpKind.U16, 11, 5), Is.EqualTo(5));
        Assert.That(LcpWidth.OverflowCount(LcpKind.U32, 11, 5), Is.Zero);
    }

    // ---- LcpGeometry

    private static long CeilDiv(long a, long b) => (a + b - 1) / b;

    [Test]
    public void Geometry_ChunkCounts_AndLengths_CoverTheTextAndTheRanksExactly()
    {
        foreach (long n in new long[] { 0, 1, 2, 3, 7, 8, 9, 100, 101, 1000 })
        {
            foreach (long r in new long[] { 1, 2, 3, 4, 7, 8, 9, 100, 101, 1000, 1001, 5000 })
            {
                var g = new LcpGeometry(n, r);
                Assert.That((g.N, g.R), Is.EqualTo((n, r)));
                Assert.That(g.TextChunks, Is.EqualTo(CeilDiv(n, r)), $"n={n} R={r}");
                Assert.That(g.RankChunks, Is.EqualTo(CeilDiv(n + 1, r)), $"n={n} R={r}");
                long text = 0, ranks = 0;
                for (long c = 0; c < g.TextChunks; c++)
                {
                    Assert.That(g.TextChunkLength(c), Is.InRange(1, r));
                    text += g.TextChunkLength(c);
                }
                for (long c = 0; c < g.RankChunks; c++)
                {
                    Assert.That(g.RankChunkLength(c), Is.InRange(1, r));
                    ranks += g.RankChunkLength(c);
                }
                Assert.That(text, Is.EqualTo(n), $"n={n} R={r}: text positions");
                Assert.That(ranks, Is.EqualTo(n + 1), $"n={n} R={r}: ranks");
            }
        }
    }

    [Test]
    public void Geometry_ChunkOf_IsTheDivisionByR_AtTheEdges()
    {
        var rng = new SplitMix64(11);
        foreach (long r in new long[] { 1, 2, 3, 7, 1000, 8_000_000, 266_667, 4_294_967_295, 4_294_967_296, 10_000_000_000 })
        {
            var g = new LcpGeometry(100, r);
            var xs = new List<uint> { 0, 1, 2, uint.MaxValue, uint.MaxValue - 1, 4_294_967_292 };
            foreach (long edge in new[] { r - 1, r, r + 1, (2 * r) - 1, 2 * r, (2 * r) + 1, 375 * r })
            {
                if (edge >= 0 && edge <= uint.MaxValue) xs.Add((uint)edge);
            }
            for (int i = 0; i < 300; i++) xs.Add((uint)rng.Next());
            foreach (uint x in xs) Assert.That((long)g.ChunkOf(x), Is.EqualTo((long)x / r), $"R={r} x={x}");
        }
    }

    [Test]
    public void Geometry_RefusesWhatItCannotRepresent()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => new LcpGeometry(-1, 8));
        Assert.Throws<ArgumentOutOfRangeException>(() => new LcpGeometry(10, 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => new LcpGeometry(10, -5));
        Assert.Throws<NotSupportedException>(() => new LcpGeometry(PackedTextFormat.MaxSymbols + 1, 8_000_000));
        var edge = new LcpGeometry(PackedTextFormat.MaxSymbols, 8_000_000);
        Assert.That(edge.TextChunks, Is.EqualTo(CeilDiv(PackedTextFormat.MaxSymbols, 8_000_000)));
        Assert.That(edge.RankChunks, Is.EqualTo(CeilDiv(PackedTextFormat.MaxSymbols + 1, 8_000_000)));
        Assert.That(edge.TextChunkLength(edge.TextChunks - 1), Is.EqualTo(PackedTextFormat.MaxSymbols - ((edge.TextChunks - 1) * 8_000_000)));
    }

    // ---- the parameter rule

    private static LcpChildOptions O(int threads = 4, long scaleN = 0, long ruleN = 0, LcpChildTuning? tuning = null, long arena = 160L << 20) => new()
    {
        ScratchDir = "s",
        Threads = threads,
        ScaleN = scaleN,
        RuleN = ruleN,
        ArenaCapBytes = arena,
        Tuning = tuning ?? new LcpChildTuning(),
    };

    [Test]
    public void Rule_ReferenceSizes_At3e9_AreTheOnesOfTheAdr()
    {
        LcpChildParameters p = LcpChildRules.Resolve(3_000_000_000, O());
        Assert.That(p.LcpChunk, Is.EqualTo(8_000_000));
        Assert.That((3_000_000_000 + p.LcpChunk - 1) / p.LcpChunk, Is.EqualTo(375), "P = 375");
        Assert.That(p.LcpWorkers, Is.EqualTo(2), "two B2 tasks fit under 160 MiB at the reference geometry");
        Assert.That(p.LcpReadBytes, Is.EqualTo(1 << 20));
        Assert.That(p.LcpTupleBytes, Is.EqualTo(32 * 1024));
        Assert.That(p.WindowBytes, Is.EqualTo(1 << 20));
        Assert.That(p.BlockEntries, Is.EqualTo(256 * 1024));
        Assert.That(p.ChildWindowBytes, Is.EqualTo(1 << 20));
        Assert.That(p.ChildSlots, Is.EqualTo(8));
        Assert.That(p.SpillFlushBytes, Is.EqualTo(16L << 20));
        Assert.That(p.Scaled, Is.False);
        Assert.That(p.Factor, Is.EqualTo(1.0));
        Assert.That(p.N, Is.EqualTo(3_000_000_000));
    }

    [Test]
    public void Rule_Unscaled_UsesTheReferenceChunk_ForEveryN()
    {
        foreach (long n in new long[] { 0, 1, 1000, 100_000_000, 3_000_000_000 })
            Assert.That(LcpChildRules.Resolve(n, O()).LcpChunk, Is.EqualTo(8_000_000), $"n={n}");
    }

    [Test]
    public void Rule_Scaled_ScalesTheChunk_ByNOverScaleN_WithAFloorOfThree()
    {
        LcpChildParameters p = LcpChildRules.Resolve(100_000_000, O(scaleN: 3_000_000_000));
        Assert.That(p.LcpChunk, Is.EqualTo(266_667), "ceil(8,000,000 x 10^8 / 3x10^9)");
        Assert.That((100_000_000 + p.LcpChunk - 1) / p.LcpChunk, Is.EqualTo(375), "P stays 375");
        Assert.That(p.Scaled, Is.True);
        Assert.That(p.Factor, Is.EqualTo(100_000_000 / 3e9).Within(1e-15));
        Assert.That(p.LcpWorkers, Is.EqualTo(2), "the workers come from the reference geometry, not from the scaled chunk");
        Assert.That(LcpChildRules.Resolve(1000, O(scaleN: 3_000_000_000)).LcpChunk, Is.EqualTo(3), "floor of 3");
        Assert.That(LcpChildRules.Resolve(1_500_000_000, O(scaleN: 3_000_000_000)).LcpChunk, Is.EqualTo(4_000_000));
    }

    [Test]
    public void Rule_RuleN_FixesTheNAtWhichTheRuleIsEvaluated()
    {
        LcpChildParameters p = LcpChildRules.Resolve(500_000_000, O(scaleN: 3_000_000_000, ruleN: 100_000_000));
        Assert.That(p.LcpChunk, Is.EqualTo(266_667));
        Assert.That(p.Factor, Is.EqualTo(100_000_000 / 3e9).Within(1e-15));
        Assert.That(p.N, Is.EqualTo(500_000_000));
        Assert.That(LcpChildRules.Resolve(500_000_000, O(ruleN: 100_000_000)).LcpChunk, Is.EqualTo(8_000_000), "RuleN is ignored without ScaleN");
    }

    [Test]
    public void Rule_Workers_FollowTheArenaCapAndTheThreads_AtTheReferenceGeometry()
    {
        long perTask = ExternalLcp.WorkerBytes(LcpChildRules.ReferenceN, new ExtLcpOptions { ScratchDir = "s", ChunkPositions = LcpChildRules.ChunkReference });
        Assert.That(perTask, Is.InRange(70L << 20, 80L << 20), "about 77 MB per task: two Phi / rank arrays of 8 M entries, the tuple buffers of 376 chunks");
        Assert.That(LcpChildRules.Resolve(1000, O(arena: 160L << 20)).LcpWorkers, Is.EqualTo(2));
        Assert.That(LcpChildRules.Resolve(1000, O(arena: 400L << 20)).LcpWorkers, Is.EqualTo(4), "capped by the threads");
        Assert.That(LcpChildRules.Resolve(1000, O(threads: 3, arena: 400L << 20)).LcpWorkers, Is.EqualTo(3));
        Assert.That(LcpChildRules.Resolve(1000, O(threads: 1, arena: 400L << 20)).LcpWorkers, Is.EqualTo(1));
        Assert.That(LcpChildRules.Resolve(1000, O(arena: 1)).LcpWorkers, Is.EqualTo(1), "at least one");
        Assert.That(LcpChildRules.Resolve(1000, O(arena: perTask)).LcpWorkers, Is.EqualTo(1));
        Assert.That(LcpChildRules.Resolve(1000, O(arena: (2 * perTask) - 1)).LcpWorkers, Is.EqualTo(1));
        Assert.That(LcpChildRules.Resolve(1000, O(arena: 2 * perTask)).LcpWorkers, Is.EqualTo(2));
    }

    [Test]
    public void Rule_ExplicitValuesWin_EveryOne()
    {
        var t = new LcpChildTuning
        {
            LcpChunk = 5,
            LcpWorkers = 3,
            LcpReadBytes = 17,
            LcpTupleBytes = 19,
            WindowBytes = 23,
            BlockEntries = 7,
            ChildWindowBytes = 64,
            ChildSlots = 5,
            SpillFlushBytes = 99,
        };
        LcpChildParameters p = LcpChildRules.Resolve(100_000_000, O(scaleN: 3_000_000_000, tuning: t));
        Assert.That((p.LcpChunk, p.LcpWorkers, p.LcpReadBytes, p.LcpTupleBytes, p.WindowBytes), Is.EqualTo((5L, 3, 17, 19, 23)));
        Assert.That((p.BlockEntries, p.ChildWindowBytes, p.ChildSlots, p.SpillFlushBytes), Is.EqualTo((7, 64, 5, 99L)));
    }

    [Test]
    public void Rule_RefusesNegativeArguments()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => LcpChildRules.Resolve(-1, O()));
        Assert.Throws<ArgumentOutOfRangeException>(() => LcpChildRules.Resolve(10, O(scaleN: -1)));
        Assert.Throws<ArgumentOutOfRangeException>(() => LcpChildRules.Resolve(10, O(ruleN: -1)));
        Assert.Throws<ArgumentOutOfRangeException>(() => LcpChildRules.Resolve(10, O(threads: 0)));
    }
}

/// <summary>The PLCP loop on one chunk of Phi (the loop of <c>LcpKernel.ComputePlcp</c> over a text of the packed layout).</summary>
[TestFixture]
public class PlcpScanTests
{
    /// <summary>Phi in text order from the suffix array (n for the smallest nonempty suffix, whose predecessor is the empty suffix).</summary>
    private static uint[] Phi(LcpChildCase c)
    {
        var phi = new uint[c.N];
        for (int rank = 1; rank <= c.N; rank++) phi[c.Sa[rank]] = (uint)c.Sa[rank - 1];
        return phi;
    }

    private static int Naive(byte[] s, int i, int j)
    {
        int k = 0;
        while (i + k < s.Length && j + k < s.Length && s[i + k] == s[j + k]) k++;
        return k;
    }

    private static void RunChunks(LcpChildCase c, int chunk, bool vectors, out uint[] plcp, out PlcpChunkResult total)
    {
        uint[] phi = Phi(c);
        plcp = new uint[c.N];
        long e8 = 0, e16 = 0, steps = 0, max = 0, first = 0;
        using var text = new PinnedText(c.Symbols, c.Bits);
        text.Reader.UseVectors = vectors;
        using var arena = new BudgetedArena(new ArenaBudget(long.MaxValue), ArenaOptions.Guarded(0xCD));
        for (int start = 0; start < c.N; start += chunk)
        {
            int count = Math.Min(chunk, c.N - start);
            Buf<uint> buf = arena.Allocate<uint>(count, zero: false, "phi");
            phi.AsSpan(start, count).CopyTo(buf.AsSpan(count));
            PlcpChunkResult r = PlcpScan.Run(text.Reader, start, buf, count);
            buf.AsSpan(count).CopyTo(plcp.AsSpan(start, count));
            e8 += r.Escapes8;
            e16 += r.Escapes16;
            steps += r.CompareSteps;
            max = Math.Max(max, r.MaxPlcp);
            first += r.FirstValue;
            arena.VerifyCanaries();
        }
        total = new PlcpChunkResult(e8, e16, steps, max, first);
    }

    private static void Check(string text, int chunk, bool vectors)
    {
        LcpChildCase c = LcpChildCase.Of(text);
        RunChunks(c, chunk, vectors, out uint[] plcp, out PlcpChunkResult total);
        uint[] phi = Phi(c);
        long e8 = 0, e16 = 0, max = 0, first = 0;
        for (int j = 0; j < c.N; j++)
        {
            int expected = phi[j] == c.N ? 0 : Naive(c.Symbols, j, (int)phi[j]);
            if (plcp[j] != expected) Assert.Fail($"'{(text.Length > 40 ? text[..40] + "..." : text)}' chunk {chunk} vectors {vectors}: PLCP[{j}] = {plcp[j]}, expected {expected}");
            if (expected >= 255) e8++;
            if (expected >= 65535) e16++;
            max = Math.Max(max, expected);
            if (j % chunk == 0) first += expected;
        }
        Assert.That((total.Escapes8, total.Escapes16, total.MaxPlcp, total.FirstValue), Is.EqualTo((e8, e16, max, first)), "counts, maximum and the restart cost");
        if (c.N > 0) Assert.That(total.CompareSteps, Is.GreaterThan(0));
    }

    [Test]
    public void Plcp_EqualsTheDefinition_ForEveryChunkSize_VectorAndScalar()
    {
        foreach (bool vectors in new[] { true, false })
        {
            foreach (string text in new[] { "a", "ab", "banana", "mississippi", "aaaaaaaaaa", "abababababab", StringFamilies.Fibonacci(200), StringFamilies.ThueMorse(150), StringFamilies.RandomDna(300, 3), "abcabcabcabcxabcabcabc" })
            {
                foreach (int chunk in new[] { 1, 2, 3, 5, 8, 13, 64, 1000 }) Check(text, chunk, vectors);
            }
        }
    }

    [Test]
    public void Plcp_CountsTheEscapes_At255And65535_Exactly()
    {
        // a^300: PLCP[j] = 299 - j for j >= 1 (the predecessor of the suffix at j is the shorter suffix at j + 1), so the values >= 255 are the positions 0..44 -> check by the definition
        foreach (string text in new[] { new string('a', 300), StringFamilies.RandomDna(1500, 4) + new string('c', 400) })
        {
            Check(text, 7, true);
            Check(text, 100000, false);
        }
        LcpChildCase c = LcpChildCase.Of(new string('a', 70_000));
        RunChunks(c, 1 << 20, true, out _, out PlcpChunkResult r);
        Assert.That(r.Escapes8, Is.EqualTo(70_000 - 255), "values 255 .. 69,999");
        Assert.That(r.Escapes16, Is.EqualTo(70_000 - 65_535), "values 65,535 .. 69,999");
        Assert.That(r.MaxPlcp, Is.EqualTo(69_999));
    }

    [Test]
    public void Plcp_ARepeatedBlock_LongerThanTheVectorStep_IsComparedAcrossSteps()
    {
        string block = StringFamilies.RandomDna(900, 5);
        string text = StringFamilies.RandomDna(2000, 6) + block + StringFamilies.RandomDna(333, 7) + block;
        foreach (bool vectors in new[] { true, false }) Check(text, 4096, vectors);
        Check(text, 500, true);
    }

    [Test]
    public void ARestartedChunk_HasTheSameValues_AsOneLongRun_ButCostsTheFirstValue()
    {
        LcpChildCase c = LcpChildCase.Of(new string('a', 500));
        RunChunks(c, 1000, true, out uint[] one, out PlcpChunkResult whole);
        RunChunks(c, 100, true, out uint[] five, out PlcpChunkResult chunked);
        Assert.That(five, Is.EqualTo(one));
        Assert.That(whole.FirstValue, Is.EqualTo(one[0]));
        Assert.That(chunked.FirstValue, Is.EqualTo(one[0] + one[100] + one[200] + one[300] + one[400]));
    }

    [Test]
    public void EmptyChunk_DoesNothing()
    {
        using var text = new PinnedText([0, 1, 0], 1);
        using var arena = new BudgetedArena(new ArenaBudget(long.MaxValue), ArenaOptions.Guarded(0xCD));
        Buf<uint> buf = arena.Allocate<uint>(1, zero: true, "phi");
        PlcpChunkResult r = PlcpScan.Run(text.Reader, 0, buf, 0);
        Assert.That(r, Is.EqualTo(new PlcpChunkResult(0, 0, 0, 0, 0)));
    }

    [Test]
    public void ACorruptPhi_NeverReadsOutsideTheText_AndStaysBelowN()
    {
        // a Phi that is garbage (entries beyond n, entries equal to the position): the loop is memory safe and every value stays within the text
        var rng = new SplitMix64(31);
        for (int round = 0; round < 50; round++)
        {
            int n = 1 + rng.NextInt(200);
            string t = StringFamilies.Random(n, 1 + rng.NextInt(3), (ulong)round + 1);
            LcpChildCase c = LcpChildCase.Of(t);
            using var text = new PinnedText(c.Symbols, c.Bits);
            using var arena = new BudgetedArena(new ArenaBudget(long.MaxValue), ArenaOptions.Guarded(0xCD));
            Buf<uint> buf = arena.Allocate<uint>(n, zero: false, "phi");
            for (int k = 0; k < n; k++) buf[k] = (uint)rng.Next();
            for (int k = 0; k < n; k += 3) buf[k] = (uint)rng.NextInt(n + 2);
            PlcpScan.Run(text.Reader, 0, buf, n);
            for (int k = 0; k < n; k++) Assert.That(buf[k], Is.LessThanOrEqualTo((uint)n), $"round {round} PLCP[{k}]");
            arena.VerifyCanaries();
        }
    }
}

/// <summary>The write side of B3: the base array in the chosen width and the overflow list.</summary>
[TestFixture]
public class LcpOutputWriterTests
{
    private TempDir _dir = null!;

    [SetUp]
    public void SetUp() => _dir = new TempDir();

    [TearDown]
    public void TearDown() => _dir.Dispose();

    private static BudgetedArena Arena() => new(new ArenaBudget(long.MaxValue), ArenaOptions.Guarded(0xCD));

    private (byte[] Base, byte[] Overflow) Write(LcpKind kind, uint[] dense, long overflowCount, int slice, int window, int step)
    {
        string bp = _dir.File($"b{kind}{slice}{window}{step}"), op = _dir.File($"o{kind}{slice}{window}{step}");
        long n = dense.Length - 1;
        using (MappedOutputFile bf = MappedOutputFile.Create(bp, (n + 1) * LcpChildRun.Width(kind)))
        using (MappedOutputFile of = MappedOutputFile.Create(op, overflowCount * 8))
        {
            using BudgetedArena arena = Arena();
            using var w = new LcpOutputWriter(arena, kind, n, overflowCount, bf.CreateWriter(window), of.CreateWriter(window), slice);
            for (int at = 0; at < dense.Length; at += step) w.Write(at, dense.AsSpan(at, Math.Min(step, dense.Length - at)));
            w.Complete(n);
        }
        return (File.ReadAllBytes(bp), File.ReadAllBytes(op));
    }

    private static uint[] Pairs(byte[] bytes) => System.Runtime.InteropServices.MemoryMarshal.Cast<byte, uint>(bytes).ToArray();

    [Test]
    public void U8_EscapesFrom255_AreAPairEach_AndTheBaseHoldsThe255()
    {
        uint[] dense = [0, 0, 1, 254, 255, 256, 7, 1_000_000, 255, 0];
        foreach (int slice in new[] { 1, 2, 3, 100 })
        {
            foreach (int step in new[] { 1, 4, 100 })
            {
                (byte[] b, byte[] o) = Write(LcpKind.U8, dense, 4, slice, 8 + slice, step);
                Assert.That(b, Is.EqualTo(new byte[] { 0, 0, 1, 254, 255, 255, 7, 255, 255, 0 }), $"slice {slice} step {step}");
                Assert.That(Pairs(o), Is.EqualTo(new uint[] { 4, 255, 5, 256, 7, 1_000_000, 8, 255 }), $"slice {slice} step {step}");
            }
        }
    }

    [Test]
    public void U16_EscapesFrom65535()
    {
        uint[] dense = [0, 0, 65_534, 65_535, 65_536, 1, 4_000_000_000, 300];
        (byte[] b, byte[] o) = Write(LcpKind.U16, dense, 3, 3, 16, 2);
        var expected = new ushort[] { 0, 0, 65_534, 65_535, 65_535, 1, 65_535, 300 };
        Assert.That(System.Runtime.InteropServices.MemoryMarshal.Cast<byte, ushort>(b).ToArray(), Is.EqualTo(expected));
        Assert.That(Pairs(o), Is.EqualTo(new uint[] { 3, 65_535, 4, 65_536, 6, 4_000_000_000 }));
    }

    [Test]
    public void U32_IsTheRawArray_AndNeverEscapes()
    {
        uint[] dense = [0, 0, 255, 65_535, uint.MaxValue - 1, 12];
        (byte[] b, byte[] o) = Write(LcpKind.U32, dense, 0, 2, 12, 3);
        Assert.That(Pairs(b), Is.EqualTo(dense));
        Assert.That(o, Is.Empty);
    }

    [Test]
    public void RanksMustContinue_WhereThePreviousCallEnded()
    {
        using var f = MappedOutputFile.Create(_dir.File("b"), 4);
        using var g = MappedOutputFile.Create(_dir.File("o"), 0);
        using BudgetedArena arena = Arena();
        using var w = new LcpOutputWriter(arena, LcpKind.U8, 3, 0, f.CreateWriter(), g.CreateWriter(), 4);
        w.Write(0, [0, 0]);
        Assert.Throws<InvalidOperationException>(() => w.Write(3, [1]), "a gap");
        Assert.Throws<InvalidOperationException>(() => w.Write(1, [0]), "a repeat");
        w.Write(2, [5, 6]);
        w.Complete(3);
    }

    [Test]
    public void Complete_NeedsEveryEntry_AndExactlyTheAnnouncedEscapes()
    {
        using (var f = MappedOutputFile.Create(_dir.File("b1"), 4))
        using (var g = MappedOutputFile.Create(_dir.File("o1"), 0))
        {
            using BudgetedArena arena = Arena();
            using var w = new LcpOutputWriter(arena, LcpKind.U8, 3, 0, f.CreateWriter(), g.CreateWriter(), 4);
            w.Write(0, [0, 0, 1]);
            Assert.Throws<InvalidOperationException>(() => w.Complete(3), "one entry short");
        }
        using (var f = MappedOutputFile.Create(_dir.File("b2"), 4))
        using (var g = MappedOutputFile.Create(_dir.File("o2"), 8))
        {
            using BudgetedArena arena = Arena();
            using var w = new LcpOutputWriter(arena, LcpKind.U8, 3, 1, f.CreateWriter(), g.CreateWriter(), 4);
            w.Write(0, [0, 0, 1, 2]);
            Assert.Throws<InvalidOperationException>(() => w.Complete(3), "an announced escape never came");
        }
    }

    [Test]
    public void MoreEscapesThanAnnounced_AreRefused()
    {
        using var f = MappedOutputFile.Create(_dir.File("b"), 6);
        using var g = MappedOutputFile.Create(_dir.File("o"), 8);
        using BudgetedArena arena = Arena();
        using var w = new LcpOutputWriter(arena, LcpKind.U8, 5, 1, f.CreateWriter(), g.CreateWriter(), 8);
        w.Write(0, [0, 0, 300]);
        Assert.Throws<InvalidOperationException>(() => w.Write(3, [400, 1, 2]));
    }

    [Test]
    public void TheRanges_MustHaveTheAnnouncedLengths()
    {
        using var f = MappedOutputFile.Create(_dir.File("b"), 100);
        using var g = MappedOutputFile.Create(_dir.File("o"), 100);
        using BudgetedArena arena = Arena();
        Assert.Throws<ArgumentException>(() => new LcpOutputWriter(arena, LcpKind.U8, 3, 0, f.CreateWriter(0, 5, 4), g.CreateWriter(0, 0, 4), 4), "base too long");
        Assert.Throws<ArgumentException>(() => new LcpOutputWriter(arena, LcpKind.U16, 3, 0, f.CreateWriter(0, 4, 4), g.CreateWriter(0, 0, 4), 4), "base too short for the width");
        Assert.Throws<ArgumentException>(() => new LcpOutputWriter(arena, LcpKind.U8, 3, 2, f.CreateWriter(0, 4, 4), g.CreateWriter(0, 8, 4), 4), "overflow too short");
        Assert.Throws<ArgumentException>(() => new LcpOutputWriter(arena, LcpKind.U8, 3, 0, f.CreateWriter(0, 4, 4), g.CreateWriter(0, 8, 4), 4), "overflow too long");
    }

    [Test]
    public void ARangeInsideALargerFile_LeavesTheRestUntouched()
    {
        string path = _dir.File("shared");
        using (MappedOutputFile f = MappedOutputFile.Create(path, 64))
        {
            using (MappedWindowWriter fill = f.CreateWriter(16)) fill.Append(Enumerable.Repeat((byte)0xEE, 64).ToArray());
            using BudgetedArena arena = Arena();
            using var g = MappedOutputFile.Create(_dir.File("ovf"), 8);
            using var w = new LcpOutputWriter(arena, LcpKind.U8, 4, 1, f.CreateWriter(10, 5, 4), g.CreateWriter(), 2);
            w.Write(0, [0, 0, 3, 300, 1]);
            w.Complete(4);
        }
        byte[] all = File.ReadAllBytes(path);
        Assert.That(all[..10], Is.All.EqualTo((byte)0xEE));
        Assert.That(all[10..15], Is.EqualTo(new byte[] { 0, 0, 3, 255, 1 }));
        Assert.That(all[15..], Is.All.EqualTo((byte)0xEE));
    }
}
