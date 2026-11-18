using System.Runtime.InteropServices;
using SuffixTree.Esa.External;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.External;

/// <summary>Records what the B3 hook delivers.</summary>
internal sealed class RecordingConsumer : ILcpChunkConsumer
{
    internal List<(long First, uint[] Values)> Calls { get; } = [];

    internal Action<long>? OnCall { get; set; }

    public void Consume(long firstRank, ReadOnlySpan<uint> lcp)
    {
        OnCall?.Invoke(firstRank);
        Calls.Add((firstRank, lcp.ToArray()));
    }
}

/// <summary>An <see cref="IReservationSource"/> that reserves normally for the first <c>succeed</c> creations and then fails like a full disk.</summary>
internal sealed class FailingAfter(int succeed) : IReservationSource
{
    private int _calls;

    public bool CanPreallocate => true;

    internal int Calls => _calls;

    public FileStream CreatePreallocated(string path, long length)
    {
        if (Interlocked.Increment(ref _calls) > succeed) throw new IOException($"no space left on device: {path}");
        return new FileStream(path, FileMode.Create, FileAccess.ReadWrite, FileShare.ReadWrite);
    }

    public long AvailableFreeBytes(string directory) => -1;
}

/// <summary>PH4 (base mode) on a given suffix-array file against the compact LCP of the product, and the B3 hook.</summary>
[TestFixture]
public class ExternalLcpTests
{
    private TempDir _dir = null!;

    [SetUp]
    public void SetUp() => _dir = new TempDir();

    [TearDown]
    public void TearDown() => _dir.Dispose();

    private string Scratch => _dir.Sub("scratch");

    private ExtLcpOptions Options(long r, int threads = 1, int workers = 0) => new()
    {
        ScratchDir = Scratch,
        ChunkPositions = r,
        Threads = threads,
        Workers = workers,
        ReadBufferBytes = 64,
        TupleBufferBytes = 48,
        WindowBytes = 128,
        ArenaGuards = ArenaOptions.Guarded(0xCD),
    };

    private (ExtLcpStats Stats, string Base, string Ovf) Run(LcpChildCase c, ExtLcpOptions o, ILcpChunkConsumer? consumer = null, string tag = "")
    {
        (string sa, string textPath) = LcpChildRun.WriteInputs(_dir, c, tag);
        string b = _dir.File($"lcp{tag}.base"), v = _dir.File($"lcp{tag}.ovf");
        using MappedTextReader text = MappedTextReader.Open(textPath, c.N, c.Bits);
        ExtLcpStats st = ExternalLcp.Build(sa, text, b, v, o, consumer);
        return (st, b, v);
    }

    private void Check(LcpChildCase c, ExtLcpStats st, string b, string v, string label)
    {
        Assert.That(st.Kind, Is.EqualTo(c.Lcp.Kind), $"{label}: kind");
        Assert.That(File.ReadAllBytes(b).AsSpan().SequenceEqual(c.LcpBaseBytes()), Is.True, $"{label}: base array ({c.Lcp.Kind})");
        Assert.That(File.ReadAllBytes(v).AsSpan().SequenceEqual(c.OverflowBytes()), Is.True, $"{label}: overflow list");
        Assert.That(st.N, Is.EqualTo(c.N), $"{label}: N");
        Assert.That(st.OverflowCount, Is.EqualTo(c.Lcp.OverflowCount), $"{label}: overflow count");
        int[] dense = c.Lcp.ToArray();
        Assert.That(st.Escapes8, Is.EqualTo(dense.Count(x => x >= 255)), $"{label}: values >= 255");
        Assert.That(st.Escapes16, Is.EqualTo(dense.Count(x => x >= 65535)), $"{label}: values >= 65,535");
        Assert.That(st.MaxPlcp, Is.EqualTo(c.MaxLcp), $"{label}: max PLCP");
        Assert.That(Directory.Exists(Scratch) ? Directory.GetFileSystemEntries(Scratch) : [], Is.Empty, $"{label}: scratch files left behind");
    }

    [Test]
    public void SmallTexts_EqualTheCompactLcp_ForEveryChunkSize_AndEveryNumberOfThreads()
    {
        foreach (string text in new[] { "a", "ab", "aa", "banana", "mississippi", "aaaaaaaaaaaaaaaa", "abababababababab", StringFamilies.Fibonacci(100), StringFamilies.ThueMorse(90), StringFamilies.RandomDna(120, 9), "abcabcabcabcabcxabc" })
        {
            LcpChildCase c = LcpChildCase.Of(text);
            foreach (long r in new long[] { 1, 2, 3, 5, 8, 13, 1000 })
            {
                foreach (int threads in new[] { 1, 2, 4 })
                {
                    (ExtLcpStats st, string b, string v) = Run(c, Options(r, threads, workers: threads));
                    Check(c, st, b, v, $"'{(text.Length > 30 ? text[..30] : text)}' R={r} threads {threads}");
                    Assert.That(st.TextChunks, Is.EqualTo((c.N + r - 1) / r));
                    Assert.That(st.RankChunks, Is.EqualTo((c.N + 1 + r - 1) / r));
                }
            }
        }
    }

    [Test]
    public void TheEmptyText_AndATextOfOne_Work()
    {
        foreach (string text in new[] { string.Empty, "x" })
        {
            LcpChildCase c = LcpChildCase.Of(text);
            (ExtLcpStats st, string b, string v) = Run(c, Options(3));
            Check(c, st, b, v, $"'{text}'");
            Assert.That(st.Kind, Is.EqualTo(LcpKind.U8));
            Assert.That(new FileInfo(b).Length, Is.EqualTo(c.N + 1));
            Assert.That(new FileInfo(v).Length, Is.Zero);
        }
    }

    [Test]
    public void EveryWidth_IsWrittenAsTheKernelWritesIt()
    {
        string u8 = LcpChildRun.RepeatRich(40_000, 3) + new string('G', 10); // few long values: U8 (+ overflow maybe)
        string u8Overflow = StringFamilies.RandomDna(20_000, 5) + StringFamilies.RandomDna(600, 11) + StringFamilies.RandomDna(100, 6) + StringFamilies.RandomDna(600, 11);
        string u16 = StringFamilies.RandomDna(20_000, 7) + StringFamilies.RandomDna(70_000, 8) + StringFamilies.RandomDna(333, 9) + StringFamilies.RandomDna(70_000, 8);
        string u32 = new string('a', 100_000); // 34,465 values at or above 65,535: more than n / 8
        var kinds = new Dictionary<LcpKind, int>();
        foreach ((string name, string text) in new[] { ("random", StringFamilies.RandomDna(30_000, 1)), ("repeat-rich", u8), ("u8 + overflow", u8Overflow), ("u16 + overflow", u16), ("u32", u32) })
        {
            LcpChildCase c = LcpChildCase.Of(text);
            kinds[c.Lcp.Kind] = kinds.GetValueOrDefault(c.Lcp.Kind) + 1;
            foreach (long r in new long[] { 7_919, 8_000_000 })
            {
                (ExtLcpStats st, string b, string v) = Run(c, Options(r, 2, 2) with { ReadBufferBytes = 4096, TupleBufferBytes = 1024, WindowBytes = 4096 });
                Check(c, st, b, v, $"{name} R={r}");
            }
            if (name == "u8 + overflow") Assert.That((c.Lcp.Kind, c.Lcp.OverflowCount > 0), Is.EqualTo((LcpKind.U8, true)));
            if (name == "u16 + overflow") Assert.That((c.Lcp.Kind, c.Lcp.OverflowCount > 0), Is.EqualTo((LcpKind.U16, true)));
            if (name == "u32") Assert.That(c.Lcp.Kind, Is.EqualTo(LcpKind.U32));
        }
        Assert.That(kinds.Keys, Is.EquivalentTo(new[] { LcpKind.U8, LcpKind.U16, LcpKind.U32 }), "all three widths are met");
    }

    [Test]
    public void TheHook_GetsEveryRankOnce_InOrder_InBlocksOfR_WithResolvedValues()
    {
        LcpChildCase c = LcpChildCase.Of(StringFamilies.RandomDna(20_000, 5) + StringFamilies.RandomDna(600, 11) + StringFamilies.RandomDna(100, 6) + StringFamilies.RandomDna(600, 11));
        Assert.That(c.Lcp.OverflowCount, Is.GreaterThan(0));
        var consumer = new RecordingConsumer();
        (ExtLcpStats st, string b, string v) = Run(c, Options(1000, 2, 2) with { ReadBufferBytes = 4096, TupleBufferBytes = 1024, WindowBytes = 4096 }, consumer);
        Check(c, st, b, v, "hook");
        long next = 0;
        uint[] expected = c.DenseLcp();
        foreach ((long first, uint[] values) in consumer.Calls)
        {
            Assert.That(first, Is.EqualTo(next), "ascending, no gap, no repeat");
            Assert.That(values.Length, Is.EqualTo(Math.Min(1000, expected.Length - next)), "one call per rank chunk");
            Assert.That(values.AsSpan().SequenceEqual(expected.AsSpan((int)first, values.Length)), Is.True, $"values of the block at {first}");
            next += values.Length;
        }
        Assert.That(next, Is.EqualTo(c.N + 1));
        Assert.That(consumer.Calls, Has.Count.EqualTo(st.RankChunks));
    }

    [Test]
    public void TheOutputs_AreCreatedAfterTheWidthIsKnown_AndOnlyThen()
    {
        LcpChildCase c = LcpChildCase.Of(StringFamilies.RandomDna(400, 5));
        string b = _dir.File("w.base"), v = _dir.File("w.ovf");
        var steps = new List<string>();
        bool existedBefore = false, existedAtFirstBlock = false;
        long lengthAtFirstBlock = -1;
        var consumer = new RecordingConsumer
        {
            OnCall = first =>
            {
                if (first != 0) return;
                existedAtFirstBlock = File.Exists(b) && File.Exists(v);
                lengthAtFirstBlock = new FileInfo(b).Length;
            },
        };
        ExtLcpOptions o = Options(50) with { Trace = s => { steps.Add(s); if (s == "b3") existedBefore = File.Exists(b) || File.Exists(v); } };
        (string sa, string textPath) = LcpChildRun.WriteInputs(_dir, c);
        using MappedTextReader text = MappedTextReader.Open(textPath, c.N, c.Bits);
        ExternalLcp.Build(sa, text, b, v, o, consumer);
        Assert.That(steps, Is.EqualTo(new[] { "b1", "b2", "b3" }));
        Assert.That(existedBefore, Is.False, "nothing exists when B3 begins");
        Assert.That(existedAtFirstBlock, Is.True, "both files exist when the first block is delivered");
        Assert.That(lengthAtFirstBlock, Is.EqualTo(c.N + 1), "the base file has its final length (U8) at the first block");
    }

    [Test]
    public void TheOutputIs_TheSame_ForEveryThreadCountWorkerCountAndChunkSize()
    {
        foreach (string text in new[] { StringFamilies.RandomDna(900, 3), new string('a', 500), StringFamilies.Fibonacci(700) })
        {
            LcpChildCase c = LcpChildCase.Of(text);
            string? reference = null;
            foreach (long r in new long[] { 3, 10, 64, 1_000_000 })
            {
                foreach ((int threads, int workers) in new[] { (1, 1), (2, 2), (4, 4), (4, 0), (3, 3) })
                {
                    (ExtLcpStats st, string b, string v) = Run(c, Options(r, threads, workers));
                    Check(c, st, b, v, $"R={r} threads {threads} workers {workers}");
                    string h = Convert.ToHexString(System.Security.Cryptography.SHA256.HashData(File.ReadAllBytes(b))) + Convert.ToHexString(System.Security.Cryptography.SHA256.HashData(File.ReadAllBytes(v)));
                    reference ??= h;
                    Assert.That(h, Is.EqualTo(reference), $"R={r} threads {threads} workers {workers}");
                }
            }
        }
    }

    [Test]
    public void TheStats_DoNotDependOnTheThreadCount_ButTheWorkerCountReportsIt()
    {
        LcpChildCase c = LcpChildCase.Of(StringFamilies.RandomDna(2000, 3));
        (ExtLcpStats one, _, _) = Run(c, Options(100, 1, 1), tag: "1");
        (ExtLcpStats four, _, _) = Run(c, Options(100, 4, 4), tag: "4");
        Assert.That((four.CompareSteps, four.RestartSymbols, four.MaxPlcp, four.Escapes8, four.Kind), Is.EqualTo((one.CompareSteps, one.RestartSymbols, one.MaxPlcp, one.Escapes8, one.Kind)));
        Assert.That((one.Workers, four.Workers), Is.EqualTo((1, 4)));
        Assert.That(one.CompareSteps, Is.GreaterThan(0));
        Assert.That(one.B1Seconds + one.B2Seconds + one.B3Seconds, Is.GreaterThanOrEqualTo(0));
    }

    [Test]
    public void TheNumberOfWorkers_IsTheSmallest_OfTheRequest_TheChunks_AndWhatTheBudgetHolds()
    {
        LcpChildCase c = LcpChildCase.Of(StringFamilies.RandomDna(600, 3));
        ExtLcpOptions o = Options(100, threads: 4);
        long per = ExternalLcp.WorkerBytes(c.N, o);
        Assert.That(Run(c, o with { Workers = 3 }, tag: "a").Stats.Workers, Is.EqualTo(3), "explicit");
        Assert.That(Run(c, o with { Workers = 0 }, tag: "b").Stats.Workers, Is.EqualTo(4), "threads, no budget");
        Assert.That(Run(c, o with { Workers = 0, Budget = new ArenaBudget((long)(2.5 * per)) }, tag: "c").Stats.Workers, Is.EqualTo(2), "what the budget holds");
        Assert.That(Run(c, o with { Workers = 0, Budget = new ArenaBudget(per) }, tag: "d").Stats.Workers, Is.EqualTo(1));
        Assert.That(Run(c, o with { Workers = 9 }, tag: "e").Stats.Workers, Is.EqualTo(6), "never more than the six text chunks");
        Assert.That(Run(c, o with { ChunkPositions = 1000, Workers = 4 }, tag: "f").Stats.Workers, Is.EqualTo(1), "one chunk, one worker");
    }

    [Test]
    public void WorkerBytes_IsWhatOneB2TaskHolds_AndTheBudgetPeakIsAtMostThat_Times_Workers()
    {
        LcpChildCase c = LcpChildCase.Of(StringFamilies.RandomDna(2000, 3));
        foreach (long r in new long[] { 50, 500, 100_000 })
        {
            var budget = new ArenaBudget(long.MaxValue);
            ExtLcpOptions o = Options(r, 1, 1) with { Budget = budget, ReadBufferBytes = 1024, TupleBufferBytes = 64, WindowBytes = 256, ArenaGuards = null };
            ExtLcpStats st = Run(c, o, tag: $"r{r}").Stats;
            long per = ExternalLcp.WorkerBytes(c.N, o);
            Assert.That(budget.LiveBytes, Is.Zero, $"R={r}: every arena is released");
            Assert.That(budget.PeakBytes, Is.EqualTo(per), $"R={r}: the peak is the B2 task, which holds exactly WorkerBytes (B1 and B3 are smaller)");
            Assert.That(st.Workers, Is.EqualTo(1));
        }
    }

    [Test]
    public void ATooSmallBudget_FailsWithTheBudgetException_AndLeavesNothing()
    {
        LcpChildCase c = LcpChildCase.Of(StringFamilies.RandomDna(600, 3));
        (string sa, string textPath) = LcpChildRun.WriteInputs(_dir, c);
        using MappedTextReader text = MappedTextReader.Open(textPath, c.N, c.Bits);
        string b = _dir.File("t.base"), v = _dir.File("t.ovf");
        var o = Options(100) with { Budget = new ArenaBudget(100) };
        Assert.Throws<ArenaBudgetExceededException>(() => ExternalLcp.Build(sa, text, b, v, o));
        Assert.That(File.Exists(b) || File.Exists(v), Is.False);
        Assert.That(Directory.Exists(Scratch) ? Directory.GetFileSystemEntries(Scratch) : [], Is.Empty);
        Assert.That(o.Budget!.LiveBytes, Is.Zero);
    }

    [Test]
    public void TheExplicitCounters_AreExactlyTheTupleVolume()
    {
        foreach (string text in new[] { StringFamilies.RandomDna(1500, 3), new string('b', 300), "banana" })
        {
            LcpChildCase c = LcpChildCase.Of(text);
            var counters = new IoCounters();
            (ExtLcpStats st, string b, string v) = Run(c, Options(64, 2, 2) with { Counters = counters });
            long n = c.N;
            Assert.That(counters.ReadBytes, Is.EqualTo((4 * (n + 1)) + (12 * n) + (8 * n)), "the suffix array (header read included), the Phi tuples, the rank tuples");
            Assert.That(counters.WriteBytes, Is.EqualTo((12 * n) + (8 * n)), "the same tuples written");
            Assert.That(counters.MappedWriteBytes, Is.EqualTo(new FileInfo(b).Length + new FileInfo(v).Length));
            Assert.That(counters.ScratchLiveBytes, Is.Zero);
            Assert.That(counters.ScratchPeakBytes, Is.InRange(12 * n, 20 * n));
            Check(c, st, b, v, text);
        }
    }

    // ---- corrupt suffix arrays

    private (string Sa, string Text, MappedTextReader Reader, LcpChildCase Case) Corrupt(Func<uint[], uint[]> damage, string text = "mississippi", int tag = 0)
    {
        LcpChildCase c = LcpChildCase.Of(text);
        string textPath = TextOracle.WriteFile(_dir, $"ctext{tag}.pk", c.Symbols, c.Bits);
        string sa = _dir.File($"csa{tag}.bin");
        FileBits.WriteUInts(sa, damage(c.Sa.Select(x => (uint)x).ToArray()));
        return (sa, textPath, MappedTextReader.Open(textPath, c.N, c.Bits), c);
    }

    private void AssertRejected(Func<uint[], uint[]> damage, string label, long r = 3, string text = "mississippi", int tag = 0)
    {
        (string sa, _, MappedTextReader reader, _) = Corrupt(damage, text, tag);
        using (reader)
        {
            string b = _dir.File($"x{tag}.base"), v = _dir.File($"x{tag}.ovf");
            var consumer = new RecordingConsumer();
            Assert.Throws<InvalidDataException>(() => ExternalLcp.Build(sa, reader, b, v, Options(r, 2, 2), consumer), label);
            Assert.That(File.Exists(b) || File.Exists(v), Is.False, $"{label}: no output is left");
            Assert.That(Directory.Exists(Scratch) ? Directory.GetFileSystemEntries(Scratch) : [], Is.Empty, $"{label}: no scratch is left");
        }
    }

    [Test]
    [CancelAfter(120_000)]
    public void ACorruptSuffixArray_EndsWithAnException_NeverAHang_AnyChunkSize()
    {
        int tag = 0;
        foreach (long r in new long[] { 1, 3, 4, 1000 })
        {
            AssertRejected(sa => { sa[0] = 0; return sa; }, $"R={r}: SA[0] != n", r, tag: tag++);
            AssertRejected(sa => { sa[0] = (uint)sa.Length; return sa; }, $"R={r}: SA[0] = n + 1", r, tag: tag++);
            AssertRejected(sa => { sa[5] = (uint)sa.Length - 1; return sa; }, $"R={r}: an entry equal to n", r, tag: tag++);
            AssertRejected(sa => { sa[5] = uint.MaxValue; return sa; }, $"R={r}: an entry of 2^32 - 1", r, tag: tag++);
            AssertRejected(sa => { sa[6] = sa[5]; return sa; }, $"R={r}: a repeated entry", r, tag: tag++);
            AssertRejected(sa => { sa[2] = sa[1]; return sa; }, $"R={r}: a repeated entry next to the empty suffix", r, tag: tag++);
            AssertRejected(sa => { sa[^1] = sa[^2]; return sa; }, $"R={r}: the last entry repeated", r, tag: tag++);
            AssertRejected(sa => sa[..^1], $"R={r}: one entry short", r, tag: tag++);
            AssertRejected(sa => [.. sa, 3u], $"R={r}: one entry too many", r, tag: tag++);
            AssertRejected(_ => [], $"R={r}: an empty file", r, tag: tag++);
            AssertRejected(sa => sa[..3], $"R={r}: a third of the file", r, tag: tag++);
            AssertRejected(sa => { Array.Fill(sa, 0u); sa[0] = (uint)sa.Length - 1; return sa; }, $"R={r}: all zero", r, tag: tag++);
        }
    }

    [Test]
    public void ARepeatedEntry_InsideOneChunk_IsCaught_EvenWhenEveryChunkHasTheRightCount()
    {
        // R covers the whole text: the chunk count is right (n tuples), only the permutation check can see the repeat
        AssertRejected(sa => { sa[4] = sa[3]; return sa; }, "one chunk", 1000, tag: 1);
        // R = 6 on n = 11: chunks of positions 0-5 and 6-10; the entry of position 2 replaced by 3: both lie in chunk 0, which still holds six tuples
        AssertRejected(sa =>
        {
            sa[Array.IndexOf(sa, 2u)] = 3;
            return sa;
        }, "a repeat inside chunk 0 with R = 6", 6, tag: 2);
    }

    [Test]
    public void ATruncatedFile_ThatLooksLikeASmallerText_IsRefusedByLength()
    {
        LcpChildCase c = LcpChildCase.Of("mississippi");
        (string sa, string textPath) = LcpChildRun.WriteInputs(_dir, c);
        using MappedTextReader text = MappedTextReader.Open(textPath, c.N, c.Bits);
        using (var f = new FileStream(sa, FileMode.Open)) f.SetLength(4 * c.N);
        InvalidDataException e = Assert.Throws<InvalidDataException>(() => ExternalLcp.Build(sa, text, _dir.File("a"), _dir.File("b"), Options(3)))!;
        Assert.That(e.Message, Does.Contain("12").And.Contain("44"), "the message names the expected entries and the bytes found");
        Assert.Throws<FileNotFoundException>(() => ExternalLcp.Build(_dir.File("missing.sa"), text, _dir.File("a"), _dir.File("b"), Options(3)));
        Assert.That(File.Exists(_dir.File("a")), Is.False);
    }

    [Test]
    public void ARejectedSuffixArray_TouchesNoExistingOutput()
    {
        LcpChildCase c = LcpChildCase.Of("mississippi");
        (string sa, string textPath) = LcpChildRun.WriteInputs(_dir, c);
        using MappedTextReader text = MappedTextReader.Open(textPath, c.N, c.Bits);
        using (var f = new FileStream(sa, FileMode.Open)) f.SetLength(8);
        string b = _dir.File("stale.base"), v = _dir.File("stale.ovf");
        File.WriteAllText(b, "stale-base");
        File.WriteAllText(v, "stale-ovf");
        Assert.Throws<InvalidDataException>(() => ExternalLcp.Build(sa, text, b, v, Options(3)));
        Assert.That(File.ReadAllText(b), Is.EqualTo("stale-base"), "the header check comes before any output is created");
        Assert.That(File.ReadAllText(v), Is.EqualTo("stale-ovf"));
    }

    [Test]
    public void ANonSortedPermutation_ProducesAnArray_OrAnException_NeverACrash()
    {
        // a permutation that is not the suffix array: the LCP of it is meaningless but every read stays inside the text and the values inside [0, n]
        var rng = new SplitMix64(77);
        for (int round = 0; round < 40; round++)
        {
            int n = 1 + rng.NextInt(150);
            string t = StringFamilies.Random(n, 1 + rng.NextInt(3), (ulong)round + 9);
            LcpChildCase c = LcpChildCase.Of(t);
            var sa = new uint[n + 1];
            sa[0] = (uint)n;
            var rest = Enumerable.Range(0, n).Select(i => (uint)i).ToArray();
            for (int i = n - 1; i > 0; i--)
            {
                int j = rng.NextInt(i + 1);
                (rest[i], rest[j]) = (rest[j], rest[i]);
            }
            Array.Copy(rest, 0, sa, 1, n);
            string saPath = _dir.File($"perm{round}.sa");
            FileBits.WriteUInts(saPath, sa);
            string textPath = TextOracle.WriteFile(_dir, $"perm{round}.pk", c.Symbols, c.Bits);
            using MappedTextReader text = MappedTextReader.Open(textPath, c.N, c.Bits);
            var consumer = new RecordingConsumer();
            ExtLcpStats st = ExternalLcp.Build(saPath, text, _dir.File($"p{round}.base"), _dir.File($"p{round}.ovf"), Options(1 + rng.NextInt(20), 2, 2), consumer);
            Assert.That(consumer.Calls.SelectMany(x => x.Values).All(x => x <= n), Is.True, $"round {round}: every value is within the text");
            Assert.That(st.Kind, Is.AnyOf(LcpKind.U8, LcpKind.U16, LcpKind.U32));
        }
    }

    // ---- failures of the environment

    [Test]
    public void AConsumerThatThrows_EndsTheBuild_AndRemovesTheOutputsAndTheScratch()
    {
        LcpChildCase c = LcpChildCase.Of(StringFamilies.RandomDna(300, 5));
        (string sa, string textPath) = LcpChildRun.WriteInputs(_dir, c);
        using MappedTextReader text = MappedTextReader.Open(textPath, c.N, c.Bits);
        string b = _dir.File("c.base"), v = _dir.File("c.ovf");
        foreach (int failAt in new[] { 0, 2 })
        {
            int calls = 0;
            var consumer = new RecordingConsumer { OnCall = _ => { if (calls++ == failAt) throw new IOException("the spill of the stack failed"); } };
            Assert.Throws<IOException>(() => ExternalLcp.Build(sa, text, b, v, Options(40, 2, 2), consumer), $"fail at call {failAt}");
            Assert.That(File.Exists(b) || File.Exists(v), Is.False, $"fail at call {failAt}: outputs");
            Assert.That(Directory.Exists(Scratch) ? Directory.GetFileSystemEntries(Scratch) : [], Is.Empty, $"fail at call {failAt}: scratch");
        }
    }

    [Test]
    public void AFullDisk_AtTheCreationOfAnOutput_IsAnIoException_AndNothingIsLeft()
    {
        LcpChildCase c = LcpChildCase.Of(StringFamilies.RandomDna(300, 5));
        (string sa, string textPath) = LcpChildRun.WriteInputs(_dir, c);
        using MappedTextReader text = MappedTextReader.Open(textPath, c.N, c.Bits);
        foreach (int succeed in new[] { 0, 1 })
        {
            string b = _dir.File($"d{succeed}.base"), v = _dir.File($"d{succeed}.ovf");
            var source = new FailingAfter(succeed);
            ExtLcpOptions o = Options(40, 2, 2) with { Files = new MappedFileOptions { Source = source } };
            Assert.Throws<IOException>(() => ExternalLcp.Build(sa, text, b, v, o), $"the creation number {succeed + 1} fails");
            Assert.That(source.Calls, Is.EqualTo(succeed + 1));
            Assert.That(File.Exists(b) || File.Exists(v), Is.False, $"{succeed} files were created before: both are removed");
            Assert.That(Directory.Exists(Scratch) ? Directory.GetFileSystemEntries(Scratch) : [], Is.Empty);
        }
    }

    [Test]
    public void ATraceThatThrows_FailsTheBuildAtThatStep_WithNothingLeft()
    {
        LcpChildCase c = LcpChildCase.Of(StringFamilies.RandomDna(300, 5));
        (string sa, string textPath) = LcpChildRun.WriteInputs(_dir, c);
        using MappedTextReader text = MappedTextReader.Open(textPath, c.N, c.Bits);
        foreach (string step in new[] { "b1", "b2", "b3" })
        {
            string b = _dir.File($"{step}.base"), v = _dir.File($"{step}.ovf");
            Assert.Throws<InvalidOperationException>(() => ExternalLcp.Build(sa, text, b, v, Options(40, 2, 2) with { Trace = s => { if (s == step) throw new InvalidOperationException(step); } }), step);
            Assert.That(File.Exists(b) || File.Exists(v), Is.False, step);
            Assert.That(Directory.Exists(Scratch) ? Directory.GetFileSystemEntries(Scratch) : [], Is.Empty, step);
        }
    }

    [Test]
    public void TheScalarCommonPrefix_GivesTheSameArray_AsTheVectorOne()
    {
        string text = StringFamilies.RandomDna(8000, 5) + StringFamilies.RandomDna(900, 12) + StringFamilies.RandomDna(77, 6) + StringFamilies.RandomDna(900, 12) + new string('A', 700);
        LcpChildCase c = LcpChildCase.Of(text);
        (string sa, string textPath) = LcpChildRun.WriteInputs(_dir, c);
        foreach (bool vectors in new[] { true, false })
        {
            using MappedTextReader reader = MappedTextReader.Open(textPath, c.N, c.Bits);
            reader.UseVectors = vectors;
            string b = _dir.File($"v{vectors}.base"), v = _dir.File($"v{vectors}.ovf");
            ExtLcpStats st = ExternalLcp.Build(sa, reader, b, v, Options(1500, 2, 2) with { ReadBufferBytes = 4096, TupleBufferBytes = 1024, WindowBytes = 4096 });
            Check(c, st, b, v, $"vectors {vectors}");
        }
    }

    [Test]
    public void TheArguments_AreChecked()
    {
        LcpChildCase c = LcpChildCase.Of("banana");
        (string sa, string textPath) = LcpChildRun.WriteInputs(_dir, c);
        using MappedTextReader text = MappedTextReader.Open(textPath, c.N, c.Bits);
        string b = _dir.File("a.base"), v = _dir.File("a.ovf");
        Assert.Throws<ArgumentNullException>(() => ExternalLcp.Build(sa, null!, b, v, Options(3)));
        Assert.Throws<ArgumentNullException>(() => ExternalLcp.Build(sa, text, b, v, null!));
        Assert.Throws<ArgumentOutOfRangeException>(() => ExternalLcp.Build(sa, text, b, v, Options(0)));
        Assert.Throws<ArgumentOutOfRangeException>(() => ExternalLcp.Build(sa, text, b, v, Options(3, threads: 0)));
        Assert.Throws<ArgumentOutOfRangeException>(() => ExternalLcp.Build(sa, text, b, v, Options(3, workers: -1)));
        Assert.Throws<ArgumentOutOfRangeException>(() => ExternalLcp.Build(sa, text, b, v, Options(3) with { ReadBufferBytes = 0 }));
        Assert.Throws<ArgumentOutOfRangeException>(() => ExternalLcp.Build(sa, text, b, v, Options(3) with { TupleBufferBytes = 0 }));
        Assert.Throws<ArgumentOutOfRangeException>(() => ExternalLcp.Build(sa, text, b, v, Options(3) with { WindowBytes = 0 }));
        Assert.That(File.Exists(b) || File.Exists(v), Is.False);
    }

    [Test]
    public void TheTuples_HaveTheDocumentedSizes()
    {
        Assert.That(Marshal.SizeOf<PhiTuple>(), Is.EqualTo(12));
        Assert.That(Marshal.SizeOf<RankTuple>(), Is.EqualTo(8));
    }
}
