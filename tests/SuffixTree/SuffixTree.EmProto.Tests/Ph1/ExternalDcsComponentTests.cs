using System.Buffers.Binary;
using System.Globalization;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;
using SuffixTree.Esa;
using SuffixTree.EmProto.Build;
using SuffixTree.EmProto.Dcs;
using SuffixTree.EmProto.Experiments;
using SuffixTree.EmProto.Harness;
using SuffixTree.EmProto.Tests.Build;
using SuffixTree.EmProto.Text;
using SuffixTree.Esa.Tests.TestSupport;
using SplitMix64 = SuffixTree.EmProto.Harness.SplitMix64;

namespace SuffixTree.EmProto.Tests.Ph1;

// Tests written AFTER the code of the external-memory PH1 (work package E3b): they pin the building blocks one by one so that the mutation audit of the new code has something to kill with besides the
// byte comparison with the in-RAM table (the oracle tests of ExternalDcsTests.cs were written first and committed red).

[TestFixture]
public class SampleGeometryTests
{
    [Test]
    public void SlotsAndReducedIndices_FollowTheDefinitions_ForManyPeriodsAndLengths()
    {
        foreach (int x in new[] { 3, 4, 7, 9, 13, 31 })
        {
            DifferenceCover cover = DifferenceCover.Create(x);
            for (int n = 0; n <= 80; n++)
            {
                var g = new SampleGeometry(n, cover);
                var sample = Enumerable.Range(0, n + 1).Where(p => cover.IndexOf[p % x] >= 0).ToList();
                Assert.That(g.Samples, Is.EqualTo(sample.Count), $"X={x} n={n}");
                Assert.That(g.Slots, Is.EqualTo(((n / x) + 1) * cover.Residues.Length));
                Assert.That(g.RedLength, Is.EqualTo(sample.Count + cover.Residues.Length));
                Assert.That(g.ClassCount.Sum(), Is.EqualTo(sample.Count));
                for (int i = 0; i < sample.Count; i++)
                {
                    Assert.That(g.SlotOf(sample[i]), Is.EqualTo(i), $"the real slots are 0..S-1 in position order (X={x} n={n})");
                    Assert.That(g.PositionOfSlot(i), Is.EqualTo(sample[i]));
                }

                // the reduced string by definition: class after class, positions r, r + X, ..., then a separator
                var red = new List<long>();
                for (int c = 0; c < cover.Residues.Length; c++)
                {
                    Assert.That(g.RedStart[c], Is.EqualTo(red.Count), $"X={x} n={n} class {c}");
                    for (long p = cover.Residues[c]; p <= n; p += x) red.Add(p);
                    red.Add(-1);
                }
                Assert.That(g.RedStart[cover.Residues.Length], Is.EqualTo(red.Count));
                for (int q = 0; q < red.Count; q++)
                {
                    bool real = g.TrySlotOfRed(q, out long slot);
                    Assert.That(real, Is.EqualTo(red[q] >= 0), $"X={x} n={n} q={q}");
                    if (real)
                    {
                        Assert.That(g.PositionOfSlot(slot), Is.EqualTo(red[q]));
                        Assert.That(g.RedIndexOfSlot(slot), Is.EqualTo(q));
                    }
                }
            }
        }
    }

    [Test]
    public void ClassesWithoutMembers_HaveOnlyTheirSeparator()
    {
        DifferenceCover cover = DifferenceCover.Create(31);
        int big = cover.Residues[^1];
        var g = new SampleGeometry(big - 1, cover);
        Assert.That(g.ClassCount[^1], Is.EqualTo(0));
        Assert.That(g.RedStart[^1] - g.RedStart[^2], Is.EqualTo(1), "an empty class is one separator");
        Assert.That(g.RedStart[^2] - g.RedStart[^3], Is.EqualTo(g.ClassCount[^3] + 1));
        Assert.Throws<ArgumentOutOfRangeException>(() => _ = new SampleGeometry(-1, cover));
    }
}

[TestFixture]
public class MappedWindowWriterTests
{
    private string _dir = string.Empty;

    [SetUp]
    public void SetUp()
    {
        _dir = Path.Combine(Path.GetTempPath(), "emproto-ph1c-" + Guid.NewGuid().ToString("N"));
        Directory.CreateDirectory(_dir);
    }

    [TearDown]
    public void TearDown() => Directory.Delete(_dir, true);

    [Test]
    public void SequentialAppends_AcrossManyWindows_ProduceTheFile_AndTheWindowsStayBounded()
    {
        string path = Path.Combine(_dir, "a.bin");
        var data = new byte[10_007];
        new SplitMix64(3).Next();
        var rng = new SplitMix64(9);
        for (int i = 0; i < data.Length; i++) data[i] = (byte)rng.Next();
        MappedIo.Reset();
        using (var w = new MappedWindowWriter(path, data.Length, 1000))
        {
            int at = 0;
            foreach (int chunk in new[] { 1, 999, 1000, 1001, 4000, 2500, 506 })
            {
                w.Append(data.AsSpan(at, chunk));
                at += chunk;
                Assert.That(w.Position, Is.EqualTo(at));
            }
            Assert.That(at, Is.EqualTo(data.Length));
            Assert.Throws<ArgumentOutOfRangeException>(() => w.Append(new byte[1]));
        }
        Assert.That(File.ReadAllBytes(path), Is.EqualTo(data));
        Assert.That(MappedIo.WriteBytes, Is.EqualTo(data.Length));
        Assert.That(MappedIo.MaxWindowBytes, Is.LessThanOrEqualTo(1000));
        Assert.That(MappedIo.Windows, Is.GreaterThanOrEqualTo(10), "10,007 bytes in windows of at most 1,000");
    }

    [Test]
    public void WriteAt_JumpsBackAndForth_AndTheUntouchedBytesStayZero()
    {
        string path = Path.Combine(_dir, "b.bin");
        using (var w = new MappedWindowWriter(path, 5000, 256))
        {
            w.WriteAt(4990, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]);
            w.WriteAt(0, [9, 9]);
            w.WriteAt(2500, new byte[600].Select((_, i) => (byte)(i % 250 + 1)).ToArray());
            Assert.Throws<ArgumentOutOfRangeException>(() => w.WriteAt(4995, new byte[6]));
            Assert.Throws<ArgumentOutOfRangeException>(() => w.WriteAt(-1, new byte[1]));
        }
        byte[] file = File.ReadAllBytes(path);
        Assert.That(file, Has.Length.EqualTo(5000));
        Assert.That(file[..2], Is.EqualTo(new byte[] { 9, 9 }));
        Assert.That(file[2..2500].All(b => b == 0), Is.True);
        Assert.That(file[2500], Is.EqualTo(1));
        Assert.That(file[3099], Is.EqualTo((byte)(599 % 250 + 1)));
        Assert.That(file[4990..], Is.EqualTo(new byte[] { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 }));
    }

    [Test]
    public void EmptyFile_NeedsNoMapping()
    {
        string path = Path.Combine(_dir, "c.bin");
        using (var w = new MappedWindowWriter(path, 0, 16))
        {
            w.Append([]);
            Assert.Throws<ArgumentOutOfRangeException>(() => w.Append(new byte[1]));
        }
        Assert.That(new FileInfo(path).Length, Is.EqualTo(0));
        Assert.Throws<ArgumentOutOfRangeException>(() => _ = new MappedWindowWriter(path, -1, 16));
        Assert.Throws<ArgumentOutOfRangeException>(() => _ = new MappedWindowWriter(path, 10, 0));
    }
}

[TestFixture]
public class BucketRunTests
{
    private string _dir = string.Empty;

    [SetUp]
    public void SetUp()
    {
        _dir = Path.Combine(Path.GetTempPath(), "emproto-ph1c-" + Guid.NewGuid().ToString("N"));
        Directory.CreateDirectory(_dir);
    }

    [TearDown]
    public void TearDown() => Directory.Delete(_dir, true);

    private static ulong[] Records(int count, int slotStep) => Enumerable.Range(0, count).Select(i => ((ulong)(i * slotStep) << 32) | (uint)(i + 1)).ToArray();

    [Test]
    public void SmallRunsStayInMemory_LargeRunsGoToAFile_AndTheFileIsTrackedAndDeleted()
    {
        ScratchTracker.Reset();
        IoCounters.Reset();
        using (BucketRun small = BucketRun.Create(Records(8, 1), _dir, "t", 64))
        {
            Assert.That(small.InMemory, Is.True);
            Assert.That(small.Count, Is.EqualTo(8));
            Assert.That(Directory.GetFiles(_dir), Is.Empty);
        }
        BucketRun big = BucketRun.Create(Records(9, 1), _dir, "t", 64);
        Assert.That(big.InMemory, Is.False, "72 bytes are above the spill size of 64");
        Assert.That(Directory.GetFiles(_dir), Has.Length.EqualTo(1));
        Assert.That(ScratchTracker.Live, Is.EqualTo(72));
        Assert.That(IoCounters.WriteBytes, Is.EqualTo(72));
        var dest = new ulong[4];
        Assert.That(big.ReadAt(2, dest), Is.EqualTo(4));
        Assert.That(dest, Is.EqualTo(Records(9, 1)[2..6]));
        Assert.That(big.ReadAt(7, dest), Is.EqualTo(2), "only two records are left");
        Assert.That(dest[..2], Is.EqualTo(Records(9, 1)[7..]));
        Assert.That(big.ReadAt(9, dest), Is.EqualTo(0));
        big.Dispose();
        Assert.That(Directory.GetFiles(_dir), Is.Empty);
        Assert.That(ScratchTracker.Live, Is.EqualTo(0));
        using BucketRun zero = BucketRun.Create(Records(1, 1), _dir, "t", 0);
        Assert.That(zero.InMemory, Is.False, "a spill size of 0 sends every run to a file (the scale-down emulation)");
    }

    [TestCase(1, 3, 100)]
    [TestCase(2, 5, 8)]
    [TestCase(100_000, 1, 16)]
    public unsafe void Drain_ScattersTheRecordsBelowTheLimit_ThroughABoundedBuffer(int spill, int slotStep, int bufferRecords)
    {
        const int count = 400;
        ulong[] records = Records(count, slotStep);
        using BucketRun run = BucketRun.Create(records, _dir, "d", spill);
        using var arena = new BudgetedArena(new ArenaBudget(1 << 20));
        Buf<ulong> buf = run.InMemory ? default : arena.Allocate<ulong>(bufferRecords, false, "t");
        var reader = new RunReader(run, buf, bufferRecords);
        int total = count * slotStep;
        var window = new uint[total + 10];
        long drained = 0;
        for (int start = 0; start < total; start += 37)
        {
            int end = Math.Min(total, start + 37);
            drained += reader.Drain((uint)end, 1000, window.AsSpan(), 0);
        }
        Assert.That(drained, Is.EqualTo(count));
        Assert.That(reader.Scattered, Is.EqualTo(count));
        for (int i = 0; i < count; i++) Assert.That(window[i * slotStep], Is.EqualTo(1000u + (uint)(i + 1)), $"record {i}");
        Assert.That(reader.Drain(uint.MaxValue, 0, window.AsSpan(), 0), Is.EqualTo(0));
    }

    [Test]
    public void Drain_StopsAtTheFirstRecordAtOrAboveTheLimit_AndRelativeToTheWindowStart()
    {
        ulong[] records = [(10UL << 32) | 1, (11UL << 32) | 2, (12UL << 32) | 3, (20UL << 32) | 4];
        using BucketRun run = BucketRun.Create(records, _dir, "e", 1000);
        var reader = new RunReader(run, default, 0);
        var window = new uint[11];
        Assert.That(reader.Drain(12, 100, window, 10), Is.EqualTo(2));
        Assert.That(window, Is.EqualTo(new uint[] { 101, 102, 0, 0, 0, 0, 0, 0, 0, 0, 0 }));
        Assert.That(reader.Drain(12, 100, window, 10), Is.EqualTo(0), "stopped at slot 12");
        Assert.That(reader.Drain(21, 100, window, 10), Is.EqualTo(2));
        Assert.That(window, Is.EqualTo(new uint[] { 101, 102, 103, 0, 0, 0, 0, 0, 0, 0, 104 }));
        Assert.That(reader.Scattered, Is.EqualTo(4));
    }
}

[TestFixture]
public class WorkPoolTests
{
    [Test]
    public void TasksPushedByWorkers_AreAllProcessed_AndTheLoopsTerminate()
    {
        var pool = new WorkPool<StrongBox<int>>();
        pool.Push(new StrongBox<int>(1));
        long sum = 0;
        var threads = Enumerable.Range(0, 4).Select(_ => new Thread(() =>
        {
            while (pool.Take(out StrongBox<int>? t))
            {
                try
                {
                    Interlocked.Add(ref sum, t.Value);
                    if (t.Value < 200)
                    {
                        pool.Push(new StrongBox<int>((t.Value * 2) + 1));
                        pool.Push(new StrongBox<int>((t.Value * 2) + 2));
                    }
                }
                finally
                {
                    pool.Done();
                }
            }
        })).ToList();
        foreach (Thread t in threads) t.Start();
        foreach (Thread t in threads) Assert.That(t.Join(TimeSpan.FromSeconds(30)), Is.True, "the pool must detect that all work is done");
        // the tree 1 -> 2k+1, 2k+2 below 200 has the nodes 1..(first level above 200)
        long expected = 0;
        var queue = new Queue<int>();
        queue.Enqueue(1);
        while (queue.Count > 0)
        {
            int v = queue.Dequeue();
            expected += v;
            if (v < 200)
            {
                queue.Enqueue((v * 2) + 1);
                queue.Enqueue((v * 2) + 2);
            }
        }
        Assert.That(sum, Is.EqualTo(expected));
        Assert.That(pool.Error, Is.Null);
    }

    [Test]
    public void AFailure_StopsEveryWorker_AndIsReported()
    {
        var pool = new WorkPool<StrongBox<int>>();
        for (int i = 0; i < 50; i++) pool.Push(new StrongBox<int>(i));
        pool.Fail(new InvalidOperationException("boom"));
        pool.Fail(new InvalidOperationException("second"));
        Assert.That(pool.Take(out _), Is.False);
        Assert.That(pool.Error!.Message, Is.EqualTo("boom"), "the first failure is kept");
    }

    [Test]
    public void LastInFirstOut_AndEmptyPoolEnds()
    {
        var pool = new WorkPool<StrongBox<int>>();
        pool.Push(new StrongBox<int>(1));
        pool.Push(new StrongBox<int>(2));
        Assert.That(pool.Take(out StrongBox<int>? a), Is.True);
        Assert.That(a!.Value, Is.EqualTo(2));
        pool.Done();
        Assert.That(pool.Take(out StrongBox<int>? b), Is.True);
        Assert.That(b!.Value, Is.EqualTo(1));
        pool.Done();
        Assert.That(pool.Take(out _), Is.False);
    }
}

[TestFixture]
public class NameSinkTests
{
    private string _dir = string.Empty;

    [SetUp]
    public void SetUp()
    {
        _dir = Path.Combine(Path.GetTempPath(), "emproto-ph1c-" + Guid.NewGuid().ToString("N"));
        Directory.CreateDirectory(_dir);
    }

    [TearDown]
    public void TearDown() => Directory.Delete(_dir, true);

    private static uint[] NamesInSlotOrder(SampleGeometry g) => Enumerable.Range(0, (int)g.Samples).Select(s => (uint)((s * 7) + 3)).ToArray();

    [Test]
    public void RankFileSink_WritesLittleEndianNamesInSlotOrder_AndZeroesTheUnusedSlots()
    {
        DifferenceCover cover = DifferenceCover.Create(7);
        var g = new SampleGeometry(40, cover);
        string path = Path.Combine(_dir, "r.bin");
        uint[] names = NamesInSlotOrder(g);
        ScratchTracker.Reset();
        using (var sink = new RankFileSink(path, g.Slots, 64, scratch: true))
        {
            sink.Accept(0, names.AsSpan(0, 10));
            sink.Accept(10, names.AsSpan(10));
        }
        Assert.That(ScratchTracker.Live, Is.EqualTo(g.Slots * 4));
        byte[] file = File.ReadAllBytes(path);
        Assert.That(file, Has.Length.EqualTo(g.Slots * 4));
        for (int s = 0; s < g.Slots; s++) Assert.That(BinaryPrimitives.ReadUInt32LittleEndian(file.AsSpan(s * 4)), Is.EqualTo(s < g.Samples ? names[s] : 0u));
        ScratchTracker.Reset();
        using (new RankFileSink(Path.Combine(_dir, "o.bin"), 5, 64, scratch: false)) Assert.That(ScratchTracker.Live, Is.EqualTo(0));
    }

    [Test]
    public void ReducedFileSink_WritesBigEndianClassMajorWithSeparatorsAndPadding()
    {
        foreach (int bufferBytes in new[] { 16, 40, 4096 })
        {
            DifferenceCover cover = DifferenceCover.Create(7);
            var g = new SampleGeometry(60, cover);
            string path = Path.Combine(_dir, $"red{bufferBytes}.bin");
            uint[] names = NamesInSlotOrder(g);
            ScratchTracker.Reset();
            using (var arena = new BudgetedArena(new ArenaBudget(1 << 20)))
            {
                var sink = new ReducedFileSink(path, g, bufferBytes, arena);
                int at = 0;
                foreach (int take in new[] { 1, 5, 7, 2 })
                {
                    sink.Accept(at, names.AsSpan(at, take));
                    at += take;
                }
                sink.Accept(at, names.AsSpan(at));
                sink.Complete();
                sink.Dispose();
            }
            Assert.That(ScratchTracker.Live, Is.EqualTo((g.RedLength * 4) + PackedText.PadBytes));
            byte[] file = File.ReadAllBytes(path);
            Assert.That(file, Has.Length.EqualTo((g.RedLength * 4) + PackedText.PadBytes));
            var expected = new uint[g.RedLength];
            for (int slot = 0; slot < g.Samples; slot++) expected[g.RedIndexOfSlot(slot)] = names[slot];
            for (int q = 0; q < g.RedLength; q++) Assert.That(BinaryPrimitives.ReadUInt32BigEndian(file.AsSpan(q * 4)), Is.EqualTo(expected[q]), $"buffer {bufferBytes} index {q}");
            Assert.That(file[^PackedText.PadBytes..].All(b => b == 0), Is.True);
            // the file is a PackedText of 32 bits per symbol
            unsafe
            {
                fixed (byte* p = file)
                {
                    using var t = new PackedText(p, g.RedLength, 32, null, false);
                    for (int q = 0; q < g.RedLength; q++) Assert.That(t.Symbol(q), Is.EqualTo(expected[q]));
                }
            }
        }
    }

    [Test]
    public void ReducedRamSink_PlacesTheNamesClassMajor()
    {
        DifferenceCover cover = DifferenceCover.Create(9);
        var g = new SampleGeometry(77, cover);
        uint[] names = NamesInSlotOrder(g);
        using var arena = new BudgetedArena(new ArenaBudget(1 << 20));
        using var sink = new ReducedRamSink(g, arena);
        sink.Accept(0, names.AsSpan(0, 13));
        sink.Accept(13, names.AsSpan(13));
        int[] red = sink.Red.AsSpan((int)g.RedLength).ToArray();
        var expected = new int[g.RedLength];
        for (int slot = 0; slot < g.Samples; slot++) expected[g.RedIndexOfSlot(slot)] = (int)names[slot];
        Assert.That(red, Is.EqualTo(expected));
    }
}

[TestFixture]
public class SaInverterTests
{
    private string _dir = string.Empty;

    [SetUp]
    public void SetUp()
    {
        _dir = Path.Combine(Path.GetTempPath(), "emproto-ph1c-" + Guid.NewGuid().ToString("N"));
        Directory.CreateDirectory(_dir);
    }

    [TearDown]
    public void TearDown() => Directory.Delete(_dir, true);

    [Test]
    public void Inversion_EqualsTheNaiveInverse_ForRandomSuffixArrays_AndManyChunkAndBufferSizes()
    {
        var rng = new SplitMix64(21);
        foreach (int x in new[] { 7, 13 })
        {
            DifferenceCover cover = DifferenceCover.Create(x);
            foreach (long n in new long[] { 0, 1, 20, 93, 400 })
            {
                var g = new SampleGeometry(n, cover);
                int reduced = (int)g.RedLength;
                // a random permutation of 0..reduced-1 stands for the suffix array of the reduced string, SA[0] = reduced
                int[] perm = Enumerable.Range(0, reduced).OrderBy(_ => rng.Next()).ToArray();
                var sa = new uint[reduced + 1];
                sa[0] = (uint)reduced;
                for (int i = 0; i < reduced; i++) sa[i + 1] = (uint)perm[i];
                string saPath = Path.Combine(_dir, "sa.bin");
                File.WriteAllBytes(saPath, MemoryMarshal.AsBytes(sa.AsSpan()).ToArray());
                var expected = new uint[g.Slots];
                for (int k = 1; k <= reduced; k++)
                {
                    if (g.TrySlotOfRed(sa[k], out long slot)) expected[slot] = (uint)k;
                }
                foreach ((long chunk, int ioBytes) in new[] { (1L, 16), (5L, 24), (g.Slots + 3, 4096), (33L, 100) })
                {
                    string rank = Path.Combine(_dir, "rank.bin");
                    ScratchTracker.Reset();
                    ExtDcsOptions o = ExtDcsHarness.Options(Path.Combine(_dir, "scratch"), 3, 3) with { ChunkSlots = chunk, IoBufferBytes = ioBytes, WindowBytes = 64 };
                    Directory.CreateDirectory(o.ScratchDir);
                    SaInverter.Invert(saPath, reduced, g, rank, true, o, "inv", new ArenaBudget(1 << 24));
                    byte[] file = File.ReadAllBytes(rank);
                    Assert.That(MemoryMarshal.Cast<byte, uint>(file).ToArray(), Is.EqualTo(expected), $"X={x} n={n} chunk={chunk} io={ioBytes}");
                    Assert.That(Directory.GetFileSystemEntries(o.ScratchDir), Is.Empty, "the tuple files are deleted");
                    Assert.That(ScratchTracker.Live, Is.EqualTo(g.Slots * 4), "only the scratch output is still tracked");
                }
            }
        }
    }
}

[TestFixture]
public class PrefixOrderTests
{
    [Test]
    public void Compare_FollowsTheXPrefixOrderWithTheEndOfTheTextSmallest()
    {
        var rng = new SplitMix64(5);
        foreach (int bits in new[] { 1, 2, 8, 32 })
        {
            foreach (int x in new[] { 1, 3, 7, 20 })
            {
                var sym = new byte[60];
                for (int i = 0; i < sym.Length; i++) sym[i] = (byte)rng.NextInt(2);
                using PackedText t = PackedText.Pack(sym, bits);
                var order = new PrefixOrder(t, x);
                for (int a = 0; a <= sym.Length; a++)
                {
                    for (int b = 0; b <= sym.Length; b++)
                    {
                        int expected = Naive(sym, x, a, b);
                        Assert.That(Math.Sign(order.Compare(a, b)), Is.EqualTo(expected), $"b={bits} X={x} a={a} b={b}");
                    }
                }
            }
        }
    }

    /// <summary>The first X symbols of the suffix compared as strings over symbol + 1 with the end of the text as 0; equal only when both suffixes have X symbols that agree (or the same position).</summary>
    private static int Naive(byte[] s, int x, int a, int b)
    {
        if (a == b) return 0;
        int la = s.Length - a, lb = s.Length - b;
        for (int i = 0; i < x; i++)
        {
            int ca = i < la ? s[a + i] + 1 : 0;
            int cb = i < lb ? s[b + i] + 1 : 0;
            if (ca != cb) return ca < cb ? -1 : 1;
            if (ca == 0) return a > b ? -1 : 1; // both ended at the same offset: impossible for a != b (different lengths), kept for completeness
        }
        return 0;
    }
}

[TestFixture]
public class E3bCommandTests
{
    private string _dir = string.Empty;

    [SetUp]
    public void SetUp()
    {
        _dir = Path.Combine(Path.GetTempPath(), "emproto-e3b-" + Guid.NewGuid().ToString("N"));
        Directory.CreateDirectory(_dir);
    }

    [TearDown]
    public void TearDown() => Directory.Delete(_dir, true);

    [TestCase("genome:300000", "K_s=105 N_base=1678 R=419 K_sort=105", "run spill 7 bytes")]
    [TestCase("a:200000", "K_s=70 N_base=1118 R=280 K_sort=70", "run spill 4 bytes")]
    public void Pipeline_Pack_Mkranks_E3b_Verify_ReportsEquality(string family, string parameters, string spill)
    {
        if (!File.Exists("/proc/self/stat")) Assert.Ignore("the report reads /proc (Linux)");
        string pk = Path.Combine(_dir, "t.pk");
        string oracle = Path.Combine(_dir, "t.oracle");
        string rank = Path.Combine(_dir, "t.rank");
        var report = new Report(Path.Combine(_dir, "rows.tsv"));
        Assert.That(E3.Pack(Options.Parse(["--family", family, "--out", pk]), report), Is.EqualTo(0));
        Assert.That(E3.MakeRanks(Options.Parse(["--text", pk, "--out", oracle, "--period", "1024"]), report), Is.EqualTo(0));
        int rc = E3b.Build(Options.Parse(["--text", pk, "--out", rank, "--scratch", Path.Combine(_dir, "scratch"), "--scale-n", "3000000000", "--verify", oracle, "--threads", "3", "--rss-cap-mb", "100000", "--case", "t"]), report);
        Assert.That(rc, Is.EqualTo(0));
        string rows = File.ReadAllText(Path.Combine(_dir, "rows.tsv"));
        Assert.That(rows, Does.Contain("rank_equals_oracle\tyes"));
        Assert.That(rows, Does.Contain("\tparameters\t" + parameters), "K_s, N_base, R and K_sort are their 3e9 reference values times n / 3e9");
        Assert.That(rows, Does.Contain(spill));
        Assert.That(rows, Does.Contain("\tnonunique_share\t"));
        Assert.That(File.ReadAllBytes(rank), Is.EqualTo(File.ReadAllBytes(oracle)));
        Assert.That(Directory.GetFileSystemEntries(Path.Combine(_dir, "scratch")), Is.Empty);
    }

    [TestCase(1001L)]
    [TestCase(2_500_000L)]
    [TestCase(8L)]
    public void PackRandom_WritesAValidPackedFile_WithMetaAndPadding(long n)
    {
        string path = Path.Combine(_dir, "r.pk");
        Assert.That(E3b.PackRandom(Options.Parse(["--n", n.ToString(CultureInfo.InvariantCulture), "--out", path, "--seed", "5"]), new Report(null)), Is.EqualTo(0));
        byte[] file = File.ReadAllBytes(path);
        Assert.That(file, Has.Length.EqualTo(PackedText.TotalBytes(n, 2)));
        Assert.That(file[^PackedText.PadBytes..].All(b => b == 0), Is.True);
        (long metaN, int bits) = E3.ReadMeta(path);
        Assert.That((metaN, bits), Is.EqualTo((n, 2)));
        unsafe
        {
            fixed (byte* p = file)
            {
                using var t = new PackedText(p, n, 2, null, false);
                var symbols = new byte[n];
                var counts = new long[4];
                for (long i = 0; i < n; i++)
                {
                    symbols[i] = (byte)t.Symbol(i);
                    counts[symbols[i]]++;
                }
                using PackedText again = PackedText.Pack(symbols, 2);
                Assert.That(new ReadOnlySpan<byte>(again.Pointer, (int)PackedText.TotalBytes(n, 2)).SequenceEqual(file), Is.True, "re-packing the symbols gives the same bytes (the unused bits of the last byte are zero)");
                if (n > 1000) Assert.That(counts.All(c => Math.Abs(c - (n / 4.0)) < n * 0.05), Is.True, "the four symbols are about equally frequent");
            }
        }
        Assert.Throws<ArgumentException>(() => E3b.PackRandom(Options.Parse(["--out", path]), new Report(null)));
    }

    [Test]
    public void TouchedArena_MakesBlocksResidentAndKeepsTheContent_UntouchedDoesNot()
    {
        if (!RssSampler.IsSupported) Assert.Ignore("needs RssAnon (Linux)");
        static long RssAnonKb()
        {
            foreach (string line in File.ReadLines("/proc/self/status"))
            {
                if (line.StartsWith("RssAnon:", StringComparison.Ordinal)) return long.Parse(line.Split(' ', StringSplitOptions.RemoveEmptyEntries)[1], CultureInfo.InvariantCulture);
            }
            return 0;
        }
        const int mb = 48;
        long before = RssAnonKb();
        using (var plain = new BudgetedArena(new ArenaBudget(1L << 30), null, touchPages: false))
        {
            _ = plain.Allocate<byte>(mb << 20, zero: true, "plain");
            Assert.That(RssAnonKb() - before, Is.LessThan(8 * 1024), "a zeroed block that is not touched is not resident");
        }
        using (var touched = new BudgetedArena(new ArenaBudget(1L << 30), null, touchPages: true))
        {
            Buf<byte> b = touched.Allocate<byte>(mb << 20, zero: true, "touched");
            Assert.That(RssAnonKb() - before, Is.GreaterThan((mb - 4) * 1024L), "every page of a touched block is resident");
            Assert.That(b.AsSpan(mb << 20).IndexOfAnyExcept((byte)0), Is.EqualTo(-1), "the content stays zero");
        }
    }

    [Test]
    public void CompareFiles_FindsTheFirstDifference_AndALengthMismatch()
    {
        string a = Path.Combine(_dir, "a"), b = Path.Combine(_dir, "b");
        var data = new byte[3_000_000];
        new SplitMix64(1).Next();
        var rng = new SplitMix64(2);
        for (int i = 0; i < data.Length; i++) data[i] = (byte)rng.Next();
        File.WriteAllBytes(a, data);
        File.WriteAllBytes(b, data);
        Assert.That(E3b.CompareFiles(a, b), Is.EqualTo((true, -1L)));
        data[2_500_001] ^= 1;
        File.WriteAllBytes(b, data);
        Assert.That(E3b.CompareFiles(a, b), Is.EqualTo((false, 2_500_001L)));
        File.WriteAllBytes(b, data[..2_000_000]);
        File.WriteAllBytes(a, data[..2_000_000]);
        Assert.That(E3b.CompareFiles(a, b), Is.EqualTo((true, -1L)));
        File.WriteAllBytes(b, data[..1_999_999]);
        Assert.That(E3b.CompareFiles(a, b), Is.EqualTo((false, 1_999_999L)));
    }

    [Test]
    public void FailedBuild_LeavesNoScratchBehind()
    {
        DifferenceCover cover = DifferenceCover.Create(31);
        using PackedText t = PackedText.Pack(new byte[5000], 1);
        string scratch = Path.Combine(_dir, "scratch");
        ExtDcsOptions o = ExtDcsHarness.Options(scratch, 10, 20, 1) with { Budget = new ArenaBudget(2000) };
        Assert.Catch<InvalidOperationException>(() => ExternalDcs.BuildRankFile(t, cover, Path.Combine(_dir, "r"), o));
        Assert.That(Directory.GetFileSystemEntries(scratch), Is.Empty);
    }
}

// Tests added after round 1 of the mutation audit (analysis/mutation/results_e3b_*.tsv): they close the survivors that were not equivalent.

[TestFixture]
public class SampleNamerAuditTests
{
    private string _dir = string.Empty;

    [SetUp]
    public void SetUp()
    {
        _dir = Path.Combine(Path.GetTempPath(), "emproto-ph1a-" + Guid.NewGuid().ToString("N"));
        Directory.CreateDirectory(_dir);
    }

    [TearDown]
    public void TearDown() => Directory.Delete(_dir, true);

    private static (NamingResult Result, ExtDcsLevelStats Stats, SampleNamer Namer, SampleGeometry Geometry, PackedText Text) Name(byte[] symbols, int bits, int x, ExtDcsOptions o)
    {
        Directory.CreateDirectory(o.ScratchDir);
        PackedText t = PackedText.Pack(symbols, bits);
        DifferenceCover cover = DifferenceCover.Create(x);
        var g = new SampleGeometry(symbols.Length, cover);
        var stats = new ExtDcsLevelStats();
        var namer = new SampleNamer(t, g, o, stats, new ArenaBudget(1 << 26), "audit");
        return (namer.Run(), stats, namer, g, t);
    }

    [Test]
    public void EverySplitMakesNonEmptyBuckets_AndEveryRunIsAFileWhenTheSpillIsZero()
    {
        var rng = new SplitMix64(31);
        for (int round = 0; round < 40; round++)
        {
            var symbols = new byte[600 + rng.NextInt(1500)];
            int family = round % 4;
            for (int i = 0; i < symbols.Length; i++) symbols[i] = family switch { 0 => (byte)rng.NextInt(2), 1 => 0, 2 => (byte)((i / 3) & 1), _ => (byte)(rng.NextInt(10) == 0 ? 1 : 0) };
            ExtDcsOptions o = ExtDcsHarness.Options(Path.Combine(_dir, "scratch"), 3 + rng.NextInt(30), 1 << 20, 1 + (round % 3)) with { ReservoirSize = 3 + rng.NextInt(20), FanOut = 2 + rng.NextInt(50), RunSpillBytes = 0, Seed = rng.Next() };
            (NamingResult r, ExtDcsLevelStats st, _, _, PackedText t) = Name(symbols, 1, 7 + rng.NextInt(20), o);
            using (t)
            using (r)
            {
                Assert.That(r.Leaves.Count, Is.EqualTo(st.Leaves), $"round {round}: a split produced an empty bucket");
                Assert.That(st.FileRuns, Is.EqualTo(st.Leaves), $"round {round}: with a spill size of 0 every bucket run is a file");
            }
        }
    }

    [Test]
    public void RootBetweenOneAndTwoTimesTheCapacity_IsSplitDirectly_WithoutAnOverflowingCollectScan()
    {
        var rng = new SplitMix64(5);
        var symbols = new byte[2000];
        for (int i = 0; i < symbols.Length; i++) symbols[i] = (byte)rng.NextInt(2);
        DifferenceCover cover = DifferenceCover.Create(31);
        long s = cover.SampleCount(symbols.Length);
        long k = (s * 2 / 3) + 1; // K < Samples <= 2 K
        ExtDcsOptions o = ExtDcsHarness.Options(Path.Combine(_dir, "scratch"), k, 1 << 20) with { ReservoirSize = 4096 };
        (NamingResult r, ExtDcsLevelStats st, _, _, PackedText t) = Name(symbols, 1, 31, o);
        using (t)
        using (r)
        {
            Assert.That(st.SplitScans, Is.EqualTo(1));
            Assert.That(st.Scans, Is.EqualTo(st.Leaves + st.SplitScans), "the root is known to be larger than K: one split scan, then one scan per leaf, nothing wasted");
        }
    }

    private sealed class CollectingSink : INameSink
    {
        internal int Calls;

        public void Accept(long firstSlot, ReadOnlySpan<uint> names) => Calls++;

        public void Dispose()
        {
        }
    }

    [Test]
    public void Merge_ReportsAMissingRecordInItsWindow_AndADuplicate()
    {
        var rng = new SplitMix64(9);
        var symbols = new byte[1500];
        for (int i = 0; i < symbols.Length; i++) symbols[i] = (byte)rng.NextInt(2);
        ExtDcsOptions o = ExtDcsHarness.Options(Path.Combine(_dir, "scratch"), 50, 1 << 20) with { RunSpillBytes = 1 << 20, ChunkSlots = 100 };
        foreach (bool duplicate in new[] { false, true })
        {
            (NamingResult r, _, SampleNamer namer, _, PackedText t) = Name(symbols, 1, 31, o);
            using (t)
            using (r)
            {
                LeafResult leaf = r.Leaves[1];
                ulong[] records = leaf.Run!.Memory.ToArray();
                ulong[] damaged = duplicate ? [.. records, records[^1]] : records[..^1];
                var broken = new LeafResult { Count = leaf.Count, Distinct = leaf.Distinct, FirstPos = leaf.FirstPos, LastPos = leaf.LastPos, FirstLen = leaf.FirstLen, LastLen = leaf.LastLen, Singles = leaf.Singles, Run = BucketRun.Create(damaged, _dir, "damaged", 1 << 20), Offset = leaf.Offset };
                r.Leaves[1] = broken;
                var ex = Assert.Throws<InvalidOperationException>(() => namer.Merge(r, new CollectingSink()))!;
                Assert.That(ex.Message, Does.Contain("merge window"), $"duplicate={duplicate}: the window check names the window with the wrong number of records ({ex.Message})");
            }
        }
    }

    [Test]
    public void MergeWindow_FollowsChunkSlots_NotTheBucketCapacity()
    {
        var rng = new SplitMix64(12);
        var symbols = new byte[5000];
        for (int i = 0; i < symbols.Length; i++) symbols[i] = (byte)rng.NextInt(4);
        DifferenceCover cover = DifferenceCover.Create(31);
        long s = cover.SampleCount(symbols.Length);
        var budget = new ArenaBudget(1 << 26);
        ExtDcsOptions o = ExtDcsHarness.Options(Path.Combine(_dir, "scratch"), 100, 1 << 20) with { RunSpillBytes = 1 << 20, ChunkSlots = 1000, Budget = budget };
        using PackedText t = PackedText.Pack(symbols, 2);
        byte[] expected = ExtDcsHarness.Oracle(t, cover);
        (byte[] got, _) = ExtDcsHarness.Run(t, cover, o, _dir);
        Assert.That(got.AsSpan().SequenceEqual(expected), Is.True);
        Assert.That(budget.PeakBytes, Is.GreaterThanOrEqualTo(s * 4), $"the merge window holds min(ChunkSlots, Samples) = {s} entries of 4 bytes; the bucket arrays of K = 100 are only {100 * 24} bytes");
    }

    [Test]
    public void ExternalSortLeaf_FollowsSortLeafCapacity_NotTheBucketCapacity()
    {
        DenseCase c = BuilderTestSupport.Dense("same", StringFamilies.AllSame(3000));
        DifferenceCover cover = DifferenceCover.Create(13);
        using PackedText t = PackedText.Pack(c.Symbols, 1);
        var budget = new ArenaBudget(1 << 26);
        ExtDcsOptions o = ExtDcsHarness.Options(Path.Combine(_dir, "scratch"), 3, 40, 1) with { SortLeafCapacity = 500, Budget = budget };
        (byte[] got, ExtDcsStats st) = ExtDcsHarness.Run(t, cover, o, _dir);
        Assert.That(st.MaxLevel, Is.GreaterThanOrEqualTo(1));
        Assert.That(got.AsSpan().SequenceEqual(ExtDcsHarness.Oracle(t, cover)), Is.True);
        Assert.That(budget.PeakBytes, Is.GreaterThanOrEqualTo(500 * 24), "the leaf sorter of the external level allocates 24 bytes for each of K_sort = 500 records");
    }

    [Test]
    public void E3bCommand_ReturnsFailure_WhenTheOracleDiffers()
    {
        if (!File.Exists("/proc/self/stat")) Assert.Ignore("the report reads /proc (Linux)");
        string pk = Path.Combine(_dir, "t.pk");
        string oracle = Path.Combine(_dir, "t.oracle");
        var report = new Report(Path.Combine(_dir, "rows.tsv"));
        Assert.That(E3.Pack(Options.Parse(["--family", "genome:100000", "--out", pk]), report), Is.EqualTo(0));
        Assert.That(E3.MakeRanks(Options.Parse(["--text", pk, "--out", oracle]), report), Is.EqualTo(0));
        byte[] bytes = File.ReadAllBytes(oracle);
        bytes[1000] ^= 1;
        File.WriteAllBytes(oracle, bytes);
        int rc = E3b.Build(Options.Parse(["--text", pk, "--out", Path.Combine(_dir, "t.rank"), "--scratch", Path.Combine(_dir, "scratch"), "--scale-n", "3000000000", "--verify", oracle, "--threads", "2", "--rss-cap-mb", "100000"]), report);
        Assert.That(rc, Is.EqualTo(1));
        Assert.That(File.ReadAllText(Path.Combine(_dir, "rows.tsv")), Does.Contain("rank_equals_oracle\tNO"));
    }
}
