using System.Runtime.InteropServices;
using SuffixTree.Esa;
using SuffixTree.EmProto.Build;
using SuffixTree.EmProto.Experiments;
using SuffixTree.EmProto.Harness;
using SuffixTree.EmProto.Ph4;
using SuffixTree.EmProto.Tests.Build;
using SuffixTree.EmProto.Text;
using SuffixTree.Esa.Tests.TestSupport;
using SplitMix64 = SuffixTree.EmProto.Harness.SplitMix64;

namespace SuffixTree.EmProto.Tests.Ph4;

// Tests written AFTER the code of PH4 (work package E4): they pin the building blocks one by one (chunk arithmetic, width rule, PLCP loop on a chunk, the spool, the output writer, the worker budget, the
// options, the commands) so that the mutation audit of the new code has something to kill with besides the oracle tests of ExternalLcpTests.cs, which were written first against a stub.

[TestFixture]
public class LcpGeometryTests
{
    [Test]
    public void ChunkOf_EqualsTheDivision_ForEveryChunkSizeAndEdgeValues()
    {
        var rng = new SplitMix64(3);
        var edges = new uint[] { 0, 1, 2, 3, 4, 5, 255, 256, 65535, 65536, int.MaxValue, (uint)int.MaxValue + 1, uint.MaxValue - 1, uint.MaxValue };
        var sizes = new List<long>(Enumerable.Range(1, 300).Select(i => (long)i)) { 1024, 4096, 65537, 1 << 20, 8_000_000, int.MaxValue, uint.MaxValue, 1L << 32, (1L << 32) + 5, 1L << 40 };
        foreach (long r in sizes)
        {
            var g = new LcpGeometry(10, r);
            foreach (uint x in edges) Assert.That(g.ChunkOf(x), Is.EqualTo((uint)(x / (ulong)r)), $"{x} / {r}");
            for (int i = 0; i < 30; i++)
            {
                uint x = (uint)rng.Next();
                Assert.That(g.ChunkOf(x), Is.EqualTo((uint)(x / (ulong)r)), $"{x} / {r}");
            }
            foreach (uint k in new uint[] { 1, 2, 3, 7 })
            {
                if (r * k > uint.MaxValue) continue;
                uint boundary = (uint)(r * k);
                Assert.That(g.ChunkOf(boundary - 1), Is.EqualTo(k - 1), $"{boundary - 1} / {r}");
                Assert.That(g.ChunkOf(boundary), Is.EqualTo(k), $"{boundary} / {r}");
            }
        }
    }

    [Test]
    public void ChunkCountsAndLengths_CoverTheTextAndTheRanksExactly()
    {
        for (long n = 0; n <= 40; n++)
        {
            for (long r = 1; r <= 45; r++)
            {
                var g = new LcpGeometry(n, r);
                Assert.That(g.TextChunks, Is.EqualTo((n + r - 1) / r), $"n={n} R={r}");
                Assert.That(g.RankChunks, Is.EqualTo((n + 1 + r - 1) / r), $"n={n} R={r}");
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
                Assert.That(text, Is.EqualTo(n), $"n={n} R={r}");
                Assert.That(ranks, Is.EqualTo(n + 1), $"n={n} R={r}");
                Assert.That(g.N, Is.EqualTo(n));
                Assert.That(g.R, Is.EqualTo(r));
            }
        }
        var big = new LcpGeometry(3_000_000_000, 8_000_000);
        Assert.That(big.TextChunks, Is.EqualTo(375), "ADR-004 section 4.2: P = n / R = 375 chunks at 3x10^9 symbols and R = 8 M");
    }

    [Test]
    public void InvalidArguments_AreRejected()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => _ = new LcpGeometry(10, 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => _ = new LcpGeometry(-1, 5));
        Assert.Throws<NotSupportedException>(() => _ = new LcpGeometry(uint.MaxValue, 5));
        Assert.DoesNotThrow(() => _ = new LcpGeometry(uint.MaxValue - 1, 1 << 20));
    }
}

[TestFixture]
public class LcpWidthTests
{
    [Test]
    public void Select_EqualsCompactLcpSelectKind_OnAGridAndAtTheBoundaries()
    {
        foreach (int n in new[] { 0, 1, 2, 3, 15, 16, 17, 31, 32, 33, 100, 255, 1000, 65535, 70000 })
        {
            var e8s = new SortedSet<long>(Enumerable.Range(0, 40).Select(i => (long)i)) { n / 16, (n / 16) + 1, n / 8, (n / 8) + 1, n };
            foreach (long e8 in e8s.Where(e => e <= n + 2))
            {
                var e16s = new SortedSet<long>(Enumerable.Range(0, 40).Select(i => (long)i)) { n / 8, (n / 8) + 1, e8 };
                foreach (long e16 in e16s.Where(e => e <= e8))
                    Assert.That(LcpWidth.Select(n, e8, e16), Is.EqualTo(CompactLcp.SelectKind(n, e8, e16)), $"n={n} e8={e8} e16={e16}");
            }
        }
        Assert.That(LcpWidth.Select(160, 10, 0), Is.EqualTo(LcpKind.U8), "16 * escapes8 == n: U8");
        Assert.That(LcpWidth.Select(160, 11, 0), Is.EqualTo(LcpKind.U16), "16 * escapes8 > n, no escapes16: U16");
        Assert.That(LcpWidth.Select(160, 30, 20), Is.EqualTo(LcpKind.U16), "8 * escapes16 == n: U16");
        Assert.That(LcpWidth.Select(160, 30, 21), Is.EqualTo(LcpKind.U32));
        Assert.That(LcpWidth.Select(4_000_000_000, 250_000_000, 0), Is.EqualTo(LcpKind.U8), "64-bit counts");
    }

    [Test]
    public void Select_RejectsMoreHighEscapesThanLowEscapes() => Assert.Throws<ArgumentException>(() => LcpWidth.Select(100, 1, 2));

    [Test]
    public void Bytes_Escape_OverflowCount_PerKind()
    {
        Assert.That(LcpWidth.Bytes(LcpKind.U8), Is.EqualTo(1));
        Assert.That(LcpWidth.Bytes(LcpKind.U16), Is.EqualTo(2));
        Assert.That(LcpWidth.Bytes(LcpKind.U32), Is.EqualTo(4));
        Assert.That(LcpWidth.Escape(LcpKind.U8), Is.EqualTo(255u));
        Assert.That(LcpWidth.Escape(LcpKind.U16), Is.EqualTo(65535u));
        Assert.That(LcpWidth.Escape(LcpKind.U32), Is.EqualTo(uint.MaxValue));
        Assert.That(LcpWidth.OverflowCount(LcpKind.U8, 7, 3), Is.EqualTo(7));
        Assert.That(LcpWidth.OverflowCount(LcpKind.U16, 7, 3), Is.EqualTo(3));
        Assert.That(LcpWidth.OverflowCount(LcpKind.U32, 7, 3), Is.Zero);
    }
}

[TestFixture]
public class PlcpScanTests
{
    private static unsafe (PlcpChunkResult Result, uint[] Plcp) Scan(PackedText t, uint[] phi, long first, int count)
    {
        using var arena = new NativeArena();
        Buf<uint> buf = arena.Allocate<uint>(count, zero: true, "phi");
        for (int i = 0; i < count; i++) buf[i] = phi[first + i];
        PlcpChunkResult r = PlcpScan.Run(t, first, buf, count);
        return (r, buf.AsSpan(count).ToArray());
    }

    private static uint[] PhiOf(uint[] sa)
    {
        int n = sa.Length - 1;
        var phi = new uint[n];
        for (int i = 1; i <= n; i++) phi[sa[i]] = sa[i - 1];
        return phi;
    }

    private static int NaiveLcp(byte[] s, int a, int b)
    {
        int l = 0;
        while (a + l < s.Length && b + l < s.Length && s[a + l] == s[b + l]) l++;
        return l;
    }

    [TestCase(1000, 4, 1, false)]
    [TestCase(1000, 4, 2, true)]
    [TestCase(1500, 2, 7, false)]
    [TestCase(777, 2, 3, true)]
    public void AnyChunking_GivesTheValuesByDefinition_AndTheCounts(int n, int sigma, int bitsExtra, bool scalar)
    {
        // a random text with a planted repeat of 700 symbols, so that values above 255 occur
        var rng = new SplitMix64((ulong)n);
        var s = new byte[n];
        for (int i = 0; i < n; i++) s[i] = (byte)rng.NextInt(sigma);
        if (n >= 1500) Array.Copy(s, 100, s, 800, 700);
        var c = new DenseCase("t", s, sigma);
        uint[] sa = c.ReferenceSa();
        uint[] phi = PhiOf(sa);
        using PackedText t = PackedText.Pack(s, c.MinBits + bitsExtra - 1);
        t.UseVectors = t.UseVectors && !scalar;
        var expected = new uint[n];
        for (int j = 0; j < n; j++) expected[j] = phi[j] == n ? 0u : (uint)NaiveLcp(s, j, (int)phi[j]);
        foreach (int chunk in new[] { 1, 2, 7, 64, 333, n })
        {
            var got = new uint[n];
            long e8 = 0, e16 = 0, firstSum = 0, max = 0;
            for (int first = 0; first < n; first += chunk)
            {
                int count = Math.Min(chunk, n - first);
                (PlcpChunkResult r, uint[] plcp) = Scan(t, phi, first, count);
                plcp.CopyTo(got, first);
                e8 += r.Escapes8;
                e16 += r.Escapes16;
                firstSum += r.FirstValue;
                max = Math.Max(max, r.MaxPlcp);
                Assert.That(r.FirstValue, Is.EqualTo(expected[first]), $"chunk {chunk} at {first}: the value at the first position is compared from scratch");
            }
            Assert.That(got, Is.EqualTo(expected), $"n={n} chunk={chunk}");
            Assert.That(e8, Is.EqualTo(expected.Count(v => v >= 255)));
            Assert.That(e16, Is.Zero);
            Assert.That(max, Is.EqualTo(expected.Max()));
            Assert.That(firstSum, Is.EqualTo(Enumerable.Range(0, (n + chunk - 1) / chunk).Sum(k => (long)expected[k * chunk])));
        }
        if (n >= 1500) Assert.That(expected.Count(v => v >= 255), Is.GreaterThan(100), "the planted repeat gives escapes");
    }

    [Test]
    public void AllSame_RestartCostIsTheLengthOfTheSuffixAtEveryChunkStart_AndTheCarryKeepsTheRestCheap()
    {
        const int n = 1000;
        var s = new byte[n];
        uint[] sa = Enumerable.Range(0, n + 1).Select(i => (uint)(n - i)).ToArray();
        uint[] phi = PhiOf(sa);
        using PackedText t = PackedText.Pack(s, 1);
        (PlcpChunkResult whole, uint[] plcp) = Scan(t, phi, 0, n);
        Assert.That(plcp, Is.EqualTo(Enumerable.Range(0, n).Select(j => (uint)(n - j - 1)).ToArray()));
        Assert.That(whole.MaxPlcp, Is.EqualTo(n - 1));
        Assert.That(whole.FirstValue, Is.EqualTo(n - 1));
        long restart = 0, steps = 0;
        for (int first = 0; first < n; first += 100)
        {
            (PlcpChunkResult r, _) = Scan(t, phi, first, 100);
            restart += r.FirstValue;
            steps += r.CompareSteps;
        }
        Assert.That(restart, Is.EqualTo(5490), "sum over the ten chunks of n - first - 1 symbols compared from scratch");
        Assert.That(steps, Is.LessThan(200), "one long compare per chunk start, then the carried length leaves every further position without a compare");
    }

    [Test]
    public void SmallestNonemptySuffix_HasValueZero_AndResetsTheCarry()
    {
        // "ba": SA = [2, 1, 0]; suffix 1 ("a") has the empty suffix as its predecessor
        byte[] s = [1, 0];
        uint[] sa = [2, 1, 0];
        using PackedText t = PackedText.Pack(s, 1);
        (PlcpChunkResult r, uint[] plcp) = Scan(t, PhiOf(sa), 0, 2);
        Assert.That(plcp, Is.EqualTo(new uint[] { 0, 0 }));
        Assert.That(r.Escapes8, Is.Zero);
        Assert.That(r.MaxPlcp, Is.Zero);
    }
}

[TestFixture]
public class SpoolTests
{
    private string _dir = string.Empty;

    [SetUp]
    public void SetUp()
    {
        _dir = Path.Combine(Path.GetTempPath(), "emproto-spool-" + Guid.NewGuid().ToString("N"));
        Directory.CreateDirectory(_dir);
        ScratchTracker.Reset();
        IoCounters.Reset();
    }

    [TearDown]
    public void TearDown() => Directory.Delete(_dir, true);

    [Test]
    public void Append_ReadBack_InOrder_WithTheBytesCounted_AndDeleteFreesThem()
    {
        using var spool = new Spool<RankTuple>(_dir, "t", 3);
        spool.Append(1, [new RankTuple(1, 10), new RankTuple(2, 20)]);
        spool.Append(1, [new RankTuple(3, 30)]);
        spool.Append(1, ReadOnlySpan<RankTuple>.Empty);
        Assert.That(spool.Entries(1), Is.EqualTo(3));
        Assert.That(spool.Entries(0), Is.Zero);
        Assert.That(spool.Chunks, Is.EqualTo(3));
        Assert.That(new FileInfo(spool.PathOf(1)).Length, Is.EqualTo(24));
        Assert.That(ScratchTracker.Live, Is.EqualTo(24));
        Assert.That(IoCounters.WriteBytes, Is.EqualTo(24));
        var buf = new RankTuple[2];
        using (SpoolReader<RankTuple> reader = spool.OpenReader(1))
        {
            Assert.That(reader.Count, Is.EqualTo(3));
            Assert.That(reader.Read(buf), Is.EqualTo(2));
            Assert.That((buf[0].Rank, buf[0].Value, buf[1].Rank, buf[1].Value), Is.EqualTo((1u, 10u, 2u, 20u)));
            Assert.That(reader.Read(buf), Is.EqualTo(1));
            Assert.That((buf[0].Rank, buf[0].Value), Is.EqualTo((3u, 30u)));
            Assert.That(reader.Read(buf), Is.Zero);
        }
        Assert.That(IoCounters.ReadBytes, Is.EqualTo(24));
        using (SpoolReader<RankTuple> empty = spool.OpenReader(2))
        {
            Assert.That(empty.Count, Is.Zero);
            Assert.That(empty.Read(buf), Is.Zero);
        }
        spool.Delete(1);
        Assert.That(File.Exists(spool.PathOf(1)), Is.False);
        Assert.That(ScratchTracker.Live, Is.Zero);
        Assert.That(ScratchTracker.Peak, Is.EqualTo(24));
        Assert.That(spool.Entries(1), Is.Zero);
        spool.Delete(1); // idempotent
        spool.Delete(0);
    }

    [Test]
    public void TupleSizes_AreTheDocumentedOnes()
    {
        unsafe
        {
            Assert.That(sizeof(PhiTuple), Is.EqualTo(12));
            Assert.That(sizeof(PhiPair), Is.EqualTo(8));
            Assert.That(sizeof(RankTuple), Is.EqualTo(8));
        }
    }

    [Test]
    public void Dispose_RemovesWhatIsLeft()
    {
        var spool = new Spool<PhiPair>(_dir, "t", 4);
        spool.Append(0, [new PhiPair(1, 2)]);
        spool.Append(3, [new PhiPair(3, 4)]);
        Assert.That(Directory.GetFiles(_dir), Has.Length.EqualTo(2));
        spool.Dispose();
        Assert.That(Directory.GetFiles(_dir), Is.Empty);
        Assert.That(ScratchTracker.Live, Is.Zero);
    }

    [Test]
    public void ConcurrentAppends_ToTheSameChunk_KeepEveryEntryExactlyOnce()
    {
        using var spool = new Spool<PhiPair>(_dir, "t", 2);
        Parallel.For(0, 8, new ParallelOptions { MaxDegreeOfParallelism = 4 }, task =>
        {
            for (uint k = 0; k < 500; k++) spool.Append(task % 2, [new PhiPair((uint)task, k), new PhiPair((uint)task, k + 1000)]);
        });
        Assert.That(spool.Entries(0) + spool.Entries(1), Is.EqualTo(8 * 500 * 2));
        var seen = new HashSet<(uint, uint)>();
        var buf = new PhiPair[64];
        for (long chunk = 0; chunk < 2; chunk++)
        {
            using SpoolReader<PhiPair> reader = spool.OpenReader(chunk);
            int got;
            while ((got = reader.Read(buf)) > 0)
            {
                for (int i = 0; i < got; i++) Assert.That(seen.Add((buf[i].Position, buf[i].Phi)), Is.True, "no entry twice");
            }
        }
        Assert.That(seen, Has.Count.EqualTo(8 * 500 * 2));
    }

    [Test]
    public void Writer_FlushesAtCapacity_KeepsTheRestUntilFlushAll_AndStaysInItsChunkRange()
    {
        using var arena = new BudgetedArena(new ArenaBudget(1 << 20));
        using var spool = new Spool<RankTuple>(_dir, "t", 10);
        var writer = new SpoolWriter<RankTuple>(spool, arena, 4, 3, 40); // chunks 4..6, 5 tuples per buffer
        Assert.That(writer.Capacity, Is.EqualTo(5));
        for (uint i = 0; i < 4; i++) writer.Add(4, new RankTuple(i, i));
        Assert.That(spool.Entries(4), Is.Zero, "below the capacity nothing is written");
        writer.Add(4, new RankTuple(4, 4));
        Assert.That(spool.Entries(4), Is.EqualTo(5), "a full buffer is flushed");
        writer.Add(4, new RankTuple(5, 5));
        writer.Add(6, new RankTuple(6, 6));
        Assert.That(spool.Entries(6), Is.Zero);
        writer.FlushAll();
        Assert.That(spool.Entries(4), Is.EqualTo(6));
        Assert.That(spool.Entries(5), Is.Zero);
        Assert.That(spool.Entries(6), Is.EqualTo(1));
        var buf = new RankTuple[8];
        using (SpoolReader<RankTuple> r = spool.OpenReader(4))
        {
            Assert.That(r.Read(buf), Is.EqualTo(6));
            Assert.That(buf.Take(6).Select(t => t.Rank), Is.EqualTo(new uint[] { 0, 1, 2, 3, 4, 5 }));
        }
        writer.FlushAll(); // nothing left
        Assert.That(spool.Entries(4), Is.EqualTo(6));
        Assert.Throws<IndexOutOfRangeException>(() => writer.Add(7, new RankTuple(1, 1)));
    }

    [Test]
    public void Writer_BufferIsAtLeastOneTuple_AndChargesTheBudget()
    {
        var budget = new ArenaBudget(1 << 20);
        using (var arena = new BudgetedArena(budget))
        {
            using var spool = new Spool<PhiTuple>(_dir, "t", 5);
            var writer = new SpoolWriter<PhiTuple>(spool, arena, 0, 5, 1);
            Assert.That(writer.Capacity, Is.EqualTo(1));
            Assert.That(budget.LiveBytes, Is.EqualTo(5 * 12));
            writer.Add(2, new PhiTuple(1, 2, 3));
            Assert.That(spool.Entries(2), Is.EqualTo(1), "a capacity of one flushes every tuple");
        }
        Assert.That(budget.LiveBytes, Is.Zero);
    }
}

[TestFixture]
public class LcpOutputWriterTests
{
    private string _dir = string.Empty;

    [SetUp]
    public void SetUp()
    {
        _dir = Path.Combine(Path.GetTempPath(), "emproto-out-" + Guid.NewGuid().ToString("N"));
        Directory.CreateDirectory(_dir);
    }

    [TearDown]
    public void TearDown() => Directory.Delete(_dir, true);

    private (byte[] Base, uint[] Pairs) Write(LcpKind kind, uint[] values, long overflow, int slice, int[] parts)
    {
        string b = Path.Combine(_dir, "b.bin"), o = Path.Combine(_dir, "o.bin");
        using (var arena = new BudgetedArena(new ArenaBudget(1 << 20)))
        using (var w = new LcpOutputWriter(arena, kind, values.Length - 1, overflow, b, o, 16, slice))
        {
            int at = 0;
            foreach (int part in parts)
            {
                w.Write(at, values.AsSpan(at, part));
                at += part;
            }
            w.Complete(values.Length - 1);
        }
        return (File.ReadAllBytes(b), MemoryMarshal.Cast<byte, uint>(File.ReadAllBytes(o)).ToArray());
    }

    [Test]
    public void U8_EscapesValuesFrom255_IncludingTheEscapeValueItself()
    {
        uint[] v = [0, 0, 5, 254, 255, 256, 7, 1000, 0];
        (byte[] b, uint[] pairs) = Write(LcpKind.U8, v, 3, 2, [4, 5]);
        Assert.That(b, Is.EqualTo(new byte[] { 0, 0, 5, 254, 255, 255, 7, 255, 0 }));
        Assert.That(pairs, Is.EqualTo(new uint[] { 4, 255, 5, 256, 7, 1000 }));
    }

    [Test]
    public void U16_EscapesValuesFrom65535_LittleEndian()
    {
        uint[] v = [0, 300, 65534, 65535, 70000, 1];
        (byte[] b, uint[] pairs) = Write(LcpKind.U16, v, 2, 4, [6]);
        Assert.That(MemoryMarshal.Cast<byte, ushort>(b).ToArray(), Is.EqualTo(new ushort[] { 0, 300, 65534, 65535, 65535, 1 }));
        Assert.That(pairs, Is.EqualTo(new uint[] { 3, 65535, 4, 70000 }));
    }

    [Test]
    public void U32_StoresEveryValue_WithoutAnOverflowFile()
    {
        uint[] v = [0, 0, 4_000_000_000, 70000, uint.MaxValue - 1];
        (byte[] b, uint[] pairs) = Write(LcpKind.U32, v, 0, 3, [2, 3]);
        Assert.That(MemoryMarshal.Cast<byte, uint>(b).ToArray(), Is.EqualTo(v));
        Assert.That(pairs, Is.Empty);
    }

    [Test]
    public void RanksMustArriveInOrder_AndTheAnnouncedEscapesMustMatch()
    {
        string b = Path.Combine(_dir, "b.bin"), o = Path.Combine(_dir, "o.bin");
        using var arena = new BudgetedArena(new ArenaBudget(1 << 20));
        using (var w = new LcpOutputWriter(arena, LcpKind.U8, 5, 1, b, o, 16, 4))
        {
            w.Write(0, [0u, 1u, 2u]);
            Assert.Throws<InvalidOperationException>(() => w.Write(4, [1u, 2u]), "a gap in the ranks");
            Assert.Throws<InvalidOperationException>(() => w.Complete(5), "only three of six entries written");
            Assert.Throws<InvalidOperationException>(() => w.Write(3, [1u, 300u, 400u]), "two escapes met, one announced");
        }
        using (var w = new LcpOutputWriter(arena, LcpKind.U8, 2, 1, b, o, 16, 4))
        {
            w.Write(0, [0u, 1u, 2u]);
            Assert.Throws<InvalidOperationException>(() => w.Complete(2), "one escape announced, none met");
        }
        using var w2 = new LcpOutputWriter(arena, LcpKind.U8, 2, 1, b, o, 16, 4);
        w2.Write(0, [0u, 1u, 255u]);
        Assert.DoesNotThrow(() => w2.Complete(2));
    }
}

[TestFixture]
public class ExternalLcpOptionsTests
{
    private string _dir = string.Empty;

    [SetUp]
    public void SetUp()
    {
        _dir = Path.Combine(Path.GetTempPath(), "emproto-opt-" + Guid.NewGuid().ToString("N"));
        Directory.CreateDirectory(_dir);
    }

    [TearDown]
    public void TearDown() => Directory.Delete(_dir, true);

    private static (PackedText Text, uint[] Sa) Small()
    {
        DenseCase c = BuilderTestSupport.Dense("dna", StringFamilies.RandomDna(6000, 5));
        return (PackedText.Pack(c.Symbols, c.MinBits), c.ReferenceSa());
    }

    [Test]
    public void InvalidOptions_AreRejected()
    {
        (PackedText t, uint[] sa) = Small();
        using (t)
        {
            string saPath = LcpHarness.WriteSa(sa, _dir);
            ExtLcpOptions ok = LcpHarness.Options(Path.Combine(_dir, "s"), 100);
            string b = Path.Combine(_dir, "b"), o = Path.Combine(_dir, "o");
            Assert.Throws<ArgumentOutOfRangeException>(() => ExternalLcp.Build(saPath, t, b, o, ok with { ChunkPositions = 0 }));
            Assert.Throws<ArgumentOutOfRangeException>(() => ExternalLcp.Build(saPath, t, b, o, ok with { Threads = 0 }));
            Assert.Throws<ArgumentOutOfRangeException>(() => ExternalLcp.Build(saPath, t, b, o, ok with { Workers = -1 }));
            Assert.Throws<ArgumentOutOfRangeException>(() => ExternalLcp.Build(saPath, t, b, o, ok with { ReadBufferBytes = 0 }));
            Assert.Throws<ArgumentOutOfRangeException>(() => ExternalLcp.Build(saPath, t, b, o, ok with { TupleBufferBytes = 0 }));
            Assert.Throws<ArgumentOutOfRangeException>(() => ExternalLcp.Build(saPath, t, b, o, ok with { WindowBytes = 0 }));
            Assert.Throws<ArgumentOutOfRangeException>(() => ExternalLcp.Build(saPath, t, b, o, ok with { RoundChunks = -1 }));
            Assert.Throws<ArgumentNullException>(() => ExternalLcp.Build(saPath, null!, b, o, ok));
            Assert.Throws<ArgumentNullException>(() => ExternalLcp.Build(saPath, t, b, o, null!));
            Assert.DoesNotThrow(() => ExternalLcp.Build(saPath, t, b, o, ok));
        }
    }

    // Added after audit round 1 (mutant B06: a chunk with too FEW tuples was accepted, because the chunk with one too many is found anyway - unless the single task meets the deficient one first and goes on
    // with stale buffers): with one task and the deficient chunk before the surplus one the check on the short chunk is the only thing that stops the build.
    [TestCase(false, true)]
    [TestCase(true, true)]
    [TestCase(false, false)]
    [TestCase(true, false)]
    public void NotAPermutation_OneTask_IsDetectedAtTheShortChunk_OrTheLongOne_WhicheverComesFirst(bool boundedMode, bool deficientFirst)
    {
        (PackedText t, uint[] sa) = Small(); // n = 6000
        using (t)
        {
            var g = new LcpGeometry(6000, 500);
            // sa[from] is overwritten by a copy of sa[j]: the chunk of the old sa[from] loses a tuple, the chunk of sa[j] gains one
            int from = 10;
            int j = 11;
            while (deficientFirst ? g.ChunkOf(sa[j]) <= g.ChunkOf(sa[from]) : g.ChunkOf(sa[j]) >= g.ChunkOf(sa[from]))
            {
                if (++j == sa.Length)
                {
                    from++;
                    j = from + 1;
                }
            }
            sa[from] = sa[j];
            ExtLcpOptions o = LcpHarness.Options(Path.Combine(_dir, "s"), 500, 1, boundedMode ? LcpMode.Bounded : LcpMode.Distributed, 0);
            string saPath = LcpHarness.WriteSa(sa, _dir);
            var ex = Assert.Throws<InvalidOperationException>(() => ExternalLcp.Build(saPath, t, Path.Combine(_dir, "b"), Path.Combine(_dir, "o"), o));
            Assert.That(ex!.Message, Does.Contain("Phi tuples, expected"), "the check on the number of tuples of a chunk, not a consequence of stale buffers");
            Assert.That(Directory.GetFileSystemEntries(Path.Combine(_dir, "s")), Is.Empty);
        }
    }

    [TestCase(false)]
    [TestCase(true)]
    public void NotAPermutation_IsDetected_WhenTheChunksDisagreeOnTheirTupleCounts(bool boundedMode)
    {
        (PackedText t, uint[] sa) = Small(); // n = 6000
        using (t)
        {
            var g = new LcpGeometry(6000, 500);
            int j = 11;
            while (g.ChunkOf(sa[j]) == g.ChunkOf(sa[10])) j++;
            sa[10] = sa[j]; // the position sa[j] appears twice, the old sa[10] never: one chunk has a tuple too many, another one too few
            ExtLcpOptions o = LcpHarness.Options(Path.Combine(_dir, "s"), 500, 2, boundedMode ? LcpMode.Bounded : LcpMode.Distributed, 4);
            string saPath = LcpHarness.WriteSa(sa, _dir);
            Assert.Throws<InvalidOperationException>(() => ExternalLcp.Build(saPath, t, Path.Combine(_dir, "b"), Path.Combine(_dir, "o"), o));
            Assert.That(Directory.GetFileSystemEntries(Path.Combine(_dir, "s")), Is.Empty, "scratch removed");
            Assert.That(File.Exists(Path.Combine(_dir, "b")), Is.False);
        }
    }

    [Test]
    public void ScratchPeak_IsTwelveBytesPerSymbolInTheBaseMode_AndEightPlusFourOverRoundsInTheBoundedMode()
    {
        (PackedText t, uint[] sa) = Small(); // n = 6000
        using (t)
        {
            ExtLcpOptions o = LcpHarness.Options(Path.Combine(_dir, "s"), 500, 4);
            ScratchTracker.Reset();
            LcpHarness.Run(t, sa, o, _dir);
            Assert.That(ScratchTracker.Peak, Is.EqualTo(12 * 6000), "all Phi tuples (12 B) are on disk when B1 ends; every chunk file is deleted before its rank tuples (8 B) are written");
            Assert.That(ScratchTracker.Live, Is.Zero);
            ScratchTracker.Reset();
            LcpHarness.Run(t, sa, o with { Mode = LcpMode.Bounded }, _dir);
            Assert.That(ScratchTracker.Peak, Is.EqualTo(8 * 6000), "one round: all Phi pairs (8 B); the PLCP file (4 B) grows as the pair files go");
            ScratchTracker.Reset();
            LcpHarness.Run(t, sa, o with { Mode = LcpMode.Bounded, RoundChunks = 6 }, _dir); // two rounds of 3000 positions
            Assert.That(ScratchTracker.Peak, Is.EqualTo((4 * 3000) + (8 * 3000)), "two rounds: the PLCP of the first half (4 B) and the pairs of the second half (8 B)");
            Assert.That(ScratchTracker.Live, Is.Zero);
        }
    }

    // Added after audit round 1 (mutant B26: the restart cost of the tasks was not summed).
    [TestCase(false)]
    [TestCase(true)]
    public void Stats_RestartSymbolsOfAllSame_AreTheSuffixLengthsAtTheChunkStarts_SummedOverTheTasks(bool boundedMode)
    {
        const int n = 1000;
        var symbols = new byte[n];
        uint[] sa = Enumerable.Range(0, n + 1).Select(i => (uint)(n - i)).ToArray();
        using PackedText t = PackedText.Pack(symbols, 1);
        ExtLcpOptions o = LcpHarness.Options(Path.Combine(_dir, "s"), 100, 3, boundedMode ? LcpMode.Bounded : LcpMode.Distributed, 0);
        LcpOutcome got = LcpHarness.Run(t, sa, o, _dir);
        Assert.That(got.Stats.RestartSymbols, Is.EqualTo(5490), "ten chunk starts compare n - start - 1 = 999, 899, ..., 99 symbols from scratch");
        Assert.That(got.Stats.MaxPlcp, Is.EqualTo(n - 1));
        Assert.That(got.Stats.Escapes8, Is.EqualTo(n - 255), "the values 255 .. 999 of a^1000");
        Assert.That(got.Stats.Escapes16, Is.Zero);
        Assert.That(got.Stats.CompareSteps, Is.GreaterThan(0));
    }

    [Test]
    public void ParallelTasks_AreLimitedByTheBudget_TheChunkCount_AndTheExplicitWorkers()
    {
        (PackedText t, uint[] sa) = Small();
        using (t)
        {
            ExtLcpOptions o = LcpHarness.Options(Path.Combine(_dir, "s"), 500, 4) with { ReadBufferBytes = 1024, TupleBufferBytes = 256 };
            long per = ExternalLcp.WorkerBytes(6000, o);
            Assert.That(per, Is.EqualTo((8 * 500) + 1020 + (13 * 256)), "Phi + rank (4 B each), 85 tuples of 12 B read, 13 rank chunks of 32 tuples of 8 B");
            Assert.That(ExternalLcp.WorkerBytes(6000, o with { Mode = LcpMode.Bounded }), Is.EqualTo((4 * 500) + (128 * 8)));
            string saPath = LcpHarness.WriteSa(sa, _dir);
            string b = Path.Combine(_dir, "b"), ov = Path.Combine(_dir, "o");
            ExtLcpStats Run(ExtLcpOptions x) => ExternalLcp.Build(saPath, t, b, ov, x);
            Assert.That(Run(o with { Budget = new ArenaBudget(per * 25 / 10) }).Workers, Is.EqualTo(2), "room for two tasks besides the other buffers");
            Assert.That(Run(o with { Budget = new ArenaBudget(per * 100) }).Workers, Is.EqualTo(4), "threads");
            Assert.That(Run(o with { Workers = 3, Budget = new ArenaBudget(per * 100) }).Workers, Is.EqualTo(3), "explicit");
            Assert.That(Run(o with { Workers = 3, Threads = 1 }).Workers, Is.EqualTo(3), "explicit beats threads");
            Assert.That(Run(o with { ChunkPositions = 3000 }).Workers, Is.EqualTo(2), "two chunks, two tasks");
            Assert.That(Run(o with { ChunkPositions = 100_000 }).Workers, Is.EqualTo(1), "one chunk");
            Assert.That(Run(o with { Threads = 1 }).Workers, Is.EqualTo(1));
            Assert.That(Run(o with { Budget = new ArenaBudget(1 << 20), Workers = 0, Threads = 4 }).Workers, Is.EqualTo(4));
        }
    }
}

[TestFixture]
public class E4CommandTests
{
    private string _dir = string.Empty;

    [SetUp]
    public void SetUp()
    {
        _dir = Path.Combine(Path.GetTempPath(), "emproto-e4-" + Guid.NewGuid().ToString("N"));
        Directory.CreateDirectory(_dir);
    }

    [TearDown]
    public void TearDown() => Directory.Delete(_dir, true);

    private (string Text, string Sa) Prepare(string family)
    {
        string pk = Path.Combine(_dir, "t.pk");
        string ranks = Path.Combine(_dir, "t.ranks");
        string sa = Path.Combine(_dir, "t.sa");
        var report = new Report(null);
        Assert.That(E3.Pack(Options.Parse(["--family", family, "--out", pk]), report), Is.EqualTo(0));
        Assert.That(E3.MakeRanks(Options.Parse(["--text", pk, "--out", ranks, "--period", "512"]), report), Is.EqualTo(0));
        Assert.That(E3.Build(Options.Parse(["--text", pk, "--ranks", ranks, "--sa", sa, "--scratch", Path.Combine(_dir, "e3scratch"), "--threads", "2", "--k-div", "500", "--budget-mb", "16", "--rss-cap-mb", "100000"]), report), Is.EqualTo(0));
        return (pk, sa);
    }

    [TestCase("genome:200000", "distributed", "")]
    [TestCase("genome:200000", "bounded", "3")]
    [TestCase("tandem:171x300", "distributed", "")]
    [TestCase("a:60000", "bounded", "2")]
    [TestCase("a:70000", "distributed", "")]
    public void Pipeline_SaFile_E4_Oracle_Verify_AndCheck(string family, string mode, string rounds)
    {
        if (!File.Exists("/proc/self/stat")) Assert.Ignore("the E4 report reads /proc (Linux)");
        (string pk, string sa) = Prepare(family);
        string rows = Path.Combine(_dir, "rows.tsv");
        var report = new Report(rows);
        string oracle = Path.Combine(_dir, "oracle");
        Assert.That(E4.Oracle(Options.Parse(["--text", pk, "--sa", sa, "--out", oracle]), report), Is.EqualTo(0));
        string[] args = ["--text", pk, "--sa", sa, "--out", Path.Combine(_dir, "out"), "--scratch", Path.Combine(_dir, "scratch"), "--scale-n", "3000000000", "--mode", mode, "--threads", "3", "--rss-cap-mb", "100000", "--verify", oracle, "--case", "t", "--io-kb", "8", "--tuple-kb", "1", "--window-kb", "4"];
        if (rounds.Length > 0) args = [.. args, "--rounds", rounds];
        Assert.That(E4.Build(Options.Parse(args), report), Is.EqualTo(0), "equal to the kernel, no scratch left, budget balanced");
        string text = File.ReadAllText(rows);
        Assert.That(text, Does.Contain("\tlcp_equals_oracle\tyes"));
        Assert.That(text, Does.Contain("\tparameters\t"));
        Assert.That(text, Does.Contain("P=375"), "the scale-down keeps the chunk count of 3x10^9 symbols");
        Assert.That(text, Does.Contain(mode == "bounded" ? "B2 workers=3" : "B2 workers=2"), "workers from the reference geometry: two fit the arena cap in the base mode, the bounded mode takes min(threads, 5)");
        if (rounds.Length > 0) Assert.That(text, Does.Contain($"rounds={rounds} "));
        Assert.That(Directory.GetFileSystemEntries(Path.Combine(_dir, "scratch")), Is.Empty);
        Assert.That(E4.CheckSample(Options.Parse(["--text", pk, "--sa", sa, "--lcp", Path.Combine(_dir, "out"), "--samples", "300"]), report), Is.EqualTo(0));
        Assert.That(File.ReadAllText(rows), Does.Contain("\tstructure_problems\t0\t"));
        Assert.That(File.ReadAllText(rows), Does.Contain("\tsample_mismatches\t0\t"));
    }

    [Test]
    public void Verify_AgainstAWrongOracle_ReturnsOne_AndSaysNo()
    {
        (string pk, string sa) = Prepare("genome:100000");
        string rows = Path.Combine(_dir, "rows.tsv");
        var report = new Report(rows);
        string oracle = Path.Combine(_dir, "oracle");
        Assert.That(E4.Oracle(Options.Parse(["--text", pk, "--sa", sa, "--out", oracle]), report), Is.EqualTo(0));
        byte[] lcp = File.ReadAllBytes(oracle + ".lcp");
        lcp[5000] ^= 1;
        File.WriteAllBytes(oracle + ".lcp", lcp);
        string[] args = ["--text", pk, "--sa", sa, "--out", Path.Combine(_dir, "out"), "--scratch", Path.Combine(_dir, "scratch"), "--threads", "2", "--rss-cap-mb", "100000", "--verify", oracle];
        Assert.That(E4.Build(Options.Parse(args), report), Is.EqualTo(1));
        Assert.That(File.ReadAllText(rows), Does.Contain("\tlcp_equals_oracle\tNO\t").And.Contain("first difference: base at byte 5000"));
        Assert.That(E4.Build(Options.Parse(["--text", pk, "--sa", sa, "--out", Path.Combine(_dir, "out"), "--scratch", Path.Combine(_dir, "scratch"), "--threads", "2", "--rss-cap-mb", "0.001"]), report), Is.EqualTo(1), "the RssAnon cap of the gate");
    }

    [Test]
    public void CheckSample_DetectsACorruptedOverflowListAndACorruptedBaseEntry()
    {
        (string pk, string sa) = Prepare("a:70000");
        var report = new Report(null);
        string outPrefix = Path.Combine(_dir, "out");
        Assert.That(E4.Build(Options.Parse(["--text", pk, "--sa", sa, "--out", outPrefix, "--scratch", Path.Combine(_dir, "scratch"), "--chunk", "1000", "--threads", "2", "--rss-cap-mb", "100000"]), report), Is.EqualTo(0));
        Assert.That(new FileInfo(outPrefix + ".ovf").Length, Is.GreaterThan(0), "a^n has escapes");
        string[] check = ["--text", pk, "--sa", sa, "--lcp", outPrefix, "--samples", "500"];
        Assert.That(E4.CheckSample(Options.Parse(check), report), Is.EqualTo(0));

        byte[] ovf = File.ReadAllBytes(outPrefix + ".ovf");
        byte[] saved = (byte[])ovf.Clone();
        ovf[^4] ^= 1; // the value of the last pair, the rank n, which is always sampled
        File.WriteAllBytes(outPrefix + ".ovf", ovf);
        Assert.That(E4.CheckSample(Options.Parse(check), report), Is.EqualTo(1), "a wrong overflow value is found by the samples");
        ovf = (byte[])saved.Clone();
        ovf[8] ^= 1; // the rank of the second pair
        File.WriteAllBytes(outPrefix + ".ovf", ovf);
        Assert.That(E4.CheckSample(Options.Parse(check), report), Is.EqualTo(1), "a pair that does not belong to its escape slot is found by the structure scan");
        ovf = (byte[])saved.Clone();
        ovf[8] ^= 2; // the rank of the second pair becomes larger than its slot (65,539 instead of 65,537): only the equality of rank and slot finds it
        File.WriteAllBytes(outPrefix + ".ovf", ovf);
        Assert.That(E4.CheckSample(Options.Parse(check), report), Is.EqualTo(1), "a pair that is ahead of its slot");
        File.WriteAllBytes(outPrefix + ".ovf", saved[..^8]);
        Assert.That(E4.CheckSample(Options.Parse(check), report), Is.EqualTo(1), "an escape slot without a pair");
        File.WriteAllBytes(outPrefix + ".ovf", saved);
        byte[] lcp = File.ReadAllBytes(outPrefix + ".lcp");
        lcp[132_000] = 0; // the entry of rank 66,000 held the escape value 65,535: now it does not, and its pair has no slot
        File.WriteAllBytes(outPrefix + ".lcp", lcp);
        Assert.That(E4.CheckSample(Options.Parse(check), report), Is.EqualTo(1));
    }

    // Added after audit round 1 (mutant E06: the closing comparison of the structure scan was not covered, because every corruption in the middle is also found by the pair that follows it).
    [Test]
    public void CheckSample_StructureScan_FindsAnUnusedLastPair()
    {
        (string pk, string sa) = Prepare("a:70000");
        string rows = Path.Combine(_dir, "rows.tsv");
        var report = new Report(rows);
        string outPrefix = Path.Combine(_dir, "out");
        Assert.That(E4.Build(Options.Parse(["--text", pk, "--sa", sa, "--out", outPrefix, "--scratch", Path.Combine(_dir, "scratch"), "--chunk", "1000", "--threads", "2", "--rss-cap-mb", "100000"]), report), Is.EqualTo(0));
        byte[] lcp = File.ReadAllBytes(outPrefix + ".lcp");
        lcp[2 * 70_000] = 0; // the entry of rank n held the escape value 65,535: the last pair has no slot any more
        lcp[(2 * 70_000) + 1] = 0;
        File.WriteAllBytes(outPrefix + ".lcp", lcp);
        File.Delete(rows);
        Assert.That(E4.CheckSample(Options.Parse(["--text", pk, "--sa", sa, "--lcp", outPrefix, "--samples", "10"]), report), Is.EqualTo(1));
        Assert.That(File.ReadAllText(rows), Does.Contain("\tstructure_problems\t1\t"), "found by the structure scan itself, not only by the sample of the last rank");

        // a pair whose value is below the escape value (mutant E09: the value of a pair was not required to reach the escape): the last pair, at the rank of the last escape slot
        byte[] saved = File.ReadAllBytes(outPrefix + ".ovf");
        lcp[2 * 70_000] = 0xFF;
        lcp[(2 * 70_000) + 1] = 0xFF; // the slot is an escape again
        File.WriteAllBytes(outPrefix + ".lcp", lcp);
        byte[] ovf = (byte[])saved.Clone();
        ovf[^4] = 100;
        ovf[^3] = ovf[^2] = ovf[^1] = 0; // value 100 < 65,535
        File.WriteAllBytes(outPrefix + ".ovf", ovf);
        File.Delete(rows);
        Assert.That(E4.CheckSample(Options.Parse(["--text", pk, "--sa", sa, "--lcp", outPrefix, "--samples", "10"]), report), Is.EqualTo(1));
        Assert.That(File.ReadAllText(rows), Does.Contain("\tstructure_problems\t1\t"), "a pair value below the escape is a structure problem");
    }

    // Added after audit round 1 (mutant E05: the exit code ignored a scratch entry left behind - the condition is now a function of its own).
    [Test]
    public void Healthy_RequiresTheCapAndNothingLeftBehind()
    {
        Assert.That(E4.Healthy(100, 256, 0, 0, 0), Is.True);
        Assert.That(E4.Healthy(256, 256, 0, 0, 0), Is.True, "the cap itself is allowed");
        Assert.That(E4.Healthy(256.001, 256, 0, 0, 0), Is.False, "above the cap");
        Assert.That(E4.Healthy(100, 256, 1, 0, 0), Is.False, "a scratch entry is left in the directory");
        Assert.That(E4.Healthy(100, 256, 0, 1, 0), Is.False, "a tracked scratch byte is left");
        Assert.That(E4.Healthy(100, 256, 0, 0, 1), Is.False, "a native byte is still charged to the budget");
        Assert.That(E4.Healthy(300, 256, 1, 1, 1), Is.False);
    }

    [Test]
    public void Commands_RejectABadModeAndMissingOptions_AndTheCliDispatchesThem()
    {
        var report = new Report(null);
        Assert.Throws<ArgumentException>(() => E4.Build(Options.Parse(["--text", "x", "--sa", "y", "--out", "z", "--scratch", "s", "--mode", "fast"]), report));
        Assert.Throws<ArgumentException>(() => E4.Build(Options.Parse(["--text", "x"]), report));
        Assert.That(Cli.Run(["e4", "--text", "x"]), Is.EqualTo(2));
        Assert.That(Cli.Run(["lcp-oracle"]), Is.EqualTo(2));
        Assert.That(Cli.Run(["check-lcp"]), Is.EqualTo(2));
    }

    [Test]
    public void ReferenceGeometry_IsThatOfThe3e9Build_AndGivesTwoTasksUnderTheArenaCap()
    {
        Assert.That(E4.ReferenceN, Is.EqualTo(3_000_000_000));
        Assert.That(E4.ChunkRef, Is.EqualTo(8_000_000));
        var o = new ExtLcpOptions { ScratchDir = "x", ChunkPositions = E4.ChunkRef };
        long per = ExternalLcp.WorkerBytes(E4.ReferenceN, o);
        Assert.That(per, Is.EqualTo((8 * 8_000_000L) + (87_381L * 12) + (376L * 32768)), "Phi + rank chunk, 87,381 tuples of 12 B read, 376 rank chunks (n + 1 ranks) of 32 KiB");
        Assert.That(160L * 1048576 / per, Is.EqualTo(2), "two B2 tasks fit the 160 MiB arena cap");
        Assert.That(160L * 1048576 / ExternalLcp.WorkerBytes(E4.ReferenceN, o with { Mode = LcpMode.Bounded }), Is.EqualTo(5), "the bounded mode holds the Phi chunk only: five tasks fit, the command takes min(threads, fit) = 4");
    }

    [Test]
    public void Cli_Usage_ListsTheNewCommands()
    {
        var sw = new StringWriter();
        TextWriter old = Console.Out;
        Console.SetOut(sw);
        try
        {
            Assert.That(Cli.Run(["--help"]), Is.EqualTo(0));
        }
        finally
        {
            Console.SetOut(old);
        }
        Assert.That(sw.ToString(), Does.Contain("e4 ").And.Contain("lcp-oracle").And.Contain("check-lcp"));
    }
}
