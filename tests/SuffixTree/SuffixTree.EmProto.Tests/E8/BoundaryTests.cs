using SuffixTree.EmProto.Build;
using SuffixTree.EmProto.Dcs;
using SuffixTree.EmProto.Experiments;
using SuffixTree.EmProto.Harness;
using SuffixTree.EmProto.Ph4;
using SuffixTree.EmProto.Ph5;
using SuffixTree.EmProto.Pipeline;
using SuffixTree.EmProto.Text;
using SuffixTree.EmProto.Tests.WholePipeline;
using SuffixTree.Esa;
using SplitMix64 = SuffixTree.Esa.Tests.TestSupport.SplitMix64;

namespace SuffixTree.EmProto.Tests.E8;

/// <summary>
/// The boundary audit of the E8 package (ADR-004 section 5a, E8): everything of the pipeline that can be exercised cheaply with synthetic values above 2^31 (and up to the 2^32 - 3 limit of the 32-bit positions): parameter arithmetic of the scale rule at n = 3x10^9,
/// chunk geometry, the fast division, the sample slots and the rank-table lookup, the packed text at symbol positions above 2^31 (bit offsets above 2^32), and the mapped writers at file offsets above 2^32 bytes (sparse files, nothing is allocated on disk).
/// What these tests cannot reach (a run over 2^31 symbols end to end) stays unverified until the owner's run (E8-RUNBOOK.md).
/// </summary>
[TestFixture]
public class BoundaryTests : TempDirFixture
{
    private const long Three = 3_000_000_000;

    // ---- the scale rule

    [Test]
    public void ScaleRules_AtTheDesignPoint_AreTheReferenceValues()
    {
        EndToEndParameters p = ScaleRules.Resolve(Three, new EndToEndOptions { InputPath = "x", OutPrefix = "o", ScratchDir = "s", ScaleN = Three });
        Assert.That(p.Factor, Is.EqualTo(1.0));
        Assert.That((p.DcsBucketCapacity, p.DcsBaseCase, p.DcsChunkSlots, p.DcsSortCapacity), Is.EqualTo((E3b.KsRef, E3b.BaseRef, E3b.ChunkRef, E3b.SortKRef)));
        Assert.That(p.DcsRunSpillBytes, Is.EqualTo(ScaleRules.RunSpillReference));
        Assert.That(p.BucketCapacity, Is.EqualTo(ScaleRules.BucketReference));
        Assert.That(p.LcpChunk, Is.EqualTo(E4.ChunkRef));
        Assert.That((Three + p.LcpChunk - 1) / p.LcpChunk, Is.EqualTo(375), "P = n / R chunks of PH4 at 3x10^9");
        Assert.That(p.LcpWorkers, Is.InRange(1, 4));
    }

    [Test]
    public void ScaleRules_AtTheProtocolSizes_ScaleByNOverThreeBillion_WithoutOverflow()
    {
        foreach (long n in new[] { 1_000_000_000L, 2_147_483_648L, 2_200_000_000L, 2_999_999_999L, 3_000_000_000L, 4_294_967_293L })
        {
            EndToEndParameters p = ScaleRules.Resolve(n, new EndToEndOptions { InputPath = "x", OutPrefix = "o", ScratchDir = "s", ScaleN = Three });
            double f = n / (double)Three;
            Assert.That(p.DcsBucketCapacity, Is.EqualTo((long)Math.Round(E3b.KsRef * f)), $"n = {n}");
            Assert.That(p.DcsBaseCase, Is.EqualTo((long)Math.Round(E3b.BaseRef * f)), $"n = {n}");
            Assert.That(p.DcsChunkSlots, Is.EqualTo((long)Math.Round(E3b.ChunkRef * f)), $"n = {n}");
            Assert.That(p.BucketCapacity, Is.EqualTo((long)((decimal)n * ScaleRules.BucketReference / Three)), $"n = {n}");
            Assert.That(p.LcpChunk, Is.EqualTo((long)Math.Ceiling(E4.ChunkRef * f)), $"n = {n}");
            Assert.That(p.DcsRunSpillBytes, Is.GreaterThan(0));
            long chunks = (n + p.LcpChunk - 1) / p.LcpChunk;
            Assert.That(chunks, Is.InRange(374, 376), $"n = {n}: the number of PH4 chunks stays 375 (the real regime of E7 series B)");
            Assert.That(p.LcpChunk, Is.LessThan(uint.MaxValue), "chunk positions fit the fast 32-bit division");
            Assert.That(new[] { p.FanOut, p.Oversample, p.SpillBytes, p.LcpReadBytes, p.LcpTupleBytes, p.BlockEntries, p.ChildWindowBytes, p.ChildSlots, p.WindowBytes, p.PrepareBlockSymbols, p.DcsFanOut, p.DcsReservoir, p.DcsIoBytes }, Has.All.GreaterThan(0));
        }
    }

    [Test]
    public void ScaleRules_RuleN_FixesTheParameters_WhateverTheLength()
    {
        EndToEndParameters a = ScaleRules.Resolve(Three, new EndToEndOptions { InputPath = "x", OutPrefix = "o", ScratchDir = "s", ScaleN = Three, RuleN = 100_000_000 });
        EndToEndParameters b = ScaleRules.Resolve(100_000_000, new EndToEndOptions { InputPath = "x", OutPrefix = "o", ScratchDir = "s", ScaleN = Three });
        Assert.That((a.BucketCapacity, a.LcpChunk, a.DcsBucketCapacity), Is.EqualTo((b.BucketCapacity, b.LcpChunk, b.DcsBucketCapacity)));
    }

    [Test]
    public void LcpWidthRule_And_Geometry_UseSixtyFourBitArithmetic_AtThreeBillion()
    {
        Assert.That(LcpWidth.Select(Three, 187_500_000, 0), Is.EqualTo(LcpKind.U8));
        Assert.That(LcpWidth.Select(Three, 187_500_001, 0), Is.EqualTo(LcpKind.U16));
        Assert.That(LcpWidth.Select(Three, 375_000_001, 375_000_001), Is.EqualTo(LcpKind.U32));
        Assert.That(LcpWidth.OverflowCount(LcpKind.U8, 5_000_000_000, 1), Is.EqualTo(5_000_000_000));
        var g = new LcpGeometry(Three, E4.ChunkRef);
        Assert.That((g.TextChunks, g.RankChunks), Is.EqualTo((375L, 376L)));
        Assert.That(g.TextChunkLength(374), Is.EqualTo(Three - (374 * E4.ChunkRef)));
        Assert.That(g.RankChunkLength(375), Is.EqualTo(1));
        Assert.That(g.ChunkOf(2_999_999_999u), Is.EqualTo(374u));
        Assert.That(g.ChunkOf(uint.MaxValue - 3), Is.EqualTo((uint)((uint.MaxValue - 3) / E4.ChunkRef)));
        Assert.Throws<NotSupportedException>(() => _ = new LcpGeometry(uint.MaxValue, E4.ChunkRef));
        var edge = new LcpGeometry(uint.MaxValue - 1, 1000);
        Assert.That(edge.RankChunks, Is.EqualTo((uint.MaxValue - 1L + 1000) / 1000));
    }

    [Test]
    public void FastMod32_IsExact_ForEveryNumeratorUpToTwoToThePowerThirtyTwo()
    {
        var rng = new SplitMix64(7);
        foreach (uint d in new uint[] { 2, 3, 7, 1024, 1025, 266_667, 4_194_304, 8_000_000, 2_147_483_649, 4_294_967_295 })
        {
            var m = new FastMod32(d);
            var values = new List<uint> { 0, 1, d - 1, d, d + 1, int.MaxValue, 1u << 31, (1u << 31) + 1, uint.MaxValue - 3, uint.MaxValue - 1, uint.MaxValue };
            for (int i = 0; i < 20_000; i++) values.Add((uint)rng.NextInt(int.MaxValue) + (uint)rng.NextInt(int.MaxValue) + (uint)rng.NextInt(3));
            foreach (uint a in values)
            {
                Assert.That(m.Div(a), Is.EqualTo(a / d), $"{a} / {d}");
                Assert.That(m.Mod(a), Is.EqualTo(a % d), $"{a} % {d}");
            }
        }
    }

    // ---- the sample slots and the rank table above 2^31

    [Test]
    public void SampleGeometry_AtThreeBillion_SlotsAndPositionsInvertEachOther()
    {
        DifferenceCover cover = DifferenceCover.Create(DifferenceCover.NearestSingerPeriod(1024));
        var g = new SampleGeometry(Three, cover);
        Assert.That(g.Samples, Is.EqualTo(cover.SampleCount(Three)));
        Assert.That(g.Samples, Is.LessThan(uint.MaxValue), "sample indices are 32-bit in the rank file and the names");
        Assert.That(g.ClassStart[g.Classes], Is.EqualTo(g.Samples));
        Assert.That(g.Slots, Is.GreaterThanOrEqualTo(g.Samples));
        var rng = new SplitMix64(3);
        for (int i = 0; i < 20_000; i++)
        {
            long p = i < 100 ? Three - i : 2_100_000_000L + (long)(rng.Next() % 900_000_000UL);
            if (cover.IndexOf[(int)(p % cover.Period)] < 0) continue;
            long slot = g.SlotOf(p);
            Assert.That(g.PositionOfSlot(slot), Is.EqualTo(p), $"p = {p}");
            Assert.That(slot, Is.LessThan(g.Slots));
            Assert.That(g.TrySlotOfRed(g.RedIndexOfSlot(slot), out long back), Is.True);
            Assert.That(back, Is.EqualTo(slot));
        }
    }

    [Test]
    public void DcsRankTable_Lookups_AboveTwoToThePowerThirtyOne_ReadTheRightSlot()
    {
        DifferenceCover cover = DifferenceCover.Create(DifferenceCover.NearestSingerPeriod(1024));
        var g = new SampleGeometry(Three, cover);
        using var block = new NativeBlock(g.Slots * 4); // zero pages, nothing is touched but the slots written below
        unsafe
        {
            uint* ranks = (uint*)block.Pointer;
            var rng = new SplitMix64(5);
            var checkedPositions = new List<(long P, uint V)>();
            for (int i = 0; i < 5_000; i++)
            {
                long p = 2_147_483_000L + (long)(rng.Next() % 852_000_000UL);
                if (cover.IndexOf[(int)(p % cover.Period)] < 0) continue;
                uint v = (uint)(i + 1) * 40_503u;
                ranks[g.SlotOf(p)] = v;
                checkedPositions.Add((p, v));
            }
            using var table = new DcsRankTable(cover, new NonOwningBlock(block), Three);
            Assert.That(table.Slots, Is.EqualTo(g.Slots));
            foreach ((long p, uint v) in checkedPositions) Assert.That(table.Rank(p), Is.EqualTo(v), $"p = {p}");
        }
    }

    // ---- the packed text above 2^31 symbols

    [TestCase(2, 4_000_000_000L)]
    [TestCase(3, 3_000_000_000L)]
    [TestCase(8, 2_300_000_000L)]
    public unsafe void PackedText_SymbolsWindowsKeysAndLce_AtPositionsAboveTwoToThePowerThirtyOne(int bits, long n)
    {
        long bytes = PackedText.TotalBytes(n, bits);
        using var block = new NativeBlock(bytes); // calloc: lazily zeroed pages
        var vectors = new PackedText(block.Pointer, n, bits, block, ownsOwner: false);
        var scalar = new PackedText(block.Pointer, n, bits, block, ownsOwner: false) { UseVectors = false };
        long a = 2_150_000_000L, b = 2_250_000_000L, tail = n - 5000;
        // plant a pattern of 300 symbols at a and at b (equal), then a different symbol after it; and a run at the tail
        var rng = new SplitMix64(1);
        int symbols = 1 << bits;
        var pattern = new byte[300];
        for (int i = 0; i < pattern.Length; i++) pattern[i] = (byte)rng.NextInt(symbols);
        Plant(block, bits, a, pattern);
        Plant(block, bits, b, pattern);
        Plant(block, bits, a + 300, [1]);
        Plant(block, bits, b + 300, [2]);
        Plant(block, bits, tail, pattern.AsSpan(0, 100).ToArray());
        for (int i = 0; i < pattern.Length; i++)
        {
            Assert.That(vectors.Symbol(a + i), Is.EqualTo((uint)pattern[i]), $"vector text, symbol {a + i}");
            Assert.That(vectors.Symbol(b + i), Is.EqualTo((uint)pattern[i]));
        }
        Assert.That(vectors.Window(a) >> (64 - bits), Is.EqualTo((ulong)pattern[0]));
        Assert.That(vectors.Key(a), Is.EqualTo(vectors.Key(b)));
        foreach (PackedText t in new[] { vectors, scalar })
        {
            long steps = 0;
            Assert.That(t.Lce(a, b, 5000, ref steps), Is.EqualTo(300), $"LCE of two positions above 2^31 (vectors {t.UseVectors})");
            Assert.That(t.CompareSuffixes(a, b, 0, ref steps, out long lce), Is.EqualTo(-1), "symbol 1 < symbol 2 after the common part");
            Assert.That(lce, Is.EqualTo(300));
            Assert.That(t.CompareSuffixes(tail, a, 0, ref steps, out long lce2), Is.Not.EqualTo(0));
            Assert.That(lce2, Is.EqualTo(100).Or.LessThan(100));
        }
        // the last symbols and the padding
        Assert.That(vectors.Symbol(n - 1), Is.EqualTo(0u));
        Assert.That(PackedText.DataBytes(n, bits) + PackedText.PadBytes, Is.EqualTo(bytes));
    }

    // ---- mapped writers and readers at file offsets above 2^32 bytes (sparse files)

    [Test]
    public void MappedWindowWriter_WritesAtOffsets_AboveTwoToThePowerThirtyTwo()
    {
        string path = Path.Combine(Run("mw"), "big.bin");
        const long length = 5_000_000_000L;
        using (var w = new MappedWindowWriter(path, length, 1 << 20))
        {
            w.WriteAt(4_900_000_000L, [1, 2, 3, 4]);
            w.WriteAt(4_294_967_295L, [9, 8]);
            w.WriteAt(length - 2, [7, 7]);
        }
        Assert.That(new FileInfo(path).Length, Is.EqualTo(length));
        using var fs = new FileStream(path, FileMode.Open, FileAccess.Read);
        var buf = new byte[4];
        fs.Position = 4_900_000_000L;
        fs.ReadExactly(buf);
        Assert.That(buf, Is.EqualTo(new byte[] { 1, 2, 3, 4 }));
        fs.Position = 4_294_967_295L;
        fs.ReadExactly(buf.AsSpan(0, 2));
        Assert.That(buf[..2], Is.EqualTo(new byte[] { 9, 8 }));
        fs.Position = length - 2;
        fs.ReadExactly(buf.AsSpan(0, 2));
        Assert.That(buf[..2], Is.EqualTo(new byte[] { 7, 7 }));
        File.Delete(path);
    }

    [Test]
    public void AlignedWindowWriter_Write32_AtEntryIndices_AboveTwoToThePowerThirty()
    {
        string path = Path.Combine(Run("aw"), "child.bin");
        const long entries = 1_300_000_000L; // 5.2 GB
        using (var w = new AlignedWindowWriter(path, 4 * entries, 1 << 20, 4, flushOnRelease: false))
        {
            w.Write32(1_250_000_000L, 0xDEADBEEF);
            w.Write32(1_073_741_824L, 0x01020304);
            w.Write32(entries - 1, 0xCAFEF00D);
        }
        using var fs = new FileStream(path, FileMode.Open, FileAccess.Read);
        Assert.That(fs.Length, Is.EqualTo(4 * entries));
        var buf = new byte[4];
        foreach ((long index, uint value) in new[] { (1_250_000_000L, 0xDEADBEEF), (1_073_741_824L, 0x01020304u), (entries - 1, 0xCAFEF00Du) })
        {
            fs.Position = 4 * index;
            fs.ReadExactly(buf);
            Assert.That(BitConverter.ToUInt32(buf), Is.EqualTo(value), $"entry {index}");
        }
        File.Delete(path);
    }

    [Test]
    public void MappedSaSink_WritesSlices_AtEntryOffsets_AboveTwoToThePowerThirty()
    {
        string path = Path.Combine(Run("sa"), "sa.bin");
        const long entries = 1_300_000_001L;
        using (var sink = new MappedSaSink(path, entries))
        {
            sink.Write(1_250_000_000L, [11, 22, 33]);
            sink.Write(entries - 2, [44, 55]);
        }
        using var fs = new FileStream(path, FileMode.Open, FileAccess.Read);
        Assert.That(fs.Length, Is.EqualTo(4 * entries));
        var buf = new byte[12];
        fs.Position = 4 * 1_250_000_000L;
        fs.ReadExactly(buf);
        Assert.That(BitConverter.ToUInt32(buf, 8), Is.EqualTo(33u));
        fs.Position = 4 * (entries - 2);
        fs.ReadExactly(buf.AsSpan(0, 8));
        Assert.That(BitConverter.ToUInt32(buf, 4), Is.EqualTo(55u));
        File.Delete(path);
    }

    [Test]
    public void MappedFile_ReadsBytes_AtOffsetsAboveTwoToThePowerThirtyTwo()
    {
        string path = Path.Combine(Run("mf"), "big.bin");
        using (var fs = new FileStream(path, FileMode.Create, FileAccess.Write))
        {
            fs.SetLength(5_368_709_120L);
            fs.Position = 5_000_000_001L;
            fs.WriteByte(0x5A);
            fs.Position = 4_294_967_296L;
            fs.WriteByte(0xA5);
        }
        unsafe
        {
            using var map = new MappedFile(path, writable: false);
            Assert.That(map.Bytes, Is.EqualTo(5_368_709_120L));
            Assert.That(map.Pointer[5_000_000_001L], Is.EqualTo(0x5A));
            Assert.That(map.Pointer[4_294_967_296L], Is.EqualTo(0xA5));
            Assert.That(map.Pointer[4_294_967_295L], Is.EqualTo(0));
        }
        File.Delete(path);
    }

    // ---- the limits of the 32-bit positions

    [Test]
    public void EndToEnd_RefusesTextsFromTwoToThePowerThirtyTwoMinusThree_AndPassesTheLengthCheckBelow()
    {
        string dir = Run("lim");
        foreach (long n in new[] { uint.MaxValue - 2L, uint.MaxValue - 1L, uint.MaxValue, uint.MaxValue + 1L, 6_200_000_000L })
        {
            string input = Path.Combine(dir, "in.txt");
            using (var fs = new FileStream(input, FileMode.Create, FileAccess.Write)) fs.SetLength(n);
            var o = new EndToEndOptions { InputPath = input, OutPrefix = Path.Combine(dir, "out"), ScratchDir = Path.Combine(dir, "scratch") };
            Assert.Throws<NotSupportedException>(() => EndToEnd.Run(o), $"n = {n}");
            File.Delete(input);
        }
        // below the limit the run starts: stop it at the first phase, before a byte of the (sparse) text is read
        string big = Path.Combine(dir, "in2.txt");
        using (var fs = new FileStream(big, FileMode.Create, FileAccess.Write)) fs.SetLength(uint.MaxValue - 3L);
        var started = new List<string>();
        var below = new EndToEndOptions
        {
            InputPath = big,
            OutPrefix = Path.Combine(dir, "out2"),
            ScratchDir = Path.Combine(dir, "scratch2"),
            Trace = (what, phase) =>
            {
                started.Add(what + " " + phase);
                throw new OperationCanceledException("stop after the length check");
            },
        };
        Assert.Throws<OperationCanceledException>(() => EndToEnd.Run(below));
        Assert.That(started, Is.EqualTo(new[] { "begin ph0_prepare" }), "n = 2^32 - 4 passes the length check");
        File.Delete(big);
    }

    [Test]
    public void ChildScan_AcceptsTwoToThePowerThirtyTwoMinusThree_AndRefusesTwoToThePowerThirtyTwoMinusTwo()
    {
        string dir = Run("cs");
        var options = new ChildBuildOptions { ScratchDir = dir, Codec = ChildCodec.Packed };
        using (ChildScan.Create(uint.MaxValue - 3L, new ChildScanOutputs(), options)) { }
        Assert.Throws<NotSupportedException>(() => ChildScan.Create(uint.MaxValue - 1L, new ChildScanOutputs(), options));
        Assert.Throws<NotSupportedException>(() => ChildScan.Create(uint.MaxValue, new ChildScanOutputs(), options));
        Assert.Throws<NotSupportedException>(() => ChildScan.Create(2_000_000_000L, new ChildScanOutputs { ChildPath = Path.Combine(dir, "c.bin") }, options with { Codec = ChildCodec.Tagged }), "the tagged codec ends at 2^30 - 1 symbols");
    }

    private static unsafe void Plant(NativeBlock block, int bits, long start, byte[] symbols)
    {
        for (int i = 0; i < symbols.Length; i++)
        {
            for (int b = 0; b < bits; b++)
            {
                long bit = ((start + i) * bits) + b;
                byte mask = (byte)(0x80 >> (int)(bit % 8));
                if (((symbols[i] >> (bits - 1 - b)) & 1) != 0) block.Pointer[bit / 8] |= mask;
                else block.Pointer[bit / 8] &= (byte)~mask;
            }
        }
    }

    /// <summary>A view of another block's memory that does not free it (the rank table disposes its owner; the test owns the block).</summary>
    private sealed unsafe class NonOwningBlock(NativeBlock inner) : IMemoryBlock
    {
        public byte* Pointer => inner.Pointer;

        public long Bytes => inner.Bytes;

        public void Dispose()
        {
        }
    }
}
