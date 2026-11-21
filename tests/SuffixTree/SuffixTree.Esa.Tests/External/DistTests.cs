using System.Runtime.InteropServices;
using SuffixTree.Esa.External;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.External;

/// <summary>
/// DIST, the one-level radix distribute-and-scatter primitive (ADR-004 section 3.3): distribute (key, value) tuples to P = ceil(slots / R) chunk files by key / R through P write buffers, then scatter every chunk into an array
/// of R entries and write it through a window. The oracle is the plain in-RAM scatter.
/// </summary>
[TestFixture]
public unsafe class DistTests
{
    private static DistOptions Opt(TempDir d, long r, int io, int read, int window, ArenaBudget? budget = null, IoCounters? counters = null)
        => new() { ScratchDir = d.Sub("scratch"), ChunkEntries = r, IoBufferBytes = io, ReadBytes = read, WindowBytes = window, Budget = budget, Counters = counters, Arena = ArenaOptions.Guarded(0xCD) };

    private static uint[] IsaOracle(uint[] sa)
    {
        var isa = new uint[sa.Length];
        for (int i = 0; i < sa.Length; i++) isa[sa[i]] = (uint)i;
        return isa;
    }

    /// <summary>Phi[SA[i]] = SA[i-1] for i = 1..n (SA has n + 1 entries, SA[0] = n).</summary>
    private static uint[] PhiOracle(uint[] sa)
    {
        var phi = new uint[sa.Length - 1];
        for (int i = 1; i < sa.Length; i++) phi[sa[i]] = sa[i - 1];
        return phi;
    }

    /// <summary>A suffix array of the text with the empty suffix first (n + 1 entries).</summary>
    private static uint[] SaOf(string text) => EsaIndex.Build(text, ChildLayoutKind.Tagged).SuffixArray.Select(x => (uint)x).ToArray();

    /// <summary>Per chunk, the number of tuples that were flushed to the file given a buffer capacity (the remainder stays in memory).</summary>
    private static long ExpectedFlushedTuples(long slots, long r, int capacity)
    {
        long flushed = 0;
        for (long start = 0; start < slots; start += r) flushed += Math.Min(r, slots - start) / capacity * capacity;
        return flushed;
    }

    // ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // ISA: contract, exhaustive small cases, properties

    [Test]
    public void Isa_EveryPermutationOfFiveOrFewer_EveryConfiguration()
    {
        using var dir = new TempDir();
        string sa = dir.File("sa.bin");
        string isa = dir.File("isa.bin");
        int runs = 0;
        for (int count = 1; count <= 5; count++)
        {
            long n = count - 1;
            (long R, int Io, int Read, int Win)[] configs =
            [
                (1, 8, 4, 1), (2, 16, 8, 3), (3, 24, 12, 7), (count, 8, 4, 16), (count + 1, 1, 4, 1 << 20), (count - 1 < 1 ? 1 : count - 1, 64, 1 << 20, 5), (1000, 4096, 1 << 20, 1 << 20),
            ];
            foreach (uint[] perm in FileBits.Permutations(count))
            {
                FileBits.WriteUInts(sa, perm);
                uint[] expected = IsaOracle(perm);
                foreach ((long r, int io, int read, int win) in configs)
                {
                    DistStats stats = Dist.BuildIsa(sa, n, isa, Opt(dir, r, io, read, win));
                    Assert.That(FileBits.ReadUInts(isa), Is.EqualTo(expected), $"perm [{string.Join(",", perm)}] R={r} io={io} read={read} win={win}");
                    Assert.That(stats.Chunks, Is.EqualTo((count + r - 1) / r));
                    Assert.That(Directory.GetFileSystemEntries(dir.Sub("scratch")), Is.Empty);
                    runs++;
                }
            }
        }
        Assert.That(runs, Is.GreaterThan(1000));
    }

    [Test]
    public void Isa_RandomPermutations_AllOptions_ExactAccounting()
    {
        var rng = new SplitMix64(2026);
        for (int round = 0; round < 60; round++)
        {
            using var dir = new TempDir();
            int n = rng.NextInt(3000);
            long rows = n + 1L;
            uint[] sa = Enumerable.Range(0, n + 1).Select(i => (uint)i).OrderBy(_ => rng.Next()).ToArray();
            string saPath = dir.File("sa.bin");
            FileBits.WriteUInts(saPath, sa);
            long r = 1 + rng.NextInt(n + 3);
            int io = 8 * (1 + rng.NextInt(40));
            int read = 4 * (1 + rng.NextInt(3000));
            int window = 1 + rng.NextInt(5000);
            var counters = new IoCounters();
            var budget = new ArenaBudget(long.MaxValue);
            DistStats stats = Dist.BuildIsa(saPath, n, dir.File("isa.bin"), Opt(dir, r, io, read, window, budget, counters));
            Assert.That(FileBits.ReadUInts(dir.File("isa.bin")), Is.EqualTo(IsaOracle(sa)), $"round {round}: n={n} R={r} io={io} read={read} window={window}");

            long chunks = (rows + r - 1) / r;
            int capacity = io / 8;
            long flushedBytes = ExpectedFlushedTuples(rows, r, capacity) * 8;
            Assert.That(stats, Is.EqualTo(new DistStats(chunks, flushedBytes, flushedBytes, stats.Reservation)), $"round {round}: the tuples that outgrew a buffer are written and read once");
            Assert.That(counters.WriteBytes, Is.EqualTo(flushedBytes));
            Assert.That(counters.ReadBytes, Is.EqualTo((4 * rows) + flushedBytes), "the source once, the chunk files once");
            Assert.That(counters.MappedWriteBytes, Is.EqualTo(4 * rows), "the target through windows");
            Assert.That(counters.MappedMaxWindowBytes, Is.EqualTo(Math.Min(window, 4 * rows)));
            Assert.That(counters.ScratchLiveBytes, Is.Zero, "every chunk file is removed");
            Assert.That(counters.ScratchPeakBytes, Is.EqualTo(flushedBytes), "all the chunk files exist together at the end of pass 1");
            Assert.That(budget.LiveBytes, Is.Zero, "the arena gives everything back");
            Assert.That(Directory.GetFileSystemEntries(dir.Sub("scratch")), Is.Empty);
            Assert.That(new FileInfo(dir.File("isa.bin")).Length, Is.EqualTo(4 * rows));
        }
    }

    [Test]
    public void Isa_Boundaries_OfTheChunkSize_AndTheBufferCapacity()
    {
        using var dir = new TempDir();
        foreach (long r in new long[] { 1, 4, 64 })
        {
            foreach (int cap in new[] { 1, 2, 3, 5 })
            {
                foreach (long rows in new[] { 1, 2, r - 1, r, r + 1, (2 * r) - 1, 2 * r, (2 * r) + 1, (3 * r * cap) + 1 })
                {
                    if (rows < 1) continue;
                    int n = (int)(rows - 1);
                    foreach (bool reversed in new[] { false, true })
                    {
                        uint[] sa = Enumerable.Range(0, n + 1).Select(i => reversed ? (uint)(n - i) : (uint)i).ToArray();
                        FileBits.WriteUInts(dir.File("sa.bin"), sa);
                        DistStats s = Dist.BuildIsa(dir.File("sa.bin"), n, dir.File("isa.bin"), Opt(dir, r, 8 * cap, 1 << 20, 1 << 20));
                        Assert.That(FileBits.ReadUInts(dir.File("isa.bin")), Is.EqualTo(IsaOracle(sa)), $"R={r} cap={cap} rows={rows} reversed={reversed}");
                        Assert.That(s.Chunks, Is.EqualTo((rows + r - 1) / r));
                        Assert.That(s.TupleBytesWritten, Is.EqualTo(ExpectedFlushedTuples(rows, r, cap) * 8));
                    }
                }
            }
        }
    }

    [Test]
    public void Isa_OfTheSuffixArraysOfTheAdversarialFamilies_IsTheInversePermutation()
    {
        using var dir = new TempDir();
        var texts = new List<(string Name, string Text)>
        {
            ("a^n", StringFamilies.AllSame(3000)),
            ("fibonacci", StringFamilies.Fibonacci(3000)),
            ("thue-morse", StringFamilies.ThueMorse(3000)),
            ("dna", StringFamilies.RandomDna(3000, 5)),
            ("periodic", StringFamilies.Periodic(3000, 7, 3)),
            ("banana", "banana"),
            ("single", "x"),
        };
        foreach ((string name, string text) in texts)
        {
            uint[] sa = SaOf(text);
            FileBits.WriteUInts(dir.File("sa.bin"), sa);
            foreach (long r in new long[] { 7, 500, 100_000 })
            {
                Dist.BuildIsa(dir.File("sa.bin"), text.Length, dir.File("isa.bin"), Opt(dir, r, 64, 100, 64));
                uint[] isa = FileBits.ReadUInts(dir.File("isa.bin"));
                Assert.That(isa, Is.EqualTo(IsaOracle(sa)), $"{name} R={r}");
                for (int i = 0; i < sa.Length; i++) Assert.That(sa[isa[i]], Is.EqualTo((uint)i), "ISA is the inverse of SA");
            }
        }
    }

    [Test]
    public void Isa_Empty_IsOneEntry()
    {
        using var dir = new TempDir();
        FileBits.WriteUInts(dir.File("sa.bin"), [0u]);
        DistStats s = Dist.BuildIsa(dir.File("sa.bin"), 0, dir.File("isa.bin"), Opt(dir, 8, 64, 64, 64));
        Assert.That(FileBits.ReadUInts(dir.File("isa.bin")), Is.EqualTo(new uint[] { 0 }));
        Assert.That(s, Is.EqualTo(new DistStats(1, 0, 0, s.Reservation)), "one tuple stays in its buffer");
    }

    // ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // Phi (the second user of the same primitive)

    [Test]
    public void Phi_EqualsThePredecessorTable_ForRealSuffixArrays_AndExhaustiveSmallTexts()
    {
        using var dir = new TempDir();
        var texts = new List<string> { "", "a", "banana", "mississippi", StringFamilies.AllSame(500), StringFamilies.Fibonacci(1000), StringFamilies.ThueMorse(1000), StringFamilies.RandomDna(2000, 9) };
        texts.AddRange(ExhaustiveStrings.All(2, 7));
        int checkedTexts = 0;
        foreach (string text in texts)
        {
            uint[] sa = SaOf(text);
            FileBits.WriteUInts(dir.File("sa.bin"), sa);
            long r = text.Length < 20 ? 3 : 333;
            Dist.BuildPhi(dir.File("sa.bin"), text.Length, dir.File("phi.bin"), Opt(dir, r, 40, 28, 12));
            Assert.That(FileBits.ReadUInts(dir.File("phi.bin")), Is.EqualTo(PhiOracle(sa)), $"'{(text.Length > 20 ? text[..20] + "..." : text)}'");
            checkedTexts++;
        }
        Assert.That(checkedTexts, Is.GreaterThan(250));
    }

    [Test]
    public void Phi_FeedsThePlcpLoop_SoItsValuesAreWhatLcpNeeds()
    {
        // LCP[i] = lcp(T[SA[i]..], T[SA[i-1]..]) = lcp(T[j..], T[Phi[j]..]) for j = SA[i]: recompute the LCP array from Phi and compare with the kernel's
        string text = StringFamilies.RandomDna(5000, 31) + StringFamilies.Fibonacci(2000);
        using var dir = new TempDir();
        EsaIndex index = EsaIndex.Build(text, ChildLayoutKind.Tagged);
        int[] lcp = index.Lcp.ToArray();
        uint[] sa = index.SuffixArray.Select(x => (uint)x).ToArray();
        FileBits.WriteUInts(dir.File("sa.bin"), sa);
        Dist.BuildPhi(dir.File("sa.bin"), text.Length, dir.File("phi.bin"), Opt(dir, 700, 64, 1 << 12, 1 << 12));
        uint[] phi = FileBits.ReadUInts(dir.File("phi.bin"));
        for (int i = 1; i < sa.Length; i++)
        {
            int j = (int)sa[i];
            int p = (int)phi[j];
            int l = 0;
            while (j + l < text.Length && p + l < text.Length && text[j + l] == text[p + l]) l++;
            Assert.That(lcp[i], Is.EqualTo(l), $"rank {i}");
        }
    }

    // ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // the generalised target

    private struct EveryThird : IDistMapper<uint>
    {
        public bool TryMap(long index, uint sourceValue, out uint key, out uint value)
        {
            key = sourceValue;
            value = (uint)(index * 10);
            return index % 3 != 0;
        }
    }

    [Test]
    public void Scatter_WithASkippingMapper_LeavesTheFillInTheUntouchedSlots_AndNeedsNoCheck()
    {
        using var dir = new TempDir();
        var keys = new uint[] { 9, 2, 4, 0, 7, 5, 3, 8 }; // distinct keys below 10
        FileBits.WriteUInts(dir.File("src.bin"), keys);
        var mapper = new EveryThird();
        var target = new DistTarget<uint>(10, 0xFFFF_FFFEu, false);
        Dist.Scatter(dir.File("src.bin"), keys.Length, ref mapper, dir.File("out.bin"), target, Opt(dir, 4, 24, 12, 8));
        var expected = Enumerable.Repeat(0xFFFF_FFFEu, 10).ToArray();
        for (int i = 0; i < keys.Length; i++)
        {
            if (i % 3 != 0) expected[keys[i]] = (uint)(i * 10);
        }
        Assert.That(FileBits.ReadUInts(dir.File("out.bin")), Is.EqualTo(expected));
    }

    private struct Widen<TV> : IDistMapper<TV>
        where TV : unmanaged
    {
        private readonly Func<long, TV> _make;

        internal Widen(Func<long, TV> make) => _make = make;

        public bool TryMap(long index, uint sourceValue, out uint key, out TV value)
        {
            key = sourceValue;
            value = _make(index);
            return true;
        }
    }

    [Test]
    public void Scatter_IsGenericInTheValue_ByteUShortULong()
    {
        using var dir = new TempDir();
        var rng = new SplitMix64(8);
        const int n = 700;
        uint[] perm = Enumerable.Range(0, n).Select(i => (uint)i).OrderBy(_ => rng.Next()).ToArray();
        FileBits.WriteUInts(dir.File("src.bin"), perm);

        var mb = new Widen<byte>(i => (byte)(i % 251));
        DistStats statsB = Dist.Scatter(dir.File("src.bin"), n, ref mb, dir.File("b.bin"), new DistTarget<byte>(n, 255, true), Opt(dir, 64, 30, 40, 100));
        Assert.That(statsB.TupleBytesWritten, Is.EqualTo(ExpectedFlushedTuples(n, 64, 30 / 5) * 5), "a tuple of a byte is 5 bytes");
        byte[] gotB = File.ReadAllBytes(dir.File("b.bin"));
        Assert.That(gotB, Has.Length.EqualTo(n));
        for (int i = 0; i < n; i++) Assert.That(gotB[perm[i]], Is.EqualTo((byte)(i % 251)));

        var mu = new Widen<ushort>(i => (ushort)(i * 3));
        Dist.Scatter(dir.File("src.bin"), n, ref mu, dir.File("u.bin"), new DistTarget<ushort>(n, ushort.MaxValue, true), Opt(dir, 100, 50, 64, 33));
        ushort[] gotU = MemoryMarshal.Cast<byte, ushort>(File.ReadAllBytes(dir.File("u.bin"))).ToArray();
        for (int i = 0; i < n; i++) Assert.That(gotU[perm[i]], Is.EqualTo((ushort)(i * 3)));

        var ml = new Widen<ulong>(i => 0x1_0000_0000UL + (ulong)i);
        DistStats statsL = Dist.Scatter(dir.File("src.bin"), n, ref ml, dir.File("l.bin"), new DistTarget<ulong>(n, ulong.MaxValue, true), Opt(dir, 128, 36, 400, 1000));
        Assert.That(statsL.TupleBytesWritten, Is.EqualTo(ExpectedFlushedTuples(n, 128, 36 / 12) * 12), "a tuple of a 64-bit value is 12 bytes");
        Assert.That(statsL.TupleBytesRead, Is.EqualTo(statsL.TupleBytesWritten));
        ulong[] gotL = MemoryMarshal.Cast<byte, ulong>(File.ReadAllBytes(dir.File("l.bin"))).ToArray();
        for (int i = 0; i < n; i++) Assert.That(gotL[perm[i]], Is.EqualTo(0x1_0000_0000UL + (ulong)i));
    }

    [Test]
    public void TheTupleLayout_IsPacked_FourPlusTheValue()
    {
        Assert.That(sizeof(DistEntry<uint>), Is.EqualTo(8));
        Assert.That(sizeof(DistEntry<ulong>), Is.EqualTo(12));
        Assert.That(sizeof(DistEntry<byte>), Is.EqualTo(5));
        var e = new DistEntry<ushort>(0xAABBCCDD, 0x1234);
        Assert.That((e.Key, e.Value), Is.EqualTo((0xAABBCCDDu, (ushort)0x1234)));
    }

    // ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // arithmetic

    [Test]
    public void ChunksFor_IsTheCeilingOfTheDivision()
    {
        Assert.That(Dist.ChunksFor(0, 5), Is.Zero);
        Assert.That(Dist.ChunksFor(1, 5), Is.EqualTo(1));
        Assert.That(Dist.ChunksFor(5, 5), Is.EqualTo(1));
        Assert.That(Dist.ChunksFor(6, 5), Is.EqualTo(2));
        Assert.That(Dist.ChunksFor(10, 5), Is.EqualTo(2));
        Assert.That(Dist.ChunksFor(11, 5), Is.EqualTo(3));
        Assert.That(Dist.ChunksFor(3_000_000_001L, 8_000_000), Is.EqualTo(376));
        Assert.That(Dist.ChunksFor(3_000_000_000L, 8_000_000), Is.EqualTo(375), "the P of the ADR at 3x10^9");
        Assert.That(Dist.ChunksFor(long.MaxValue, 1), Is.EqualTo(long.MaxValue), "no overflow of the rounding");
        Assert.That(Dist.ChunksFor(long.MaxValue, long.MaxValue), Is.EqualTo(1));
        Assert.Throws<ArgumentOutOfRangeException>(() => Dist.ChunksFor(1, 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => Dist.ChunksFor(-1, 5));
    }

    [Test]
    public void ChunkDivider_EqualsTheIntegerDivision_ForEveryKindOfDivisor_AndTheEdgesOfTheKeyRange()
    {
        var rng = new SplitMix64(99);
        long[] divisors = [1, 2, 3, 5, 7, 8, 10, 255, 256, 257, 65_535, 65_536, 266_667, 8_000_000, int.MaxValue, 1L << 31, uint.MaxValue - 1L, uint.MaxValue, 1L << 32, (1L << 32) + 5, long.MaxValue];
        foreach (long r in divisors)
        {
            var d = new ChunkDivider(r);
            var keys = new List<uint> { 0, 1, uint.MaxValue, uint.MaxValue - 1 };
            foreach (long m in new[] { r - 1, r, r + 1, (2 * r) - 1, 2 * r, (2 * r) + 1, 3 * r }.Where(m => m is >= 0 and <= uint.MaxValue)) keys.Add((uint)m);
            for (int i = 0; i < 20_000; i++) keys.Add((uint)rng.Next());
            foreach (uint key in keys) Assert.That(d.Of(key), Is.EqualTo(r > uint.MaxValue ? 0u : (uint)(key / r)), $"R={r} key={key}");
        }
        Assert.Throws<ArgumentOutOfRangeException>(() => _ = new ChunkDivider(0));
        Assert.Throws<ArgumentOutOfRangeException>(() => _ = new ChunkDivider(-5));
    }

    // ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // memory

    [Test]
    public void TheArena_PeaksAtTheDocumentedBlocks_AndTheCapIsEnforcedAtTheBoundary()
    {
        using var dir = new TempDir();
        const int n = 999; // rows 1000
        var rng = new SplitMix64(3);
        uint[] sa = Enumerable.Range(0, n + 1).Select(i => (uint)i).OrderBy(_ => rng.Next()).ToArray();
        FileBits.WriteUInts(dir.File("sa.bin"), sa);
        // R = 300 -> 4 chunks; io 64 B -> 8 tuples per buffer; read 400 B -> 100 source entries; chunk array 300 x 4 B; chunk read 8 x 8 B
        long pass1 = (4L * 8 * 8) + (100 * 4);
        long pass2 = (4L * 8 * 8) + (300 * 4) + (8 * 8);
        long expectedPeak = Math.Max(pass1, pass2);
        var budget = new ArenaBudget(long.MaxValue);
        Dist.BuildIsa(dir.File("sa.bin"), n, dir.File("isa.bin"), Opt(dir, 300, 64, 400, 1 << 20, budget));
        Assert.That(budget.PeakBytes, Is.EqualTo(expectedPeak));
        Assert.That(budget.LiveBytes, Is.Zero);

        var exact = new ArenaBudget(expectedPeak);
        Assert.DoesNotThrow(() => Dist.BuildIsa(dir.File("sa.bin"), n, dir.File("isa.bin"), Opt(dir, 300, 64, 400, 1 << 20, exact)));
        Assert.That(FileBits.ReadUInts(dir.File("isa.bin")), Is.EqualTo(IsaOracle(sa)));
        var tight = new ArenaBudget(expectedPeak - 1);
        Assert.Throws<ArenaBudgetExceededException>(() => Dist.BuildIsa(dir.File("sa.bin"), n, dir.File("isa.bin"), Opt(dir, 300, 64, 400, 1 << 20, tight)));
        Assert.That(tight.LiveBytes, Is.Zero, "a refused build leaves nothing charged");
        Assert.That(Directory.GetFileSystemEntries(dir.Sub("scratch")), Is.Empty, "and nothing on disk");
        var tiny = new ArenaBudget(10);
        Assert.Throws<ArenaBudgetExceededException>(() => Dist.BuildIsa(dir.File("sa.bin"), n, dir.File("isa.bin"), Opt(dir, 300, 64, 400, 1 << 20, tiny)));
    }

    [Test]
    public void Phi_RefusesASuffixArrayWhoseFirstEntryIsNotN_OrWhoseLaterEntryIsNOrMore()
    {
        using var dir = new TempDir();
        var budget = new ArenaBudget(long.MaxValue);
        var counters = new IoCounters();
        foreach (uint[] sa in new[] { new uint[] { 2, 0, 1, 3 }, new uint[] { 3, 0, 3, 1 }, new uint[] { 3, 0, 1, 9 }, new uint[] { 2, 0, 1, 2 }, new uint[] { 0, 0, 1, 2 } })
        {
            FileBits.WriteUInts(dir.File("sa.bin"), sa);
            Assert.Throws<InvalidDataException>(() => Dist.BuildPhi(dir.File("sa.bin"), 3, dir.File("phi.bin"), Opt(dir, 2, 8, 8, 8, budget, counters)), string.Join(",", sa));
            AssertClean(dir, budget, counters);
        }
        // the first entry is the only thing wrong with [2, 0, 1, 2]: the other three are positions below n
        InvalidDataException first = Assert.Throws<InvalidDataException>(() => { FileBits.WriteUInts(dir.File("f.bin"), [2u, 0, 1, 2]); Dist.BuildPhi(dir.File("f.bin"), 3, dir.File("phi.bin"), Opt(dir, 2, 8, 8, 8)); })!;
        Assert.That(first.Message, Does.Contain("SA[0]"));
        FileBits.WriteUInts(dir.File("ok.bin"), [3u, 2, 0, 1]);
        Dist.BuildPhi(dir.File("ok.bin"), 3, dir.File("phi.bin"), Opt(dir, 2, 8, 8, 8));
        Assert.That(FileBits.ReadUInts(dir.File("phi.bin")), Is.EqualTo(new uint[] { 2, 0, 3 }), "Phi[SA[i]] = SA[i-1]: Phi[0] = 2, Phi[1] = 0, Phi[2] = 3");
    }

    [Test]
    public void TheLimitOfThe32BitPositions_IsExactlyAtTwoToThe32MinusTwo()
    {
        using var dir = new TempDir();
        FileBits.WriteUInts(dir.File("sa.bin"), [0u]);
        foreach (long n in new[] { uint.MaxValue - 2L, uint.MaxValue - 1L, uint.MaxValue, uint.MaxValue + 7L })
        {
            Assert.Throws<NotSupportedException>(() => Dist.BuildIsa(dir.File("sa.bin"), n, dir.File("isa.bin"), Opt(dir, 4, 8, 8, 8)), $"ISA n={n}");
            Assert.Throws<NotSupportedException>(() => Dist.BuildPhi(dir.File("sa.bin"), n, dir.File("phi.bin"), Opt(dir, 4, 8, 8, 8)), $"Phi n={n}");
        }
        // one below the limit the call is accepted and fails on the length of the source file instead
        Assert.Throws<InvalidDataException>(() => Dist.BuildIsa(dir.File("sa.bin"), uint.MaxValue - 3L, dir.File("isa.bin"), Opt(dir, 4, 8, 8, 8)));
        Assert.Throws<InvalidDataException>(() => Dist.BuildPhi(dir.File("sa.bin"), uint.MaxValue - 3L, dir.File("phi.bin"), Opt(dir, 4, 8, 8, 8)));
    }

    [Test]
    public void EnsureKeysFit_AllowsUpToTwoToThe32Slots()
    {
        Assert.DoesNotThrow(() => Dist.EnsureKeysFit(0));
        Assert.DoesNotThrow(() => Dist.EnsureKeysFit(1L << 32));
        Assert.Throws<NotSupportedException>(() => Dist.EnsureKeysFit((1L << 32) + 1));
        Assert.Throws<NotSupportedException>(() => Dist.EnsureKeysFit(long.MaxValue));
    }

    [Test]
    public void ScatterInto_PlacesTuplesBySlot_AndRefusesADamagedChunk()
    {
        var checkedTarget = new DistTarget<uint>(20, uint.MaxValue, true);
        var plainTarget = new DistTarget<uint>(20, uint.MaxValue, false);
        uint[] array = Enumerable.Repeat(uint.MaxValue, 5).ToArray(); // the chunk of the slots 5..9
        Dist.ScatterInto<uint>([new DistEntry<uint>(5, 100), new DistEntry<uint>(9, 200), new DistEntry<uint>(7, 300)], array, 5, checkedTarget);
        Assert.That(array, Is.EqualTo(new uint[] { 100, uint.MaxValue, 300, uint.MaxValue, 200 }));
        Assert.Throws<InvalidDataException>(() => Dist.ScatterInto<uint>([new DistEntry<uint>(4, 1)], array, 5, checkedTarget), "a key below the chunk");
        Assert.Throws<InvalidDataException>(() => Dist.ScatterInto<uint>([new DistEntry<uint>(10, 1)], array, 5, checkedTarget), "a key just above the chunk");
        Assert.Throws<InvalidDataException>(() => Dist.ScatterInto<uint>([new DistEntry<uint>(uint.MaxValue, 1)], array, 5, checkedTarget), "a key far above");
        Assert.Throws<InvalidDataException>(() => Dist.ScatterInto<uint>([new DistEntry<uint>(5, 1)], array, 5, checkedTarget), "slot 5 was written: a repeat");
        Assert.Throws<InvalidDataException>(() => Dist.ScatterInto<uint>([new DistEntry<uint>(6, 1), new DistEntry<uint>(6, 2)], array, 5, checkedTarget), "a repeat inside one run");
        array[1] = uint.MaxValue;
        Dist.ScatterInto<uint>([new DistEntry<uint>(5, 111), new DistEntry<uint>(5, 222)], array, 5, plainTarget);
        Assert.That(array[0], Is.EqualTo(222u), "without the check the last tuple wins");
        Assert.Throws<InvalidDataException>(() => Dist.ScatterInto<uint>([new DistEntry<uint>(4, 1)], array, 5, plainTarget), "a damaged chunk is refused with or without the check");
        Dist.ScatterInto<uint>([], array, 5, checkedTarget); // an empty run is nothing
    }

    [Test]
    public void ScratchIsReleasedChunkByChunk_WhilePassTwoRuns_NotAtTheEnd()
    {
        // ADR-004 Q1: a chunk file is deleted as soon as its chunk is in memory. A poller records (live scratch bytes, bytes written to the target) in order: after the peak (the end of pass 1, all the chunk files exist)
        // the scratch falls in steps while the target is still being written. (Disposing the scratch at the end of the call also deletes the files, but only after the last chunk was written.)
        using var dir = new TempDir();
        const int n = 1_500_000;
        var sa = new uint[n + 1];
        for (int i = 0; i <= n; i++) sa[i] = (uint)((i * 7919L) % (n + 1));
        FileBits.WriteUInts(dir.File("sa.bin"), sa);
        var counters = new IoCounters();
        var series = new List<(long Scratch, long Mapped)>();
        using var stop = new ManualResetEventSlim(false);
        var poller = new Thread(() =>
        {
            (long, long) last = (-1, -1);
            while (!stop.IsSet)
            {
                (long, long) now = (counters.ScratchLiveBytes, counters.MappedWriteBytes);
                if (now == last) continue;
                series.Add(now);
                last = now;
            }
        });
        poller.Start();
        DistStats stats = Dist.BuildIsa(dir.File("sa.bin"), n, dir.File("isa.bin"), Opt(dir, 5_000, 4096, 1 << 20, 1 << 20, null, counters));
        stop.Set();
        poller.Join();
        long peak = stats.TupleBytesWritten;
        long targetBytes = 4L * (n + 1);
        Assert.That(peak, Is.GreaterThan(0));
        int atPeak = series.FindIndex(x => x.Scratch == peak);
        Assert.That(atPeak, Is.GreaterThanOrEqualTo(0), "the poller saw the peak");
        long fallingWhileWriting = series.Skip(atPeak + 1).Count(x => x.Scratch > 0 && x.Scratch < peak && x.Mapped < targetBytes);
        Assert.That(fallingWhileWriting, Is.GreaterThanOrEqualTo(3), $"the live scratch bytes fell in steps while the target was still being written (series of {series.Count} values)");
    }

    [Test]
    public void TheDefaultOptions_AreTheAdrValues_AndTheWindowIsAtLeastOneMiB()
    {
        var o = new DistOptions();
        Assert.That(o.ChunkEntries, Is.EqualTo(8_000_000));
        Assert.That(o.IoBufferBytes, Is.EqualTo(64 * 1024));
        Assert.That(o.ReadBytes, Is.EqualTo(1 << 20));
        Assert.That(o.WindowBytes, Is.GreaterThanOrEqualTo(1 << 20));
        Assert.That(o.Budget, Is.Null);
        Assert.That(o.Counters, Is.Null);
    }

    [Test]
    public void WithDefaultOptionsExceptTheScratch_TheTargetIsWrittenThroughWindowsOfOneMiB()
    {
        using var dir = new TempDir();
        const int n = 700_000; // 2.8 MB of ISA
        uint[] sa = Enumerable.Range(0, n + 1).Select(i => (uint)((i * 7919L) % (n + 1))).ToArray();
        FileBits.WriteUInts(dir.File("sa.bin"), sa);
        var c = new IoCounters();
        Dist.BuildIsa(dir.File("sa.bin"), n, dir.File("isa.bin"), new DistOptions { ScratchDir = dir.Sub("scratch"), Counters = c });
        Assert.That(FileBits.ReadUInts(dir.File("isa.bin")), Is.EqualTo(IsaOracle(sa)));
        Assert.That(c.MappedMaxWindowBytes, Is.EqualTo(1 << 20));
    }

    // ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // failure modes: never a hang, never silent damage, nothing left behind

    private static void AssertClean(TempDir dir, ArenaBudget budget, IoCounters counters)
    {
        Assert.That(budget.LiveBytes, Is.Zero, "the arena was disposed");
        Assert.That(counters.ScratchLiveBytes, Is.Zero, "the scratch accounting is back at zero");
        if (Directory.Exists(dir.Sub("scratch"))) Assert.That(Directory.GetFileSystemEntries(dir.Sub("scratch")), Is.Empty, "no chunk file is left");
    }

    [Test]
    public void NotAPermutation_ARepeat_AMiss_AnOutOfRangeKey_AreRefused_AndCleanedUp()
    {
        using var dir = new TempDir();
        foreach (long r in new long[] { 1, 2, 3, 100 })
        {
            foreach (int io in new[] { 8, 24, 4096 })
            {
                (string Name, uint[] Sa)[] cases =
                [
                    ("repeat", [3, 1, 1, 2, 0]),
                    ("repeat-first", [0, 0, 2, 3, 4]),
                    ("repeat-last", [0, 1, 2, 3, 3]),
                    ("beyond-n", [0, 1, 2, 3, 9]),
                    ("beyond-n-by-one", [0, 1, 2, 3, 5]),
                    ("huge", [0, 1, 2, 3, uint.MaxValue]),
                ];
                foreach ((string name, uint[] sa) in cases)
                {
                    FileBits.WriteUInts(dir.File("sa.bin"), sa);
                    var budget = new ArenaBudget(long.MaxValue);
                    var counters = new IoCounters();
                    Assert.Throws<InvalidDataException>(() => Dist.BuildIsa(dir.File("sa.bin"), 4, dir.File("isa.bin"), Opt(dir, r, io, 8, 8, budget, counters)), $"{name} R={r} io={io}");
                    AssertClean(dir, budget, counters);
                }
            }
        }
    }

    [Test]
    public void AMissingPosition_IsReported_WithItsNumber()
    {
        using var dir = new TempDir();
        FileBits.WriteUInts(dir.File("sa.bin"), [0u, 1, 1, 3]); // 2 never occurs, 1 occurs twice: whichever is found first is reported
        InvalidDataException ex = Assert.Throws<InvalidDataException>(() => Dist.BuildIsa(dir.File("sa.bin"), 3, dir.File("isa.bin"), Opt(dir, 100, 64, 64, 64)))!;
        Assert.That(ex.Message, Does.Contain("1").Or.Contain("2"));
        FileBits.WriteUInts(dir.File("sa2.bin"), [0u, 1, 3, 3]); // 2 misses, 3 repeats
        Assert.Throws<InvalidDataException>(() => Dist.BuildIsa(dir.File("sa2.bin"), 3, dir.File("isa.bin"), Opt(dir, 2, 8, 64, 64)));
        // a pure miss cannot be built from a file with n + 1 entries all below n + 1 without a repeat, so the miss path is exercised through a skipping mapper below
    }

    private struct SkipAll : IDistMapper<uint>
    {
        public bool TryMap(long index, uint sourceValue, out uint key, out uint value)
        {
            key = 0;
            value = 0;
            return index == 0;
        }
    }

    [Test]
    public void ExactlyOnce_ReportsAMissingSlot()
    {
        using var dir = new TempDir();
        FileBits.WriteUInts(dir.File("src.bin"), [0u, 1, 2]);
        var mapper = new SkipAll();
        var budget = new ArenaBudget(long.MaxValue);
        var counters = new IoCounters();
        InvalidDataException ex = Assert.Throws<InvalidDataException>(() => Dist.Scatter(dir.File("src.bin"), 3, ref mapper, dir.File("out.bin"), new DistTarget<uint>(3, uint.MaxValue, true), Opt(dir, 2, 8, 8, 8, budget, counters)))!;
        Assert.That(ex.Message, Does.Contain("1"), "slot 1 is the first never written (slot 0 was)");
        AssertClean(dir, budget, counters);
        // the same input without the check is fine
        Dist.Scatter(dir.File("src.bin"), 3, ref mapper, dir.File("out.bin"), new DistTarget<uint>(3, 7u, false), Opt(dir, 2, 8, 8, 8));
        Assert.That(FileBits.ReadUInts(dir.File("out.bin")), Is.EqualTo(new uint[] { 0, 7, 7 }));
    }

    private struct SkipKey(uint skip) : IDistMapper<uint>
    {
        public bool TryMap(long index, uint sourceValue, out uint key, out uint value)
        {
            key = sourceValue;
            value = (uint)index;
            return sourceValue != skip;
        }
    }

    [Test]
    public void ExactlyOnce_ReportsAMissingSlot_AtTheFirstAndLastSlotOfEveryChunk()
    {
        using var dir = new TempDir();
        const int n = 23; // slots 0..22, R = 5: chunks start at 0, 5, 10, 15, 20
        uint[] perm = Enumerable.Range(0, n).Select(i => (uint)((i * 7) % n)).ToArray();
        FileBits.WriteUInts(dir.File("src.bin"), perm);
        foreach (uint missing in new uint[] { 0, 4, 5, 9, 10, 15, 19, 20, 22 })
        {
            var mapper = new SkipKey(missing);
            var budget = new ArenaBudget(long.MaxValue);
            var counters = new IoCounters();
            InvalidDataException ex = Assert.Throws<InvalidDataException>(() => Dist.Scatter(dir.File("src.bin"), n, ref mapper, dir.File("out.bin"), new DistTarget<uint>(n, uint.MaxValue, true), Opt(dir, 5, 16, 16, 16, budget, counters)), $"slot {missing}")!;
            Assert.That(ex.Message, Does.Contain($"slot {missing} is never written"));
            AssertClean(dir, budget, counters);
        }
    }

    private struct FixedKey(uint key) : IDistMapper<uint>
    {
        public bool TryMap(long index, uint sourceValue, out uint key1, out uint value)
        {
            key1 = key;
            value = 1;
            return true;
        }
    }

    [Test]
    public void AKeyAtOrBeyondTheSlotCount_IsRefused_BeforeItCanLandInAWrongBuffer()
    {
        using var dir = new TempDir();
        FileBits.WriteUInts(dir.File("src.bin"), [0u]);
        foreach (uint key in new uint[] { 5, 6, 1000, uint.MaxValue })
        {
            var mapper = new FixedKey(key);
            var budget = new ArenaBudget(long.MaxValue);
            var counters = new IoCounters();
            InvalidDataException ex = Assert.Throws<InvalidDataException>(() => Dist.Scatter(dir.File("src.bin"), 1, ref mapper, dir.File("out.bin"), new DistTarget<uint>(5, 0, false), Opt(dir, 2, 8, 8, 8, budget, counters)), $"key {key}")!;
            Assert.That(ex.Message, Does.Contain("maps to key").And.Contain("5 slots"), $"key {key}: refused at the mapper, before it can reach a chunk");
            AssertClean(dir, budget, counters);
        }
        var last = new FixedKey(4);
        Dist.Scatter(dir.File("src.bin"), 1, ref last, dir.File("out.bin"), new DistTarget<uint>(5, 0, false), Opt(dir, 2, 8, 8, 8));
        Assert.That(FileBits.ReadUInts(dir.File("out.bin")), Is.EqualTo(new uint[] { 0, 0, 0, 0, 1 }), "slots - 1 is the last valid key");
    }

    private struct Boom : IDistMapper<uint>
    {
        public bool TryMap(long index, uint sourceValue, out uint key, out uint value)
        {
            if (index == 5) throw new NotSupportedException("mapper failed");
            key = sourceValue;
            value = 0;
            return true;
        }
    }

    [Test]
    public void AnExceptionOfTheMapper_PropagatesAfterTheCleanUp()
    {
        using var dir = new TempDir();
        FileBits.WriteUInts(dir.File("src.bin"), Enumerable.Range(0, 20).Select(i => (uint)i).ToArray());
        var mapper = new Boom();
        var budget = new ArenaBudget(long.MaxValue);
        var counters = new IoCounters();
        Assert.Throws<NotSupportedException>(() => Dist.Scatter(dir.File("src.bin"), 20, ref mapper, dir.File("out.bin"), new DistTarget<uint>(20, 0, false), Opt(dir, 4, 8, 16, 8, budget, counters)));
        AssertClean(dir, budget, counters);
    }

    [Test]
    public void ASourceOfTheWrongLength_IsRefused_ShortLongOrMissing()
    {
        using var dir = new TempDir();
        FileBits.WriteUInts(dir.File("sa.bin"), [0u, 1, 2, 3]);
        Assert.Throws<InvalidDataException>(() => Dist.BuildIsa(dir.File("sa.bin"), 4, dir.File("isa.bin"), Opt(dir, 4, 8, 8, 8)), "n = 4 needs 5 entries");
        Assert.Throws<InvalidDataException>(() => Dist.BuildIsa(dir.File("sa.bin"), 2, dir.File("isa.bin"), Opt(dir, 4, 8, 8, 8)), "n = 2 needs 3 entries");
        Assert.Throws<InvalidDataException>(() => Dist.BuildPhi(dir.File("sa.bin"), 4, dir.File("phi.bin"), Opt(dir, 4, 8, 8, 8)));
        Assert.Throws<FileNotFoundException>(() => Dist.BuildIsa(dir.File("nothing.bin"), 3, dir.File("isa.bin"), Opt(dir, 4, 8, 8, 8)));
        File.WriteAllBytes(dir.File("odd.bin"), new byte[15]);
        Assert.Throws<InvalidDataException>(() => Dist.BuildIsa(dir.File("odd.bin"), 3, dir.File("isa.bin"), Opt(dir, 4, 8, 8, 8)), "not a multiple of 4");
    }

    [Test]
    public void BadArguments_AreRefused()
    {
        using var dir = new TempDir();
        FileBits.WriteUInts(dir.File("sa.bin"), [0u]);
        Assert.Throws<ArgumentNullException>(() => Dist.BuildIsa(dir.File("sa.bin"), 0, dir.File("isa.bin"), null!));
        Assert.Throws<ArgumentOutOfRangeException>(() => Dist.BuildIsa(dir.File("sa.bin"), -1, dir.File("isa.bin"), Opt(dir, 4, 8, 8, 8)));
        Assert.Throws<ArgumentOutOfRangeException>(() => Dist.BuildIsa(dir.File("sa.bin"), 0, dir.File("isa.bin"), Opt(dir, 0, 8, 8, 8)));
        Assert.Throws<ArgumentOutOfRangeException>(() => Dist.BuildIsa(dir.File("sa.bin"), 0, dir.File("isa.bin"), Opt(dir, 4, 0, 8, 8)));
        Assert.Throws<ArgumentOutOfRangeException>(() => Dist.BuildIsa(dir.File("sa.bin"), 0, dir.File("isa.bin"), Opt(dir, 4, 8, 3, 8)));
        Assert.Throws<ArgumentOutOfRangeException>(() => Dist.BuildIsa(dir.File("sa.bin"), 0, dir.File("isa.bin"), Opt(dir, 4, 8, 8, 0)));
        Assert.Throws<NotSupportedException>(() => Dist.BuildIsa(dir.File("sa.bin"), 4_294_967_294L, dir.File("isa.bin"), Opt(dir, 4, 8, 8, 8)), "positions are 32-bit: n + 1 must fit");
        var m = new IsaMapper(0);
        Assert.Throws<ArgumentOutOfRangeException>(() => Dist.Scatter(dir.File("sa.bin"), 1, ref m, dir.File("o.bin"), new DistTarget<uint>(-1, 0, false), Opt(dir, 4, 8, 8, 8)));
        Assert.Throws<ArgumentOutOfRangeException>(() => Dist.Scatter(dir.File("sa.bin"), -1, ref m, dir.File("o.bin"), new DistTarget<uint>(1, 0, false), Opt(dir, 4, 8, 8, 8)));
    }

    [Test]
    public void DiskErrors_AreSurfaced_AndCleanedUp()
    {
        using var dir = new TempDir();
        FileBits.WriteUInts(dir.File("sa.bin"), Enumerable.Range(0, 50).Select(i => (uint)(49 - i)).ToArray());
        // the target is in a directory that does not exist: the writer fails after pass 1 has written chunk files
        var budget = new ArenaBudget(long.MaxValue);
        var counters = new IoCounters();
        Assert.Throws<DirectoryNotFoundException>(() => Dist.BuildIsa(dir.File("sa.bin"), 49, dir.File("no/such/isa.bin"), Opt(dir, 10, 16, 16, 16, budget, counters)));
        AssertClean(dir, budget, counters);
        // the scratch directory is a file
        File.WriteAllText(dir.File("blocker"), "x");
        var o = Opt(dir, 10, 16, 16, 16) with { ScratchDir = dir.File("blocker") };
        Assert.Throws<IOException>(() => Dist.BuildIsa(dir.File("sa.bin"), 49, dir.File("isa.bin"), o));
    }

    [Test]
    public void TheScratchDirectory_IsCreatedWhenMissing_AndKeepsForeignFiles()
    {
        using var dir = new TempDir();
        FileBits.WriteUInts(dir.File("sa.bin"), [2u, 0, 1]);
        string scratch = dir.Sub("deep");
        Directory.CreateDirectory(scratch);
        File.WriteAllText(Path.Combine(scratch, "keep.txt"), "mine");
        var o = Opt(dir, 1, 8, 8, 8) with { ScratchDir = scratch };
        Dist.BuildIsa(dir.File("sa.bin"), 2, dir.File("isa.bin"), o);
        Assert.That(Directory.GetFileSystemEntries(scratch).Select(Path.GetFileName), Is.EqualTo(new[] { "keep.txt" }), "only its own chunk files are removed");
        string created = dir.Sub("new/nested/scratch");
        Dist.BuildIsa(dir.File("sa.bin"), 2, dir.File("isa.bin"), Opt(dir, 1, 8, 8, 8) with { ScratchDir = created });
        Assert.That(Directory.Exists(created), Is.True);
        Assert.That(FileBits.ReadUInts(dir.File("isa.bin")), Is.EqualTo(IsaOracle([2u, 0, 1])));
    }
}
