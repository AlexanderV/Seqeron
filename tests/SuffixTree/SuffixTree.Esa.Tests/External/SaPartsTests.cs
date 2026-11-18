using System.Buffers.Binary;
using System.Runtime.InteropServices;
using SuffixTree.Esa.External;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.External;

/// <summary>The difference cover of PH1 (work package 1b.7d): the cover property of every period that can be asked for, proved by exhaustion in the test itself, and the tables the comparator reads.</summary>
[TestFixture]
public class DifferenceCoverTests
{
    /// <summary>Every residue difference (b - a) mod X, the definition of a difference cover, computed here and not by the class under test.</summary>
    private static bool Covers(int[] residues, int x)
    {
        var seen = new bool[x];
        foreach (int a in residues)
        {
            foreach (int b in residues) seen[(((b - a) % x) + x) % x] = true;
        }
        return seen.All(s => s);
    }

    [Test]
    [CancelAfter(600_000)]
    public void EveryPeriodFrom2To4096_IsACover_ByExhaustion_WithSortedDistinctResidues()
    {
        for (int x = 2; x <= 4096; x++)
        {
            DifferenceCover cover = DifferenceCover.Create(x);
            Assert.That(cover.Period, Is.EqualTo(x));
            int[] r = cover.Residues;
            Assert.That(r, Is.Ordered.Ascending.And.Unique, $"X = {x}");
            Assert.That(r[0], Is.GreaterThanOrEqualTo(0));
            Assert.That(r[^1], Is.LessThan(x));
            Assert.That(Covers(r, x), Is.True, $"X = {x}: not a difference cover");
            Assert.That(cover.Verify(), Is.True, $"X = {x}: Verify disagrees");
        }
    }

    [Test]
    public void TheSingerPeriods_ArePerfectDifferenceSets_OfQPlusOneResidues()
    {
        foreach (int q in new[] { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43 })
        {
            int x = (q * q) + q + 1;
            DifferenceCover cover = DifferenceCover.Create(x);
            Assert.That(cover.Kind, Is.EqualTo(CoverKind.Singer), $"q = {q}");
            Assert.That(cover.Residues, Has.Length.EqualTo(q + 1));
            var count = new int[x];
            foreach (int a in cover.Residues)
            {
                foreach (int b in cover.Residues)
                {
                    if (a != b) count[(((b - a) % x) + x) % x]++;
                }
            }
            for (int d = 1; d < x; d++) Assert.That(count[d], Is.EqualTo(1), $"q = {q}: difference {d} occurs {count[d]} times");
            Assert.That(cover.IsPerfect(), Is.True);
        }
    }

    [Test]
    public void AnImperfectCover_IsNotPerfect_AndAPrunedOneIsKnownAsSuch()
    {
        DifferenceCover cover = DifferenceCover.Create(9);
        Assert.That(cover.Kind, Is.EqualTo(CoverKind.Pruned));
        Assert.That(cover.IsPerfect(), Is.False);
        Assert.That(DifferenceCover.Create(1024).Kind, Is.EqualTo(CoverKind.Pruned));
    }

    [Test]
    public void TheNominalPeriodsOfTheAdr_MapToTheNearestSingerPeriod()
    {
        Assert.That(DifferenceCover.NearestSingerPeriod(512), Is.EqualTo(553));
        Assert.That(DifferenceCover.NearestSingerPeriod(1024), Is.EqualTo(993));
        Assert.That(DifferenceCover.NearestSingerPeriod(2048), Is.EqualTo(1893));
        foreach (int x in new[] { 7, 13, 31, 57, 133, 183, 993 }) Assert.That(DifferenceCover.NearestSingerPeriod(x), Is.EqualTo(x));
        Assert.That(DifferenceCover.NearestSingerPeriod(2), Is.EqualTo(7), "the smallest Singer period is the fallback below it");
        Assert.That(DifferenceCover.NearestSingerPeriod(8), Is.EqualTo(7));
        Assert.That(DifferenceCover.NearestSingerPeriod(11), Is.EqualTo(13));
    }

    [Test]
    public void TheDefaultPeriod_IsTheSingerSetOf993_With32Residues()
    {
        DifferenceCover cover = DifferenceCover.Create(DifferenceCover.NearestSingerPeriod(1024));
        Assert.That(cover.Period, Is.EqualTo(993));
        Assert.That(cover.Kind, Is.EqualTo(CoverKind.Singer));
        Assert.That(cover.Residues, Has.Length.EqualTo(32));
        Assert.That(cover.Density, Is.EqualTo(32.0 / 993).Within(1e-12));
        Assert.That(cover.PeriodBits, Is.EqualTo(10));
    }

    [Test]
    public void PeriodBits_IsLog2OfThePeriodRoundedUp()
    {
        foreach ((int x, int bits) in new[] { (2, 1), (3, 2), (4, 2), (5, 3), (7, 3), (8, 3), (9, 4), (16, 4), (17, 5), (993, 10), (1024, 10), (1025, 11) })
            Assert.That(DifferenceCover.Create(x).PeriodBits, Is.EqualTo(bits), $"X = {x}");
    }

    [Test]
    public void Witness_GivesForEveryDifferenceAnOffsetThatLandsBothSuffixesInTheSample()
    {
        foreach (int x in new[] { 2, 3, 7, 9, 13, 16, 21, 25, 31, 64, 100, 257, 993 })
        {
            DifferenceCover cover = DifferenceCover.Create(x);
            var member = new bool[x];
            foreach (int r in cover.Residues) member[r] = true;
            for (int d = 0; d < x; d++)
            {
                int a = cover.Witness[d];
                Assert.That(a, Is.InRange(0, x - 1), $"X = {x}, d = {d}: no witness");
                Assert.That(member[a] && member[(a + d) % x], Is.True, $"X = {x}, d = {d}: a = {a} and a + d are not both in D");
            }
        }
    }

    [Test]
    public void IndexOf_IsTheInverseOfResidues_AndMinusOneElsewhere()
    {
        foreach (int x in new[] { 7, 16, 31, 100 })
        {
            DifferenceCover cover = DifferenceCover.Create(x);
            Assert.That(cover.IndexOf, Has.Length.EqualTo(x));
            for (int i = 0; i < cover.Residues.Length; i++) Assert.That(cover.IndexOf[cover.Residues[i]], Is.EqualTo(i));
            Assert.That(cover.IndexOf.Count(v => v >= 0), Is.EqualTo(cover.Residues.Length));
            Assert.That(cover.IndexOf.Where(v => v < 0).Distinct(), Is.EqualTo(new[] { -1 }));
        }
    }

    [Test]
    public void SampleCount_IsTheNumberOfPositionsUpToN_WithAResidueInD()
    {
        foreach (int x in new[] { 7, 13, 16, 31, 64 })
        {
            DifferenceCover cover = DifferenceCover.Create(x);
            var member = new bool[x];
            foreach (int r in cover.Residues) member[r] = true;
            for (int n = 0; n <= 300; n++)
            {
                long brute = Enumerable.Range(0, n + 1).Count(p => member[p % x]);
                Assert.That(cover.SampleCount(n), Is.EqualTo(brute), $"X = {x}, n = {n}");
            }
        }
        Assert.That(DifferenceCover.Create(993).SampleCount(3_000_000_000), Is.EqualTo((3_000_000_000 / 993 * 32) + DifferenceCover.Create(993).Residues.Count(r => r <= 3_000_000_000 % 993)));
    }

    [Test]
    public void TheCover_IsDeterministic_AndTheArgumentsAreChecked()
    {
        foreach (int x in new[] { 9, 16, 100, 1000 })
        {
            int[] first = DifferenceCover.Create(x).Residues;
            int[] second = DifferenceCover.Create(x).Residues;
            Assert.That(first, Is.EqualTo(second));
        }
        Assert.Throws<ArgumentOutOfRangeException>(() => DifferenceCover.Create(1));
        Assert.Throws<ArgumentOutOfRangeException>(() => DifferenceCover.Create(0));
        Assert.Throws<ArgumentOutOfRangeException>(() => DifferenceCover.Create(-5));
    }
}

/// <summary>The arithmetic of the sample of one text: which positions are in it, where a position lives in the rank file and in the reduced string.</summary>
[TestFixture]
public class SampleGeometryTests
{
    [Test]
    public void SlotsAndReducedIndices_AreBijections_ForEveryNUpTo130_AndSeveralPeriods()
    {
        foreach (int x in new[] { 7, 9, 13, 16 })
        {
            DifferenceCover cover = DifferenceCover.Create(x);
            int m = cover.Residues.Length;
            for (int n = 0; n <= 130; n++)
            {
                var g = new SampleGeometry(n, cover);
                var positions = new List<long>();
                for (long p = 0; p <= n; p++)
                {
                    if (cover.IndexOf[(int)(p % x)] >= 0) positions.Add(p);
                }
                Assert.That(g.Samples, Is.EqualTo(positions.Count), $"X = {x}, n = {n}");
                Assert.That(g.Slots, Is.EqualTo(((n / x) + 1L) * m));
                Assert.That(g.Classes, Is.EqualTo(m));
                Assert.That(g.Period, Is.EqualTo(x));
                Assert.That(g.N, Is.EqualTo(n));
                Assert.That(g.Cover, Is.SameAs(cover));
                Assert.That(g.RedLength, Is.EqualTo(positions.Count + m));
                // the slot of a sample position is block-major, the real slots are 0 .. Samples - 1, and PositionOfSlot inverts SlotOf
                var slots = positions.Select(g.SlotOf).ToList();
                Assert.That(slots.Order(), Is.EqualTo(Enumerable.Range(0, positions.Count).Select(i => (long)i)), $"X = {x}, n = {n}: the real slots are 0 .. Samples - 1");
                Assert.That(slots, Is.Ordered.Ascending, "block-major slots increase with the position");
                foreach (long p in positions) Assert.That(g.PositionOfSlot(g.SlotOf(p)), Is.EqualTo(p));
                // the reduced string is class-major with one separator per class; TrySlotOfRed is false exactly on the separators
                var reducedIndex = new HashSet<long>();
                long lastOfClass = -1;
                for (int c = 0; c < m; c++)
                {
                    long count = positions.Count(p => cover.IndexOf[(int)(p % x)] == c);
                    Assert.That(g.ClassCount[c], Is.EqualTo(count));
                    Assert.That(g.ClassStart[c + 1] - g.ClassStart[c], Is.EqualTo(count));
                    Assert.That(g.RedStart[c], Is.EqualTo(g.ClassStart[c] + c));
                    Assert.That(g.RedStart[c + 1], Is.EqualTo(g.RedStart[c] + count + 1), "a class and its separator");
                    lastOfClass = g.RedStart[c] + count;
                    Assert.That(g.TrySlotOfRed(lastOfClass, out _), Is.False, $"the separator of class {c}");
                }
                Assert.That(lastOfClass + 1, Is.EqualTo(g.RedLength));
                foreach (long p in positions)
                {
                    long slot = g.SlotOf(p);
                    long red = g.RedIndexOfSlot(slot);
                    Assert.That(reducedIndex.Add(red), Is.True, "distinct reduced indices");
                    Assert.That(g.TrySlotOfRed(red, out long back), Is.True);
                    Assert.That(back, Is.EqualTo(slot), $"X = {x}, n = {n}, p = {p}");
                }
                Assert.That(reducedIndex, Has.Count.EqualTo(positions.Count));
                int separators = 0;
                for (long q = 0; q < g.RedLength; q++)
                {
                    if (!g.TrySlotOfRed(q, out long slot))
                    {
                        separators++;
                        Assert.That(slot, Is.EqualTo(-1));
                    }
                    else
                    {
                        Assert.That(reducedIndex.Contains(q), Is.True);
                    }
                }
                Assert.That(separators, Is.EqualTo(m));
            }
        }
    }

    [Test]
    public void WithinAClass_TheReducedIndexIncreasesWithThePosition_TheOrderOfTheChain()
    {
        DifferenceCover cover = DifferenceCover.Create(13);
        var g = new SampleGeometry(200, cover);
        for (int c = 0; c < cover.Residues.Length; c++)
        {
            long previous = -1;
            for (long p = cover.Residues[c]; p <= 200; p += 13)
            {
                long red = g.RedIndexOfSlot(g.SlotOf(p));
                Assert.That(red, Is.EqualTo(previous < 0 ? g.RedStart[c] : previous + 1));
                previous = red;
            }
        }
    }

    [Test]
    public void ARealSize_HasTheExpectedCounts()
    {
        DifferenceCover cover = DifferenceCover.Create(993);
        var g = new SampleGeometry(3_000_000_000, cover);
        Assert.That(g.Slots, Is.EqualTo((3_000_000_000 / 993 + 1) * 32L));
        Assert.That(g.Samples, Is.LessThanOrEqualTo(g.Slots));
        Assert.That(g.Slots - g.Samples, Is.LessThan(32));
        Assert.That(g.RedLength, Is.EqualTo(g.Samples + 32));
    }

    [Test]
    public void ANegativeLength_IsRefused()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => _ = new SampleGeometry(-1, DifferenceCover.Create(7)));
    }
}

/// <summary>The reader of the reduced strings (32 bits per symbol): the same layout and the same answers as the byte-symbol reader, on the symbols of the recursion.</summary>
[TestFixture]
public unsafe class ReducedTextReaderTests
{
    private TempDir _dir = null!;

    [SetUp]
    public void SetUp() => _dir = new TempDir();

    [TearDown]
    public void TearDown() => _dir.Dispose();

    /// <summary>The file of the reduced string: 32-bit big-endian symbols and the 64 zero bytes of padding.</summary>
    private static byte[] File32(int[] symbols)
    {
        var bytes = new byte[(symbols.Length * 4) + 64];
        for (int i = 0; i < symbols.Length; i++) BinaryPrimitives.WriteInt32BigEndian(bytes.AsSpan(i * 4), symbols[i]);
        return bytes;
    }

    private static int[] RandomSymbols(int n, int sigma, ulong seed, int repeat = 0)
    {
        var rng = new SplitMix64(seed);
        var s = new int[n];
        for (int i = 0; i < n; i++) s[i] = sigma == int.MaxValue ? (int)(rng.Next() >> 33) : rng.NextInt(sigma);
        if (repeat > 0)
        {
            for (int i = repeat; i < n; i++)
            {
                if ((i / repeat) % 3 != 0) s[i] = s[i - repeat];
            }
        }
        return s;
    }

    private MappedTextReader Open(int[] symbols)
    {
        string path = _dir.File("red-" + Guid.NewGuid().ToString("N")[..6]);
        File.WriteAllBytes(path, File32(symbols));
        return MappedTextReader.OpenReduced(path, symbols.Length);
    }

    [Test]
    public void TheShape_IsThirtyTwoBitsPerSymbol_WithTwoSymbolsInAKey()
    {
        using MappedTextReader r = Open(RandomSymbols(50, 100, 1));
        Assert.That(r.Length, Is.EqualTo(50));
        Assert.That(r.Bits, Is.EqualTo(32));
        Assert.That(r.KeySymbols, Is.EqualTo(2));
        Assert.That(r.KeyShift, Is.Zero);
        Assert.That(r.VectorSymbols, Is.EqualTo(8));
    }

    [Test]
    public void SymbolsAndKeys_AreTheOracle_AtEveryPosition_IncludingTheLargestValuesAndThePadding()
    {
        foreach (int sigma in new[] { 3, 1000, int.MaxValue })
        {
            int[] s = RandomSymbols(300, sigma, 7);
            s[0] = int.MaxValue;
            s[^1] = int.MaxValue;
            using MappedTextReader r = Open(s);
            for (int i = -0; i < s.Length + 3; i++)
            {
                uint sym = i < s.Length ? (uint)s[i] : 0;
                uint next = i + 1 < s.Length ? (uint)s[i + 1] : 0;
                Assert.That(r.Symbol(i), Is.EqualTo(sym), $"sigma {sigma}, symbol {i}");
                Assert.That(r.Key(i), Is.EqualTo(((ulong)sym << 32) | next), $"sigma {sigma}, key {i}");
            }
        }
    }

    [Test]
    public void TheCommonPrefix_AndTheSuffixComparison_AreTheOracle_OnTheVectorAndTheScalarPath()
    {
        int[] s = RandomSymbols(700, 3, 11, repeat: 90);
        s[300] = 77;
        using MappedTextReader r = Open(s);
        var rng = new SplitMix64(5);
        for (int round = 0; round < 3000; round++)
        {
            int i = rng.NextInt(s.Length + 1), j = rng.NextInt(s.Length + 1);
            long max = rng.NextInt(s.Length + 40);
            int expected = 0;
            int cap = (int)Math.Min(max, s.Length - Math.Max(i, j));
            while (expected < cap && s[i + expected] == s[j + expected]) expected++;
            long steps = 0;
            foreach (bool vectors in new[] { true, false })
            {
                r.UseVectors = vectors;
                Assert.That(r.Lce(i, j, max, ref steps), Is.EqualTo(expected), $"lce({i}, {j}, {max}) vectors {vectors}");
                int cmp = r.CompareSuffixes(i, j, 0, ref steps, out long lce);
                int want = Math.Sign(((ReadOnlySpan<int>)s.AsSpan(Math.Min(i, s.Length))).SequenceCompareTo(s.AsSpan(Math.Min(j, s.Length))));
                Assert.That(cmp, Is.EqualTo(want), $"compare({i}, {j}) vectors {vectors}");
                if (i != j) Assert.That(lce, Is.EqualTo(Math.Min(TextLce(s, i, j), s.Length - Math.Max(i, j))));
            }
        }
    }

    private static int TextLce(int[] s, int i, int j)
    {
        int k = 0;
        while (i + k < s.Length && j + k < s.Length && s[i + k] == s[j + k]) k++;
        return k;
    }

    [Test]
    public void FromMemory_GivesTheSameReader()
    {
        int[] s = RandomSymbols(120, 5, 3);
        byte[] bytes = File32(s);
        fixed (byte* p = bytes)
        {
            using MappedTextReader r = MappedTextReader.FromMemoryReduced(p, bytes.Length, s.Length);
            for (int i = 0; i < s.Length; i++) Assert.That(r.Symbol(i), Is.EqualTo((uint)s[i]));
            Assert.That(r.Bits, Is.EqualTo(32));
        }
    }

    [Test]
    public void ADamagedFile_IsInvalidData_NotAnArgumentError()
    {
        int[] s = RandomSymbols(40, 9, 2);
        byte[] good = File32(s);
        string path = _dir.File("bad");
        foreach ((string what, byte[] bytes) in new (string, byte[])[]
        {
            ("short", good[..^1]),
            ("long", [.. good, 0]),
            ("padding", [.. good[..^1], 1]),
            ("empty", []),
        })
        {
            File.WriteAllBytes(path, bytes);
            Assert.Throws<InvalidDataException>(() => MappedTextReader.OpenReduced(path, s.Length), what);
        }
        File.WriteAllBytes(path, good);
        Assert.Throws<InvalidDataException>(() => MappedTextReader.OpenReduced(path, s.Length + 1), "a longer text than the file holds");
        Assert.Throws<InvalidDataException>(() => MappedTextReader.OpenReduced(path, -1));
    }

    [Test]
    public void TheByteReader_StillRefusesAWideSymbol_AndAWideReaderCannotUnpackBytes()
    {
        int[] s = RandomSymbols(40, 9, 2);
        string path = _dir.File("w");
        File.WriteAllBytes(path, File32(s));
        Assert.Throws<InvalidDataException>(() => MappedTextReader.Open(path, s.Length, 32), "the file codec has 1 to 8 bits");
        using MappedTextReader r = MappedTextReader.OpenReduced(path, s.Length);
        Assert.Throws<NotSupportedException>(() => r.ReadSymbols(0, new byte[4]));
    }

    [Test]
    public void ADisposedWideReader_Throws_AndReleasesItsFile()
    {
        int[] s = RandomSymbols(40, 9, 2);
        string path = _dir.File("d");
        File.WriteAllBytes(path, File32(s));
        MappedTextReader r = MappedTextReader.OpenReduced(path, s.Length);
        Assert.That(ProcMaps.Contains(path) ?? true, Is.True);
        r.Dispose();
        Assert.That(ProcMaps.Contains(path) ?? false, Is.False);
        Assert.Throws<ObjectDisposedException>(() => r.Symbol(3 + 100));
    }
}

/// <summary>The remainder by the period without a division, and the read-only rank file that PH3 looks the sample ranks up in.</summary>
[TestFixture]
public class RankTableTests
{
    private TempDir _dir = null!;

    [SetUp]
    public void SetUp() => _dir = new TempDir();

    [TearDown]
    public void TearDown() => _dir.Dispose();

    [Test]
    public void FastMod32_IsDivisionAndRemainder_ForEveryDivisorUpTo2100_AndTheExtremesOfTheDomain()
    {
        var rng = new SplitMix64(3);
        var probes = new List<uint> { 0, 1, 2, 3, 4, uint.MaxValue, uint.MaxValue - 1, 1u << 31, (1u << 31) - 1, 3_000_000_000u, 4_294_967_291u };
        for (int i = 0; i < 400; i++) probes.Add((uint)rng.Next());
        for (uint d = 2; d <= 2100; d++)
        {
            var f = new FastMod32(d);
            Assert.That(f.Divisor, Is.EqualTo(d));
            foreach (uint a in probes)
            {
                Assert.That(f.Div(a), Is.EqualTo(a / d), $"{a} / {d}");
                Assert.That(f.Mod(a), Is.EqualTo(a % d), $"{a} % {d}");
            }
            for (uint a = 0; a < 3 * d; a++)
            {
                Assert.That(f.Div(a), Is.EqualTo(a / d));
                Assert.That(f.Mod(a), Is.EqualTo(a % d));
            }
        }
        foreach (uint d in new[] { 65_535u, 65_536u, 1_000_003u, (1u << 31) - 1, 1u << 31, uint.MaxValue })
        {
            var f = new FastMod32(d);
            foreach (uint a in probes)
            {
                Assert.That(f.Div(a), Is.EqualTo(a / d));
                Assert.That(f.Mod(a), Is.EqualTo(a % d));
            }
        }
    }

    [Test]
    public void FastMod32_RefusesADivisorBelowTwo()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => _ = new FastMod32(0));
        Assert.Throws<ArgumentOutOfRangeException>(() => _ = new FastMod32(1));
    }

    private string WriteRanks(SaCase c, DifferenceCover cover)
    {
        uint[] ranks = DcsOracle.ExpectedRanks(c.Symbols, cover, out _);
        string path = _dir.File("ranks-" + Guid.NewGuid().ToString("N")[..6]);
        File.WriteAllBytes(path, MemoryMarshal.AsBytes(ranks.AsSpan()).ToArray());
        return path;
    }

    [Test]
    public void Rank_IsTheValueOfTheSlotOfThePosition_AtEverySamplePosition()
    {
        foreach (int x in new[] { 7, 16, 31 })
        {
            DifferenceCover cover = DifferenceCover.Create(x);
            foreach (string text in new[] { StringFamilies.RandomDna(300, 4), StringFamilies.Fibonacci(200), "a", string.Empty, "mississippi" })
            {
                SaCase c = SaCase.Of(text);
                uint[] expected = DcsOracle.ExpectedRanks(c.Symbols, cover, out _);
                using RankTable t = RankTable.Open(WriteRanks(c, cover), cover, c.N);
                Assert.That(t.N, Is.EqualTo(c.N));
                Assert.That(t.Slots, Is.EqualTo(expected.Length));
                Assert.That(t.Cover, Is.SameAs(cover));
                Assert.That(t.Mod.Divisor, Is.EqualTo((uint)x));
                for (int p = 0; p <= c.N; p++)
                {
                    int idx = cover.IndexOf[p % x];
                    if (idx < 0) continue;
                    Assert.That(t.Rank(p), Is.EqualTo(expected[((long)(p / x) * cover.Residues.Length) + idx]), $"X = {x}, p = {p}");
                }
                t.Verify();
            }
        }
    }

    [Test]
    public void ARankFileOfAnotherLength_IsInvalidData()
    {
        DifferenceCover cover = DifferenceCover.Create(7);
        SaCase c = SaCase.Of(StringFamilies.RandomDna(100, 4));
        string good = WriteRanks(c, cover);
        byte[] bytes = File.ReadAllBytes(good);
        string path = _dir.File("odd");
        foreach (byte[] wrong in new[] { bytes[..^4], [.. bytes, 0, 0, 0, 0], bytes[..^1], [] })
        {
            File.WriteAllBytes(path, wrong);
            Assert.Throws<InvalidDataException>(() => RankTable.Open(path, cover, c.N));
        }
        Assert.Throws<InvalidDataException>(() => RankTable.Open(good, cover, c.N + 7), "the file belongs to a shorter text");
    }

    [Test]
    public void Verify_FindsARepeatedRank_AZeroRank_AndAnUsedUnusedSlot()
    {
        DifferenceCover cover = DifferenceCover.Create(7);
        SaCase c = SaCase.Of(StringFamilies.RandomDna(100, 4));
        uint[] ranks = DcsOracle.ExpectedRanks(c.Symbols, cover, out _);
        var g = new SampleGeometry(c.N, cover);
        string path = _dir.File("v");
        foreach ((string what, Action<uint[]> damage) in new (string, Action<uint[]>)[]
        {
            ("repeated", r => r[3] = r[4]),
            ("zero", r => r[5] = 0),
            ("too large", r => r[5] = (uint)(g.Samples + cover.Residues.Length + 1)),
            ("unused slot", r => r[^1] = 9),
        })
        {
            uint[] copy = [.. ranks];
            if (what == "unused slot" && g.Slots == g.Samples) continue;
            damage(copy);
            File.WriteAllBytes(path, MemoryMarshal.AsBytes(copy.AsSpan()).ToArray());
            using RankTable t = RankTable.Open(path, cover, c.N);
            Assert.Throws<InvalidDataException>(() => t.Verify(), what);
        }
    }

    [Test]
    public void ADisposedTable_ReleasesItsMapping()
    {
        DifferenceCover cover = DifferenceCover.Create(7);
        SaCase c = SaCase.Of("banana");
        string path = WriteRanks(c, cover);
        RankTable t = RankTable.Open(path, cover, c.N);
        Assert.That(ProcMaps.Contains(path) ?? true, Is.True);
        t.Dispose();
        t.Dispose();
        Assert.That(ProcMaps.Contains(path) ?? false, Is.False);
        Assert.Throws<ObjectDisposedException>(t.Verify);
    }
}

/// <summary>The generator of the splitters and the work list of the bucket tasks.</summary>
[TestFixture]
public class SaPlumbingTests
{
    [Test]
    public void SeededRandom_IsTheSplitMix64Sequence_AndNextIntIsTheModuloRule()
    {
        foreach (ulong seed in new ulong[] { 0, 1, 0x9E3779B97F4A7C15UL, ulong.MaxValue })
        {
            var a = new SeededRandom(seed);
            var b = new SplitMix64(seed);
            for (int i = 0; i < 1000; i++) Assert.That(a.Next(), Is.EqualTo(b.Next()), $"seed {seed}, value {i}");
            var c = new SeededRandom(seed);
            var d = new SplitMix64(seed);
            foreach (int bound in new[] { 1, 2, 3, 10, 1000, int.MaxValue })
                for (int i = 0; i < 50; i++) Assert.That(c.NextInt(bound), Is.EqualTo(d.NextInt(bound)));
        }
        Assert.Throws<ArgumentOutOfRangeException>(() => new SeededRandom(1).NextInt(0));
        Assert.Throws<ArgumentOutOfRangeException>(() => new SeededRandom(1).NextInt(-3));
    }

    [Test]
    public void WorkPool_IsLastInFirstOut_AndEndsWhenTheListIsEmptyAndNoTaskIsRunning()
    {
        var pool = new WorkPool<string>();
        pool.Push("a");
        pool.Push("b");
        pool.Push("c");
        Assert.That(pool.Take(out string t1), Is.True);
        Assert.That(t1, Is.EqualTo("c"));
        pool.Done();
        Assert.That(pool.Take(out string t2), Is.True);
        Assert.That(t2, Is.EqualTo("b"));
        pool.Done();
        Assert.That(pool.Take(out string t3), Is.True);
        Assert.That(t3, Is.EqualTo("a"));
        pool.Done();
        Assert.That(pool.Take(out _), Is.False);
        Assert.That(pool.Error, Is.Null);
    }

    [Test]
    public void WorkPool_ATaskThatPushesChildren_KeepsTheOtherThreadsWaiting_UntilTheLastTaskIsDone()
    {
        foreach (int threads in new[] { 1, 2, 4, 8 })
        {
            var pool = new WorkPool<int[]>();
            pool.Push([0, 12]);
            long sum = 0;
            long count = 0;
            var workers = Enumerable.Range(0, threads).Select(_ => new Thread(() =>
            {
                while (pool.Take(out int[] task))
                {
                    try
                    {
                        Interlocked.Increment(ref count);
                        Interlocked.Add(ref sum, task[0]);
                        Thread.SpinWait(2000);
                        if (task[1] > 0)
                        {
                            pool.Push([task[0] + 1, task[1] - 1]);
                            pool.Push([task[0] + 2, task[1] - 1]);
                        }
                    }
                    finally
                    {
                        pool.Done();
                    }
                }
            })).ToList();
            workers.ForEach(w => w.Start());
            workers.ForEach(w => Assert.That(w.Join(60_000), Is.True, "a worker did not end"));
            Assert.That(count, Is.EqualTo((1L << 13) - 1), $"{threads} threads: every task of the binary tree ran exactly once");
            Assert.That(sum, Is.GreaterThan(0));
        }
    }

    [Test]
    public void WorkPool_AFailureEndsEveryTake_KeepsTheFirstError_AndDropsThePendingTasks()
    {
        var pool = new WorkPool<string>();
        pool.Push("a");
        pool.Push("b");
        Assert.That(pool.Take(out _), Is.True);
        var first = new IOException("first");
        pool.Fail(first);
        pool.Fail(new InvalidOperationException("second"));
        Assert.That(pool.Error, Is.SameAs(first));
        pool.Done();
        Assert.That(pool.Take(out _), Is.False, "the failed pool hands out nothing more");
        Assert.That(pool.Take(out _), Is.False);
    }

    [Test]
    public void WorkPool_ABlockedTake_IsWokenByAFailure()
    {
        var pool = new WorkPool<string>();
        pool.Push("a");
        Assert.That(pool.Take(out _), Is.True); // one task is running: a second taker has to wait
        bool? got = null;
        var waiter = new Thread(() => got = pool.Take(out _));
        waiter.Start();
        using (var pause = new ManualResetEventSlim()) pause.Wait(100);
        Assert.That(got, Is.Null, "the second taker waits");
        pool.Fail(new IOException("x"));
        Assert.That(waiter.Join(10_000), Is.True);
        Assert.That(got, Is.False);
    }
}

/// <summary>The runs of the sample buckets (memory or scratch file), the ledger of the scratch bytes, and the reader that scatters a run into the name window.</summary>
[TestFixture]
public class BucketRunTests
{
    private TempDir _dir = null!;

    [SetUp]
    public void SetUp() => _dir = new TempDir();

    [TearDown]
    public void TearDown() => _dir.Dispose();

    private static ulong Rec(uint slot, uint name) => ((ulong)slot << 32) | name;

    [Test]
    public void ARunOfAtMostTheSpillSize_StaysInMemory_AndALargerOneGoesToAFile()
    {
        var counters = new IoCounters();
        var ledger = new ScratchLedger(counters);
        var store = new RunStore(_dir.Path, "t", 16, ledger, counters);
        using BucketRun small = store.Create([Rec(1, 1), Rec(2, 2)]);
        Assert.That(small.InMemory, Is.True);
        Assert.That(small.Count, Is.EqualTo(2));
        Assert.That(small.Memory.ToArray(), Is.EqualTo(new[] { Rec(1, 1), Rec(2, 2) }));
        Assert.That(ledger.Live, Is.Zero);
        Assert.That(Directory.GetFiles(_dir.Path), Is.Empty);
        BucketRun big = store.Create([Rec(1, 1), Rec(2, 2), Rec(3, 3)]);
        Assert.That(big.InMemory, Is.False);
        Assert.That(ledger.Live, Is.EqualTo(24));
        Assert.That(counters.ScratchLiveBytes, Is.EqualTo(24));
        Assert.That(counters.WriteBytes, Is.EqualTo(24));
        string file = Directory.GetFiles(_dir.Path).Single();
        Assert.That(new FileInfo(file).Length, Is.EqualTo(24));
        Assert.That(Path.GetFileName(file), Does.StartWith("t-run-").And.EndWith(".bin"));
        big.Dispose();
        big.Dispose();
        Assert.That(File.Exists(file), Is.False);
        Assert.That(ledger.Live, Is.Zero);
        Assert.That(counters.ScratchLiveBytes, Is.Zero);
        Assert.That(counters.ScratchPeakBytes, Is.EqualTo(24));
    }

    [Test]
    public void TheFilesOfOneStore_AreNumberedFromOne_AndTwoStoresNeverShareAName()
    {
        var ledger = new ScratchLedger(null);
        var a = new RunStore(_dir.Path, "a", 0, ledger, null);
        var b = new RunStore(_dir.Path, "b", 0, ledger, null);
        using BucketRun r1 = a.Create([Rec(1, 1)]);
        using BucketRun r2 = a.Create([Rec(2, 1)]);
        using BucketRun r3 = b.Create([Rec(3, 1)]);
        Assert.That(Directory.GetFiles(_dir.Path).Select(Path.GetFileName).Order(), Is.EqualTo(new[] { "a-run-00000001.bin", "a-run-00000002.bin", "b-run-00000001.bin" }));
    }

    [Test]
    public void ReadAt_ReadsTheRecordsFromTheIndex_NeverBeyondTheEnd_AndCountsTheBytes()
    {
        var counters = new IoCounters();
        var store = new RunStore(_dir.Path, "t", 0, new ScratchLedger(counters), counters);
        using BucketRun run = store.Create([Rec(1, 11), Rec(2, 12), Rec(3, 13), Rec(4, 14), Rec(5, 15)]);
        var buffer = new ulong[8];
        Assert.That(run.ReadAt(0, buffer), Is.EqualTo(5), "fewer than the buffer holds");
        Assert.That(buffer[..5], Is.EqualTo(new[] { Rec(1, 11), Rec(2, 12), Rec(3, 13), Rec(4, 14), Rec(5, 15) }));
        Assert.That(run.ReadAt(3, buffer.AsSpan(0, 8)), Is.EqualTo(2));
        Assert.That(buffer[..2], Is.EqualTo(new[] { Rec(4, 14), Rec(5, 15) }));
        Assert.That(run.ReadAt(1, buffer.AsSpan(0, 2)), Is.EqualTo(2), "bounded by the buffer");
        Assert.That(buffer[..2], Is.EqualTo(new[] { Rec(2, 12), Rec(3, 13) }));
        Assert.That(run.ReadAt(5, buffer), Is.Zero);
        Assert.That(counters.ReadBytes, Is.EqualTo(8 * (5 + 2 + 2)));
    }

    [Test]
    public void AFailedWrite_GivesTheBytesBack_AndLeavesNoFile()
    {
        var counters = new IoCounters();
        var ledger = new ScratchLedger(counters);
        var store = new RunStore(_dir.Sub("missing"), "t", 0, ledger, counters);
        Assert.Throws(Is.InstanceOf<IOException>(), () => store.Create([Rec(1, 1), Rec(2, 2)]));
        Assert.That(ledger.Live, Is.Zero);
        Assert.That(counters.ScratchLiveBytes, Is.Zero);
        Assert.That(counters.ScratchPeakBytes, Is.EqualTo(16), "the bytes were taken with the reservation");
    }

    [Test]
    public void TheLedger_KnowsWhatItTook_AndGivesAllOfItBack()
    {
        var counters = new IoCounters();
        counters.AddScratch(1000); // another stage of the same counters
        var ledger = new ScratchLedger(counters);
        ledger.Add(40);
        ledger.Add(2);
        ledger.Remove(10);
        Assert.That(ledger.Live, Is.EqualTo(32));
        Assert.That(counters.ScratchLiveBytes, Is.EqualTo(1032));
        ledger.ReleaseAll();
        Assert.That(ledger.Live, Is.Zero);
        Assert.That(counters.ScratchLiveBytes, Is.EqualTo(1000), "what another stage holds is not touched");
        Assert.That(counters.ScratchPeakBytes, Is.EqualTo(1042));
        Assert.DoesNotThrow(() => new ScratchLedger(null).Add(5));
    }

    private static void Drain(BucketRun run, bool inMemory)
    {
        using var arena = new BudgetedArena(new ArenaBudget(long.MaxValue), ArenaOptions.Guarded(0xCD));
        int records = inMemory ? 1 : 2;
        Buf<ulong> buf = inMemory ? default : arena.Allocate<ulong>(records, zero: false, "test");
        var reader = new RunReader(run, buf, records);
        var window = new uint[4];
        // slots 3, 5, 6, 9 with the names 1, 2, 3, 1; the window covers the slots 3..6, the next one 7..10
        Assert.That(reader.Drain(7, 10, window, 3), Is.EqualTo(3));
        Assert.That(window[0], Is.EqualTo(11u));
        Assert.That(window[2], Is.EqualTo(12u));
        Assert.That(window[3], Is.EqualTo(13u));
        Assert.That(reader.Scattered, Is.EqualTo(3));
        var next = new uint[4];
        Assert.That(reader.Drain(11, 10, next, 7), Is.EqualTo(1));
        Assert.That(next[2], Is.EqualTo(11u));
        Assert.That(reader.Scattered, Is.EqualTo(4));
        Assert.That(reader.Drain(100, 10, next, 7), Is.Zero, "the run is used up");
    }

    [Test]
    public void ARunReader_ScattersTheRecordsBelowTheLimit_WithTheOffset_FromMemoryAndFromAFile()
    {
        var store = new RunStore(_dir.Path, "t", 0, new ScratchLedger(null), null);
        var records = new[] { Rec(3, 1), Rec(5, 2), Rec(6, 3), Rec(9, 1) };
        using BucketRun file = store.Create(records);
        Assert.That(file.InMemory, Is.False);
        Drain(file, inMemory: false);
        var memoryStore = new RunStore(_dir.Path, "m", 1 << 20, new ScratchLedger(null), null);
        using BucketRun memory = memoryStore.Create(records);
        Assert.That(memory.InMemory, Is.True);
        Drain(memory, inMemory: true);
    }
}

/// <summary>The inverse permutation of the suffix array of a reduced string, written straight into the rank file by DIST.</summary>
[TestFixture]
public class RankInverterTests
{
    private TempDir _dir = null!;

    [SetUp]
    public void SetUp() => _dir = new TempDir();

    [TearDown]
    public void TearDown() => _dir.Dispose();

    /// <summary>A reduced string of the geometry: the names (1..sigma) at the sample slots, a 0 after every class, and its suffix array with the empty suffix first.</summary>
    private static (int[] Red, uint[] Sa) Reduced(SampleGeometry g, int sigma, ulong seed)
    {
        var rng = new SplitMix64(seed);
        var red = new int[g.RedLength];
        for (long q = 0; q < red.Length; q++) red[q] = g.TrySlotOfRed(q, out _) ? 1 + rng.NextInt(sigma) : 0;
        int[] order = Enumerable.Range(0, red.Length).ToArray();
        Array.Sort(order, (a, b) => Math.Sign(red.AsSpan(a).SequenceCompareTo(red.AsSpan(b))));
        var sa = new uint[red.Length + 1];
        sa[0] = (uint)red.Length;
        for (int i = 0; i < order.Length; i++) sa[i + 1] = (uint)order[i];
        return (red, sa);
    }

    private (string Sa, string Rank, DistOptions Options) Files(uint[] sa, int chunk = 7)
    {
        string path = _dir.File("sa.bin");
        FileBits.WriteUInts(path, sa);
        return (path, _dir.File("rank.bin"), new DistOptions { ScratchDir = _dir.Sub("s"), ChunkEntries = chunk, IoBufferBytes = 24, ReadBytes = 16, WindowBytes = 32 });
    }

    [Test]
    public void TheRankOfASamplePosition_IsItsIndexInTheSuffixArray_AndTheUnusedSlotsStayZero()
    {
        foreach (int x in new[] { 7, 13 })
        {
            DifferenceCover cover = DifferenceCover.Create(x);
            foreach (long n in new long[] { 1, 20, 61, 100 })
            {
                var g = new SampleGeometry(n, cover);
                (int[] red, uint[] sa) = Reduced(g, 3, (ulong)n);
                (string saPath, string rankPath, DistOptions o) = Files(sa);
                RankInverter.Invert(saPath, g, rankPath, o);
                uint[] got = FileBits.ReadUInts(rankPath);
                Assert.That(got, Has.Length.EqualTo(g.Slots));
                var rankOf = new uint[red.Length];
                for (int i = 1; i < sa.Length; i++) rankOf[sa[i]] = (uint)i;
                for (long slot = 0; slot < g.Slots; slot++)
                {
                    uint expected = slot < g.Samples ? rankOf[g.RedIndexOfSlot(slot)] : 0;
                    Assert.That(got[slot], Is.EqualTo(expected), $"X = {x}, n = {n}, slot {slot}");
                }
                File.Delete(rankPath);
                if (Directory.Exists(_dir.Sub("s"))) Directory.Delete(_dir.Sub("s"), true);
            }
        }
    }

    [Test]
    public void ASuffixArrayThatIsWrong_IsInvalidData_WhateverTheWay()
    {
        DifferenceCover cover = DifferenceCover.Create(7);
        var g = new SampleGeometry(40, cover);
        (_, uint[] good) = Reduced(g, 3, 5);
        int sampleEntry = Array.FindIndex(good, 1, q => g.TrySlotOfRed(q, out _));
        int separatorEntry = Array.FindIndex(good, 1, q => !g.TrySlotOfRed(q, out _));
        var cases = new List<(string What, uint[] Sa)>
        {
            ("an entry beyond the reduced string", [.. good.Take(3), (uint)g.RedLength + 4, .. good.Skip(4)]),
            ("short", good[..^1]),
            ("long", [.. good, 0u]),
            ("a sample position that is missing", [.. good.Select((v, i) => i == sampleEntry ? good[separatorEntry] : v)]),
        };
        foreach ((string what, uint[] sa) in cases)
        {
            (string saPath, string rankPath, DistOptions o) = Files(sa);
            Assert.Throws<InvalidDataException>(() => RankInverter.Invert(saPath, g, rankPath, o), what);
            if (Directory.Exists(_dir.Sub("s"))) Assert.That(Directory.GetFileSystemEntries(_dir.Sub("s")), Is.Empty, $"{what}: no chunk file is left");
        }
    }

    [Test]
    public void TheMapper_SkipsTheEmptySuffixAndTheSeparators_AndCountsWhatItEmits()
    {
        DifferenceCover cover = DifferenceCover.Create(7);
        var g = new SampleGeometry(30, cover);
        var mapper = new SampleRankMapper(g);
        Assert.That(mapper.TryMap(0, (uint)g.RedLength, out _, out _), Is.False, "the empty suffix");
        long separator = g.RedStart[1] - 1;
        Assert.That(mapper.TryMap(5, (uint)separator, out _, out _), Is.False, "a separator");
        Assert.That(mapper.Emitted, Is.Zero);
        long q = g.RedIndexOfSlot(4);
        Assert.That(mapper.TryMap(9, (uint)q, out uint key, out uint value), Is.True);
        Assert.That((key, value), Is.EqualTo((4u, 9u)), "the slot is the key and the index in the suffix array is the rank");
        Assert.That(mapper.Emitted, Is.EqualTo(1));
        Assert.Throws<InvalidDataException>(() => mapper.TryMap(3, (uint)g.RedLength, out _, out _));
    }
}
