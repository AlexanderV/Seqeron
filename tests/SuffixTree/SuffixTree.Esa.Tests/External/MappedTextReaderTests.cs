using System.Runtime.Intrinsics.X86;
using SuffixTree.Esa.External;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.External;

/// <summary>
/// The reader of the packed text (TEXT codec 4) against the in-RAM oracle of <see cref="TextOracle"/> (work package 1b.7b, gate G1): every position and every key of small and exhaustive texts for b = 1..8, the alignment
/// cases, past-the-end reads, bulk reads, the common prefix with the vector and the scalar path, positions and offsets above 2^31 and 2^32 on sparse files, damaged input.
/// </summary>
[TestFixture]
public class MappedTextReaderTests
{
    private static readonly int[] FirstClass = [2, 3, 4, 8];
    private static readonly int[] AllBits = [1, 2, 3, 4, 5, 6, 7, 8];

    private static void CheckEverything(PinnedText t, byte[] s, int bits, string label)
    {
        MappedTextReader r = t.Reader;
        int n = s.Length;
        Assert.That(r.Length, Is.EqualTo(n), label);
        Assert.That(r.Bits, Is.EqualTo(bits), label);
        Assert.That(r.KeySymbols, Is.EqualTo(64 / bits), label);
        Assert.That(r.KeyShift, Is.EqualTo(64 - (64 / bits * bits)), label);
        Assert.That(r.VectorSymbols, Is.EqualTo(256 / bits), label);
        for (long i = 0; i <= n + 2; i++)
        {
            uint want = i < n ? s[i] : 0u;
            Assert.That(r.Symbol(i), Is.EqualTo(want), $"{label}: Symbol({i})");
            Assert.That(r.Window(i), Is.EqualTo(TextOracle.Window(s, i, bits)), $"{label}: Window({i})");
            Assert.That(r.Key(i), Is.EqualTo(TextOracle.Key(s, i, bits)), $"{label}: Key({i})");
        }
    }

    // ---- the contract on small texts

    [Test]
    public void EveryPosition_Symbol_Window_Key_EqualTheOracle_ForEveryWidth()
    {
        foreach (int bits in AllBits)
        {
            for (int n = 0; n <= 90; n++)
            {
                byte[] s = TextOracle.Random(n, 1 << bits, (ulong)((bits * 1000) + n));
                using var t = new PinnedText(s, bits);
                CheckEverything(t, s, bits, $"b={bits} n={n}");
            }
        }
    }

    [Test]
    public void ExhaustiveSmallTexts_EveryPosition_EveryKey_EveryCommonPrefix_EveryComparison()
    {
        (int Bits, int Sigma, int MaxLength)[] cases = [(1, 2, 11), (2, 4, 5), (3, 8, 3), (4, 16, 3), (8, 256, 2), (3, 5, 4)];
        foreach ((int bits, int sigma, int maxLength) in cases)
        {
            foreach (byte[] s in TextOracle.AllTexts(sigma, maxLength))
            {
                using var t = new PinnedText(s, bits);
                CheckEverything(t, s, bits, $"b={bits} text=[{string.Join(',', s)}]");
                foreach (bool vectors in new[] { false, true })
                {
                    t.Reader.UseVectors = vectors;
                    for (long i = 0; i <= s.Length; i++)
                    {
                        for (long j = 0; j <= s.Length; j++)
                        {
                            long max = Math.Max(0, s.Length - Math.Max(i, j));
                            long steps = 0;
                            long want = TextOracle.Lce(s, i, j, max);
                            Assert.That(t.Reader.Lce(i, j, max, ref steps), Is.EqualTo(want), $"b={bits} [{string.Join(',', s)}] Lce({i},{j}) vectors={vectors}");
                            if (i != j) Assert.That(Math.Sign(t.Reader.CompareSuffixes(i, j, 0, ref steps, out long lce)), Is.EqualTo(TextOracle.CompareSuffixes(s, i, j)), $"Compare({i},{j})");
                        }
                    }
                }
            }
        }
    }

    [Test]
    public void KeyOrder_IsTheLexicographicOrderOfTheFirstSymbols_AndMonotoneInTheSuffixOrder()
    {
        foreach (int bits in FirstClass)
        {
            byte[] s = TextOracle.Random(300, 1 << bits, 77 + (ulong)bits, repeatEvery: 40);
            using var t = new PinnedText(s, bits);
            int ks = 64 / bits;
            for (int i = 0; i < s.Length; i += 3)
            {
                for (int j = 0; j < s.Length; j += 5)
                {
                    int cmp = 0;
                    for (int k = 0; k < ks && cmp == 0; k++) cmp = (i + k < s.Length ? s[i + k] : 0).CompareTo(j + k < s.Length ? s[j + k] : 0);
                    Assert.That(t.Reader.Key(i).CompareTo(t.Reader.Key(j)), Is.EqualTo(cmp), $"b={bits} i={i} j={j}");
                    if (TextOracle.CompareSuffixes(s, i, j) < 0) Assert.That(t.Reader.Key(i), Is.LessThanOrEqualTo(t.Reader.Key(j)), "suffix(i) < suffix(j) implies Key(i) <= Key(j)");
                }
            }
        }
    }

    // ---- alignment: every bit offset, word and page boundaries, the last symbols, past the end

    [Test]
    public void EveryBitOffset_EveryWidth_AroundWordAndPageBoundaries_AndAtTheEnd()
    {
        foreach (int bits in AllBits)
        {
            int n = 70_000;
            byte[] s = TextOracle.Random(n, 1 << bits, 4242 + (ulong)bits);
            using var t = new PinnedText(s, bits);
            var centres = new List<long> { 0, 4096L * 8 / bits, 8192L * 8 / bits, 65536L * 8 / bits, (1L << 16) };
            for (int c = 0; c < centres.Count; c++)
            {
                for (long i = centres[c] - 80; i <= centres[c] + 80; i++)
                {
                    if (i < 0 || i >= n) continue;
                    Assert.That(t.Reader.Window(i), Is.EqualTo(TextOracle.Window(s, i, bits)), $"b={bits} Window({i})");
                    Assert.That(t.Reader.Key(i), Is.EqualTo(TextOracle.Key(s, i, bits)), $"b={bits} Key({i})");
                    Assert.That(t.Reader.Symbol(i), Is.EqualTo((uint)s[i]), $"b={bits} Symbol({i})");
                }
            }
            for (long i = n - 150; i <= n + 200; i++)
            {
                Assert.That(t.Reader.Key(i), Is.EqualTo(TextOracle.Key(s, i, bits)), $"b={bits} tail Key({i})");
                Assert.That(t.Reader.Window(i), Is.EqualTo(TextOracle.Window(s, i, bits)), $"b={bits} tail Window({i})");
            }
        }
    }

    [Test]
    public void PastTheEnd_AnyPosition_ReadsTheZeroPadding_AndNothingOutOfTheMappingIsRead()
    {
        foreach (int bits in new[] { 2, 3, 8 })
        {
            byte[] s = TextOracle.Random(100, 1 << bits, 6);
            using var t = new PinnedText(s, bits);
            foreach (long i in new[] { 100L, 101, 164, 165, 1_000, int.MaxValue, 1L << 32, 1L << 40, long.MaxValue / 8, long.MaxValue })
            {
                Assert.That(t.Reader.Symbol(i), Is.Zero, $"b={bits} Symbol({i})");
                Assert.That(t.Reader.Key(i), Is.Zero, $"b={bits} Key({i})");
                Assert.That(t.Reader.Window(i), Is.Zero, $"b={bits} Window({i})");
            }
        }
    }

    [Test]
    public void NegativePositions_AreRefused()
    {
        using var t = new PinnedText(TextOracle.Random(100, 4, 6), 2);
        foreach (long i in new[] { -1L, -2, long.MinValue, int.MinValue })
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => t.Reader.Symbol(i), $"Symbol({i})");
            Assert.Throws<ArgumentOutOfRangeException>(() => t.Reader.Key(i), $"Key({i})");
            Assert.Throws<ArgumentOutOfRangeException>(() => t.Reader.Window(i), $"Window({i})");
        }
        long steps = 0;
        Assert.Throws<ArgumentOutOfRangeException>(() => t.Reader.Lce(-1, 0, 5, ref steps));
        Assert.Throws<ArgumentOutOfRangeException>(() => t.Reader.Lce(0, -1, 5, ref steps));
        Assert.Throws<ArgumentOutOfRangeException>(() => t.Reader.Lce(0, 1, -1, ref steps));
        Assert.Throws<ArgumentOutOfRangeException>(() => t.Reader.CompareSuffixes(-1, 0, 0, ref steps, out _));
        Assert.Throws<ArgumentOutOfRangeException>(() => t.Reader.CompareSuffixes(0, -1, 0, ref steps, out _));
        Assert.Throws<ArgumentOutOfRangeException>(() => t.Reader.CompareSuffixes(0, 1, -1, ref steps, out _));
    }

    // ---- bulk access

    [Test]
    public void ReadSymbols_EveryRange_EqualsTheOracle_ForEveryWidth()
    {
        foreach (int bits in AllBits)
        {
            byte[] s = TextOracle.Random(70, 1 << bits, 31 + (ulong)bits);
            using var t = new PinnedText(s, bits);
            for (int start = 0; start <= s.Length; start++)
            {
                for (int len = 0; start + len <= s.Length; len++)
                {
                    var got = new byte[len];
                    t.Reader.ReadSymbols(start, got);
                    Assert.That(got, Is.EqualTo(s[start..(start + len)]), $"b={bits} [{start},{start + len})");
                }
            }
        }
    }

    [Test]
    public void ReadSymbols_LargeRandomWindows_AndTheWholeText()
    {
        foreach (int bits in FirstClass)
        {
            int n = 1_000_003;
            byte[] s = TextOracle.Random(n, 1 << bits, 55 + (ulong)bits);
            using var t = new PinnedText(s, bits);
            var whole = new byte[n];
            t.Reader.ReadSymbols(0, whole);
            Assert.That(whole, Is.EqualTo(s), $"b={bits} the whole text");
            var rng = new SplitMix64(9);
            for (int round = 0; round < 300; round++)
            {
                int start = rng.NextInt(n);
                int len = rng.NextInt(Math.Min(5000, n - start) + 1);
                var got = new byte[len];
                t.Reader.ReadSymbols(start, got);
                Assert.That(got.AsSpan().SequenceEqual(s.AsSpan(start, len)), Is.True, $"b={bits} [{start},{start + len})");
            }
        }
    }

    [Test]
    public void ReadSymbols_RefusesRangesOutsideTheText()
    {
        using var t = new PinnedText(TextOracle.Random(20, 4, 1), 2);
        Assert.Throws<ArgumentOutOfRangeException>(() => t.Reader.ReadSymbols(-1, new byte[3]));
        Assert.Throws<ArgumentOutOfRangeException>(() => t.Reader.ReadSymbols(19, new byte[2]));
        Assert.Throws<ArgumentOutOfRangeException>(() => t.Reader.ReadSymbols(21, []));
        Assert.Throws<ArgumentOutOfRangeException>(() => t.Reader.ReadSymbols(long.MaxValue, new byte[1]));
        t.Reader.ReadSymbols(20, []);
        t.Reader.ReadSymbols(19, new byte[1]);
    }

    // ---- common prefix: vector path, scalar path

    [Test]
    public void Lce_ScalarAndVector_EqualTheNaiveCount_OnRepeatRichTexts()
    {
        foreach (int bits in AllBits)
        {
            int n = 3000;
            byte[] s = TextOracle.Random(n, 1 << bits, 1000 + (ulong)bits, repeatEvery: 211);
            using var t = new PinnedText(s, bits);
            var rng = new SplitMix64(5 + (ulong)bits);
            for (int round = 0; round < 800; round++)
            {
                int i = rng.NextInt(n);
                int j = rng.NextInt(n);
                if (round % 5 == 0) j = Math.Min(n - 1, i + (211 * (1 + rng.NextInt(3))));
                if (round % 7 == 0) i = n - 1 - rng.NextInt(40);
                long max = n - Math.Max(i, j);
                long cap = rng.NextInt(4) == 0 ? max : Math.Min(max, 1 + rng.NextInt(900));
                long want = TextOracle.Lce(s, i, j, cap);
                foreach (bool vectors in new[] { false, true })
                {
                    t.Reader.UseVectors = vectors;
                    long steps = 0;
                    Assert.That(t.Reader.Lce(i, j, cap, ref steps), Is.EqualTo(want), $"b={bits} i={i} j={j} cap={cap} vectors={vectors}");
                }
            }
        }
    }

    [Test]
    public void Lce_ClampsToTheTextAndNeverReadsPastIt()
    {
        foreach (int bits in new[] { 2, 3, 8 })
        {
            byte[] s = new byte[500]; // all the same: every comparison runs to the end
            using var t = new PinnedText(s, bits);
            foreach (bool vectors in new[] { false, true })
            {
                t.Reader.UseVectors = vectors;
                long steps = 0;
                Assert.That(t.Reader.Lce(0, 1, 1_000_000, ref steps), Is.EqualTo(499), $"b={bits} the cap is clamped to n - max(i, j)");
                Assert.That(t.Reader.Lce(10, 490, long.MaxValue, ref steps), Is.EqualTo(10));
                Assert.That(t.Reader.Lce(500, 3, 100, ref steps), Is.Zero, "i = n: the empty suffix");
                Assert.That(t.Reader.Lce(3, 501, 100, ref steps), Is.Zero, "j beyond n");
                Assert.That(t.Reader.Lce(1L << 40, 1L << 41, 100, ref steps), Is.Zero);
                Assert.That(t.Reader.Lce(0, 0, 500, ref steps), Is.EqualTo(500), "a suffix with itself");
                Assert.That(t.Reader.Lce(0, 1, 0, ref steps), Is.Zero);
            }
        }
    }

    [Test]
    public void VectorPath_TakesFewerSteps_AndTheFlagIsGuardedByTheHardware()
    {
        var s = new byte[100_000];
        using var t = new PinnedText(s, 2);
        Assert.That(t.Reader.UseVectors, Is.EqualTo(Avx2.IsSupported), "the default");
        t.Reader.UseVectors = true;
        Assert.That(t.Reader.UseVectors, Is.EqualTo(Avx2.IsSupported), "never on without AVX2");
        long vector = 0;
        long a = t.Reader.Lce(0, 1, 99_999, ref vector);
        t.Reader.UseVectors = false;
        Assert.That(t.Reader.UseVectors, Is.False);
        long scalar = 0;
        long b = t.Reader.Lce(0, 1, 99_999, ref scalar);
        Assert.That(a, Is.EqualTo(99_999));
        Assert.That(b, Is.EqualTo(a), "identical results on both paths (the scalar one is the one that runs on arm64)");
        Assert.That(scalar, Is.EqualTo((99_999 + 31) / 32), "the scalar path compares 64 / b = 32 symbols per step");
        if (Avx2.IsSupported) Assert.That(vector, Is.EqualTo((99_999 / 128) + 1), "the vector path compares 256 / b = 128 symbols per step, and one scalar step for the last 31 symbols");
        else Assert.That(vector, Is.EqualTo(scalar));
    }

    [Test]
    public void TheStepCounts_AreThoseOfTheCostModel_AtTheEdgesOfTheVectorLoop()
    {
        var s = new byte[2_000];
        using var t = new PinnedText(s, 2);
        foreach (bool vectors in new[] { false, true })
        {
            t.Reader.UseVectors = vectors;
            bool v = t.Reader.UseVectors;
            foreach ((long max, long scalarSteps, long vectorSteps) in new[] { (0L, 0L, 0L), (1L, 1L, 1L), (31L, 1L, 1L), (32L, 1L, 1L), (33L, 2L, 2L), (127L, 4L, 4L), (128L, 4L, 1L), (129L, 5L, 2L), (256L, 8L, 2L), (1_000L, 32L, 11L) })
            {
                long steps = 0;
                Assert.That(t.Reader.Lce(0, 1, max, ref steps), Is.EqualTo(max));
                Assert.That(steps, Is.EqualTo(v ? vectorSteps : scalarSteps), $"max={max} vectors={v}");
            }
        }
    }

    [Test]
    public void CompareSuffixes_PositionsAtOrBeyondTheEnd_AreTheEmptySuffix()
    {
        byte[] s = TextOracle.Random(50, 4, 3);
        using var t = new PinnedText(s, 2);
        long steps = 0;
        Assert.That(t.Reader.CompareSuffixes(50, 3, 0, ref steps, out long lce), Is.EqualTo(-1), "the empty suffix is the smallest");
        Assert.That(lce, Is.Zero);
        Assert.That(t.Reader.CompareSuffixes(3, 55, 0, ref steps, out lce), Is.EqualTo(1));
        Assert.That(lce, Is.Zero);
        Assert.That(t.Reader.CompareSuffixes(50, 90, 0, ref steps, out lce), Is.Zero, "two empty suffixes are the same");
        Assert.That(t.Reader.CompareSuffixes(1L << 40, 60, 5, ref steps, out lce), Is.Zero);
        Assert.That(lce, Is.Zero);
        Assert.That(t.Reader.CompareSuffixes(0, 50, 10, ref steps, out lce), Is.EqualTo(1));
        Assert.That(lce, Is.Zero, "nothing is common with the empty suffix, whatever skip says");
    }

    [Test]
    public void CompareSuffixes_EqualsTheNaiveOrder_WithSkip_AndReportsTheCommonPrefix()
    {
        foreach (int bits in new[] { 1, 2, 3, 8 })
        {
            int n = 400;
            byte[] s = TextOracle.Random(n, 1 << bits, 31 + (ulong)bits, repeatEvery: 17);
            using var t = new PinnedText(s, bits);
            foreach (bool vectors in new[] { false, true })
            {
                t.Reader.UseVectors = vectors;
                for (int i = 0; i <= n; i += 7)
                {
                    for (int j = 0; j <= n; j += 11)
                    {
                        if (i == j) continue;
                        long steps = 0;
                        int got = t.Reader.CompareSuffixes(i, j, 0, ref steps, out long lce);
                        Assert.That(Math.Sign(got), Is.EqualTo(TextOracle.CompareSuffixes(s, i, j)), $"b={bits} i={i} j={j}");
                        Assert.That(lce, Is.EqualTo(TextOracle.Lce(s, i, j, n)), $"b={bits} lce({i},{j})");
                        long common = TextOracle.Lce(s, i, j, n);
                        int skipped = t.Reader.CompareSuffixes(i, j, common / 2, ref steps, out long lce2);
                        Assert.That(Math.Sign(skipped), Is.EqualTo(Math.Sign(got)), $"skip: i={i} j={j}");
                        Assert.That(lce2, Is.EqualTo(lce), "the skipped symbols are counted in the result");
                        int big = t.Reader.CompareSuffixes(i, j, 1_000_000, ref steps, out long lce3);
                        Assert.That(lce3, Is.EqualTo(n - Math.Max(i, j)), "a skip beyond the shorter suffix is clamped to it");
                        Assert.That(big, Is.EqualTo(i > j ? -1 : 1), "the shorter suffix is the smaller one");
                    }
                }
                long st = 0;
                Assert.That(t.Reader.CompareSuffixes(5, 5, 0, ref st, out long same), Is.Zero, "a suffix with itself");
                Assert.That(same, Is.EqualTo(n - 5));
            }
        }
    }

    // ---- larger texts, random positions, concurrency

    [Test]
    public void RandomPositions_OnLargeTexts_EqualTheOracle()
    {
        foreach (int bits in AllBits)
        {
            int n = 2_000_003;
            byte[] s = TextOracle.Random(n, 1 << bits, 600 + (ulong)bits, repeatEvery: 100_003);
            using var t = new PinnedText(s, bits);
            var rng = new SplitMix64(17 + (ulong)bits);
            for (int round = 0; round < 20_000; round++)
            {
                long i = rng.NextInt(n + 70);
                Assert.That(t.Reader.Key(i), Is.EqualTo(TextOracle.Key(s, i, bits)), $"b={bits} Key({i})");
                Assert.That(t.Reader.Symbol(i), Is.EqualTo(i < n ? (uint)s[i] : 0u), $"b={bits} Symbol({i})");
                if (round % 4 == 0) Assert.That(t.Reader.Window(i), Is.EqualTo(TextOracle.Window(s, i, bits)), $"b={bits} Window({i})");
            }
            for (int round = 0; round < 2_000; round++)
            {
                int i = rng.NextInt(n);
                int j = round % 2 == 0 ? Math.Min(n - 1, i + 100_003) : rng.NextInt(n);
                long cap = n - Math.Max(i, j);
                long want = TextOracle.Lce(s, i, j, cap);
                foreach (bool vectors in new[] { false, true })
                {
                    t.Reader.UseVectors = vectors;
                    long steps = 0;
                    Assert.That(t.Reader.Lce(i, j, cap, ref steps), Is.EqualTo(want), $"b={bits} Lce({i},{j}) vectors={vectors}");
                }
            }
        }
    }

    [Test]
    public void ConcurrentReaders_AgreeWithTheOracle()
    {
        int n = 1_000_003;
        byte[] s = TextOracle.Random(n, 8, 3);
        using var t = new PinnedText(s, 3);
        var failures = new List<string>();
        Parallel.For(0, 8, new ParallelOptions { MaxDegreeOfParallelism = 8 }, worker =>
        {
            var rng = new SplitMix64((ulong)worker + 1);
            for (int round = 0; round < 60_000; round++)
            {
                long i = rng.NextInt(n + 30);
                if (t.Reader.Key(i) != TextOracle.Key(s, i, 3))
                {
                    lock (failures) failures.Add($"worker {worker} Key({i})");
                    return;
                }
            }
        });
        Assert.That(failures, Is.Empty);
    }

    // ---- through files: the whole-file open, a section of a larger file, sparse files above 2^31 and 2^32

    [Test]
    public void OpenedFromAFile_ReadsTheSame_AndOpenRegionFindsASectionInsideALargerFile()
    {
        using var dir = new TempDir();
        foreach (int bits in FirstClass)
        {
            byte[] s = TextOracle.Random(5003, 1 << bits, 70 + (ulong)bits);
            byte[] packed = TextOracle.Pack(s, bits);
            string whole = dir.File($"w{bits}.pk");
            File.WriteAllBytes(whole, packed);
            using (MappedTextReader r = MappedTextReader.Open(whole, s.Length, bits))
            {
                for (int i = 0; i < s.Length + 3; i++) Assert.That(r.Key(i), Is.EqualTo(TextOracle.Key(s, i, bits)), $"b={bits} Key({i})");
            }
            string container = dir.File($"c{bits}.bin");
            const int front = 4096 + 64;
            var file = new byte[front + packed.Length + 4096];
            Array.Fill(file, (byte)0xEE); // whatever surrounds the section must not matter
            packed.CopyTo(file, front);
            File.WriteAllBytes(container, file);
            using MappedTextReader sec = MappedTextReader.OpenRegion(container, front, packed.Length, s.Length, bits);
            for (int i = 0; i < s.Length + 70; i++) Assert.That(sec.Key(i), Is.EqualTo(TextOracle.Key(s, i, bits)), $"b={bits} section Key({i})");
            var all = new byte[s.Length];
            sec.ReadSymbols(0, all);
            Assert.That(all, Is.EqualTo(s));
        }
    }

    [Test]
    public void SparseFiles_PositionsAboveTwoToThe31AndTwoToThe32_AreReadThroughTheMapping()
    {
        using var dir = new TempDir();
        (int Bits, long N)[] cases = [(2, 4_400_000_000L), (3, 3_000_000_000L), (4, 4_300_000_000L), (8, 2_300_000_000L)];
        foreach ((int bits, long n) in cases)
        {
            string path = dir.File($"big{bits}.pk");
            long data = PackedTextFormat.DataBytes(n, bits);
            using (var fs = new FileStream(path, FileMode.Create, FileAccess.Write)) fs.SetLength(PackedTextFormat.TotalBytes(n, bits));
            var planted = new SortedDictionary<long, byte>();
            var rng = new SplitMix64((ulong)bits);
            foreach (long centre in new[] { (1L << 31) - 3, 1L << 31, (1L << 32) - 3, 1L << 32, 2_150_000_000L, n - 40 })
            {
                if (centre + 33 > n) continue;
                for (int k = 0; k < 33; k++)
                {
                    byte sym = (byte)(1 + rng.NextInt((1 << bits) - 1));
                    planted[centre + k] = sym;
                    TextOracle.Plant(path, 0, bits, centre + k, sym);
                }
            }
            long Sym(long p) => planted.TryGetValue(p, out byte v) ? v : (byte)0;
            using MappedTextReader r = MappedTextReader.Open(path, n, bits);
            Assert.That(r.Length, Is.EqualTo(n));
            foreach (long p in planted.Keys)
            {
                Assert.That(r.Symbol(p), Is.EqualTo((uint)Sym(p)), $"b={bits} Symbol({p})");
                int ks = 64 / bits;
                ulong key = 0;
                for (int k = 0; k < ks; k++) key = (key << bits) | (p + k < n ? (ulong)Sym(p + k) : 0UL);
                Assert.That(r.Key(p), Is.EqualTo(key), $"b={bits} Key({p})");
                Assert.That(r.Window(p) >> (64 - bits), Is.EqualTo((ulong)Sym(p)));
            }
            long lo = planted.Keys.First(p => p >= (1L << 31) - 3);
            var buf = new byte[40];
            r.ReadSymbols(lo - 3, buf);
            for (int k = 0; k < buf.Length; k++) Assert.That(buf[k], Is.EqualTo(Sym(lo - 3 + k)), $"b={bits} ReadSymbols at {lo - 3 + k}");
            Assert.That(r.Symbol(n), Is.Zero);
            Assert.That(r.Key(n + 5), Is.Zero);
            Assert.That(data, Is.GreaterThan(1L << 30));
            long steps = 0;
            long a = (1L << 31) + 0, b = (1L << 31) + 100;
            Assert.That(r.Lce(a, b, 20, ref steps), Is.EqualTo(TextOracle_SparseLce(Sym, a, b, 20)));
        }
    }

    private static long TextOracle_SparseLce(Func<long, long> sym, long a, long b, long max)
    {
        long k = 0;
        while (k < max && sym(a + k) == sym(b + k)) k++;
        return k;
    }

    [Test]
    public void ASection_AtAFileOffsetAboveTwoToThe32_IsOpenedAndRead()
    {
        using var dir = new TempDir();
        string path = dir.File("far.bin");
        const long offset = 5_000_000_000L;
        byte[] s = TextOracle.Random(2_000, 8, 4);
        byte[] packed = TextOracle.Pack(s, 3);
        using (var fs = new FileStream(path, FileMode.Create, FileAccess.Write))
        {
            fs.SetLength(offset + packed.Length + 4096);
            fs.Position = offset;
            fs.Write(packed);
        }
        using MappedTextReader r = MappedTextReader.OpenRegion(path, offset, packed.Length, s.Length, 3);
        for (int i = 0; i < s.Length + 3; i++) Assert.That(r.Key(i), Is.EqualTo(TextOracle.Key(s, i, 3)), $"Key({i})");
        Assert.Throws<InvalidDataException>(() => MappedTextReader.OpenRegion(path, offset + 4097, packed.Length, s.Length, 3), "a region that starts later runs past the end of the file");
    }

    // ---- verification and damaged input

    [Test]
    public void Verify_AcceptsAnIntactText_ForEveryWidthAndSigma()
    {
        foreach (int bits in AllBits)
        {
            foreach (int sigma in new[] { (1 << (bits - 1)) + 1, 1 << bits })
            {
                byte[] s = TextOracle.Random(20_000, sigma, (ulong)((bits * 100) + sigma));
                using var t = new PinnedText(s, bits);
                PackedTextInfo info = TextOracle.Info(s, bits, sigma);
                t.Reader.Verify(info);
                t.Reader.Verify(info, checkHash: false);
            }
        }
        using var empty = new PinnedText([], 1);
        empty.Reader.Verify(new PackedTextInfo(0, 1, [], [], Convert.ToHexStringLower(System.Security.Cryptography.SHA256.HashData([]))));
    }

    [Test]
    public void Verify_ASymbolAboveSigma_IsReported()
    {
        byte[] s = TextOracle.Random(5_000, 7, 1);
        s[4_321] = 7;
        using var t = new PinnedText(s, 3);
        PackedTextInfo info = TextOracle.Info(s, 3, 8) with { Letters = [.. Enumerable.Range('A', 7).Select(x => (byte)x)], Counts = new long[7] };
        InvalidDataException? e = Assert.Throws<InvalidDataException>(() => t.Reader.Verify(info));
        Assert.That(e!.Message, Does.Contain("symbol").And.Contain("7").And.Contain("4321").IgnoreCase);
    }

    [Test]
    public void Verify_ACountThatDiffers_ALengthThatDiffers_AWidthThatDiffers_AHashThatDiffers_AreReported()
    {
        byte[] s = TextOracle.Random(5_000, 8, 2);
        using var t = new PinnedText(s, 3);
        PackedTextInfo good = TextOracle.Info(s, 3, 8);
        long[] counts = [.. good.Counts];
        counts[2]++;
        counts[3]--;
        InvalidDataException? e = Assert.Throws<InvalidDataException>(() => t.Reader.Verify(good with { Counts = counts }));
        Assert.That(e!.Message, Does.Contain("count").IgnoreCase);
        e = Assert.Throws<InvalidDataException>(() => t.Reader.Verify(good with { N = 4_999 }));
        Assert.That(e!.Message, Does.Contain("length").IgnoreCase);
        e = Assert.Throws<InvalidDataException>(() => t.Reader.Verify(good with { Bits = 4 }));
        Assert.That(e!.Message, Does.Contain("bits").IgnoreCase);
        e = Assert.Throws<InvalidDataException>(() => t.Reader.Verify(good with { Sha256 = new string('0', 64) }));
        Assert.That(e!.Message, Does.Contain("SHA-256").IgnoreCase);
        t.Reader.Verify(good with { Sha256 = new string('0', 64) }, checkHash: false);
        t.Reader.Verify(good with { Sha256 = good.Sha256.ToUpperInvariant() });
        Assert.Throws<InvalidDataException>(() => t.Reader.Verify(good with { Letters = good.Letters[..7] }), "the alphabet and the counts must agree");
        e = Assert.Throws<InvalidDataException>(() => t.Reader.Verify(good with { Counts = [.. good.Counts, 0L] }), "one count too many, every symbol below sigma");
        Assert.That(e!.Message, Does.Contain("letters").And.Contain("counts"));
    }

    [Test]
    public void AWrongWidth_IsCaught_ByTheLengthOrByTheSymbolCounts()
    {
        using var dir = new TempDir();
        foreach (int trueBits in new[] { 2, 3, 4, 8 })
        {
            foreach (int openBits in new[] { 1, 2, 3, 4, 5, 8 })
            {
                if (openBits == trueBits) continue;
                foreach (int n in new[] { 1, 5, 8, 13, 100, 1001 })
                {
                    byte[] s = TextOracle.Random(n, 1 << trueBits, (ulong)(trueBits * 31 + openBits * 7 + n));
                    s[0] = (byte)((1 << trueBits) - 1); // not all zero: a text of zeros is valid at every width of the same size
                    string path = TextOracle.WriteFile(dir, $"x{trueBits}_{openBits}_{n}.pk", s, trueBits);
                    bool caught = false;
                    try
                    {
                        using MappedTextReader r = MappedTextReader.Open(path, n, openBits);
                        r.Verify(TextOracle.Info(s, openBits, 1 << trueBits));
                    }
                    catch (InvalidDataException)
                    {
                        caught = true;
                    }
                    Assert.That(caught, Is.True, $"file packed with b={trueBits} opened as b={openBits}, n={n}");
                }
            }
        }
    }

    [Test]
    public void DamagedFiles_AreRefusedWithAClearException_NeverReadOutOfTheMapping()
    {
        using var dir = new TempDir();
        byte[] s = TextOracle.Random(1001, 8, 5);
        byte[] packed = TextOracle.Pack(s, 3); // 1001 * 3 = 3003 bits: 5 unused bits in the last data byte
        long data = PackedTextFormat.DataBytes(1001, 3);
        string Make(string name, byte[] bytes)
        {
            string p = dir.File(name);
            File.WriteAllBytes(p, bytes);
            return p;
        }

        // truncated by one byte, by the whole padding, to the data only, to nothing, one byte
        foreach (int keep in new[] { packed.Length - 1, packed.Length - 64, (int)data, 1, 0 })
        {
            string p = Make($"t{keep}.pk", packed[..keep]);
            InvalidDataException? e = Assert.Throws<InvalidDataException>(() => MappedTextReader.Open(p, 1001, 3).Dispose(), $"kept {keep}");
            Assert.That(e!.Message, Does.Contain("truncated").Or.Contain("length").IgnoreCase, $"kept {keep}");
        }

        // longer than the layout says: a whole-file open wants the exact length
        string longer = Make("longer.pk", [.. packed, 0]);
        Assert.Throws<InvalidDataException>(() => MappedTextReader.Open(longer, 1001, 3).Dispose());
        MappedTextReader.OpenRegion(longer, 0, packed.Length, 1001, 3).Dispose(); // the same bytes as a section of a larger file are fine
        Assert.Throws<InvalidDataException>(() => MappedTextReader.OpenRegion(longer, 0, packed.Length + 1, 1001, 3).Dispose(), "a section whose LENGTH differs from the layout");
        Assert.Throws<InvalidDataException>(() => MappedTextReader.OpenRegion(longer, 0, packed.Length - 1, 1001, 3).Dispose());
        Assert.Throws<InvalidDataException>(() => MappedTextReader.OpenRegion(longer, 2, packed.Length, 1001, 3).Dispose(), "past the end of the file");
        Assert.Throws<InvalidDataException>(() => MappedTextReader.OpenRegion(longer, -1, packed.Length, 1001, 3).Dispose());

        // unused bits of the last data byte
        byte[] unused = (byte[])packed.Clone();
        unused[data - 1] |= 0x01;
        InvalidDataException? bad = Assert.Throws<InvalidDataException>(() => MappedTextReader.Open(Make("unused.pk", unused), 1001, 3).Dispose());
        Assert.That(bad!.Message, Does.Contain("unused").IgnoreCase);

        // padding bytes
        foreach (int k in new[] { 0, 1, 31, 62, 63 })
        {
            byte[] pad = (byte[])packed.Clone();
            pad[data + k] = 1;
            bad = Assert.Throws<InvalidDataException>(() => MappedTextReader.Open(Make($"pad{k}.pk", pad), 1001, 3).Dispose(), $"padding byte {k}");
            Assert.That(bad!.Message, Does.Contain("padding").IgnoreCase);
        }

        // widths and lengths that cannot be right
        string zeros = Make("zeros.pk", new byte[64]);
        Assert.Throws<InvalidDataException>(() => MappedTextReader.Open(zeros, 5, 0).Dispose(), "a file of 64 bytes would fit zero bits per symbol: the width is refused, not divided by");
        Assert.Throws<InvalidDataException>(() => MappedTextReader.Open(zeros, 5, -1).Dispose());
        foreach (int bits in new[] { 0, -1, 9, 16, 33, int.MinValue })
        {
            Assert.Throws<InvalidDataException>(() => MappedTextReader.Open(Make("b.pk", packed), 1001, bits).Dispose(), $"bits {bits}");
        }
        foreach (long n in new[] { -1L, long.MinValue, long.MaxValue, (long.MaxValue / 8) + 1, 1L << 61 })
        {
            Assert.Throws<InvalidDataException>(() => MappedTextReader.Open(Make("n.pk", packed), n, 3).Dispose(), $"n {n}");
        }

        // not a file, no file
        Assert.Throws<FileNotFoundException>(() => MappedTextReader.Open(dir.File("missing.pk"), 0, 1).Dispose());
        Assert.Throws<DirectoryNotFoundException>(() => MappedTextReader.Open(dir.File("no/such/dir.pk"), 0, 1).Dispose());
        Assert.That(() => MappedTextReader.Open(dir.Path, 0, 1).Dispose(), Throws.InstanceOf<IOException>().Or.InstanceOf<UnauthorizedAccessException>());
        Assert.Throws<ArgumentNullException>(() => MappedTextReader.Open(null!, 0, 1).Dispose());
    }

    [Test]
    public void FromMemory_ChecksTheLengthAndThePadding_LikeTheFileOpen()
    {
        byte[] s = TextOracle.Random(100, 8, 5);
        byte[] packed = TextOracle.Pack(s, 3);
        using (var ok = new PinnedText(packed, 100, 3)) Assert.That(ok.Reader.Length, Is.EqualTo(100));
        Assert.Throws<InvalidDataException>(() => new PinnedText(packed[..^1], 100, 3).Dispose());
        InvalidDataException? longer = Assert.Throws<InvalidDataException>(() => new PinnedText([.. packed, 0], 100, 3).Dispose());
        Assert.That(longer!.Message, Does.Contain("length"));
        byte[] pad = (byte[])packed.Clone();
        pad[^1] = 7;
        Assert.Throws<InvalidDataException>(() => new PinnedText(pad, 100, 3).Dispose());
        byte[] unused = (byte[])packed.Clone();
        unused[(int)PackedTextFormat.DataBytes(100, 3) - 1] |= 1;
        Assert.Throws<InvalidDataException>(() => new PinnedText(unused, 100, 3).Dispose());
    }

    [Test]
    public void AfterDispose_EveryReadThrowsObjectDisposed_NothingIsRead_AndDisposeIsIdempotent()
    {
        using var dir = new TempDir();
        byte[] s = TextOracle.Random(100, 4, 5);
        MappedTextReader r = MappedTextReader.Open(TextOracle.WriteFile(dir, "d.pk", s, 2), 100, 2);
        Assert.That(r.Key(3), Is.EqualTo(TextOracle.Key(s, 3, 2)));
        bool? mapped = ProcMaps.Contains(dir.File("d.pk"));
        if (mapped is not null) Assert.That(mapped, Is.True, "the file is mapped while the reader is open");
        r.Dispose();
        r.Dispose();
        if (mapped is not null) Assert.That(ProcMaps.Contains(dir.File("d.pk")), Is.False, "Dispose unmaps the view");
        long steps = 0;
        Assert.Throws<ObjectDisposedException>(() => r.Symbol(3));
        Assert.Throws<ObjectDisposedException>(() => r.Key(3));
        Assert.Throws<ObjectDisposedException>(() => r.Window(3));
        Assert.Throws<ObjectDisposedException>(() => r.Symbol(1000), "even past the end");
        Assert.Throws<ObjectDisposedException>(() => r.ReadSymbols(0, new byte[3]));
        Assert.Throws<ObjectDisposedException>(() => r.Lce(0, 1, 5, ref steps));
        Assert.Throws<ObjectDisposedException>(() => r.CompareSuffixes(0, 1, 0, ref steps, out _));
        Assert.Throws<ObjectDisposedException>(() => r.Verify(TextOracle.Info(s, 2, 4)));
        Assert.That(r.Length, Is.EqualTo(100), "the facts stay readable");
        Assert.That(r.Bits, Is.EqualTo(2));
        using var exclusive = new FileStream(dir.File("d.pk"), FileMode.Open, FileAccess.ReadWrite, FileShare.None);
    }

    [Test]
    public void AnEmptyText_HasAnEmptyDataPartAndTheZeroPadding()
    {
        using var dir = new TempDir();
        string p = dir.File("e.pk");
        File.WriteAllBytes(p, new byte[64]);
        using MappedTextReader r = MappedTextReader.Open(p, 0, 1);
        Assert.That(r.Length, Is.Zero);
        Assert.That(r.Key(0), Is.Zero);
        Assert.That(r.Symbol(0), Is.Zero);
        r.ReadSymbols(0, []);
        long steps = 0;
        Assert.That(r.Lce(0, 0, 0, ref steps), Is.Zero);
    }
}
