using SuffixTree.Esa.Tests.External;
using SuffixTree.Esa.Tests.Fuzz;
using SuffixTree.Esa.Tests.Large;

namespace SuffixTree.Esa.Tests.Views;

/// <summary>
/// The LCP pyramid of work package S2.3 with 64-bit ranks and bounds, over every form of <see cref="CompactLcp"/> (managed, mapped at width 4, mapped at width 5, U8 / U16 / U32 / U64): every query against brute force, the 32-bit
/// members that the navigation pack calls are the same queries, and a block minimum that the 32-bit levels cannot hold is <see cref="NotSupportedException"/> (a bridge to be removed by S2.5), never a wrapped number.
/// </summary>
[TestFixture]
public unsafe class LcpPyramidWideTests
{
    private static long[] Values(int n, int seed, int maxSmall, int bigEvery, long big)
    {
        var v = new long[n];
        ulong x = (ulong)seed * 0x9E3779B97F4A7C15UL;
        for (int i = 0; i < n; i++)
        {
            x = (x ^ (x >> 29)) * 0xBF58476D1CE4E5B9UL + 1;
            v[i] = i < 2 ? 0 : (long)((x >> 33) % (ulong)maxSmall);
            if (bigEvery > 0 && i > 2 && i % bigEvery == 0) v[i] = big + (i % 11);
        }
        return v;
    }

    private static long BruteMin(long[] v, long lo, long hi) => v.Skip((int)lo).Take((int)(hi - lo + 1)).Min();

    private static long BrutePrev(long[] v, long pos, long bound)
    {
        for (long j = pos; j >= 0; j--) if (v[j] < bound) return j;
        return -1;
    }

    private static long BruteNext(long[] v, long pos, long bound)
    {
        for (long j = pos; j < v.Length; j++) if (v[j] < bound) return j;
        return v.Length;
    }

    private static void CheckAll(LcpPyramid pyramid, long[] v, string context)
    {
        Assert.That(pyramid.Length, Is.EqualTo(v.Length), context);
        long[] bounds = [0, 1, 2, 3, 7, 20, 254, 255, 256, 300, 65_534, 65_535, 65_536, 70_000, 1L << 31, 1L << 33];
        for (long p = 0; p < v.Length; p += 17)
        {
            Assert.That(pyramid.ValueAt(p), Is.EqualTo(v[p]), $"{context} ValueAt({p})");
            foreach (long bound in bounds)
            {
                Assert.That(pyramid.PrevLess(p, bound), Is.EqualTo(BrutePrev(v, p, bound)), $"{context} PrevLess({p}, {bound})");
                Assert.That(pyramid.NextLess(p, bound), Is.EqualTo(BruteNext(v, p, bound)), $"{context} NextLess({p}, {bound})");
            }
        }
        Assert.That(pyramid.PrevLess(-1, 5), Is.EqualTo(-1));
        Assert.That(pyramid.NextLess(v.Length, 5), Is.EqualTo(v.Length));
        Assert.That(pyramid.NextLess(long.MaxValue, 5), Is.EqualTo(v.Length), "a position beyond the array");
        for (long lo = 0; lo < v.Length; lo += 211)
        {
            foreach (long len in new long[] { 1, 2, 31, 32, 33, 500, 2000 })
            {
                long hi = Math.Min(v.Length - 1, lo + len - 1);
                Assert.That(pyramid.Min(lo, hi), Is.EqualTo(BruteMin(v, lo, hi)), $"{context} Min({lo}, {hi})");
            }
        }
        Assert.That(pyramid.Min(0, v.Length - 1), Is.EqualTo(BruteMin(v, 0, v.Length - 1)));
        Assert.That(() => pyramid.ValueAt(v.Length), Throws.TypeOf<ArgumentOutOfRangeException>());
        Assert.That(() => pyramid.ValueAt(-1L), Throws.TypeOf<ArgumentOutOfRangeException>());
        Assert.That(() => pyramid.ValueAt(1L << 40), Throws.TypeOf<ArgumentOutOfRangeException>());
        Assert.That(() => pyramid.Min(5, 4), Throws.TypeOf<ArgumentOutOfRangeException>());
        Assert.That(() => pyramid.Min(0, v.Length), Throws.TypeOf<ArgumentOutOfRangeException>());
        Assert.That(() => pyramid.Min(-1L, 3), Throws.TypeOf<ArgumentOutOfRangeException>());
        Assert.That(() => pyramid.PrevLess(v.Length, 3), Throws.TypeOf<ArgumentOutOfRangeException>());
        Assert.That(() => pyramid.NextLess(-1, 3), Throws.TypeOf<ArgumentOutOfRangeException>());

        // the borders of a range: the previous and the next entry below the bound, and their values (-1 for no entry)
        foreach (long a in new long[] { 0, 1, 2, 17, 500, v.Length - 1 })
        {
            foreach (long len in new long[] { 1, 5, 40 })
            {
                long b = Math.Min(v.Length - 1, a + len - 1);
                foreach (long t in new long[] { 0, 1, 3, 7, 20, 300, 1L << 33 })
                {
                    pyramid.Borders(a, b, t, out long left, out long right, out long leftValue, out long rightValue);
                    long expectedLeft = BrutePrev(v, a, t), expectedRight = BruteNext(v, b + 1, t);
                    Assert.That((left, right), Is.EqualTo((expectedLeft, expectedRight)), $"{context} Borders({a}, {b}, {t})");
                    Assert.That(leftValue, Is.EqualTo(expectedLeft >= 0 ? v[expectedLeft] : -1), $"{context} Borders({a}, {b}, {t}) left value");
                    Assert.That(rightValue, Is.EqualTo(expectedRight < v.Length ? v[expectedRight] : -1), $"{context} Borders({a}, {b}, {t}) right value");
                }
            }
        }
        // the 32-bit members are the same queries, and -1 (no entry) is an answer, not a damaged number
        foreach (int bound in new[] { 0, 1, 3, 7, 20, 300, 70_000 })
        {
            for (int p = 0; p < v.Length; p += 29)
            {
                Assert.That(pyramid.PrevLess(p, bound), Is.EqualTo((int)BrutePrev(v, p, bound)), $"{context} int PrevLess({p}, {bound})");
                Assert.That(pyramid.NextLess(p, bound), Is.EqualTo((int)BruteNext(v, p, bound)), $"{context} int NextLess({p}, {bound})");
            }
            foreach (int a in new[] { 0, 3, 100, v.Length - 1 })
            {
                int b = Math.Min(v.Length - 1, a + 6);
                pyramid.Borders(a, b, bound, out int l, out int r, out int lv, out int rv);
                int el = (int)BrutePrev(v, a, bound), er = (int)BruteNext(v, b + 1, bound);
                Assert.That((l, r), Is.EqualTo((el, er)), $"{context} int Borders({a}, {b}, {bound})");
                Assert.That(lv, Is.EqualTo(el >= 0 ? (int)v[el] : -1));
                Assert.That(rv, Is.EqualTo(er < v.Length ? (int)v[er] : -1));
            }
        }
        Assert.That(pyramid.PrevLess(0, 0), Is.EqualTo(-1), "no entry below 0");
        Assert.That(() => pyramid.Borders(5L, 4L, 3L, out _, out _, out _, out _), Throws.TypeOf<ArgumentOutOfRangeException>());
        Assert.That(() => pyramid.Borders(5L, v.Length, 3L, out _, out _, out _, out _), Throws.TypeOf<ArgumentOutOfRangeException>());
        Assert.That(() => pyramid.Borders(v.Length, v.Length, 3L, out _, out _, out _, out _), Throws.TypeOf<ArgumentOutOfRangeException>());
    }

    [Test]
    public void EveryFormOfTheLcpArray_GivesTheBruteForceAnswers([Values("U8", "U16", "U32", "U64")] string kindName, [Values(4, 5)] int width, [Values(true, false)] bool isChecked)
    {
        var kind = Enum.Parse<LcpKind>(kindName);
        long big = kind == LcpKind.U8 ? 300 : kind == LcpKind.U16 ? 70_000 : 2_500; // a U32 / U64 value above the text is damage to the Checked reader
        long[] v = Values(3_000, 1, 40, 53, big);
        (byte[] baseBytes, List<(long Rank, long Value)> pairs) = ViewImage.EncodeLcp(v, kind);
        byte[] pairBytes = ViewImage.Encode(pairs.SelectMany(p => new[] { (ulong)p.Rank, (ulong)p.Value }), width, width == 5 ? 3 : 0);
        using var b = new GuardedImage(baseBytes, flush: true);
        using var p = new GuardedImage(pairBytes, flush: true);
        CompactLcp lcp = CompactLcp.FromMappedWide(kind, v.Length, b.Data, pairs.Count == 0 ? null : p.Data, pairs.Count, width, isChecked);
        using LcpPyramid pyramid = LcpPyramid.Build(lcp, 32);
        CheckAll(pyramid, v, $"{kind} w{width} checked={isChecked}");
    }

    [Test]
    public void TheManagedAndTheLegacyMappedForms_GiveTheSameAnswers_ThroughTheLongMembers([Values("U8", "U16", "U32")] string kindName, [Values(2, 4, 32)] int fanout)
    {
        var kind = Enum.Parse<LcpKind>(kindName);
        long big = kind == LcpKind.U8 ? 400 : 80_000;
        long[] v = Values(2_500, 2, 30, 61, big);
        if (kind == LcpKind.U32) v = [.. v.Select(x => Math.Min(x, 2_000))];
        using LcpPyramid pyramid = LcpPyramid.Build(CompactLcp.FromDense([.. v.Select(x => (int)x)], kind), fanout);
        CheckAll(pyramid, v, $"managed {kind} F{fanout}");
        (byte[] baseBytes, List<(long Rank, long Value)> pairs) = ViewImage.EncodeLcp(v, kind);
        uint[] words = [.. pairs.SelectMany(p => new[] { (uint)p.Rank, (uint)p.Value })];
        fixed (byte* bp = baseBytes)
        fixed (uint* wp = words)
        {
            CompactLcp mapped = CompactLcp.FromMapped(kind, v.Length, bp, words.Length == 0 ? null : wp, pairs.Count, isChecked: false);
            using LcpPyramid legacy = LcpPyramid.Build(mapped, fanout);
            CheckAll(legacy, v, $"mapped {kind} F{fanout}");
        }
    }

    [Test]
    public void TheIntegerMembers_AreTheSameQueries()
    {
        long[] v = Values(1_500, 3, 25, 97, 400);
        using LcpPyramid pyramid = LcpPyramid.Build(CompactLcp.FromDense([.. v.Select(x => (int)x)]), 8);
        for (int p = 0; p < v.Length; p += 13)
        {
            Assert.That(pyramid.ValueAt(p), Is.EqualTo((int)v[p]));
            Assert.That(pyramid.PrevLess(p, 6), Is.EqualTo((int)BrutePrev(v, p, 6)));
            Assert.That(pyramid.NextLess(p, 6), Is.EqualTo((int)BruteNext(v, p, 6)));
            Assert.That(pyramid.PreviousSmaller(p), Is.EqualTo((int)BrutePrev(v, p - 1, v[p])));
            Assert.That(pyramid.NextSmaller(p), Is.EqualTo((int)BruteNext(v, p + 1, v[p])));
        }
        pyramid.Borders(500, 520, 5, out int left, out int right, out int leftValue, out int rightValue);
        Assert.That((left, right), Is.EqualTo(((int)BrutePrev(v, 500, 5), (int)BruteNext(v, 521, 5))));
        Assert.That(leftValue, Is.EqualTo(left >= 0 ? (int)v[left] : -1));
        Assert.That(rightValue, Is.EqualTo(right < v.Length ? (int)v[right] : -1));
        pyramid.Borders(500L, 520L, 5L, out long l2, out long r2, out long lv2, out long rv2);
        Assert.That((l2, r2, lv2, rv2), Is.EqualTo((left, right, leftValue, rightValue)));
    }

    [Test]
    public void ABlockMinimumThatThe32BitLevelsCannotHold_IsNotSupported_AndAHugeValueThatIsNoMinimumIsFine()
    {
        // U64: every value above int.MaxValue: the minimum of the first block is 2^31 + 5: the levels are int
        long[] v = [.. Enumerable.Range(0, 200).Select(i => i < 2 ? 0L : (1L << 31) + 5 + i)];
        (byte[] baseBytes, _) = ViewImage.EncodeLcp(v, LcpKind.U64);
        using var b = new GuardedImage(baseBytes, flush: true);
        CompactLcp all = CompactLcp.FromMappedWide(LcpKind.U64, v.Length, b.Data, null, 0, 5, isChecked: false);
        Assert.That(() => LcpPyramid.Build(all, 4), Throws.TypeOf<NotSupportedException>().With.Message.Contains("32-bit"));
        // one huge entry per block among small ones: the minima are small, the huge value is read through the base reader
        long[] w = [.. Enumerable.Range(0, 200).Select(i => i % 4 == 3 ? (1L << 33) + i : i % 5)];
        (byte[] baseW, _) = ViewImage.EncodeLcp(w, LcpKind.U64);
        using var bw = new GuardedImage(baseW, flush: true);
        CompactLcp some = CompactLcp.FromMappedWide(LcpKind.U64, w.Length, bw.Data, null, 0, 5, isChecked: false);
        using LcpPyramid pyramid = LcpPyramid.Build(some, 4);
        Assert.That(pyramid.ValueAt(3L), Is.EqualTo((1L << 33) + 3));
        Assert.That(() => pyramid.ValueAt(3), Throws.TypeOf<InvalidDataException>(), "the int member says that the value is no int");
        Assert.That(pyramid.Min(0, 199), Is.EqualTo(0));
    }

    [Test]
    public void ALcpArrayWhoseLevelsHoldMoreThanInt32MinimaButLessThanUInt32_IsRefusedByTheCountOfTheLevels_BeforeAnythingIsAllocated()
    {
        // Orchestrator verification of S2.3: a^n at 2^33 with fanout 4 has 2^31 + 2^29 ... = about 2.86 x 10^9 block minima, between int.MaxValue and uint.MaxValue; the 32-bit table of levels cannot index them and the pyramid says so
        // by name (a limit of 4 x 10^9 instead of 2 x 10^9 would let the count wrap to a negative int and fail later and elsewhere). The check is on sizes: nothing is allocated and no value is read, so the test is instant.
        var f = new PowerOfA(1L << 33);
        CompactLcp lcp = CompactLcp.FromReaders(LcpKind.U64, new PosGen<W8>(f, GenSection.LcpBase), new PosGen<W5>(f, GenSection.LcpPairsOfPowerOfA, long.MaxValue), isChecked: false);
        Assert.That(() => LcpPyramid.Build(lcp, 4), Throws.TypeOf<NotSupportedException>().With.Message.Contains("block minima"));
    }

    [Test]
    public void AnEscapeThatTheListDoesNotHold_IsInvalidData_ThroughThePyramid()
    {
        long[] v = [.. Enumerable.Range(0, 120).Select(i => i == 60 ? 700L : 2L)];
        (byte[] baseBytes, _) = ViewImage.EncodeLcp(v, LcpKind.U8);
        using var b = new GuardedImage(baseBytes, flush: true);
        // the base array has the escape at rank 60 and the list is empty: the file is damaged
        Assert.That(() => LcpPyramid.Build(CompactLcp.FromMappedWide(LcpKind.U8, v.Length, b.Data, null, 0, 5, true), 4), Throws.TypeOf<InvalidDataException>());
    }

    [Test]
    public void ALoneMinimum_IsFoundAtEveryPositionOfTheBlocks_ByEveryQuery([Values(2, 4, 8)] int fanout, [Values(4, 5)] int width)
    {
        // one entry of 1 among entries of 100, the 1 at every position in turn: a range holds it or it does not, whatever the block borders, so a partial block that is read one entry short or one entry late, a level that
        // is entered at the wrong block, shows as the wrong minimum. Every range, every position, both directions.
        const int n = 70;
        for (int p = 0; p < n; p++)
        {
            long[] v = [.. Enumerable.Range(0, n).Select(i => i == p ? 1L : 100L)];
            (byte[] baseBytes, _) = ViewImage.EncodeLcp(v, LcpKind.U8);
            using var b = new GuardedImage(baseBytes, flush: true);
            using LcpPyramid pyramid = LcpPyramid.Build(CompactLcp.FromMappedWide(LcpKind.U8, n, b.Data, null, 0, width, isChecked: true), fanout);
            for (long lo = 0; lo < n; lo++)
            {
                for (long hi = lo; hi < n; hi++)
                {
                    long expected = lo <= p && p <= hi ? 1 : 100;
                    long actual = pyramid.Min(lo, hi);
                    if (actual != expected) Assert.Fail($"Min({lo}, {hi}) with the minimum at {p}, fanout {fanout}: {actual}, expected {expected}");
                }
                long prev = pyramid.PrevLess(lo, 50);
                if (prev != (lo >= p ? p : -1)) Assert.Fail($"PrevLess({lo}, 50) with the minimum at {p}, fanout {fanout}: {prev}");
                long next = pyramid.NextLess(lo, 50);
                if (next != (lo <= p ? p : n)) Assert.Fail($"NextLess({lo}, 50) with the minimum at {p}, fanout {fanout}: {next}");
            }
        }
    }

    [Test]
    public void TheManagedPyramid_RefusesARankThatIsNoInt_WithTheMessageOfThePyramid()
    {
        using LcpPyramid pyramid = LcpPyramid.Build(CompactLcp.FromDense([0, 0, 1, 2, 0, 3, 1, 1]), 4);
        Assert.That(() => pyramid.ValueAt(int.MaxValue + 1L), Throws.TypeOf<ArgumentOutOfRangeException>().With.Message.Contains("a rank outside the LCP array"));
        Assert.That(() => pyramid.ValueAt(-1L), Throws.TypeOf<ArgumentOutOfRangeException>().With.Message.Contains("a rank outside the LCP array"));
        Assert.That(() => pyramid.Min(int.MaxValue + 1L, 5), Throws.TypeOf<ArgumentOutOfRangeException>().With.Message.Contains("a rank outside the LCP array"));
        Assert.That(() => pyramid.Min(0L, 8L), Throws.TypeOf<ArgumentOutOfRangeException>().With.Property("ParamName").EqualTo("hi"), "the right end at the length");
        Assert.That(pyramid.ValueAt(5L), Is.EqualTo(3));
        Assert.That(pyramid.Min(2L, 5L), Is.EqualTo(0));
        Assert.That(pyramid.Min(0L, 7L), Is.EqualTo(0));
    }

    [Test]
    public void ABlockMinimumOfExactlyIntMaxValue_FitsTheLevels()
    {
        long[] v = [.. Enumerable.Repeat((long)int.MaxValue, 200)];
        (byte[] baseBytes, _) = ViewImage.EncodeLcp(v, LcpKind.U64);
        using var b = new GuardedImage(baseBytes, flush: true);
        using LcpPyramid pyramid = LcpPyramid.Build(CompactLcp.FromMappedWide(LcpKind.U64, v.Length, b.Data, null, 0, 4, isChecked: false), 4);
        Assert.That(pyramid.Min(0, 199), Is.EqualTo(int.MaxValue));
        Assert.That(pyramid.Min(3, 150), Is.EqualTo(int.MaxValue));
        Assert.That(pyramid.PrevLess(199L, int.MaxValue), Is.EqualTo(-1));
        Assert.That(pyramid.PrevLess(199, int.MaxValue), Is.EqualTo(-1), "no entry is an answer, not a damaged number");
        Assert.That(pyramid.NextLess(0L, (long)int.MaxValue + 1), Is.EqualTo(0));
    }

    [Test]
    public void ALengthWhoseBlockMinimaDoNotFitTheLevels_IsNotSupported_BeforeAnythingIsRead()
    {
        // 2^32 + 10 entries in blocks of 2: more than 2^31 block minima at the first level; the check precedes the allocation and the scan (the generator holds no memory, a scan would not end)
        var f = new PowerOfA((1L << 32) + 10);
        CompactLcp lcp = CompactLcp.FromReaders(LcpKind.U32, new PosGen<W4>(f, GenSection.Synthetic), new PosGen<W4>(f, GenSection.LcpPairsOfPowerOfA, long.MaxValue), isChecked: false);
        Assert.That(() => LcpPyramid.Build(lcp, 2), Throws.TypeOf<NotSupportedException>().With.Message.Contains("block minima"));
    }
}
