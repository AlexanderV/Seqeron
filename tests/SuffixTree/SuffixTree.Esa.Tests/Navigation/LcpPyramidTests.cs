using SuffixTree.Esa.Bench.Links;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.Navigation;

/// <summary>
/// The production range-minimum / previous-smaller / next-smaller helper (work package 1b.4). Oracles: brute force over the dense array; the prototype <see cref="LcpRangeHelper"/> of 1b.4a as a
/// second implementation. Arrays are stored in all three representations of <see cref="CompactLcp"/> with the values that matter for the escape-aware pointer reader (254, 255, 256, 65,534,
/// 65,535, 65,536): a bound above the escape needs the true value of an escaped entry, a bound at or below it does not.
/// </summary>
[TestFixture]
public class LcpPyramidTests
{
    private static int BruteMin(int[] v, int lo, int hi)
    {
        int m = int.MaxValue;
        for (int j = lo; j <= hi; j++) m = Math.Min(m, v[j]);
        return m;
    }

    private static int BrutePrevLess(int[] v, int pos, int bound)
    {
        for (int j = pos; j >= 0; j--)
            if (v[j] < bound) return j;
        return -1;
    }

    private static int BruteNextLess(int[] v, int pos, int bound)
    {
        for (int j = pos; j < v.Length; j++)
            if (v[j] < bound) return j;
        return v.Length;
    }

    private static void Fail(string context, string what, int actual, int expected) => Assert.Fail($"{context}: {what} gave {actual}, expected {expected}");

    /// <summary>Every query on one array against brute force (and the prototype when it is given); bounds are the values around every distinct value of the array. Plain comparisons: the exhaustive sweeps run millions of queries.</summary>
    internal static void CheckAll(int[] v, LcpKind kind, int fanout, bool prototype, string context)
    {
        CompactLcp compact = CompactLcp.FromDense(v, kind);
        using LcpPyramid pyramid = LcpPyramid.Build(compact, fanout);
        LcpRangeHelper? proto = prototype ? new LcpRangeHelper(compact, fanout) : null;
        string ctx = $"{context} {kind} F{fanout} [{(v.Length <= 12 ? string.Join(",", v) : v.Length + " entries")}]";
        if (pyramid.Length != v.Length) Fail(ctx, "Length", pyramid.Length, v.Length);
        if (pyramid.Fanout != fanout) Fail(ctx, "Fanout", pyramid.Fanout, fanout);
        var bounds = new SortedSet<int> { 0, 1, 2, int.MaxValue };
        foreach (int x in v.Distinct())
        {
            bounds.Add(x);
            bounds.Add(x + 1);
            if (x > 0) bounds.Add(x - 1);
        }
        // long arrays: a spread of bounds and positions (the brute force is O(n) per query)
        int[] boundList = bounds.ToArray();
        if (boundList.Length > 24) boundList = Enumerable.Range(0, 24).Select(i => boundList[i * (boundList.Length - 1) / 23]).Distinct().ToArray();
        int posStep = Math.Max(1, v.Length / 48);
        for (int p = 0; p < v.Length; p++)
        {
            if (v.Length > 96 && p % posStep != 0 && p >= 3 && p < v.Length - 3 && (p & 511) != 0 && (p & 511) != 511) continue;
            int got = pyramid.ValueAt(p);
            if (got != v[p]) Fail(ctx, $"ValueAt({p})", got, v[p]);
            foreach (int bound in boundList)
            {
                got = pyramid.PrevLess(p, bound);
                int want = BrutePrevLess(v, p, bound);
                if (got != want) Fail(ctx, $"PrevLess({p}, {bound})", got, want);
                if (proto is not null && got != proto.PrevLess(p, bound)) Fail(ctx, $"PrevLess({p}, {bound}) against the prototype", got, proto.PrevLess(p, bound));
                got = pyramid.NextLess(p, bound);
                want = BruteNextLess(v, p, bound);
                if (got != want) Fail(ctx, $"NextLess({p}, {bound})", got, want);
                if (proto is not null && got != proto.NextLess(p, bound)) Fail(ctx, $"NextLess({p}, {bound}) against the prototype", got, proto.NextLess(p, bound));
            }
            got = pyramid.PreviousSmaller(p);
            if (got != BrutePrevLess(v, p - 1, v[p])) Fail(ctx, $"PreviousSmaller({p})", got, BrutePrevLess(v, p - 1, v[p]));
            got = pyramid.NextSmaller(p);
            if (got != BruteNextLess(v, p + 1, v[p])) Fail(ctx, $"NextSmaller({p})", got, BruteNextLess(v, p + 1, v[p]));
        }
        int step = v.Length > 60 ? v.Length / 30 : 1;
        for (int lo = 0; lo < v.Length; lo += step)
            for (int hi = lo; hi < v.Length; hi += step)
            {
                int got = pyramid.Min(lo, hi);
                int want = BruteMin(v, lo, hi);
                if (got != want) Fail(ctx, $"Min({lo}, {hi})", got, want);
                if (proto is not null && got != proto.Min(lo, hi)) Fail(ctx, $"Min({lo}, {hi}) against the prototype", got, proto.Min(lo, hi));
            }
        if (v.Length > 0 && pyramid.Min(0, v.Length - 1) != BruteMin(v, 0, v.Length - 1)) Fail(ctx, "Min(whole)", pyramid.Min(0, v.Length - 1), BruteMin(v, 0, v.Length - 1));
    }

    internal static IEnumerable<int[]> AllArrays(int[] palette, int maxLength)
    {
        for (int len = 1; len <= maxLength; len++)
        {
            var idx = new int[len];
            while (true)
            {
                yield return idx.Select(i => palette[i]).ToArray();
                int pos = len - 1;
                while (pos >= 0 && idx[pos] == palette.Length - 1) idx[pos--] = 0;
                if (pos < 0) break;
                idx[pos]++;
            }
        }
    }

    [Test]
    public void AllBinaryArrays_UpTo10_AllKindsAndFanouts()
    {
        foreach (int[] v in AllArrays([0, 1], 10))
            foreach (int fanout in new[] { 2, 4, 8 })
                CheckAll(v, LcpKind.U8, fanout, false, "binary");
    }

    [Test]
    public void AllTernaryArrays_UpTo7_WithPrototype()
    {
        foreach (int[] v in AllArrays([0, 1, 2], 7))
            CheckAll(v, LcpKind.U32, 2, true, "ternary");
    }

    [Test]
    public void EscapeBoundaryValues_U8_AllArraysUpTo4()
    {
        foreach (int[] v in AllArrays([0, 1, 254, 255, 256, 300], 4))
            foreach (int fanout in new[] { 2, 4 })
                CheckAll(v, LcpKind.U8, fanout, false, "u8 escapes");
    }

    [Test]
    public void EscapeBoundaryValues_U16_AllArraysUpTo4()
    {
        foreach (int[] v in AllArrays([0, 1, 65534, 65535, 65536, 70000], 4))
            foreach (int fanout in new[] { 2, 4 })
                CheckAll(v, LcpKind.U16, fanout, false, "u16 escapes");
    }

    [Test]
    public void EscapeBoundaryValues_AllKindsHoldAnyValue()
    {
        foreach (int[] v in AllArrays([0, 3, 255, 65535, 100000], 4))
            foreach (LcpKind kind in Enum.GetValues<LcpKind>())
                CheckAll(v, kind, 2, true, "mixed");
    }

    [Test]
    public void RandomArrays_AllKinds_AllFanouts_WithPrototype()
    {
        var rng = new SplitMix64(0x1B4B10);
        int[] ranges = [2, 5, 40, 300, 70000];
        for (int rep = 0; rep < 60; rep++)
        {
            int len = 1 + rng.NextInt(rep < 40 ? 120 : 1500);
            int range = ranges[rng.NextInt(ranges.Length)];
            var v = new int[len];
            for (int i = 0; i < len; i++) v[i] = rng.NextInt(range);
            // sprinkle escape-relevant values
            for (int k = 0; k < len / 10; k++) v[rng.NextInt(len)] = new[] { 254, 255, 256, 65534, 65535, 65536 }[rng.NextInt(6)];
            LcpKind kind = Enum.GetValues<LcpKind>()[rep % 3];
            int fanout = new[] { 2, 4, 16, 32 }[rng.NextInt(4)];
            CheckAll(v, kind, fanout, len <= 400, $"random #{rep} len {len} range {range}");
        }
    }

    [Test]
    public void StructuredArrays()
    {
        foreach (int n in new[] { 1, 2, 31, 32, 33, 64, 65, 1023, 1024, 1025, 2000 })
        {
            var flat = new int[n];
            Array.Fill(flat, 7);
            var up = Enumerable.Range(0, n).ToArray();
            var down = Enumerable.Range(0, n).Reverse().ToArray();
            var saw = Enumerable.Range(0, n).Select(i => i % 13).ToArray();
            var dip = Enumerable.Repeat(9, n).ToArray();
            dip[n / 2] = 0;
            var bump = new int[n];
            bump[n / 3] = 500;
            foreach (var (name, v) in new[] { ("flat", flat), ("up", up), ("down", down), ("saw", saw), ("dip", dip), ("bump", bump) })
                foreach (int fanout in new[] { 2, 32 })
                    CheckAll(v, LcpKind.U8, fanout, false, $"{name} n {n}");
        }
    }

    [Test]
    public void LargeArray_FourLevels_AgainstBruteForceOnSamples()
    {
        var rng = new SplitMix64(0x1B4B11);
        const int Length = 40_000;
        var v = new int[Length];
        for (int i = 0; i < Length; i++) v[i] = 30 + rng.NextInt(1000);
        foreach (int i in new[] { 5, 777, 12_345, 31_999, 39_999 }) v[i] = rng.NextInt(30);
        foreach (LcpKind kind in Enum.GetValues<LcpKind>())
        {
            using LcpPyramid pyramid = LcpPyramid.Build(CompactLcp.FromDense(v, kind));
            Assert.That(pyramid.Fanout, Is.EqualTo(32));
            Assert.That(pyramid.Levels, Is.EqualTo(4), "40,000 entries: 1,250 / 40 / 2 / 1 blocks");
            for (int k = 0; k < 400; k++)
            {
                int p = rng.NextInt(Length), bound = rng.NextInt(1100);
                Assert.That(pyramid.PrevLess(p, bound), Is.EqualTo(BrutePrevLess(v, p, bound)), $"{kind} PrevLess({p}, {bound})");
                Assert.That(pyramid.NextLess(p, bound), Is.EqualTo(BruteNextLess(v, p, bound)), $"{kind} NextLess({p}, {bound})");
                int lo = rng.NextInt(Length), hi = Math.Min(Length - 1, lo + rng.NextInt(Length));
                Assert.That(pyramid.Min(lo, hi), Is.EqualTo(BruteMin(v, lo, hi)), $"{kind} Min({lo}, {hi})");
            }
        }
    }

    [Test]
    public void Shape_LevelsAndBytes()
    {
        static LcpPyramid Of(int length, int fanout)
            => LcpPyramid.Build(CompactLcp.FromDense(new int[length]), fanout);

        using (LcpPyramid p = Of(1, 32)) Assert.That((p.Levels, p.LevelBytes), Is.EqualTo((0, 0L)), "one entry: nothing to summarise");
        using (LcpPyramid p = Of(2, 32)) Assert.That((p.Levels, p.LevelBytes), Is.EqualTo((1, 4L)));
        using (LcpPyramid p = Of(32, 32)) Assert.That((p.Levels, p.LevelBytes), Is.EqualTo((1, 4L)));
        using (LcpPyramid p = Of(33, 32)) Assert.That((p.Levels, p.LevelBytes), Is.EqualTo((2, 12L)), "2 blocks, then 1");
        using (LcpPyramid p = Of(1025, 32)) Assert.That((p.Levels, p.LevelBytes), Is.EqualTo((3, 4L * (33 + 2 + 1))));
        using (LcpPyramid p = Of(8, 2)) Assert.That((p.Levels, p.LevelBytes), Is.EqualTo((3, 4L * (4 + 2 + 1))));
        using LcpPyramid big = Of(1_000_001, 32);
        Assert.That(big.LevelBytes / 1_000_001.0, Is.InRange(0.12, 0.14), "about 4 / (F - 1) bytes per entry");
        Assert.That(LcpPyramid.DefaultFanout, Is.EqualTo(32));
    }

    [Test]
    public void Contract_ArgumentChecks()
    {
        CompactLcp lcp = CompactLcp.FromDense(new[] { 0, 0, 1, 2, 0 });
        Assert.Throws<ArgumentNullException>(() => LcpPyramid.Build(null!));
        foreach (int bad in new[] { -2, 0, 1, 3, 6, 33 })
            Assert.Throws<ArgumentOutOfRangeException>(() => LcpPyramid.Build(lcp, bad), $"fanout {bad}");
        using LcpPyramid p = LcpPyramid.Build(lcp, 2);
        Assert.That(p.PrevLess(-1, 5), Is.EqualTo(-1), "nothing at or before a negative position");
        Assert.That(p.NextLess(5, 5), Is.EqualTo(5), "nothing at or after the end");
        Assert.Throws<ArgumentOutOfRangeException>(() => p.PrevLess(5, 1));
        Assert.Throws<ArgumentOutOfRangeException>(() => p.NextLess(-1, 1));
        Assert.Throws<ArgumentOutOfRangeException>(() => p.Min(3, 2));
        Assert.Throws<ArgumentOutOfRangeException>(() => p.Min(-1, 2));
        Assert.Throws<ArgumentOutOfRangeException>(() => p.Min(0, 5));
        Assert.Throws<ArgumentOutOfRangeException>(() => p.ValueAt(5));
        Assert.Throws<ArgumentOutOfRangeException>(() => p.ValueAt(-1));
        Assert.DoesNotThrow(p.Dispose);
        Assert.DoesNotThrow(p.Dispose);
    }

    [Test]
    public void LcpOfARealIndex_FibonacciAndRepeats_AgainstPrototype()
    {
        foreach (string text in new[] { StringFamilies.Fibonacci(3000), StringFamilies.AllSame(1200), StringFamilies.ThueMorse(2000), StringFamilies.RandomDna(4000, 2) })
        {
            EsaIndex index = EsaIndex.Build(text);
            int[] dense = index.Lcp.ToArray();
            foreach (int fanout in new[] { 4, 32 })
                CheckAll(dense.Length > 1500 ? dense[..1500] : dense, index.Lcp.Kind, fanout, false, $"real lcp {text.Length}");
            using LcpPyramid pyramid = LcpPyramid.Build(index.Lcp);
            var proto = new LcpRangeHelper(index.Lcp, 32);
            var rng = new SplitMix64(11);
            for (int k = 0; k < 300; k++)
            {
                int p = rng.NextInt(dense.Length), bound = rng.NextInt(Math.Max(2, dense.Max() + 2));
                Assert.That(pyramid.PrevLess(p, bound), Is.EqualTo(proto.PrevLess(p, bound)));
                Assert.That(pyramid.NextLess(p, bound), Is.EqualTo(proto.NextLess(p, bound)));
            }
        }
    }
}


/// <summary>The long enumerations of <see cref="LcpPyramidTests"/>: every array over the escape-relevant palettes up to length 6 and binary arrays up to length 13. Heavy; explicit, category <c>Exhaustive</c>.</summary>
[TestFixture]
[Category("Exhaustive")]
[Explicit("Heavy enumeration; run with --filter Category=Exhaustive (see PROGRESS-S1.md, Log 1b.4)")]
public class LcpPyramidExhaustiveTests
{
    [Test]
    public void AllBinaryArrays_UpTo13()
    {
        foreach (int[] v in LcpPyramidTests.AllArrays([0, 1], 13))
            foreach (int fanout in new[] { 2, 4, 8 })
                LcpPyramidTests.CheckAll(v, LcpKind.U8, fanout, false, "binary");
    }

    [Test]
    public void AllTernaryArrays_UpTo9_Fanout2()
    {
        foreach (int[] v in LcpPyramidTests.AllArrays([0, 1, 2], 9))
            LcpPyramidTests.CheckAll(v, LcpKind.U32, 2, false, "ternary");
    }

    [Test]
    public void EscapePalettes_U8_AllArraysUpTo6()
    {
        foreach (int[] v in LcpPyramidTests.AllArrays([0, 1, 254, 255, 256, 300], 6))
            LcpPyramidTests.CheckAll(v, LcpKind.U8, 2, false, "u8 escapes");
    }

    [Test]
    public void EscapePalettes_U16_AllArraysUpTo6()
    {
        foreach (int[] v in LcpPyramidTests.AllArrays([0, 1, 65534, 65535, 65536, 70000], 6))
            LcpPyramidTests.CheckAll(v, LcpKind.U16, 2, false, "u16 escapes");
    }
}
