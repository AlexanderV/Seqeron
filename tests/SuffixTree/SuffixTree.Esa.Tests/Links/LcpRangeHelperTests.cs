using SuffixTree.Esa.Bench.Links;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.Links;

/// <summary>
/// The RMQ + PSV/NSV helper of link kernel (c) (work package 1b.4a, prototype P1) against brute force on the dense array: range minimum with the leftmost position,
/// previous / next smaller value (strict) and the value-parameterised PrevLess / NextLess, on random, all-equal, strictly monotone, saw-tooth and escape-boundary arrays,
/// for fanouts 2, 4, 16 and 64 (small fanouts make a short array a deep pyramid).
/// </summary>
[TestFixture]
public class LcpRangeHelperTests
{
    private static readonly int[] Fanouts = [2, 4, 16, 64];

    private static int BruteMin(int[] a, int lo, int hi)
    {
        int m = int.MaxValue;
        for (int i = lo; i <= hi; i++) m = Math.Min(m, a[i]);
        return m;
    }

    private static int BruteArgMin(int[] a, int lo, int hi)
    {
        int best = lo;
        for (int i = lo + 1; i <= hi; i++)
            if (a[i] < a[best]) best = i;
        return best;
    }

    private static int BrutePrevLess(int[] a, int pos, int bound)
    {
        for (int j = pos; j >= 0; j--)
            if (a[j] < bound) return j;
        return -1;
    }

    private static int BruteNextLess(int[] a, int pos, int bound)
    {
        for (int j = pos; j < a.Length; j++)
            if (a[j] < bound) return j;
        return a.Length;
    }

    /// <summary>Checks every query of the helper on <paramref name="a"/>: all (lo, hi) pairs when the array is short, a seeded sample otherwise; all positions and the interesting bounds.</summary>
    private static void CheckArray(int[] a, string name, bool allPairs = true, int sampledPairs = 3000, ulong seed = 1)
    {
        foreach (LcpKind kind in new[] { LcpKind.U32, LcpKind.U8, LcpKind.U16 })
        {
            CompactLcp lcp = CompactLcp.FromDense(a, kind);
            foreach (int fanout in Fanouts)
            {
                var h = new LcpRangeHelper(lcp, fanout);
                if (h.Length != a.Length) Assert.Fail($"{name}: Length {h.Length} != {a.Length} (fanout {fanout})");
                string ctx = $"{name} kind {kind} fanout {fanout}";

                if (allPairs)
                {
                    for (int lo = 0; lo < a.Length; lo++)
                        for (int hi = lo; hi < a.Length; hi++)
                            CheckRange(h, a, lo, hi, ctx);
                }
                else
                {
                    var rng = new SplitMix64(seed ^ (ulong)fanout);
                    for (int k = 0; k < sampledPairs; k++)
                    {
                        int lo = rng.NextInt(a.Length);
                        int hi = lo + rng.NextInt(Math.Min(a.Length - lo, k % 3 == 0 ? 40 : a.Length) );
                        CheckRange(h, a, lo, hi, ctx);
                    }
                    CheckRange(h, a, 0, a.Length - 1, ctx);
                }

                int[] bounds = a.Distinct().Concat([0, 1, int.MaxValue - 1]).Concat(a.Select(v => v + 1)).Distinct().OrderBy(v => v).Take(60).ToArray();
                IEnumerable<int> positions = allPairs ? Enumerable.Range(0, a.Length) : Enumerable.Range(0, 400).Select(i => (int)((long)i * (a.Length - 1) / 399));
                foreach (int pos in positions)
                {
                    foreach (int bound in bounds)
                    {
                        int pl = h.PrevLess(pos, bound);
                        if (pl != BrutePrevLess(a, pos, bound)) Assert.Fail($"{ctx}: PrevLess({pos}, {bound}) = {pl}, expected {BrutePrevLess(a, pos, bound)}");
                        int nl = h.NextLess(pos, bound);
                        if (nl != BruteNextLess(a, pos, bound)) Assert.Fail($"{ctx}: NextLess({pos}, {bound}) = {nl}, expected {BruteNextLess(a, pos, bound)}");
                    }
                    int psv = h.PreviousSmaller(pos);
                    if (psv != BrutePrevLess(a, pos - 1, a[pos])) Assert.Fail($"{ctx}: PreviousSmaller({pos}) = {psv}, expected {BrutePrevLess(a, pos - 1, a[pos])}");
                    int nsv = h.NextSmaller(pos);
                    if (nsv != BruteNextLess(a, pos + 1, a[pos])) Assert.Fail($"{ctx}: NextSmaller({pos}) = {nsv}, expected {BruteNextLess(a, pos + 1, a[pos])}");
                }
            }
        }
    }

    private static void CheckRange(LcpRangeHelper h, int[] a, int lo, int hi, string ctx)
    {
        int min = h.Min(lo, hi);
        if (min != BruteMin(a, lo, hi)) Assert.Fail($"{ctx}: Min({lo}, {hi}) = {min}, expected {BruteMin(a, lo, hi)}");
        int arg = h.ArgMin(lo, hi);
        if (arg != BruteArgMin(a, lo, hi)) Assert.Fail($"{ctx}: ArgMin({lo}, {hi}) = {arg}, expected {BruteArgMin(a, lo, hi)}");
    }

    [Test]
    public void AllArraysOverThreeValues_UpToLength8()
    {
        for (int len = 1; len <= 8; len++)
        {
            int total = (int)Math.Pow(3, len);
            for (int code = 0; code < total; code++)
            {
                var a = new int[len];
                for (int i = 0, c = code; i < len; i++, c /= 3) a[i] = c % 3;
                CheckArray(a, $"ternary-{len}-{code}");
            }
        }
    }

    [Test]
    public void AllBinaryArrays_UpToLength14()
    {
        for (int len = 9; len <= 14; len++)
            for (int code = 0; code < 1 << len; code++)
            {
                var a = new int[len];
                for (int i = 0; i < len; i++) a[i] = (code >> i) & 1;
                CheckArray(a, $"binary-{len}-{code}", allPairs: len <= 11, sampledPairs: 500);
            }
    }

    [Test]
    public void RandomArrays_ManyValueRanges()
    {
        var rng = new SplitMix64(0xBEEF);
        foreach (int range in new[] { 2, 3, 5, 16, 100, 300, 70000 })
            for (int rep = 0; rep < 12; rep++)
            {
                int len = 1 + rng.NextInt(range <= 5 ? 300 : 150);
                var a = new int[len];
                for (int i = 0; i < len; i++) a[i] = rng.NextInt(range);
                CheckArray(a, $"random-r{range}-{rep}-n{len}", allPairs: len <= 90, sampledPairs: 1500, seed: (ulong)rep);
            }
    }

    [Test]
    public void AllEqual_StrictlyMonotone_SawTooth_SingleDip()
    {
        foreach (int len in new[] { 1, 2, 3, 15, 16, 17, 63, 64, 65, 255, 256, 257, 1000 })
        {
            CheckArray(Enumerable.Repeat(7, len).ToArray(), $"all-equal-{len}", allPairs: len <= 70, sampledPairs: 1200);
            CheckArray(Enumerable.Range(0, len).ToArray(), $"increasing-{len}", allPairs: len <= 70, sampledPairs: 1200);
            CheckArray(Enumerable.Range(0, len).Reverse().ToArray(), $"decreasing-{len}", allPairs: len <= 70, sampledPairs: 1200);
            CheckArray(Enumerable.Range(0, len).Select(i => i % 5).ToArray(), $"saw5-{len}", allPairs: len <= 70, sampledPairs: 1200);
            CheckArray(Enumerable.Range(0, len).Select(i => i % 17 * 31).ToArray(), $"saw17-{len}", allPairs: len <= 70, sampledPairs: 1200);
            CheckArray(Enumerable.Range(0, len).Select(i => i == len / 2 ? 0 : 5).ToArray(), $"dip-{len}", allPairs: len <= 70, sampledPairs: 1200);
            CheckArray(Enumerable.Range(0, len).Select(i => i == len / 2 ? 9 : 5).ToArray(), $"bump-{len}", allPairs: len <= 70, sampledPairs: 1200);
            CheckArray(Enumerable.Range(0, len).Select(i => (i & 1) == 0 ? 3 : 4).ToArray(), $"zigzag-{len}", allPairs: len <= 70, sampledPairs: 1200);
        }
    }

    [Test]
    public void ValuesAroundTheEscapeBoundaries()
    {
        foreach (int centre in new[] { 255, 65535 })
        {
            var rng = new SplitMix64((ulong)centre);
            for (int rep = 0; rep < 6; rep++)
            {
                int len = 40 + rng.NextInt(120);
                var a = new int[len];
                for (int i = 0; i < len; i++) a[i] = centre - 2 + rng.NextInt(5);
                CheckArray(a, $"escape-{centre}-{rep}", allPairs: false, sampledPairs: 800, seed: (ulong)rep);
            }
        }
    }

    [Test]
    public void LargeArrays_SampledQueries()
    {
        var rng = new SplitMix64(0x1B4A);
        foreach (int len in new[] { 5000, 40000 })
        {
            var a = new int[len];
            for (int i = 0; i < len; i++) a[i] = rng.NextInt(1000) < 995 ? rng.NextInt(30) : 100 + rng.NextInt(2000);
            CheckArray(a, $"large-mixed-{len}", allPairs: false, sampledPairs: 2500);
            CheckArray(Enumerable.Range(0, len).Select(i => i).ToArray(), $"large-increasing-{len}", allPairs: false, sampledPairs: 2500);
            CheckArray(new int[len], $"large-zero-{len}", allPairs: false, sampledPairs: 2500);
        }
    }

    [Test]
    public void PayloadBytes_AreFourBytesPerBlockMinimumOfEveryLevel()
    {
        var lcp = CompactLcp.FromDense(new int[1000 + 1]);
        foreach (int fanout in Fanouts)
        {
            var h = new LcpRangeHelper(lcp, fanout);
            long expected = 0;
            int size = lcp.Length;
            int levels = 0;
            while (size > 1)
            {
                size = (size + fanout - 1) / fanout;
                expected += 4L * size;
                levels++;
            }
            Assert.That(h.PayloadBytes, Is.EqualTo(expected), $"fanout {fanout}");
            Assert.That(h.Levels, Is.EqualTo(levels), $"fanout {fanout}");
            Assert.That(h.Fanout, Is.EqualTo(fanout));
        }
    }

    [Test]
    public void Fanout_MustBeAPowerOfTwoAtLeastTwo()
    {
        var lcp = CompactLcp.FromDense(new[] { 0, 1, 2 });
        foreach (int bad in new[] { 0, 1, 3, 6, -4 })
            Assert.Throws<ArgumentOutOfRangeException>(() => _ = new LcpRangeHelper(lcp, bad), $"fanout {bad}");
        Assert.Throws<ArgumentNullException>(() => _ = new LcpRangeHelper(null!, 16));
    }

    [Test]
    public void StepCounting_CountsEveryEntryRead_AndZeroForNoCount()
    {
        var a = new int[2000];
        for (int i = 0; i < a.Length; i++) a[i] = 10 + (i % 7);
        a[3] = 0;
        var h = new LcpRangeHelper(CompactLcp.FromDense(a), 16);
        var counted = new CountSteps();
        Assert.That(h.PrevLess(1900, 1, ref counted), Is.EqualTo(3));
        Assert.That(counted.Steps, Is.GreaterThan(0));
        Assert.That(counted.Steps, Is.LessThan(200), "the pyramid skips the 1,900 entries in between");
        var none = new NoCount();
        Assert.That(h.NextLess(4, 1, ref none), Is.EqualTo(a.Length));
        Assert.That(none.Steps, Is.EqualTo(0));
    }
}
