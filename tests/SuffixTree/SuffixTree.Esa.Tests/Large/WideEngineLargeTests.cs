using SuffixTree.Esa.External;

namespace SuffixTree.Esa.Tests.Large;

/// <summary>
/// Work package S2.4, the generator-backed operation matrix of the long engine (ADR-005 section 3.7 layer 1; category <c>Large</c>, explicit): a text of 2^31 + 10, 2^32 + 10 or 2^33 symbols behind <see cref="PosGen{TWidth}"/> readers (no
/// memory), the long descent against the closed forms of S2.1 (<see cref="ClosedFormText"/>): the count of <c>a^k</c> for k = 1 .. 10,000 and for k near n, the factors of the text at and around 2^31 and 2^32, the pages of
/// <c>Locate</c> that straddle them, the statistics and the distinct substrings by one pass over the generated LCP array, and the boundaries of the stored width (n = 2^31 - 2 .. 2^32 - 3). The long runs of <c>a^k</c> (k near n) walk a chain of k nodes: one
/// step is a few reads, so they take a minute each (the readers of <see cref="PowerReaders"/> compute <c>a^n</c> in a few instructions).
/// </summary>
[TestFixture]
[Category("Large")]
[Explicit(ExplicitReason)]
public class WideEngineLargeTests : LargeFixtureBase
{
    private const long TwoTo31 = 1L << 31;
    private const long TwoTo32 = 1L << 32;
    private const long TwoTo33 = 1L << 33;

    /// <summary>The sizes of the matrix of the brief: the first text above 2^31, a text above 2^32, 2^33.</summary>
    internal static readonly long[] Sizes = [TwoTo31 + 10, TwoTo32 + 10, TwoTo33];

    /// <summary>The boundaries of the stored width and of the 32-bit views.</summary>
    internal static readonly long[] Boundaries = [TwoTo31 - 2, TwoTo31 - 1, TwoTo31, TwoTo32 - 4, TwoTo32 - 3];

    /// <summary>The count of <c>a^k</c> in <c>a^n</c> by the long descent over the cheap readers of <c>a^n</c>: the interval and the number of suffixes that start with it.</summary>
    private static (bool Found, long Lo, long Hi) RunOfA(long n, long k)
    {
        var text = new PowerReaders.Text(n);
        var cells = new WideCells<CheckedAccess, PowerReaders.Child, PowerReaders.Lcp>(new PowerReaders.Child(n), new PowerReaders.Lcp(n), Lcp(n), tagged: false, n);
        bool found = WideSearch.Run<CheckedAccess, PowerReaders.Sa, WideCells<CheckedAccess, PowerReaders.Child, PowerReaders.Lcp>, PowerReaders.Text, RunPattern>(new PowerReaders.Sa(n), cells, text, n, new RunPattern(k), out long lo, out long hi);
        return (found, lo, hi);
    }

    private static CompactLcp Lcp(long n) => GenDescent.Lcp64(new PowerOfA(n));

    // ---------------------------------------------------------------------------------------------- the cheap readers against the closed form

    [Test]
    public void ThePowerReaders_AreTheClosedFormOfPowerOfA([ValueSource(nameof(AllSizes))] long n)
    {
        var f = new PowerOfA(n);
        var sa = new PowerReaders.Sa(n);
        var child = new PowerReaders.Child(n);
        var lcp = new PowerReaders.Lcp(n);
        foreach (long r in GeneratedViewLargeTests.Ranks(n))
        {
            Assert.That((long)sa.Get(r), Is.EqualTo(f.SuffixAt(r)), $"SA[{r}]");
            Assert.That((long)child.Get(r), Is.EqualTo(f.PackedCell(r)), $"CHILD[{r}]");
            Assert.That((long)lcp.Get(r), Is.EqualTo(f.Lcp(r)), $"LCP[{r}]");
        }
    }

    public static IEnumerable<long> AllSizes() => Sizes.Concat(Boundaries);

    // ---------------------------------------------------------------------------------------------- the count of a^k

    /// <summary>k = 1 .. 10,000 at every size: the count n - k + 1 and the interval [k..n], found by the long descent over the cheap readers.</summary>
    [Test]
    public void TheCountOfARun_ForK1To10000_IsTheClosedForm([ValueSource(nameof(AllSizes))] long n)
    {
        for (long k = 1; k <= 10_000; k++)
        {
            (bool found, long lo, long hi) = RunOfA(n, k);
            Assert.That(found, Is.True, $"a^{k} in a^{n}");
            Assert.That((lo, hi, hi - lo + 1), Is.EqualTo((k, n, n - k + 1)), $"a^{k} in a^{n}");
        }
        Assert.That(RunOfA(n, 1_000_000).Found, Is.True);
    }

    /// <summary>The same descent over the general generated readers (<see cref="PosGen{TWidth}"/>) at both widths, for k up to 2,000: the cheap readers and the general ones agree with the closed form.</summary>
    [Test]
    public void TheCountOfARun_OverTheGeneralGeneratedReaders_IsTheClosedForm([ValueSource(nameof(Sizes))] long n)
    {
        var f = new PowerOfA(n);
        for (long k = 1; k <= 2_000; k += 7)
        {
            RankRange expected = f.RangeOfRun(k);
            Assert.That(GenDescent.Width5(f, new RunPattern(k), out long lo, out long hi), Is.True);
            Assert.That(new RankRange(lo, hi).SameAs(expected), Is.True, $"a^{k} in a^{n}: {lo}..{hi}");
            if (n <= TwoTo32 - 4 || n == TwoTo31 + 10)
            {
                Assert.That(GenDescent.Width4(f, new RunPattern(k), out long lo4, out long hi4), Is.True);
                Assert.That((lo4, hi4), Is.EqualTo((lo, hi)), $"a^{k} in a^{n}: the 4-byte readers");
            }
        }
    }

    /// <summary>
    /// k near n: the descent walks the chain of nodes a^1 .. a^k (one step per level), so the matched length, the depth of the node and the rank pass 2^31 and 2^32 inside the loop. One minute each at 2^31 + 10, a few at 2^33.
    /// </summary>
    [TestCase(TwoTo31 - 2, 0)]
    [TestCase(TwoTo31 + 10, 0)]
    [TestCase(TwoTo31 + 10, 1)]
    [TestCase(TwoTo31 + 10, 2)]
    [TestCase(TwoTo32 + 10, 1)]
    [TestCase(TwoTo33, 1)]
    public void TheCountOfARun_NearN_IsTheClosedForm_AfterAChainOfMoreThan2To31Steps(long n, int back)
    {
        long k = n - back;
        (bool found, long lo, long hi) = RunOfA(n, k);
        Assert.That(found, Is.True, $"a^{k} in a^{n}");
        Assert.That((lo, hi, hi - lo + 1), Is.EqualTo((k, n, n - k + 1)), $"a^{k} in a^{n}: a^(n - 1) occurs twice, a^n once");
    }

    /// <summary>A run one symbol longer than the text: the walk goes down the whole chain of n nodes and the last leaf says no (the remainder is longer than the suffix).</summary>
    [TestCase(TwoTo31 - 2)]
    [TestCase(TwoTo31 - 1)]
    [TestCase(TwoTo31 + 10)]
    public void ARunLongerThanTheText_IsAbsent_AfterTheWholeChain_NearN(long n)
    {
        (bool found, long lo, long hi) = RunOfA(n, n + 1);
        Assert.That(found, Is.False, $"a^{n + 1} in a^{n}");
        Assert.That((lo, hi), Is.EqualTo((0L, -1L)));
    }

    // ---------------------------------------------------------------------------------------------- the factors of a periodic text at the boundaries

    public static IEnumerable<ClosedFormText> PeriodicForms()
    {
        foreach (long n in new[] { TwoTo31 + 10, TwoTo32 + 10, TwoTo33 })
        {
            yield return new PeriodicAkB(n, 1);
            yield return new PeriodicAkB(n, 3);
            yield return new PeriodicAkB(n, 20);
        }
    }

    [Test]
    public void TheFactors_OfAPeriodicText_AroundTheBoundaries2To31And2To32_AreFound([ValueSource(nameof(PeriodicForms))] ClosedFormText f)
    {
        long n = f.Length;
        foreach (long start in new long[] { 0, 1, 7, TwoTo31 - 9, TwoTo31 - 1, TwoTo31, TwoTo31 + 3, TwoTo32 - 9, TwoTo32 - 1, TwoTo32, TwoTo32 + 1, n / 2, n - 700 })
        {
            foreach (long length in new long[] { 1, 2, 5, 17, 100, 600 })
            {
                if (start + length > n) continue;
                RankRange expected = f.RangeOfFactor(start, length);
                Assert.That(GenDescent.Width5(f, new FactorPattern(f, start, length), out long lo, out long hi), Is.True, $"{f} factor [{start}, +{length})");
                Assert.That(new RankRange(lo, hi).SameAs(expected), Is.True, $"{f} factor [{start}, +{length}): {lo}..{hi}, the closed form says {expected}");
            }
        }
    }

    [Test]
    public void AnAbsentPattern_OfAPeriodicText_IsAbsent_AtEverySize([ValueSource(nameof(PeriodicForms))] ClosedFormText f)
    {
        // "bb" never occurs in (a^k b)*, a letter that is not in the text is no match, and a factor with one symbol changed is absent whenever the closed form says so
        Assert.That(GenDescent.Width5(f, new DensePattern([1, 1]), out long lo, out long hi), Is.False);
        Assert.That((lo, hi), Is.EqualTo((0L, -1L)));
        Assert.That(GenDescent.Width5(f, new DensePattern([0, 0, 2]), out _, out _), Is.False);
        long start = TwoTo31 - 30;
        var chars = new int[40];
        for (int i = 0; i < chars.Length; i++) chars[i] = f.CodeAt(start + i);
        chars[20] = 1 - chars[20];
        string mutated = new([.. chars.Select(c => ClosedFormText.LetterOf(c))]);
        RankRange expected = f.RangeOf(mutated);
        bool found = GenDescent.Width5(f, new DensePattern([.. chars]), out long mlo, out long mhi);
        Assert.That(found, Is.EqualTo(!expected.IsEmpty), $"{f} '{mutated}'");
        if (found) Assert.That(new RankRange(mlo, mhi).SameAs(expected), Is.True);
    }

    private readonly struct DensePattern(int[] symbols) : IWidePattern
    {
        public long Length => symbols.Length;

        public int this[long index] => symbols[index];
    }

    // ---------------------------------------------------------------------------------------------- the boundaries of the stored width

    [Test]
    public void TheDescent_AtTheBoundariesOfTheStoredWidth_HasNoOffByOne([ValueSource(nameof(Boundaries))] long n)
    {
        foreach (ClosedFormText f in new ClosedFormText[] { new PowerOfA(n), new PeriodicAkB(n, 1), new PeriodicAkB(n, 4) })
        {
            foreach (long start in new long[] { 0, 1, n - 12, n - 8, n - 5, n - 1 })
            {
                foreach (long length in new long[] { 1, 2, 3, 5, 11 })
                {
                    if (start < 0 || start + length > n) continue;
                    RankRange expected = f.RangeOfFactor(start, length);
                    var pat = new FactorPattern(f, start, length);
                    Assert.That(GenDescent.Width5(f, pat, out long lo5, out long hi5), Is.True, $"{f} [{start}, +{length}) at width 5");
                    Assert.That(new RankRange(lo5, hi5).SameAs(expected), Is.True, $"{f} [{start}, +{length}): {lo5}..{hi5}, the closed form says {expected}");
                    Assert.That(GenDescent.Width4(f, pat, out long lo4, out long hi4), Is.True, $"{f} [{start}, +{length}) at width 4");
                    Assert.That((lo4, hi4), Is.EqualTo((lo5, hi5)), $"{f} [{start}, +{length}): the widths agree");
                }
            }
        }
    }

    // ---------------------------------------------------------------------------------------------- the pages of Locate

    [Test]
    public void ThePagesOfLocate_AroundTheBoundaries_ReadTheSuffixArrayAtRanksAbove2To31([ValueSource(nameof(AllSizes))] long n)
    {
        var f = new PowerOfA(n);
        var sa = new SaView<PosGen<W5>>(new PosGen<W5>(f, GenSection.Sa));
        long lo = 4, hi = n; // the node a^4: the ranks 4..n
        foreach (long centre in new long[] { 10, TwoTo31, TwoTo32, n - 20 })
        {
            if (centre < 12 || centre > n - 10) continue;
            foreach (int capacity in new[] { 1, 16, 33 })
            {
                long skip = centre - capacity / 2 - lo;
                int count = EsaIndex.PageOf(n, lo, hi, skip, capacity, out long start);
                int expectedCount = (int)Math.Min(capacity, n - (lo + skip) + 1); // the last page of the node is shorter
                Assert.That((count, start), Is.EqualTo((expectedCount, lo + skip)), $"page at {centre}, capacity {capacity}");
                var page = new long[count];
                sa.CopyTo(start, page);
                for (int i = 0; i < count; i++) Assert.That(page[i], Is.EqualTo(n - (start + i)), $"page at {centre}: element {i} is the suffix at the rank {start + i}");
            }
        }
        // the page that holds the last ranks, and the one beyond them
        Assert.That(EsaIndex.PageOf(n, lo, hi, n - lo - 2, 10, out long tail), Is.EqualTo(3));
        Assert.That(tail, Is.EqualTo(n - 2));
        Assert.That(EsaIndex.PageOf(n, lo, hi, n - lo + 1, 10, out _), Is.EqualTo(0));
    }

    // ---------------------------------------------------------------------------------------------- the statistics by one pass, and the distinct substrings

    private static CompactLcp GeneratedLcp(ClosedFormText f) =>
        CompactLcp.FromReaders(LcpKind.U32, new PosGen<W4>(f, GenSection.LcpBase), new PosGen<W4>(f, GenSection.LcpPairsOfPowerOfA, long.MaxValue), isChecked: true);

    /// <summary>
    /// One pass over 2^31 + 2^20 generated LCP values of <c>u^2049</c> for the de Bruijn word of order 10 over 4 letters (the periodic text whose scan stack stays short: 2,049 + 10 open intervals; <c>(ab)^m</c> and <c>a^n</c> nest n / 2 and n
    /// deep): the header statistics by the wide scan, against the closed form.
    /// </summary>
    [Test]
    public void TheScanOfTheLcpArray_OfADeBruijnRepeat_GivesTheStatisticsOfTheHeader_At2To31()
    {
        var f = new DeBruijnRepeat(4, 10, 2049);
        Assert.That(f.Length, Is.GreaterThan(TwoTo31));
        EsaStatistics st = LcpIntervalScan.StatisticsOf(GeneratedLcp(f));
        Assert.That((st.Length, st.NodeCount, st.InternalCount, st.LeafCount, st.MaxDepth, st.LrsDepth, st.LrsRank), Is.EqualTo((f.Length, f.NodeCount, f.InternalCount, f.Length, f.Length, f.LrsDepth, f.LrsRank)));
        Assert.That(st.StackPeak, Is.InRange(2049, 2049 + 12), "the open intervals of u^m: one per period, and the levels of the word");
    }

    /// <summary>
    /// One pass (about two minutes) over the generated LCP array of <c>a^n</c>, n = 2^32 + 10: the distinct substrings of a^n are n of them, and the term n (n + 1) / 2 is above 2^63 (it wraps in 64 bits: the sum is taken in 128).
    /// </summary>
    [Test]
    public void TheDistinctSubstrings_OfPowerOfA_AreCounted_At2To32Plus10_WhereNTimesNPlus1WrapsIn64Bits()
    {
        var f = new PowerOfA(TwoTo32 + 10);
        Assert.That(EsaIndex.CountDistinctSubstrings(GenDescent.Lcp64(f), f.Length), Is.EqualTo(f.Length), "a^n has n distinct substrings: one of every length");
    }

    /// <summary>One pass (about a minute) over the generated LCP array of <c>a^n</c>, n = 2^31 + 10: one distinct substring of every length, the histogram of the LCP values clamped at 6.</summary>
    [Test]
    public void TheDistinctSubstringsByLength_OfPowerOfA_AreCounted_At2To31Plus10()
    {
        var f = new PowerOfA(TwoTo31 + 10);
        Assert.That(EsaIndex.CountDistinctSubstringsByLength(GenDescent.Lcp64(f), f.Length, 6), Is.EqualTo(new long[] { 1, 1, 1, 1, 1, 1, 1 }));
    }
}
