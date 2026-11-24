using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.Large;

/// <summary>
/// Work package S2.1, deliverable A: the closed forms of <see cref="PowerOfA"/> (<c>a^n</c>), <see cref="PeriodicAkB"/> (<c>(a^k b)^m</c> truncated to any n, <c>(ab)^m</c> the case k = 1) are checked against three oracles that
/// do not use them: the definitions (naive suffix array, naive LCP, the lcp-interval tree by definition, naive counting; n up to 2,000), the in-memory index of S1 (n up to 10^6) and the periodicity of the text itself
/// (n = 2^31 + 10, 2^32 + 10, 2^33, 2^40). The default set: fast, no disk.
/// </summary>
[TestFixture]
public class ClosedFormSelfTests
{
    private const long TwoTo31 = 1L << 31;
    private const long TwoTo32 = 1L << 32;
    private const long TwoTo33 = 1L << 33;
    private const long TwoTo40 = 1L << 40;

    public static IEnumerable<TestCaseData> SmallShapes()
    {
        for (int n = 0; n <= 70; n++) yield return Case($"a^{n}", () => new PowerOfA(n));
        foreach (int n in new[] { 100, 255, 256, 257, 500, 1000, 2000 }) yield return Case($"a^{n}", () => new PowerOfA(n));
        for (int k = 1; k <= 8; k++)
            for (int n = 0; n <= 70; n++)
            {
                int kk = k, nn = n;
                yield return Case($"(a^{k} b)* n={n}", () => new PeriodicAkB(nn, kk));
            }
        foreach (int k in new[] { 1, 2, 3, 5, 9, 17 })
            foreach (int n in new[] { 127, 128, 129, 255, 257, 511, 1000, 1999, 2000 })
            {
                int kk = k, nn = n;
                yield return Case($"(a^{k} b)* n={n}", () => new PeriodicAkB(nn, kk));
            }
        for (int m = 0; m <= 20; m++)
        {
            int mm = m;
            yield return Case($"(ab)^{m}", () => PeriodicAkB.AbRepeat(mm));
        }
    }

    private static TestCaseData Case(string name, Func<ClosedFormText> make) => new TestCaseData(make).SetArgDisplayNames(name);

    [TestCaseSource(nameof(SmallShapes))]
    public void EqualsTheDefinitions_UpTo2000(Func<ClosedFormText> make) => ClosedFormChecks.AgainstBruteForce(make());

    public static IEnumerable<TestCaseData> LargerShapes()
    {
        foreach (long n in new long[] { 100_000, 1_000_000 }) yield return Case($"a^{n}", () => new PowerOfA(n));
        foreach (long n in new long[] { 99_999, 100_000, 100_001, 1_000_000 })
            foreach (long k in new long[] { 1, 2, 3, 7, 100, 1000 })
            {
                long nn = n, kk = k;
                yield return Case($"(a^{k} b)* n={n}", () => new PeriodicAkB(nn, kk));
            }
        yield return Case("(ab)^500000", () => PeriodicAkB.AbRepeat(500_000));
    }

    [TestCaseSource(nameof(LargerShapes))]
    public void EqualsTheS1Index_UpTo1e6(Func<ClosedFormText> make)
    {
        ClosedFormText f = make();
        var patterns = new List<string> { string.Empty, "a", "b", "aa", "ab", "ba", "bb", "c", "aaa", "aab", "baa", "abab", "aaaaaaaaaa" };
        string text = f.Materialize();
        foreach (long len in new long[] { 1, 2, 3, 5, 8, 13, 100, 999, f.Length - 1, f.Length, f.Length + 1 })
            if (len >= 0 && len <= f.Length + 1) patterns.Add(new string('a', (int)len));
        var rng = new SplitMix64((ulong)f.Length);
        for (int c = 0; c < 40 && f.Length > 0; c++)
        {
            int s = rng.NextInt((int)f.Length);
            int l = 1 + rng.NextInt((int)Math.Min(f.Length - s, c % 5 == 0 ? 3000 : 30));
            patterns.Add(text.Substring(s, l));
        }
        if (f.Length > 0) patterns.Add(text[1..]);
        ClosedFormChecks.AgainstTheS1Index(f, patterns);
    }

    public static IEnumerable<TestCaseData> PackedShapes()
    {
        foreach (long n in new long[] { 1, 2, 3, 10, 100, 4_321, 20_000 }) yield return Case($"a^{n}", () => new PowerOfA(n));
        foreach (long k in new long[] { 1, 2, 3, 6, 50 })
            foreach (long n in new long[] { 7, 64, 1_001, 20_000 })
            {
                long nn = n, kk = k;
                yield return Case($"(a^{k} b)* n={n}", () => new PeriodicAkB(nn, kk));
            }
    }

    [TestCaseSource(nameof(PackedShapes))]
    public void PackedChildCells_EqualTheCellsOfTheS1Table(Func<ClosedFormText> make) => ClosedFormChecks.PackedCellsAgainstTheS1Index(make());

    [Test]
    public void TheTextIsTheOneOfTheStringFamilies()
    {
        foreach (int n in new[] { 0, 1, 2, 5, 40, 1000 })
        {
            Assert.That(new PowerOfA(n).Materialize(), Is.EqualTo(StringFamilies.AllSame(n)), $"a^{n}");
            Assert.That(PeriodicAkB.AbRepeat(n).Materialize(), Is.EqualTo(StringFamilies.AbRepeat(2 * n)), $"(ab)^{n}");
            foreach (int k in new[] { 1, 2, 3, 10 })
                Assert.That(new PeriodicAkB(n, k).Materialize(), Is.EqualTo(StringFamilies.AkB(n, k)), $"(a^{k} b)* n={n}");
        }
    }

    [Test]
    public void ArgumentsOutsideTheTextAreRefused()
    {
        foreach (ClosedFormText f in new ClosedFormText[] { new PowerOfA(10), new PeriodicAkB(10, 2) })
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => f.SuffixAt(-1), f.Name);
            Assert.Throws<ArgumentOutOfRangeException>(() => f.SuffixAt(11), f.Name);
            Assert.Throws<ArgumentOutOfRangeException>(() => f.RankOf(11), f.Name);
            Assert.Throws<ArgumentOutOfRangeException>(() => f.Lcp(11), f.Name);
            Assert.Throws<ArgumentOutOfRangeException>(() => f.CodeAt(10), f.Name);
            Assert.Throws<ArgumentOutOfRangeException>(() => f.CodeAt(-1), f.Name);
        }
        Assert.Throws<ArgumentOutOfRangeException>(() => _ = new PowerOfA(-1));
        Assert.Throws<ArgumentOutOfRangeException>(() => _ = new PeriodicAkB(5, 0));
    }

    // ---------------------------------------------------------------------------------------- huge texts

    public static IEnumerable<TestCaseData> HugeShapes()
    {
        foreach (long n in new[] { TwoTo31 + 10, TwoTo32 + 10, TwoTo33, TwoTo40 })
        {
            long nn = n;
            yield return Case($"a^{n}", () => new PowerOfA(nn));
        }
        // (ab)^m and (a^k b)^m at the sizes of the plan; n not a multiple of the period too
        foreach (long n in new[] { TwoTo31 + 10, TwoTo32 + 10, TwoTo33, TwoTo40 })
            foreach (long k in new long[] { 1, 2, 3, 7, 12 })
            {
                long nn = n, kk = k;
                yield return Case($"(a^{k} b)* n={n}", () => new PeriodicAkB(nn, kk));
            }
        yield return Case("(ab)^(2^32 + 5)", () => PeriodicAkB.AbRepeat(TwoTo32 + 5));
        yield return Case("(a^5 b)^(2^31) exactly", () => new PeriodicAkB(6 * TwoTo31, 5));
    }

    [TestCaseSource(nameof(HugeShapes))]
    public void TheIdentitiesHold_AtTwoToThe31AndAbove(Func<ClosedFormText> make)
    {
        ClosedFormText f = make();
        long n = f.Length;
        string who = f.ToString();
        var ranks = new SortedSet<long>();
        foreach (long center in new long[] { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, TwoTo31, TwoTo32, TwoTo31 + 10, n / 2, n - 1 })
            for (long d = -4; d <= 4; d++)
                if (center + d >= 0 && center + d <= n) ranks.Add(center + d);
        for (long q = 1; q < f.Period + 3 && q < 40; q++)
            foreach (long basis in new[] { n / (q + 1), (n / 3) + q, n - (q * 7) }) if (basis >= 0 && basis <= n) ranks.Add(basis);
        var rng = new SplitMix64((ulong)n ^ 0xC10EDUL);
        for (int c = 0; c < 400; c++) ranks.Add((long)(rng.Next() % (ulong)(n + 1)));
        ranks.Add(n);

        long previous = -1;
        foreach (long r in ranks)
        {
            long position = f.SuffixAt(r);
            if (position < 0 || position > n) Assert.Fail($"{who}: SuffixAt({r}) = {position} is outside [0, {n}]");
            if (f.RankOf(position) != r) Assert.Fail($"{who}: RankOf(SuffixAt({r}) = {position}) = {f.RankOf(position)}, not {r}");
            // order and LCP of the neighbours, from the periodicity of the text
            if (r >= 2)
            {
                long before = f.SuffixAt(r - 1);
                if (PeriodicOracle.Compare(f, before, position) >= 0) Assert.Fail($"{who}: the suffixes at the ranks {r - 1} and {r} ({before}, {position}) are not in ascending order");
                long lcp = PeriodicOracle.Lcp(f, before, position);
                if (f.Lcp(r) != lcp) Assert.Fail($"{who}: Lcp({r}) = {f.Lcp(r)}, the periodicity says {lcp}");
            }
            if (r <= 1 && f.Lcp(r) != 0) Assert.Fail($"{who}: Lcp({r}) must be 0");
            // the child table: the relations between a stored field and the LCP values at the boundaries
            long down = f.Down(r), up = f.Up(r), next = f.Next(r);
            if (down != 0 && (down <= r || down > n || (r > 0 && f.Lcp(down) <= f.Lcp(r)))) Assert.Fail($"{who}: Down({r}) = {down} is not an l-index of an interval that starts at {r}");
            if (up != 0 && (up > r || r == n || f.Lcp(up) <= f.Lcp(r + 1))) Assert.Fail($"{who}: Up({r}) = {up} is not the first l-index of a non-last child that ends at {r}");
            if (next != 0 && (next <= r || next > n || f.Lcp(next) != f.Lcp(r))) Assert.Fail($"{who}: Next({r}) = {next} is not the next l-index of {r}");
            if (previous >= 0 && r == previous) Assert.Fail("the sample must be distinct");
            previous = r;
        }
        Assert.That(f.RankOf(n), Is.EqualTo(0), $"{who}: RankOf(n)");
        Assert.That(f.SuffixAt(0), Is.EqualTo(n), $"{who}: SuffixAt(0)");
        Assert.That(f.Up(n), Is.EqualTo(0), $"{who}: nothing ends a non-last child at n");

        // positions: the inverse in the other direction
        foreach (long p in new long[] { 0, 1, 2, 3, TwoTo31 - 1, TwoTo31, TwoTo31 + 1, TwoTo32 - 1, TwoTo32, TwoTo32 + 1, n / 2, n - 3, n - 2, n - 1, n })
        {
            if (p < 0 || p > n) continue;
            long rank = f.RankOf(p);
            if (rank < 0 || rank > n || f.SuffixAt(rank) != p) Assert.Fail($"{who}: SuffixAt(RankOf({p}) = {rank}) = {(rank >= 0 && rank <= n ? f.SuffixAt(rank) : -1)}");
        }

        // the sum of the LCP array against the number of distinct substrings
        Int128 all = (Int128)n * (n + 1) / 2;
        Int128 distinct = f is PowerOfA ? n : PeriodicOracle.DistinctSubstringsOfAkB(n, f.Period - 1);
        if (all - f.LcpSum != distinct) Assert.Fail($"{who}: n (n + 1) / 2 - sum LCP is {all - f.LcpSum}, the number of distinct substrings is {distinct}");

        // the statistics, each from first principles. The longest repeat of a periodic text is the text without its first period (a^n: without one symbol) and occurs one period earlier;
        // the internal nodes are the non-empty right-special factors: for a^n the a^j with j < n; for (a^k b)* the a^j with j < k and the repeated suffixes of the text (n - period of them), the a^j that are suffixes counted once
        long t = n % f.Period;
        long expectedLrs = f is PowerOfA ? n - 1 : n - f.Period;
        long expectedInternal = f is PowerOfA ? n - 1 : n - f.Period + (f.Period - 2) - Math.Min(t, f.Period - 2);
        Assert.That(f.LrsDepth, Is.EqualTo(expectedLrs), $"{who}: LrsDepth");
        long other = f.LrsPosition - f.Period >= 0 ? f.LrsPosition - f.Period : f.LrsPosition + f.Period;
        Assert.That(PeriodicOracle.Lcp(f, f.LrsPosition, other), Is.GreaterThanOrEqualTo(f.LrsDepth), $"{who}: the longest repeat occurs again one period away");
        Assert.That(f.InternalCount, Is.EqualTo(expectedInternal), $"{who}: InternalCount (the non-empty right-special factors)");
        Assert.That(f.NodeCount, Is.EqualTo(1 + f.InternalCount + n + 1), $"{who}: NodeCount");
        Assert.That(f.MaxDepth, Is.EqualTo(n), $"{who}: MaxDepth");

        // patterns: the run of a's, factors of the text, the ends of the range
        var runs = new List<long> { 1, 2, 3, 4, 5, 10, 100, 12345, TwoTo31 - 1, TwoTo31, TwoTo31 + 1, TwoTo32 - 1, TwoTo32 + 1, n - 1, n, n + 1 };
        foreach (long len in runs)
        {
            RankRange range = f.RangeOfRun(len);
            long expected = RunCount(f, len);
            if (range.Count != expected) Assert.Fail($"{who}: RangeOfRun({len}) has {range.Count} ranks, the phase count says {expected}");
            if (!range.IsEmpty) CheckRangeEnds(f, range, 0, len); // a run that occurs at all occurs at position 0
        }
        for (int c = 0; c < 60; c++)
        {
            long len = c % 3 == 0 ? 1 + (long)(rng.Next() % (ulong)Math.Min(n, 2 * f.Period + 3)) : 1 + (long)(rng.Next() % (ulong)n);
            long start = (long)(rng.Next() % (ulong)(n - len + 1));
            RankRange range = f.RangeOfFactor(start, len);
            if (range.IsEmpty) Assert.Fail($"{who}: the factor ({start}, {len}) of the text has no occurrence");
            long expected = FactorCount(f, start, len);
            if (range.Count != expected) Assert.Fail($"{who}: RangeOfFactor({start}, {len}) has {range.Count} ranks, the phase count says {expected}");
            CheckRangeEnds(f, range, start, len);
        }
    }

    [Test]
    public void ExplicitPatterns_HaveTheCountOfThePhaseSum_AtTwoToThe33()
    {
        foreach (ClosedFormText f in new ClosedFormText[] { new PowerOfA(TwoTo33), new PeriodicAkB(TwoTo33, 1), new PeriodicAkB(TwoTo33 + 3, 3), new PeriodicAkB(TwoTo32 + 10, 2) })
        {
            var patterns = new List<string> { "a", "b", "aa", "ab", "ba", "bb", "aaa", "aab", "aba", "baa", "bab", "abab", "baba", "aaaa", "aaab", "abaa", "baab", "c", "ac", "cb", string.Empty };
            string period = new string(Enumerable.Range(0, (int)f.Period + 6).Select(i => ClosedFormText.LetterOf(f.InfiniteCodeAt(i))).ToArray());
            for (int s = 0; s < 6; s++) patterns.Add(period.Substring(s, Math.Min(period.Length - s, (int)f.Period + 2)));
            foreach (string w in patterns)
            {
                long expected = PeriodicOracle.CountOccurrences(f, w);
                RankRange range = f.RangeOf(w);
                if (range.Count != expected) Assert.Fail($"{f}: RangeOf('{w}') has {range.Count} ranks, the phase count says {expected}");
                if (!range.IsEmpty && w.Length > 0) CheckRangeEnds(f, range, w);
            }
        }
    }

    // ---------------------------------------------------------------------------------------- helpers

    /// <summary>The ranks of a pattern are exactly <c>[first..last]</c>: the suffixes at both ends begin with it and the neighbours do not.</summary>
    private static void CheckRangeEnds(ClosedFormText f, RankRange range, string w)
    {
        bool Begins(long position)
        {
            if (position + w.Length > f.Length) return false;
            for (int x = 0; x < w.Length; x++)
                if (ClosedFormText.LetterOf(f.InfiniteCodeAt(position + x)) != w[x]) return false;
            return true;
        }
        if (!Begins(f.SuffixAt(range.First)) || !Begins(f.SuffixAt(range.Last))) Assert.Fail($"{f}: the ends of {range} do not begin with '{w}'");
        if (range.First > 1 && Begins(f.SuffixAt(range.First - 1))) Assert.Fail($"{f}: the rank before {range} begins with '{w}'");
        if (range.Last < f.Length && Begins(f.SuffixAt(range.Last + 1))) Assert.Fail($"{f}: the rank after {range} begins with '{w}'");
    }

    /// <summary>The same for a pattern that is a factor of the text (<c>text[start..start+len)</c>), by the longest common prefix of the periodic text.</summary>
    private static void CheckRangeEnds(ClosedFormText f, RankRange range, long start, long len)
    {
        bool Begins(long position) => position + len <= f.Length && PeriodicOracle.Lcp(f, position, start) >= len;
        if (!Begins(f.SuffixAt(range.First)) || !Begins(f.SuffixAt(range.Last))) Assert.Fail($"{f}: the ends of {range} do not begin with the factor ({start}, {len})");
        if (range.First > 1 && Begins(f.SuffixAt(range.First - 1))) Assert.Fail($"{f}: the rank before {range} begins with the factor ({start}, {len})");
        if (range.Last < f.Length && Begins(f.SuffixAt(range.Last + 1))) Assert.Fail($"{f}: the rank after {range} begins with the factor ({start}, {len})");
    }

    /// <summary>The occurrences of <c>a^len</c>, phase by phase: it fits at a phase when the first min(len, period + 1) symbols there are a's (a full period of a's means the word is a^n).</summary>
    private static long RunCount(ClosedFormText f, long len)
    {
        long n = f.Length, p = f.Period, total = 0;
        for (long phase = 0; phase < p; phase++)
        {
            bool all = true;
            for (long x = 0; x < Math.Min(len, p + 1) && all; x++) all = f.InfiniteCodeAt(phase + x) == 0;
            if (all && n - len - phase >= 0) total += ((n - len - phase) / p) + 1;
        }
        return total;
    }

    /// <summary>The occurrences of the factor <c>text[start..start+len)</c>, phase by phase.</summary>
    private static long FactorCount(ClosedFormText f, long start, long len)
    {
        long n = f.Length, p = f.Period, total = 0, startPhase = start % p;
        for (long phase = 0; phase < p; phase++)
        {
            bool match = true;
            if (phase != startPhase)
                for (long x = 0; x < Math.Min(len, p + 1) && match; x++) match = f.InfiniteCodeAt(phase + x) == f.InfiniteCodeAt(startPhase + x);
            if (match && n - len - phase >= 0) total += ((n - len - phase) / p) + 1;
        }
        return total;
    }
}
