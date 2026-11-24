using SuffixTree.Esa.Tests.External;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.Large;

/// <summary>
/// What a closed form is compared with (work package S2.1). Three independent oracles, none of which uses the closed form it checks:
/// <list type="number">
/// <item>the definitions (<see cref="NaiveSuffixArray"/>, <see cref="NaiveLcp"/>, <see cref="BruteLcpTree"/>, <see cref="BruteForce"/>): <see cref="AgainstBruteForce"/>, n up to a few thousand;</item>
/// <item>the in-memory index of S1 (<see cref="EsaIndex"/>, the reference child table): <see cref="AgainstTheS1Index"/>, n up to 10^6;</item>
/// <item>the periodicity of the text itself (<see cref="PeriodicOracle"/>: two suffixes of the same phase are prefixes of each other, two of different phases differ within one period): the identities at n = 2^31 and above.</item>
/// </list>
/// </summary>
internal static class ClosedFormChecks
{
    /// <summary>The closed form against the definitions: every array, the statistics, the child table and the patterns, for a text of at most a few thousand symbols.</summary>
    internal static void AgainstBruteForce(ClosedFormText f)
    {
        string text = f.Materialize();
        int n = text.Length;
        string who = f.ToString();
        int[] sa = NaiveSuffixArray.Build(text);
        int[] lcp = NaiveLcp.Build(text, sa);
        var isa = new int[n + 1];
        for (int r = 0; r <= n; r++) isa[sa[r]] = r;

        for (int i = 0; i < n; i++) Same(who, $"CodeAt({i})", text[i] - 'a', f.CodeAt(i));
        Same(who, "Sigma", text.Distinct().Count(), f.Sigma);
        Same(who, "RankOf(n)", 0, f.RankOf(n));
        for (int r = 0; r <= n; r++)
        {
            Same(who, $"SuffixAt({r})", sa[r], f.SuffixAt(r));
            Same(who, $"RankOf({r})", isa[r], f.RankOf(r));
            Same(who, $"Lcp({r})", lcp[r], f.Lcp(r));
        }
        SameWide(who, "LcpSum", lcp.Sum(x => (long)x), f.LcpSum);

        // the child table, from the intervals of the definition
        List<BruteInterval> intervals = BruteLcpTree.Intervals(lcp);
        var down = new long[n + 1];
        var up = new long[n + 1];
        var next = new long[n + 1];
        foreach (BruteInterval iv in intervals)
        {
            long first = iv.LIndices[0];
            if (iv.LastChild) down[iv.Lb] = first;
            else up[iv.Rb] = first;
            for (int t = 0; t + 1 < iv.LIndices.Count; t++) next[iv.LIndices[t]] = iv.LIndices[t + 1];
        }
        for (int i = 0; i <= n; i++)
        {
            Same(who, $"Down({i})", down[i], f.Down(i));
            Same(who, $"Up({i})", up[i], f.Up(i));
            Same(who, $"Next({i})", next[i], f.Next(i));
            Same(who, $"FirstLIndex({i}, last child)", down[i], f.FirstLIndex(i, 0, true));
            Same(who, $"FirstLIndex({i}, not last)", up[i], f.FirstLIndex(0, i, false));
            Same(who, $"NextLIndex({i})", next[i] == 0 ? -1 : next[i], f.NextLIndex(i, 0));
        }

        // statistics
        Same(who, "InternalCount", intervals.Count(iv => iv.Ell > 0), f.InternalCount);
        Same(who, "NodeCount", n == 0 ? 1 : 1 + intervals.Count(iv => iv.Ell > 0) + n + 1, f.NodeCount);
        Same(who, "MaxDepth", n, f.MaxDepth);
        int bestEll = 0, bestLb = -1, bestRb = int.MaxValue;
        foreach (BruteInterval iv in intervals)
        {
            if (iv.Ell > bestEll || (iv.Ell == bestEll && iv.Ell > 0 && iv.Rb < bestRb))
            {
                bestEll = iv.Ell;
                bestLb = iv.Lb;
                bestRb = iv.Rb;
            }
        }
        Same(who, "LrsDepth", bestEll, f.LrsDepth);
        Same(who, "LrsRank", bestEll == 0 ? -1 : bestLb, f.LrsRank);
        Same(who, "LrsPosition", bestEll == 0 ? -1 : sa[bestLb], f.LrsPosition);

        // patterns: every distinct short substring, some long ones, the absent ones
        var patterns = new HashSet<string>(StringComparer.Ordinal) { string.Empty, "c", "ac", "ca", "b", "a", "ab", "ba", "bb", "aa" };
        for (int i = 0; i < n; i++)
            for (int l = 1; l <= Math.Min(8, n - i); l++) patterns.Add(text.Substring(i, l));
        if (n > 0)
        {
            patterns.Add(text);
            patterns.Add(text + "a");
            patterns.Add(text + "b");
            patterns.Add("a" + text);
            patterns.Add("b" + text);
            patterns.Add(text[1..]);
            patterns.Add(text[..^1]);
            patterns.Add(text[(n / 2)..]);
            patterns.Add(text[..(n / 2)]);
        }
        foreach (string w in patterns)
        {
            RankRange expected = ExpectedRange(text, isa, w);
            RankRange actual = f.RangeOf(w);
            if (!expected.SameAs(actual)) Assert.Fail($"{who}: RangeOf('{Shorten(w)}') expected {expected}, got {actual}");
            Same(who, $"CountOf('{Shorten(w)}')", w.Length == 0 ? n : expected.Count, f.CountOf(w));
        }
        for (int len = 0; len <= n + 1; len++)
        {
            RankRange expected = ExpectedRange(text, isa, new string('a', len));
            RankRange actual = f.RangeOfRun(len);
            if (!expected.SameAs(actual)) Assert.Fail($"{who}: RangeOfRun({len}) expected {expected}, got {actual}");
        }
        if (n <= 40)
        {
            for (int s = 0; s <= n; s++)
            {
                for (int l = 0; s + l <= n; l++)
                {
                    RankRange expected = ExpectedRange(text, isa, text.Substring(s, l));
                    RankRange actual = f.RangeOfFactor(s, l);
                    if (!expected.SameAs(actual)) Assert.Fail($"{who}: RangeOfFactor({s}, {l}) expected {expected}, got {actual}");
                }
            }
        }
        else
        {
            var rng = new SplitMix64(0x5EEDUL + (ulong)n);
            for (int c = 0; c < 400; c++)
            {
                int s = rng.NextInt(n + 1);
                int l = rng.NextInt(Math.Min(n - s, c % 4 == 0 ? n : 40) + 1);
                RankRange expected = ExpectedRange(text, isa, text.Substring(s, l));
                RankRange actual = f.RangeOfFactor(s, l);
                if (!expected.SameAs(actual)) Assert.Fail($"{who}: RangeOfFactor({s}, {l}) expected {expected}, got {actual}");
            }
        }
    }

    /// <summary>The closed form against the in-memory index of S1 (reference child table), for a text of up to 10^6 symbols.</summary>
    internal static void AgainstTheS1Index(ClosedFormText f, IEnumerable<string> patterns)
    {
        string text = f.Materialize();
        int n = text.Length;
        string who = f.ToString();
        EsaIndex index = EsaIndex.Build(text, ChildLayoutKind.Reference3);
        int[] sa = index.SuffixArray;
        CompactLcp lcp = index.Lcp;
        for (int r = 0; r <= n; r++)
        {
            if (f.SuffixAt(r) != sa[r]) Assert.Fail($"{who}: SuffixAt({r}) is {f.SuffixAt(r)}, the S1 index has {sa[r]}");
            if (f.Lcp(r) != lcp[r]) Assert.Fail($"{who}: Lcp({r}) is {f.Lcp(r)}, the S1 index has {lcp[r]}");
        }
        var isa = new int[n + 1];
        for (int r = 0; r <= n; r++) isa[sa[r]] = r;
        for (int p = 0; p <= n; p++)
            if (f.RankOf(p) != isa[p]) Assert.Fail($"{who}: RankOf({p}) is {f.RankOf(p)}, the S1 index has {isa[p]}");
        ChildTable child = index.Child;
        for (int i = 0; i <= n; i++)
        {
            if (f.Down(i) != child.FirstLIndex(i, 0, true)) Assert.Fail($"{who}: Down({i}) is {f.Down(i)}, the S1 child table has {child.FirstLIndex(i, 0, true)}");
            if (f.Up(i) != child.FirstLIndex(0, i, false)) Assert.Fail($"{who}: Up({i}) is {f.Up(i)}, the S1 child table has {child.FirstLIndex(0, i, false)}");
            if (f.NextLIndex(i, 0) != child.NextLIndex(i, 0)) Assert.Fail($"{who}: NextLIndex({i}) is {f.NextLIndex(i, 0)}, the S1 child table has {child.NextLIndex(i, 0)}");
        }
        EsaStatistics s = index.Statistics;
        Same(who, "InternalCount", s.InternalCount, f.InternalCount);
        Same(who, "NodeCount", s.NodeCount, f.NodeCount);
        Same(who, "MaxDepth", s.MaxDepth, f.MaxDepth);
        Same(who, "LrsDepth", s.LrsDepth, f.LrsDepth);
        Same(who, "LrsRank", s.LrsDepth == 0 ? -1 : s.LrsRank, f.LrsRank);
        SameWide(who, "LcpSum", lcp.ToArray().Sum(x => (long)x), f.LcpSum);
        foreach (string w in patterns)
        {
            RankRange range = f.RangeOf(w);
            long count = index.CountOccurrences(w);
            Same(who, $"CountOf('{Shorten(w)}')", count, f.CountOf(w));
            if (range.IsEmpty || w.Length == 0) continue;
            // the occurrences in suffix-array order are SA[first..last]
            IReadOnlyList<int> occurrences = index.FindAllOccurrences(w);
            if (occurrences.Count != range.Count) Assert.Fail($"{who}: '{Shorten(w)}' has {occurrences.Count} occurrences in the S1 index, the closed form says {range.Count}");
            if (occurrences[0] != f.SuffixAt(range.First) || occurrences[^1] != f.SuffixAt(range.Last))
                Assert.Fail($"{who}: '{Shorten(w)}' first / last occurrence {occurrences[0]} / {occurrences[^1]} in the S1 index, the closed form says {f.SuffixAt(range.First)} / {f.SuffixAt(range.Last)}");
        }
    }

    /// <summary>The cell of CHILD codec 1 (packed) of a closed form against the raw cells of the in-memory index.</summary>
    internal static void PackedCellsAgainstTheS1Index(ClosedFormText f)
    {
        LcpChildCase c = LcpChildCase.Of(f.Materialize());
        uint[] packed = c.Packed;
        for (int r = 0; r <= f.Length; r++)
            if (f.PackedCell(r) != packed[r]) Assert.Fail($"{f}: PackedCell({r}) is {f.PackedCell(r)}, the packed child table of S1 has {packed[r]}");
    }

    /// <summary>The ranks of the occurrences of <paramref name="w"/> in the text of the definitions.</summary>
    private static RankRange ExpectedRange(string text, int[] isa, string w)
    {
        if (w.Length == 0) return text.Length == 0 ? RankRange.Empty : new RankRange(1, text.Length);
        List<int> positions = BruteForce.Occurrences(text, w);
        if (positions.Count == 0) return RankRange.Empty;
        int lo = positions.Min(p => isa[p]);
        int hi = positions.Max(p => isa[p]);
        Assert.That(hi - lo + 1, Is.EqualTo(positions.Count), "the suffixes that start with a pattern are consecutive in the suffix array (the definitions disagree with themselves)");
        return new RankRange(lo, hi);
    }

    private static void Same(string who, string what, long expected, long actual)
    {
        if (expected != actual) Assert.Fail($"{who}: {what} expected {expected}, got {actual}");
    }

    private static void SameWide(string who, string what, Int128 expected, Int128 actual)
    {
        if (expected != actual) Assert.Fail($"{who}: {what} expected {expected}, got {actual}");
    }

    private static string Shorten(string w) => w.Length <= 24 ? w : w[..12] + "..." + w[^8..] + $" ({w.Length})";
}

/// <summary>
/// The oracle of the huge texts: facts that follow from the definition of a periodic text alone, with no table and no closed form of a family. Two suffixes whose starts differ by a multiple of the period are prefixes
/// of each other (the shorter is a prefix of the longer); two suffixes of different phases differ within one period (the period word is primitive).
/// </summary>
internal static class PeriodicOracle
{
    /// <summary>The longest common prefix of the suffixes at <paramref name="i"/> and <paramref name="j"/> of the text of <paramref name="f"/>.</summary>
    internal static long Lcp(ClosedFormText f, long i, long j)
    {
        long n = f.Length;
        if (i == j) return n - i;
        long hi = Math.Max(i, j);
        if (hi >= n) return 0;
        if ((i - j) % f.Period == 0) return n - hi;
        for (long x = 0; x <= f.Period && hi + x < n; x++)
            if (f.InfiniteCodeAt(i + x) != f.InfiniteCodeAt(j + x)) return x;
        return Math.Min(n - hi, f.Period);
    }

    /// <summary>The order of two suffixes (-1, 0, 1): a proper prefix is smaller (the end of the text is smaller than every symbol).</summary>
    internal static int Compare(ClosedFormText f, long i, long j)
    {
        if (i == j) return 0;
        long l = Lcp(f, i, j);
        long n = f.Length;
        if (i + l == n) return -1;
        if (j + l == n) return 1;
        return f.InfiniteCodeAt(i + l).CompareTo(f.InfiniteCodeAt(j + l));
    }

    /// <summary>The number of occurrences of the pattern in the prefix of length n of the periodic word, counted phase by phase: the positions <c>phase + c * period</c> that leave room for the pattern.</summary>
    internal static long CountOccurrences(ClosedFormText f, ReadOnlySpan<char> w)
    {
        long n = f.Length, p = f.Period, len = w.Length;
        if (len == 0) return n;
        long total = 0;
        for (long phase = 0; phase < p; phase++)
        {
            bool matches = true;
            for (int x = 0; x < len && matches; x++) matches = ClosedFormText.LetterOf(f.InfiniteCodeAt(phase + x)) == w[x];
            if (matches && n - len - phase >= 0) total += ((n - len - phase) / p) + 1;
        }
        return total;
    }

    /// <summary>
    /// The number of distinct non-empty substrings of the prefix of length n of the periodic word of <c>a^k b</c> (n at least twice the period): for each length l the distinct factors of that length are
    /// <c>min(l + 1, period, n - l + 1)</c>. With the sum of the LCP array it must be <c>n (n + 1) / 2</c>: <c>D = n (n + 1) / 2 - sum LCP</c>.
    /// </summary>
    internal static Int128 DistinctSubstringsOfAkB(long n, long k)
    {
        long p = k + 1;
        if (n < 2 * p) throw new ArgumentOutOfRangeException(nameof(n), n, "the formula needs n >= 2 * period");
        Int128 head = ((Int128)k * (k + 1) / 2) + k;                       // lengths 1..k: l + 1 factors
        Int128 body = (Int128)p * (n - p - k);                             // lengths k + 1 .. n - p: the period
        Int128 tail = (Int128)p * (p + 1) / 2;                             // lengths n - p + 1 .. n: n - l + 1 windows
        return head + body + tail;
    }
}
