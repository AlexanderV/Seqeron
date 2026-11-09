using SuffixTree.Esa.Tests.TestSupport;
using SuffixTree.Persistent;
using UkkonenTree = global::SuffixTree.SuffixTree;

namespace SuffixTree.Esa.Tests.Search;

/// <summary>
/// Search on the in-memory ESA (work package 1b.2): pattern descent, <c>Contains</c>, <c>CountOccurrences</c>, <c>FindAllOccurrences</c>, lazy locate. Oracles: brute force (L0),
/// the Ukkonen tree (L2) and the persistent v6 tree (L3). Occurrence lists are compared as sorted sets (the order is engine order, <c>ISuffixTreeSearch.cs:49-54</c>).
/// </summary>
[TestFixture]
public class EsaSearchTests
{
    private static readonly ChildLayoutKind[] Layouts = Enum.GetValues<ChildLayoutKind>();

    private static IEnumerable<string> DistinctSubstrings(string text, int maxLength = int.MaxValue)
    {
        var seen = new HashSet<string>(StringComparer.Ordinal);
        for (int i = 0; i < text.Length; i++)
            for (int l = 1; l <= Math.Min(maxLength, text.Length - i); l++)
                if (seen.Add(text.Substring(i, l))) yield return text.Substring(i, l);
    }

    /// <summary>Patterns that are NOT (necessarily) substrings: one substituted symbol, a symbol outside the text, extension at both ends, longer than the text.</summary>
    private static IEnumerable<string> NearMisses(string text, SplitMix64 rng, int count)
    {
        char[] alphabet = text.Distinct().ToArray();
        for (int k = 0; k < count && text.Length > 0; k++)
        {
            int len = 1 + rng.NextInt(Math.Min(12, text.Length));
            int start = rng.NextInt(text.Length - len + 1);
            char[] pattern = text.Substring(start, len).ToCharArray();
            pattern[rng.NextInt(len)] = alphabet[rng.NextInt(alphabet.Length)];
            yield return new string(pattern);
            pattern[rng.NextInt(len)] = (char)0x4E00;       // a code unit the text does not contain (unless it is CJK text)
            yield return new string(pattern);
            yield return text.Substring(start, len) + alphabet[rng.NextInt(alphabet.Length)];
            yield return alphabet[rng.NextInt(alphabet.Length)] + text.Substring(start, len);
        }
        yield return text + text;
        yield return text + "~";
        yield return "~" + text;
        yield return "￾";
        yield return "\0";
    }

    private static void AssertSameAnswer(EsaIndex esa, UkkonenTree ukkonen, PersistentSuffixTree? v6, string text, string pattern, string context)
    {
        List<int> expected = BruteForce.Occurrences(text, pattern);
        string where = $"{context} pattern '{pattern}'";
        Assert.That(esa.Contains(pattern.AsSpan()), Is.EqualTo(expected.Count > 0), $"Contains {where}");
        Assert.That(esa.CountOccurrences(pattern.AsSpan()), Is.EqualTo(expected.Count), $"Count {where}");
        Assert.That(esa.FindAllOccurrences(pattern.AsSpan()).OrderBy(x => x), Is.EqualTo(expected), $"FindAll {where}");
        Assert.That(ukkonen.CountOccurrences(pattern), Is.EqualTo(expected.Count), $"Ukkonen count {where}");
        Assert.That(ukkonen.FindAllOccurrences(pattern).OrderBy(x => x), Is.EqualTo(expected), $"Ukkonen FindAll {where}");
        if (v6 is not null)
        {
            Assert.That(v6.CountOccurrences(pattern), Is.EqualTo(expected.Count), $"v6 count {where}");
            Assert.That(v6.FindAllOccurrences(pattern).OrderBy(x => x), Is.EqualTo(expected), $"v6 FindAll {where}");
            Assert.That(v6.Contains(pattern), Is.EqualTo(expected.Count > 0), $"v6 Contains {where}");
        }
    }

    private static void AssertEverySubstring(string text, string context, bool withV6, int nearMisses = 20)
    {
        UkkonenTree ukkonen = UkkonenTree.Build(text);
        using PersistentSuffixTree? v6 = withV6 ? PersistentSuffixTreeFactory.CreatePersistent(new StringTextSource(text)) : null;
        foreach (ChildLayoutKind layout in Layouts)
        {
            EsaIndex esa = EsaIndex.Build(text, layout);
            foreach (string pattern in DistinctSubstrings(text))
                AssertSameAnswer(esa, ukkonen, layout == ChildLayoutKind.Packed ? v6 : null, text, pattern, $"{context} {layout}");
            var rng = new SplitMix64((ulong)text.Length * 7919 + 13);
            foreach (string pattern in NearMisses(text, rng, nearMisses))
                AssertSameAnswer(esa, ukkonen, layout == ChildLayoutKind.Packed ? v6 : null, text, pattern, $"{context} {layout} (absent?)");
        }
    }

    [Test]
    public void Banana_EverySubstringAndAbsentPatterns()
    {
        foreach (string s in new[] { "banana", "mississippi", "abracadabra", "a", "aa", "ab" })
            AssertEverySubstring(s, $"'{s}'", withV6: true);
    }

    [Test]
    public void EmptyPatternSemantics_AsThePersistentTree()
    {
        // PersistentSuffixTree.Search.cs:10,18 (Contains true), :35-36 (count n), :58-60 (all start positions 0..n-1)
        foreach (string s in new[] { "", "a", "banana", "aaaa" })
        {
            EsaIndex esa = EsaIndex.Build(s);
            using PersistentSuffixTree v6 = PersistentSuffixTreeFactory.CreatePersistent(new StringTextSource(s));
            Assert.That(esa.Contains(""), Is.True, $"'{s}'");
            Assert.That(esa.Contains(ReadOnlySpan<char>.Empty), Is.True, $"'{s}'");
            Assert.That(esa.CountOccurrences(""), Is.EqualTo(s.Length), $"'{s}'");
            Assert.That(esa.CountOccurrences(ReadOnlySpan<char>.Empty), Is.EqualTo(v6.CountOccurrences("")), $"'{s}'");
            Assert.That(esa.FindAllOccurrences(""), Is.EqualTo(Enumerable.Range(0, s.Length)), $"'{s}'");
            Assert.That(esa.FindAllOccurrences(ReadOnlySpan<char>.Empty), Is.EqualTo(v6.FindAllOccurrences("")), $"'{s}'");
        }
    }

    [Test]
    public void EmptyText_AnyNonEmptyPatternIsAbsent()
    {
        EsaIndex esa = EsaIndex.Build("");
        foreach (string p in new[] { "a", "ab", "\0", "~" })
        {
            Assert.That(esa.Contains(p), Is.False, p);
            Assert.That(esa.CountOccurrences(p), Is.EqualTo(0), p);
            Assert.That(esa.FindAllOccurrences(p), Is.Empty, p);
            Assert.That(esa.TryMatch(p, out _, out _), Is.False, p);
        }
    }

    [Test]
    public void NullPattern_Throws()
    {
        EsaIndex esa = EsaIndex.Build("banana");
        Assert.That(() => esa.Contains((string)null!), Throws.TypeOf<ArgumentNullException>());
        Assert.That(() => esa.CountOccurrences((string)null!), Throws.TypeOf<ArgumentNullException>());
        Assert.That(() => esa.FindAllOccurrences((string)null!), Throws.TypeOf<ArgumentNullException>());
    }

    [Test]
    public void AllSmallStrings_EverySubstring()
    {
        foreach (var (sigma, maxLen) in new[] { (2, 8), (3, 5), (4, 4) })
            foreach (string s in ExhaustiveStrings.All(sigma, maxLen).Where(s => s.Length > 0))
            {
                UkkonenTree ukkonen = UkkonenTree.Build(s);
                EsaIndex esa = EsaIndex.Build(s);
                foreach (string pattern in DistinctSubstrings(s).Concat(new[] { s + "a", "c" + s, "z" }))
                    AssertSameAnswer(esa, ukkonen, null, s, pattern, $"'{s}'");
            }
    }

    [Test]
    public void AdversarialFamilies_EverySubstring_UpTo200Characters()
    {
        var rng = new SplitMix64(0x1B20);
        for (int i = 0; i < 6; i++)
        {
            int n = 40 + rng.NextInt(161);
            ulong seed = rng.Next();
            foreach (var (name, text) in StringFamilies.Adversarial(n, seed).Where((_, k) => (k + i) % 6 == 0))
                AssertEverySubstring(text, $"{name} n {n}", withV6: i == 0);
        }
    }

    [Test]
    public void EdgeCases_Unicode_AndWideAlphabets()
    {
        foreach (string s in new[]
        {
            "\uD800a𐀀\uDBFF", "\0￿\0￿\0", StringFamilies.Bmp(150, 3), StringFamilies.SurrogatePairs(120, 5),
            StringFamilies.LoneSurrogates(150, 7), StringFamilies.NulAndFfff(150, 9), StringFamilies.AllDistinct(200, 11),
            StringFamilies.Random(180, 300, 13),    // sigma > 256: the ushort path
        })
            AssertEverySubstring(s, $"unicode n {s.Length}", withV6: false, nearMisses: 40);
    }

    [Test]
    public void CodeUnitsAroundTheAlphabetBounds_AreAbsent()
    {
        // the symbol map ends at the largest code unit of the text: the next code unit, the largest possible one and the smallest one must be absent, not out of range
        foreach (string text in new[] { "banana", "\u0001\u0002ab", "zzz\uFFFE", "a\uFFFF" })
        {
            EsaIndex esa = EsaIndex.Build(text);
            char max = text.Max();
            foreach (char c in new[] { (char)(max + 1), char.MaxValue, '\0', (char)(max - 1) })
            {
                bool present = text.Contains(c, StringComparison.Ordinal);
                Assert.That(esa.Contains(c.ToString()), Is.EqualTo(present), $"'{text}' / U+{(int)c:X4}");
                Assert.That(esa.CountOccurrences(c.ToString()), Is.EqualTo(text.Count(x => x == c)), $"'{text}' / U+{(int)c:X4}");
                Assert.That(esa.Contains(text.Substring(0, 1) + c), Is.EqualTo(text.Contains(text.Substring(0, 1) + c, StringComparison.Ordinal)), $"'{text}' / U+{(int)c:X4} second");
            }
        }
    }

    [Test]
    public void PatternsBeyondTheStackBuffer_AreHandled()
    {
        string text = StringFamilies.RandomDna(3000, 5) + new string('A', 700) + StringFamilies.RandomDna(100, 6);
        EsaIndex esa = EsaIndex.Build(text);
        UkkonenTree ukkonen = UkkonenTree.Build(text);
        foreach (int len in new[] { 511, 512, 513, 700, 1000, 3800 })
        {
            string pattern = text.Substring(10, Math.Min(len, text.Length - 10));
            AssertSameAnswer(esa, ukkonen, null, text, pattern, $"len {len}");
            AssertSameAnswer(esa, ukkonen, null, text, pattern + "A", $"len {len} + A");
        }
        AssertSameAnswer(esa, ukkonen, null, text, new string('A', 650), "A^650");
        AssertSameAnswer(esa, ukkonen, null, text, new string('A', 701), "A^701 (absent)");
    }

    [Test]
    public void TryMatch_ReturnsTheRankIntervalOfTheNonEmptySuffixes()
    {
        string text = "banana";
        EsaIndex esa = EsaIndex.Build(text);
        Assert.That(esa.TryMatch("ana", out int lo, out int hi), Is.True);
        Assert.That((hi - lo + 1), Is.EqualTo(2));
        Assert.That(lo, Is.GreaterThanOrEqualTo(1));
        for (int r = lo; r <= hi; r++) Assert.That(text.AsSpan(esa.SuffixArray[r]).StartsWith("ana"), Is.True);
        Assert.That(esa.TryMatch("a", out lo, out hi), Is.True);
        Assert.That(hi - lo + 1, Is.EqualTo(3));
        Assert.That(() => esa.TryMatch("", out _, out _), Throws.TypeOf<ArgumentException>());
    }

    [Test]
    public void Locate_PagesCoverEveryPositionOnce_ForAllWidthsOfPage()
    {
        string text = StringFamilies.Periodic(300, 3, 4);
        EsaIndex esa = EsaIndex.Build(text);
        foreach (string pattern in DistinctSubstrings(text, 4))
        {
            Assert.That(esa.TryMatch(pattern, out int lo, out int hi), Is.True);
            List<int> expected = BruteForce.Occurrences(text, pattern);
            foreach (int pageSize in new[] { 1, 2, 5, 1000 })
            {
                var collected = new List<int>();
                var page = new int[pageSize];
                long skip = 0;
                for (int written; (written = esa.Locate(lo, hi, skip, page)) > 0; skip += written)
                    collected.AddRange(page.Take(written));
                Assert.That(collected.OrderBy(x => x), Is.EqualTo(expected), $"'{pattern}' page {pageSize}");

                var wide = new long[pageSize];
                var collectedWide = new List<long>();
                skip = 0;
                for (int written; (written = esa.Locate(lo, hi, skip, wide)) > 0; skip += written)
                    collectedWide.AddRange(wide.Take(written));
                Assert.That(collectedWide.OrderBy(x => x), Is.EqualTo(expected.Select(x => (long)x)), $"'{pattern}' wide page {pageSize}");
            }
        }
    }

    [Test]
    public void Locate_RejectsBadArguments()
    {
        EsaIndex esa = EsaIndex.Build("banana");
        Assert.That(esa.TryMatch("an", out int lo, out int hi), Is.True);
        Assert.That(() => esa.Locate(lo, hi, -1, new int[2]), Throws.TypeOf<ArgumentOutOfRangeException>());
        Assert.That(() => esa.Locate(hi, lo, 0, new int[2]), Throws.TypeOf<ArgumentOutOfRangeException>());
        Assert.That(() => esa.Locate(lo, esa.Length + 1, 0, new int[2]), Throws.TypeOf<ArgumentOutOfRangeException>());
        Assert.That(esa.Locate(lo, hi, 5, new int[2]), Is.EqualTo(0));
        Assert.That(esa.Locate(lo, hi, 0, Span<int>.Empty), Is.EqualTo(0));
    }

    [Test]
    public void LargerTexts_RandomSubstringsAndAbsentPatterns()
    {
        foreach (var (name, text) in new[]
        {
            ("dna", StringFamilies.RandomDna(20_000, 21)),
            ("protein", StringFamilies.RandomProtein(8_000, 22)),
            ("fibonacci", StringFamilies.Fibonacci(6_000)),
            ("abab", StringFamilies.AbRepeat(5_000)),
            ("aaaa", StringFamilies.AllSame(3_000)),
        })
        {
            UkkonenTree ukkonen = UkkonenTree.Build(text);
            foreach (ChildLayoutKind layout in Layouts)
            {
                EsaIndex esa = EsaIndex.Build(text, layout);
                var rng = new SplitMix64(0x1B21);
                for (int q = 0; q < 300; q++)
                {
                    int len = 1 + rng.NextInt(40);
                    int start = rng.NextInt(text.Length - len + 1);
                    string pattern = text.Substring(start, len);
                    Assert.That(esa.CountOccurrences(pattern), Is.EqualTo(ukkonen.CountOccurrences(pattern)), $"{name} {layout} '{pattern}'");
                    string absent = pattern[..^1] + (char)(pattern[^1] ^ 1);
                    Assert.That(esa.CountOccurrences(absent), Is.EqualTo(ukkonen.CountOccurrences(absent)), $"{name} {layout} '{absent}'");
                }
            }
        }
    }

    [Test]
    public void ConcurrentReaders_AgreeWithTheSequentialAnswers()
    {
        string text = StringFamilies.RandomDna(5_000, 31);
        EsaIndex esa = EsaIndex.Build(text);
        var rng = new SplitMix64(33);
        string[] patterns = Enumerable.Range(0, 400).Select(_ => text.Substring(rng.NextInt(4_980), 1 + rng.NextInt(15))).ToArray();
        int[] expected = patterns.Select(p => esa.CountOccurrences(p)).ToArray();
        var errors = new System.Collections.Concurrent.ConcurrentBag<string>();
        Parallel.For(0, 4, new ParallelOptions { MaxDegreeOfParallelism = 4 }, _ =>
        {
            for (int round = 0; round < 5; round++)
                for (int k = 0; k < patterns.Length; k++)
                    if (esa.CountOccurrences(patterns[k]) != expected[k]) errors.Add(patterns[k]);
        });
        Assert.That(errors, Is.Empty);
    }
}
