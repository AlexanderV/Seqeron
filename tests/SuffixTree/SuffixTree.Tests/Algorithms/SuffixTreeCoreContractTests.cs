using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;

namespace SuffixTree.Tests.Algorithms
{
    /// <summary>
    /// SUFFIXTREE-CORE review (2026-09, B05): documented contracts of LCS / FindAll-LCS /
    /// exact-match anchors / distinct-substring counting, locked against independent brute force
    /// (values cross-checked with a Python brute-force enumeration).
    /// </summary>
    [TestFixture]
    [Category("Algorithms")]
    public class SuffixTreeCoreContractTests
    {
        // ---------------- LCS ----------------

        [Test]
        public void Lcs_GusfieldXabxacExample_IsAbxa()
        {
            // Gusfield 1997 running example strings; brute force: "abxa" (text 1, other 3).
            var info = SuffixTree.Build("xabxac").LongestCommonSubstringInfo("abcabxabcd");
            Assert.That(info, Is.EqualTo(("abxa", 1, 3)));
        }

        [Test]
        public void LcsInfo_Tie_ReturnsSubstringFirstFoundInOther()
        {
            // "ab" vs "bandana": both "a" and "b" have maximal length 1; "b" occurs first in other.
            var info = SuffixTree.Build("ab").LongestCommonSubstringInfo("bandana");
            Assert.That(info, Is.EqualTo(("b", 1, 0)));
        }

        [Test]
        public void FindAllLcs_PositionsAreAscendingAndDistinct()
        {
            var (s, inText, inOther) = SuffixTree.Build("abracadabra").FindAllLongestCommonSubstrings("xaxa");
            Assert.Multiple(() =>
            {
                Assert.That(s, Is.EqualTo("a"));
                Assert.That(inText, Is.EqualTo(new[] { 0, 3, 5, 7, 10 }));
                Assert.That(inOther, Is.EqualTo(new[] { 1, 3 }));
            });
        }

        [Test]
        public void FindAllLcs_Tie_ReportsOnlyCanonicalSubstring()
        {
            // "ab" vs "bandana": canonical is "b"; "a" (same length) is not reported.
            var (s, inText, inOther) = SuffixTree.Build("ab").FindAllLongestCommonSubstrings("bandana");
            Assert.Multiple(() =>
            {
                Assert.That(s, Is.EqualTo("b"));
                Assert.That(inText, Is.EqualTo(new[] { 1 }));
                Assert.That(inOther, Is.EqualTo(new[] { 0 }));
            });
        }

        [Test]
        public void Lcs_RandomInputs_TieBreakAndPositionsMatchBruteForce([Values(0, 1, 2, 3, 4, 5, 6, 7)] int seed)
        {
            var rng = new Random(seed);
            for (int iter = 0; iter < 40; iter++)
            {
                string text = Random(rng, rng.Next(1, 30), seed % 2 == 0 ? "AC" : "ACGT");
                string other = Random(rng, rng.Next(1, 30), "ACGT");
                var tree = SuffixTree.Build(text);

                // Brute force: longest substring of other contained in text, earliest end in other.
                string expected = string.Empty;
                for (int end = 0; end < other.Length; end++)
                    for (int start = 0; start <= end; start++)
                    {
                        int len = end - start + 1;
                        if (len > expected.Length && text.Contains(other.Substring(start, len), StringComparison.Ordinal))
                            expected = other.Substring(start, len);
                    }

                var info = tree.LongestCommonSubstringInfo(other);
                var all = tree.FindAllLongestCommonSubstrings(other);
                string ctx = $"text={text} other={other}";
                Assert.That(info.Substring, Is.EqualTo(expected), ctx);
                Assert.That(all.Substring, Is.EqualTo(expected), ctx);
                if (expected.Length == 0) continue;

                Assert.That(info.PositionInOther, Is.EqualTo(other.IndexOf(expected, StringComparison.Ordinal)), ctx);
                Assert.That(text.Substring(info.PositionInText, expected.Length), Is.EqualTo(expected), ctx);
                Assert.That(all.PositionsInText, Is.EqualTo(Occurrences(text, expected)), ctx);
                Assert.That(all.PositionsInOther, Is.EqualTo(Occurrences(other, expected)), ctx);
            }
        }

        // ---------------- Exact-match anchors ----------------

        [Test]
        public void Anchors_AdjacentRunsMayOverlapInQuery()
        {
            // Matching statistics of "ababa" vs "aba": 1,2,3,2,3 → two runs (ms >= 3) → two anchors
            // that share query position 2 (documented; brute-force verified).
            var anchors = SuffixTree.Build("aba").FindExactMatchAnchors("ababa", 3);
            Assert.That(anchors, Is.EqualTo(new[] { (0, 0, 3), (0, 2, 3) }));
        }

        [Test]
        public void Anchors_RandomInputs_AreRunPeaksOfMatchingStatistics_AndMems([Values(0, 1, 2, 3, 4, 5)] int seed)
        {
            var rng = new Random(100 + seed);
            for (int iter = 0; iter < 40; iter++)
            {
                string text = Random(rng, rng.Next(1, 35), "ACGT");
                string query = Random(rng, rng.Next(1, 35), "ACGT");
                int minLength = rng.Next(1, 5);
                var tree = SuffixTree.Build(text);
                var anchors = tree.FindExactMatchAnchors(query, minLength);

                // Independent matching statistics: ms[i] = longest suffix of query[0..i] occurring in text.
                var ms = new int[query.Length];
                for (int i = 0; i < query.Length; i++)
                    for (int s = 0; s <= i; s++)
                        if (text.Contains(query.Substring(s, i - s + 1), StringComparison.Ordinal))
                        {
                            ms[i] = i - s + 1;
                            break;
                        }

                var expected = new List<(int Q, int L)>();
                for (int i = 0; i < query.Length;)
                {
                    if (ms[i] < minLength) { i++; continue; }
                    int peak = i;
                    int j = i;
                    for (; j < query.Length && ms[j] >= minLength; j++)
                        if (ms[j] > ms[peak]) peak = j;
                    expected.Add((peak - ms[peak] + 1, ms[peak]));
                    i = j;
                }

                string ctx = $"text={text} query={query} min={minLength}";
                Assert.That(anchors.Select(a => (a.PositionInQuery, a.Length)).ToList(), Is.EqualTo(expected), ctx);
                foreach (var (t, q, l) in anchors)
                {
                    Assert.That(text.Substring(t, l), Is.EqualTo(query.Substring(q, l)), ctx);
                    bool leftExtends = t > 0 && q > 0 && text[t - 1] == query[q - 1];
                    bool rightExtends = t + l < text.Length && q + l < query.Length && text[t + l] == query[q + l];
                    Assert.That(leftExtends || rightExtends, Is.False, $"{ctx}: anchor ({t},{q},{l}) is not a MEM");
                }
            }
        }

        // ---------------- Distinct substrings (Troyanskaya 2002) ----------------

        [TestCase("a$b$", 9L, new long[] { 1, 3, 3, 2, 1 })]
        [TestCase("\0a\0a", 7L, new long[] { 1, 2, 2, 2, 1 })]
        [TestCase("ATATAT", 11L, new long[] { 1, 2, 2, 2, 2, 2, 1 })]
        [TestCase("héllo héllo", 50L, new long[] { 1, 5, 6, 6, 6, 6, 6, 5, 4, 3, 2, 1 })]
        [TestCase("\U0001F600a\U0001F600b", 18L, new long[] { 1, 4, 4, 4, 3, 2, 1 })] // UTF-16 code units
        public void DistinctSubstrings_TerminatorLikeAndNonAsciiText(string text, long total, long[] byLength)
        {
            var tree = SuffixTree.Build(text);
            Assert.Multiple(() =>
            {
                Assert.That(tree.CountDistinctSubstrings(), Is.EqualTo(total));
                Assert.That(tree.CountDistinctSubstringsByLength(text.Length), Is.EqualTo(byLength));
                Assert.That(SuffixTreeAlgorithms.CountDistinctSubstrings(tree), Is.EqualTo(total));
                Assert.That(SuffixTreeAlgorithms.CountDistinctSubstringsByLength(tree, text.Length), Is.EqualTo(byLength));
            });
        }

        [Test]
        public void Lrs_TerminatorLikeText_IsMaximalRepeat()
        {
            Assert.Multiple(() =>
            {
                Assert.That(SuffixTree.Build("a$b$a$").LongestRepeatedSubstring(), Is.EqualTo("a$"));
                Assert.That(SuffixTree.Build("\0a\0a").LongestRepeatedSubstring(), Is.EqualTo("\0a"));
            });
        }

        private static string Random(Random rng, int n, string alphabet)
        {
            var c = new char[n];
            for (int i = 0; i < n; i++) c[i] = alphabet[rng.Next(alphabet.Length)];
            return new string(c);
        }

        private static List<int> Occurrences(string s, string p)
        {
            var result = new List<int>();
            for (int i = 0; i + p.Length <= s.Length; i++)
                if (string.CompareOrdinal(s, i, p, 0, p.Length) == 0) result.Add(i);
            return result;
        }
    }
}
