using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;

namespace SuffixTree.Tests.Algorithms
{
    /// <summary>
    /// FindAllDistinctLongestCommonSubstrings — every tied longest common substring with all positions
    /// (matching statistics, Gusfield 1997 §7.8). Literal expectations come from a Python brute force
    /// (all substrings of the query, ordered by first query occurrence).
    /// </summary>
    [TestFixture]
    [Category("Algorithms")]
    public class AllDistinctLongestCommonSubstringsTests
    {
        private static (string, int[], int[])[] All(ISuffixTreeAnalysis tree, string other)
            => tree.FindAllDistinctLongestCommonSubstrings(other)
                .Select(r => (r.Substring, r.PositionsInText.ToArray(), r.PositionsInOther.ToArray())).ToArray();

        [Test]
        public void TwoTiedSubstrings_BothReported()
        {
            var all = All(SuffixTree.Build("abcxyzabc"), "xyzqabcpxyz");
            Assert.That(all, Is.EqualTo(new[]
            {
                ("xyz", new[] { 3 }, new[] { 0, 8 }),
                ("abc", new[] { 0, 6 }, new[] { 4 }),
            }));
        }

        [Test]
        public void FourTies_GattacaTagacca()
        {
            var all = All(SuffixTree.Build("GATTACA"), "TAGACCA");
            Assert.That(all, Is.EqualTo(new[]
            {
                ("TA", new[] { 3 }, new[] { 0 }),
                ("GA", new[] { 0 }, new[] { 2 }),
                ("AC", new[] { 4 }, new[] { 3 }),
                ("CA", new[] { 5 }, new[] { 5 }),
            }));
        }

        [Test]
        public void FirstEntry_EqualsCanonicalFindAllLongestCommonSubstrings()
        {
            var tree = SuffixTree.Build("abcxyzabc");
            var canonical = tree.FindAllLongestCommonSubstrings("xyzqabcpxyz");
            var first = tree.FindAllDistinctLongestCommonSubstrings("xyzqabcpxyz")[0];
            Assert.That(first.Substring, Is.EqualTo(canonical.Substring));
            Assert.That(first.PositionsInText, Is.EqualTo(canonical.PositionsInText));
            Assert.That(first.PositionsInOther, Is.EqualTo(canonical.PositionsInOther));
        }

        [Test]
        public void NoCommonCharacter_Empty()
        {
            Assert.That(SuffixTree.Build("AAAA").FindAllDistinctLongestCommonSubstrings("CCC"), Is.Empty);
            Assert.That(SuffixTree.Build("").FindAllDistinctLongestCommonSubstrings("CCC"), Is.Empty);
            Assert.That(SuffixTree.Build("ACGT").FindAllDistinctLongestCommonSubstrings(""), Is.Empty);
            Assert.Throws<ArgumentNullException>(() => SuffixTree.Build("ACGT").FindAllDistinctLongestCommonSubstrings(null!));
        }

        [Test]
        public void RandomPairs_EqualBruteForce_AndDefaultImplementation([Values(0, 1)] int seed)
        {
            var rng = new Random(20261002 + seed);
            for (int iter = 0; iter < 150; iter++)
            {
                string alphabet = iter % 3 == 0 ? "AB" : "ACGT";
                string text = Random(rng, rng.Next(0, 50), alphabet);
                string other = Random(rng, rng.Next(0, 50), alphabet);
                var tree = SuffixTree.Build(text);
                var expected = BruteForce(text, other);
                string ctx = $"seed {seed} iter {iter}: text={text} other={other}";
                Assert.That(All(tree, other), Is.EqualTo(expected), ctx);
                Assert.That(All(new ExternalAnalysis(tree), other), Is.EqualTo(expected), ctx + " default");
            }
        }

        private static string Random(Random rng, int length, string alphabet)
        {
            var c = new char[length];
            for (int i = 0; i < length; i++) c[i] = alphabet[rng.Next(alphabet.Length)];
            return new string(c);
        }

        internal static (string, int[], int[])[] BruteForce(string text, string other)
        {
            int best = 0;
            for (int i = 0; i < other.Length; i++)
                for (int j = i + best + 1; j <= other.Length && text.Contains(other[i..j], StringComparison.Ordinal); j++)
                    best = j - i;
            if (best == 0) return Array.Empty<(string, int[], int[])>();

            var order = new List<string>();
            var inOther = new Dictionary<string, List<int>>();
            for (int q = 0; q + best <= other.Length; q++)
            {
                string s = other.Substring(q, best);
                if (!text.Contains(s, StringComparison.Ordinal)) continue;
                if (!inOther.TryGetValue(s, out var list)) { inOther[s] = list = new List<int>(); order.Add(s); }
                list.Add(q);
            }
            return order.Select(s => (s,
                Enumerable.Range(0, text.Length - best + 1).Where(r => string.CompareOrdinal(text, r, s, 0, best) == 0).ToArray(),
                inOther[s].ToArray())).ToArray();
        }

        /// <summary>An implementer outside the library: only the default interface body is available.</summary>
        private sealed class ExternalAnalysis(SuffixTree inner) : ISuffixTreeSearch, ISuffixTreeAnalysis
        {
            public ITextSource Text => inner.Text;
            public int NodeCount => inner.NodeCount;
            public int LeafCount => inner.LeafCount;
            public int MaxDepth => inner.MaxDepth;
            public bool IsEmpty => inner.IsEmpty;
            public bool Contains(string value) => inner.Contains(value);
            public bool Contains(ReadOnlySpan<char> value) => inner.Contains(value);
            public IReadOnlyList<int> FindAllOccurrences(string pattern) => inner.FindAllOccurrences(pattern);
            public IReadOnlyList<int> FindAllOccurrences(ReadOnlySpan<char> pattern) => inner.FindAllOccurrences(pattern);
            public int CountOccurrences(string pattern) => inner.CountOccurrences(pattern);
            public int CountOccurrences(ReadOnlySpan<char> pattern) => inner.CountOccurrences(pattern);
            public string LongestRepeatedSubstring() => inner.LongestRepeatedSubstring();
            public IReadOnlyList<string> GetAllSuffixes() => inner.GetAllSuffixes();
            public IEnumerable<string> EnumerateSuffixes() => inner.EnumerateSuffixes();
            public string LongestCommonSubstring(string other) => inner.LongestCommonSubstring(other);
            public string LongestCommonSubstring(ReadOnlySpan<char> other) => inner.LongestCommonSubstring(other);
            public (string Substring, int PositionInText, int PositionInOther) LongestCommonSubstringInfo(string other)
                => inner.LongestCommonSubstringInfo(other);
            public (string Substring, IReadOnlyList<int> PositionsInText, IReadOnlyList<int> PositionsInOther) FindAllLongestCommonSubstrings(string other)
                => inner.FindAllLongestCommonSubstrings(other);
            public IReadOnlyList<(int PositionInText, int PositionInQuery, int Length)> FindExactMatchAnchors(string query, int minLength)
                => inner.FindExactMatchAnchors(query, minLength);
        }
    }
}
