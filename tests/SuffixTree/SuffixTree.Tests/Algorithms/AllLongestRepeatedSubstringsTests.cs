using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;

namespace SuffixTree.Tests.Algorithms
{
    /// <summary>
    /// FindAllLongestRepeatedSubstrings (B05 follow-up): every distinct longest repeated substring with all
    /// start positions ascending, ordered by first occurrence. Expected values are from an independent
    /// O(n²) brute force (Python and the C# oracle below); definition: Gusfield 1997 §7.1.
    /// </summary>
    [TestFixture]
    [Category("Algorithms")]
    public class AllLongestRepeatedSubstringsTests
    {
        private static IEnumerable<TestCaseData> Literature()
        {
            yield return new TestCaseData("banana", new[] { ("ana", new[] { 1, 3 }) });
            yield return new TestCaseData("mississippi", new[] { ("issi", new[] { 1, 4 }) });
            yield return new TestCaseData("abracadabra", new[] { ("abra", new[] { 0, 7 }) });
            yield return new TestCaseData("ATATAT", new[] { ("ATAT", new[] { 0, 2 }) });
            yield return new TestCaseData("aaaa", new[] { ("aaa", new[] { 0, 1 }) });
            // Ties: all distinct maximal repeats are reported.
            yield return new TestCaseData("aabb", new[] { ("a", new[] { 0, 1 }), ("b", new[] { 2, 3 }) });
            yield return new TestCaseData("abcabcxyzxyz", new[] { ("abc", new[] { 0, 3 }), ("xyz", new[] { 6, 9 }) });
            yield return new TestCaseData("GATTACAGATTTACA", new[] { ("TTACA", new[] { 2, 10 }) });
            yield return new TestCaseData("xyxzyzxz", new[] { ("xz", new[] { 2, 6 }) });
            yield return new TestCaseData("abcdbcdacdab", new[] { ("bcd", new[] { 1, 4 }), ("cda", new[] { 5, 8 }) });
            yield return new TestCaseData("abcxbcaxcab", new[] { ("ab", new[] { 0, 9 }), ("bc", new[] { 1, 4 }), ("ca", new[] { 5, 8 }) });
            yield return new TestCaseData("xyzzyxxz", new[] { ("x", new[] { 0, 5, 6 }), ("y", new[] { 1, 4 }), ("z", new[] { 2, 3, 7 }) });
        }

        [TestCaseSource(nameof(Literature))]
        public void KnownTexts(string text, (string Substring, int[] Positions)[] expected)
        {
            var actual = SuffixTree.Build(text).FindAllLongestRepeatedSubstrings();
            AssertEqual(expected.Select(e => (e.Substring, (IReadOnlyList<int>)e.Positions)).ToList(), actual, text);
        }

        [TestCase("")]
        [TestCase("a")]
        [TestCase("abcd")]
        public void NoRepeat_ReturnsEmpty(string text)
        {
            Assert.That(SuffixTree.Build(text).FindAllLongestRepeatedSubstrings(), Is.Empty);
        }

        [Test]
        public void LongestRepeatedSubstring_Unchanged_AndIsOneOfTheTies()
        {
            foreach (var text in new[] { "banana", "aabb", "abcabcxyzxyz", "GATTACAGATTTACA", "xyxzyzxz" })
            {
                var tree = SuffixTree.Build(text);
                var all = tree.FindAllLongestRepeatedSubstrings();
                string lrs = tree.LongestRepeatedSubstring();
                Assert.That(all.Select(a => a.Substring), Does.Contain(lrs), text);
                Assert.That(all.All(a => a.Substring.Length == lrs.Length), text);
            }
        }

        [Test]
        public void InterfaceDefaultAndSharedAlgorithm_AgreeWithClassMethod()
        {
            var tree = SuffixTree.Build("abcabcxyzxyz");
            ISuffixTree asInterface = tree;
            AssertEqual(tree.FindAllLongestRepeatedSubstrings(), asInterface.FindAllLongestRepeatedSubstrings(), "interface");
            AssertEqual(tree.FindAllLongestRepeatedSubstrings(), SuffixTreeAlgorithms.FindAllLongestRepeatedSubstrings(tree), "shared");
            Assert.Throws<ArgumentNullException>(() => SuffixTreeAlgorithms.FindAllLongestRepeatedSubstrings(null!));
        }

        [Test]
        public void RandomTexts_MatchBruteForce([Values(0, 1, 2, 3)] int seed)
        {
            var rng = new Random(20_260_930 + seed);
            string[] alphabets = { "ab", "ACGT", "abc", "é中😀".Substring(0, 3) };
            for (int round = 0; round < 100; round++)
            {
                string alphabet = alphabets[round % alphabets.Length];
                int length = rng.Next(0, 120);
                string text = new(Enumerable.Range(0, length).Select(_ => alphabet[rng.Next(alphabet.Length)]).ToArray());
                AssertEqual(BruteForce(text), SuffixTree.Build(text).FindAllLongestRepeatedSubstrings(), $"seed {seed} round {round} '{text}'");
            }
        }

        [Test]
        public void Properties_EveryResultRepeatsAndNothingLongerRepeats()
        {
            var rng = new Random(424_242);
            for (int round = 0; round < 60; round++)
            {
                string text = new(Enumerable.Range(0, rng.Next(2, 300)).Select(_ => "ACG"[rng.Next(3)]).ToArray());
                var tree = SuffixTree.Build(text);
                var all = tree.FindAllLongestRepeatedSubstrings();
                if (all.Count == 0) continue;
                int length = all[0].Substring.Length;
                Assert.That(all.Select(a => a.Substring).Distinct().Count(), Is.EqualTo(all.Count), "distinct");
                foreach (var (substring, positions) in all)
                {
                    Assert.That(positions.Count, Is.GreaterThanOrEqualTo(2));
                    Assert.That(positions, Is.Ordered.Ascending);
                    Assert.That(positions, Is.EqualTo(tree.FindAllOccurrences(substring).OrderBy(p => p)));
                }
                Assert.That(all.Select(a => a.Positions[0]), Is.Ordered.Ascending, "first-occurrence order");
                // No repeated substring of length L + 1 exists.
                var seen = new HashSet<string>();
                for (int i = 0; i + length + 1 <= text.Length; i++)
                    Assert.That(seen.Add(text.Substring(i, length + 1)), Is.True, "longer repeat");
            }
        }

        internal static List<(string Substring, IReadOnlyList<int> Positions)> BruteForce(string text)
        {
            for (int length = text.Length - 1; length >= 1; length--)
            {
                var map = new Dictionary<string, List<int>>(StringComparer.Ordinal);
                for (int i = 0; i + length <= text.Length; i++)
                {
                    string s = text.Substring(i, length);
                    if (!map.TryGetValue(s, out var list)) map[s] = list = new List<int>();
                    list.Add(i);
                }
                var repeats = map.Where(kv => kv.Value.Count >= 2)
                    .Select(kv => (kv.Key, (IReadOnlyList<int>)kv.Value)).OrderBy(t => t.Item2[0]).ToList();
                if (repeats.Count > 0) return repeats;
            }
            return new List<(string, IReadOnlyList<int>)>();
        }

        internal static void AssertEqual(
            IReadOnlyList<(string Substring, IReadOnlyList<int> Positions)> expected,
            IReadOnlyList<(string Substring, IReadOnlyList<int> Positions)> actual, string label)
        {
            Assert.That(actual.Select(a => a.Substring), Is.EqualTo(expected.Select(e => e.Substring)), label + " substrings");
            for (int i = 0; i < expected.Count; i++)
                Assert.That(actual[i].Positions, Is.EqualTo(expected[i].Positions), $"{label} positions of '{expected[i].Substring}'");
        }
    }
}
