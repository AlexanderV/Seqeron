using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;

namespace SuffixTree.Tests.Algorithms
{
    /// <summary>
    /// Longest common substring of k strings / k-common substring problem (Gusfield 1997 §7.6) over a
    /// generalized suffix tree. Literal expectations: Rosalind LCSM sample (any of AC, CA, TA is accepted
    /// there; the library returns all of them) and a Python brute force over all substrings.
    /// </summary>
    [TestFixture]
    [Category("Algorithms")]
    public class KCommonSubstringTests
    {
        private static readonly string[] RosalindLcsm = { "GATTACA", "TAGACCA", "ATACA" };

        [Test]
        public void RosalindLcsmSample_AllLongestCommonSubstrings()
            => Assert.That(SuffixTree.FindLongestCommonSubstrings(RosalindLcsm), Is.EqualTo(new[] { "AC", "CA", "TA" }));

        [Test]
        public void RosalindLcsmSample_SupportTwo()
            => Assert.That(SuffixTree.FindLongestCommonSubstrings(RosalindLcsm, 2), Is.EqualTo(new[] { "TACA" }));

        [Test]
        public void RosalindLcsmSample_LengthsBySupport()
            => Assert.That(SuffixTree.LongestCommonSubstringLengthsBySupport(RosalindLcsm), Is.EqualTo(new[] { 7, 7, 4, 2 }));

        [Test]
        public void SupportOne_IsLongestText()
        {
            Assert.That(SuffixTree.FindLongestCommonSubstrings(new[] { "AB", "CDE", "XYZ", "AB" }, 1), Is.EqualTo(new[] { "CDE", "XYZ" }));
            Assert.That(SuffixTree.FindLongestCommonSubstrings(new[] { "", "" }, 1), Is.Empty);
        }

        [Test]
        public void SingleText_IsItself()
            => Assert.That(SuffixTree.FindLongestCommonSubstrings(new[] { "banana" }), Is.EqualTo(new[] { "banana" }));

        [Test]
        public void NothingShared_OrEmptyText_Empty()
        {
            Assert.That(SuffixTree.FindLongestCommonSubstrings(new[] { "AAA", "CCC" }), Is.Empty);
            Assert.That(SuffixTree.FindLongestCommonSubstrings(new[] { "ACGT", "", "ACGT" }), Is.Empty);
            Assert.That(SuffixTree.FindLongestCommonSubstrings(new[] { "ACGT", "", "ACGT" }, 2), Is.EqualTo(new[] { "ACGT" }));
        }

        [Test]
        public void IdenticalTexts_AndRepeatsWithinOneText()
        {
            Assert.That(SuffixTree.FindLongestCommonSubstrings(new[] { "ACGT", "ACGT", "ACGT" }), Is.EqualTo(new[] { "ACGT" }));
            // "AAAA" repeated inside text 0 must not count as support from two texts.
            Assert.That(SuffixTree.FindLongestCommonSubstrings(new[] { "AAAAAAA", "CAAC", "GGG" }, 2), Is.EqualTo(new[] { "AA" }));
        }

        [Test]
        public void InputsUsingPrivateUseCharacters_GetOtherSeparators()
        {
            string[] texts = { "ab", "xab", "ab" };
            Assert.That(SuffixTree.FindLongestCommonSubstrings(texts), Is.EqualTo(new[] { "ab" }));
            var (text, owner, starts) = SuffixTreeAlgorithms.BuildGeneralizedText(texts);
            Assert.That(starts, Is.EqualTo(new[] { 0, 6, 12 }));
            Assert.That(owner.Count(o => o < 0), Is.EqualTo(3));
            var separators = Enumerable.Range(0, text.Length).Where(i => owner[i] < 0).Select(i => text[i]).ToArray();
            Assert.That(separators.Distinct().Count(), Is.EqualTo(3));
            Assert.That(separators.Any(c => texts.Any(t => t.Contains(c))), Is.False);
        }

        [Test]
        public void Guards()
        {
            Assert.Throws<ArgumentNullException>(() => SuffixTree.FindLongestCommonSubstrings(null!));
            Assert.Throws<ArgumentException>(() => SuffixTree.FindLongestCommonSubstrings(Array.Empty<string>()));
            Assert.Throws<ArgumentException>(() => SuffixTree.FindLongestCommonSubstrings(new[] { "A", null! }));
            Assert.Throws<ArgumentOutOfRangeException>(() => SuffixTree.FindLongestCommonSubstrings(new[] { "A", "B" }, 0));
            Assert.Throws<ArgumentOutOfRangeException>(() => SuffixTree.FindLongestCommonSubstrings(new[] { "A", "B" }, 3));
            Assert.Throws<ArgumentNullException>(() => SuffixTreeAlgorithms.FindLongestCommonSubstrings(new[] { "A" }, 1, null!));
        }

        [Test]
        public void RandomSets_EqualBruteForce([Values(0, 1, 2)] int seed)
        {
            var rng = new Random(20261003 + seed);
            for (int iter = 0; iter < 60; iter++)
            {
                string alphabet = (iter % 3) switch { 0 => "AC", 1 => "ACG", _ => "ACGT" };
                int k = rng.Next(1, 7);
                var texts = Enumerable.Range(0, k)
                    .Select(_ => new string(Enumerable.Range(0, rng.Next(0, 30)).Select(_ => alphabet[rng.Next(alphabet.Length)]).ToArray()))
                    .ToArray();
                var (lengths, bySupport) = BruteForce(texts);
                string ctx = $"seed {seed} iter {iter}: {string.Join(",", texts)}";
                Assert.That(SuffixTree.LongestCommonSubstringLengthsBySupport(texts), Is.EqualTo(lengths), ctx);
                for (int q = 1; q <= k; q++)
                    Assert.That(SuffixTree.FindLongestCommonSubstrings(texts, q), Is.EqualTo(bySupport[q]), $"{ctx} q={q}");
            }
        }

        internal static (int[] Lengths, string[][] BySupport) BruteForce(IReadOnlyList<string> texts)
        {
            var support = new Dictionary<string, HashSet<int>>(StringComparer.Ordinal);
            for (int t = 0; t < texts.Count; t++)
                for (int i = 0; i < texts[t].Length; i++)
                    for (int j = i + 1; j <= texts[t].Length; j++)
                    {
                        string s = texts[t][i..j];
                        if (!support.TryGetValue(s, out var set)) support[s] = set = new HashSet<int>();
                        set.Add(t);
                    }

            int k = texts.Count;
            var lengths = new int[k + 1];
            var bySupport = new string[k + 1][];
            for (int q = 1; q <= k; q++)
            {
                lengths[q] = support.Where(kv => kv.Value.Count >= q).Select(kv => kv.Key.Length).DefaultIfEmpty(0).Max();
                int l = lengths[q];
                bySupport[q] = support.Where(kv => kv.Value.Count >= q && kv.Key.Length == l && l > 0)
                    .Select(kv => kv.Key).OrderBy(s => s, StringComparer.Ordinal).ToArray();
            }
            lengths[0] = lengths[1];
            return (lengths, bySupport);
        }
    }
}
