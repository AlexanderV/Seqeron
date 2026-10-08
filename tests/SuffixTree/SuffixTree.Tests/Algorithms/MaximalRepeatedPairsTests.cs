using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;

namespace SuffixTree.Tests.Algorithms
{
    /// <summary>
    /// FindMaximalRepeatedPairs — maximal repeated pairs (Gusfield 1997 §7.12), forward strand.
    ///
    /// Expected values marked "MUMmer 3.23" were produced by the MUMmer 3.23 binary
    /// (Ubuntu package mummer 3.23+dfsg-8) with <c>repeat-match -f -n L g.fa</c>, converted from
    /// MUMmer's 1-based (Start1, Start2, Length) to 0-based and ordered by (FirstPosition, SecondPosition)
    /// (repeat-match prints in suffix-tree traversal order). repeat-match lower-cases its input and
    /// matches every character, including N, with itself (the library default).
    /// </summary>
    [TestFixture]
    [Category("Algorithms")]
    public class MaximalRepeatedPairsTests
    {
        private static (int, int, int)[] Pairs(string text, int minLength, Func<char, bool>? unique = null)
            => SuffixTree.Build(text).FindMaximalRepeatedPairs(minLength, unique)
                .Select(p => (p.FirstPosition, p.SecondPosition, p.Length)).ToArray();

        #region MUMmer 3.23 literal outputs

        // repeat-match -f -n 1|2|3: 2 10 4 / 6 10 3 / 2 6 3 (Gusfield's xabcyabcwabcyz)
        [TestCase(1)]
        [TestCase(2)]
        [TestCase(3)]
        public void GusfieldExample_MatchesRepeatMatch(int minLength)
        {
            Assert.That(Pairs("xabcyabcwabcyz", minLength),
                Is.EqualTo(new[] { (1, 5, 3), (1, 9, 4), (5, 9, 3) }));
        }

        // repeat-match -f -n 4: 2 10 4
        [Test]
        public void GusfieldExample_MinLength4_MatchesRepeatMatch()
            => Assert.That(Pairs("xabcyabcwabcyz", 4), Is.EqualTo(new[] { (1, 9, 4) }));

        // repeat-match -f -n 3: 1 5 8 / 1 9 4 (overlapping tandem copies)
        [Test]
        public void TandemRepeat_OverlappingCopies_MatchesRepeatMatch()
            => Assert.That(Pairs("acgtacgtacgt", 3), Is.EqualTo(new[] { (0, 4, 8), (0, 8, 4) }));

        // repeat-match -f -n 1: 1 2 4 / 1 3 3 / 1 4 2 / 1 5 1
        [Test]
        public void Homopolymer_MatchesRepeatMatch()
            => Assert.That(Pairs("aaaaa", 1), Is.EqualTo(new[] { (0, 1, 4), (0, 2, 3), (0, 3, 2), (0, 4, 1) }));

        // repeat-match -f -n 2: 6 14 2 / 1 8 6;  -n 3: 1 8 6
        [Test]
        public void Gattaca_MatchesRepeatMatch()
        {
            Assert.That(Pairs("gattacagattacca", 2), Is.EqualTo(new[] { (0, 7, 6), (5, 13, 2) }));
            Assert.That(Pairs("gattacagattacca", 3), Is.EqualTo(new[] { (0, 7, 6) }));
        }

        // repeat-match -f -n 1 (19 lines): 6 14 2 / 1 8 6 / 12 15 1 / 5 15 1 / 9 15 1 / 2 15 1 / 7 12 1 /
        // 5 7 1 / 7 9 1 / 2 7 1 / 9 12 1 / 2 12 1 / 5 9 1 / 2 5 1 / 10 11 1 / 3 11 1 / 4 10 1 / 3 4 1 / 13 14 1
        [Test]
        public void Gattaca_MinLength1_MatchesRepeatMatch()
        {
            var mummer = new[]
            {
                (6, 14, 2), (1, 8, 6), (12, 15, 1), (5, 15, 1), (9, 15, 1), (2, 15, 1), (7, 12, 1),
                (5, 7, 1), (7, 9, 1), (2, 7, 1), (9, 12, 1), (2, 12, 1), (5, 9, 1), (2, 5, 1),
                (10, 11, 1), (3, 11, 1), (4, 10, 1), (3, 4, 1), (13, 14, 1),
            };
            var expected = mummer.Select(m => (m.Item1 - 1, m.Item2 - 1, m.Item3)).OrderBy(m => m.Item1).ThenBy(m => m.Item2);
            Assert.That(Pairs("gattacagattacca", 1), Is.EqualTo(expected));
        }

        // repeat-match -f -n 3: 1 6 9 / 1 11 4 — N matches N by default (repeat-match has no N handling)
        [Test]
        public void NotUnique_NMatchesItself_MatchesRepeatMatch()
            => Assert.That(Pairs("ACGTNACGTNACGT", 3), Is.EqualTo(new[] { (0, 5, 9), (0, 10, 4) }));

        #endregion

        #region Unique symbols

        // Brute-force definition with N as a mismatch: the 9-long N-spanning pair splits into three ACGT pairs.
        [Test]
        public void UniqueSymbol_NeverMatches()
        {
            Assert.That(Pairs("ACGTNACGTNACGT", 3, c => c == 'N'),
                Is.EqualTo(new[] { (0, 5, 4), (0, 10, 4), (5, 10, 4) }));
            Assert.That(Pairs("ACGNNACGNNACG", 2, c => c == 'N'),
                Is.EqualTo(new[] { (0, 5, 3), (0, 10, 3), (5, 10, 3) }));
            Assert.That(Pairs("AC$AC$AC", 2, c => c == '$'),
                Is.EqualTo(new[] { (0, 3, 2), (3, 6, 2), (0, 6, 2) }.OrderBy(p => p.Item1).ThenBy(p => p.Item2)));
        }

        [Test]
        public void UniqueSymbol_PrecedingUniqueCharacter_IsLeftMaximal()
        {
            // "NAC" / "NAC": the N's never match, so (1, 4, 2) is left-maximal although both are preceded by N.
            Assert.That(Pairs("NACNAC", 2, c => c == 'N'), Is.EqualTo(new[] { (1, 4, 2) }));
            Assert.That(Pairs("NACNAC", 2), Is.EqualTo(new[] { (0, 3, 3) }));
        }

        [Test]
        public void AllUnique_NoPairs()
            => Assert.That(Pairs("NNNNNN", 1, c => c == 'N'), Is.Empty);

        #endregion

        #region Guards and edge cases

        [Test]
        public void Guards()
        {
            var tree = SuffixTree.Build("ACGT");
            Assert.Throws<ArgumentOutOfRangeException>(() => tree.FindMaximalRepeatedPairs(0));
            Assert.Throws<ArgumentNullException>(() => SuffixTreeAlgorithms.FindMaximalRepeatedPairs(null!, 1));
            Assert.That(SuffixTree.Build("").FindMaximalRepeatedPairs(1), Is.Empty);
            Assert.That(SuffixTree.Build("A").FindMaximalRepeatedPairs(1), Is.Empty);
            Assert.That(SuffixTree.Build("AA").FindMaximalRepeatedPairs(2), Is.Empty);
            Assert.That(SuffixTree.Build("ACGT").FindMaximalRepeatedPairs(1), Is.Empty);
        }

        [Test]
        public void InterfaceDefault_EqualsSharedAlgorithm()
        {
            var tree = SuffixTree.Build("abracadabra");
            ISuffixTree view = tree;
            Assert.That(view.FindMaximalRepeatedPairs(1), Is.EqualTo(SuffixTreeAlgorithms.FindMaximalRepeatedPairs(tree, 1)));
            Assert.That(Pairs("abracadabra", 4), Is.EqualTo(new[] { (0, 7, 4) }));
        }

        [Test]
        public void NonAsciiAndSurrogates_TreatedAsPlainCodeUnits()
        {
            string text = "é😀xé😀y";
            Assert.That(Pairs(text, 2), Is.EqualTo(BruteForce(text, 2, null)));
        }

        #endregion

        #region Seeded property test vs brute force

        [Test]
        public void RandomTexts_EqualBruteForceDefinition([Values(0, 1, 2)] int seed)
        {
            var rng = new Random(20261001 + seed);
            Func<char, bool>?[] predicates = { null, c => c == 'N', c => "ACGT".IndexOf(c) < 0, c => c is 'N' or '$' };
            for (int iter = 0; iter < 80; iter++)
            {
                string alphabet = (iter % 4) switch { 0 => "AB", 1 => "ACG", 2 => "ACGTN", _ => "ACGTNN$" };
                var chars = new char[rng.Next(0, 70)];
                for (int i = 0; i < chars.Length; i++) chars[i] = alphabet[rng.Next(alphabet.Length)];
                string text = new(chars);
                int minLength = rng.Next(1, 6);
                var unique = predicates[rng.Next(predicates.Length)];
                Assert.That(Pairs(text, minLength, unique), Is.EqualTo(BruteForce(text, minLength, unique)),
                    $"seed {seed} iter {iter}: text={text} L={minLength}");
            }
        }

        /// <summary>O(n³) definition: every i &lt; j, left-maximal, extended to the first mismatch / unique symbol.</summary>
        internal static (int, int, int)[] BruteForce(string t, int minLength, Func<char, bool>? unique)
        {
            bool U(char c) => unique != null && unique(c);
            var result = new List<(int, int, int)>();
            for (int i = 0; i < t.Length; i++)
            {
                for (int j = i + 1; j < t.Length; j++)
                {
                    if (i > 0 && !U(t[i - 1]) && t[i - 1] == t[j - 1])
                        continue;
                    int l = 0;
                    while (j + l < t.Length && t[i + l] == t[j + l] && !U(t[i + l]))
                        l++;
                    if (l >= minLength)
                        result.Add((i, j, l));
                }
            }
            return result.ToArray();
        }

        #endregion
    }
}
