using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;

namespace SuffixTree.Tests.Core
{
    /// <summary>
    /// Tests for fundamental suffix tree invariants from theory.
    /// Every test here is a mathematical property that MUST hold for any correct
    /// suffix tree. If a test fails, the implementation is wrong — not the test.
    /// References: Ukkonen (1995), Gusfield (1997), Weiner (1973).
    /// </summary>
    [TestFixture]
    public class InvariantTests
    {
        #region Leaf Count Invariant

        /// <summary>
        /// Theorem (Gusfield §5.2): A suffix tree for a string of length n
        /// (with unique terminator) has exactly n leaves.
        /// </summary>
        [Test]
        [TestCase("a")]
        [TestCase("ab")]
        [TestCase("banana")]
        [TestCase("mississippi")]
        [TestCase("abracadabra")]
        [TestCase("aaaaaaa")]
        [TestCase("abababab")]
        public void LeafCount_EqualsTextLength(string text)
        {
            var tree = SuffixTree.Build(text);
            Assert.That(tree.LeafCount, Is.EqualTo(text.Length),
                $"Suffix tree for \"{text}\" must have exactly {text.Length} leaves");
        }

        #endregion

        #region Node Count Bounds

        /// <summary>
        /// Theorem (Gusfield §5.2): A suffix tree for a string of length n has
        /// at most 2n nodes (including root). NodeCount >= LeafCount + 1 (root).
        /// </summary>
        [Test]
        [TestCase("a")]
        [TestCase("ab")]
        [TestCase("banana")]
        [TestCase("mississippi")]
        [TestCase("abracadabra")]
        [TestCase("aaaaaaa")]
        public void NodeCount_WithinTheoreticalBounds(string text)
        {
            var tree = SuffixTree.Build(text);

            Assert.Multiple(() =>
            {
                Assert.That(tree.NodeCount, Is.GreaterThanOrEqualTo(tree.LeafCount + 1),
                    "NodeCount must include at least all leaves + root");
                Assert.That(tree.NodeCount, Is.LessThanOrEqualTo(2 * text.Length + 1),
                    "NodeCount cannot exceed 2n+1 (Gusfield bound)");
            });
        }

        #endregion

        #region Suffix Containment Invariant

        /// <summary>
        /// Fundamental invariant: Every suffix of the text is a path from
        /// root to a leaf. Therefore Contains(suffix) must return true.
        /// </summary>
        [Test]
        [TestCase("banana")]
        [TestCase("mississippi")]
        [TestCase("abracadabra")]
        [TestCase("xabxac")]
        [TestCase("abcabxabcd")]
        public void AllSuffixes_ExistAsSubstrings(string text)
        {
            var tree = SuffixTree.Build(text);

            for (int i = 0; i < text.Length; i++)
            {
                string suffix = text.Substring(i);
                Assert.That(tree.Contains(suffix), Is.True,
                    $"Suffix '{suffix}' starting at position {i} not found");
            }
        }

        #endregion

        #region Count/Positions Consistency

        /// <summary>
        /// Invariant: CountOccurrences (= leaf count under matched node) must equal
        /// |FindAllOccurrences| (= number of leaves collected via DFS).
        /// </summary>
        [Test]
        [TestCase("banana", "a")]
        [TestCase("banana", "an")]
        [TestCase("banana", "ana")]
        [TestCase("mississippi", "issi")]
        [TestCase("abracadabra", "abra")]
        [TestCase("aaaaaa", "aa")]
        public void CountOccurrences_EqualsPositionsCount(string text, string pattern)
        {
            var tree = SuffixTree.Build(text);

            int count = tree.CountOccurrences(pattern);
            var positions = tree.FindAllOccurrences(pattern);

            Assert.That(count, Is.EqualTo(positions.Count),
                $"CountOccurrences={count} but FindAllOccurrences returned {positions.Count} positions");
        }

        /// <summary>
        /// Invariant: Each position from FindAllOccurrences must be a valid occurrence —
        /// the substring at that position must match the pattern exactly.
        /// </summary>
        [Test]
        [TestCase("banana", "ana")]
        [TestCase("mississippi", "issi")]
        [TestCase("aaaaaa", "aa")]
        [TestCase("abcabcabc", "abc")]
        public void AllPositions_AreValidOccurrences(string text, string pattern)
        {
            var tree = SuffixTree.Build(text);
            var positions = tree.FindAllOccurrences(pattern);

            foreach (var pos in positions)
            {
                Assert.That(pos, Is.GreaterThanOrEqualTo(0));
                Assert.That(pos + pattern.Length, Is.LessThanOrEqualTo(text.Length));
                Assert.That(text.Substring(pos, pattern.Length), Is.EqualTo(pattern),
                    $"Position {pos} does not contain pattern '{pattern}'");
            }
        }

        #endregion

        #region Statistics Immutability

        /// <summary>
        /// NodeCount and LeafCount are computed once at construction.
        /// They must not change regardless of subsequent query operations.
        /// </summary>
        [Test]
        public void Stats_DoNotChangeAfterQueries()
        {
            var st = SuffixTree.Build("banana");

            var nodesBefore = st.NodeCount;
            var leafBefore = st.LeafCount;

            // Perform various operations that touch different code paths
            _ = st.Contains("ana");
            _ = st.FindAllOccurrences("a").ToList();
            _ = st.LongestRepeatedSubstring();
            _ = st.LongestCommonSubstring("other");
            _ = st.GetAllSuffixes();
            _ = st.CountOccurrences("an");

            Assert.Multiple(() =>
            {
                Assert.That(st.NodeCount, Is.EqualTo(nodesBefore));
                Assert.That(st.LeafCount, Is.EqualTo(leafBefore));
            });
        }

        #endregion

        #region MaxDepth Invariant

        /// <summary>
        /// MaxDepth equals the length of the longest suffix, which is text.Length.
        /// </summary>
        [Test]
        [TestCase("a")]
        [TestCase("banana")]
        [TestCase("mississippi")]
        [TestCase("aaaaaaa")]
        public void MaxDepth_EqualsTextLength(string text)
        {
            var tree = SuffixTree.Build(text);
            Assert.That(tree.MaxDepth, Is.EqualTo(text.Length));
        }

        #endregion

        #region Random Invariant Verification

        /// <summary>
        /// Verify all invariants on random input — property-based testing.
        /// </summary>
        [Test]
        [Repeat(10)]
        public void RandomInput_AllInvariantsHold()
        {
            var random = new Random();
            int length = random.Next(10, 100);
            string text = GenerateRandomString(random, length, "abcd");

            var tree = SuffixTree.Build(text);

            // Invariant 1: Correct leaf count
            Assert.That(tree.LeafCount, Is.EqualTo(text.Length));

            // Invariant 2: Node count bounds
            Assert.That(tree.NodeCount, Is.LessThanOrEqualTo(2 * text.Length + 1));

            // Invariant 3: All suffixes exist
            for (int i = 0; i < text.Length; i++)
            {
                Assert.That(tree.Contains(text.Substring(i)), Is.True);
            }

            // Invariant 4: Suffixes are correctly enumerated (lexicographic order)
            var suffixes = tree.GetAllSuffixes();
            var expected = Enumerable.Range(0, text.Length)
                .Select(i => text.Substring(i))
                .OrderBy(s => s, StringComparer.Ordinal)
                .ToList();
            Assert.That(suffixes, Is.EqualTo(expected));

            // Invariant 5: Count matches positions for random pattern
            string pattern = text.Substring(random.Next(text.Length / 2),
                Math.Min(3, text.Length / 2));
            Assert.That(tree.CountOccurrences(pattern),
                Is.EqualTo(tree.FindAllOccurrences(pattern).Count));
        }

        private static string GenerateRandomString(Random random, int length, string alphabet)
        {
            return new string(Enumerable.Range(0, length)
                .Select(_ => alphabet[random.Next(alphabet.Length)])
                .ToArray());
        }

        #endregion
    }
}
