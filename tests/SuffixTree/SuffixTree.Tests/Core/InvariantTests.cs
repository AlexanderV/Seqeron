using System;
using System.Linq;
using NUnit.Framework;

namespace SuffixTree.Tests.Core
{
    /// <summary>
    /// Tests for fundamental suffix tree invariants from theory.
    /// References: Ukkonen (1995), Gusfield (1997).
    /// </summary>
    [TestFixture]
    public class InvariantTests
    {
        #region Suffix Containment Invariant

        /// <summary>
        /// Invariant: Every suffix of the text can be found in the tree.
        /// This is the fundamental property of suffix trees.
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
        /// Invariant: CountOccurrences must equal the number of positions from FindAllOccurrences.
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
        /// Invariant: Each position from FindAllOccurrences must be a valid occurrence.
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

        #region Random Invariant Verification

        /// <summary>
        /// Verify all invariants on random input.
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

            // Invariant 2: All suffixes exist
            for (int i = 0; i < text.Length; i++)
            {
                Assert.That(tree.Contains(text.Substring(i)), Is.True);
            }

            // Invariant 3: Suffixes are correctly enumerated
            var suffixes = tree.GetAllSuffixes();
            var expected = Enumerable.Range(0, text.Length)
                .Select(i => text.Substring(i))
                .OrderBy(s => s, StringComparer.Ordinal)
                .ToList();
            Assert.That(suffixes, Is.EqualTo(expected));

            // Invariant 4: Count matches positions
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
