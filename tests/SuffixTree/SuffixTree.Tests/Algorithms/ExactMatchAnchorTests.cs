using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;

namespace SuffixTree.Tests.Algorithms
{
    /// <summary>
    /// Tests for FindExactMatchAnchors — the O(n+m) suffix-link streaming
    /// algorithm that finds right-maximal exact matches (anchors) between
    /// the tree's text and a query.
    ///
    /// Theory: An exact-match anchor at (posInText, posInQuery, length) means
    /// text[posInText..posInText+length) == query[posInQuery..posInQuery+length).
    /// The algorithm emits the peak (longest) match within each contiguous run
    /// where matchLen >= minLength. Anchors are non-overlapping in the query.
    /// </summary>
    [TestFixture]
    [Category("Algorithms")]
    public class ExactMatchAnchorTests
    {
        private static readonly int[] RandomSeeds = Enumerable.Range(0, 5).ToArray();

        #region Null and Edge Cases

        [Test]
        public void FindAnchors_NullQuery_ThrowsArgumentNullException()
        {
            var tree = SuffixTree.Build("banana");
            Assert.Throws<ArgumentNullException>(() => tree.FindExactMatchAnchors(null!, 3));
        }

        [Test]
        public void FindAnchors_EmptyQuery_ReturnsEmpty()
        {
            var tree = SuffixTree.Build("banana");
            var result = tree.FindExactMatchAnchors("", 1);
            Assert.That(result, Is.Empty);
        }

        [Test]
        public void FindAnchors_EmptyTree_ReturnsEmpty()
        {
            var tree = SuffixTree.Build("");
            var result = tree.FindExactMatchAnchors("abc", 1);
            Assert.That(result, Is.Empty);
        }

        [Test]
        public void FindAnchors_MinLengthZero_ReturnsEmpty()
        {
            var tree = SuffixTree.Build("banana");
            var result = tree.FindExactMatchAnchors("banana", 0);
            Assert.That(result, Is.Empty);
        }

        [Test]
        public void FindAnchors_MinLengthNegative_ReturnsEmpty()
        {
            var tree = SuffixTree.Build("banana");
            var result = tree.FindExactMatchAnchors("banana", -1);
            Assert.That(result, Is.Empty);
        }

        #endregion

        #region Basic Anchor Detection

        [Test]
        public void FindAnchors_IdenticalStrings_ReturnsOneAnchorCoveringAll()
        {
            var tree = SuffixTree.Build("abcdef");
            var result = tree.FindExactMatchAnchors("abcdef", 3);

            Assert.That(result.Count, Is.GreaterThanOrEqualTo(1));
            // The single anchor should cover the entire string
            var anchor = result[0];
            Assert.That(anchor.Length, Is.EqualTo(6));
            Assert.That(anchor.PositionInQuery, Is.EqualTo(0));
            Assert.That(anchor.PositionInText, Is.EqualTo(0));
        }

        [Test]
        public void FindAnchors_NoCommonSubstring_ReturnsEmpty()
        {
            var tree = SuffixTree.Build("aaaaa");
            var result = tree.FindExactMatchAnchors("bbbbb", 1);
            Assert.That(result, Is.Empty);
        }

        [Test]
        public void FindAnchors_QuerySubstringOfText_FindsIt()
        {
            var tree = SuffixTree.Build("abcdefghij");
            var result = tree.FindExactMatchAnchors("defgh", 3);

            Assert.That(result.Count, Is.GreaterThanOrEqualTo(1));
            // Must find "defgh" (length 5) starting at position 3 in text
            var anchor = result[0];
            Assert.That(anchor.Length, Is.EqualTo(5));
            Assert.That(anchor.PositionInText, Is.EqualTo(3));
            Assert.That(anchor.PositionInQuery, Is.EqualTo(0));
        }

        #endregion

        #region Multiple Anchors

        [Test]
        public void FindAnchors_TwoDisjointMatches_FindsBoth()
        {
            // Text: "abcXYZdef" — query "abcPPPdef" shares "abc" and "def"
            var tree = SuffixTree.Build("abcXYZdef");
            var result = tree.FindExactMatchAnchors("abcPPPdef", 3);

            Assert.That(result.Count, Is.EqualTo(2));

            // First anchor: "abc"
            Assert.That(result[0].Length, Is.GreaterThanOrEqualTo(3));
            ValidateAnchor("abcXYZdef", "abcPPPdef", result[0]);

            // Second anchor: "def"
            Assert.That(result[1].Length, Is.GreaterThanOrEqualTo(3));
            ValidateAnchor("abcXYZdef", "abcPPPdef", result[1]);
        }

        [Test]
        public void FindAnchors_MinLengthFiltersShortMatches()
        {
            // "abXXcdXXef" vs "abYYcdYYef" — "ab", "cd", "ef" are length 2
            var tree = SuffixTree.Build("abXXcdXXef");
            var shortResult = tree.FindExactMatchAnchors("abYYcdYYef", 2);
            var longResult = tree.FindExactMatchAnchors("abYYcdYYef", 3);

            Assert.That(shortResult.Count, Is.GreaterThan(0), "minLength=2 should find anchors");
            Assert.That(longResult.Count, Is.EqualTo(0), "minLength=3 should filter out 2-char matches");
        }

        #endregion

        #region Peak Tracking

        [Test]
        public void FindAnchors_ReportsPeakNotAllPositions()
        {
            // For "abcabc" with query "abcabc", the entire string is one contiguous match
            var tree = SuffixTree.Build("abcabc");
            var result = tree.FindExactMatchAnchors("abcabc", 2);

            // Should be a single anchor covering the peak match
            Assert.That(result.Count, Is.EqualTo(1));
            Assert.That(result[0].Length, Is.EqualTo(6));
        }

        [Test]
        public void FindAnchors_TrailingMatch_Emitted()
        {
            // Trailing match at end of query should still be emitted
            var tree = SuffixTree.Build("XXXabc");
            var result = tree.FindExactMatchAnchors("YYYabc", 3);

            Assert.That(result.Count, Is.GreaterThanOrEqualTo(1));
            // The trailing "abc" must appear
            var lastAnchor = result[result.Count - 1];
            Assert.That(lastAnchor.Length, Is.GreaterThanOrEqualTo(3));
            ValidateAnchor("XXXabc", "YYYabc", lastAnchor);
        }

        #endregion

        #region Anchor Validity Invariant

        /// <summary>
        /// Every reported anchor must satisfy the fundamental contract:
        /// text[posInText..posInText+len) == query[posInQuery..posInQuery+len)
        /// </summary>
        [Test]
        [TestCase("banana", "bananarama", 3)]
        [TestCase("mississippi", "missouri", 3)]
        [TestCase("abcabxabcd", "abcabcabc", 2)]
        [TestCase("abracadabra", "abracadabra", 4)]
        [TestCase("aaaaabbbbb", "bbbbbaaaaa", 3)]
        public void FindAnchors_AllAnchorsAreValidExactMatches(string text, string query, int minLen)
        {
            var tree = SuffixTree.Build(text);
            var result = tree.FindExactMatchAnchors(query, minLen);

            foreach (var anchor in result)
            {
                ValidateAnchor(text, query, anchor);
                Assert.That(anchor.Length, Is.GreaterThanOrEqualTo(minLen),
                    $"Anchor length {anchor.Length} below minLength {minLen}");
            }
        }

        [Test]
        [TestCaseSource(nameof(RandomSeeds))]
        public void FindAnchors_RandomInput_AllAnchorsValid(int seed)
        {
            var random = new Random(seed);
            int textLen = random.Next(50, 200);
            int queryLen = random.Next(50, 200);
            string text = GenerateRandomString(random, textLen, "abcd");
            string query = GenerateRandomString(random, queryLen, "abcd");

            var tree = SuffixTree.Build(text);
            var result = tree.FindExactMatchAnchors(query, 3);

            foreach (var anchor in result)
            {
                ValidateAnchor(text, query, anchor);
            }
        }

        #endregion

        #region Suffix Link Traversal

        [Test]
        public void FindAnchors_RequiresSuffixLinkBacktrack()
        {
            // "xabxac" — suffix links from "xab..." to "ab..." are critical
            var tree = SuffixTree.Build("xabxac");
            var result = tree.FindExactMatchAnchors("xabYac", 2);

            // Should find at least "xab" and "ac" as anchors
            Assert.That(result.Count, Is.GreaterThanOrEqualTo(1));
            foreach (var anchor in result)
                ValidateAnchor("xabxac", "xabYac", anchor);
        }

        [Test]
        public void FindAnchors_DeepSuffixLinkChain()
        {
            // "abcabxabcd" has deep suffix links: abc->bc->c
            var tree = SuffixTree.Build("abcabxabcd");
            var result = tree.FindExactMatchAnchors("abcabcabc", 2);

            foreach (var anchor in result)
                ValidateAnchor("abcabxabcd", "abcabcabc", anchor);
        }

        #endregion

        #region Helpers

        private static void ValidateAnchor(string text, string query,
            (int PositionInText, int PositionInQuery, int Length) anchor)
        {
            Assert.That(anchor.PositionInText, Is.GreaterThanOrEqualTo(0));
            Assert.That(anchor.PositionInQuery, Is.GreaterThanOrEqualTo(0));
            Assert.That(anchor.Length, Is.GreaterThan(0));
            Assert.That(anchor.PositionInText + anchor.Length, Is.LessThanOrEqualTo(text.Length),
                $"Anchor overflows text: pos={anchor.PositionInText}, len={anchor.Length}, textLen={text.Length}");
            Assert.That(anchor.PositionInQuery + anchor.Length, Is.LessThanOrEqualTo(query.Length),
                $"Anchor overflows query: pos={anchor.PositionInQuery}, len={anchor.Length}, queryLen={query.Length}");

            string fromText = text.Substring(anchor.PositionInText, anchor.Length);
            string fromQuery = query.Substring(anchor.PositionInQuery, anchor.Length);
            Assert.That(fromText, Is.EqualTo(fromQuery),
                $"Anchor mismatch: text[{anchor.PositionInText}..+{anchor.Length}]=\"{fromText}\" != " +
                $"query[{anchor.PositionInQuery}..+{anchor.Length}]=\"{fromQuery}\"");
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

