using System;
using System.Linq;
using NUnit.Framework;

namespace SuffixTree.Tests.Regression
{
    /// <summary>
    /// Regression tests for previously found bugs.
    /// Each test documents a specific bug that was fixed.
    /// These tests should NOT be deleted when refactoring - they prevent regressions.
    /// </summary>
    [TestFixture]
    public class RegressionTests
    {
        #region LCS Backtracking (Bug: incorrect LCS when mismatch in mid-edge)

        /// <summary>
        /// Regression: LCS algorithm must properly backtrack when
        /// matching fails in the middle of an edge.
        /// Bug: Algorithm returned wrong LCS because edge offset wasn't reset.
        /// </summary>
        [Test]
        public void LCS_BacktrackMidEdge_FindsCorrect()
        {
            // Tree: "abcdef" has edge "abcdef" from root
            // Pattern: "abxyz" matches "ab" then fails at 'x'
            // Must backtrack to find "ab"
            var st = SuffixTree.Build("abcdef");

            var result = st.LongestCommonSubstring("abxyz");

            Assert.That(result, Is.EqualTo("ab"));
        }

        /// <summary>
        /// Regression: LCS must handle multiple suffix link jumps correctly.
        /// Bug: State was lost after multiple consecutive suffix link traversals.
        /// </summary>
        [Test]
        public void LCS_MultipleSuffixLinkJumps_Works()
        {
            var st = SuffixTree.Build("aaaaaab");

            var result = st.LongestCommonSubstring("aaaaxaaaaa");

            // "aaaaa" appears in both (5 a's in "aaaaaab" and in "aaaaxaaaaa")
            Assert.That(result, Is.EqualTo("aaaaa"));
        }

        /// <summary>
        /// Regression: LCS edge offset must be reset after backtracking.
        /// Bug: edgeOffset variable retained stale value causing incorrect matching.
        /// </summary>
        [Test]
        public void LCS_EdgeOffsetResetAfterBacktrack_Works()
        {
            var st = SuffixTree.Build("xyzabcdef");

            var result = st.LongestCommonSubstring("abcXdef");

            // Should find "abc" or "def" (both length 3)
            Assert.That(result.Length, Is.EqualTo(3));
        }

        #endregion

        #region TERMINATOR_KEY Handling (Bug: null char conflict)

        /// <summary>
        /// Regression: TERMINATOR_KEY (-1) must not conflict with null character (0).
        /// Bug: Original implementation used '\0' as terminator, causing issues
        /// when text contained null characters.
        /// </summary>
        [Test]
        public void TerminatorKey_DoesNotConflictWithNullChar()
        {
            var text = "ab\0cd";
            var st = SuffixTree.Build(text);

            Assert.Multiple(() =>
            {
                Assert.That(st.Contains("ab\0cd"), Is.True);
                Assert.That(st.Contains("\0"), Is.True);
                Assert.That(st.Contains("b\0c"), Is.True);
                Assert.That(st.GetAllSuffixes().Count, Is.EqualTo(5));
            });
        }

        #endregion

        #region Consistency Invariants

        /// <summary>
        /// Regression: CountOccurrences must always equal FindAllOccurrences().Count.
        /// Bug: Precomputed leaf counts were incorrect for some tree structures.
        /// </summary>
        [Test]
        public void Count_EqualsFind_ForAllPatterns()
        {
            var text = "abracadabra";
            var st = SuffixTree.Build(text);

            // Test all possible substrings
            for (int i = 0; i < text.Length; i++)
            {
                for (int len = 1; len <= text.Length - i; len++)
                {
                    var pattern = text.Substring(i, len);
                    var count = st.CountOccurrences(pattern);
                    var findCount = st.FindAllOccurrences(pattern).Count;

                    Assert.That(count, Is.EqualTo(findCount),
                        $"Mismatch for pattern '{pattern}'");
                }
            }
        }

        /// <summary>
        /// Regression: FindAllLongestCommonSubstrings and LongestCommonSubstring
        /// must return the same substring.
        /// Bug: FindAll variant used different traversal logic.
        /// </summary>
        [Test]
        public void FindAllLCS_MatchesSingleLCS()
        {
            var testCases = new[]
            {
                ("hello world", "world hello"),
                ("ABAB", "BABA"),
                ("abcdef", "xyzabc"),
            };

            foreach (var (text, other) in testCases)
            {
                var st = SuffixTree.Build(text);
                var singleLcs = st.LongestCommonSubstring(other);
                var (allLcsSubstring, _, _) = st.FindAllLongestCommonSubstrings(other);

                Assert.That(allLcsSubstring, Is.EqualTo(singleLcs),
                    $"FindAll substring should match single LCS for '{text}' vs '{other}'");
            }
        }

        /// <summary>
        /// Regression: Positions returned by FindAllOccurrences must be accurate.
        /// Bug: Leaf suffix indices were off by terminator position.
        /// </summary>
        [Test]
        public void FindAll_PositionsAreAccurate()
        {
            var text = "abcabcabc";
            var st = SuffixTree.Build(text);
            var pattern = "abc";

            var positions = st.FindAllOccurrences(pattern).ToList();

            foreach (var pos in positions)
            {
                var extracted = text.Substring(pos, pattern.Length);
                Assert.That(extracted, Is.EqualTo(pattern),
                    $"Position {pos} does not yield pattern");
            }
        }

        #endregion
    }
}
