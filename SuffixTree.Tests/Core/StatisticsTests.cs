using System;
using System.Linq;
using NUnit.Framework;

namespace SuffixTree.Tests.Core
{
    /// <summary>
    /// Tests for tree statistics: NodeCount, LeafCount.
    /// </summary>
    [TestFixture]
    public class StatisticsTests
    {
        #region Empty Tree

        [Test]
        public void Stats_EmptyTree_HasMinimalCounts()
        {
            var st = SuffixTree.Build("");

            Assert.Multiple(() =>
            {
                // Root always exists
                Assert.That(st.NodeCount, Is.GreaterThanOrEqualTo(1));
                Assert.That(st.LeafCount, Is.GreaterThanOrEqualTo(0));
            });
        }

        #endregion

        #region Single Character

        [Test]
        public void Stats_SingleChar_BasicStructure()
        {
            var st = SuffixTree.Build("a");

            Assert.Multiple(() =>
            {
                // Root + one leaf minimum
                Assert.That(st.NodeCount, Is.GreaterThanOrEqualTo(2));
                Assert.That(st.LeafCount, Is.GreaterThanOrEqualTo(1));
            });
        }

        #endregion

        #region Consistency

        [Test]
        public void Stats_NodeCountIncludesLeaves()
        {
            var testCases = new[] { "a", "ab", "abc", "banana", "mississippi" };

            foreach (var text in testCases)
            {
                var st = SuffixTree.Build(text);

                // NodeCount should include leaves
                Assert.That(st.NodeCount, Is.GreaterThanOrEqualTo(st.LeafCount),
                    $"NodeCount should be >= LeafCount for '{text}'");
            }
        }

        [Test]
        public void Stats_LeafCountEqualsTextLength()
        {
            var testCases = new[] { "a", "ab", "abc", "banana", "mississippi" };

            foreach (var text in testCases)
            {
                var st = SuffixTree.Build(text);

                // One leaf per suffix
                Assert.That(st.LeafCount, Is.EqualTo(text.Length),
                    $"LeafCount should equal text length for '{text}'");
            }
        }

        #endregion

        #region Growth Pattern

        [Test]
        public void Stats_RepeatingText_HasNodes()
        {
            var stUnique = SuffixTree.Build("abcdefghij");
            var stRepeating = SuffixTree.Build("ababababab");

            // Both trees have nodes > 0
            Assert.Multiple(() =>
            {
                Assert.That(stUnique.NodeCount, Is.GreaterThan(0));
                Assert.That(stRepeating.NodeCount, Is.GreaterThan(0));
            });
        }

        #endregion

        #region Known Structure

        [Test]
        public void Stats_AllSameCharacter_ReasonableNodeCount()
        {
            var st = SuffixTree.Build("aaaaa");

            // For "aaaa...", tree is valid and has some structure
            // NodeCount should be reasonable (depends on implementation)
            Assert.Multiple(() =>
            {
                Assert.That(st.NodeCount, Is.GreaterThan(0));
                Assert.That(st.LeafCount, Is.EqualTo(5));
            });
        }

        #endregion

        #region Immutability

        [Test]
        public void Stats_DoNotChangeAfterQueries()
        {
            var st = SuffixTree.Build("banana");

            var nodesBefore = st.NodeCount;
            var leafBefore = st.LeafCount;

            // Perform various operations
            _ = st.Contains("ana");
            _ = st.FindAllOccurrences("a").ToList();
            _ = st.LongestRepeatedSubstring();
            _ = st.LongestCommonSubstring("other");

            Assert.Multiple(() =>
            {
                Assert.That(st.NodeCount, Is.EqualTo(nodesBefore));
                Assert.That(st.LeafCount, Is.EqualTo(leafBefore));
            });
        }

        #endregion
    }
}
