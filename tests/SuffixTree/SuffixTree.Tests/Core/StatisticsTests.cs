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
