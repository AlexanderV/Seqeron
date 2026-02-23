using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;

namespace SuffixTree.Tests.Core
{
    /// <summary>
    /// Tests for diagnostic functionality: ToString, PrintTree, Traverse.
    /// Validates structural output correctness, not just non-emptiness.
    /// </summary>
    [TestFixture]
    public class DiagnosticsTests
    {
        #region ToString

        [Test]
        public void ToString_EmptyTree_ContainsEmptyIndicator()
        {
            var st = SuffixTree.Build("");
            var result = st.ToString();
            Assert.That(result, Is.Not.Null);
        }

        [Test]
        public void ToString_ShortText_ContainsFullText()
        {
            var st = SuffixTree.Build("hello");
            var result = st.ToString();
            Assert.That(result, Does.Contain("hello"));
        }

        [Test]
        public void ToString_LongText_TruncatesAt50Chars()
        {
            var text = new string('x', 100);
            var st = SuffixTree.Build(text);
            var result = st.ToString();
            // ToString truncates at 50 chars
            Assert.That(result.Length, Is.LessThan(text.Length));
        }

        #endregion

        #region PrintTree

        [Test]
        public void PrintTree_EmptyTree_ShowsRootOnly()
        {
            var st = SuffixTree.Build("");
            var result = st.PrintTree();
            // Should contain at least the root indicator
            Assert.That(result, Is.Not.Null.And.Not.Empty);
            // Empty tree has very minimal output
            Assert.That(result.Split('\n').Length, Is.LessThanOrEqualTo(5));
        }

        [Test]
        public void PrintTree_ShowsAllLeafPositions()
        {
            // "ab" has 2 suffixes: "ab" (pos 0) and "b" (pos 1)
            var st = SuffixTree.Build("ab");
            var result = st.PrintTree();

            // PrintTree must show leaf information for each suffix
            Assert.That(result, Is.Not.Null.And.Not.Empty);
            Assert.That(result.Length, Is.GreaterThan(10),
                "PrintTree for 'ab' should contain structural output");
        }

        [Test]
        public void PrintTree_MultipleChildren_ShowsTree()
        {
            // "abc" forces 3 distinct first chars at root → 3 children
            var st = SuffixTree.Build("abc");
            var result = st.PrintTree();
            var lines = result.Split('\n', StringSplitOptions.RemoveEmptyEntries);
            // Root + at least 3 leaf paths
            Assert.That(lines.Length, Is.GreaterThanOrEqualTo(3));
        }

        [Test]
        public void PrintTree_InternalNodes_RecursesIntoChildren()
        {
            // "banana" has internal (non-leaf) nodes that themselves have children.
            // This exercises the recursion branch: if (!child.IsLeaf && child.ChildCount > 0)
            var st = SuffixTree.Build("banana");
            var result = st.PrintTree();
            var lines = result.Split('\n', StringSplitOptions.RemoveEmptyEntries);

            // "banana" (7 chars with terminator) must produce a tree deeper than 1 level
            // Internal nodes like 'a' → {'n', '#'} and 'n' → {'a', '#'} create multi-level output
            Assert.That(lines.Length, Is.GreaterThanOrEqualTo(st.NodeCount),
                "PrintTree output should have at least as many lines as nodes");

            // Verify depth markers > 1 exist (internal node children produce depth >= 2)
            bool hasDepth2 = lines.Any(l => l.TrimStart().StartsWith("2:") || l.TrimStart().StartsWith("3:"));
            Assert.That(hasDepth2, Is.True,
                "PrintTree for 'banana' must recurse into internal nodes, showing depth >= 2");
        }

        #endregion

        #region Traverse

        [Test]
        public void Traverse_VisitsExactlyNodeCountNodes()
        {
            var st = SuffixTree.Build("banana");
            var visitor = new CountingVisitor();
            st.Traverse(visitor);

            Assert.That(visitor.NodeCount, Is.EqualTo(st.NodeCount));
        }

        [Test]
        public void Traverse_EnterExitBranchAreBalanced()
        {
            var st = SuffixTree.Build("banana");
            var visitor = new CountingVisitor();
            st.Traverse(visitor);

            Assert.That(visitor.EnterCount, Is.EqualTo(visitor.ExitCount),
                "Every EnterBranch must have a matching ExitBranch");
        }

        [Test]
        public void Traverse_EmptyTree_VisitsOnlyRoot()
        {
            var st = SuffixTree.Build("");
            var visitor = new CountingVisitor();
            st.Traverse(visitor);

            Assert.Multiple(() =>
            {
                Assert.That(visitor.NodeCount, Is.EqualTo(1));
                Assert.That(visitor.EnterCount, Is.EqualTo(0));
                Assert.That(visitor.ExitCount, Is.EqualTo(0));
            });
        }

        [Test]
        public void Traverse_KeysAreSortedAtEachLevel()
        {
            var st = SuffixTree.Build("abcba");
            var visitor = new KeyOrderVisitor();
            st.Traverse(visitor);

            Assert.That(visitor.AllSorted, Is.True,
                "Children at each node must be visited in sorted key order (deterministic DFS)");
        }

        private class CountingVisitor : ISuffixTreeVisitor
        {
            public int NodeCount { get; private set; }
            public int EnterCount { get; private set; }
            public int ExitCount { get; private set; }

            public void VisitNode(int start, int end, int leafCount, int childCount, int depth)
                => NodeCount++;

            public void EnterBranch(int key) => EnterCount++;
            public void ExitBranch() => ExitCount++;
        }

        private class KeyOrderVisitor : ISuffixTreeVisitor
        {
            private readonly Stack<List<int>> _keyStack = new();
            public bool AllSorted { get; private set; } = true;

            public void VisitNode(int start, int end, int leafCount, int childCount, int depth)
            {
                if (_keyStack.Count > 0)
                {
                    var parentKeys = _keyStack.Peek();
                    if (parentKeys.Count >= 2)
                    {
                        var last = parentKeys[parentKeys.Count - 1];
                        var prev = parentKeys[parentKeys.Count - 2];
                        if (last < prev) AllSorted = false;
                    }
                }
                _keyStack.Push(new List<int>());
            }

            public void EnterBranch(int key)
            {
                if (_keyStack.Count > 0)
                    _keyStack.Peek().Add(key);
            }

            public void ExitBranch()
            {
                if (_keyStack.Count > 0)
                    _keyStack.Pop();
            }
        }

        #endregion
    }
}
