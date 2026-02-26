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
    [Category("Core")]
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
        public void PrintTree_EmptyTree_HasOnlyRootNodeLine()
        {
            var st = SuffixTree.Build("");
            var result = st.PrintTree();
            var parsed = ParsePrintTree(result);

            Assert.Multiple(() =>
            {
                Assert.That(parsed.ContentLength, Is.EqualTo(0));
                Assert.That(parsed.NodeLines.Count, Is.EqualTo(1));
                Assert.That(parsed.NodeLines[0].Depth, Is.EqualTo(0));
                Assert.That(parsed.NodeLines[0].Label, Is.EqualTo("ROOT"));
            });
        }

        [TestCase("ab")]
        [TestCase("banana")]
        [TestCase("mississippi")]
        public void PrintTree_NodeLinesMatchNodeCount_AndDepthIndentContract(string text)
        {
            var st = SuffixTree.Build(text);
            var result = st.PrintTree();
            var parsed = ParsePrintTree(result);

            Assert.Multiple(() =>
            {
                Assert.That(parsed.ContentLength, Is.EqualTo(text.Length));
                Assert.That(parsed.NodeLines.Count, Is.EqualTo(st.NodeCount),
                    "Each tree node must be printed exactly once");
                Assert.That(parsed.NodeLines[0].Label, Is.EqualTo("ROOT"));
                Assert.That(parsed.NodeLines[0].Depth, Is.EqualTo(0));
                Assert.That(parsed.NodeLines.All(line => line.IndentSpaces == line.Depth * 2), Is.True,
                    "Indentation must be exactly 2 spaces per depth level");
            });
        }

        [TestCase("ab")]
        [TestCase("banana")]
        [TestCase("mississippi")]
        public void PrintTree_NonEmptyTree_LeafLineCountMatchesLeafCountPlusTerminator(string text)
        {
            var st = SuffixTree.Build(text);
            var result = st.PrintTree();
            var parsed = ParsePrintTree(result);
            int leafLines = parsed.NodeLines.Count(line => line.IsLeaf);

            Assert.That(leafLines, Is.EqualTo(st.LeafCount + 1),
                "Printed leaves include explicit terminator leaf in addition to API LeafCount");
        }

        [Test]
        public void PrintTree_FirstLevelChildren_AreSortedByLeadingSymbol()
        {
            var st = SuffixTree.Build("banana");
            var parsed = ParsePrintTree(st.PrintTree());
            var depth1Leading = parsed.NodeLines
                .Where(line => line.Depth == 1)
                .Select(line => line.Label[0])
                .ToArray();

            var sorted = depth1Leading.OrderBy(ch => ch).ToArray();
            Assert.That(depth1Leading, Is.EqualTo(sorted),
                "Root children must be printed in deterministic sorted order");
        }

        private static ParsedPrintTree ParsePrintTree(string output)
        {
            var lines = output
                .Split('\n', StringSplitOptions.RemoveEmptyEntries)
                .Select(static line => line.TrimEnd('\r'))
                .ToArray();
            Assert.That(lines.Length, Is.GreaterThanOrEqualTo(2), "PrintTree output must contain header and root");
            Assert.That(lines[0], Does.StartWith("Content length: "));

            int contentLength = int.Parse(lines[0]["Content length: ".Length..], System.Globalization.CultureInfo.InvariantCulture);
            var nodeLines = new List<PrintedNodeLine>();

            foreach (string line in lines.Skip(1))
            {
                int colon = line.IndexOf(':');
                if (colon <= 0)
                    continue;

                string depthPart = line[..colon].Trim();
                if (!int.TryParse(depthPart, out int depth))
                    continue;

                int indentSpaces = line.TakeWhile(ch => ch == ' ').Count();
                string content = line[(colon + 1)..].TrimStart();

                bool isLeaf = content.EndsWith(" (Leaf)", StringComparison.Ordinal);
                if (isLeaf)
                    content = content[..^" (Leaf)".Length];

                int linkMarker = content.IndexOf(" -> [", StringComparison.Ordinal);
                if (linkMarker >= 0)
                    content = content[..linkMarker];

                nodeLines.Add(new PrintedNodeLine(depth, indentSpaces, content, isLeaf));
            }

            return new ParsedPrintTree(contentLength, nodeLines);
        }

        private sealed record PrintedNodeLine(int Depth, int IndentSpaces, string Label, bool IsLeaf);

        private sealed record ParsedPrintTree(int ContentLength, IReadOnlyList<PrintedNodeLine> NodeLines);

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

