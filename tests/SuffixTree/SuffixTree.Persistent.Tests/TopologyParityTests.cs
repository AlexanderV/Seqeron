using System;
using System.Collections.Generic;
using NUnit.Framework;
using SuffixTree;

namespace SuffixTree.Persistent.Tests
{
    [TestFixture]
    public class TopologyParityTests
    {
        [TestCase("banana")]
        [TestCase("mississippi")]
        [TestCase("aaaaaaaa")]
        [TestCase("abcde")]
        [TestCase("abacaba")]
        [TestCase("")]
        public void Topology_MatchesReference(string text)
        {
            var reference = global::SuffixTree.SuffixTree.Build(text);
            using (var persistent = PersistentSuffixTreeFactory.Create(new StringTextSource(text)) as IDisposable)
            {
                var st = (ISuffixTree)persistent!;

                var refNodes = new List<NodeInfo>();
                var perNodes = new List<NodeInfo>();

                reference.Traverse(new TopologyCollector(refNodes));
                st.Traverse(new TopologyCollector(perNodes));

                Assert.That(perNodes.Count, Is.EqualTo(refNodes.Count), "Node count mismatch");

                for (int i = 0; i < refNodes.Count; i++)
                {
                    Assert.That(perNodes[i], Is.EqualTo(refNodes[i]), $"Node structure mismatch at index {i}");
                }
            }
        }

        private class NodeInfo : IEquatable<NodeInfo>
        {
            public int Start { get; }
            public int End { get; }
            public int LeafCount { get; }
            public int ChildCount { get; }
            public List<int> ChildrenKeys { get; }

            public NodeInfo(int start, int end, int leafCount, int childCount, List<int> childrenKeys)
            {
                Start = start;
                End = end;
                LeafCount = leafCount;
                ChildCount = childCount;
                ChildrenKeys = childrenKeys;
            }

            public bool Equals(NodeInfo? other)
            {
                if (other == null) return false;
                return Start == other.Start &&
                       End == other.End &&
                       LeafCount == other.LeafCount &&
                       ChildCount == other.ChildCount &&
                       ChildrenKeys.SequenceEqual(other.ChildrenKeys);
            }

            public override bool Equals(object? obj) => Equals(obj as NodeInfo);
            public override int GetHashCode() => HashCode.Combine(Start, End, LeafCount, ChildCount);

            public override string ToString() =>
                $"Node(S={Start}, E={End}, L={LeafCount}, C={ChildCount}, Keys=[{string.Join(",", ChildrenKeys)}])";
        }

        private class TopologyCollector : ISuffixTreeVisitor
        {
            private readonly List<NodeInfo> _nodes;
            private readonly Stack<List<int>> _childrenStack = new();

            public TopologyCollector(List<NodeInfo> nodes) => _nodes = nodes;

            public void VisitNode(int startIndex, int endIndex, int leafCount, int childCount, int depth)
            {
                var children = new List<int>();
                _nodes.Add(new NodeInfo(startIndex, endIndex, leafCount, childCount, children));
                _childrenStack.Push(children);
            }

            public void EnterBranch(int key)
            {
                if (_childrenStack.Count > 0)
                {
                    _childrenStack.Peek().Add(key);
                }
            }

            public void ExitBranch()
            {
                if (_childrenStack.Count > 0)
                {
                    _childrenStack.Pop();
                }
            }
        }
    }
}
