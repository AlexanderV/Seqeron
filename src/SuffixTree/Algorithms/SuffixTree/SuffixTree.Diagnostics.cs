using System;
using System.Collections.Generic;
using System.Text;

namespace SuffixTree;

public partial class SuffixTree
{
    private (int nodeCount, int leafCount, int maxDepth, SuffixTreeNode? deepestInternalNode, int maxInternalDepth) CalculateTreeStatisticsInternal()
    {
        if (_text == null || _text.Length == 0) return (1, 0, 0, null, 0);

        var stack = new Stack<(SuffixTreeNode Node, int Depth, bool Visited)>();
        stack.Push((_root, 0, false));

        var childBuffer = new List<SuffixTreeNode>(8);

        int nodeCount = 0;
        int leafCount = 0;
        int maxDepth = 0;
        int maxInternalDepth = 0;
        SuffixTreeNode? deepestInternalNode = null;

        while (stack.Count > 0)
        {
            var (node, depth, visited) = stack.Pop();

            if (visited)
            {
                nodeCount++;

                int nodeFullDepth = depth + (node == _root ? 0 : LengthOf(node));
                if (node.IsLeaf)
                {
                    node.LeafCount = 1;
                    leafCount++;
                    if (nodeFullDepth > maxDepth) maxDepth = nodeFullDepth;
                }
                else
                {
                    int sum = 0;
                    node.GetChildren(childBuffer);
                    foreach (var child in childBuffer)
                        sum += child.LeafCount;
                    node.LeafCount = sum;

                    // Check if this is the deepest internal node for LRS
                    if (nodeFullDepth > maxInternalDepth)
                    {
                        maxInternalDepth = nodeFullDepth;
                        deepestInternalNode = node;
                    }
                }
            }
            else
            {
                stack.Push((node, depth, true));

                if (!node.IsLeaf)
                {
                    int nodeDepth = depth + (node == _root ? 0 : LengthOf(node));
                    node.GetChildren(childBuffer);
                    for (int k = childBuffer.Count - 1; k >= 0; k--)
                        stack.Push((childBuffer[k], nodeDepth, false));
                }
            }
        }

        return (
            nodeCount,
            leafCount > 0 ? leafCount - 1 : 0,
            maxDepth > 0 ? maxDepth - 1 : 0,
            deepestInternalNode,
            maxInternalDepth
        );
    }

    /// <inheritdoc />
    public string PrintTree()
    {
        int estimatedNodes = Math.Max(1, _text!.Length * 2);
        var sb = new StringBuilder(Math.Max(256, estimatedNodes * 50));
        var ci = System.Globalization.CultureInfo.InvariantCulture;
        sb.Append(ci, $"Content length: {_text!.Length}").AppendLine();
        sb.AppendLine();

        var stack = new Stack<(SuffixTreeNode Node, int Depth, int ChildIndex, List<SuffixTreeNode> SortedChildren)>();
        var keyBuffer = new List<int>(8);

        sb.Append(ci, $"0:ROOT").AppendLine();

        var rootChildren = GetSortedChildren(_root, keyBuffer);
        if (rootChildren.Count > 0)
            stack.Push((_root, 0, 0, rootChildren));

        while (stack.Count > 0)
        {
            var (node, depth, childIndex, sortedChildren) = stack.Pop();

            if (childIndex < sortedChildren.Count)
            {
                stack.Push((node, depth, childIndex + 1, sortedChildren));

                var child = sortedChildren[childIndex];
                int childDepth = depth + 1;

                var nodeLabel = LabelOf(child);
                var leafMark = child.IsLeaf ? " (Leaf)" : "";
                var linkMark = child.SuffixLink != null && child.SuffixLink != _root && !child.IsLeaf
                    ? $" -> [{FirstCharOf(child.SuffixLink)}]"
                    : "";
                sb.Append(' ', childDepth * 2);
                sb.Append(ci, $"{childDepth}: {nodeLabel}{leafMark}{linkMark}").AppendLine();

                if (!child.IsLeaf && child.ChildCount > 0)
                {
                    var grandChildren = GetSortedChildren(child, keyBuffer);
                    stack.Push((child, childDepth, 0, grandChildren));
                }
            }
        }

        return sb.ToString();
    }

    private static List<SuffixTreeNode> GetSortedChildren(SuffixTreeNode node, List<int> keyBuffer)
    {
        node.GetKeys(keyBuffer);
        keyBuffer.Sort();
        var result = new List<SuffixTreeNode>(keyBuffer.Count);
        foreach (var key in keyBuffer)
        {
            if (node.TryGetChild(key, out var child))
                result.Add(child!);
        }
        return result;
    }

    /// <summary>
    /// Validates tree integrity by checking suffix links.
    /// Useful for debugging algorithm changes.
    /// </summary>
    internal void ValidateSuffixLinks()
    {
        var stack = new Stack<SuffixTreeNode>();
        stack.Push(_root);
        var visited = new HashSet<SuffixTreeNode>();
        var buffer = new List<SuffixTreeNode>();

        while (stack.Count > 0)
        {
            var node = stack.Pop();
            if (!visited.Add(node)) continue;

            if (!node.IsLeaf && node != _root)
            {
                if (node.SuffixLink == null)
                    throw new InvalidOperationException($"Node at depth {node.DepthFromRoot} missing suffix link");
            }

            node.GetChildren(buffer);
            foreach (var child in buffer)
                stack.Push(child);
        }
    }
    private string LabelOf(SuffixTreeNode node)
    {
        if (node == _root) return "ROOT";
        var sb = new StringBuilder();
        int len = LengthOf(node);
        for (int i = 0; i < len; i++)
        {
            int s = GetSymbolAt(node.Start + i);
            if (s == TERMINATOR_KEY)
            {
                sb.Append('#');
                break;
            }
            sb.Append((char)s);
        }
        return sb.ToString();
    }

    /// <inheritdoc/>
    public void Traverse(ISuffixTreeVisitor visitor)
    {
        ArgumentNullException.ThrowIfNull(visitor);
        visitor.VisitNode(_root.Start, _root.End, _root.LeafCount, _root.ChildCount, 0);
        if (_root.IsLeaf)
            return;

        var rootKeys = new List<int>(_root.ChildCount);
        _root.GetKeys(rootKeys);
        rootKeys.Sort();

        var stack = new Stack<(SuffixTreeNode Node, int Depth, List<int> Keys, int Index, bool IsRoot)>();
        stack.Push((_root, 0, rootKeys, 0, true));

        while (stack.Count > 0)
        {
            var (node, depth, keys, index, isRoot) = stack.Pop();

            if (index >= keys.Count)
            {
                if (!isRoot)
                    visitor.ExitBranch();
                continue;
            }

            int key = keys[index];
            stack.Push((node, depth, keys, index + 1, isRoot));

            if (!node.TryGetChild(key, out var child) || child is null)
                continue;

            visitor.EnterBranch(key);

            int childDepth = depth + (node == _root ? 0 : LengthOf(node));
            visitor.VisitNode(child.Start, child.End, child.LeafCount, child.ChildCount, childDepth);

            if (child.IsLeaf)
            {
                visitor.ExitBranch();
                continue;
            }

            var childKeys = new List<int>(child.ChildCount);
            child.GetKeys(childKeys);
            childKeys.Sort();
            stack.Push((child, childDepth, childKeys, 0, false));
        }
    }
}
