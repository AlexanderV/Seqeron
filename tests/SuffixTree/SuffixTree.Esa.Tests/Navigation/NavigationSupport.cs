using SuffixTree.Esa.Tests.TestSupport;
using UkkonenTree = global::SuffixTree.SuffixTree;

namespace SuffixTree.Esa.Tests.Navigation;

/// <summary>One node of a tree walk in preorder, in a form that does not depend on the engine: the path label, the leaf flag, the leaf count, the child count and the key of the branch that led to it.</summary>
internal sealed record TreeRecord(string Path, bool IsLeaf, int LeafCount, int ChildCount, int EnterKey);

/// <summary>The shape of the Ukkonen tree (public <c>Traverse</c> events) as <see cref="TreeRecord"/>s: the second oracle of the navigator walk (the first is the definition).</summary>
internal static class UkkonenShape
{
    private sealed class Visitor(string text) : ISuffixTreeVisitor
    {
        private readonly Stack<string> _paths = new();
        private int _key;
        private bool _root = true;

        internal List<TreeRecord> Records { get; } = [];

        public void VisitNode(int startIndex, int endIndex, int leafCount, int childCount, int depth)
        {
            if (_root)
            {
                _root = false;
                _paths.Push("");
                Records.Add(new TreeRecord("", childCount == 0, leafCount, childCount, 0));
                return;
            }
            string edge = endIndex < 0 ? text.Substring(Math.Min(startIndex, text.Length)) : text.Substring(startIndex, endIndex - startIndex);
            string path = _paths.Peek() + edge;
            _paths.Push(path);
            Records.Add(new TreeRecord(path, childCount == 0, leafCount, childCount, _key));
        }

        public void EnterBranch(int key) => _key = key;

        public void ExitBranch() => _paths.Pop();
    }

    internal static List<TreeRecord> Of(string text)
    {
        var visitor = new Visitor(text);
        UkkonenTree.Build(text).Traverse(visitor);
        return visitor.Records;
    }
}

/// <summary>
/// Walks the ESA tree through <see cref="EsaNavigator"/> only (children by first edge symbol, terminator first) and checks, node by node, everything the generic algorithms rely on:
/// parent depth, edge length and symbols, <c>TryGetChild</c> for every present and some absent keys, <c>FindAnyLeafPosition</c> and <c>CollectLeaves</c>. Returns the preorder
/// <see cref="TreeRecord"/>s for the comparison with the Ukkonen tree.
/// </summary>
internal static class NavigatorWalk
{
    internal static List<TreeRecord> Walk(EsaIndex index, bool bruteLeaves, string context)
    {
        string text = index.Text;
        int n = text.Length;
        EsaNavigator nav = index.CreateNavigator();
        var records = new List<TreeRecord>();
        EsaNode root = nav.Root;
        Assert.That(nav.IsRoot(root), $"{context}: Root is the root");
        Assert.That(nav.IsNull(root), Is.False, $"{context}: Root is not null");
        Assert.That(nav.IsNull(nav.NullNode), $"{context}: NullNode is null");
        Assert.That(nav.IsRoot(nav.NullNode), Is.False, $"{context}: NullNode is not the root");
        Assert.That(nav.LengthOf(root), Is.EqualTo(0), $"{context}: the root has no edge");
        Assert.That(root.Lb, Is.EqualTo(0), context);
        Assert.That(root.Rb, Is.EqualTo(n), context);
        if (n == 0)
        {
            records.Add(new TreeRecord("", true, 0, 0, 0));
            var none = new List<EsaNode> { root };
            nav.GetChildren(root, none);
            Assert.That(none, Is.Empty, $"{context}: the root of the empty text has no children and GetChildren replaces the list");
            Assert.That(nav.TryGetChild(root, -1, out _), Is.False, $"{context}: no child in the empty tree");
            return records;
        }

        var probeKeys = text.Distinct().Select(c => (int)c).Append(-1).Append(text.Max() + 1).ToList();
        var stack = new Stack<(EsaNode Node, int Key, string Path)>();
        stack.Push((root, 0, ""));
        var children = new List<EsaNode>();
        while (stack.Count > 0)
        {
            var (node, key, path) = stack.Pop();
            nav.GetChildren(node, children);
            var kids = children.OrderBy(c => nav.GetEdgeSymbol(c, 0)).ToList();
            bool leaf = node.IsLeaf;
            Assert.That(kids.Count == 0, Is.EqualTo(leaf), $"{context}: node '{path}' leaf <=> no children");
            records.Add(new TreeRecord(path, leaf, node.Rb - node.Lb + 1, kids.Count, key));
            CheckNode(nav, text, node, path, bruteLeaves, context);
            if (leaf) continue;

            var seen = new HashSet<int>();
            foreach (EsaNode child in kids)
            {
                Assert.That(child.ParentDepth, Is.EqualTo(node.Depth), $"{context}: parent depth of a child of '{path}'");
                Assert.That(child.Depth, Is.GreaterThan(node.Depth), $"{context}: child of '{path}' is deeper");
                Assert.That(nav.LengthOf(child), Is.EqualTo(child.Depth - node.Depth), context);
                Assert.That(child.Lb, Is.GreaterThanOrEqualTo(node.Lb), context);
                Assert.That(child.Rb, Is.LessThanOrEqualTo(node.Rb), context);
                int firstSymbol = nav.GetEdgeSymbol(child, 0);
                Assert.That(seen.Add(firstSymbol), $"{context}: two children of '{path}' share the key {firstSymbol}");
                string childPath = child.IsLeaf ? text.Substring(Math.Min(child.Pos, n)) : text.Substring(child.Pos, child.Depth);
                Assert.That(childPath.StartsWith(path, StringComparison.Ordinal), $"{context}: '{childPath}' extends '{path}'");
                Assert.That(nav.TryGetChild(node, firstSymbol, out EsaNode found), $"{context}: TryGetChild('{path}', {firstSymbol})");
                Assert.That(found, Is.EqualTo(child), $"{context}: TryGetChild returns the child that GetChildren lists");
            }
            foreach (int probe in probeKeys)
            {
                if (seen.Contains(probe)) continue;
                Assert.That(nav.TryGetChild(node, probe, out EsaNode absent), Is.False, $"{context}: '{path}' has no child with key {probe}");
                Assert.That(nav.IsNull(absent), $"{context}: a failed TryGetChild yields the null node");
            }
            for (int i = kids.Count - 1; i >= 0; i--)
            {
                EsaNode child = kids[i];
                string childPath = child.IsLeaf ? text.Substring(Math.Min(child.Pos, n)) : text.Substring(child.Pos, child.Depth);
                stack.Push((child, nav.GetEdgeSymbol(child, 0), childPath));
            }
        }
        return records;
    }

    private static void CheckNode(EsaNavigator nav, string text, EsaNode node, string path, bool bruteLeaves, string context)
    {
        int n = text.Length;
        bool isRoot = node.Depth == 0;
        Assert.That(nav.IsRoot(node), Is.EqualTo(isRoot), $"{context}: IsRoot('{path}')");
        Assert.That(nav.IsNull(node), Is.False, context);
        // edge symbols: the path label past the parent depth, then the terminator for a leaf
        for (int o = 0; o < node.Depth - node.ParentDepth; o++)
        {
            int at = node.ParentDepth + o;
            int expected = at < path.Length ? path[at] : -1;
            Assert.That(nav.GetEdgeSymbol(node, o), Is.EqualTo(expected), $"{context}: edge symbol {o} of '{path}'");
        }
        if (node.IsLeaf) Assert.That(nav.GetEdgeSymbol(node, node.Depth - node.ParentDepth - 1), Is.EqualTo(-1), $"{context}: a leaf edge ends with the terminator");
        else if (!isRoot) Assert.That(nav.GetEdgeSymbol(node, node.Depth - node.ParentDepth - 1), Is.GreaterThanOrEqualTo(0), $"{context}: an internal edge never holds the terminator");

        var leaves = new List<int> { -7 };
        nav.CollectLeaves(node, node.ParentDepth, leaves);
        Assert.That(leaves[0], Is.EqualTo(-7), $"{context}: CollectLeaves appends");
        leaves.RemoveAt(0);
        int expectedCount = node.IsLeaf ? (path.Length > 0 ? 1 : 0) : (isRoot ? n : node.Rb - node.Lb + 1);
        Assert.That(leaves.Count, Is.EqualTo(expectedCount), $"{context}: CollectLeaves('{path}') count");
        Assert.That(leaves.Distinct().Count(), Is.EqualTo(leaves.Count), $"{context}: leaf positions are distinct");
        foreach (int p in leaves)
        {
            Assert.That((uint)p, Is.LessThan((uint)n), context);
            Assert.That(string.CompareOrdinal(text, p, path, 0, path.Length) == 0 && (!node.IsLeaf || p + path.Length == n), $"{context}: leaf position {p} of '{path}'");
        }
        if (bruteLeaves && !node.IsLeaf)
            Assert.That(leaves.Order().ToList(), Is.EqualTo(BruteForce.Occurrences(text, path)), $"{context}: leaves of '{path}' against brute force");

        int any = nav.FindAnyLeafPosition(node, node.ParentDepth);
        if (expectedCount == 0) Assert.That(any, Is.EqualTo(-1), $"{context}: no real suffix below '{path}'");
        else Assert.That(leaves, Does.Contain(any), $"{context}: FindAnyLeafPosition('{path}') is one of the leaves");
    }

    /// <summary>The walk of <paramref name="index"/> and the Ukkonen shape of the same text must be equal record by record.</summary>
    internal static void AssertEqualsUkkonen(EsaIndex index, bool bruteLeaves, string context)
    {
        List<TreeRecord> esa = Walk(index, bruteLeaves, context);
        List<TreeRecord> ukk = UkkonenShape.Of(index.Text);
        Assert.That(esa.Count, Is.EqualTo(ukk.Count), $"{context}: node count");
        for (int i = 0; i < esa.Count; i++)
            if (!esa[i].Equals(ukk[i]))
                Assert.Fail($"{context}: node {i}: ESA {esa[i]} against Ukkonen {ukk[i]}");
    }
}
