using System.Collections.Generic;
using System.Runtime.CompilerServices;

namespace SuffixTree;

public partial class SuffixTree
{
    /// <summary>
    /// Struct navigator for the in-memory suffix tree.
    /// Implements <see cref="ISuffixTreeNavigator{TNode}"/> to enable
    /// shared algorithm dispatch with zero overhead (JIT specialization).
    /// </summary>
    internal struct SuffixTreeNavigator : ISuffixTreeNavigator<SuffixTreeNode>
    {
        private readonly SuffixTree _tree;

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public SuffixTreeNavigator(SuffixTree tree) => _tree = tree;

        public ITextSource Text
        {
            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            get => _tree._text;
        }

        public SuffixTreeNode Root
        {
            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            get => _tree._root;
        }

        public SuffixTreeNode NullNode
        {
            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            get => null!;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public bool IsNull(SuffixTreeNode node) => node == null;

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public bool IsRoot(SuffixTreeNode node) => node == _tree._root;

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public int GetEdgeSymbol(SuffixTreeNode node, int offset)
            => _tree.GetSymbolAt(node.Start + offset);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public int LengthOf(SuffixTreeNode node) => _tree.LengthOf(node);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public int GetNodeDepth(SuffixTreeNode node) => _tree.GetNodeDepth(node);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public int GetDepthFromRoot(SuffixTreeNode node) => node.DepthFromRoot;

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public SuffixTreeNode GetSuffixLink(SuffixTreeNode node)
            => node.SuffixLink ?? _tree._root;

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public bool TryGetChild(SuffixTreeNode node, int key, out SuffixTreeNode child)
        {
            bool found = node.TryGetChild(key, out var c);
            child = c!;
            return found && c != null;
        }

        public void CollectLeaves(SuffixTreeNode node, int depth, List<int> results)
            => _tree.CollectLeaves(node, depth, results);

        public int FindAnyLeafPosition(SuffixTreeNode node)
            => _tree.FindAnyLeafPosition(node);
    }
}
