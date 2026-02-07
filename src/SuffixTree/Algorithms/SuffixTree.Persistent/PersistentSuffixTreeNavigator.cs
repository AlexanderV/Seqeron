using System.Collections.Generic;
using System.Runtime.CompilerServices;

namespace SuffixTree.Persistent;

/// <summary>
/// Struct navigator for the persistent (memory-mapped) suffix tree.
/// Implements <see cref="ISuffixTreeNavigator{TNode}"/> to enable
/// shared algorithm dispatch with zero overhead (JIT specialization).
/// </summary>
internal struct PersistentSuffixTreeNavigator : ISuffixTreeNavigator<PersistentSuffixTreeNode>
{
    private readonly IStorageProvider _storage;
    private readonly long _rootOffset;
    private readonly ITextSource _textSource;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public PersistentSuffixTreeNavigator(IStorageProvider storage, long rootOffset, ITextSource textSource)
    {
        _storage = storage;
        _rootOffset = rootOffset;
        _textSource = textSource;
    }

    public ITextSource Text
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get => _textSource;
    }

    public PersistentSuffixTreeNode Root
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get => new PersistentSuffixTreeNode(_storage, _rootOffset);
    }

    public PersistentSuffixTreeNode NullNode
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get => PersistentSuffixTreeNode.Null(_storage);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public bool IsNull(PersistentSuffixTreeNode node) => node.IsNull;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public bool IsRoot(PersistentSuffixTreeNode node) => node.Offset == _rootOffset;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public int GetEdgeSymbol(PersistentSuffixTreeNode node, int offset)
    {
        int index = (int)(node.Start + (uint)offset);
        return index >= _textSource.Length ? -1 : _textSource[index];
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public int LengthOf(PersistentSuffixTreeNode node)
    {
        uint end = node.End == PersistentConstants.BOUNDLESS
            ? (uint)(_textSource.Length + 1)
            : node.End;
        return (int)(end - node.Start);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public int GetNodeDepth(PersistentSuffixTreeNode node)
    {
        if (node.Offset == _rootOffset) return 0;
        return (int)node.DepthFromRoot + LengthOf(node);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public int GetDepthFromRoot(PersistentSuffixTreeNode node) => (int)node.DepthFromRoot;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public PersistentSuffixTreeNode GetSuffixLink(PersistentSuffixTreeNode node)
    {
        long suffixLink = node.SuffixLink;
        return suffixLink != PersistentConstants.NULL_OFFSET
            ? new PersistentSuffixTreeNode(_storage, suffixLink)
            : new PersistentSuffixTreeNode(_storage, _rootOffset);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public bool TryGetChild(PersistentSuffixTreeNode node, int key, out PersistentSuffixTreeNode child)
        => node.TryGetChild((uint)key, out child);

    public void CollectLeaves(PersistentSuffixTreeNode node, int depth, List<int> results)
    {
        var stack = new Stack<(PersistentSuffixTreeNode Node, int Depth)>();
        stack.Push((node, depth));

        while (stack.Count > 0)
        {
            var (current, currentDepth) = stack.Pop();

            if (current.IsLeaf)
            {
                int suffixLength = currentDepth + LengthOf(current);
                int startPosition = (_textSource.Length + 1) - suffixLength;
                if (startPosition < _textSource.Length)
                    results.Add(startPosition);
                continue;
            }

            int childDepth = currentDepth + LengthOf(current);
            int childCount = current.ChildCount;
            long arrayBase = current.ChildrenHead;
            for (int ci = 0; ci < childCount; ci++)
            {
                var entry = new PersistentChildEntry(_storage, arrayBase + (long)ci * PersistentConstants.CHILD_ENTRY_SIZE);
                stack.Push((new PersistentSuffixTreeNode(_storage, entry.ChildNodeOffset), childDepth));
            }
        }
    }

    public int FindAnyLeafPosition(PersistentSuffixTreeNode node)
    {
        var current = node;
        while (!current.IsLeaf)
        {
            int childCount = current.ChildCount;
            long arrayBase = current.ChildrenHead;
            PersistentSuffixTreeNode bestChild = PersistentSuffixTreeNode.Null(_storage);
            for (int ci = 0; ci < childCount; ci++)
            {
                var entry = new PersistentChildEntry(_storage, arrayBase + (long)ci * PersistentConstants.CHILD_ENTRY_SIZE);
                bestChild = new PersistentSuffixTreeNode(_storage, entry.ChildNodeOffset);
                if (entry.Key != PersistentConstants.TERMINATOR_KEY)
                    break;
            }
            if (bestChild.IsNull) return -1;
            current = bestChild;
        }

        int leafDepth = GetNodeDepth(current);
        int pos = (_textSource.Length + 1) - leafDepth;
        return (pos >= 0 && pos < _textSource.Length) ? pos : -1;
    }
}
