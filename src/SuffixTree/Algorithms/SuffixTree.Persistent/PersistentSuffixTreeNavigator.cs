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
    private readonly HybridLayout _hybrid;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public PersistentSuffixTreeNavigator(IStorageProvider storage, long rootOffset, ITextSource textSource,
        HybridLayout hybrid)
    {
        _storage = storage;
        _rootOffset = rootOffset;
        _textSource = textSource;
        _hybrid = hybrid;
    }

    /// <summary>Returns the correct layout for a node at the given offset. Delegates to <see cref="HybridLayout"/>.</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private NodeLayout LayoutForOffset(long offset) => _hybrid.LayoutForOffset(offset);

    /// <summary>Resolves a jump-table entry if the offset falls in the jump zone. Delegates to <see cref="HybridLayout"/>.</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private long ResolveJump(long offset) => _hybrid.ResolveJump(offset);

    /// <summary>Creates a node with the correct zone layout. Delegates to <see cref="HybridLayout"/>.</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private PersistentSuffixTreeNode NodeAt(long offset) => _hybrid.NodeAt(offset);

    /// <summary>Reads child-array info, handling jumped/large-entry arrays. Delegates to <see cref="HybridLayout"/>.</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private (long ArrayBase, NodeLayout EntryLayout, int Count) ReadChildArrayInfo(PersistentSuffixTreeNode parent)
        => _hybrid.ReadChildArrayInfo(parent);

    public ITextSource Text
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get => _textSource;
    }

    public PersistentSuffixTreeNode Root
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get => NodeAt(_rootOffset);
    }

    public PersistentSuffixTreeNode NullNode
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get => PersistentSuffixTreeNode.Null(_storage, _hybrid.Layout);
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
        long raw = node.SuffixLink;
        if (raw == PersistentConstants.NULL_OFFSET)
            return NodeAt(_rootOffset);
        long resolved = ResolveJump(raw);
        return NodeAt(resolved);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public bool TryGetChild(PersistentSuffixTreeNode node, int key, out PersistentSuffixTreeNode child)
    {
        var (arrayBase, entryLayout, count) = ReadChildArrayInfo(node);
        if (count == 0) { child = PersistentSuffixTreeNode.Null(_storage, _hybrid.Layout); return false; }

        int lo = 0, hi = count - 1;
        int signedKey = key;

        while (lo <= hi)
        {
            int mid = lo + ((hi - lo) >> 1);
            long entryOffset = arrayBase + (long)mid * entryLayout.ChildEntrySize;
            int midKey = (int)_storage.ReadUInt32(entryOffset + NodeLayout.ChildOffsetKey);

            if (midKey == signedKey)
            {
                long childOffset = entryLayout.ReadOffset(_storage, entryOffset + NodeLayout.ChildOffsetNode);
                child = NodeAt(childOffset);
                return true;
            }
            if (midKey < signedKey) lo = mid + 1;
            else hi = mid - 1;
        }
        child = PersistentSuffixTreeNode.Null(_storage, _hybrid.Layout);
        return false;
    }

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
            var (arrayBase, entryLayout, childCount) = ReadChildArrayInfo(current);
            for (int ci = 0; ci < childCount; ci++)
            {
                long entryOffset = arrayBase + (long)ci * entryLayout.ChildEntrySize;
                long childOffset = entryLayout.ReadOffset(_storage, entryOffset + NodeLayout.ChildOffsetNode);
                stack.Push((NodeAt(childOffset), childDepth));
            }
        }
    }

    public int FindAnyLeafPosition(PersistentSuffixTreeNode node)
    {
        var current = node;
        while (!current.IsLeaf)
        {
            var (arrayBase, entryLayout, childCount) = ReadChildArrayInfo(current);
            PersistentSuffixTreeNode bestChild = PersistentSuffixTreeNode.Null(_storage, _hybrid.Layout);
            for (int ci = 0; ci < childCount; ci++)
            {
                long entryOffset = arrayBase + (long)ci * entryLayout.ChildEntrySize;
                int key = (int)_storage.ReadUInt32(entryOffset + NodeLayout.ChildOffsetKey);
                long childOffset = entryLayout.ReadOffset(_storage, entryOffset + NodeLayout.ChildOffsetNode);
                bestChild = NodeAt(childOffset);
                if ((uint)key != PersistentConstants.TERMINATOR_KEY)
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
