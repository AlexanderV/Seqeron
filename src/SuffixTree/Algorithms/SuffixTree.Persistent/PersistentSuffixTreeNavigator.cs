using System.Runtime.CompilerServices;

namespace SuffixTree.Persistent;

/// <summary>
/// Struct navigator for the persistent (memory-mapped) suffix tree.
/// Implements <see cref="ISuffixTreeNavigator{TNode}"/> to enable
/// shared algorithm dispatch with zero overhead (JIT specialization).
/// </summary>
internal unsafe struct PersistentSuffixTreeNavigator : ISuffixTreeNavigator<PersistentSuffixTreeNode>
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
        if (_hybrid.TryGetChildFast(node.Offset, (uint)key, out long childOffset))
        {
            child = NodeAt(childOffset);
            return true;
        }
        child = PersistentSuffixTreeNode.Null(_storage, _hybrid.Layout);
        return false;
    }

    public void CollectLeaves(PersistentSuffixTreeNode node, int depth, List<int> results)
    {
        byte* ptr = _hybrid.Ptr;
        if (ptr != null)
        {
            // Fast path: all reads via direct pointer — no interface dispatch or bounds checks
            CollectLeavesUnsafe(ptr, node.Offset, depth, results);
            return;
        }

        // Safe fallback for non-MMF storage
        CollectLeavesSequential(node, depth, results);
    }

    /// <summary>
    /// High-performance leaf collection using direct MMF pointer reads.
    /// Pure DFS — eliminates IStorageProvider interface dispatch, bounds checks,
    /// and node struct creation. No extra allocations beyond the explicit stack.
    /// </summary>
    private void CollectLeavesUnsafe(byte* ptr, long nodeOffset, int depth, List<int> results)
    {
        int textLen = _textSource.Length;
        int textLenPlus1 = textLen + 1;
        long transOff = _hybrid.TransitionOffset;

        var stack = new Stack<(long Offset, int Depth)>(256);
        stack.Push((nodeOffset, depth));

        while (stack.Count > 0)
        {
            var (off, d) = stack.Pop();

            uint nodeStart = *(uint*)(ptr + off);
            uint nodeEnd = *(uint*)(ptr + off + 4);

            if (nodeEnd == PersistentConstants.BOUNDLESS)
            {
                int edgeLen = textLenPlus1 - (int)nodeStart;
                int pos = textLenPlus1 - (d + edgeLen);
                if (pos < textLen)
                    results.Add(pos);
                continue;
            }

            int childDepth = d + (int)(nodeEnd - nodeStart);

            bool isLarge = transOff >= 0 && off >= transOff;
            int rawCount = *(int*)(ptr + off + (isLarge ? 28 : 20));
            bool isJumped = (rawCount & unchecked((int)0x80000000)) != 0;
            int count = isJumped ? (rawCount & 0x7FFFFFFF) : rawCount;

            long head;
            if (isLarge)
                head = *(long*)(ptr + off + 20);
            else
            {
                uint rawHead = *(uint*)(ptr + off + 16);
                head = rawHead == uint.MaxValue ? PersistentConstants.NULL_OFFSET : (long)rawHead;
            }

            long arrayBase;
            int entrySize;
            bool entryIs64;
            if (isJumped)
            {
                arrayBase = *(long*)(ptr + head);
                entrySize = 12;
                entryIs64 = true;
            }
            else
            {
                arrayBase = head;
                entrySize = isLarge ? 12 : 8;
                entryIs64 = isLarge;
            }

            for (int ci = count - 1; ci >= 0; ci--)
            {
                byte* entry = ptr + arrayBase + (long)ci * entrySize;
                long childOff = entryIs64 ? *(long*)(entry + 4) : (long)*(uint*)(entry + 4);
                stack.Push((childOff, childDepth));
            }
        }
    }

    private void CollectLeavesSequential(PersistentSuffixTreeNode node, int depth, List<int> results)
    {
        var stack = new Stack<(long NodeOffset, int Depth)>();
        stack.Push((node.Offset, depth));

        while (stack.Count > 0)
        {
            var (currentOffset, currentDepth) = stack.Pop();
            var current = NodeAt(currentOffset);

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
            for (int ci = childCount - 1; ci >= 0; ci--)
            {
                long entryOffset = arrayBase + (long)ci * entryLayout.ChildEntrySize;
                long childOffset = entryLayout.ReadOffset(_storage, entryOffset + NodeLayout.ChildOffsetNode);
                stack.Push((childOffset, childDepth));
            }
        }
    }

    public int FindAnyLeafPosition(PersistentSuffixTreeNode node, int depthFromRoot)
    {
        var current = node;
        int depth = depthFromRoot;
        while (!current.IsLeaf)
        {
            depth += LengthOf(current);
            var (arrayBase, entryLayout, childCount) = ReadChildArrayInfo(current);
            PersistentSuffixTreeNode bestChild = PersistentSuffixTreeNode.Null(_storage, _hybrid.Layout);
            for (int ci = 0; ci < childCount; ci++)
            {
                long entryOffset = arrayBase + (long)ci * entryLayout.ChildEntrySize;
                int key = (int)_hybrid.ReadUInt32Fast(entryOffset + NodeLayout.ChildOffsetKey);
                long childOffset = entryLayout.OffsetIs64Bit
                    ? _hybrid.ReadInt64Fast(entryOffset + NodeLayout.ChildOffsetNode)
                    : (long)_hybrid.ReadUInt32Fast(entryOffset + NodeLayout.ChildOffsetNode);
                bestChild = NodeAt(childOffset);
                if ((uint)key != PersistentConstants.TERMINATOR_KEY)
                    break;
            }
            if (bestChild.IsNull) return -1;
            current = bestChild;
        }

        depth += LengthOf(current);
        int pos = (_textSource.Length + 1) - depth;
        return (pos >= 0 && pos < _textSource.Length) ? pos : -1;
    }
}
