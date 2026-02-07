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
    private readonly NodeLayout _layout;
    private readonly long _transitionOffset;
    private readonly long _jumpTableStart;
    private readonly long _jumpTableEnd;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public PersistentSuffixTreeNavigator(IStorageProvider storage, long rootOffset, ITextSource textSource,
        NodeLayout layout, long transitionOffset = -1, long jumpTableStart = -1, long jumpTableEnd = -1)
    {
        _storage = storage;
        _rootOffset = rootOffset;
        _textSource = textSource;
        _layout = layout;
        _transitionOffset = transitionOffset;
        _jumpTableStart = jumpTableStart;
        _jumpTableEnd = jumpTableEnd;
    }

    /// <summary>Returns the correct layout for a node at the given offset.</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private NodeLayout LayoutForOffset(long offset)
    {
        if (_transitionOffset < 0) return _layout;
        return offset < _transitionOffset ? NodeLayout.Compact : NodeLayout.Large;
    }

    /// <summary>Resolves a jump-table entry if the offset falls in the jump zone.</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private long ResolveJump(long offset)
    {
        if (_jumpTableStart >= 0 && offset >= _jumpTableStart && offset < _jumpTableEnd)
            return _storage.ReadInt64(offset);
        return offset;
    }

    /// <summary>Creates a node with the correct zone layout.</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private PersistentSuffixTreeNode NodeAt(long offset)
        => new PersistentSuffixTreeNode(_storage, offset, LayoutForOffset(offset));

    /// <summary>Reads child-array info, handling jumped/large-entry arrays.</summary>
    private (long ArrayBase, NodeLayout EntryLayout, int Count) ReadChildArrayInfo(PersistentSuffixTreeNode parent)
    {
        int rawCount = parent.ChildCount;
        long head = parent.ChildrenHead;
        bool isJumped = (rawCount & unchecked((int)0x80000000)) != 0;
        int count = isJumped ? (rawCount & 0x7FFFFFFF) : rawCount;

        if (isJumped)
        {
            long realArrayOffset = _storage.ReadInt64(head);
            return (realArrayOffset, NodeLayout.Large, count);
        }
        return (head, LayoutForOffset(parent.Offset), count);
    }

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
        get => PersistentSuffixTreeNode.Null(_storage, _layout);
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
        if (count == 0) { child = PersistentSuffixTreeNode.Null(_storage, _layout); return false; }

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
        child = PersistentSuffixTreeNode.Null(_storage, _layout);
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
            PersistentSuffixTreeNode bestChild = PersistentSuffixTreeNode.Null(_storage, _layout);
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
