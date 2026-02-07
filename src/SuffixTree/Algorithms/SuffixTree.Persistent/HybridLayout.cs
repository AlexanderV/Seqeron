using System.Runtime.CompilerServices;

namespace SuffixTree.Persistent;

/// <summary>
/// Encapsulates hybrid v5 layout resolution logic: zone detection, jump-table
/// dereferencing, and child-array info reading.
/// <para>
/// Shared by <see cref="PersistentSuffixTree"/> and
/// <see cref="PersistentSuffixTreeNavigator"/> to eliminate code duplication.
/// </para>
/// </summary>
internal readonly struct HybridLayout
{
    private readonly IStorageProvider _storage;
    private readonly NodeLayout _layout;
    private readonly long _transitionOffset;
    private readonly long _jumpTableStart;
    private readonly long _jumpTableEnd;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public HybridLayout(IStorageProvider storage, NodeLayout layout,
        long transitionOffset, long jumpTableStart, long jumpTableEnd)
    {
        _storage = storage;
        _layout = layout;
        _transitionOffset = transitionOffset;
        _jumpTableStart = jumpTableStart;
        _jumpTableEnd = jumpTableEnd;
    }

    /// <summary>The base layout (Compact for most trees).</summary>
    public NodeLayout Layout
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get => _layout;
    }

    /// <summary>Whether this layout uses the hybrid v5 format with dual zones.</summary>
    public bool IsHybrid
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get => _transitionOffset >= 0;
    }

    /// <summary>
    /// Returns the correct <see cref="NodeLayout"/> for a node at the given offset.
    /// For pure Compact/Large trees, always returns the base layout.
    /// For hybrid v5 trees, returns Compact for offsets below the transition and Large above.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public NodeLayout LayoutForOffset(long offset)
    {
        if (_transitionOffset < 0)
            return _layout;
        return offset < _transitionOffset ? NodeLayout.Compact : NodeLayout.Large;
    }

    /// <summary>
    /// Resolves an offset that might be a jump-table entry.
    /// If the offset falls within [jumpTableStart, jumpTableEnd), reads the int64 target
    /// from the jump entry and returns it. Otherwise returns the offset unchanged.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public long ResolveJump(long offset)
    {
        if (_jumpTableStart >= 0
            && offset >= _jumpTableStart
            && offset < _jumpTableEnd)
        {
            return _storage.ReadInt64(offset);
        }
        return offset;
    }

    /// <summary>
    /// Creates a <see cref="PersistentSuffixTreeNode"/> with the correct layout for its zone.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public PersistentSuffixTreeNode NodeAt(long offset)
        => new PersistentSuffixTreeNode(_storage, offset, LayoutForOffset(offset));

    /// <summary>
    /// Reads child information from a parent node, handling hybrid jump entries.
    /// Returns the real child-array offset, child-entry layout, and actual child count.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public (long ArrayBase, NodeLayout EntryLayout, int Count) ReadChildArrayInfo(PersistentSuffixTreeNode parent)
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
}
