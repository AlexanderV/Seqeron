using System.Diagnostics;
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
internal readonly unsafe struct HybridLayout
{
    private readonly IStorageProvider _storage;
    private readonly NodeLayout _layout;
    private readonly long _transitionOffset;
    private readonly long _jumpTableStart;
    private readonly long _jumpTableEnd;
    private readonly byte* _ptr; // non-null when backed by MappedFileStorageProvider

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public HybridLayout(IStorageProvider storage, NodeLayout layout,
        long transitionOffset, long jumpTableStart, long jumpTableEnd)
    {
        _storage = storage;
        _layout = layout;
        _transitionOffset = transitionOffset;
        _jumpTableStart = jumpTableStart;
        _jumpTableEnd = jumpTableEnd;
        _ptr = storage is MappedFileStorageProvider mmf ? mmf.RawPointer : null;
    }

    /// <summary>The base layout (Compact for most trees).</summary>
    public NodeLayout Layout
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get => _layout;
    }

    /// <summary>Transition offset (compact/large boundary), or -1 for single-format.</summary>
    public long TransitionOffset
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get => _transitionOffset;
    }

    /// <summary>Start of contiguous jump table, or -1 for single-format.</summary>
    public long JumpTableStart
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get => _jumpTableStart;
    }

    /// <summary>End of jump table, or -1 for single-format.</summary>
    public long JumpTableEnd
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get => _jumpTableEnd;
    }

    /// <summary>Whether this layout uses the hybrid v5 format with dual zones.</summary>
    public bool IsHybrid
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get => _transitionOffset >= 0;
    }

    /// <summary>Raw MMF pointer for unchecked fast-path reads (null if not MMF-backed).</summary>
    internal byte* Ptr
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get => _ptr;
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
            if (_ptr != null) return *(long*)(_ptr + offset);
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
            long realArrayOffset = _ptr != null ? *(long*)(_ptr + head) : _storage.ReadInt64(head);
            return (realArrayOffset, NodeLayout.Large, count);
        }

        return (head, LayoutForOffset(parent.Offset), count);
    }

    // ──────────────── Unchecked fast-path for query hot loops ────────────────

    /// <summary>
    /// Reads child-array info directly from the MMF pointer, bypassing the node struct
    /// and all IStorageProvider bounds checks. Falls back to safe path when not MMF-backed.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    internal (long ArrayBase, NodeLayout EntryLayout, int Count) ReadChildArrayInfoFast(long nodeOffset)
    {
        if (_ptr == null) return ReadChildArrayInfo(NodeAt(nodeOffset));

        NodeLayout layout = LayoutForOffset(nodeOffset);
        int rawCount = *(int*)(_ptr + nodeOffset + layout.OffsetChildCount);

        long head;
        if (layout.OffsetIs64Bit)
            head = *(long*)(_ptr + nodeOffset + layout.OffsetChildrenHead);
        else
        {
            uint rawHead = *(uint*)(_ptr + nodeOffset + layout.OffsetChildrenHead);
            head = rawHead == uint.MaxValue ? PersistentConstants.NULL_OFFSET : (long)rawHead;
        }

        bool isJumped = (rawCount & unchecked((int)0x80000000)) != 0;
        int count = isJumped ? (rawCount & 0x7FFFFFFF) : rawCount;

        if (isJumped)
        {
            long realArrayOffset = *(long*)(_ptr + head);
            return (realArrayOffset, NodeLayout.Large, count);
        }

        return (head, layout, count);
    }

    /// <summary>Read uint32 directly from MMF or via storage.</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    internal uint ReadUInt32Fast(long offset)
        => _ptr != null ? *(uint*)(_ptr + offset) : _storage.ReadUInt32(offset);

    /// <summary>Read int64 directly from MMF or via storage.</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    internal long ReadInt64Fast(long offset)
        => _ptr != null ? *(long*)(_ptr + offset) : _storage.ReadInt64(offset);

    /// <summary>
    /// Binary search a sorted child array and return the child node offset if found.
    /// Uses direct pointer access when available.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    internal bool TryGetChildFast(long nodeOffset, uint key, out long childNodeOffset)
    {
        var (arrayBase, entryLayout, count) = ReadChildArrayInfoFast(nodeOffset);
        if (count == 0) { childNodeOffset = PersistentConstants.NULL_OFFSET; return false; }

        int lo = 0, hi = count - 1;
        int signedKey = (int)key;
        int entrySize = entryLayout.ChildEntrySize;
        bool is64 = entryLayout.OffsetIs64Bit;

        if (_ptr != null)
        {
            byte* p = _ptr;
            while (lo <= hi)
            {
                int mid = lo + ((hi - lo) >> 1);
                byte* entry = p + arrayBase + (long)mid * entrySize;
                int midKey = (int)*(uint*)entry;
#if DEBUG
                if (mid > 0)
                {
                    int prevKey = (int)*(uint*)(p + arrayBase + (long)(mid - 1) * entrySize);
                    Debug.Assert(prevKey <= midKey, $"Child array keys not sorted: prev={prevKey}, mid={midKey}");
                }
#endif

                if (midKey == signedKey)
                {
                    childNodeOffset = is64 ? *(long*)(entry + 4) : (long)*(uint*)(entry + 4);
                    return true;
                }
                if (midKey < signedKey) lo = mid + 1;
                else hi = mid - 1;
            }
        }
        else
        {
            while (lo <= hi)
            {
                int mid = lo + ((hi - lo) >> 1);
                long entryOffset = arrayBase + (long)mid * entrySize;
                int midKey = (int)_storage.ReadUInt32(entryOffset + NodeLayout.ChildOffsetKey);
#if DEBUG
                if (mid > 0)
                {
                    int prevKey = (int)_storage.ReadUInt32(arrayBase + (long)(mid - 1) * entrySize + NodeLayout.ChildOffsetKey);
                    Debug.Assert(prevKey <= midKey, $"Child array keys not sorted (safe path): prev={prevKey}, mid={midKey}");
                }
#endif

                if (midKey == signedKey)
                {
                    childNodeOffset = entryLayout.ReadOffset(_storage, entryOffset + NodeLayout.ChildOffsetNode);
                    return true;
                }
                if (midKey < signedKey) lo = mid + 1;
                else hi = mid - 1;
            }
        }
        childNodeOffset = PersistentConstants.NULL_OFFSET;
        return false;
    }
}
