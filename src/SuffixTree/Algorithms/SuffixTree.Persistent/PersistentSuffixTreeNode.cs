namespace SuffixTree.Persistent;

/// <summary>
/// A handle to a node in the persistent storage.
/// Represents a fixed-size block whose layout is described by <see cref="NodeLayout"/>.
/// </summary>
public readonly struct PersistentSuffixTreeNode
{
    private readonly IStorageProvider _storage;
    private readonly long _offset;
    private readonly NodeLayout _layout;

    public PersistentSuffixTreeNode(IStorageProvider storage, long offset, NodeLayout layout)
    {
        _storage = storage;
        _offset = offset;
        _layout = layout;
    }

    public long Offset => _offset;
    public bool IsNull => _offset == PersistentConstants.NULL_OFFSET;
    internal NodeLayout Layout => _layout;

    public uint Start
    {
        get => _storage.ReadUInt32(_offset + PersistentConstants.OFFSET_START);
        internal set => _storage.WriteUInt32(_offset + PersistentConstants.OFFSET_START, value);
    }

    public uint End
    {
        get => _storage.ReadUInt32(_offset + PersistentConstants.OFFSET_END);
        internal set => _storage.WriteUInt32(_offset + PersistentConstants.OFFSET_END, value);
    }

    internal long SuffixLink
    {
        get => _layout.ReadOffset(_storage, _offset + _layout.OffsetSuffixLink);
        set => _layout.WriteOffset(_storage, _offset + _layout.OffsetSuffixLink, value);
    }

    public uint DepthFromRoot
    {
        get => _storage.ReadUInt32(_offset + _layout.OffsetDepth);
        internal set => _storage.WriteUInt32(_offset + _layout.OffsetDepth, value);
    }

    public uint LeafCount
    {
        get => _storage.ReadUInt32(_offset + _layout.OffsetLeafCount);
        internal set => _storage.WriteUInt32(_offset + _layout.OffsetLeafCount, value);
    }

    internal long ChildrenHead
    {
        get => _layout.ReadOffset(_storage, _offset + _layout.OffsetChildrenHead);
        set => _layout.WriteOffset(_storage, _offset + _layout.OffsetChildrenHead, value);
    }

    internal int ChildCount
    {
        get => _storage.ReadInt32(_offset + _layout.OffsetChildCount);
        set => _storage.WriteInt32(_offset + _layout.OffsetChildCount, value);
    }

    public bool IsLeaf => End == PersistentConstants.BOUNDLESS;

    internal bool TryGetChild(uint key, out PersistentSuffixTreeNode child)
    {
        int rawCount = ChildCount;

        // Hybrid v5 trees store a jumped-flag in the high bit of ChildCount.
        // A jumped child array requires hybrid context (transition offset,
        // jump table bounds) that this struct does not carry.  Throw early
        // instead of silently returning wrong results.
        if ((rawCount & unchecked((int)0x80000000)) != 0)
            throw new InvalidOperationException(
                "Cannot resolve children of a hybrid-jumped node through " +
                "PersistentSuffixTreeNode.TryGetChild. Use ISuffixTree methods instead.");

        int count = rawCount;
        if (count == 0) { child = Null(_storage, _layout); return false; }

        long arrayBase = ChildrenHead;
        int lo = 0, hi = count - 1;
        int signedKey = (int)key;

        while (lo <= hi)
        {
            int mid = lo + ((hi - lo) >> 1);
            long entryOffset = arrayBase + (long)mid * _layout.ChildEntrySize;
            int midKey = (int)_storage.ReadUInt32(entryOffset + NodeLayout.ChildOffsetKey);

            if (midKey == signedKey)
            {
                long childOffset = _layout.ReadOffset(_storage, entryOffset + NodeLayout.ChildOffsetNode);
                child = new PersistentSuffixTreeNode(_storage, childOffset, _layout);
                return true;
            }
            if (midKey < signedKey) lo = mid + 1;
            else hi = mid - 1;
        }
        child = Null(_storage, _layout);
        return false;
    }

    public bool HasChildren => (ChildCount & 0x7FFFFFFF) > 0;

    public static PersistentSuffixTreeNode Null(IStorageProvider storage, NodeLayout layout)
        => new(storage, PersistentConstants.NULL_OFFSET, layout);
}
