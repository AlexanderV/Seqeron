namespace SuffixTree.Persistent;

/// <summary>
/// A handle to a node in the persistent storage.
/// Represents a fixed-size block in the storage provider.
/// </summary>
public readonly struct PersistentSuffixTreeNode
{
    private readonly IStorageProvider _storage;
    private readonly long _offset;

    public PersistentSuffixTreeNode(IStorageProvider storage, long offset)
    {
        _storage = storage;
        _offset = offset;
    }

    public long Offset => _offset;
    public bool IsNull => _offset == PersistentConstants.NULL_OFFSET;

    public uint Start
    {
        get => _storage.ReadUInt32(_offset + PersistentConstants.OFFSET_START);
        set => _storage.WriteUInt32(_offset + PersistentConstants.OFFSET_START, value);
    }

    public uint End
    {
        get => _storage.ReadUInt32(_offset + PersistentConstants.OFFSET_END);
        set => _storage.WriteUInt32(_offset + PersistentConstants.OFFSET_END, value);
    }

    public long SuffixLink
    {
        get => _storage.ReadInt64(_offset + PersistentConstants.OFFSET_SUFFIX_LINK);
        set => _storage.WriteInt64(_offset + PersistentConstants.OFFSET_SUFFIX_LINK, value);
    }

    public uint DepthFromRoot
    {
        get => _storage.ReadUInt32(_offset + PersistentConstants.OFFSET_DEPTH);
        set => _storage.WriteUInt32(_offset + PersistentConstants.OFFSET_DEPTH, value);
    }

    public uint LeafCount
    {
        get => _storage.ReadUInt32(_offset + PersistentConstants.OFFSET_LEAF_COUNT);
        set => _storage.WriteUInt32(_offset + PersistentConstants.OFFSET_LEAF_COUNT, value);
    }

    public long ChildrenHead
    {
        get => _storage.ReadInt64(_offset + PersistentConstants.OFFSET_CHILDREN_HEAD);
        set => _storage.WriteInt64(_offset + PersistentConstants.OFFSET_CHILDREN_HEAD, value);
    }

    public int ChildCount
    {
        get => _storage.ReadInt32(_offset + PersistentConstants.OFFSET_CHILD_COUNT);
        set => _storage.WriteInt32(_offset + PersistentConstants.OFFSET_CHILD_COUNT, value);
    }

    public bool IsLeaf => End == PersistentConstants.BOUNDLESS;

    public bool TryGetChild(uint key, out PersistentSuffixTreeNode child)
    {
        int count = ChildCount;
        if (count == 0) { child = Null(_storage); return false; }

        long arrayBase = ChildrenHead;
        int lo = 0, hi = count - 1;
        int signedKey = (int)key;

        while (lo <= hi)
        {
            int mid = lo + ((hi - lo) >> 1);
            long entryOffset = arrayBase + (long)mid * PersistentConstants.CHILD_ENTRY_SIZE;
            int midKey = (int)_storage.ReadUInt32(entryOffset + PersistentConstants.CHILD_OFFSET_KEY);

            if (midKey == signedKey)
            {
                long childOffset = _storage.ReadInt64(entryOffset + PersistentConstants.CHILD_OFFSET_NODE);
                child = new PersistentSuffixTreeNode(_storage, childOffset);
                return true;
            }
            if (midKey < signedKey) lo = mid + 1;
            else hi = mid - 1;
        }
        child = Null(_storage);
        return false;
    }

    public bool HasChildren => ChildCount > 0;

    public static PersistentSuffixTreeNode Null(IStorageProvider storage)
        => new(storage, PersistentConstants.NULL_OFFSET);
}
