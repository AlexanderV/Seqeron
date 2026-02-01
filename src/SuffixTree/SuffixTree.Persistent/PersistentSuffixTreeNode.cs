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

    public int Start
    {
        get => _storage.ReadInt32(_offset + PersistentConstants.OFFSET_START);
        set => _storage.WriteInt32(_offset + PersistentConstants.OFFSET_START, value);
    }

    public int End
    {
        get => _storage.ReadInt32(_offset + PersistentConstants.OFFSET_END);
        set => _storage.WriteInt32(_offset + PersistentConstants.OFFSET_END, value);
    }

    public long SuffixLink
    {
        get => _storage.ReadInt64(_offset + PersistentConstants.OFFSET_SUFFIX_LINK);
        set => _storage.WriteInt64(_offset + PersistentConstants.OFFSET_SUFFIX_LINK, value);
    }

    public int DepthFromRoot
    {
        get => _storage.ReadInt32(_offset + PersistentConstants.OFFSET_DEPTH);
        set => _storage.WriteInt32(_offset + PersistentConstants.OFFSET_DEPTH, value);
    }

    public int LeafCount
    {
        get => _storage.ReadInt32(_offset + PersistentConstants.OFFSET_LEAF_COUNT);
        set => _storage.WriteInt32(_offset + PersistentConstants.OFFSET_LEAF_COUNT, value);
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

    public bool TryGetChild(int key, out PersistentSuffixTreeNode child)
    {
        long currentOffset = ChildrenHead;
        while (currentOffset != PersistentConstants.NULL_OFFSET)
        {
            var entry = new PersistentChildEntry(_storage, currentOffset);
            if (entry.Key == key)
            {
                child = new PersistentSuffixTreeNode(_storage, entry.ChildNodeOffset);
                return true;
            }
            currentOffset = entry.NextEntryOffset;
        }
        child = Null(_storage);
        return false;
    }

    public void SetChild(int key, PersistentSuffixTreeNode child)
    {
        long currentOffset = ChildrenHead;
        long tailOffset = PersistentConstants.NULL_OFFSET;

        while (currentOffset != PersistentConstants.NULL_OFFSET)
        {
            var entry = new PersistentChildEntry(_storage, currentOffset);
            if (entry.Key == key)
            {
                entry.ChildNodeOffset = child.Offset;
                return;
            }
            tailOffset = currentOffset;
            currentOffset = entry.NextEntryOffset;
        }

        // Not found, create new entry
        long newEntryOffset = _storage.Allocate(PersistentConstants.CHILD_ENTRY_SIZE);
        var newEntry = new PersistentChildEntry(_storage, newEntryOffset);
        newEntry.Key = key;
        newEntry.ChildNodeOffset = child.Offset;
        newEntry.NextEntryOffset = PersistentConstants.NULL_OFFSET;

        if (tailOffset == PersistentConstants.NULL_OFFSET)
        {
            ChildrenHead = newEntryOffset;
        }
        else
        {
            var tailEntry = new PersistentChildEntry(_storage, tailOffset);
            tailEntry.NextEntryOffset = newEntryOffset;
        }

        ChildCount++;
    }

    public bool HasChildren => ChildrenHead != PersistentConstants.NULL_OFFSET;

    public static PersistentSuffixTreeNode Null(IStorageProvider storage) 
        => new(storage, PersistentConstants.NULL_OFFSET);
}
