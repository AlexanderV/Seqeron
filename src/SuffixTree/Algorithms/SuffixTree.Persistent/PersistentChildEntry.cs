namespace SuffixTree.Persistent;

/// <summary>
/// A handle to a child entry in the persistent storage.
/// Represents an entry in the sorted contiguous children array of a node.
/// </summary>
public readonly struct PersistentChildEntry
{
    private readonly IStorageProvider _storage;
    private readonly long _offset;

    public PersistentChildEntry(IStorageProvider storage, long offset)
    {
        _storage = storage;
        _offset = offset;
    }

    public long Offset => _offset;
    public bool IsNull => _offset == PersistentConstants.NULL_OFFSET;

    public uint Key
    {
        get => _storage.ReadUInt32(_offset + PersistentConstants.CHILD_OFFSET_KEY);
        set => _storage.WriteUInt32(_offset + PersistentConstants.CHILD_OFFSET_KEY, value);
    }

    public long ChildNodeOffset
    {
        get => _storage.ReadInt64(_offset + PersistentConstants.CHILD_OFFSET_NODE);
        set => _storage.WriteInt64(_offset + PersistentConstants.CHILD_OFFSET_NODE, value);
    }

    public static PersistentChildEntry Null(IStorageProvider storage)
        => new(storage, PersistentConstants.NULL_OFFSET);
}
