namespace SuffixTree.Persistent;

/// <summary>
/// A handle to a child entry in the persistent storage.
/// Represents an entry in the sorted contiguous children array of a node.
/// </summary>
public readonly struct PersistentChildEntry
{
    private readonly IStorageProvider _storage;
    private readonly long _offset;
    private readonly NodeLayout _layout;

    public PersistentChildEntry(IStorageProvider storage, long offset, NodeLayout layout)
    {
        _storage = storage;
        _offset = offset;
        _layout = layout;
    }

    public long Offset => _offset;
    public bool IsNull => _offset == PersistentConstants.NULL_OFFSET;

    public uint Key
    {
        get => _storage.ReadUInt32(_offset + NodeLayout.ChildOffsetKey);
        set => _storage.WriteUInt32(_offset + NodeLayout.ChildOffsetKey, value);
    }

    public long ChildNodeOffset
    {
        get => _layout.ReadOffset(_storage, _offset + NodeLayout.ChildOffsetNode);
        set => _layout.WriteOffset(_storage, _offset + NodeLayout.ChildOffsetNode, value);
    }

    public static PersistentChildEntry Null(IStorageProvider storage, NodeLayout layout)
        => new(storage, PersistentConstants.NULL_OFFSET, layout);
}
