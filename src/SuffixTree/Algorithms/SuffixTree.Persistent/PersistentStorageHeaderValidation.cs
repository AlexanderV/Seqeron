namespace SuffixTree.Persistent;

internal readonly struct PersistentStorageHeader
{
    public PersistentStorageHeader(
        long magic,
        int version,
        long recordedSize,
        long rootOffset,
        long textOffset,
        int textLength,
        int flags,
        long transitionOffset,
        long jumpTableStart,
        long jumpTableEnd,
        long deepestNodeRawOffset,
        int lrsDepth,
        int baseNodeSize)
    {
        Magic = magic;
        Version = version;
        RecordedSize = recordedSize;
        RootOffset = rootOffset;
        TextOffset = textOffset;
        TextLength = textLength;
        Flags = flags;
        TransitionOffset = transitionOffset;
        JumpTableStart = jumpTableStart;
        JumpTableEnd = jumpTableEnd;
        DeepestNodeRawOffset = deepestNodeRawOffset;
        LrsDepth = lrsDepth;
        BaseNodeSize = baseNodeSize;
    }

    public long Magic { get; }
    public int Version { get; }
    public long RecordedSize { get; }
    public long RootOffset { get; }
    public long TextOffset { get; }
    public int TextLength { get; }
    public int Flags { get; }
    public long TransitionOffset { get; }
    public long JumpTableStart { get; }
    public long JumpTableEnd { get; }
    public long DeepestNodeRawOffset { get; }
    public int LrsDepth { get; }
    public int BaseNodeSize { get; }
}

internal readonly struct PersistentStorageMetadata
{
    public PersistentStorageMetadata(
        NodeLayout layout,
        long rootOffset,
        long transitionOffset,
        long jumpTableStart,
        long jumpTableEnd,
        long deepestNodeOffset,
        int lrsDepth)
    {
        Layout = layout;
        RootOffset = rootOffset;
        TransitionOffset = transitionOffset;
        JumpTableStart = jumpTableStart;
        JumpTableEnd = jumpTableEnd;
        DeepestNodeOffset = deepestNodeOffset;
        LrsDepth = lrsDepth;
    }

    public NodeLayout Layout { get; }
    public long RootOffset { get; }
    public long TransitionOffset { get; }
    public long JumpTableStart { get; }
    public long JumpTableEnd { get; }
    public long DeepestNodeOffset { get; }
    public int LrsDepth { get; }
}

internal static class PersistentStorageHeaderReader
{
    public static PersistentStorageHeader Read(IStorageProvider storage)
    {
        if (storage.Size < PersistentConstants.HEADER_SIZE_V6)
        {
            throw new InvalidOperationException(
                $"Invalid storage format: storage size {storage.Size} is smaller than v6 header size {PersistentConstants.HEADER_SIZE_V6}.");
        }

        try
        {
            return new PersistentStorageHeader(
                magic: storage.ReadInt64(PersistentConstants.HEADER_OFFSET_MAGIC),
                version: storage.ReadInt32(PersistentConstants.HEADER_OFFSET_VERSION),
                recordedSize: storage.ReadInt64(PersistentConstants.HEADER_OFFSET_SIZE),
                rootOffset: storage.ReadInt64(PersistentConstants.HEADER_OFFSET_ROOT),
                textOffset: storage.ReadInt64(PersistentConstants.HEADER_OFFSET_TEXT_OFF),
                textLength: storage.ReadInt32(PersistentConstants.HEADER_OFFSET_TEXT_LEN),
                flags: storage.ReadInt32(PersistentConstants.HEADER_OFFSET_FLAGS),
                transitionOffset: storage.ReadInt64(PersistentConstants.HEADER_OFFSET_TRANSITION),
                jumpTableStart: storage.ReadInt64(PersistentConstants.HEADER_OFFSET_JUMP_START),
                jumpTableEnd: storage.ReadInt64(PersistentConstants.HEADER_OFFSET_JUMP_END),
                deepestNodeRawOffset: storage.ReadInt64(PersistentConstants.HEADER_OFFSET_DEEPEST_NODE),
                lrsDepth: storage.ReadInt32(PersistentConstants.HEADER_OFFSET_LRS_DEPTH),
                baseNodeSize: storage.ReadInt32(PersistentConstants.HEADER_OFFSET_BASE_NODE_SIZE));
        }
        catch (ArgumentOutOfRangeException ex)
        {
            throw new InvalidOperationException("Invalid storage format: header is truncated.", ex);
        }
    }
}

internal static class PersistentStorageHeaderValidator
{
    public static PersistentStorageMetadata Validate(PersistentStorageHeader header, long storageSize)
    {
        if (header.Magic != PersistentConstants.MAGIC_NUMBER)
            throw new InvalidOperationException("Invalid storage format: Magic number mismatch.");

        NodeLayout layout = NodeLayout.ForVersion(header.Version);
        if (header.BaseNodeSize == NodeLayout.Large.NodeSize)
            layout = NodeLayout.Large;

        int headerSize = PersistentConstants.HEADER_SIZE_V6;

        if (header.RecordedSize > 0 && header.RecordedSize != storageSize)
            throw new InvalidOperationException(
                $"Invalid storage format: header size {header.RecordedSize} does not match actual storage size {storageSize}. File may be truncated or corrupted.");

        if (header.RootOffset < headerSize || header.RootOffset >= storageSize)
            throw new InvalidOperationException(
                $"Invalid storage format: root offset {header.RootOffset} is outside valid range [{headerSize}, {storageSize}).");

        if (header.TransitionOffset >= 0)
        {
            if (header.TransitionOffset < headerSize || header.TransitionOffset > storageSize)
                throw new InvalidOperationException(
                    $"Invalid storage format: transition offset {header.TransitionOffset} is outside valid range [{headerSize}, {storageSize}].");

            if (header.JumpTableEnd < header.JumpTableStart)
                throw new InvalidOperationException(
                    $"Invalid storage format: jump table end ({header.JumpTableEnd}) < start ({header.JumpTableStart}).");

            if (header.JumpTableStart < headerSize || header.JumpTableEnd > storageSize)
                throw new InvalidOperationException(
                    $"Invalid storage format: jump table [{header.JumpTableStart}, {header.JumpTableEnd}) is outside valid storage range.");
        }

        bool isAscii = (header.Flags & PersistentConstants.FLAG_TEXT_ASCII) != 0;
        int bytesPerChar = isAscii ? 1 : 2;

        if (header.TextLength < 0)
            throw new InvalidOperationException(
                $"Invalid storage format: text length {header.TextLength} is negative.");

        if (header.TextOffset < headerSize || header.TextOffset >= storageSize)
            throw new InvalidOperationException(
                $"Invalid storage format: text offset {header.TextOffset} is outside valid range [{headerSize}, {storageSize}).");

        long textEnd = header.TextOffset + (long)header.TextLength * bytesPerChar;
        if (textEnd > storageSize)
            throw new InvalidOperationException(
                $"Invalid storage format: text region [{header.TextOffset}, {textEnd}) exceeds storage size {storageSize}.");

        long deepestNodeOffset = PersistentConstants.NULL_OFFSET;
        if (header.DeepestNodeRawOffset != 0 && header.DeepestNodeRawOffset != PersistentConstants.NULL_OFFSET)
        {
            if (header.DeepestNodeRawOffset < headerSize || header.DeepestNodeRawOffset >= storageSize)
                throw new InvalidOperationException(
                    $"Invalid storage format: deepest internal node offset {header.DeepestNodeRawOffset} is outside valid range [{headerSize}, {storageSize}).");
            deepestNodeOffset = header.DeepestNodeRawOffset;
        }

        return new PersistentStorageMetadata(
            layout,
            header.RootOffset,
            header.TransitionOffset,
            header.JumpTableStart,
            header.JumpTableEnd,
            deepestNodeOffset,
            header.LrsDepth);
    }
}
