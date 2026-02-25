using System.Buffers;
using System.Text;

namespace SuffixTree.Persistent;

public partial class PersistentSuffixTreeBuilder
{
    private void FinalizeTree()
    {
        if (IsHybrid)
        {
            _progress?.Report(("Jump table materialization", 91.0));
            MaterializeJumpTable();
        }

        // Sequential multi-pass finalize: replaces DFS to eliminate random I/O.
        // Disposes child store and depth store internally when no longer needed.
        SequentialFinalize(_rootOffset);

        _progress?.Report(("DFS: Serializing text", 98.0));
        var (textOffset, textIsAscii) = SerializeTextPayload();
        WriteHeader(textOffset, textIsAscii);
        ClearDeferredJumpMetadata();
    }

    private (long TextOffset, bool TextIsAscii) SerializeTextPayload()
    {
        bool textIsAscii = IsTextAscii();
        int bytesPerChar = textIsAscii ? 1 : 2;
        long textByteLen = (long)_text.Length * bytesPerChar;
        long textOffset = _storage.Allocate((int)textByteLen);
        const int ChunkChars = 4096;
        byte[] chunkBuf = ArrayPool<byte>.Shared.Rent(ChunkChars * 2);
        try
        {
            int written = 0;
            while (written < _text.Length)
            {
                int remaining = _text.Length - written;
                int chunkLen = remaining < ChunkChars ? remaining : ChunkChars;
                var charSpan = _text.Slice(written, chunkLen);
                int byteCount = textIsAscii
                    ? Encoding.ASCII.GetBytes(charSpan, chunkBuf.AsSpan())
                    : Encoding.Unicode.GetBytes(charSpan, chunkBuf.AsSpan());
                _storage.WriteBytes(textOffset + (long)written * bytesPerChar, chunkBuf, 0, byteCount);
                written += chunkLen;
            }
        }
        finally
        {
            ArrayPool<byte>.Shared.Return(chunkBuf);
        }

        return (textOffset, textIsAscii);
    }

    private void WriteHeader(long textOffset, bool textIsAscii)
    {
        _storage.WriteInt64(PersistentConstants.HEADER_OFFSET_MAGIC, PersistentConstants.MAGIC_NUMBER);
        _storage.WriteInt32(PersistentConstants.HEADER_OFFSET_VERSION, 6);
        _storage.WriteInt64(PersistentConstants.HEADER_OFFSET_ROOT, _rootOffset);
        _storage.WriteInt64(PersistentConstants.HEADER_OFFSET_TEXT_OFF, textOffset);
        _storage.WriteInt32(PersistentConstants.HEADER_OFFSET_TEXT_LEN, _text.Length);
        _storage.WriteInt32(PersistentConstants.HEADER_OFFSET_NODE_COUNT, _nodeCount);
        _storage.WriteInt32(PersistentConstants.HEADER_OFFSET_FLAGS, textIsAscii ? PersistentConstants.FLAG_TEXT_ASCII : 0);
        _storage.WriteInt64(PersistentConstants.HEADER_OFFSET_SIZE, _storage.Size);
        _storage.WriteInt64(PersistentConstants.HEADER_OFFSET_TRANSITION, IsHybrid ? _transitionOffset : -1);
        _storage.WriteInt64(PersistentConstants.HEADER_OFFSET_JUMP_START, IsHybrid ? _jumpTableStart : -1);
        _storage.WriteInt64(PersistentConstants.HEADER_OFFSET_JUMP_END, IsHybrid ? _jumpTableEnd : -1);
        _storage.WriteInt64(PersistentConstants.HEADER_OFFSET_DEEPEST_NODE, _deepestInternalNodeOffset);
        _storage.WriteInt32(PersistentConstants.HEADER_OFFSET_LRS_DEPTH, _lrsDepth);
        _storage.WriteInt32(PersistentConstants.HEADER_OFFSET_BASE_NODE_SIZE, _initialLayout.NodeSize);
    }

    private void ClearDeferredJumpMetadata()
    {
        CrossZoneSuffixLinkCount = _deferredSuffixLinks.Count;
        _deferredSuffixLinks.Clear();
        _crossZoneSuffixLinks.Clear();
    }

    private void MaterializeJumpTable()
    {
        int suffixLinkJumps = _deferredSuffixLinks.Count;
        long compactEnd = _compactNodesEnd >= 0
            ? _compactNodesEnd
            : Math.Min(_transitionOffset, _maxNodeEndOffset);
        var mmfMain = _mmfStorage;

        if (TryMaterializeReservedJumpTable(suffixLinkJumps, compactEnd, mmfMain))
            return;

        MaterializeFallbackJumpTable(suffixLinkJumps, compactEnd, mmfMain);
    }

    // Single-pass path when reserve is pre-allocated (genome-scale builds).
    private bool TryMaterializeReservedJumpTable(
        int suffixLinkJumps,
        long compactEnd,
        MappedFileStorageProvider? mmfMain)
    {
        if (_jumpTableReserveStart < 0)
            return false;

        long reserveCapacity = _jumpTableReserveEnd - _jumpTableReserveStart;
        _jumpTableStart = _jumpTableReserveStart;
        int slotIndex = WriteSuffixLinkJumpEntries(suffixLinkJumps, _jumpTableStart);

        for (long offset = HeaderSize; offset < compactEnd; offset += _initialLayout.NodeSize)
        {
            if (ReadNodeEnd(offset, mmfMain) == PersistentConstants.BOUNDLESS)
                continue;

            long jumpEntryOffset = _jumpTableStart + (long)slotIndex * 8;
            if (jumpEntryOffset + 8 > _jumpTableReserveEnd)
            {
                throw new InvalidOperationException(
                    $"Jump table requires more than {reserveCapacity:N0} reserved bytes. " +
                    "The input exceeds the jump table reserve capacity.");
            }

            WriteChildArrayJumpSlot(offset, jumpEntryOffset, mmfMain);
            slotIndex++;
        }

        SetJumpTableRange(slotIndex);
        return true;
    }

    // Fallback two-pass path for non-reserved mode (tiny compact limits in tests).
    private void MaterializeFallbackJumpTable(
        int suffixLinkJumps,
        long compactEnd,
        MappedFileStorageProvider? mmfMain)
    {
        int childArrayJumps = CountCompactInternalNodes(compactEnd, mmfMain);
        int totalEntries = suffixLinkJumps + childArrayJumps;
        if (totalEntries == 0)
        {
            ResetJumpTableRange();
            return;
        }

        int tableSize = totalEntries * 8;
        _jumpTableStart = _storage.Allocate(tableSize);
        _jumpTableEnd = _jumpTableStart + tableSize;
        int slotIdx = WriteSuffixLinkJumpEntries(suffixLinkJumps, _jumpTableStart);

        for (long offset = HeaderSize; offset < compactEnd; offset += _initialLayout.NodeSize)
        {
            if (ReadNodeEnd(offset, mmfMain) == PersistentConstants.BOUNDLESS)
                continue;

            long jumpEntryOffset = _jumpTableStart + (long)slotIdx * 8;
            WriteChildArrayJumpSlot(offset, jumpEntryOffset, mmfMain);
            slotIdx++;
        }
    }

    private int CountCompactInternalNodes(long compactEnd, MappedFileStorageProvider? mmfMain)
    {
        int internalCount = 0;
        for (long offset = HeaderSize; offset < compactEnd; offset += _initialLayout.NodeSize)
        {
            if (ReadNodeEnd(offset, mmfMain) != PersistentConstants.BOUNDLESS)
                internalCount++;
        }

        return internalCount;
    }

    private int WriteSuffixLinkJumpEntries(int suffixLinkJumps, long jumpTableStart)
    {
        int slotIdx = 0;
        for (int i = 0; i < suffixLinkJumps; i++, slotIdx++)
        {
            var (compactNodeOffset, largeTargetOffset) = _deferredSuffixLinks[i];
            long jumpEntryOffset = jumpTableStart + (long)slotIdx * 8;
            _storage.WriteInt64(jumpEntryOffset, largeTargetOffset);

            var compactLayout = LayoutOf(compactNodeOffset);
            var node = new PersistentSuffixTreeNode(_storage, compactNodeOffset, compactLayout);
            node.SuffixLink = jumpEntryOffset;
        }

        return slotIdx;
    }

    private uint ReadNodeEnd(long offset, MappedFileStorageProvider? mmfMain)
        => mmfMain != null
            ? mmfMain.ReadUInt32Unchecked(offset + 4)
            : _storage.ReadUInt32(offset + 4);

    private void WriteChildArrayJumpSlot(long nodeOffset, long jumpEntryOffset, MappedFileStorageProvider? mmfMain)
    {
        if (mmfMain != null)
            mmfMain.WriteUInt32Unchecked(nodeOffset + NodeLayout.Compact.OffsetLeafCount, (uint)jumpEntryOffset);
        else
            _storage.WriteUInt32(nodeOffset + NodeLayout.Compact.OffsetLeafCount, (uint)jumpEntryOffset);
    }

    private void SetJumpTableRange(int totalEntries)
    {
        if (totalEntries == 0)
        {
            ResetJumpTableRange();
            return;
        }

        _jumpTableEnd = _jumpTableStart + totalEntries * 8L;
    }

    private void ResetJumpTableRange()
    {
        _jumpTableStart = -1;
        _jumpTableEnd = -1;
    }
}
