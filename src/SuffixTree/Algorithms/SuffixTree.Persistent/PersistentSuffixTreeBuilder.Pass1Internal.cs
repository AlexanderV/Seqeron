namespace SuffixTree.Persistent;

#pragma warning disable CA1001 // _childStore and _depthStore ownership finalized in FinalizeTree/SequentialFinalize
public partial class PersistentSuffixTreeBuilder
#pragma warning restore CA1001
{
    /// <summary>
    /// Pass 1 helper: process a single internal node. Reads child linked list,
    /// sorts children, writes child array, updates node header, builds parent graph in temp.
    /// </summary>
    private unsafe void Pass1Internal(
        long off, int nodeIdx, NodeLayout layout, bool useMmfUnchecked,
        MappedFileStorageProvider? mmfMain, byte* childBase, byte* tp,
        uint* cKeys, long* cOffs, int cBufMax,
        ref long batchPos, ref long batchStart, ref long batchEnd, int batchSize)
    {
        bool isCompact = !layout.OffsetIs64Bit;

        // Read jump slot from LeafCount field (placed by MaterializeJumpTable)
        long jumpSlotOffset = -1;
        bool hasJumpSlot = false;
        if (isCompact && IsHybrid)
        {
            uint jumpOff = (useMmfUnchecked && mmfMain != null)
                ? mmfMain.ReadUInt32Unchecked(off + NodeLayout.Compact.OffsetLeafCount)
                : _storage.ReadUInt32(off + NodeLayout.Compact.OffsetLeafCount);
            if (jumpOff > 0)
            {
                jumpSlotOffset = (long)jumpOff;
                hasJumpSlot = true;
            }
        }

        // Read ChildrenHead (linked list head in child store)
        long headIndex;
        if (useMmfUnchecked && mmfMain != null)
        {
            uint raw = mmfMain.ReadUInt32Unchecked(off + NodeLayout.Compact.OffsetChildrenHead);
            headIndex = raw == uint.MaxValue ? PersistentConstants.NULL_OFFSET : (long)raw;
        }
        else
        {
            headIndex = new PersistentSuffixTreeNode(_storage, off, layout).ChildrenHead;
        }

        // Walk child linked list → collect (key, childOffset)
        int childCount = 0;
        if (headIndex != PersistentConstants.NULL_OFFSET)
        {
            long ci = headIndex;
            if (childBase != null)
            {
                while (ci >= 0 && childCount < cBufMax)
                {
                    byte* ce = childBase + ci * CHILD_ENTRY_SIZE;
                    cKeys[childCount] = *(uint*)(ce + CE_OFF_KEY);
                    cOffs[childCount] = *(long*)(ce + CE_OFF_CHILD);
                    childCount++;
                    ci = *(int*)(ce + CE_OFF_NEXT);
                }
            }
            else
            {
                while (ci >= 0 && childCount < cBufMax)
                {
                    _childStoreAdapter.ReadEntry((int)ci, out uint entryKey, out int nextIndex, out long entryChildOffset);
                    cKeys[childCount] = entryKey;
                    cOffs[childCount] = entryChildOffset;
                    childCount++;
                    ci = nextIndex;
                }
            }
        }

        // Set parent index for each child in temp file
        for (int c = 0; c < childCount; c++)
        {
            int childIdx = NodeIndex(cOffs[c]);
            *(int*)(tp + (long)childIdx * TEMP_ENTRY + TEMP_OFF_PARENT) = nodeIdx;
        }
        // Set remaining count for this node (leaf count stays 0, computed in Pass 2)
        *(int*)(tp + (long)nodeIdx * TEMP_ENTRY + TEMP_OFF_REMAINING) = childCount;

        // Sort children by key (signed comparison: terminator -1 sorts first)
        for (int i = 1; i < childCount; i++)
        {
            uint ik = cKeys[i];
            long io = cOffs[i];
            int j = i - 1;
            while (j >= 0 && (int)cKeys[j] > (int)ik)
            {
                cKeys[j + 1] = cKeys[j];
                cOffs[j + 1] = cOffs[j];
                j--;
            }
            cKeys[j + 1] = ik;
            cOffs[j + 1] = io;
        }

        // Write sorted child array (batch-allocated)
        NodeLayout arrayLayout = hasJumpSlot ? NodeLayout.Large : layout;
        int entrySize = arrayLayout.ChildEntrySize;
        int totalBytes = checked(childCount * entrySize);

        if (batchPos + totalBytes > batchEnd)
        {
            int allocSz = Math.Max(batchSize, totalBytes);
            batchStart = _storage.Allocate(allocSz);
            batchEnd = batchStart + allocSz;
            batchPos = batchStart;
        }
        long arrayOffset = batchPos;
        batchPos += totalBytes;

        // Direct unsafe write to MMF when possible
        var mmf = _mmfStorage;
        if (mmf != null && arrayOffset + totalBytes <= mmf.Size)
        {
            byte* dst = mmf.RawPointer + arrayOffset;
            if (arrayLayout.OffsetIs64Bit)
            {
                for (int c = 0; c < childCount; c++)
                {
                    *(uint*)dst = cKeys[c];
                    *(long*)(dst + 4) = cOffs[c];
                    dst += 12;
                }
            }
            else
            {
                for (int c = 0; c < childCount; c++)
                {
                    *(uint*)dst = cKeys[c];
                    *(uint*)(dst + 4) = (uint)cOffs[c];
                    dst += 8;
                }
            }
        }
        else
        {
            for (int c = 0; c < childCount; c++)
            {
                long eOff = arrayOffset + (long)c * entrySize;
                _storage.WriteUInt32(eOff, cKeys[c]);
                if (arrayLayout.OffsetIs64Bit)
                    _storage.WriteInt64(eOff + 4, cOffs[c]);
                else
                    _storage.WriteUInt32(eOff + 4, (uint)cOffs[c]);
            }
        }

        // Update node header (ChildrenHead + ChildCount)
        int cCountRaw = childCount | (hasJumpSlot ? unchecked((int)0x80000000) : 0);

        if (hasJumpSlot)
        {
            _storage.WriteInt64(jumpSlotOffset, arrayOffset);
            if (useMmfUnchecked && mmfMain != null)
            {
                mmfMain.WriteUInt32Unchecked(off + (uint)layout.OffsetChildrenHead, (uint)jumpSlotOffset);
                mmfMain.WriteInt32Unchecked(off + layout.OffsetChildCount, cCountRaw);
            }
            else
            {
                var node = new PersistentSuffixTreeNode(_storage, off, layout);
                node.ChildrenHead = jumpSlotOffset;
                node.ChildCount = cCountRaw;
            }
        }
        else
        {
            if (useMmfUnchecked && mmfMain != null)
            {
                if (!layout.OffsetIs64Bit)
                    mmfMain.WriteUInt32Unchecked(off + (uint)layout.OffsetChildrenHead, (uint)arrayOffset);
                else
                    mmfMain.WriteInt64Unchecked(off + layout.OffsetChildrenHead, arrayOffset);
                mmfMain.WriteInt32Unchecked(off + layout.OffsetChildCount, cCountRaw);
            }
            else
            {
                var node = new PersistentSuffixTreeNode(_storage, off, layout);
                node.ChildrenHead = arrayOffset;
                node.ChildCount = cCountRaw;
            }
        }
    }
}
