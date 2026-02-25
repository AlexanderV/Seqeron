namespace SuffixTree.Persistent;

#pragma warning disable CA1001 // _childStore and _depthStore ownership finalized in FinalizeTree/SequentialFinalize
public partial class PersistentSuffixTreeBuilder
#pragma warning restore CA1001
{
    private bool BuilderTryGetChild(long nodeOffset, uint key, out long childOffset)
    {
        // Read ChildrenHead — fast-path for MMF
        long headIndex;
        if (_mmfStorage != null && (_transitionOffset < 0 || nodeOffset < _transitionOffset))
        {
            // Compact: ChildrenHead is uint32 at nodeOffset + NodeLayout.Compact.OffsetChildrenHead
            uint raw = _mmfStorage.ReadUInt32Unchecked(nodeOffset + NodeLayout.Compact.OffsetChildrenHead);
            headIndex = raw == uint.MaxValue ? PersistentConstants.NULL_OFFSET : (long)raw;
        }
        else if (_mmfStorage != null) // Large: ChildrenHead is int64 at nodeOffset + NodeLayout.Large.OffsetChildrenHead
        {
            headIndex = _mmfStorage.ReadInt64Unchecked(nodeOffset + NodeLayout.Large.OffsetChildrenHead);
        }
        else
        {
            var node = new PersistentSuffixTreeNode(_storage, nodeOffset, LayoutOf(nodeOffset));
            headIndex = node.ChildrenHead;
        }

        if (headIndex == PersistentConstants.NULL_OFFSET)
        {
            childOffset = PersistentConstants.NULL_OFFSET;
            return false;
        }

        // Walk linked list in child store — fast-path for MMF child store
        if (_mmfChildStore != null)
        {
            int idx = (int)headIndex;
            while (idx >= 0)
            {
                long entryOff = (long)idx * CHILD_ENTRY_SIZE;
                uint entryKey = _mmfChildStore.ReadUInt32Unchecked(entryOff + CE_OFF_KEY);
                if (entryKey == key)
                {
                    childOffset = _mmfChildStore.ReadInt64Unchecked(entryOff + CE_OFF_CHILD);
                    return true;
                }
                idx = _mmfChildStore.ReadInt32Unchecked(entryOff + CE_OFF_NEXT);
            }
        }
        else
        {
            int idx = (int)headIndex;
            while (idx >= 0)
            {
                long entryOff = (long)idx * CHILD_ENTRY_SIZE;
                uint entryKey = _childStore.ReadUInt32(entryOff + CE_OFF_KEY);
                if (entryKey == key)
                {
                    childOffset = _childStore.ReadInt64(entryOff + CE_OFF_CHILD);
                    return true;
                }
                idx = _childStore.ReadInt32(entryOff + CE_OFF_NEXT);
            }
        }

        childOffset = PersistentConstants.NULL_OFFSET;
        return false;
    }

    private void BuilderSetChild(long nodeOffset, uint key, long childOffset)
    {
        // Read ChildrenHead — fast-path for MMF
        bool isCompactMmf = _mmfStorage != null && (_transitionOffset < 0 || nodeOffset < _transitionOffset);
        bool isLargeMmf = !isCompactMmf && _mmfStorage != null;
        long headIndex;
        if (isCompactMmf)
        {
            uint raw = _mmfStorage!.ReadUInt32Unchecked(nodeOffset + NodeLayout.Compact.OffsetChildrenHead);
            headIndex = raw == uint.MaxValue ? PersistentConstants.NULL_OFFSET : (long)raw;
        }
        else if (isLargeMmf) // Large: ChildrenHead is int64 at nodeOffset + NodeLayout.Large.OffsetChildrenHead
        {
            headIndex = _mmfStorage!.ReadInt64Unchecked(nodeOffset + NodeLayout.Large.OffsetChildrenHead);
        }
        else
        {
            var nodeLayout = LayoutOf(nodeOffset);
            var node = new PersistentSuffixTreeNode(_storage, nodeOffset, nodeLayout);
            headIndex = node.ChildrenHead;
        }

        // Check if key already exists (edge split replaces child pointer)
        if (headIndex != PersistentConstants.NULL_OFFSET)
        {
            if (_mmfChildStore != null)
            {
                int idx = (int)headIndex;
                while (idx >= 0)
                {
                    long entryOff = (long)idx * CHILD_ENTRY_SIZE;
                    uint entryKey = _mmfChildStore.ReadUInt32Unchecked(entryOff + CE_OFF_KEY);
                    if (entryKey == key)
                    {
                        _mmfChildStore.WriteInt64Unchecked(entryOff + CE_OFF_CHILD, childOffset);
                        return;
                    }
                    idx = _mmfChildStore.ReadInt32Unchecked(entryOff + CE_OFF_NEXT);
                }
            }
            else
            {
                int idx = (int)headIndex;
                while (idx >= 0)
                {
                    long entryOff = (long)idx * CHILD_ENTRY_SIZE;
                    uint entryKey = _childStore.ReadUInt32(entryOff + CE_OFF_KEY);
                    if (entryKey == key)
                    {
                        _childStore.WriteInt64(entryOff + CE_OFF_CHILD, childOffset);
                        return;
                    }
                    idx = _childStore.ReadInt32(entryOff + CE_OFF_NEXT);
                }
            }
        }

        // New child — allocate entry and prepend to linked list
        long newOff = _childStore.Allocate(CHILD_ENTRY_SIZE);
        int newIndex = (int)(newOff / CHILD_ENTRY_SIZE);

        if (_mmfChildStore != null)
        {
            _mmfChildStore.WriteUInt32Unchecked(newOff + CE_OFF_KEY, key);
            _mmfChildStore.WriteInt32Unchecked(newOff + CE_OFF_NEXT, headIndex == PersistentConstants.NULL_OFFSET ? -1 : (int)headIndex);
            _mmfChildStore.WriteInt64Unchecked(newOff + CE_OFF_CHILD, childOffset);
        }
        else
        {
            _childStore.WriteUInt32(newOff + CE_OFF_KEY, key);
            _childStore.WriteInt32(newOff + CE_OFF_NEXT, headIndex == PersistentConstants.NULL_OFFSET ? -1 : (int)headIndex);
            _childStore.WriteInt64(newOff + CE_OFF_CHILD, childOffset);
        }

        // Update node's head pointer to the new entry index
        if (isCompactMmf)
        {
            _mmfStorage!.WriteUInt32Unchecked(nodeOffset + NodeLayout.Compact.OffsetChildrenHead, (uint)newIndex);
        }
        else if (isLargeMmf) // Large: ChildrenHead is int64 at nodeOffset + NodeLayout.Large.OffsetChildrenHead
        {
            _mmfStorage!.WriteInt64Unchecked(nodeOffset + NodeLayout.Large.OffsetChildrenHead, (long)newIndex);
        }
        else
        {
            var nodeLayout = LayoutOf(nodeOffset);
            var node = new PersistentSuffixTreeNode(_storage, nodeOffset, nodeLayout);
            node.ChildrenHead = newIndex;
        }
    }
}
