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
        int idx = (int)headIndex;
        while (idx >= 0)
        {
            _childStoreAdapter.ReadEntry(idx, out uint entryKey, out int nextIndex, out long entryChildOffset);
            if (entryKey == key)
            {
                childOffset = entryChildOffset;
                return true;
            }
            idx = nextIndex;
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
            int idx = (int)headIndex;
            while (idx >= 0)
            {
                _childStoreAdapter.ReadEntry(idx, out uint entryKey, out int nextIndex, out _);
                if (entryKey == key)
                {
                    _childStoreAdapter.WriteChildOffset(idx, childOffset);
                    return;
                }
                idx = nextIndex;
            }
        }

        // New child — allocate entry and prepend to linked list
        int newIndex = _childStoreAdapter.AddEntry(
            headIndex == PersistentConstants.NULL_OFFSET ? -1 : (int)headIndex,
            key,
            childOffset);

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
