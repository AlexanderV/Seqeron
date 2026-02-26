namespace SuffixTree.Persistent;

public partial class PersistentSuffixTreeBuilder
{
    private void ExtendTree(uint key)
    {
        _position++;
        _remainder++;
        _lastCreatedInternalNodeOffset = PersistentConstants.NULL_OFFSET;

        while (_remainder > 0)
        {
            if (_activeLength == 0)
                _activeEdgeIndex = _position;

            uint activeEdgeKey = GetSymbolAt(_activeEdgeIndex);

            // Fast-path: read active node depth via unchecked MMF (Start/End at offsets 0/4 for both layouts)
            uint activeNodeDepth;
            if (_mmfStorage != null)
            {
                uint aStart = ReadStartUnchecked(_activeNodeOffset);
                uint aEnd = ReadEndUnchecked(_activeNodeOffset);
                int aEdgeLen = (int)((aEnd == PersistentConstants.BOUNDLESS ? (uint)(_position + 1) : aEnd) - aStart);
                activeNodeDepth = GetBuildDepth(_activeNodeOffset) + (uint)aEdgeLen;
            }
            else
            {
                var activeNode = new PersistentSuffixTreeNode(_storage, _activeNodeOffset, LayoutOf(_activeNodeOffset));
                activeNodeDepth = GetNodeDepth(activeNode);
            }

            if (!BuilderTryGetChild(_activeNodeOffset, activeEdgeKey, out var nextChildOffset))
            {
                var leafOffset = CreateNode((uint)_position, PersistentConstants.BOUNDLESS, activeNodeDepth);
                BuilderSetChild(_activeNodeOffset, activeEdgeKey, leafOffset);
                AddSuffixLink(_activeNodeOffset);
            }
            else
            {
                // Read nextChild Start/End via unchecked path (offsets 0/4 same for both layouts)
                uint ncStart, ncEnd;
                if (_mmfStorage != null)
                {
                    ncStart = ReadStartUnchecked(nextChildOffset);
                    ncEnd = ReadEndUnchecked(nextChildOffset);
                }
                else
                {
                    var nextChild = new PersistentSuffixTreeNode(_storage, nextChildOffset, LayoutOf(nextChildOffset));
                    ncStart = nextChild.Start;
                    ncEnd = nextChild.End;
                }

                int edgeLen = (int)((ncEnd == PersistentConstants.BOUNDLESS ? (uint)(_position + 1) : ncEnd) - ncStart);
                if (_activeLength >= edgeLen)
                {
                    _activeEdgeIndex += edgeLen;
                    _activeLength -= edgeLen;
                    _activeNodeOffset = nextChildOffset;
                    continue;
                }

                if (GetSymbolAt((int)(ncStart + (uint)_activeLength)) == key)
                {
                    _activeLength++;
                    AddSuffixLink(_activeNodeOffset);
                    break;
                }

                // Split edge
                uint nextChildDFR = GetBuildDepth(nextChildOffset);
                long splitOffset = CreateNode(ncStart, ncStart + (uint)_activeLength, nextChildDFR);

                // Read split node's edge length for depth calc (offsets 0/4 same for both layouts)
                uint splitStart, splitEnd;
                if (_mmfStorage != null)
                {
                    splitStart = ReadStartUnchecked(splitOffset);
                    splitEnd = ReadEndUnchecked(splitOffset);
                }
                else
                {
                    var split = new PersistentSuffixTreeNode(_storage, splitOffset, LayoutOf(splitOffset));
                    splitStart = split.Start;
                    splitEnd = split.End;
                }

                BuilderSetChild(_activeNodeOffset, activeEdgeKey, splitOffset);

                int splitEdgeLen = (int)((splitEnd == PersistentConstants.BOUNDLESS ? (uint)(_position + 1) : splitEnd) - splitStart);
                uint splitEndDepth = GetBuildDepth(splitOffset) + (uint)splitEdgeLen;
                long leafOffset = CreateNode((uint)_position, PersistentConstants.BOUNDLESS, splitEndDepth);
                BuilderSetChild(splitOffset, key, leafOffset);

                // Update the original child's Start (offset 0 same for both layouts)
                uint newStart = ncStart + (uint)_activeLength;
                if (_mmfStorage != null)
                    WriteStartUnchecked(nextChildOffset, newStart);
                else
                    new PersistentSuffixTreeNode(_storage, nextChildOffset, LayoutOf(nextChildOffset)).Start = newStart;

                SetBuildDepth(nextChildOffset, splitEndDepth);
                BuilderSetChild(splitOffset, GetSymbolAt((int)newStart), nextChildOffset);

                AddSuffixLink(splitOffset);
            }

            _remainder--;
            if (_activeNodeOffset == _rootOffset && _activeLength > 0)
            {
                _activeLength--;
                _activeEdgeIndex = _position - _remainder + 1;
            }
            else if (_activeNodeOffset != _rootOffset)
            {
                long suffLink;
                if (_crossZoneSuffixLinks.TryGetValue(_activeNodeOffset, out long resolved))
                    suffLink = resolved;
                else if (IsCompactMmf(_activeNodeOffset))
                    suffLink = ReadSuffixLinkUnchecked(_activeNodeOffset);
                else if (_mmfStorage != null) // large zone: SuffixLink is int64 at offset 8
                    suffLink = _mmfStorage.ReadInt64Unchecked(_activeNodeOffset + 8);
                else
                {
                    var node = new PersistentSuffixTreeNode(_storage, _activeNodeOffset, LayoutOf(_activeNodeOffset));
                    suffLink = node.SuffixLink;
                }
                _activeNodeOffset = suffLink != PersistentConstants.NULL_OFFSET ? suffLink : _rootOffset;
            }
        }
    }

    /// <summary>
    /// Allocates storage with automatic compact→large transition.
    /// When the Compact layout's address space is exhausted, switches to Large
    /// and continues allocating in the same storage.
    /// </summary>
    private long AllocateChecked(int size)
    {
        if (!_layout.OffsetIs64Bit)
        {
            // Still in compact mode — check if this allocation would overflow.
            // For large compact zones, transition early to leave room for a jump
            // table reserve block whose offsets are still compact-addressable.
            long currentSize = _storage.Size;
            long jumpReserveBytes = Math.Max(
                _compactOffsetLimit / JUMP_TABLE_RESERVE_FRACTION,
                JUMP_TABLE_MIN_RESERVE);
            bool largeCompactZone = _compactOffsetLimit > jumpReserveBytes * 4;
            long effectiveLimit = largeCompactZone
                ? _compactOffsetLimit - jumpReserveBytes
                : _compactOffsetLimit;

            if ((currentSize + size) > effectiveLimit)
            {
                // Record where compact nodes end (before any reserve)
                _compactNodesEnd = currentSize;

                // Reserve jump table block within the compact address space
                if (largeCompactZone)
                {
                    _jumpTableReserveStart = _storage.Allocate((int)jumpReserveBytes);
                    _jumpTableReserveEnd = _jumpTableReserveStart + jumpReserveBytes;
                }

                // Transition to Large layout
                _transitionOffset = _storage.Size;
                _layout = NodeLayout.Large;
                // The caller passed the Compact node size; now that _layout
                // has switched to Large/SlimLarge we must allocate the larger block.
                size = _layout.NodeSize;
            }
        }
        return _storage.Allocate(size);
    }

    private long CreateNode(uint start, uint end, uint depthFromRoot)
    {
        _nodeCount++;
        long offset = AllocateChecked(_layout.NodeSize);
        _maxNodeEndOffset = offset + _layout.NodeSize;

        // Fast-path: Compact layout on MMF — bulk-write all 24 bytes in one shot
        // instead of 6 individual IStorageProvider calls (each with ThrowIfDisposed+CheckBounds).
        if (_mmfStorage != null && !_layout.OffsetIs64Bit)
        {
            unsafe
            {
                byte* p = _mmfStorage.RawPointer + offset;
                *(uint*)(p + 0) = start;                   // Start
                *(uint*)(p + 4) = end;                     // End
                *(uint*)(p + 8) = uint.MaxValue;           // SuffixLink = NULL (compact sentinel)
                *(uint*)(p + 12) = 0;                      // LeafCount = 0
                *(uint*)(p + 16) = uint.MaxValue;          // ChildrenHead = NULL (compact sentinel)
                *(int*)(p + 20) = 0;                       // ChildCount = 0
            }
        }
        else if (_mmfStorage != null) // Large layout on MMF — bulk-write all 32 bytes
        {
            unsafe
            {
                byte* p = _mmfStorage.RawPointer + offset;
                *(uint*)(p + 0) = start;                   // Start
                *(uint*)(p + 4) = end;                     // End
                *(long*)(p + 8) = PersistentConstants.NULL_OFFSET; // SuffixLink (int64)
                *(uint*)(p + 16) = 0;                      // LeafCount
                *(long*)(p + 20) = PersistentConstants.NULL_OFFSET; // ChildrenHead (int64)
                *(int*)(p + 28) = 0;                       // ChildCount
            }
        }
        else
        {
            var node = new PersistentSuffixTreeNode(_storage, offset, _layout);
            node.Start = start;
            node.End = end;
            node.SuffixLink = PersistentConstants.NULL_OFFSET;
            node.ChildrenHead = PersistentConstants.NULL_OFFSET;
            node.LeafCount = 0;
        }

        SetBuildDepth(offset, depthFromRoot);
        return offset;
    }

    private void AddSuffixLink(long nodeOffset)
    {
        if (_lastCreatedInternalNodeOffset != PersistentConstants.NULL_OFFSET)
        {
            bool sourceIsCompact = _transitionOffset < 0 || _lastCreatedInternalNodeOffset < _transitionOffset;
            bool targetInLargeZone = _transitionOffset >= 0 && nodeOffset >= _transitionOffset;

            if (sourceIsCompact && targetInLargeZone)
            {
                // Cross-zone link: 64-bit offset can't fit in a compact 32-bit field.
                if (IsCompactMmf(_lastCreatedInternalNodeOffset))
                    WriteSuffixLinkUnchecked(_lastCreatedInternalNodeOffset, PersistentConstants.NULL_OFFSET);
                else
                    new PersistentSuffixTreeNode(_storage, _lastCreatedInternalNodeOffset, LayoutOf(_lastCreatedInternalNodeOffset)).SuffixLink = PersistentConstants.NULL_OFFSET;

                _deferredSuffixLinks.Add((_lastCreatedInternalNodeOffset, nodeOffset));
                _crossZoneSuffixLinks[_lastCreatedInternalNodeOffset] = nodeOffset;
            }
            else
            {
                if (IsCompactMmf(_lastCreatedInternalNodeOffset))
                    WriteSuffixLinkUnchecked(_lastCreatedInternalNodeOffset, nodeOffset);
                else if (_mmfStorage != null) // large zone: SuffixLink is int64 at offset 8
                    _mmfStorage.WriteInt64Unchecked(_lastCreatedInternalNodeOffset + 8, nodeOffset);
                else
                    new PersistentSuffixTreeNode(_storage, _lastCreatedInternalNodeOffset, LayoutOf(_lastCreatedInternalNodeOffset)).SuffixLink = nodeOffset;
            }
        }
        _lastCreatedInternalNodeOffset = nodeOffset;
    }

    /// <summary>
    /// Returns the layout a node was written with, determined by its offset relative
    /// to the transition boundary. No dictionary lookup — O(1), zero allocations.
    /// </summary>
    [System.Runtime.CompilerServices.MethodImpl(System.Runtime.CompilerServices.MethodImplOptions.AggressiveInlining)]
    private NodeLayout LayoutOf(long nodeOffset)
        => (_transitionOffset >= 0 && nodeOffset >= _transitionOffset)
            ? NodeLayout.Large
            : _initialLayout;

    /// <summary>Whether the node is in compact zone backed by MMF — enables unchecked direct pointer access.</summary>
    [System.Runtime.CompilerServices.MethodImpl(System.Runtime.CompilerServices.MethodImplOptions.AggressiveInlining)]
    private bool IsCompactMmf(long nodeOffset)
        => _mmfStorage != null && (_transitionOffset < 0 || nodeOffset < _transitionOffset);

    // ──── Unchecked compact-node field readers (bypass interface dispatch + bounds checks) ────

    [System.Runtime.CompilerServices.MethodImpl(System.Runtime.CompilerServices.MethodImplOptions.AggressiveInlining)]
    private uint ReadStartUnchecked(long nodeOffset) => _mmfStorage!.ReadUInt32Unchecked(nodeOffset + 0);

    [System.Runtime.CompilerServices.MethodImpl(System.Runtime.CompilerServices.MethodImplOptions.AggressiveInlining)]
    private uint ReadEndUnchecked(long nodeOffset) => _mmfStorage!.ReadUInt32Unchecked(nodeOffset + 4);

    [System.Runtime.CompilerServices.MethodImpl(System.Runtime.CompilerServices.MethodImplOptions.AggressiveInlining)]
    private long ReadSuffixLinkUnchecked(long nodeOffset)
    {
        uint raw = _mmfStorage!.ReadUInt32Unchecked(nodeOffset + 8);
        return raw == uint.MaxValue ? PersistentConstants.NULL_OFFSET : (long)raw;
    }

    [System.Runtime.CompilerServices.MethodImpl(System.Runtime.CompilerServices.MethodImplOptions.AggressiveInlining)]
    private void WriteSuffixLinkUnchecked(long nodeOffset, long value)
        => _mmfStorage!.WriteUInt32Unchecked(nodeOffset + 8,
            value == PersistentConstants.NULL_OFFSET ? uint.MaxValue : (uint)value);

    [System.Runtime.CompilerServices.MethodImpl(System.Runtime.CompilerServices.MethodImplOptions.AggressiveInlining)]
    private void WriteStartUnchecked(long nodeOffset, uint value)
        => _mmfStorage!.WriteUInt32Unchecked(nodeOffset + 0, value);

    [System.Runtime.CompilerServices.MethodImpl(System.Runtime.CompilerServices.MethodImplOptions.AggressiveInlining)]
    private uint GetSymbolAt(int index)
    {
        if (index > _position) return PersistentConstants.TERMINATOR_KEY;
        var raw = _rawString;
        if (raw != null)
            return index < raw.Length ? (uint)raw[index] : PersistentConstants.TERMINATOR_KEY;
        return (index < _text.Length) ? (uint)_text[index] : PersistentConstants.TERMINATOR_KEY;
    }

    private int LengthOf(PersistentSuffixTreeNode node)
        => (int)((node.End == PersistentConstants.BOUNDLESS ? (uint)(_position + 1) : node.End) - node.Start);

    private uint GetNodeDepth(PersistentSuffixTreeNode node)
        => GetBuildDepth(node.Offset) + (uint)LengthOf(node);

    /// <summary>
    /// Checks whether every character in the input text is pure ASCII (0-127).
    /// When true, the builder stores text as 1 byte/char, halving disk usage.
    /// Uses SIMD-accelerated <c>System.Text.Ascii.IsValid</c> when available.
    /// </summary>
    private bool IsTextAscii()
    {
        if (_rawString != null)
            return System.Text.Ascii.IsValid(_rawString);
        for (int i = 0; i < _text.Length; i++)
            if (_text[i] > 127) return false;
        return true;
    }

    // ──────────────── Off-heap depth helpers ────────────────

    /// <summary>Maps a node storage offset to a sequential 0-based index.</summary>
    [System.Runtime.CompilerServices.MethodImpl(System.Runtime.CompilerServices.MethodImplOptions.AggressiveInlining)]
    private int NodeIndex(long nodeOffset)
    {
        int hdr = HeaderSize;
        if (_transitionOffset >= 0 && nodeOffset >= _transitionOffset)
        {
            // Use _compactNodesEnd (excludes jump reserve) for accurate compact node count
            long compactBound = _compactNodesEnd >= 0 ? _compactNodesEnd : _transitionOffset;
            int compactCount = (int)((compactBound - hdr) / _initialLayout.NodeSize);
            return compactCount + (int)((nodeOffset - _transitionOffset) / NodeLayout.Large.NodeSize);
        }
        return (int)((nodeOffset - hdr) / _initialLayout.NodeSize);
    }

    /// <summary>Reads DepthFromRoot for a node from off-heap depth store.</summary>
    [System.Runtime.CompilerServices.MethodImpl(System.Runtime.CompilerServices.MethodImplOptions.AggressiveInlining)]
    private uint GetBuildDepth(long nodeOffset)
    {
        long off = (long)NodeIndex(nodeOffset) * 4;
        return _depthStoreAdapter.ReadUInt32(off);
    }

    /// <summary>Writes DepthFromRoot for a node to off-heap depth store.</summary>
    [System.Runtime.CompilerServices.MethodImpl(System.Runtime.CompilerServices.MethodImplOptions.AggressiveInlining)]
    private void SetBuildDepth(long nodeOffset, uint value)
    {
        int idx = NodeIndex(nodeOffset);
        long off = (long)idx * 4;
        _depthStoreAdapter.EnsureSizeAtLeast(off + 4);
        _depthStoreAdapter.WriteUInt32(off, value);
    }
}
