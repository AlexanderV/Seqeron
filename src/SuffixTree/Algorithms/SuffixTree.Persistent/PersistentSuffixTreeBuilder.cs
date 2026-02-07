using System;
using System.Buffers;
using System.Buffers.Binary;
using System.Collections.Generic;
using System.Text;

namespace SuffixTree.Persistent;

/// <summary>
/// Handles construction of a persistent suffix tree using Ukkonen's algorithm.
/// Writes nodes and child entries directly to the storage provider.
/// <para>
/// <b>Hybrid continuation</b>: construction starts with the Compact (28-byte, uint32)
/// layout. When the compact address space (~4 GB) is exhausted, the builder
/// seamlessly switches to the Large (40-byte, int64) layout and continues
/// in the same storage without restarting. A small jump table at the
/// compact/large boundary allows compact-zone nodes to reference large-zone
/// addresses through indirection entries.
/// </para>
/// </summary>
public class PersistentSuffixTreeBuilder
{
    private readonly IStorageProvider _storage;
    private readonly NodeLayout _initialLayout;  // layout before any transition (always the compact variant)
    private NodeLayout _layout;                // switches Compact → Large on overflow
    private long _compactOffsetLimit = NodeLayout.CompactMaxOffset;
    private long _rootOffset;
    private ITextSource _text = new StringTextSource(string.Empty);
    private int _nodeCount = 0;

    // Hybrid continuation state
    private long _transitionOffset = -1;       // -1 = no transition (pure Compact or pure Large)
    private long _jumpTableStart = -1;         // start of contiguous jump table
    private long _jumpTableEnd = -1;           // end of the jump table area
    private readonly List<long> _jumpEntries = new();                     // allocated jump slot offsets

    // Algorithm state
    private long _activeNodeOffset;
    private int _activeEdgeIndex;
    private int _activeLength;
    private int _remainder;
    private int _position = -1;
    private long _lastCreatedInternalNodeOffset = PersistentConstants.NULL_OFFSET;
    private readonly Dictionary<long, List<(uint Key, long ChildOffset)>> _children = new();
    private bool _built;

    // Deferred cross-zone suffix links: compact-node → large-zone target
    // Written at FinalizeTree when the jump table is materialized.
    private readonly List<(long CompactNodeOffset, long LargeTargetOffset)> _deferredSuffixLinks = new();

    // Cross-zone suffix link resolution for builder: compact-node offset →
    // correct large-zone target (avoids uint32 truncation when reading back
    // from storage during the Ukkonen walk).  Only populated when transition occurs.
    private readonly Dictionary<long, long> _crossZoneSuffixLinks = new();

    // Pre-computed during CalculateLeafCount — the internal node with maximum depth.
    private long _deepestInternalNodeOffset = PersistentConstants.NULL_OFFSET;

    /// <summary>
    /// Override the compact address-space limit. Used by <see cref="PersistentSuffixTreeFactory"/>
    /// and tests to control the Compact → Large promotion threshold.
    /// </summary>
    internal long CompactOffsetLimit
    {
        get => _compactOffsetLimit;
        set => _compactOffsetLimit = value;
    }

    /// <summary>The offset where the compact→large transition occurred, or -1.</summary>
    internal long TransitionOffset => _transitionOffset;

    /// <summary>The first offset of the contiguous jump table, or -1 if no transition.</summary>
    internal long JumpTableStart => _jumpTableStart;

    /// <summary>The first offset after the jump table, or -1 if no transition.</summary>
    internal long JumpTableEnd => _jumpTableEnd;

    /// <summary>Whether a compact→large transition occurred during the build.</summary>
    internal bool IsHybrid => _transitionOffset >= 0;

    /// <summary>
    /// Offset of the deepest internal node (maximum DepthFromRoot + edge length).
    /// Computed during <see cref="Build"/> at zero extra cost (piggybacks on leaf-count traversal).
    /// Returns <see cref="PersistentConstants.NULL_OFFSET"/> if Build has not been called.
    /// </summary>
    internal long DeepestInternalNodeOffset => _deepestInternalNodeOffset;

    public PersistentSuffixTreeBuilder(IStorageProvider storage, NodeLayout? layout = null)
    {
        _storage = storage;
        _layout = layout ?? NodeLayout.Compact;
        _initialLayout = _layout;

        // Allocate header — use the larger v5 header size to leave room for hybrid fields.
        // For pure Compact/Large builds the extra 16 bytes are unused but harmless.
        _storage.Allocate(PersistentConstants.HEADER_SIZE_V5);

        _rootOffset = _storage.Allocate(_layout.NodeSize);
        _nodeCount = 1;
        var root = new PersistentSuffixTreeNode(_storage, _rootOffset, _layout);
        root.Start = 0;
        root.End = 0;
        root.SuffixLink = PersistentConstants.NULL_OFFSET;
        root.ChildrenHead = PersistentConstants.NULL_OFFSET;

        _activeNodeOffset = _rootOffset;
    }

    public long Build(ITextSource text)
    {
        ArgumentNullException.ThrowIfNull(text);
        if (_built)
            throw new InvalidOperationException("Build() has already been called. Create a new builder instance to build another tree.");
        _built = true;

        _text = text;
        if (text.Length > 0)
        {
            for (int i = 0; i < text.Length; i++)
            {
                ExtendTree((uint)text[i]);
            }
            ExtendTree(PersistentConstants.TERMINATOR_KEY);
        }

        FinalizeTree();
        return _rootOffset;
    }

    private void ExtendTree(uint key)
    {
        _position++;
        _remainder++;
        _lastCreatedInternalNodeOffset = PersistentConstants.NULL_OFFSET;

        while (_remainder > 0)
        {
            if (_activeLength == 0)
                _activeEdgeIndex = _position;

            var activeLayout = LayoutOf(_activeNodeOffset);
            var activeNode = new PersistentSuffixTreeNode(_storage, _activeNodeOffset, activeLayout);
            uint activeEdgeKey = GetSymbolAt(_activeEdgeIndex);

            if (!BuilderTryGetChild(_activeNodeOffset, activeEdgeKey, out var nextChildOffset))
            {
                var leafOffset = CreateNode((uint)_position, PersistentConstants.BOUNDLESS, GetNodeDepth(activeNode));
                BuilderSetChild(_activeNodeOffset, activeEdgeKey, leafOffset);
                AddSuffixLink(_activeNodeOffset);
            }
            else
            {
                var nextChildLayout = LayoutOf(nextChildOffset);
                var nextChild = new PersistentSuffixTreeNode(_storage, nextChildOffset, nextChildLayout);
                int edgeLen = LengthOf(nextChild);
                if (_activeLength >= edgeLen)
                {
                    _activeEdgeIndex += edgeLen;
                    _activeLength -= edgeLen;
                    _activeNodeOffset = nextChildOffset;
                    continue;
                }

                if (GetSymbolAt((int)(nextChild.Start + (uint)_activeLength)) == key)
                {
                    _activeLength++;
                    AddSuffixLink(_activeNodeOffset);
                    break;
                }

                // Split edge — new split node uses CURRENT layout (may be Large if transitioned)
                long splitOffset = CreateNode(nextChild.Start, nextChild.Start + (uint)_activeLength, nextChild.DepthFromRoot);
                var splitLayout = LayoutOf(splitOffset);
                var split = new PersistentSuffixTreeNode(_storage, splitOffset, splitLayout);
                BuilderSetChild(_activeNodeOffset, activeEdgeKey, splitOffset);

                long leafOffset = CreateNode((uint)_position, PersistentConstants.BOUNDLESS, split.DepthFromRoot + (uint)LengthOf(split));
                BuilderSetChild(splitOffset, key, leafOffset);

                // Update the original child's metadata using its original layout
                nextChild.Start += (uint)_activeLength;
                nextChild.DepthFromRoot = split.DepthFromRoot + (uint)LengthOf(split);
                BuilderSetChild(splitOffset, GetSymbolAt((int)nextChild.Start), nextChildOffset);

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
                // Read suffix link — use in-memory dictionary for cross-zone
                // links to avoid uint32 truncation on compact nodes (C11).
                long suffLink;
                if (_crossZoneSuffixLinks.TryGetValue(_activeNodeOffset, out long resolved))
                    suffLink = resolved;
                else
                {
                    var nodeLayout = LayoutOf(_activeNodeOffset);
                    var node = new PersistentSuffixTreeNode(_storage, _activeNodeOffset, nodeLayout);
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
            // Still in compact mode — check if this allocation would overflow
            long currentSize = _storage.Size;
            if ((currentSize + size) > _compactOffsetLimit)
            {
                // Transition to Large layout
                _transitionOffset = currentSize;
                _layout = NodeLayout.Large;
                // The caller passed the Compact node size; now that _layout
                // has switched to Large we must allocate the larger block.
                size = _layout.NodeSize;
            }
        }
        return _storage.Allocate(size);
    }

    private long CreateNode(uint start, uint end, uint depthFromRoot)
    {
        _nodeCount++;
        long offset = AllocateChecked(_layout.NodeSize);
        var node = new PersistentSuffixTreeNode(_storage, offset, _layout);
        node.Start = start;
        node.End = end;
        node.DepthFromRoot = depthFromRoot;
        node.SuffixLink = PersistentConstants.NULL_OFFSET;
        node.ChildrenHead = PersistentConstants.NULL_OFFSET;
        node.LeafCount = 0;
        return offset;
    }

    private void AddSuffixLink(long nodeOffset)
    {
        if (_lastCreatedInternalNodeOffset != PersistentConstants.NULL_OFFSET)
        {
            var sourceLayout = LayoutOf(_lastCreatedInternalNodeOffset);
            bool sourceIsCompact = !sourceLayout.OffsetIs64Bit;
            bool targetInLargeZone = _transitionOffset >= 0 && nodeOffset >= _transitionOffset;

            // Always write the suffix link directly so the Ukkonen algorithm
            // can follow it during the remaining build steps.
            var lastNode = new PersistentSuffixTreeNode(_storage, _lastCreatedInternalNodeOffset, sourceLayout);
            lastNode.SuffixLink = nodeOffset;

            if (sourceIsCompact && targetInLargeZone)
            {
                // Also record for finalization: the reader needs a jump-table
                // entry so that a compact node's uint32 SuffixLink field can
                // be redirected to an int64 slot when the offset exceeds 32 bits.
                _deferredSuffixLinks.Add((_lastCreatedInternalNodeOffset, nodeOffset));

                // Keep the correct target in memory so ExtendTree can follow
                // cross-zone suffix links accurately even if the uint32 write
                // truncated the value (C11 safety guard).
                _crossZoneSuffixLinks[_lastCreatedInternalNodeOffset] = nodeOffset;
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

    private uint GetSymbolAt(int index)
    {
        if (index > _position) return PersistentConstants.TERMINATOR_KEY;
        return (index < _text.Length) ? (uint)_text[index] : PersistentConstants.TERMINATOR_KEY;
    }

    private int LengthOf(PersistentSuffixTreeNode node)
        => (int)((node.End == PersistentConstants.BOUNDLESS ? (uint)(_position + 1) : node.End) - node.Start);

    private uint GetNodeDepth(PersistentSuffixTreeNode node)
        => node.DepthFromRoot + (uint)LengthOf(node);

    // ──────────────── Finalize ────────────────

    private void FinalizeTree()
    {
        // Calculate leaf counts (zone-aware)
        CalculateLeafCount(_rootOffset);

        // Materialize a contiguous jump table for ALL cross-zone references:
        // 1. Deferred suffix links (compact node → large zone target)
        // 2. Child array jumps (compact parent with large-zone children)
        if (IsHybrid)
            MaterializeJumpTable();

        // Write sorted children arrays to storage (zone-aware)
        WriteChildrenArrays();

        // Store text in storage for true persistence (chunked write, no full-string copy)
        long textByteLen = (long)_text.Length * 2;
        if (textByteLen > int.MaxValue)
            throw new InvalidOperationException(
                $"Text length {_text.Length} exceeds maximum serializable size ({int.MaxValue / 2} characters).");
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
                int byteCount = Encoding.Unicode.GetBytes(charSpan, chunkBuf.AsSpan());
                _storage.WriteBytes(textOffset + (long)written * 2, chunkBuf, 0, byteCount);
                written += chunkLen;
            }
        }
        finally
        {
            ArrayPool<byte>.Shared.Return(chunkBuf);
        }

        // Write Header
        int version = IsHybrid ? 5 : _layout.Version;
        _storage.WriteInt64(PersistentConstants.HEADER_OFFSET_MAGIC, PersistentConstants.MAGIC_NUMBER);
        _storage.WriteInt32(PersistentConstants.HEADER_OFFSET_VERSION, version);
        _storage.WriteInt64(PersistentConstants.HEADER_OFFSET_ROOT, _rootOffset);
        _storage.WriteInt64(PersistentConstants.HEADER_OFFSET_TEXT_OFF, textOffset);
        _storage.WriteInt32(PersistentConstants.HEADER_OFFSET_TEXT_LEN, _text.Length);
        _storage.WriteInt32(PersistentConstants.HEADER_OFFSET_NODE_COUNT, _nodeCount);
        _storage.WriteInt64(PersistentConstants.HEADER_OFFSET_SIZE, _storage.Size);

        if (IsHybrid)
        {
            _storage.WriteInt64(PersistentConstants.HEADER_OFFSET_TRANSITION, _transitionOffset);
            _storage.WriteInt64(PersistentConstants.HEADER_OFFSET_JUMP_START, _jumpTableStart);
            _storage.WriteInt64(PersistentConstants.HEADER_OFFSET_JUMP_END, _jumpTableEnd);
        }

        // P17: Persist deepest internal node offset for O(1) LRS on Load
        _storage.WriteInt64(PersistentConstants.HEADER_OFFSET_DEEPEST_NODE, _deepestInternalNodeOffset);
    }

    // ──────────────── Jump table ────────────────

    /// <summary>
    /// Pre-calculates and allocates a contiguous jump table for ALL cross-zone references.
    /// Must be called BEFORE WriteChildrenArrays so all jump slots are available.
    /// </summary>
    private void MaterializeJumpTable()
    {
        // Count suffix-link jump entries
        int suffixLinkJumps = _deferredSuffixLinks.Count;

        // Count child-array jump entries: compact parents with any child in large zone
        int childArrayJumps = 0;
        var compactParentsNeedingJump = new List<long>();
        foreach (var (nodeOffset, childList) in _children)
        {
            var parentLayout = LayoutOf(nodeOffset);
            if (!parentLayout.OffsetIs64Bit)
            {
                foreach (var (_, childOffset) in childList)
                {
                    if (childOffset >= _transitionOffset)
                    {
                        childArrayJumps++;
                        compactParentsNeedingJump.Add(nodeOffset);
                        break;
                    }
                }
            }
        }

        int totalEntries = suffixLinkJumps + childArrayJumps;
        if (totalEntries == 0)
        {
            _jumpTableStart = -1;
            _jumpTableEnd = -1;
            return;
        }

        // Allocate one contiguous block for all entries
        int tableSize = totalEntries * 8;
        _jumpTableStart = _storage.Allocate(tableSize);
        _jumpTableEnd = _jumpTableStart + tableSize;

        // Fill suffix-link entries first
        int slotIndex = 0;
        for (int i = 0; i < suffixLinkJumps; i++, slotIndex++)
        {
            var (compactNodeOffset, largeTargetOffset) = _deferredSuffixLinks[i];
            long jumpEntryOffset = _jumpTableStart + (long)slotIndex * 8;
            _jumpEntries.Add(jumpEntryOffset);

            _storage.WriteInt64(jumpEntryOffset, largeTargetOffset);

            // Point the compact node's SuffixLink to the jump entry
            var compactLayout = LayoutOf(compactNodeOffset);
            var node = new PersistentSuffixTreeNode(_storage, compactNodeOffset, compactLayout);
            node.SuffixLink = jumpEntryOffset;
        }

        // Pre-allocate slots for child-array jumps (filled later in WriteChildrenArrays)
        _childArrayJumpSlots = new Dictionary<long, long>(childArrayJumps);
        for (int i = 0; i < compactParentsNeedingJump.Count; i++, slotIndex++)
        {
            long jumpEntryOffset = _jumpTableStart + (long)slotIndex * 8;
            _jumpEntries.Add(jumpEntryOffset);
            _childArrayJumpSlots[compactParentsNeedingJump[i]] = jumpEntryOffset;
        }
    }

    // Map from compact parent offset → pre-allocated jump entry offset for child arrays
    private Dictionary<long, long>? _childArrayJumpSlots;

    // ──────────────── Leaf count ────────────────

    private void CalculateLeafCount(long rootOffset)
    {
        var workStack = new Stack<long>();
        var resultStack = new Stack<long>();

        workStack.Push(rootOffset);

        while (workStack.Count > 0)
        {
            long offset = workStack.Pop();
            resultStack.Push(offset);

            if (_children.TryGetValue(offset, out var childList))
            {
                foreach (var (_, childOffset) in childList)
                    workStack.Push(childOffset);
            }
        }

        // Track deepest internal node during the bottom-up pass
        long deepestOffset = rootOffset;
        int maxDepth = 0;

        while (resultStack.Count > 0)
        {
            long offset = resultStack.Pop();
            var layout = LayoutOf(offset);
            var node = new PersistentSuffixTreeNode(_storage, offset, layout);

            if (node.IsLeaf)
            {
                node.LeafCount = 1;
            }
            else
            {
                uint totalLeaves = 0;
                if (_children.TryGetValue(offset, out var childList))
                {
                    foreach (var (_, childOffset) in childList)
                    {
                        var childLayout = LayoutOf(childOffset);
                        var child = new PersistentSuffixTreeNode(_storage, childOffset, childLayout);
                        totalLeaves += child.LeafCount;
                    }
                }
                node.LeafCount = totalLeaves;

                // Track deepest internal node (by DepthFromRoot + edge length)
                int nodeDepth = (int)node.DepthFromRoot + LengthOf(node);
                if (nodeDepth > maxDepth)
                {
                    maxDepth = nodeDepth;
                    deepestOffset = offset;
                }
            }
        }

        _deepestInternalNodeOffset = deepestOffset;
    }

    // ──────────────── Builder child management (in-memory) ────────────────

    private bool BuilderTryGetChild(long nodeOffset, uint key, out long childOffset)
    {
        if (_children.TryGetValue(nodeOffset, out var childList))
        {
            foreach (var entry in childList)
            {
                if (entry.Key == key)
                {
                    childOffset = entry.ChildOffset;
                    return true;
                }
            }
        }
        childOffset = PersistentConstants.NULL_OFFSET;
        return false;
    }

    private void BuilderSetChild(long nodeOffset, uint key, long childOffset)
    {
        if (!_children.TryGetValue(nodeOffset, out var childList))
        {
            childList = new List<(uint Key, long ChildOffset)>();
            _children[nodeOffset] = childList;
        }

        for (int i = 0; i < childList.Count; i++)
        {
            if (childList[i].Key == key)
            {
                childList[i] = (key, childOffset);
                return;
            }
        }

        childList.Add((key, childOffset));
    }

    // ──────────────── Write children arrays ────────────────

    private void WriteChildrenArrays()
    {
        foreach (var (nodeOffset, childList) in _children)
        {
            // Sort by key using signed comparison (terminator=-1 first)
            childList.Sort((a, b) => ((int)a.Key).CompareTo((int)b.Key));

            var parentLayout = LayoutOf(nodeOffset);
            bool parentIsCompact = !parentLayout.OffsetIs64Bit;

            // Check if this compact parent has a pre-allocated jump slot
            bool hasJumpSlot = parentIsCompact && _childArrayJumpSlots != null
                && _childArrayJumpSlots.ContainsKey(nodeOffset);

            // If a compact parent has large-zone children, the child array must use
            // Large entry format (12 bytes) to hold int64 offsets.
            NodeLayout arrayLayout;
            if (hasJumpSlot)
                arrayLayout = NodeLayout.Large; // 12-byte entries for int64 child offsets
            else
                arrayLayout = parentLayout;

            int count = childList.Count;
            int entrySize = arrayLayout.ChildEntrySize;
            int totalBytes = checked(count * entrySize);
            long arrayOffset = _storage.Allocate(totalBytes);

            // Serialize all entries into a single buffer, then batch-write
            byte[] buf = ArrayPool<byte>.Shared.Rent(totalBytes);
            try
            {
                for (int i = 0; i < count; i++)
                {
                    int off = checked(i * entrySize);
                    BinaryPrimitives.WriteUInt32LittleEndian(buf.AsSpan(off, 4), childList[i].Key);
                    if (arrayLayout.OffsetIs64Bit)
                        BinaryPrimitives.WriteInt64LittleEndian(buf.AsSpan(off + 4, 8), childList[i].ChildOffset);
                    else
                        BinaryPrimitives.WriteUInt32LittleEndian(buf.AsSpan(off + 4, 4), (uint)childList[i].ChildOffset);
                }
                _storage.WriteBytes(arrayOffset, buf, 0, totalBytes);
            }
            finally
            {
                ArrayPool<byte>.Shared.Return(buf);
            }

            var node = new PersistentSuffixTreeNode(_storage, nodeOffset, parentLayout);

            if (hasJumpSlot)
            {
                // Write the child array offset into the pre-allocated jump entry
                long jumpOffset = _childArrayJumpSlots![nodeOffset];
                _storage.WriteInt64(jumpOffset, arrayOffset);
                node.ChildrenHead = jumpOffset;
                // Store child count with the high bit set to signal "jump + large entries"
                node.ChildCount = count | unchecked((int)0x80000000);
            }
            else
            {
                node.ChildrenHead = arrayOffset;
                node.ChildCount = count;
            }
        }
    }
}
