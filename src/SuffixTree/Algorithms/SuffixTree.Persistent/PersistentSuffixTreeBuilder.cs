using System.Buffers;
using System.Buffers.Binary;
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
/// <para>
/// <b>Off-heap children</b>: parent→child edges are stored in a separate
/// <see cref="IStorageProvider"/> (<c>_childStore</c>) as linked lists instead
/// of a managed <c>Dictionary</c>. This keeps the managed heap near-zero
/// during genome-scale builds where the dictionary would exceed 20 GB.
/// </para>
/// </summary>
#pragma warning disable CA1001 // _childStore disposed in FinalizeTree (owned) or by caller (external)
public class PersistentSuffixTreeBuilder
#pragma warning restore CA1001
{
    private readonly IStorageProvider _storage;
    private readonly MappedFileStorageProvider? _mmfStorage;   // non-null when storage is MMF — enables unchecked fast-path
    private readonly MappedFileStorageProvider? _mmfChildStore; // non-null when child storage is MMF
    private readonly NodeLayout _initialLayout;  // layout before any transition (always the compact variant)
    private NodeLayout _layout;                // switches Compact → Large on overflow
    private long _compactOffsetLimit = NodeLayout.CompactMaxOffset;
    private readonly long _rootOffset;
    private ITextSource _text = new StringTextSource(string.Empty);
    private string? _rawString;
    private int _nodeCount;

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
    private bool _built;

    // Off-heap child storage: linked lists in a separate IStorageProvider.
    // Each entry: [Key:uint4, NextIndex:int4, ChildOffset:long8] = 16 bytes.
    // Node.ChildrenHead stores the entry INDEX (not byte offset) during build.
    private readonly IStorageProvider _childStore;
    private readonly bool _ownsChildStore;
    private const int CHILD_ENTRY_SIZE = 16;
    private const int CE_OFF_KEY = 0;
    private const int CE_OFF_NEXT = 4;
    private const int CE_OFF_CHILD = 8;

    // Tracks the end of the last allocated node (used by MaterializeJumpTable
    // to iterate compact-zone nodes without needing _children.Keys).
    private long _maxNodeEndOffset;

    // Deferred cross-zone suffix links: compact-node → large-zone target
    // Written at FinalizeTree when the jump table is materialized.
    private readonly List<(long CompactNodeOffset, long LargeTargetOffset)> _deferredSuffixLinks = new();

    // Cross-zone suffix link resolution for builder: compact-node offset →
    // correct large-zone target (avoids uint32 truncation when reading back
    // from storage during the Ukkonen walk).  Only populated when transition occurs.
    private readonly Dictionary<long, long> _crossZoneSuffixLinks = new();

    // Pre-computed during finalization DFS — the internal node with maximum depth.
    private long _deepestInternalNodeOffset = PersistentConstants.NULL_OFFSET;

    // Off-heap depth storage for v6 Slim layouts (no DepthFromRoot in node).
    // Indexed by sequential node number: depthOffset = NodeIndex(nodeOffset) * 4.
    private readonly IStorageProvider? _depthStore;
    private readonly MappedFileStorageProvider? _mmfDepthStore; // unchecked fast-path
    private readonly bool _ownsDepthStore;

    // Pre-computed LRS total depth (v6 Slim): DepthFromRoot + EdgeLen of deepest internal node.
    private int _lrsDepth;

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

    /// <summary>Number of cross-zone suffix links deferred to the jump table during build.</summary>
    internal int CrossZoneSuffixLinkCount { get; private set; }

    /// <summary>
    /// Offset of the deepest internal node (maximum DepthFromRoot + edge length).
    /// Computed during <see cref="Build"/> at zero extra cost (piggybacks on leaf-count traversal).
    /// Returns <see cref="PersistentConstants.NULL_OFFSET"/> if Build has not been called.
    /// </summary>
    internal long DeepestInternalNodeOffset => _deepestInternalNodeOffset;

    /// <summary>
    /// Pre-computed LRS total depth for v6 Slim layouts. 0 if not slim.
    /// </summary>
    internal int LrsDepth => _lrsDepth;

    /// <summary>
    /// Header size used by this builder's layout.
    /// </summary>
    private int HeaderSize => PersistentConstants.HEADER_SIZE_V6;

    /// <summary>Initializes a new builder with the specified storage provider and optional layout.</summary>
    /// <param name="storage">Main tree storage (the MMF or heap that holds the final tree).</param>
    /// <param name="layout">Initial node layout (default: Compact).</param>
    /// <param name="childStorage">
    /// Optional separate storage for build-time child linked lists.
    /// For genome-scale MMF builds, pass a <see cref="MappedFileStorageProvider"/> backed
    /// by a temp file so that the ~8 GB of child-edge data stays off the managed heap.
    /// When <c>null</c>, an internal <see cref="HeapStorageProvider"/> is created (fine for small trees and tests).
    /// </param>
    /// <param name="depthStorage">
    /// Optional separate storage for build-time DepthFromRoot values.
    /// When <c>null</c>, an internal <see cref="HeapStorageProvider"/> is created.
    /// </param>
    public PersistentSuffixTreeBuilder(IStorageProvider storage, NodeLayout? layout = null,
                                        IStorageProvider? childStorage = null,
                                        IStorageProvider? depthStorage = null)
    {
        _storage = storage;
        _mmfStorage = storage as MappedFileStorageProvider;
        _layout = layout ?? NodeLayout.Compact;
        _initialLayout = _layout;

        if (childStorage != null)
        {
            _childStore = childStorage;
            _ownsChildStore = false;
        }
        else
        {
            _childStore = new HeapStorageProvider();
            _ownsChildStore = true;
        }
        _mmfChildStore = _childStore as MappedFileStorageProvider;

        // Off-heap depth storage (DepthFromRoot not stored in node)
        if (depthStorage != null)
        {
            _depthStore = depthStorage;
            _ownsDepthStore = false;
        }
        else
        {
            _depthStore = new HeapStorageProvider();
            _ownsDepthStore = true;
        }
        _mmfDepthStore = _depthStore as MappedFileStorageProvider;

        // Allocate header (88 bytes).
        _storage.Allocate(HeaderSize);

        _rootOffset = _storage.Allocate(_layout.NodeSize);
        _maxNodeEndOffset = _rootOffset + _layout.NodeSize;
        _nodeCount = 1;
        var root = new PersistentSuffixTreeNode(_storage, _rootOffset, _layout);
        root.Start = 0;
        root.End = 0;
        root.SuffixLink = PersistentConstants.NULL_OFFSET;
        root.ChildrenHead = PersistentConstants.NULL_OFFSET;

        // Root depth = 0 (default from zero-init for both node field and depth store)
        _depthStore.Allocate(4); // index 0 for root, value = 0

        _activeNodeOffset = _rootOffset;
    }

    /// <summary>Builds the suffix tree from the given text source using Ukkonen's algorithm.</summary>
    public long Build(ITextSource text)
    {
        ArgumentNullException.ThrowIfNull(text);
        if (_built)
            throw new InvalidOperationException("Build() has already been called. Create a new builder instance to build another tree.");
        _built = true;

        _text = text;
        _rawString = (text as StringTextSource)?.Value;
        if (text.Length > 0)
        {
            // Direct string indexing eliminates ITextSource virtual dispatch
            var raw = _rawString;
            if (raw != null)
            {
                for (int i = 0; i < raw.Length; i++)
                    ExtendTree((uint)raw[i]);
            }
            else
            {
                for (int i = 0; i < text.Length; i++)
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

            uint activeEdgeKey = GetSymbolAt(_activeEdgeIndex);

            // Fast-path: read active node depth via unchecked compact MMF
            uint activeNodeDepth;
            if (IsCompactMmf(_activeNodeOffset))
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
                // Read nextChild Start/End via unchecked path if possible
                uint ncStart, ncEnd;
                if (IsCompactMmf(nextChildOffset))
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

                // Read split node's edge length for depth calc
                uint splitStart, splitEnd;
                if (IsCompactMmf(splitOffset))
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

                // Update the original child's Start
                uint newStart = ncStart + (uint)_activeLength;
                if (IsCompactMmf(nextChildOffset))
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
            // Still in compact mode — check if this allocation would overflow
            long currentSize = _storage.Size;
            if ((currentSize + size) > _compactOffsetLimit)
            {
                // Transition to Large layout
                _transitionOffset = currentSize;
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
            int compactCount = (int)((_transitionOffset - hdr) / _initialLayout.NodeSize);
            return compactCount + (int)((nodeOffset - _transitionOffset) / NodeLayout.Large.NodeSize);
        }
        return (int)((nodeOffset - hdr) / _initialLayout.NodeSize);
    }

    /// <summary>Reads DepthFromRoot for a node from off-heap depth store.</summary>
    [System.Runtime.CompilerServices.MethodImpl(System.Runtime.CompilerServices.MethodImplOptions.AggressiveInlining)]
    private uint GetBuildDepth(long nodeOffset)
    {
        long off = (long)NodeIndex(nodeOffset) * 4;
        if (_mmfDepthStore != null)
            return _mmfDepthStore.ReadUInt32Unchecked(off);
        return _depthStore!.ReadUInt32(off);
    }

    /// <summary>Writes DepthFromRoot for a node to off-heap depth store.</summary>
    [System.Runtime.CompilerServices.MethodImpl(System.Runtime.CompilerServices.MethodImplOptions.AggressiveInlining)]
    private void SetBuildDepth(long nodeOffset, uint value)
    {
        int idx = NodeIndex(nodeOffset);
        long off = (long)idx * 4;
        while (_depthStore!.Size < off + 4)
            _depthStore.Allocate(4);
        if (_mmfDepthStore != null)
            _mmfDepthStore.WriteUInt32Unchecked(off, value);
        else
            _depthStore.WriteUInt32(off, value);
    }

    // ──────────────── Finalize ────────────────

    private void FinalizeTree()
    {
        // Materialize a contiguous jump table for ALL cross-zone references
        // BEFORE the combined DFS so child-array jump slots are available.
        if (IsHybrid)
            MaterializeJumpTable();

        // Release depth storage before the DFS — depth is tracked inline
        // in the DFS stack, so we no longer need the off-heap depth data.
        // This frees the temp MMF file early for genome-scale builds.
        if (_ownsDepthStore)
            _depthStore?.Dispose();

        // Single post-order DFS: compute leaf counts, find deepest internal node,
        // and write sorted children arrays — all in one traversal.
        CalculateLeafCountAndWriteChildrenArrays(_rootOffset);

        // Release temp child storage — no longer needed after child arrays are written
        if (_ownsChildStore)
            _childStore.Dispose();

        // Store text in storage for true persistence (chunked write, no full-string copy).
        // Use ASCII (1 byte/char) when all characters fit in 0-127, else UTF-16 (2 bytes/char).
        bool textIsAscii = IsTextAscii();
        int bytesPerChar = textIsAscii ? 1 : 2;
        long textByteLen = (long)_text.Length * bytesPerChar;
        if (textByteLen > int.MaxValue)
            throw new InvalidOperationException(
                $"Text length {_text.Length} exceeds maximum serializable size.");
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
                int byteCount;
                if (textIsAscii)
                {
                    // Write 1 byte per char (ASCII fast path)
                    for (int i = 0; i < chunkLen; i++)
                        chunkBuf[i] = (byte)charSpan[i];
                    byteCount = chunkLen;
                }
                else
                {
                    byteCount = Encoding.Unicode.GetBytes(charSpan, chunkBuf.AsSpan());
                }
                _storage.WriteBytes(textOffset + (long)written * bytesPerChar, chunkBuf, 0, byteCount);
                written += chunkLen;
            }
        }
        finally
        {
            ArrayPool<byte>.Shared.Return(chunkBuf);
        }

        // Write v6 Header (88 bytes)
        _storage.WriteInt64(PersistentConstants.HEADER_OFFSET_MAGIC, PersistentConstants.MAGIC_NUMBER);
        _storage.WriteInt32(PersistentConstants.HEADER_OFFSET_VERSION, 6);
        _storage.WriteInt64(PersistentConstants.HEADER_OFFSET_ROOT, _rootOffset);
        _storage.WriteInt64(PersistentConstants.HEADER_OFFSET_TEXT_OFF, textOffset);
        _storage.WriteInt32(PersistentConstants.HEADER_OFFSET_TEXT_LEN, _text.Length);
        _storage.WriteInt32(PersistentConstants.HEADER_OFFSET_NODE_COUNT, _nodeCount);
        _storage.WriteInt32(PersistentConstants.HEADER_OFFSET_FLAGS,
            textIsAscii ? PersistentConstants.FLAG_TEXT_ASCII : 0);
        _storage.WriteInt64(PersistentConstants.HEADER_OFFSET_SIZE, _storage.Size);

        // Hybrid fields — set to -1 for non-hybrid trees
        _storage.WriteInt64(PersistentConstants.HEADER_OFFSET_TRANSITION,
            IsHybrid ? _transitionOffset : -1);
        _storage.WriteInt64(PersistentConstants.HEADER_OFFSET_JUMP_START,
            IsHybrid ? _jumpTableStart : -1);
        _storage.WriteInt64(PersistentConstants.HEADER_OFFSET_JUMP_END,
            IsHybrid ? _jumpTableEnd : -1);

        // Deepest internal node offset for O(1) LRS on Load
        _storage.WriteInt64(PersistentConstants.HEADER_OFFSET_DEEPEST_NODE, _deepestInternalNodeOffset);

        // LRS total depth
        _storage.WriteInt32(PersistentConstants.HEADER_OFFSET_LRS_DEPTH, _lrsDepth);

        // Base layout node size (for Load() to distinguish Compact vs Large base)
        _storage.WriteInt32(PersistentConstants.HEADER_OFFSET_BASE_NODE_SIZE, _initialLayout.NodeSize);

        // Note: _depthStore was already disposed before the DFS (depth tracked inline).

        // S20: Release in-memory build collections — the tree is fully serialized,
        // keeping these alive wastes memory if the builder instance is held.
        CrossZoneSuffixLinkCount = _deferredSuffixLinks.Count;
        _deferredSuffixLinks.Clear();
        _crossZoneSuffixLinks.Clear();
        _jumpEntries.Clear();
        _childArrayJumpSlots = null;
    }

    // ──────────────── Jump table ────────────────

    /// <summary>
    /// Pre-calculates and allocates a contiguous jump table for ALL cross-zone references.
    /// Must be called BEFORE the combined leaf-count/child-array DFS so all jump slots are available.
    /// </summary>
    private void MaterializeJumpTable()
    {
        // Count suffix-link jump entries
        int suffixLinkJumps = _deferredSuffixLinks.Count;

        // Count child-array jump entries: compact parents with any child in large zone
        int childArrayJumps = 0;
        var compactParentsNeedingJump = new List<long>();

        long compactEnd = Math.Min(_transitionOffset, _maxNodeEndOffset);
        for (long offset = HeaderSize; offset < compactEnd; offset += _initialLayout.NodeSize)
        {
            if (ChildListHasLargeZoneTarget(offset))
            {
                childArrayJumps++;
                compactParentsNeedingJump.Add(offset);
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

        // Pre-allocate slots for child-array jumps (filled later in the combined DFS)
        _childArrayJumpSlots = new Dictionary<long, long>(childArrayJumps);
        for (int i = 0; i < compactParentsNeedingJump.Count; i++, slotIndex++)
        {
            long jumpEntryOffset = _jumpTableStart + (long)slotIndex * 8;
            _jumpEntries.Add(jumpEntryOffset);
            _childArrayJumpSlots[compactParentsNeedingJump[i]] = jumpEntryOffset;
        }
    }

    /// <summary>Checks whether any child of the node at <paramref name="nodeOffset"/> is in the large zone.</summary>
    private bool ChildListHasLargeZoneTarget(long nodeOffset)
    {
        var node = new PersistentSuffixTreeNode(_storage, nodeOffset, LayoutOf(nodeOffset));
        long headIndex = node.ChildrenHead;
        if (headIndex == PersistentConstants.NULL_OFFSET) return false;

        int idx = (int)headIndex;
        while (idx >= 0)
        {
            long entryOff = (long)idx * CHILD_ENTRY_SIZE;
            long childOffset = _childStore.ReadInt64(entryOff + CE_OFF_CHILD);
            if (childOffset >= _transitionOffset)
                return true;
            idx = _childStore.ReadInt32(entryOff + CE_OFF_NEXT);
        }
        return false;
    }

    // Map from compact parent offset → pre-allocated jump entry offset for child arrays
    private Dictionary<long, long>? _childArrayJumpSlots;

    // ──────────────── Combined leaf-count + child-array pass ────────────────

    /// <summary>
    /// Single post-order DFS that (a) computes LeafCount for every node,
    /// (b) finds the deepest internal node, and (c) collects, sorts, and
    /// writes the final child arrays — eliminating the previous two-pass
    /// approach (CalculateLeafCount + WriteChildrenArrays).
    /// <para>
    /// Depth is tracked inline in the DFS stack (no _depthStore reads during
    /// finalization), which avoids O(N) random reads from a temp MMF.
    /// Peak stack size is O(branching_factor × tree_depth), typically a few MB
    /// for genome-scale data (DNA alphabet 5, depth ~50K → ~250K entries).
    /// </para>
    /// </summary>
    private void CalculateLeafCountAndWriteChildrenArrays(long rootOffset)
    {
        // Stack entry: (nodeOffset, childrenPushed, depthFromRoot).
        // depthFromRoot = parent's depthFromRoot + parent's edge length.
        var stack = new Stack<(long Offset, bool ChildrenPushed, int DepthFromRoot)>();
        stack.Push((rootOffset, false, 0));

        long deepestOffset = rootOffset;
        int maxDepth = 0;

        // Reusable buffer for collecting children (cleared each use)
        var childBuf = new List<(uint Key, long ChildOffset)>();

        // Cache MMF references to avoid repeated null-checks in inner loop
        var mmfMain = _mmfStorage;
        var mmfChild = _mmfChildStore;

        while (stack.Count > 0)
        {
            var (offset, childrenPushed, depthFromRoot) = stack.Pop();
            var layout = LayoutOf(offset);
            bool isCompact = mmfMain != null && !layout.OffsetIs64Bit;

            // Read End to check IsLeaf
            uint nodeEnd = isCompact
                ? mmfMain!.ReadUInt32Unchecked(offset + 4)
                : new PersistentSuffixTreeNode(_storage, offset, layout).End;

            if (nodeEnd == PersistentConstants.BOUNDLESS) // IsLeaf
            {
                // Write LeafCount = 1
                if (isCompact)
                    mmfMain!.WriteUInt32Unchecked(offset + 12, 1); // LeafCount offset=12 in Compact
                else
                    new PersistentSuffixTreeNode(_storage, offset, layout).LeafCount = 1;
                continue;
            }

            if (childrenPushed)
            {
                // ── Post-order: all children already processed ──

                // 1. Collect children from _childStore and sum leaf counts in one pass
                childBuf.Clear();
                uint totalLeaves = 0;

                long headIndex;
                if (isCompact)
                {
                    uint raw = mmfMain!.ReadUInt32Unchecked(offset + 16);
                    headIndex = raw == uint.MaxValue ? PersistentConstants.NULL_OFFSET : (long)raw;
                }
                else
                {
                    headIndex = new PersistentSuffixTreeNode(_storage, offset, layout).ChildrenHead;
                }

                if (headIndex != PersistentConstants.NULL_OFFSET)
                {
                    int idx = (int)headIndex;
                    if (mmfChild != null)
                    {
                        while (idx >= 0)
                        {
                            long entryOff = (long)idx * CHILD_ENTRY_SIZE;
                            uint key = mmfChild.ReadUInt32Unchecked(entryOff + CE_OFF_KEY);
                            long childOffset = mmfChild.ReadInt64Unchecked(entryOff + CE_OFF_CHILD);

                            // Read child's LeafCount
                            var childLayout = LayoutOf(childOffset);
                            bool childCompact = mmfMain != null && !childLayout.OffsetIs64Bit;
                            uint childLeaves = childCompact
                                ? mmfMain!.ReadUInt32Unchecked(childOffset + 12)
                                : new PersistentSuffixTreeNode(_storage, childOffset, childLayout).LeafCount;
                            totalLeaves += childLeaves;

                            childBuf.Add((key, childOffset));
                            idx = mmfChild.ReadInt32Unchecked(entryOff + CE_OFF_NEXT);
                        }
                    }
                    else
                    {
                        while (idx >= 0)
                        {
                            long entryOff = (long)idx * CHILD_ENTRY_SIZE;
                            uint key = _childStore.ReadUInt32(entryOff + CE_OFF_KEY);
                            long childOffset = _childStore.ReadInt64(entryOff + CE_OFF_CHILD);

                            var childLayout = LayoutOf(childOffset);
                            var child = new PersistentSuffixTreeNode(_storage, childOffset, childLayout);
                            totalLeaves += child.LeafCount;

                            childBuf.Add((key, childOffset));
                            idx = _childStore.ReadInt32(entryOff + CE_OFF_NEXT);
                        }
                    }
                }

                // Write LeafCount
                if (isCompact)
                    mmfMain!.WriteUInt32Unchecked(offset + 12, totalLeaves);
                else
                    new PersistentSuffixTreeNode(_storage, offset, layout).LeafCount = totalLeaves;

                // 2. Track deepest internal node (by DepthFromRoot + edge length)
                uint nodeStart = isCompact
                    ? mmfMain!.ReadUInt32Unchecked(offset + 0)
                    : new PersistentSuffixTreeNode(_storage, offset, layout).Start;
                int edgeLenHere = (int)((nodeEnd == PersistentConstants.BOUNDLESS ? (uint)(_position + 1) : nodeEnd) - nodeStart);
                int nodeDepth = depthFromRoot + edgeLenHere;
                if (nodeDepth > maxDepth)
                {
                    maxDepth = nodeDepth;
                    deepestOffset = offset;
                }

                // 3. Sort by key using signed comparison (terminator=-1 first)
                childBuf.Sort((a, b) => ((int)a.Key).CompareTo((int)b.Key));

                // 4. Write sorted child array to main storage (skip if node has no children, e.g. empty-text root)
                int count = childBuf.Count;
                if (count > 0)
                {
                    bool parentIsCompact = !layout.OffsetIs64Bit;
                    bool hasJumpSlot = parentIsCompact && _childArrayJumpSlots != null
                        && _childArrayJumpSlots.ContainsKey(offset);

                    NodeLayout arrayLayout = hasJumpSlot ? NodeLayout.Large : layout;

                    int entrySize = arrayLayout.ChildEntrySize;
                    int totalBytes = checked(count * entrySize);
                    long arrayOffset = _storage.Allocate(totalBytes);

                    byte[] buf = ArrayPool<byte>.Shared.Rent(totalBytes);
                    try
                    {
                        for (int i = 0; i < count; i++)
                        {
                            int off = checked(i * entrySize);
                            BinaryPrimitives.WriteUInt32LittleEndian(buf.AsSpan(off, 4), childBuf[i].Key);
                            if (arrayLayout.OffsetIs64Bit)
                                BinaryPrimitives.WriteInt64LittleEndian(buf.AsSpan(off + 4, 8), childBuf[i].ChildOffset);
                            else
                                BinaryPrimitives.WriteUInt32LittleEndian(buf.AsSpan(off + 4, 4), (uint)childBuf[i].ChildOffset);
                        }
                        _storage.WriteBytes(arrayOffset, buf, 0, totalBytes);
                    }
                    finally
                    {
                        ArrayPool<byte>.Shared.Return(buf);
                    }

                    if (hasJumpSlot)
                    {
                        long jumpOffset = _childArrayJumpSlots![offset];
                        _storage.WriteInt64(jumpOffset, arrayOffset);
                        if (isCompact)
                        {
                            mmfMain!.WriteUInt32Unchecked(offset + 16, (uint)jumpOffset); // ChildrenHead
                            mmfMain!.WriteInt32Unchecked(offset + 20, count | unchecked((int)0x80000000)); // ChildCount
                        }
                        else
                        {
                            var node2 = new PersistentSuffixTreeNode(_storage, offset, layout);
                            node2.ChildrenHead = jumpOffset;
                            node2.ChildCount = count | unchecked((int)0x80000000);
                        }
                    }
                    else
                    {
                        if (isCompact)
                        {
                            mmfMain!.WriteUInt32Unchecked(offset + 16, (uint)arrayOffset); // ChildrenHead
                            mmfMain!.WriteInt32Unchecked(offset + 20, count);               // ChildCount
                        }
                        else
                        {
                            var node2 = new PersistentSuffixTreeNode(_storage, offset, layout);
                            node2.ChildrenHead = arrayOffset;
                            node2.ChildCount = count;
                        }
                    }
                }
            }
            else
            {
                // First visit: re-push self as processed, then push all children
                uint nodeStart = isCompact
                    ? mmfMain!.ReadUInt32Unchecked(offset + 0)
                    : new PersistentSuffixTreeNode(_storage, offset, layout).Start;
                int edgeLen = (int)((nodeEnd == PersistentConstants.BOUNDLESS ? (uint)(_position + 1) : nodeEnd) - nodeStart);
                int childDepth = depthFromRoot + edgeLen;
                stack.Push((offset, true, depthFromRoot));

                long headIndex;
                if (isCompact)
                {
                    uint raw = mmfMain!.ReadUInt32Unchecked(offset + 16);
                    headIndex = raw == uint.MaxValue ? PersistentConstants.NULL_OFFSET : (long)raw;
                }
                else
                {
                    headIndex = new PersistentSuffixTreeNode(_storage, offset, layout).ChildrenHead;
                }

                if (headIndex != PersistentConstants.NULL_OFFSET)
                {
                    int idx = (int)headIndex;
                    if (mmfChild != null)
                    {
                        while (idx >= 0)
                        {
                            long entryOff = (long)idx * CHILD_ENTRY_SIZE;
                            long childOffset = mmfChild.ReadInt64Unchecked(entryOff + CE_OFF_CHILD);
                            stack.Push((childOffset, false, childDepth));
                            idx = mmfChild.ReadInt32Unchecked(entryOff + CE_OFF_NEXT);
                        }
                    }
                    else
                    {
                        while (idx >= 0)
                        {
                            long entryOff = (long)idx * CHILD_ENTRY_SIZE;
                            long childOffset = _childStore.ReadInt64(entryOff + CE_OFF_CHILD);
                            stack.Push((childOffset, false, childDepth));
                            idx = _childStore.ReadInt32(entryOff + CE_OFF_NEXT);
                        }
                    }
                }
            }
        }

        _deepestInternalNodeOffset = deepestOffset;
        _lrsDepth = maxDepth;
    }

    // ──────────────── Builder child management (off-heap linked lists) ────────────────

    private bool BuilderTryGetChild(long nodeOffset, uint key, out long childOffset)
    {
        // Read ChildrenHead — fast-path for compact MMF
        long headIndex;
        if (_mmfStorage != null && (_transitionOffset < 0 || nodeOffset < _transitionOffset))
        {
            // Compact: ChildrenHead is uint32 at nodeOffset + 16
            uint raw = _mmfStorage.ReadUInt32Unchecked(nodeOffset + 16);
            headIndex = raw == uint.MaxValue ? PersistentConstants.NULL_OFFSET : (long)raw;
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
        // Read ChildrenHead — fast-path for compact MMF
        bool isCompactMmf = _mmfStorage != null && (_transitionOffset < 0 || nodeOffset < _transitionOffset);
        long headIndex;
        if (isCompactMmf)
        {
            uint raw = _mmfStorage!.ReadUInt32Unchecked(nodeOffset + 16);
            headIndex = raw == uint.MaxValue ? PersistentConstants.NULL_OFFSET : (long)raw;
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
            _mmfStorage!.WriteUInt32Unchecked(nodeOffset + 16, (uint)newIndex);
        }
        else
        {
            var nodeLayout = LayoutOf(nodeOffset);
            var node = new PersistentSuffixTreeNode(_storage, nodeOffset, nodeLayout);
            node.ChildrenHead = newIndex;
        }
    }

}
