using System.Buffers;
using System.Text;

namespace SuffixTree.Persistent;

/// <summary>
/// Handles construction of a persistent suffix tree using Ukkonen's algorithm.
/// Writes nodes and child entries directly to the storage provider.
/// <para>
/// <b>Hybrid continuation</b>: construction starts with the Compact v6 layout
/// (24-byte nodes, uint32 offsets). When the compact address space (~4 GB) is exhausted,
/// the builder seamlessly switches to the Large v6 layout (32-byte nodes, int64 offsets) and continues
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
    // Jump offsets for child arrays are stored directly in each compact internal
    // node's LeafCount field during MaterializeJumpTable (zero extra RAM).

    // Jump table reserve: block pre-allocated within the compact zone at transition time.
    // This ensures jump entry offsets fit in uint32 fields of compact-zone nodes.
    private long _compactNodesEnd = -1;          // offset after last compact node (before jump reserve)
    private long _jumpTableReserveStart = -1;    // start of pre-reserved jump table block
    private long _jumpTableReserveEnd = -1;      // end of pre-reserved jump table block
    private const int JUMP_TABLE_MIN_RESERVE = 16 * 1024 * 1024; // 16 MB floor for small trees
    private const int JUMP_TABLE_RESERVE_FRACTION = 4;            // reserve 1/4 (25%) of compact zone
    // NOTE: 25% is needed because ALL compact internal nodes require jump slots —
    // child arrays are allocated in the large zone (>4 GB) during the DFS phase,
    // and their offsets cannot fit in the compact uint32 ChildrenHead field.

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

    // Temp file entry layout for sequential finalize (per node, 12 bytes).
    private const int TEMP_ENTRY = 12;
    private const int TEMP_OFF_PARENT = 0;     // int32: parent node index (-1 = root)
    private const int TEMP_OFF_LEAFCOUNT = 4;  // uint32: leaf count
    private const int TEMP_OFF_REMAINING = 8;  // int32: remaining children for bottom-up

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

    // Progress reporting
    private IProgress<(string Stage, double Percent)>? _progress;

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
    /// <param name="text">The text source to index.</param>
    /// <param name="progress">Optional progress reporter: (Stage, Percent 0–100).</param>
    public long Build(ITextSource text, IProgress<(string Stage, double Percent)>? progress = null)
    {
        ArgumentNullException.ThrowIfNull(text);
        if (_built)
            throw new InvalidOperationException("Build() has already been called. Create a new builder instance to build another tree.");
        _built = true;

        _text = text;
        _rawString = (text as StringTextSource)?.Value;
        _progress = progress;
        if (text.Length > 0)
        {
            int total = text.Length;
            int reportInterval = Math.Max(total / 200, 1); // ~0.5% granularity
            // Direct string indexing eliminates ITextSource virtual dispatch
            var raw = _rawString;
            if (raw != null)
            {
                for (int i = 0; i < raw.Length; i++)
                {
                    ExtendTree((uint)raw[i]);
                    if (i % reportInterval == 0)
                        progress?.Report(("Ukkonen build", (double)i / total * 90.0));
                }
            }
            else
            {
                for (int i = 0; i < text.Length; i++)
                {
                    ExtendTree((uint)text[i]);
                    if (i % reportInterval == 0)
                        progress?.Report(("Ukkonen build", (double)i / total * 90.0));
                }
            }
            ExtendTree(PersistentConstants.TERMINATOR_KEY);
            progress?.Report(("Ukkonen build", 90.0));
        }

        progress?.Report(("Finalizing tree", 91.0));
        FinalizeTree();
        progress?.Report(("Build complete", 100.0));
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
        if (IsHybrid)
        {
            _progress?.Report(("Jump table materialization", 91.0));
            MaterializeJumpTable();
        }

        // Sequential multi-pass finalize: replaces DFS to eliminate random I/O.
        // Disposes child store and depth store internally when no longer needed.
        SequentialFinalize(_rootOffset);

        _progress?.Report(("DFS: Serializing text", 98.0));

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

        CrossZoneSuffixLinkCount = _deferredSuffixLinks.Count;
        _deferredSuffixLinks.Clear();
        _crossZoneSuffixLinks.Clear();
    }
    private void MaterializeJumpTable()
    {
        // Count suffix-link jump entries
        int suffixLinkJumps = _deferredSuffixLinks.Count;

        // Iterate only over actual compact nodes (exclude jump reserve block)
        long compactEnd = _compactNodesEnd >= 0
            ? _compactNodesEnd
            : Math.Min(_transitionOffset, _maxNodeEndOffset);
        var mmfMain = _mmfStorage;

        // ── Single-pass path when reserve is pre-allocated (genome-scale builds) ──
        // We know the reserve capacity and can write slots directly while scanning.
        if (_jumpTableReserveStart >= 0)
        {
            long reserveCapacity = _jumpTableReserveEnd - _jumpTableReserveStart;
            _jumpTableStart = _jumpTableReserveStart;

            // Fill suffix-link entries first
            int slotIndex = 0;
            for (int i = 0; i < suffixLinkJumps; i++, slotIndex++)
            {
                var (compactNodeOffset, largeTargetOffset) = _deferredSuffixLinks[i];
                long jumpEntryOffset = _jumpTableStart + (long)slotIndex * 8;
                _storage.WriteInt64(jumpEntryOffset, largeTargetOffset);

                var compactLayout = LayoutOf(compactNodeOffset);
                var node = new PersistentSuffixTreeNode(_storage, compactNodeOffset, compactLayout);
                node.SuffixLink = jumpEntryOffset;
            }

            // Single pass: write child-array jump slots while scanning compact nodes
            for (long offset = HeaderSize; offset < compactEnd; offset += _initialLayout.NodeSize)
            {
                uint nodeEnd = mmfMain != null
                    ? mmfMain.ReadUInt32Unchecked(offset + 4)
                    : _storage.ReadUInt32(offset + 4);
                if (nodeEnd != PersistentConstants.BOUNDLESS)
                {
                    long jumpEntryOffset = _jumpTableStart + (long)slotIndex * 8;
                    if (jumpEntryOffset + 8 > _jumpTableReserveEnd)
                        throw new InvalidOperationException(
                            $"Jump table requires more than {reserveCapacity:N0} reserved bytes. " +
                            "The input exceeds the jump table reserve capacity.");
                    if (mmfMain != null)
                        mmfMain.WriteUInt32Unchecked(offset + NodeLayout.Compact.OffsetLeafCount, (uint)jumpEntryOffset);
                    else
                        _storage.WriteUInt32(offset + NodeLayout.Compact.OffsetLeafCount, (uint)jumpEntryOffset);
                    slotIndex++;
                }
            }

            int totalEntries = slotIndex;
            if (totalEntries == 0)
            {
                _jumpTableStart = -1;
                _jumpTableEnd = -1;
            }
            else
            {
                _jumpTableEnd = _jumpTableStart + totalEntries * 8;
            }
            return;
        }

        // ── Fallback: two-pass for non-reserved (tiny compact limits in tests) ──

        // Pass 1: count compact internal nodes
        int childArrayJumps = 0;
        for (long offset = HeaderSize; offset < compactEnd; offset += _initialLayout.NodeSize)
        {
            uint nodeEnd = mmfMain != null
                ? mmfMain.ReadUInt32Unchecked(offset + 4)
                : _storage.ReadUInt32(offset + 4);
            if (nodeEnd != PersistentConstants.BOUNDLESS)
                childArrayJumps++;
        }

        int totalEntriesFallback = suffixLinkJumps + childArrayJumps;
        if (totalEntriesFallback == 0)
        {
            _jumpTableStart = -1;
            _jumpTableEnd = -1;
            return;
        }

        int tableSize = totalEntriesFallback * 8;
        _jumpTableStart = _storage.Allocate(tableSize);
        _jumpTableEnd = _jumpTableStart + tableSize;

        // Fill suffix-link entries
        int slotIdx = 0;
        for (int i = 0; i < suffixLinkJumps; i++, slotIdx++)
        {
            var (compactNodeOffset, largeTargetOffset) = _deferredSuffixLinks[i];
            long jumpEntryOffset = _jumpTableStart + (long)slotIdx * 8;
            _storage.WriteInt64(jumpEntryOffset, largeTargetOffset);

            var compactLayout = LayoutOf(compactNodeOffset);
            var node = new PersistentSuffixTreeNode(_storage, compactNodeOffset, compactLayout);
            node.SuffixLink = jumpEntryOffset;
        }

        // Pass 2: write child-array jump slots
        for (long offset = HeaderSize; offset < compactEnd; offset += _initialLayout.NodeSize)
        {
            uint nodeEnd = mmfMain != null
                ? mmfMain.ReadUInt32Unchecked(offset + 4)
                : _storage.ReadUInt32(offset + 4);
            if (nodeEnd != PersistentConstants.BOUNDLESS)
            {
                long jumpEntryOffset = _jumpTableStart + (long)slotIdx * 8;
                if (mmfMain != null)
                    mmfMain.WriteUInt32Unchecked(offset + NodeLayout.Compact.OffsetLeafCount, (uint)jumpEntryOffset);
                else
                    _storage.WriteUInt32(offset + NodeLayout.Compact.OffsetLeafCount, (uint)jumpEntryOffset);
                slotIdx++;
            }
        }
    }

    // ──────────────── Sequential multi-pass finalize ────────────────

    /// <summary>
    /// Replaces the old DFS-based finalization with sequential disk passes to
    /// eliminate random I/O that caused 13+ minute DFS phases on genome-scale builds.
    ///
    /// <b>Pass 1</b> — Sequential scan of main storage (19 GB):
    ///   For each node in allocation order: read child linked list from child store,
    ///   sort children, write sorted child array (batch-allocated), update
    ///   ChildrenHead/ChildCount in node header. Build parent-graph in temp file.
    ///
    /// <b>Pass 2</b> — Bottom-up leaf count propagation (temp files only, ~7 GB):
    ///   BFS from leaves → parents using temp file + sequential queue file.
    ///   Zero main storage I/O.
    ///
    /// <b>Pass 3</b> — Sequential write-back (19 GB):
    ///   Write LeafCount from temp file. Track deepest internal node via depth store.
    ///
    /// Disk: ~7 GB temp files (auto-deleted). RAM: near-zero (just OS page cache).
    /// </summary>
    private unsafe void SequentialFinalize(long rootOffset)
    {
        var mmfMain = _mmfStorage;
        byte* childBase = _mmfChildStore != null ? _mmfChildStore.RawPointer : null;
        int nodeCount = _nodeCount;
        int hdr = HeaderSize;

        string tempNodePath = Path.Combine(Path.GetTempPath(), $"sfx_{Guid.NewGuid():N}_node.tmp");
        string tempQueuePath = Path.Combine(Path.GetTempPath(), $"sfx_{Guid.NewGuid():N}_queue.tmp");

        try
        {
            // ── Create temp file for node graph metadata ──
            long tempNodeCap = (long)nodeCount * TEMP_ENTRY;
            using var tempNodeMmf = new MappedFileStorageProvider(tempNodePath, Math.Max(tempNodeCap, 65536));
            tempNodeMmf.SetSize(tempNodeCap);
            byte* tp = tempNodeMmf.RawPointer;

            // Initialize parent = -1 (leafCount=0 and remaining=0 are zero-filled by OS)
            for (long i = 0; i < nodeCount; i++)
                *(int*)(tp + i * TEMP_ENTRY + TEMP_OFF_PARENT) = -1;

            // ═══════════ PASS 1: Sequential scan → child arrays + parent graph ═══════════
            _progress?.Report(("Pass 1: child arrays", 92.0));

            // Prefetch child store into OS page cache before linked list walks.
            // This converts random page faults into sequential read-ahead.
            if (_mmfChildStore != null)
                _mmfChildStore.PrefetchForBuild();

            // Batch allocation for child arrays
            long estimatedChildArrayBytes = (long)nodeCount * 10;
            int batchSize = (int)Math.Min(estimatedChildArrayBytes, int.MaxValue - 65536);
            batchSize = Math.Max(batchSize, 4096);
            long batchStart = _storage.Allocate(batchSize);
            long batchEnd = batchStart + batchSize;
            long batchPos = batchStart;

            // Reusable child buffer (max 256 children per node)
            const int CHILD_BUF_MAX = 256;
            uint* cKeys = stackalloc uint[CHILD_BUF_MAX];
            long* cOffs = stackalloc long[CHILD_BUF_MAX];

            long compactEnd = IsHybrid
                ? (_compactNodesEnd >= 0 ? _compactNodesEnd : _transitionOffset)
                : _maxNodeEndOffset;
            int compactNodeSize = _initialLayout.NodeSize;
            int totalNodes = nodeCount;
            int reportInterval = Math.Max(totalNodes / 200, 1);
            int reportCountdown = reportInterval;

            int idx = 0;

            // Compact zone — sequential scan
            for (long off = hdr; off < compactEnd; off += compactNodeSize, idx++)
            {
                uint nodeStart, nodeEnd;
                if (mmfMain != null)
                {
                    nodeStart = mmfMain.ReadUInt32Unchecked(off);
                    nodeEnd = mmfMain.ReadUInt32Unchecked(off + 4);
                }
                else
                {
                    nodeStart = _storage.ReadUInt32(off);
                    nodeEnd = _storage.ReadUInt32(off + 4);
                }

                if (nodeEnd == PersistentConstants.BOUNDLESS)
                {
                    // Leaf: leafCount = 1, remaining = 0 (already zero)
                    *(uint*)(tp + (long)idx * TEMP_ENTRY + TEMP_OFF_LEAFCOUNT) = 1;
                }
                else
                {
                    // Internal node: process children
                    Pass1Internal(off, idx, _initialLayout, /* useMmfUnchecked */ mmfMain != null,
                        mmfMain, childBase, tp, cKeys, cOffs, CHILD_BUF_MAX,
                        ref batchPos, ref batchStart, ref batchEnd, batchSize);
                }

                if (--reportCountdown <= 0)
                {
                    reportCountdown = reportInterval;
                    _progress?.Report(("Pass 1: child arrays", 92.0 + (double)idx / totalNodes * 2.0));
                }
            }

            // Large zone — sequential scan (hybrid only)
            if (IsHybrid)
            {
                int largeNodeSize = NodeLayout.Large.NodeSize;
                for (long off = _transitionOffset; off < _maxNodeEndOffset; off += largeNodeSize, idx++)
                {
                    uint nodeEnd = _storage.ReadUInt32(off + 4);

                    if (nodeEnd == PersistentConstants.BOUNDLESS)
                    {
                        *(uint*)(tp + (long)idx * TEMP_ENTRY + TEMP_OFF_LEAFCOUNT) = 1;
                    }
                    else
                    {
                        Pass1Internal(off, idx, NodeLayout.Large, /* useMmfUnchecked */ false,
                            mmfMain, childBase, tp, cKeys, cOffs, CHILD_BUF_MAX,
                            ref batchPos, ref batchStart, ref batchEnd, batchSize);
                    }

                    if (--reportCountdown <= 0)
                    {
                        reportCountdown = reportInterval;
                        _progress?.Report(("Pass 1: child arrays", 92.0 + (double)idx / totalNodes * 2.0));
                    }
                }
            }

            // Dispose child store — free its page cache for subsequent passes
            childBase = null;
            if (_ownsChildStore) _childStore.Dispose();

            // ═══════════ PASS 2: Bottom-up leaf count propagation ═══════════
            _progress?.Report(("Pass 2: leaf counts", 94.0));

            long queueCap = (long)nodeCount * 4;
            using var tempQueueMmf = new MappedFileStorageProvider(tempQueuePath, Math.Max(queueCap, 65536));
            tempQueueMmf.SetSize(queueCap);
            byte* qp = tempQueueMmf.RawPointer;

            // Enqueue all leaves (remaining == 0 AND leafCount == 1)
            long qWrite = 0;
            for (int i = 0; i < nodeCount; i++)
            {
                byte* e = tp + (long)i * TEMP_ENTRY;
                if (*(int*)(e + TEMP_OFF_REMAINING) == 0 && *(uint*)(e + TEMP_OFF_LEAFCOUNT) > 0)
                {
                    *(int*)(qp + qWrite) = i;
                    qWrite += 4;
                }
            }

            // BFS: propagate leaf counts bottom-up
            long qRead = 0;
            while (qRead < qWrite)
            {
                int ni = *(int*)(qp + qRead);
                qRead += 4;

                byte* e = tp + (long)ni * TEMP_ENTRY;
                int parentIdx = *(int*)(e + TEMP_OFF_PARENT);
                if (parentIdx < 0) continue; // root — no parent

                uint myLeaves = *(uint*)(e + TEMP_OFF_LEAFCOUNT);
                byte* pe = tp + (long)parentIdx * TEMP_ENTRY;
                *(uint*)(pe + TEMP_OFF_LEAFCOUNT) += myLeaves;
                int rem = *(int*)(pe + TEMP_OFF_REMAINING) - 1;
                *(int*)(pe + TEMP_OFF_REMAINING) = rem;

                if (rem == 0)
                {
                    *(int*)(qp + qWrite) = parentIdx;
                    qWrite += 4;
                }
            }

            // ═══════════ PASS 3: Write LeafCount + find deepest internal node ═══════════
            _progress?.Report(("Pass 3: write-back", 96.0));

            long deepestOffset = rootOffset;
            int maxLrsDepth = 0;
            idx = 0;
            reportCountdown = reportInterval;

            // Compact zone
            for (long off = hdr; off < compactEnd; off += compactNodeSize, idx++)
            {
                uint lc = *(uint*)(tp + (long)idx * TEMP_ENTRY + TEMP_OFF_LEAFCOUNT);
                if (mmfMain != null)
                    mmfMain.WriteUInt32Unchecked(off + (uint)_initialLayout.OffsetLeafCount, lc);
                else
                    _storage.WriteUInt32(off + _initialLayout.OffsetLeafCount, lc);

                // Track deepest internal node (read Start/End sequentially)
                uint nodeEnd = mmfMain != null
                    ? mmfMain.ReadUInt32Unchecked(off + 4)
                    : _storage.ReadUInt32(off + 4);
                if (nodeEnd != PersistentConstants.BOUNDLESS)
                {
                    uint nodeStart = mmfMain != null
                        ? mmfMain.ReadUInt32Unchecked(off)
                        : _storage.ReadUInt32(off);
                    int edgeLen = (int)(nodeEnd - nodeStart);
                    uint dfr = GetBuildDepth(off);
                    int sd = (int)dfr + edgeLen;
                    if (sd > maxLrsDepth) { maxLrsDepth = sd; deepestOffset = off; }
                }

                if (--reportCountdown <= 0)
                {
                    reportCountdown = reportInterval;
                    _progress?.Report(("Pass 3: write-back", 96.0 + (double)idx / totalNodes * 1.5));
                }
            }

            // Large zone
            if (IsHybrid)
            {
                int largeNodeSize = NodeLayout.Large.NodeSize;
                for (long off = _transitionOffset; off < _maxNodeEndOffset; off += largeNodeSize, idx++)
                {
                    uint lc = *(uint*)(tp + (long)idx * TEMP_ENTRY + TEMP_OFF_LEAFCOUNT);
                    _storage.WriteUInt32(off + NodeLayout.Large.OffsetLeafCount, lc);

                    uint nodeEnd = _storage.ReadUInt32(off + 4);
                    if (nodeEnd != PersistentConstants.BOUNDLESS)
                    {
                        uint nodeStart = _storage.ReadUInt32(off);
                        int edgeLen = (int)(nodeEnd - nodeStart);
                        uint dfr = GetBuildDepth(off);
                        int sd = (int)dfr + edgeLen;
                        if (sd > maxLrsDepth) { maxLrsDepth = sd; deepestOffset = off; }
                    }

                    if (--reportCountdown <= 0)
                    {
                        reportCountdown = reportInterval;
                        _progress?.Report(("Pass 3: write-back", 96.0 + (double)idx / totalNodes * 1.5));
                    }
                }
            }

            _deepestInternalNodeOffset = deepestOffset;
            _lrsDepth = maxLrsDepth;

            // Dispose depth store after Pass 3 (needed for deepest node tracking)
            if (_ownsDepthStore) _depthStore?.Dispose();
        }
        finally
        {
            try { File.Delete(tempNodePath); } catch (IOException) { }
            try { File.Delete(tempQueuePath); } catch (IOException) { }
        }
    }

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
                    long entryOff = ci * CHILD_ENTRY_SIZE;
                    cKeys[childCount] = _childStore.ReadUInt32(entryOff + CE_OFF_KEY);
                    cOffs[childCount] = _childStore.ReadInt64(entryOff + CE_OFF_CHILD);
                    childCount++;
                    ci = _childStore.ReadInt32(entryOff + CE_OFF_NEXT);
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

    // ──────────────── Builder child management (off-heap linked lists) ────────────────

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
