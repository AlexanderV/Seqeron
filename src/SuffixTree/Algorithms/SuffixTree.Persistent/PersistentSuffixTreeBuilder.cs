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

        _progress?.Report(("DFS: Topology extraction", 92.0));

        CombinedFinalizeDFS(_rootOffset);

        if (_ownsDepthStore) _depthStore?.Dispose();
        if (_ownsChildStore) _childStore.Dispose();

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

        // Pass 1: count compact internal nodes that need child-array jump slots.
        // ALL compact internal nodes need them because child arrays are allocated
        // in the large zone (>4 GB) where offsets exceed uint32.
        int childArrayJumps = 0;
        for (long offset = HeaderSize; offset < compactEnd; offset += _initialLayout.NodeSize)
        {
            uint nodeEnd = mmfMain != null
                ? mmfMain.ReadUInt32Unchecked(offset + 4)
                : _storage.ReadUInt32(offset + 4);
            if (nodeEnd != PersistentConstants.BOUNDLESS)
                childArrayJumps++;
        }

        int totalEntries = suffixLinkJumps + childArrayJumps;
        if (totalEntries == 0)
        {
            _jumpTableStart = -1;
            _jumpTableEnd = -1;
            return;
        }

        // Use pre-reserved space within the compact zone if available,
        // otherwise allocate normally (works for tests with tiny compact limits).
        int tableSize = totalEntries * 8;
        if (_jumpTableReserveStart >= 0)
        {
            long reserveCapacity = _jumpTableReserveEnd - _jumpTableReserveStart;
            if (tableSize > reserveCapacity)
                throw new InvalidOperationException(
                    $"Jump table requires {tableSize:N0} bytes ({totalEntries:N0} entries) " +
                    $"but only {reserveCapacity:N0} bytes were reserved. " +
                    "The input exceeds the jump table reserve capacity.");
            _jumpTableStart = _jumpTableReserveStart;
            _jumpTableEnd = _jumpTableReserveStart + tableSize;
        }
        else
        {
            _jumpTableStart = _storage.Allocate(tableSize);
            _jumpTableEnd = _jumpTableStart + tableSize;
        }

        // Fill suffix-link entries first
        int slotIndex = 0;
        for (int i = 0; i < suffixLinkJumps; i++, slotIndex++)
        {
            var (compactNodeOffset, largeTargetOffset) = _deferredSuffixLinks[i];
            long jumpEntryOffset = _jumpTableStart + (long)slotIndex * 8;
            _storage.WriteInt64(jumpEntryOffset, largeTargetOffset);

            // Point the compact node's SuffixLink to the jump entry
            var compactLayout = LayoutOf(compactNodeOffset);
            var node = new PersistentSuffixTreeNode(_storage, compactNodeOffset, compactLayout);
            node.SuffixLink = jumpEntryOffset;
        }

        // Pass 2: write child-array jump entry offsets directly into each
        // compact internal node's LeafCount field (uint32 at offset+12).
        // The jump table lives in the compact zone so offsets always fit in uint32.
        // This replaces the old 1.4 GB flat long[] array — ZERO extra RAM.
        // The DFS will read the jump offset from LeafCount before overwriting
        // it with the real leaf count.
        for (long offset = HeaderSize; offset < compactEnd; offset += _initialLayout.NodeSize)
        {
            uint nodeEnd = mmfMain != null
                ? mmfMain.ReadUInt32Unchecked(offset + 4)
                : _storage.ReadUInt32(offset + 4);
            if (nodeEnd != PersistentConstants.BOUNDLESS)
            {
                long jumpEntryOffset = _jumpTableStart + (long)slotIndex * 8;
                // Store jump offset in LeafCount field (temporary, overwritten by DFS later)
                if (mmfMain != null)
                    mmfMain.WriteUInt32Unchecked(offset + 12, (uint)jumpEntryOffset);
                else
                    _storage.WriteUInt32(offset + 12, (uint)jumpEntryOffset);
                slotIndex++;
            }
        }
    }

    // ──────────────── Combined finalize DFS ────────────────

    /// <summary>
    /// Single post-order DFS that does ALL finalization work in one pass:
    ///   • Computes LeafCount for every node
    ///   • Writes sorted child arrays (batch-allocated — ~30 allocs, not 175M)
    ///   • Tracks the deepest internal node (depth accumulated on stack — no depth store)
    ///   • Reads jump offsets from LeafCount field, then overwrites with the real count
    ///
    /// Children are cached in a native flat buffer during the first visit so that
    /// the child store is read only ONCE per internal node (not 2-3 times).
    /// LeafCounts are written directly into main storage — no intermediate 1.4 GB array.
    ///
    /// Memory: DFS stack ~5-20 MB + child cache ~3 MB, zero large arrays.
    /// I/O: single pass over the tree — minimum possible page touches.
    /// </summary>
    private unsafe void CombinedFinalizeDFS(long rootOffset)
    {
        var mmfMain = _mmfStorage;
        var mmfChild = _mmfChildStore;

        // ── Native child cache ──
        // During first visit we read (key, childOffset) from child store and cache them
        // in a flat native buffer. Post-order visit reads from cache — zero child store reads.
        // Buffer holds entries for nodes currently on the DFS path: O(Height × Branching).
        // For DNA genome: ~50K depth × 5 branching × 12 bytes = ~3 MB.
        const int CHILD_CACHE_ENTRY = 12; // uint key (4) + long childOffset (8)
        int childCacheCapacity = 512 * 1024; // initial entries
        byte* childCache = (byte*)System.Runtime.InteropServices.NativeMemory.Alloc(
            (nuint)childCacheCapacity * CHILD_CACHE_ENTRY);
        int childCacheUsed = 0; // high-water mark of entries currently in use

        // ── DFS stack ──
        // Entry: (Offset, ChildCacheStart, ChildCount, DepthAtEdgeStart)
        //   ChildCacheStart = -1 means first visit (children not yet pushed)
        //   ChildCacheStart >= 0 means post-order visit, children cached at [start..start+count)
        var stack = new System.Collections.Generic.Stack<(long Offset, int ChildCacheStart, int ChildCount, int DepthAtEdgeStart)>(4096);
        stack.Push((rootOffset, -1, 0, 0));

        // ── Batch allocation state for child arrays ──
        const int MAX_BATCH = 256 * 1024 * 1024;
        int batchSize = Math.Min(MAX_BATCH, Math.Max(4096, _nodeCount * 48));
        long batchStart = _storage.Allocate(batchSize);
        long batchEnd = batchStart + batchSize;
        long batchPos = batchStart;

        // ── Reusable buffers ──
        var childBuf = new System.Collections.Generic.List<(uint Key, long ChildOffset)>(8);
        byte[] writeBuf = System.Buffers.ArrayPool<byte>.Shared.Rent(256);

        long deepestOffset = rootOffset;
        int maxLrsDepth = 0;
        int nodesProcessed = 0;
        int totalNodes = _nodeCount;
        int reportInterval = Math.Max(totalNodes / 200, 1);

        try
        {
            while (stack.Count > 0)
            {
                var (offset, cacheStart, cacheCount, depthAtEdge) = stack.Pop();
                var layout = LayoutOf(offset);
                bool useMmf = mmfMain != null && !layout.OffsetIs64Bit;

                // ── Read Start / End ──
                uint nodeStart, nodeEnd;
                if (useMmf)
                {
                    nodeStart = mmfMain!.ReadUInt32Unchecked(offset);
                    nodeEnd = mmfMain!.ReadUInt32Unchecked(offset + 4);
                }
                else
                {
                    var n = new PersistentSuffixTreeNode(_storage, offset, layout);
                    nodeStart = n.Start;
                    nodeEnd = n.End;
                }

                // ── Leaf ──
                if (nodeEnd == PersistentConstants.BOUNDLESS)
                {
                    // Write LeafCount = 1 directly to main storage
                    if (useMmf)
                        mmfMain!.WriteUInt32Unchecked(offset + (uint)layout.OffsetLeafCount, 1);
                    else
                        new PersistentSuffixTreeNode(_storage, offset, layout).LeafCount = 1;

                    nodesProcessed++;
                    if (nodesProcessed % reportInterval == 0)
                        _progress?.Report(("DFS: Topology extraction", 92.0 + (double)nodesProcessed / totalNodes * 5.0));
                    continue;
                }

                int edgeLen = (int)(nodeEnd - nodeStart);
                int stringDepth = depthAtEdge + edgeLen;

                if (cacheStart == -1)
                {
                    // ══════ First visit: read children from child store, cache them, push for DFS ══════

                    long headIndex;
                    if (useMmf)
                    {
                        uint raw = mmfMain!.ReadUInt32Unchecked(offset + 16);
                        headIndex = raw == uint.MaxValue ? PersistentConstants.NULL_OFFSET : (long)raw;
                    }
                    else
                    {
                        headIndex = new PersistentSuffixTreeNode(_storage, offset, layout).ChildrenHead;
                    }

                    // Count and cache children — THE ONLY child store read for this node
                    int myStart = childCacheUsed;
                    int myCount = 0;

                    if (headIndex != PersistentConstants.NULL_OFFSET)
                    {
                        long ci = headIndex;
                        if (mmfChild != null)
                        {
                            while (ci >= 0)
                            {
                                long entryOff = ci * CHILD_ENTRY_SIZE;
                                uint key = mmfChild.ReadUInt32Unchecked(entryOff + CE_OFF_KEY);
                                long childOff = mmfChild.ReadInt64Unchecked(entryOff + CE_OFF_CHILD);

                                // Grow cache if needed
                                if (childCacheUsed >= childCacheCapacity)
                                {
                                    int newCap = childCacheCapacity * 2;
                                    byte* newBuf = (byte*)System.Runtime.InteropServices.NativeMemory.Alloc(
                                        (nuint)newCap * CHILD_CACHE_ENTRY);
                                    System.Buffer.MemoryCopy(childCache, newBuf,
                                        (long)newCap * CHILD_CACHE_ENTRY,
                                        (long)childCacheCapacity * CHILD_CACHE_ENTRY);
                                    System.Runtime.InteropServices.NativeMemory.Free(childCache);
                                    childCache = newBuf;
                                    childCacheCapacity = newCap;
                                }

                                byte* cEntry = childCache + (long)childCacheUsed * CHILD_CACHE_ENTRY;
                                *(uint*)cEntry = key;
                                *(long*)(cEntry + 4) = childOff;
                                childCacheUsed++;
                                myCount++;

                                ci = mmfChild.ReadInt32Unchecked(entryOff + CE_OFF_NEXT);
                            }
                        }
                        else
                        {
                            while (ci >= 0)
                            {
                                long entryOff = ci * CHILD_ENTRY_SIZE;
                                uint key = _childStore.ReadUInt32(entryOff + CE_OFF_KEY);
                                long childOff = _childStore.ReadInt64(entryOff + CE_OFF_CHILD);

                                // Grow cache if needed
                                if (childCacheUsed >= childCacheCapacity)
                                {
                                    int newCap = childCacheCapacity * 2;
                                    byte* newBuf = (byte*)System.Runtime.InteropServices.NativeMemory.Alloc(
                                        (nuint)newCap * CHILD_CACHE_ENTRY);
                                    System.Buffer.MemoryCopy(childCache, newBuf,
                                        (long)newCap * CHILD_CACHE_ENTRY,
                                        (long)childCacheCapacity * CHILD_CACHE_ENTRY);
                                    System.Runtime.InteropServices.NativeMemory.Free(childCache);
                                    childCache = newBuf;
                                    childCacheCapacity = newCap;
                                }

                                byte* cEntry = childCache + (long)childCacheUsed * CHILD_CACHE_ENTRY;
                                *(uint*)cEntry = key;
                                *(long*)(cEntry + 4) = childOff;
                                childCacheUsed++;
                                myCount++;

                                ci = _childStore.ReadInt32(entryOff + CE_OFF_NEXT);
                            }
                        }
                    }

                    // Push post-order marker with cached children info
                    stack.Push((offset, myStart, myCount, depthAtEdge));

                    // Push children for DFS traversal (read from cache, not child store)
                    for (int c = myStart; c < myStart + myCount; c++)
                    {
                        byte* cEntry = childCache + (long)c * CHILD_CACHE_ENTRY;
                        long childOff = *(long*)(cEntry + 4);
                        stack.Push((childOff, -1, 0, stringDepth));
                    }
                }
                else
                {
                    // ══════ Post-order visit: all children are finalized ══════

                    nodesProcessed++;
                    if (nodesProcessed % reportInterval == 0)
                        _progress?.Report(("DFS: Topology extraction", 92.0 + (double)nodesProcessed / totalNodes * 5.0));

                    // 1. Track deepest internal node
                    if (stringDepth > maxLrsDepth)
                    {
                        maxLrsDepth = stringDepth;
                        deepestOffset = offset;
                    }

                    // 2. Read jump slot offset from LeafCount (put there by MaterializeJumpTable).
                    bool isCompact = !layout.OffsetIs64Bit;
                    long jumpSlotOffset = -1;
                    bool hasJumpSlot = false;
                    if (isCompact && IsHybrid)
                    {
                        uint jumpOff = useMmf
                            ? mmfMain!.ReadUInt32Unchecked(offset + 12)
                            : _storage.ReadUInt32(offset + 12);
                        if (jumpOff > 0)
                        {
                            jumpSlotOffset = (long)jumpOff;
                            hasJumpSlot = true;
                        }
                    }

                    // 3. Collect children from cache + sum leaf counts from main storage
                    childBuf.Clear();
                    uint totalLeaves = 0;

                    for (int c = cacheStart; c < cacheStart + cacheCount; c++)
                    {
                        byte* cEntry = childCache + (long)c * CHILD_CACHE_ENTRY;
                        uint key = *(uint*)cEntry;
                        long childOff = *(long*)(cEntry + 4);
                        childBuf.Add((key, childOff));

                        // Read leaf count from child node's storage (already written by child's post-order)
                        var cl = LayoutOf(childOff);
                        bool cUseMmf = mmfMain != null && !cl.OffsetIs64Bit;
                        totalLeaves += cUseMmf
                            ? mmfMain!.ReadUInt32Unchecked(childOff + (uint)cl.OffsetLeafCount)
                            : new PersistentSuffixTreeNode(_storage, childOff, cl).LeafCount;
                    }

                    // Release cache entries (all children processed, allow reuse for sibling subtrees)
                    childCacheUsed = cacheStart;

                    // 4. Sort children by key
                    childBuf.Sort((a, b) => ((int)a.Key).CompareTo((int)b.Key));
                    int count = childBuf.Count;

                    // 5. Write sorted child array (batch-allocated)
                    long arrayOffset = PersistentConstants.NULL_OFFSET;
                    if (count > 0)
                    {
                        NodeLayout arrayLayout = hasJumpSlot ? NodeLayout.Large : layout;
                        int entrySize = arrayLayout.ChildEntrySize;
                        int totalBytes = checked(count * entrySize);

                        // Bump-allocate from current batch
                        if (batchPos + totalBytes > batchEnd)
                        {
                            int allocSz = Math.Max(batchSize, totalBytes);
                            batchStart = _storage.Allocate(allocSz);
                            batchEnd = batchStart + allocSz;
                            batchPos = batchStart;
                        }
                        arrayOffset = batchPos;
                        batchPos += totalBytes;

                        // Serialize
                        if (totalBytes > writeBuf.Length)
                        {
                            byte[] newBuf = System.Buffers.ArrayPool<byte>.Shared.Rent(totalBytes);
                            SerializeAndWriteChildArray(childBuf, arrayLayout, newBuf, totalBytes, arrayOffset);
                            System.Buffers.ArrayPool<byte>.Shared.Return(newBuf);
                        }
                        else
                        {
                            SerializeAndWriteChildArray(childBuf, arrayLayout, writeBuf, totalBytes, arrayOffset);
                        }
                    }

                    // 6. Update node header
                    int cCountRaw = count | (hasJumpSlot ? unchecked((int)0x80000000) : 0);

                    if (hasJumpSlot)
                    {
                        _storage.WriteInt64(jumpSlotOffset, arrayOffset);
                        if (useMmf)
                        {
                            mmfMain!.WriteUInt32Unchecked(offset + (uint)layout.OffsetLeafCount, totalLeaves);
                            if (!layout.OffsetIs64Bit)
                            {
                                mmfMain!.WriteUInt32Unchecked(offset + (uint)layout.OffsetChildrenHead, (uint)jumpSlotOffset);
                                mmfMain!.WriteInt32Unchecked(offset + layout.OffsetChildCount, cCountRaw);
                            }
                            else
                            {
                                mmfMain!.WriteInt64Unchecked(offset + layout.OffsetChildrenHead, jumpSlotOffset);
                                mmfMain!.WriteInt32Unchecked(offset + layout.OffsetChildCount, cCountRaw);
                            }
                        }
                        else
                        {
                            var node = new PersistentSuffixTreeNode(_storage, offset, layout);
                            node.LeafCount = totalLeaves;
                            node.ChildrenHead = jumpSlotOffset;
                            node.ChildCount = cCountRaw;
                        }
                    }
                    else
                    {
                        if (useMmf)
                        {
                            mmfMain!.WriteUInt32Unchecked(offset + (uint)layout.OffsetLeafCount, totalLeaves);
                            if (!layout.OffsetIs64Bit)
                            {
                                mmfMain!.WriteUInt32Unchecked(offset + (uint)layout.OffsetChildrenHead, (uint)arrayOffset);
                                mmfMain!.WriteInt32Unchecked(offset + layout.OffsetChildCount, cCountRaw);
                            }
                            else
                            {
                                mmfMain!.WriteInt64Unchecked(offset + layout.OffsetChildrenHead, arrayOffset);
                                mmfMain!.WriteInt32Unchecked(offset + layout.OffsetChildCount, cCountRaw);
                            }
                        }
                        else
                        {
                            var node = new PersistentSuffixTreeNode(_storage, offset, layout);
                            node.LeafCount = totalLeaves;
                            node.ChildrenHead = arrayOffset;
                            node.ChildCount = cCountRaw;
                        }
                    }
                }
            }
        }
        finally
        {
            System.Buffers.ArrayPool<byte>.Shared.Return(writeBuf);
            System.Runtime.InteropServices.NativeMemory.Free(childCache);
        }

        _deepestInternalNodeOffset = deepestOffset;
        _lrsDepth = maxLrsDepth;
    }

    /// <summary>Serializes child entries into a buffer and writes to storage.</summary>
    private void SerializeAndWriteChildArray(
        List<(uint Key, long ChildOffset)> children,
        NodeLayout arrayLayout, byte[] buf, int totalBytes, long arrayOffset)
    {
        int entrySize = arrayLayout.ChildEntrySize;
        for (int i = 0; i < children.Count; i++)
        {
            int off = checked(i * entrySize);
            BinaryPrimitives.WriteUInt32LittleEndian(buf.AsSpan(off, 4), children[i].Key);
            if (arrayLayout.OffsetIs64Bit)
                BinaryPrimitives.WriteInt64LittleEndian(buf.AsSpan(off + 4, 8), children[i].ChildOffset);
            else
                BinaryPrimitives.WriteUInt32LittleEndian(buf.AsSpan(off + 4, 4), (uint)children[i].ChildOffset);
        }
        _storage.WriteBytes(arrayOffset, buf, 0, totalBytes);
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
