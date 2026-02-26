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
public partial class PersistentSuffixTreeBuilder
#pragma warning restore CA1001
{
    private readonly IStorageProvider _storage;
    private readonly MappedFileStorageProvider? _mmfStorage;   // non-null when storage is MMF — enables unchecked fast-path
    private readonly MappedFileStorageProvider? _mmfChildStore; // non-null when child storage is MMF
    private readonly ChildStoreAdapter _childStoreAdapter;
    private readonly DepthStoreAdapter _depthStoreAdapter;
    private readonly NodeLayout _initialLayout;  // layout before any transition (always the compact variant)
    private NodeLayout _layout;                // switches Compact → Large on overflow
    private long _compactOffsetLimit = NodeLayout.CompactMaxOffset;
    private readonly long _rootOffset;
    private ITextSource _text = new StringTextSource(string.Empty);
    private string? _rawString;
    private int _nodeCount;
    private readonly HybridTransitionState _hybridTransitionState = new();
    private readonly UkkonenRuntimeState _ukkonenRuntimeState = new();

    // Hybrid continuation state
    // Jump offsets for child arrays are stored directly in each compact internal
    // node's LeafCount field during MaterializeJumpTable (zero extra RAM).

    // Jump table reserve: block pre-allocated within the compact zone at transition time.
    // This ensures jump entry offsets fit in uint32 fields of compact-zone nodes.
    private const int JUMP_TABLE_MIN_RESERVE = 16 * 1024 * 1024; // 16 MB floor for small trees
    private const int JUMP_TABLE_RESERVE_FRACTION = 4;            // reserve 1/4 (25%) of compact zone
    // NOTE: 25% is needed because ALL compact internal nodes require jump slots —
    // child arrays are allocated in the large zone (>4 GB) during the DFS phase,
    // and their offsets cannot fit in the compact uint32 ChildrenHead field.

    // Algorithm state
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

    // Accessors preserving the existing field names while storing state in dedicated containers.
    private long _transitionOffset
    {
        get => _hybridTransitionState.TransitionOffset;
        set => _hybridTransitionState.TransitionOffset = value;
    }

    private long _jumpTableStart
    {
        get => _hybridTransitionState.JumpTableStart;
        set => _hybridTransitionState.JumpTableStart = value;
    }

    private long _jumpTableEnd
    {
        get => _hybridTransitionState.JumpTableEnd;
        set => _hybridTransitionState.JumpTableEnd = value;
    }

    private long _compactNodesEnd
    {
        get => _hybridTransitionState.CompactNodesEnd;
        set => _hybridTransitionState.CompactNodesEnd = value;
    }

    private long _jumpTableReserveStart
    {
        get => _hybridTransitionState.JumpTableReserveStart;
        set => _hybridTransitionState.JumpTableReserveStart = value;
    }

    private long _jumpTableReserveEnd
    {
        get => _hybridTransitionState.JumpTableReserveEnd;
        set => _hybridTransitionState.JumpTableReserveEnd = value;
    }

    private long _activeNodeOffset
    {
        get => _ukkonenRuntimeState.ActiveNodeOffset;
        set => _ukkonenRuntimeState.ActiveNodeOffset = value;
    }

    private int _activeEdgeIndex
    {
        get => _ukkonenRuntimeState.ActiveEdgeIndex;
        set => _ukkonenRuntimeState.ActiveEdgeIndex = value;
    }

    private int _activeLength
    {
        get => _ukkonenRuntimeState.ActiveLength;
        set => _ukkonenRuntimeState.ActiveLength = value;
    }

    private int _remainder
    {
        get => _ukkonenRuntimeState.Remainder;
        set => _ukkonenRuntimeState.Remainder = value;
    }

    private int _position
    {
        get => _ukkonenRuntimeState.Position;
        set => _ukkonenRuntimeState.Position = value;
    }

    private long _lastCreatedInternalNodeOffset
    {
        get => _ukkonenRuntimeState.LastCreatedInternalNodeOffset;
        set => _ukkonenRuntimeState.LastCreatedInternalNodeOffset = value;
    }

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
        _childStoreAdapter = new ChildStoreAdapter(_childStore, _mmfChildStore);

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
        _depthStoreAdapter = new DepthStoreAdapter(_depthStore, _mmfDepthStore);

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
}
