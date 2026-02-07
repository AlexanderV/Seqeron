namespace SuffixTree.Persistent;

public static class PersistentConstants
{
    // ──── Legacy v3 (Large) layout sizes — superseded by NodeLayout ────
    // Kept for reference only. All runtime code uses NodeLayout.NodeSize / .ChildEntrySize.

    /// <summary>Node size for v3 (Large) format. Use <see cref="NodeLayout.NodeSize"/> instead.</summary>
    [Obsolete("Use NodeLayout.NodeSize", error: true)]
    public const int NODE_SIZE = 40;

    /// <summary>Child entry size for v3 (Large) format. Use <see cref="NodeLayout.ChildEntrySize"/> instead.</summary>
    [Obsolete("Use NodeLayout.ChildEntrySize", error: true)]
    public const int CHILD_ENTRY_SIZE = 12;

    public const int HEADER_SIZE = 48;

    /// <summary>
    /// Extended header size for Hybrid v5 format (80 bytes).
    /// Extra fields beyond base 48: TRANSITION (int64) + JUMP_START (int64)
    /// + JUMP_END (int64) + DEEPEST_NODE (int64).
    /// </summary>
    public const int HEADER_SIZE_V5 = 80;

    public const uint BOUNDLESS = uint.MaxValue;
    public const uint TERMINATOR_KEY = uint.MaxValue;  // Same bit pattern as -1 when cast to int
    public const long NULL_OFFSET = -1L;
    public const long MAGIC_NUMBER = 0x5452454558494646L; // "SUFFIXTR"

    /// <summary>
    /// Binary storage format version for the v3 (Large) layout.
    /// Current default is v4 (Compact). Use <see cref="NodeLayout.Version"/> instead.
    /// Independent of <see cref="SuffixTreeSerializer"/> logical format (v2).
    /// </summary>
    [Obsolete("Use NodeLayout.Version — format is auto-selected", error: true)]
    public const int CURRENT_VERSION = 3;

    // Header Layout (48 bytes, packed, naturally aligned)
    //  0-7  : MAGIC      (int64)
    //  8-11 : VERSION    (int32)
    // 12-15 : TEXT_LEN   (int32)
    // 16-23 : ROOT       (int64)
    // 24-31 : TEXT_OFF   (int64)
    // 32-35 : NODE_COUNT (int32)
    // 36-39 : reserved
    // 40-47 : SIZE       (int64)
    public const int HEADER_OFFSET_MAGIC = 0;
    public const int HEADER_OFFSET_VERSION = 8;
    public const int HEADER_OFFSET_TEXT_LEN = 12;
    public const int HEADER_OFFSET_ROOT = 16;
    public const int HEADER_OFFSET_TEXT_OFF = 24;
    public const int HEADER_OFFSET_NODE_COUNT = 32;
    public const int HEADER_OFFSET_SIZE = 40;

    // ──── Additional v5 (Hybrid) header fields ────
    // 48-55 : TRANSITION  (int64) — first offset of large-zone nodes
    // 56-63 : JUMP_START  (int64) — first offset of the contiguous jump table
    // 64-71 : JUMP_END    (int64) — first offset after the jump table
    /// <summary>Offset of the compact→large transition boundary in a v5 header.</summary>
    public const int HEADER_OFFSET_TRANSITION = 48;
    /// <summary>Start of the contiguous jump table in a v5 header.</summary>
    public const int HEADER_OFFSET_JUMP_START = 56;
    /// <summary>End of the contiguous jump table in a v5 header.</summary>
    public const int HEADER_OFFSET_JUMP_END = 64;

    /// <summary>Offset of the pre-computed deepest internal node (for O(1) LRS). 8 bytes.</summary>
    public const int HEADER_OFFSET_DEEPEST_NODE = 72;

    // Invariant node field offsets (same in both Compact and Large layouts)
    public const int OFFSET_START = 0;
    public const int OFFSET_END = 4;

    // ──── Legacy v3 (Large) node field offsets — superseded by NodeLayout ────

    /// <summary>SuffixLink offset for v3 (Large) format. Use <see cref="NodeLayout.OffsetSuffixLink"/> instead.</summary>
    [Obsolete("Use NodeLayout.OffsetSuffixLink", error: true)]
    public const int OFFSET_SUFFIX_LINK = 8;

    /// <summary>Depth offset for v3 (Large) format. Use <see cref="NodeLayout.OffsetDepth"/> instead.</summary>
    [Obsolete("Use NodeLayout.OffsetDepth", error: true)]
    public const int OFFSET_DEPTH = 16;

    /// <summary>LeafCount offset for v3 (Large) format. Use <see cref="NodeLayout.OffsetLeafCount"/> instead.</summary>
    [Obsolete("Use NodeLayout.OffsetLeafCount", error: true)]
    public const int OFFSET_LEAF_COUNT = 20;

    /// <summary>ChildrenHead offset for v3 (Large) format. Use <see cref="NodeLayout.OffsetChildrenHead"/> instead.</summary>
    [Obsolete("Use NodeLayout.OffsetChildrenHead", error: true)]
    public const int OFFSET_CHILDREN_HEAD = 24;

    /// <summary>ChildCount offset for v3 (Large) format. Use <see cref="NodeLayout.OffsetChildCount"/> instead.</summary>
    [Obsolete("Use NodeLayout.OffsetChildCount", error: true)]
    public const int OFFSET_CHILD_COUNT = 32;

    // Child Entry Layout Offsets — use NodeLayout.ChildOffsetKey / ChildOffsetNode instead
    [Obsolete("Use NodeLayout.ChildOffsetKey", error: true)]
    public const int CHILD_OFFSET_KEY = 0;

    [Obsolete("Use NodeLayout.ChildOffsetNode", error: true)]
    public const int CHILD_OFFSET_NODE = 4;
}
