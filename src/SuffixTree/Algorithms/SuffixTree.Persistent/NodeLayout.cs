using System.Runtime.CompilerServices;

namespace SuffixTree.Persistent;

/// <summary>
/// Describes the binary layout of nodes and child entries in persistent storage.
/// Two singleton instances for the v6 format:
/// <list type="bullet">
///   <item><see cref="Compact"/> — 24-byte nodes, 32-bit offsets, 8-byte child entries.</item>
///   <item><see cref="Large"/> — 32-byte nodes, 64-bit offsets, 12-byte child entries.</item>
/// </list>
/// <para>
/// Thread-safe. Instances are immutable and shared across all trees of the same format.
/// </para>
/// </summary>
public sealed class NodeLayout
{
    // ──────────────── Singleton instances ────────────────

    /// <summary>
    /// Compact format (v6): 32-bit offsets, 24-byte nodes, 8-byte child entries.
    /// <para>Max file ~4 GB → ~72 M characters.</para>
    /// </summary>
    /// <remarks>
    /// <b>Node layout (24 bytes):</b>
    /// <code>
    ///  0- 3 : Start           (uint32)
    ///  4- 7 : End             (uint32)
    ///  8-11 : SuffixLink      (uint32)
    /// 12-15 : LeafCount       (uint32)
    /// 16-19 : ChildrenHead    (uint32)
    /// 20-23 : ChildCount      (int32)
    /// </code>
    /// <b>Child entry (8 bytes):</b> Key(4) + ChildNodeOffset(4)
    /// </remarks>
    public static readonly NodeLayout Compact = new(
        version: 6,
        nodeSize: 24,
        childEntrySize: 8,
        offsetIs64Bit: false,
        offsetSuffixLink: 8,
        offsetLeafCount: 12,
        offsetChildrenHead: 16,
        offsetChildCount: 20);

    /// <summary>
    /// Large format (v6): 64-bit offsets, 32-byte nodes, 12-byte child entries.
    /// Used as the large-zone layout in hybrid trees when the compact address space is exhausted.
    /// <para>No practical size limit.</para>
    /// </summary>
    /// <remarks>
    /// <b>Node layout (32 bytes):</b>
    /// <code>
    ///  0- 3 : Start           (uint32)
    ///  4- 7 : End             (uint32)
    ///  8-15 : SuffixLink      (int64)
    /// 16-19 : LeafCount       (uint32)
    /// 20-27 : ChildrenHead    (int64)
    /// 28-31 : ChildCount      (int32)
    /// </code>
    /// <b>Child entry (12 bytes):</b> Key(4) + ChildNodeOffset(8)
    /// </remarks>
    public static readonly NodeLayout Large = new(
        version: 6,
        nodeSize: 32,
        childEntrySize: 12,
        offsetIs64Bit: true,
        offsetSuffixLink: 8,
        offsetLeafCount: 16,
        offsetChildrenHead: 20,
        offsetChildCount: 28);

    // ──────────────── Fields ────────────────

    /// <summary>Sentinel stored in uint32 offset fields to represent "null".</summary>
    private const uint COMPACT_NULL = uint.MaxValue; // 0xFFFFFFFF

    /// <summary>Format version written to the file header.</summary>
    public int Version { get; }

    /// <summary>Size of a single node block in bytes.</summary>
    public int NodeSize { get; }

    /// <summary>Size of a single child entry in bytes.</summary>
    public int ChildEntrySize { get; }

    /// <summary>True if offset fields (SuffixLink, ChildrenHead, ChildNodeOffset) are 64-bit.</summary>
    public bool OffsetIs64Bit { get; }

    // Node field offsets (Start=0, End=4 are invariant across formats)
    /// <summary>Offset of the SuffixLink field within a node.</summary>
    public int OffsetSuffixLink { get; }
    /// <summary>Offset of the LeafCount field within a node.</summary>
    public int OffsetLeafCount { get; }
    /// <summary>Offset of the ChildrenHead field within a node.</summary>
    public int OffsetChildrenHead { get; }
    /// <summary>Offset of the ChildCount field within a node.</summary>
    public int OffsetChildCount { get; }

    // Child entry field offsets (invariant)
    /// <summary>Offset of the Key field within a child entry (always 0).</summary>
    public const int ChildOffsetKey = 0;
    /// <summary>Offset of the ChildNodeOffset field within a child entry (always 4).</summary>
    public const int ChildOffsetNode = 4;

    // ──────────────── Constructor ────────────────

    private NodeLayout(int version, int nodeSize, int childEntrySize, bool offsetIs64Bit,
        int offsetSuffixLink, int offsetLeafCount,
        int offsetChildrenHead, int offsetChildCount)
    {
        Version = version;
        NodeSize = nodeSize;
        ChildEntrySize = childEntrySize;
        OffsetIs64Bit = offsetIs64Bit;
        OffsetSuffixLink = offsetSuffixLink;
        OffsetLeafCount = offsetLeafCount;
        OffsetChildrenHead = offsetChildrenHead;
        OffsetChildCount = offsetChildCount;
    }

    // ──────────────── Offset read/write helpers ────────────────

    /// <summary>
    /// Reads an offset field (SuffixLink, ChildrenHead, or ChildNodeOffset)
    /// from storage. Handles 32-bit ↔ 64-bit and NULL translation.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public long ReadOffset(IStorageProvider storage, long position)
    {
        if (OffsetIs64Bit)
            return storage.ReadInt64(position);
        uint raw = storage.ReadUInt32(position);
        return raw == COMPACT_NULL ? PersistentConstants.NULL_OFFSET : (long)raw;
    }

    /// <summary>
    /// Writes an offset field to storage. Handles long → uint32 truncation
    /// and NULL translation for compact format.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public void WriteOffset(IStorageProvider storage, long position, long value)
    {
        if (OffsetIs64Bit)
        {
            storage.WriteInt64(position, value);
        }
        else
        {
            if (value != PersistentConstants.NULL_OFFSET && (value < 0 || value > CompactMaxOffset))
                throw new InvalidOperationException(
                    $"WriteOffset: value {value} exceeds Compact format limit ({CompactMaxOffset}). " +
                    "Use Large layout or hybrid continuation for offsets > 4 GB.");
            storage.WriteUInt32(position,
                value == PersistentConstants.NULL_OFFSET ? COMPACT_NULL : (uint)value);
        }
    }

    // ──────────────── Compact offset limit ────────────────

    /// <summary>
    /// Maximum byte offset addressable by the Compact (uint32) format.
    /// <para>
    /// <c>uint.MaxValue</c> is reserved as the NULL sentinel (<see cref="COMPACT_NULL"/>),
    /// so the highest valid offset is <c>uint.MaxValue − 1 = 0xFFFF_FFFE</c> (≈ 4.29 GB).
    /// </para>
    /// </summary>
    internal const long CompactMaxOffset = (long)uint.MaxValue - 1; // 0xFFFF_FFFE

    // ──────────────── Version lookup ────────────────

    /// <summary>
    /// Returns the <see cref="NodeLayout"/> for a given format version read from the header.
    /// Only v6 is supported.
    /// </summary>
    public static NodeLayout ForVersion(int version) => version switch
    {
        6 => Compact, // hybrid v6: base layout is Compact; large zone uses Large
        _ => throw new InvalidOperationException(
            $"Unsupported storage format version: {version}. Expected 6.")
    };
}
