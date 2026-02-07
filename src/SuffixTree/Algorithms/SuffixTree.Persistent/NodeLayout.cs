using System;
using System.Runtime.CompilerServices;

namespace SuffixTree.Persistent;

/// <summary>
/// Describes the binary layout of nodes and child entries in persistent storage.
/// Two singleton instances: <see cref="Compact"/> (32-bit offsets, v4) and
/// <see cref="Large"/> (64-bit offsets, v3).
/// <para>
/// Thread-safe. Instances are immutable and shared across all trees of the same format.
/// </para>
/// </summary>
public sealed class NodeLayout
{
    // ──────────────── Singleton instances ────────────────

    /// <summary>
    /// Compact format (v4): 32-bit offsets, 28-byte nodes, 8-byte child entries.
    /// <para>Max file ~4 GB → ~58 M characters.</para>
    /// </summary>
    /// <remarks>
    /// <b>Node layout (28 bytes):</b>
    /// <code>
    ///  0- 3 : Start           (uint32)
    ///  4- 7 : End             (uint32)
    ///  8-11 : SuffixLink      (uint32)   ← was int64
    /// 12-15 : DepthFromRoot   (uint32)
    /// 16-19 : LeafCount       (uint32)
    /// 20-23 : ChildrenHead    (uint32)   ← was int64
    /// 24-27 : ChildCount      (int32)
    /// </code>
    /// <b>Child entry (8 bytes):</b> Key(4) + ChildNodeOffset(4)
    /// </remarks>
    public static readonly NodeLayout Compact = new(
        version: 4,
        nodeSize: 28,
        childEntrySize: 8,
        offsetIs64Bit: false,
        offsetSuffixLink: 8,
        offsetDepth: 12,
        offsetLeafCount: 16,
        offsetChildrenHead: 20,
        offsetChildCount: 24);

    /// <summary>
    /// Large format (v3): 64-bit offsets, 40-byte nodes, 12-byte child entries.
    /// <para>No practical size limit.</para>
    /// </summary>
    /// <remarks>
    /// <b>Node layout (40 bytes):</b>
    /// <code>
    ///  0- 3 : Start           (uint32)
    ///  4- 7 : End             (uint32)
    ///  8-15 : SuffixLink      (int64)
    /// 16-19 : DepthFromRoot   (uint32)
    /// 20-23 : LeafCount       (uint32)
    /// 24-31 : ChildrenHead    (int64)
    /// 32-35 : ChildCount      (int32)
    /// 36-39 : (padding)
    /// </code>
    /// <b>Child entry (12 bytes):</b> Key(4) + ChildNodeOffset(8)
    /// </remarks>
    public static readonly NodeLayout Large = new(
        version: 3,
        nodeSize: 40,
        childEntrySize: 12,
        offsetIs64Bit: true,
        offsetSuffixLink: 8,
        offsetDepth: 16,
        offsetLeafCount: 20,
        offsetChildrenHead: 24,
        offsetChildCount: 32);

    // ──────────────── Fields ────────────────

    /// <summary>Sentinel stored in uint32 offset fields to represent "null".</summary>
    private const uint COMPACT_NULL = uint.MaxValue; // 0xFFFFFFFF

    /// <summary>Format version written to the file header.</summary>
    public readonly int Version;

    /// <summary>Size of a single node block in bytes.</summary>
    public readonly int NodeSize;

    /// <summary>Size of a single child entry in bytes.</summary>
    public readonly int ChildEntrySize;

    /// <summary>True if offset fields (SuffixLink, ChildrenHead, ChildNodeOffset) are 64-bit.</summary>
    public readonly bool OffsetIs64Bit;

    // Node field offsets (Start=0, End=4 are invariant across formats)
    /// <summary>Offset of the SuffixLink field within a node.</summary>
    public readonly int OffsetSuffixLink;
    /// <summary>Offset of the DepthFromRoot field within a node.</summary>
    public readonly int OffsetDepth;
    /// <summary>Offset of the LeafCount field within a node.</summary>
    public readonly int OffsetLeafCount;
    /// <summary>Offset of the ChildrenHead field within a node.</summary>
    public readonly int OffsetChildrenHead;
    /// <summary>Offset of the ChildCount field within a node.</summary>
    public readonly int OffsetChildCount;

    // Child entry field offsets (invariant)
    /// <summary>Offset of the Key field within a child entry (always 0).</summary>
    public const int ChildOffsetKey = 0;
    /// <summary>Offset of the ChildNodeOffset field within a child entry (always 4).</summary>
    public const int ChildOffsetNode = 4;

    // ──────────────── Constructor ────────────────

    private NodeLayout(int version, int nodeSize, int childEntrySize, bool offsetIs64Bit,
        int offsetSuffixLink, int offsetDepth, int offsetLeafCount,
        int offsetChildrenHead, int offsetChildCount)
    {
        Version = version;
        NodeSize = nodeSize;
        ChildEntrySize = childEntrySize;
        OffsetIs64Bit = offsetIs64Bit;
        OffsetSuffixLink = offsetSuffixLink;
        OffsetDepth = offsetDepth;
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
            storage.WriteUInt32(position,
                value == PersistentConstants.NULL_OFFSET ? COMPACT_NULL : (uint)value);
        }
    }

    // ──────────────── Version lookup ────────────────

    /// <summary>
    /// Returns the <see cref="NodeLayout"/> for a given format version read from the header.
    /// </summary>
    public static NodeLayout ForVersion(int version) => version switch
    {
        3 => Large,
        4 => Compact,
        _ => throw new InvalidOperationException(
            $"Unsupported storage format version: {version}. Expected 3 (Large) or 4 (Compact).")
    };

    /// <summary>
    /// Returns the <see cref="NodeLayout"/> for a given <see cref="StorageFormat"/> enum.
    /// </summary>
    public static NodeLayout ForFormat(StorageFormat format) => format switch
    {
        StorageFormat.Compact => Compact,
        StorageFormat.Large => Large,
        _ => throw new ArgumentOutOfRangeException(nameof(format), format, "Unknown storage format.")
    };
}
