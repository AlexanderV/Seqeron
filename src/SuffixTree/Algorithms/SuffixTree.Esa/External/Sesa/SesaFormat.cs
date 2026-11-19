using System.Buffers.Binary;

namespace SuffixTree.Esa.External;

/// <summary>One entry of the section table (ADR-003 section 3, 64 bytes).</summary>
/// <param name="Kind">The section kind (ADR-003 section 4).</param>
/// <param name="Flags">Bit 0 <see cref="SesaFormat.FlagRequired"/>, bit 1 <see cref="SesaFormat.FlagHasCrc"/>.</param>
/// <param name="Offset">Absolute file offset of the payload.</param>
/// <param name="Length">Payload bytes (padding excluded).</param>
/// <param name="Count">Number of elements.</param>
/// <param name="ElemWidth">Bytes per element (0 for a payload that is not an array).</param>
/// <param name="Codec">Kind specific.</param>
/// <param name="Param0">Kind specific.</param>
/// <param name="Param1">Kind specific.</param>
/// <param name="Crc">CRC-32C of the payload when <see cref="SesaFormat.FlagHasCrc"/> is set.</param>
internal sealed record SesaSection(int Kind, int Flags, long Offset, long Length, long Count, int ElemWidth, int Codec, long Param0, long Param1, uint Crc);

/// <summary>The header of a SESA file (ADR-003 section 3, 128 bytes) except the magic and the CRC, which the encoder adds.</summary>
internal sealed record SesaHeader
{
    internal long TextLength { get; init; }

    internal long NodeCount { get; init; }

    internal long FileSize { get; init; }

    internal long LrsDepth { get; init; }

    /// <summary>The rank of the first occurrence of the longest repeated substring, -1 when there is none.</summary>
    internal long LrsRank { get; init; } = -1;

    internal int PosWidth { get; init; } = 4;

    internal int Sigma { get; init; }

    /// <summary>The 32 bytes of the SHA-256 of the code units as bytes.</summary>
    internal byte[] TextSha256 { get; init; } = new byte[32];

    internal long BuildFlags { get; init; } = SesaFormat.BuildFlagSentinelSa;
}

/// <summary>The constants and the encoder of the SESA header page (ADR-003 sections 2-4 as amended on 2026-11-18).</summary>
internal static class SesaFormat
{
    /// <summary>Bytes of the header.</summary>
    internal const int HeaderBytes = 128;

    /// <summary>Bytes of the header page: the header, the section table and zero padding.</summary>
    internal const int PageBytes = 4096;

    /// <summary>Bytes of one section entry.</summary>
    internal const int EntryBytes = 64;

    /// <summary>Entries the page holds.</summary>
    internal const int MaxSections = 62;

    internal const int VersionMajor = 1;

    internal const int VersionMinor = 0;

    /// <summary>Stored little endian: the bytes <c>04 03 02 01</c>.</summary>
    internal const uint EndianProbe = 0x01020304;

    /// <summary>Header offset of <c>HEADER_CRC32C</c>.</summary>
    internal const int CrcOffset = 124;

    /// <summary>Sections of at least this many bytes start at a multiple of <see cref="LargeAlignment"/>.</summary>
    internal const long LargeSectionBytes = 1 << 20;

    internal const int SmallAlignment = 64;

    internal const int LargeAlignment = 4096;

    internal const int KindText = 1;
    internal const int KindAlphabet = 2;
    internal const int KindSa = 3;
    internal const int KindLcp = 4;
    internal const int KindLcpOverflow = 5;
    internal const int KindChild = 6;
    internal const int KindBuildInfo = 11;
    internal const int KindPsi = 12;
    internal const int KindCounts = 13;

    internal const int FlagRequired = 1;
    internal const int FlagHasCrc = 2;

    /// <summary>Bit 0 of <c>BUILD_FLAGS</c>: the suffix array holds the sentinel suffix first (<c>SA[0] = n</c>).</summary>
    internal const long BuildFlagSentinelSa = 1;

    /// <summary>The magic: <c>89 'S' 'E' 'S' 'A' CR LF SUB</c>.</summary>
    internal static ReadOnlySpan<byte> Magic => [0x89, 0x53, 0x45, 0x53, 0x41, 0x0D, 0x0A, 0x1A];

    /// <summary>The next offset at or after <paramref name="offset"/> where a section of <paramref name="length"/> bytes may start (64-byte aligned; 4096-byte aligned from 1 MiB).</summary>
    internal static long AlignedOffset(long offset, long length)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(offset);
        ArgumentOutOfRangeException.ThrowIfNegative(length);
        long a = length >= LargeSectionBytes ? LargeAlignment : SmallAlignment;
        return (offset + a - 1) / a * a;
    }

    /// <summary>
    /// The header page: the header, the section table in the order given, zero padding to <see cref="PageBytes"/>. The CRC covers the bytes <c>[0, 124)</c> and the table, computed with the magic in place;
    /// with <paramref name="withMagic"/> false the first eight bytes of the returned page are zero (the state in which the page is written first), the CRC is the same.
    /// </summary>
    internal static byte[] EncodeHeaderPage(SesaHeader header, IReadOnlyList<SesaSection> sections, bool withMagic)
    {
        ArgumentNullException.ThrowIfNull(header);
        ArgumentNullException.ThrowIfNull(sections);
        if (sections.Count > MaxSections) throw new ArgumentException($"the header page holds {MaxSections} sections, {sections.Count} were given", nameof(sections));
        if (header.TextSha256 is not { Length: 32 }) throw new ArgumentException("TEXT_SHA256 is 32 bytes", nameof(header));
        var page = new byte[PageBytes];
        Magic.CopyTo(page);
        BinaryPrimitives.WriteUInt16LittleEndian(page.AsSpan(8), VersionMajor);
        BinaryPrimitives.WriteUInt16LittleEndian(page.AsSpan(10), VersionMinor);
        BinaryPrimitives.WriteUInt32LittleEndian(page.AsSpan(12), EndianProbe);
        BinaryPrimitives.WriteInt64LittleEndian(page.AsSpan(16), header.TextLength);
        BinaryPrimitives.WriteInt64LittleEndian(page.AsSpan(24), header.NodeCount);
        BinaryPrimitives.WriteInt64LittleEndian(page.AsSpan(32), header.FileSize);
        BinaryPrimitives.WriteInt64LittleEndian(page.AsSpan(40), HeaderBytes);
        BinaryPrimitives.WriteInt32LittleEndian(page.AsSpan(48), sections.Count);
        BinaryPrimitives.WriteInt32LittleEndian(page.AsSpan(52), EntryBytes);
        BinaryPrimitives.WriteInt64LittleEndian(page.AsSpan(56), header.LrsDepth);
        BinaryPrimitives.WriteInt64LittleEndian(page.AsSpan(64), header.LrsRank);
        BinaryPrimitives.WriteInt32LittleEndian(page.AsSpan(72), header.PosWidth);
        BinaryPrimitives.WriteInt32LittleEndian(page.AsSpan(76), header.Sigma);
        header.TextSha256.CopyTo(page.AsSpan(80));
        BinaryPrimitives.WriteInt64LittleEndian(page.AsSpan(112), header.BuildFlags);
        for (int i = 0; i < sections.Count; i++)
        {
            SesaSection s = sections[i];
            Span<byte> e = page.AsSpan(HeaderBytes + (EntryBytes * i), EntryBytes);
            BinaryPrimitives.WriteInt32LittleEndian(e, s.Kind);
            BinaryPrimitives.WriteInt32LittleEndian(e[4..], s.Flags);
            BinaryPrimitives.WriteInt64LittleEndian(e[8..], s.Offset);
            BinaryPrimitives.WriteInt64LittleEndian(e[16..], s.Length);
            BinaryPrimitives.WriteInt64LittleEndian(e[24..], s.Count);
            BinaryPrimitives.WriteInt32LittleEndian(e[32..], s.ElemWidth);
            BinaryPrimitives.WriteInt32LittleEndian(e[36..], s.Codec);
            BinaryPrimitives.WriteInt64LittleEndian(e[40..], s.Param0);
            BinaryPrimitives.WriteInt64LittleEndian(e[48..], s.Param1);
            BinaryPrimitives.WriteUInt32LittleEndian(e[56..], s.Crc);
        }
        uint crc = Crc32C.Append(Crc32C.Append(0, page.AsSpan(0, CrcOffset)), page.AsSpan(HeaderBytes, EntryBytes * sections.Count));
        BinaryPrimitives.WriteUInt32LittleEndian(page.AsSpan(CrcOffset), crc);
        if (!withMagic) page.AsSpan(0, Magic.Length).Clear();
        return page;
    }
}
