using System.Buffers.Binary;
using System.Runtime.InteropServices;
using System.Text;
using SuffixTree.Esa.External;
using SuffixTree.Esa.Tests.External;

namespace SuffixTree.Esa.Tests.Large;

/// <summary>
/// The test-only re-encoder of work package S2.1, deliverable C (ADR-005 section 3.7 layer 3, section 6.1): a SESA container of POS_WIDTH 4 becomes a container of POS_WIDTH 5 with the same content. The position sections
/// (SA, CHILD, PSI, ISA, BUCKET) get 5-byte little-endian entries and 3 bytes of zero padding that belong to the payload; LCP_OVERFLOW pairs get 10 bytes (two 5-byte fields) and the same padding; a tagged CHILD (codec 2)
/// becomes a packed one (codec 1; there is no tagged cell at width 5); every other section is copied byte for byte; the sections are laid out again in file order (section 6.1), every length, offset, CRC-32C and the header CRC
/// are recomputed, the version word says 1.1. The output is a container that a 1.0 reader (the S1 head) refuses at <c>POS_WIDTH</c>, and that <see cref="Sesa5Image"/> reads back to the same logical arrays.
/// </summary>
internal static class Sesa5Encoder
{
    private const int Padding = 3;
    private const ulong Mask40 = 0xFF_FFFF_FFFFUL;
    private const uint ChildValueMask = (1u << 30) - 1;

    /// <summary>The position sections: kinds whose entries are positions of <c>POS_WIDTH</c> bytes.</summary>
    internal static bool IsPositionKind(int kind) => kind is SesaFormat.KindSa or SesaFormat.KindChild or SesaFormat.KindPsi or Sesa5Image.KindIsa or Sesa5Image.KindBucket;

    /// <summary>Re-encodes the bytes of a container of POS_WIDTH 4.</summary>
    /// <exception cref="InvalidDataException">The input is not a container of version 1.0 and POS_WIDTH 4 (the header CRC and FILE_SIZE are checked: the encoder never makes a container out of damage).</exception>
    /// <exception cref="NotSupportedException">The input has a section of a kind that section 6.1 does not describe.</exception>
    internal static byte[] Encode(byte[] container)
    {
        ArgumentNullException.ThrowIfNull(container);
        if (container.Length < SesaFormat.PageBytes) throw new InvalidDataException($"the input has {container.Length} bytes, shorter than the header page");
        SesaOracle o = SesaOracle.FromBytes(container);
        if (!o.Magic.AsSpan().SequenceEqual(SesaFormat.Magic)) throw new InvalidDataException("the input is not a SESA container: the magic differs");
        if ((o.VersionMajor, o.VersionMinor) != (1, 0)) throw new InvalidDataException($"the input has version {o.VersionMajor}.{o.VersionMinor}, the encoder takes 1.0");
        if (o.PosWidth != 4) throw new InvalidDataException($"the input has POS_WIDTH {o.PosWidth}, the encoder takes 4");
        if (o.FileSize != container.Length) throw new InvalidDataException("FILE_SIZE is not the length of the input");
        if (o.HeaderCrc != o.ExpectedHeaderCrc()) throw new InvalidDataException("the header CRC-32C of the input is wrong");

        var entries = new List<SesaSection>();
        var payloads = new List<byte[]>();
        long offset = SesaFormat.PageBytes;
        foreach (OracleSection s in o.Sections.OrderBy(x => x.Offset))
        {
            (byte[] payload, int elemWidth, int codec) = Convert(o, s);
            long at = SesaFormat.AlignedOffset(offset, payload.Length);
            entries.Add(new SesaSection(s.Kind, s.Flags, at, payload.Length, s.Count, elemWidth, codec, s.Param0, s.Param1, s.HasCrc ? Crc32C.Compute(payload) : 0));
            payloads.Add(payload);
            offset = at + payload.Length;
        }

        var header = new SesaHeader
        {
            TextLength = o.TextLen,
            NodeCount = o.NodeCount,
            FileSize = offset,
            LrsDepth = o.LrsDepth,
            LrsRank = o.LrsRank,
            PosWidth = 5,
            Sigma = o.Sigma,
            TextSha256 = o.TextSha256,
            BuildFlags = o.BuildFlags,
        };
        byte[] page = SesaFormat.EncodeHeaderPage(header, entries, withMagic: true);
        BinaryPrimitives.WriteUInt16LittleEndian(page.AsSpan(10), 1);
        SparseContainer.ResealHeader(page);
        var file = new byte[offset];
        page.CopyTo(file, 0);
        for (int i = 0; i < entries.Count; i++) payloads[i].CopyTo(file, entries[i].Offset);
        return file;
    }

    /// <summary>Re-encodes the container at <paramref name="source"/> into <paramref name="target"/> (a path that does not exist).</summary>
    internal static void EncodeFile(string source, string target) => File.WriteAllBytes(target, Encode(File.ReadAllBytes(source)));

    /// <summary>The payload, <c>ELEM_WIDTH</c> and codec of a section at width 5.</summary>
    private static (byte[] Payload, int ElemWidth, int Codec) Convert(SesaOracle o, OracleSection s)
    {
        byte[] payload = o.Payload(s);
        if (s.Kind == SesaFormat.KindLcpOverflow)
        {
            if (s.ElemWidth != 8) throw new InvalidDataException($"LCP_OVERFLOW has ELEM_WIDTH {s.ElemWidth}, a 1.0 container has 8");
            ReadOnlySpan<uint> pairs = MemoryMarshal.Cast<byte, uint>(payload);
            return (PackPositions(pairs), 10, s.Codec); // a pair is two consecutive fields: the flat sequence of 5-byte entries is the layout of section 6.1
        }
        if (IsPositionKind(s.Kind))
        {
            if (s.ElemWidth != 4) throw new InvalidDataException($"section {s.Kind} has ELEM_WIDTH {s.ElemWidth}, a 1.0 container has 4");
            uint[] cells = [.. MemoryMarshal.Cast<byte, uint>(payload)];
            int codec = s.Codec;
            if (s.Kind == SesaFormat.KindChild && s.Codec == SesaFormat.ChildCodecTagged)
            {
                cells = TaggedToPacked(cells);
                codec = SesaFormat.ChildCodecPacked;
            }
            return (PackPositions(cells), 5, codec);
        }
        if (s.Kind == SesaFormat.KindBuildInfo)
        {
            string text = Encoding.UTF8.GetString(payload).Replace("format=SESA 1.0", "format=SESA 1.1", StringComparison.Ordinal).Replace("child=layout A, tagged", "child=layout A, packed", StringComparison.Ordinal);
            return (Encoding.UTF8.GetBytes(text), s.ElemWidth, s.Codec);
        }
        if (s.Kind is SesaFormat.KindText or SesaFormat.KindAlphabet or SesaFormat.KindCounts or SesaFormat.KindLcp) return (payload, s.ElemWidth, s.Codec);
        throw new NotSupportedException($"the encoder knows no layout for the section kind {s.Kind}");
    }

    /// <summary>
    /// A tagged CHILD table as a packed one: kind 0 (empty) is 0; kind 1 (next l-index) with the value <c>k + 1</c> is 0 (the packed table leaves that cell empty); every other cell is its value with the kind bits dropped.
    /// </summary>
    internal static uint[] TaggedToPacked(uint[] tagged)
    {
        var packed = new uint[tagged.Length];
        for (int k = 0; k < tagged.Length; k++)
        {
            uint kind = tagged[k] >> 30, value = tagged[k] & ChildValueMask;
            packed[k] = kind == 0 || (kind == 1 && value == (uint)k + 1) ? 0 : value;
        }
        return packed;
    }

    /// <summary>The 4-byte values of a 1.0 section as 5-byte little-endian entries followed by 3 zero bytes.</summary>
    private static byte[] PackPositions(ReadOnlySpan<uint> values)
    {
        var wide = new ulong[values.Length];
        for (int i = 0; i < values.Length; i++) wide[i] = values[i];
        return PackPositions(wide);
    }

    /// <summary>The values as 5-byte little-endian entries followed by the 3 zero bytes of padding.</summary>
    /// <exception cref="ArgumentOutOfRangeException">A value is 2^40 or more.</exception>
    internal static byte[] PackPositions(ReadOnlySpan<ulong> values)
    {
        var bytes = new byte[(5 * values.Length) + Padding];
        for (int i = 0; i < values.Length; i++)
        {
            if (values[i] > Mask40) throw new ArgumentOutOfRangeException(nameof(values), values[i], "a 5-byte position is below 2^40");
            for (int j = 0; j < 5; j++) bytes[(5 * i) + j] = (byte)(values[i] >> (8 * j));
        }
        return bytes;
    }
}
