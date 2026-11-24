using System.Buffers.Binary;
using SuffixTree.Esa.External;
using SuffixTree.Esa.Tests.External;

namespace SuffixTree.Esa.Tests.Large;

/// <summary>
/// The test-side reader of a SESA container of POS_WIDTH 5 (work package S2.1, deliverable C): it parses the bytes by the rules of ADR-005 section 6.1 and by nothing else (no product type decodes; the encoder is not used), and
/// gives the logical arrays back. Every violation is an <see cref="InvalidDataException"/>. A position is read as section 6.1 says: one unaligned 8-byte little-endian load and a mask, which is why the padding is checked.
/// </summary>
internal sealed class Sesa5Image
{
    /// <summary>Section kind 7 (<c>BUCKET</c>).</summary>
    internal const int KindBucket = 7;

    /// <summary>Section kind 8 (<c>ISA</c>).</summary>
    internal const int KindIsa = 8;

    private const ulong Mask40 = 0xFF_FFFF_FFFFUL;

    private Sesa5Image(byte[] file, SesaOracle oracle)
    {
        File = file;
        Oracle = oracle;
    }

    internal byte[] File { get; }

    internal SesaOracle Oracle { get; }

    /// <summary>The sections in the order of the table.</summary>
    internal IReadOnlyList<OracleSection> Sections => Oracle.Sections;

    internal long N => Oracle.TextLen;

    internal int VersionMinor => Oracle.VersionMinor;

    internal OracleSection? Find(int kind) => Sections.FirstOrDefault(s => s.Kind == kind);

    internal OracleSection Get(int kind) => Find(kind) ?? throw new InvalidDataException($"the container has no section of kind {kind}");

    /// <summary>The payload of a section (the padding included).</summary>
    internal ReadOnlySpan<byte> Payload(OracleSection s) => File.AsSpan((int)s.Offset, (int)s.Length);

    /// <summary>The entry <paramref name="index"/> of a section of 5-byte entries that starts at <paramref name="sectionOffset"/>: the 8-byte little-endian load at <c>5 index</c> and the mask.</summary>
    internal static ulong Pos40At(byte[] file, long sectionOffset, long index)
    {
        ReadOnlySpan<byte> load = file.AsSpan((int)(sectionOffset + (5 * index)), 8);
        return BinaryPrimitives.ReadUInt64LittleEndian(load) & Mask40;
    }

    /// <summary>Every entry of a position section.</summary>
    internal ulong[] Positions(int kind)
    {
        OracleSection s = Get(kind);
        var values = new ulong[s.Count];
        for (long i = 0; i < s.Count; i++) values[i] = Pos40At(File, s.Offset, i);
        return values;
    }

    /// <summary>The pairs of LCP_OVERFLOW: both fields are 5-byte positions, the second at the pair offset + 5.</summary>
    internal (ulong Rank, ulong Value)[] Overflow()
    {
        OracleSection? s = Find(SesaFormat.KindLcpOverflow);
        if (s is null) return [];
        var pairs = new (ulong, ulong)[s.Count];
        for (long i = 0; i < s.Count; i++)
        {
            ulong rank = BinaryPrimitives.ReadUInt64LittleEndian(File.AsSpan((int)(s.Offset + (10 * i)), 8)) & Mask40;
            ulong value = BinaryPrimitives.ReadUInt64LittleEndian(File.AsSpan((int)(s.Offset + (10 * i) + 5), 8)) & Mask40;
            pairs[i] = (rank, value);
        }
        return pairs;
    }

    /// <summary>The LCP array of the container, escapes resolved through LCP_OVERFLOW (codecs 1 - 4).</summary>
    internal long[] Lcp()
    {
        OracleSection s = Get(SesaFormat.KindLcp);
        ReadOnlySpan<byte> bytes = Payload(s);
        var lcp = new long[s.Count];
        for (long i = 0; i < s.Count; i++)
        {
            lcp[i] = s.Codec switch
            {
                1 => bytes[(int)i],
                2 => BinaryPrimitives.ReadUInt16LittleEndian(bytes[(int)(2 * i)..]),
                3 => BinaryPrimitives.ReadUInt32LittleEndian(bytes[(int)(4 * i)..]),
                4 => (long)BinaryPrimitives.ReadUInt64LittleEndian(bytes[(int)(8 * i)..]),
                _ => throw new InvalidDataException($"LCP codec {s.Codec}"),
            };
        }
        long escape = s.Codec == 1 ? 255 : s.Codec == 2 ? 65_535 : -1;
        foreach ((ulong rank, ulong value) in Overflow())
        {
            if (lcp[(long)rank] != escape) throw new InvalidDataException($"LCP_OVERFLOW rank {rank} is not an escape entry");
            lcp[(long)rank] = (long)value;
        }
        return lcp;
    }

    /// <summary>Parses and validates a container of POS_WIDTH 5 (section 6.1).</summary>
    internal static Sesa5Image Parse(byte[] file)
    {
        ArgumentNullException.ThrowIfNull(file);
        if (file.Length < 4096) throw new InvalidDataException($"truncated: {file.Length} bytes");
        var o = SesaOracle.FromBytes(file);
        if (!o.Magic.AsSpan().SequenceEqual(SesaFormat.Magic)) throw new InvalidDataException("not a SESA container");
        if (o.VersionMajor != 1) throw new InvalidDataException($"major version {o.VersionMajor}");
        if (o.PosWidth != 5) throw new InvalidDataException($"POS_WIDTH is {o.PosWidth}, not 5");
        if (o.VersionMinor != 1) throw new InvalidDataException($"POS_WIDTH 5 with VERSION_MINOR {o.VersionMinor}: a 1.1 feature in a container that says 1.0");
        if (o.EndianProbe != 0x01020304) throw new InvalidDataException("endian probe");
        if (o.TableOffset != 128 || o.EntrySize != 64 || o.SectionCount is < 1 or > 62) throw new InvalidDataException("section table");
        if (o.HeaderCrc != o.ExpectedHeaderCrc()) throw new InvalidDataException("HEADER_CRC32C");
        if (o.HeaderReserved != 0) throw new InvalidDataException("reserved header word");
        if (o.FileSize != file.Length) throw new InvalidDataException($"FILE_SIZE {o.FileSize} but the file has {file.Length} bytes");
        if (o.TextLen is < 0 or >= 1L << 40) throw new InvalidDataException($"TEXT_LEN {o.TextLen}");
        if (o.Page.AsSpan(128 + (64 * o.SectionCount)).IndexOfAnyExcept((byte)0) >= 0) throw new InvalidDataException("the header page is not zero after the table");

        long end = 4096;
        foreach (OracleSection s in o.Sections.OrderBy(x => x.Offset))
        {
            string what = $"section {s.Kind}";
            if (s.Reserved != 0) throw new InvalidDataException($"{what}: reserved word");
            if (s.Offset < end) throw new InvalidDataException($"{what}: overlaps the one before it");
            if (s.Offset + s.Length > file.Length) throw new InvalidDataException($"{what}: outside the file");
            if (s.Offset % (s.Length >= 1 << 20 ? 4096 : 64) != 0) throw new InvalidDataException($"{what}: alignment");
            if (file.AsSpan((int)end, (int)(s.Offset - end)).IndexOfAnyExcept((byte)0) >= 0) throw new InvalidDataException($"{what}: the gap before it is not zero");
            end = s.Offset + s.Length;
            if (s.HasCrc && s.Crc != Crc32C.Compute(file.AsSpan((int)s.Offset, (int)s.Length))) throw new InvalidDataException($"{what}: CRC-32C");
            if (!s.HasCrc && s.Crc != 0) throw new InvalidDataException($"{what}: CRC without HAS_CRC");
            bool positions = Sesa5Encoder.IsPositionKind(s.Kind), pairs = s.Kind == SesaFormat.KindLcpOverflow;
            if (positions || pairs)
            {
                int width = positions ? 5 : 10;
                if (s.ElemWidth != width) throw new InvalidDataException($"{what}: ELEM_WIDTH {s.ElemWidth}, expected {width}");
                if (s.Length != (width * s.Count) + 3) throw new InvalidDataException($"{what}: LENGTH {s.Length}, expected {(width * s.Count) + 3} (the 3 bytes of padding belong to the payload)");
                if (file.AsSpan((int)(s.Offset + (width * s.Count)), 3).IndexOfAnyExcept((byte)0) >= 0) throw new InvalidDataException($"{what}: the padding is not zero");
            }
            if (s.Kind == SesaFormat.KindChild && s.Codec != 1) throw new InvalidDataException($"CHILD codec {s.Codec} does not exist at POS_WIDTH 5");
        }
        if (end != file.Length) throw new InvalidDataException("the last section does not end at FILE_SIZE");
        return new Sesa5Image(file, o);
    }
}
