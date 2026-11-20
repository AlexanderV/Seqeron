using System.Runtime.CompilerServices;

namespace SuffixTree.Esa.External;

/// <summary>What the header page of a container says, once it has been validated at the <see cref="EsaValidationLevel.Header"/> level.</summary>
/// <param name="Header">The header fields.</param>
/// <param name="Sections">The section table in the order of the file's table.</param>
/// <param name="Text">TEXT (codec 4).</param>
/// <param name="Alphabet">ALPHABET.</param>
/// <param name="Sa">SA.</param>
/// <param name="Lcp">LCP.</param>
/// <param name="Child">CHILD (codec 1, layout A packed).</param>
/// <param name="Overflow">LCP_OVERFLOW when the LCP has escapes.</param>
/// <param name="Psi">PSI when present.</param>
/// <param name="Counts">COUNTS when present.</param>
/// <param name="Info">BUILD_INFO when present.</param>
/// <param name="LcpKind">The width of the LCP entries.</param>
/// <param name="Bits">Bits per packed symbol.</param>
/// <param name="Letters">The alphabet, <c>Letters[symbol]</c> the code unit.</param>
/// <param name="Counts64">The symbol counts of COUNTS (null when the section is absent).</param>
/// <param name="HeaderCrc">The CRC-32C stored in the header.</param>
internal sealed record SesaParsed(
    SesaHeader Header,
    IReadOnlyList<SesaSection> Sections,
    SesaSection Text,
    SesaSection Alphabet,
    SesaSection Sa,
    SesaSection Lcp,
    SesaSection Child,
    SesaSection? Overflow,
    SesaSection? Psi,
    SesaSection? Counts,
    SesaSection? Info,
    LcpKind LcpKind,
    int Bits,
    char[] Letters,
    ulong[]? Counts64,
    uint HeaderCrc);

/// <summary>
/// The header level of the validation (ADR-003 section 6, as amended on 2026-11-19): the header page and the section table of a container are read from memory and checked, in O(sections), without touching a payload. Every
/// violation is an <see cref="InvalidDataException"/> that names the field, except a version above the supported one and a feature of a valid container that this reader does not have (a position width of 8, a codec other than the
/// ones it knows, an unknown section marked REQUIRED), which are <see cref="NotSupportedException"/>. Nothing here can read outside the file length it is given.
/// </summary>
internal static unsafe class SesaHeaderReader
{
    private const long MaxTextLength = uint.MaxValue - 3L;
    private const long MaxBuildInfoBytes = 1 << 20;

    /// <summary>
    /// The CHILD codec against the length of the text (work package 1b.9a): codec 1 (packed) and codec 2 (tagged, ranks in 30 bits) exist. STUB of the red commit.
    /// </summary>
    /// <exception cref="NotSupportedException">Another codec.</exception>
    /// <exception cref="InvalidDataException">Codec 2 for a text whose ranks do not fit 30 bits.</exception>
    internal static void CheckChildCodec(int codec, long textLength) => throw new NotImplementedException($"1b.9a: not written yet ({codec}, {textLength})");

    /// <summary>Parses and validates the header page at <paramref name="p"/> of a file of <paramref name="length"/> bytes.</summary>
    internal static SesaParsed Parse(byte* p, long length)
    {
        if (!BitConverter.IsLittleEndian) throw new NotSupportedException("the SESA container is little endian and is read in place; a big-endian host cannot map it");
        if (length < SesaFormat.PageBytes) throw new InvalidDataException($"truncated: the file has {length} bytes, shorter than the header page of {SesaFormat.PageBytes}");
        if (!new ReadOnlySpan<byte>(p, 8).SequenceEqual(SesaFormat.Magic)) throw new InvalidDataException("not a SESA container: the magic bytes differ");

        int major = U16(p, 8), minor = U16(p, 10);
        if (major > SesaFormat.VersionMajor) throw new NotSupportedException($"the container has format version {major}.{minor}, this reader supports major version {SesaFormat.VersionMajor}");
        if (major < SesaFormat.VersionMajor) throw new InvalidDataException($"the container has format version {major}.{minor}: major version {major} was never written");
        if (U32(p, 12) != SesaFormat.EndianProbe) throw new InvalidDataException("the endian probe is not 01 02 03 04 stored little endian");

        long tableOffset = I64(p, 40);
        int count = I32(p, 48);
        int entrySize = I32(p, 52);
        if (tableOffset != SesaFormat.HeaderBytes) throw new InvalidDataException($"SECTION_TABLE_OFFSET is {tableOffset}, not {SesaFormat.HeaderBytes}");
        if (entrySize != SesaFormat.EntryBytes) throw new InvalidDataException($"SECTION_ENTRY_SIZE is {entrySize}, not {SesaFormat.EntryBytes}");
        if (count is < 1 or > SesaFormat.MaxSections) throw new InvalidDataException($"SECTION_COUNT is {count}: the header page holds 1 to {SesaFormat.MaxSections} sections");

        uint stored = U32(p, SesaFormat.CrcOffset);
        uint crc = Crc32C.Append(Crc32C.Append(0, p, SesaFormat.CrcOffset), p + SesaFormat.HeaderBytes, (long)SesaFormat.EntryBytes * count);
        if (crc != stored) throw new InvalidDataException($"the header CRC-32C is {stored:X8}, the header page and the section table have {crc:X8}");
        if (U32(p, 120) != 0) throw new InvalidDataException("the reserved header word is not zero");
        int tableEnd = SesaFormat.HeaderBytes + (SesaFormat.EntryBytes * count);
        if (new ReadOnlySpan<byte>(p + tableEnd, SesaFormat.PageBytes - tableEnd).IndexOfAnyExcept((byte)0) >= 0) throw new InvalidDataException("the header page is not zero after the section table");

        long n = I64(p, 16);
        long nodeCount = I64(p, 24);
        long fileSize = I64(p, 32);
        long lrsDepth = I64(p, 56);
        long lrsRank = I64(p, 64);
        int posWidth = I32(p, 72);
        int sigma = I32(p, 76);
        long buildFlags = I64(p, 112);
        if (fileSize != length) throw new InvalidDataException($"FILE_SIZE is {fileSize} but the file has {length} bytes");
        if (posWidth == 8) throw new NotSupportedException("POS_WIDTH 8 (64-bit positions) is not supported by this reader");
        if (posWidth != 4) throw new InvalidDataException($"POS_WIDTH is {posWidth}: 4 or 8");
        if (n < 0 || n > MaxTextLength) throw new InvalidDataException($"TEXT_LEN is {n}: a text of at most {MaxTextLength} symbols fits 4-byte positions");
        if ((buildFlags & SesaFormat.BuildFlagSentinelSa) == 0) throw new NotSupportedException("BUILD_FLAGS says the suffix array has no sentinel entry; this reader needs SA[0] = n");
        if (sigma < 0 || sigma > (n == 0 ? 0 : Math.Min(n, 65536))) throw new InvalidDataException($"SIGMA is {sigma} for a text of {n} symbols");
        CheckStatistics(n, nodeCount, lrsDepth, lrsRank);

        var sections = new SesaSection[count];
        var byKind = new SesaSection?[SesaFormat.KindCounts + 1];
        for (int i = 0; i < count; i++)
        {
            SesaSection s = Entry(p + SesaFormat.HeaderBytes + (SesaFormat.EntryBytes * i), i, length);
            sections[i] = s;
            if (s.Kind is >= 1 and <= SesaFormat.KindCounts)
            {
                if (byKind[s.Kind] is not null) throw new InvalidDataException($"the section kind {s.Kind} occurs twice");
                byKind[s.Kind] = s;
            }
            else if ((s.Flags & SesaFormat.FlagRequired) != 0)
            {
                throw new NotSupportedException($"the section kind {s.Kind} is marked REQUIRED and this reader does not know it");
            }
        }
        CheckLayout(sections, length);

        // kinds that the format defines but this reader does not implement: ignored when optional, refused when REQUIRED
        foreach (int kind in new[] { 7, 8, 9, 10 })
        {
            if (byKind[kind] is { } s && (s.Flags & SesaFormat.FlagRequired) != 0) throw new NotSupportedException($"the section kind {kind} is marked REQUIRED and this reader does not implement it");
        }

        SesaSection text = Required(byKind, SesaFormat.KindText, "TEXT");
        SesaSection alphabet = Required(byKind, SesaFormat.KindAlphabet, "ALPHABET");
        SesaSection sa = Required(byKind, SesaFormat.KindSa, "SA");
        SesaSection lcp = Required(byKind, SesaFormat.KindLcp, "LCP");
        SesaSection child = Required(byKind, SesaFormat.KindChild, "CHILD");
        long rows = n + 1;

        if (text.Codec != PackedTextFormat.CodecId) throw new NotSupportedException($"TEXT codec {text.Codec} is not supported (this reader reads codec {PackedTextFormat.CodecId}, the packed text)");
        if (sigma > 1 << PackedTextFormat.MaxBits) throw new InvalidDataException($"SIGMA is {sigma}, more than the {1 << PackedTextFormat.MaxBits} symbols of the packed text");
        int bits = PackedTextFormat.BitsFor(sigma);
        if (text.Param0 != bits || text.Param1 != 0 || text.ElemWidth != 0 || text.Count != n)
            throw new InvalidDataException($"the TEXT entry (COUNT {text.Count}, ELEM_WIDTH {text.ElemWidth}, PARAM0 {text.Param0}, PARAM1 {text.Param1}) does not describe {n} symbols of {bits} bits (SIGMA {sigma})");
        if (text.Length != PackedTextFormat.TotalBytes(n, bits)) throw new InvalidDataException($"TEXT has LENGTH {text.Length}, {n} symbols of {bits} bits need {PackedTextFormat.TotalBytes(n, bits)}");

        CheckArray(alphabet, "ALPHABET", sigma, 2, 1);
        char[] letters = new char[sigma];
        var letterWords = new ReadOnlySpan<ushort>(p + alphabet.Offset, sigma);
        for (int i = 0; i < sigma; i++)
        {
            if (letterWords[i] > byte.MaxValue) throw new InvalidDataException($"ALPHABET letter {i} is {letterWords[i]}: the packed text holds code units below 256 only");
            if (i > 0 && letterWords[i] <= letterWords[i - 1]) throw new InvalidDataException($"ALPHABET is not strictly ascending at letter {i}");
            letters[i] = (char)letterWords[i];
        }

        CheckArray(sa, "SA", rows, 4, 1);
        CheckArray(child, "CHILD", rows, 4, 1);
        CheckRequired(sa, "SA");
        CheckRequired(child, "CHILD");
        CheckRequired(text, "TEXT");
        CheckRequired(alphabet, "ALPHABET");
        CheckRequired(lcp, "LCP");

        if (lcp.Codec is < 1 or > 3) throw new NotSupportedException($"LCP codec {lcp.Codec} is not supported (1: 8-bit, 2: 16-bit, 3: 32-bit entries)");
        LcpKind lcpKind = lcp.Codec switch { 1 => LcpKind.U8, 2 => LcpKind.U16, _ => LcpKind.U32 };
        CheckArray(lcp, "LCP", rows, lcp.Codec switch { 1 => 1, 2 => 2, _ => 4 }, lcp.Codec);
        if (lcp.Param1 != 0) throw new InvalidDataException("LCP PARAM1 is not zero");
        long escapes = lcp.Param0;
        SesaSection? overflow = byKind[SesaFormat.KindLcpOverflow];
        if (escapes < 0 || escapes > rows || (lcpKind == LcpKind.U32 && escapes != 0)) throw new InvalidDataException($"LCP PARAM0 (the number of overflow entries) is {escapes}");
        if ((escapes > 0) != (overflow is not null)) throw new InvalidDataException(escapes > 0 ? "the LCP has overflow entries but there is no LCP_OVERFLOW section" : "there is an LCP_OVERFLOW section but the LCP has no overflow entries");
        if (overflow is not null)
        {
            CheckArray(overflow, "LCP_OVERFLOW", escapes, 2 * posWidth, 1);
            CheckRequired(overflow, "LCP_OVERFLOW");
        }

        SesaSection? psi = byKind[SesaFormat.KindPsi];
        if (psi is not null)
        {
            CheckArray(psi, "PSI", rows, posWidth, 1);
            if (psi.Param0 < 0 || psi.Param0 > n || psi.Param1 != 0) throw new InvalidDataException($"PSI PARAM0 (the primary row) is {psi.Param0} and PARAM1 {psi.Param1} for a text of {n} symbols");
        }

        ulong[]? counts = null;
        SesaSection? countsSection = byKind[SesaFormat.KindCounts];
        if (countsSection is not null)
        {
            CheckArray(countsSection, "COUNTS", sigma, 8, 1);
            counts = new ulong[sigma];
            new ReadOnlySpan<ulong>(p + countsSection.Offset, sigma).CopyTo(counts);
            ulong sum = 0;
            foreach (ulong c in counts)
            {
                if (c == 0) throw new InvalidDataException("COUNTS has a symbol that does not occur: the alphabet lists only the code units of the text");
                if (c > (ulong)n - sum) throw new InvalidDataException($"COUNTS sums to more than the {n} symbols of the text");
                sum += c;
            }
            if (sum != (ulong)n) throw new InvalidDataException($"COUNTS sums to {sum}, not to the {n} symbols of the text");
        }

        SesaSection? info = byKind[SesaFormat.KindBuildInfo];
        if (info is not null && (info.ElemWidth != 0 || info.Length > MaxBuildInfoBytes || info.Count < 0 || info.Count > info.Length))
            throw new InvalidDataException("the BUILD_INFO entry is malformed");

        var header = new SesaHeader
        {
            TextLength = n,
            NodeCount = nodeCount,
            FileSize = fileSize,
            LrsDepth = lrsDepth,
            LrsRank = lrsRank,
            PosWidth = posWidth,
            Sigma = sigma,
            TextSha256 = new ReadOnlySpan<byte>(p + 80, 32).ToArray(),
            BuildFlags = buildFlags,
        };
        return new SesaParsed(header, sections, text, alphabet, sa, lcp, child, overflow, psi, countsSection, info, lcpKind, bits, letters, counts, stored);
    }

    private static void CheckStatistics(long n, long nodeCount, long lrsDepth, long lrsRank)
    {
        if (n == 0)
        {
            if (nodeCount != 1 || lrsDepth != 0 || lrsRank != -1) throw new InvalidDataException("the statistics of an empty text are NODE_COUNT 1, LRS_DEPTH 0, LRS_RANK -1");
            return;
        }
        if (nodeCount < n + 2 || nodeCount > (2 * n) + 1) throw new InvalidDataException($"NODE_COUNT is {nodeCount}: a tree of a text of {n} symbols has {n + 2} to {(2 * n) + 1} nodes");
        if (lrsDepth < 0 || lrsDepth > n - 1) throw new InvalidDataException($"LRS_DEPTH is {lrsDepth} for a text of {n} symbols");
        if (lrsDepth == 0 ? lrsRank != -1 : lrsRank < 0 || lrsRank > n - 1) throw new InvalidDataException($"LRS_RANK is {lrsRank} for LRS_DEPTH {lrsDepth} and a text of {n} symbols");
    }

    private static SesaSection Required(SesaSection?[] byKind, int kind, string name) => byKind[kind] ?? throw new InvalidDataException($"the required section {name} is missing");

    private static void CheckRequired(SesaSection s, string name)
    {
        if ((s.Flags & SesaFormat.FlagRequired) == 0) throw new InvalidDataException($"the section {name} is not marked REQUIRED");
    }

    private static void CheckArray(SesaSection s, string name, long count, int width, int codec)
    {
        if (s.Codec != codec) throw new NotSupportedException($"the {name} codec is {s.Codec}; this reader reads codec {codec}");
        if (s.Count != count || s.ElemWidth != width) throw new InvalidDataException($"the {name} entry has COUNT {s.Count} and ELEM_WIDTH {s.ElemWidth}, expected {count} and {width}");
        if (s.Length != count * width) throw new InvalidDataException($"{name} has LENGTH {s.Length}, COUNT * ELEM_WIDTH is {count * width}");
    }

    private static SesaSection Entry(byte* e, int index, long fileLength)
    {
        if (U32(e, 60) != 0) throw new InvalidDataException($"the reserved word of section entry {index} is not zero");
        int flags = I32(e, 4);
        if ((flags & ~(SesaFormat.FlagRequired | SesaFormat.FlagHasCrc)) != 0) throw new InvalidDataException($"section entry {index} has the undefined flag bits {flags:X}");
        var s = new SesaSection(I32(e, 0), flags, I64(e, 8), I64(e, 16), I64(e, 24), I32(e, 32), I32(e, 36), I64(e, 40), I64(e, 48), U32(e, 56));
        if (s.Offset < SesaFormat.PageBytes) throw new InvalidDataException($"section {s.Kind} starts at offset {s.Offset}, inside the header page");
        if (s.Length < 0 || s.Count < 0 || s.ElemWidth < 0 || s.Length > fileLength || s.Offset > fileLength - s.Length) throw new InvalidDataException($"section {s.Kind} [{s.Offset}, {s.Offset + s.Length}) lies outside the file of {fileLength} bytes");
        long alignment = s.Length >= SesaFormat.LargeSectionBytes ? SesaFormat.LargeAlignment : SesaFormat.SmallAlignment;
        if (s.Offset % alignment != 0) throw new InvalidDataException($"section {s.Kind} starts at {s.Offset}, which is not a multiple of {alignment}");
        if ((s.Flags & SesaFormat.FlagHasCrc) == 0 && s.Crc != 0) throw new InvalidDataException($"section {s.Kind} has no CRC but its CRC field is {s.Crc:X8}");
        return s;
    }

    /// <summary>The sections do not overlap, and the last one ends at the end of the file.</summary>
    private static void CheckLayout(SesaSection[] sections, long length)
    {
        SesaSection[] byOffset = (SesaSection[])sections.Clone();
        Array.Sort(byOffset, static (a, b) => a.Offset != b.Offset ? a.Offset.CompareTo(b.Offset) : a.Length.CompareTo(b.Length)); // no LINQ here: the first open of a process pays the JIT of everything it touches
        long end = SesaFormat.PageBytes;
        foreach (SesaSection s in byOffset)
        {
            if (s.Offset < end) throw new InvalidDataException($"section {s.Kind} at offset {s.Offset} overlaps the section before it, which ends at {end}");
            end = s.Offset + s.Length;
        }
        if (end != length) throw new InvalidDataException($"the last section ends at {end} but the file has {length} bytes: FILE_SIZE must be the end of the last section");
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static int U16(byte* p, int offset) => Unsafe.ReadUnaligned<ushort>(p + offset);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static uint U32(byte* p, int offset) => Unsafe.ReadUnaligned<uint>(p + offset);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static int I32(byte* p, int offset) => Unsafe.ReadUnaligned<int>(p + offset);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static long I64(byte* p, int offset) => Unsafe.ReadUnaligned<long>(p + offset);
}
