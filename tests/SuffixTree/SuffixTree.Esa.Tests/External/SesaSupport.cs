using System.Buffers.Binary;
using System.Runtime.InteropServices;
using System.Security.Cryptography;
using System.Text;
using SuffixTree.Esa.External;

namespace SuffixTree.Esa.Tests.External;

/// <summary>One entry of the section table, decoded by <see cref="SesaOracle"/> straight from the bytes of ADR-003 section 3.</summary>
internal sealed record OracleSection(int Index, int Kind, int Flags, long Offset, long Length, long Count, int ElemWidth, int Codec, long Param0, long Param1, uint Crc, uint Reserved)
{
    internal bool Required => (Flags & 1) != 0;

    internal bool HasCrc => (Flags & 2) != 0;
}

/// <summary>
/// An independent reader of a SESA file, written from the text of ADR-003 section 3 (as amended on 2026-11-18) and nothing else: no type of the product is used to decode (the CRC-32C is the bit-by-bit definition, not
/// a table and not an instruction). Reads the whole file into memory, so it is for the small files of the tests.
/// </summary>
internal sealed class SesaOracle
{
    private SesaOracle(byte[] file)
    {
        File = file;
        Page = file.Length >= 4096 ? file[..4096] : file;
    }

    internal byte[] File { get; }

    internal byte[] Page { get; }

    internal byte[] Magic => Page[..8];

    internal int VersionMajor => BinaryPrimitives.ReadUInt16LittleEndian(Page.AsSpan(8));

    internal int VersionMinor => BinaryPrimitives.ReadUInt16LittleEndian(Page.AsSpan(10));

    internal uint EndianProbe => BinaryPrimitives.ReadUInt32LittleEndian(Page.AsSpan(12));

    internal long TextLen => BinaryPrimitives.ReadInt64LittleEndian(Page.AsSpan(16));

    internal long NodeCount => BinaryPrimitives.ReadInt64LittleEndian(Page.AsSpan(24));

    internal long FileSize => BinaryPrimitives.ReadInt64LittleEndian(Page.AsSpan(32));

    internal long TableOffset => BinaryPrimitives.ReadInt64LittleEndian(Page.AsSpan(40));

    internal int SectionCount => BinaryPrimitives.ReadInt32LittleEndian(Page.AsSpan(48));

    internal int EntrySize => BinaryPrimitives.ReadInt32LittleEndian(Page.AsSpan(52));

    internal long LrsDepth => BinaryPrimitives.ReadInt64LittleEndian(Page.AsSpan(56));

    internal long LrsRank => BinaryPrimitives.ReadInt64LittleEndian(Page.AsSpan(64));

    internal int PosWidth => BinaryPrimitives.ReadInt32LittleEndian(Page.AsSpan(72));

    internal int Sigma => BinaryPrimitives.ReadInt32LittleEndian(Page.AsSpan(76));

    internal byte[] TextSha256 => Page[80..112];

    internal long BuildFlags => BinaryPrimitives.ReadInt64LittleEndian(Page.AsSpan(112));

    internal uint HeaderReserved => BinaryPrimitives.ReadUInt32LittleEndian(Page.AsSpan(120));

    internal uint HeaderCrc => BinaryPrimitives.ReadUInt32LittleEndian(Page.AsSpan(124));

    internal List<OracleSection> Sections
    {
        get
        {
            var list = new List<OracleSection>();
            for (int i = 0; i < SectionCount; i++)
            {
                ReadOnlySpan<byte> e = Page.AsSpan(128 + (64 * i), 64);
                list.Add(new OracleSection(
                    i,
                    BinaryPrimitives.ReadInt32LittleEndian(e),
                    BinaryPrimitives.ReadInt32LittleEndian(e[4..]),
                    BinaryPrimitives.ReadInt64LittleEndian(e[8..]),
                    BinaryPrimitives.ReadInt64LittleEndian(e[16..]),
                    BinaryPrimitives.ReadInt64LittleEndian(e[24..]),
                    BinaryPrimitives.ReadInt32LittleEndian(e[32..]),
                    BinaryPrimitives.ReadInt32LittleEndian(e[36..]),
                    BinaryPrimitives.ReadInt64LittleEndian(e[40..]),
                    BinaryPrimitives.ReadInt64LittleEndian(e[48..]),
                    BinaryPrimitives.ReadUInt32LittleEndian(e[56..]),
                    BinaryPrimitives.ReadUInt32LittleEndian(e[60..])));
            }
            return list;
        }
    }

    internal static SesaOracle Read(string path) => new(System.IO.File.ReadAllBytes(path));

    internal static SesaOracle FromBytes(byte[] file) => new(file);

    internal OracleSection? Find(int kind) => Sections.FirstOrDefault(s => s.Kind == kind);

    internal OracleSection Get(int kind) => Find(kind) ?? throw new InvalidOperationException($"the file has no section of kind {kind}");

    internal byte[] Payload(OracleSection s) => File[(int)s.Offset..(int)(s.Offset + s.Length)];

    internal byte[] Payload(int kind) => Payload(Get(kind));

    /// <summary>The CRC-32C of the definition (polynomial 0x82F63B78 reflected, initial and final xor 0xFFFFFFFF), one bit at a time.</summary>
    internal static uint Crc32cBitwise(ReadOnlySpan<byte> data, uint seed = 0)
    {
        uint crc = ~seed;
        foreach (byte b in data)
        {
            crc ^= b;
            for (int k = 0; k < 8; k++) crc = (crc & 1) != 0 ? (crc >> 1) ^ 0x82F63B78u : crc >> 1;
        }
        return ~crc;
    }

    /// <summary>The header CRC the text of ADR-003 defines: over the bytes <c>[0, 124)</c> and the whole section table (the CRC field is not in either).</summary>
    internal uint ExpectedHeaderCrc()
    {
        uint crc = Crc32cBitwise(Page.AsSpan(0, 124));
        return Crc32cBitwise(Page.AsSpan(128, 64 * SectionCount), crc);
    }

    /// <summary>The element width of the section kinds whose payload is an array of whole elements (the others have 0).</summary>
    private static bool IsArray(OracleSection s) => !(s.Kind == 1 && s.Codec == 4) && s.Kind != 11;

    /// <summary>Everything that ADR-003 sections 3, 4 and 6 (level <c>Header</c> plus the CRCs of level <c>Sections</c>) say about a file, checked from the bytes; the first violation fails the test.</summary>
    internal void AssertConforms(string label, bool psiExpected)
    {
        Assert.That(File.Length, Is.GreaterThanOrEqualTo(4096), $"{label}: file shorter than the header page");
        Assert.That(Magic, Is.EqualTo(new byte[] { 0x89, 0x53, 0x45, 0x53, 0x41, 0x0D, 0x0A, 0x1A }), $"{label}: magic");
        Assert.That((VersionMajor, VersionMinor), Is.EqualTo((1, 0)), $"{label}: version");
        Assert.That(Page.AsSpan(12, 4).ToArray(), Is.EqualTo(new byte[] { 4, 3, 2, 1 }), $"{label}: endian probe bytes");
        Assert.That(TableOffset, Is.EqualTo(128), $"{label}: SECTION_TABLE_OFFSET");
        Assert.That(EntrySize, Is.EqualTo(64), $"{label}: SECTION_ENTRY_SIZE");
        Assert.That(SectionCount, Is.InRange(1, 62), $"{label}: SECTION_COUNT");
        Assert.That(FileSize, Is.EqualTo(File.Length), $"{label}: FILE_SIZE");
        Assert.That(PosWidth, Is.EqualTo(4), $"{label}: POS_WIDTH");
        Assert.That(BuildFlags, Is.EqualTo(1), $"{label}: BUILD_FLAGS (sentinel SA)");
        Assert.That(HeaderReserved, Is.Zero, $"{label}: reserved header word");
        Assert.That(HeaderCrc, Is.EqualTo(ExpectedHeaderCrc()), $"{label}: HEADER_CRC32C");
        Assert.That(Page.AsSpan(128 + (64 * SectionCount)).IndexOfAnyExcept((byte)0), Is.EqualTo(-1), $"{label}: the page after the table is zero");
        List<OracleSection> sections = Sections;
        long end = 4096;
        var covered = new bool[File.Length];
        foreach (OracleSection s in sections)
        {
            string what = $"{label}: section {s.Kind}";
            Assert.That(s.Reserved, Is.Zero, $"{what}: reserved");
            Assert.That(s.Offset % 64, Is.Zero, $"{what}: 64-byte alignment");
            if (s.Length >= (1 << 20)) Assert.That(s.Offset % 4096, Is.Zero, $"{what}: 4096-byte alignment of a section of 1 MiB or more");
            Assert.That(s.Offset, Is.GreaterThanOrEqualTo(4096), $"{what}: inside the header page");
            Assert.That(s.Offset + s.Length, Is.LessThanOrEqualTo(File.Length), $"{what}: outside the file");
            Assert.That(s.Offset, Is.GreaterThanOrEqualTo(end), $"{what}: overlaps the previous section or is out of file order");
            end = s.Offset + s.Length;
            for (long b = s.Offset; b < s.Offset + s.Length; b++) covered[b] = true;
            Assert.That(s.HasCrc, Is.True, $"{what}: HAS_CRC");
            Assert.That(s.Crc, Is.EqualTo(Crc32cBitwise(Payload(s))), $"{what}: CRC32C");
            if (IsArray(s)) Assert.That(s.Count * s.ElemWidth, Is.EqualTo(s.Length), $"{what}: COUNT * ELEM_WIDTH");
        }
        for (long b = 4096; b < File.Length; b++)
        {
            if (!covered[b]) Assert.That(File[b], Is.Zero, $"{label}: gap byte {b} between sections is not zero");
        }
        Assert.That(end, Is.EqualTo(File.Length), $"{label}: the file ends with its last section");
        Assert.That(sections.Select(s => s.Kind).Distinct().Count(), Is.EqualTo(sections.Count), $"{label}: a kind occurs twice");

        OracleSection text = Get(1), alphabet = Get(2), sa = Get(3), lcp = Get(4), child = Get(6), counts = Get(13), info = Get(11);
        Assert.That(text.Required && alphabet.Required && sa.Required && lcp.Required && child.Required, Is.True, $"{label}: required sections carry REQUIRED");
        Assert.That(counts.Required || info.Required, Is.False, $"{label}: COUNTS and BUILD_INFO are written with REQUIRED clear");
        Assert.That((text.Count, text.ElemWidth, text.Codec, text.Param1), Is.EqualTo((TextLen, 0, 4, 0L)), $"{label}: TEXT entry");
        Assert.That(text.Param0, Is.EqualTo(Math.Max(1, Sigma <= 1 ? 0 : (int)Math.Ceiling(Math.Log2(Sigma)))), $"{label}: TEXT b = max(1, ceil(log2 SIGMA))");
        Assert.That(text.Length, Is.EqualTo(((text.Count * text.Param0) + 7) / 8 + 64), $"{label}: TEXT LENGTH");
        Assert.That(1L << (int)text.Param0, Is.GreaterThanOrEqualTo(Sigma), $"{label}: 2^b >= SIGMA");
        Assert.That(alphabet.Count, Is.EqualTo(Sigma), $"{label}: ALPHABET COUNT");
        Assert.That(alphabet.ElemWidth, Is.EqualTo(2), $"{label}: ALPHABET width");
        ushort[] letters = MemoryMarshal.Cast<byte, ushort>(Payload(alphabet)).ToArray();
        Assert.That(letters.OrderBy(l => l).Distinct().Count(), Is.EqualTo(letters.Length), $"{label}: ALPHABET sorted distinct");
        Assert.That(letters.SequenceEqual(letters.OrderBy(l => l)), Is.True, $"{label}: ALPHABET sorted");
        Assert.That((counts.Count, counts.ElemWidth), Is.EqualTo(((long)Sigma, 8)), $"{label}: COUNTS entry");
        Assert.That(MemoryMarshal.Cast<byte, ulong>(Payload(counts)).ToArray().Aggregate(0UL, (a, c) => a + c), Is.EqualTo((ulong)TextLen), $"{label}: COUNTS sum to n");
        Assert.That((sa.Count, sa.ElemWidth, sa.Codec), Is.EqualTo((TextLen + 1, 4, 1)), $"{label}: SA entry");
        Assert.That((child.Count, child.ElemWidth, child.Codec), Is.EqualTo((TextLen + 1, 4, 1)), $"{label}: CHILD entry (layout A, packed)");
        Assert.That(lcp.Count, Is.EqualTo(TextLen + 1), $"{label}: LCP COUNT");
        Assert.That(lcp.Codec, Is.InRange(1, 3), $"{label}: LCP codec");
        Assert.That(lcp.ElemWidth, Is.EqualTo(lcp.Codec switch { 1 => 1, 2 => 2, _ => 4 }), $"{label}: LCP width");
        OracleSection? overflow = Find(5);
        Assert.That(overflow is not null, Is.EqualTo(lcp.Param0 > 0), $"{label}: LCP_OVERFLOW iff PARAM0 > 0");
        if (overflow is not null)
        {
            Assert.That((overflow.Count, overflow.ElemWidth, overflow.Required), Is.EqualTo((lcp.Param0, 8, true)), $"{label}: LCP_OVERFLOW entry");
        }
        OracleSection? psi = Find(12);
        Assert.That(psi is not null, Is.EqualTo(psiExpected), $"{label}: PSI present");
        if (psi is not null)
        {
            Assert.That((psi.Count, psi.ElemWidth, psi.Codec, psi.Required, psi.Param1), Is.EqualTo((TextLen + 1, 4, 1, false, 0L)), $"{label}: PSI entry");
        }
        Assert.That((info.Codec, info.ElemWidth), Is.EqualTo((1, 0)), $"{label}: BUILD_INFO entry");
        string[] lines = Encoding.UTF8.GetString(Payload(info)).Split('\n', StringSplitOptions.RemoveEmptyEntries);
        Assert.That(info.Count, Is.EqualTo(lines.Length), $"{label}: BUILD_INFO COUNT = lines");
        Assert.That(lines.All(l => l.Contains('=', StringComparison.Ordinal)), Is.True, $"{label}: BUILD_INFO key=value lines");
    }
}

/// <summary>One container built from a string and the in-RAM oracles of everything in it.</summary>
internal static class SesaRun
{
    internal sealed record Outcome(SesaBuildResult Result, string Input, string Target, string OutDir, SesaOracle Oracle, LcpChildCase Case);

    /// <summary>Writes the text as one byte per character (every test text is below U+0100) and returns the path.</summary>
    internal static string WriteInput(TempDir dir, string text, string name = "input.txt")
    {
        string path = dir.File(name);
        File.WriteAllBytes(path, Encoding.Latin1.GetBytes(text));
        return path;
    }

    /// <summary>Tiny parameters that rotate with <paramref name="v"/> in every stage, so that every multi-bucket, multi-chunk, spill, window and pass path runs on texts of 10^1-10^4 symbols.</summary>
    internal static SesaBuildOptions Tiny(int v, int threads = 1, bool psi = true) => new()
    {
        Threads = threads,
        IncludePsi = psi,
        Period = SaTiny.Periods[v % SaTiny.Periods.Length],
        ExactPeriod = true,
        SaTuning = SaTiny.Tuning(v),
        LcpTuning = LcpChildRun.TinyTuning(v),
        Ph6Tuning = new SesaPh6Tuning { WindowBytes = 4 + (4 * (v % 7)), SymbolsPerPass = 1 + (v % 5), ReadBytes = 4 + (4 * (v % 9)) },
        TextBlockSymbols = 8 * (1 + (v % 4)),
        TextWindowBytes = 1 + (v * 29 % 120),
        Seed = (ulong)(v + 1),
        ArenaGuards = ArenaOptions.Guarded(0xCD),
        Fsync = false,
    };

    internal static Outcome Build(TempDir dir, string text, SesaBuildOptions options, string tag = "")
    {
        string input = WriteInput(dir, text, $"input{tag}.txt");
        string outDir = dir.Sub($"out{tag}");
        Directory.CreateDirectory(outDir);
        string target = Path.Combine(outDir, "index.sesa");
        SesaBuildResult r = SesaBuilder.Build(input, target, options);
        return new Outcome(r, input, target, outDir, SesaOracle.Read(target), LcpChildCase.Of(text));
    }

    /// <summary>The PSI array from the definition: <c>PSI[i] = ISA[SA[i] + 1]</c>, 0 for the empty suffix (and for the suffix of the last symbol, whose successor is the empty suffix at rank 0).</summary>
    internal static uint[] NaivePsi(ReadOnlySpan<int> sa)
    {
        int rows = sa.Length;
        var isa = new int[rows + 1];
        for (int i = 0; i < rows; i++) isa[sa[i]] = i;
        var psi = new uint[rows];
        for (int i = 1; i < rows; i++) psi[i] = (uint)(sa[i] + 1 <= rows - 1 ? isa[sa[i] + 1] : 0);
        return psi;
    }

    /// <summary>The row of the suffix 0.</summary>
    internal static int PrimaryRow(ReadOnlySpan<int> sa) => sa.IndexOf(0);

    /// <summary>Compares every section of the file with what the in-RAM index of the text holds, and the header with its statistics.</summary>
    internal static void Compare(Outcome o, string label, bool psi)
    {
        LcpChildCase c = o.Case;
        SesaOracle f = o.Oracle;
        f.AssertConforms(label, psi);
        int n = c.N;
        Assert.That(f.TextLen, Is.EqualTo(n), $"{label}: TEXT_LEN");
        Assert.That(f.Sigma, Is.EqualTo(c.Sigma), $"{label}: SIGMA");
        Assert.That(f.NodeCount, Is.EqualTo(c.Statistics.NodeCount), $"{label}: NODE_COUNT");
        Assert.That(f.LrsDepth, Is.EqualTo(c.Statistics.LrsDepth), $"{label}: LRS_DEPTH");
        Assert.That(f.LrsRank, Is.EqualTo(c.Statistics.LrsDepth == 0 ? -1 : c.Statistics.LrsRank), $"{label}: LRS_RANK");
        Assert.That(f.TextSha256, Is.EqualTo(SHA256.HashData(Encoding.Latin1.GetBytes(c.Text))), $"{label}: TEXT_SHA256");
        byte[] packed = TextOracle.Pack(c.Symbols, c.Bits);
        Assert.That(f.Payload(1), Is.EqualTo(packed), $"{label}: TEXT payload");
        Assert.That(f.Get(1).Param0, Is.EqualTo(c.Bits), $"{label}: TEXT bits");
        ushort[] letters = [.. c.Text.Distinct().OrderBy(ch => ch).Select(ch => (ushort)ch)];
        Assert.That(MemoryMarshal.Cast<byte, ushort>(f.Payload(2)).ToArray(), Is.EqualTo(letters), $"{label}: ALPHABET");
        var counts = new ulong[c.Sigma];
        foreach (byte s in c.Symbols) counts[s]++;
        Assert.That(MemoryMarshal.Cast<byte, ulong>(f.Payload(13)).ToArray(), Is.EqualTo(counts), $"{label}: COUNTS");
        Assert.That(MemoryMarshal.Cast<byte, int>(f.Payload(3)).ToArray(), Is.EqualTo(c.Sa), $"{label}: SA");
        OracleSection lcp = f.Get(4);
        Assert.That(lcp.Codec, Is.EqualTo(c.Lcp.Kind switch { LcpKind.U8 => 1, LcpKind.U16 => 2, _ => 3 }), $"{label}: LCP codec");
        Assert.That(lcp.Param0, Is.EqualTo(c.Lcp.OverflowCount), $"{label}: LCP PARAM0 (overflow entries)");
        Assert.That(f.Payload(4).AsSpan().SequenceEqual(c.LcpBaseBytes()), Is.True, $"{label}: LCP base array ({c.Lcp.Kind})");
        if (c.Lcp.OverflowCount > 0) Assert.That(f.Payload(5).AsSpan().SequenceEqual(c.OverflowBytes()), Is.True, $"{label}: LCP_OVERFLOW");
        Assert.That(MemoryMarshal.Cast<byte, uint>(f.Payload(6)).ToArray(), Is.EqualTo(c.Packed), $"{label}: CHILD");
        if (psi)
        {
            OracleSection p = f.Get(12);
            Assert.That(p.Param0, Is.EqualTo(PrimaryRow(c.Sa)), $"{label}: PSI PARAM0 (the row of the suffix 0)");
            Assert.That(MemoryMarshal.Cast<byte, uint>(f.Payload(12)).ToArray(), Is.EqualTo(NaivePsi(c.Sa)), $"{label}: PSI");
        }
        else
        {
            Assert.That(f.Find(12), Is.Null, $"{label}: no PSI section");
        }
        SesaBuildResult r = o.Result;
        Assert.That(r.FileSize, Is.EqualTo(f.File.Length), $"{label}: result FILE_SIZE");
        Assert.That(r.N, Is.EqualTo(n), $"{label}: result n");
        Assert.That(r.HasPsi, Is.EqualTo(psi), $"{label}: result HasPsi");
        Assert.That(r.Sections.Select(s => (s.Kind, s.Offset, s.Length, s.Crc)), Is.EqualTo(f.Sections.Select(s => (s.Kind, s.Offset, s.Length, s.Crc))), $"{label}: result table = file table");
    }

    internal static string Sha(string path) => Convert.ToHexString(SHA256.HashData(File.ReadAllBytes(path)));
}
