using System.Buffers.Binary;
using System.Runtime.InteropServices;
using System.Security.Cryptography;
using System.Text;
using SuffixTree.Esa.External;

namespace SuffixTree.Esa.Tests.External;

/// <summary>
/// Containers for the tests of the reader (work package 1b.7f): <see cref="Image"/> writes a SESA file from the in-RAM oracles of <see cref="LcpChildCase"/> without the end-to-end builder (microseconds for a tiny text, so thousands
/// of files can be read back), <see cref="Reseal"/> recomputes every checksum of a damaged image so that only the structure checks can see the damage, and <see cref="Open"/> opens a path at a level and an access mode.
/// </summary>
internal static class SesaTestFile
{
    /// <summary>The bytes of the container of the text (a text of code units below 256): the layout of <see cref="SesaLayout"/>, the table in file order, the flags of the builder.</summary>
    internal static byte[] Image(LcpChildCase c, bool psi = true, bool counts = true, bool buildInfo = true)
    {
        PackedTextInfo info = TextOracle.Info(c.Symbols, c.Bits, c.Sigma);
        SesaPlan plan = SesaLayout.Plan(info, psi);
        SesaTail tail = SesaLayout.PlaceTail(plan, c.Lcp.Kind, c.Lcp.OverflowCount);
        var file = new byte[tail.FinalLength];
        int[] sa = c.Sa;
        void Put(long offset, ReadOnlySpan<byte> payload) => payload.CopyTo(file.AsSpan((int)offset));

        byte[] text = TextOracle.Pack(c.Symbols, c.Bits);
        ushort[] letters = [.. c.Text.Distinct().OrderBy(ch => ch).Select(ch => (ushort)ch)];
        var symbolCounts = new ulong[c.Sigma];
        foreach (byte s in c.Symbols) symbolCounts[s]++;
        Put(plan.TextOffset, text);
        Put(plan.AlphabetOffset, MemoryMarshal.AsBytes(letters.AsSpan()));
        if (counts) Put(plan.CountsOffset, MemoryMarshal.AsBytes(symbolCounts.AsSpan()));
        Put(plan.SaOffset, MemoryMarshal.AsBytes(sa.AsSpan()));
        Put(plan.ChildOffset, MemoryMarshal.AsBytes(c.Packed.AsSpan()));
        uint[] psiArray = SesaRun.NaivePsi(sa);
        if (psi) Put(plan.PsiOffset, MemoryMarshal.AsBytes(psiArray.AsSpan()));
        string infoText = SesaBuilder.BuildInfoText(info, c.Lcp.Kind, c.Lcp.OverflowCount, psi);
        byte[] infoBytes = Encoding.UTF8.GetBytes(infoText);
        if (buildInfo) Put(plan.BuildInfoOffset, infoBytes);
        byte[] lcpBase = c.LcpBaseBytes();
        Put(plan.LcpOffset, lcpBase);
        byte[] overflow = c.OverflowBytes();
        if (c.Lcp.OverflowCount > 0) Put(tail.OverflowOffset, overflow);

        const int required = SesaFormat.FlagRequired | SesaFormat.FlagHasCrc;
        const int optional = SesaFormat.FlagHasCrc;
        long rows = c.N + 1;
        var sections = new List<SesaSection>
        {
            new(SesaFormat.KindText, required, plan.TextOffset, plan.TextLength, c.N, 0, 4, c.Bits, 0, Crc32C.Compute(text)),
            new(SesaFormat.KindAlphabet, required, plan.AlphabetOffset, plan.AlphabetLength, c.Sigma, 2, 1, 0, 0, Crc32C.Compute(MemoryMarshal.AsBytes(letters.AsSpan()))),
        };
        if (counts) sections.Add(new SesaSection(SesaFormat.KindCounts, optional, plan.CountsOffset, plan.CountsLength, c.Sigma, 8, 1, 0, 0, Crc32C.Compute(MemoryMarshal.AsBytes(symbolCounts.AsSpan()))));
        sections.Add(new SesaSection(SesaFormat.KindSa, required, plan.SaOffset, plan.SaLength, rows, 4, 1, 0, 0, Crc32C.Compute(MemoryMarshal.AsBytes(sa.AsSpan()))));
        sections.Add(new SesaSection(SesaFormat.KindChild, required, plan.ChildOffset, plan.ChildLength, rows, 4, 1, 0, 0, Crc32C.Compute(MemoryMarshal.AsBytes(c.Packed.AsSpan()))));
        if (psi) sections.Add(new SesaSection(SesaFormat.KindPsi, optional, plan.PsiOffset, plan.PsiLength, rows, 4, 1, SesaRun.PrimaryRow(sa), 0, Crc32C.Compute(MemoryMarshal.AsBytes(psiArray.AsSpan()))));
        if (buildInfo) sections.Add(new SesaSection(SesaFormat.KindBuildInfo, optional, plan.BuildInfoOffset, infoBytes.Length, infoText.Count(ch => ch == '\n'), 0, 1, 0, 0, Crc32C.Compute(infoBytes)));
        sections.Add(new SesaSection(SesaFormat.KindLcp, required, plan.LcpOffset, tail.LcpLength, rows, tail.LcpWidth, c.Lcp.Kind switch { LcpKind.U8 => 1, LcpKind.U16 => 2, _ => 3 }, c.Lcp.OverflowCount, 0, Crc32C.Compute(lcpBase)));
        if (c.Lcp.OverflowCount > 0) sections.Add(new SesaSection(SesaFormat.KindLcpOverflow, required, tail.OverflowOffset, tail.OverflowLength, c.Lcp.OverflowCount, 8, 1, 0, 0, Crc32C.Compute(overflow)));

        var header = new SesaHeader
        {
            TextLength = c.N,
            NodeCount = c.Statistics.NodeCount,
            FileSize = file.Length,
            LrsDepth = c.Statistics.LrsDepth,
            LrsRank = c.Statistics.LrsDepth == 0 ? -1 : c.Statistics.LrsRank,
            PosWidth = 4,
            Sigma = c.Sigma,
            TextSha256 = SHA256.HashData(Encoding.Latin1.GetBytes(c.Text)),
        };
        SesaFormat.EncodeHeaderPage(header, sections, withMagic: true).CopyTo(file.AsSpan());
        return file;
    }

    /// <summary>Writes the container of <paramref name="text"/> and returns its path.</summary>
    internal static string Write(TempDir dir, string text, string name = "t.sesa", bool psi = true)
    {
        string path = dir.File(name);
        File.WriteAllBytes(path, Image(LcpChildCase.Of(text), psi));
        return path;
    }

    /// <summary>Recomputes the CRC-32C of every section that carries one and then the header CRC, in place: the damage that was done to the body is now invisible to the checksums.</summary>
    internal static void Reseal(byte[] image)
    {
        SesaOracle o = SesaOracle.FromBytes(image);
        foreach (OracleSection s in o.Sections)
        {
            if (!s.HasCrc) continue;
            uint crc = Crc32C.Compute(image.AsSpan((int)s.Offset, (int)s.Length));
            BinaryPrimitives.WriteUInt32LittleEndian(image.AsSpan(128 + (64 * s.Index) + 56), crc);
        }
        ResealHeader(image);
    }

    /// <summary>Recomputes the header CRC only: the CRC-32C of bytes [0, 124) followed by the section table (the four CRC bytes are in neither).</summary>
    internal static void ResealHeader(byte[] image)
    {
        int count = BinaryPrimitives.ReadInt32LittleEndian(image.AsSpan(48));
        uint crc = Crc32C.Append(Crc32C.Append(0, image.AsSpan(0, 124)), image.AsSpan(128, 64 * count));
        BinaryPrimitives.WriteUInt32LittleEndian(image.AsSpan(124), crc);
    }

    internal static OracleSection Section(byte[] image, int kind) => SesaOracle.FromBytes(image).Get(kind);

    /// <summary>Opens the container (the caller disposes the tree).</summary>
    internal static EsaSuffixTree Open(string path, EsaValidationLevel level = EsaValidationLevel.Header, EsaAccessMode access = EsaAccessMode.Checked, bool sha = false)
        => EsaSuffixTreeFactory.Open(path, new EsaOpenOptions { Validation = level, Access = access, VerifyTextSha256 = sha });

    /// <summary>Writes the image to a file of its own and opens it.</summary>
    internal static EsaSuffixTree OpenImage(TempDir dir, byte[] image, EsaValidationLevel level = EsaValidationLevel.Header, EsaAccessMode access = EsaAccessMode.Checked, bool sha = false, string name = "x.sesa")
    {
        string path = dir.File(name);
        File.WriteAllBytes(path, image);
        return Open(path, level, access, sha);
    }

    /// <summary>The tree of a text through the file: written, opened (mapped), and the path removed at once (a mapping outlives the name on Linux), so that thousands of trees leave nothing behind.</summary>
    internal static EsaSuffixTree FileTree(string text, bool psi = true, EsaAccessMode access = EsaAccessMode.Checked, EsaValidationLevel level = EsaValidationLevel.Header)
    {
        string path = Path.Combine(Path.GetTempPath(), "esa-rd-" + Guid.NewGuid().ToString("N") + ".sesa");
        File.WriteAllBytes(path, Image(LcpChildCase.Of(text), psi));
        try
        {
            return Open(path, level, access);
        }
        finally
        {
            File.Delete(path);
        }
    }

    /// <summary>True when every code unit of the text fits a byte (what the container's packed text can hold).</summary>
    internal static bool IsByteText(string text) => text.All(c => c < 256);
}
