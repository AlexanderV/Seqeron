// SPDX-License-Identifier: LGPL-2.1-or-later
// Port of the EMBOSS 6.6.0 AJAX / NUCLEUS libraries (library source files carry the GNU Lesser General Public License, version 2.1 or later).
// Part of Seqeron.Genomics.Emboss.SequenceFormats, licensed under the GNU Lesser General Public License, version 2.1 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
#pragma warning disable S1121, S3358, S2325, S1066, S3776, S107 // the writer keeps the control flow of ajseqwrite.c / ajseqbam.c

using System.Text;

namespace Seqeron.Genomics.Emboss.SequenceFormats;

using Seqeron.Genomics.Infrastructure;

// seqWriteBam (ajseqwrite.c) and the BGZF output of ajseqbam.c (ajSeqBamHeaderWrite, ajSeqBamWrite, ajSeqBamBgzfWrite / Flush / Close, bamDeflateBlock) as the installed library writes them: the BAM
// bytes of the sequences are kept until the output is closed; the blocks of 64 KiB are deflated by zlib (EmbossZlibDeflate). A5-8c3h, B09 WP59, F272.
public static partial class EmbossSeqWrite
{
    private const int BgzfBlockSize = 64 * 1024, BgzfHeaderLength = 18, BgzfFooterLength = 8;

    // bam_nt16_table
    private static readonly byte[] BamNt16 = BuildNt16();

    private static byte[] BuildNt16()
    {
        var t = new byte[256];
        Array.Fill(t, (byte)15);
        string up = "=ACMGRSVTWYHKDBN";
        for (int i = 0; i < up.Length; i++)
        {
            t[up[i]] = (byte)i;
            t[char.ToLowerInvariant(up[i])] = (byte)i;
        }

        t['='] = 0;
        return t;
    }

    private static void BamInt(List<byte> b, int v)
    {
        b.Add((byte)v);
        b.Add((byte)(v >> 8));
        b.Add((byte)(v >> 16));
        b.Add((byte)(v >> 24));
    }

    // seqWriteBam: the first sequence of the run makes the file (header "@HD VN:1.3 SO:unsorted"), every sequence is an unmapped record
    private static void WriteBam(EmbossSeqWriteRun run, string name, string seq, float[]? accuracy)
    {
        if (run.HeaderCount == 0)
        {
            run.Bam = [];
            const string text = "@HD\tVN:1.3\tSO:unsorted\n";
            run.Bam.AddRange("BAM\u0001"u8);
            BamInt(run.Bam, text.Length);
            run.Bam.AddRange(Encoding.Latin1.GetBytes(text));
            BamInt(run.Bam, 0);
        }
        else if (run.Count == 0 || run.Bam is null)
        {
            // -ossingle: the later sequences go to the BGZF file of the first one, whose file the library closed
            throw new InvalidOperationException("A BAM file of -ossingle: the later sequences are written to a file the library has closed (the reference program reads and writes freed memory).");
        }

        var b = run.Bam;
        int ilen = seq.Length;
        int lQname = 1 + name.Length;
        var qual = new byte[ilen];
        if (accuracy is not null)
        {
            for (int i = 0; i < ilen; i++)
                qual[i] = unchecked((byte)(int)AccuracyAt(accuracy, i));
        }
        else
        {
            Array.Fill(qual, (byte)('"' - 33));
        }

        int slen = (ilen + 1) / 2;
        int dataLen = lQname + slen + ilen;
        BamInt(b, dataLen + 32);                         // block_len
        BamInt(b, -1);                                   // tid
        BamInt(b, -1);                                   // pos
        BamInt(b, lQname);                               // bin 0, qual 0, l_qname
        BamInt(b, 0x0004 << 16);                         // flag 4 (unmapped), n_cigar 0
        BamInt(b, ilen);                                 // l_qseq
        BamInt(b, -1);                                   // mtid
        BamInt(b, 0);                                    // mpos
        BamInt(b, 0);                                    // isize
        foreach (char c in name)
            b.Add((byte)c);
        b.Add(0);
        var packed = new byte[slen];
        for (int i = 0; i < ilen; i++)
            packed[i / 2] |= (byte)(BamNt16[seq[i] & 0xFF] << (4 * (1 - i % 2)));
        b.AddRange(packed);
        b.AddRange(qual);
    }

    private static uint[]? _crcTable;

    private static uint Crc32(ReadOnlySpan<byte> data)
    {
        if (_crcTable is null)
        {
            var t = new uint[256];
            for (uint n = 0; n < 256; n++)
            {
                uint c = n;
                for (int k = 0; k < 8; k++)
                    c = (c & 1) != 0 ? 0xEDB88320u ^ (c >> 1) : c >> 1;
                t[n] = c;
            }

            _crcTable = t;
        }

        uint crc = 0xFFFFFFFFu;
        foreach (byte x in data)
            crc = _crcTable[(crc ^ x) & 0xFF] ^ (crc >> 8);
        return crc ^ 0xFFFFFFFFu;
    }

    // ajSeqBamBgzfClose: the saved bytes are deflated in blocks of 64 KiB, then an empty block ends the file
    public static string? FinishBam(EmbossSeqWriteRun run)
    {
        if (run.Bam is null)
            return null;
        byte[] all = [.. run.Bam];
        var o = new StringBuilder();
        int pos = 0, forced = 0;
        while (pos < all.Length)
        {
            int want = forced > 0 ? forced : Math.Min(BgzfBlockSize, all.Length - pos);
            int len = want;
            byte[] comp;
            while (true)
            {
                comp = EmbossZlibDeflate.Compress(all.AsSpan(pos, len));
                if (comp.Length < BgzfBlockSize - BgzfHeaderLength - BgzfFooterLength)
                    break;
                len -= 1024;                          // bamDeflateBlock: the input does not compress enough, the amount is reduced until it fits
                if (len <= 0)
                    throw new InvalidOperationException("bgzf: input reduction failed");
            }

            AppendBgzfBlock(o, comp, Crc32(all.AsSpan(pos, len)), len);
            forced = want - len;                      // the rest of the block is flushed as a block of its own
            pos += len;
        }

        AppendBgzfBlock(o, EmbossZlibDeflate.Compress([]), 0, 0);
        run.Bam = null;
        return o.ToString();
    }

    private static void AppendBgzfBlock(StringBuilder o, byte[] deflated, uint crc, int inputLength)
    {
        int total = deflated.Length + BgzfHeaderLength + BgzfFooterLength;
        byte[] head = [31, 139, 8, 4, 0, 0, 0, 0, 0, 255, 6, 0, 66, 67, 2, 0, (byte)(total - 1), (byte)((total - 1) >> 8)];
        foreach (byte x in head)
            o.Append((char)x);
        foreach (byte x in deflated)
            o.Append((char)x);
        for (int i = 0; i < 4; i++)
            o.Append((char)((crc >> (8 * i)) & 0xFF));
        for (int i = 0; i < 4; i++)
            o.Append((char)(((uint)inputLength >> (8 * i)) & 0xFF));
    }
}
