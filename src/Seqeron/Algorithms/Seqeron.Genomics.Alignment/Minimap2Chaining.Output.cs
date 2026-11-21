using System.Globalization;
using System.Text;

#pragma warning disable S127, S3358

namespace Seqeron.Genomics.Alignment;

/// <summary>A read of a query file: the name (the first word of the header), the sequence (U converted to T as the program does), the base qualities of a FASTQ record and the comment (the rest of the header).</summary>
/// <param name="Name">The read name.</param>
/// <param name="Sequence">The read sequence.</param>
/// <param name="Quality">The quality string, null for FASTA.</param>
/// <param name="Comment">The header text after the name, null when there is none.</param>
public sealed record MinimapQuery(string Name, string Sequence, string? Quality = null, string? Comment = null);

/// <summary>How the <c>cs</c> tag is written (<c>--cs[=short|long]</c>).</summary>
public enum MinimapCsMode
{
    /// <summary>No <c>cs</c> tag.</summary>
    None,
    /// <summary><c>--cs</c> / <c>--cs=short</c>: matches as <c>:N</c>.</summary>
    Short,
    /// <summary><c>--cs=long</c> / <c>-S</c>: matches as <c>=ACGT</c>.</summary>
    Long,
}

/// <summary>The output options of <c>minimap2</c> (<c>-a -c --cs --MD -Y -L -R -Q -y --paf-no-hit --sam-hit-only --secondary=no --secondary-seq</c>).</summary>
public sealed record MinimapOutputOptions
{
    /// <summary><c>-a</c>: SAM instead of PAF.</summary>
    public bool Sam { get; init; }

    /// <summary><c>-c</c> (<c>MM_F_OUT_CG</c>): the <c>cg</c> tag in PAF.</summary>
    public bool CigarTag { get; init; }

    /// <summary><c>--cs[=short|long]</c> / <c>-S</c>.</summary>
    public MinimapCsMode Cs { get; init; }

    /// <summary><c>--MD</c>: the <c>MD</c> tag.</summary>
    public bool Md { get; init; }

    /// <summary><c>-Y</c>: soft-clip supplementary and secondary alignments (and write their full sequence).</summary>
    public bool SoftClip { get; init; }

    /// <summary><c>-L</c>: a CIGAR with more than 65535 operations goes to the <c>CG</c> tag.</summary>
    public bool LongCigar { get; init; }

    /// <summary><c>--secondary-seq</c>: write the (hard-clipped) sequence of secondary alignments.</summary>
    public bool SecondarySeq { get; init; }

    /// <summary><c>-Q</c>: no base qualities in SAM.</summary>
    public bool NoQuality { get; init; }

    /// <summary><c>-y</c>: append the comment of the read.</summary>
    public bool CopyComment { get; init; }

    /// <summary><c>--paf-no-hit</c>: a PAF line for a read without a hit.</summary>
    public bool PafNoHit { get; init; }

    /// <summary><c>--sam-hit-only</c>: no SAM line for a read without a hit.</summary>
    public bool SamHitOnly { get; init; }

    /// <summary><c>--secondary=no</c> (<c>MM_F_NO_PRINT_2ND</c>): secondary hits are not printed.</summary>
    public bool NoSecondary { get; init; }

    /// <summary><c>-R</c>: the read group line (<c>@RG\tID:foo</c>); the tag <c>RG:Z:foo</c> is added to every record.</summary>
    public string? ReadGroup { get; init; }
}

public static partial class Minimap2Chaining
{
    private static readonly byte[] CompTable = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,84,86,71,72,69,70,67,68,73,74,77,76,75,78,79,80,81,89,83,65,65,66,87,88,82,90,91,92,93,94,95,96,116,118,103,104,101,102,99,100,105,106,109,108,107,110,111,112,113,121,115,97,97,98,119,120,114,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,204,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,237,238,239,240,241,242,243,244,245,246,247,248,249,250,251,252,253,254,255];

    // sam_write_sq: the sequence as is, or reversed (and complemented with the IUPAC table when comp)
    private static string SamSq(string s, bool rev, bool comp)
    {
        if (!rev)
            return s;
        var c = new char[s.Length];
        for (int i = 0; i < c.Length; i++)
        {
            int ch = s[s.Length - 1 - i];
            c[i] = (char)(ch < 128 && comp ? CompTable[ch] : ch);
        }

        return new string(c);
    }

    private static int QnameLength(string s) => s.Length >= 3 && s[^1] >= '0' && s[^1] <= '9' && s[^2] == '/' ? s.Length - 2 : s.Length;

    /// <summary><c>mm_qname_same</c>: the names are equal after the suffix <c>/N</c> is removed.</summary>
    internal static bool QnameSame(string a, string b)
    {
        int l1 = QnameLength(a), l2 = QnameLength(b);
        return l1 == l2 && string.CompareOrdinal(a, 0, b, 0, l1) == 0;
    }

    private static string C4(double v) => v == 0.0 ? "0" : CFormat.Fixed(v, 4);

    /// <summary>The read group id and line of <c>-R</c> (<c>sam_write_rg_line</c>): the line must start with <c>@RG</c>, contain no TAB (use the two characters <c>\\t</c>) and have an <c>ID:</c> of at most 255 characters; null when the option is invalid (the program prints an error and stops).</summary>
    public static (string Id, string Line)? ParseReadGroup(string? readGroup) => ParseReadGroupChecked(readGroup, out _);

    /// <summary><see cref="ParseReadGroup"/> with the error text of the program (<c>[ERROR] ...</c>, without line feed) when the line is invalid.</summary>
    public static (string Id, string Line)? ParseReadGroupChecked(string? readGroup, out string? error)
    {
        error = null;
        if (readGroup is null)
            return ("", "");
        if (!readGroup.StartsWith("@RG", StringComparison.Ordinal))
        {
            error = "[ERROR] the read group line is not started with @RG";
            return null;
        }

        if (readGroup.Contains('\t'))
        {
            error = "[ERROR] the read group line contained literal <tab> characters -- replace with escaped tabs: \\t";
            return null;
        }

        var sb = new StringBuilder();
        for (int i = 0; i < readGroup.Length; i++)
        {
            if (readGroup[i] == '\\')
            {
                i++;
                if (i < readGroup.Length && readGroup[i] == 't')
                    sb.Append('\t');
                else if (i < readGroup.Length && readGroup[i] == '\\')
                    sb.Append('\\');
            }
            else
            {
                sb.Append(readGroup[i]);
            }
        }

        string line = sb.ToString();
        int p = line.IndexOf("\tID:", StringComparison.Ordinal);
        if (p < 0)
        {
            error = "[ERROR] no ID within the read group line";
            return null;
        }

        p += 4;
        int q = p;
        while (q < line.Length && line[q] != '\t' && line[q] != '\n')
            q++;
        if (q - p + 1 > 256)
        {
            error = "[ERROR] @RG:ID is longer than 255 characters";
            return null;
        }

        return (line[p..q], line);
    }

    /// <summary>The SAM header of <c>minimap2 -a</c> (<c>mm_write_sam_hdr</c>): <c>@HD</c>, the <c>@SQ</c> lines, the <c>-R</c> line and <c>@PG</c> (<c>VN:2.26-r1175</c> and the command line <c>CL:minimap2 ARGS</c> when arguments are given).</summary>
    public static string FormatSamHeader(IReadOnlyList<string?> targetNames, IReadOnlyList<int> targetLengths, string? readGroup = null, IReadOnlyList<string>? commandLine = null, bool withSequences = true)
    {
        var sb = new StringBuilder();
        sb.Append("@HD\tVN:1.6\tSO:unsorted\tGO:query\n");
        if (withSequences)
        {
            for (int i = 0; i < targetNames.Count; i++)
                sb.Append("@SQ\tSN:").Append(targetNames[i] ?? i.ToString(CultureInfo.InvariantCulture)).Append("\tLN:").Append(targetLengths[i].ToString(CultureInfo.InvariantCulture)).Append('\n');
        }

        var rg = ParseReadGroup(readGroup);
        if (rg is not null && readGroup is not null)
            sb.Append(rg.Value.Line).Append('\n');
        sb.Append("@PG\tID:minimap2\tPN:minimap2\tVN:2.26-r1175");
        if (commandLine is { Count: > 0 })
            sb.Append("\tCL:minimap2 ").Append(string.Join(' ', commandLine));
        sb.Append('\n');
        return sb.ToString();
    }

    private static string NtChars => "ACGTN";

    // the target / query bases of an aligned hit as write_cs_or_MD builds them
    private static (byte[] Q, byte[] T) AlignedBases(MinimapHit h, string query, IReadOnlyList<string> targets, bool qstrand)
    {
        int ql = h.QueryEnd - h.QueryStart, tl = h.TargetEnd - h.TargetStart;
        var q = new byte[ql];
        var t = new byte[tl];
        string tgt = targets[h.TargetId];
        if (qstrand)
        {
            var tt = h.Reverse ? MmIndexSeqRev(tgt, h.TargetStart, h.TargetEnd) : MmIndexSeq(tgt, h.TargetStart, h.TargetEnd);
            Array.Copy(tt, t, tl);
            for (int i = h.QueryStart; i < h.QueryEnd; i++)
                q[i - h.QueryStart] = (byte)Nt4(query[i]);
        }
        else
        {
            var tt = MmIndexSeq(tgt, h.TargetStart, h.TargetEnd);
            Array.Copy(tt, t, tl);
            if (!h.Reverse)
            {
                for (int i = h.QueryStart; i < h.QueryEnd; i++)
                    q[i - h.QueryStart] = (byte)Nt4(query[i]);
            }
            else
            {
                for (int i = h.QueryStart; i < h.QueryEnd; i++)
                {
                    int c = Nt4(query[i]);
                    q[h.QueryEnd - i - 1] = c >= 4 ? (byte)4 : (byte)(3 - c);
                }
            }
        }

        return (q, t);
    }

    private static byte[] MmIndexSeq(string s, int st, int en)
    {
        en = Math.Min(en, s.Length);
        var r = new byte[Math.Max(0, en - st)];
        for (int i = st; i < en; i++)
            r[i - st] = (byte)Nt4(s[i]);
        return r;
    }

    private static byte[] MmIndexSeqRev(string s, int st, int en)
    {
        en = Math.Min(en, s.Length);
        int st1 = s.Length - en, en1 = s.Length - st;
        var r = new byte[Math.Max(0, en - st)];
        for (int i = st1; i < en1; i++)
        {
            int c = Nt4(s[i]);
            r[en1 - i - 1] = c < 4 ? (byte)(3 - c) : (byte)c;
        }

        return r;
    }

    // write_cs_core
    private static void WriteCs(StringBuilder s, byte[] tseq, byte[] qseq, MinimapAlignment p, bool noIden)
    {
        int qOff = 0, tOff = 0;
        var tmp = new StringBuilder();
        foreach (var (opc, len) in p.Cigar)
        {
            if (opc is 'M' or '=' or 'X')
            {
                tmp.Clear();
                for (int j = 0; j < len; ++j)
                {
                    if (qseq[qOff + j] != tseq[tOff + j])
                    {
                        if (tmp.Length > 0)
                        {
                            if (!noIden)
                                s.Append('=').Append(tmp);
                            else
                                s.Append(':').Append(tmp.Length.ToString(CultureInfo.InvariantCulture));
                            tmp.Clear();
                        }

                        s.Append('*').Append("acgtn"[tseq[tOff + j]]).Append("acgtn"[qseq[qOff + j]]);
                    }
                    else
                    {
                        tmp.Append(NtChars[qseq[qOff + j]]);
                    }
                }

                if (tmp.Length > 0)
                {
                    if (!noIden)
                        s.Append('=').Append(tmp);
                    else
                        s.Append(':').Append(tmp.Length.ToString(CultureInfo.InvariantCulture));
                }

                qOff += len;
                tOff += len;
            }
            else if (opc == 'I')
            {
                s.Append('+');
                for (int j = 0; j < len; ++j)
                    s.Append("acgtn"[qseq[qOff + j]]);
                qOff += len;
            }
            else if (opc == 'D')
            {
                s.Append('-');
                for (int j = 0; j < len; ++j)
                    s.Append("acgtn"[tseq[tOff + j]]);
                tOff += len;
            }
            else
            {
                s.Append('~').Append("acgtn"[tseq[tOff]]).Append("acgtn"[tseq[tOff + 1]]).Append(len.ToString(CultureInfo.InvariantCulture)).Append("acgtn"[tseq[tOff + len - 2]]).Append("acgtn"[tseq[tOff + len - 1]]);
                tOff += len;
            }
        }
    }

    // write_MD_core
    private static void WriteMd(StringBuilder s, byte[] tseq, byte[] qseq, MinimapAlignment p)
    {
        int qOff = 0, tOff = 0, lMd = 0;
        foreach (var (opc, len) in p.Cigar)
        {
            if (opc is 'M' or '=' or 'X')
            {
                for (int j = 0; j < len; ++j)
                {
                    if (qseq[qOff + j] != tseq[tOff + j])
                    {
                        s.Append(lMd.ToString(CultureInfo.InvariantCulture)).Append(NtChars[tseq[tOff + j]]);
                        lMd = 0;
                    }
                    else
                    {
                        ++lMd;
                    }
                }

                qOff += len;
                tOff += len;
            }
            else if (opc == 'I')
            {
                qOff += len;
            }
            else if (opc == 'D')
            {
                s.Append(lMd.ToString(CultureInfo.InvariantCulture)).Append('^');
                for (int j = 0; j < len; ++j)
                    s.Append(NtChars[tseq[tOff + j]]);
                lMd = 0;
                tOff += len;
            }
            else if (opc == 'N')
            {
                tOff += len;
            }
        }

        if (lMd > 0)
            s.Append(lMd.ToString(CultureInfo.InvariantCulture));
    }

    // write_cs_or_MD with the tag prefix
    private static void WriteCsOrMd(StringBuilder s, MinimapHit h, string query, IReadOnlyList<string> targets, MinimapOutputOptions oo, bool qstrand, bool isMd, bool writeTag)
    {
        if (h.Alignment is null)
            return;
        var (q, t) = AlignedBases(h, query, targets, qstrand);
        if (isMd)
        {
            if (writeTag)
                s.Append("\tMD:Z:");
            WriteMd(s, t, q, h.Alignment);
        }
        else
        {
            if (writeTag)
                s.Append("\tcs:Z:");
            WriteCs(s, t, q, h.Alignment, oo.Cs != MinimapCsMode.Long);
        }
    }

    /// <summary><c>mm_gen_cs</c> / <c>mm_gen_MD</c>: the <c>cs</c> string (short or long) or the <c>MD</c> string of an aligned hit (without the tag).</summary>
    /// <param name="hit">An aligned hit (<see cref="MinimapHit.Alignment"/> not null).</param>
    /// <param name="query">The query the hit belongs to.</param>
    /// <param name="targets">The target sequences.</param>
    /// <param name="md"><c>true</c> for MD, <c>false</c> for cs.</param>
    /// <param name="longCs"><c>true</c> for <c>--cs=long</c>.</param>
    /// <param name="queryStrand"><c>--qstrand</c>.</param>
    public static string GenerateCsOrMd(MinimapHit hit, string query, IReadOnlyList<string> targets, bool md = false, bool longCs = false, bool queryStrand = false)
    {
        ArgumentNullException.ThrowIfNull(hit);
        var sb = new StringBuilder();
        WriteCsOrMd(sb, hit, query, targets, new MinimapOutputOptions { Cs = longCs ? MinimapCsMode.Long : MinimapCsMode.Short }, queryStrand, md, false);
        return sb.ToString();
    }

    private static string CigarText(MinimapAlignment p)
    {
        var sb = new StringBuilder();
        foreach (var (op, len) in p.Cigar)
            sb.Append(len.ToString(CultureInfo.InvariantCulture)).Append(op);
        return sb.ToString();
    }

    // mm_event_identity of a hit
    private static double HitIdentity(MinimapHit h)
    {
        var p = h.Alignment!;
        int nGap = 0, nGapo = 0;
        foreach (var (op, len) in p.Cigar)
            if (op is 'I' or 'D')
            {
                ++nGapo;
                nGap += len;
            }

        return (double)h.MatchLength / (h.BlockLength + p.AmbiguousBases - nGap + nGapo);
    }

    // write_tags
    private static void WriteTags(StringBuilder s, MinimapHit r)
    {
        char type = r.Id == r.Parent ? (r.IsInversion ? 'I' : 'P') : (r.IsInversion ? 'i' : 'S');
        if (r.Alignment is { } p)
        {
            s.Append("\tNM:i:").Append(r.BlockLength - r.MatchLength + p.AmbiguousBases).Append("\tms:i:").Append(p.DpMax).Append("\tAS:i:").Append(p.DpScore).Append("\tnn:i:").Append(p.AmbiguousBases);
            if (p.TranscriptStrand is 1 or 2)
                s.Append("\tts:A:").Append("?+-?"[p.TranscriptStrand]);
        }

        s.Append("\ttp:A:").Append(type).Append("\tcm:i:").Append(r.AnchorCount).Append("\ts1:i:").Append(r.Score);
        if (r.Parent == r.Id)
            s.Append("\ts2:i:").Append(r.SecondaryScore);
        if (r.Alignment is not null)
        {
            double div = 1.0 - HitIdentity(r);
            s.Append("\tde:f:").Append(C4(div));
        }
        else if (r.Divergence >= 0f && r.Divergence <= 1f)
        {
            s.Append("\tdv:f:").Append(r.Divergence == 0f ? "0" : CFormat.Fixed(r.Divergence, 4));
        }

        if (r.Split != 0)
            s.Append("\tzd:i:").Append(r.Split);
    }

    private static string FormatPafLine(MinimapQuery t, MinimapHit? r, int repLen, IReadOnlyList<string?> targetNames, IReadOnlyList<int> targetLengths, IReadOnlyList<string> targets, MinimapOutputOptions oo, bool qstrand)
    {
        var s = new StringBuilder();
        if (r is null)
        {
            s.Append(t.Name).Append('\t').Append(t.Sequence.Length).Append("\t0\t0\t*\t*\t0\t0\t0\t0\t0\t0");
            if (repLen >= 0)
                s.Append("\trl:i:").Append(repLen);
            return s.Append('\n').ToString();
        }

        int len = targetLengths[r.TargetId];
        s.Append(t.Name).Append('\t').Append(t.Sequence.Length).Append('\t').Append(r.QueryStart).Append('\t').Append(r.QueryEnd).Append('\t').Append(r.Reverse ? '-' : '+').Append('\t');
        s.Append(r.TargetId < targetNames.Count && targetNames[r.TargetId] is { } name ? name : r.TargetId.ToString(CultureInfo.InvariantCulture));
        s.Append('\t').Append(len);
        if (qstrand && r.Reverse)
            s.Append('\t').Append(len - r.TargetEnd).Append('\t').Append(len - r.TargetStart);
        else
            s.Append('\t').Append(r.TargetStart).Append('\t').Append(r.TargetEnd);
        s.Append('\t').Append(r.MatchLength).Append('\t').Append(r.BlockLength).Append('\t').Append(r.MapQ);
        WriteTags(s, r);
        if (repLen >= 0)
            s.Append("\trl:i:").Append(repLen);
        if (r.Alignment is not null && oo.CigarTag)
            s.Append("\tcg:Z:").Append(CigarText(r.Alignment));
        if (r.Alignment is not null && (oo.Cs != MinimapCsMode.None || oo.Md))
            WriteCsOrMd(s, r, t.Sequence, targets, oo, qstrand, oo.Md, true);
        if (oo.CopyComment && t.Comment is not null)
            s.Append('\t').Append(t.Comment);
        return s.Append('\n').ToString();
    }

    private static MinimapHit? GetSamPri(IReadOnlyList<MinimapHit> regs)
    {
        foreach (var r in regs)
            if (r.IsSamPrimary)
                return r;
        return null;
    }

    private static void WriteSamCigar(StringBuilder s, int samFlag, bool inTag, int qlen, MinimapHit r, MinimapOutputOptions oo)
    {
        if (r.Alignment is null)
        {
            s.Append('*');
            return;
        }

        int clip0 = r.Reverse ? qlen - r.QueryEnd : r.QueryStart;
        int clip1 = r.Reverse ? r.QueryStart : qlen - r.QueryEnd;
        bool hard = ((samFlag & 0x800) != 0 || ((samFlag & 0x100) != 0 && oo.SecondarySeq)) && !oo.SoftClip;
        var ops = "MIDNSHP=XB";
        if (inTag)
        {
            int clipChar = hard ? 5 : 4;
            s.Append("\tCG:B:I");
            if (clip0 != 0)
                s.Append(',').Append(((uint)clip0 << 4 | (uint)clipChar).ToString(CultureInfo.InvariantCulture));
            foreach (var (op, len) in r.Alignment.Cigar)
                s.Append(',').Append(((uint)len << 4 | (uint)ops.IndexOf(op)).ToString(CultureInfo.InvariantCulture));
            if (clip1 != 0)
                s.Append(',').Append(((uint)clip1 << 4 | (uint)clipChar).ToString(CultureInfo.InvariantCulture));
        }
        else
        {
            char clipChar = hard ? 'H' : 'S';
            if (clip0 != 0)
                s.Append(clip0).Append(clipChar);
            s.Append(CigarText(r.Alignment));
            if (clip1 != 0)
                s.Append(clip1).Append(clipChar);
        }
    }

    // mm_write_sam3: the SAM line of hit regIdx (-1: unmapped) of segment segIdx
    private static string FormatSamLine(IReadOnlyList<MinimapQuery> queries, int segIdx, int regIdx, IReadOnlyList<IReadOnlyList<MinimapHit>> regss, int repLen, IReadOnlyList<string?> targetNames, IReadOnlyList<string> targets, MinimapOutputOptions oo, string rgId)
    {
        const int maxBamCigarOp = 65535;
        int nSeg = queries.Count;
        var t = queries[segIdx];
        var regs = regss[segIdx];
        int nRegs = regs.Count;
        bool cigarInTag = false;
        int thisRid = -1, thisPos = -1;
        MinimapHit? r = nRegs > 0 && regIdx < nRegs && regIdx >= 0 ? regs[regIdx] : null;
        MinimapHit? rPrev = null, rNext = null;
        string TName(int rid) => targetNames[rid] ?? rid.ToString(CultureInfo.InvariantCulture);
        if (nSeg > 1)
        {
            int nextSid = (segIdx + 1) % nSeg;
            rNext = GetSamPri(regss[nextSid]);
            if (nSeg > 2)
            {
                for (int i = 1; i <= nSeg - 1; ++i)
                {
                    int prevSid = (segIdx + nSeg - i) % nSeg;
                    if (regss[prevSid].Count > 0)
                    {
                        rPrev = GetSamPri(regss[prevSid]);
                        break;
                    }
                }
            }
            else
            {
                rPrev = rNext;
            }
        }

        var s = new StringBuilder();
        string qname = nSeg > 1 ? t.Name[..QnameLength(t.Name)] : t.Name;
        s.Append(qname);
        int flag = nSeg > 1 ? 0x1 : 0x0;
        if (r is null)
        {
            flag |= 0x4;
        }
        else
        {
            if (r.Reverse)
                flag |= 0x10;
            if (r.Parent != r.Id)
                flag |= 0x100;
            else if (!r.IsSamPrimary)
                flag |= 0x800;
        }

        if (nSeg > 1)
        {
            if (r is not null && r.IsProperFragment)
                flag |= 0x2;
            if (segIdx == 0)
                flag |= 0x40;
            else if (segIdx == nSeg - 1)
                flag |= 0x80;
            if (rNext is null)
                flag |= 0x8;
            else if (rNext.Reverse)
                flag |= 0x20;
        }

        s.Append('\t').Append(flag);
        int qlen = t.Sequence.Length;
        if (r is null)
        {
            if (rPrev is not null)
            {
                thisRid = rPrev.TargetId;
                thisPos = rPrev.TargetStart;
                s.Append('\t').Append(TName(thisRid)).Append('\t').Append(thisPos + 1).Append("\t0\t*");
            }
            else
            {
                s.Append("\t*\t0\t0\t*");
            }
        }
        else
        {
            thisRid = r.TargetId;
            thisPos = r.TargetStart;
            s.Append('\t').Append(TName(r.TargetId)).Append('\t').Append(r.TargetStart + 1).Append('\t').Append(r.MapQ).Append('\t');
            if (oo.LongCigar && r.Alignment is not null && r.Alignment.Cigar.Count > maxBamCigarOp - 2)
            {
                int nCigar = r.Alignment.Cigar.Count;
                if (r.QueryStart != 0)
                    ++nCigar;
                if (r.QueryEnd != qlen)
                    ++nCigar;
                if (nCigar > maxBamCigarOp)
                    cigarInTag = true;
            }

            if (cigarInTag)
            {
                int slen;
                if ((flag & 0x900) == 0 || oo.SoftClip)
                    slen = qlen;
                else if ((flag & 0x100) != 0 && !oo.SecondarySeq)
                    slen = 0;
                else
                    slen = r.QueryEnd - r.QueryStart;
                s.Append(slen).Append('S').Append(r.TargetEnd - r.TargetStart).Append('N');
            }
            else
            {
                WriteSamCigar(s, flag, false, qlen, r, oo);
            }
        }

        if (nSeg > 1)
        {
            int tlen = 0;
            if (thisRid >= 0 && rNext is not null)
            {
                if (thisRid == rNext.TargetId)
                {
                    if (r is not null)
                    {
                        int thisPos5 = r.Reverse ? r.TargetEnd - 1 : thisPos;
                        int nextPos5 = rNext.Reverse ? rNext.TargetEnd - 1 : rNext.TargetStart;
                        tlen = nextPos5 - thisPos5;
                    }

                    s.Append("\t=\t");
                }
                else
                {
                    s.Append('\t').Append(TName(rNext.TargetId)).Append('\t');
                }

                s.Append(rNext.TargetStart + 1).Append('\t');
            }
            else if (rNext is not null)
            {
                s.Append('\t').Append(TName(rNext.TargetId)).Append('\t').Append(rNext.TargetStart + 1).Append('\t');
            }
            else if (thisRid >= 0)
            {
                s.Append("\t=\t").Append(thisPos + 1).Append('\t');
            }
            else
            {
                s.Append("\t*\t0\t");
            }

            if (tlen > 0)
                ++tlen;
            else if (tlen < 0)
                --tlen;
            s.Append(tlen).Append('\t');
        }
        else
        {
            s.Append("\t*\t0\t0\t");
        }

        if (r is null)
        {
            s.Append(t.Sequence).Append('\t').Append(t.Quality ?? "*");
        }
        else if ((flag & 0x900) == 0 || oo.SoftClip)
        {
            s.Append(SamSq(t.Sequence, r.Reverse, r.Reverse)).Append('\t').Append(t.Quality is null ? "*" : SamSq(t.Quality, r.Reverse, false));
        }
        else if ((flag & 0x100) != 0 && !oo.SecondarySeq)
        {
            s.Append("*\t*");
        }
        else
        {
            int ql = r.QueryEnd - r.QueryStart;
            s.Append(SamSq(t.Sequence.Substring(r.QueryStart, ql), r.Reverse, r.Reverse)).Append('\t').Append(t.Quality is null ? "*" : SamSq(t.Quality.Substring(r.QueryStart, ql), r.Reverse, false));
        }

        if (rgId.Length > 0)
            s.Append("\tRG:Z:").Append(rgId);
        if (nSeg > 2)
            s.Append("\tFI:i:").Append(segIdx);
        if (r is not null)
        {
            WriteTags(s, r);
            if (r.Parent == r.Id && r.Alignment is not null && nRegs > 1)
            {
                int nSa = 0;
                for (int i = 0; i < nRegs; ++i)
                    if (i != regIdx && regs[i].Parent == regs[i].Id && regs[i].Alignment is not null)
                        ++nSa;
                if (nSa > 0)
                {
                    s.Append("\tSA:Z:");
                    for (int i = 0; i < nRegs; ++i)
                    {
                        var q = regs[i];
                        int lM, lI = 0, lD = 0, clip5, clip3;
                        if (i == regIdx || q.Parent != q.Id || q.Alignment is null)
                            continue;
                        if (q.QueryEnd - q.QueryStart < q.TargetEnd - q.TargetStart)
                        {
                            lM = q.QueryEnd - q.QueryStart;
                            lD = (q.TargetEnd - q.TargetStart) - lM;
                        }
                        else
                        {
                            lM = q.TargetEnd - q.TargetStart;
                            lI = (q.QueryEnd - q.QueryStart) - lM;
                        }

                        clip5 = q.Reverse ? qlen - q.QueryEnd : q.QueryStart;
                        clip3 = q.Reverse ? q.QueryStart : qlen - q.QueryEnd;
                        s.Append(TName(q.TargetId)).Append(',').Append(q.TargetStart + 1).Append(',').Append(q.Reverse ? '-' : '+').Append(',');
                        if (clip5 != 0)
                            s.Append(clip5).Append('S');
                        if (lM != 0)
                            s.Append(lM).Append('M');
                        if (lI != 0)
                            s.Append(lI).Append('I');
                        if (lD != 0)
                            s.Append(lD).Append('D');
                        if (clip3 != 0)
                            s.Append(clip3).Append('S');
                        s.Append(',').Append(q.MapQ).Append(',').Append(q.BlockLength - q.MatchLength + q.Alignment.AmbiguousBases).Append(';');
                    }
                }
            }

            if (r.Alignment is not null && (oo.Cs != MinimapCsMode.None || oo.Md))
                WriteCsOrMd(s, r, t.Sequence, targets, oo, false, oo.Md, true);
            if (cigarInTag)
                WriteSamCigar(s, flag, true, qlen, r, oo);
        }

        if (repLen >= 0)
            s.Append("\trl:i:").Append(repLen);
        if (oo.CopyComment && t.Comment is not null)
            s.Append('\t').Append(t.Comment);
        return s.Append('\n').ToString();
    }

    /// <summary>
    /// The output lines of one fragment (the worker step of <c>mm_map_file_frag</c>): for every read of the fragment its hits as PAF lines (<c>mm_write_paf3</c>: the 12 columns and the tags <c>NM ms AS nn ts tp cm s1 s2 de|dv zd rl cg cs MD</c>) or SAM lines (<c>mm_write_sam3</c>: flags, mate fields,
    /// <c>SA</c>, <c>RG</c>, <c>FI</c>, the tags), secondary hits skipped with <c>--secondary=no</c>, a line for an unmapped read with <c>--paf-no-hit</c> / (SAM) unless <c>--sam-hit-only</c>.
    /// </summary>
    /// <param name="queries">The reads of the fragment (1, or 2 for paired reads).</param>
    /// <param name="result">The result of <see cref="MapFragment"/> (or of <see cref="Map"/> as a one-element fragment).</param>
    /// <param name="targetNames">The target names (null entries print the index).</param>
    /// <param name="targetLengths">The target lengths.</param>
    /// <param name="targets">The target sequences (for <c>cs</c> / <c>MD</c>).</param>
    /// <param name="output">The output options.</param>
    /// <param name="queryStrand"><c>--qstrand</c>.</param>
    public static string FormatFragment(IReadOnlyList<MinimapQuery> queries, MinimapFragmentResult result, IReadOnlyList<string?> targetNames, IReadOnlyList<int> targetLengths, IReadOnlyList<string> targets, MinimapOutputOptions? output = null, bool queryStrand = false)
    {
        ArgumentNullException.ThrowIfNull(queries);
        ArgumentNullException.ThrowIfNull(result);
        var oo = output ?? new MinimapOutputOptions();
        var rg = ParseReadGroup(oo.ReadGroup) ?? ("", "");          // an invalid -R line leaves the read group empty (the program stops after the header in SAM mode)
        var regss = result.Segments.Select(sg => sg.Hits).ToArray();
        var sb = new StringBuilder();
        for (int i = 0; i < queries.Count; i++)
        {
            var seg = result.Segments[i];
            if (seg.Hits.Count > 0)
            {
                for (int j = 0; j < seg.Hits.Count; j++)
                {
                    var r = seg.Hits[j];
                    if (oo.NoSecondary && r.Id != r.Parent)
                        continue;
                    sb.Append(oo.Sam ? FormatSamLine(queries, i, j, regss, seg.RepeatLength, targetNames, targets, oo, rg.Id) : FormatPafLine(queries[i], r, seg.RepeatLength, targetNames, targetLengths, targets, oo, queryStrand));
                }
            }
            else if (oo.PafNoHit || (oo.Sam && !oo.SamHitOnly))
            {
                sb.Append(oo.Sam ? FormatSamLine(queries, i, -1, regss, seg.RepeatLength, targetNames, targets, oo, rg.Id) : FormatPafLine(queries[i], null, seg.RepeatLength, targetNames, targetLengths, targets, oo, queryStrand));
            }
        }

        return sb.ToString();
    }
}
