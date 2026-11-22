#pragma warning disable S1121, S3358, S2325, S1066, S3776, S107, S1854, S1481, S2589, S3626 // the writers keep the control flow of ajseqwrite.c

using System.Globalization;
using System.Text;

namespace Seqeron.Genomics.Infrastructure;

// The sequence writers of the formats that write a set of sequences when the output is closed (seqWriteMsf, seqWriteClustal, seqWriteSelex, seqWritePhylip, seqWritePhylipnon, seqWriteHennig86, seqWriteMega, seqWriteMeganon,
// seqWriteNexus, seqWriteNexusnon, seqWriteJackknifer, seqWriteJackknifernon, seqWriteTreecon) and the list they read (seqWriteListAppend). A5-8c3h, B09 F266.
internal static partial class EmbossSeqWrite
{
    /// <summary>The writers of this port that keep the sequences and write them when the output is closed.</summary>
    internal static bool IsSetWriter(EmbossSeqOutFormat f) => f.Save;

    /// <summary>
    /// With <c>-ossingle</c> a set writer writes the list at once, to the output object whose file is not open: the formatted print of the library ignores a null file, the functions that write a line or a newline do not, and the
    /// conversions <c>%F</c> (msf) dereference it. The writers that call them (selex, msf, hennig86, mega, nexus, nexusnon always; clustal for a length above 60, phylip and jackknifernon for a length above 50) kill the process by a
    /// segmentation fault; the others (phylipnon, meganon, jackknifer, treecon) write nothing.
    /// </summary>
    internal static bool WritesToNullFileDies(EmbossSeqOutFormat f, int maxLength) => f.Writer switch
    {
        "Msf" or "Selex" or "Hennig86" or "Mega" or "Nexus" or "Nexusnon" => true,
        "Clustal" => maxLength > 60,
        "Phylip" or "Jackknifernon" => maxLength > 50,
        _ => false,
    };

    /// <summary>A program that writes new sequences (extractfeat, coderet, union) gives a set format the sequence to keep: true when it was kept (written when the output is closed with <see cref="CloseSet"/>).</summary>
    internal static bool Save(EmbossSeqOutFormat f, EmbossSeqWriteRun run, string name, string desc, string type, string seq)
    {
        if (!f.Save)
            return false;
        if (type.Length > 0)
            run.OutType = type;
        ListAppend(run, new EmbossSeqSaved { Name = name, Desc = desc, Type = type, Seq = seq }, "");
        return true;
    }

    /// <summary><c>ajSeqoutClose</c> of a set format: the saved sequences are written.</summary>
    internal static void CloseSet(StringBuilder o, EmbossSeqOutFormat? f, EmbossSeqWriteRun run, int count)
    {
        if (f is { Save: true } && run.Saved.Count > 0)
        {
            run.Count = count;
            WriteSaved(o, f, run);
        }

        if (f is { Writer: "Bam" } && FinishBam(run) is { } bam)
            o.Append(bam);
    }

    /// <summary>The feature format that a sequence format writes into its file (<c>ajFeattabOutNewCSF</c> of the writer).</summary>
    internal static string LocalFeatureFormat(EmbossSeqOutFormat f, bool protein = false) => f.Writer switch
    {
        "Experiment" when protein => "swiss",                          // seqWriteExperiment writes a protein with seqWriteSwiss
        "Embl" or "Emblnew" or "Experiment" => "embl",
        "Genbank" => "genbank",
        "Genpept" => "genpept",
        "Refseq" => "refseq",
        "Refseqp" => "refseqp",
        "Swiss" or "Swissnew" => "swiss",
        "Nbrf" => "pir",
        "Gff2" => "gff2",
        "Gff3" => "gff3",
        _ => "gff3",
    };

    /// <summary>The name of the writer function in the warning <c>seqWriteEmbl features output failed UFO</c>.</summary>
    internal static string LocalFeatureWriterName(EmbossSeqOutFormat f, bool protein = false) => f.Writer switch
    {
        "Experiment" when protein => "seqWriteSwiss",
        "Embl" or "Emblnew" or "Experiment" => "seqWriteEmbl",
        "Swiss" or "Swissnew" => "seqWriteSwiss",
        _ => "seqWrite" + f.Writer,
    };

    private static string Today(EmbossSeqWriteRun run) => run.RunDate.ToString("dd/MM/yy", CultureInfo.InvariantCulture);

    // seqGapSL: the characters of " .~-" become the gap character, the ends of the sequence the pad character (when there is one), and the sequence is padded to a length
    private static string GapSL(string seq, char gapc, char padc, int ilen)
    {
        const string test = " .~-";
        char[] cp = seq.ToCharArray();
        for (int i = 0; i < cp.Length; i++)
        {
            if (test.Contains(cp[i], StringComparison.Ordinal))
                cp[i] = gapc;
        }

        char endc = gapc;
        if (padc != 0)
        {
            endc = padc;
            if (cp.Length > 0 && test.Contains(cp[0], StringComparison.Ordinal))
            {
                int k = 0;
                while (k < cp.Length && test.Contains(cp[k], StringComparison.Ordinal))
                    cp[k++] = padc;
                if (k >= cp.Length)
                    throw new InvalidOperationException("The sequence has only gaps: the library pads past the end of the string (the reference program overwrites memory it does not own).");
            }

            if (cp.Length > 0 && test.Contains(cp[^1], StringComparison.Ordinal))
            {
                for (int i = cp.Length - 1; i > 0 && test.Contains(cp[i], StringComparison.Ordinal); i--)
                    cp[i] = padc;
            }
        }

        var r = new StringBuilder().Append(cp);
        while (r.Length < ilen)
            r.Append(endc);
        return r.ToString();
    }

    private static string GapStr(string seq, char gapc) => GapSL(seq, gapc, '\0', 0);

    // ajStrAssignSubS(pos1, pos2) with the positions of ajCvtSposToPos / ajCvtSposToPosStart: a position past the end is the last character
    private static string SubS(string s, int pos1, int pos2)
    {
        int len = s.Length;
        if (len == 0)
            return "";
        int ibegin = pos1 < 0 ? len + pos1 : pos1;
        if (ibegin >= len)
            ibegin = len - 1;
        if (ibegin < 0)
            ibegin = 0;
        int iend = pos2 < 0 ? len + pos2 : pos2;
        if (iend >= len)
            iend = len - 1;
        if (iend < ibegin)
            iend = ibegin;
        return s.Substring(ibegin, iend - ibegin + 1);
    }

    // ajStrMaskIdent: the characters of the string that equal those of the first one are masked
    private static string MaskIdent(string s, string first, char mask)
    {
        char[] cp = s.ToCharArray();
        for (int i = 0; i < first.Length; i++)
        {
            if (i >= cp.Length)
                break;
            if (cp[i] == first[i])
                cp[i] = mask;
        }

        return new string(cp);
    }

    // ajStrFmtBlock: a blank after every len characters
    private static string FmtBlock(string s, int len)
    {
        var b = new StringBuilder();
        for (int i = 0; i < s.Length; i++)
        {
            if (i > 0 && i % len == 0)
                b.Append(' ');
            b.Append(s[i]);
        }

        return b.ToString();
    }

    // %-W.PS
    private static string Left(string s, int width, int prec)
    {
        if (s.Length > prec)
            s = s[..prec];
        return s.Length >= width ? s : s + new string(' ', width - s.Length);
    }

    private static string Right(string s, int width) => s.Length >= width ? s : new string(' ', width - s.Length) + s;

    /// <summary>
    /// <c>seqWriteListAppend</c>: the sequence (a copy cut to its range) is added to the list that the set formats write when the output is closed; a sequence with no name takes the name of the entry or <c>EMBOSS_001</c> ... (a counter of the library
    /// that lasts for the program).
    /// </summary>
    internal static void ListAppend(EmbossSeqWriteRun run, EmbossSeqSaved item, string entryname)
    {
        if (item.Name.Length == 0)
        {
            if (entryname.Length > 0)
            {
                if (run.NameMultiCount != 0)
                    item.Name = entryname + "_" + (++run.NameMultiCount).ToString("000", CultureInfo.InvariantCulture);
                else
                {
                    item.Name = entryname;
                    ++run.NameMultiCount;
                }
            }
            else
            {
                item.Name = "EMBOSS_" + (++run.NameMultiCount).ToString("000", CultureInfo.InvariantCulture);
            }
        }

        run.Saved.Add(item);
    }

    /// <summary>Writes the saved sequences in the format of a set writer (the call of <c>ajSeqoutClose</c> and, with <c>-ossingle</c>, of every <c>seqWriteListAppend</c>).</summary>
    internal static void WriteSaved(StringBuilder o, EmbossSeqOutFormat format, EmbossSeqWriteRun run)
    {
        var seqs = run.Saved;
        if (seqs.Count == 0)
            return;
        switch (format.Writer)
        {
            case "Msf": WriteMsf(o, run); break;
            case "Clustal": WriteClustal(o, run); break;
            case "Selex": WriteSelex(o, run); break;
            case "Phylip": WritePhylip(o, run, true); break;
            case "Phylipnon": WritePhylip(o, run, false); break;
            case "Hennig86": WriteHennig86(o, run); break;
            case "Mega": WriteMega(o, run, true); break;
            case "Meganon": WriteMega(o, run, false); break;
            case "Nexus": WriteNexus(o, run, true); break;
            case "Nexusnon": WriteNexus(o, run, false); break;
            case "Jackknifer": WriteJackknifer(o, run, true); break;
            case "Jackknifernon": WriteJackknifer(o, run, false); break;
            case "Treecon": WriteTreecon(o, run); break;
            default: throw new NotSupportedException($"The sequence output format '{format.Name}' is not written by this port.");
        }
    }

    private static int MaxLen(List<EmbossSeqSaved> seqs)
    {
        int ilen = 0;
        foreach (var s in seqs)
        {
            if (ilen < s.Seq.Length)
                ilen = s.Seq.Length;
        }

        return ilen;
    }

    private static void WriteHennig86(StringBuilder o, EmbossSeqWriteRun run)
    {
        var seqs = run.Saved;
        int ilen = MaxLen(seqs);
        o.Append("xread\n");
        o.Append("' Written by EMBOSS ").Append(Today(run)).Append(" '\n");
        o.Append(ilen.ToString(CultureInfo.InvariantCulture)).Append(' ').Append(seqs.Count.ToString(CultureInfo.InvariantCulture)).Append('\n');
        foreach (var seq in seqs)
        {
            var b = new StringBuilder();
            foreach (char c in seq.Seq)
            {
                b.Append(c switch { 'A' or 'a' => '0', 'T' or 't' or 'U' or 'u' => '1', 'G' or 'g' => '2', 'C' or 'c' => '3', _ => '?' });
            }

            o.Append(seq.Name).Append('\n');
            if (b.Length > 0)
                o.Append(b).Append('\n');                              // ajWritelineNewline of an empty string writes nothing (fwrite returns 0: no newline either)
        }

        o.Append(";\n");
    }

    private static void WriteTreecon(StringBuilder o, EmbossSeqWriteRun run)
    {
        var seqs = run.Saved;
        o.Append(MaxLen(seqs).ToString(CultureInfo.InvariantCulture)).Append('\n');
        foreach (var seq in seqs)
            o.Append(seq.Name).Append('\n').Append(GapStr(seq.Seq, '-')).Append('\n');
    }

    private static void WriteJackknifer(StringBuilder o, EmbossSeqWriteRun run, bool interleaved)
    {
        var seqs = run.Saved;
        int ilen = MaxLen(seqs);
        const int wid = 50;
        o.Append("' Written by EMBOSS ").Append(Today(run)).Append(" \n");
        if (interleaved)
        {
            for (int ipos = 1; ipos <= ilen; ipos += wid)
            {
                int iend = Math.Min(ipos + wid - 1, ilen);
                foreach (var seq in seqs)
                    o.Append(Left("(" + seq.Name + ")", 20, 20)).Append(' ').Append(GapStr(SubS(seq.Seq, ipos - 1, iend - 1), '-')).Append('\n');
            }
        }
        else
        {
            foreach (var seq in seqs)
            {
                for (int ipos = 1; ipos <= ilen; ipos += wid)
                {
                    int iend = Math.Min(ipos + wid - 1, ilen);
                    string sseq = GapStr(SubS(seq.Seq, ipos - 1, iend - 1), '-');
                    if (ipos == 1)
                        o.Append(Left("(" + seq.Name + ")", 20, 20)).Append(' ').Append(sseq).Append('\n');
                    else
                        o.Append(sseq).Append('\n');
                }
            }
        }

        o.Append(";\n");
    }

    private static void WriteNexus(StringBuilder o, EmbossSeqWriteRun run, bool interleaved)
    {
        var seqs = run.Saved;
        int ilen = MaxLen(seqs);
        if (interleaved)
        {
            foreach (var seq in seqs)
                seq.Seq = GapSL(seq.Seq, '-', '-', ilen);
        }

        bool prot = run.OutType.Length > 0 && run.OutType[0] == 'P';
        o.Append("#NEXUS\n");
        o.Append("[TITLE: Written by EMBOSS ").Append(Today(run)).Append("]\n\n");
        o.Append("begin data;\n");
        o.Append(CultureInfo.InvariantCulture, $"dimensions ntax={seqs.Count} nchar={ilen};\n");
        string inter = interleaved ? "interleave " : "";
        o.Append(prot ? $"format {inter}datatype=protein missing=X gap=-;\n" : $"format {inter}datatype=DNA missing=N gap=-;\n");
        o.Append('\n');
        o.Append("matrix\n");
        if (interleaved)
        {
            const int wid = 50;
            for (int ipos = 1; ipos <= ilen; ipos += wid)
            {
                int iend = Math.Min(ipos + wid - 1, ilen);
                if (ipos > 1)
                    o.Append('\n');
                foreach (var seq in seqs)
                    o.Append(Left(seq.Name, 20, 20)).Append(' ').Append(GapStr(SubS(seq.Seq, ipos - 1, iend - 1), '-')).Append('\n');
            }
        }
        else
        {
            foreach (var seq in seqs)
                o.Append(seq.Name).Append('\n').Append(GapStr(seq.Seq, '-')).Append('\n');
        }

        o.Append(";\n\n");
        o.Append("end;\n");
        o.Append("begin assumptions;\n");
        o.Append("options deftype=unord;\n");
        o.Append("end;\n");
    }

    private static void WriteMega(StringBuilder o, EmbossSeqWriteRun run, bool interleaved)
    {
        var seqs = run.Saved;
        int isize = seqs.Count;
        var seqfirst = seqs[0];
        int namewidth = seqfirst.Name.Length;
        int ilen = 0;
        bool onedesc, onegene;
        string genestr = "";
        if (isize == 1)
        {
            onegene = false;
            onedesc = false;
        }
        else
        {
            onegene = true;
            onedesc = true;
            ilen = seqfirst.Seq.Length;
            if (seqfirst.Desc.Length == 0)
                onedesc = false;
            if (seqfirst.Gene is not null)
                genestr = seqfirst.Gene;
            else
                onegene = false;
            for (int i = 1; i < isize; i++)
            {
                var seq = seqs[i];
                if (seq.Name.Length > namewidth)
                    namewidth = seq.Name.Length;
                if (ilen < seq.Seq.Length)
                    ilen = seq.Seq.Length;
                if (onedesc && seq.Desc != seqfirst.Desc)
                    onedesc = false;
                if (onegene && seq.Gene is not null && genestr != seq.Gene)
                    onegene = false;
            }
        }

        if (namewidth > 40)
            namewidth = 40;
        if (!interleaved)
        {
            foreach (var seq in seqs)
            {
                if (ilen < seq.Seq.Length)
                    ilen = seq.Seq.Length;
                if (seq.Desc != seqfirst.Desc)
                    onedesc = false;
            }
        }

        bool prot = seqfirst.Type.Length > 0 && seqfirst.Type[0] == 'P';
        if (run.Count == isize)
        {
            o.Append("#mega\n");
            if (onedesc)
                o.Append("!Title: ").Append(seqfirst.Desc).Append(";\n");
            else
                o.Append("!Title: Written by EMBOSS ").Append(Today(run)).Append(";\n");
            o.Append("!Format\n");
            if (interleaved)
                o.Append(prot ? "    DataType=Protein DataFormat=Interleaved\n" : "    DataType=Nucleotide DataFormat=Interleaved\n");
            else
                o.Append(prot ? "    DataType=Protein\n" : "    DataType=Nucleotide\n");
            o.Append("    Identical=. Indel=- Missing=?\n");
            if (!prot)
                o.Append("    CodeTable=Standard\n");
            o.Append("    ;\n\n");
        }

        if (interleaved)
        {
            o.Append('\n').Append('\n');
            if (onegene)
                o.Append("!Gene=").Append(genestr).Append(";\n");
            const int wid = 50;
            for (int ipos = 1; ipos <= ilen; ipos += wid)
            {
                int iend = Math.Min(ipos + wid - 1, ilen);
                string sseqfirst = GapStr(SubS(seqfirst.Seq, ipos - 1, iend - 1), '-');
                o.Append('\n');
                for (int i = 0; i < isize; i++)
                {
                    var seq = seqs[i];
                    string sseq = GapStr(SubS(seq.Seq, ipos - 1, iend - 1), '-');
                    if (i > 0)
                        sseq = MaskIdent(sseq, sseqfirst, '.');
                    if (!onedesc && ipos == 1 && seq.Desc.Length > 0)
                        o.Append('[').Append(seq.Desc).Append("]\n");
                    o.Append('#').Append(Left(seq.Name, namewidth, namewidth)).Append(' ').Append(sseq).Append('\n');
                }
            }

            o.Append('\n');
        }
        else
        {
            string sseqfirst = GapStr(seqfirst.Seq, '-');
            if (onegene)
                o.Append("!Gene=").Append(genestr).Append(";\n");
            for (int i = 0; i < isize; i++)
            {
                var seq = seqs[i];
                string sseq = GapStr(seq.Seq, '-');
                if (i > 0)
                    sseq = MaskIdent(sseq, sseqfirst, '.');
                if (!onedesc && seq.Desc.Length > 0)
                    o.Append('[').Append(seq.Desc).Append("]\n");
                o.Append('#').Append(Left(seq.Name, namewidth, namewidth)).Append('\n').Append(sseq).Append('\n');
            }
        }
    }

    private static void WritePhylip(StringBuilder o, EmbossSeqWriteRun run, bool interleaved)
    {
        var seqs = run.Saved;
        int ilen = MaxLen(seqs);
        string Padded(EmbossSeqSaved s) => s.Seq.Length < ilen ? s.Seq + new string('-', ilen - s.Seq.Length) : s.Seq;
        if (interleaved)
        {
            o.Append(CultureInfo.InvariantCulture, $" {seqs.Count} {ilen}\n");
            int iend = 0;
            for (int ipos = 1; ipos <= ilen; ipos += 50)
            {
                iend = Math.Min(ipos + 49, ilen);
                foreach (var seq in seqs)
                {
                    string sseq = FmtBlock(GapStr(SubS(Padded(seq), ipos - 1, iend - 1), '-'), 10);
                    if (ipos == 1)
                        o.Append(Left(seq.Name, 10, 10)).Append(sseq).Append('\n');
                    else
                        o.Append(new string(' ', 10)).Append(sseq).Append('\n');
                }

                if (iend < ilen)
                    o.Append('\n');
            }
        }
        else
        {
            o.Append(CultureInfo.InvariantCulture, $"{seqs.Count} {ilen}\n");
            int iend = 0;
            foreach (var seq in seqs)
            {
                string tstr = Padded(seq);
                for (int ipos = 1; ipos <= ilen; ipos += 50)
                {
                    iend = Math.Min(ipos + 49, ilen);
                    string sseq = FmtBlock(GapStr(SubS(tstr, ipos - 1, iend - 1), '-'), 10);
                    if (ipos == 1)
                        o.Append(Left(seq.Name, 10, 10)).Append(sseq).Append('\n');
                    else
                        o.Append(new string(' ', 10)).Append(sseq).Append('\n');
                }

                if (iend < ilen)
                    o.Append('\n');
            }
        }
    }

    private static void WriteClustal(StringBuilder o, EmbossSeqWriteRun run)
    {
        var seqs = run.Saved;
        int ilen = MaxLen(seqs);
        foreach (var seq in seqs)
        {
            if (ilen > seq.Seq.Length)
                seq.Seq += new string('-', ilen - seq.Seq.Length);
        }

        o.Append("CLUSTAL W (1.83) multiple sequence alignment\n");
        o.Append("\n\n");
        int iwidth = 60;
        for (int ipos = 1; ipos <= ilen; ipos += 60)
        {
            int iend = ipos + 60 - 1;
            if (iend > ilen)
            {
                iend = ilen;
                iwidth = ilen - ipos + 1;
            }

            foreach (var seq in seqs)
                o.Append(Left(seq.Name, 15, 15)).Append(' ').Append(GapStr(SubS(seq.Seq, ipos - 1, iend - 1), '-')).Append('\n');
            o.Append(new string(' ', 15)).Append(' ').Append(new string(' ', iwidth)).Append('\n');
            if (iend < ilen)
                o.Append('\n');
        }
    }

    private static void WriteSelex(StringBuilder o, EmbossSeqWriteRun run)
    {
        var seqs = run.Saved;
        int n = seqs.Count, nlen = 0, slen = 0, len = 0;
        var names = new string[n];
        var aseqs = new string[n];
        for (int i = 0; i < n; i++)
        {
            names[i] = seqs[i].Name;
            if ((len = names[i].Length) > nlen)
                nlen = len;
            if ((len = seqs[i].Seq.Length) > slen)
                slen = len;
            aseqs[i] = seqs[i].Seq;
        }

        for (int i = 0; i < n; i++)
        {
            names[i] += new string(' ', nlen - names[i].Length);
            aseqs[i] += new string(' ', slen - seqs[i].Seq.Length);
            o.Append("#=SQ ").Append(names[i]).Append(' ').Append(seqs[i].Weight.ToString("F2", CultureInfo.InvariantCulture)).Append(" - - 0..0:0 ");
            if (seqs[i].Desc.Length > 0)
                o.Append(seqs[i].Desc).Append('\n');
            else
                o.Append("-\n");
        }

        o.Append('\n');
        for (int i = 0; i < slen; i += 50)
        {
            for (int j = 0; j < n; j++)
            {
                string p = aseqs[j][i..];
                if (i + 50 >= len)
                    o.Append(names[j]).Append(' ').Append(p).Append('\n');
                else
                    o.Append(names[j]).Append(' ').Append(Left(p, 50, 50)).Append('\n');
            }

            if (i + 50 < len)
                o.Append('\n');
        }
    }

    private static void WriteMsf(StringBuilder o, EmbossSeqWriteRun run)
    {
        var seqs = run.Saved;
        int isize = seqs.Count, ilen = 0, maxnamelen = 10;
        foreach (var seq in seqs)
        {
            if (ilen < seq.Seq.Length)
                ilen = seq.Seq.Length;
            if (seq.Name.Length > maxnamelen)
                maxnamelen = seq.Name.Length;
        }

        int checktot = 0;
        foreach (var seq in seqs)
        {
            seq.Seq = GapSL(seq.Seq, '.', '~', ilen);
            checktot += CheckGcg(seq.Seq);
            checktot %= 10000;
        }

        string type = run.OutType.Length > 0 ? run.OutType : EmbossSequenceReader.IsNucleotide(seqs[0].Seq) ? "N" : "P";               // ajSeqType of the first sequence
        if (type.Length > 0 && type[0] == 'P')
        {
            o.Append("!!AA_MULTIPLE_ALIGNMENT 1.0\n\n");
            o.Append(CultureInfo.InvariantCulture, $"  {run.FileName} MSF:  {ilen} Type: P {Today(run)} CompCheck: {checktot,4} ..\n\n");
        }
        else
        {
            o.Append("!!NA_MULTIPLE_ALIGNMENT 1.0\n\n");
            o.Append(CultureInfo.InvariantCulture, $"  {run.FileName} MSF: {ilen} Type: N {Today(run)} CompCheck: {checktot,4} ..\n\n");
        }

        foreach (var seq in seqs)
            o.Append("  Name: ").Append(Left(seq.Name, maxnamelen, int.MaxValue)).Append(CultureInfo.InvariantCulture, $" Len: {seq.Seq.Length}  Check: {CheckGcg(seq.Seq),4} Weight: ").Append(seq.Weight.ToString("F2", CultureInfo.InvariantCulture)).Append('\n');
        o.Append("\n//\n\n");
        for (int ipos = 1; ipos <= ilen; ipos += 50)
        {
            int iend = Math.Min(ipos + 49, ilen);
            string sqbeg = ipos.ToString(CultureInfo.InvariantCulture), sqend = iend.ToString(CultureInfo.InvariantCulture);
            if (iend == ilen)
            {
                int igap = iend - ipos - sqbeg.Length;
                if (unchecked((uint)igap) >= (uint)sqend.Length)
                {
                    // "%*s %S %*S": a negative width (the unsigned gap that wrapped) left-justifies
                    string tail = igap >= 0 ? Right(sqend, igap) : Left(sqend, -igap, int.MaxValue);
                    o.Append(Right(" ", maxnamelen)).Append(' ').Append(sqbeg).Append(' ').Append(tail).Append('\n');
                }
                else
                {
                    o.Append("           ").Append(sqbeg).Append('\n');
                }
            }
            else
            {
                o.Append("           ").Append(Left(sqbeg, 25, int.MaxValue)).Append(Right(sqend, 25)).Append('\n');
            }

            foreach (var seq in seqs)
                o.Append(Left(seq.Name, maxnamelen, int.MaxValue)).Append(' ').Append(SubS(seq.Seq, ipos - 1, iend - 1)).Append('\n');
            o.Append('\n');
        }
    }
}
