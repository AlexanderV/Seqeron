#pragma warning disable S1121, S3358, S2325, S1066, S3776, S107, S1854, S1481, S2589, S3626, S2692, S1643 // the writers keep the control flow of ajseqwrite.c

using System.Globalization;
using System.Text;

namespace Seqeron.Genomics.Infrastructure;

// seqWriteSwiss (swold, swissold) and seqWriteSwissnew (swiss, sw, swissprot, uniprot, swnew): the Swiss-Prot / UniProt entry of a protein sequence. A5-8c3h, B09 F270.
internal static partial class EmbossSeqWrite
{
    private static readonly ulong[] Crc64Table = MakeCrc64Table();

    private static ulong[] MakeCrc64Table()
    {
        var t = new ulong[256];
        for (int i = 0; i < 256; i++)
        {
            ulong v = (ulong)i;
            for (int j = 0; j < 8; j++)
                v = (v & 1) != 0 ? 0xd800000000000000UL ^ (v >> 1) : v >> 1;
            t[i] = v;
        }

        return t;
    }

    // ajMathCrc64
    private static ulong Crc64(string seq)
    {
        ulong crc = 0;
        foreach (char c in seq)
            crc = Crc64Table[(crc ^ (ulong)(long)(sbyte)c) & 0xff] ^ (crc >> 8);
        return crc;
    }

    private static readonly double[] AaMass =
    [
        89.0894, 132.6068, 121.1494, 133.0992, 147.1261, 165.1872, 75.0625, 155.1517, 131.1700, 131.1700, 146.1847, 131.1700, 149.2032, 132.1144, 255.3124, 115.1273, 146.1413, 174.1981, 105.0888, 119.1157, 168.0494, 117.1432,
        204.2238, 136.8966, 181.1866, 146.6337,
    ];

    // ajSeqstrCalcMolwt
    private static float MolWeight(string seq)
    {
        double mw = 18.01057;
        foreach (char c in seq)
        {
            if (c == '-')
                continue;
            int i = char.ToUpperInvariant(c) - 'A';
            if (i is > 25 or < 0)
                i = 'X' - 'A';
            mw += AaMass[i] - 18.01057;
        }

        return (float)mw;
    }

    private static string SwissSq(string seq, bool nw)
    {
        ulong crc = Crc64(seq);
        int mw = (int)(0.5 + MolWeight(seq));
        string hi = ((crc >> 32) & 0xffffffff).ToString("X8", CultureInfo.InvariantCulture), lo = (crc & 0xffffffff).ToString("X8", CultureInfo.InvariantCulture);
        return nw ? string.Create(CultureInfo.InvariantCulture, $"SQ   SEQUENCE   {seq.Length} AA;  {mw} MW;  {hi}{lo} CRC64;\n") : string.Create(CultureInfo.InvariantCulture, $"SQ   SEQUENCE {seq.Length,5} AA; {mw,6} MW;  {hi}{lo} CRC64;\n");
    }

    private static void SwissSequence(StringBuilder o, string seq)
    {
        var sf = NewFormat(seq.Length);
        sf.EndStr = "\n//";
        sf.Tab = 4;
        sf.Spacer = 11;
        sf.Width = 60;
        WriteSeq(o, seq, sf);
    }

    private static void SwissSource(StringBuilder o, EmbossSeqOut s, int taxlimit)
    {
        if (s.Tax.Length > 0)
        {
            string tmp = s.Tax + (s.Taxcommon.Length > 0 ? " (" + s.Taxcommon + ")" : "") + ".";
            WrapLines(o, "OS   ", tmp, 75);
        }
    }

    // seqWriteSwiss
    private static void WriteSwissOld(StringBuilder o, EmbossSeqOut s, EmbossSeqWriteRun run)
    {
        if (s.Type.Length > 0 && s.Type[0] == 'N')
        {
            WriteEmblOld(o, s, run);
            return;
        }

        o.Append("ID   ").Append(Left(s.Name, 10, int.MaxValue)).Append(CultureInfo.InvariantCulture, $"     STANDARD;      PRT; {s.Seq.Length,5} AA.\n");
        if (s.Acclist.Count > 0)
            ListLines(o, "AC   ", s.Acclist, true);
        if (s.Desc.Length > 0)
            o.Append("DE   ").Append(s.Desc).Append('\n');
        SwissSource(o, s, 79);
        if (s.Taxlist.Count > 1)
            ListLines(o, "OC   ", s.Taxlist, false);
        if (s.Keylist.Count > 0)
            ListLines(o, "KW   ", s.Keylist, false);
        if (run.UfoLocal)
            o.Append(s.FeatureText);
        o.Append(SwissSq(s.Seq, false));
        SwissSequence(o, s.Seq);
    }

    private static void SubDescLines(StringBuilder o, EmbossSeqin.SubDesc sub, string first, string rest, string cont)
    {
        string alttext = "AltName:";
        if (sub.Name.Length > 0)
        {
            o.Append(first).Append(alttext).Append(" Full=").Append(sub.Name).Append(";\n");
            alttext = "        ";
        }

        foreach (var c in sub.Inn)
        {
            o.Append(first).Append(alttext).Append(" INN=").Append(c).Append(";\n");
            alttext = "        ";
        }

        foreach (var c in sub.Short)
            o.Append(cont).Append("Short=").Append(c).Append(";\n");
        foreach (var c in sub.Ec)
            o.Append(cont).Append("EC=").Append(c).Append(";\n");
        foreach (var c in sub.Allergen)
            o.Append(cont).Append("Allergen=").Append(c).Append(";\n");
        foreach (var c in sub.Biotech)
            o.Append(cont).Append("Biotech=").Append(c).Append(";\n");
        foreach (var c in sub.Cdantigen)
            o.Append(cont).Append("CD_antigen=").Append(c).Append(";\n");
        _ = rest;
    }

    // the Includes and Contains blocks (the Contains block writes the names of its SubNames with a different indentation)
    private static void IncludedBlock(StringBuilder o, EmbossSeqin.Desc d, string title, bool contains)
    {
        o.Append("DE   ").Append(title).Append('\n');
        o.Append("DE     RecName: Full=").Append(d.Name).Append(";\n");
        foreach (var c in d.Short)
            o.Append("DE              Short=").Append(c).Append(";\n");
        foreach (var c in d.Ec)
            o.Append("DE              EC=").Append(c).Append(";\n");
        foreach (var sub in d.AltNames)
            SubDescLines(o, sub, "DE     ", "", "DE              ");
        foreach (var sub in d.SubNames)
        {
            o.Append("DE     SubName: Full=").Append(sub.Name).Append(";\n");
            string ind = contains ? "DE            " : "DE              ";
            foreach (var c in sub.Short)
                o.Append(ind).Append("Short=").Append(c).Append(";\n");
            foreach (var c in sub.Ec)
                o.Append(ind).Append("EC=").Append(c).Append(";\n");
        }
    }

    // seqWriteSwissnew
    private static void WriteSwissNew(StringBuilder o, EmbossSeqOut s, EmbossSeqWriteRun run)
    {
        if (s.Type.Length > 0 && s.Type[0] == 'N')
        {
            WriteEmblOld(o, s, run);
            return;
        }

        if (s.Name.IndexOf('_', StringComparison.Ordinal) > 0)
            o.Append("ID   ").Append(Left(s.Name, 19, int.MaxValue)).Append(CultureInfo.InvariantCulture, $"     Reviewed;    {s.Seq.Length,8} AA.\n");
        else
            o.Append("ID   ").Append(Left(s.Name, 19, int.MaxValue)).Append(CultureInfo.InvariantCulture, $"     Unreviewed;    {s.Seq.Length,8} AA.\n");
        if (s.Acclist.Count > 0)
        {
            int ilen = 0;
            foreach (string cur in s.Acclist)
            {
                if (ilen + cur.Length > 73)
                {
                    o.Append(";\n");
                    ilen = 0;
                }

                if (ilen == 0)
                {
                    o.Append("AC   ");
                    ilen = 6;
                }
                else
                {
                    o.Append("; ");
                    ilen += 2;
                }

                o.Append(cur);
                ilen += cur.Length;
            }

            o.Append(";\n");
        }

        if (s.Date is { } date)
        {
            if (date.CreDate is { } cre)
                o.Append("DT   ").Append(EmbossSeqDate.Format(cre)).Append(", integrated into ").Append(EmbossSeqDate.Str(date.CreVer)).Append(".\n");
            if (date.SeqDate is { } seqd)
                o.Append("DT   ").Append(EmbossSeqDate.Format(seqd)).Append(", sequence version ").Append(EmbossSeqDate.Str(date.SeqVer)).Append(".\n");
            if (date.ModDate is { } mod)
                o.Append("DT   ").Append(EmbossSeqDate.Format(mod)).Append(", entry version ").Append(EmbossSeqDate.Str(date.ModVer)).Append(".\n");
        }

        if (s.Fulldesc is { } fd && fd.Name.Length > 0)
        {
            o.Append("DE   RecName: Full=").Append(fd.Name).Append(";\n");
            foreach (var c in fd.Short)
                o.Append("DE            Short=").Append(c).Append(";\n");
            foreach (var c in fd.Ec)
                o.Append("DE            EC=").Append(c).Append(";\n");
            foreach (var sub in fd.AltNames)
                SubDescLines(o, sub, "DE   ", "", "DE            ");
            foreach (var sub in fd.SubNames)
            {
                o.Append("DE     SubName: Full=").Append(sub.Name).Append(";\n");
                foreach (var c in sub.Short)
                    o.Append("DE              Short=").Append(c).Append(";\n");
                foreach (var c in sub.Ec)
                    o.Append("DE              EC=").Append(c).Append(";\n");
            }

            foreach (var d in fd.Includes)
                IncludedBlock(o, d, "Includes:", false);
            foreach (var d in fd.Contains)
                IncludedBlock(o, d, "Contains:", true);
            if (fd.Fragments != 0 || fd.Precursor)
            {
                o.Append("DE   Flags:");
                if (fd.Fragments == 1)
                    o.Append(" Fragment;");
                if (fd.Fragments == 2)
                    o.Append(" Fragments;");
                if (fd.Precursor)
                    o.Append(" Precursor;");
                o.Append('\n');
            }
        }
        else if (s.Desc.Length > 0)
        {
            o.Append("DE   ").Append(s.Desc).Append('\n');
        }

        if (s.Genelist.Count > 0)
        {
            string tmp = "";
            foreach (var g in s.Genelist)
            {
                if (tmp.Length > 0)
                    o.Append("GN   and\n");
                tmp = "";
                if (g.Name.Length > 0)
                    tmp += " Name=" + g.Name + ";";
                if (g.Synonyms.Length > 0)
                    tmp += " Synonyms=" + g.Synonyms + ";";
                if (g.Oln.Length > 0)
                    tmp += " OrderedLocusNames=" + g.Oln + ";";
                if (g.Orf.Length > 0)
                    tmp += " ORFNames=" + g.Orf + ";";
                if (tmp.Length > 0)
                    o.Append("GN  ").Append(tmp).Append('\n');
            }
        }

        SwissSource(o, s, 74);
        if (s.Organelle.Length > 0)
            o.Append("OG   ").Append(s.Organelle).Append(".\n");
        if (s.Taxlist.Count > 1)
        {
            int ilen = 0;
            foreach (string cur in s.Taxlist)
            {
                if (ilen + cur.Length >= 74)
                {
                    o.Append(";\n");
                    ilen = 0;
                }

                if (ilen == 0)
                {
                    o.Append("OC   ");
                    ilen = 6;
                }
                else
                {
                    o.Append("; ");
                    ilen += 2;
                }

                o.Append(cur);
                ilen += cur.Length;
            }

            o.Append(".\n");
        }

        if (s.Taxid.Length > 0)
            o.Append("OX   NCBI_TaxID=").Append(s.Taxid).Append(";\n");
        foreach (var r in s.Reflist)
        {
            o.Append(CultureInfo.InvariantCulture, $"RN   [{r.Number}]\n");
            if (r.Position.Length > 0)
                WrapLines(o, "RP   ", r.Position, 70);
            if (r.Comment.Length > 0)
                WrapLines(o, "RC   ", r.Comment, 70);
            if (r.Xref.Length > 0)
                WrapLines(o, "RX   ", r.Xref, 70);
            if (r.Groupname.Length > 0)
                WrapLines(o, "RG   ", r.Groupname, 70);
            if (r.Authors.Length > 0)
            {
                foreach (string line in Parts(EmbossStrFmt.WrapAt(r.Authors + ";", 70, ',')))
                    o.Append("RA   ").Append(line).Append('\n');
            }

            if (r.Title.Length > 0)
                WrapLines(o, "RT   ", "\"" + r.Title + "\";", 70);
            if (r.Location.Length > 0)
                WrapLines(o, "RL   ", r.Location + ".", 70);
        }

        foreach (string cmt in s.Cmtlist)
        {
            foreach (string line in Parts(cmt))
                o.Append(line == " " ? "CC   \n" : "CC   " + line + "\n");
        }

        foreach (var x in s.Xreflist)
        {
            if (x.Type != 1)
                continue;
            if (x.Quatid.Length > 0)
                o.Append("DR   ").Append(x.Db).Append("; ").Append(x.Id).Append("; ").Append(x.Secid).Append("; ").Append(x.Terid).Append("; ").Append(x.Quatid).Append(".\n");
            else if (x.Terid.Length > 0)
                o.Append("DR   ").Append(x.Db).Append("; ").Append(x.Id).Append("; ").Append(x.Secid).Append("; ").Append(x.Terid).Append(".\n");
            else if (x.Secid.Length > 0)
                o.Append("DR   ").Append(x.Db).Append("; ").Append(x.Id).Append("; ").Append(x.Secid).Append(".\n");
            else
                o.Append("DR   ").Append(x.Db).Append("; ").Append(x.Id).Append(".\n");
        }

        if (s.Evidence.Length > 0)
            o.Append("PE   ").Append(s.Evidence).Append('\n');
        if (s.Keylist.Count > 0)
        {
            int ilen = 0;
            foreach (string cur in s.Keylist)
            {
                if (ilen + cur.Length >= 74)
                {
                    o.Append(";\n");
                    ilen = 0;
                }

                if (ilen == 0)
                {
                    o.Append("KW   ");
                    ilen = 6;
                }
                else
                {
                    o.Append("; ");
                    ilen += 2;
                }

                o.Append(cur);
                ilen += cur.Length;
            }

            o.Append(".\n");
        }

        if (run.UfoLocal)
            o.Append(s.FeatureText);
        o.Append(SwissSq(s.Seq, true));
        SwissSequence(o, s.Seq);
    }
}
