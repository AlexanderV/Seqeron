#pragma warning disable S1121, S3358, S2325, S1066, S3776, S107, S1854, S1481, S2589, S3626 // the writers keep the control flow of ajseqwrite.c

using System.Globalization;
using System.Text;

namespace Seqeron.Genomics.Infrastructure;

// The sequence writers of the entry formats that print the annotation of the sequence (references, keywords, comments, taxonomy, dates): embl (emblold, emold), emblnew (embl, em), genbank, genpept, refseq, refseqp,
// swiss (swold), swissnew (swiss, sw, swissprot, uniprot). A5-8c3h, B09 F268.
internal static partial class EmbossSeqWrite
{
    // ajStrParseC: the non-empty parts between newlines
    private static IEnumerable<string> Parts(string s) => s.Split('\n', StringSplitOptions.RemoveEmptyEntries);

    private static void WrapLines(StringBuilder o, string prefix, string text, int width)
    {
        foreach (string line in Parts(EmbossStrFmt.Wrap(text, width)))
            o.Append(prefix).Append(line).Append('\n');
    }

    // the accession, keyword and classification lines: the items are joined by "; " and a line is closed when the next item does not fit
    private static void ListLines(StringBuilder o, string prefix, List<string> list, bool accession)
    {
        int ilen = 0;
        foreach (string cur in list)
        {
            if (accession ? ilen + cur.Length > 79 : ilen + cur.Length >= 79)
            {
                o.Append(";\n");
                ilen = 0;
            }

            if (ilen == 0)
            {
                o.Append(prefix);
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

        o.Append(accession ? ";\n" : ".\n");
    }

    private static string BaseCounts(string seq)
    {
        int a = 0, c = 0, g = 0, t = 0;
        foreach (char ch in seq)
        {
            switch (ch)
            {
                case 'A' or 'a': a++; break;
                case 'C' or 'c': c++; break;
                case 'G' or 'g': g++; break;
                case 'T' or 't' or 'U' or 'u': t++; break;
            }
        }

        return string.Create(CultureInfo.InvariantCulture, $"SQ   Sequence {seq.Length} BP; {a} A; {c} C; {g} G; {t} T; {seq.Length - a - c - g - t} other;\n");
    }

    private static void EmblSequence(StringBuilder o, string seq)
    {
        var sf = NewFormat(seq.Length);
        sf.EndStr = "\n//";
        sf.Tab = 4;
        sf.Spacer = 11;
        sf.Width = 60;
        sf.NumRight = true;
        sf.NumWidth = 9;
        sf.NumJust = true;
        WriteSeq(o, seq, sf);
    }

    private static void EmblTaxonomy(StringBuilder o, EmbossSeqOut s)
    {
        if (s.Tax.Length > 0)
            o.Append("OS   ").Append(s.Tax).Append(s.Taxcommon.Length > 0 ? " (" + s.Taxcommon + ")" : "").Append('\n');
        if (s.Taxlist.Count > 1)
            ListLines(o, "OC   ", s.Taxlist, false);
    }

    // seqWriteEmbl (the format before release 87: no XX lines)
    private static void WriteEmblOld(StringBuilder o, EmbossSeqOut s, EmbossSeqWriteRun run)
    {
        if (s.IsProtein)
            throw new NotSupportedException("The protein sequences of the embl writers are written by the swiss writer.");
        o.Append("ID   ").Append(Left(s.Name, 10, int.MaxValue)).Append(CultureInfo.InvariantCulture, $" standard; DNA; UNC; {s.Seq.Length} BP.\n");
        if (s.Acclist.Count > 0)
            ListLines(o, "AC   ", s.Acclist, true);
        if (s.Sv.Length > 0)
            o.Append("SV   ").Append(s.Sv).Append('\n');
        if (s.Desc.Length > 0)
            WrapLines(o, "DE   ", s.Desc, 75);
        if (s.Keylist.Count > 0)
            ListLines(o, "KW   ", s.Keylist, false);
        EmblTaxonomy(o, s);
        if (run.UfoLocal)
            o.Append(s.FeatureText);
        o.Append(BaseCounts(s.Seq));
        EmblSequence(o, s.Seq);
    }

    // seqWriteEmblnew
    private static void WriteEmblNew(StringBuilder o, EmbossSeqOut s, EmbossSeqWriteRun run)
    {
        if (s.IsProtein)
            throw new NotSupportedException("The protein sequences of the embl writers are written by the swiss writer.");
        string svstr = "1";
        if (s.Sv.Length > 0)
        {
            svstr = s.Sv;
            int ifind = svstr.IndexOf('.', StringComparison.Ordinal);
            if (ifind >= 0)
                svstr = svstr[(ifind + 1)..];
        }

        string idstr = s.Acc.Length > 0 ? s.Acc : s.Name;
        o.Append("ID   ").Append(idstr).Append("; SV ").Append(svstr).Append("; ").Append(s.Circular ? "circular" : "linear").Append(';');
        o.Append(' ').Append(EmbossSeqTables.MolEmbl(s.Molecule)).Append(';');
        o.Append(s.Class.Length > 0 ? " " + s.Class + ";" : " STD;");
        o.Append(s.Division.Length > 0 ? " " + s.Division + ";" : " UNC;");
        o.Append(CultureInfo.InvariantCulture, $" {s.Seq.Length} BP.\nXX\n");
        if (s.Acclist.Count > 0)
        {
            ListLines(o, "AC   ", s.Acclist, true);
            o.Append("XX\n");
        }

        if (s.Date is { Exists: true } date)
        {
            if (date.CreDate is { } cre)
                o.Append("DT   ").Append(EmbossSeqDate.Format(cre)).Append(" (Rel. ").Append(EmbossSeqDate.Str(date.CreRel)).Append(", Created)\n");
            else if (date.ModDate is { } mod0)
                o.Append("DT   ").Append(EmbossSeqDate.Format(mod0)).Append(" (Rel. ").Append(EmbossSeqDate.Str(date.ModRel)).Append(", Created)\n");
            if (date.ModDate is { } mod)
                o.Append("DT   ").Append(EmbossSeqDate.Format(mod)).Append(" (Rel. ").Append(EmbossSeqDate.Str(date.ModRel)).Append(", Last updated, Version ").Append(EmbossSeqDate.Str(date.ModVer)).Append(")\n");
            o.Append("XX\n");
        }

        if (s.Desc.Length > 0)
        {
            WrapLines(o, "DE   ", s.Desc, 75);
            o.Append("XX\n");
        }

        if (s.Keylist.Count > 0)
        {
            ListLines(o, "KW   ", s.Keylist, false);
            o.Append("XX\n");
        }

        EmblTaxonomy(o, s);
        if (s.Organelle.Length > 0)
            o.Append("OG   ").Append(s.Organelle).Append('\n');
        if (s.Tax.Length > 0 || s.Organelle.Length > 0 || s.Taxlist.Count > 1)
            o.Append("XX\n");
        foreach (var r in s.Reflist)
        {
            o.Append(CultureInfo.InvariantCulture, $"RN   [{r.Number}]\n");
            if (r.Comment.Length > 0)
                WrapLines(o, "RC   ", r.Comment, 75);
            if (r.Position.Length > 0)
                WrapLines(o, "RP   ", r.Position, 75);
            if (r.Xref.Length > 0)
                WrapLines(o, "RX   ", r.Xref + ".", 75);
            if (r.Groupname.Length > 0)
                WrapLines(o, "RG   ", r.Groupname, 75);
            if (r.Authors.Length > 0)
            {
                foreach (string line in Parts(EmbossStrFmt.WrapAt(r.Authors + ";", 75, ',')))
                    o.Append("RA   ").Append(line).Append('\n');
            }

            if (r.Title.Length > 0)
                WrapLines(o, "RT   ", "\"" + r.Title + "\";", 75);
            else
                o.Append("RT   ;\n");
            if (r.Location.Length > 0)
                WrapLines(o, "RL   ", r.Location + ".", 75);
            o.Append("XX\n");
        }

        if (s.Xreflist.Count > 0)
        {
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

            o.Append("XX\n");
        }

        foreach (string cmt in s.Cmtlist)
        {
            foreach (string line in Parts(EmbossStrFmt.WrapAt(cmt, 75, ',')))
                o.Append(line == " " ? "CC   \n" : "CC   " + line + "\n");
            o.Append("XX\n");
        }

        if (run.UfoLocal)
        {
            o.Append(s.FeatureText);
            o.Append("XX\n");
        }

        o.Append(BaseCounts(s.Seq));
        EmblSequence(o, s.Seq);
    }

    // seqWriteGenbank (kind 0), seqWriteGenpept (1), seqWriteRefseq (2), seqWriteRefseqp (3): the writers differ in the LOCUS line, the source line, the unit of the positions of a reference and the layout of the sequence
    private static void WriteGenbankLike(StringBuilder o, EmbossSeqOut s, EmbossSeqWriteRun run, int kind)
    {
        bool refseqp = kind == 3;
        if (refseqp)
            o.Append("LOCUS       ").Append(Left(s.Name, 9, int.MaxValue)).Append(' ').Append(Right(s.Seq.Length.ToString(CultureInfo.InvariantCulture), 7)).Append(" aa                   ");
        else
            o.Append("LOCUS       ").Append(Left(s.Name, 17, int.MaxValue)).Append(' ').Append(Right(s.Seq.Length.ToString(CultureInfo.InvariantCulture), 10)).Append(" bp   ");
        if (!refseqp)
        {
            o.Append(' ').Append(Left(EmbossSeqTables.MolGb(s.Molecule), 7, int.MaxValue));
            o.Append(' ').Append(Left(s.Circular ? "circular" : "linear", 8, int.MaxValue));
        }

        if (s.Division.Length > 0)
            o.Append(' ').Append(Left(EmbossSeqTables.DivGb(s.Division), 3, int.MaxValue));
        else if (kind == 0 && s.Class.Length > 0)
            o.Append(' ').Append(Left(EmbossSeqTables.ClsGb(s.Class), 3, int.MaxValue));
        else
            o.Append(" UNC");
        string pad = refseqp ? "       " : " ";
        if (s.Date is { } date)
        {
            if (date.ModDate is { } mod)
                o.Append(pad).Append(EmbossSeqDate.Format(mod));
            else if (date.CreDate is { } cre)
                o.Append(pad).Append(EmbossSeqDate.Format(cre));
        }
        else
        {
            o.Append(pad).Append(EmbossSeqDate.Format(run.RunDate));
        }

        o.Append('\n');
        if (s.Desc.Length > 0)
        {
            string tmp = s.Desc.EndsWith('.') ? s.Desc : s.Desc + ".";
            bool first = true;
            foreach (string line in Parts(EmbossStrFmt.Wrap(tmp, 67)))
            {
                o.Append(first ? "DEFINITION  " : "            ").Append(line).Append('\n');
                first = false;
            }
        }

        if (s.Acclist.Count > 0)
        {
            int ilen = 0;
            foreach (string cur in s.Acclist)
            {
                if (ilen == 0)
                {
                    o.Append("ACCESSION   ");
                    ilen = 11;
                }

                if (ilen + cur.Length > 79)
                {
                    o.Append("\n            ");
                    ilen = 11;
                }

                if (ilen > 11)
                    o.Append(' ');
                ilen += 1;
                o.Append(cur);
                ilen += cur.Length;
            }

            if (ilen > 0)
                o.Append('\n');
        }

        if (s.Sv.Length > 0)
        {
            if (s.Gi.Length > 0)
                o.Append("VERSION     ").Append(s.Sv).Append("  GI:").Append(s.Gi).Append('\n');
            else
                o.Append("VERSION     ").Append(s.Sv).Append('\n');
        }

        if (s.Keylist.Count > 0)
        {
            int ilen = 0;
            foreach (string cur in s.Keylist)
            {
                if (ilen == 0)
                {
                    o.Append("KEYWORDS    ");
                    ilen = 11;
                }

                if (ilen + cur.Length >= 79)
                {
                    o.Append(";\n            ");
                    ilen = 11;
                }

                if (ilen > 11)
                    o.Append("; ");
                ilen += 2;
                o.Append(cur);
                ilen += cur.Length;
            }

            o.Append(".\n");
        }

        if (s.Tax.Length > 0)
        {
            if (kind == 0 && s.Taxcommon.Length > 0)
                o.Append("SOURCE      ").Append(s.Tax).Append(" (").Append(s.Taxcommon).Append(")\n");
            else
                o.Append("SOURCE      ").Append(s.Tax).Append('\n');
            o.Append("  ORGANISM  ").Append(s.Tax).Append('\n');
            if (s.Taxlist.Count > 0)
            {
                int ilen = 0;
                foreach (string cur in s.Taxlist)
                {
                    if (ilen + cur.Length >= 79)
                    {
                        o.Append(";\n");
                        ilen = 0;
                    }

                    if (ilen == 0)
                    {
                        o.Append("            ");
                        ilen = 12;
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
        }

        foreach (var r in s.Reflist)
        {
            o.Append(CultureInfo.InvariantCulture, $"REFERENCE   {r.Number}");
            if (r.Position.Length > 0)
                o.Append("  (").Append(refseqp ? "residues " : "bases ").Append(r.Position.Replace("-", " to ", StringComparison.Ordinal)).Append(')');
            o.Append('\n');
            if (r.Authors.Length > 0)
            {
                bool first = true;
                foreach (string line in Parts(EmbossStrFmt.WrapAt(r.FmtAuthorsGb(), 68, ',')))
                {
                    o.Append(first ? "  AUTHORS   " : "            ").Append(line).Append('\n');
                    first = false;
                }
            }

            string title = r.FmtTitleGb();
            if (title.Length > 0)
            {
                bool first = true;
                foreach (string line in Parts(EmbossStrFmt.Wrap(title, 68)))
                {
                    o.Append(first ? "  TITLE     " : "            ").Append(line).Append('\n');
                    first = false;
                }
            }

            if (r.Location.Length > 0)
            {
                bool first = true;
                foreach (string line in Parts(EmbossStrFmt.Wrap(r.FmtLocationGb(), 68)))
                {
                    o.Append(first ? "  JOURNAL   " : "            ").Append(line).Append('\n');
                    first = false;
                }
            }

            if (r.Xref.Length > 0)
            {
                foreach (string line in Parts(EmbossStrFmt.Wrap(r.Xref, 75)))
                {
                    if (line.StartsWith("PUBMED; ", StringComparison.Ordinal))
                        o.Append("  PUBMED    ").Append(line[8..]).Append('\n');
                }
            }
        }

        bool firstcmt = true;
        foreach (string cmt in s.Cmtlist)
        {
            bool first = true;
            foreach (string line in Parts(EmbossStrFmt.WrapAt(cmt, 68, ',')))
            {
                if (first)
                {
                    if (firstcmt)
                    {
                        firstcmt = false;
                        o.Append("COMMENT     ").Append(line).Append('\n');
                    }
                    else
                    {
                        o.Append("            \n");
                        o.Append("            ").Append(line).Append('\n');
                    }

                    first = false;
                }
                else
                {
                    o.Append("            ").Append(line).Append('\n');
                }
            }
        }

        if (run.UfoLocal)
            o.Append(s.FeatureText);
        o.Append(refseqp ? "ORIGIN      \n" : "ORIGIN\n");
        var sf = NewFormat(s.Seq.Length);
        sf.EndStr = "\n//";
        sf.Tab = kind is 1 or 2 ? 1 : 0;
        sf.Spacer = 11;
        sf.Width = 60;
        sf.NumLeft = true;
        sf.NumWidth = kind is 1 or 2 ? 8 : 9;
        WriteSeq(o, s.Seq, sf);
    }
}
