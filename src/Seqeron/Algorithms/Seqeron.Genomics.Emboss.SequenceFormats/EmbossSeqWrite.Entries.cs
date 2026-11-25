// SPDX-License-Identifier: LGPL-2.1-or-later
// Port of the EMBOSS 6.6.0 AJAX / NUCLEUS libraries (library source files carry the GNU Lesser General Public License, version 2.1 or later).
// Part of Seqeron.Genomics.Emboss.SequenceFormats, licensed under the GNU Lesser General Public License, version 2.1 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
#pragma warning disable S1121, S3358, S2325, S1066, S3776, S107, S1854, S1481, S2589, S3626 // the writers keep the control flow of ajseqwrite.c

using System.Globalization;
using System.Text;

namespace Seqeron.Genomics.Emboss.SequenceFormats;

using Seqeron.Genomics.Infrastructure;

// The sequence writers of the entry formats that print the annotation of the sequence (references, keywords, comments, taxonomy, dates): embl (emblold, emold), emblnew (embl, em), genbank, genpept, refseq, refseqp,
// swiss (swold), swissnew (swiss, sw, swissprot, uniprot). A5-8c3h, B09 F268.
public static partial class EmbossSeqWrite
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
        {
            WriteSwissOld(o, s, run);                                    // seqWriteEmbl / seqWriteEmblnew: a protein sequence is written by seqWriteSwiss
            return;
        }

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
        {
            WriteSwissOld(o, s, run);                                    // seqWriteEmbl / seqWriteEmblnew: a protein sequence is written by seqWriteSwiss
            return;
        }

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

    // seqWriteExperiment: the Staden experiment file (EMBL like lines with EX for the description and AV for the accuracy values); a protein is written by the Swiss-Prot writer
    private static void WriteExperiment(StringBuilder o, EmbossSeqOut s, EmbossSeqWriteRun run)
    {
        if (s.IsProtein)
        {
            WriteSwissOld(o, s, run);
            return;
        }

        o.Append("ID   ").Append(Left(s.Name, 10, int.MaxValue)).Append(CultureInfo.InvariantCulture, $" standard; DNA; UNC; {s.Seq.Length} BP.\n");
        if (s.Acclist.Count > 0)
            ListLines(o, "AC   ", s.Acclist, true);
        if (s.Sv.Length > 0)
            o.Append("SV   ").Append(s.Sv).Append('\n');
        if (s.Desc.Length > 0)
            o.Append("EX   ").Append(s.Desc).Append('\n');
        if (s.Keylist.Count > 0)
            ListLines(o, "KW   ", s.Keylist, false);
        if (s.Tax.Length > 0)
            o.Append("OS   ").Append(s.Tax).Append('\n');
        if (s.Taxlist.Count > 0)
            ListLines(o, "OC   ", s.Taxlist, false);
        if (run.UfoLocal)
            o.Append(s.FeatureText);
        if (s.Accuracy is not null)
        {
            int ilen = s.Seq.Length;
            for (int i = 0; i < ilen; i += 20)
            {
                o.Append("AV  ");
                int jend = Math.Min(i + 20, ilen);
                for (int j = i; j < jend; j++)
                    o.Append(' ').Append(((int)AccuracyAt(s.Accuracy, j)).ToString(CultureInfo.InvariantCulture).PadLeft(2));
                o.Append('\n');
            }
        }

        o.Append(BaseCounts(s.Seq));
        var sf = NewFormat(s.Seq.Length);
        sf.EndStr = "\n//";
        sf.Tab = 4;
        sf.Spacer = 11;
        sf.Width = 60;
        sf.NumRight = true;
        sf.NumWidth = 9;
        sf.NumJust = true;
        WriteSeq(o, s.Seq, sf);
    }

    private static readonly string[] XrefTypes = ["unknown", "DRline", "db_xref", "ECnumber", "DescTag", "TaxId", "RXline"];

    private static string Yes(bool b) => b ? "Yes" : "No";

    // seqWriteDebug: a trace of the output object
    private static void WriteDebug(StringBuilder o, EmbossSeqOut s, EmbossSeqWriteRun run)
    {
        var d = run.Debug ?? throw new NotSupportedException("The program that writes the sequences gives no information for the debug writer.");
        var c = CultureInfo.InvariantCulture;
        o.Append("Sequence output trace\n=====================\n\n");
        o.Append("  Name: '").Append(DbName(s.Name, "")).Append("'\n");
        o.Append("  Accession: '").Append(s.Acc).Append("'\n");
        if (s.Acclist.Count > 0)
        {
            o.Append(c, $"  Acclist: ({s.Acclist.Count})");
            foreach (string a in s.Acclist)
                o.Append(' ').Append(a).Append('\n');
            o.Append('\n');
        }

        o.Append("  SeqVersion: '").Append(s.Sv).Append("'\n");
        o.Append("  GenInfo Id: '").Append(s.Gi).Append("'\n");
        o.Append("  Description: '").Append(s.Desc).Append("'\n");
        if (s.Keylist.Count > 0)
        {
            o.Append(c, $"  Keywordlist: ({s.Keylist.Count})\n");
            foreach (string k in s.Keylist)
                o.Append("    '").Append(k).Append("'\n");
        }

        o.Append("  Taxonomy: '").Append(s.Tax).Append("'\n");
        o.Append("  Taxcommon: '").Append(s.Taxcommon).Append("'\n");
        o.Append("  TaxId: '").Append(s.Taxid).Append("'\n");
        o.Append("  Organelle: '").Append(s.Organelle).Append("'\n");
        if (s.Taxlist.Count > 0)
        {
            o.Append(c, $"  Taxlist: ({s.Taxlist.Count})\n");
            foreach (string k in s.Taxlist)
                o.Append("    '").Append(k).Append("'\n");
        }

        o.Append("  Type: '").Append(s.Type).Append("'\n");
        o.Append("  Output type: '").Append(d.OutputType).Append("'\n");
        o.Append("  Molecule: '").Append(s.Molecule).Append("'\n");
        o.Append("  Class: '").Append(s.Class).Append("'\n");
        o.Append("  Division: '").Append(s.Division).Append("'\n");
        o.Append("  Database: '").Append(s.Db).Append("'\n");
        o.Append("  Set database: '").Append(s.Setdb).Append("'\n");
        o.Append("  Output database: '").Append(run.Setoutdb).Append("'\n");
        o.Append("  Full name: ''\n");
        if (s.Date is { } date)
        {
            if (date.CreDate is { } cre)
                o.Append("  Created: '").Append(EmbossSeqDate.Format(cre)).Append("' Rel. '").Append(EmbossSeqDate.Str(date.CreRel)).Append("' Ver. '").Append(EmbossSeqDate.Str(date.CreVer)).Append("'\n");
            if (date.ModDate is { } mod)
                o.Append("  Modified: '").Append(EmbossSeqDate.Format(mod)).Append("' Rel. '").Append(EmbossSeqDate.Str(date.ModRel)).Append("' Ver. '").Append(EmbossSeqDate.Str(date.ModVer)).Append("'\n");
            if (date.SeqDate is { } sd)
                o.Append("  Seq modified: '").Append(EmbossSeqDate.Format(sd)).Append("' Rel. '").Append(EmbossSeqDate.Str(date.SeqRel)).Append("' Ver. '").Append(EmbossSeqDate.Str(date.SeqVer)).Append("'\n");
        }

        if (s.Cmtlist.Count > 0)
        {
            o.Append(c, $"  Commentlist: ({s.Cmtlist.Count})\n");
            foreach (string k in s.Cmtlist)
                o.Append("    '").Append(k).Append("'\n");
        }

        if (s.Xreflist.Count > 0)
        {
            o.Append(c, $"  Xreflist: ({s.Xreflist.Count})\n");
            int i = 0;
            foreach (var x in s.Xreflist)
            {
                if (x.Type < XrefTypes.Length)
                    o.Append(c, $"{i++,3} Type: '{x.Type}' {XrefTypes[x.Type]}\n");
                else
                    o.Append(c, $"{i++,3} Type: '{x.Type}' ..BAD..\n");
                o.Append("       Db: '").Append(x.Db).Append("'\n");
                o.Append("       Id: '").Append(x.Id).Append("'\n");
                if (x.Secid.Length > 0)
                    o.Append("    Secid: '").Append(x.Secid).Append("'\n");
                if (x.Terid.Length > 0)
                    o.Append("    Terid: '").Append(x.Terid).Append("'\n");
                if (x.Quatid.Length > 0)
                    o.Append("   Quatid: '").Append(x.Quatid).Append("'\n");
                if (x.Start != 0 || x.End != 0)
                    o.Append(c, $"    Range: {x.Start}..{x.End}\n");
            }
        }

        if (s.Reflist.Count > 0)
        {
            o.Append(c, $"  Citationlist: ({s.Reflist.Count})\n");
            foreach (var r in s.Reflist)
            {
                o.Append(c, $"    Number: {r.Number}\n");
                if (r.Position.Length > 0)
                    o.Append("      Position: '").Append(r.Position).Append("'\n");
                if (r.Groupname.Length > 0)
                    o.Append("      Groupname: '").Append(r.Groupname).Append("'\n");
                if (r.Authors.Length > 0)
                    o.Append("      Authors: '").Append(r.Authors).Append("'\n");
                if (r.Title.Length > 0)
                    o.Append("      Title: '").Append(r.Title).Append("'\n");
                if (r.Comment.Length > 0)
                    o.Append("      Comment: '").Append(r.Comment).Append("'\n");
                if (r.Xref.Length > 0)
                    o.Append("      Xref: '").Append(r.Xref).Append("'\n");
                if (r.Location.Length > 0)
                    o.Append("      Location: '").Append(r.Location).Append("'\n");
                if (r.Loctype.Length > 0)
                    o.Append("      Loctype: '").Append(r.Loctype).Append("'\n");
            }
        }

        o.Append("  Usa: '").Append(d.Usa).Append("'\n");
        o.Append("  Ufo: '").Append(d.Ufo).Append("'\n");
        o.Append("  Input format: '").Append(d.InputFormat).Append("'\n");
        o.Append("  Output format: '").Append(d.OutputFormat).Append("'\n");
        o.Append("  Filename: '").Append(d.Filename).Append("'\n");
        o.Append("  Output directory: '").Append(d.Directory).Append("'\n");
        o.Append("  Entryname: '").Append(d.Entryname).Append("'\n");
        o.Append("  File name: '").Append(d.FileNameOpen).Append("'\n");
        o.Append("  Known file name: '<null>'\n");
        o.Append("  Extension: '").Append(d.Extension).Append("'\n");
        o.Append("  Single: '").Append(Yes(d.Single)).Append("'\n");
        o.Append("  Features: '").Append(Yes(d.Features)).Append("'\n");
        o.Append("  Rev: '").Append(Yes(s.Rev)).Append("'\n");
        o.Append("  Circular: '").Append(Yes(s.Circular)).Append("'\n");
        o.Append(c, $"  Offset: '{s.Offset}'\n");
        o.Append(c, $"  Count: '{run.HeaderCount}'\n");
        o.Append("  Documentation:...\n\n");
        o.Append("  Features filename: '<null>'\n");
        o.Append("  Features format: '<null>'\n");
        o.Append(c, $"  Feature table size: '{s.FeatureCount}'\n");
        if (s.Accuracy is not null)
        {
            int ilen = s.Seq.Length;
            o.Append("  Accuracy: \n");
            for (int i = 0; i < ilen; i += 20)
            {
                o.Append("    ");
                int jend = Math.Min(i + 20, ilen);
                for (int j = i; j < jend; j++)
                    o.Append(' ').Append(((int)AccuracyAt(s.Accuracy, j)).ToString(c).PadLeft(2));
                o.Append('\n');
            }
        }

        var sf = NewFormat(s.Seq.Length);
        sf.NumRight = true;
        sf.NumLeft = true;
        sf.NumJust = true;
        sf.Tab = 1;
        sf.Spacer = 11;
        sf.Width = 50;
        WriteSeq(o, s.Seq, sf);
    }
}
