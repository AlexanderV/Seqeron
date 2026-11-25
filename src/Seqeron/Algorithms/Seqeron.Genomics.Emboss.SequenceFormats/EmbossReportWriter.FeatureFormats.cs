// SPDX-License-Identifier: LGPL-2.1-or-later
// Port of the EMBOSS 6.6.0 AJAX / NUCLEUS libraries (library source files carry the GNU Lesser General Public License, version 2.1 or later).
// Part of Seqeron.Genomics.Emboss.SequenceFormats, licensed under the GNU Lesser General Public License, version 2.1 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using System.Globalization;
using System.Text;

namespace Seqeron.Genomics.Emboss.SequenceFormats;

using Seqeron.Genomics.Alignment;

// The feature-table formats of ajfeatwrite.c (embl, genbank, pir, swiss, dasgff) and the debug trace of ajreport.c / ajfeat.c.
public static partial class EmbossReportWriter
{
    // ------------------------------------------------------------------------------------ embl / genbank ----

    // featTagEmblWrapC: single spaces only, then wrap at the margin on the last space (a word longer than the width is cut)
    private static void TagWrap(string value, int margin, string prefix, StringBuilder sb)
    {
        string pval = CollapseWhite(value);
        int width = margin - prefix.Length;
        int len = pval.Length;
        int k;
        for (int i = 0; i < len; i += k)
        {
            int last = i + width - 1;
            string chunk;
            int j;
            if (last + 1 >= len)
            {
                chunk = pval[i..];
                j = 0;
            }
            else if (pval[last + 1] == ' ')
            {
                chunk = pval.Substring(i, last - i + 1);
                j = last + 1 - i;
            }
            else
            {
                chunk = pval.Substring(i, last - i + 1);
                j = chunk.LastIndexOf(' ');
            }

            if (j < 1)
            {
                j = chunk.Length;
                k = j;
            }
            else
                k = j + 1;
            sb.Append(prefix).Append(chunk, 0, Math.Min(j, chunk.Length)).Append('\n');
        }
    }

    // ajStrRemoveWhiteExcess: leading / trailing white space removed, inner runs reduced to one blank
    private static string CollapseWhite(string s)
    {
        var sb = new StringBuilder();
        bool pending = false;
        foreach (char c in s)
        {
            if (c is ' ' or '\t' or '\n' or '\r' or '\f' or '\v')
            {
                pending = sb.Length > 0;
                continue;
            }

            if (pending)
                sb.Append(' ');
            pending = false;
            sb.Append(c);
        }

        return sb.ToString();
    }

    // featLocEmblWrapC: no white space; split after the last comma that fits (a location without comma is cut at the width)
    private static void LocWrap(string loc, int margin, string prefix, string preftyp, StringBuilder sb)
    {
        int width = margin - prefix.Length;
        string ploc = new([.. loc.Where(c => !char.IsWhiteSpace(c))]);
        int len = ploc.Length;
        int k;
        for (int i = 0; i < len; i += k)
        {
            int last = i + width - 1;
            string chunk;
            int j;
            if (last + 1 >= len)
            {
                chunk = ploc[i..];
                j = 0;
            }
            else
            {
                chunk = ploc.Substring(i, last - i + 1);
                j = chunk.LastIndexOf(',');
            }

            if (j < 1)
            {
                j = chunk.Length;
                k = j;
            }
            else
            {
                j++;
                k = j;
            }

            sb.Append(i == 0 ? preftyp : prefix).Append(chunk, 0, Math.Min(j, chunk.Length)).Append('\n');
        }
    }

    // featTagQuoteEmbl: the value in double quotes, an inner quote doubled
    private static string QuoteEmbl(string v) => "\"" + v.Replace("\"", "\"\"", StringComparison.Ordinal) + "\"";

    // the EMBL tag vocabulary (Etags.embl) for the tags of the programs: note is QTEXT, rpt_type LIMITED
    private static readonly string[] RptTypes = ["TANDEM", "INVERTED", "FLANKING", "TERMINAL", "DIRECT", "DISPERSED", "OTHER"];

    private static void Embl(State st, EmbossReport t, bool isEmbl)
    {
        var sb = st.Out;
        if (t.Protein)
            throw new EmbossReportException($"Report format '{FormatName(st)}' not supported for protein sequences");
        sb.Append(isEmbl ? "FH   Key             Location/Qualifiers\nFH\n" : "FEATURES             Location/Qualifiers\n");
        string lead = isEmbl ? "FT" : "  ";
        int margin = isEmbl ? 80 : 79;
        string preftag = lead + "                   ";
        foreach (var f in t.Features)
        {
            var ti = TypeOf(f);
            string typ = ti.EmblName ?? "misc_feature";
            string preftyp = lead + "   " + (typ.Length > 15 ? typ[..15] : typ).PadRight(15) + " ";
            string pos = Int(f.Start) + (f.Between ? "^" : "..") + Int(f.End);
            string loc = f.Strand == '-' ? "complement(" + pos + ")" : pos;
            LocWrap(loc, margin, preftag, preftyp, sb);
            if (ti.EmblName is null)
                TagWrap($"/note=\"*Type {ti.Id} {ti.Name}\"\n", margin, preftag, sb);
            foreach (var tag in f.Tags)
            {
                string name = tag.Name.ToLowerInvariant();
                string line;
                if (name == "rpt_type")
                {
                    string val = tag.Value;
                    foreach (string lim in RptTypes)
                        if (string.Equals(lim, val, StringComparison.OrdinalIgnoreCase))
                        {
                            val = lim;
                            break;
                        }

                    line = "/rpt_type=" + val + "\n";
                }
                else if (name == "note")
                    line = "/note=" + QuoteEmbl(tag.Value) + "\n";
                else if (name == "replace")
                    line = ReplaceTag(tag.Value, true);
                else
                    throw new ArgumentException($"The feature tag '{tag.Name}' is not in the tag vocabulary of the report layer (note, rpt_type, replace).");
                TagWrap(line, margin, preftag, sb);
                st.Leftover = CollapseWhite(line);
            }
        }
    }

    // the QSPECIAL tag replace (featTagSpecialAllReplace: white space removed, then bases and "/" "-" only): the quoted value, or the note "replace: value" the writers make of a bad one
    private static string ReplaceTag(string value, bool embl)
    {
        string v = new([.. value.Where(c => !char.IsWhiteSpace(c))]);
        if (v.All(c => "abcdghkmnrstuvwxyABCDGHKMNRSTUVWXY/-".Contains(c, StringComparison.Ordinal)))
            return "/replace=" + QuoteEmbl(v) + "\n";
        return embl ? "/note=\"replace: " + QuoteEmbl(v) + "\"" : "";
    }

    // ----------------------------------------------------------------------------------------------- pir ----

    private static void Pir(State st, EmbossReport t)
    {
        var sb = st.Out;
        if (!t.Protein)
            throw new EmbossReportException($"Report format '{FormatName(st)}' not supported for nucleotide sequences");
        foreach (var f in t.Features)
        {
            var ti = TypeOf(f);
            string location = Int(f.Start) + (f.End != f.Start ? "-" + Int(f.End) : "");
            string outtyp = ti.PirType.Length > 0 ? ti.PirType : throw new ArgumentException($"The feature type '{f.Type}' has no PIR type.");
            // featTypePirOut: when the external PIR type differs from the internal one the internal type name is appended to the first note
            var tags = f.Tags.ToList();
            int note = tags.FindIndex(g => g.Name.Equals("note", StringComparison.OrdinalIgnoreCase));
            if (ti.Id == "SO:0000839")
            {
                // diffseq's generic feature of a protein table: its PIR type Region is the external type, nothing is appended
            }
            else if (note >= 0)
                tags[note] = tags[note] with { Value = tags[note].Value + ", " + ti.Name };
            else
                tags.Insert(0, new EmbossFeatureTag("note", ti.Name));
            sb.Append("F;").Append(location).Append('/').Append(outtyp).Append(':');
            var comment = new StringBuilder();
            foreach (var tag in tags)
            {
                string name = tag.Name.ToLowerInvariant();
                if (name == "comment")
                {
                    comment.Append(" #").Append(tag.Value);
                    continue;
                }

                sb.Append(st.Leftover).Append(' ').Append(tag.Value);
                st.Leftover = "";
            }

            sb.Append(comment).Append('\n');
        }
    }

    // --------------------------------------------------------------------------------------------- swiss ----

    private static void Swiss(State st, EmbossReport t)
    {
        var sb = st.Out;
        if (!t.Protein)
            throw new EmbossReportException($"Report format '{FormatName(st)}' not supported for nucleotide sequences");
        foreach (var f in t.Features)
        {
            var ti = TypeOf(f);
            string typ = ti.SwissType.Length > 0 ? ti.SwissType : throw new ArgumentException($"The feature type '{f.Type}' has no Swiss-Prot type.");
            var line = new StringBuilder("FT   ");
            line.Append((typ.Length > 8 ? typ[..8] : typ).PadRight(8)).Append(' ');
            string from = Int(f.Start), to = Int(f.End);
            line.Append((from.Length > 6 ? from[..6] : from).PadLeft(6)).Append(' ').Append((to.Length > 6 ? to[..6] : to).PadLeft(6));
            int i = 0;
            foreach (var tag in f.Tags)
            {
                line.Append(i++ > 0 ? " " : "       ");
                if (!tag.Name.Equals("note", StringComparison.OrdinalIgnoreCase))
                    throw new ArgumentException($"The feature tag '{tag.Name}' is not in the tag vocabulary of the report layer (note).");
                line.Append(tag.Value).Append('.');                          // TEXT: simple text, wrap at space
            }

            SwissWrap(line.ToString(), 75, "FT                                ", sb);
            st.Leftover = "";
        }
    }

    // ---------------------------------------------------------------------------------------------- dasgff ----

    private static void DasGff(State st, EmbossReport t)
    {
        var sb = st.Out;
        if (!st.Opened)
        {
            st.Opened = true;
            sb.Append("<?xml version=\"1.0\" standalone=\"no\"?>\n<DASGFF>\n  <GFF version=\"1.0\" href=\"url\">\n");
        }

        sb.Append("    <SEGMENT id=\"").Append(t.SeqId).Append("\" start=\"").Append(Int(1 + t.Offset + t.Start)).Append("\" stop=\"").Append(Int(1 + t.Offset + t.End)).Append("\"\n");
        sb.Append("                version=\"0.0\">\n");
        int group = 0;
        foreach (var f in t.Features)
        {
            group++;
            var ti = TypeOf(f);
            sb.Append("      <FEATURE id=\"").Append(t.SeqId).Append('.').Append(f.Source).Append('.').Append(group).Append("\">\n");
            sb.Append("         <TYPE id=\"").Append(ti.Id).Append("\" category=\"").Append(ti.DasCategory).Append("\" reference=\"no\" subparts=\"no\">").Append(ti.Gff3Name).Append("</TYPE>\n");
            if (f.Source.Length > 0)
                sb.Append("         <METHOD id=\"").Append(f.Source).Append("\"></METHOD>\n");
            else
                sb.Append("         <METHOD></METHOD>\n");
            sb.Append("         <START>").Append(Int(f.Start)).Append("</START>\n");
            sb.Append("         <END>").Append(Int(f.End)).Append("</END>\n");
            if ((float)f.Score != 0f)
                sb.Append("         <SCORE>").Append(CFormatFixed6(f.Score)).Append("</SCORE>\n");
            if (f.Strand != '\0')
                sb.Append("         <ORIENTATION>").Append(f.Strand).Append("</ORIENTATION>\n");
            else
                sb.Append("         <ORIENTATION>0</ORIENTATION>\n");
            if (f.Frame == 0)
                sb.Append("         <PHASE>-</PHASE>\n");
            else
                sb.Append("         <PHASE>").Append(f.Frame is > 0 and <= 3 ? (char)('0' + f.Frame - 1) : '.').Append("</PHASE>\n");
            foreach (var tag in f.Tags)
            {
                string name = tag.Name.ToLowerInvariant();
                sb.Append("         <NOTE>");
                string v = tag.Value;
                string note;
                if (name == "note" && v.Length > 0 && v[0] == '*')
                {
                    v = v[1..];
                    int sp = v.IndexOf(' ');
                    string w = sp < 0 ? v : v[..sp];
                    string rest = sp < 0 ? "" : v[(sp + 1)..].TrimStart();
                    note = w + ":" + rest;
                }
                else
                    note = name + ":" + v;
                sb.Append(note);
                st.Leftover = note;
                sb.Append("</NOTE>\n");
            }

            sb.Append("      </FEATURE>\n");
        }

        sb.Append("    </SEGMENT>\n");
    }

    private static string CFormatFixed6(double v) => Infrastructure.CFormat.Fixed((float)v, 6);

    // ---------------------------------------------------------------------------------------------- debug ----

    // reportWriteTrace + ajFeattablePrint
    private static void DebugReport(State st, EmbossReport t)
    {
        var sb = st.Out;
        var o = st.O;
        Header(st, t, t.Sequence);
        sb.Append("Trace output\n");
        sb.Append("Name: '").Append(o.Name).Append("'\n");
        sb.Append("Type: '").Append(st.Type).Append("'\n");
        sb.Append("Formatstr: '").Append(FormatName(st)).Append("'\n");
        sb.Append("Feattable (internal): 0 features\n");
        sb.Append("FtQuery: '<null>'\n");
        sb.Append("Extension: '").Append(o.Extension).Append("'\n");
        sb.Append("File: '").Append(o.ReportFile).Append("'\n");
        sb.Append("Tagnames: ").Append(o.Tags.Count).Append('\n');
        sb.Append("Tagprints: ").Append(o.Tags.Count).Append('\n');
        sb.Append("Tagtypes: ").Append(o.Tags.Count).Append('\n');
        sb.Append("FileNames: ").Append(st.O.AdditionalFiles.Count).Append("\nFileTypes: ").Append(st.O.AdditionalFiles.Count).Append('\n');
        sb.Append("Precision: ").Append(o.Precision).Append('\n');
        sb.Append("Showacc: ").Append(YesNo(o.ShowAccession)).Append('\n');
        sb.Append("Showdes: ").Append(YesNo(o.ShowDescription)).Append('\n');
        sb.Append("Showusa: ").Append(YesNo(st.ShowUsa)).Append('\n');
        sb.Append("Showscore: ").Append(YesNo(o.ShowScore)).Append('\n');
        sb.Append("Showstrand: ").Append(YesNo(o.ShowStrand)).Append('\n');
        sb.Append("Multi: ").Append(o.Multi ? "Yes" : "No").Append("\nMintags: 0\n");
        sb.Append("CountSeq: ").Append(st.CountSeq).Append('\n');
        sb.Append("CountHit: ").Append(st.CountHit).Append('\n');
        sb.Append("TotHits: ").Append(st.TotHits).Append('\n');
        sb.Append("MaxHitAll: ").Append(o.MaxAll).Append('\n');
        sb.Append("MaxHitSeq: ").Append(o.MaxSeq).Append('\n');
        sb.Append("Format: ").Append(st.FormatIndex).Append('\n');
        if (t.FeatureTable is not null)
        {
            // ajFeattablePrint of the full table
            sb.Append("\nFeattable (external): ").Append(t.FeatureTable.Features.Count).Append(" features\n");
            sb.Append(t.FeatureTable.Debug());
            Tail(st, t, t.Sequence);
            return;
        }

        sb.Append("\nFeattable (external): ").Append(t.Features.Count).Append(" features\n");
        sb.Append("  Format: 0 (").Append(t.TableFormat).Append(")  Start: ").Append(t.Start).Append('\n');
        sb.Append("  End: ").Append(t.End).Append('\n');
        sb.Append("  Len: ").Append(t.TableLength).Append('\n');
        sb.Append("  Offset: ").Append(t.Offset).Append('\n');
        sb.Append("  Groups: ").Append(t.Groups ?? t.Features.Count).Append('\n');
        int n = 0;
        foreach (var f in t.Features)
        {
            n++;
            sb.Append("\n  Feature ").Append(n).Append('\n');
            sb.Append("    Source: '").Append(f.Source).Append("'\n");
            sb.Append("    Type: '").Append(TypeOf(f).Id).Append("'\n");
            sb.Append("    Score: ").Append(Infrastructure.CFormat.Fixed((float)f.Score, 6)).Append('\n');
            sb.Append("    Protein: ").Append(YesNo(f.Protein ?? t.Protein)).Append('\n');
            sb.Append("    Strand: '").Append(f.Strand).Append("'\n");
            sb.Append("    Start: ").Append(f.Start).Append('\n');
            sb.Append("    End: ").Append(f.End).Append('\n');
            sb.Append("    Start2: 0\n    End2: 0\n    Remote: '<null>'\n    Label: '<null>'\n");
            sb.Append("    Frame: ").Append(f.Frame).Append('\n');
            sb.Append("    Exon: 0\n");
            sb.Append("    Group: ").Append(n).Append('\n');
            sb.Append(f.Between ? "    Flags: 8\n      BETWEEN_SEQ\n" : "    Flags: 0\n");
            sb.Append("    Tags: ").Append(f.Tags.Count).Append(" tags\n");
            int j = 0;
            foreach (var tag in f.Tags)
                sb.Append("      Tag ").Append((++j).ToString(CultureInfo.InvariantCulture).PadLeft(3)).Append(' ').Append(tag.Name).Append(" : '").Append(tag.Value).Append("'\n");
        }

        Tail(st, t, t.Sequence);
    }

    private static string YesNo(bool b) => b ? "Yes" : "No";
}
