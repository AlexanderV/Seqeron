// SPDX-License-Identifier: LGPL-2.1-or-later
// Port of the EMBOSS 6.6.0 AJAX / NUCLEUS libraries (library source files carry the GNU Lesser General Public License, version 2.1 or later).
// Part of Seqeron.Genomics.Emboss.SequenceFormats, licensed under the GNU Lesser General Public License, version 2.1 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using System.Globalization;
using System.Text;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Emboss.SequenceFormats;

using Seqeron.Genomics.Alignment;


/// <summary>The static string <c>featoutStr</c> of <c>ajfeatwrite.c</c> that the feature writers of one program run share: what the last tag of one writer left in it is printed in front of the first tag of the next writer that appends to it (gff2, pir).</summary>
public sealed class EmbossFeatureWriteState
{
    /// <summary>The current content of the shared string.</summary>
    public string Leftover { get; set; } = "";
}

// The FEATURE OUTPUT formats of ajfeatwrite.c (ajFeattableWrite and the ajFeattableWrite* of the output format table) for the feature tables the programs write with -aoutfeat / -boutfeat (diffseq's conflicts, wordmatch's matches).
public static partial class EmbossReportWriter
{
    /// <summary>The names of the feature output formats (<c>featOutFormatDef</c>): name, whether it is valid for nucleotide and for protein sequences.</summary>
    private static readonly (string Name, string Writer, bool Nuc, bool Protein)[] FeatureFormatTable =
    [
        ("unknown", "gff3", false, false), ("gff", "gff3", true, true), ("gff2", "gff2", true, true), ("gff3", "gff3", true, true), ("bed", "bed", true, true), ("embl", "embl", true, false), ("em", "embl", true, false),
        ("genbank", "genbank", true, false), ("gb", "genbank", true, false), ("ddbj", "genbank", true, false), ("refseq", "genbank", true, false), ("refseqp", "refseqp", false, true), ("pir", "pir", false, true),
        ("nbrf", "pir", false, true), ("swiss", "swiss", false, true), ("sw", "swiss", false, true), ("swissprot", "swiss", false, true), ("uniprot", "swiss", false, true), ("dasgff", "dasgff", true, true),
        ("draw", "draw", true, true), ("debug", "debug", true, true),
    ];

    /// <summary><c>featFindOutFormatS</c>: the index of the format in the output format table (the name is compared in lower case), -1 when there is none.</summary>
    public static int FeatureFormatIndex(string name)
    {
        ArgumentNullException.ThrowIfNull(name);
        if (name.Length == 0)
            return -1;
        string lower = name.ToLowerInvariant();
        for (int i = 0; i < FeatureFormatTable.Length; i++)
            if (FeatureFormatTable[i].Name == lower)
                return i;
        return -1;
    }

    /// <summary><c>featOutFormat[i].Name</c>.</summary>
    public static string FeatureFormatName(int index) => FeatureFormatTable[index].Name;

    /// <summary>
    /// <c>ajFeattableWrite</c> (<paramref name="tablesBefore"/> = the tables already written to the same output file, <c>ftout->Count</c>: dasgff writes its XML header for the first table only): the text of the feature table in the output format at <paramref name="formatIndex"/> (see <see cref="FeatureFormatIndex"/>), or null when the program refuses the format for the kind of sequence
    /// ("Feature format 'embl' not supported for protein sequence features": nothing is written).
    /// </summary>
    public static string? WriteFeatureTable(EmbossReport table, int formatIndex, DateTime date, string program = "diffseq", EmbossFeatureWriteState? state = null, int tablesBefore = 0)
    {
        ArgumentNullException.ThrowIfNull(table);
        var (name, writer, nuc, protein) = FeatureFormatTable[formatIndex];
        if (!(table.Protein ? protein : nuc))
            return null;
        var options = new EmbossReportOptions { Program = program, RunDate = date, Tags = table.Tags ?? [], Format = EmbossReportFormat.Gff };
        state ??= new EmbossFeatureWriteState();
        var sb = new StringBuilder();
        switch (writer)
        {
            case "gff3":
                return Write(table, EmbossReportFormat.Gff, date, options);
            case "embl":
                return Output(options with { Format = EmbossReportFormat.Embl }, table, state, (st, t) => Embl(st, t, true));
            case "genbank":
                return Output(options with { Format = EmbossReportFormat.Genbank }, table, state, (st, t) => Embl(st, t, false));
            case "pir":
                return Output(options with { Format = EmbossReportFormat.Pir }, table, state, Pir);
            case "swiss":
                return Output(options with { Format = EmbossReportFormat.Swiss }, table, state, Swiss);
            case "dasgff":
                return Output(options with { Format = EmbossReportFormat.DasGff }, table, state, DasGff, tablesBefore > 0);          // the closing </GFF></DASGFF> is written by the report only

            case "gff2":
                Gff2(sb, table, date, state);
                return sb.ToString();
            case "bed":
                Bed(sb, table);
                return sb.ToString();
            case "draw":
                Draw(sb, table);
                return sb.ToString();
            case "debug":
                FeatureTablePrint(sb, table);
                return sb.ToString();
            case "refseqp":
                RefseqP(sb, table, state);
                return sb.ToString();
            default:
                throw new InvalidOperationException($"No writer for the feature output format '{name}'.");
        }
    }

    private static string Output(EmbossReportOptions options, EmbossReport table, EmbossFeatureWriteState state, Action<State, EmbossReport> body, bool opened = false)
    {
        var st = new State { O = options, FormatIndex = EmbossReportFormats.IndexOf(options.Format), Leftover = state.Leftover, Opened = opened };
        body(st, table);
        state.Leftover = st.Leftover;
        return st.Out.ToString();
    }

    // ----------------------------------------------------------------------------------------------- gff2 ----

    // featTagQuoteGff2: every inner quote escaped with a backslash, then the value in quotes
    private static string QuoteGff2(string v)
    {
        string q = v.Replace("\"", "\\\"", StringComparison.Ordinal);
        if (q.Length == 0 || q[0] != '"')
            q = "\"" + q;
        if (q.Length == 1 || q[^1] != '"')
            q += "\"";
        return q;
    }

    private static void Gff2(StringBuilder sb, EmbossReport t, DateTime date, EmbossFeatureWriteState state)
    {
        sb.Append("##gff-version 2.0\n");
        sb.Append("##date ").Append(EmbossDates.Gff(date)).Append('\n');
        string seqid = t.SeqId;
        if (t.TableType == "N")
            sb.Append("##Type DNA ").Append(seqid).Append('\n');
        else if (t.TableType == "P")
            sb.Append("##Type Protein ").Append(seqid).Append('\n');
        else
            sb.Append("##Type unknown <").Append(t.TableType).Append("> ").Append(seqid).Append('\n');
        int group = 0;
        foreach (var f in t.Features)
        {
            group++;
            var ti = TypeOf(f);
            char strand = FeatStrand(f, t.Protein);
            int a = f.Start, b = f.End;
            if (strand == '-' && f.End < f.Start)
            {
                a = f.End;
                b = f.Start;
            }

            sb.Append(seqid).Append('\t').Append(f.Source).Append('\t').Append(ti.Gff2Name ?? ti.Gff3Name).Append('\t').Append(Int(a)).Append('\t').Append(Int(b)).Append('\t');
            sb.Append((float)f.Score != 0f ? CFormat.General((float)f.Score, 4) : ".");
            sb.Append('\t').Append(strand).Append('\t').Append(f.Frame is > 0 and <= 3 ? (char)('0' + f.Frame - 1) : '.').Append('\t');
            string flagdata = f.Between ? "0x8" : "";
            sb.Append("Sequence \"").Append(seqid).Append('.').Append(group).Append('"');
            if (flagdata.Length > 0)
                sb.Append(" ; FeatFlags \"").Append(flagdata).Append('"');
            foreach (var tag in f.Tags)
            {
                string name = tag.Name.ToLowerInvariant();
                string text;
                if (name == "replace")
                {
                    string v = new([.. tag.Value.Where(c => !char.IsWhiteSpace(c))]);
                    text = v.All(c => "abcdghkmnrstuvwxyABCDGHKMNRSTUVWXY/-".Contains(c, StringComparison.Ordinal)) ? "replace " + QuoteGff2(v) : "note \"replace: " + QuoteGff2(v) + "\"";
                }
                else
                    text = name + " " + QuoteGff2(tag.Value);
                sb.Append(" ; ").Append(state.Leftover).Append(text);
                state.Leftover = "";
            }

            sb.Append('\n');
        }
    }

    // ----------------------------------------------------------------------------------------------- bed ----

    private static void Bed(StringBuilder sb, EmbossReport t)
    {
        // ajFeattableGetScorerange: the scores are never used (the range stays 0 .. 0 for a score of 0 or more: the program assigns curmax = 0)
        sb.Append("track name=").Append(t.SeqId).Append(" description=\"EMBOSS:6.6.0.0\" useScore=0\n");
        foreach (var f in t.Features)
        {
            var ti = TypeOf(f);
            char strand = FeatStrand(f, t.Protein);
            int istart, iend;                                           // ajuint printed with %d
            if (strand == '-' && f.End < f.Start)
            {
                iend = f.Start - 1;
                istart = f.End - 1;
            }
            else
            {
                istart = f.Start - 1;
                iend = f.End - 1;
            }

            sb.Append(t.SeqId).Append('\t').Append(istart).Append('\t').Append(iend).Append('\t').Append(ti.Gff2Name ?? ti.Gff3Name).Append("\t0\t").Append(strand).Append('\t').Append(istart).Append('\t').Append(iend).Append("\t0");
            sb.Append("\t1\t0\t").Append(Int(f.End - f.Start + 1)).Append('\n');
        }
    }

    // ---------------------------------------------------------------------------------------------- draw ----

    private static void Draw(StringBuilder sb, EmbossReport t)
    {
        sb.Append("Start ").Append(Int(t.Start != 0 ? t.Start : 1)).Append('\n');
        sb.Append("End   ").Append(Int(t.End != 0 ? t.End : unchecked((int)t.TableLength))).Append("\n\ngroup\n");
        foreach (var f in t.Features)
        {
            var ti = TypeOf(f);
            sb.Append("label\nRange ").Append(Int(f.Start)).Append(' ').Append(Int(f.End)).Append(" | | 9 H\n").Append(t.Protein ? ti.Name : ti.NucName).Append("\nendlabel\n\n");
        }

        sb.Append("endgroup\n");
    }

    // ---------------------------------------------------------------------------------------------- debug ----

    // ajFeattablePrint
    private static void FeatureTablePrint(StringBuilder sb, EmbossReport t)
    {
        sb.Append("  Format: 0 (").Append(t.TableFormat).Append(")  Start: ").Append(Int(t.Start)).Append('\n');
        sb.Append("  End: ").Append(Int(t.End)).Append('\n');
        sb.Append("  Len: ").Append(Int(t.TableLength)).Append('\n');
        sb.Append("  Offset: ").Append(Int(t.Offset)).Append('\n');
        sb.Append("  Groups: ").Append(Int(t.Groups ?? t.Features.Count)).Append('\n');
        int n = 0;
        foreach (var f in t.Features)
        {
            n++;
            sb.Append("\n  Feature ").Append(n).Append('\n');
            sb.Append("    Source: '").Append(f.Source).Append("'\n");
            sb.Append("    Type: '").Append(TypeOf(f).Id).Append("'\n");
            sb.Append("    Score: ").Append(CFormat.Fixed((float)f.Score, 6)).Append('\n');
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
    }

    // --------------------------------------------------------------------------------------------- refseqp ----

    // feattableWriteRefseqp: the GenBank layout of a protein table (no complement; the type of the refseqp vocabulary is misc_feature for a type it does not know)
    private static void RefseqP(StringBuilder sb, EmbossReport t, EmbossFeatureWriteState state)
    {
        sb.Append("FEATURES             Location/Qualifiers\n");
        const string lead = "  ";
        string preftag = lead + "                   ";
        foreach (var f in t.Features)
        {
            string preftyp = lead + "   " + TypeOf(f).RefseqpName.PadRight(15) + " ";
            string loc = Int(f.Start) + (f.Between ? "^" : "..") + Int(f.End);
            LocWrap(loc, 79, preftag, preftyp, sb);
            foreach (var tag in f.Tags)
            {
                if (!tag.Name.Equals("note", StringComparison.OrdinalIgnoreCase))
                    throw new ArgumentException($"The feature tag '{tag.Name}' is not in the tag vocabulary of the feature output layer (note).");
                string line = "/note=" + QuoteEmbl(tag.Value) + "\n";
                TagWrap(line, 79, preftag, sb);
                state.Leftover = CollapseWhite(line);
            }
        }
    }
}
