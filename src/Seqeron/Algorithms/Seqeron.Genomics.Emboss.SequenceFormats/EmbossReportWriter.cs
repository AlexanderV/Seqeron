// SPDX-License-Identifier: LGPL-2.1-or-later
// Port of the EMBOSS 6.6.0 AJAX / NUCLEUS libraries (library source files carry the GNU Lesser General Public License, version 2.1 or later).
// Part of Seqeron.Genomics.Emboss.SequenceFormats, licensed under the GNU Lesser General Public License, version 2.1 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using System.Globalization;
using System.Text;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Emboss.SequenceFormats;

using Seqeron.Genomics.Alignment;

/// <summary>The date layouts of the EMBOSS 6.6 time module (<c>ajTimeRefTodayFmt</c>) shared by the alignment and report writers.</summary>
public static class EmbossDates
{
    private static readonly string[] Days = ["Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat"];
    private static readonly string[] Months = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"];

    /// <summary><c>ajTimeRefTodayFmt("report")</c>: <c>%a %e %b %Y %H:%M:%S</c>.</summary>
    public static string Report(DateTime d)
        => string.Create(CultureInfo.InvariantCulture, $"{Days[(int)d.DayOfWeek]} {d.Day,2} {Months[d.Month - 1]} {d.Year:D4} {d.Hour:D2}:{d.Minute:D2}:{d.Second:D2}");

    /// <summary><c>ajTimeRefTodayFmt("GFF")</c>: <c>%Y-%m-%d</c>.</summary>
    public static string Gff(DateTime d) => d.ToString("yyyy-MM-dd", CultureInfo.InvariantCulture);
}

/// <summary>The program refused the report (<c>ajErr</c> of <c>ajreport.c</c> / <c>ajfeatwrite.c</c>: an unsupported format for the sequence type, a format needing more tags than the program has).</summary>
public sealed class EmbossReportException : InvalidOperationException
{
    /// <summary>Creates the exception with the text the program prints after <c>Error:</c>.</summary>
    public EmbossReportException(string message) : base(message) { }

    /// <summary>The <c>-origfile</c> text the program had already written when the report was refused (etandem, equicktandem), else null.</summary>
    public string? OriginalFile { get; init; }
}

/// <summary>
/// The ONE writer of the EMBOSS feature-report files (<c>ajreport.c</c> <c>reportFormat[]</c>: every <c>-rformat</c> of the reports of fuzznuc, fuzzpro, fuzztran, etandem and equicktandem; <c>ajfeatwrite.c</c> for the feature-table
/// formats gff / gff3 / embl / genbank / pir / swiss / dasgff), checked against the installed EMBOSS 6.6.0 byte for byte on the files of the five programs (TestData/Alignment/emboss_report_formats_oracle.jsonl.gz,
/// scripts/alignment/generate_emboss_report_formats_oracle.py; audit A5-8b, B09 F184). One report file = a header (first sequence), per sequence header + body + tail, a footer; <see cref="EmbossReportOptions"/> carries the report options
/// <c>-raccshow -rdesshow -rscoreshow -rstrandshow -rusashow -rmaxall -rmaxseq</c>. The feature vocabulary (type names per format) covers the types these programs write.
/// </summary>
public static partial class EmbossReportWriter
{
    /// <summary>The end of <c>##sequence-region</c> of a table whose length is not set (fuzzpro, fuzztran): the largest coordinate read as an unsigned 32-bit integer (<c>ajFeattableGetLen</c>), so a table with negative reverse-frame coordinates names a number near 2^32.</summary>
    public static long UnsignedRegionEnd(IEnumerable<EmbossFeature> features)
    {
        long end = 0;
        foreach (var f in features)
            end = Math.Max(end, Math.Max(unchecked((uint)f.Start), unchecked((uint)f.End)));
        return end;
    }

    // ------------------------------------------------------------------ the single-sequence entry point of F176 ----

    /// <summary>Writes the report of ONE sequence in <paramref name="format"/> with the default options of the five programs (no command line, <paramref name="date"/> as the run date).</summary>
    /// <exception cref="ArgumentNullException"><paramref name="report"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="format"/> is not a defined <see cref="EmbossReportFormat"/>.</exception>
    /// <exception cref="ArgumentException">A field of <paramref name="report"/> or of one of its features is null (audit A6-4).</exception>
    /// <exception cref="EmbossReportException">The format is not written for the sequence type.</exception>
    public static string Write(EmbossReport report, EmbossReportFormat format, DateTime date, EmbossReportOptions? options = null)
    {
        ArgumentNullException.ThrowIfNull(report);
        if (!Enum.IsDefined(format))
            throw new ArgumentOutOfRangeException(nameof(format), format, "Unknown report format.");
        var o = (options ?? new EmbossReportOptions { Program = report.Features is { Count: > 0 } ? report.Features[0]?.Source ?? "" : "", Tags = report.Tags ?? [] }) with { Format = format, RunDate = date };
        Validate(o, [report], nameof(report));
        if (report.Features.Count == 0 && !report.WriteWhenEmpty && !o.WriteEmptyTables)         // fuzznuc, fuzzpro and fuzztran hand the report only a table with a hit
            return string.Empty;
        return WriteFile(o, [report]);
    }

    // ------------------------------------------------------------------------------------- the report file ----

    private sealed class State
    {
        public required EmbossReportOptions O;
        public StringBuilder Out = new();
        public int CountSeq, CountHit, TotHits;
        public bool MaxLimit;
        public bool ShowUsa;
        public string Type = "";
        public string SubTail = "";
        public int FormatIndex;
        public bool Opened;

        /// <summary>The feature tables written to this report file so far (<c>ftout-&gt;Count</c>: dasgff writes its XML header for the first one only).</summary>
        public int FeatureTables;
        public string Leftover = "";          // the static featoutStr of ajfeatwrite.c between two writers (feature output)
    }

    private static void Validate(EmbossReportOptions options, IReadOnlyList<EmbossReport> tables, string tablesParam = "tables")
    {
        if (options is null || options.Program is null || options.CommandLine is null || options.Tags is null || options.Header is null || options.ReportFile is null)
            throw new ArgumentException("The options need a Program, a CommandLine, Tags, a Header and a ReportFile.", nameof(options));
        if (options.CommandLine.Any(c => c is null) || options.Tags.Any(t => t is null || t.Type is null || t.Name is null || t.Print is null))
            throw new ArgumentException("The options' CommandLine or Tags contain null.", nameof(options));
        ArgumentNullException.ThrowIfNull(tables);
        for (int t = 0; t < tables.Count; t++)
        {
            var r = tables[t] ?? throw new ArgumentException($"Table {t} is null.", tablesParam);
            if (r.Sequence is null || r.Sequence.Name is null || r.Sequence.Usa is null || r.Sequence.Accession is null || r.Sequence.Description is null || r.Sequence.Residues is null || r.TableType is null)
                throw new ArgumentException($"The sequence of table {t} has a null field.", tablesParam);
            if (r.Features is null)
                throw new ArgumentException("The report's Features list is null.", tablesParam);
            for (int i = 0; i < r.Features.Count; i++)
            {
                var f = r.Features[i];
                if (f is null)
                    throw new ArgumentException($"Feature {i} is null.", tablesParam);
                if (f.Type is null || f.Source is null || f.Tags is null || f.Tags.Any(g => g is null || g.Name is null || g.Value is null))
                    throw new ArgumentException($"Feature {i} has a null Type / Source / Tags or a null tag field.", tablesParam);
            }
        }
    }

    /// <summary>
    /// Writes the whole report file the program writes: <c>ajReportWrite</c> for each table in turn (the program passes the tables that have a feature, etandem / equicktandem every table), then <c>ajReportClose</c>
    /// (<see cref="EmbossReportOptions.WriteFooter"/>); <paramref name="totals"/> names the sequences read and their total length for the footer once it is known which table stopped the run (-1 = none). Tables after the <c>-rmaxall</c> limit are not written (the program stops reading sequences at the first one it cannot write).
    /// </summary>
    /// <exception cref="EmbossReportException">A format that needs more tags than the program has (<c>diffseq</c>), or a format not written for the sequence type (embl / genbank for proteins, pir / swiss for nucleotides): the program writes an empty file.</exception>
    public static string WriteFile(EmbossReportOptions options, IReadOnlyList<EmbossReport> tables, Func<int, (long Sequences, long Length)>? totals = null)
    {
        ArgumentNullException.ThrowIfNull(options);
        Validate(options, tables);
        if (!Enum.IsDefined(options.Format))
            throw new ArgumentOutOfRangeException(nameof(options), options.Format, "Unknown report format.");
        var st = new State { O = options, FormatIndex = EmbossReportFormats.IndexOf(options.Format), Leftover = options.FeatureLeftover };
        var flags = EmbossReportFormats.Flags(options.Format);
        if (flags.MinTags > options.Tags.Count)
            throw new EmbossReportException($"Report format '{FormatName(st)}' needs {flags.MinTags} tags, has only {options.Tags.Count}");
        int stoppedAt = -1;
        for (int i = 0; i < tables.Count; i++)
        {
            if (!WriteTable(st, tables[i]))
            {
                stoppedAt = i;
                break;
            }
        }

        if (st.O.Format == EmbossReportFormat.DasGff && st.Opened)
            st.Out.Append("  </GFF>\n</DASGFF>\n");
        if (totals is not null)
        {
            var (sequences, length) = totals(stoppedAt);
            st.O = st.O with { TotalSequences = sequences, TotalLength = length };
        }

        if (flags.Header && st.O.WriteFooter)
        {
            if (st.CountSeq == 0)
                Header(st, null, null);
            Tail(st, null, null);
        }

        return st.Out.ToString();
    }

    private static string FormatName(State st) => st.O.FormatText ?? EmbossReportFormats.Name(st.O.Format);

    // ajReportWrite
    private static bool WriteTable(State st, EmbossReport t)
    {
        var o = st.O;
        var flags = EmbossReportFormats.Flags(o.Format);
        bool tableProt = t.Protein || t.TypeUnknown, tableNuc = !t.Protein || t.TypeUnknown;
        bool isnuc = t.Sequence.Nucleotide && tableNuc, isprot = !t.Sequence.Nucleotide && tableProt;
        bool ok = (isnuc && flags.Nuc) || (isprot && flags.Prot);
        if (!ok)
        {
            if (!isnuc && !isprot)
                throw new EmbossReportException($"Report format '{FormatName(st)}' failed: unknown sequence type");
            throw new EmbossReportException($"Report format '{FormatName(st)}' not supported for {(isnuc ? "nucleotide" : "protein")} sequences");
        }

        int maxreport = 0;
        if (o.MaxAll != 0)
        {
            maxreport = o.MaxAll - st.CountHit;
            if (maxreport <= 0)
            {
                st.MaxLimit = true;
                return false;
            }
        }

        if (o.MaxSeq != 0)
            maxreport = o.MaxAll != 0 ? Math.Min(maxreport, o.MaxSeq) : o.MaxSeq;
        if (st.Type.Length == 0)                                       // ajReportSetType: the first table's type letter
            st.Type = t.TableType.Length > 0 ? t.TableType : DefaultType(t);
        var features = t.Features;
        if (maxreport != 0 && maxreport < features.Count)
        {
            st.SubTail = AppendTail(st.SubTail, $"HitLimit: {maxreport}/{features.Count}");
            var limited = t with { Features = [.. features.Take(maxreport)], Groups = t.Groups ?? features.Count, TableFormat = "", FeatureTable = t.FeatureTable?.NewFtableLimit(maxreport) };
            Body(st, limited);
            st.CountHit += maxreport;
            st.TotHits += features.Count;
        }
        else
        {
            Body(st, t);
            st.CountHit += features.Count;
            st.TotHits += features.Count;
        }

        st.CountSeq++;
        return true;
    }

    private static string DefaultType(EmbossReport t)
    {
        if (t.Sequence.Nucleotide)
            return "N";
        return "P";
    }

    // ajReportAppendSubtailS: a second line goes below the first
    private static string AppendTail(string old, string text) => old.Length == 0 ? text : old + "\n" + text;

    // --------------------------------------------------------------------------------------------- helpers ----

    private static string Int(long v) => v.ToString(CultureInfo.InvariantCulture);

    // ajFmt "%*S" / "%-*S": a negative width left-justifies
    private static string Pad(string s, int width, bool left)
    {
        if (width < 0)
        {
            width = -width;
            left = true;
        }

        return left ? s.PadRight(width) : s.PadLeft(width);
    }

    // ajStrAssignSubS: positions 0-based inclusive, negative counted from the end, both clipped
    private static string SubS(string s, long pos1, long pos2)
    {
        long len = s.Length;
        if (len == 0)
            return "";
        long b = pos1 < 0 ? len + pos1 : pos1;
        if (b >= len)
            b = len - 1;
        if (b < 0)
            b = 0;
        long e = pos2 < 0 ? len + pos2 : pos2;
        if (e >= len)
            e = len - 1;
        if (e < b)
            e = b;
        return s.Substring((int)b, (int)(e - b + 1));
    }

    private static string RevComp(string s) => EmbossInputSequence.ReverseComplement(s);

    // the sequence text the writers use: ajSeqstrReverse of it when the sequence is reversed
    private static string SeqAlias(EmbossReportSequence s) => s.Reversed ? RevComp(s.Residues) : s.Residues;


    private static char FeatStrand(EmbossFeature f, bool protein)
    {
        if (f.Protein ?? protein)
            return '.';
        return f.Strand == '-' ? '-' : '+';
    }

    // the feature types of the programs: internal SO id, the EMBOSS name (ajFeatGetType / ajFeatTypeNuc / ajFeatTypeProt), the GFF3 name (also the DAS one), the DAS category
    private sealed record TypeInfo(string Id, string Name, string NucName, string Gff3Name, string DasCategory, string? EmblName, string PirType, string SwissType, string? Gff2Name = null, string RefseqpName = "misc_feature");

    private static readonly Dictionary<string, TypeInfo> Types = new(StringComparer.Ordinal)
    {
        ["SO:0000714"] = new("SO:0000714", "nucleotide_motif", "nucleotide_motif", "nucleotide_motif", "other", null, "", ""),
        ["SO:0001067"] = new("SO:0001067", "polypeptide_motif", "misc_feature", "polypeptide_motif", "other", null, "Region", "REGION"),
        ["SO:0000657"] = new("SO:0000657", "repeat_region", "repeat_region", "repeat_region", "repeat", "repeat_region", "", ""),
        ["SO:0000110"] = new("SO:0000110", "misc_feature", "misc_feature", "sequence_feature", "other", "misc_feature", "", "", "located_sequence_feature"),
        ["SO:0005836"] = new("SO:0005836", "misc_signal", "misc_signal", "regulatory_region", "transcription", "misc_signal", "", ""),                                     // marscan's MRS: the feature key of EMBL, regulatory_region in GFF3
        ["SO:0000839"] = new("SO:0000839", "polypeptide_region", "misc_feature", "sequence_feature", "other", null, "Region", "REGION"),          // diffseq's generic feature in a protein table
        ["polypeptide_region"] = new("SO:0000839", "polypeptide_region", "misc_feature", "polypeptide_region", "other", null, "Region", "REGION", null, "Region"),   // wordmatch's misc_feature on a protein table
        ["conflict"] = new("SO:0001085", "sequence_conflict", "conflict", "sequence_conflict", "other", "conflict", "Region", "CONFLICT"),         // diffseq's feature files (-aoutfeat / -boutfeat)
    };

    private static TypeInfo TypeOf(EmbossFeature f)
    {
        if (!f.FromTable && Types.TryGetValue(f.Type, out var t))
            return t;
        if (!f.FromTable)
        {
            foreach (var v in Types.Values)
                if (string.Equals(v.Name, f.Type, StringComparison.Ordinal))
                    return v;
        }

        // a type of the feature library (the tables of featreport and twofeat): its name in the vocabulary of the sequence kind (ajFeatGetType)
        var vocab = f.Protein == true ? EmbossFeatureVocab.Protein : EmbossFeatureVocab.Dna;
        string name = EmbossFeatureVocab.TypeExternal(vocab.Types, f.Type);
        string nucName = EmbossFeatureVocab.TypeExternal(EmbossFeatureVocab.Dna.Types, f.Type);          // ajFeatTypeNuc
        return new TypeInfo(f.Type, name, nucName, name, "other", null, "", "");
    }

    private static string TypeName(EmbossFeature f) => TypeOf(f).Name;

    // the regions and draw reports name diffseq's polypeptide_region by its feature-table name misc_feature (the tagseq report keeps the internal name)
    private static string RegionTypeName(EmbossFeature f) => !f.FromTable && TypeOf(f).Id == "SO:0000839" ? "misc_feature" : TypeName(f);

    // ajFeatGetNoteS / ajFeatGetNoteCI: the first note "*name value" whose text after '*' starts (case-insensitively) with name; the value after "*name " ("" when none); false when that note has another character after the name
    private static bool GetNote(EmbossFeature f, string name, out string value, int count = 0)
    {
        value = "";
        int n = 0;
        foreach (var tag in f.Tags)
        {
            if (!string.Equals(tag.Name, "note", StringComparison.OrdinalIgnoreCase))
                continue;
            string v = tag.Value;
            if (v.Length == 0 || v[0] != '*')
                continue;
            if (v.Length - 1 >= name.Length && string.Compare(v, 1, name, 0, name.Length, StringComparison.OrdinalIgnoreCase) == 0)
            {
                n++;
                if (n >= count)
                {
                    if (v.Length > name.Length + 1)
                    {
                        if (v[name.Length + 1] != ' ')
                            return false;
                        value = v[(name.Length + 2)..];
                    }

                    return true;
                }
            }
        }

        return false;
    }

    // ------------------------------------------------------------------------ ajReportWriteHeader / Tail ----

    private static string ExchangeBlock(string text)
    {
        string t = text.Replace("\n", "\n# ", StringComparison.Ordinal).TrimEnd(' ');
        return t;
    }

    private static void Header(State st, EmbossReport? t, EmbossReportSequence? seq)
    {
        var o = st.O;
        var sb = st.Out;
        if (st.CountSeq == 0)
        {
            sb.Append("########################################\n");
            sb.Append("# Program: ").Append(o.Program).Append('\n');
            sb.Append("# Rundate: ").Append(EmbossDates.Report(o.RunDate)).Append('\n');
            sb.Append("# Commandline: ").Append(o.Program).Append('\n');
            if (o.CommandLine.Count > 0)
                sb.Append("#    ").Append(string.Join("\n#    ", o.CommandLine)).Append('\n');
            if (o.MaxAll != 0)
                sb.Append("# Max_hits_total: ").Append(o.MaxAll).Append('\n');
            if (o.MaxSeq != 0)
                sb.Append("# Max_hits_sequence: ").Append(o.MaxSeq).Append('\n');
            sb.Append("# Report_format: ").Append(FormatName(st)).Append('\n');
            sb.Append("# Report_file: ").Append(o.ReportFile).Append('\n');
            if (o.AdditionalFiles.Count > 0)
            {
                sb.Append("# Additional_files: ").Append(o.AdditionalFiles.Count).Append('\n');
                for (int i = 0; i < o.AdditionalFiles.Count; i++)
                    sb.Append("# ").Append(i + 1).Append(": ").Append(o.AdditionalFiles[i].Name).Append(" (").Append(o.AdditionalFiles[i].Description).Append(")\n");
            }

            sb.Append("########################################\n\n");
        }

        sb.Append("#=======================================\n#\n");
        if (seq is not null)
        {
            string name = st.ShowUsa ? seq.Usa : seq.Name;
            if (seq.Reversed)
                sb.Append("# Sequence: ").Append(name).Append("     from: ").Append(Int(seq.Begin + seq.Offset)).Append("   to: ").Append(Int(seq.End + seq.Offset)).Append("   (Reversed)\n");
            else
                sb.Append("# Sequence: ").Append(name).Append("     from: ").Append(Int(seq.Begin + seq.Offset)).Append("   to: ").Append(Int(seq.End + seq.Offset)).Append('\n');
            if (o.ShowAccession)
                sb.Append("# Accession: ").Append(seq.Accession).Append('\n');
            if (o.ShowDescription)
                sb.Append("# Description: ").Append(seq.Description).Append('\n');
        }

        if (t is not null)
            sb.Append("# HitCount: ").Append(t.Features.Count).Append('\n');
        if (o.Header.Length > 0)
            AppendBlock(sb, o.Header);
        sb.Append("#=======================================\n\n");
    }

    private static void AppendBlock(StringBuilder sb, string text)
    {
        string tmp = ExchangeBlock(text);
        sb.Append("#\n");
        sb.Append("# ").Append(tmp);
        if (!tmp.EndsWith("\n#", StringComparison.Ordinal))
            sb.Append("\n#");
        sb.Append('\n');
    }

    private static void Tail(State st, EmbossReport? t, EmbossReportSequence? seq)
    {
        var sb = st.Out;
        sb.Append("\n#---------------------------------------\n");
        if (t is not null)
        {
            if (st.SubTail.Length > 0)
            {
                AppendBlock(sb, st.SubTail);
                st.SubTail = "";
            }

            if (st.O.Tail.Length > 0)
                AppendBlock(sb, st.O.Tail);
        }
        else
        {
            var o = st.O;
            if (o.TotalSequences != 0)
                sb.Append("# Total_sequences: ").Append(o.TotalSequences).Append('\n');
            if (o.TotalLength != 0)
                sb.Append("# Total_length: ").Append(o.TotalLength).Append('\n');
            sb.Append("# Reported_sequences: ").Append(st.CountSeq).Append('\n');
            sb.Append("# Reported_hitcount: ").Append(st.CountHit).Append('\n');
            if (st.CountHit < st.TotHits)
                sb.Append("# Unreported_hitcount: ").Append(st.TotHits - st.CountHit).Append('\n');
            if (o.MaxAll != 0)
                sb.Append("# Max_hitcount: ").Append(o.MaxAll).Append('\n');
            if (st.MaxLimit)
                sb.Append("# Maxhits_stop: Yes\n");
        }

        sb.Append("#---------------------------------------\n");
    }

    // reportWriteDiffseq: the differences between two sequences (the program's seven tags start, end, length, name, sequence, first_feature, second_feature are notes "*name value" of the feature)
    private static void DiffSeq(State st, EmbossReport t)
    {
        var sb = st.Out;
        string seqalias = SeqAlias(t.Sequence);
        string name = st.ShowUsa ? t.Sequence.Usa : t.Sequence.Name;
        Header(st, t, t.Sequence);
        foreach (var f in t.Features)
        {
            char strand = FeatStrand(f, t.Protein);
            long istart = f.Start, iend = f.End;
            long ilen = iend - istart + 1;
            if (ilen == 2 && f.Between)
                ilen = 0;
            // reportWriteDiffseq keeps the feature's start and end in ajuint: a start before the sequence (-1) is the huge positive number the string cut then clips to the last character
            string subseq = SubS(seqalias, unchecked((uint)(f.Start - 1)), unchecked((uint)(f.End - 1)));
            if (strand == '-')
                subseq = RevComp(subseq);
            uint jstart = GetNote(f, "start", out var v) ? NoteUInt(v) : 0;
            uint jend = GetNote(f, "end", out v) ? NoteUInt(v) : 0;
            uint jlen = GetNote(f, "length", out v) ? NoteUInt(v) : 0;           // printed with %d: a wrapped negative number reads as negative again
            char jstrand = '+';
            if (GetNote(f, "strand", out v))
                jstrand = v.Length > 0 ? v[0] : '\0';
            if (!GetNote(f, "name", out string jname))
                jname = "";
            if (ilen > 0)
            {
                sb.Append('\n').Append(name).Append(' ').Append(Int(istart)).Append('-').Append(Int(iend)).Append(" Length: ").Append(Int(ilen)).Append(strand == '-' ? " (Reversed)\n" : "\n");
                DiffSeqFeatures(sb, f, "first_feature");
                sb.Append("Sequence: ").Append(subseq).Append('\n');
            }
            else
            {
                sb.Append('\n').Append(name).Append(' ').Append(Int(istart)).Append(" Length: ").Append(Int(ilen)).Append('\n');
                DiffSeqFeatures(sb, f, "first_feature");
                sb.Append("Sequence: \n");
            }

            if (!GetNote(f, "sequence", out subseq))
                subseq = "";
            if (jlen > 0)
            {
                sb.Append("Sequence: ").Append(subseq).Append('\n');
                DiffSeqFeatures(sb, f, "second_feature");
                sb.Append(jname).Append(' ').Append(Int(unchecked((int)jstart))).Append('-').Append(Int(unchecked((int)jend))).Append(" Length: ").Append(Int(unchecked((int)jlen))).Append(jstrand == '-' ? " (Reversed)\n" : "\n");
            }
            else
            {
                sb.Append("Sequence: \n");
                DiffSeqFeatures(sb, f, "second_feature");
                sb.Append(jname).Append(' ').Append(Int(unchecked((int)jstart))).Append(" Length: ").Append(Int(unchecked((int)jlen))).Append(jstrand == '-' ? " (Reversed)\n" : "\n");
            }
        }

        Tail(st, t, t.Sequence);
    }

    // while(ajFeatGetNoteCI(feature, tag, n++, &text)) "Feature: %S"
    private static void DiffSeqFeatures(StringBuilder sb, EmbossFeature f, string tag)
    {
        for (int n = 1; GetNote(f, tag, out string text, n); n++)
            sb.Append("Feature: ").Append(text).Append('\n');
    }

    // ajStrToUint: the leading unsigned decimal number (0 when none)
    private static uint NoteUInt(string v)
    {
        string t = v.TrimStart();
        bool negative = t.StartsWith('-');
        if (negative)
            t = t[1..];
        int e = 0;
        while (e < t.Length && char.IsAsciiDigit(t[e]))
            e++;
        if (e == 0 || !uint.TryParse(t.AsSpan(0, e), NumberStyles.None, CultureInfo.InvariantCulture, out uint r))
            return 0;
        return negative ? unchecked(0u - r) : r;
    }

    private static string CharName(State st) => st.Type.Length > 0 && (st.Type[0] == 'n' || st.Type[0] == 'N') ? "bases" : "residues";

    private static string Score(State st, EmbossFeature f) => CFormat.Fixed((float)f.Score, st.O.Precision);

    private sealed record Tag(string Type, string Name, string Print, int Size);

    private static List<Tag> GetTags(State st) => [.. st.O.Tags.Select(t => new Tag(t.Type, t.Name, t.Print, Math.Max(6, t.Print.Length)))];

    // -------------------------------------------------------------------------------------- the bodies ----

    private static void Body(State st, EmbossReport t)
    {
        st.ShowUsa = st.ShowUsa || st.O.ShowUsa;
        switch (st.O.Format)
        {
            case EmbossReportFormat.Gff:
            case EmbossReportFormat.Gff3:
                if (t.FeatureTable is not null)
                    TableBody(st, t, "gff", "ajReportWriteGff features output failed format: 'gff'");
                else
                    Gff(st, t);
                break;
            case EmbossReportFormat.Excel:
                Excel(st, t);
                break;
            case EmbossReportFormat.Embl:
                if (t.FeatureTable is not null)
                    TableBody(st, t, "embl", "ajReportWriteEmbl features output failed format: 'embl'");
                else
                    Embl(st, t, true);
                break;
            case EmbossReportFormat.Genbank:
                if (t.FeatureTable is not null)
                    TableBody(st, t, "genbank", "ajReportWriteGenbank features output failed format: 'genbank'");
                else
                    Embl(st, t, false);
                break;
            case EmbossReportFormat.Pir:
                if (t.FeatureTable is not null)
                    TableBody(st, t, "pir", "ajReportWritePir features output failed format: 'pir'");
                else
                    Pir(st, t);
                break;
            case EmbossReportFormat.Swiss:
                if (t.FeatureTable is not null)
                    TableBody(st, t, "swiss", "ajReportWriteSwiss features output failed format: 'swiss'");
                else
                    Swiss(st, t);
                break;
            case EmbossReportFormat.DasGff:
                if (t.FeatureTable is not null)
                {
                    st.Opened = true;                                      // the end of the XML is written once, when the report is closed
                    TableBody(st, t, "dasgff", "ajReportWriteDasgff features output failed format: 'dasgff'");
                }
                else
                    DasGff(st, t);
                break;
            case EmbossReportFormat.Debug:
                DebugReport(st, t);
                break;
            case EmbossReportFormat.ListFile:
                st.ShowUsa = true;                                         // "so we get a usable USA"
                ListFile(st, t);
                break;
            case EmbossReportFormat.DbMotif:
                DbMotif(st, t);
                break;
            case EmbossReportFormat.Draw:
                Draw(st, t);
                break;
            case EmbossReportFormat.Restrict:
                Restrict(st, t);
                break;
            case EmbossReportFormat.FeatTable:
                FeatTable(st, t);
                break;
            case EmbossReportFormat.DiffSeq:
                DiffSeq(st, t);
                break;
            case EmbossReportFormat.Motif:
                Motif(st, t);
                break;
            case EmbossReportFormat.NameTable:
                TableLike(st, t, true, false);
                break;
            case EmbossReportFormat.Regions:
                Regions(st, t);
                break;
            case EmbossReportFormat.SeqTable:
                SeqTable(st, t);
                break;
            case EmbossReportFormat.Simple:
                SrsFlags(st, t, false);
                break;
            case EmbossReportFormat.Srs:
                SrsFlags(st, t, true);
                break;
            case EmbossReportFormat.Summary:
                if (t.Features.Count > 0)
                    Header(st, t, t.Sequence);
                break;
            case EmbossReportFormat.Table:
                TableLike(st, t, false, false);
                break;
            case EmbossReportFormat.TagSeq:
                TagSeq(st, t);
                break;
            default:
                throw new EmbossReportException($"Report format '{FormatName(st)}' is not written for these programs.");
        }
    }

    // reportWriteGff / Embl / Genbank / Pir / Swiss / Dasgff of a full feature table: ajFeattabOutNewCSF(format, ...) and ajFeattableWrite, a failure is a warning
    private static void TableBody(State st, EmbossReport t, string format, string failure)
    {
        var context = st.O.Context ?? new EmbossFeatureReadContext();
        string? text = WriteFeatureTable(t.FeatureTable!, format, st.O.RunDate, st.O.Program, context, new EmbossFeatureWriteState(), st.FeatureTables);
        if (text is null)
        {
            context.Warn(failure);
        }
        else
        {
            st.Out.Append(text);
            st.FeatureTables++;
        }
    }

    // reportWriteExcel (no header, no tail)
    private static void Excel(State st, EmbossReport t)
    {
        var sb = st.Out;
        var tags = GetTags(st);
        sb.Append(st.O.ShowScore ? "SeqName\tStart\tEnd\tScore\tStrand" : "SeqName\tStart\tEnd\tStrand");
        foreach (var g in tags)
            sb.Append('\t').Append(g.Print);
        sb.Append('\n');
        foreach (var f in t.Features)
        {
            char strand = f.Strand == '-' ? '-' : '+';
            string name = st.ShowUsa ? t.Sequence.Usa : t.Sequence.Name;
            if (st.O.ShowScore)
                sb.Append(name).Append('\t').Append(Int(f.Start)).Append('\t').Append(Int(f.End)).Append('\t').Append(Score(st, f)).Append('\t').Append(strand);
            else
                sb.Append(name).Append('\t').Append(Int(f.Start)).Append('\t').Append(Int(f.End)).Append('\t').Append(strand);
            foreach (var g in tags)
                sb.Append('\t').Append(GetNote(f, g.Name, out var v) ? v : ".");
            sb.Append('\n');
        }
    }

    // reportWriteSrsFlags (simple / srs)
    private static void SrsFlags(State st, EmbossReport t, bool withSeq)
    {
        var sb = st.Out;
        Header(st, t, t.Sequence);
        var tags = GetTags(st);
        string seq = t.Sequence.Residues;
        int ift = 0;
        foreach (var f in t.Features)
        {
            uint istart = unchecked((uint)f.Start), iend = unchecked((uint)f.End);
            int ilen = unchecked((int)(iend - istart + 1));
            string subseq = SubS(seq, unchecked((long)(istart - 1)), unchecked((long)(iend - 1)));
            if (ift > 0)
                sb.Append('\n');
            ift++;
            sb.Append("Feature: ").Append(ift).Append('\n');
            sb.Append("Name: ").Append(st.ShowUsa ? t.Sequence.Usa : t.Sequence.Name).Append('\n');
            sb.Append("Start: ").Append(Int(unchecked((int)istart))).Append('\n');
            sb.Append("End: ").Append(Int(unchecked((int)iend))).Append('\n');
            sb.Append("Length: ").Append(Int(ilen)).Append('\n');
            if (withSeq)
                sb.Append("Sequence: ").Append(subseq).Append('\n');
            sb.Append("Score: ").Append(Score(st, f)).Append('\n');
            sb.Append("Strand: ").Append(FeatStrand(f, t.Protein)).Append('\n');
            if (f.Frame != 0)
                sb.Append("Frame: ").Append(f.Frame).Append('\n');
            foreach (var g in tags)
                sb.Append(g.Print).Append(": ").Append(GetNote(f, g.Name, out var v) ? v : ".").Append('\n');
        }

        Tail(st, t, t.Sequence);
    }

    // reportWriteTable / reportWriteNameTable
    private static void TableLike(State st, EmbossReport t, bool nameColumn, bool unused)
    {
        _ = unused;
        var sb = st.Out;
        bool dostrand = st.O.ShowStrand && t.Sequence.Nucleotide;
        Header(st, t, t.Sequence);
        var tags = GetTags(st);
        string seq = t.Sequence.Residues;
        foreach (var f in t.Features)
            for (int j = 0; j < tags.Count; j++)
                if (GetNote(f, tags[j].Name, out var v))
                    tags[j] = tags[j] with { Size = Math.Max(tags[j].Size, v.Length) };
        if (nameColumn)
            sb.Append(Pad("USA", 20, true)).Append(' ').Append("Start".PadLeft(7)).Append(' ').Append("End".PadLeft(7));
        else
            sb.Append("Start".PadLeft(7)).Append(' ').Append("End".PadLeft(7));
        if (dostrand)
            sb.Append(' ').Append("Strand".PadLeft(7));
        if (st.O.ShowScore)
            sb.Append(' ').Append("Score".PadLeft(7));
        foreach (var g in tags)
            sb.Append(' ').Append(Pad(g.Print, g.Size, g.Type == "str"));
        sb.Append('\n');
        foreach (var f in t.Features)
        {
            int istart = f.Start, iend = f.End;
            if (f.Between)                                      // reportWriteTable / NameTable: a between-feature prints its end first
                (istart, iend) = (iend, istart);
            if (nameColumn)
                sb.Append(Pad(st.ShowUsa ? t.Sequence.Usa : t.Sequence.Name, 20, true)).Append(' ').Append(Int(istart).PadLeft(7)).Append(' ').Append(Int(iend).PadLeft(7));
            else
                sb.Append(Int(istart).PadLeft(7)).Append(' ').Append(Int(iend).PadLeft(7));
            if (dostrand)
                sb.Append("       ").Append(FeatStrand(f, t.Protein));
            if (st.O.ShowScore)
                sb.Append(' ').Append(Score(st, f).PadLeft(7));
            foreach (var g in tags)
            {
                string v = GetNote(f, g.Name, out var val) ? val : ".";
                sb.Append(' ').Append(Pad(v, g.Size, g.Type == "str"));
            }

            sb.Append('\n');
        }

        _ = seq;
        Tail(st, t, t.Sequence);
    }

    // reportWriteSeqTable
    private static void SeqTable(State st, EmbossReport t)
    {
        var sb = st.Out;
        string seqalias = SeqAlias(t.Sequence);
        bool dostrand = st.O.ShowStrand && t.Sequence.Nucleotide;
        Header(st, t, t.Sequence);
        var tags = GetTags(st);
        foreach (var f in t.Features)
            for (int j = 0; j < tags.Count; j++)
                if (GetNote(f, tags[j].Name, out var v))
                    tags[j] = tags[j] with { Size = Math.Max(tags[j].Size, v.Length) };
        sb.Append("Start".PadLeft(7)).Append(' ').Append("End".PadLeft(7));
        if (dostrand)
            sb.Append(' ').Append("Strand".PadLeft(7));
        foreach (var g in tags)
            sb.Append(' ').Append(Pad(g.Print, g.Size, g.Type == "str"));
        sb.Append(" Sequence\n");
        foreach (var f in t.Features)
        {
            uint istart = unchecked((uint)f.Start), iend = unchecked((uint)f.End);
            if (f.Between)                                      // reportWriteSeqTable: a between-feature prints its end first
                (istart, iend) = (iend, istart);
            string subseq = f.Between ? "" : SubS(seqalias, unchecked((long)(istart - 1)), unchecked((long)(iend - 1)));
            if (f.Strand == '-')
                subseq = RevComp(subseq);
            sb.Append(Int(unchecked((int)istart)).PadLeft(7)).Append(' ').Append(Int(unchecked((int)iend)).PadLeft(7));
            if (dostrand)
                sb.Append("       ").Append(FeatStrand(f, t.Protein));
            foreach (var g in tags)
            {
                string v = GetNote(f, g.Name, out var val) ? val : ".";
                sb.Append(' ').Append(Pad(v, g.Size, g.Type == "str"));
            }

            sb.Append(' ').Append(subseq).Append('\n');
        }

        Tail(st, t, t.Sequence);
    }

    // reportWriteRegions
    private static void Regions(State st, EmbossReport t)
    {
        var sb = st.Out;
        Header(st, t, t.Sequence);
        var tags = GetTags(st);
        foreach (var f in t.Features)
        {
            uint istart = unchecked((uint)f.Start), iend = unchecked((uint)f.End);
            int ilen = unchecked((int)(iend - istart + 1));
            sb.Append(RegionTypeName(f)).Append(" from ").Append(Int(unchecked((int)istart))).Append(" to ").Append(Int(unchecked((int)iend))).Append(" (").Append(Int(ilen)).Append(' ').Append(CharName(st)).Append(")\n");
            if (st.O.ShowScore)
                sb.Append("   Max score: ").Append(Score(st, f));
            if (tags.Count > 0)
            {
                sb.Append(" (");
                for (int j = 0; j < tags.Count; j++)
                    if (GetNote(f, tags[j].Name, out var v))
                    {
                        if (j > 0)
                            sb.Append(", ");
                        sb.Append(tags[j].Print).Append(": ").Append(v);
                    }

                sb.Append(')');
            }

            sb.Append("\n\n");
        }

        Tail(st, t, t.Sequence);
    }

    // reportWriteListFile
    private static void ListFile(State st, EmbossReport t)
    {
        var sb = st.Out;
        Header(st, t, t.Sequence);
        foreach (var f in t.Features)
        {
            uint istart = unchecked((uint)f.Start), iend = unchecked((uint)f.End);
            var range = new StringBuilder("[");
            if (istart != 0)
                range.Append(Int(unchecked((int)istart)));
            range.Append(':');
            if (iend != 0)
                range.Append(Int(unchecked((int)iend)));
            if (f.Strand == '-')
                range.Append(":r");
            range.Append(']');
            string name = t.Sequence.Usa;
            if (name.Length > 0 && name[^1] == ']')
            {
                int ipos = name.LastIndexOf('[');
                name = name.Remove(ipos);
            }

            sb.Append(name);
            if (range.Length > 3)
                sb.Append(range);
            sb.Append('\n');
        }

        Tail(st, t, t.Sequence);
    }

    // reportWriteDraw
    private static void Draw(State st, EmbossReport t)
    {
        var sb = st.Out;
        Header(st, t, t.Sequence);
        sb.Append("Start ").Append(Int(t.Sequence.Begin + t.Sequence.Offset)).Append('\n');
        sb.Append("End   ").Append(Int(t.Sequence.End + t.Sequence.Offset)).Append('\n');
        sb.Append('\n');
        sb.Append("group\n");
        foreach (var f in t.Features)
        {
            sb.Append("label\n");
            sb.Append("Range ").Append(Int(f.Start)).Append(' ').Append(Int(f.End)).Append(" | | 9 H\n");
            sb.Append(RegionTypeName(f)).Append('\n');
            sb.Append("endlabel\n");
            sb.Append('\n');
        }

        sb.Append("endgroup\n");
        Tail(st, t, t.Sequence);
    }

    // reportWriteRestrict: without the 5prime / 3prime / enzyme tags of the restriction programs the groups are empty lines
    private static void Restrict(State st, EmbossReport t)
    {
        var sb = st.Out;
        Header(st, t, t.Sequence);
        var tags = GetTags(st);
        int jenz = tags.FindIndex(g => string.Equals(g.Name, "enzyme", StringComparison.OrdinalIgnoreCase));
        int j5 = tags.FindIndex(g => string.Equals(g.Name, "5prime", StringComparison.OrdinalIgnoreCase));
        int j3 = tags.FindIndex(g => string.Equals(g.Name, "3prime", StringComparison.OrdinalIgnoreCase));
        sb.Append("Start ").Append(Int(t.Sequence.Begin + t.Sequence.Offset)).Append('\n');
        sb.Append("End   ").Append(Int(t.Sequence.End + t.Sequence.Offset)).Append("\n\n");
        sb.Append("group\nForward\n");
        foreach (var f in t.Features)
        {
            if (j5 >= 0)
            {
                sb.Append("label\n");
                GetNote(f, tags[j5].Name, out var v);
                sb.Append("Tick ").Append(Int(ToUint(v))).Append(" 8\n");
                if (jenz >= 0)
                {
                    GetNote(f, tags[jenz].Name, out var e);
                    sb.Append(e).Append('\n');
                }
                else
                    sb.Append("Enz\n");
                sb.Append("endlabel\n");
            }

            sb.Append('\n');
        }

        sb.Append("endgroup\ngroup\nReverse\n");
        foreach (var f in t.Features)
        {
            if (j3 >= 0)
            {
                sb.Append("label\n");
                GetNote(f, tags[j3].Name, out var v);
                sb.Append("Tick ").Append(Int(ToUint(v))).Append(" 3\n");
                if (jenz >= 0)
                {
                    GetNote(f, tags[jenz].Name, out var e);
                    sb.Append(e).Append('\n');
                }
                else
                    sb.Append("Enz\n");
                sb.Append("endlabel\n");
            }

            sb.Append('\n');
        }

        sb.Append("endgroup\n");
        Tail(st, t, t.Sequence);
    }

    // ajStrToUint: leading decimal digits (0 when not a number)
    private static uint ToUint(string s)
        => uint.TryParse(s.Trim(), NumberStyles.None, CultureInfo.InvariantCulture, out var v) ? v : 0;

    // reportWriteFeatTable
    private static void FeatTable(State st, EmbossReport t)
    {
        var sb = st.Out;
        var tags = GetTags(st);
        foreach (var f in t.Features)
        {
            string type = RegionTypeName(f);
            string range = f.Between ? $"{Int(f.Start)}^{Int(f.End)}" : $"{Int(f.Start)}..{Int(f.End)}";
            string loc = f.Strand == '-' ? $"complement({range})" : range;
            sb.Append("FT   ").Append((type.Length > 15 ? type[..15] : type).PadRight(15)).Append(' ').Append(loc).Append('\n');
            foreach (var g in tags)
                if (GetNote(f, g.Name, out var v))
                    sb.Append("FT                   /").Append(g.Print).Append("=\"").Append(v).Append("\"\n");
        }
    }

    // reportWriteDbMotif
    private static void DbMotif(State st, EmbossReport t)
    {
        var sb = st.Out;
        string seqalias = SeqAlias(t.Sequence);
        Header(st, t, t.Sequence);
        var tags = GetTags(st);
        uint seqLen = (uint)t.Sequence.Length;
        foreach (var f in t.Features)
        {
            char strand = f.Strand == '-' ? '-' : '+';
            uint istart = unchecked((uint)f.Start), iend = unchecked((uint)f.End);
            int ilen = unchecked((int)(iend - istart + 1));
            uint jstart = Math.Max(6u, istart);
            jstart -= 6;
            uint jend = Math.Min(seqLen, unchecked(iend + 4));
            string spaces = Spaces(unchecked(istart + 1), iend);
            string subseq = SubS(seqalias, jstart, jend);
            sb.Append("Length = ").Append(Int(ilen)).Append('\n');
            sb.Append("Start = position ").Append(Int(unchecked((int)istart))).Append(" of sequence\n");
            sb.Append("End = position ").Append(Int(unchecked((int)iend))).Append(" of sequence\n\n");
            if (strand == '-')
                sb.Append("Strand = Reversed\n");
            foreach (var g in tags)
                if (GetNote(f, g.Name, out var v))
                    sb.Append(g.Print).Append(" = ").Append(v).Append('\n');
            sb.Append('\n');
            sb.Append(subseq).Append('\n');
            int w1 = unchecked((int)(istart - jstart - 1));
            if (istart == iend)
            {
                sb.Append(Pad("", w1, false)).Append("|\n");
                sb.Append(Pad(Int(unchecked((int)istart)), (int)Math.Min(6u, istart), false)).Append("\n\n");
            }
            else
            {
                sb.Append(Pad("", w1, false)).Append('|').Append(spaces).Append("|\n");
                sb.Append(Pad(Int(unchecked((int)istart)), (int)Math.Min(6u, istart), false)).Append(spaces).Append(Int(unchecked((int)iend))).Append("\n\n");
            }
        }

        Tail(st, t, t.Sequence);
    }

    // the blanks of "for(j = from; j < to; j++) append ' '" over unsigned positions; the program writes about 2^32 of them for a feature with a negative coordinate (fuzztran, reverse frames) and does not finish in practice
    private static string Spaces(uint from, uint to)
    {
        if (to <= from)
            return "";
        uint count = to - from;
        if (count > 50_000_000)
            throw new NotSupportedException($"The program pads this feature with {count} blanks (a negative coordinate read as unsigned); it does not finish in practice, the report is not reproduced.");
        return new string(' ', (int)count);
    }

    // reportWriteMotif
    private static void Motif(State st, EmbossReport t)
    {
        var sb = st.Out;
        string seqalias = SeqAlias(t.Sequence);
        Header(st, t, t.Sequence);
        var tags = GetTags(st);
        int jstar = -1, jplus = -1;
        for (int j = 0; j < tags.Count; j++)
        {
            if (string.Equals(tags[j].Name, "pos", StringComparison.OrdinalIgnoreCase))
            {
                jstar = j;
                sb.Append(tags[j].Print).Append(" at \"*\"\n");
            }

            if (string.Equals(tags[j].Name, "pos2", StringComparison.OrdinalIgnoreCase))
            {
                jplus = j;
                sb.Append(tags[j].Print).Append(" at \"+\"\n");
            }
        }

        if (jstar >= 0 || jplus >= 0)
            sb.Append('\n');
        int i = 0;
        var tmpRes = new MotifTmpRes();
        foreach (var f in t.Features)
        {
            uint istart = unchecked((uint)f.Start), iend = unchecked((uint)f.End);
            int ilen = unchecked((int)(iend - istart + 1));
            string subseq = SubS(seqalias, unchecked((long)(istart - 1)), unchecked((long)(iend - 1)));
            tmpRes.Assign(unchecked(iend - istart + 1));                // ajStrAssignResC(&tmpstr, ilen, "")
            i++;
            sb.Append('(').Append(i).Append(')');
            if (st.O.ShowScore)
                sb.Append(" Score ").Append(Score(st, f));
            sb.Append(" length ").Append(Int(ilen)).Append(" at ").Append(CharName(st)).Append(' ').Append(Int(unchecked((int)istart))).Append("->").Append(Int(unchecked((int)iend)));
            sb.Append('\n');
            uint jposmin = 0, jposmax = 0, jpos = 0, jpos2 = 0;
            if (jstar >= 0 && GetNote(f, tags[jstar].Name, out var sv))
            {
                jpos = ToUint(sv);
                jposmin = jposmax = jpos;
            }

            if (jplus >= 0 && GetNote(f, tags[jplus].Name, out var pv))
            {
                jpos2 = ToUint(pv);
                if (jpos2 < jpos)
                    jposmin = jpos2;
                else
                    jposmax = jpos2;
            }

            if (jposmin > iend)
                jposmin = iend;
            if (jposmax > iend)
                jposmax = jposmin;
            var tmp = new StringBuilder();
            uint jj;
            for (jj = istart; jj < jposmin; jj++)
                tmp.Append(' ');
            if (jpos == jj)
                tmp.Append('*');
            else if (jpos2 == jj)
                tmp.Append('+');
            if (jposmax > jposmin)
            {
                for (jj++; jj < jposmax; jj++)
                    tmp.Append(' ');
                if (jpos == jj)
                    tmp.Append('*');
                else if (jpos2 == jj)
                    tmp.Append('+');
            }

            sb.Append("           ").Append(tmp).Append('\n');
            sb.Append(" Sequence: ").Append(subseq).Append('\n');
            string between = Spaces(unchecked(istart + 1), iend);
            if (istart == iend)
            {
                sb.Append("           |\n");
                sb.Append(Int(unchecked((int)istart)).PadLeft(12)).Append('\n');
            }
            else
            {
                sb.Append("           |").Append(between).Append("|\n");
                sb.Append(Int(unchecked((int)istart)).PadLeft(12)).Append(between).Append(Int(unchecked((int)iend))).Append('\n');
            }

            tmpRes.Assign(unchecked(iend - istart + 1));                // the second ajStrAssignResC(&tmpstr, ilen, "") of the feature
            for (int j = 0; j < tags.Count; j++)
                if (GetNote(f, tags[j].Name, out var v))
                    sb.Append(' ').Append(tags[j].Print).Append(": ").Append(v).Append('\n');
            sb.Append('\n');
        }

        Tail(st, t, t.Sequence);
    }

    /// <summary>
    /// The reserved size of <c>tmpstr</c> of <c>reportWriteMotif</c> (<c>ajStrAssignResC</c> / <c>ajStrSetResRound</c>): the feature whose end is before its start has a length of about 2^32; the rounding of such a size wraps to 0, the reserved size
    /// becomes 0, and the next size of 512 or more doubles it (<c>trysize+=trysize</c>) for ever: the program never stops (measured on 6.6.0: a table with a second feature after one whose location is <c>complement(join(a,b,c,d))</c> with overlapping parts).
    /// </summary>
    private sealed class MotifTmpRes
    {
        private bool _null = true;
        private uint _res;

        private static uint Round(uint i, uint v)
        {
            uint m = i % v;
            return m == 0 ? i : unchecked(i + (v - m));
        }

        public void Assign(uint size)
        {
            uint isize = size == 0 ? 1 : size;
            if (_null)
            {
                _null = false;
                uint r = Round(isize, 32);
                _res = r == 0 ? 32 : Round(r, 16);
                return;
            }

            if (_res >= isize)
                return;
            if (isize >= 512)
            {
                if (_res == 0)
                    throw new NotSupportedException("The reference program never stops on this input: the reserved size of its work string is 0 and ajStrSetResRound doubles it for ever (a feature of a negative length in the motif report).");
                ulong trysize = _res;
                while (trysize < isize)
                    trysize += trysize;
                _res = Round(unchecked((uint)trysize), 512);
            }
            else
            {
                _res = Round(isize, 32);
            }
        }
    }

    // reportWriteTagseq
    private static void TagSeq(State st, EmbossReport t)
    {
        var sb = st.Out;
        string seqalias = SeqAlias(t.Sequence);
        var s = t.Sequence;
        uint seqlen = (uint)s.Length;
        uint seqbeg = unchecked((uint)(s.Begin + s.Offset));
        uint seqend = unchecked((uint)(s.End + s.Offset));
        Header(st, t, s);
        var keys = new List<string>();
        var markup = new Dictionary<string, char[]>(StringComparer.Ordinal);
        bool isnuc = t.TableType != "P" || s.Nucleotide;        // ajFeattableIsNuc of the table or ajSeqIsNuc of the sequence: a reset table (fuzzpro, 2nd sequence on) counts as nucleotide
        var tags = GetTags(st);
        bool isgarnier = tags.Any(g => g.Name == "garnier");
        int margin = 5;
        string[] garnier = ["helix", "sheet", "turns", "coil"];
        if (isgarnier)
        {
            foreach (var gk in garnier)
            {
                markup[gk] = Enumerable.Repeat(' ', (int)seqlen).ToArray();
                keys.Add(gk);
                margin = Math.Max(margin, gk.Length);
            }
        }
        else
        {
            foreach (var f in t.Features)
                if (!markup.ContainsKey(f.Type))
                {
                    markup[f.Type] = Enumerable.Repeat(' ', (int)seqlen).ToArray();
                    keys.Add(f.Type);
                    margin = Math.Max(margin, TagName(f, isnuc).Length);
                }
        }

        var numbers = new StringBuilder("    .").Append(Int(10).PadLeft(5));
        uint istop = seqend > 9 ? seqend - 9 : 0;
        for (uint i = 10; i < istop; i += 10)
            numbers.Append("    .").Append(Int(i + 10).PadLeft(5));
        foreach (var f in t.Features)
        {
            char strand = FeatStrand(f, t.Protein);
            int istart = f.Start, iend = f.End;
            uint ilen = unchecked((uint)(iend - istart + 1));
            char[] m = markup[f.Type];
            PasteCountK(m, istart - 1, strand, ilen);
        }

        uint iinit = unchecked(seqbeg - 1);
        for (uint i = iinit; i < seqend; i += 50)
        {
            uint ilast = Math.Min(i + 50 - 1, seqend - 1);
            uint jlast = Math.Min(i + 50 - 1, (uint)numbers.Length - 1);
            string sub = SubS(numbers.ToString(), i, jlast);
            sb.Append(Pad("", margin, false)).Append(' ').Append(sub).Append('\n');
            sb.Append(Pad("", margin, false)).Append(' ').Append(SubS(seqalias, i, ilast)).Append('\n');
            foreach (var key in keys)
            {
                string line = SubS(new string(markup[key]), i, ilast);
                string label = key;
                if (!isgarnier)
                    label = TypeNameOfKey(t, key, isnuc);
                sb.Append(Pad(label, margin, false)).Append(' ').Append(line).Append('\n');
            }
        }

        Tail(st, t, s);
    }

    private static string TagName(EmbossFeature f, bool nuc) => nuc ? TypeOf(f).NucName : TypeOf(f).Name;

    private static string TypeNameOfKey(EmbossReport t, string key, bool nuc)
    {
        foreach (var f in t.Features)
            if (f.Type == key)
                return TagName(f, nuc);
        return key;
    }

    // ajStrPasteCountK
    private static void PasteCountK(char[] str, long pos, char chr, ulong num)
    {
        long len = str.Length;
        long b = pos < 0 ? len + pos : pos;
        if (b >= len)
            b = len - 1;
        if (b < 0)
            b = 0;
        ulong iend = (ulong)b + num;
        if (iend > (ulong)len)
            iend = (ulong)(len - 1);
        for (long i = b; i < (long)iend; i++)
            str[i] = chr;
    }

    // -------------------------------------------------------------------------- feature-table formats ----

    // ajFeattableWriteGff3 / featDumpGff3
    private static void Gff(State st, EmbossReport t)
    {
        var sb = st.Out;
        var o = st.O;
        sb.Append("##gff-version 3\n");
        sb.Append("##sequence-region ").Append(t.SeqId).Append(" 1 ").Append(t.TableLength != 0 ? t.TableLength : UnsignedRegionEnd(t.Features)).Append('\n');
        sb.Append("#!Date ").Append(EmbossDates.Gff(o.RunDate)).Append('\n');
        sb.Append("#!Type ").Append(t.TableType == "P" ? "Protein" : "DNA").Append('\n');
        sb.Append("#!Source-version EMBOSS 6.6.0.0\n");
        int group = 0;
        foreach (var f in t.Features)
        {
            group++;
            var ti = TypeOf(f);
            char frameChar = f.Frame is > 0 and <= 3 ? (char)('0' + f.Frame - 1) : '.';
            int tmpbeg = f.Start, tmpend = f.End;
            if (FeatStrand(f, t.Protein) == '-' && unchecked((uint)f.End) < unchecked((uint)f.Start))        // compared as unsigned
            {
                tmpbeg = f.End;
                tmpend = f.Start;
            }

            sb.Append(t.SeqId).Append('\t').Append(f.Source).Append('\t').Append(ti.Gff3Name).Append('\t').Append(Int(tmpbeg)).Append('\t').Append(Int(tmpend)).Append('\t');
            if ((float)f.Score != 0f)
                sb.Append(CFormat.General((float)f.Score, 4));
            else
                sb.Append('.');
            sb.Append('\t').Append(FeatStrand(f, t.Protein)).Append('\t').Append(frameChar).Append("\tID=").Append(t.SeqId).Append('.').Append(group);
            if (f.Between)
                sb.Append(";featflags=between:true");
            foreach (var tag in f.Tags)
            {
                string name = tag.Name.ToLowerInvariant();
                sb.Append(';').Append(name).Append('=').Append(PercentEncode(tag.Value));
            }

            sb.Append('\n');
        }
    }

    // ajStrFmtPercentEncodeC(pval, ";=%&,\t")
    private static string PercentEncode(string v)
    {
        var sb = new StringBuilder();
        foreach (char c in v)
        {
            if (c is ';' or '=' or '%' or '&' or ',' or '\t')
                sb.Append('%').Append(((int)c).ToString("X2", CultureInfo.InvariantCulture));
            else
                sb.Append(c);
        }

        return sb.ToString();
    }
}
