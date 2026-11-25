// SPDX-License-Identifier: LGPL-2.1-or-later
// Port of the EMBOSS 6.6.0 AJAX / NUCLEUS libraries (library source files carry the GNU Lesser General Public License, version 2.1 or later).
// Part of Seqeron.Genomics.Emboss.SequenceFormats, licensed under the GNU Lesser General Public License, version 2.1 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Emboss.SequenceFormats;

using Seqeron.Genomics.Alignment;

/// <summary>
/// The EMBOSS 6.6 report formats (<c>-rformat</c>, the <c>reportFormat[]</c> table of <c>ajreport.c</c>) the report layer writes. <see cref="DiffSeq"/> needs seven report tags (only the <c>diffseq</c> program has them): the five
/// programs of this layer refuse it exactly as the installed programs do. <see cref="Gff3"/> is the alias of <see cref="Gff"/>.
/// </summary>
public enum EmbossReportFormat
{
    /// <summary>GFF version 3 (<c>-rformat gff</c>).</summary>
    Gff = 0,
    /// <summary>Tab-separated columns with a header line (<c>-rformat excel</c>).</summary>
    Excel = 1,
    /// <summary>Simple report (<c>simple</c>): one block of <c>Name: value</c> lines per feature.</summary>
    Simple,
    /// <summary>EMBL feature table (<c>embl</c>; nucleotide features only).</summary>
    Embl,
    /// <summary>GenBank feature table (<c>genbank</c>; nucleotide features only).</summary>
    Genbank,
    /// <summary>Alias of <see cref="Gff"/> (<c>gff3</c>); writes the same text, the report header names the format as typed.</summary>
    Gff3,
    /// <summary>PIR feature lines (<c>pir</c>; protein features only).</summary>
    Pir,
    /// <summary>Swiss-Prot feature lines (<c>swiss</c>; protein features only).</summary>
    Swiss,
    /// <summary>DAS GFF (<c>dasgff</c>, XML).</summary>
    DasGff,
    /// <summary>Debugging trace of the full internal data content (<c>debug</c>).</summary>
    Debug,
    /// <summary>EMBOSS list file of sequence USAs with ranges (<c>listfile</c>).</summary>
    ListFile,
    /// <summary>Motif database hits (<c>dbmotif</c>).</summary>
    DbMotif,
    /// <summary>Differences between a pair of sequences (<c>diffseq</c>); needs seven report tags.</summary>
    DiffSeq,
    /// <summary>lindna / cirdna input (<c>draw</c>).</summary>
    Draw,
    /// <summary>lindna / cirdna input for restriction sites (<c>restrict</c>).</summary>
    Restrict,
    /// <summary>EMBL format feature table with internal tags (<c>feattable</c>).</summary>
    FeatTable,
    /// <summary>Motif report (<c>motif</c>).</summary>
    Motif,
    /// <summary>Simple table with the sequence name (<c>nametable</c>).</summary>
    NameTable,
    /// <summary>Annotated sequence regions (<c>regions</c>).</summary>
    Regions,
    /// <summary>Simple table with the sequence on each line (<c>seqtable</c>, the default of fuzznuc and fuzzpro).</summary>
    SeqTable,
    /// <summary>Simple report format for SRS (<c>srs</c>).</summary>
    Srs,
    /// <summary>Summary of hits (<c>summary</c>): the headers only.</summary>
    Summary,
    /// <summary>Simple table (<c>table</c>, the default of fuzztran, etandem and equicktandem).</summary>
    Table,
    /// <summary>Sequence with features marked below (<c>tagseq</c>).</summary>
    TagSeq,
}

/// <summary>The names, aliases and table positions of the <see cref="EmbossReportFormat"/> members (<c>ajReportFindFormat</c>: case-insensitive).</summary>
public static class EmbossReportFormats
{
    // the order of reportFormat[] (ajreport.c): the position is what the debug report prints as "Format:"
    private static readonly (string Name, EmbossReportFormat Format, bool Nuc, bool Prot, bool Header, int MinTags)[] Table =
    [
        ("simple", EmbossReportFormat.Simple, true, true, true, 0),
        ("embl", EmbossReportFormat.Embl, true, false, false, 0),
        ("genbank", EmbossReportFormat.Genbank, true, false, false, 0),
        ("gff", EmbossReportFormat.Gff, true, true, false, 0),
        ("gff3", EmbossReportFormat.Gff3, true, true, false, 0),
        ("pir", EmbossReportFormat.Pir, false, true, false, 0),
        ("swiss", EmbossReportFormat.Swiss, false, true, false, 0),
        ("dasgff", EmbossReportFormat.DasGff, true, true, false, 0),
        ("debug", EmbossReportFormat.Debug, true, true, true, 0),
        ("listfile", EmbossReportFormat.ListFile, true, true, true, 0),
        ("dbmotif", EmbossReportFormat.DbMotif, true, true, true, 0),
        ("diffseq", EmbossReportFormat.DiffSeq, true, true, true, 7),
        ("draw", EmbossReportFormat.Draw, true, true, true, 0),
        ("restrict", EmbossReportFormat.Restrict, true, true, true, 0),
        ("excel", EmbossReportFormat.Excel, true, true, false, 0),
        ("feattable", EmbossReportFormat.FeatTable, true, true, false, 0),
        ("motif", EmbossReportFormat.Motif, true, true, true, 0),
        ("nametable", EmbossReportFormat.NameTable, true, true, true, 0),
        ("regions", EmbossReportFormat.Regions, true, true, true, 0),
        ("seqtable", EmbossReportFormat.SeqTable, true, true, true, 0),
        ("srs", EmbossReportFormat.Srs, true, true, true, 0),
        ("summary", EmbossReportFormat.Summary, true, true, true, 0),
        ("table", EmbossReportFormat.Table, true, true, true, 0),
        ("tagseq", EmbossReportFormat.TagSeq, true, true, true, 0),
    ];

    /// <summary>The format names in the order of the program's list.</summary>
    public static IReadOnlyList<string> Names => [.. Table.Select(t => t.Name)];

    /// <summary>The format of <paramref name="name"/> (case-insensitive, <c>ajReportFindFormat</c>); false for an unknown name.</summary>
    public static bool TryParse(string? name, out EmbossReportFormat format)
    {
        format = default;
        if (string.IsNullOrEmpty(name))
            return false;
        foreach (var t in Table)
            if (string.Equals(t.Name, name, StringComparison.OrdinalIgnoreCase))
            {
                format = t.Format;
                return true;
            }

        return false;
    }

    /// <summary>The name the program prints for <paramref name="format"/> (<c>Report_format:</c> of the header).</summary>
    public static string Name(EmbossReportFormat format) => Table[IndexOf(format)].Name;

    internal static int IndexOf(EmbossReportFormat format)
    {
        for (int i = 0; i < Table.Length; i++)
            if (Table[i].Format == format)
                return i;
        throw new ArgumentOutOfRangeException(nameof(format), format, "Unknown report format.");
    }

    internal static (bool Nuc, bool Prot, bool Header, int MinTags) Flags(EmbossReportFormat format)
    {
        var t = Table[IndexOf(format)];
        return (t.Nuc, t.Prot, t.Header, t.MinTags);
    }
}

/// <summary>One tag of a feature (<c>AjPTagval</c>): <c>note</c> tags carrying <c>*name value</c> are how the programs store their own values.</summary>
/// <param name="Name">The tag name (<c>note</c>, <c>rpt_type</c>).</param>
/// <param name="Value">The value.</param>
public sealed record EmbossFeatureTag(string Name, string Value);

/// <summary>
/// A feature of an EMBOSS feature table in the form the programs store it (<c>AjPFeature</c>): the report writers read these fields only.
/// </summary>
/// <param name="Type">The internal type (a Sequence Ontology id: <c>SO:0000714</c> nucleotide_motif, <c>SO:0001067</c> polypeptide_motif, <c>SO:0000657</c> repeat_region, <c>SO:0000110</c> sequence_feature).</param>
/// <param name="Source">The program name.</param>
/// <param name="Start">First coordinate as stored (may be negative: fuzztran's reverse frames).</param>
/// <param name="End">Second coordinate as stored.</param>
/// <param name="Score">The score.</param>
/// <param name="Strand">The strand: <c>'\0'</c> forward, <c>'-'</c> reverse (the stored character; <c>ajFeatGetStrand</c> reports <c>+</c> / <c>-</c>).</param>
/// <param name="Frame">The frame (0 = none).</param>
/// <param name="Tags">The tags in the order they were added.</param>
/// <param name="Protein">The feature's own protein flag (<c>feature->Protein</c>) when it differs from its table's (diffseq's <c>ajFeatNewII</c> features in a protein table are not protein features: strand <c>+</c>); null = the table's.</param>
/// <param name="Between">The feature lies between two positions (<c>AJFEATFLAG_BETWEEN_SEQ</c>, <c>ajFeatNewBetween</c>: start = the position, end = position + 1; embl / genbank write <c>203^204</c>, the GFF has <c>featflags=between:true</c>).</param>
public sealed record EmbossFeature(string Type, string Source, int Start, int End, double Score, char Strand, int Frame, IReadOnlyList<EmbossFeatureTag> Tags, bool Between = false, bool? Protein = null)
{
    /// <summary>The feature comes from a feature table of the feature library (featreport, twofeat): its type is an internal type of the library, named by the vocabulary and not by the small table of the programs that make their own features.</summary>
    internal bool FromTable { get; init; }
}

/// <summary>An additional output file named in the report header (<c>ajReportAddFileF</c>: <c># 1: a.gff (Feature file for first sequence)</c>).</summary>
/// <param name="Name">The file name.</param>
/// <param name="Description">The description in brackets.</param>
public sealed record EmbossReportFile(string Name, string Description);

/// <summary>One report tag of the program (<c>taglist</c> of the ACD report: <c>str:pat=Pattern</c>).</summary>
/// <param name="Type">The value type (<c>str</c>, <c>int</c>, <c>float</c>): <c>str</c> columns are left-justified.</param>
/// <param name="Name">The tag name looked up in the features (<c>pat</c>).</param>
/// <param name="Print">The column name printed (<c>Pattern</c>).</param>
public sealed record EmbossReportTag(string Type, string Name, string Print);

/// <summary>
/// What the report writers read from the sequence of a feature table (<c>AjPSeq</c> / <c>AjPFeattable</c>): names, the text the positions refer to and the range state after <c>ajSeqTrim</c>.
/// </summary>
/// <param name="Name">The sequence name (<c>ajSeqGetNameS</c>).</param>
/// <param name="Usa">The USA (<c>fasta::file:name</c> for a FASTA file read by name; with the range <c>[b:e]</c> / <c>[b:e:r]</c> when one was given).</param>
/// <param name="Accession">The accession (empty when the header has none).</param>
/// <param name="Description">The description.</param>
/// <param name="Residues">The sequence text as the report reads it (<c>ajSeqGetSeqS</c>: case as read).</param>
/// <param name="Begin">The range start (<c>ajSeqGetBegin</c>).</param>
/// <param name="End">The range end (<c>ajSeqGetEnd</c>).</param>
/// <param name="Offset">The offset added to the range (<c>ajSeqGetOffset</c>).</param>
/// <param name="Reversed">The sequence is reverse-complemented (<c>ajSeqIsReversedTrue</c>).</param>
/// <param name="Nucleotide">The sequence is nucleotide (<c>ajSeqIsNuc</c>).</param>
public sealed record EmbossReportSequence(string Name, string Usa, string Accession, string Description, string Residues, int Begin, int End, int Offset, bool Reversed, bool Nucleotide)
{
    /// <summary>The sequence length (<c>ajSeqGetLen</c>).</summary>
    public int Length => Residues.Length;
}

/// <summary>The feature table of one sequence as a program hands it to the report (<c>ajReportWrite(report, ftable, seq)</c>).</summary>
/// <param name="Sequence">The sequence.</param>
/// <param name="Features">The features in table order.</param>
/// <param name="Protein">The table is a protein table (<c>ajFeattableNewProt</c>), else nucleotide.</param>
/// <param name="TableType">The table's type text (<c>Feattab-&gt;Type</c>: <c>N</c>, <c>P</c>, or empty for a reset table: fuzzpro resets the table of its second and later sequences, which then reads as <c>p</c> in the GFF header).</param>
/// <param name="Start">The table's start (the debug report prints it).</param>
/// <param name="End">The table's end.</param>
/// <param name="Offset">The table's offset.</param>
/// <param name="TableLength">The table's length (<c>Len</c>; 0 = not set: the GFF header then names the largest coordinate, read as unsigned).</param>
/// <param name="Groups">The groups of the table (<c>Groups</c> of the debug report); null = the number of features.</param>
/// <param name="Tags">The report tags of the program that made the table (<see cref="EmbossReportWriter.Write"/> uses them when no options are given).</param>
/// <param name="WriteWhenEmpty">The program writes the table even without a feature (etandem, equicktandem); used by <see cref="EmbossReportWriter.Write"/>.</param>
/// <param name="TableFormat">The table's format string as the debug report prints it (<c>&lt;null&gt;</c> for a table made by the program, empty for the limited copy <c>-rmaxall</c> / <c>-rmaxseq</c> writes).</param>
public sealed record EmbossReport(EmbossReportSequence Sequence, IReadOnlyList<EmbossFeature> Features, bool Protein, string TableType = "", int Start = 0, int End = 0, int Offset = 0, long TableLength = 0, int? Groups = null, string TableFormat = "<null>", IReadOnlyList<EmbossReportTag>? Tags = null, bool WriteWhenEmpty = false)
{
    /// <summary>The sequence name (what the old single-sequence layer called <c>SequenceName</c>).</summary>
    public string SequenceName => Sequence.Name;

    /// <summary>
    /// The full feature table the report was made from (sub-features, flags, GFF3 tags; <c>ajReportWrite</c> of <c>featreport</c> and <c>twofeat</c>): the report formats that are feature output formats (gff, gff3, embl, genbank, pir, swiss, dasgff) write it with the feature writers
    /// of the library (<c>ajFeattableWrite</c>) instead of the flat <see cref="Features"/>.
    /// </summary>
    internal EmbossSeqFeatureTable? FeatureTable { get; init; }

    /// <summary>The table has no type (<c>ajFeattableIsNuc</c> and <c>ajFeattableIsProt</c> are both true): it can be reported for a sequence of either kind.</summary>
    public bool TypeUnknown { get; init; }

    /// <summary>The name of the feature table when it is not the name of the sequence (<c>thys-&gt;Seqid</c>: the table of a program that clears it between sequences has none after the first): the feature-table writers print it.</summary>
    public string? TableSeqId { get; init; }

    /// <summary>The name the feature-table writers (gff3, gff2, bed, dasgff) print for the table.</summary>
    internal string SeqId => TableSeqId ?? Sequence.Name;
}

/// <summary>The options and context of a report file (<c>AjPReport</c>: the ACD report parameter plus the program's header text).</summary>
public sealed record EmbossReportOptions
{
    /// <summary>The report is <c>multiple</c> (the ACD attribute of the report item: <c>Multi:</c> of the debug format).</summary>
    public bool Multi { get; init; }

    /// <summary>The program name (<c># Program:</c>, the GFF source).</summary>
    public string Program { get; init; } = "";

    /// <summary>The command-line entries after the program name, one per header line (<c>[-sequence] s.fa</c>, <c>-pattern ACGT</c>).</summary>
    public IReadOnlyList<string> CommandLine { get; init; } = [];

    /// <summary>The run date (header <c>Rundate</c>, <c>#!Date</c>).</summary>
    public DateTime RunDate { get; init; } = new(2000, 1, 1, 0, 0, 0, DateTimeKind.Unspecified);

    /// <summary>The report file name (<c>Report_file:</c>).</summary>
    public string ReportFile { get; init; } = "stdout";

    /// <summary>The format (<c>-rformat</c>).</summary>
    public EmbossReportFormat Format { get; init; } = EmbossReportFormat.Gff;

    /// <summary>The format name as typed (printed in <c>Report_format:</c>; null = the canonical name).</summary>
    public string? FormatText { get; init; }

    /// <summary>The report tags of the program (the ACD <c>taglist</c>).</summary>
    public IReadOnlyList<EmbossReportTag> Tags { get; init; } = [];

    /// <summary>The header text of the program (<c>ajReportSetHeaderS</c>: pattern documentation, parameters), lines separated by <c>\n</c>.</summary>
    public string Header { get; init; } = "";

    /// <summary>The tail text of the program (<c>ajReportSetTailS</c>; printed in the tail of the report of each sequence after the hit-limit text), lines separated by <c>\n</c>.</summary>
    public string Tail { get; init; } = "";

    /// <summary>The additional output files (<c>ajReportAddFileF</c>; <c># Additional_files: n</c> in the file header).</summary>
    public IReadOnlyList<EmbossReportFile> AdditionalFiles { get; init; } = [];

    /// <summary><c>-rscoreshow</c> (default true).</summary>
    public bool ShowScore { get; init; } = true;

    /// <summary><c>-rstrandshow</c> (default true).</summary>
    public bool ShowStrand { get; init; } = true;

    /// <summary><c>-raccshow</c>: <c># Accession:</c> in the sequence header.</summary>
    public bool ShowAccession { get; init; }

    /// <summary><c>-rdesshow</c>: <c># Description:</c>.</summary>
    public bool ShowDescription { get; init; }

    /// <summary><c>-rusashow</c>: the USA instead of the name.</summary>
    public bool ShowUsa { get; init; }

    /// <summary><c>-rmaxall</c>: maximum total hits to report (0 = no limit).</summary>
    public int MaxAll { get; init; }

    /// <summary><c>-rmaxseq</c>: maximum hits to report for one sequence (0 = no limit).</summary>
    public int MaxSeq { get; init; }

    /// <summary>The decimals of the scores (<c>precision</c> of the ACD report; 0 for the five programs).</summary>
    public int Precision { get; init; }

    /// <summary>The number of sequences read (<c>ajReportSetSeqstats</c>: <c># Total_sequences</c>); 0 = not set (etandem / equicktandem).</summary>
    public long TotalSequences { get; init; }

    /// <summary>The total length of the sequences read (<c># Total_length</c>); 0 = not set.</summary>
    public long TotalLength { get; init; }

    /// <summary>The program closes the report (<c>ajReportClose</c>: the <c># Total_sequences ...</c> footer); etandem / equicktandem do not (only the end of the file is written by the file close).</summary>
    public bool WriteFooter { get; init; } = true;

    /// <summary>What the feature writers left in the static string of <c>ajfeatwrite.c</c> before the report is written (a PIR report prints it in front of its first note; diffseq writes its feature files first).</summary>
    public string FeatureLeftover { get; init; } = "";

    /// <summary>The program writes the report even when the feature table is empty (etandem, equicktandem); fuzznuc, fuzzpro and fuzztran write only tables with a feature.</summary>
    public bool WriteEmptyTables { get; init; }

    /// <summary>The state of the feature library of the run: the messages the feature writers print (the report formats that are feature formats of a full table).</summary>
    internal EmbossFeatureReadContext? Context { get; init; }

    /// <summary>The ACD qualifier <c>-rformat</c>'s default is not relevant here: the report type letter of <c>ajReportSetType</c> comes from the first table (set by the layer).</summary>
    public string Name { get; init; } = "";

    /// <summary>The report's file extension text printed by the debug report (<c>Extension:</c>).</summary>
    public string Extension { get; init; } = "";
}
