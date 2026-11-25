// SPDX-License-Identifier: GPL-2.0-or-later
// This file is part of Seqeron.Mcp.EmbossSequenceFormats, licensed under the terms stated in LICENSE (this directory).
using System.ComponentModel;
using System.Text;
using ModelContextProtocol.Server;
using Seqeron.Mcp.Alignment.Tools;
using static Seqeron.Mcp.Alignment.Tools.AlignmentTools;

namespace Seqeron.Mcp.EmbossSequenceFormats.Tools;

public partial class EmbossSequenceFormatsTools
{
    /// <summary>
    /// Makes the EMBOSS alignment formats (pair, markx10, msf, ...) available to the format parameter of the core alignment tools of the MIT Alignment server that run in this process
    /// (<see cref="AlignmentTools.ExternalFormatWriter"/>): the MIT server itself does not contain the EMBOSS writers.
    /// </summary>
    public static void InstallAlignmentFormats() =>
        AlignmentTools.ExternalFormatWriter = (format, rows, names, score, offsets, offends, program) =>
            WriteAlignment(rows, format, names, null, offsets, offends, score, program: program).Text;

    [McpServerTool(Name = "write_alignment", Title = "Alignment — Write an alignment file (EMBOSS -aformat)", ReadOnly = true)]
    [Description("Writes an alignment (gapped rows) in an EMBOSS 6.6.0 -aformat output format, byte for byte as the installed programs print it: PAIRWISE formats (exactly two rows; header with program, run date, command line, matrix, gap penalties and Length / Identity / Similarity / Gaps / Score): pair, simple, srspair, srs, markx0, markx1, markx2, markx3, markx10, match, score, sam, tcoffee, unknown, multiple; unknown, multiple, simple, srs and tcoffee also take THREE OR MORE rows (the consensus of EMBOSS alignConsStats under the rows, as aligncopy writes them; the tcoffee library lists every pair of rows; 1 040 / 1 040 committed alignments of 3 - 8 rows equal the installed program), while the pair-markup formats refuse more than two rows with the program's message (Alignment format pair specifies at most 2 sequences, alignment has 3); SEQUENCE-SET formats (any number of rows): fasta, a2m, msf, clustal, mega, meganon, nexus, paup, nexusnon, paupnon, phylip, phylipnon, selex, treecon. The pair formats take the score, the matrix (decides the similarity marks and count), the gap penalties and the program name (needle, water, stretcher, matcher, supermatcher: stretcher prints integer scores). VERIFICATION: 3 303 needle / water / stretcher files (29 formats x 114 alignment sets, faked run dates) and 1 309 seqret files (11 sequence-set formats, 1–9 rows) are identical to the installed EMBOSS 6.6.0. Offsets (residues of the original sequences outside the aligned part) are given per row for the local alignments.")]
    public static WriteAlignmentResult WriteAlignment(
        [Description("The aligned rows ('-' = gap); two for the pairwise formats, one or more for the sequence-set formats.")] string[] alignedSequences,
        [Description("The EMBOSS -aformat name: pair, simple, srspair, srs, markx0, markx1, markx2, markx3, markx10, match, score, sam, tcoffee, unknown, multiple, fasta, a2m, msf, clustal, mega, meganon, nexus, paup, nexusnon, paupnon, phylip, phylipnon, selex, treecon.")] string format = "srspair",
        [Description("Sequence names (the first word of the FASTA header); default seq1, seq2, ...")] string[]? names = null,
        [Description("Optional descriptions (one per row; shown by the FASTA writer and with showDescription).")] string[]? descriptions = null,
        [Description("Residues of each original sequence before the aligned part (water: the local start - 1); default 0.")] int[]? offsets = null,
        [Description("Residues of each original sequence after the aligned part; default 0.")] int[]? offends = null,
        [Description("The alignment score (pairwise formats).")] double score = 0,
        [Description("Substitution matrix deciding the Similarity count and the marks (EMBOSS name or NCBI name). Omit: EDNAFULL for nucleotide rows, EBLOSUM62 otherwise.")] string? matrix = null,
        [Description("Gap open penalty printed in the header (Gap_penalty); omit for no line.")] double? gapOpen = null,
        [Description("Gap extension penalty printed in the header (Extend_penalty); omit for no line.")] double? gapExtend = null,
        [Description("The program that produced the alignment: needle, water, stretcher, matcher, supermatcher (header, SAM @PG; stretcher = integer scores, -aglobal Y; water = -aglobal N).")] string program = "needle",
        [Description("-awidth: residues per line (0 = 50).")] int width = 0,
        [Description("Run date (ISO 8601, e.g. 2026-11-19T21:43:53); default: now (UTC).")] string? runDate = null,
        [Description("The command line words printed in the header, one per line (e.g. '-asequence a.fa').")] string[]? commandLine = null,
        [Description("-aaccshow: print (accession) in the header.")] bool showAccession = false,
        [Description("-adesshow: print the description in the header.")] bool showDescription = false,
        [Description("'emboss' (default): the EMBOSS 6.6 writers named above. 'biopython': the Biopython 1.88 Bio.Align writers fasta, a2m, clustal, phylip, nexus, stockholm (byte for byte as Bio.Align.write; the descriptions are written, the Clustal header says Biopython 1.88).")] string dialect = "emboss")
    {
        ArgumentNullException.ThrowIfNull(alignedSequences);
        if (string.Equals(dialect?.Trim(), "biopython", StringComparison.OrdinalIgnoreCase))
        {
            if (!global::Seqeron.Genomics.Alignment.AlignmentFormats.TryParseFormat(format, out var bio))
                throw new ArgumentException($"Unknown Biopython alignment format '{format}'. Supported: fasta, a2m, clustal, phylip, nexus, stockholm.", nameof(format));
            if (alignedSequences.Length == 0)
                throw new ArgumentException("At least one aligned row is required.", nameof(alignedSequences));
            return new WriteAlignmentResult(WriteBiopythonFormat(bio, alignedSequences, names, descriptions), global::Seqeron.Genomics.Alignment.AlignmentFormats.FormatName(bio), "biopython");
        }

        if (!string.IsNullOrWhiteSpace(dialect) && !string.Equals(dialect.Trim(), "emboss", StringComparison.OrdinalIgnoreCase))
            throw new ArgumentException($"dialect must be 'emboss' or 'biopython', got '{dialect}'.", nameof(dialect));
        if (!EmbossAlignmentWriter.TryParseFormat(format, out var fmt))
            throw new ArgumentException($"Unknown EMBOSS alignment format '{format}'. Supported: unknown, fasta, a2m, msf, clustal, mega, meganon, nexus, paup, nexusnon, paupnon, phylip, phylipnon, selex, treecon, markx0, markx1, markx2, markx3, markx10, match, multiple, pair, simple, sam, score, srs, srspair, tcoffee.", nameof(format));
        if (alignedSequences.Length == 0)
            throw new ArgumentException("At least one aligned row is required.", nameof(alignedSequences));
        string prog = string.IsNullOrWhiteSpace(program) ? "needle" : program.Trim().ToLowerInvariant();
        var date = runDate is null
            ? DateTime.UtcNow
            : DateTime.Parse(runDate, CultureInfo.InvariantCulture, DateTimeStyles.AdjustToUniversal | DateTimeStyles.AssumeUniversal);
        var sub = ResolveMatrix(matrix);
        var options = new EmbossAlignOptions
        {
            Program = prog,
            CommandLine = commandLine ?? [],
            RunDate = date,
            Width = width,
            Global = prog != "water",
            IntegerScores = prog == "stretcher",
            StatisticsSetByProgram = prog == "stretcher",
            ShowAccession = showAccession,
            ShowDescription = showDescription,
            Matrix = sub,
            GapPenalty = gapOpen,
            ExtendPenalty = gapExtend,
        };
        var rows = new List<EmbossAlignedSequence>();
        for (int i = 0; i < alignedSequences.Length; i++)
        {
            string name = names is not null && i < names.Length && !string.IsNullOrWhiteSpace(names[i]) ? names[i] : "seq" + (i + 1).ToString(CultureInfo.InvariantCulture);
            rows.Add(new EmbossAlignedSequence(
                name,
                alignedSequences[i] ?? throw new ArgumentException($"Row {i} is null.", nameof(alignedSequences)),
                offsets is not null && i < offsets.Length ? offsets[i] : 0,
                offends is not null && i < offends.Length ? offends[i] : 0,
                Description: descriptions is not null && i < descriptions.Length ? descriptions[i] ?? "" : "",
                Usa: name));
        }

        string text;
        if (EmbossAlignmentWriter.IsPairFormat(fmt) || fmt is EmbossAlignFormat.Simple or EmbossAlignFormat.Srs
            or EmbossAlignFormat.Multiple or EmbossAlignFormat.Unknown or EmbossAlignFormat.TCoffee)
        {
            if (rows.Count > 2)
            {
                var multi = new EmbossAlignmentWriter(fmt, options);
                multi.Write(new EmbossMultiAlignment(rows));
                text = multi.Close();
                return new WriteAlignmentResult(text, EmbossAlignmentWriter.FormatName(fmt), "emboss");
            }

            if (rows.Count != 2)
                throw new ArgumentException($"The pairwise format '{format}' takes exactly two aligned rows (got {rows.Count}); the sequence-set formats take any number.", nameof(alignedSequences));
            text = EmbossAlignmentWriter.Format(fmt, [new EmbossPairAlignment(rows[0], rows[1], score)], options);
        }
        else
        {
            text = EmbossAlignmentWriter.Format(fmt, rows, options);
        }

        return new WriteAlignmentResult(text, EmbossAlignmentWriter.FormatName(fmt), "emboss");
    }
}
