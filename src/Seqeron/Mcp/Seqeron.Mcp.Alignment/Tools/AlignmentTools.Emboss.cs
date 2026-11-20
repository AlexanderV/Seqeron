using System.ComponentModel;
using System.Globalization;
using ModelContextProtocol.Server;

namespace Seqeron.Mcp.Alignment.Tools;

// B09 A2-27 / A2-28: the EMBOSS alignment output formats (-aformat) and the stretcher program; delegate to EmbossAlignmentWriter and SequenceAligner.StretcherAlign.
public partial class AlignmentTools
{
    [McpServerTool(Name = "stretcher_align", Title = "Pairwise — Global Alignment, Myers–Miller (EMBOSS stretcher)", ReadOnly = true)]
    [Description("EMBOSS 6.6 stretcher: global alignment by the Myers–Miller linear-space algorithm (Myers & Miller 1988) with affine gaps and INTEGER scores; a gap of length k costs gapOpen + (k-1)·gapExtend (the needle/EMBOSS convention; EMBOSS -gapopen / -gapextend, positive penalties), end gaps are charged like any gap, residues are upper-cased and a residue that is not in the matrix scores 0 (EMBOSS 'bad character'). The rows are the rows of the installed stretcher 6.6.0 (its midpoint tie rules differ from global_align_affine; the score is the same): VERIFICATION: rows and score identical on 689 / 689 committed random DNA / IUPAC / protein pairs (1–260 residues, gap penalties 0–20 / 0–6, 10 matrices). Default matrix EDNAFULL (nucleotide) or EBLOSUM62 (protein) with gapOpen 16 / gapExtend 4 (nucleotide) or 12 / 2 (protein). Use write_alignment (program 'stretcher') for the EMBOSS output files.")]
    public static AlignmentResultDto StretcherAlign(
        [Description("First sequence (-asequence).")] string sequence1,
        [Description("Second sequence (-bsequence).")] string sequence2,
        [Description("Optional EMBOSS data matrix (-datafile): EBLOSUM62, EPAM250, EDNAFULL, ... or the NCBI name. Omit for EDNAFULL (nucleotide) / EBLOSUM62 (protein).")] string? matrix = null,
        [Description("Gap opening penalty (-gapopen, integer >= 0): the cost of a gap of length 1. Omit for 16 (nucleotide) / 12 (protein).")] int? gapOpen = null,
        [Description("Gap extension penalty (-gapextend, integer >= 0): the cost of each further gap position. Omit for 4 (nucleotide) / 2 (protein).")] int? gapExtend = null,
        [Description(FormatParamDescription)] string? format = null,
        [Description(NamesParamDescription)] string[]? names = null,
        CancellationToken cancellationToken = default)
    {
        return FormatPair(StretcherAlignCore(sequence1, sequence2, matrix, gapOpen, gapExtend, cancellationToken), format, sequence1.Length, sequence2.Length, names);
    }

    private static AlignmentResultDto StretcherAlignCore(
        [Description("First sequence (-asequence).")] string sequence1,
        [Description("Second sequence (-bsequence).")] string sequence2,
        [Description("Optional EMBOSS data matrix (-datafile): EBLOSUM62, EPAM250, EDNAFULL, ... or the NCBI name. Omit for EDNAFULL (nucleotide) / EBLOSUM62 (protein).")] string? matrix = null,
        [Description("Gap opening penalty (-gapopen, integer >= 0): the cost of a gap of length 1. Omit for 16 (nucleotide) / 12 (protein).")] int? gapOpen = null,
        [Description("Gap extension penalty (-gapextend, integer >= 0): the cost of each further gap position. Omit for 4 (nucleotide) / 2 (protein).")] int? gapExtend = null,
        CancellationToken cancellationToken = default)
    {
        using var cancelScope = global::Seqeron.Genomics.Alignment.SequenceAligner.BeginCancellationScope(cancellationToken);
        RequireSequence(sequence1, nameof(sequence1));
        RequireSequence(sequence2, nameof(sequence2));
        return ToAlignmentDto(global::Seqeron.Genomics.Alignment.SequenceAligner.StretcherAlign(sequence1, sequence2, ResolveMatrix(matrix), gapOpen, gapExtend));
    }

    [McpServerTool(Name = "write_alignment", Title = "Alignment — Write an alignment file (EMBOSS -aformat)", ReadOnly = true)]
    [Description("Writes an alignment (gapped rows) in an EMBOSS 6.6.0 -aformat output format, byte for byte as the installed programs print it: PAIRWISE formats (exactly two rows; header with program, run date, command line, matrix, gap penalties and Length / Identity / Similarity / Gaps / Score): pair, simple, srspair, srs, markx0, markx1, markx2, markx3, markx10, match, score, sam, tcoffee, unknown, multiple; SEQUENCE-SET formats (any number of rows): fasta, a2m, msf, clustal, mega, meganon, nexus, paup, nexusnon, paupnon, phylip, phylipnon, selex, treecon. The pair formats take the score, the matrix (decides the similarity marks and count), the gap penalties and the program name (needle, water, stretcher, matcher, supermatcher: stretcher prints integer scores). VERIFICATION: 3 303 needle / water / stretcher files (29 formats x 114 alignment sets, faked run dates) and 1 309 seqret files (11 sequence-set formats, 1–9 rows) are identical to the installed EMBOSS 6.6.0. Offsets (residues of the original sequences outside the aligned part) are given per row for the local alignments.")]
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
        if (!global::Seqeron.Genomics.Alignment.EmbossAlignmentWriter.TryParseFormat(format, out var fmt))
            throw new ArgumentException($"Unknown EMBOSS alignment format '{format}'. Supported: unknown, fasta, a2m, msf, clustal, mega, meganon, nexus, paup, nexusnon, paupnon, phylip, phylipnon, selex, treecon, markx0, markx1, markx2, markx3, markx10, match, multiple, pair, simple, sam, score, srs, srspair, tcoffee.", nameof(format));
        if (alignedSequences.Length == 0)
            throw new ArgumentException("At least one aligned row is required.", nameof(alignedSequences));
        string prog = string.IsNullOrWhiteSpace(program) ? "needle" : program.Trim().ToLowerInvariant();
        var date = runDate is null
            ? DateTime.UtcNow
            : DateTime.Parse(runDate, CultureInfo.InvariantCulture, DateTimeStyles.AdjustToUniversal | DateTimeStyles.AssumeUniversal);
        var sub = ResolveMatrix(matrix);
        var options = new global::Seqeron.Genomics.Infrastructure.EmbossAlignOptions
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
        var rows = new List<global::Seqeron.Genomics.Infrastructure.EmbossAlignedSequence>();
        for (int i = 0; i < alignedSequences.Length; i++)
        {
            string name = names is not null && i < names.Length && !string.IsNullOrWhiteSpace(names[i]) ? names[i] : "seq" + (i + 1).ToString(CultureInfo.InvariantCulture);
            rows.Add(new global::Seqeron.Genomics.Infrastructure.EmbossAlignedSequence(
                name,
                alignedSequences[i] ?? throw new ArgumentException($"Row {i} is null.", nameof(alignedSequences)),
                offsets is not null && i < offsets.Length ? offsets[i] : 0,
                offends is not null && i < offends.Length ? offends[i] : 0,
                Description: descriptions is not null && i < descriptions.Length ? descriptions[i] ?? "" : "",
                Usa: name));
        }

        string text;
        if (global::Seqeron.Genomics.Alignment.EmbossAlignmentWriter.IsPairFormat(fmt) || fmt is global::Seqeron.Genomics.Infrastructure.EmbossAlignFormat.Simple or global::Seqeron.Genomics.Infrastructure.EmbossAlignFormat.Srs
            or global::Seqeron.Genomics.Infrastructure.EmbossAlignFormat.Multiple or global::Seqeron.Genomics.Infrastructure.EmbossAlignFormat.Unknown or global::Seqeron.Genomics.Infrastructure.EmbossAlignFormat.TCoffee)
        {
            if (rows.Count != 2)
                throw new ArgumentException($"The pairwise format '{format}' takes exactly two aligned rows (got {rows.Count}); the sequence-set formats take any number.", nameof(alignedSequences));
            text = global::Seqeron.Genomics.Alignment.EmbossAlignmentWriter.Format(fmt, [new global::Seqeron.Genomics.Infrastructure.EmbossPairAlignment(rows[0], rows[1], score)], options);
        }
        else
        {
            text = global::Seqeron.Genomics.Alignment.EmbossAlignmentWriter.Format(fmt, rows, options);
        }

        return new WriteAlignmentResult(text, global::Seqeron.Genomics.Alignment.EmbossAlignmentWriter.FormatName(fmt), "emboss");
    }

    [McpServerTool(Name = "dual_affine_align", Title = "Pairwise — Dual Affine Gaps (ksw2 / minimap2 style)", ReadOnly = true)]
    [Description("Pairwise alignment with the TWO-PIECE (dual) affine gap model of ksw2 / minimap2 (ksw2_extd2): a gap of length k scores max(gapOpen1 + k*gapExtend1, gapOpen2 + k*gapExtend2) (a short-gap and a long-gap regime; all four gap numbers are scores <= 0, ksw2's q + k*e cost is gapOpen = -q, gapExtend = -e; defaults = minimap2 -A 2 -B 4 -O 4,24 -E 2,1), Gotoh's O(m n) program with two gap states per direction, mode 'global' (end to end) or 'local' (Smith-Waterman), optional diagonal band (ksw2 w: only the cells |i - j| <= band; a global alignment needs |m - n| <= band). VERIFICATION: 1 600 committed cases (400 global, 400 local, 400 banded global, 400 banded local; scores equal Biopython 1.88 PairwiseAligner with the callable gap score and an independent Waterman-Smith-Beyer program, every alignment re-scores to its score). Z-drop / end-bonus extension alignment is not included (B09.md A2-29b).")]
    public static AlignmentResultDto DualAffineAlign(
        [Description("First sequence.")] string sequence1,
        [Description("Second sequence.")] string sequence2,
        [Description("'global' (default) or 'local'.")] string mode = "global",
        [Description("Score of an identical pair (minimap2 -A).")] int match = 2,
        [Description("Score of a different pair (minimap2 -B is the penalty: -4).")] int mismatch = -4,
        [Description("Short regime: gap open score (<= 0; minimap2 -O first value as -4).")] int gapOpen1 = -4,
        [Description("Short regime: gap extension score per position (<= 0; -E first value as -2).")] int gapExtend1 = -2,
        [Description("Long regime: gap open score (<= 0; -O second value as -24).")] int gapOpen2 = -24,
        [Description("Long regime: gap extension score per position (<= 0; -E second value as -1).")] int gapExtend2 = -1,
        [Description("Band half-width (cells with |i - j| <= band); -1 = full matrix.")] int band = -1,
        [Description(FormatParamDescription)] string? format = null,
        [Description(NamesParamDescription)] string[]? names = null,
        CancellationToken cancellationToken = default)
    {
        return FormatPair(DualAffineAlignCore(sequence1, sequence2, mode, match, mismatch, gapOpen1, gapExtend1, gapOpen2, gapExtend2, band, cancellationToken), format, sequence1.Length, sequence2.Length, names);
    }

    private static AlignmentResultDto DualAffineAlignCore(
        [Description("First sequence.")] string sequence1,
        [Description("Second sequence.")] string sequence2,
        [Description("'global' (default) or 'local'.")] string mode = "global",
        [Description("Score of an identical pair (minimap2 -A).")] int match = 2,
        [Description("Score of a different pair (minimap2 -B is the penalty: -4).")] int mismatch = -4,
        [Description("Short regime: gap open score (<= 0; minimap2 -O first value as -4).")] int gapOpen1 = -4,
        [Description("Short regime: gap extension score per position (<= 0; -E first value as -2).")] int gapExtend1 = -2,
        [Description("Long regime: gap open score (<= 0; -O second value as -24).")] int gapOpen2 = -24,
        [Description("Long regime: gap extension score per position (<= 0; -E second value as -1).")] int gapExtend2 = -1,
        [Description("Band half-width (cells with |i - j| <= band); -1 = full matrix.")] int band = -1,
        CancellationToken cancellationToken = default)
    {
        using var cancelScope = global::Seqeron.Genomics.Alignment.SequenceAligner.BeginCancellationScope(cancellationToken);
        RequireSequence(sequence1, nameof(sequence1));
        RequireSequence(sequence2, nameof(sequence2));
        var scoring = new global::Seqeron.Genomics.Alignment.DualAffineScoring(match, mismatch, gapOpen1, gapExtend1, gapOpen2, gapExtend2);
        return (mode ?? "global").Trim().ToLowerInvariant() switch
        {
            "global" => ToAlignmentDto(global::Seqeron.Genomics.Alignment.SequenceAligner.GlobalAlignDualAffine(sequence1, sequence2, scoring, band)),
            "local" => ToAlignmentDto(global::Seqeron.Genomics.Alignment.SequenceAligner.LocalAlignDualAffine(sequence1, sequence2, scoring, band)),
            _ => throw new ArgumentException($"mode must be 'global' or 'local', got '{mode}'.", nameof(mode)),
        };
    }
}
