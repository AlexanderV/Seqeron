// SPDX-License-Identifier: GPL-2.0-or-later
// This file is part of Seqeron.Mcp.EmbossPairwiseAlignment, licensed under the terms stated in LICENSE (this directory).
using System.ComponentModel;
using System.Text;
using ModelContextProtocol.Server;
using Seqeron.Mcp.Alignment.Tools;
using static Seqeron.Mcp.Alignment.Tools.AlignmentTools;

namespace Seqeron.Mcp.EmbossPairwiseAlignment.Tools;

public partial class EmbossPairwiseAlignmentTools
{
    [McpServerTool(Name = "stretcher_align", Title = "Pairwise — Global Alignment, Myers–Miller (EMBOSS stretcher)", ReadOnly = true)]
    [Description("EMBOSS 6.6 stretcher: global alignment by the Myers–Miller linear-space algorithm (Myers & Miller 1988) with affine gaps and INTEGER scores; a gap of length k costs gapOpen + (k-1)·gapExtend (the needle/EMBOSS convention; EMBOSS -gapopen / -gapextend, positive penalties), end gaps are charged like any gap, residues are upper-cased and a residue that is not in the matrix scores 0 (EMBOSS 'bad character'). The rows are the rows of the installed stretcher 6.6.0 (its midpoint tie rules differ from global_align_affine; the score is the same): VERIFICATION: rows and score identical on 689 / 689 committed random DNA / IUPAC / protein pairs (1–260 residues, gap penalties 0–20 / 0–6, 10 matrices). Default matrix EDNAFULL (nucleotide) or EBLOSUM62 (protein) with gapOpen 16 / gapExtend 4 (nucleotide) or 12 / 2 (protein). Use write_alignment (program 'stretcher') for the EMBOSS output files. Input sequences are read as EMBOSS reads them (A4-3, one shared reader): only letters and '*' are kept; gap characters (- . ~ ?), '#', '+', '_', digits, blanks and every other character are dropped first, and positions, lengths and scores count the cleaned sequence (e.g. needle of AC-GTTAGC.TT* and ACGTTAGCTT: score 58, length 11, as the program).")]
    public static AlignmentResultDto StretcherAlign(
        [Description("First sequence (-asequence).")] string sequence1,
        [Description("Second sequence (-bsequence).")] string sequence2,
        [Description("Optional EMBOSS data matrix (-datafile): EBLOSUM62, EPAM250, EDNAFULL, ... or the NCBI name. Omit for EDNAFULL (nucleotide) / EBLOSUM62 (protein).")] string? matrix = null,
        [Description("Gap opening penalty (-gapopen, integer >= 0): the cost of a gap of length 1. Omit for 16 (nucleotide) / 12 (protein).")] int? gapOpen = null,
        [Description("Gap extension penalty (-gapextend, integer >= 0): the cost of each further gap position. Omit for 4 (nucleotide) / 2 (protein).")] int? gapExtend = null,
        [Description(FormatParamDescription)] string? format = null,
        [Description(NamesParamDescription)] string[]? names = null,
        CancellationToken cancellationToken = default,
        [Description(SequenceQualifiersDescription)] string? sequenceQualifiers = null)
    {
        if (sequenceQualifiers is not null)
        {
            using var qualifiedScope = global::Seqeron.Genomics.Alignment.SequenceAligner.BeginCancellationScope(cancellationToken);
            var (first, second) = QualifiedPair("stretcher", sequence1, sequence2, sequenceQualifiers, names);
            var result = EmbossPairwisePrograms.StretcherProgram(first, second, ResolveMatrix(matrix), gapOpen, gapExtend);
            bool protein = result.MatrixName != "EDNAFULL" && ResolveMatrix(matrix) is null;
            return QualifiedPairResult("stretcher", result, format, gapOpen ?? (protein ? 12 : 16), gapExtend ?? (protein ? 2 : 4), true, true);
        }

        return FormatPair(StretcherAlignCore(sequence1, sequence2, matrix, gapOpen, gapExtend, cancellationToken), format, EmbossLength(sequence1), EmbossLength(sequence2), names, "stretcher");
    }
    [McpServerTool(Name = "needle_align", Title = "Pairwise — Global Alignment (EMBOSS needle, with end-gap penalties)", ReadOnly = true)]
    [Description("EMBOSS 6.6 needle: global alignment (Needleman-Wunsch with Gotoh's three matrices) in the arithmetic and with the conventions of the program: a gap of length k costs gapOpen + (k-1)*gapExtend (positive penalties, EMBOSS -gapopen / -gapextend, rounded to 16 mantissa bits as EMBOSS does), scores are single-precision floats, and end gaps are FREE unless endWeight (-endweight): then an end gap of length k costs endOpen + (k-1)*endExtend (-endopen / -endextend). Rows (ties included), score and the printed score text are those of the installed needle 6.6.0. Default matrix EDNAFULL (both sequences nucleotide) or EBLOSUM62. exactScore is the single-precision score of the program, score its rounding. VERIFICATION: rows and printed score identical on 887 / 887 committed random DNA / IUPAC / protein pairs (1-200 residues, integer / dyadic / non-dyadic gap penalties, 445 without and 442 with end-gap weights, 10 matrices). Sequences whose matrices would exceed 2^27 = 134 217 728 cells are refused (the program keeps three float matrices and a compass matrix, 16 bytes per cell): use stretcher_align (linear space). Use format / write_alignment (program 'needle') for the EMBOSS output files. Input sequences are read as EMBOSS reads them (A4-3, one shared reader): only letters and '*' are kept; gap characters (- . ~ ?), '#', '+', '_', digits, blanks and every other character are dropped first, and positions, lengths and scores count the cleaned sequence (e.g. needle of AC-GTTAGC.TT* and ACGTTAGCTT: score 58, length 11, as the program).")]
    public static AlignmentResultDto NeedleAlign(
        [Description("First sequence (-asequence), at least one residue.")] string sequence1,
        [Description("Second sequence (-bsequence), at least one residue.")] string sequence2,
        [Description("Optional EMBOSS data matrix (-datafile): EBLOSUM62, EPAM250, EDNAFULL, ... or the NCBI name. Omit for EDNAFULL (nucleotide) / EBLOSUM62 (protein).")] string? matrix = null,
        [Description("Gap opening penalty (-gapopen, 1 - 100): the cost of a gap of length 1.")] double gapOpen = 10,
        [Description("Gap extension penalty (-gapextend, 0 - 10): the cost of each further gap position.")] double gapExtend = 0.5,
        [Description("-endweight: apply the end-gap penalties (default false: end gaps are free).")] bool endWeight = false,
        [Description("-endopen (1 - 100): cost of an end gap of length 1 (used only with endWeight).")] double endOpen = 10,
        [Description("-endextend (0 - 10): cost of each further end-gap position (used only with endWeight).")] double endExtend = 0.5,
        [Description(FormatParamDescription)] string? format = null,
        [Description(NamesParamDescription)] string[]? names = null,
        CancellationToken cancellationToken = default,
        [Description(SequenceQualifiersDescription)] string? sequenceQualifiers = null)
    {
        using var cancelScope = global::Seqeron.Genomics.Alignment.SequenceAligner.BeginCancellationScope(cancellationToken);
        if (sequenceQualifiers is not null)
        {
            var (first, second) = QualifiedPair("needle", sequence1, sequence2, sequenceQualifiers, names);
            return QualifiedPairResult("needle", EmbossPairwisePrograms.NeedleProgram(first, second, ResolveMatrix(matrix), gapOpen, gapExtend, endWeight, endOpen, endExtend), format, gapOpen, gapExtend, true, false);
        }

        RequireSequence(sequence1, nameof(sequence1));
        RequireSequence(sequence2, nameof(sequence2));
        var r = global::Seqeron.Genomics.Emboss.PairwiseAlignment.EmbossPairwisePrograms.NeedleAlign(sequence1, sequence2, ResolveMatrix(matrix), gapOpen, gapExtend, endWeight, endOpen, endExtend);
        return FormatPair(ToAlignmentDto(r), format, EmbossLength(sequence1), EmbossLength(sequence2), names);
    }
    [McpServerTool(Name = "water_align", Title = "Pairwise — Local Alignment (EMBOSS water)", ReadOnly = true)]
    [Description("EMBOSS 6.6 water: local alignment (Smith-Waterman with Gotoh's affine gaps) as the program computes and reports it: a gap of length k costs gapOpen + (k-1)*gapExtend (positive penalties, EMBOSS -gapopen / -gapextend, rounded to 16 mantissa bits), single-precision scores, the cell of the best score is the LAST one scanned within the float epsilon of the maximum, and the walk-back follows the program's own rules; rows (ties included), start / end coordinates and score are those of the installed water 6.6.0. startPosition1 / startPosition2 are the 0-based first aligned residues. Needs sequences of at least two residues (shorter input returns the empty alignment, as the program's maximum is undefined for a one-row matrix); more than 2^28 matrix cells are refused. Default matrix EDNAFULL (both sequences nucleotide) or EBLOSUM62. VERIFICATION: rows, start coordinates and printed score identical on 860 / 860 committed random DNA / IUPAC / protein pairs (10 matrices). Use format / write_alignment (program 'water') for the EMBOSS output files. Input sequences are read as EMBOSS reads them (A4-3, one shared reader): only letters and '*' are kept; gap characters (- . ~ ?), '#', '+', '_', digits, blanks and every other character are dropped first, and positions, lengths and scores count the cleaned sequence (e.g. needle of AC-GTTAGC.TT* and ACGTTAGCTT: score 58, length 11, as the program).")]
    public static AlignmentResultDto WaterAlign(
        [Description("First sequence (-asequence).")] string sequence1,
        [Description("Second sequence (-bsequence).")] string sequence2,
        [Description("Optional EMBOSS data matrix (-datafile): EBLOSUM62, EPAM250, EDNAFULL, ... or the NCBI name. Omit for EDNAFULL (nucleotide) / EBLOSUM62 (protein).")] string? matrix = null,
        [Description("Gap opening penalty (-gapopen, 1 - 100): the cost of a gap of length 1.")] double gapOpen = 10,
        [Description("Gap extension penalty (-gapextend, 0 - 10): the cost of each further gap position.")] double gapExtend = 0.5,
        [Description(FormatParamDescription)] string? format = null,
        [Description(NamesParamDescription)] string[]? names = null,
        CancellationToken cancellationToken = default,
        [Description(SequenceQualifiersDescription)] string? sequenceQualifiers = null,
        [Description(SignificanceParamDescription + " EMBOSS gap costs G / E must be integers: BLAST existence G - E, extension E (not with sequenceQualifiers).")] bool significance = false,
        [Description(DbLengthParamDescription)] long? dbLength = null)
    {
        using var cancelScope = global::Seqeron.Genomics.Alignment.SequenceAligner.BeginCancellationScope(cancellationToken);
        if (sequenceQualifiers is not null)
        {
            var (first, second) = QualifiedPair("water", sequence1, sequence2, sequenceQualifiers, names);
            return QualifiedPairResult("water", EmbossPairwisePrograms.WaterProgram(first, second, ResolveMatrix(matrix), gapOpen, gapExtend), format, gapOpen, gapExtend, false, false);
        }

        RequireSequence(sequence1, nameof(sequence1));
        RequireSequence(sequence2, nameof(sequence2));
        var r = global::Seqeron.Genomics.Emboss.PairwiseAlignment.EmbossPairwisePrograms.WaterAlign(sequence1, sequence2, ResolveMatrix(matrix), gapOpen, gapExtend);
        var dto = ToAlignmentDto(r);
        if (significance)
        {
            if (sequenceQualifiers is not null)
                throw new ArgumentException("significance cannot be combined with sequenceQualifiers.", nameof(significance));
            string used = matrix ?? (EmbossSequenceReader.IsNucleotide(EmbossSequenceReader.Clean(sequence1)) && EmbossSequenceReader.IsNucleotide(EmbossSequenceReader.Clean(sequence2)) ? "EDNAFULL" : "EBLOSUM62");
            dto = WithSignificance(dto, true, used, 0, 0, -(gapOpen - gapExtend), -gapExtend, EmbossLength(sequence1), EmbossLength(sequence2), dbLength, 1);
        }

        return FormatPair(dto, format, EmbossLength(sequence1), EmbossLength(sequence2), names);
    }
}
