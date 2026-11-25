using System.ComponentModel;
using ModelContextProtocol.Server;

namespace Seqeron.Mcp.Alignment.Tools;

/// <summary>
/// MCP tools for sequence alignment, approximate pattern matching, and
/// genome assembly operations.
/// </summary>
[McpServerToolType]
public partial class AlignmentTools
{
    public const string GapModelDescription =
        "'linear' (default): every gap position costs gapExtend and gapOpen is ignored; totalScore is the column sum of pairs (gap-gap columns 0). " +
        "'affine': gap runs cost gapOpen + k * gapExtend; the pairwise / profile alignments use affine (Gotoh) dynamic programming and totalScore is the sum of pairs with " +
        "natural affine gap costs (Altschul 1989: each pair of rows is projected out of the alignment, gap-gap columns skipped, every gap run of the projection pays gapOpen once). " +
        "For two sequences the affine result is an optimal pairwise alignment; for three or more sequences every method is a heuristic.";

    public const string GapOpenDescription =
        "Gap-open penalty. Used only with gapModel='affine' (a gap of length k then costs gapOpen + k * gapExtend, Gotoh); ignored by the default linear model. May be fractional only with a substitution matrix (type='protein' or matrix).";

    public const string GapExtendDescription =
        "Gap penalty per gap position (linear model: every gap position; affine model: the per-position part of gapOpen + k * gapExtend). May be fractional (e.g. -0.5) only with a substitution matrix (type='protein' or matrix).";

    public const string MsaTypeDescription =
        "Sequence type: 'dna' (default; A, C, G, T scored with match/mismatch) or 'protein' (amino-acid sequences scored with the substitution matrix, default BLOSUM62; match/mismatch are then ignored). " +
        "The protein path shares the dynamic programs, merge, profile and objective code of the DNA path (no second implementation).";

    public const string MsaMatrixDescription =
        "Optional substitution matrix: any built-in protein matrix by EMBOSS or NCBI name (BLOSUM62, PAM250, EBLOSUM45, EPAM120, ...) for type='protein', a nucleotide matrix (EDNAFULL, EDNAMAT, EDNASIMPLE, ENUC.4.4) for type='dna'. When given (or type='protein'), match/mismatch are ignored, residues are scored by the published table " +
        "(case-insensitive, any symbol of the matrix), gapOpen/gapExtend may be fractional (EMBOSS needle numbers gapopen G / gapextend E = gapOpen -(G - E), gapExtend -E) and the result also carries the real-valued exactTotalScore " +
        "(totalScore = rounded). Verified against independent Python implementations of every method, Biopython 1.88 pairwise optima and exhaustive enumeration of tiny sets.";

    // Utility holder for static MCP tools; never instantiated (S1118).
    private AlignmentTools() { }

    #region SequenceAligner

    [McpServerTool(Name = "global_align", Title = "Pairwise — Global Alignment (Needleman–Wunsch)", ReadOnly = true)]
    [Description("Performs global pairwise alignment using the Needleman–Wunsch dynamic programming algorithm (aligns full length of both sequences end-to-end). Linear gap model: every gap position costs gapExtend; gapOpen is IGNORED by this tool (use global_align_affine for gapOpen + k·gapExtend). With matrix = BLOSUM62 | PAM250 | EDNAFULL the alignment is scored by that published residue-pair table (protein or IUPAC DNA, any symbol of the matrix) instead of match/mismatch and gapOpen/gapExtend may be fractional (EMBOSS needle/water 6.6.0 = Biopython PairwiseAligner = parasail on >1000 oracle rows; EMBOSS -gapopen G -gapextend E = gapOpen -(G - E), gapExtend -E); exactScore is then the real-valued score (score = rounded). Memory: inputs above 25,000,000 dynamic-programming cells ((len1+1)*(len2+1), e.g. two 5 kb sequences) are aligned automatically in exact LINEAR SPACE (Hirschberg 1975 / Myers-Miller 1988 / Huang-Miller 1991-style locate-then-solve) - the returned score, rows and coordinates are identical to the full-matrix result (ties included), only the running time grows (about 1.3-2x) and the memory stays O(len1 + len2) instead of O(len1 * len2) (20 kb x 20 kb affine: 111 MB peak instead of 4.8 GB).")]
    public static AlignmentResultDto GlobalAlign(
        [Description("First sequence (DNA, uppercased internally).")] string sequence1,
        [Description("Second sequence.")] string sequence2,
        [Description("Score for matching base.")] int match = 1,
        [Description("Score for mismatching base.")] int mismatch = -1,
        [Description("Gap-open score (ignored by this linear-gap tool; kept for schema compatibility). May be fractional when matrix is given.")] double gapOpen = -2,
        [Description("Gap score charged per gap position (<= 0). May be fractional (e.g. -0.5) when matrix is given.")] double gapExtend = -1,
        [Description(MatrixParamDescription)] string? matrix = null,
        [Description(PairFormatParamDescription)] string? format = null,
        [Description(NamesParamDescription)] string[]? names = null,
        CancellationToken cancellationToken = default,
        [Description(TrimEndGapsParamDescription)] bool trimEndGaps = false)
    {
        return FormatPair(GlobalAlignCore(sequence1, sequence2, match, mismatch, gapOpen, gapExtend, matrix, cancellationToken), format, sequence1.Length, sequence2.Length, names, sequence1: sequence1, sequence2: sequence2, trimEndGaps: trimEndGaps);
    }

    private static AlignmentResultDto GlobalAlignCore(
        [Description("First sequence (DNA, uppercased internally).")] string sequence1,
        [Description("Second sequence.")] string sequence2,
        [Description("Score for matching base.")] int match = 1,
        [Description("Score for mismatching base.")] int mismatch = -1,
        [Description("Gap-open score (ignored by this linear-gap tool; kept for schema compatibility). May be fractional when matrix is given.")] double gapOpen = -2,
        [Description("Gap score charged per gap position (<= 0). May be fractional (e.g. -0.5) when matrix is given.")] double gapExtend = -1,
        [Description(MatrixParamDescription)] string? matrix = null,
        CancellationToken cancellationToken = default)
    {
        using var cancelScope = global::Seqeron.Genomics.Alignment.SequenceAligner.BeginCancellationScope(cancellationToken);
        RequireSequence(sequence1, nameof(sequence1));
        RequireSequence(sequence2, nameof(sequence2));

        var sub = ResolveMatrix(matrix);
        if (sub is not null)
            return ToAlignmentDto(global::Seqeron.Genomics.Alignment.SequenceAligner.GlobalAlign(sequence1, sequence2, new global::Seqeron.Genomics.Infrastructure.MatrixScoring(sub, gapOpen, gapExtend)));

        var scoring = new global::Seqeron.Genomics.Infrastructure.ScoringMatrix(match, mismatch, RequireInteger(gapOpen, nameof(gapOpen)), RequireInteger(gapExtend, nameof(gapExtend)));
        var r = global::Seqeron.Genomics.Alignment.SequenceAligner.GlobalAlign(sequence1, sequence2, scoring);
        return ToAlignmentDto(r);
    }

    [McpServerTool(Name = "local_align", Title = "Pairwise — Local Alignment (Smith–Waterman)", ReadOnly = true)]
    [Description("Performs local pairwise alignment using the Smith–Waterman algorithm (finds best-scoring substring alignment, with zero floor). Linear gap model: every gap position costs gapExtend; gapOpen is IGNORED by this tool (use local_align_affine for gapOpen + k·gapExtend). With matrix = BLOSUM62 | PAM250 | EDNAFULL the alignment is scored by that published residue-pair table (protein or IUPAC DNA, any symbol of the matrix) instead of match/mismatch and gapOpen/gapExtend may be fractional (EMBOSS needle/water 6.6.0 = Biopython PairwiseAligner = parasail on >1000 oracle rows; EMBOSS -gapopen G -gapextend E = gapOpen -(G - E), gapExtend -E); exactScore is then the real-valued score (score = rounded). Memory: inputs above 25,000,000 dynamic-programming cells ((len1+1)*(len2+1), e.g. two 5 kb sequences) are aligned automatically in exact LINEAR SPACE (Hirschberg 1975 / Myers-Miller 1988 / Huang-Miller 1991-style locate-then-solve) - the returned score, rows and coordinates are identical to the full-matrix result (ties included), only the running time grows (about 1.3-2x) and the memory stays O(len1 + len2) instead of O(len1 * len2) (20 kb x 20 kb affine: 111 MB peak instead of 4.8 GB).")]
    public static AlignmentResultDto LocalAlign(
        [Description("First sequence.")] string sequence1,
        [Description("Second sequence.")] string sequence2,
        [Description("Score for matching base.")] int match = 1,
        [Description("Score for mismatching base.")] int mismatch = -1,
        [Description("Gap-open score (ignored by this linear-gap tool; kept for schema compatibility). May be fractional when matrix is given.")] double gapOpen = -2,
        [Description("Gap score charged per gap position (<= 0). May be fractional (e.g. -0.5) when matrix is given.")] double gapExtend = -1,
        [Description(MatrixParamDescription)] string? matrix = null,
        [Description(PairFormatParamDescription)] string? format = null,
        [Description(NamesParamDescription)] string[]? names = null,
        CancellationToken cancellationToken = default,
        [Description(SignificanceParamDescription + " This linear-gap tool has gap existence 0 and extension -gapExtend (the linear rows of the blastn tables).")] bool significance = false,
        [Description(DbLengthParamDescription)] long? dbLength = null,
        [Description(DbSequencesParamDescription)] long dbSequences = 1)
    {
        var dto = WithSignificance(LocalAlignCore(sequence1, sequence2, match, mismatch, gapOpen, gapExtend, matrix, cancellationToken), significance, matrix, match, mismatch, 0, gapExtend, sequence1.Length, sequence2.Length, dbLength, dbSequences);
        return FormatPair(dto, format, sequence1.Length, sequence2.Length, names, sequence1: sequence1, sequence2: sequence2);
    }

    private static AlignmentResultDto LocalAlignCore(
        [Description("First sequence.")] string sequence1,
        [Description("Second sequence.")] string sequence2,
        [Description("Score for matching base.")] int match = 1,
        [Description("Score for mismatching base.")] int mismatch = -1,
        [Description("Gap-open score (ignored by this linear-gap tool; kept for schema compatibility). May be fractional when matrix is given.")] double gapOpen = -2,
        [Description("Gap score charged per gap position (<= 0). May be fractional (e.g. -0.5) when matrix is given.")] double gapExtend = -1,
        [Description(MatrixParamDescription)] string? matrix = null,
        CancellationToken cancellationToken = default)
    {
        using var cancelScope = global::Seqeron.Genomics.Alignment.SequenceAligner.BeginCancellationScope(cancellationToken);
        RequireSequence(sequence1, nameof(sequence1));
        RequireSequence(sequence2, nameof(sequence2));

        var sub = ResolveMatrix(matrix);
        if (sub is not null)
            return ToAlignmentDto(global::Seqeron.Genomics.Alignment.SequenceAligner.LocalAlign(sequence1, sequence2, new global::Seqeron.Genomics.Infrastructure.MatrixScoring(sub, gapOpen, gapExtend)));

        var scoring = new global::Seqeron.Genomics.Infrastructure.ScoringMatrix(match, mismatch, RequireInteger(gapOpen, nameof(gapOpen)), RequireInteger(gapExtend, nameof(gapExtend)));
        var r = global::Seqeron.Genomics.Alignment.SequenceAligner.LocalAlign(sequence1, sequence2, scoring);
        return ToAlignmentDto(r);
    }

    [McpServerTool(Name = "semi_global_align", Title = "Pairwise — Semi-Global (Fitting) Alignment", ReadOnly = true)]
    [Description("Performs semi-global (fitting / glocal) pairwise alignment with free end gaps in sequence2 by default (endGaps = 'Fitting'; useful for fitting a shorter query into a longer reference); the optional endGaps parameter selects any of the 16 end-gap policies (overlap / suffix-prefix in both orientations, sequence2 inside sequence1, all four ends free = EMBOSS needle's default -endweight N, none = global; Biopython PairwiseAligner left/right_insertion/deletion_score = 0, parasail sg_* variants). Linear gap model: every gap position costs gapExtend; gapOpen is IGNORED by this tool (use semi_global_align_affine for gapOpen + k·gapExtend). With the default endGaps startPosition2/endPosition2 report the fitted interval of sequence2 (0-based inclusive; the free flanks lie outside it; an empty interval has endPosition2 = startPosition2 - 1) and startPosition1/endPosition1 span the whole query; for any endGaps they are the consumed intervals of each sequence. With matrix = BLOSUM62 | PAM250 | EDNAFULL the alignment is scored by that published residue-pair table (protein or IUPAC DNA, any symbol of the matrix) instead of match/mismatch and gapOpen/gapExtend may be fractional (EMBOSS needle/water 6.6.0 = Biopython PairwiseAligner = parasail on >1000 oracle rows; EMBOSS -gapopen G -gapextend E = gapOpen -(G - E), gapExtend -E); exactScore is then the real-valued score (score = rounded). Matrix mode accepts protein sequences (no DNA validation). Memory: inputs above 25,000,000 dynamic-programming cells ((len1+1)*(len2+1), e.g. two 5 kb sequences) are aligned automatically in exact LINEAR SPACE (Hirschberg 1975 / Myers-Miller 1988 / Huang-Miller 1991-style locate-then-solve) - the returned score, rows and coordinates are identical to the full-matrix result (ties included), only the running time grows (about 1.3-2x) and the memory stays O(len1 + len2) instead of O(len1 * len2) (20 kb x 20 kb affine: 111 MB peak instead of 4.8 GB).")]
    public static AlignmentResultDto SemiGlobalAlign(
        [Description("Query sequence (typically shorter; any letters - N, IUPAC, protein - uppercased internally).")] string sequence1,
        [Description("Reference sequence (typically longer; any letters, uppercased internally).")] string sequence2,
        [Description("Score for matching base.")] int match = 1,
        [Description("Score for mismatching base.")] int mismatch = -1,
        [Description("Gap-open score (ignored by this linear-gap tool; kept for schema compatibility). May be fractional when matrix is given.")] double gapOpen = -2,
        [Description("Gap score charged per gap position (<= 0). May be fractional (e.g. -0.5) when matrix is given.")] double gapExtend = -1,
        [Description(MatrixParamDescription)] string? matrix = null,
        [Description(EndGapsParamDescription)] string? endGaps = null,
        [Description(PairFormatParamDescription)] string? format = null,
        [Description(NamesParamDescription)] string[]? names = null,
        CancellationToken cancellationToken = default,
        [Description(TrimEndGapsParamDescription)] bool trimEndGaps = false)
    {
        return FormatPair(SemiGlobalAlignCore(sequence1, sequence2, match, mismatch, gapOpen, gapExtend, matrix, endGaps, cancellationToken), format, sequence1.Length, sequence2.Length, names, sequence1: sequence1, sequence2: sequence2, queryFirst: true, trimEndGaps: trimEndGaps);
    }

    private static AlignmentResultDto SemiGlobalAlignCore(
        [Description("Query sequence (typically shorter; any letters - N, IUPAC, protein - uppercased internally).")] string sequence1,
        [Description("Reference sequence (typically longer; any letters, uppercased internally).")] string sequence2,
        [Description("Score for matching base.")] int match = 1,
        [Description("Score for mismatching base.")] int mismatch = -1,
        [Description("Gap-open score (ignored by this linear-gap tool; kept for schema compatibility). May be fractional when matrix is given.")] double gapOpen = -2,
        [Description("Gap score charged per gap position (<= 0). May be fractional (e.g. -0.5) when matrix is given.")] double gapExtend = -1,
        [Description(MatrixParamDescription)] string? matrix = null,
        [Description(EndGapsParamDescription)] string? endGaps = null,
        CancellationToken cancellationToken = default)
    {
        using var cancelScope = global::Seqeron.Genomics.Alignment.SequenceAligner.BeginCancellationScope(cancellationToken);
        var ends = ResolveEndGaps(endGaps);
        var sub = ResolveMatrix(matrix);
        if (sub is not null)
            return ToAlignmentDto(global::Seqeron.Genomics.Alignment.SequenceAligner.SemiGlobalAlign(RequireSequence(sequence1, nameof(sequence1)), RequireSequence(sequence2, nameof(sequence2)),
                new global::Seqeron.Genomics.Infrastructure.MatrixScoring(sub, gapOpen, gapExtend), ends));

        // R14-2: any string (N, IUPAC, protein letters) like global_align / local_align: the string overload of the library, not the strict ACGT DnaSequence
        var scoring = new global::Seqeron.Genomics.Infrastructure.ScoringMatrix(match, mismatch, RequireInteger(gapOpen, nameof(gapOpen)), RequireInteger(gapExtend, nameof(gapExtend)));
        var r = global::Seqeron.Genomics.Alignment.SequenceAligner.SemiGlobalAlign(RequireSequence(sequence1, nameof(sequence1)), RequireSequence(sequence2, nameof(sequence2)), scoring, ends);
        return ToAlignmentDto(r);
    }

    [McpServerTool(Name = "alignment_statistics", Title = "Alignment — Statistics", ReadOnly = true)]
    [Description("Computes identity / similarity / gap statistics for a pairwise alignment (EMBOSS needle convention: every percentage is over the full alignment length, gap columns included; a column is a gap column when either row has '-'; residues compared case-insensitively). Similarity counts columns with a POSITIVE substitution score (BLAST positives): with matrix = BLOSUM62 | PAM250 | EDNAFULL the published residue-pair table is used, otherwise the scalar match/mismatch model (similarity = identity for match > 0 >= mismatch). Also returns positives and gapOpens (number of gap runs, BLAST gapopen) and the BLAST+ report values: blastPositives (identities + non-identical positives), the INTEGER percentages of the BLAST pairwise report (blastIdentityPercent, blastPositivesPercent, blastGapPercent: round half up, never 100 unless every column matches; NCBI GetPercentMatch), blastPident (tabular pident, 3 decimals) and blastReport (' Identities = a/L (p%), [Positives = b/L (q%),] Gaps = g/L (r%)'; Positives only with a protein matrix).")]
    public static AlignmentStatisticsDto AlignmentStatistics(
        [Description("Aligned representation of sequence1 (with '-' gaps).")] string alignedSequence1,
        [Description("Aligned representation of sequence2 (with '-' gaps).")] string alignedSequence2,
        [Description("Optional substitution matrix deciding which non-gap columns are similar: any built-in matrix by EMBOSS or NCBI name - BLOSUM62, PAM250, EBLOSUM30 ... EBLOSUM90, EPAM10 ... EPAM500 (protein) or EDNAFULL, EDNAMAT, EDNASIMPLE, ENUC.4.4 (IUPAC DNA). Omit to use match/mismatch.")] string? matrix = null,
        [Description("Scalar model (used when matrix is omitted): score of an identical pair; the pair is 'similar' iff this is > 0.")] int match = 1,
        [Description("Scalar model (used when matrix is omitted): score of a non-identical pair; the pair is 'similar' iff this is > 0.")] int mismatch = -1)
    {
        var ar = ToPairwise(alignedSequence1, alignedSequence2);
        var sub = ResolveMatrix(matrix);
        var s = sub is null
            ? global::Seqeron.Genomics.Alignment.SequenceAligner.CalculateStatistics(
                ar, new global::Seqeron.Genomics.Infrastructure.ScoringMatrix(match, mismatch, 0, 0))
            : global::Seqeron.Genomics.Alignment.SequenceAligner.CalculateStatistics(ar, sub);
        return new AlignmentStatisticsDto(
            s.Matches, s.Mismatches, s.Gaps, s.AlignmentLength,
            s.Identity, s.Similarity, s.GapPercent, s.Positives, s.GapOpens,
            s.BlastPositives, s.BlastPercentIdentity, s.BlastPercentPositives, s.BlastPercentGaps,
            s.BlastTabularPidentText, s.ToBlastString(protein: sub is not null && !sub.IsNucleotide));
    }

    [McpServerTool(Name = "format_alignment", Title = "Alignment — Pretty Format", ReadOnly = true)]
    [Description("Renders a pairwise alignment as three-line blocks in the EMBOSS needle/water (srspair) style: the markup line has '|' for identical columns, ':' for similar columns (positive substitution score), '.' for mismatches and ' ' for gap columns; wrapped to a configurable line width. matrix (BLOSUM62 | PAM250 | EDNAFULL) or the scalar match/mismatch decide which mismatches are similar.")]
    public static FormatAlignmentResult FormatAlignment(
        [Description("Aligned representation of sequence1.")] string alignedSequence1,
        [Description("Aligned representation of sequence2.")] string alignedSequence2,
        [Description("Maximum characters per line in the formatted output (>= 1).")] int lineWidth = 60,
        [Description("Optional substitution matrix deciding which non-gap columns are similar: any built-in matrix by EMBOSS or NCBI name (BLOSUM62, PAM250, EBLOSUM45, EPAM120, EDNAFULL, ...). Omit to use match/mismatch.")] string? matrix = null,
        [Description("Scalar model (used when matrix is omitted): score of an identical pair.")] int match = 1,
        [Description("Scalar model (used when matrix is omitted): score of a non-identical pair; > 0 marks mismatches ':' (similar).")] int mismatch = -1)
    {
        if (string.IsNullOrEmpty(alignedSequence1))
            throw new ArgumentException("Aligned sequence cannot be null or empty.", nameof(alignedSequence1));
        if (string.IsNullOrEmpty(alignedSequence2))
            throw new ArgumentException("Aligned sequence cannot be null or empty.", nameof(alignedSequence2));
        if (lineWidth < 1)
            throw new ArgumentOutOfRangeException(nameof(lineWidth), "lineWidth must be >= 1.");

        var ar = ToPairwise(alignedSequence1, alignedSequence2);
        var sub = ResolveMatrix(matrix);
        var formatted = sub is null
            ? global::Seqeron.Genomics.Alignment.SequenceAligner.FormatAlignment(
                ar, lineWidth, new global::Seqeron.Genomics.Infrastructure.ScoringMatrix(match, mismatch, 0, 0))
            : global::Seqeron.Genomics.Alignment.SequenceAligner.FormatAlignment(ar, lineWidth, sub);
        return new FormatAlignmentResult(formatted);
    }

    private static global::Seqeron.Genomics.Infrastructure.AlignmentResult ToPairwise(string alignedSequence1, string alignedSequence2)
    {
        if (string.IsNullOrEmpty(alignedSequence1))
            throw new ArgumentException("Aligned sequence cannot be null or empty.", nameof(alignedSequence1));
        if (string.IsNullOrEmpty(alignedSequence2))
            throw new ArgumentException("Aligned sequence cannot be null or empty.", nameof(alignedSequence2));
        if (alignedSequence1.Length != alignedSequence2.Length)
            throw new ArgumentException("Aligned sequences must have equal length.", nameof(alignedSequence2));

        return new global::Seqeron.Genomics.Infrastructure.AlignmentResult(
            alignedSequence1, alignedSequence2, 0,
            global::Seqeron.Genomics.Infrastructure.AlignmentType.Global,
            0, 0, 0, 0);
    }

    private const string MatrixParamDescription =
        "Optional substitution matrix: any built-in matrix by EMBOSS, NCBI or Biopython name - BLOSUM62, PAM250, EBLOSUM30 ... EBLOSUM90, EPAM10 ... EPAM500, BLASTP (BLOSUM62 with the residues U = selenocysteine, O = pyrrolysine and J = leucine / isoleucine, as blastp scores them), GONNET1992, JONES, DAYHOFF, MDM78, BENNER6 / 22 / 74, JOHNSON, RISLER, GENETIC, FENG, LEVIN, MCLACHLAN, RAO, STR (protein; the fractional ones are scored exactly) or EDNAFULL, EDNAMAT, EDNASIMPLE, ENUC.4.4, BLASTN, MEGABLAST, HOXD70, TRANS (IUPAC DNA); or THE TEXT of a matrix instead of a name: several lines in the NCBI / EMBOSS / Biopython data-file format (a header line of symbols, then one row per symbol: the symbol and its scores; # lines are comments; scores may be decimals). When given, match/mismatch are ignored, residues are scored by the published table (any protein or IUPAC symbol of the matrix; case-insensitive) and gapOpen/gapExtend may be fractional (the EMBOSS needle/water models: gapOpen = -(G - E), gapExtend = -E for -gapopen G -gapextend E). The result then also carries the real-valued exactScore.";

    /// <summary>
    /// The pairwise tools DELEGATE to the library, whose contract (F39) accepts an empty sequence (the all-gap alignment of the other sequence;
    /// no top-N / empty co-optimal alignments for the local tools) - the tool only refuses a missing (null) value.
    /// </summary>
    public static string RequireSequence(string sequence, string name) =>
        sequence ?? throw new ArgumentException("Sequence is required (null); an empty string is accepted and aligns to all gaps.", name);

    public static int RequireInteger(double value, string name) =>
        value == Math.Round(value) && Math.Abs(value) <= int.MaxValue
            ? (int)value
            : throw new ArgumentException($"{name} = {value} is not an integer: fractional gap scores need a substitution matrix (parameter 'matrix').", name);

    private const string EndGapsParamDescription =
        "Which end gaps are free (default 'Fitting' = the classic semi-global: free leading/trailing sequence2 residues). Comma-separated, case-insensitive names of the flags FreeStart1 (leading sequence1 residues free), FreeEnd1 (trailing sequence1), FreeStart2 (leading sequence2), FreeEnd2 (trailing sequence2) or of a preset: None (= global, all end gaps penalised; result tagged Global), Fitting (FreeStart2+FreeEnd2; sequence1 inside sequence2; parasail sg_dx, Biopython end_insertion_score=0), FittingSequence2 (FreeStart1+FreeEnd1; sequence2 inside sequence1; sg_qx), OverlapSuffix1Prefix2 (FreeStart1+FreeEnd2: a suffix of sequence1 overlaps a prefix of sequence2; sg_qb_de), OverlapSuffix2Prefix1 (FreeEnd1+FreeStart2; sg_qe_db), All (all four free; parasail sg, EMBOSS needle's default -endweight N). Example: 'FreeStart1,FreeEnd2'. startPosition1/endPosition1 and startPosition2/endPosition2 are the consumed intervals of each sequence (the free flanks lie outside them; an empty interval has end = start - 1); a free flank is shown against '-' in the aligned rows.";

    private static global::Seqeron.Genomics.Infrastructure.EndGaps ResolveEndGaps(string? value)
    {
        if (string.IsNullOrWhiteSpace(value))
            return global::Seqeron.Genomics.Infrastructure.EndGaps.Fitting;
        var ends = global::Seqeron.Genomics.Infrastructure.EndGaps.None;
        foreach (var token in value.Split(new[] { ',', '|', '+', ' ' }, StringSplitOptions.RemoveEmptyEntries))
        {
            if (char.IsDigit(token[0]) || !Enum.TryParse<global::Seqeron.Genomics.Infrastructure.EndGaps>(token, ignoreCase: true, out var flag))
                throw new ArgumentException(
                    $"Unknown endGaps '{token}'. Use None, Fitting, FittingSequence2, OverlapSuffix1Prefix2, OverlapSuffix2Prefix1, All or a comma-separated list of FreeStart1, FreeEnd1, FreeStart2, FreeEnd2.", nameof(value));
            ends |= flag;
        }
        return ends;
    }

    private static global::Seqeron.Genomics.Infrastructure.SubstitutionMatrix? ResolveMatrix(string? name)
    {
        if (string.IsNullOrWhiteSpace(name))
            return null;
        // R14-8: a text with a line break is the matrix itself (NCBI / EMBOSS / Biopython format), not a name
        if (name.Contains('\n') || name.Contains('\r'))
            return global::Seqeron.Genomics.Infrastructure.SubstitutionMatrix.Parse(name, "USER");
        return global::Seqeron.Genomics.Infrastructure.SubstitutionMatrix.TryGet(name, out var m)
            ? m
            : throw new ArgumentException($"Unknown matrix '{name}'. Supported: the EMBOSS data matrices (EBLOSUM30 ... EBLOSUM90, EBLOSUM62-12, EBLOSUMN, EPAM10 ... EPAM500, EDNAFULL, EDNAMAT, EDNASIMPLE, ENUC.4.4) and their NCBI names (BLOSUM62, PAM250, BLOSUM45, PAM30, NUC.4.4, ...), the Biopython matrices BLASTN, BLASTP, FENG, HOXD70, LEVIN, MCLACHLAN, MEGABLAST, RAO, STR, TRANS and the fractional BENNER6 / 22 / 74, DAYHOFF, GENETIC, GONNET1992, JOHNSON, JONES, MDM78, RISLER, or the text of a matrix (several lines).", nameof(name));
    }

    [McpServerTool(Name = "global_align_affine", Title = "Pairwise — Global Alignment, Affine Gaps (Gotoh)", ReadOnly = true)]
    [Description("Global pairwise alignment with affine gap costs (Gotoh 1982): a gap of length k scores gapOpen + k·gapExtend (the BLAST+ / Biopython open_gap_score = gapOpen + gapExtend convention), so BLAST scoring (match 2, mismatch -3, gapOpen -5, gapExtend -2) can be applied. Border gaps are charged an opening each (Flouri et al. 2015). Equals Biopython PairwiseAligner(mode='global') and parasail nw_trace. With matrix = BLOSUM62 | PAM250 | EDNAFULL the alignment is scored by that published residue-pair table (protein or IUPAC DNA, any symbol of the matrix) instead of match/mismatch and gapOpen/gapExtend may be fractional (EMBOSS needle/water 6.6.0 = Biopython PairwiseAligner = parasail on >1000 oracle rows; EMBOSS -gapopen G -gapextend E = gapOpen -(G - E), gapExtend -E); exactScore is then the real-valued score (score = rounded). Memory: inputs above 25,000,000 dynamic-programming cells ((len1+1)*(len2+1), e.g. two 5 kb sequences) are aligned automatically in exact LINEAR SPACE (Hirschberg 1975 / Myers-Miller 1988 / Huang-Miller 1991-style locate-then-solve) - the returned score, rows and coordinates are identical to the full-matrix result (ties included), only the running time grows (about 1.3-2x) and the memory stays O(len1 + len2) instead of O(len1 * len2) (20 kb x 20 kb affine: 111 MB peak instead of 4.8 GB).")]
    public static AlignmentResultDto GlobalAlignAffine(
        [Description("First sequence (DNA, uppercased internally).")] string sequence1,
        [Description("Second sequence.")] string sequence2,
        [Description("Score for matching base.")] int match = 1,
        [Description("Score for mismatching base.")] int mismatch = -1,
        [Description("Gap-open score added once per gap run (negative = penalty). May be fractional when matrix is given; EMBOSS gapopen G / gapextend E = gapOpen -(G - E), gapExtend -E.")] double gapOpen = -2,
        [Description("Gap score per gap position, in addition to gapOpen for the first one. May be fractional (e.g. -0.5) when matrix is given.")] double gapExtend = -1,
        [Description(MatrixParamDescription)] string? matrix = null,
        [Description(PairFormatParamDescription)] string? format = null,
        [Description(NamesParamDescription)] string[]? names = null,
        CancellationToken cancellationToken = default,
        [Description(TrimEndGapsParamDescription)] bool trimEndGaps = false)
    {
        return FormatPair(GlobalAlignAffineCore(sequence1, sequence2, match, mismatch, gapOpen, gapExtend, matrix, cancellationToken), format, sequence1.Length, sequence2.Length, names, sequence1: sequence1, sequence2: sequence2, trimEndGaps: trimEndGaps);
    }

    private static AlignmentResultDto GlobalAlignAffineCore(
        [Description("First sequence (DNA, uppercased internally).")] string sequence1,
        [Description("Second sequence.")] string sequence2,
        [Description("Score for matching base.")] int match = 1,
        [Description("Score for mismatching base.")] int mismatch = -1,
        [Description("Gap-open score added once per gap run (negative = penalty). May be fractional when matrix is given; EMBOSS gapopen G / gapextend E = gapOpen -(G - E), gapExtend -E.")] double gapOpen = -2,
        [Description("Gap score per gap position, in addition to gapOpen for the first one. May be fractional (e.g. -0.5) when matrix is given.")] double gapExtend = -1,
        [Description(MatrixParamDescription)] string? matrix = null,
        CancellationToken cancellationToken = default)
    {
        using var cancelScope = global::Seqeron.Genomics.Alignment.SequenceAligner.BeginCancellationScope(cancellationToken);
        RequireSequence(sequence1, nameof(sequence1));
        RequireSequence(sequence2, nameof(sequence2));

        var sub = ResolveMatrix(matrix);
        if (sub is not null)
            return ToAlignmentDto(global::Seqeron.Genomics.Alignment.SequenceAligner.GlobalAlignAffine(sequence1, sequence2, new global::Seqeron.Genomics.Infrastructure.MatrixScoring(sub, gapOpen, gapExtend)));

        var scoring = new global::Seqeron.Genomics.Infrastructure.ScoringMatrix(match, mismatch, RequireInteger(gapOpen, nameof(gapOpen)), RequireInteger(gapExtend, nameof(gapExtend)));
        var r = global::Seqeron.Genomics.Alignment.SequenceAligner.GlobalAlignAffine(sequence1, sequence2, scoring);
        return ToAlignmentDto(r);
    }

    [McpServerTool(Name = "local_align_affine", Title = "Pairwise — Local Alignment, Affine Gaps (Smith–Waterman–Gotoh)", ReadOnly = true)]
    [Description("Local pairwise alignment with affine gap costs (Smith–Waterman 1981 with Gotoh 1982 gaps): a gap of length k scores gapOpen + k·gapExtend; zero floor. Equals Biopython PairwiseAligner(mode='local') and parasail sw_trace. With matrix = BLOSUM62 | PAM250 | EDNAFULL the alignment is scored by that published residue-pair table (protein or IUPAC DNA, any symbol of the matrix) instead of match/mismatch and gapOpen/gapExtend may be fractional (EMBOSS needle/water 6.6.0 = Biopython PairwiseAligner = parasail on >1000 oracle rows; EMBOSS -gapopen G -gapextend E = gapOpen -(G - E), gapExtend -E); exactScore is then the real-valued score (score = rounded). Memory: inputs above 25,000,000 dynamic-programming cells ((len1+1)*(len2+1), e.g. two 5 kb sequences) are aligned automatically in exact LINEAR SPACE (Hirschberg 1975 / Myers-Miller 1988 / Huang-Miller 1991-style locate-then-solve) - the returned score, rows and coordinates are identical to the full-matrix result (ties included), only the running time grows (about 1.3-2x) and the memory stays O(len1 + len2) instead of O(len1 * len2) (20 kb x 20 kb affine: 111 MB peak instead of 4.8 GB).")]
    public static AlignmentResultDto LocalAlignAffine(
        [Description("First sequence.")] string sequence1,
        [Description("Second sequence.")] string sequence2,
        [Description("Score for matching base.")] int match = 1,
        [Description("Score for mismatching base.")] int mismatch = -1,
        [Description("Gap-open score added once per gap run (negative = penalty). May be fractional when matrix is given; EMBOSS gapopen G / gapextend E = gapOpen -(G - E), gapExtend -E.")] double gapOpen = -2,
        [Description("Gap score per gap position, in addition to gapOpen for the first one. May be fractional (e.g. -0.5) when matrix is given.")] double gapExtend = -1,
        [Description(MatrixParamDescription)] string? matrix = null,
        [Description(PairFormatParamDescription)] string? format = null,
        [Description(NamesParamDescription)] string[]? names = null,
        CancellationToken cancellationToken = default,
        [Description(SignificanceParamDescription)] bool significance = false,
        [Description(DbLengthParamDescription)] long? dbLength = null,
        [Description(DbSequencesParamDescription)] long dbSequences = 1)
    {
        var dto = WithSignificance(LocalAlignAffineCore(sequence1, sequence2, match, mismatch, gapOpen, gapExtend, matrix, cancellationToken), significance, matrix, match, mismatch, gapOpen, gapExtend, sequence1.Length, sequence2.Length, dbLength, dbSequences);
        return FormatPair(dto, format, sequence1.Length, sequence2.Length, names, sequence1: sequence1, sequence2: sequence2);
    }

    private static AlignmentResultDto LocalAlignAffineCore(
        [Description("First sequence.")] string sequence1,
        [Description("Second sequence.")] string sequence2,
        [Description("Score for matching base.")] int match = 1,
        [Description("Score for mismatching base.")] int mismatch = -1,
        [Description("Gap-open score added once per gap run (negative = penalty). May be fractional when matrix is given; EMBOSS gapopen G / gapextend E = gapOpen -(G - E), gapExtend -E.")] double gapOpen = -2,
        [Description("Gap score per gap position, in addition to gapOpen for the first one. May be fractional (e.g. -0.5) when matrix is given.")] double gapExtend = -1,
        [Description(MatrixParamDescription)] string? matrix = null,
        CancellationToken cancellationToken = default)
    {
        using var cancelScope = global::Seqeron.Genomics.Alignment.SequenceAligner.BeginCancellationScope(cancellationToken);
        RequireSequence(sequence1, nameof(sequence1));
        RequireSequence(sequence2, nameof(sequence2));

        var sub = ResolveMatrix(matrix);
        if (sub is not null)
            return ToAlignmentDto(global::Seqeron.Genomics.Alignment.SequenceAligner.LocalAlignAffine(sequence1, sequence2, new global::Seqeron.Genomics.Infrastructure.MatrixScoring(sub, gapOpen, gapExtend)));

        var scoring = new global::Seqeron.Genomics.Infrastructure.ScoringMatrix(match, mismatch, RequireInteger(gapOpen, nameof(gapOpen)), RequireInteger(gapExtend, nameof(gapExtend)));
        var r = global::Seqeron.Genomics.Alignment.SequenceAligner.LocalAlignAffine(sequence1, sequence2, scoring);
        return ToAlignmentDto(r);
    }

    [McpServerTool(Name = "local_align_top_n", Title = "Pairwise — Top-N Non-Intersecting Local Alignments (Waterman–Eggert)", ReadOnly = true)]
    [Description("Waterman–Eggert (1987, J Mol Biol 197:723) declumping: the best local alignment, then the best local alignment that shares no ALIGNED RESIDUE PAIR with it, then the best sharing no pair with either, ... up to maxAlignments alignments of score >= minScore (and always > 0) - the 'alternatives' of EMBOSS matcher -alternatives N (Huang-Miller SIM; exact same definition, verified round by round on 2 902 reported alignments) and FASTA lalign. An aligned pair is a diagonal column (a_i opposite b_j); gap columns are not pairs, so two reported alignments may cross through a gap column but never align the same two residues. Alignment 1 is exactly local_align (affine = false) / local_align_affine (affine = true); scores are non-increasing; ties between co-optimal alignments are broken, with tieOrder = WatermanEggert (default), as local_align does (first maximal cell in row-major order, traceback diagonal > up > left, M > X > Y), so on tie-rich inputs a later alignment can differ from matcher's (its score is always the true maximum of the remaining problem); tieOrder = Sim reproduces EMBOSS matcher's own tie order exactly (the Huang-Miller SIM program of matcher.c ported line by line: K-vertex LIST, forward pass by start point, locate / small-pass recomputation, linear-space Diff; 1 138 oracle rows equal alignment by alignment, 59 more correctly abort) - it needs integer scores (or a matrix with integer entries) and, like the program, depends on maxAlignments (the program's LIST size); the program aborts when fewer than maxAlignments alignments exist, this tool then returns the alignments found. affine = false: linear gaps (every gap position costs gapExtend, gapOpen ignored); affine = true: a gap of length k costs gapOpen + k*gapExtend (EMBOSS -gapopen G -gapextend E = gapOpen -(G - E), gapExtend -E). With matrix = BLOSUM62 | PAM250 | EDNAFULL residues are scored by that table, gaps may be fractional and minScore/exactScore are real scores. Returns alignments (same fields as local_align) in the order found; fewer than maxAlignments when nothing positive >= minScore is left (empty list when no pair scores > 0). Memory: up to 25,000,000 cells ((len1+1)*(len2+1)) the original recomputation scheme (matrices kept, only cells whose value can change are recomputed after each report; O(len1*len2) memory); ABOVE that (e.g. two 5-20 kb sequences) the same alignments - identical rows, scores, coordinates, order and ties - are computed in exact LINEAR SPACE (O(len1 + len2) memory per alignment: one run of the Hirschberg / Myers-Miller / Huang-Miller locate-then-solve local program per reported alignment with the aligned pairs of the earlier ones forbidden, the SIM / EMBOSS matcher scheme), so the running time grows with maxAlignments (about one extra local alignment per reported alignment) while the memory does not (20 kb x 20 kb: tens of MB instead of 1.6 GB); the linear-space program needs non-positive gap scores (a positive gap score on such an input is rejected with a clear error). With tieOrder = Sim (EMBOSS matcher's own tie order) the sequences are also read as the matcher program reads them (A4-3): only letters and '*' are kept; gap characters (- . ~ ?), '#', '+', '_', digits, blanks and every other character are dropped first and positions count the cleaned sequence. With tieOrder = Sim, maxAlignments = 2147483647 means every alignment (as in the Waterman-Eggert order); the matcher candidate list is bounded at 4 194 304 start points (SequenceAligner.TopNSimMaxVertices, ArgumentException).")]
    public static LocalAlignmentListResult LocalAlignTopN(
        [Description("First sequence.")] string sequence1,
        [Description("Second sequence.")] string sequence2,
        [Description("Maximum number of alignments to report (>= 1).")] int maxAlignments = 10,
        [Description("Smallest accepted score (inclusive); alignments scoring <= 0 are never reported. Real number (fractional with matrix).")] double minScore = 0,
        [Description("false = linear gap model (gapExtend per gap position); true = affine (gapOpen + k*gapExtend).")] bool affine = false,
        [Description("Score for matching base (ignored with matrix).")] int match = 1,
        [Description("Score for mismatching base (ignored with matrix).")] int mismatch = -1,
        [Description("Gap-open score (affine = true only; ignored by the linear model). May be fractional when matrix is given.")] double gapOpen = -2,
        [Description("Gap score charged per gap position (<= 0). May be fractional when matrix is given.")] double gapExtend = -1,
        [Description(MatrixParamDescription)] string? matrix = null,
        [Description("Which co-optimal alignment is reported on a tie: WatermanEggert (default; the fixed local_align order) or Sim (EMBOSS matcher's own order, exact; Sim also reads the sequences as the matcher program does: only letters and '*' are kept, gap characters, digits, blanks and every other character are dropped first and positions count the cleaned sequence).")] string tieOrder = "WatermanEggert",
        CancellationToken cancellationToken = default)
    {
        using var cancelScope = global::Seqeron.Genomics.Alignment.SequenceAligner.BeginCancellationScope(cancellationToken);
        RequireSequence(sequence1, nameof(sequence1));
        RequireSequence(sequence2, nameof(sequence2));
        if (!Enum.TryParse<global::Seqeron.Genomics.Alignment.TopNTieOrder>(tieOrder, ignoreCase: true, out var order) || !Enum.IsDefined(order))
            throw new ArgumentException($"tieOrder must be WatermanEggert or Sim, got '{tieOrder}'.", nameof(tieOrder));
        var sub = ResolveMatrix(matrix);
        if (sub is not null)
        {
            var ms = new global::Seqeron.Genomics.Infrastructure.MatrixScoring(sub, gapOpen, gapExtend);
            return ToAlignmentListDto(affine
                ? global::Seqeron.Genomics.Alignment.SequenceAligner.LocalAlignAffineTopN(sequence1, sequence2, ms, maxAlignments, minScore, order)
                : global::Seqeron.Genomics.Alignment.SequenceAligner.LocalAlignTopN(sequence1, sequence2, ms, maxAlignments, minScore, order));
        }

        var scoring = new global::Seqeron.Genomics.Infrastructure.ScoringMatrix(match, mismatch, RequireInteger(gapOpen, nameof(gapOpen)), RequireInteger(gapExtend, nameof(gapExtend)));
        int min = (int)Math.Clamp(Math.Ceiling(minScore), int.MinValue, int.MaxValue);
        return ToAlignmentListDto(affine
            ? global::Seqeron.Genomics.Alignment.SequenceAligner.LocalAlignAffineTopN(sequence1, sequence2, scoring, maxAlignments, min, order)
            : global::Seqeron.Genomics.Alignment.SequenceAligner.LocalAlignTopN(sequence1, sequence2, scoring, maxAlignments, min, order));
    }

    [McpServerTool(Name = "local_align_co_optimal", Title = "Pairwise — All Co-Optimal Local Alignments", ReadOnly = true)]
    [Description("Enumerates ALL distinct local alignments whose score equals the Smith–Waterman optimum (the set Biopython PairwiseAligner(mode='local').align(a, b) iterates; equal on 1 345 oracle pairs), capped at maxAlignments because the number of co-optimal alignments grows exponentially with repeats (tandem repeats, long gap-placement ties). Every alignment ends with an aligned pair, begins right after a zero-score cell, and is not a zero-score extension of a shorter optimal alignment (Biopython rule). Order is deterministic: end cell in row-major order, then diagonal > up > left (affine: M > X > Y) at each branch, so the first alignment equals local_align / local_align_affine. affine = false: linear gaps (gapExtend per position, gapOpen ignored); affine = true: gapOpen + k*gapExtend. With matrix = BLOSUM62 | PAM250 | EDNAFULL residues are scored by that table and gaps may be fractional (exactScore is the real score). Returns alignments (same fields as local_align); empty list when no pair scores > 0. O(len1*len2) memory: inputs above 25,000,000 cells ((len1+1)*(len2+1)) are rejected with a clear error.")]
    public static LocalAlignmentListResult LocalAlignCoOptimal(
        [Description("First sequence.")] string sequence1,
        [Description("Second sequence.")] string sequence2,
        [Description("Cap on the number of alignments returned (>= 1).")] int maxAlignments = 100,
        [Description("false = linear gap model (gapExtend per gap position); true = affine (gapOpen + k*gapExtend).")] bool affine = false,
        [Description("Score for matching base (ignored with matrix).")] int match = 1,
        [Description("Score for mismatching base (ignored with matrix).")] int mismatch = -1,
        [Description("Gap-open score (affine = true only; ignored by the linear model). May be fractional when matrix is given.")] double gapOpen = -2,
        [Description("Gap score charged per gap position (<= 0). May be fractional when matrix is given.")] double gapExtend = -1,
        [Description(MatrixParamDescription)] string? matrix = null,
        CancellationToken cancellationToken = default)
    {
        using var cancelScope = global::Seqeron.Genomics.Alignment.SequenceAligner.BeginCancellationScope(cancellationToken);
        RequireSequence(sequence1, nameof(sequence1));
        RequireSequence(sequence2, nameof(sequence2));
        var sub = ResolveMatrix(matrix);
        if (sub is not null)
        {
            var ms = new global::Seqeron.Genomics.Infrastructure.MatrixScoring(sub, gapOpen, gapExtend);
            return ToAlignmentListDto(affine
                ? global::Seqeron.Genomics.Alignment.SequenceAligner.LocalAlignAffineCoOptimal(sequence1, sequence2, ms, maxAlignments)
                : global::Seqeron.Genomics.Alignment.SequenceAligner.LocalAlignCoOptimal(sequence1, sequence2, ms, maxAlignments));
        }

        var scoring = new global::Seqeron.Genomics.Infrastructure.ScoringMatrix(match, mismatch, RequireInteger(gapOpen, nameof(gapOpen)), RequireInteger(gapExtend, nameof(gapExtend)));
        return ToAlignmentListDto(affine
            ? global::Seqeron.Genomics.Alignment.SequenceAligner.LocalAlignAffineCoOptimal(sequence1, sequence2, scoring, maxAlignments)
            : global::Seqeron.Genomics.Alignment.SequenceAligner.LocalAlignCoOptimal(sequence1, sequence2, scoring, maxAlignments));
    }

    [McpServerTool(Name = "global_align_co_optimal", Title = "Pairwise — All Co-Optimal Global Alignments", ReadOnly = true)]
    [Description("Enumerates ALL distinct co-optimal GLOBAL (Needleman–Wunsch / Gotoh) alignments - the set Biopython PairwiseAligner(mode='global').align(a, b) iterates; equal as a set on 1 920 oracle pairs covering all 16 end-gap policies, linear and affine gaps, scalar and BLOSUM62 scoring with fractional gaps, plus exhaustive enumeration of every column sequence of tiny inputs - capped at maxAlignments because the number of co-optimal alignments grows exponentially with ties (tandem repeats, long gap-placement ties). An alignment is a gapped column sequence; distinct gapped rows are listed once; an alignment is optimal when its score (substitution columns + gaps; leading/trailing gap columns of a FREE end cost 0) is maximal. Order is deterministic: end cell (last row ascending, then last column ascending when end gaps are free), then diagonal > up > left (affine: M > X > Y) at each branch of the traceback, so the first alignment equals global_align / global_align_affine / semi_global_align(_affine) with the same endGaps. affine = false: linear gaps (every gap position costs gapExtend, gapOpen ignored); affine = true: a gap of length k scores gapOpen + k*gapExtend (a gap in one row directly followed by a gap in the other row is two runs). endGaps (default 'None' = global; any of the 16 end-gap policies, as semi_global_align) selects which end gaps are free; startPosition/endPosition are the consumed intervals of each alignment. With matrix = BLOSUM62 | PAM250 | EDNAFULL residues are scored by that table and gaps may be fractional (exactScore is the real score). An empty sequence has the single all-gap alignment. O(len1*len2) memory: inputs above 25,000,000 cells ((len1+1)*(len2+1)) are rejected with a clear error (use global_align, which falls back to linear space, when one optimal alignment is enough).")]
    public static LocalAlignmentListResult GlobalAlignCoOptimal(
        [Description("First sequence.")] string sequence1,
        [Description("Second sequence.")] string sequence2,
        [Description("Cap on the number of alignments returned (>= 1).")] int maxAlignments = 100,
        [Description("false = linear gap model (gapExtend per gap position); true = affine (gapOpen + k*gapExtend).")] bool affine = false,
        [Description("Score for matching base (ignored with matrix).")] int match = 1,
        [Description("Score for mismatching base (ignored with matrix).")] int mismatch = -1,
        [Description("Gap-open score (affine = true only; ignored by the linear model). May be fractional when matrix is given.")] double gapOpen = -2,
        [Description("Gap score charged per gap position (<= 0). May be fractional when matrix is given.")] double gapExtend = -1,
        [Description(MatrixParamDescription)] string? matrix = null,
        [Description("Which end gaps are free; default 'None' = global (all end gaps penalised). " + EndGapsParamDescription)] string? endGaps = null,
        CancellationToken cancellationToken = default)
    {
        using var cancelScope = global::Seqeron.Genomics.Alignment.SequenceAligner.BeginCancellationScope(cancellationToken);
        RequireSequence(sequence1, nameof(sequence1));
        RequireSequence(sequence2, nameof(sequence2));
        var ends = string.IsNullOrWhiteSpace(endGaps) ? global::Seqeron.Genomics.Infrastructure.EndGaps.None : ResolveEndGaps(endGaps);
        var sub = ResolveMatrix(matrix);
        if (sub is not null)
        {
            var ms = new global::Seqeron.Genomics.Infrastructure.MatrixScoring(sub, gapOpen, gapExtend);
            return ToAlignmentListDto(affine
                ? global::Seqeron.Genomics.Alignment.SequenceAligner.SemiGlobalAlignAffineCoOptimal(sequence1, sequence2, ms, ends, maxAlignments)
                : global::Seqeron.Genomics.Alignment.SequenceAligner.SemiGlobalAlignCoOptimal(sequence1, sequence2, ms, ends, maxAlignments));
        }

        var scoring = new global::Seqeron.Genomics.Infrastructure.ScoringMatrix(match, mismatch, RequireInteger(gapOpen, nameof(gapOpen)), RequireInteger(gapExtend, nameof(gapExtend)));
        return ToAlignmentListDto(affine
            ? global::Seqeron.Genomics.Alignment.SequenceAligner.SemiGlobalAlignAffineCoOptimal(sequence1, sequence2, scoring, ends, maxAlignments)
            : global::Seqeron.Genomics.Alignment.SequenceAligner.SemiGlobalAlignCoOptimal(sequence1, sequence2, scoring, ends, maxAlignments));
    }

    [McpServerTool(Name = "semi_global_align_affine", Title = "Pairwise — Semi-Global (Fitting) Alignment, Affine Gaps", ReadOnly = true)]
    [Description("Semi-global (fitting) pairwise alignment with affine gap costs: by default sequence1 is aligned end to end inside sequence2 (free leading/trailing sequence2 residues, endGaps = 'Fitting'); the optional endGaps parameter selects any of the 16 end-gap policies (overlap / suffix-prefix in both orientations, sequence2 inside sequence1, all four ends free = EMBOSS needle's default -endweight N, none = global); a gap of length k scores gapOpen + k·gapExtend, a free end gap 0. Equals Biopython PairwiseAligner(mode='global', end_insertion_score=0 for the default; left/right_insertion/deletion_score = 0 in general) and parasail sg_dx (sg_* variants). With the default endGaps startPosition2/endPosition2 report the fitted interval of sequence2 (0-based inclusive; the free flanks lie outside it; an empty interval has endPosition2 = startPosition2 - 1) and startPosition1/endPosition1 span the whole query; for any endGaps they are the consumed intervals of each sequence. With matrix = BLOSUM62 | PAM250 | EDNAFULL the alignment is scored by that published residue-pair table (protein or IUPAC DNA, any symbol of the matrix) instead of match/mismatch and gapOpen/gapExtend may be fractional (EMBOSS needle/water 6.6.0 = Biopython PairwiseAligner = parasail on >1000 oracle rows; EMBOSS -gapopen G -gapextend E = gapOpen -(G - E), gapExtend -E); exactScore is then the real-valued score (score = rounded). Matrix mode accepts protein sequences (no DNA validation). Memory: inputs above 25,000,000 dynamic-programming cells ((len1+1)*(len2+1), e.g. two 5 kb sequences) are aligned automatically in exact LINEAR SPACE (Hirschberg 1975 / Myers-Miller 1988 / Huang-Miller 1991-style locate-then-solve) - the returned score, rows and coordinates are identical to the full-matrix result (ties included), only the running time grows (about 1.3-2x) and the memory stays O(len1 + len2) instead of O(len1 * len2) (20 kb x 20 kb affine: 111 MB peak instead of 4.8 GB).")]
    public static AlignmentResultDto SemiGlobalAlignAffine(
        [Description("Query sequence (typically shorter; any letters - N, IUPAC, protein - uppercased internally).")] string sequence1,
        [Description("Reference sequence (typically longer; any letters, uppercased internally).")] string sequence2,
        [Description("Score for matching base.")] int match = 1,
        [Description("Score for mismatching base.")] int mismatch = -1,
        [Description("Gap-open score added once per gap run (negative = penalty). May be fractional when matrix is given; EMBOSS gapopen G / gapextend E = gapOpen -(G - E), gapExtend -E.")] double gapOpen = -2,
        [Description("Gap score per gap position, in addition to gapOpen for the first one. May be fractional (e.g. -0.5) when matrix is given.")] double gapExtend = -1,
        [Description(MatrixParamDescription)] string? matrix = null,
        [Description(EndGapsParamDescription)] string? endGaps = null,
        [Description(PairFormatParamDescription)] string? format = null,
        [Description(NamesParamDescription)] string[]? names = null,
        CancellationToken cancellationToken = default,
        [Description(TrimEndGapsParamDescription)] bool trimEndGaps = false)
    {
        return FormatPair(SemiGlobalAlignAffineCore(sequence1, sequence2, match, mismatch, gapOpen, gapExtend, matrix, endGaps, cancellationToken), format, sequence1.Length, sequence2.Length, names, sequence1: sequence1, sequence2: sequence2, queryFirst: true, trimEndGaps: trimEndGaps);
    }

    private static AlignmentResultDto SemiGlobalAlignAffineCore(
        [Description("Query sequence (typically shorter; any letters - N, IUPAC, protein - uppercased internally).")] string sequence1,
        [Description("Reference sequence (typically longer; any letters, uppercased internally).")] string sequence2,
        [Description("Score for matching base.")] int match = 1,
        [Description("Score for mismatching base.")] int mismatch = -1,
        [Description("Gap-open score added once per gap run (negative = penalty). May be fractional when matrix is given; EMBOSS gapopen G / gapextend E = gapOpen -(G - E), gapExtend -E.")] double gapOpen = -2,
        [Description("Gap score per gap position, in addition to gapOpen for the first one. May be fractional (e.g. -0.5) when matrix is given.")] double gapExtend = -1,
        [Description(MatrixParamDescription)] string? matrix = null,
        [Description(EndGapsParamDescription)] string? endGaps = null,
        CancellationToken cancellationToken = default)
    {
        using var cancelScope = global::Seqeron.Genomics.Alignment.SequenceAligner.BeginCancellationScope(cancellationToken);
        var ends = ResolveEndGaps(endGaps);
        var sub = ResolveMatrix(matrix);
        if (sub is not null)
            return ToAlignmentDto(global::Seqeron.Genomics.Alignment.SequenceAligner.SemiGlobalAlignAffine(RequireSequence(sequence1, nameof(sequence1)), RequireSequence(sequence2, nameof(sequence2)),
                new global::Seqeron.Genomics.Infrastructure.MatrixScoring(sub, gapOpen, gapExtend), ends));

        // R14-2: any string (N, IUPAC, protein letters) like global_align / local_align: the string overload of the library, not the strict ACGT DnaSequence
        var scoring = new global::Seqeron.Genomics.Infrastructure.ScoringMatrix(match, mismatch, RequireInteger(gapOpen, nameof(gapOpen)), RequireInteger(gapExtend, nameof(gapExtend)));
        var r = global::Seqeron.Genomics.Alignment.SequenceAligner.SemiGlobalAlignAffine(RequireSequence(sequence1, nameof(sequence1)), RequireSequence(sequence2, nameof(sequence2)), scoring, ends);
        return ToAlignmentDto(r);
    }

    [McpServerTool(Name = "general_gap_align", Title = "Pairwise — General Gap-Cost Alignment (Waterman–Smith–Beyer)", ReadOnly = true)]
    [Description("Pairwise alignment with an ARBITRARY gap-cost function (Waterman–Smith–Beyer 1976): the score of a gap of length k is gapScores[k-1], so concave (logarithmic, square-root), convex (quadratic), capped, step-wise or any other gap model can be used instead of linear/affine. Every maximal gap run is charged once (a run is never split), the model of Biopython PairwiseAligner with a callable gap_score (1.88, global and local: equal on 960 oracle pairs with 8 gap-function families). mode = global (Needleman–Wunsch), local (Smith–Waterman) or semi_global (free end gaps chosen by endGaps; default Fitting). gapScores are scores (<= 0 = penalty) for k = 1 .. max(length1, length2) (at least that many values; extra values are ignored); a value of -Infinity is not expressible in JSON, use a very large penalty. O(m·n·(m+n)) time: practical to a few hundred residues; a table that is linear or affine (gapScores[k-1] = w1 + (k-1)·e) is recognised and aligned with the O(m·n) linear/affine programs (same score). With matrix = BLOSUM62 | PAM250 | EDNAFULL residues are scored by that table instead of match/mismatch (protein or IUPAC DNA); exactScore is the real-valued score (score = rounded). Coordinates as global_align / local_align / semi_global_align. Limit: a non-affine gap function needs O(m*n) memory (33 bytes per cell) and O(m*n*(m+n)) time, so (len1+1)*(len2+1) above 4,000,000 cells (about 2,000 x 2,000) is rejected with a clear error; linear and affine functions are delegated to the linear-space programs and have no limit.")]
    public static AlignmentResultDto GeneralGapAlign(
        [Description("First sequence (any letters - N, IUPAC, protein - uppercased internally).")] string sequence1,
        [Description("Second sequence (any letters, uppercased internally).")] string sequence2,
        [Description("Gap scores: gapScores[k-1] is the score of a gap of length k (negative = penalty), for k = 1 .. max(length1, length2) at least.")] double[] gapScores,
        [Description("Alignment mode: global, local or semi_global.")] string mode = "global",
        [Description("Score for matching base (ignored with matrix).")] int match = 1,
        [Description("Score for mismatching base (ignored with matrix).")] int mismatch = -1,
        [Description(MatrixParamDescription)] string? matrix = null,
        [Description("semi_global only: which end gaps are free, see semi_global_align (default 'Fitting').")] string? endGaps = null,
        [Description(PairFormatParamDescription)] string? format = null,
        [Description(NamesParamDescription)] string[]? names = null,
        CancellationToken cancellationToken = default,
        [Description(TrimEndGapsParamDescription)] bool trimEndGaps = false)
    {
        return FormatPair(GeneralGapAlignCore(sequence1, sequence2, gapScores, mode, match, mismatch, matrix, endGaps, cancellationToken), format, sequence1.Length, sequence2.Length, names, sequence1: sequence1, sequence2: sequence2, queryFirst: mode.Trim().StartsWith("semi", StringComparison.OrdinalIgnoreCase), trimEndGaps: trimEndGaps);
    }

    private static AlignmentResultDto GeneralGapAlignCore(
        [Description("First sequence (any letters - N, IUPAC, protein - uppercased internally).")] string sequence1,
        [Description("Second sequence (any letters, uppercased internally).")] string sequence2,
        [Description("Gap scores: gapScores[k-1] is the score of a gap of length k (negative = penalty), for k = 1 .. max(length1, length2) at least.")] double[] gapScores,
        [Description("Alignment mode: global, local or semi_global.")] string mode = "global",
        [Description("Score for matching base (ignored with matrix).")] int match = 1,
        [Description("Score for mismatching base (ignored with matrix).")] int mismatch = -1,
        [Description(MatrixParamDescription)] string? matrix = null,
        [Description("semi_global only: which end gaps are free, see semi_global_align (default 'Fitting').")] string? endGaps = null,
        CancellationToken cancellationToken = default)
    {
        using var cancelScope = global::Seqeron.Genomics.Alignment.SequenceAligner.BeginCancellationScope(cancellationToken);
        ArgumentNullException.ThrowIfNull(gapScores);
        RequireSequence(sequence1, nameof(sequence1));
        RequireSequence(sequence2, nameof(sequence2));
        int needed = Math.Max(sequence1.Length, sequence2.Length);
        if (gapScores.Length < needed)
            throw new ArgumentException($"gapScores needs at least max(length1, length2) = {needed} values, got {gapScores.Length}.", nameof(gapScores));
        foreach (var g in gapScores)
            if (double.IsNaN(g) || double.IsInfinity(g))
                throw new ArgumentException("gapScores must be finite numbers.", nameof(gapScores));
        Func<int, double> w = k => gapScores[k - 1];

        var ends = ResolveEndGaps(endGaps);
        var kind = (mode ?? "global").Trim().ToLowerInvariant().Replace("-", "_");
        if (kind is not ("global" or "local" or "semi_global" or "semiglobal"))
            throw new ArgumentException($"Unknown mode '{mode}'. Use global, local or semi_global.", nameof(mode));

        var sub = ResolveMatrix(matrix);
        if (sub is not null)
            return ToAlignmentDto(kind switch
            {
                "global" => global::Seqeron.Genomics.Alignment.SequenceAligner.GlobalAlignGeneralGap(sequence1, sequence2, w, sub),
                "local" => global::Seqeron.Genomics.Alignment.SequenceAligner.LocalAlignGeneralGap(sequence1, sequence2, w, sub),
                _ => global::Seqeron.Genomics.Alignment.SequenceAligner.SemiGlobalAlignGeneralGap(sequence1, sequence2, w, sub, ends),
            });

        // R14-2: any string (N, IUPAC, protein letters), like global_align / local_align
        var scoring = new global::Seqeron.Genomics.Infrastructure.ScoringMatrix(match, mismatch, 0, 0);
        return ToAlignmentDto(kind switch
        {
            "global" => global::Seqeron.Genomics.Alignment.SequenceAligner.GlobalAlignGeneralGap(sequence1, sequence2, w, scoring),
            "local" => global::Seqeron.Genomics.Alignment.SequenceAligner.LocalAlignGeneralGap(sequence1, sequence2, w, scoring),
            _ => global::Seqeron.Genomics.Alignment.SequenceAligner.SemiGlobalAlignGeneralGap(sequence1, sequence2, w, scoring, ends),
        });
    }

    [McpServerTool(Name = "multiple_align_muscle", Title = "MSA — MUSCLE 3.8 (port of the program)", ReadOnly = true)]
    [Description("MUSCLE 3.8.1551's multiple alignment (Edgar 2004, NAR 32:1792; muscle -in x.fa -out y.fa), ported step by step from the program's source: letters upper-cased, gaps / white space dropped, letters outside the alphabet replaced by the wildcard (X protein, N nucleotide), a leading M of a protein aligned as X and restored (methionine hack); first guide tree from a k-mer distance (distance1: kmer6_6 over the six MAFFT residue groups for proteins, kmer4_6 for nucleotides, kmer20_3, kmer20_4, kbit20_3, pwkimura, pwscoredist) clustered by UPGMB (0.1 x average + 0.9 x minimum linkage; cluster1: upgma, upgmamin, upgmamax, neighborjoining) and rooted (root1 pseudo / midlongestspan / minavgleafdist), ClustalW tree weights (weight1: none, henikoff, henikoffpb, gsc, clustalw, threeway), profile-profile dynamic programming with position-specific gap open / close penalties and half-open terminal gaps (termGaps half / full / ext), the log-expectation score (VTML200 likelihood ratios, center -0.52, gap open -2.9) for proteins or MUSCLE's nucleotide sum-of-pairs score (BLASTZ matrix, gap open -400) for DNA / RNA (scoreModel sp = PAM200, sv = VTML240 for proteins, a matrix file replaces the matrix of the SP score); maxIters >= 2 rebuilds the tree from the distances of the draft (distance2 pctidkimura / pctidlog / edit / scoredist, cluster2, root2) and realigns only the subtrees whose leaf sets changed (maxTrees); maxIters >= 3 adds maxIters - 2 rounds of restricted-partition refinement (every guide-tree edge splits the sequences in two, the two profiles are re-aligned and kept when the objective score objScore sp / dp / xp / ps / spf / spm improves, inside the blocks between anchor columns unless anchors = false). The program's own default is maxIters 8 (params.cpp g_uMaxIters = 8; verified with the installed muscle3, whose -help text 'default 16' is stale: a plain run writes the rows of -maxiters 8, not 16); the default here is 2 (draft + tree refinement), pass maxIters = 8 for the program's behaviour. Options accepted but without effect, because the 3.8.1551 source compiles them out: diags / diags1 / diags2 / diagLength / diagBreak / diagMargin, hydro / hydroFactor / tomHydro, gapOpen2 / gapExtend2 / gapAmbig; maxMb / maxHours (resource limits). useTree = a rooted binary Newick tree with edge lengths whose labels are the names (skips the tree refinement); clusterOnly / computeWeights stop early. Rows are returned in input order (treeOrder = the order muscle writes) with the files the program can write (tree1, tree2, distance matrices, weights). totalScore is the natural affine sum of pairs (BLOSUM62 10 / 0.2 for proteins, EDNAFULL 10 / 0.5 for nucleotides), not MUSCLE's objective. VERIFICATION: maxIters 1 and 2 with the option set of the first Docker round (le / sp / sv, gapOpen / gapExtend / center, distance1 kmer6_6, cluster1 / cluster2 upgma*, sueff, maxTrees 0) equal the MUSCLE 3.8.1551 binary built from the same source on 420 / 420 user-run sets (rows and output order). Everything else (maxIters >= 3, objScore, weights, the other distances, neighborjoining and roots, anchors, termGaps, brenner, useTree, the tree / distance / weight files) is verified against the same binary on 630 further user-run sets (committed oracle muscle3_r2_oracle.tsv): rows in the written order, files and the program's refusals equal; the user matrix (matrix) is verified for proteins with a '*' column (the file without it and nucleotides: round 3).")]
    public static MuscleRunResultDto MultipleAlignMuscle(
        [Description("Sequences to align (at least one, none empty): letters, '-', '.' and white space only.")] string[] sequences,
        [Description("'auto' (default; MUSCLE's guess: nucleotide when >= 95 % of the first 100 residues are A, C, G, T, U or N), 'protein', 'dna' or 'rna' (MUSCLE -seqtype).")] string type = "auto",
        [Description("MUSCLE -maxiters (at least 1): 1 = progressive draft, 2 (default) = plus tree refinement, 3 and above = plus restricted-partition refinement iterations (the program's default is 8: pass 8 for a plain 'muscle -in x.fa'; the -help text 'default 16' of 3.8.1551 is stale).")] int maxIters = 2,
        [Description("Profile-profile score (-le / -sp / -sv / -spn): 'default' (LE for proteins, SPN for nucleotides), 'le', 'sp' (PAM200), 'sv' (VTML240), 'spn'. A score of the wrong alphabet is replaced like the program does (SPN for nucleotides, LE for proteins).")] string scoreModel = "default",
        [Description("MUSCLE -gapopen (a negative score; default -2.9 LE, -1439 SP, -300 SV, -400 SPN DNA, -420 RNA).")] double? gapOpen = null,
        [Description("MUSCLE -gapextend (default 0).")] double? gapExtend = null,
        [Description("MUSCLE -center (default -0.52 LE, 0 otherwise, -300 RNA).")] double? center = null,
        [Description("MUSCLE -distance1: 'default' (kmer6_6 protein, kmer4_6 nucleotide), 'kmer6_6', 'kmer4_6' (nucleotides only), 'kmer20_3', 'kmer20_4', 'kbit20_3', 'pwkimura', 'pwscoredist' (protein only). kmer20_3 / kmer20_4 / kbit20_3 need sequences of at least 4 residues (the program reads outside its buffer for shorter ones).")] string distance1 = "default",
        [Description("MUSCLE -cluster1 (first guide tree): 'upgmb' (default), 'upgma', 'upgmamin', 'upgmamax' or 'neighborjoining'.")] string cluster1 = "upgmb",
        [Description("MUSCLE -cluster2 (tree of iteration 2): 'upgmb' (default), 'upgma', 'upgmamin', 'upgmamax' or 'neighborjoining'.")] string cluster2 = "upgmb",
        [Description("MUSCLE -sueff (default 0.1): weight of the average linkage in UPGMB.")] float sueff = 0.1f,
        [Description("MUSCLE -maxtrees (default 1): 0 switches the tree refinement of iteration 2 off; more than 1 repeats it only for a weighting other than ClustalW.")] int maxTrees = 1,
        [Description("true: rows in the order muscle writes them (depth-first leaf order of the final guide tree, -group) instead of the input order.")] bool treeOrder = false,
        [Description("MUSCLE -distance2 (distance of the aligned pairs of iteration 2): 'default' (pctidkimura), 'pctidkimura', 'pctidlog', 'edit', 'scoredist' (protein only).")] string distance2 = "default",
        [Description("MUSCLE -root1: 'pseudo' (default: the root of the clustering), 'midlongestspan' or 'minavgleafdist'.")] string root1 = "pseudo",
        [Description("MUSCLE -root2 (tree of iteration 2): 'pseudo', 'midlongestspan' or 'minavgleafdist'.")] string root2 = "pseudo",
        [Description("MUSCLE -weight1 (sequence weights of the progressive stage): 'clustalw' (default), 'none', 'henikoff', 'henikoffpb', 'gsc', 'threeway'.")] string weight1 = "clustalw",
        [Description("MUSCLE -weight2 (sequence weights of the refinement): same values as weight1.")] string weight2 = "clustalw",
        [Description("MUSCLE -objscore (refinement objective): 'spm' (default: xp up to 100 sequences, spf above), 'sp', 'dp', 'xp', 'ps' (LE score only), 'spf'.")] string objScore = "spm",
        [Description("MUSCLE -termgaps: 'half' (default), 'full', 'ext'.")] string termGaps = "half",
        [Description("MUSCLE -anchors (true, default) / -noanchors (false): refine inside the blocks between anchor columns.")] bool anchors = true,
        [Description("MUSCLE -anchorspacing (default 32).")] int anchorSpacing = 32,
        [Description("MUSCLE -smoothwindow (odd, default 7).")] int smoothWindow = 7,
        [Description("MUSCLE -smoothscoreceil (default 3 LE, 200 SP, 90 SV, 999 SPN).")] double? smoothScoreCeil = null,
        [Description("MUSCLE -minbestcolscore (default 2 LE, 300 SP, 130 SV, 90 SPN).")] double? minBestColScore = null,
        [Description("MUSCLE -minsmoothscore (default 1 LE, 125 SP, 40 SV, 90 SPN).")] double? minSmoothScore = null,
        [Description("MUSCLE -brenner: build the root alignment by aligning every sequence to the root profile.")] bool brenner = false,
        [Description("MUSCLE -diags1 (no effect: compiled out of 3.8.1551; -diags sets diags1 and diags2).")] bool diags1 = false,
        [Description("MUSCLE -diags2 (no effect).")] bool diags2 = false,
        [Description("MUSCLE -diaglength (default 24; no effect).")] int diagLength = 24,
        [Description("MUSCLE -diagbreak (default 1; no effect).")] int diagBreak = 1,
        [Description("MUSCLE -diagmargin (default 5; no effect).")] int diagMargin = 5,
        [Description("MUSCLE -hydro (default 5; no effect).")] int hydro = 5,
        [Description("MUSCLE -hydrofactor (default 1.2; no effect).")] double hydroFactor = 1.2,
        [Description("MUSCLE -tomhydro (no effect).")] bool tomHydro = false,
        [Description("MUSCLE -gapopen2 (no effect).")] double? gapOpen2 = null,
        [Description("MUSCLE -gapextend2 (no effect).")] double? gapExtend2 = null,
        [Description("MUSCLE -gapambig (no effect).")] double? gapAmbig = null,
        [Description("MUSCLE -maxmb (no effect on the result).")] int? maxMb = null,
        [Description("MUSCLE -maxhours (no effect on the result).")] double? maxHours = null,
        [Description("MUSCLE -matrix: the text of a score matrix file (NCBI format) replacing the matrix of the SP score. For DNA/RNA the program reads the file under the nucleotide alphabet (the A/C/G/T entries of the matrix) and then runs with the amino-acid SP score.")] string? matrix = null,
        [Description("MUSCLE -usetree: a rooted binary Newick tree with edge lengths, leaf labels = names; skips the first tree and the tree refinement.")] string? useTree = null,
        [Description("FASTA labels of the sequences (default s0, s1, ...): the leaf names of useTree and of the returned trees.")] string[]? names = null,
        [Description("MUSCLE -clusteronly with -tree1: return only the first guide tree (tree1), no alignment.")] bool clusterOnly = false,
        [Description("MUSCLE -computeweights: stop after the progressive alignment and return the weights file instead of an alignment.")] bool computeWeights = false,
        [Description("MUSCLE -stable: the 3.8.1551 program refuses it (\"-stable not supported in this version of muscle\"), so true is rejected; treeOrder = false (default) already returns the rows in input order.")] bool stable = false,
        CancellationToken cancellationToken = default)
    {
        using var cancelScope = global::Seqeron.Genomics.Alignment.SequenceAligner.BeginCancellationScope(cancellationToken);
        var options = BuildMuscleOptions(type, maxIters, scoreModel, gapOpen, gapExtend, center, distance1, distance2, cluster1, cluster2, root1, root2, sueff, maxTrees, weight1, weight2, objScore, termGaps,
            anchors, anchorSpacing, smoothWindow, smoothScoreCeil, minBestColScore, minSmoothScore) with
        {
            TreeOrder = treeOrder,
            Brenner = brenner,
            Diags1 = diags1,
            Diags2 = diags2,
            DiagLength = diagLength,
            DiagBreak = diagBreak,
            DiagMargin = diagMargin,
            Hydro = hydro,
            HydroFactor = hydroFactor,
            TomHydro = tomHydro,
            GapOpen2 = gapOpen2,
            GapExtend2 = gapExtend2,
            GapAmbig = gapAmbig,
            MaxMb = maxMb,
            MaxHours = maxHours,
            UserMatrix = matrix,
            UserTree = useTree,
            SequenceNames = names,
            ClusterOnly = clusterOnly,
            ComputeWeights = computeWeights,
            Stable = stable,
        };
        var r = WrapNotSupported(() => global::Seqeron.Genomics.Alignment.SequenceAligner.MuscleRun(sequences, options));
        var a = r.Alignment;
        return new MuscleRunResultDto(a.AlignedSequences, a.Consensus, a.TotalScore, a.ExactTotalScore, r.Names.ToArray(), r.OutputOrder.ToArray(), r.Tree1, r.Tree2, r.DistanceMatrix1, r.DistanceMatrix2, r.WeightsFile, r.ReferenceBehaviourUndefined);
    }

    [McpServerTool(Name = "muscle_msa_command", Title = "MSA — MUSCLE 3.8 command lines on alignments (-profile, -refine, -spscore, -ppscore, -maketree, -refinew, -sw)", ReadOnly = true)]
    [Description("The other command lines of MUSCLE 3.8.1551 that work on existing alignments, ported from the source: mode 'profile' (muscle -profile -in1 a -in2 b: profile-profile alignment of two alignments, rows of alignment followed by those of alignment2), 'refine' (-refine: maxIters rounds of restricted-partition refinement of an alignment, the program's default is 8), 'spscore' (-spscore: the sum-of-pairs objective score of an alignment, 'SP=%.4g'), 'ppscore' (-ppscore: the profile-profile score of two alignments of equal length; nucleotides only), 'maketree' (-maketree: the guide tree of an alignment as a Newick file text; the program never parses -distance2 in this mode, the distance is always pctidkimura), 'refinew' (-refinew: every window of refineWindow columns re-aligned from scratch with the full pipeline; nucleotides only), 'sw' (-sw: the local Smith-Waterman alignment of two profiles, residues outside the local alignment in lower case and padded with '.'; nucleotides only; needs distance1 pctidkimura when a profile has two or more sequences). For proteins the program never initialises the score matrix in ppscore / refinew / sw and crashes, so those modes are rejected for proteins; mode 'profdb' (-profdb) is always rejected: the program never initialises the alphabet there and crashes. The sequence weights are weight1 (default ClustalW weights of a guide tree built from the alignment with cluster2 / distance2 / root1 or root2). VERIFICATION STATUS: implemented from the program's source and verified against the MUSCLE 3.8.1551 binary on the second Docker round (committed oracle muscle3_r2_oracle.tsv: rows, trees, scores and refusals equal on every recorded set; protein ppscore / sw / refinew and profdb are refused because the program crashes).")]
    public static MuscleCommandResultDto MuscleMsaCommand(
        [Description("Mode: 'profile', 'refine', 'spscore', 'ppscore', 'maketree', 'refinew', 'sw' or 'profdb' (always rejected, the program crashes).")] string mode,
        [Description("The (first) alignment: equal-length rows, '-' = gap.")] string[] alignment,
        [Description("The second alignment (profile, ppscore, sw).")] string[]? alignment2 = null,
        [Description("'auto' (MUSCLE's guess from the alignment), 'protein', 'dna' or 'rna'.")] string type = "auto",
        [Description("MUSCLE -maxiters for refine / refinew (default 8, the program's).")] int maxIters = 8,
        [Description("'default', 'le', 'sp', 'sv', 'spn' (profile-profile score).")] string scoreModel = "default",
        [Description("MUSCLE -gapopen.")] double? gapOpen = null,
        [Description("MUSCLE -gapextend.")] double? gapExtend = null,
        [Description("MUSCLE -center.")] double? center = null,
        [Description("MUSCLE -distance1 (tree of the sw profiles; also the first tree of refinew).")] string distance1 = "default",
        [Description("MUSCLE -distance2 (guide tree of an alignment): 'default' (pctidkimura), 'pctidlog', 'edit', 'scoredist'.")] string distance2 = "default",
        [Description("MUSCLE -cluster1.")] string cluster1 = "upgmb",
        [Description("MUSCLE -cluster2 (guide tree of an alignment): 'upgmb', 'upgma', 'upgmamin', 'upgmamax', 'neighborjoining'.")] string cluster2 = "upgmb",
        [Description("MUSCLE -root1 (profile / ppscore / sw trees).")] string root1 = "pseudo",
        [Description("MUSCLE -root2 (refine / spscore / maketree tree).")] string root2 = "pseudo",
        [Description("MUSCLE -weight1 (sequence weights): 'clustalw' (default), 'none', 'henikoff', 'henikoffpb', 'gsc', 'threeway'.")] string weight1 = "clustalw",
        [Description("MUSCLE -weight2.")] string weight2 = "clustalw",
        [Description("MUSCLE -objscore (refinement objective).")] string objScore = "spm",
        [Description("MUSCLE -termgaps.")] string termGaps = "half",
        [Description("MUSCLE -anchors (true) / -noanchors (false).")] bool anchors = true,
        [Description("MUSCLE -refinewindow (refinew, default 200).")] int refineWindow = 200,
        [Description("MUSCLE -fromwindow (refinew).")] int fromWindow = 0,
        [Description("MUSCLE -towindow (refinew; 0 = the last window).")] int toWindow = 0,
        [Description("MUSCLE -windowoffset (refinew).")] int windowOffset = 0,
        [Description("FASTA labels of the rows (default s0, s1, ...), the leaf names of the maketree tree.")] string[]? names = null,
        CancellationToken cancellationToken = default)
    {
        using var cancelScope = global::Seqeron.Genomics.Alignment.SequenceAligner.BeginCancellationScope(cancellationToken);
        var options = BuildMuscleOptions(type, maxIters, scoreModel, gapOpen, gapExtend, center, distance1, distance2, cluster1, cluster2, root1, root2, 0.1f, 1, weight1, weight2, objScore, termGaps,
            anchors, 32, 7, null, null, null) with
        {
            RefineWindow = refineWindow,
            FromWindow = fromWindow,
            ToWindow = toWindow,
            WindowOffset = windowOffset,
            SequenceNames = names,
        };
        if (alignment == null || alignment.Length == 0)
            throw new ArgumentException("alignment must contain at least one row.", nameof(alignment));
        string[] Second() => alignment2 is { Length: > 0 } ? alignment2 : throw new ArgumentException("alignment2 is required for this mode.");
        switch ((mode ?? "").Trim().ToLowerInvariant())
        {
            case "profile":
                return new MuscleCommandResultDto(global::Seqeron.Genomics.Alignment.SequenceAligner.MuscleProfileAlign(alignment, Second(), options), null, null);
            case "refine":
                return new MuscleCommandResultDto(global::Seqeron.Genomics.Alignment.SequenceAligner.MuscleRefine(alignment, options), null, null);
            case "spscore":
                return new MuscleCommandResultDto([], global::Seqeron.Genomics.Alignment.SequenceAligner.MuscleScoreAlignment(alignment, options), null);
            case "ppscore":
                return new MuscleCommandResultDto([], WrapNotSupported(() => global::Seqeron.Genomics.Alignment.SequenceAligner.MuscleProfileScore(alignment, Second(), options)), null);
            case "maketree":
                return new MuscleCommandResultDto([], null, global::Seqeron.Genomics.Alignment.SequenceAligner.MuscleGuideTree(alignment, options));
            case "refinew":
                return new MuscleCommandResultDto(WrapNotSupported(() => global::Seqeron.Genomics.Alignment.SequenceAligner.MuscleRefineWindows(alignment, options)), null, null);
            case "sw":
                return new MuscleCommandResultDto(WrapNotSupported(() => global::Seqeron.Genomics.Alignment.SequenceAligner.MuscleLocalAlign(alignment, Second(), options)), null, null);
            case "profdb":
                return new MuscleCommandResultDto(WrapNotSupported(() => global::Seqeron.Genomics.Alignment.SequenceAligner.MuscleProfileDatabase(alignment, Second(), options)), null, null);
            default:
                throw new ArgumentException("mode must be 'profile', 'refine', 'spscore', 'ppscore', 'maketree', 'refinew', 'sw' or 'profdb'.", nameof(mode));
        }
    }

    private static T WrapNotSupported<T>(Func<T> run)
    {
        try
        {
            return run();
        }
        catch (NotSupportedException ex)
        {
            throw new ArgumentException(ex.Message, ex);
        }
    }

    private static global::Seqeron.Genomics.Infrastructure.MuscleOptions BuildMuscleOptions(string? type, int maxIters, string? scoreModel, double? gapOpen, double? gapExtend, double? center,
        string? distance1, string? distance2, string? cluster1, string? cluster2, string? root1, string? root2, float sueff, int maxTrees, string? weight1, string? weight2, string? objScore,
        string? termGaps, bool anchors, int anchorSpacing, int smoothWindow, double? smoothScoreCeil, double? minBestColScore, double? minSmoothScore)
    {
        static T Pick<T>(string? v, string name, string dflt, params (string Key, T Value)[] map)
        {
            string key = (v ?? dflt).Trim().ToLowerInvariant();
            foreach (var (k, value) in map)
                if (k == key)
                    return value;
            throw new ArgumentException($"must be one of: {string.Join(", ", map.Select(m => m.Key))}.", name);
        }

        var cl = new (string, global::Seqeron.Genomics.Infrastructure.MuscleClusterMethod)[]
        {
            ("upgmb", global::Seqeron.Genomics.Infrastructure.MuscleClusterMethod.Upgmb), ("upgma", global::Seqeron.Genomics.Infrastructure.MuscleClusterMethod.Upgma),
            ("upgmamin", global::Seqeron.Genomics.Infrastructure.MuscleClusterMethod.UpgmaMin), ("upgmamax", global::Seqeron.Genomics.Infrastructure.MuscleClusterMethod.UpgmaMax),
            ("neighborjoining", global::Seqeron.Genomics.Infrastructure.MuscleClusterMethod.NeighborJoining),
        };
        var roots = new (string, global::Seqeron.Genomics.Infrastructure.MuscleRoot)[]
        {
            ("pseudo", global::Seqeron.Genomics.Infrastructure.MuscleRoot.Pseudo), ("midlongestspan", global::Seqeron.Genomics.Infrastructure.MuscleRoot.MidLongestSpan),
            ("minavgleafdist", global::Seqeron.Genomics.Infrastructure.MuscleRoot.MinAvgLeafDist),
        };
        var weights = new (string, global::Seqeron.Genomics.Infrastructure.MuscleSeqWeight)[]
        {
            ("clustalw", global::Seqeron.Genomics.Infrastructure.MuscleSeqWeight.ClustalW), ("none", global::Seqeron.Genomics.Infrastructure.MuscleSeqWeight.None),
            ("henikoff", global::Seqeron.Genomics.Infrastructure.MuscleSeqWeight.Henikoff), ("henikoffpb", global::Seqeron.Genomics.Infrastructure.MuscleSeqWeight.HenikoffPb),
            ("gsc", global::Seqeron.Genomics.Infrastructure.MuscleSeqWeight.Gsc), ("threeway", global::Seqeron.Genomics.Infrastructure.MuscleSeqWeight.ThreeWay),
        };
        var distances = new (string, global::Seqeron.Genomics.Infrastructure.MuscleDistance)[]
        {
            ("default", global::Seqeron.Genomics.Infrastructure.MuscleDistance.Default), ("kmer6_6", global::Seqeron.Genomics.Infrastructure.MuscleDistance.Kmer6_6),
            ("kmer4_6", global::Seqeron.Genomics.Infrastructure.MuscleDistance.Kmer4_6), ("kmer20_3", global::Seqeron.Genomics.Infrastructure.MuscleDistance.Kmer20_3),
            ("kmer20_4", global::Seqeron.Genomics.Infrastructure.MuscleDistance.Kmer20_4), ("kbit20_3", global::Seqeron.Genomics.Infrastructure.MuscleDistance.Kbit20_3),
            ("pctidkimura", global::Seqeron.Genomics.Infrastructure.MuscleDistance.PctIdKimura), ("pctid_kimura", global::Seqeron.Genomics.Infrastructure.MuscleDistance.PctIdKimura),
            ("pctidlog", global::Seqeron.Genomics.Infrastructure.MuscleDistance.PctIdLog), ("pctid_log", global::Seqeron.Genomics.Infrastructure.MuscleDistance.PctIdLog),
            ("pwkimura", global::Seqeron.Genomics.Infrastructure.MuscleDistance.PwKimura), ("pwscoredist", global::Seqeron.Genomics.Infrastructure.MuscleDistance.PwScoreDist),
            ("scoredist", global::Seqeron.Genomics.Infrastructure.MuscleDistance.ScoreDist), ("edit", global::Seqeron.Genomics.Infrastructure.MuscleDistance.Edit),
        };
        return new global::Seqeron.Genomics.Infrastructure.MuscleOptions
        {
            SequenceType = Pick(type, nameof(type), "auto", ("auto", global::Seqeron.Genomics.Infrastructure.MuscleSequenceType.Auto), ("protein", global::Seqeron.Genomics.Infrastructure.MuscleSequenceType.Protein),
                ("dna", global::Seqeron.Genomics.Infrastructure.MuscleSequenceType.Dna), ("rna", global::Seqeron.Genomics.Infrastructure.MuscleSequenceType.Rna)),
            MaxIterations = maxIters,
            ScoreModel = Pick(scoreModel, nameof(scoreModel), "default", ("default", global::Seqeron.Genomics.Infrastructure.MuscleScoreModel.Default), ("le", global::Seqeron.Genomics.Infrastructure.MuscleScoreModel.LogExpectation),
                ("sp", global::Seqeron.Genomics.Infrastructure.MuscleScoreModel.SumOfPairs), ("sv", global::Seqeron.Genomics.Infrastructure.MuscleScoreModel.SumOfPairsVtml),
                ("spn", global::Seqeron.Genomics.Infrastructure.MuscleScoreModel.SumOfPairsNucleotide)),
            GapOpen = gapOpen,
            GapExtend = gapExtend,
            Center = center,
            Distance1 = Pick(distance1, nameof(distance1), "default", distances),
            Distance2 = Pick(distance2, nameof(distance2), "default", distances),
            Cluster1 = Pick(cluster1, nameof(cluster1), "upgmb", cl),
            Cluster2 = Pick(cluster2, nameof(cluster2), "upgmb", cl),
            Root1 = Pick(root1, nameof(root1), "pseudo", roots),
            Root2 = Pick(root2, nameof(root2), "pseudo", roots),
            Sueff = sueff,
            MaxTrees = maxTrees,
            Weight1 = Pick(weight1, nameof(weight1), "clustalw", weights),
            Weight2 = Pick(weight2, nameof(weight2), "clustalw", weights),
            ObjectiveScore = Pick(objScore, nameof(objScore), "spm", ("spm", global::Seqeron.Genomics.Infrastructure.MuscleObjectiveScore.Spm), ("sp", global::Seqeron.Genomics.Infrastructure.MuscleObjectiveScore.Sp),
                ("dp", global::Seqeron.Genomics.Infrastructure.MuscleObjectiveScore.Dp), ("xp", global::Seqeron.Genomics.Infrastructure.MuscleObjectiveScore.Xp), ("ps", global::Seqeron.Genomics.Infrastructure.MuscleObjectiveScore.Ps),
                ("spf", global::Seqeron.Genomics.Infrastructure.MuscleObjectiveScore.Spf)),
            TermGaps = Pick(termGaps, nameof(termGaps), "half", ("half", global::Seqeron.Genomics.Infrastructure.MuscleTermGaps.Half), ("full", global::Seqeron.Genomics.Infrastructure.MuscleTermGaps.Full),
                ("ext", global::Seqeron.Genomics.Infrastructure.MuscleTermGaps.Ext)),
            Anchors = anchors,
            AnchorSpacing = anchorSpacing,
            SmoothWindow = smoothWindow,
            SmoothScoreCeil = smoothScoreCeil,
            MinBestColScore = minBestColScore,
            MinSmoothScore = minSmoothScore,
        };
    }

    #endregion

    #region ApproximateMatcher

    [McpServerTool(Name = "find_with_mismatches", Title = "Approximate — Find with Mismatches (Hamming)", ReadOnly = true)]
    [Description("Finds all occurrences of pattern inside sequence allowing up to maxMismatches substitutions (Hamming-style, fixed-length window). Reports each match position, matched substring, mismatch count, and 0-based mismatch positions.")]
    public static ApproximateMatchListResult FindWithMismatches(
        [Description("Sequence to search in.")] string sequence,
        [Description("Pattern to find.")] string pattern,
        [Description("Maximum number of allowed mismatches (>= 0).")] int maxMismatches)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty.", nameof(sequence));
        if (string.IsNullOrEmpty(pattern))
            throw new ArgumentException("Pattern cannot be null or empty.", nameof(pattern));
        if (maxMismatches < 0)
            throw new ArgumentOutOfRangeException(nameof(maxMismatches), "maxMismatches must be >= 0.");

        var items = new List<ApproximateMatchDto>();
        foreach (var m in global::Seqeron.Genomics.Alignment.ApproximateMatcher.FindWithMismatches(
            sequence, pattern, maxMismatches))
        {
            items.Add(ToApproximateMatchDto(m));
        }
        return new ApproximateMatchListResult(items.ToArray());
    }

    [McpServerTool(Name = "find_with_edits", Title = "Approximate — Find with Edits (Levenshtein)", ReadOnly = true)]
    [Description("Finds all approximate matches of pattern in sequence allowing up to maxEdits Levenshtein edits (insertions, deletions, substitutions) with variable-length sliding windows of length pattern.Length ± maxEdits.")]
    public static ApproximateMatchListResult FindWithEdits(
        [Description("Sequence to search in.")] string sequence,
        [Description("Pattern to find.")] string pattern,
        [Description("Maximum allowed edit distance (>= 0).")] int maxEdits)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty.", nameof(sequence));
        if (string.IsNullOrEmpty(pattern))
            throw new ArgumentException("Pattern cannot be null or empty.", nameof(pattern));
        if (maxEdits < 0)
            throw new ArgumentOutOfRangeException(nameof(maxEdits), "maxEdits must be >= 0.");

        var items = new List<ApproximateMatchDto>();
        foreach (var m in global::Seqeron.Genomics.Alignment.ApproximateMatcher.FindWithEdits(
            sequence, pattern, maxEdits))
        {
            items.Add(ToApproximateMatchDto(m));
        }
        return new ApproximateMatchListResult(items.ToArray());
    }

    [McpServerTool(Name = "find_edit_end_positions", Title = "Approximate — Edit End Positions (Sellers)", ReadOnly = true)]
    [Description("Sellers (1980) k-differences search: reports every 0-based end position j in sequence at which some substring ending at j is within maxEdits Levenshtein edits of pattern, with that minimum distance (case-insensitive; Myers bit-parallel engine). Ordered by increasing end position. Optional insertionCost/deletionCost/substitutionCost (default 1) switch to the weighted Sellers DP (pattern = s1, text window = s2, rapidfuzz weights convention); maxEdits is then the maximum weighted cost.")]
    public static EditEndPositionsResult FindEditEndPositions(
        [Description("Sequence to search in.")] string sequence,
        [Description("Pattern to find.")] string pattern,
        [Description("Maximum allowed edit distance (>= 0); the maximum weighted cost when costs are given.")] int maxEdits,
        [Description("Cost of a text character absent from the pattern (insertion into the pattern, >= 0; default 1).")] int insertionCost = 1,
        [Description("Cost of a pattern character absent from the text (deletion from the pattern, >= 0; default 1).")] int deletionCost = 1,
        [Description("Cost of a substitution (>= 0; default 1).")] int substitutionCost = 1)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty.", nameof(sequence));
        if (string.IsNullOrEmpty(pattern))
            throw new ArgumentException("Pattern cannot be null or empty.", nameof(pattern));
        if (maxEdits < 0)
            throw new ArgumentOutOfRangeException(nameof(maxEdits), "maxEdits must be >= 0.");

        var costs = ToEditCosts(insertionCost, deletionCost, substitutionCost);
        var positions = costs == global::Seqeron.Genomics.Alignment.EditCosts.Unit
            ? global::Seqeron.Genomics.Alignment.ApproximateMatcher.FindEditEndPositions(sequence, pattern, maxEdits)
            : global::Seqeron.Genomics.Alignment.ApproximateMatcher.FindEditEndPositions(sequence, pattern, maxEdits, costs);
        var items = positions
            .Select(e => new EditEndPositionItem(e.EndPosition, e.Distance))
            .ToArray();
        return new EditEndPositionsResult(items);
    }

    [McpServerTool(Name = "edit_alignment", Title = "Approximate — Edit (Levenshtein) Alignment", ReadOnly = true)]
    [Description("Optimal global unit-cost (Levenshtein) alignment of query against target in edlib's convention: operations '=' match, 'X' mismatch, 'I' query character absent from the target, 'D' target character absent from the query; returns the distance, extended and standard CIGAR, gapped strings and substitution positions. Case-sensitive. linearSpace=true uses Hirschberg's O(m+n)-space algorithm (same distance, possibly a different co-optimal path). Optional insertionCost/deletionCost/substitutionCost (default 1) give a weighted alignment with rapidfuzz Levenshtein weights semantics (query = s1, target = s2: 'D' costs insertionCost, 'I' deletionCost, 'X' substitutionCost); distance is then the weighted cost.")]
    public static EditAlignmentDto EditAlignment(
        [Description("Query sequence (alignment rows).")] string query,
        [Description("Target sequence (alignment columns).")] string target,
        [Description("Use Hirschberg's linear-space algorithm instead of the full-matrix diagonal-first traceback (default false).")] bool linearSpace = false,
        [Description("Cost of inserting a target character, i.e. a 'D' column (>= 0; default 1).")] int insertionCost = 1,
        [Description("Cost of deleting a query character, i.e. an 'I' column (>= 0; default 1).")] int deletionCost = 1,
        [Description("Cost of a substitution, i.e. an 'X' column (>= 0; default 1).")] int substitutionCost = 1)
    {
        if (string.IsNullOrEmpty(query))
            throw new ArgumentException("Query cannot be null or empty.", nameof(query));
        if (string.IsNullOrEmpty(target))
            throw new ArgumentException("Target cannot be null or empty.", nameof(target));

        var costs = ToEditCosts(insertionCost, deletionCost, substitutionCost);
        bool unit = costs == global::Seqeron.Genomics.Alignment.EditCosts.Unit;
        var a = (linearSpace, unit) switch
        {
            (true, true) => global::Seqeron.Genomics.Alignment.ApproximateMatcher.GetEditAlignmentLinearSpace(query, target),
            (false, true) => global::Seqeron.Genomics.Alignment.ApproximateMatcher.GetEditAlignment(query, target),
            (true, false) => global::Seqeron.Genomics.Alignment.ApproximateMatcher.GetEditAlignmentLinearSpace(query, target, costs),
            (false, false) => global::Seqeron.Genomics.Alignment.ApproximateMatcher.GetEditAlignment(query, target, costs),
        };
        return new EditAlignmentDto(
            a.Distance, a.Operations, a.Cigar, a.StandardCigar,
            a.AlignedQuery, a.AlignedTarget, a.SubstitutionPositions.ToArray(), a.HasIndels);
    }

    private static global::Seqeron.Genomics.Alignment.EditCosts ToEditCosts(int insertionCost, int deletionCost, int substitutionCost)
    {
        if (insertionCost < 0)
            throw new ArgumentOutOfRangeException(nameof(insertionCost), "insertionCost must be >= 0.");
        if (deletionCost < 0)
            throw new ArgumentOutOfRangeException(nameof(deletionCost), "deletionCost must be >= 0.");
        if (substitutionCost < 0)
            throw new ArgumentOutOfRangeException(nameof(substitutionCost), "substitutionCost must be >= 0.");
        return new global::Seqeron.Genomics.Alignment.EditCosts(insertionCost, deletionCost, substitutionCost);
    }

    [McpServerTool(Name = "damerau_levenshtein_distance", Title = "Approximate — Damerau–Levenshtein Distance", ReadOnly = true)]
    [Description("Edit distance with adjacent transpositions. variant 'unrestricted' (default) is the true Damerau–Levenshtein metric (Lowrance & Wagner 1975; DL(CA,ABC)=2); 'osa' is the optimal string alignment (restricted) distance where no substring is edited twice (OSA(CA,ABC)=3). Case-sensitive. Optional insertionCost/deletionCost/substitutionCost/transpositionCost (default 1) give the weighted distance (sequence1 = s1 transformed into sequence2: insertion adds a sequence2 character, deletion removes a sequence1 character); 'unrestricted' requires 2*transpositionCost >= insertionCost + deletionCost (Lowrance-Wagner exactness condition).")]
    public static DamerauDistanceResult DamerauLevenshteinDistance(
        [Description("First sequence.")] string sequence1,
        [Description("Second sequence.")] string sequence2,
        [Description("'unrestricted' (true Damerau–Levenshtein, default) or 'osa' (optimal string alignment).")] string variant = "unrestricted",
        [Description("Cost of inserting a sequence2 character (>= 0; default 1).")] int insertionCost = 1,
        [Description("Cost of deleting a sequence1 character (>= 0; default 1).")] int deletionCost = 1,
        [Description("Cost of a substitution (>= 0; default 1).")] int substitutionCost = 1,
        [Description("Cost of swapping two adjacent characters (>= 0; default 1; 'unrestricted' needs 2*transpositionCost >= insertionCost + deletionCost).")] int transpositionCost = 1)
    {
        if (string.IsNullOrEmpty(sequence1))
            throw new ArgumentException("Sequence cannot be null or empty.", nameof(sequence1));
        if (string.IsNullOrEmpty(sequence2))
            throw new ArgumentException("Sequence cannot be null or empty.", nameof(sequence2));

        string v = ParseDamerauVariant(variant);
        var costs = ToDamerauCosts(insertionCost, deletionCost, substitutionCost, transpositionCost, v);
        bool unit = costs == global::Seqeron.Genomics.Alignment.DamerauCosts.Unit;
        int distance = (v, unit) switch
        {
            ("unrestricted", true) => global::Seqeron.Genomics.Alignment.ApproximateMatcher.DamerauLevenshteinDistance(sequence1, sequence2),
            ("unrestricted", false) => global::Seqeron.Genomics.Alignment.ApproximateMatcher.DamerauLevenshteinDistance(sequence1, sequence2, costs),
            (_, true) => global::Seqeron.Genomics.Alignment.ApproximateMatcher.OptimalStringAlignmentDistance(sequence1, sequence2),
            (_, false) => global::Seqeron.Genomics.Alignment.ApproximateMatcher.OptimalStringAlignmentDistance(sequence1, sequence2, costs),
        };
        return new DamerauDistanceResult(distance, v);
    }

    [McpServerTool(Name = "damerau_alignment", Title = "Approximate — Damerau–Levenshtein Alignment (Edit Script)", ReadOnly = true)]
    [Description("Optimal transposition-aware edit script turning sequence1 into sequence2 (Lowrance & Wagner 1975 trace). variant 'unrestricted' (default, true Damerau–Levenshtein: a transposition block a_k..a_i -> b_l..b_j deletes the characters between the swapped pair and inserts b_(l+1..j-1) between them) or 'osa' (adjacent swaps only). Returns the distance (summed cost), a compact script ('=' match, 'X' substitution, 'I' sequence1 character deleted, 'D' sequence2 character inserted, 'T' transposition followed by one 'i' per character deleted and one 'd' per character inserted inside the block) and the operations with 0-based positions and costs. Optional weighted costs as in damerau_levenshtein_distance. Case-sensitive. O(m*n) time and space.")]
    public static DamerauAlignmentDto DamerauAlignment(
        [Description("Source sequence (s1).")] string sequence1,
        [Description("Target sequence (s2).")] string sequence2,
        [Description("'unrestricted' (true Damerau–Levenshtein, default) or 'osa' (optimal string alignment).")] string variant = "unrestricted",
        [Description("Cost of inserting a sequence2 character (>= 0; default 1).")] int insertionCost = 1,
        [Description("Cost of deleting a sequence1 character (>= 0; default 1).")] int deletionCost = 1,
        [Description("Cost of a substitution (>= 0; default 1).")] int substitutionCost = 1,
        [Description("Cost of swapping two adjacent characters (>= 0; default 1; 'unrestricted' needs 2*transpositionCost >= insertionCost + deletionCost).")] int transpositionCost = 1)
    {
        if (string.IsNullOrEmpty(sequence1))
            throw new ArgumentException("Sequence cannot be null or empty.", nameof(sequence1));
        if (string.IsNullOrEmpty(sequence2))
            throw new ArgumentException("Sequence cannot be null or empty.", nameof(sequence2));

        string v = ParseDamerauVariant(variant);
        var costs = ToDamerauCosts(insertionCost, deletionCost, substitutionCost, transpositionCost, v);
        var a = v == "unrestricted"
            ? global::Seqeron.Genomics.Alignment.ApproximateMatcher.GetDamerauLevenshteinAlignment(sequence1, sequence2, costs)
            : global::Seqeron.Genomics.Alignment.ApproximateMatcher.GetOptimalStringAlignment(sequence1, sequence2, costs);
        var ops = a.Operations
            .Select(o => new DamerauOperationDto(
                o.Kind.ToString().ToLowerInvariant(), o.SourcePosition, o.SourceLength, o.TargetPosition, o.TargetLength, o.Cost))
            .ToArray();
        return new DamerauAlignmentDto(a.Distance, v, a.Script, a.TranspositionCount, ops);
    }

    private static string ParseDamerauVariant(string variant)
    {
        string v = (variant ?? string.Empty).Trim().ToLowerInvariant();
        if (v != "unrestricted" && v != "osa")
            throw new ArgumentException("Variant must be 'unrestricted' or 'osa'.", nameof(variant));
        return v;
    }

    private static global::Seqeron.Genomics.Alignment.DamerauCosts ToDamerauCosts(
        int insertionCost, int deletionCost, int substitutionCost, int transpositionCost, string variant)
    {
        var edit = ToEditCosts(insertionCost, deletionCost, substitutionCost);
        if (transpositionCost < 0)
            throw new ArgumentOutOfRangeException(nameof(transpositionCost), "transpositionCost must be >= 0.");
        if (variant == "unrestricted" && 2L * transpositionCost < (long)insertionCost + deletionCost)
            throw new ArgumentException(
                "variant 'unrestricted' requires 2*transpositionCost >= insertionCost + deletionCost (Lowrance & Wagner 1975); use variant 'osa'.",
                nameof(transpositionCost));
        return new global::Seqeron.Genomics.Alignment.DamerauCosts(edit, transpositionCost);
    }

    [McpServerTool(Name = "find_best_match", Title = "Approximate — Find Best (Minimum Hamming)", ReadOnly = true)]
    [Description("Returns the single best (minimum-Hamming-distance) fixed-length window of pattern inside sequence. Stops early on perfect (distance=0) match.")]
    public static FindBestMatchResult FindBestMatch(
        [Description("Sequence to search in.")] string sequence,
        [Description("Pattern to find.")] string pattern)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty.", nameof(sequence));
        if (string.IsNullOrEmpty(pattern))
            throw new ArgumentException("Pattern cannot be null or empty.", nameof(pattern));

        var m = global::Seqeron.Genomics.Alignment.ApproximateMatcher.FindBestMatch(
            sequence, pattern);
        if (m is null)
            return new FindBestMatchResult(null);
        return new FindBestMatchResult(ToApproximateMatchDto(m.Value));
    }

    [McpServerTool(Name = "frequent_kmers_with_mismatches", Title = "K-mers — Frequent with Mismatches", ReadOnly = true)]
    [Description("Finds the most-frequent k-mers within sequence allowing up to d mismatches (counts each k-mer plus all its DNA neighbors at Hamming distance ≤ d). Returns all k-mers tied for the maximum count.")]
    public static FrequentKmersResult FrequentKmersWithMismatches(
        [Description("Sequence to analyze.")] string sequence,
        [Description("K-mer length (> 0).")] int k,
        [Description("Maximum mismatches in neighborhood (>= 0).")] int d)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty.", nameof(sequence));
        if (k <= 0)
            throw new ArgumentOutOfRangeException(nameof(k), "k must be > 0.");
        if (d < 0)
            throw new ArgumentOutOfRangeException(nameof(d), "d must be >= 0.");

        var items = new List<FrequentKmerItem>();
        foreach (var t in global::Seqeron.Genomics.Alignment.ApproximateMatcher.FindFrequentKmersWithMismatches(
            sequence, k, d))
        {
            items.Add(new FrequentKmerItem(t.Kmer, t.Count));
        }
        return new FrequentKmersResult(items.ToArray());
    }

    [McpServerTool(Name = "frequent_kmers_with_mismatches_and_revcomp", Title = "K-mers — Frequent with Mismatches + Reverse Complements", ReadOnly = true)]
    [Description("Frequent words with mismatches and reverse complements (ROSALIND BA1J): all DNA k-mers P maximising Count_d(sequence, P) + Count_d(sequence, reverseComplement(P)), ties included; the set is closed under reverse complement. Items sorted by k-mer (ordinal).")]
    public static FrequentKmersResult FrequentKmersWithMismatchesAndRevcomp(
        [Description("Sequence to analyze.")] string sequence,
        [Description("K-mer length (> 0).")] int k,
        [Description("Maximum mismatches in neighborhood (>= 0).")] int d)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty.", nameof(sequence));
        if (k <= 0)
            throw new ArgumentOutOfRangeException(nameof(k), "k must be > 0.");
        if (d < 0)
            throw new ArgumentOutOfRangeException(nameof(d), "d must be >= 0.");

        // The library leaves the order unspecified; sort for a deterministic tool output.
        var items = global::Seqeron.Genomics.Alignment.ApproximateMatcher
            .FindFrequentKmersWithMismatchesAndReverseComplements(sequence, k, d)
            .OrderBy(t => t.Kmer, StringComparer.Ordinal)
            .Select(t => new FrequentKmerItem(t.Kmer, t.Count))
            .ToArray();
        return new FrequentKmersResult(items);
    }

    #endregion

    #region SequenceAssembler

    [McpServerTool(Name = "assemble_olc", Title = "Assembly — Overlap-Layout-Consensus", ReadOnly = true)]
    [Description("Assembles input reads into contigs using the overlap–layout–consensus approach: detects pairwise suffix–prefix overlaps above thresholds, builds the overlap graph, lays out contigs, and reports N50 / longest / total length.")]
    public static AssemblyResultDto AssembleOlc(
        [Description("Input sequence reads.")] string[] reads,
        [Description("Minimum overlap length to accept (bp).")] int minOverlap = 20,
        [Description("Minimum identity ratio for an overlap to be accepted (0.0-1.0).")] double minIdentity = 0.9,
        [Description("K-mer size (unused by OLC path; part of shared assembly parameters).")] int kmerSize = 31,
        [Description("Minimum contig length to keep (bp).")] int minContigLength = 100)
    {
        var p = new global::Seqeron.Genomics.Alignment.SequenceAssembler.AssemblyParameters(
            minOverlap, minIdentity, kmerSize, minContigLength);
        var r = global::Seqeron.Genomics.Alignment.SequenceAssembler.AssembleOLC(
            reads ?? Array.Empty<string>(), p);
        return ToAssemblyDto(r);
    }

    [McpServerTool(Name = "assemble_de_bruijn", Title = "Assembly — de Bruijn Graph", ReadOnly = true)]
    [Description("Assembles input reads using a de Bruijn graph: shreds reads into k-mers, builds the (k-1)-mer node graph, traces non-branching paths into contigs, then reports assembly statistics.")]
    public static AssemblyResultDto AssembleDeBruijn(
        [Description("Input sequence reads.")] string[] reads,
        [Description("Minimum overlap length (shared parameter; not the operative knob here).")] int minOverlap = 20,
        [Description("Minimum identity (shared parameter).")] double minIdentity = 0.9,
        [Description("K-mer size used to build the de Bruijn graph.")] int kmerSize = 31,
        [Description("Minimum contig length to keep (bp).")] int minContigLength = 100)
    {
        var p = new global::Seqeron.Genomics.Alignment.SequenceAssembler.AssemblyParameters(
            minOverlap, minIdentity, kmerSize, minContigLength);
        var r = global::Seqeron.Genomics.Alignment.SequenceAssembler.AssembleDeBruijn(
            reads ?? Array.Empty<string>(), p);
        return ToAssemblyDto(r);
    }

    [McpServerTool(Name = "find_all_overlaps", Title = "Assembly — All Pairwise Overlaps", ReadOnly = true)]
    [Description("Computes all suffix-of-i / prefix-of-j overlaps between read pairs above minOverlap length and minIdentity ratio.")]
    public static OverlapListResult FindAllOverlaps(
        [Description("Input sequence reads.")] string[] reads,
        [Description("Minimum overlap length (bp).")] int minOverlap = 20,
        [Description("Minimum identity ratio (0.0-1.0).")] double minIdentity = 0.9)
    {
        var src = global::Seqeron.Genomics.Alignment.SequenceAssembler.FindAllOverlaps(
            reads ?? Array.Empty<string>(), minOverlap, minIdentity);
        var arr = new OverlapDto[src.Count];
        for (int i = 0; i < src.Count; i++)
        {
            var o = src[i];
            arr[i] = new OverlapDto(o.ReadIndex1, o.ReadIndex2, o.OverlapLength, o.Position1, o.Position2);
        }
        return new OverlapListResult(arr);
    }

    [McpServerTool(Name = "find_overlap", Title = "Assembly — Pairwise Suffix–Prefix Overlap", ReadOnly = true)]
    [Description("Detects the longest suffix-of-sequence1 / prefix-of-sequence2 overlap satisfying minOverlap and minIdentity thresholds. Returns null if no qualifying overlap exists.")]
    public static FindOverlapResult FindOverlap(
        [Description("First sequence (suffix candidate).")] string sequence1,
        [Description("Second sequence (prefix candidate).")] string sequence2,
        [Description("Minimum overlap length (bp).")] int minOverlap = 20,
        [Description("Minimum identity ratio (0.0-1.0).")] double minIdentity = 0.9)
    {
        var t = global::Seqeron.Genomics.Alignment.SequenceAssembler.FindOverlap(
            sequence1 ?? string.Empty, sequence2 ?? string.Empty, minOverlap, minIdentity);
        if (t is null)
            return new FindOverlapResult(null);
        var v = t.Value;
        return new FindOverlapResult(new FindOverlapInfo(v.length, v.pos1, v.pos2));
    }

    [McpServerTool(Name = "sequence_identity", Title = "Sequence — Percent Identity (Equal Length)", ReadOnly = true)]
    [Description("Computes percent-identity (fraction in [0,1]) between two equal-length sequences via case-insensitive position-by-position comparison. Returns 0 for unequal-length inputs.")]
    public static IdentityResult SequenceIdentity(
        [Description("First sequence.")] string sequence1,
        [Description("Second sequence (must equal sequence1.Length).")] string sequence2)
    {
        var s1 = sequence1 ?? string.Empty;
        var s2 = sequence2 ?? string.Empty;
        var v = global::Seqeron.Genomics.Alignment.SequenceAssembler.CalculateIdentity(s1, s2);
        return new IdentityResult(v);
    }

    [McpServerTool(Name = "assembly_stats", Title = "Assembly — N50 / Length Stats", ReadOnly = true)]
    [Description("Computes assembly statistics (N50, longest contig, total length, assembled-reads accounting) for a precomputed list of contigs.")]
    public static AssemblyResultDto AssemblyStats(
        [Description("Assembled contig sequences.")] string[] contigs,
        [Description("Total number of input reads (used to fill TotalReads).")] int totalReads)
    {
        var r = global::Seqeron.Genomics.Alignment.SequenceAssembler.CalculateStats(
            contigs ?? Array.Empty<string>(), totalReads);
        return ToAssemblyDto(r);
    }

    [McpServerTool(Name = "merge_contigs", Title = "Assembly — Merge Contigs by Overlap", ReadOnly = true)]
    [Description("Concatenates two contigs, collapsing the specified suffix/prefix overlap of overlapLength bases. Invalid overlapLength falls back to plain concatenation.")]
    public static MergedContigResult MergeContigs(
        [Description("First contig (suffix donor).")] string contig1,
        [Description("Second contig (prefix donor).")] string contig2,
        [Description("Length of overlap to collapse.")] int overlapLength)
    {
        var s = global::Seqeron.Genomics.Alignment.SequenceAssembler.MergeContigs(
            contig1 ?? string.Empty, contig2 ?? string.Empty, overlapLength);
        return new MergedContigResult(s);
    }

    [McpServerTool(Name = "scaffold_contigs", Title = "Assembly — Scaffold by Paired Links", ReadOnly = true)]
    [Description("Builds scaffolds by joining contigs using paired-end link records, inserting gapCharacter between linked contigs to span the indicated gap size.")]
    public static ScaffoldsResult ScaffoldContigs(
        [Description("Contig sequences.")] string[] contigs,
        [Description("Paired-end links between contigs (by index).")] ContigLinkInput[] links,
        [Description("Single-character gap filler (e.g. \"N\").")] string gapCharacter = "N")
    {
        if (string.IsNullOrEmpty(gapCharacter) || gapCharacter.Length != 1)
            throw new ArgumentException("gapCharacter must be exactly one character.", nameof(gapCharacter));

        var srcLinks = new List<(int contig1, int contig2, int gapSize)>(links?.Length ?? 0);
        if (links != null)
        {
            foreach (var l in links)
                srcLinks.Add((l.Contig1, l.Contig2, l.GapSize));
        }

        var scaffolds = global::Seqeron.Genomics.Alignment.SequenceAssembler.Scaffold(
            contigs ?? Array.Empty<string>(), srcLinks, gapCharacter[0]);

        var arr = new string[scaffolds.Count];
        for (int i = 0; i < scaffolds.Count; i++)
            arr[i] = scaffolds[i];
        return new ScaffoldsResult(arr);
    }

    [McpServerTool(Name = "calculate_coverage", Title = "Assembly — Per-Base Coverage Depth", ReadOnly = true)]
    [Description("Maps each read to its best-matching position on the reference (≥ minOverlap matching bases) and returns per-base coverage depth as an integer array of length reference.Length.")]
    public static CoverageResult CalculateCoverage(
        [Description("Reference sequence.")] string reference,
        [Description("Reads to map onto reference.")] string[] reads,
        [Description("Minimum matching bases required to map a read.")] int minOverlap = 20)
    {
        var arr = global::Seqeron.Genomics.Alignment.SequenceAssembler.CalculateCoverage(
            reference ?? string.Empty, reads ?? Array.Empty<string>(), minOverlap);
        return new CoverageResult(arr);
    }

    [McpServerTool(Name = "compute_consensus", Title = "Assembly — Column-Majority Consensus", ReadOnly = true)]
    [Description("Builds a consensus sequence from a list of pre-aligned reads (same length, '-' / 'N' ignored) by majority vote per column. Columns with no informative bases yield 'N'.")]
    public static ConsensusResult ComputeConsensus(
        [Description("Pre-aligned reads of identical length.")] string[] alignedReads)
    {
        if (alignedReads != null && alignedReads.Length > 1)
        {
            int len = alignedReads[0]?.Length ?? 0;
            for (int i = 1; i < alignedReads.Length; i++)
            {
                if ((alignedReads[i]?.Length ?? 0) != len)
                    throw new ArgumentException(
                        $"All aligned reads must have the same length; alignedReads[{i}] differs.",
                        nameof(alignedReads));
            }
        }

        var s = global::Seqeron.Genomics.Alignment.SequenceAssembler.ComputeConsensus(
            alignedReads ?? Array.Empty<string>());
        return new ConsensusResult(s);
    }

    [McpServerTool(Name = "quality_trim_reads", Title = "Reads — Quality Trim (Phred+33)", ReadOnly = true)]
    [Description("Trims Phred+33 quality-encoded reads from both ends, dropping bases whose decoded quality is below minQuality. Reads shorter than minLength after trimming are filtered out.")]
    public static TrimmedReadsResult QualityTrimReads(
        [Description("Reads, each with sequence and matching Phred+33 quality string.")] QualityReadInput[] reads,
        [Description("Minimum acceptable quality score (Phred).")] int minQuality = 20,
        [Description("Minimum length after trimming for a read to be kept.")] int minLength = 50)
    {
        var input = new List<(string sequence, string quality)>(reads?.Length ?? 0);
        if (reads != null)
        {
            for (int i = 0; i < reads.Length; i++)
            {
                var r = reads[i];
                if (r is null)
                    throw new ArgumentException($"reads[{i}] is null.", nameof(reads));
                var seq = r.Sequence ?? string.Empty;
                var qual = r.Quality ?? string.Empty;
                if (seq.Length != qual.Length)
                    throw new ArgumentException(
                        $"reads[{i}]: sequence and quality must be the same length.", nameof(reads));
                input.Add((seq, qual));
            }
        }

        var trimmed = global::Seqeron.Genomics.Alignment.SequenceAssembler.QualityTrimReads(
            input, minQuality, minLength);

        var arr = new string[trimmed.Count];
        for (int i = 0; i < trimmed.Count; i++)
            arr[i] = trimmed[i];
        return new TrimmedReadsResult(arr);
    }

    [McpServerTool(Name = "error_correct_reads", Title = "Reads — k-mer Frequency Error Correction", ReadOnly = true)]
    [Description("Corrects single-base errors in reads using k-mer frequency: any k-mer occurring fewer than minKmerFrequency times is corrected by substituting its middle base if a higher-frequency k-mer can be obtained.")]
    public static CorrectedReadsResult ErrorCorrectReads(
        [Description("Input reads.")] string[] reads,
        [Description("K-mer size used to build the frequency table.")] int kmerSize = 21,
        [Description("Minimum k-mer frequency considered \"trusted\".")] int minKmerFrequency = 3)
    {
        var corrected = global::Seqeron.Genomics.Alignment.SequenceAssembler.ErrorCorrectReads(
            reads ?? Array.Empty<string>(), kmerSize, minKmerFrequency);

        var arr = new string[corrected.Count];
        for (int i = 0; i < corrected.Count; i++)
            arr[i] = corrected[i];
        return new CorrectedReadsResult(arr);
    }

    #endregion

    #region Mappers

    private static AlignmentResultDto ToAlignmentDto(global::Seqeron.Genomics.Infrastructure.AlignmentResult r) =>
        new(
            r.AlignedSequence1,
            r.AlignedSequence2,
            r.Score,
            r.AlignmentType.ToString(),
            r.StartPosition1,
            r.StartPosition2,
            r.EndPosition1,
            r.EndPosition2,
            r.ExactScore);

    private static LocalAlignmentListResult ToAlignmentListDto(IReadOnlyList<global::Seqeron.Genomics.Infrastructure.AlignmentResult> list)
    {
        var items = new AlignmentResultDto[list.Count];
        for (int i = 0; i < items.Length; i++)
            items[i] = ToAlignmentDto(list[i]);
        return new LocalAlignmentListResult(items, items.Length);
    }

    private static AssemblyResultDto ToAssemblyDto(global::Seqeron.Genomics.Alignment.SequenceAssembler.AssemblyResult r)
    {
        var arr = new string[r.Contigs.Count];
        for (int i = 0; i < r.Contigs.Count; i++)
            arr[i] = r.Contigs[i];
        return new AssemblyResultDto(arr, r.TotalReads, r.AssembledReads, r.N50, r.LongestContig, r.TotalLength);
    }

    private static ApproximateMatchDto ToApproximateMatchDto(
        global::Seqeron.Genomics.Alignment.ApproximateMatchResult m)
    {
        var positions = new int[m.MismatchPositions.Count];
        for (int i = 0; i < m.MismatchPositions.Count; i++)
            positions[i] = m.MismatchPositions[i];
        return new ApproximateMatchDto(
            m.Position,
            m.MatchedSequence,
            m.Distance,
            positions,
            m.MismatchType.ToString());
    }

    #endregion


    [McpServerTool(Name = "multiple_align_exact", Title = "MSA — Exact Sum-of-Pairs Optimal Alignment (Carrillo–Lipman)", ReadOnly = true)]
    [Description("EXACT sum-of-pairs optimal multiple alignment: no other alignment of the sequences has a larger totalScore under the chosen gap model. gapModel='linear' (default): totalScore = column sum of pairs (residue pair match/mismatch, residue against gap gapExtend, gap against gap 0; gapOpen ignored). gapModel='affine': natural affine sum of pairs (Altschul 1989: every pair of rows is projected out of the alignment, gap-gap columns skipped, every gap run of the projection pays gapOpen + k*gapExtend) - the objective of the heuristic multiple_align* tools, which are all heuristics for three or more sequences (their totalScore never exceeds this tool's; the star / progressive / iterative aligners reach it on a large share of small related sets). Two sequences = global_align / global_align_affine; three or more = the k-dimensional dynamic program over the 2^k - 1 neighbours of every cell (Altschul & Lipman 1989; for k = 3 the classic 7-neighbour 3-D recurrence) with, for affine gaps, one value per gap-history state of every cell (the exact state-per-pair program, 13 / 75 / 541 states for k = 3 / 4 / 5), pruned by the Carrillo-Lipman (1988) bound: a cell is dropped when the sum over all pairs of the best pairwise alignment through it (forward + backward Gotoh tables) cannot reach the sum-of-pairs score of the best heuristic alignment (star, progressive or iterative); carrilloLipman=false evaluates the whole lattice (identical result, used as the oracle). The problem is NP-complete in the number of sequences: at most 8 sequences (linear) / 6 (affine), and the dynamic program is stopped with a clear error when it would hold more than maxStateCells values (cells x gap states; default 50 000 000, about 200 MB) - similar sequences are cheap (a thin corridor), unrelated ones explode. Verified equal to exhaustive enumeration of every alignment of tiny sets and to two independent Python exact programs on 804 random sets (k = 2..5), both objectives. Returns the alignment rows (input order), majority consensus, the optimal totalScore, lowerBound (score of the best heuristic seed), latticeCells (product of length+1), visitedCells and stateValues (size actually solved). PROTEIN / MATRIX: type='protein' (or matrix = BLOSUM62 | PAM250 | EDNAFULL) runs the same lattice program with the substitution matrix and fractional gap scores (at most 8 / 6 sequences); totalScore is then rounded and exactTotalScore is the real optimum (scoreScale = the integer denominator of lowerBound). Verified on 308 protein / IUPAC sets (k = 2 = Biopython, k = 3, 4 = lattice DP, tiny sets = exhaustive enumeration).")]
    public static ExactMultipleAlignmentResultDto MultipleAlignExact(
        [Description("Sequences to align (at least one): DNA (A, C, G, T) by default, amino acids with type='protein'.")] string[] sequences,
        [Description("Score for matching base.")] int match = 1,
        [Description("Score for mismatching base.")] int mismatch = -1,
        [Description(GapOpenDescription)] double gapOpen = -2,
        [Description(GapExtendDescription)] double gapExtend = -1,
        [Description("'linear' (default; column sum of pairs, gapOpen ignored) or 'affine' (natural affine sum of pairs, gap of length k = gapOpen + k*gapExtend).")] string gapModel = "linear",
        [Description("true (default): prune the dynamic-programming lattice with the Carrillo-Lipman bound; false: evaluate every cell (same result, exponentially slower).")] bool carrilloLipman = true,
        [Description("Complexity guard: largest number of dynamic-programming values (surviving cells x gap states) the exact program may hold; above it the call fails with a clear error instead of exhausting memory.")] long maxStateCells = 50_000_000,
        [Description(MsaTypeDescription)] string type = "dna",
        [Description(MsaMatrixDescription)] string? matrix = null,
        CancellationToken cancellationToken = default)
    {
        using var cancelScope = global::Seqeron.Genomics.Alignment.SequenceAligner.BeginCancellationScope(cancellationToken);
        var (sub, _) = ResolveMsaMatrix(type, matrix);
        var options = new global::Seqeron.Genomics.Infrastructure.MsaOptions { GapModel = ParseExactGapModel(gapModel), ExactCarrilloLipman = carrilloLipman, ExactMaxStateCells = maxStateCells };
        if (sequences == null || sequences.Length == 0)
        {
            var e = global::Seqeron.Genomics.Infrastructure.MultipleAlignmentResult.Empty;
            return new ExactMultipleAlignmentResultDto(e.AlignedSequences, e.Consensus, 0, 0, 1, 0, 0);
        }

        if (sub is not null)
        {
            var m = global::Seqeron.Genomics.Alignment.SequenceAligner.MultipleAlignExactWithStatistics(
                sequences, new global::Seqeron.Genomics.Infrastructure.MatrixScoring(sub, gapOpen, gapExtend), options);
            return new ExactMultipleAlignmentResultDto(m.Alignment.AlignedSequences, m.Alignment.Consensus, m.Alignment.TotalScore, m.LowerBound, m.LatticeCells, m.VisitedCells, m.StateValues,
                m.Alignment.ExactTotalScore, m.ScoreScale);
        }

        ArgumentNullException.ThrowIfNull(sequences);
        var dnaList = ParseDna(sequences);
        var scoring = new global::Seqeron.Genomics.Infrastructure.ScoringMatrix(match, mismatch, RequireInteger(gapOpen, nameof(gapOpen)), RequireInteger(gapExtend, nameof(gapExtend)));
        var r = global::Seqeron.Genomics.Alignment.SequenceAligner.MultipleAlignExactWithStatistics(dnaList, scoring, options);
        return new ExactMultipleAlignmentResultDto(r.Alignment.AlignedSequences, r.Alignment.Consensus, r.Alignment.TotalScore, r.LowerBound, r.LatticeCells, r.VisitedCells, r.StateValues,
            r.Alignment.ExactTotalScore, 1);
    }

    private static global::Seqeron.Genomics.Infrastructure.MsaGapModel ParseExactGapModel(string? gapModel) => (gapModel ?? "linear").Trim().ToLowerInvariant() switch
    {
        "linear" => global::Seqeron.Genomics.Infrastructure.MsaGapModel.Linear,
        "affine" => global::Seqeron.Genomics.Infrastructure.MsaGapModel.Affine,
        _ => throw new ArgumentException("gapModel must be 'linear' or 'affine'.", nameof(gapModel)),
    };

    public static List<global::Seqeron.Genomics.Core.DnaSequence> ParseDna(string[] sequences)
    {
        ArgumentNullException.ThrowIfNull(sequences);
        var dnaList = new List<global::Seqeron.Genomics.Core.DnaSequence>(sequences.Length);
        for (int i = 0; i < sequences.Length; i++)
        {
            if (!global::Seqeron.Genomics.Core.DnaSequence.TryCreate(sequences[i], out var d))
                throw new ArgumentException($"Invalid DNA at sequences[{i}].", nameof(sequences));
            dnaList.Add(d!);
        }

        return dnaList;
    }

    /// <summary>
    /// Scoring of the multiple-alignment tools: <c>null</c> = the scalar match/mismatch model on <c>DnaSequence</c>; otherwise the
    /// substitution matrix (type='protein' defaults to BLOSUM62; type='dna' accepts only EDNAFULL).
    /// </summary>
    public static (global::Seqeron.Genomics.Infrastructure.SubstitutionMatrix? Matrix, bool Protein) ResolveMsaMatrix(string? type, string? matrix)
    {
        var kind = (type ?? "dna").Trim().ToLowerInvariant();
        if (kind is not ("dna" or "protein"))
            throw new ArgumentException("type must be 'dna' or 'protein'.", nameof(type));

        var sub = ResolveMatrix(matrix);
        bool edna = sub is not null && sub.IsNucleotide;
        if (kind == "protein")
        {
            if (edna)
                throw new ArgumentException("matrix " + sub!.Name + " is a nucleotide matrix; use type='dna' with it, or a protein matrix (BLOSUM62, PAM250, ...) for type='protein'.", nameof(matrix));
            return (sub ?? global::Seqeron.Genomics.Infrastructure.SubstitutionMatrix.Blosum62, true);
        }

        if (sub is not null && !edna)
            throw new ArgumentException("type='dna' accepts only a nucleotide matrix (EDNAFULL, EDNAMAT, EDNASIMPLE, ENUC.4.4); use type='protein' for a protein matrix (BLOSUM62, PAM250, ...).", nameof(matrix));
        return (sub, false);
    }

    [McpServerTool(Name = "multiple_align", Title = "MSA — Center-Star Alignment (Gusfield)", ReadOnly = true)]
    [Description("Center-star multiple sequence alignment (Gusfield 1993) with exact pairwise alignments: the center is the sequence with the highest total optimal pairwise global score to the others, every other sequence is aligned to it with the canonical global aligner (Needleman–Wunsch, or Gotoh affine with gapModel='affine') and the pairwise alignments are merged column-wise (once a gap, always a gap). Default linear gap model (every gap position costs gapExtend; gapOpen is ignored); gapModel='affine' charges gapOpen + k·gapExtend per gap run and reports the natural affine sum-of-pairs score. For two sequences the result is an optimal pairwise global alignment under either model. Returns the aligned rows, majority consensus and sum-of-pairs score (gap-gap columns score 0). PROTEIN / MATRIX: type='protein' (or matrix = BLOSUM62 | PAM250 | EDNAFULL) scores residue pairs with the published substitution table and real-valued (fractional) gap scores through the same code (string / ProteinSequence overloads of SequenceAligner.MultipleAlign); the result carries exactTotalScore. Equal row for row to an independent Python implementation on 238 protein / IUPAC sets (linear and affine) and, for two sequences, to the Biopython 1.88 optimum.")]
    public static MultipleAlignmentResultDto MultipleAlign(
        [Description("Sequences to align (at least one): DNA (A, C, G, T) by default, amino acids with type='protein'.")] string[] sequences,
        [Description("Score for matching base.")] int match = 1,
        [Description("Score for mismatching base.")] int mismatch = -1,
        [Description(GapOpenDescription)] double gapOpen = -2,
        [Description(GapExtendDescription)] double gapExtend = -1,
        [Description(GapModelDescription)] string gapModel = "linear",
        [Description(MsaTypeDescription)] string type = "dna",
        [Description(MsaMatrixDescription)] string? matrix = null,
        [Description(FormatParamDescription)] string? format = null,
        [Description(NamesParamDescription)] string[]? names = null,
        CancellationToken cancellationToken = default)
    {
        return FormatMultiple(MultipleAlignCore(sequences, match, mismatch, gapOpen, gapExtend, gapModel, type, matrix, cancellationToken), format, names);
    }

    [McpServerTool(Name = "multiple_align_anchored", Title = "MSA — Anchored Center-Star (approximate, suffix-tree seeds)", ReadOnly = true)]
    [Description("Approximate, faster center-star MSA: picks the center by 4-mer count cosine similarity, builds a suffix tree on it and aligns every other sequence with exact-match anchors (MUMmer-style seed, weighted colinear chain, Needleman–Wunsch — or Gotoh affine with gapModel='affine' — on the gaps between anchors). The pairwise alignments are NOT guaranteed optimal, so the sum-of-pairs score can be lower than multiple_align's. Default linear gap model (gapOpen ignored); gapModel='affine' charges gapOpen + k·gapExtend per gap run and reports the natural affine sum-of-pairs score. PROTEIN / MATRIX: type='protein' (or matrix) scores anchors as the sum of matrix(c, c) and fills the gaps between anchors with the matrix aligner (any alphabet); rows equal an independent Python implementation (anchors = suffix-tree matching statistics, colinear chain, Needleman-Wunsch / Gotoh between anchors).")]
    public static MultipleAlignmentResultDto MultipleAlignAnchored(
        [Description("Sequences to align (at least one): DNA (A, C, G, T) by default, amino acids with type='protein'.")] string[] sequences,
        [Description("Score for matching base.")] int match = 1,
        [Description("Score for mismatching base.")] int mismatch = -1,
        [Description(GapOpenDescription)] double gapOpen = -2,
        [Description(GapExtendDescription)] double gapExtend = -1,
        [Description(GapModelDescription)] string gapModel = "linear",
        [Description(MsaTypeDescription)] string type = "dna",
        [Description(MsaMatrixDescription)] string? matrix = null,
        [Description(FormatParamDescription)] string? format = null,
        [Description(NamesParamDescription)] string[]? names = null,
        CancellationToken cancellationToken = default)
    {
        return FormatMultiple(MultipleAlignAnchoredCore(sequences, match, mismatch, gapOpen, gapExtend, gapModel, type, matrix, cancellationToken), format, names);
    }

    public static MultipleAlignmentResultDto FormatMultiple(MultipleAlignmentResultDto d, string? format, string[]? names)
    {
        if (string.IsNullOrWhiteSpace(format))
            return d;
        double score = d.ExactTotalScore != 0 ? d.ExactTotalScore : d.TotalScore;
        return d with { Formatted = FormatRowsForOutput(format, d.AlignedSequences, names, score) };
    }

    public static MultipleAlignmentResultDto MultipleAlignCore(
        [Description("Sequences to align (at least one): DNA (A, C, G, T) by default, amino acids with type='protein'.")] string[] sequences,
        [Description("Score for matching base.")] int match = 1,
        [Description("Score for mismatching base.")] int mismatch = -1,
        [Description(GapOpenDescription)] double gapOpen = -2,
        [Description(GapExtendDescription)] double gapExtend = -1,
        [Description(GapModelDescription)] string gapModel = "linear",
        [Description(MsaTypeDescription)] string type = "dna",
        [Description(MsaMatrixDescription)] string? matrix = null,
        CancellationToken cancellationToken = default)
    {
        using var cancelScope = global::Seqeron.Genomics.Alignment.SequenceAligner.BeginCancellationScope(cancellationToken);
        return RunMsa(sequences, match, mismatch, gapOpen, gapExtend, type, matrix,
            (seqs, scoring) => global::Seqeron.Genomics.Alignment.SequenceAligner.MultipleAlign(seqs, scoring, ParseMsaOptions(null, gapModel)),
            (seqs, scoring) => global::Seqeron.Genomics.Alignment.SequenceAligner.MultipleAlign(seqs, scoring, ParseMsaOptions(null, gapModel)));
    }

    public static MultipleAlignmentResultDto MultipleAlignAnchoredCore(
        [Description("Sequences to align (at least one): DNA (A, C, G, T) by default, amino acids with type='protein'.")] string[] sequences,
        [Description("Score for matching base.")] int match = 1,
        [Description("Score for mismatching base.")] int mismatch = -1,
        [Description(GapOpenDescription)] double gapOpen = -2,
        [Description(GapExtendDescription)] double gapExtend = -1,
        [Description(GapModelDescription)] string gapModel = "linear",
        [Description(MsaTypeDescription)] string type = "dna",
        [Description(MsaMatrixDescription)] string? matrix = null,
        CancellationToken cancellationToken = default)
    {
        using var cancelScope = global::Seqeron.Genomics.Alignment.SequenceAligner.BeginCancellationScope(cancellationToken);
        return RunMsa(sequences, match, mismatch, gapOpen, gapExtend, type, matrix,
            (seqs, scoring) => global::Seqeron.Genomics.Alignment.SequenceAligner.MultipleAlignAnchored(seqs, scoring, ParseMsaOptions(null, gapModel)),
            (seqs, scoring) => global::Seqeron.Genomics.Alignment.SequenceAligner.MultipleAlignAnchored(seqs, scoring, ParseMsaOptions(null, gapModel)));
    }

    public static MultipleAlignmentResultDto RunMsa(
        string[] sequences, int match, int mismatch, int gapOpen, int gapExtend,
        Func<IEnumerable<global::Seqeron.Genomics.Core.DnaSequence>, global::Seqeron.Genomics.Infrastructure.ScoringMatrix?, global::Seqeron.Genomics.Infrastructure.MultipleAlignmentResult> align) =>
        RunMsa(sequences, match, mismatch, gapOpen, gapExtend, "dna", null, align, null, programPort: true);

    public static MultipleAlignmentResultDto RunMsa(
        string[] sequences, int match, int mismatch, double gapOpen, double gapExtend, string? type, string? matrix,
        Func<IEnumerable<global::Seqeron.Genomics.Core.DnaSequence>, global::Seqeron.Genomics.Infrastructure.ScoringMatrix?, global::Seqeron.Genomics.Infrastructure.MultipleAlignmentResult> align,
        Func<IEnumerable<string>, global::Seqeron.Genomics.Infrastructure.MatrixScoring, global::Seqeron.Genomics.Infrastructure.MultipleAlignmentResult>? alignMatrix,
        bool programPort = false)
    {
        var (sub, _) = ResolveMsaMatrix(type, matrix);
        // an empty list: the Seqeron aligners return an empty result, the program ports refuse it as their programs do (audit A3-4: the port throws)
        if (!programPort && (sequences == null || sequences.Length == 0))
        {
            var e = global::Seqeron.Genomics.Infrastructure.MultipleAlignmentResult.Empty;
            return new MultipleAlignmentResultDto(e.AlignedSequences, e.Consensus, e.TotalScore, e.ExactTotalScore);
        }

        if (sub is not null)
        {
            var m = alignMatrix!(sequences, new global::Seqeron.Genomics.Infrastructure.MatrixScoring(sub, gapOpen, gapExtend));
            return new MultipleAlignmentResultDto(m.AlignedSequences, m.Consensus, m.TotalScore, m.ExactTotalScore);
        }

        var dnaList = ParseDna(sequences);
        var scoring = new global::Seqeron.Genomics.Infrastructure.ScoringMatrix(match, mismatch, RequireInteger(gapOpen, nameof(gapOpen)), RequireInteger(gapExtend, nameof(gapExtend)));
        var r = align(dnaList, scoring);
        return new MultipleAlignmentResultDto(r.AlignedSequences, r.Consensus, r.TotalScore, r.ExactTotalScore);
    }

    public static global::Seqeron.Genomics.Infrastructure.MsaOptions ParseMsaOptions(string? distance, string? gapModel) => new() { GapModel = ParseExactGapModel(gapModel) };
}