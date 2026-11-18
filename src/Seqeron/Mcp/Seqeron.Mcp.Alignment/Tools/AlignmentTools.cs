using System.ComponentModel;
using ModelContextProtocol.Server;

namespace Seqeron.Mcp.Alignment.Tools;

/// <summary>
/// MCP tools for sequence alignment, approximate pattern matching, and
/// genome assembly operations.
/// </summary>
[McpServerToolType]
public class AlignmentTools
{
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
        [Description(MatrixParamDescription)] string? matrix = null)
    {
        if (string.IsNullOrEmpty(sequence1))
            throw new ArgumentException("Sequence cannot be null or empty.", nameof(sequence1));
        if (string.IsNullOrEmpty(sequence2))
            throw new ArgumentException("Sequence cannot be null or empty.", nameof(sequence2));

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
        [Description(MatrixParamDescription)] string? matrix = null)
    {
        if (string.IsNullOrEmpty(sequence1))
            throw new ArgumentException("Sequence cannot be null or empty.", nameof(sequence1));
        if (string.IsNullOrEmpty(sequence2))
            throw new ArgumentException("Sequence cannot be null or empty.", nameof(sequence2));

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
        [Description("Query sequence (typically shorter).")] string sequence1,
        [Description("Reference sequence (typically longer).")] string sequence2,
        [Description("Score for matching base.")] int match = 1,
        [Description("Score for mismatching base.")] int mismatch = -1,
        [Description("Gap-open score (ignored by this linear-gap tool; kept for schema compatibility). May be fractional when matrix is given.")] double gapOpen = -2,
        [Description("Gap score charged per gap position (<= 0). May be fractional (e.g. -0.5) when matrix is given.")] double gapExtend = -1,
        [Description(MatrixParamDescription)] string? matrix = null,
        [Description(EndGapsParamDescription)] string? endGaps = null)
    {
        var ends = ResolveEndGaps(endGaps);
        var sub = ResolveMatrix(matrix);
        if (sub is not null)
            return ToAlignmentDto(global::Seqeron.Genomics.Alignment.SequenceAligner.SemiGlobalAlign(RequireSequence(sequence1, nameof(sequence1)), RequireSequence(sequence2, nameof(sequence2)),
                new global::Seqeron.Genomics.Infrastructure.MatrixScoring(sub, gapOpen, gapExtend), ends));

        if (!global::Seqeron.Genomics.Core.DnaSequence.TryCreate(sequence1, out var d1))
            throw new ArgumentException("Invalid DNA in sequence1.", nameof(sequence1));
        if (!global::Seqeron.Genomics.Core.DnaSequence.TryCreate(sequence2, out var d2))
            throw new ArgumentException("Invalid DNA in sequence2.", nameof(sequence2));

        var scoring = new global::Seqeron.Genomics.Infrastructure.ScoringMatrix(match, mismatch, RequireInteger(gapOpen, nameof(gapOpen)), RequireInteger(gapExtend, nameof(gapExtend)));
        var r = global::Seqeron.Genomics.Alignment.SequenceAligner.SemiGlobalAlign(d1!, d2!, scoring, ends);
        return ToAlignmentDto(r);
    }

    [McpServerTool(Name = "alignment_statistics", Title = "Alignment — Statistics", ReadOnly = true)]
    [Description("Computes identity / similarity / gap statistics for a pairwise alignment (EMBOSS needle convention: every percentage is over the full alignment length, gap columns included; a column is a gap column when either row has '-'; residues compared case-insensitively). Similarity counts columns with a POSITIVE substitution score (BLAST positives): with matrix = BLOSUM62 | PAM250 | EDNAFULL the published residue-pair table is used, otherwise the scalar match/mismatch model (similarity = identity for match > 0 >= mismatch). Also returns positives and gapOpens (number of gap runs, BLAST gapopen) and the BLAST+ report values: blastPositives (identities + non-identical positives), the INTEGER percentages of the BLAST pairwise report (blastIdentityPercent, blastPositivesPercent, blastGapPercent: round half up, never 100 unless every column matches; NCBI GetPercentMatch), blastPident (tabular pident, 3 decimals) and blastReport (' Identities = a/L (p%), [Positives = b/L (q%),] Gaps = g/L (r%)'; Positives only with a protein matrix).")]
    public static AlignmentStatisticsDto AlignmentStatistics(
        [Description("Aligned representation of sequence1 (with '-' gaps).")] string alignedSequence1,
        [Description("Aligned representation of sequence2 (with '-' gaps).")] string alignedSequence2,
        [Description("Optional substitution matrix deciding which non-gap columns are similar: BLOSUM62, PAM250 (protein) or EDNAFULL (IUPAC DNA). Omit to use match/mismatch.")] string? matrix = null,
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
            s.BlastTabularPidentText, s.ToBlastString(protein: sub is not null && !ReferenceEquals(sub, global::Seqeron.Genomics.Infrastructure.SubstitutionMatrix.EdnaFull)));
    }

    [McpServerTool(Name = "format_alignment", Title = "Alignment — Pretty Format", ReadOnly = true)]
    [Description("Renders a pairwise alignment as three-line blocks in the EMBOSS needle/water (srspair) style: the markup line has '|' for identical columns, ':' for similar columns (positive substitution score), '.' for mismatches and ' ' for gap columns; wrapped to a configurable line width. matrix (BLOSUM62 | PAM250 | EDNAFULL) or the scalar match/mismatch decide which mismatches are similar.")]
    public static FormatAlignmentResult FormatAlignment(
        [Description("Aligned representation of sequence1.")] string alignedSequence1,
        [Description("Aligned representation of sequence2.")] string alignedSequence2,
        [Description("Maximum characters per line in the formatted output (>= 1).")] int lineWidth = 60,
        [Description("Optional substitution matrix deciding which non-gap columns are similar: BLOSUM62, PAM250 or EDNAFULL. Omit to use match/mismatch.")] string? matrix = null,
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
        "Optional substitution matrix: BLOSUM62, PAM250 (protein) or EDNAFULL (IUPAC DNA). When given, match/mismatch are ignored, residues are scored by the published table (any protein or IUPAC symbol of the matrix; case-insensitive) and gapOpen/gapExtend may be fractional (the EMBOSS needle/water models: gapOpen = -(G - E), gapExtend = -E for -gapopen G -gapextend E). The result then also carries the real-valued exactScore.";

    private static string RequireSequence(string sequence, string name) =>
        string.IsNullOrEmpty(sequence) ? throw new ArgumentException("Sequence cannot be null or empty.", name) : sequence;

    private static int RequireInteger(double value, string name) =>
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
        return name.Trim().ToUpperInvariant() switch
        {
            "BLOSUM62" or "EBLOSUM62" => global::Seqeron.Genomics.Infrastructure.SubstitutionMatrix.Blosum62,
            "PAM250" or "EPAM250" => global::Seqeron.Genomics.Infrastructure.SubstitutionMatrix.Pam250,
            "EDNAFULL" or "NUC.4.4" or "NUC44" => global::Seqeron.Genomics.Infrastructure.SubstitutionMatrix.EdnaFull,
            _ => throw new ArgumentException($"Unknown matrix '{name}'. Supported: BLOSUM62, PAM250, EDNAFULL.", nameof(name))
        };
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
        [Description(MatrixParamDescription)] string? matrix = null)
    {
        if (string.IsNullOrEmpty(sequence1))
            throw new ArgumentException("Sequence cannot be null or empty.", nameof(sequence1));
        if (string.IsNullOrEmpty(sequence2))
            throw new ArgumentException("Sequence cannot be null or empty.", nameof(sequence2));

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
        [Description(MatrixParamDescription)] string? matrix = null)
    {
        if (string.IsNullOrEmpty(sequence1))
            throw new ArgumentException("Sequence cannot be null or empty.", nameof(sequence1));
        if (string.IsNullOrEmpty(sequence2))
            throw new ArgumentException("Sequence cannot be null or empty.", nameof(sequence2));

        var sub = ResolveMatrix(matrix);
        if (sub is not null)
            return ToAlignmentDto(global::Seqeron.Genomics.Alignment.SequenceAligner.LocalAlignAffine(sequence1, sequence2, new global::Seqeron.Genomics.Infrastructure.MatrixScoring(sub, gapOpen, gapExtend)));

        var scoring = new global::Seqeron.Genomics.Infrastructure.ScoringMatrix(match, mismatch, RequireInteger(gapOpen, nameof(gapOpen)), RequireInteger(gapExtend, nameof(gapExtend)));
        var r = global::Seqeron.Genomics.Alignment.SequenceAligner.LocalAlignAffine(sequence1, sequence2, scoring);
        return ToAlignmentDto(r);
    }

    [McpServerTool(Name = "local_align_top_n", Title = "Pairwise — Top-N Non-Intersecting Local Alignments (Waterman–Eggert)", ReadOnly = true)]
    [Description("Waterman–Eggert (1987, J Mol Biol 197:723) declumping: the best local alignment, then the best local alignment that shares no ALIGNED RESIDUE PAIR with it, then the best sharing no pair with either, ... up to maxAlignments alignments of score >= minScore (and always > 0) - the 'alternatives' of EMBOSS matcher -alternatives N (Huang-Miller SIM; exact same definition, verified round by round on 2 902 reported alignments) and FASTA lalign. An aligned pair is a diagonal column (a_i opposite b_j); gap columns are not pairs, so two reported alignments may cross through a gap column but never align the same two residues. Alignment 1 is exactly local_align (affine = false) / local_align_affine (affine = true); scores are non-increasing; ties between co-optimal alignments are broken as local_align does (first maximal cell in row-major order, traceback diagonal > up > left, M > X > Y), so on tie-rich inputs a later alignment can differ from matcher's (its score is always the true maximum of the remaining problem). affine = false: linear gaps (every gap position costs gapExtend, gapOpen ignored); affine = true: a gap of length k costs gapOpen + k*gapExtend (EMBOSS -gapopen G -gapextend E = gapOpen -(G - E), gapExtend -E). With matrix = BLOSUM62 | PAM250 | EDNAFULL residues are scored by that table, gaps may be fractional and minScore/exactScore are real scores. Returns alignments (same fields as local_align) in the order found; fewer than maxAlignments when nothing positive >= minScore is left (empty list when no pair scores > 0). Implementation: the original recomputation scheme (matrices kept, only cells whose value can change are recomputed after each report); O(len1*len2) memory, inputs above 25,000,000 cells ((len1+1)*(len2+1)) are rejected with a clear error.")]
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
        [Description(MatrixParamDescription)] string? matrix = null)
    {
        RequireSequence(sequence1, nameof(sequence1));
        RequireSequence(sequence2, nameof(sequence2));
        var sub = ResolveMatrix(matrix);
        if (sub is not null)
        {
            var ms = new global::Seqeron.Genomics.Infrastructure.MatrixScoring(sub, gapOpen, gapExtend);
            return ToAlignmentListDto(affine
                ? global::Seqeron.Genomics.Alignment.SequenceAligner.LocalAlignAffineTopN(sequence1, sequence2, ms, maxAlignments, minScore)
                : global::Seqeron.Genomics.Alignment.SequenceAligner.LocalAlignTopN(sequence1, sequence2, ms, maxAlignments, minScore));
        }

        var scoring = new global::Seqeron.Genomics.Infrastructure.ScoringMatrix(match, mismatch, RequireInteger(gapOpen, nameof(gapOpen)), RequireInteger(gapExtend, nameof(gapExtend)));
        int min = (int)Math.Clamp(Math.Ceiling(minScore), int.MinValue, int.MaxValue);
        return ToAlignmentListDto(affine
            ? global::Seqeron.Genomics.Alignment.SequenceAligner.LocalAlignAffineTopN(sequence1, sequence2, scoring, maxAlignments, min)
            : global::Seqeron.Genomics.Alignment.SequenceAligner.LocalAlignTopN(sequence1, sequence2, scoring, maxAlignments, min));
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
        [Description(MatrixParamDescription)] string? matrix = null)
    {
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

    [McpServerTool(Name = "semi_global_align_affine", Title = "Pairwise — Semi-Global (Fitting) Alignment, Affine Gaps", ReadOnly = true)]
    [Description("Semi-global (fitting) pairwise alignment with affine gap costs: by default sequence1 is aligned end to end inside sequence2 (free leading/trailing sequence2 residues, endGaps = 'Fitting'); the optional endGaps parameter selects any of the 16 end-gap policies (overlap / suffix-prefix in both orientations, sequence2 inside sequence1, all four ends free = EMBOSS needle's default -endweight N, none = global); a gap of length k scores gapOpen + k·gapExtend, a free end gap 0. Equals Biopython PairwiseAligner(mode='global', end_insertion_score=0 for the default; left/right_insertion/deletion_score = 0 in general) and parasail sg_dx (sg_* variants). With the default endGaps startPosition2/endPosition2 report the fitted interval of sequence2 (0-based inclusive; the free flanks lie outside it; an empty interval has endPosition2 = startPosition2 - 1) and startPosition1/endPosition1 span the whole query; for any endGaps they are the consumed intervals of each sequence. With matrix = BLOSUM62 | PAM250 | EDNAFULL the alignment is scored by that published residue-pair table (protein or IUPAC DNA, any symbol of the matrix) instead of match/mismatch and gapOpen/gapExtend may be fractional (EMBOSS needle/water 6.6.0 = Biopython PairwiseAligner = parasail on >1000 oracle rows; EMBOSS -gapopen G -gapextend E = gapOpen -(G - E), gapExtend -E); exactScore is then the real-valued score (score = rounded). Matrix mode accepts protein sequences (no DNA validation). Memory: inputs above 25,000,000 dynamic-programming cells ((len1+1)*(len2+1), e.g. two 5 kb sequences) are aligned automatically in exact LINEAR SPACE (Hirschberg 1975 / Myers-Miller 1988 / Huang-Miller 1991-style locate-then-solve) - the returned score, rows and coordinates are identical to the full-matrix result (ties included), only the running time grows (about 1.3-2x) and the memory stays O(len1 + len2) instead of O(len1 * len2) (20 kb x 20 kb affine: 111 MB peak instead of 4.8 GB).")]
    public static AlignmentResultDto SemiGlobalAlignAffine(
        [Description("Query sequence (typically shorter).")] string sequence1,
        [Description("Reference sequence (typically longer).")] string sequence2,
        [Description("Score for matching base.")] int match = 1,
        [Description("Score for mismatching base.")] int mismatch = -1,
        [Description("Gap-open score added once per gap run (negative = penalty). May be fractional when matrix is given; EMBOSS gapopen G / gapextend E = gapOpen -(G - E), gapExtend -E.")] double gapOpen = -2,
        [Description("Gap score per gap position, in addition to gapOpen for the first one. May be fractional (e.g. -0.5) when matrix is given.")] double gapExtend = -1,
        [Description(MatrixParamDescription)] string? matrix = null,
        [Description(EndGapsParamDescription)] string? endGaps = null)
    {
        var ends = ResolveEndGaps(endGaps);
        var sub = ResolveMatrix(matrix);
        if (sub is not null)
            return ToAlignmentDto(global::Seqeron.Genomics.Alignment.SequenceAligner.SemiGlobalAlignAffine(RequireSequence(sequence1, nameof(sequence1)), RequireSequence(sequence2, nameof(sequence2)),
                new global::Seqeron.Genomics.Infrastructure.MatrixScoring(sub, gapOpen, gapExtend), ends));

        if (!global::Seqeron.Genomics.Core.DnaSequence.TryCreate(sequence1, out var d1))
            throw new ArgumentException("Invalid DNA in sequence1.", nameof(sequence1));
        if (!global::Seqeron.Genomics.Core.DnaSequence.TryCreate(sequence2, out var d2))
            throw new ArgumentException("Invalid DNA in sequence2.", nameof(sequence2));

        var scoring = new global::Seqeron.Genomics.Infrastructure.ScoringMatrix(match, mismatch, RequireInteger(gapOpen, nameof(gapOpen)), RequireInteger(gapExtend, nameof(gapExtend)));
        var r = global::Seqeron.Genomics.Alignment.SequenceAligner.SemiGlobalAlignAffine(d1!, d2!, scoring, ends);
        return ToAlignmentDto(r);
    }

    [McpServerTool(Name = "general_gap_align", Title = "Pairwise — General Gap-Cost Alignment (Waterman–Smith–Beyer)", ReadOnly = true)]
    [Description("Pairwise alignment with an ARBITRARY gap-cost function (Waterman–Smith–Beyer 1976): the score of a gap of length k is gapScores[k-1], so concave (logarithmic, square-root), convex (quadratic), capped, step-wise or any other gap model can be used instead of linear/affine. Every maximal gap run is charged once (a run is never split), the model of Biopython PairwiseAligner with a callable gap_score (1.88, global and local: equal on 960 oracle pairs with 8 gap-function families). mode = global (Needleman–Wunsch), local (Smith–Waterman) or semi_global (free end gaps chosen by endGaps; default Fitting). gapScores are scores (<= 0 = penalty) for k = 1 .. max(length1, length2) (at least that many values; extra values are ignored); a value of -Infinity is not expressible in JSON, use a very large penalty. O(m·n·(m+n)) time: practical to a few hundred residues; a table that is linear or affine (gapScores[k-1] = w1 + (k-1)·e) is recognised and aligned with the O(m·n) linear/affine programs (same score). With matrix = BLOSUM62 | PAM250 | EDNAFULL residues are scored by that table instead of match/mismatch (protein or IUPAC DNA); exactScore is the real-valued score (score = rounded). Coordinates as global_align / local_align / semi_global_align. Limit: a non-affine gap function needs O(m*n) memory (33 bytes per cell) and O(m*n*(m+n)) time, so (len1+1)*(len2+1) above 4,000,000 cells (about 2,000 x 2,000) is rejected with a clear error; linear and affine functions are delegated to the linear-space programs and have no limit.")]
    public static AlignmentResultDto GeneralGapAlign(
        [Description("First sequence (DNA unless matrix is given).")] string sequence1,
        [Description("Second sequence (DNA unless matrix is given).")] string sequence2,
        [Description("Gap scores: gapScores[k-1] is the score of a gap of length k (negative = penalty), for k = 1 .. max(length1, length2) at least.")] double[] gapScores,
        [Description("Alignment mode: global, local or semi_global.")] string mode = "global",
        [Description("Score for matching base (ignored with matrix).")] int match = 1,
        [Description("Score for mismatching base (ignored with matrix).")] int mismatch = -1,
        [Description(MatrixParamDescription)] string? matrix = null,
        [Description("semi_global only: which end gaps are free, see semi_global_align (default 'Fitting').")] string? endGaps = null)
    {
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

        if (!global::Seqeron.Genomics.Core.DnaSequence.TryCreate(sequence1, out var d1))
            throw new ArgumentException("Invalid DNA in sequence1.", nameof(sequence1));
        if (!global::Seqeron.Genomics.Core.DnaSequence.TryCreate(sequence2, out var d2))
            throw new ArgumentException("Invalid DNA in sequence2.", nameof(sequence2));

        var scoring = new global::Seqeron.Genomics.Infrastructure.ScoringMatrix(match, mismatch, 0, 0);
        return ToAlignmentDto(kind switch
        {
            "global" => global::Seqeron.Genomics.Alignment.SequenceAligner.GlobalAlignGeneralGap(d1!, d2!, w, scoring),
            "local" => global::Seqeron.Genomics.Alignment.SequenceAligner.LocalAlignGeneralGap(d1!, d2!, w, scoring),
            _ => global::Seqeron.Genomics.Alignment.SequenceAligner.SemiGlobalAlignGeneralGap(d1!, d2!, w, scoring, ends),
        });
    }

    private const string GapOpenDescription =
        "Gap-open penalty. Used only with gapModel='affine' (a gap of length k then costs gapOpen + k * gapExtend, Gotoh); ignored by the default linear model.";

    private const string GapExtendDescription =
        "Gap penalty per gap position (linear model: every gap position; affine model: the per-position part of gapOpen + k * gapExtend).";

    private const string GapModelDescription =
        "'linear' (default): every gap position costs gapExtend and gapOpen is ignored; totalScore is the column sum of pairs (gap-gap columns 0). " +
        "'affine': gap runs cost gapOpen + k * gapExtend; the pairwise / profile alignments use affine (Gotoh) dynamic programming and totalScore is the sum of pairs with " +
        "natural affine gap costs (Altschul 1989: each pair of rows is projected out of the alignment, gap-gap columns skipped, every gap run of the projection pays gapOpen once). " +
        "For two sequences the affine result is an optimal pairwise alignment; for three or more sequences every method is a heuristic.";

    private const string DistanceModeDescription =
        "Guide-tree distance: 'linear' (default) = 1 - identities / min(length) of the repository's global alignment under the given scores (linear NW, or Gotoh when gapModel='affine'); " +
        "'clustalw' = the same formula on ClustalW 2.1's own pairwise stage (affine local Smith-Waterman-Gotoh alignment with its IUB matrix x100 and gap penalties 15 / 6.66, " +
        "the 'slow' full-DP mode of clustalw -ALIGN; the match/mismatch/gap parameters are then used only for the profile alignment).";

    private const string GuideTreeDescription =
        "Guide-tree method: 'upgma' (default; Sokal & Michener 1958, the Feng-Doolittle / CLUSTAL-1 choice) or 'nj' (neighbour-joining, Saitou & Nei 1987 - " +
        "ClustalW 2.1's default tree, T-Coffee/Clustal-style; the unrooted NJ tree is rooted and ordered exactly as clustalw does: Tree::reRoot, left subtree first). " +
        "With distance='clustalw' and guideTree='nj' the tree equals the .dnd file of clustalw -ALIGN. Both are built by the canonical phylogenetics builders.";

    [McpServerTool(Name = "multiple_align", Title = "MSA — Center-Star Alignment (Gusfield)", ReadOnly = true)]
    [Description("Center-star multiple sequence alignment (Gusfield 1993) with exact pairwise alignments: the center is the sequence with the highest total optimal pairwise global score to the others, every other sequence is aligned to it with the canonical global aligner (Needleman–Wunsch, or Gotoh affine with gapModel='affine') and the pairwise alignments are merged column-wise (once a gap, always a gap). Default linear gap model (every gap position costs gapExtend; gapOpen is ignored); gapModel='affine' charges gapOpen + k·gapExtend per gap run and reports the natural affine sum-of-pairs score. For two sequences the result is an optimal pairwise global alignment under either model. Returns the aligned rows, majority consensus and sum-of-pairs score (gap-gap columns score 0).")]
    public static MultipleAlignmentResultDto MultipleAlign(
        [Description("DNA sequences to align (at least one).")] string[] sequences,
        [Description("Score for matching base.")] int match = 1,
        [Description("Score for mismatching base.")] int mismatch = -1,
        [Description(GapOpenDescription)] int gapOpen = -2,
        [Description(GapExtendDescription)] int gapExtend = -1,
        [Description(GapModelDescription)] string gapModel = "linear") =>
        RunMsa(sequences, match, mismatch, gapOpen, gapExtend,
            (seqs, scoring) => global::Seqeron.Genomics.Alignment.SequenceAligner.MultipleAlign(seqs, scoring, ParseMsaOptions(null, gapModel)));

    [McpServerTool(Name = "multiple_align_anchored", Title = "MSA — Anchored Center-Star (approximate, suffix-tree seeds)", ReadOnly = true)]
    [Description("Approximate, faster center-star MSA: picks the center by 4-mer count cosine similarity, builds a suffix tree on it and aligns every other sequence with exact-match anchors (MUMmer-style seed, weighted colinear chain, Needleman–Wunsch — or Gotoh affine with gapModel='affine' — on the gaps between anchors). The pairwise alignments are NOT guaranteed optimal, so the sum-of-pairs score can be lower than multiple_align's. Default linear gap model (gapOpen ignored); gapModel='affine' charges gapOpen + k·gapExtend per gap run and reports the natural affine sum-of-pairs score.")]
    public static MultipleAlignmentResultDto MultipleAlignAnchored(
        [Description("DNA sequences to align (at least one).")] string[] sequences,
        [Description("Score for matching base.")] int match = 1,
        [Description("Score for mismatching base.")] int mismatch = -1,
        [Description(GapOpenDescription)] int gapOpen = -2,
        [Description(GapExtendDescription)] int gapExtend = -1,
        [Description(GapModelDescription)] string gapModel = "linear") =>
        RunMsa(sequences, match, mismatch, gapOpen, gapExtend,
            (seqs, scoring) => global::Seqeron.Genomics.Alignment.SequenceAligner.MultipleAlignAnchored(seqs, scoring, ParseMsaOptions(null, gapModel)));

    [McpServerTool(Name = "multiple_align_progressive", Title = "MSA — Progressive Guide-Tree Alignment (Clustal-style)", ReadOnly = true)]
    [Description("Guide-tree progressive MSA (Feng & Doolittle 1987 / Clustal-style): the ClustalW pairwise distance of every pairwise alignment, d = 1 − (identical aligned residue pairs) / (length of the shorter sequence) (Thompson, Higgins & Gibson 1994; pairalign.c tmat = (100 − 100·count/MIN(len1, len2)) / 100), gives the distance matrix, a guide tree (UPGMA by default, guideTree='nj' = ClustalW's neighbour-joining tree) orders the merges, and profiles are aligned by profile–profile dynamic programming with sum-of-pairs column scores (once a gap, always a gap). Default linear gap model (gapOpen ignored); gapModel='affine' uses a Gotoh three-state profile DP (occupancy-scaled gapOpen/gapExtend, ClustalW/MUSCLE style) and reports the natural affine sum-of-pairs score. distance='clustalw' takes the guide-tree distances from ClustalW 2.1's own affine pairwise stage instead (percent identity of clustalw -ALIGN). Returns rows, majority consensus and sum-of-pairs score.")]
    public static MultipleAlignmentResultDto MultipleAlignProgressive(
        [Description("DNA sequences to align (at least one).")] string[] sequences,
        [Description("Score for matching base.")] int match = 1,
        [Description("Score for mismatching base.")] int mismatch = -1,
        [Description(GapOpenDescription)] int gapOpen = -2,
        [Description(GapExtendDescription)] int gapExtend = -1,
        [Description(DistanceModeDescription)] string distance = "linear",
        [Description(GapModelDescription)] string gapModel = "linear",
        [Description(GuideTreeDescription)] string guideTree = "upgma") =>
        RunMsa(sequences, match, mismatch, gapOpen, gapExtend,
            (seqs, scoring) => global::Seqeron.Genomics.Alignment.SequenceAligner.MultipleAlignProgressive(seqs, scoring, ParseMsaOptions(distance, gapModel, guideTree)));

    [McpServerTool(Name = "multiple_align_iterative", Title = "MSA — Progressive + Iterative Refinement (MUSCLE-style)", ReadOnly = true)]
    [Description("Progressive MSA followed by MUSCLE-style tree-dependent restricted-partitioning refinement (Edgar 2004, stage 3): each guide-tree edge splits the alignment into two profiles that are re-aligned and kept only if the sum-of-pairs score strictly improves; repeated until convergence or maxIterations passes. The result's sum-of-pairs score is never below multiple_align_progressive's. Default linear gap model (gapOpen ignored); gapModel='affine' refines under the natural affine sum-of-pairs objective with the affine profile DP. distance='clustalw' selects ClustalW 2.1's pairwise stage for the guide tree and guideTree='nj' neighbour-joining (see multiple_align_progressive); the refinement partitions are the edges of that tree.")]
    public static MultipleAlignmentResultDto MultipleAlignIterative(
        [Description("DNA sequences to align (at least one).")] string[] sequences,
        [Description("Score for matching base.")] int match = 1,
        [Description("Score for mismatching base.")] int mismatch = -1,
        [Description(GapOpenDescription)] int gapOpen = -2,
        [Description(GapExtendDescription)] int gapExtend = -1,
        [Description("Upper bound on full refinement passes (>= 1).")] int maxIterations = 16,
        [Description(DistanceModeDescription)] string distance = "linear",
        [Description(GapModelDescription)] string gapModel = "linear",
        [Description(GuideTreeDescription)] string guideTree = "upgma")
    {
        if (maxIterations < 1)
            throw new ArgumentOutOfRangeException(nameof(maxIterations), "maxIterations must be >= 1.");

        var options = ParseMsaOptions(distance, gapModel, guideTree);
        return RunMsa(sequences, match, mismatch, gapOpen, gapExtend,
            (seqs, scoring) => global::Seqeron.Genomics.Alignment.SequenceAligner.MultipleAlignIterative(seqs, scoring, options, maxIterations));
    }

    [McpServerTool(Name = "multiple_align_consistency", Title = "MSA — Consistency-Based Alignment (T-Coffee)", ReadOnly = true)]
    [Description("Consistency-based MSA (T-Coffee; Notredame, Higgins & Heringa 2000): a primary library of global (Needleman–Wunsch) and local (Smith–Waterman) pairwise alignments weighted by percent identity, extended by the triplet consistency transformation (direct weight + Σ min(W1, W2) over third sequences), then progressively aligned along a guide tree (UPGMA by default; guideTree='nj' = neighbour-joining as in the T-Coffee paper) with the extended-library weights as position-specific scores and zero gap penalties (T-Coffee default). The scoring parameters seed the pairwise library (gapModel='affine': Gotoh global and Smith–Waterman–Gotoh local alignments with gapOpen + k·gapExtend) and define totalScore, the sum of pairs of the result (linear column score, or natural affine costs with gapModel='affine'). distance='clustalw' takes the guide-tree distances from ClustalW 2.1's own pairwise stage.")]
    public static MultipleAlignmentResultDto MultipleAlignConsistency(
        [Description("DNA sequences to align (at least one).")] string[] sequences,
        [Description("Score for matching base (pairwise library and reported score).")] int match = 1,
        [Description("Score for mismatching base.")] int mismatch = -1,
        [Description(GapOpenDescription)] int gapOpen = -2,
        [Description(GapExtendDescription)] int gapExtend = -1,
        [Description(DistanceModeDescription)] string distance = "linear",
        [Description(GapModelDescription)] string gapModel = "linear",
        [Description(GuideTreeDescription)] string guideTree = "upgma") =>
        RunMsa(sequences, match, mismatch, gapOpen, gapExtend,
            (seqs, scoring) => global::Seqeron.Genomics.Alignment.SequenceAligner.MultipleAlignConsistency(seqs, scoring, ParseMsaOptions(distance, gapModel, guideTree)));

    private static global::Seqeron.Genomics.Infrastructure.MsaOptions ParseMsaOptions(string? distance, string? gapModel, string? guideTree = null) => new()
    {
        GuideTree = (guideTree ?? "upgma").Trim().ToLowerInvariant() switch
        {
            "upgma" => global::Seqeron.Genomics.Infrastructure.MsaGuideTree.Upgma,
            "nj" or "neighbor_joining" or "neighbour_joining" or "neighborjoining" => global::Seqeron.Genomics.Infrastructure.MsaGuideTree.NeighborJoining,
            _ => throw new ArgumentException("guideTree must be 'upgma' or 'nj'.", nameof(guideTree)),
        },
        Distance = (distance ?? "linear").Trim().ToLowerInvariant() switch
        {
            "linear" => global::Seqeron.Genomics.Infrastructure.MsaDistanceMode.Linear,
            "clustalw" or "clustalw_affine" or "clustalwaffine" => global::Seqeron.Genomics.Infrastructure.MsaDistanceMode.ClustalWAffine,
            _ => throw new ArgumentException("distance must be 'linear' or 'clustalw'.", nameof(distance)),
        },
        GapModel = (gapModel ?? "linear").Trim().ToLowerInvariant() switch
        {
            "linear" => global::Seqeron.Genomics.Infrastructure.MsaGapModel.Linear,
            "affine" => global::Seqeron.Genomics.Infrastructure.MsaGapModel.Affine,
            _ => throw new ArgumentException("gapModel must be 'linear' or 'affine'.", nameof(gapModel)),
        },
    };

    private static MultipleAlignmentResultDto RunMsa(
        string[] sequences, int match, int mismatch, int gapOpen, int gapExtend,
        Func<IEnumerable<global::Seqeron.Genomics.Core.DnaSequence>, global::Seqeron.Genomics.Infrastructure.ScoringMatrix?, global::Seqeron.Genomics.Infrastructure.MultipleAlignmentResult> align)
    {
        if (sequences == null || sequences.Length == 0)
        {
            var e = global::Seqeron.Genomics.Infrastructure.MultipleAlignmentResult.Empty;
            return new MultipleAlignmentResultDto(e.AlignedSequences, e.Consensus, e.TotalScore);
        }

        var dnaList = new List<global::Seqeron.Genomics.Core.DnaSequence>(sequences.Length);
        for (int i = 0; i < sequences.Length; i++)
        {
            if (!global::Seqeron.Genomics.Core.DnaSequence.TryCreate(sequences[i], out var d))
                throw new ArgumentException($"Invalid DNA at sequences[{i}].", nameof(sequences));
            dnaList.Add(d!);
        }

        var scoring = new global::Seqeron.Genomics.Infrastructure.ScoringMatrix(match, mismatch, gapOpen, gapExtend);
        var r = align(dnaList, scoring);
        return new MultipleAlignmentResultDto(r.AlignedSequences, r.Consensus, r.TotalScore);
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
}
