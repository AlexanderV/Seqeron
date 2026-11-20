using System.ComponentModel;
using ModelContextProtocol.Server;

namespace Seqeron.Mcp.Alignment.Tools;

// B09 A3-10: per-direction / per-end gap scores (Biopython PairwiseAligner's twelve gap-score properties); delegates to SequenceAligner.{Global,SemiGlobal,Local}Align[Affine](..., GapScoreProfile).
public partial class AlignmentTools
{
    [McpServerTool(Name = "align_gap_profile", Title = "Pairwise — Alignment with per-direction / per-end gap scores", ReadOnly = true)]
    [Description("Pairwise alignment whose gap scores differ by direction and by position: Biopython 1.88 PairwiseAligner's twelve properties {open,extend}_{internal,left,right}_{insertion,deletion}_score (insertion = a gap in sequence1 = the target, deletion = a gap in sequence2 = the query; the deprecated names target_left_open_gap_score, query_right_extend_gap_score, ... and every aggregate such as open_gap_score, end_extend_gap_score, left_gap_score, insertion_score are accepted too, more specific names override general ones). gapScores = name -> value, Biopython convention (an open score is the score of the FIRST gap position, extend the score of each further one; all <= 0 up to Open + Extend <= 0). A gap run that begins the alignment is a left (start) run, one that finishes it a right (end) run, every other run is internal; unnamed slots use the interior gapOpen / gapExtend of this tool (Seqeron convention: a gap of length k scores gapOpen + k*gapExtend, linear when affine = false: k*gapExtend). A left / right score of 0 is a free overhang (the semi-global family); endGaps (comma-separated flags FreeStart1, FreeEnd1, FreeStart2, FreeEnd2, or the names Fitting, Overlap..., All) makes the matching slots free unless gapScores names them. mode = global (default), semi_global (the result is tagged SemiGlobal when an end is free) or local (a local alignment never starts or ends in a gap: only the internal scores apply, as in Biopython's local mode). affine = false is the linear program (the per-position score is the extend value; open is ignored). With matrix = BLOSUM62 | PAM250 | EDNAFULL ... residues are scored by that table and every gap score may be fractional (exactScore = the real score). Coordinates: a free start / end run is an overhang outside startPosition / endPosition, a priced one is inside. No new dynamic program: the existing linear / Gotoh / Smith-Waterman cores take the gap scores as a function of the matrix position. VERIFICATION: Biopython 1.88 PairwiseAligner score on 1000 random pairs (scalar, EDNAFULL, BLOSUM62; global / local; linear / affine; free, equal and independently priced ends) equal on 998; the two others are a Biopython linear Needleman-Wunsch corner (its first column of the last row is charged the right score nA*right_deletion_extend, its own Gotoh routine the left score). Inputs above 25,000,000 dynamic-programming cells are refused (the linear-space program has one gap model).")]
    public static AlignmentResultDto AlignGapProfile(
        [Description("First sequence (Biopython target); uppercased internally.")] string sequence1,
        [Description("Second sequence (Biopython query).")] string sequence2,
        [Description("global (default) | semi_global | local (case-insensitive).")] string mode = "global",
        [Description("true (default) = affine gaps (Gotoh), false = linear (per-position extend score).")] bool affine = true,
        [Description("Score for matching base (ignored when matrix is given).")] int match = 1,
        [Description("Score for mismatching base (ignored when matrix is given).")] int mismatch = -1,
        [Description("Interior gap-open score for slots not named in gapScores (Seqeron convention: added once per run; <= 0). May be fractional when matrix is given.")] double gapOpen = -2,
        [Description("Interior gap score per position for slots not named in gapScores (<= 0). May be fractional when matrix is given.")] double gapExtend = -1,
        [Description("Biopython gap-score properties, name -> value (e.g. {\"open_left_insertion_score\": 0, \"extend_left_insertion_score\": 0, \"end_deletion_score\": -3, \"open_gap_score\": -5, \"extend_gap_score\": -1}); see the description for the accepted names.")] Dictionary<string, double>? gapScores = null,
        [Description("Optional free end gaps (comma-separated FreeStart1, FreeEnd1, FreeStart2, FreeEnd2 or None, Fitting, FittingSequence2, OverlapSuffix1Prefix2, OverlapSuffix2Prefix1, All); slots named in gapScores keep their score.")] string? endGaps = null,
        [Description(MatrixParamDescription)] string? matrix = null,
        CancellationToken cancellationToken = default)
    {
        using var cancelScope = global::Seqeron.Genomics.Alignment.SequenceAligner.BeginCancellationScope(cancellationToken);
        RequireSequence(sequence1, nameof(sequence1));
        RequireSequence(sequence2, nameof(sequence2));
        string which = (mode ?? "global").Trim().ToLowerInvariant().Replace("-", "_");
        if (which is not ("global" or "semi_global" or "local"))
            throw new ArgumentException("mode must be global, semi_global or local.", nameof(mode));

        var profile = BuildGapProfile(gapScores, endGaps);
        var sub = ResolveMatrix(matrix);
        global::Seqeron.Genomics.Infrastructure.AlignmentResult r;
        if (sub is not null)
        {
            var ms = new global::Seqeron.Genomics.Infrastructure.MatrixScoring(sub, gapOpen, gapExtend);
            r = (which, affine) switch
            {
                ("global", true) => global::Seqeron.Genomics.Alignment.SequenceAligner.GlobalAlignAffine(sequence1, sequence2, ms, profile),
                ("global", false) => global::Seqeron.Genomics.Alignment.SequenceAligner.GlobalAlign(sequence1, sequence2, ms, profile),
                ("semi_global", true) => global::Seqeron.Genomics.Alignment.SequenceAligner.SemiGlobalAlignAffine(sequence1, sequence2, ms, profile),
                ("semi_global", false) => global::Seqeron.Genomics.Alignment.SequenceAligner.SemiGlobalAlign(sequence1, sequence2, ms, profile),
                (_, true) => global::Seqeron.Genomics.Alignment.SequenceAligner.LocalAlignAffine(sequence1, sequence2, ms, profile),
                _ => global::Seqeron.Genomics.Alignment.SequenceAligner.LocalAlign(sequence1, sequence2, ms, profile),
            };
        }
        else
        {
            var sc = new global::Seqeron.Genomics.Infrastructure.ScoringMatrix(match, mismatch, RequireInteger(gapOpen, nameof(gapOpen)), RequireInteger(gapExtend, nameof(gapExtend)));
            r = (which, affine) switch
            {
                ("global", true) => global::Seqeron.Genomics.Alignment.SequenceAligner.GlobalAlignAffine(sequence1, sequence2, sc, profile),
                ("global", false) => global::Seqeron.Genomics.Alignment.SequenceAligner.GlobalAlign(sequence1, sequence2, sc, profile),
                ("semi_global", true) => global::Seqeron.Genomics.Alignment.SequenceAligner.SemiGlobalAlignAffine(sequence1, sequence2, sc, profile),
                ("semi_global", false) => global::Seqeron.Genomics.Alignment.SequenceAligner.SemiGlobalAlign(sequence1, sequence2, sc, profile),
                (_, true) => global::Seqeron.Genomics.Alignment.SequenceAligner.LocalAlignAffine(sequence1, sequence2, sc, profile),
                _ => global::Seqeron.Genomics.Alignment.SequenceAligner.LocalAlign(sequence1, sequence2, sc, profile),
            };
        }

        return ToAlignmentDto(r);
    }

    private static global::Seqeron.Genomics.Infrastructure.GapScoreProfile BuildGapProfile(Dictionary<string, double>? gapScores, string? endGaps)
    {
        var named = gapScores is null || gapScores.Count == 0
            ? new global::Seqeron.Genomics.Infrastructure.GapScoreProfile()
            : global::Seqeron.Genomics.Infrastructure.GapScoreProfile.FromBiopython(gapScores);
        if (string.IsNullOrWhiteSpace(endGaps))
            return named;
        var free = global::Seqeron.Genomics.Infrastructure.GapScoreProfile.FromEndGaps(ResolveEndGapsOrNone(endGaps));
        return named with
        {
            Start1 = named.Start1 ?? free.Start1,
            End1 = named.End1 ?? free.End1,
            Start2 = named.Start2 ?? free.Start2,
            End2 = named.End2 ?? free.End2,
        };
    }

    private static global::Seqeron.Genomics.Infrastructure.EndGaps ResolveEndGapsOrNone(string endGaps) =>
        endGaps.Trim().Equals("None", StringComparison.OrdinalIgnoreCase)
            ? global::Seqeron.Genomics.Infrastructure.EndGaps.None
            : ResolveEndGaps(endGaps);
}
