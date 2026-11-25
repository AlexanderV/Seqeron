using System.ComponentModel;
using System.Globalization;
using ModelContextProtocol.Server;

namespace Seqeron.Mcp.Alignment.Tools;

// B09 A2-28: dual-affine (ksw2 / minimap2 style) pairwise alignment.
public partial class AlignmentTools
{

    [McpServerTool(Name = "dual_affine_align", Title = "Pairwise — Dual Affine Gaps (ksw2 / minimap2 style)", ReadOnly = true)]
    [Description("Pairwise alignment with the TWO-PIECE (dual) affine gap model of ksw2 / minimap2 (ksw2_extd2): a gap of length k scores max(gapOpen1 + k*gapExtend1, gapOpen2 + k*gapExtend2) (a short-gap and a long-gap regime; all four gap numbers are scores <= 0, ksw2's q + k*e cost is gapOpen = -q, gapExtend = -e; defaults = minimap2 -A 2 -B 4 -O 4,24 -E 2,1), Gotoh's O(m n) program with two gap states per direction, mode 'global' (end to end) or 'local' (Smith-Waterman), optional diagonal band (ksw2 w: only the cells |i - j| <= band; a global alignment needs |m - n| <= band). VERIFICATION: 1 600 committed cases (400 global, 400 local, 400 banded global, 400 banded local; scores equal Biopython 1.88 PairwiseAligner with the callable gap score and an independent Waterman-Smith-Beyer program, every alignment re-scores to its score). Z-drop / end-bonus extension alignment is the separate tool ksw2_extend.")]
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
        [Description(PairFormatParamDescription)] string? format = null,
        [Description(NamesParamDescription)] string[]? names = null,
        CancellationToken cancellationToken = default,
        [Description(TrimEndGapsParamDescription)] bool trimEndGaps = false)
    {
        return FormatPair(DualAffineAlignCore(sequence1, sequence2, mode, match, mismatch, gapOpen1, gapExtend1, gapOpen2, gapExtend2, band, cancellationToken), format, sequence1.Length, sequence2.Length, names, sequence1: sequence1, sequence2: sequence2, trimEndGaps: trimEndGaps);
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
