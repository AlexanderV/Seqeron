// SPDX-License-Identifier: GPL-2.0-or-later
// This file is part of Seqeron.Mcp.EmbossPairwiseAlignment, licensed under the terms stated in LICENSE (this directory).
using System.ComponentModel;
using ModelContextProtocol.Server;
using Seqeron.Mcp.Alignment.Tools;
using static Seqeron.Mcp.Alignment.Tools.AlignmentTools;

namespace Seqeron.Mcp.EmbossPairwiseAlignment.Tools;

public partial class EmbossPairwiseAlignmentTools
{
    [McpServerTool(Name = "matcher_align", Title = "Pairwise — EMBOSS matcher: Top-N Non-Intersecting Local Alignments (Huang-Miller SIM)", ReadOnly = true)]
    [Description("Waterman–Eggert (1987, J Mol Biol 197:723) declumping: the best local alignment, then the best local alignment that shares no ALIGNED RESIDUE PAIR with it, then the best sharing no pair with either, ... up to maxAlignments alignments of score >= minScore (and always > 0) - the 'alternatives' of EMBOSS matcher -alternatives N (Huang-Miller SIM; exact same definition, verified round by round on 2 902 reported alignments) and FASTA lalign. An aligned pair is a diagonal column (a_i opposite b_j); gap columns are not pairs, so two reported alignments may cross through a gap column but never align the same two residues. Alignment 1 is exactly local_align (affine = false) / local_align_affine (affine = true); scores are non-increasing; ties between co-optimal alignments are broken, with tieOrder = WatermanEggert (default), as local_align does (first maximal cell in row-major order, traceback diagonal > up > left, M > X > Y), so on tie-rich inputs a later alignment can differ from matcher's (its score is always the true maximum of the remaining problem); tieOrder = Sim reproduces EMBOSS matcher's own tie order exactly (the Huang-Miller SIM program of matcher.c ported line by line: K-vertex LIST, forward pass by start point, locate / small-pass recomputation, linear-space Diff; 1 138 oracle rows equal alignment by alignment, 59 more correctly abort) - it needs integer scores (or a matrix with integer entries) and, like the program, depends on maxAlignments (the program's LIST size); the program aborts when fewer than maxAlignments alignments exist, this tool then returns the alignments found. affine = false: linear gaps (every gap position costs gapExtend, gapOpen ignored); affine = true: a gap of length k costs gapOpen + k*gapExtend (EMBOSS -gapopen G -gapextend E = gapOpen -(G - E), gapExtend -E). With matrix = BLOSUM62 | PAM250 | EDNAFULL residues are scored by that table, gaps may be fractional and minScore/exactScore are real scores. Returns alignments (same fields as local_align) in the order found; fewer than maxAlignments when nothing positive >= minScore is left (empty list when no pair scores > 0). Memory: up to 25,000,000 cells ((len1+1)*(len2+1)) the original recomputation scheme (matrices kept, only cells whose value can change are recomputed after each report; O(len1*len2) memory); ABOVE that (e.g. two 5-20 kb sequences) the same alignments - identical rows, scores, coordinates, order and ties - are computed in exact LINEAR SPACE (O(len1 + len2) memory per alignment: one run of the Hirschberg / Myers-Miller / Huang-Miller locate-then-solve local program per reported alignment with the aligned pairs of the earlier ones forbidden, the SIM / EMBOSS matcher scheme), so the running time grows with maxAlignments (about one extra local alignment per reported alignment) while the memory does not (20 kb x 20 kb: tens of MB instead of 1.6 GB); the linear-space program needs non-positive gap scores (a positive gap score on such an input is rejected with a clear error). With tieOrder = Sim (EMBOSS matcher's own tie order) the sequences are also read as the matcher program reads them (A4-3): only letters and '*' are kept; gap characters (- . ~ ?), '#', '+', '_', digits, blanks and every other character are dropped first and positions count the cleaned sequence. With tieOrder = Sim, maxAlignments = 2147483647 means every alignment (as in the Waterman-Eggert order); the matcher candidate list is bounded at 4 194 304 start points (SequenceAligner.TopNSimMaxVertices, ArgumentException).")]
    public static LocalAlignmentListResult MatcherAlign(
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
        [Description("Which co-optimal alignment is reported on a tie: WatermanEggert (default; the fixed local_align order) or Sim (EMBOSS matcher's own order, exact; Sim also reads the sequences as the matcher program does: only letters and '*' are kept, gap characters, digits, blanks and every other character are dropped first and positions count the cleaned sequence).")] string tieOrder = "Sim",
        CancellationToken cancellationToken = default)
    {
        using var cancelScope = global::Seqeron.Genomics.Alignment.SequenceAligner.BeginCancellationScope(cancellationToken);
        RequireSequence(sequence1, nameof(sequence1));
        RequireSequence(sequence2, nameof(sequence2));
        if (!Enum.TryParse<global::Seqeron.Genomics.Emboss.PairwiseAlignment.TopNTieOrder>(tieOrder, ignoreCase: true, out var order) || !Enum.IsDefined(order))
            throw new ArgumentException($"tieOrder must be WatermanEggert or Sim, got '{tieOrder}'.", nameof(tieOrder));
        var sub = ResolveMatrix(matrix);
        if (sub is not null)
        {
            var ms = new global::Seqeron.Genomics.Infrastructure.MatrixScoring(sub, gapOpen, gapExtend);
            return ToAlignmentListDto(affine
                ? global::Seqeron.Genomics.Emboss.PairwiseAlignment.EmbossPairwisePrograms.LocalAlignAffineTopN(sequence1, sequence2, ms, maxAlignments, minScore, order)
                : global::Seqeron.Genomics.Emboss.PairwiseAlignment.EmbossPairwisePrograms.LocalAlignTopN(sequence1, sequence2, ms, maxAlignments, minScore, order));
        }

        var scoring = new global::Seqeron.Genomics.Infrastructure.ScoringMatrix(match, mismatch, RequireInteger(gapOpen, nameof(gapOpen)), RequireInteger(gapExtend, nameof(gapExtend)));
        int min = (int)Math.Clamp(Math.Ceiling(minScore), int.MinValue, int.MaxValue);
        return ToAlignmentListDto(affine
            ? global::Seqeron.Genomics.Emboss.PairwiseAlignment.EmbossPairwisePrograms.LocalAlignAffineTopN(sequence1, sequence2, scoring, maxAlignments, min, order)
            : global::Seqeron.Genomics.Emboss.PairwiseAlignment.EmbossPairwisePrograms.LocalAlignTopN(sequence1, sequence2, scoring, maxAlignments, min, order));
    }
}
