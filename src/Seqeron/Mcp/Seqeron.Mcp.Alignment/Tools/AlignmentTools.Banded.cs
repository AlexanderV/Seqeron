using System.ComponentModel;
using ModelContextProtocol.Server;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Mcp.Alignment.Tools;

// B09 R14-12 (F305): banded (Ukkonen / Chao-Pearson-Miller) and X-drop alignment on the core aligners.
public partial class AlignmentTools
{
    [McpServerTool(Name = "banded_align", Title = "Pairwise — Banded / X-drop Alignment (Ukkonen, Chao-Pearson-Miller)", ReadOnly = true)]
    [Description("Banded pairwise alignment: only the cells of the diagonals diagonalLow .. diagonalHigh (j - i) of the dynamic program are computed, O(length1 x band) time and memory instead of O(length1 x length2) (30 000 x 30 000 residues with affine gaps: about 1 s instead of 98 s in the Debug build). mode global (default), local or semi_global (endGaps as in semi_global_align), linear (default) or affine gaps (gapOpen + k x gapExtend), match / mismatch or a matrix as in global_align. The alignment is the one of the full program (the same rows, score and coordinates) whenever the full alignment lies inside the band; otherwise the best alignment inside the band (never better than the full optimum). GLOBAL: bandWidth w gives the diagonals from min(0, n - m) - w to max(0, n - m) + w; exact = true (the default when no bandWidth is given) doubles the band from initialBandWidth until every path that leaves it is proven to score strictly less (Ukkonen's bound; the result is then the full program's alignment) and reports proven. LOCAL / SEMI-GLOBAL: give diagonalLow and diagonalHigh (a semi-global band must contain the diagonal 0). xDrop prunes cells more than X below the best score seen (BLAST gapped extension; a heuristic, equal to the banded result when X is at least the largest fall of the score along the optimal path). VERIFICATION: equal to GlobalAlign / GlobalAlignAffine / LocalAlign[Affine] / SemiGlobalAlign[Affine] field for field with a band that covers the matrix and with the full alignment inside a narrow band (seeded random pairs, linear / affine, scalar / BLOSUM62).")]
    public static BandedAlignResult BandedAlignPair(
        [Description("First sequence.")] string sequence1,
        [Description("Second sequence.")] string sequence2,
        [Description("global (default), local or semi_global.")] string mode = "global",
        [Description("true: affine gaps (gapOpen + k x gapExtend); false (default): linear gaps (every position costs gapExtend).")] bool affine = false,
        [Description("Score for matching base.")] int match = 1,
        [Description("Score for mismatching base.")] int mismatch = -1,
        [Description("Gap-open score (affine only; <= 0). May be fractional when matrix is given.")] double gapOpen = -2,
        [Description("Gap score per gap position (<= 0). May be fractional when matrix is given.")] double gapExtend = -1,
        [Description(MatrixParamDescription)] string? matrix = null,
        [Description("global: the half-width w of the band beyond the corridor between the diagonal 0 and the end diagonal n - m (omit for the exact mode).")] int? bandWidth = null,
        [Description("global: double the band until the optimum is proven (default true when bandWidth is omitted).")] bool? exact = null,
        [Description("exact mode: the first band half-width (default 16).")] int initialBandWidth = 16,
        [Description("local / semi_global: the lowest diagonal j - i of the band.")] int? diagonalLow = null,
        [Description("local / semi_global: the highest diagonal j - i of the band.")] int? diagonalHigh = null,
        [Description(EndGapsParamDescription)] string? endGaps = null,
        [Description("X-drop pruning: cells more than xDrop below the best score so far are dropped (heuristic; omit for none).")] int? xDrop = null,
        CancellationToken cancellationToken = default)
    {
        using var cancelScope = global::Seqeron.Genomics.Alignment.SequenceAligner.BeginCancellationScope(cancellationToken);
        RequireSequence(sequence1, nameof(sequence1));
        RequireSequence(sequence2, nameof(sequence2));
        var sub = ResolveMatrix(matrix);
        ScoringMatrix? scalar = null;
        MatrixScoring? matrixScoring = null;
        if (sub is not null)
            matrixScoring = new MatrixScoring(sub, affine ? gapOpen : 0, gapExtend);
        else
            scalar = new ScoringMatrix(match, mismatch, affine ? RequireInteger(gapOpen, nameof(gapOpen)) : 0, RequireInteger(gapExtend, nameof(gapExtend)));
        string kind = (mode ?? "global").Trim().ToLowerInvariant().Replace("-", "_", StringComparison.Ordinal);
        if (kind is not ("global" or "local" or "semi_global" or "semiglobal"))
            throw new ArgumentException($"Unknown mode '{mode}': global, local or semi_global.", nameof(mode));
        if (kind == "global")
        {
            bool runExact = exact ?? bandWidth is null;
            if (runExact && bandWidth is null)
            {
                var r = matrixScoring is not null
                    ? SequenceAligner.GlobalAlignBandedExact(sequence1, sequence2, matrixScoring, affine, initialBandWidth)
                    : SequenceAligner.GlobalAlignBandedExact(sequence1, sequence2, scalar, affine, initialBandWidth);
                int d = sequence2.Length - sequence1.Length;
                return new BandedAlignResult(ToAlignmentDto(r.Alignment), Math.Min(0, d) - r.BandWidth, Math.Max(0, d) + r.BandWidth, r.Proven, r.Rounds, r.Cells);
            }

            int w = bandWidth ?? throw new ArgumentException("bandWidth is required when exact is false.", nameof(bandWidth));
            if (w < 0)
                throw new ArgumentOutOfRangeException(nameof(bandWidth), "bandWidth must not be negative.");
            int dd = sequence2.Length - sequence1.Length;
            var options = new BandedAlignOptions { Type = AlignmentType.Global, Affine = affine, Scoring = scalar, MatrixScoring = matrixScoring, DiagonalLow = Math.Min(0, dd) - w, DiagonalHigh = Math.Max(0, dd) + w, XDrop = xDrop };
            return new BandedAlignResult(ToAlignmentDto(SequenceAligner.BandedAlign(sequence1, sequence2, options)), options.DiagonalLow, options.DiagonalHigh, null, 1, 0);
        }

        if (diagonalLow is null || diagonalHigh is null)
            throw new ArgumentException("diagonalLow and diagonalHigh are required for local and semi_global alignment.", nameof(diagonalLow));
        bool local = kind == "local";
        var o = new BandedAlignOptions
        {
            Type = local ? AlignmentType.Local : AlignmentType.SemiGlobal,
            Affine = affine,
            Scoring = scalar,
            MatrixScoring = matrixScoring,
            EndGaps = local ? EndGaps.None : ResolveEndGaps(endGaps),
            DiagonalLow = diagonalLow.Value,
            DiagonalHigh = diagonalHigh.Value,
            XDrop = xDrop,
        };
        return new BandedAlignResult(ToAlignmentDto(SequenceAligner.BandedAlign(sequence1, sequence2, o)), o.DiagonalLow, o.DiagonalHigh, null, 1, 0);
    }
}
