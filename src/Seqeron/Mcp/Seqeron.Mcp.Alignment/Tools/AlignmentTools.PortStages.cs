using System.ComponentModel;
using ModelContextProtocol.Server;

namespace Seqeron.Mcp.Alignment.Tools;

// B09 A2-36: tools for the public stage ports that had no tool (the multiple_align_* tools run only whole pipelines): ClustalW pairwise stages,
// T-Coffee pair method / lalign / primary library / guide tree / k-tuple similarity, and the sum-of-pairs score of an alignment.
// Every tool delegates to SequenceAligner; MUSCLE's guide tree of an alignment is muscle_msa_command (mode 'maketree').
public partial class AlignmentTools
{

    [McpServerTool(Name = "sum_of_pairs_score", Title = "MSA — Sum-of-Pairs Score of an Alignment", ReadOnly = true)]
    [Description("Sum-of-pairs score of a multiple alignment (equal-length rows, '-' = gap). gapModel 'linear' (default): column sum of pairs - residue pairs match / mismatch, residue against gap gapExtend, gap against gap 0 (the totalScore of the linear aligners); 'affine': the natural affine objective (Altschul 1989) - for every pair of rows the gap-gap columns are removed and every maximal gap run costs gapOpen + length x gapExtend. Scores: match / mismatch / gapOpen / gapExtend (integers; negative = penalty); with matrix (BLOSUM62, PAM250, EDNAFULL or any EMBOSS matrix) the residue pairs are scored by that table (protein / IUPAC rows, residues case-insensitive) and the gap scores may be fractional - the score is then a real number. This is the objective the exact (Carrillo-Lipman) and affine aligners optimise; MUSCLE's own SP score is muscle_msa_command mode 'spscore'.")]
    public static SumOfPairsScoreResult SumOfPairsScore(
        [Description("Equal-length aligned rows, '-' = gap.")] string[] alignedRows,
        [Description("Score for a matching residue pair (ignored with matrix).")] int match = 1,
        [Description("Score for a mismatching residue pair (ignored with matrix).")] int mismatch = -1,
        [Description("Gap-open score added once per gap run (affine model only; negative = penalty; integer unless matrix is given).")] double gapOpen = -2,
        [Description("Gap score per gap position (negative = penalty; integer unless matrix is given).")] double gapExtend = -1,
        [Description("Gap model: 'linear' (default) or 'affine'.")] string gapModel = "linear",
        [Description(MatrixParamDescription)] string? matrix = null)
    {
        if (alignedRows is null || alignedRows.Length == 0)
            throw new ArgumentException("alignedRows must contain at least one row.", nameof(alignedRows));
        var model = (gapModel ?? "linear").Trim().ToLowerInvariant() switch
        {
            "linear" => global::Seqeron.Genomics.Infrastructure.MsaGapModel.Linear,
            "affine" => global::Seqeron.Genomics.Infrastructure.MsaGapModel.Affine,
            _ => throw new ArgumentException("gapModel must be 'linear' or 'affine'.", nameof(gapModel)),
        };
        var table = ResolveMatrix(matrix);
        if (table is not null)
        {
            double score = global::Seqeron.Genomics.Alignment.SequenceAligner.MatrixSumOfPairsScore(
                alignedRows, new global::Seqeron.Genomics.Infrastructure.MatrixScoring(table, gapOpen, gapExtend), model);
            return new SumOfPairsScoreResult(score, model.ToString(), $"matrix {matrix}");
        }

        int total = global::Seqeron.Genomics.Alignment.SequenceAligner.SumOfPairsScore(
            alignedRows,
            new global::Seqeron.Genomics.Infrastructure.ScoringMatrix(match, mismatch, RequireInteger(gapOpen, nameof(gapOpen)), RequireInteger(gapExtend, nameof(gapExtend))),
            model);
        return new SumOfPairsScoreResult(total, model.ToString(), $"match {match} mismatch {mismatch}");
    }
}
