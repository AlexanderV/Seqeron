using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Alignment;

public static partial class SequenceAligner
{
    #region Multiple alignment with affine gaps (natural sum-of-pairs objective, affine profile-profile DP)

    // Affine-gap multiple sequence alignment (audit A1-10).
    //
    // OBJECTIVE (exact). The sum-of-pairs score with NATURAL affine gap costs (Altschul 1989, J Theor Biol 138:297, "gap scores
    // should be defined as the sum of gap scores for the implied pairwise alignments"): every pair of rows is projected out of the
    // alignment, columns that are gaps in both rows are deleted, substitution columns score Match/Mismatch, and every maximal run of
    // columns of the projection in which the same row has a gap pays GapOpen + length * GapExtend (a gap in row i directly followed by a
    // gap in row j are two runs). SumOfPairsScore computes exactly this for any set of rows, O(k^2 L); for k = 2 it is the Gotoh score
    // that GlobalAlignAffine optimises, so a pairwise optimum is an MSA optimum.
    //
    // OPTIMISATION (heuristic for k >= 3, exact for k = 2).
    //  * Star methods: every sequence is aligned to the center with the exact affine pairwise aligner (see MultipleAlign).
    //  * Profile methods (progressive / iterative): a Gotoh three-state dynamic program over profile columns. The natural cost of a
    //    cross pair depends on the whole column history of both rows (a gap run of a pair continues across columns in which the other
    //    row is also gapped), which a three-state DP cannot track; the DP therefore uses the standard position-specific profile
    //    approximation (ClustalW / MUSCLE: gap penalties scaled by the number of residues of the column):
    //      M(i, j)   column i of A against column j of B  = sum over cross pairs of Match/Mismatch for residue pairs, GapExtend for
    //                a residue against an existing gap (an internal gap of the other profile: extension only), 0 for gap against gap;
    //      X(i, ·)   column i of A against a NEW all-gap column of B (cost per cross pair with a residue in column i: GapExtend, plus
    //                GapOpen when the run of new gap columns starts: r_i * |B| * (GapOpen + GapExtend) at the first column of the run,
    //                r_i * |B| * GapExtend for every further column), Y symmetric.
    //    For two single sequences (r = 1, |B| = 1) this is exactly the Gotoh recurrence of GlobalAlignAffine (same border
    //    initialisation of Flouri et al. 2015, same M > X > Y / own-state tie order), so the two-sequence result equals GlobalAlignAffine
    //    row for row; for profiles it is a heuristic whose result is scored with the exact natural objective. The iterative method keeps a
    //    re-alignment only if that exact objective strictly improves.


    public static MsaGapModel ResolveGapModel(MsaOptions? options)
    {
        var model = options?.GapModel ?? MsaGapModel.Linear;
        if (!Enum.IsDefined(model))
            throw new ArgumentOutOfRangeException(nameof(options), model, "Unknown MSA gap model.");
        return model;
    }

    /// <summary>
    /// Sum-of-pairs score of a multiple alignment. <see cref="MsaGapModel.Linear"/>: column sum of pairs (residue pairs Match/Mismatch,
    /// residue against gap <see cref="ScoringMatrix.GapExtend"/>, gap against gap 0) — the <c>TotalScore</c> of the linear aligners.
    /// <see cref="MsaGapModel.Affine"/>: the natural affine objective (Altschul 1989): for every pair of rows the gap–gap columns are
    /// removed and every maximal run of gap columns of the projection costs <see cref="ScoringMatrix.GapOpen"/> + length·<see cref="ScoringMatrix.GapExtend"/>.
    /// </summary>
    /// <param name="alignedRows">Equal-length rows; '-' is the gap.</param>
    /// <param name="scoring">Scores (default <see cref="SimpleDna"/>).</param>
    /// <param name="gapModel">Gap model of the objective.</param>
    /// <exception cref="ArgumentNullException">Rows are null.</exception>
    /// <exception cref="ArgumentException">Rows have different lengths.</exception>
    public static int SumOfPairsScore(IReadOnlyList<string> alignedRows, ScoringMatrix? scoring = null, MsaGapModel gapModel = MsaGapModel.Linear)
    {
        ArgumentNullException.ThrowIfNull(alignedRows);
        if (!Enum.IsDefined(gapModel))
            throw new ArgumentOutOfRangeException(nameof(gapModel), gapModel, "Unknown MSA gap model.");
        foreach (var row in alignedRows)
        {
            ArgumentNullException.ThrowIfNull(row);
            if (row.Length != alignedRows[0].Length)
                throw new ArgumentException("All aligned rows must have the same length.", nameof(alignedRows));
        }

        return ComputeSumOfPairsScore(alignedRows, scoring ?? SimpleDna, gapModel);
    }

    /// <summary>
    /// <see cref="SumOfPairsScore(IReadOnlyList{string}, ScoringMatrix?, MsaGapModel)"/> under a substitution matrix and real-valued (possibly
    /// fractional) gap scores: the exact sum-of-pairs of protein / IUPAC rows, returned as a real number (residues compared case-insensitively
    /// against the matrix alphabet, '-' = gap).
    /// </summary>
    /// <param name="alignedRows">Equal-length rows; '-' is the gap.</param>
    /// <param name="scoring">Substitution matrix and gap scores (<see cref="MatrixScoring.GapOpen"/> unused by the linear model).</param>
    /// <param name="gapModel">Gap model of the objective.</param>
    /// <exception cref="ArgumentNullException">Rows or scoring are null.</exception>
    /// <exception cref="ArgumentException">Rows have different lengths, or a residue is not in the matrix alphabet.</exception>
    public static double MatrixSumOfPairsScore(IReadOnlyList<string> alignedRows, MatrixScoring scoring, MsaGapModel gapModel = MsaGapModel.Linear)
    {
        ArgumentNullException.ThrowIfNull(alignedRows);
        ArgumentNullException.ThrowIfNull(scoring);
        if (!Enum.IsDefined(gapModel))
            throw new ArgumentOutOfRangeException(nameof(gapModel), gapModel, "Unknown MSA gap model.");
        var rows = new string[alignedRows.Count];
        for (int i = 0; i < rows.Length; i++)
        {
            ArgumentNullException.ThrowIfNull(alignedRows[i]);
            if (alignedRows[i].Length != alignedRows[0].Length)
                throw new ArgumentException("All aligned rows must have the same length.", nameof(alignedRows));
            rows[i] = alignedRows[i].ToUpperInvariant();
        }

        var ps = PairScoring.Create(scoring, rows.Select(r => r.Replace("-", "")).ToArray());
        return (double)SumOfPairsScoreLong(rows, ps, gapModel) / ps.Scale;
    }

    private static long AffineSumOfPairsScoreCore(IReadOnlyList<string> rows, PairScoring scoring)
    {
        long total = 0;
        for (int i = 0; i < rows.Count; i++)
            for (int j = i + 1; j < rows.Count; j++)
                total += NaturalAffinePairScore(rows[i], rows[j], scoring);

        return total;
    }

    /// <summary>Natural affine score of the projection of two rows (gap–gap columns skipped, one opening per maximal run).</summary>
    private static long NaturalAffinePairScore(string a, string b, PairScoring scoring)
    {
        long score = 0;
        int run = 0; // 0 = last projected column was residue/residue (or none), 1 = gap in a, 2 = gap in b
        int length = Math.Min(a.Length, b.Length);
        for (int c = 0; c < length; c++)
        {
            char x = a[c], y = b[c];
            if (x == '-' && y == '-')
                continue;

            if (x == '-' || y == '-')
            {
                int kind = x == '-' ? 1 : 2;
                if (kind != run)
                    score += scoring.GapOpen;
                score += scoring.GapExtend;
                run = kind;
            }
            else
            {
                score += scoring.Sub(x, y);
                run = 0;
            }
        }

        return score;
    }


    #endregion
}
