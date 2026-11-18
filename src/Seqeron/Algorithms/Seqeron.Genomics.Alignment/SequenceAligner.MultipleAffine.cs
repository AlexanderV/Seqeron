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

    /// <summary>
    /// Largest profile-profile cell count of the affine Gotoh program (one traceback byte per cell and O(width) work rows).
    /// </summary>
    public const long AffineProfileMatrixMaxCells = 200_000_000;

    private const long ProfileNegInf = long.MinValue / 4;

    private static MsaGapModel ResolveGapModel(MsaOptions? options)
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
        return (double)ComputeSumOfPairsScore(rows, ps, gapModel) / ps.Scale;
    }

    private static int AffineSumOfPairsScoreCore(IReadOnlyList<string> rows, PairScoring scoring)
    {
        long total = 0;
        for (int i = 0; i < rows.Count; i++)
            for (int j = i + 1; j < rows.Count; j++)
                total += NaturalAffinePairScore(rows[i], rows[j], scoring);

        return checked((int)total);
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

    private static void RequireAffineProfileMatrixFits(int m, int n)
    {
        long cells = (long)(m + 1) * (n + 1);
        if (cells > AffineProfileMatrixMaxCells)
            throw new ArgumentException(
                $"The affine profile-profile dynamic program needs one traceback byte per cell: {m} x {n} alignment columns = {cells:N0} cells exceeds the limit of {AffineProfileMatrixMaxCells:N0}. " +
                "Use MultipleAlign (exact center-star, linear-space pairwise stage) for long sequences.");
    }

    /// <summary>
    /// Aligns two profiles with the affine (Gotoh three-state) recurrences over columns and the occupancy-scaled sum-of-pairs
    /// column scores described in the region comment. Existing columns are never edited ("once a gap, always a gap"); all arithmetic is
    /// exact (integers), ties are broken M &gt; own gap state &gt; other gap state exactly as in <c>GlobalAlignAffine</c>.
    /// </summary>
    private static Profile AlignProfilesAffine(Profile p1, Profile p2, PairScoring scoring)
    {
        int m = p1.Width;
        int n = p2.Width;
        RequireAffineProfileMatrixFits(m, n);

        int rows1 = p1.Rows.Count;
        int rows2 = p2.Rows.Count;
        var alphabet = new List<char>();
        foreach (var profile in new[] { p1, p2 })
            foreach (var row in profile.Rows)
                foreach (char ch in row)
                    if (ch != '-' && !alphabet.Contains(ch))
                        alphabet.Add(ch);
        var (counts1, gaps1) = ProfileColumnCounts(p1, alphabet);
        var (counts2, gaps2) = ProfileColumnCounts(p2, alphabet);
        int sym = alphabet.Count;

        // Residue-pair substitution table over the profile alphabet.
        var pair = new long[sym, sym];
        for (int c = 0; c < sym; c++)
            for (int d = 0; d < sym; d++)
                pair[c, d] = scoring.Sub(alphabet[c], alphabet[d]);

        long open = scoring.GapOpen;
        long ext = scoring.GapExtend;

        // Per-column residue counts (occupancy).
        var res1 = new int[m];
        for (int i = 0; i < m; i++) res1[i] = rows1 - gaps1[i];
        var res2 = new int[n];
        for (int j = 0; j < n; j++) res2[j] = rows2 - gaps2[j];

        long Sub(int i, int j)
        {
            long s = 0;
            for (int c = 0; c < sym; c++)
            {
                int ca = counts1[i * sym + c];
                if (ca == 0) continue;
                for (int d = 0; d < sym; d++)
                {
                    int cb = counts2[j * sym + d];
                    if (cb != 0) s += (long)ca * cb * pair[c, d];
                }
            }

            // Residue against an existing gap: extension only; gap against gap: neutral.
            s += ext * ((long)gaps1[i] * res2[j] + (long)res1[i] * gaps2[j]);
            return s;
        }

        // Traceback pointers: bits 0-1 predecessor state of M, bits 2-3 of X, bits 4-5 of Y (0 = M, 1 = X, 2 = Y).
        var tb = new byte[(long)(m + 1) * (n + 1)];
        long Cell(int i, int j) => (long)i * (n + 1) + j;

        var prevM = new long[n + 1];
        var prevX = new long[n + 1];
        var prevY = new long[n + 1];
        var curM = new long[n + 1];
        var curX = new long[n + 1];
        var curY = new long[n + 1];

        for (int i = 0; i <= m; i++)
        {
            for (int j = 0; j <= n; j++)
            {
                if (i == 0 && j == 0)
                {
                    curM[0] = 0;
                    curX[0] = ProfileNegInf;
                    curY[0] = ProfileNegInf;
                    continue;
                }

                byte pointers = 0;

                if (i > 0 && j > 0)
                {
                    // M(i,j) = best of M/X/Y(i-1,j-1) + column score; ties M > X > Y.
                    long bm = prevM[j - 1];
                    int from = 0;
                    if (prevX[j - 1] > bm) { bm = prevX[j - 1]; from = 1; }
                    if (prevY[j - 1] > bm) { bm = prevY[j - 1]; from = 2; }
                    curM[j] = bm <= ProfileNegInf ? ProfileNegInf : bm + Sub(i - 1, j - 1);
                    pointers |= (byte)from;
                }
                else
                {
                    curM[j] = ProfileNegInf;
                }

                if (i > 0)
                {
                    // X(i,j): column i of A against a new gap column of B.
                    long residues = (long)res1[i - 1] * rows2;
                    long e = ext * residues;
                    long oe = (open + ext) * residues;
                    long cm = prevM[j] + oe;
                    long cx = prevX[j] + e;
                    long cy = prevY[j] + oe;
                    long bx = cm;
                    int from = 0;
                    if (cx > bx) { bx = cx; from = 1; }
                    if (cy > bx) { bx = cy; from = 2; }
                    curX[j] = bx < ProfileNegInf ? ProfileNegInf : bx;
                    pointers |= (byte)(from << 2);
                }
                else
                {
                    curX[j] = ProfileNegInf;
                }

                if (j > 0)
                {
                    // Y(i,j): column j of B against a new gap column of A; priority M > Y (own) > X.
                    long residues = (long)res2[j - 1] * rows1;
                    long e = ext * residues;
                    long oe = (open + ext) * residues;
                    long cm = curM[j - 1] + oe;
                    long cy = curY[j - 1] + e;
                    long cx = curX[j - 1] + oe;
                    long by = cm;
                    int from = 0;
                    if (cy > by) { by = cy; from = 2; }
                    if (cx > by) { by = cx; from = 1; }
                    curY[j] = by < ProfileNegInf ? ProfileNegInf : by;
                    pointers |= (byte)(from << 4);
                }
                else
                {
                    curY[j] = ProfileNegInf;
                }

                tb[Cell(i, j)] = pointers;
            }

            (prevM, curM) = (curM, prevM);
            (prevX, curX) = (curX, prevX);
            (prevY, curY) = (curY, prevY);
        }

        // After the last swap the row i = m values are in prev*.
        int state = 0;
        long best = prevM[n];
        if (prevX[n] > best) { best = prevX[n]; state = 1; }
        if (prevY[n] > best) state = 2;

        var cols1 = new List<int>(m + n);
        var cols2 = new List<int>(m + n);
        int ii = m, jj = n;
        while (ii > 0 || jj > 0)
        {
            byte pointers = tb[Cell(ii, jj)];
            if (state == 0)
            {
                cols1.Add(ii - 1);
                cols2.Add(jj - 1);
                state = pointers & 3;
                ii--; jj--;
            }
            else if (state == 1)
            {
                cols1.Add(ii - 1);
                cols2.Add(-1);
                state = (pointers >> 2) & 3;
                ii--;
            }
            else
            {
                cols1.Add(-1);
                cols2.Add(jj - 1);
                state = (pointers >> 4) & 3;
                jj--;
            }
        }

        cols1.Reverse();
        cols2.Reverse();

        var merged = new Profile();
        merged.SequenceIndices.AddRange(p1.SequenceIndices);
        merged.SequenceIndices.AddRange(p2.SequenceIndices);
        AppendReprojectedRows(merged, p1, cols1);
        AppendReprojectedRows(merged, p2, cols2);
        return merged;
    }

    /// <summary>
    /// Per-column residue counts of a profile over the given alphabet: <c>counts[col * alphabet + symbol]</c> and the number of gap
    /// characters per column.
    /// </summary>
    private static (int[] Counts, int[] Gaps) ProfileColumnCounts(Profile profile, List<char> alphabet)
    {
        int width = profile.Width;
        int sym = alphabet.Count;
        var counts = new int[width * sym];
        var gaps = new int[width];
        foreach (var row in profile.Rows)
        {
            for (int col = 0; col < width; col++)
            {
                char ch = row[col];
                if (ch == '-')
                    gaps[col]++;
                else
                    counts[col * sym + alphabet.IndexOf(ch)]++;
            }
        }

        return (counts, gaps);
    }

    #endregion
}
