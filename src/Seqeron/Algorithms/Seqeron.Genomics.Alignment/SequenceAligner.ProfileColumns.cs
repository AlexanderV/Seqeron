namespace Seqeron.Genomics.Alignment;

public static partial class SequenceAligner
{
    #region Column compositions of a profile (B09 WP78, R14-16, F362)

    // The profile-profile column score of the progressive family is the sum over all cross pairs of rows of s(a, b) (residue pair), GapExtend (residue against gap) or 0 (gap against gap), divided by the number of pairs (the row-by-row sum of the former ProfileColumnScore). The sum
    // only depends on how many rows hold each symbol in the two columns: Σ_x Σ_y c1(x) · c2(y) · s(x, y) + GapExtend · (gaps1 · residues2 + residues1 · gaps2). All terms are integers (exact in a double at any size that
    // fits the integer score model), so the quotient is the very double of the row-by-row sum, at O(nonzero symbols) per cell instead of O(rows1 · rows2).
    private sealed class ProfileColumnTable
    {
        public ProfileColumnTable(Profile profile, Dictionary<char, int> symbolIndex)
        {
            Rows = profile.Rows.Count;
            int width = profile.Width;
            int k = symbolIndex.Count;
            Counts = new int[width][];
            Residues = new int[width];
            for (int col = 0; col < width; col++)
            {
                var counts = new int[k];
                int residues = 0;
                foreach (var row in profile.Rows)
                {
                    char c = row[col];
                    if (c == '-')
                        continue;
                    counts[symbolIndex[c]]++;
                    residues++;
                }

                Counts[col] = counts;
                Residues[col] = residues;
            }
        }

        public int Rows { get; }

        /// <summary>Per column: the number of rows that hold each symbol.</summary>
        public int[][] Counts { get; }

        /// <summary>Per column: the number of rows that hold a residue (not a gap).</summary>
        public int[] Residues { get; }

        /// <summary>Per column of this profile and symbol y: Σ_x count(x) · s(x, y) (filled by <c>Project</c>).</summary>
        public long[][]? Projection { get; private set; }

        /// <summary>Per column: the symbols with a non-zero count.</summary>
        public int[][]? NonZero { get; private set; }

        public void Project(int[][] substitution)
        {
            int k = substitution.Length;
            Projection = new long[Counts.Length][];
            NonZero = new int[Counts.Length][];
            for (int col = 0; col < Counts.Length; col++)
            {
                var counts = Counts[col];
                var nonZero = new List<int>();
                for (int x = 0; x < k; x++)
                {
                    if (counts[x] != 0)
                        nonZero.Add(x);
                }

                var projection = new long[k];
                foreach (int x in nonZero)
                {
                    for (int y = 0; y < k; y++)
                        projection[y] += (long)counts[x] * substitution[x][y];
                }

                Projection[col] = projection;
                NonZero[col] = nonZero.ToArray();
            }
        }
    }

    /// <summary>The profile-profile score of column <paramref name="col1"/> of <paramref name="t1"/> against column <paramref name="col2"/> of <paramref name="t2"/> from their compositions (the double of the row-by-row sum divided by the number of pairs).</summary>
    private static double ProfileColumnScoreFast(ProfileColumnTable t1, int col1, ProfileColumnTable t2, int col2, int gapExtend)
    {
        var projection = t1.Projection![col1];
        var counts2 = t2.Counts[col2];
        long total = 0;
        foreach (int y in t2.NonZero![col2])
            total += projection[y] * counts2[y];
        long gaps1 = t1.Rows - t1.Residues[col1];
        long gaps2 = t2.Rows - t2.Residues[col2];
        total += (long)gapExtend * (gaps1 * t2.Residues[col2] + (long)t1.Residues[col1] * gaps2);
        return (double)total / ((double)t1.Rows * t2.Rows);
    }

    /// <summary>The score of column <paramref name="col"/> of <paramref name="t"/> against an all-gap column of a profile of <paramref name="otherRows"/> rows (residue against gap costs GapExtend, gap against gap nothing).</summary>
    private static double ProfileColumnVsGapFast(ProfileColumnTable t, int col, int otherRows, int gapExtend)
        => (double)((long)gapExtend * otherRows * t.Residues[col]) / ((double)t.Rows * otherRows);

    #endregion
}
