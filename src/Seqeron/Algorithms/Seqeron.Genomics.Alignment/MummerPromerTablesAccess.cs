namespace Seqeron.Genomics.Alignment;

internal static partial class MummerPromerTables
{
    /// <summary>
    /// <c>MATCH_SCORE [matrix] [i] [j]</c> read as the C array <c>int [4][26][26]</c> is laid out in memory, so that an index outside 0 - 25 (the character <c>*</c> of a translation is 23 below <c>A</c>) reads the neighbouring row exactly as
    /// <c>show-aligns</c> does (the flat position is always inside the four matrices for a matrix type of 1 - 3).
    /// </summary>
    public static int FlatMatchScore(int matrix, int i, int j)
    {
        int flat = ((matrix * 26) + i) * 26 + j;
        int block = flat / 676, within = flat % 676;
        if (block == 0)
        {
            int row = within / 26, col = within % 26;
            return row == col && row is 0 or 2 or 6 or 19 ? 3 : -7;   // A C G T
        }

        return Blosum[block - 1][within];
    }
}
