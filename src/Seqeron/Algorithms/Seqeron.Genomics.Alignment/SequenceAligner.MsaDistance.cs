namespace Seqeron.Genomics.Alignment;

public static partial class SequenceAligner
{
    #region Guide-tree distances of the progressive family: exact identity without a traceback, k-tuple estimate for large inputs (B09 WP78, R14-16, F362)

    /// <summary>
    /// The number of identical aligned residue pairs of the optimal linear-gap global alignment that <c>GlobalAlignCore</c> returns (traceback ties diagonal &gt; up &gt; left), computed by the forward
    /// recurrence alone in O(n) memory: at every cell the move that the traceback would take from that cell is the one that produced the maximum, with the same preference, so the identity count of the preferred
    /// path to a cell is its predecessor's count plus one for an identical diagonal step. The count of the cell (m, n) is the count of the traceback path; no matrix, no aligned strings. A residue pair counts when the two
    /// characters are equal and not '-' (as <see cref="PairwiseIdentityDistance"/> always counted them).
    /// </summary>
    internal static int LinearIdentityCount(string a, string b, PairScoring scoring)
    {
        int m = a.Length;
        int n = b.Length;
        scoring.EnsureRange(m, n);
        AlignmentCancellation.Poll(default);

        var upGap = new int[n + 1];
        for (int j = 0; j <= n; j++)
            upGap[j] = scoring.XExt(j, n);

        var prevScore = new int[n + 1];
        var prevCount = new int[n + 1];
        var curScore = new int[n + 1];
        var curCount = new int[n + 1];
        int firstRowGap = scoring.YExt(0, m);
        for (int j = 0; j <= n; j++)
            prevScore[j] = j * firstRowGap;

        for (int i = 1; i <= m; i++)
        {
            if (i % 100 == 0)
                AlignmentCancellation.Poll(default);

            char c1 = a[i - 1];
            int leftGap = scoring.YExt(i, m);
            curScore[0] = i * upGap[0];
            curCount[0] = 0;
            for (int j = 1; j <= n; j++)
            {
                char c2 = b[j - 1];
                int diag = prevScore[j - 1] + scoring.Sub(c1, c2);
                int up = prevScore[j] + upGap[j];
                int left = curScore[j - 1] + leftGap;
                int best = Math.Max(diag, Math.Max(up, left));
                curScore[j] = best;
                if (best == diag)
                    curCount[j] = prevCount[j - 1] + (c1 != '-' && c1 == c2 ? 1 : 0);
                else
                    curCount[j] = best == up ? prevCount[j] : curCount[j - 1];
            }

            (prevScore, curScore) = (curScore, prevScore);
            (prevCount, curCount) = (curCount, prevCount);
        }

        return prevCount[n];
    }

    /// <summary>
    /// Above this many dynamic-programming cells (the sum of len1 · len2 over all pairs of sequences) the all-pairs distance stage of the progressive family uses the k-tuple estimate
    /// (<see cref="KTupleDistance"/>) instead of one global alignment per pair; see <see cref="MsaOptions.FastDistanceCells"/>.
    /// </summary>
    internal const long DefaultFastDistanceCells = 400_000_000; // the default of MsaOptions.FastDistanceCells

    internal const int KTupleDnaLength = 8;
    internal const int KTupleProteinLength = 4;

    /// <summary>The sorted k-mer codes of a sequence (5 bits per character, upper case) for the k-tuple distance; null when it is shorter than 2 k.</summary>
    internal static ulong[]? KTupleCodes(string s, int k)
    {
        if (s.Length < 2 * k)
            return null;
        var codes = new ulong[s.Length - k + 1];
        ulong mask = (1UL << (5 * k)) - 1;
        ulong window = 0;
        for (int i = 0; i < s.Length; i++)
        {
            window = ((window << 5) | ((uint)char.ToUpperInvariant(s[i]) & 31u)) & mask;
            if (i >= k - 1)
                codes[i - k + 1] = window;
        }

        Array.Sort(codes);
        return codes;
    }

    /// <summary>
    /// k-tuple distance of two sequences: the fraction F of the k-mers of the shorter sequence that the longer one shares (multiset intersection, MUSCLE's k-mer similarity, Edgar 2004, Nucleic Acids Res 32:1792)
    /// turned into the identity scale of the ClustalW distance: a k-mer survives with probability p<sup>k</sup> when a fraction p of the residues is conserved, so p = F<sup>1/k</sup> and d = 1 − p
    /// (k = 8 for nucleotide scoring, 4 for a protein matrix; sorted code arrays, one merge per pair, O(len1 + len2)).
    /// </summary>
    internal static double KTupleDistance(ulong[] a, ulong[] b, int lengthA, int lengthB, int k)
    {
        int shared = 0, i = 0, j = 0;
        while (i < a.Length && j < b.Length)
        {
            if (a[i] < b[j])
            {
                i++;
            }
            else if (a[i] > b[j])
            {
                j++;
            }
            else
            {
                ulong value = a[i];
                int runA = 0, runB = 0;
                while (i < a.Length && a[i] == value)
                {
                    i++;
                    runA++;
                }

                while (j < b.Length && b[j] == value)
                {
                    j++;
                    runB++;
                }

                shared += Math.Min(runA, runB);
            }
        }

        int windows = Math.Min(lengthA, lengthB) - k + 1;
        double fraction = Math.Min(1.0, (double)shared / windows);
        return 1.0 - Math.Pow(fraction, 1.0 / k);
    }

    /// <summary>True when the all-pairs alignment work (the sum of len1 · len2 over the pairs) is above <paramref name="fastDistanceCells"/> (0 or negative: never).</summary>
    private static bool UseKTupleDistances(string[] rawSeqs, long fastDistanceCells)
    {
        if (fastDistanceCells <= 0)
            return false;
        double cells = 0;
        long total = 0, squares = 0;
        foreach (var s in rawSeqs)
        {
            total += s.Length;
            squares += (long)s.Length * s.Length;
        }

        cells = ((double)total * total - squares) / 2.0;
        return cells > fastDistanceCells;
    }

    #endregion
}
