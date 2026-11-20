namespace Seqeron.Genomics.Alignment;

// =====================================================================================================================
// The FASTA heuristics of the T-Coffee 13.41 program (audit A2-2; util_dp_fasta_nw.c of the Debian source package t-coffee
// 13.45.0.4846264+really13.41.0.28bdc39+dfsg-1):
//
//   evaluate_diagonals_with_ktup  : the two sequences (lower-cased majority consensus of the rows) are hashed by k-tuples of identical letters; every
//                                   pair of equal k-tuples scores the diagonal (position2 - position1) with the sum of the matrix costs of its ktup
//                                   residue pairs; the diagonals are then sorted by score (stable: glibc qsort merges).
//   flag_diagonals / extract_N_diag: the diagonals the banded DP may use: z-score threshold (`-diag_threshold` 0 -> max z-score / 25) or the N best
//                                   (`-ndiag`), always the corner diagonals 1, l1 + l2 - 1 and the main diagonal l2 (very_fast: 20 best with a window of 5).
//   make_fasta_gotoh_pair_wise    : Gotoh DP restricted to the flagged diagonals (the arrays hold the program's initial "gap" scores outside the band).
//
// Used by fast_pair (-method fast_pair), -distance_matrix_mode very_fast (-quicktree) and geometric, and -dp_mode fasta_pair_wise / cfasta_pair_wise /
// very_fast_pair_wise on profiles. All of it is integer arithmetic (the z-score threshold uses the program's int / double / float mix).
// =====================================================================================================================

public static partial class SequenceAligner
{
    private const int TcUndefinedDiagonal = -99999999;

    /// <summary>
    /// <c>evaluate_diagonals_with_ktup</c>: score and hit count of every diagonal (number d = position2 - position1 + l1, 0-based positions; rows 1..l1+l2-1 sorted by score, ascending,
    /// stable; row 0 stays {0, 0, 0}). <paramref name="cost"/> is the program's <c>get_dp_cost</c> on 0-based profile columns.
    /// </summary>
    internal static int[][] TcDiagonalsByKtup(string consensus1, string consensus2, int ktup, Func<int, int, int> cost)
    {
        int l1 = consensus1.Length, l2 = consensus2.Length;
        RequireTcPairDpCells("evaluate_diagonals_with_ktup", "keeps a 3-int row (about 48 bytes) per diagonal (len1 + len2 + 1 of them)", 48, (long)l1 + l2 + 1,
            string.Create(System.Globalization.CultureInfo.InvariantCulture, $"{l1:N0} x {l2:N0} positions"));
        int nDiag = l1 + l2 - 1;
        var diag = new int[nDiag + 2][];
        for (int a = 0; a < nDiag + 2; a++) diag[a] = new[] { a, 0, 0 };
        if (ktup >= 1)
        {
            var h1 = new Dictionary<string, List<int>>(StringComparer.Ordinal);
            for (int p = 0; p + ktup <= l1; p++)
            {
                string key = consensus1.Substring(p, ktup);
                if (!h1.TryGetValue(key, out var list)) h1[key] = list = new List<int>();
                list.Add(p);
            }

            var h2 = new Dictionary<string, List<int>>(StringComparer.Ordinal);
            for (int p = 0; p + ktup <= l2; p++)
            {
                string key = consensus2.Substring(p, ktup);
                if (!h2.TryGetValue(key, out var list)) h2[key] = list = new List<int>();
                list.Add(p);
            }

            foreach (var (key, offsets1) in h1)
            {
                if (!h2.TryGetValue(key, out var offsets2)) continue;
                foreach (int p1 in offsets1)
                {
                    AlignmentCancellation.ThrowIfRequested();       // A5-2: one poll per k-tuple occurrence (a repeat gives offsets1 x offsets2 pairs)
                    foreach (int p2 in offsets2)
                    {
                        var row = diag[p2 - p1 + l1];
                        row[2]++;
                        int score = 0;
                        for (int d = 0; d < ktup; d++) score += cost(p1 + d, p2 + d);
                        row[1] += score;
                    }
                }
            }
        }

        TcSortDiagonals(diag, nDiag);
        return diag;
    }

    /// <summary><c>sort_int(diag + 1, ..., 1, 0, n_diag - 1)</c>: the rows 1..n_diag by score, ascending (glibc qsort = stable merge sort).</summary>
    internal static void TcSortDiagonals(int[][] diag, int nDiag)
    {
        var rows = new int[nDiag][];
        Array.Copy(diag, 1, rows, 0, nDiag);
        var sorted = rows.OrderBy(r => r[1]).ToArray();
        Array.Copy(sorted, 0, diag, 1, nDiag);
    }

    /// <summary>
    /// <c>flag_diagonals(l1, l2, sorted_diag, T, 0)</c> with the z-score rule: <c>mean</c> over the n_diag + 1 rows (row 0 included), <c>sd</c> = (int) max(sqrt(sum of squared deviations of the
    /// rows from (int) mean) / n, 1); T = 0 means (best z-score) / 25 (stored in a float); the diagonals from the best downwards while z &gt; T, plus 1, l2 and l1 + l2 - 1. Returns the program's list: [0] = count, [1..] ascending.
    /// </summary>
    internal static int[] TcFlagDiagonals(int l1, int l2, int[][] diag, float threshold)
    {
        int nDiag = l1 + l2 - 1;
        int len = nDiag + 1;
        double sum = 0;
        for (int a = 0; a < len; a++) sum += diag[a][1];
        double mean = sum / len;
        int imean = (int)mean;
        double c = 0;
        for (int a = 0; a < len; a++)
        {
            int v = diag[a][1];
            if (imean - v != 0) c += ((double)imean - v) * ((double)imean - v);
        }

        c = Math.Sqrt(c) / len;
        double sd = (int)Math.Max(c, 1);
        float t = threshold;
        if (t == 0) t = (float)(((double)diag[nDiag][1] - mean) / sd / 25);

        var slopes = new int[nDiag + 1];
        for (int a = nDiag; a > 0; a--)
        {
            if (((double)diag[a][1] - mean) / sd > t) slopes[Math.Min(nDiag, Math.Max(1, diag[a][0]))] = 1;
            else break;
        }

        slopes[1] = 1;
        slopes[l1 + l2 - 1] = 1;
        slopes[l2] = 1;
        var list = new int[l1 + l2 + 1];
        for (int a = 0; a <= l1 + l2 - 1; a++)
            if (slopes[a] != 0) list[++list[0]] = a;
        return list;
    }

    /// <summary><c>extract_N_diag(l1, l2, sorted_diag, n_chosen_diag, window)</c>: the best n_chosen_diag diagonals (each widened by window on both sides) plus 1, l1 + l2 - 1 and l1, l2 (widened as well).</summary>
    internal static int[] TcExtractNDiagonals(int l1, int l2, int[][] diag, int nChosen, int window)
    {
        int nDiag = l1 + l2 - 1;
        var slopes = new int[nDiag + 1];
        for (int a = nDiag; a > 0 && a > nDiag - nChosen; a--)
        {
            int cur = diag[a][0];
            int up = Math.Max(1, cur - window), low = Math.Min(nDiag, cur + window);
            for (int b = up; b <= low; b++) slopes[b] = 1;
        }

        foreach (int centre in new[] { 1, l1 + l2 - 1, l1, l2 })
        {
            int up = Math.Max(1, centre - window), low = Math.Min(nDiag, centre + window);
            for (int a = up; a <= low; a++) slopes[a] = 1;
        }

        var list = new int[l1 + l2 + 1];
        for (int a = 0; a <= l1 + l2 - 1; a++)
            if (slopes[a] != 0) list[++list[0]] = a;
        return list;
    }

    /// <summary>
    /// <c>make_fasta_gotoh_pair_wise</c>: the Gotoh DP over the flagged diagonals (maximising) and its traceback; <paramref name="cost"/> is <c>get_dp_cost</c> on 0-based columns, gop / gep are the
    /// penalties already multiplied by SCORE_K, <paramref name="tgMode"/> the terminal gap rule (0 = open + extend, 1 = extend only, 2 = none). Returns the presence masks of the two rows (true = column holds a residue)
    /// and the score of the last cell.
    /// </summary>
    internal static (List<bool> Present1, List<bool> Present2, int Score) TcMakeFastaGotoh(int l1, int l2, int[] diag, Func<int, int, int> cost, int gop, int gep, int tgMode)
    {
        int nDiag = diag[0];
        int rows = l1 + l2 + 1, cols = nDiag + 2; // the program allocates l1 + l2 + 1 rows: the traceback may step through zero rows below the matrix (diag[nDiag + 1] reads 0)
        RequireTcPairDpCells("fasta_pair_wise (the banded Gotoh of fast_pair)", "keeps five int matrices of (len1 + len2 + 1) rows x (flagged diagonals + 2) columns", 20, ((long)l1 + l2 + 1) * ((long)nDiag + 2),
            string.Create(System.Globalization.CultureInfo.InvariantCulture, $"{l1:N0} x {l2:N0} positions with {nDiag:N0} flagged diagonals"), "use the slow_pair or lalign_id_pair primary library (linear space), or -ndiag / -diag_threshold to flag fewer diagonals");
        var C = new int[rows][]; var D = new int[rows][]; var I = new int[rows][]; var LD = new int[rows][]; var trace = new int[rows][];
        for (int r = 0; r < rows; r++)
        {
            if ((r & 255) == 0) AlignmentCancellation.ThrowIfRequested();
            C[r] = new int[cols]; D[r] = new int[cols]; I[r] = new int[cols]; LD[r] = new int[cols]; trace[r] = new int[cols];
        }

        int termGop = tgMode == 0 ? gop : 0, termGep = tgMode == 2 ? 0 : gep;
        int j;
        for (j = 1; j <= nDiag; j++)
        {
            int lGop = termGop;
            int lGep = termGep;
            if (diag[j] - l1 < 0)
            {
                trace[0][j] = TcUndefinedDiagonal;
                continue;
            }

            C[0][j] = (diag[j] - l1) * lGep + lGop;
            D[0][j] = (diag[j] - l1) * lGep + lGop + gop;
        }

        D[0][j] = D[0][j - 1] + gep;

        int t = termGop;
        int score = 0, lastI = 0, lastJ = 0;
        for (int i = 1; i <= l1; i++)
        {
            AlignmentCancellation.ThrowIfRequested();       // A5-2: one poll per row
            C[i][0] = C[i][nDiag + 1] = t = t + termGep;
            I[i][0] = D[i][nDiag + 1] = t + gop;
            for (j = 1; j <= nDiag; j++)
            {
                C[i][j] = C[i][0];
                D[i][j] = I[i][j] = I[i][0];
            }

            int eg = 0;
            for (j = 1; j <= nDiag; j++)
            {
                int posJ = diag[j] - l1 + i;
                if (posJ <= 0 || posJ > l2)
                {
                    trace[i][j] = TcUndefinedDiagonal;
                    continue;
                }

                int sub = cost(i - 1, posJ - 1);

                int lGop = i == l1 ? termGop : gop;
                int lGep = i == l1 ? termGep : gep;
                int len = j == 1 ? 0 : diag[j] - diag[j - 1];
                if (I[i][j - 1] > C[i][j - 1] + lGop) eg++;
                else eg = 1;
                I[i][j] = Math.Max(I[i][j - 1], C[i][j - 1] + lGop) + len * lGep;

                lGop = posJ == l2 ? termGop : gop;
                lGep = posJ == l2 ? termGep : gep;
                len = j == nDiag ? 0 : diag[j + 1] - diag[j];
                int deltaI = i - len > 0 ? i - len : 0;
                if (D[deltaI][j + 1] > C[deltaI][j + 1] + lGop) LD[i][j] = LD[deltaI][j + 1] + 1;
                else LD[i][j] = 1;
                D[i][j] = Math.Max(D[deltaI][j + 1], C[deltaI][j + 1] + lGop) + len * lGep;

                int su = C[i - 1][j] + sub, inn = I[i][j], de = D[i][j], tr;
                if (su >= inn && su >= de) { score = su; tr = 0; }
                else if (inn >= de) { score = inn; tr = -eg; }
                else { score = de; tr = LD[i][j]; }
                trace[i][j] = tr;
                C[i][j] = score;
                lastI = i;
                lastJ = j;
            }
        }

        int ii = lastI, jj = lastJ;
        var al0 = new List<bool>(); var al1 = new List<bool>();
        int match1 = 0, match2 = 0;
        while (!(match1 == l1 && match2 == l2))
        {
            int k;
            if (match1 == l1)
            {
                int len = l2 - match2;
                for (int a = 0; a < len; a++) { al0.Add(false); al1.Add(true); match2++; }
                break;
            }
            else if (match2 == l2)
            {
                int len = l1 - match1;
                for (int a = 0; a < len; a++) { al0.Add(true); al1.Add(false); match1++; }
                break;
            }
            else k = trace[ii][jj];

            if (k == TcUndefinedDiagonal)
                throw new TCoffeeAbortedException("segmentation fault in make_fasta_gotoh_pair_wise: the traceback enters a cell outside the band (trace = UNDEFINED indexes the diagonal list with 99999999)", 11);
            if (k == 0)
            {
                al0.Add(true); al1.Add(true);
                ii--; match1++; match2++;
            }
            else if (k > 0)
            {
                int len = diag[jj + k] - diag[jj];
                for (int a = 0; a < len; a++)
                {
                    if (match1 == l1) break;
                    al0.Add(true); al1.Add(false);
                    match1++;
                }

                ii -= len;
                jj += k;
            }
            else
            {
                k = -k;
                int len = diag[jj] - diag[jj - k];
                for (int a = 0; a < len; a++)
                {
                    if (match2 == l2) break;
                    al0.Add(false); al1.Add(true);
                    match2++;
                }

                jj -= k;
            }
        }

        al0.Reverse(); al1.Reverse();
        return (al0, al1, score);
    }

    /// <summary>Applies the presence masks of <see cref="TcMakeFastaGotoh"/> to two sequences.</summary>
    private static (string, string) TcMasksToRows(string a, string b, List<bool> m0, List<bool> m1)
    {
        var ra = new char[m0.Count]; var rb = new char[m0.Count];
        int pa = 0, pb = 0;
        for (int c = 0; c < m0.Count; c++)
        {
            ra[c] = m0[c] ? a[pa++] : '-';
            rb[c] = m1[c] ? b[pb++] : '-';
        }

        return (new string(ra), new string(rb));
    }

    /// <summary>
    /// <c>fasta_gotoh_pair_wise</c> of two single sequences with the program's <c>slow_get_dp_cost</c> (matrix score x SCORE_K minus SCORE_K x nomatch) and the k-tuple diagonals of the lower-cased letters.
    /// </summary>
    internal static (string RowA, string RowB) TcFastaPairRows(string a, string b, TcMatrix matrix, int gop, int gep, int ktup, int nomatch, int fastaStep, float diagonalThreshold, int tgMode)
    {
        int l1 = a.Length, l2 = b.Length;
        Func<int, int, int> cost = (i, j) => matrix.Score(a[i], b[j]) * TcScoreK - TcScoreK * nomatch;
        var diag = TcDiagonalsByKtup(a.ToLowerInvariant(), b.ToLowerInvariant(), ktup, cost);
        int[] list = fastaStep == 0 ? TcFlagDiagonals(l1, l2, diag, diagonalThreshold) : TcExtractNDiagonals(l1, l2, diag, fastaStep, 0);
        var (m0, m1, _) = TcMakeFastaGotoh(l1, l2, list, cost, gop * TcScoreK, gep * TcScoreK, tgMode);
        return TcMasksToRows(a, b, m0, m1);
    }
}
