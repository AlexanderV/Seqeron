namespace Seqeron.Genomics.Alignment;

// =====================================================================================================================
// The dynamic-programming modes of the progressive stage of the T-Coffee 13.41 program other than the default sparse linked DP (audit A2-5):
// `-dp_mode gotoh_pair_wise | myers_miller_pair_wise | fasta_pair_wise | cfasta_pair_wise | very_fast_pair_wise` (+ the aliases `ssec_pair_wise`, `default`,
// `gotoh_pair_wise_test`) together with the library-driven column scores of `-extend_mode` (choose_extension_mode, evaluate.c):
//
//   very_fast_triplet (default)  residue_pair_extended_list + fast_get_dp_cost_2 : the (at most -max_n_pair) pairs of sequences of the two profile columns with the
//                                highest identity of the distance matrix, the extended (triplet) score of every pair averaged over the pairs with two residues, minus SCORE_K x nomatch
//   slow_triplet                 residue_pair_extended_list + slow_get_dp_cost   : every pair of sequences
//   triplet                      residue_pair_extended_list + get_dp_cost        : slow_get_dp_cost, nomatch subtracted a second time
//   fast_triplet                 residue_pair_extended_list + fast_get_dp_cost   : hash of the library partners of the two columns, normalised by max_ext_value
//
// util_dp_gotoh_nw.c gotoh_pair_wise (insertions / deletions in the first and last row / column and in the first row / column pay the terminal gap rule `-tg_mode`),
// util_dp_mm_nw.c myers_miller_pair_wise (the Myers-Miller divide and conquer of slow_pair on profile columns) and util_dp_fasta_nw.c (the banded heuristics, see
// SequenceAligner.TCoffeeFasta.cs) run on those scores. The default linked_pair_wise has its own sparse DP (SequenceAligner.TCoffee.cs).
// =====================================================================================================================

/// <summary><c>-dp_mode</c> of the program: the dynamic program the progressive stage aligns two profiles with.</summary>
public enum TCoffeeDpMode
{
    /// <summary><c>linked_pair_wise</c> (default): the sparse DP over the scored column pairs of the library (<c>cl2pair_list_ext</c>).</summary>
    Linked = 0,

    /// <summary><c>gotoh_pair_wise</c> (also <c>ssec_pair_wise</c>): affine-gap Needleman-Wunsch over the column scores of <c>-extend_mode</c>.</summary>
    Gotoh = 1,

    /// <summary><c>gotoh_pair_wise_test</c>: <c>gotoh_pair_wise</c> with the gap penalties forced to open 0 / extend -1 (x SCORE_K for the open only) and the tie order insertion, substitution, deletion.</summary>
    GotohTest = 2,

    /// <summary><c>myers_miller_pair_wise</c> (also <c>default</c>): the linear-space Myers-Miller DP on the column scores.</summary>
    MyersMiller = 3,

    /// <summary><c>fasta_pair_wise</c>: the Gotoh DP restricted to the diagonals with a z-score above the threshold (<c>-diag_threshold</c>, 0 = best z / 25) or the <c>-ndiag</c> best ones; diagonals scored by the library.</summary>
    Fasta = 4,

    /// <summary><c>cfasta_pair_wise</c>: <c>fasta_pair_wise</c> repeated with 20 + n x step diagonals until the score stops changing.</summary>
    CFasta = 5,

    /// <summary><c>very_fast_pair_wise</c>: the 20 best diagonals with a shadow of 5 on each side.</summary>
    VeryFast = 6,

    /// <summary><c>biphasic_pair_wise</c>: five states, two affine gap models (opening x 2 / extension / 2 and opening / 2 / extension x 2) on the column scores.</summary>
    Biphasic = 7,

    /// <summary><c>glocal_pair_wise</c>: four states (match, deletion, insertion, an unmatched pair that costs the extension) on the column scores.</summary>
    Glocal = 8,

    /// <summary>Any other name: the program refuses it (see <see cref="TCoffeeProgramOptions.DpModeText"/>).</summary>
    Other = 9,

    /// <summary><c>gotoh_pair_wise_lgp</c>: the three-state Gotoh DP with position specific gap opening / closing penalties derived from the gaps of the two profiles (<c>aln2local_penalties</c>).</summary>
    GotohLgp = 10,

    /// <summary><c>ktup_pair_wise</c> (LOCAL): scores the k-tuple similarity of the two groups and leaves every row as it is (no alignment is produced).</summary>
    Ktup = 11,

    /// <summary><c>precomputed_pair_wise</c>: the identity of two already aligned groups; the program stops when the first rows of the groups differ in length, otherwise the rows stay as they are.</summary>
    Precomputed = 12,

    /// <summary><c>gotoh_pair_wise_sw</c> (LOCAL): the Smith-Waterman variant of the Gotoh DP; the rows keep only the columns of the best local alignment.</summary>
    GotohSw = 13,

    /// <summary>
    /// <c>test_pair_wise</c> (a debugging routine): no alignment, every row of the first group is followed by <c>l1 - n</c> gaps and every row of the second group by <c>l0 - n</c> gaps with <c>n = l0 / 2</c> (<c>l0 &lt; 5</c>) or 5, so that the two blocks of rows
    /// lie side by side from column 0 (the groups overlap). The routine writes the longer rows into buffers sized for the longest input sequence (a heap overflow: the program dies on part of the sets, see <see cref="TCoffeeAbortedException.HeapDependent"/>)
    /// and calls <c>generate_null</c> with a negative size when the second group has fewer than <c>n</c> columns (a null allocation: deterministic segmentation fault).
    /// </summary>
    Test = 14,
}

/// <summary><c>-extend_mode</c> of the program (<c>choose_extension_mode</c>): how the library is used to score a pair of profile columns.</summary>
public enum TCoffeeExtendMode
{
    /// <summary><c>very_fast_triplet</c> (default): <c>residue_pair_extended_list</c> and <c>fast_get_dp_cost_2</c> (the <c>-max_n_pair</c> most similar pairs of sequences).</summary>
    VeryFastTriplet = 0,

    /// <summary><c>fast_triplet</c>: <c>fast_get_dp_cost</c> (hashes of the library partners of the two columns).</summary>
    FastTriplet = 1,

    /// <summary><c>slow_triplet</c>: <c>slow_get_dp_cost</c> (every pair of sequences).</summary>
    SlowTriplet = 2,

    /// <summary><c>triplet</c>: <c>get_dp_cost</c> = <c>slow_get_dp_cost</c> with the nomatch penalty subtracted twice.</summary>
    Triplet = 3,

    /// <summary>
    /// <c>relative_triplet</c>: <c>fast_get_dp_cost_2</c> with <c>residue_pair_relative_extended_list</c>, which sums the weights of the library list of the FIRST residue at the positions of the second residue's list: when that list is longer
    /// the read runs past the array (undefined behaviour of the program: the port raises <see cref="NotSupportedException"/> there and equals the program on every set where no read leaves the array).
    /// </summary>
    RelativeTriplet = 4,

    /// <summary>Any other name: the program refuses or dies (see <see cref="TCoffeeProgramOptions.ExtendModeText"/>); no effect on the default linked DP of the names the program accepts.</summary>
    Other = 5,
}

public static partial class SequenceAligner
{
    private const int TcForbidden = -99999999;

    /// <summary>The parameters of the program's CL that the generic dynamic programs read.</summary>
    private sealed record TcDpSettings(TCoffeeDpMode Mode, TCoffeeExtendMode Extend, int Gop, int Gep, int TgMode, int Ktup, int FastaStep, float DiagonalThreshold, int MaxNPair, int Normalise);

    /// <summary>
    /// The constraint list as the generic dynamic programs see it: the (relaxed, re-weighted) library, the identity matrix of the distance mode (<c>CL->DM->similarity_matrix</c>, <c>null</c> = no distance matrix),
    /// <c>CL->nomatch</c> (rescaled) and <c>CL->max_ext_value</c>.
    /// </summary>
    private sealed class TcDpContext
    {
        public readonly TcIndex Index;
        public readonly int[,]? Sim;
        public readonly TcDpSettings Settings;
        public readonly int NoMatch;
        public readonly int MaxExt;
        private readonly int[][] _hasch;
        private readonly int[][] _hasch1, _hasch2, _n1, _n2, _in1, _in2;

        public TcDpContext(TcIndex index, int[,]? sim, TcDpSettings settings, int noMatch, int maxExt)
        {
            Index = index;
            Sim = sim;
            Settings = settings;
            NoMatch = noMatch;
            MaxExt = maxExt;
            int[][] Make() => index.Len.Select(l => new int[l + 1]).ToArray();
            _hasch = Make();
            _hasch1 = Make(); _hasch2 = Make(); _n1 = Make(); _n2 = Make(); _in1 = Make(); _in2 = Make();
        }

        /// <summary><c>residue_pair_extended_list</c> (s1 &lt; s2): the triplet-extended score of the residue pair, normalised by the sum of the weights of the two residues (scaled by SCORE_K), double arithmetic as in the original.</summary>
        public int ResiduePairExtended(int s1, int r1, int s2, int r2)
        {
            if (r1 <= 0 || r2 <= 0) return 0;
            double score = 0, maxScore = 0, maxVal = 0;
            var own = Index.Edges[s1][r1];
            _hasch[s1][r1] = TcForbidden;
            foreach (var e in own)
            {
                _hasch[e.Seq][e.Res] = e.Weight;
                maxScore += e.Weight;
            }

            foreach (var x in Index.Edges[s2][r2])
            {
                int h = _hasch[x.Seq][x.Res];
                if (h != 0)
                {
                    if (h == TcForbidden)
                    {
                        score += x.Weight;
                        maxScore += x.Weight;
                    }
                    else
                    {
                        double delta = Math.Min(h, x.Weight);
                        score += delta;
                        maxScore -= h;
                        maxScore += delta;
                        maxVal = Math.Max(maxVal, delta);
                    }
                }
                else maxScore += x.Weight;
            }

            maxScore -= _hasch[s2][r2];
            foreach (var e in own) _hasch[e.Seq][e.Res] = 0;
            _hasch[s1][r1] = _hasch[s2][r2] = 0;

            if (maxScore == 0) score = 0;
            else if (Settings.Normalise != 0)
            {
                score = ((score * Settings.Normalise) / maxScore) * TcScoreK;
                if (maxVal > Settings.Normalise) score *= maxVal / (double)Settings.Normalise;
            }

            return (int)score;
        }

        /// <summary>The column score <c>get_dp_cost</c> of the profile columns <paramref name="col1"/> of <paramref name="g1"/> and <paramref name="col2"/> of <paramref name="g2"/> for <c>-extend_mode</c>.</summary>
        /// <summary><c>residue_pair_relative_extended_list</c> (s1 &lt; s2): int arithmetic; the weights of the list of s1 are read at the positions of the list of s2 (a position past the end of the list is the program's undefined read: <see cref="NotSupportedException"/>).</summary>
        public int ResiduePairRelative(int s1, int r1, int s2, int r2)
        {
            if (r1 <= 0 || r2 <= 0) return 0;
            int score = 0, total = 0;
            var own = Index.Edges[s1][r1];
            _hasch[s1][r1] = 100000;
            foreach (var e in own)
            {
                _hasch[e.Seq][e.Res] = e.Weight;
                total += e.Weight;
            }

            var other = Index.Edges[s2][r2];
            for (int a = 0; a < other.Count; a++)
            {
                var x = other[a];
                if (a >= own.Count)
                    throw new NotSupportedException("-extend_mode relative_triplet: residue_pair_relative_extended_list reads the weights of the library list of the first residue at the positions of the longer list of the second residue, past the end of its array (undefined behaviour of the program: heap contents); the result cannot be reproduced.");
                total += own[a].Weight;
                if (_hasch[x.Seq][x.Res] != 0) score += 2 * Math.Min(_hasch[x.Seq][x.Res], x.Weight);
            }

            if (total == 0) throw new TCoffeeAbortedException("floating point exception in residue_pair_relative_extended_list: integer division by the total weight 0", 8);
            score = (Settings.Normalise * score) / total;
            foreach (var e in own) _hasch[e.Seq][e.Res] = 0;
            _hasch[s1][r1] = _hasch[s2][r2] = 0;
            return score * TcScoreK;
        }

        public Func<int, int, int> ColumnScore(TcGroup g1, TcGroup g2)
        {
            var l1 = g1.Order.OrderBy(x => x).ToArray();
            var l2 = g2.Order.OrderBy(x => x).ToArray();
            switch (Settings.Extend)
            {
                case TCoffeeExtendMode.SlowTriplet:
                    return (c1, c2) => SlowCost(g1, l1, c1, g2, l2, c2);
                case TCoffeeExtendMode.Triplet:
                    return (c1, c2) => SlowCost(g1, l1, c1, g2, l2, c2) - TcScoreK * NoMatch;
                case TCoffeeExtendMode.FastTriplet:
                    return (c1, c2) => FastCost(g1, l1, c1, g2, l2, c2);
                default:
                    {
                        // fast_get_dp_cost_2: the pairs of sequences sorted by the identity of the distance matrix (qsort: stable), visited from the most similar one
                        var pairs = new List<(int S1, int S2, int Sim)>();
                        foreach (int a in l1)
                            foreach (int b in l2)
                                pairs.Add((a, b, Sim == null ? 0 : Sim[a, b]));
                        var sorted = pairs.OrderBy(p => p.Sim).ToArray();
                        return (c1, c2) => FastCost2(g1, c1, g2, c2, sorted);
                    }
            }
        }

        private int SlowCost(TcGroup g1, int[] l1, int col1, TcGroup g2, int[] l2, int col2)
        {
            double score = 0;
            int resRes = 0;
            foreach (int a in l1)
                foreach (int b in l2)
                {
                    int rs1 = a, rs2 = b, r1 = g1.Rows[a][col1], r2 = g2.Rows[b][col2];
                    if (rs1 > rs2) { (rs1, rs2) = (rs2, rs1); (r1, r2) = (r2, r1); }
                    if (r1 > 0 && r2 > 0)
                    {
                        resRes++;
                        score += ResiduePairExtended(rs1, r1, rs2, r2);
                    }
                }

            score = resRes == 0 ? 0 : score / resRes;
            return (int)(score - TcScoreK * NoMatch);
        }

        private int FastCost2(TcGroup g1, int col1, TcGroup g2, int col2, (int S1, int S2, int Sim)[] sorted)
        {
            bool relative = Settings.Extend == TCoffeeExtendMode.RelativeTriplet;
            double score = 0;
            int resRes = 0;
            for (int a = sorted.Length - 1; a >= 0; a--)
            {
                int rs1 = sorted[a].S1, rs2 = sorted[a].S2;
                int r1 = g1.Rows[rs1][col1], r2 = g2.Rows[rs2][col2];
                if (r1 > 0 && r2 > 0) resRes++;
                if (rs1 > rs2) { (rs1, rs2) = (rs2, rs1); (r1, r2) = (r2, r1); }
                score += relative ? ResiduePairRelative(rs1, r1, rs2, r2) : ResiduePairExtended(rs1, r1, rs2, r2);
                if (resRes >= Settings.MaxNPair && Settings.MaxNPair != 0) break;
            }

            score = resRes == 0 ? 0 : score / resRes;
            return (int)(score - TcScoreK * NoMatch);
        }

        /// <summary><c>fast_get_dp_cost</c> (static hashes of the library partners of the columns, literal port including the cross-clearing of the shared targets).</summary>
        private int FastCost(TcGroup g1, int[] l1, int col1, TcGroup g2, int[] l2, int col2)
        {
            double score = 0;
            int nGap1 = 0, nGap2 = 0;
            foreach (int s in l1)
            {
                int r1 = g1.Rows[s][col1];
                if (r1 <= 0) nGap1++;
                else
                {
                    _in1[s][r1] = 1;
                    foreach (var e in Index.Edges[s][r1])
                    {
                        _hasch1[e.Seq][e.Res] += e.Weight;
                        _n1[e.Seq][e.Res]++;
                    }
                }
            }

            foreach (int s in l2)
            {
                int r2 = g2.Rows[s][col2];
                if (r2 <= 0) nGap2++;
                else
                {
                    _in2[s][r2] = 1;
                    foreach (var e in Index.Edges[s][r2])
                    {
                        _hasch2[e.Seq][e.Res] += e.Weight;
                        _n2[e.Seq][e.Res]++;
                    }
                }
            }

            if (l2.Length < l1.Length)
            {
                foreach (int s in l2)
                {
                    int r2 = g2.Rows[s][col2];
                    if (r2 <= 0) continue;
                    foreach (var e in Index.Edges[s][r2])
                    {
                        int ts = e.Seq, tr = e.Res;
                        if (_hasch2[ts][tr] != 0 && _hasch1[ts][tr] != 0 && !(_in1[ts][tr] != 0 || _in2[ts][tr] != 0))
                            score += Math.Min((long)_hasch2[ts][tr] * _n1[ts][tr], (long)_hasch1[ts][tr] * _n2[ts][tr]);
                        else if (_hasch2[ts][tr] != 0 && _in1[ts][tr] != 0)
                            score += (long)_hasch2[ts][tr] * (_n1[ts][tr] + 1);
                        else if (_hasch1[ts][tr] != 0 && _in2[ts][tr] != 0)
                            score += (long)_hasch1[ts][tr] * (_n2[ts][tr] + 1);
                        _hasch2[ts][tr] = 0;
                        _n2[ts][tr] = 0;
                    }

                    _hasch2[s][r2] = 0;
                    _in2[s][r2] = 0;
                }

                foreach (int s in l1)
                {
                    int r1 = g1.Rows[s][col1];
                    if (r1 <= 0) continue;
                    _in1[s][r1] = 0;
                    _hasch1[s][r1] = 0;
                    foreach (var e in Index.Edges[s][r1])
                    {
                        _hasch1[e.Seq][e.Res] = 0;
                        _n1[e.Seq][e.Res] = 0;
                    }
                }
            }
            else
            {
                foreach (int s in l1)
                {
                    int r1 = g1.Rows[s][col1];
                    if (r1 <= 0) continue;
                    foreach (var e in Index.Edges[s][r1])
                    {
                        int ts = e.Seq, tr = e.Res;
                        if (_hasch1[ts][tr] != 0 && _hasch2[ts][tr] != 0 && !(_in2[ts][tr] != 0 || _in1[ts][tr] != 0))
                            score += Math.Min((long)_hasch1[ts][tr] * _n2[ts][tr], (long)_hasch2[ts][tr] * _n1[ts][tr]);
                        else if (_hasch1[ts][tr] != 0 && _in2[ts][tr] != 0)
                            score += (long)_hasch1[ts][tr] * (_n2[ts][tr] + 1);
                        else if (_hasch2[ts][tr] != 0 && _in1[ts][tr] != 0)
                            score += (long)_hasch2[ts][tr] * (_n1[ts][tr] + 1);
                        _hasch1[ts][tr] = 0;
                        _n1[ts][tr] = 0;
                    }

                    _hasch1[s][r1] = 0;
                    _in1[s][r1] = 0;
                }

                foreach (int s in l2)
                {
                    int r2 = g2.Rows[s][col2];
                    if (r2 <= 0) continue;
                    _in2[s][r2] = 0;
                    _hasch1[s][r2] = 0;
                    foreach (var e in Index.Edges[s][r2])
                    {
                        _hasch2[e.Seq][e.Res] = 0;
                        _n2[e.Seq][e.Res] = 0;
                    }
                }
            }

            score = (score * TcScoreK) / ((l1.Length - nGap1) * (l2.Length - nGap2));
            if (Settings.Normalise != 0) score = (score * Settings.Normalise) / MaxExt;
            return (int)(score - TcScoreK * NoMatch);
        }
    }

    // ----- gotoh_pair_wise ------------------------------------------------------------------------------------------------

    /// <summary><c>gotoh_pair_wise</c> (and <c>gotoh_pair_wise_test</c>) on column scores: returns the presence masks of the two groups and the score of the last cell.</summary>
    private static (List<bool> Present1, List<bool> Present2, int Score) TcGotohProfile(int l1, int l2, Func<int, int, int> cost, int gopIn, int gepIn, int tgMode, bool test)
    {
        int gop = gopIn * TcScoreK, gep = gepIn * TcScoreK;
        if (test) { gop = 0; gep = -1; }
        int len = Math.Max(l1, l2) + 1;
        var cc = new int[len]; var dd = new int[len]; var ddg = new int[len];
        var trace = new int[l1 + 1][];
        for (int i = 0; i <= l1; i++) trace[i] = new int[l2 + 1];
        trace[0][0] = 1;
        for (int j = 1; j <= l2; j++) trace[0][j] = -1;

        int termGop = tgMode == 0 ? gop : 0, termGep = tgMode == 2 ? 0 : gep;
        int t = termGop;
        for (int j = 1; j <= l2; j++)
        {
            cc[j] = t = t + termGep;
            dd[j] = t + gop;
        }

        t = termGop;
        int c = 0, score = 0;
        for (int i = 1; i <= l1; i++)
        {
            var tr = trace[i];
            int s = cc[0];
            cc[0] = c = t = t + termGep;
            int e = t + gop;
            tr[0] = 1;
            int eg = 0;
            for (int j = 1; j <= l2; j++)
            {
                int sub = cost(i - 1, j - 1);
                bool edgeRow = i == l1 || i == 1;
                int lGop = edgeRow ? termGop : gop, lGep = edgeRow ? termGep : gep;
                if (e > c + lGop) eg++;
                else eg = 1;
                e = Math.Max(e, c + lGop) + lGep;
                bool edgeCol = j == l2 || j == 1;
                lGop = edgeCol ? termGop : gop;
                lGep = edgeCol ? termGep : gep;
                if (dd[j] > cc[j] + lGop) ddg[j]++;
                else ddg[j] = 1;
                dd[j] = Math.Max(dd[j], cc[j] + lGop) + lGep;

                int fop;
                if (test)
                {
                    if (e >= dd[j] && e >= s + sub) { c = e; fop = -1; }
                    else if (s + sub >= dd[j]) { c = s + sub; fop = 0; }
                    else { c = dd[j]; fop = 1; }
                }
                else
                {
                    // best_int(e, s + sub, dd[j]) takes the first strictly greater value; the tie rules then prefer the substitution over the insertion
                    c = e; fop = 0;
                    if (s + sub > c) { c = s + sub; fop = 1; }
                    if (dd[j] > c) { c = dd[j]; fop = 2; }
                    if (fop == 0 && s + sub == e) fop = 1;
                    else if (fop == 2 && s + sub == dd[j]) fop = 1;
                    fop -= 1;
                }

                s = cc[j];
                cc[j] = c;
                if (fop < 0) tr[j] = fop * eg;
                else if (fop > 0) tr[j] = fop * ddg[j];
                else tr[j] = 0;
            }

            score = c;
        }

        int ii = l1, jj = l2;
        var al0 = new List<bool>(); var al1 = new List<bool>();
        while (ii >= 0 && jj >= 0 && ii + jj != 0)
        {
            int k;
            if (ii == 0) k = -1;
            else if (jj == 0) k = 1;
            else k = trace[ii][jj];

            if (k == 0)
            {
                al0.Add(true); al1.Add(true);
                ii--; jj--;
            }
            else if (k > 0)
            {
                for (int a = 0; a < k; a++) { al0.Add(true); al1.Add(false); ii--; }
            }
            else
            {
                for (int a = 0; a > k; a--) { al0.Add(false); al1.Add(true); jj--; }
            }
        }

        al0.Reverse(); al1.Reverse();
        return (al0, al1, score);
    }

    // ----- myers_miller_pair_wise on profile columns ------------------------------------------------------------------------

    private static (List<bool> Present1, List<bool> Present2, int Score) TcMyersMillerProfile(int l1, int l2, Func<int, int, int> cost, int gop, int gep)
    {
        var st = new MyersMillerState(l1, l2, cost, gop, gep);
        int score = st.Diff(0, l1, 0, l2, 0, 0);
        var al0 = new List<bool>(); var al1 = new List<bool>();
        int i = 0, j = 0, idx = 0;
        while (!(i == l1 && j == l2))
        {
            int op = st.Script[idx++];
            if (op == 0)
            {
                al0.Add(true); al1.Add(true);
                i++; j++;
            }
            else if (op > 0)
            {
                for (int a = 0; a < op; a++) { al0.Add(false); al1.Add(true); j++; }
            }
            else
            {
                for (int a = 0; a < -op; a++) { al0.Add(true); al1.Add(false); i++; }
            }
        }

        return (al0, al1, score);
    }

    // ----- the fasta family on profiles ---------------------------------------------------------------------------------------

    /// <summary><c>evaluate_diagonals_with_clist</c>: every library entry of a residue of a sequence of group 1 to a sequence of group 2 adds the column score of (residue, partner) - the residue NUMBERS serve as column indices - to its diagonal.</summary>
    private static int[][] TcDiagonalsByLibrary(TcDpContext ctx, int[] l1ids, int[] l2ids, int l1, int l2, Func<int, int, int> cost)
    {
        int nDiag = l1 + l2 - 1;
        var diag = new int[nDiag + 2][];
        for (int a = 0; a < nDiag + 2; a++) diag[a] = new[] { a, 0, 0 };
        foreach (int s1 in l1ids)
            foreach (int s2 in l2ids)
                for (int r1 = 1; r1 <= ctx.Index.Len[s1]; r1++)
                    foreach (var e in ctx.Index.Edges[s1][r1])
                        if (e.Seq == s2) diag[e.Res - r1 + l1][1] += cost(r1 - 1, e.Res - 1);
        TcSortDiagonals(diag, nDiag);
        return diag;
    }

    /// <summary>One merge of two profiles with the generic dynamic program of <c>ctx.Settings.Mode</c>: the masks of the two groups and the score of the program's return value.</summary>
    private static (List<bool> Present1, List<bool> Present2, int Score) TcProfileDp(TcDpContext ctx, TcGroup g0, TcGroup g1)
    {
        var st = ctx.Settings;
        int l1 = g0.Columns, l2 = g1.Columns;
        var cost = ctx.ColumnScore(g0, g1);
        var ids0 = g0.Order.OrderBy(x => x).ToArray();
        var ids1 = g1.Order.OrderBy(x => x).ToArray();
        switch (st.Mode)
        {
            case TCoffeeDpMode.Gotoh:
            case TCoffeeDpMode.GotohTest:
                return TcGotohProfile(l1, l2, cost, st.Gop, st.Gep, st.TgMode, st.Mode == TCoffeeDpMode.GotohTest);
            case TCoffeeDpMode.GotohLgp:
                return TcGotohLgpMasks(g0, g1, cost, st.Gop, st.Gep);
            case TCoffeeDpMode.MyersMiller:
                return TcMyersMillerProfile(l1, l2, cost, st.Gop, st.Gep);
            case TCoffeeDpMode.Biphasic:
                return TcBiphasicMasks(l1, l2, cost, st.Gop, st.Gep);
            case TCoffeeDpMode.Glocal:
                return TcGlocalMasks(l1, l2, cost, st.Gop, st.Gep);
            case TCoffeeDpMode.Fasta:
                {
                    var tot = TcDiagonalsByLibrary(ctx, ids0, ids1, l1, l2, cost);
                    int[] list = st.FastaStep == 0 ? TcFlagDiagonals(l1, l2, tot, st.DiagonalThreshold) : TcExtractNDiagonals(l1, l2, tot, st.FastaStep, 0);
                    return TcMakeFastaGotoh(l1, l2, list, cost, st.Gop * TcScoreK, st.Gep * TcScoreK, st.TgMode);
                }

            case TCoffeeDpMode.VeryFast:
                {
                    var tot = TcDiagonalsByLibrary(ctx, ids0, ids1, l1, l2, cost);
                    int[] list = TcExtractNDiagonals(l1, l2, tot, 20, 5);
                    return TcMakeFastaGotoh(l1, l2, list, cost, st.Gop * TcScoreK, st.Gep * TcScoreK, st.TgMode);
                }

            case TCoffeeDpMode.CFasta:
                {
                    var tot = TcDiagonalsByLibrary(ctx, ids0, ids1, l1, l2, cost);
                    int step = st.FastaStep != 0 ? st.FastaStep : Math.Max((int)Math.Log(Math.Max(Math.Min(l1, l2), 1)), 20);
                    int maxChosen = l1 + l2 - 1;
                    int chosen = Math.Min(20 + step, maxChosen);
                    var result = TcMakeFastaGotoh(l1, l2, TcExtractNDiagonals(l1, l2, tot, chosen, 0), cost, st.Gop * TcScoreK, st.Gep * TcScoreK, st.TgMode);
                    int score = result.Score, newScore = 0;
                    while (newScore != score && chosen < maxChosen)
                    {
                        score = newScore;
                        chosen = Math.Min(chosen + step, maxChosen);
                        result = TcMakeFastaGotoh(l1, l2, TcExtractNDiagonals(l1, l2, tot, chosen, 0), cost, st.Gop * TcScoreK, st.Gep * TcScoreK, st.TgMode);
                        newScore = result.Score;
                    }

                    return (result.Present1, result.Present2, newScore);
                }

            default:
                throw new InvalidOperationException($"The dynamic program {st.Mode} has no generic implementation.");
        }
    }

    /// <summary><c>rec_local_tree_aln</c> with a generic dynamic program: left subtree, right subtree, then <c>profile_pair_wise (L, R)</c> (the left group is profile 1, no size-based swap).</summary>
    private static TcGroup TcAlignNodeDp(TCoffeeTreeNode node, IReadOnlyList<string> seqs, TcDpContext ctx)
    {
        if (node.Left is null)
        {
            var g = new TcGroup();
            g.Order.Add(node.LeafIndex);
            g.Rows[node.LeafIndex] = Enumerable.Range(1, seqs[node.LeafIndex].Length).ToArray();
            return g;
        }

        var left = TcAlignNodeDp(node.Left, seqs, ctx);
        var right = TcAlignNodeDp(node.Right!, seqs, ctx);
        if (ctx.Settings.Mode is TCoffeeDpMode.Ktup or TCoffeeDpMode.Precomputed or TCoffeeDpMode.GotohSw or TCoffeeDpMode.Test) return TcMergeWithoutGlobalAlignment(ctx, left, right);
        var (m0, m1, _) = TcProfileDp(ctx, left, right);
        var merged = new TcGroup();
        merged.Order.AddRange(left.Order);
        merged.Order.AddRange(right.Order);
        TcApply(left, m0, merged);
        TcApply(right, m1, merged);
        return merged;
    }

    /// <summary>
    /// The dynamic programs that do not return a global alignment of the two groups: <c>ktup_pair_wise</c> only computes a k-tuple score (the rows stay as they are, even of unequal length),
    /// <c>precomputed_pair_wise</c> computes the identity of two aligned groups and stops the program ("improper use of the function precomputed pairwise", exit status 1) when the first rows of the two groups differ in length,
    /// <c>gotoh_pair_wise_sw</c> keeps only the columns of the best local alignment (the rows of BOTH groups lose everything outside it; no positive cell empties every row, and a group with empty rows leaves the pair untouched: <c>empty_pair_wise</c> of a LOCAL mode).
    /// </summary>
    /// <summary>Diagnostics for the tests: when set, every merge of <see cref="TcMergeWithoutGlobalAlignment"/> appends (mode, columns of the first group, columns of the second group, row length written, longest input sequence).</summary>
    [ThreadStatic]
    internal static List<(TCoffeeDpMode Mode, int L0, int L1, int Written, int MaxLen)>? TcMergeTrace;

    private static TcGroup TcMergeWithoutGlobalAlignment(TcDpContext ctx, TcGroup left, TcGroup right)
    {
        var merged = new TcGroup();
        merged.Order.AddRange(left.Order);
        merged.Order.AddRange(right.Order);
        var mode = ctx.Settings.Mode;
        if (mode == TCoffeeDpMode.Precomputed && left.Columns != right.Columns)
            throw new TCoffeeAbortedException("ERROR: improper use of the function precomputed pairwise:[FATAL:T-COFFEE] (precomputed_pair_wise needs two groups whose first rows have the same length)");
        if (mode == TCoffeeDpMode.Test)
        {
            int l0 = left.Columns, l1 = right.Columns;
            int n = l0 < 5 ? l0 / 2 : 5;
            if (l1 < n) throw new TCoffeeAbortedException("segmentation fault in test_pair_wise: generate_null (l1 - n) is called with a negative size, vcalloc (size + 1 <= 0) returns NULL and the terminator is written through it (exit status 11)", 11);
            int pad0 = l1 - n, pad1 = l0 - n;
            TcMergeTrace?.Add((mode, l0, l1, l0 + pad0, ctx.Index.Len.Max()));
            foreach (int sq in left.Order) merged.Rows[sq] = left.Rows[sq].Concat(new int[pad0]).ToArray();
            foreach (int sq in right.Order) merged.Rows[sq] = right.Rows[sq].Concat(new int[pad1]).ToArray();
            return merged;
        }

        if (mode == TCoffeeDpMode.GotohSw && left.Columns > 0 && right.Columns > 0)
        {
            var sw = TcGotohSwMasks(left.Columns, right.Columns, ctx.ColumnScore(left, right), ctx.Settings.Gop, ctx.Settings.Gep, ctx.NoMatch);
            TcMergeTrace?.Add((mode, left.Columns, right.Columns, sw?.Al0.Count ?? -1, ctx.Index.Len.Max()));
            foreach (var (g, which) in new[] { (left, 0), (right, 1) })
                foreach (int sq in g.Order)
                {
                    if (sw is null) { merged.Rows[sq] = Array.Empty<int>(); continue; }
                    var mask = which == 0 ? sw.Value.Al0 : sw.Value.Al1;
                    int start = which == 0 ? sw.Value.Start0 : sw.Value.Start1;
                    var src = g.Rows[sq];
                    var dst = new int[mask.Count];
                    int p = start;
                    for (int c = 0; c < mask.Count; c++) dst[c] = mask[c] ? src[p++] : 0;
                    merged.Rows[sq] = dst;
                }

            return merged;
        }

        foreach (var g in new[] { left, right })
            foreach (int sq in g.Order) merged.Rows[sq] = g.Rows[sq];
        return merged;
    }

    /// <summary><c>pair_wise</c> of the two single sequences a and b (<c>-distance_matrix_mode slow | clustalw</c> with a generic dp mode): rows and the score of the dynamic program.</summary>
    private static (string RowA, string RowB, int Score) TcDpPairRows(IReadOnlyList<string> seqs, int a, int b, TcDpContext ctx)
    {
        var ga = new TcGroup(); ga.Order.Add(a); ga.Rows[a] = Enumerable.Range(1, seqs[a].Length).ToArray();
        var gb = new TcGroup(); gb.Order.Add(b); gb.Rows[b] = Enumerable.Range(1, seqs[b].Length).ToArray();
        var (m0, m1, score) = TcProfileDp(ctx, ga, gb);
        var (ra, rb) = TcMasksToRows(seqs[a], seqs[b], m0, m1);
        return (ra, rb, score);
    }

    // ----- what the program does with the modes that are not ported ------------------------------------------------------------

    private static readonly string[] TcExtendModesWithoutPort = { "pc", "mixt", "quadruplet", "g_coffee", "g_coffee_quadruplets", "rna", "rna1", "rna2", "rna3", "rna4", "test" };

    /// <summary>
    /// The refusals of <c>choose_extension_mode</c> / <c>initialize_scoring_scheme</c> (an unknown <c>-extend_mode</c> or <c>matrix</c> stops the program for any number of sequences, exit status 1; <c>ssp</c> / <c>tm</c> die in the first evaluation (exit status 11) and
    /// <c>curvature</c> stops (exit status 1) as soon as there are two sequences, whatever the dp mode is) and of <c>get_pair_wise_function</c> / the dynamic programs named by <c>-dp_mode</c> (measured on five sets of 2-6 protein / DNA sequences each; the names whose outcome depends on the input are not
    /// listed and raise <see cref="NotSupportedException"/>).
    /// </summary>
    private static void TcRefuseModes(TCoffeeProgramOptions o, int k)
    {
        if (o.ExtendMode == TCoffeeExtendMode.Other)
        {
            string m = o.ExtendModeText ?? string.Empty;
            if (m == "rna0") { /* residue_pair_extended_list + slow_get_dp_cost: handled as slow_triplet by the caller */ }
            else if (m == "matrix") throw new TCoffeeAbortedException("ERROR: use -in <mat_name> to set the matrix to be used [FATAL] (choose_extension_mode: -extend_mode matrix needs CL->M)");
            else if (m is "ssp" or "tm")
            {
                if (k >= 2) throw new TCoffeeAbortedException($"segmentation fault in initialize_scoring_scheme: -extend_mode {m} reads secondary-structure data that does not exist", 11);
            }
            else if (m == "curvature")
            {
                if (k >= 2) throw new TCoffeeAbortedException("ERROR: FORCED EXIT (NON INTERACTIVE MODE) in evaluate_curvature_score (initialize_scoring_scheme)");
            }
            else if (!TcExtendModesWithoutPort.Contains(m))
                throw new TCoffeeAbortedException($"ERROR: {m} is an unknown extend_mode[FATAL:T-COFFEE] (choose_extension_mode)");
        }

        if (k < 2 || o.DpMode != TCoffeeDpMode.Other) return;
        string dp = o.DpModeText ?? string.Empty;
        switch (dp)
        {
            case "idscore_pair_wise" or "subop1_pair_wise" or "subop2_pair_wise" or "proba_pair_wise" or "pavie_pair_wise" or "domain_pair_wise" or "fasta_cdna_pair_wise" or "cfasta_cdna_pair_wise":
                throw new TCoffeeAbortedException($"segmentation fault in {dp}: the function needs data the progressive stage does not provide (exit status 11 on every measured set)", 11);
            case "gotoh_pair_wise_lalign":
                throw new TCoffeeAbortedException("gotoh_pair_wise_lalign: the program aborts (SIGABRT, exit status 6) on every measured set", 6);
            case "seq_viterbi_pair_wise" or "hh_pair_wise":
                throw new TCoffeeAbortedException($"{dp}: the program stops with a FAILURE status (exit status 1) on every measured set");
            case "cfasta_sw_pair_wise":
                throw new TCoffeeAbortedException("cfasta_gotoh_pair_wise_sw not implemented yet (exit status 1)");
            default:
                if (!TcDpModesWithoutPort.Contains(dp)) throw new TCoffeeAbortedException($"[{dp}] is an unknown mode for dp_mode[FATAL] (get_pair_wise_function)");
                throw new NotSupportedException($"-dp_mode {dp} runs in the program but is not ported yet.");
        }
    }

    private static readonly string[] TcDpModesWithoutPort =
    {
        "viterbi_pair_wise", "viterbiL_pair_wise", "viterbiD_pair_wise", "procoffee_pair_wise", "co_pair_wise",
        "linked_pair_wise_collapse", "sim_pair_wise_lalign",
    };
}
