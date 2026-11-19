namespace Seqeron.Genomics.Alignment;

// =====================================================================================================================
// The distance matrix, guide tree modes and sequence weights of the T-Coffee 13.41 program (audit A1-13d, B09).
//
//   seq2distance_matrix (util_constraints_list.c) fills for every pair (a < b) an identity `id` (similarity_matrix, 0..100) and a `score` (score_similarity_matrix):
//     ktup    id = ktup_dist_mat similarity (default)                                   score = id x SCORE_K (10)
//     fast    id = get_seq_sim of the Myers-Miller alignment (idmat / blosum62mt)      score = id x 10
//     idscore id = idscore_pairseq (affine DP, own gap rule)                           score = id            (NOT scaled: a quirk, all distances are then near 100)
//     no      id = 100                                                                  score = 1000
//   make_tree: distances = normalize_array(sim_array2dist_array(score, 1000), 1000, 100) (integer division) and then
//     nj    : slow_nj_tree / fast_nj_tree on distance / 100, reroot of the trichotomy (as the default tree)
//     upgma : int_dist2upgma_tree: merge of the lowest integer distance (first on ties, rows a < b), the new row is the integer mean (x + y) / 2.
//   weight_constraint_list(seq_weight = t_coffee): compute_t_coffee_weight on `id` (float arithmetic), re_weight_constraint_list multiplies every directed entry by
//   MIN(w[s1], w[s2]) in float and truncates.
// =====================================================================================================================

public static partial class SequenceAligner
{
    private const int TcMaxIdentity = 100;

    private sealed class TcDistances
    {
        public int[,] Sim = new int[0, 0];
        public int[,] Score = new int[0, 0];
    }

    /// <summary><c>is_in_same_group_aa(r1, r2, 0, NULL, mode)</c> of the distance measure: with groups made from the matrix <c>idmat</c> (a library method ran first) or one group per letter.</summary>
    private static bool TcSameGroupMode(char a, char b, bool matrixGroups)
    {
        if (matrixGroups) return TcSameGroup(a, b);
        return char.IsAsciiLetter(a) && char.IsAsciiLetter(b) && char.ToUpperInvariant(a) == char.ToUpperInvariant(b);
    }

    /// <summary><c>get_seq_sim(row1, row2, "-", "idmat_sim1" | ...)</c> for the measures of <c>-distance_matrix_sim_mode</c>.</summary>
    private static int TcSeqSimMeasure(string r1, string r2, string simMode, bool matrixGroups)
    {
        string mode = simMode.TrimStart('_');
        string measure = "";
        bool hasMeasure = false;
        int us = mode.IndexOf('_');
        if (us >= 0)
        {
            measure = mode[(us + 1)..];
            mode = mode[..us];
            hasMeasure = true;
        }

        if (mode.Contains("idscore", StringComparison.Ordinal))
            throw new NotSupportedException("A distance_matrix_sim_mode containing 'idscore' is not supported.");
        if (mode != "idmat" && !mode.Contains("sim", StringComparison.Ordinal))
            throw new NotSupportedException($"distance_matrix_sim_mode '{simMode}' is not supported (the residue groups of other matrices are not ported); use idmat_sim1/sim2/sim3/gap1/logid.");

        int pos0 = 0, pos1 = 0, pos2 = 0, sim = 0, gap = 0;
        for (int a = 0; a < r1.Length; a++)
        {
            bool p1 = r1[a] != '-', p2 = r2[a] != '-';
            if (p1) pos1++;
            if (p2) pos2++;
            if (p1 && p2)
            {
                pos0++;
                if (TcSameGroupMode(r1[a], r2[a], matrixGroups)) sim++;
            }
            else if (p1 != p2) gap++;
        }

        if (mode.Contains("cov", StringComparison.Ordinal)) return pos0 + gap == 0 ? 0 : pos0 * TcMaxIdentity / (pos0 + gap);
        if (!hasMeasure || measure == "sim1" || measure == "sim") return pos0 == 0 ? 0 : sim * TcMaxIdentity / pos0;
        if (measure == "sim2") return pos1 == 0 || pos2 == 0 ? 0 : sim * TcMaxIdentity / Math.Min(pos1, pos2);
        if (measure == "sim3") return pos1 == 0 || pos2 == 0 ? 0 : sim * TcMaxIdentity / Math.Max(pos1, pos2);
        if (measure == "gap1") return TcMaxIdentity - (r1.Length == 0 ? TcMaxIdentity : gap * TcMaxIdentity / r1.Length);
        if (measure == "logid") return TcLogId(pos0, sim);
        if (!mode.Contains("sim", StringComparison.Ordinal) || pos0 == 0) return 0;
        return sim * TcMaxIdentity / pos0;
    }

    /// <summary><c>logid_score(len = pos0, sim)</c> (note the program's argument order: its parameters are named <c>(sim, len)</c> and called with (pos0, sim)).</summary>
    private static int TcLogId(int sim, int len)
    {
        if (len == 0) return (int)(0.33 * TcMaxIdentity);
        float score = (float)sim / len;
        score = score > 0.9f ? 1.0f : (float)-Math.Log10(1.0 - score);
        score *= TcMaxIdentity;
        return (int)score;
    }

    /// <summary><c>idscore_pairseq</c>: the affine DP of util_dp_gotoh_nw.c that returns the identity directly (letters lower-cased, terminal gaps free of the opening penalty).</summary>
    private static int TcIdscorePair(string seq1, string seq2, int gop, int gep, TcMatrix matrix, string compMode)
    {
        string s1 = seq1.ToLowerInvariant(), s2 = seq2.ToLowerInvariant();
        int l1 = s1.Length, l2 = s2.Length;
        // rows: [0]/[1] score previous / current, [2]/[3] identities, [4]/[5] aligned columns
        var I = new int[6][]; var D = new int[6][]; var M = new int[6][];
        for (int r = 0; r < 6; r++) { I[r] = new int[l2 + 1]; D[r] = new int[l2 + 1]; M[r] = new int[l2 + 1]; }
        for (int j = 0; j <= l2; j++) { D[0][j] = gep * j; M[0][j] = 2 * gep * j; D[4][j] = 0; }
        for (int i = 1; i <= l1; i++)
        {
            I[1][0] = i * gep;
            M[1][0] = 2 * i * gep;
            for (int j = 1; j <= l2; j++)
            {
                int score = matrix.Score(char.ToUpperInvariant(s1[i - 1]), char.ToUpperInvariant(s2[j - 1]));
                int id = s1[i - 1] == s2[j - 1] ? 1 : 0;
                int igop = i == l1 || j == l2 ? 0 : gop;

                if (D[0][j] + gep > M[0][j] + igop + gep) { D[1][j] = D[0][j] + gep; D[3][j] = D[2][j]; D[5][j] = D[4][j]; }
                else { D[1][j] = M[0][j] + igop + gep; D[3][j] = M[2][j]; D[5][j] = M[4][j]; }

                if (I[1][j - 1] + gep > M[1][j - 1] + igop + gep) { I[1][j] = I[1][j - 1] + gep; I[3][j] = I[3][j - 1]; I[5][j] = I[5][j - 1]; }
                else { I[1][j] = M[1][j - 1] + igop + gep; I[3][j] = M[3][j - 1]; I[5][j] = M[5][j - 1]; }

                int match = M[0][j - 1] + score;
                if (I[1][j] > match && I[1][j] > D[1][j]) { M[1][j] = I[1][j]; M[3][j] = I[3][j]; M[5][j] = I[5][j]; }
                else if (D[1][j] > match) { M[1][j] = D[1][j]; M[3][j] = D[3][j]; M[5][j] = D[5][j]; }
                else { M[1][j] = match; M[3][j] = M[2][j - 1] + id; M[5][j] = M[4][j - 1] + 1; }
            }

            for (int r = 0; r < 6; r += 2)
            {
                (I[r], I[r + 1]) = (I[r + 1], I[r]);
                (D[r], D[r + 1]) = (D[r + 1], D[r]);
                (M[r], M[r + 1]) = (M[r + 1], M[r]);
            }
        }

        int l;
        if (compMode.Contains("sim2", StringComparison.Ordinal)) { l = Math.Min(l1, l2); return l == 0 ? 0 : M[2][l2] * 100 / l; }
        if (compMode.Contains("sim3", StringComparison.Ordinal)) { l = Math.Max(l1, l2); return l == 0 ? 0 : M[2][l2] * 100 / l; }
        if (compMode.Contains("cov", StringComparison.Ordinal)) { l = Math.Max(l1, l2); return l == 0 ? 0 : M[4][l2] * 100 / l; }
        l = M[4][l2];
        return l == 0 ? 0 : M[2][l2] * 100 / l;
    }

    /// <summary>The distance matrix of the program (<c>seq2distance_matrix</c>) for the chosen <c>-distance_matrix_mode</c>.</summary>
    private static TcDistances TcDistanceMatrix(IReadOnlyList<string> seqs, TcSeqKind kind, TCoffeeProgramOptions o, bool matrixGroups)
    {
        int n = seqs.Count;
        var dm = new TcDistances { Sim = new int[n, n], Score = new int[n, n] };
        bool protein = kind == TcSeqKind.Protein;
        if (o.DistanceMatrixMode == TCoffeeDistanceMode.Ktup)
        {
            var ktup = TCoffeeKtupSimilarity(seqs, protein);
            for (int a = 0; a < n; a++)
                for (int b = 0; b < n; b++)
                {
                    dm.Sim[a, b] = ktup[a, b];
                    dm.Score[a, b] = a == b ? 0 : ktup[a, b] * TcScoreK;
                }

            return dm;
        }

        TcMatrix mat = protein ? TcMatrix.Blosum62 : TcMatrix.IdMat;
        int gop = protein ? (int)(mat.AverageMismatch() * 10) : -10;
        for (int a = 0; a < n; a++)
        {
            dm.Sim[a, a] = TcMaxIdentity;
            for (int b = a + 1; b < n; b++)
            {
                int id, score;
                switch (o.DistanceMatrixMode)
                {
                    case TCoffeeDistanceMode.Fast:
                        {
                            var (ra, rb) = TCoffeeSlowPairAlign(seqs[a], seqs[b], protein);
                            id = TcSeqSimMeasure(ra, rb, o.DistanceMatrixSimMode, matrixGroups);
                            score = id * TcScoreK;
                            break;
                        }

                    case TCoffeeDistanceMode.Idscore:
                        id = score = TcIdscorePair(seqs[a], seqs[b], gop, -1, mat, o.DistanceMatrixSimMode);
                        break;
                    case TCoffeeDistanceMode.No:
                        id = TcMaxIdentity;
                        score = id * TcScoreK;
                        break;
                    default:
                        throw new NotSupportedException($"distance_matrix_mode {o.DistanceMatrixMode} is not supported.");
                }

                dm.Sim[a, b] = dm.Sim[b, a] = id;
                dm.Score[a, b] = dm.Score[b, a] = score;
            }
        }

        return dm;
    }

    /// <summary>make_tree: <c>normalize_array(sim_array2dist_array(score, 1000), 1000, 100)</c> followed by NJ (+ reroot) or UPGMA.</summary>
    internal static TCoffeeTreeNode TcGuideTreeFor(IReadOnlyList<string> seqs, TcSeqKind kind, TCoffeeProgramOptions o, bool matrixGroups)
    {
        int k = seqs.Count;
        if (k == 1) return TCoffeeTreeNode.Leaf(0);
        if (k == 2) return TCoffeeTreeNode.Join(TCoffeeTreeNode.Leaf(0), TCoffeeTreeNode.Leaf(1));
        var dm = TcDistanceMatrix(seqs, kind, o, matrixGroups);
        var dist = new int[k, k];
        for (int a = 0; a < k; a++)
            for (int b = 0; b < k; b++)
                dist[a, b] = (TcMaxIdentity * TcScoreK - dm.Score[a, b]) * 100 / (TcMaxIdentity * TcScoreK);

        if (o.TreeMode == TCoffeeTreeMode.Upgma) return TcUpgmaTree(dist, k);
        var d = new double[k, k];
        for (int a = 0; a < k; a++)
            for (int b = 0; b < k; b++)
                d[a, b] = a == b ? 0 : (double)dist[a, b] / TcMaxIdentity;
        var unrooted = k < 100 ? TcNeighborJoining(d, k) : TcNeighborJoiningFast(d, k);
        return ToTCoffeeTree(ToGuideNode(ClustalWRootNodes(unrooted, k).Root));
    }

    /// <summary><c>int_dist2upgma_tree</c>: the pair with the smallest integer distance (first in row-major order on ties) is merged, the new row is the integer mean.</summary>
    private static TCoffeeTreeNode TcUpgmaTree(int[,] mat, int n)
    {
        var nodes = new TCoffeeTreeNode[n];
        for (int a = 0; a < n; a++) nodes[a] = TCoffeeTreeNode.Leaf(a);
        var used = new bool[n];
        TCoffeeTreeNode root = nodes[0];
        for (int left = n; left > 1; left--)
        {
            int bestA = 0, bestB = 0, best = 0;
            bool set = false;
            for (int a = 0; a < n - 1; a++)
            {
                if (used[a]) continue;
                for (int b = a + 1; b < n; b++)
                {
                    if (used[b]) continue;
                    int w = mat[a, b];
                    if (!set || w < best) { best = w; bestA = a; bestB = b; set = true; }
                }
            }

            for (int a = 0; a < n; a++)
                if (!used[a]) mat[bestA, a] = mat[a, bestA] = (mat[bestB, a] + mat[bestA, a]) / 2;
            used[bestB] = true;
            // read_tree turns the written (left, right) into (right, left): the walk of the alignment starts with the merged-in (higher index) subtree
            root = TCoffeeTreeNode.Join(nodes[bestB], nodes[bestA]);
            nodes[bestA] = root;
        }

        return root;
    }

    /// <summary><c>compute_t_coffee_weight</c> on the identities of the distance matrix (single-precision arithmetic as in the original).</summary>
    private static float[] TcTCoffeeSequenceWeights(int[,] sim, int n)
    {
        var w = new float[n];
        Array.Fill(w, 1f);
        for (int a = 0; a < n - 1; a++)
            for (int b = a + 1; b < n; b++)
            {
                float d = (float)sim[a, b] / TcMaxIdentity * 10;
                float p = (float)Math.Pow(d, 3);
                w[a] += p;
                w[b] += p;
            }

        float total = 0;
        for (int b = 0; b < n; b++)
        {
            w[b] = 2 / w[b];
            total += w[b];
        }

        for (int b = 0; b < n; b++) w[b] = w[b] * ((float)n / total);
        return w;
    }
}
