using Seqeron.Genomics.Phylogenetics;
namespace Seqeron.Genomics.Alignment;

// =====================================================================================================================
// The default GUIDE TREE of the T-Coffee 13.41 program (audit A1-13b, F66): `-tree_mode nj -distance_matrix_mode ktup`.
//
//   distances  : ktup_dist_mat (util_dp_fasta_nw.c, adapted from MAFFT 5): residues are coded into 4 (DNA: A, G, T/U, C) or 6 (protein:
//                agjopst / ilmv / bdenqz / hkr / fwy / c) classes, k = min(shortest length, 6) tuples are counted and the multiset intersection
//                of two sequences is compared with the tuple count of the FIRST one: pscore = (int)((s0 - common) / s0 * 3 * 10.0 + 0.5),
//                similarity = 100 - min(pscore(i,j), pscore(j,i)); the tree distance is (100 - similarity) / 100.
//   tree       : `slow_nj_tree` (util_make_tree.c, ClustalW 1.x neighbour joining in double arithmetic with the `MY_EPS` quirk: a candidate whose
//                value is lower by at most 1000 DBL_EPSILON replaces the pair without replacing the minimum), printed with %7.5f, read back
//                as floats; the trichotomy at the root is turned into a binary tree by `reroot` (tree_util.c) exactly as ClustalW does
//                (shared port `ClustalWRootNodes`); two sequences get the dummy tree (a, b).
//   order      : post-order, left subtree first (`rec_local_tree_aln`), the profile with fewer sequences x columns is aligned as profile 1.
// =====================================================================================================================

public static partial class SequenceAligner
{
    private const double TcMyEps = 1000 * 2.2204460492503131e-16; // 1000 * DBL_EPSILON

    /// <summary>
    /// <c>ktup_dist_mat</c>: the pair-wise similarity (0..100) the program derives its guide tree from (k-tuple identity of the MAFFT-style reduced alphabet).
    /// </summary>
    /// <param name="sequences">The sequences.</param>
    /// <param name="protein"><c>true</c> = the six-group protein alphabet, <c>false</c> = DNA/RNA (A, G, T/U, C).</param>
    /// <returns>Symmetric matrix of integer similarities (the diagonal is 100).</returns>
    public static int[,] TCoffeeKtupSimilarity(IReadOnlyList<string> sequences, bool protein = false)
    {
        ArgumentNullException.ThrowIfNull(sequences);
        int n = sequences.Count;
        var pscore = new int[n, n];
        if (n == 0) return pscore;

        string[] groups = protein
            ? new[] { "agjopstAGJOPST", "ilmvILMV", "bdenqzBDENQZ", "hkrHKR", "fwyFWY", "cC" }
            : new[] { "Aa", "Gg", "TtUu", "Cc" };
        var code = new int[256];
        for (int g = 0; g < groups.Length; g++)
            foreach (char c in groups[g]) code[c] = g;
        int ng = groups.Length;

        int minl = sequences.Min(s => s.Length);
        int ktup = Math.Min(minl, 6);
        if (ktup < 1) throw new ArgumentException("Empty sequences are not supported.", nameof(sequences));

        int tsize = (int)Math.Pow(ng, ktup);
        var tuples = new int[n][];
        for (int i = 0; i < n; i++)
        {
            string s = sequences[i];
            var pt = new List<int>();
            for (int p = 0; p + ktup <= s.Length; p++)
            {
                int key = 0;
                for (int a = 0; a < ktup; a++) key = key * ng + code[s[p + a] & 0xFF];
                pt.Add(key);
            }

            tuples[i] = pt.ToArray();
        }

        var mtx = new double[n, n];
        var table = new int[tsize];
        var memo = new int[tsize];
        for (int i = 0; i < n; i++)
        {
            Array.Clear(table);
            foreach (int p in tuples[i]) table[p]++;
            for (int j = i; j < n; j++)
            {
                int value = 0;
                foreach (int p in tuples[j])
                {
                    int tmp = memo[p]++;
                    if (tmp < table[p]) value++;
                }

                foreach (int p in tuples[j]) memo[p] = 0;
                mtx[i, j] = value;
            }
        }

        for (int i = 0; i < n; i++)
        {
            double score0 = mtx[i, i];
            for (int j = 0; j < n; j++)
                pscore[i, j] = (int)((score0 - mtx[Math.Min(i, j), Math.Max(i, j)]) / score0 * 3 * 10.0 + 0.5);
        }

        for (int i = 0; i < n - 1; i++)
            for (int j = i + 1; j < n; j++)
                pscore[i, j] = pscore[j, i] = 100 - Math.Min(pscore[i, j], pscore[j, i]);
        for (int i = 0; i < n; i++) pscore[i, i] = 100;
        return pscore;
    }

    /// <summary>
    /// <c>slow_nj_tree</c> (util_make_tree.c): the unrooted neighbour-joining tree of the program for <c>k</c> &gt;= 3 sequences; the three subtrees of the last
    /// cycle are the children of the returned node, the branch length of every node is stored in its <see cref="PhylogeneticAnalyzer.PhyloNode.BranchLength"/>.
    /// </summary>
    internal static PhylogeneticAnalyzer.PhyloNode TcNeighborJoining(double[,] distance, int k)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(k, 3);
        var tmat = new double[k + 1, k + 1];
        for (int a = 1; a <= k; a++)
            for (int b = 1; b <= k; b++)
                tmat[a, b] = a == b ? 0.0 : distance[a - 1, b - 1];

        var av = new double[k + 1];
        var tkill = new int[k + 1];
        var node = new PhylogeneticAnalyzer.PhyloNode[k + 1];
        for (int i = 1; i <= k; i++)
            node[i] = new PhylogeneticAnalyzer.PhyloNode((i - 1).ToString(System.Globalization.CultureInfo.InvariantCulture))
            {
                Taxa = new List<string> { (i - 1).ToString(System.Globalization.CultureInfo.InvariantCulture) },
            };

        double fnseqs = k, fnseqs2 = 0;
        int mini = 0, minj = 0;
        for (int nc = 1; nc <= k - 3; ++nc)
        {
            double sumd = 0.0;
            for (int j = 2; j <= k; ++j)
                for (int i = 1; i < j; ++i)
                {
                    tmat[j, i] = tmat[i, j];
                    sumd += tmat[i, j];
                }

            double tmin = 99999.0;
            for (int jj = 2; jj <= k; ++jj)
            {
                if (tkill[jj] == 1) continue;
                for (int ii = 1; ii < jj; ++ii)
                {
                    if (tkill[ii] == 1) continue;
                    double diq = 0.0, djq = 0.0;
                    for (int i = 1; i <= k; ++i)
                    {
                        diq += tmat[i, ii];
                        djq += tmat[i, jj];
                    }

                    double dij = tmat[ii, jj];
                    double d2r = diq + djq - (2.0 * dij);
                    double dr = sumd - dij - d2r;
                    fnseqs2 = fnseqs - 2.0;
                    double total = d2r + fnseqs2 * dij + dr * 2.0;
                    total /= 2.0 * fnseqs2;
                    if (total < tmin)
                    {
                        if (tmin - total > TcMyEps) tmin = total;
                        mini = ii;
                        minj = jj;
                    }
                }
            }

            double dio = 0.0, djo = 0.0;
            for (int i = 1; i <= k; ++i)
            {
                dio += tmat[i, mini];
                djo += tmat[i, minj];
            }

            double dmin = tmat[mini, minj];
            dio = (dio - dmin) / fnseqs2;
            djo = (djo - dmin) / fnseqs2;
            double bi = (dmin + dio - djo) * 0.5;
            double bj = dmin - bi;
            bi -= av[mini];
            bj -= av[minj];
            if (Math.Abs(bi) < 0.0001) bi = 0.0;
            if (Math.Abs(bj) < 0.0001) bj = 0.0;

            node[mini].BranchLength = bi;
            node[minj].BranchLength = bj;
            var joined = new PhylogeneticAnalyzer.PhyloNode
            {
                Children = new List<PhylogeneticAnalyzer.PhyloNode> { node[mini], node[minj] },
                Taxa = node[mini].Taxa.Concat(node[minj].Taxa).ToList(),
            };

            if (dmin <= 0.0) dmin = 0.000001;
            av[mini] = dmin * 0.5;
            fnseqs -= 1.0;
            tkill[minj] = 1;
            for (int j = 1; j <= k; ++j)
            {
                if (tkill[j] == 1) continue;
                double da = (tmat[mini, j] + tmat[minj, j]) * 0.5;
                if (mini - j < 0) tmat[mini, j] = da;
                if (mini - j > 0) tmat[j, mini] = da;
            }

            for (int j = 1; j <= k; ++j) tmat[minj, j] = tmat[j, minj] = 0.0;
            node[mini] = joined;
        }

        var l = new List<int>();
        for (int i = 1; i <= k; ++i)
            if (tkill[i] != 1) l.Add(i);
        int l1 = l[0], l2 = l[1], l3 = l[2];
        double b1 = (tmat[l1, l2] + tmat[l1, l3] - tmat[l2, l3]) * 0.5;
        double b2 = tmat[l1, l2] - b1;
        double b3 = tmat[l1, l3] - b1;
        double[] branch = { b1 - av[l1], b2 - av[l2], b3 - av[l3] };
        for (int i = 0; i < 3; i++)
            if (Math.Abs(branch[i]) < 0.0001) branch[i] = 0.0;
        node[l1].BranchLength = branch[0];
        node[l2].BranchLength = branch[1];
        node[l3].BranchLength = branch[2];
        return new PhylogeneticAnalyzer.PhyloNode
        {
            Children = new List<PhylogeneticAnalyzer.PhyloNode> { node[l1], node[l2], node[l3] },
            Taxa = node[l1].Taxa.Concat(node[l2].Taxa).Concat(node[l3].Taxa).ToList(),
        };
    }

    /// <summary>
    /// <c>fast_nj_tree</c> (util_make_tree.c; the program switches to it for 100 or more sequences): the same neighbour joining with running row / column sums and the
    /// candidate value <c>fnseqs2 dij + 2 sumd - diq - djq</c> compared with a strict <c>&lt;</c> (no <c>MY_EPS</c>).
    /// </summary>
    internal static PhylogeneticAnalyzer.PhyloNode TcNeighborJoiningFast(double[,] distance, int k)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(k, 3);
        var tmat = new double[k + 1, k + 1];
        for (int a = 1; a <= k; a++)
            for (int b = 1; b <= k; b++)
                tmat[a, b] = a == b ? 0.0 : distance[a - 1, b - 1];

        var av = new double[k + 1];
        var sumCols = new double[k + 1];
        var sumRows = new double[k + 1];
        var join = new double[k + 1];
        var valid = Enumerable.Range(1, k).ToList();
        var node = new PhylogeneticAnalyzer.PhyloNode[k + 1];
        for (int i = 1; i <= k; i++)
            node[i] = new PhylogeneticAnalyzer.PhyloNode((i - 1).ToString(System.Globalization.CultureInfo.InvariantCulture))
            {
                Taxa = new List<string> { (i - 1).ToString(System.Globalization.CultureInfo.InvariantCulture) },
            };

        double fnseqs = k, fnseqs2 = 0;
        for (int j = 1; j <= k; j++)
        {
            double tmp = 0.0;
            for (int i = 1; i < j; i++) tmp += tmat[i, j];
            sumRows[j] = tmp;
            tmp = 0.0;
            for (int i = j + 1; i <= k; i++) tmp += tmat[j, i];
            sumCols[j] = tmp;
        }

        for (int nc = 1; nc <= k - 3; ++nc)
        {
            double sumd = 0.0;
            foreach (int j in valid) sumd += sumCols[j];
            fnseqs2 = fnseqs - 2.0;
            double tmin = 99999.0 * 2.0 * fnseqs2;
            int mini = 0, minj = 0;
            foreach (int jj in valid)
            {
                foreach (int ii in valid)
                {
                    if (ii >= jj) break;
                    double diq = sumCols[ii] + sumRows[ii];
                    double djq = sumCols[jj] + sumRows[jj];
                    double dij = tmat[ii, jj];
                    double total = fnseqs2 * dij + 2.0 * sumd - diq - djq;
                    if (total < tmin)
                    {
                        tmin = total;
                        mini = ii;
                        minj = jj;
                    }
                }
            }

            double dio = sumCols[mini] + sumRows[mini];
            double djo = sumCols[minj] + sumRows[minj];
            double dmin = tmat[mini, minj];
            dio = (dio - dmin) / fnseqs2;
            djo = (djo - dmin) / fnseqs2;
            double bi = (dmin + dio - djo) * 0.5;
            double bj = dmin - bi;
            bi -= av[mini];
            bj -= av[minj];
            if (Math.Abs(bi) < 0.0001) bi = 0.0;
            if (Math.Abs(bj) < 0.0001) bj = 0.0;

            node[mini].BranchLength = bi;
            node[minj].BranchLength = bj;
            var joined = new PhylogeneticAnalyzer.PhyloNode
            {
                Children = new List<PhylogeneticAnalyzer.PhyloNode> { node[mini], node[minj] },
                Taxa = node[mini].Taxa.Concat(node[minj].Taxa).ToList(),
            };

            if (dmin <= 0.0) dmin = 0.000001;
            av[mini] = dmin * 0.5;
            fnseqs -= 1.0;
            valid.Remove(minj);

            foreach (int j in valid)
            {
                double tmpDi = 0.0, tmpDj = 0.0;
                if (j < minj)
                {
                    tmpDj = tmat[j, minj];
                    sumCols[j] -= tmpDj;
                }
                else if (j > minj)
                {
                    tmpDj = tmat[minj, j];
                    sumRows[j] -= tmpDj;
                }

                if (j < mini)
                {
                    tmpDi = tmat[j, mini];
                    sumCols[j] -= tmpDi;
                }
                else if (j > mini)
                {
                    tmpDi = tmat[mini, j];
                    sumRows[j] -= tmpDi;
                }

                join[j] = (tmpDj + tmpDi) * 0.5;
            }

            foreach (int j in valid)
            {
                if (j < mini)
                {
                    tmat[j, mini] = join[j];
                    sumCols[j] += join[j];
                }
                else if (j > mini)
                {
                    tmat[mini, j] = join[j];
                    sumRows[j] += join[j];
                }
            }

            sumCols[mini] = sumRows[mini] = 0.0;
            double da = 0.0;
            foreach (int j in valid)
            {
                if (j >= mini) break;
                da += join[j];
            }

            sumRows[mini] = da;
            da = 0.0;
            foreach (int j in valid)
            {
                if (j <= mini) continue;
                da += join[j];
            }

            sumCols[mini] = da;
            sumCols[minj] = sumRows[minj] = 0.0;
            for (int j = 1; j <= k; ++j) tmat[minj, j] = tmat[j, minj] = join[j] = 0.0;
            node[mini] = joined;
        }

        int l1 = valid[0], l2 = valid[1], l3 = valid[2];
        double b1 = (tmat[l1, l2] + tmat[l1, l3] - tmat[l2, l3]) * 0.5;
        double b2 = tmat[l1, l2] - b1;
        double b3 = tmat[l1, l3] - b1;
        double[] branch = { b1 - av[l1], b2 - av[l2], b3 - av[l3] };
        for (int i = 0; i < 3; i++)
            if (Math.Abs(branch[i]) < 0.0001) branch[i] = 0.0;
        node[l1].BranchLength = branch[0];
        node[l2].BranchLength = branch[1];
        node[l3].BranchLength = branch[2];
        return new PhylogeneticAnalyzer.PhyloNode
        {
            Children = new List<PhylogeneticAnalyzer.PhyloNode> { node[l1], node[l2], node[l3] },
            Taxa = node[l1].Taxa.Concat(node[l2].Taxa).Concat(node[l3].Taxa).ToList(),
        };
    }

    /// <summary>
    /// The guide tree the program builds for <paramref name="sequences"/> with its defaults (<c>-tree_mode nj -distance_matrix_mode ktup</c>):
    /// the k-tuple distances, the neighbour-joining tree and the rooting of the trichotomy; the merge order of
    /// <see cref="MultipleAlignTCoffee"/> is the post-order of the returned tree (left subtree first).
    /// </summary>
    public static TCoffeeTreeNode TCoffeeGuideTree(IReadOnlyList<string> sequences, bool? protein = null)
    {
        ArgumentNullException.ThrowIfNull(sequences);
        int k = sequences.Count;
        if (k == 0) throw new ArgumentException("At least one sequence is required.", nameof(sequences));
        if (k == 1) return TCoffeeTreeNode.Leaf(0);
        if (k == 2) return TCoffeeTreeNode.Join(TCoffeeTreeNode.Leaf(0), TCoffeeTreeNode.Leaf(1));

        var sim = TCoffeeKtupSimilarity(sequences, protein ?? TcIsProtein(sequences));
        var d = new double[k, k];
        for (int a = 0; a < k; a++)
            for (int b = 0; b < k; b++)
                d[a, b] = a == b ? 0 : (100 - sim[a, b]) / 100.0; // sim_array2dist_array + normalize_array + / MAXID
        var unrooted = k < 100 ? TcNeighborJoining(d, k) : TcNeighborJoiningFast(d, k);
        return ToTCoffeeTree(ToGuideNode(ClustalWRootNodes(unrooted, k).Root));
    }
}
