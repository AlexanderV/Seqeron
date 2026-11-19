using System.Globalization;

namespace Seqeron.Genomics.Alignment;

// =====================================================================================================================
// T-Coffee progressive stage of the reference implementation (audit A1-13).
//
// The published T-Coffee method (Notredame, Higgins & Heringa 2000) is a pipeline: primary library -> extension ->
// progressive alignment. The reference program (T-Coffee 13.41, cbcrg/tcoffee; the Debian source package t-coffee
// 13.45.0.4846264+really13.41.0.28bdc39+dfsg-1 was read) implements the last two stages differently from the paper's
// formula (sum of extended weights over every residue pair of a column):
//
//   * LIBRARY RELAXATION (`-relax_lib 1`, default; fork_relax_constraint_list in util_lib/util_constraints_list.c): every
//     directed library edge (s1,r1)->(s2,r2) is re-weighted by the triplet consistency
//         w'(s1,r1 -> s2,r2) = trunc( ( w(s2,r2 -> s1,r1) + sum_{(s3,r3) common neighbour} min(w(s1,r1 -> s3,r3), w(s2,r2 -> s3,r3)) )
//                                      / min(deg(s1,r1), deg(s2,r2)) )
//     (double accumulation, integer truncation; the LAST residue of every sequence is skipped by the loop bound `r1 < len`),
//     then edges with weight <= 10 are removed (filter_constraint_list, strict ">" threshold).
//
//   * PROGRESSIVE ALIGNMENT (dp_mode linked_pair_wise, util_dp_gotoh_nw.c): the two profiles of a guide-tree node are
//     aligned by a SPARSE dynamic programme over a list of scored profile-column pairs. Column (p1,p2) is scored by
//     cl2pair_list_ext: for every residue of column p1-1 of profile 1 and every library edge to a residue outside profile 1,
//     the direct weight and every triplet min(w1, w2) landing in a column of profile 2 are added (float arithmetic, the sum
//     divided by the number of such edges, scaled by NORM_F = 1000, kept only when it carries > 1 % of the column's mass and
//     is not in the last row/column); segments are joined with zero-score caps (cl2diag_cap) and the borders (cl2list_borders)
//     are added; the DP (list2linked_pair_wise) is a three-state (M, I, J) maximum-weight chain with gap open = gap extend = 0
//     (`-gapopen 0 -gapext 0` is what the program prints for this mode: "MANUAL PENALTIES: gapopen=0 gapext=0"), ties
//     resolved by the comparison order of the original (M >= I, J ; I >= J ; previous diagonal / row / column entry).
//     The smaller profile (sequences x columns) is always profile 1.
//
// This file is an exact port of those functions for an EXPLICIT primary library and an explicit binary guide tree. What is not
// ported (and therefore not claimed) is how the program builds the primary library (slow_pair = Myers-Miller global alignment
// with the dna_idmat matrix, lalign_id_pair = SIM top-10 local alignments, or the default proba_pair pair-HMM) and its
// default guide tree (k-tuple distances, NJ): see docs/Evidence/ALIGN-MULTI-001-Evidence.md (A1-13).
// =====================================================================================================================

/// <summary>
/// A T-Coffee primary library: weighted residue pairs between sequences (<c>TC_LIB_FORMAT_01</c>). Weights are the program's
/// integers (a pairwise alignment contributes <c>10 x percent identity</c> to each of its aligned pairs, so 1000 is the maximum).
/// </summary>
public sealed class TCoffeeLibrary
{
    private readonly List<(int SeqA, int PosA, int SeqB, int PosB, int Weight)> _pairs = new();

    /// <summary>Number of sequences the library refers to.</summary>
    public int SequenceCount { get; }

    /// <summary>The pairs in insertion order: 0-based sequence indices, 0-based residue positions.</summary>
    public IReadOnlyList<(int SeqA, int PosA, int SeqB, int PosB, int Weight)> Pairs => _pairs;

    /// <summary>Creates an empty library for <paramref name="sequenceCount"/> sequences.</summary>
    public TCoffeeLibrary(int sequenceCount)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(sequenceCount);
        SequenceCount = sequenceCount;
    }

    /// <summary>Adds the pair (<paramref name="seqA"/>[<paramref name="posA"/>], <paramref name="seqB"/>[<paramref name="posB"/>]) with <paramref name="weight"/> (all 0-based).</summary>
    public void Add(int seqA, int posA, int seqB, int posB, int weight)
    {
        if (seqA < 0 || seqA >= SequenceCount) throw new ArgumentOutOfRangeException(nameof(seqA));
        if (seqB < 0 || seqB >= SequenceCount) throw new ArgumentOutOfRangeException(nameof(seqB));
        if (seqA == seqB) throw new ArgumentException("A library pair joins two different sequences.");
        ArgumentOutOfRangeException.ThrowIfNegative(posA);
        ArgumentOutOfRangeException.ThrowIfNegative(posB);
        _pairs.Add((seqA, posA, seqB, posB, weight));
    }
}

/// <summary>A node of an explicit binary guide tree for <see cref="SequenceAligner.MultipleAlignTCoffee"/>: a leaf (sequence index) or a join of two subtrees (left, right).</summary>
public sealed record TCoffeeTreeNode
{
    /// <summary>0-based sequence index of a leaf, -1 for an internal node.</summary>
    public int LeafIndex { get; init; } = -1;

    /// <summary>Left subtree (first profile of the merge); null for a leaf.</summary>
    public TCoffeeTreeNode? Left { get; init; }

    /// <summary>Right subtree (second profile of the merge); null for a leaf.</summary>
    public TCoffeeTreeNode? Right { get; init; }

    /// <summary>A leaf for sequence <paramref name="index"/>.</summary>
    public static TCoffeeTreeNode Leaf(int index) => new() { LeafIndex = index };

    /// <summary>The join of two subtrees (left is the first profile of the merge).</summary>
    public static TCoffeeTreeNode Join(TCoffeeTreeNode left, TCoffeeTreeNode right) => new() { Left = left, Right = right };
}

public static partial class SequenceAligner
{
    // T-Coffee constants (coffee_defines.h)
    private const int TcNormF = 1000;
    private const int TcUndefined = -99999999;
    private const int TcStart = 64;
    private const int TcTerm = 65;
    private const int TcRelaxThreshold = 10;

    /// <summary>
    /// Progressive alignment exactly as the T-Coffee 13.41 program performs it for an explicit primary library and an explicit
    /// binary guide tree (library relaxation, sparse linked-pair-wise DP with zero gap penalties; see the file header).
    /// For every library/tree/sequence set the rows equal those of <c>t_coffee -seq s.fa -lib lib.tc_lib -usetree t.dnd
    /// [-relax_lib 0|1] -output fasta_aln -n_core 1</c> (verified on 1 500 random DNA sets, see the Evidence file).
    /// </summary>
    /// <param name="sequences">The sequences (non-empty strings).</param>
    /// <param name="library">Primary library (weights as in a <c>TC_LIB_FORMAT_01</c> file; <c>10 x percent identity</c> for the program's own methods).</param>
    /// <param name="tree">Binary guide tree over all sequence indices.</param>
    /// <param name="relaxLibrary"><c>true</c> = <c>-relax_lib 1</c> (the program's default), <c>false</c> = <c>-relax_lib 0</c>.</param>
    /// <returns>Aligned rows in input order.</returns>
    public static string[] MultipleAlignTCoffee(
        IReadOnlyList<string> sequences, TCoffeeLibrary library, TCoffeeTreeNode tree, bool relaxLibrary = true)
    {
        ArgumentNullException.ThrowIfNull(sequences);
        ArgumentNullException.ThrowIfNull(library);
        ArgumentNullException.ThrowIfNull(tree);
        int k = sequences.Count;
        if (library.SequenceCount != k)
            throw new ArgumentException("The library was created for a different number of sequences.", nameof(library));
        for (int s = 0; s < k; s++)
        {
            ArgumentNullException.ThrowIfNull(sequences[s]);
            if (sequences[s].Length == 0) throw new ArgumentException("Empty sequences are not supported.", nameof(sequences));
        }

        var leaves = new List<int>();
        CollectTcLeaves(tree, leaves);
        if (leaves.Count != k || leaves.Distinct().Count() != k || leaves.Any(l => l < 0 || l >= k))
            throw new ArgumentException("The guide tree must contain every sequence index exactly once.", nameof(tree));

        var index = new TcIndex(sequences.Select(s => s.Length).ToArray());
        foreach (var (a, pa, b, pb, w) in library.Pairs)
        {
            if (pa >= sequences[a].Length || pb >= sequences[b].Length)
                throw new ArgumentException("A library pair lies outside its sequence.", nameof(library));
            index.Add(a, pa + 1, b, pb + 1, w);
        }

        return TcAlignWithIndex(sequences, index, tree, relaxLibrary);
    }

    /// <summary>The progressive stage on a prepared library index (relaxation, sparse DP along the tree); rows in input order.</summary>
    private static string[] TcAlignWithIndex(IReadOnlyList<string> sequences, TcIndex index, TCoffeeTreeNode tree, bool relaxLibrary)
        => TcAlignDetailed(sequences, index, tree, relaxLibrary ? 1 : 0, 0, null).Rows;

    /// <summary>
    /// The program's order of operations on the library (t_coffee.c L5109-5172): <c>-filter_lib</c> (drop weights &lt;= T), <c>-relax_lib</c> times the relaxation,
    /// the sequence weights (<c>-seq_weight</c>), then the progressive alignment. Returns the rows in input order and the order of the sequences in the program's
    /// own alignment (leaf order of the tree walk, the smaller profile first), which <c>-outorder aligned</c> writes.
    /// </summary>
    private static (string[] Rows, List<int> AlignedOrder) TcAlignDetailed(
        IReadOnlyList<string> sequences, TcIndex index, TCoffeeTreeNode? tree, int relaxCount, int filterLib, float[]? seqWeights, TcPenalties? pen = null, int doNormalise = 1000,
        Func<TcIndex, TcPenalties, (TCoffeeTreeNode Tree, float[]? Weights, int[,]? Sim)>? lateStage = null, TcDpSettings? dp = null)
    {
        int k = sequences.Count;
        if (filterLib > 0) index.Filter(filterLib);
        for (int a = 0; a < relaxCount; a++) index.Relax();
        // evaluate_constraint_list_reference runs after the relaxation and again after the re-weighting: with -do_normalise (1000) the nomatch penalty is rescaled by the
        // largest sum of the weights of a residue's entries each time (integer division)
        var penalties = pen ?? TcPenalties.Default;
        penalties = penalties with { NoMatch = index.NormaliseNoMatch(penalties.NoMatch, doNormalise) };
        // the distance matrix (guide tree and t_coffee sequence weights) is computed here, from the relaxed library, before the library is re-weighted (t_coffee.c L5156-5172)
        int[,]? sim = null;
        if (lateStage != null) (tree, seqWeights, sim) = lateStage(index, penalties);
        if (tree is null) throw new InvalidOperationException("No guide tree.");
        if (seqWeights != null)
        {
            index.Reweight(seqWeights);
            penalties = penalties with { NoMatch = index.NormaliseNoMatch(penalties.NoMatch, doNormalise) };
        }

        var groups = dp is { Mode: not TCoffeeDpMode.Linked }
            ? TcAlignNodeDp(tree, sequences, new TcDpContext(index, sim, dp, penalties.NoMatch, (int)Math.Min(int.MaxValue, index.MaxExtValue())))
            : TcAlignNode(tree, sequences, index, penalties);
        var rows = new string[k];
        for (int s = 0; s < k; s++)
        {
            var cols = groups.Rows[s];
            var chars = new char[cols.Length];
            for (int c = 0; c < cols.Length; c++) chars[c] = cols[c] > 0 ? sequences[s][cols[c] - 1] : '-';
            rows[s] = new string(chars);
        }

        return (rows, groups.Order);
    }

    /// <summary>
    /// <see cref="MsaConsistencyEngine.TCoffeeLinked"/>: this repository's primary library (weights x 10, the program's scale) and
    /// guide tree through the program's progressive stage.
    /// </summary>
    private static MultipleAlignmentResult BuildTCoffeeLinkedResult(
        string[] seqs, Dictionary<ResiduePairKey, double> primary, ProgressiveGuideNode root, PairScoring scoring, MsaGapModel gapModel, float[]? seqWeights = null)
    {
        var library = new TCoffeeLibrary(seqs.Length);
        foreach (var (key, weight) in primary.OrderBy(p => p.Key.SeqA).ThenBy(p => p.Key.SeqB).ThenBy(p => p.Key.PosA).ThenBy(p => p.Key.PosB))
            library.Add(key.SeqA, key.PosA, key.SeqB, key.PosB, (int)(weight * 10));

        string[] rows;
        if (seqWeights is null)
            rows = MultipleAlignTCoffee(seqs, library, ToTCoffeeTree(root));
        else
        {
            // MultipleAlignTCoffee with the program's -seq_weight step (re_weight_constraint_list after the relaxation)
            var index = new TcIndex(seqs.Select(s => s.Length).ToArray());
            foreach (var (a, pa, b, pb, w) in library.Pairs)
                index.Add(a, pa + 1, b, pb + 1, w);
            rows = TcAlignDetailed(seqs, index, ToTCoffeeTree(root), 1, 0, seqWeights).Rows;
        }

        int maxLen = rows.Max(r => r.Length);
        return new MultipleAlignmentResult(rows, BuildConsensus(rows.ToList(), maxLen, scoring.IsMatrix), ComputeSumOfPairsScore(rows.ToList(), scoring, gapModel));
    }

    private static TCoffeeTreeNode ToTCoffeeTree(ProgressiveGuideNode node) =>
        node.IsLeaf ? TCoffeeTreeNode.Leaf(node.LeafIndex) : TCoffeeTreeNode.Join(ToTCoffeeTree(node.Left!), ToTCoffeeTree(node.Right!));

    private static void CollectTcLeaves(TCoffeeTreeNode node, List<int> leaves)
    {
        if (node.Left is null && node.Right is null) { leaves.Add(node.LeafIndex); return; }
        if (node.Left is null || node.Right is null) throw new ArgumentException("The guide tree must be binary.");
        CollectTcLeaves(node.Left, leaves);
        CollectTcLeaves(node.Right, leaves);
    }

    // ----- library ---------------------------------------------------------------------------------------------

    private sealed class TcEdge
    {
        public int Seq, Res, Weight;
    }

    private sealed class TcIndex
    {
        public readonly int[] Len;
        public readonly List<TcEdge>[][] Edges; // [seq][residue 1..len]

        public TcIndex(int[] len)
        {
            Len = len;
            Edges = new List<TcEdge>[len.Length][];
            for (int s = 0; s < len.Length; s++)
            {
                Edges[s] = new List<TcEdge>[len[s] + 1];
                for (int r = 0; r <= len[s]; r++) Edges[s][r] = new List<TcEdge>();
            }
        }

        /// <summary>Number of directed entries (<c>CL->ne</c> after a filter / relaxation).</summary>
        public int EdgeCount
        {
            get
            {
                int n = 0;
                for (int s = 0; s < Len.Length; s++)
                    for (int r = 1; r <= Len[s]; r++) n += Edges[s][r].Count;
                return n;
            }
        }

        public void Add(int s1, int r1, int s2, int r2, int w)
        {
            Edges[s1][r1].Add(new TcEdge { Seq = s2, Res = r2, Weight = w });
            Edges[s2][r2].Add(new TcEdge { Seq = s1, Res = r1, Weight = w });
        }

        /// <summary>The index of the program's own (possibly asymmetric) constraint list: the directed entries exactly as its residue lists hold them.</summary>
        public static TcIndex FromConstraintList(TcConstraintList cl, int[] len)
        {
            var index = new TcIndex(len);
            for (int s = 0; s < len.Length; s++)
                for (int r = 0; r <= len[s]; r++)
                    foreach (var e in cl.Entries(s, r))
                        index.Edges[s][r].Add(new TcEdge { Seq = e.Seq2, Res = e.Res2, Weight = e.Weight });
            return index;
        }

        /// <summary><c>filter_constraint_list(CL, WE, T)</c>: keeps the entries whose weight is greater than T.</summary>
        public void Filter(int threshold)
        {
            for (int s = 0; s < Len.Length; s++)
                for (int r = 1; r <= Len[s]; r++)
                    Edges[s][r].RemoveAll(e => e.Weight <= threshold);
        }

        /// <summary><c>evaluate_constraint_list_reference</c>: max_ext_value = the largest sum of the weights of the entries that end or start at one residue (every directed entry counts for both ends), at least 1;
        /// <c>nomatch</c> becomes <c>nomatch x normalise / max_ext_value</c> when normalising.</summary>
        public int NormaliseNoMatch(int noMatch, int normalise)
        {
            if (normalise == 0) return noMatch;
            return (int)((long)noMatch * normalise / MaxExtValue());
        }

        /// <summary><c>CL->max_ext_value</c> of <c>evaluate_constraint_list_reference</c>: the largest sum of the weights of the entries that end or start at one residue (every directed entry counts for both ends), at least 1.</summary>
        public long MaxExtValue()
        {
            var sum = new long[Len.Length][];
            for (int s = 0; s < Len.Length; s++) sum[s] = new long[Len[s] + 1];
            for (int s = 0; s < Len.Length; s++)
                for (int r = 1; r <= Len[s]; r++)
                    foreach (var e in Edges[s][r])
                    {
                        sum[s][r] += e.Weight;
                        sum[e.Seq][e.Res] += e.Weight;
                    }

            long max = 1;
            foreach (var row in sum) foreach (long v in row) max = Math.Max(max, v);
            return max;
        }

        /// <summary><c>re_weight_constraint_list</c>: every directed entry (s, r) -> (t, r2) is multiplied (float product, truncated to int) by MIN(weight[s], weight[t]).</summary>
        public void Reweight(float[] weight)
        {
            for (int s = 0; s < Len.Length; s++)
                for (int r = 1; r <= Len[s]; r++)
                    foreach (var e in Edges[s][r])
                    {
                        float m = Math.Min(weight[s], weight[e.Seq]);
                        e.Weight = (int)((float)e.Weight * m);
                    }
        }

        /// <summary>fork_relax_constraint_list (normalisation_mode 0) followed by filter_constraint_list(WE, 10).</summary>
        public void Relax()
        {
            int k = Len.Length;
            var update = new int[k][][];
            var hasch = new int[k][];
            for (int s = 0; s < k; s++) hasch[s] = new int[Len[s] + 1];

            for (int s1 = 0; s1 < k; s1++)
            {
                update[s1] = new int[Len[s1] + 1][];
                for (int r1 = 1; r1 < Len[s1]; r1++) // r1 < len: the last residue is not relaxed
                {
                    var own = Edges[s1][r1];
                    int norm1 = own.Count;
                    foreach (var e in own) hasch[e.Seq][e.Res] = e.Weight;
                    var outW = new int[own.Count];
                    for (int a = 0; a < own.Count; a++)
                    {
                        double score = 0;
                        int norm2 = 0;
                        var partner = Edges[own[a].Seq][own[a].Res];
                        foreach (var x in partner)
                        {
                            if (x.Seq == s1 && x.Res == r1) score += (float)x.Weight;
                            else if (hasch[x.Seq][x.Res] != 0)
                                score += Math.Min((float)hasch[x.Seq][x.Res], (float)x.Weight);
                            norm2++;
                        }

                        int norm = Math.Min(norm1, norm2);
                        score = norm != 0 ? score / norm : 0;
                        outW[a] = (int)score;
                    }

                    update[s1][r1] = outW;
                    foreach (var e in own) hasch[e.Seq][e.Res] = 0;
                }
            }

            for (int s1 = 0; s1 < k; s1++)
                for (int r1 = 1; r1 < Len[s1]; r1++)
                {
                    var own = Edges[s1][r1];
                    for (int a = 0; a < own.Count; a++) own[a].Weight = update[s1][r1][a];
                }

            for (int s1 = 0; s1 < k; s1++)
                for (int r1 = 1; r1 <= Len[s1]; r1++)
                    Edges[s1][r1].RemoveAll(e => e.Weight <= TcRelaxThreshold);
        }
    }

    // ----- progressive alignment ------------------------------------------------------------------------------------

    private sealed class TcGroup
    {
        public readonly List<int> Order = new();             // sequence indices in leaf order
        public readonly Dictionary<int, int[]> Rows = new(); // sequence -> residue (1-based) per column, 0 = gap
        public int Columns => Rows[Order.Min()].Length; // strlen (A->seq_al[l_s[0][0]]): the lowest sequence index of the group (lseq is ascending)
    }

    private static TcGroup TcAlignNode(TCoffeeTreeNode node, IReadOnlyList<string> seqs, TcIndex lib, TcPenalties pen)
    {
        if (node.Left is null)
        {
            var g = new TcGroup();
            g.Order.Add(node.LeafIndex);
            g.Rows[node.LeafIndex] = Enumerable.Range(1, seqs[node.LeafIndex].Length).ToArray();
            return g;
        }

        var left = TcAlignNode(node.Left, seqs, lib, pen);
        var right = TcAlignNode(node.Right!, seqs, lib, pen);
        return TcMerge(left, right, lib, pen);
    }

    /// <summary>The gap model of the sparse DP (<c>CL->gop</c>, <c>CL->gep</c>, <c>CL->nomatch</c>): the program's default is 0 / 0 / 0.</summary>
    internal readonly record struct TcPenalties(int GapOpen, int GapExtend, int NoMatch)
    {
        public static TcPenalties Default => new(0, 0, 0);
    }

    private sealed class TcCell
    {
        public int I, J, Diag, Score, Flag;
    }

    private static TcGroup TcMerge(TcGroup left, TcGroup right, TcIndex lib, TcPenalties pen)
    {
        // linked_pair_wise: the group with the smaller (sequences x columns) is profile 1
        var g0 = left; var g1 = right;
        if ((long)g0.Order.Count * g0.Columns > (long)g1.Order.Count * g1.Columns) (g0, g1) = (g1, g0);
        int l1 = g0.Columns, l2 = g1.Columns;

        var list = TcPairListExt(g0, g1, l1, l2, lib);
        var (al0, al1, _) = TcSparseDp(list, l1, l2, pen);

        var merged = new TcGroup();
        merged.Order.AddRange(g0.Order);
        merged.Order.AddRange(g1.Order);
        TcApply(g0, al0, merged);
        TcApply(g1, al1, merged);
        return merged;
    }

    private static void TcApply(TcGroup g, List<bool> al, TcGroup into)
    {
        foreach (var s in g.Order)
        {
            var src = g.Rows[s];
            var dst = new int[al.Count];
            int p = 0;
            for (int c = 0; c < al.Count; c++) dst[c] = al[c] ? src[p++] : 0;
            into.Rows[s] = dst;
        }
    }

    /// <summary>cl2pair_list_ext: scored profile-column pairs (float arithmetic as in the original).</summary>
    private static List<TcCell> TcPairListExt(TcGroup g0, TcGroup g1, int l1, int l2, TcIndex lib)
    {
        var in1 = new HashSet<int>(g0.Order);
        var in2 = new HashSet<int>(g1.Order);
        var inv = new Dictionary<(int, int), int>();
        foreach (var s in g1.Order)
        {
            var row = g1.Rows[s];
            for (int c = 0; c < row.Length; c++) if (row[c] > 0) inv[(s, row[c])] = c + 1;
        }

        // profile_pair_wise walks the group's lseq of the tree node: recode_tree fills it with the sequence indices in ASCENDING order (not the order of the tree walk), and the
        // float sums of the pair list depend on that order (cell (14, 20) of the oracle set of F-A2-7: 21.999998 against 22.000002)
        var sortedG0 = g0.Order.OrderBy(x => x).ToArray();
        var norm = new float[l1 + 1];
        var usedScore = new float[l2 + 1];
        var usedCount = new int[l2 + 1];
        var usedList = new List<int>();
        var list = new List<TcCell>();
        const float filter = 0.01f;

        for (int p1 = 1; p1 <= l1; p1++)
        {
            float tot = 0;
            usedList.Clear();
            foreach (var s in sortedG0)
            {
                int r = g0.Rows[s][p1 - 1];
                if (r <= 0) continue;
                foreach (var e in lib.Edges[s][r])
                {
                    if (in1.Contains(e.Seq)) continue;
                    norm[p1] += 1;
                    // b == 0 is the edge itself, then the edges of (e.Seq, e.Res)
                    TcCount(e.Seq, e.Res, e.Weight);
                    foreach (var x in lib.Edges[e.Seq][e.Res]) TcCount(x.Seq, x.Res, x.Weight);

                    void TcCount(int s2, int r2, int w2)
                    {
                        if (!in2.Contains(s2)) return;
                        int p2 = inv[(s2, r2)];
                        float score = Math.Min((float)e.Weight / TcNormF, (float)w2 / TcNormF);
                        if (usedCount[p2] == 0 && score > 0) usedList.Add(p2);
                        tot += score;
                        usedScore[p2] += score;
                        usedCount[p2]++;
                    }
                }
            }

            foreach (int p2 in usedList)
            {
                float nscore = usedScore[p2] / tot;
                float score = usedScore[p2];
                usedScore[p2] = 0; usedCount[p2] = 0;
                if (nscore > filter && p1 != l1 && p2 != l2)
                {
                    score = norm[p1] > 0 ? score / norm[p1] : 0;
                    score *= TcNormF;
                    // the original writes the value with "%f" to a temporary file and reads it back as a float
                    float back = float.Parse(((double)score).ToString("F6", CultureInfo.InvariantCulture), CultureInfo.InvariantCulture);
                    list.Add(new TcCell { I = p1, J = p2, Diag = l1 - p1 + p2, Score = (int)back });
                }
            }
        }

        return list;
    }

    /// <summary>cl2diag_cap + cl2list_borders + list2nodup_list + list2linked_pair_wise; returns the two alignment masks (true = residue column).</summary>
    // M >= I, J ; else I >= J ; else J (the original's comparison order)
    private static int TcBestState(long m, long i, long j)
    {
        if (m >= i && m >= j) return 0;
        return i >= j ? 1 : 2;
    }

    private static int CompareTcIJ(TcCell a, TcCell b)
    {
        int c = a.I.CompareTo(b.I);
        return c != 0 ? c : a.J.CompareTo(b.J);
    }

    private static int CompareTcJI(TcCell a, TcCell b)
    {
        int c = a.J.CompareTo(b.J);
        return c != 0 ? c : a.I.CompareTo(b.I);
    }

    private static (List<bool>, List<bool>, int) TcSparseDp(List<TcCell> pairList, int l1, int l2, TcPenalties pen)
    {
        var cells = pairList;

        // ---- cl2diag_cap -------------------------------------------------------------------------------------
        cells.Sort((a, b) =>
        {
            int c = a.Diag.CompareTo(b.Diag);
            return c != 0 ? c : a.I.CompareTo(b.I);
        });
        int n = cells.Count;
        var ll = new List<TcCell>();
        for (int a = 0; a < n; a++)
        {
            var e = cells[a];
            if (e.Score == 0) continue;
            int i = e.I, j = e.J;
            int pi = a == 0 ? -10 : cells[a - 1].I, pj = a == 0 ? -10 : cells[a - 1].J;
            int ni = a == n - 1 ? -10 : cells[a + 1].I, nj = a == n - 1 ? -10 : cells[a + 1].J;

            // a segment start: not on the border, not sharing a row/column with the previous entry, previous entry not diagonal-adjacent
            if (i != 0 && j != 0 && i != pi && j != pj && (i - pi != 1 || j - pj != 1)) { ll.Add(e); e.Flag = TcStart; }

            // a segment end (the same test against the next entry); an entry can be both, the last flag wins (it stays twice in `ll`)
            if (i != l1 && j != l2 && i != ni && j != nj && (ni - i != 1 || nj - j != 1)) { ll.Add(e); e.Flag = TcTerm; }
        }

        ll.Sort(CompareTcIJ);
        var extra = new List<TcCell>();
        for (int a = 0; a < ll.Count; a++)
        {
            if (ll[a].Flag != TcTerm) continue;
            int ti = ll[a].I, tj = ll[a].J, ci = ti, nl = 0, bestD = -1, bestS = -1;
            for (int b = a + 1; b < ll.Count && nl < 100; b++)
            {
                if (ll[b].Flag != TcStart) continue;
                int si = ll[b].I, sj = ll[b].J;
                if (si > ci) { nl++; ci = si; }
                int d = Math.Min(si - ti, sj - tj);
                if (d > 0 && (bestD == -1 || bestD > d)) { bestD = d; bestS = b; }
            }

            if (bestD == -1) continue;
            int sI = ll[bestS].I, sJ = ll[bestS].J;
            for (int i = ti, j = tj; i <= sI && j <= sJ; i++, j++) extra.Add(new TcCell { I = i, J = j, Diag = l1 - i + j });
            for (int i = sI, j = sJ; i >= ti && j >= tj; i--, j--) extra.Add(new TcCell { I = i, J = j, Diag = l1 - i + j });
        }

        foreach (var e in ll) e.Flag = 0;
        cells.AddRange(extra);

        // ---- cl2list_borders -----------------------------------------------------------------------------------
        for (int p1 = 0; p1 <= l1; p1++)
        {
            // the border cells carry the score gep x SCORE_K x position (addE of cl2list_borders)
            if (p1 == 0 || p1 == l1)
                for (int p2 = 0; p2 <= l2; p2++) cells.Add(new TcCell { I = p1, J = p2, Diag = l1 - p1 + p2, Score = pen.GapExtend * TcScoreK * p2 });
            else
            {
                cells.Add(new TcCell { I = p1, J = 0, Diag = l1 - p1, Score = pen.GapExtend * TcScoreK * p1 });
                cells.Add(new TcCell { I = p1, J = l2, Diag = l1 - p1 + l2, Score = pen.GapExtend * TcScoreK * p1 });
            }
        }

        // ---- list2nodup_list ------------------------------------------------------------------------------------
        cells.Sort(CompareTcIJ);
        var list = new List<TcCell>(cells.Count);
        foreach (var e in cells)
        {
            if (list.Count > 0 && list[^1].I == e.I && list[^1].J == e.J) list[^1].Score = Math.Max(list[^1].Score, e.Score);
            else list.Add(new TcCell { I = e.I, J = e.J, Diag = e.Diag, Score = e.Score });
        }

        n = list.Count;

        // ---- list2linked_pair_wise ------------------------------------------------------------------------------
        var jOrder = Enumerable.Range(0, n).ToArray();
        Array.Sort(jOrder, (x, y) => CompareTcJI(list[x], list[y]));
        var dOrder = Enumerable.Range(0, n).ToArray();
        Array.Sort(dOrder, (x, y) =>
        {
            int c = list[x].Diag.CompareTo(list[y].Diag);
            return c != 0 ? c : CompareTcIJ(list[x], list[y]);
        });
        var rankJ = new int[n];
        var rankD = new int[n];
        for (int p = 0; p < n; p++) { rankJ[jOrder[p]] = p; rankD[dOrder[p]] = p; }

        // states 0 = M, 1 = I, 2 = J: score, previous entry, residues consumed in profile 1 / profile 2 (fields 2 / 3), predecessor state
        var sc = new long[3][]; var prv = new int[3][]; var f2 = new int[3][]; var f3 = new int[3][]; var typ = new int[3][];
        for (int t = 0; t < 3; t++) { sc[t] = new long[n]; prv[t] = new int[n]; f2[t] = new int[n]; f3[t] = new int[n]; typ[t] = new int[n]; }
        const int M = 0, I = 1, J = 2;
        int start = -1;

        for (int a = 0; a < n; a++)
        {
            int i = list[a].I, j = list[a].J;
            if (i == l1 && j == l2) start = a;
            else if (i == 0 || j == 0)
            {
                sc[M][a] = -1000000;
                if (j == 0) { sc[J][a] = -10000000; sc[I][a] = (long)pen.GapExtend * i; }
                else { sc[I][a] = -10000000; sc[J][a] = (long)pen.GapExtend * j; }
                for (int t = 0; t < 3; t++) { prv[t][a] = -1; f2[t][a] = i; f3[t][a] = j; }
                continue;
            }

            int pi = jOrder[rankJ[a] - 1];
            int pj = a - 1;
            int ij = dOrder[rankD[a] - 1];

            // a gap that touches the last row or column pays no opening penalty
            int gop = i == l1 || j == l2 ? 0 : pen.GapOpen;
            int deltaI = i - list[pi].I, deltaJ = j - list[pj].J;
            sc[I][a] = Math.Max(sc[I][pi], sc[M][pi] + gop) + (long)deltaI * pen.GapExtend;
            prv[I][a] = pi; f2[I][a] = deltaI; f3[I][a] = 0;
            typ[I][a] = sc[I][pi] >= sc[M][pi] + gop ? I : M;

            sc[J][a] = Math.Max(sc[J][pj], sc[M][pj] + gop) + (long)deltaJ * pen.GapExtend;
            prv[J][a] = pj; f2[J][a] = 0; f3[J][a] = deltaJ;
            typ[J][a] = sc[J][pj] >= sc[M][pj] + gop ? J : M;

            int ls = i - list[ij].I;
            if (a > 1 && ls == j - list[ij].J)
            {
                sc[M][a] = Math.Max(sc[M][ij], Math.Max(sc[I][ij], sc[J][ij])) + list[a].Score - (long)ls * pen.NoMatch;
                prv[M][a] = ij; f2[M][a] = f3[M][a] = ls;
                typ[M][a] = TcBestState(sc[M][ij], sc[I][ij], sc[J][ij]);
            }
            else
            {
                sc[M][a] = TcUndefined; prv[M][a] = -1;
            }
        }

        // traceback
        int cur = start;
        int linkedScore = (int)Math.Max(sc[M][cur], Math.Max(sc[I][cur], sc[J][cur])); // list2linked_pair_wise: score = MAX3(MM, MI, MJ) of the end cell (a long returned as int)
        int state = TcBestState(sc[M][cur], sc[I][cur], sc[J][cur]);
        int ri = l1, rj = l2;
        var al0 = new List<bool>(); var al1 = new List<bool>();
        while (!(ri == 0 && rj == 0))
        {
            int rows1 = f2[state][cur], rows2 = f3[state][cur];
            int l = Math.Max(rows1, rows2);
            if (ri == 0) { while (rj > 0) { al0.Add(false); al1.Add(true); rj--; } }
            else if (rj == 0) { while (ri > 0) { al0.Add(true); al1.Add(false); ri--; } }
            else if (l == 0) throw new InvalidOperationException("T-Coffee DP traceback reached an empty step.");
            else
            {
                for (int b = 0; b < l; b++)
                {
                    if (rows1 != 0) { al0.Add(true); ri--; } else al0.Add(false);
                    if (rows2 != 0) { al1.Add(true); rj--; } else al1.Add(false);
                }

                int next = prv[state][cur];
                state = typ[state][cur];
                cur = next;
            }
        }

        al0.Reverse(); al1.Reverse();
        return (al0, al1, linkedScore);
    }

    /// <summary>
    /// <c>linked_pair_wise</c> of the two single sequences a and b (the distance modes slow / cscore / clustalw): the sparse DP on the library <paramref name="lib"/> (the shorter sequence is profile 1),
    /// the aligned rows (as residue masks for a and b) and the DP score the program returns.
    /// </summary>
    private static (string RowA, string RowB, int Score) TcLinkedPairRows(IReadOnlyList<string> seqs, int a, int b, TcIndex lib, TcPenalties pen)
    {
        var ga = new TcGroup(); ga.Order.Add(a); ga.Rows[a] = Enumerable.Range(1, seqs[a].Length).ToArray();
        var gb = new TcGroup(); gb.Order.Add(b); gb.Rows[b] = Enumerable.Range(1, seqs[b].Length).ToArray();
        var (g0, g1) = (ga, gb);
        if (g0.Columns > g1.Columns) (g0, g1) = (g1, g0);
        var list = TcPairListExt(g0, g1, g0.Columns, g1.Columns, lib);
        var (al0, al1, score) = TcSparseDp(list, g0.Columns, g1.Columns, pen);
        var (r0, r1) = TcMasksToRows(seqs[g0.Order[0]], seqs[g1.Order[0]], al0, al1);
        return g0 == ga ? (r0, r1, score) : (r1, r0, score);
    }
}
