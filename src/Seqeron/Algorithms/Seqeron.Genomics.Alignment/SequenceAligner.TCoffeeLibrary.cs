using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Alignment;

// =====================================================================================================================
// The T-Coffee 13.41 PRIMARY LIBRARY as the program's constraint list holds it (audit A1-13b, F66).
//
// A library is not simply the union of the methods' pair alignments. The program keeps for every residue (s, r) a SORTED list of directed
// entries (s2, r2, weight, count) (`residue_index`, util_constraints_list.c) and builds it with forked processes:
//
//   fork_read_n_constraint_list : one child per METHOD, started from the library merged so far; the child returns (`dump`) its whole list when it
//                                 grew (`CL->ne > in`) and the parent merges it (`undump` = add_entry2list2: equal key -> weight = MAX, count += count;
//                                 new key -> inserted at its sorted place).
//   fork_subset_produce_list    : inside a method child, one further child runs all pair jobs (row-major pairs, measured: the result does not depend on
//                                 `-n_core`, 2, 3 and 4 cores equal 1 core on 150 sets); at its start the child records the length of every residue
//                                 list (`freeze_constraint_list`) and at its end dumps ONLY the entries at list positions >= that length
//                                 (`safe_dump_constraint_list`, `b = freeze[s1][r1]`) - positions, not "the entries it added". The first method
//                                 starts from empty lists, so nothing is lost; for the second and later methods an entry the child inserted BEFORE
//                                 an older one is not returned (its place is below the recorded length), an old entry that moved into the tail is returned
//                                 again, and updates of old entries are lost unless they sit in the tail. The in-memory list therefore becomes
//                                 asymmetric (the directed entry of one side can be missing); `-out_lib` prints the s1 &lt; s2 side only.
//
// This is the reason why `t_coffee -method slow_pair,lalign_id_pair` is not the union of the two single-method libraries. The model below performs
// exactly these operations on the sorted lists (a binary-search insertion is order independent, so a sorted insertion is the same operation).
// =====================================================================================================================

public static partial class SequenceAligner
{
    private const int TcMaxId = 100;
    private const string TcIdentityLetters = "ABCDEFGHIKLMNPQRSTUVWXYZ";

    /// <summary><c>is_in_same_group_aa(r1, r2, 0, NULL, "")</c>: groups of the identity matrix (every letter of the 23-letter alphabet alone, T and U together); '-' is ignored by the caller.</summary>
    private static bool TcSameGroup(char a, char b)
    {
        a = char.ToUpperInvariant(a);
        b = char.ToUpperInvariant(b);
        if (TcIdentityLetters.IndexOf(a) < 0 || TcIdentityLetters.IndexOf(b) < 0) return false;
        if (a == 'U') a = 'T';
        if (b == 'U') b = 'T';
        return a == b;
    }

    /// <summary><c>get_seq_sim(s1, s2, "-", NULL)</c>: 100 x identical columns / columns with a residue in both rows (integer division), 0 when there is none.</summary>
    private static int TcSeqSim(string r1, string r2)
    {
        int pos0 = 0, sim = 0;
        for (int a = 0; a < r1.Length; a++)
        {
            if (r1[a] == '-' || r2[a] == '-') continue;
            pos0++;
            if (TcSameGroup(r1[a], r2[a])) sim++;
        }

        return pos0 == 0 ? 0 : sim * 100 / pos0;
    }

    /// <summary><c>fix_seq_seq</c> for one sequence: index[r] = coordinate (1-based) in <paramref name="full"/> of the r-th (1-based) residue of <paramref name="segment"/>.</summary>
    private static int[] TcFixIndex(string segment, string full)
    {
        var index = new int[segment.Length + 1];
        if (string.Equals(segment, full, StringComparison.OrdinalIgnoreCase))
        {
            for (int r = 1; r <= segment.Length; r++) index[r] = r;
            return index;
        }

        // align_two_sequences(..., "idmat" | "blosum62mt", -4, -1): the matrix follows the type of the sequences built by aln2seq, which is never set -> idmat
        var (ra, rb) = TcMyersMiller(segment.ToLowerInvariant(), full.ToLowerInvariant(), TcMatrix.IdMat, -4, -1);
        int nr0 = 0, nr1 = 0;
        for (int c = 0; c < ra.Length; c++)
        {
            bool g0 = ra[c] == '-', g1 = rb[c] == '-';
            if (!g0) nr0++;
            if (!g1) nr1++;
            if (!g0 && !g1) index[nr0] = nr1;
        }

        return index;
    }

    /// <summary>One directed library entry in the program's representation (residues 1-based, 0 = "no residue").</summary>
    internal readonly record struct TcEntry(int Seq1, int Res1, int Seq2, int Res2, int Weight, int Count);

    /// <summary>The program's constraint list: per residue a list sorted by (sequence, residue) of directed entries.</summary>
    internal sealed class TcConstraintList
    {
        private readonly int[] _len;
        private readonly List<(int Seq2, int Res2, int Weight, int Count)>[][] _lists;

        public int EntryCount { get; private set; }

        public TcConstraintList(int[] lengths)
        {
            _len = lengths;
            _lists = new List<(int, int, int, int)>[lengths.Length][];
            for (int s = 0; s < lengths.Length; s++)
            {
                _lists[s] = new List<(int, int, int, int)>[lengths[s] + 1];
                for (int r = 0; r <= lengths[s]; r++) _lists[s][r] = new List<(int, int, int, int)>();
            }
        }

        public TcConstraintList Clone()
        {
            var c = new TcConstraintList(_len) { EntryCount = EntryCount };
            for (int s = 0; s < _len.Length; s++)
                for (int r = 0; r <= _len[s]; r++)
                    c._lists[s][r].AddRange(_lists[s][r]);
            return c;
        }

        /// <summary><c>add_entry2list2</c>: one directed entry; an equal key keeps the larger weight and adds the counts.</summary>
        public void AddDirected(int fromSeq, int fromRes, int toSeq, int toRes, int weight, int count)
        {
            var list = _lists[fromSeq][fromRes];
            int s2 = toSeq, r2 = toRes;
            int lo = 0, hi = list.Count;
            while (lo < hi)
            {
                int mid = (lo + hi) >> 1;
                var e = list[mid];
                if (e.Seq2 < s2 || (e.Seq2 == s2 && e.Res2 < r2)) lo = mid + 1;
                else hi = mid;
            }

            if (lo < list.Count && list[lo].Seq2 == s2 && list[lo].Res2 == r2)
            {
                var old = list[lo];
                list[lo] = (s2, r2, Math.Max(weight, old.Weight), old.Count + count);
            }
            else
            {
                list.Insert(lo, (s2, r2, weight, count));
                EntryCount++;
            }
        }

        /// <summary><c>add_entry2list</c>: the mirror entry first, then the entry.</summary>
        public void AddBoth(int s1, int r1, int s2, int r2, int weight, int count)
        {
            AddDirected(fromSeq: s2, fromRes: r2, toSeq: s1, toRes: r1, weight: weight, count: count);
            AddDirected(s1, r1, s2, r2, weight, count);
        }

        /// <summary><c>freeze_constraint_list</c>: the current length of every residue list.</summary>
        public int[][] Freeze()
        {
            var f = new int[_len.Length][];
            for (int s = 0; s < _len.Length; s++)
            {
                f[s] = new int[_len[s] + 1];
                for (int r = 1; r <= _len[s]; r++) f[s][r] = _lists[s][r].Count;
            }

            return f;
        }

        /// <summary><c>safe_dump_constraint_list</c>: the entries at list positions &gt;= the frozen length (all when <paramref name="frozen"/> is null) with a residue &gt; 0.</summary>
        public List<TcEntry> Dump(int[][]? frozen)
        {
            var d = new List<TcEntry>();
            if (EntryCount == 0) return d;
            for (int s = 0; s < _len.Length; s++)
                for (int r = 1; r <= _len[s]; r++)
                {
                    var list = _lists[s][r];
                    for (int b = frozen?[s][r] ?? 0; b < list.Count; b++)
                    {
                        var e = list[b];
                        if (e.Res2 <= 0) continue;
                        d.Add(new TcEntry(s, r, e.Seq2, e.Res2, e.Weight, e.Count));
                    }
                }

            return d;
        }

        public void Undump(IEnumerable<TcEntry> entries)
        {
            foreach (var e in entries) AddDirected(e.Seq1, e.Res1, e.Seq2, e.Res2, e.Weight, e.Count);
        }

        /// <summary>The directed entries of residue (s, r) in list order.</summary>
        public IReadOnlyList<(int Seq2, int Res2, int Weight, int Count)> Entries(int s, int r) => _lists[s][r];

        public int Length(int s) => _len[s];
    }

    /// <summary>The entries one pair job (x, y) adds, in the order of <c>aln2constraint_list</c> (alignments of the chain, columns left to right).</summary>
    private static List<(int S1, int R1, int S2, int R2, int W)> TcPairJobEntries(
        IReadOnlyList<string> seqs, int x, int y, TCoffeePairMethod method, int lalignTop, TcSeqKind kind, string? matrix, TCoffeeMethodParameters? mp)
    {
        bool protein = kind == TcSeqKind.Protein;
        if (method == TCoffeePairMethod.ProbaPair)
            return TcProbaPairEntries(seqs[x], seqs[y], kind, matrix).Select(e => (x, e.R1, y, e.R2, e.Weight)).ToList();

        var rows = new List<(string, string)>();
        switch (method)
        {
            case TCoffeePairMethod.SlowPair:
                rows.Add(TCoffeeSlowPairAlign(seqs[x], seqs[y], protein, matrix, mp));
                break;
            case TCoffeePairMethod.LalignIdPair:
                {
                    var alns = TCoffeeLalignTopAlignments(seqs[x], seqs[y], lalignTop, protein, mp);
                    if (alns.Count == 0)
                    {
                        // sim_pair_wise_lalign found nothing: the program then turns the padded (unaligned) pair into the library
                        int len = Math.Max(seqs[x].Length, seqs[y].Length);
                        rows.Add((seqs[x].PadRight(len, '-'), seqs[y].PadRight(len, '-')));
                    }

                    foreach (var a in alns) rows.Add((a.RowA, a.RowB));
                    break;
                }
            case TCoffeePairMethod.ClustalwPair:
                {
                    var cw = protein
                        ? MultipleAlignClustalW(new[] { seqs[x], seqs[y] }, new ClustalWProteinOptions())
                        : MultipleAlignClustalW(new[] { seqs[x], seqs[y] });
                    rows.Add((cw.AlignedSequences[0], cw.AlignedSequences[1]));
                    break;
                }

            default:
                throw new NotSupportedException($"The {method} library is not available here.");
        }

        var entries = new List<(int, int, int, int, int)>();
        foreach (var (rowA, rowB) in rows)
        {
            int weight = TcSeqSim(rowA, rowB);
            if (weight == 0) weight = 1;
            weight *= 1000 / TcMaxId;
            int[] fa = TcFixIndex(rowA.Replace("-", ""), seqs[x]);
            int[] fb = TcFixIndex(rowB.Replace("-", ""), seqs[y]);
            int na = 0, nb = 0;
            for (int c = 0; c < rowA.Length; c++)
            {
                bool ga = rowA[c] == '-', gb = rowB[c] == '-';
                if (!ga) na++;
                if (!gb) nb++;
                if (!ga && !gb) entries.Add((x, fa[na], y, fb[nb], weight));
            }
        }

        return entries;
    }

    /// <summary>
    /// What the child process of one METHOD returns (<c>fork_read_n_constraint_list</c> child): its start list (the library merged so far) is cloned, the pair
    /// jobs run in a grand-child whose tail is merged back, and the whole list is returned when it grew.
    /// </summary>
    private static List<TcEntry>? TcRunMethod(
        TcConstraintList start, IReadOnlyList<string> seqs, List<(int X, int Y)> pairs, TCoffeePairMethod method, int lalignTop, TcSeqKind kind, string? matrix, TCoffeeMethodParameters? mp)
    {
        var methodList = start.Clone();
        int before = methodList.EntryCount;
        var child = methodList.Clone();
        var frozen = child.Freeze();
        foreach (var (x, y) in pairs)
            foreach (var (s1, r1, s2, r2, w) in TcPairJobEntries(seqs, x, y, method, lalignTop, kind, matrix, mp))
                child.AddBoth(s1, r1, s2, r2, w, 1);
        methodList.Undump(child.Dump(frozen));
        return methodList.EntryCount > before ? methodList.Dump(null) : null;
    }

    /// <summary>
    /// The program's constraint list after the pair-wise methods <paramref name="methods"/> were read in this order
    /// (see the file header for the process structure that is reproduced). A method child starts from the list merged so far.
    /// </summary>
    internal static TcConstraintList TcBuildConstraintList(
        IReadOnlyList<string> seqs, IReadOnlyList<TCoffeePairMethod> methods, int lalignTop, TcSeqKind kind, string? matrix = null, IReadOnlyDictionary<TCoffeePairMethod, TCoffeeMethodParameters>? methodParameters = null)
    {
        var lengths = seqs.Select(s => s.Length).ToArray();
        var main = new TcConstraintList(lengths);
        var pairs = new List<(int X, int Y)>();
        for (int x = 0; x < seqs.Count; x++)
            for (int y = x + 1; y < seqs.Count; y++)
                pairs.Add((x, y));

        // a method named twice is read once
        foreach (var method in methods.Distinct())
        {
            var dump = TcRunMethod(main, seqs, pairs, method, lalignTop, kind, matrix, methodParameters != null && methodParameters.TryGetValue(method, out var mp) ? mp : null);
            if (dump != null) main.Undump(dump);
        }

        return main;
    }

    /// <summary>
    /// The primary library the T-Coffee 13.41 program builds for <paramref name="sequences"/> with the given pair-wise methods
    /// (<c>t_coffee -seq s.fa -method m1,m2 -lib_only -out_lib x.tc_lib</c>, independent of <c>-n_core</c>): pairs (printed side: sequence index of the first
    /// sequence smaller than that of the second) with weight 10 x percent identity of the alignment they come from, and the number of
    /// times the pair was added. With two or more methods the program's own merge loses and duplicates entries (see the file header);
    /// the result is what it writes.
    /// </summary>
    /// <param name="sequences">The sequences.</param>
    /// <param name="methods">The pair-wise methods in command-line order.</param>
    /// <param name="lalignTopAlignments"><c>-lalign_n_top</c> (default 10).</param>
    /// <param name="protein"><c>null</c> = decide from the letters (as the program does for the parameters it selects).</param>
    public static IReadOnlyList<TCoffeeLibraryEntry> TCoffeePrimaryLibraryEntries(
        IReadOnlyList<string> sequences, IEnumerable<TCoffeePairMethod> methods, int lalignTopAlignments = 10, bool? protein = null)
    {
        var cl = TcBuildChecked(sequences, methods, lalignTopAlignments, TcKindOf(sequences, protein));
        var list = new List<TCoffeeLibraryEntry>();
        for (int s1 = 0; s1 < sequences.Count; s1++)
            for (int r1 = 1; r1 <= cl.Length(s1); r1++)
                foreach (var e in cl.Entries(s1, r1))
                    if (s1 < e.Seq2 && e.Res2 > 0)
                        list.Add(new TCoffeeLibraryEntry(s1, r1 - 1, e.Seq2, e.Res2 - 1, e.Weight, e.Count));
        return list.OrderBy(e => e.SeqA).ThenBy(e => e.SeqB).ThenBy(e => e.PosA).ThenBy(e => e.PosB).ToList();
    }

    private static TcConstraintList TcBuildChecked(IReadOnlyList<string> sequences, IEnumerable<TCoffeePairMethod> methods, int lalignTop, TcSeqKind kind, string? matrix = null, IReadOnlyDictionary<TCoffeePairMethod, TCoffeeMethodParameters>? methodParameters = null)
    {
        ArgumentNullException.ThrowIfNull(sequences);
        ArgumentNullException.ThrowIfNull(methods);
        foreach (var s in sequences)
            if (string.IsNullOrEmpty(s)) throw new ArgumentException("Empty sequences are not supported.", nameof(sequences));
        return TcBuildConstraintList(sequences, methods.ToList(), lalignTop, kind, matrix, methodParameters);
    }
}

/// <summary><c>-tree_mode</c> of the program.</summary>
public enum TCoffeeTreeMode
{
    /// <summary>Neighbour joining (ClustalW 1.x routine) with the trichotomy at the root turned into a binary tree; default.</summary>
    Nj = 0,

    /// <summary>UPGMA (<c>upgma_merge</c>) on the integer distances.</summary>
    Upgma = 1,
}

/// <summary><c>-distance_matrix_mode</c> of the program (the pair-wise distances of the guide tree and of the sequence weights).</summary>
public enum TCoffeeDistanceMode
{
    /// <summary><c>ktup</c> (default): k-tuple identity of the reduced alphabet (<c>ktup_dist_mat</c>).</summary>
    Ktup = 0,

    /// <summary><c>fast</c>: Myers-Miller alignment (identity matrix, DNA -10 / -1; BLOSUM62 protein) and its identity.</summary>
    Fast = 1,

    /// <summary><c>idscore</c>: the identity of <c>idscore_pairseq</c> (affine DP; note the program does not scale this score by 10 like the others).</summary>
    Idscore = 2,

    /// <summary><c>no</c>: every pair has identity 100 (all distances 0).</summary>
    No = 3,
}

/// <summary>Options of <see cref="SequenceAligner.MultipleAlignTCoffeeProgram"/> = the command line of the program they mirror.</summary>
public sealed partial record TCoffeeProgramOptions
{
    /// <summary><c>-method m1,m2</c>: the primary library methods in command-line order. Default = <c>proba_pair</c>, the program's own default of 13.41 (<c>t_coffee seqs.fa</c>).</summary>
    public IReadOnlyList<TCoffeePairMethod> Methods { get; init; } = new[] { TCoffeePairMethod.ProbaPair };

    /// <summary><c>-lalign_n_top</c>: alignments reported by <c>lalign_id_pair</c> per pair (default 10).</summary>
    public int LalignTopAlignments { get; init; } = 10;

    /// <summary><c>-relax_lib 1</c> (default) or <c>0</c>; <see cref="RelaxLibIterations"/> overrides it.</summary>
    public bool RelaxLibrary { get; init; } = true;

    /// <summary><c>-relax_lib N</c>: the number of times the library is relaxed (<c>null</c> = 1 when <see cref="RelaxLibrary"/>, else 0).</summary>
    public int? RelaxLibIterations { get; init; }

    /// <summary><c>-filter_lib T</c>: before the relaxation every library entry with a weight of at most T is removed (0 = off).</summary>
    public int FilterLib { get; init; }

    /// <summary><c>-tree_mode</c>: <see cref="TCoffeeTreeMode.Nj"/> (default) or <see cref="TCoffeeTreeMode.Upgma"/>.</summary>
    public TCoffeeTreeMode TreeMode { get; init; } = TCoffeeTreeMode.Nj;

    /// <summary><c>-distance_matrix_mode</c> (guide tree and sequence weights): the default is <see cref="TCoffeeDistanceMode.Ktup"/>.</summary>
    public TCoffeeDistanceMode DistanceMatrixMode { get; init; } = TCoffeeDistanceMode.Ktup;

    /// <summary><c>-distance_matrix_sim_mode</c>: the identity measure of the <c>fast</c> / <c>idscore</c> distances (<c>idmat_sim1</c> default; <c>sim1</c>, <c>sim2</c>, <c>sim3</c>, <c>cov</c>, <c>gap1</c>, <c>logid</c>, optionally with an <c>idmat_</c> prefix).</summary>
    public string DistanceMatrixSimMode { get; init; } = "idmat_sim1";

    /// <summary>The <c>method@EP@KEY@value</c> overrides (GOP, GEP, MATRIX) of <c>slow_pair</c> / <c>lalign_id_pair</c> by method; the first occurrence of a method counts, as the program reads a method once.</summary>
    public IReadOnlyDictionary<TCoffeePairMethod, TCoffeeMethodParameters>? MethodParameters { get; init; }

    /// <summary>
    /// <c>-matrix NAME</c> (any matrix of the program, e.g. blosum50mt, pam250mt): the PROTEIN <c>slow_pair</c> matrix (and its gap open from the matrix) and the emission probabilities of
    /// <c>proba_pair</c> (<c>get_emitPairs</c>: p[a][b] = s[a] s[b] exp(score x ln2 / 2), normalised); <c>lalign_id_pair</c>, <c>slow_pair</c> of DNA / RNA and the distance modes keep their own matrices. <c>null</c> = default.
    /// </summary>
    public string? Matrix { get; init; }

    /// <summary><c>-seq_weight t_coffee</c> weights every sequence by its similarity to the others (the library weights are multiplied by the smaller of the two sequence weights); default none.</summary>
    public bool TCoffeeSequenceWeights { get; init; }

    /// <summary><c>-seq_weight file</c>: explicit weights per sequence in input order (a file's weights; sequences missing from it have weight 1). Overrides <see cref="TCoffeeSequenceWeights"/>.</summary>
    public IReadOnlyList<float>? SequenceWeights { get; init; }

    /// <summary><c>-gapopen</c> (default 0): gap opening penalty of the progressive DP, in library weight units (best match 1000; a gap touching the last row or column is free).</summary>
    public int GapOpen { get; init; }

    /// <summary><c>-gapext</c> (default 0): gap extension penalty per residue of the progressive DP.</summary>
    public int GapExtend { get; init; }

    /// <summary><c>-cosmetic_penalty</c> (default 0): the gap opening penalty used when neither <c>-gapopen</c> nor <c>-gapext</c> is given.</summary>
    public int CosmeticPenalty { get; init; }

    /// <summary><c>-nomatch</c> (default 0): penalty per aligned column of the progressive DP (rescaled by <see cref="DoNormalise"/> / the largest residue weight sum).</summary>
    public int NoMatch { get; init; }

    /// <summary><c>-do_normalise</c> (default 1000): the factor of the <see cref="NoMatch"/> rescaling; 0 = no rescaling.</summary>
    public int DoNormalise { get; init; } = 1000;

    /// <summary><c>-outorder aligned</c>: the rows come in the program's internal order (the sequences sorted by <see cref="Names"/>, byte order) instead of the input order (default).</summary>
    public bool OutOrderAligned { get; init; }

    /// <summary>
    /// The sequence names (in input order). The program (<c>-inorder aligned</c>, its default) sorts the sequences by name (byte order) BEFORE building the library and the tree,
    /// which changes tie-breaking; with <c>null</c> the names are taken to sort in the input order (as q000, q001, ...). Ignored when <see cref="InOrderInput"/>.
    /// </summary>
    public IReadOnlyList<string>? Names { get; init; }

    /// <summary><c>-inorder input</c>: keep the input order instead of sorting the sequences by name.</summary>
    public bool InOrderInput { get; init; }

    /// <summary><c>true</c> = protein parameters, <c>false</c> = DNA/RNA, <c>null</c> = decided from the letters.</summary>
    public bool? Protein { get; init; }
}

public static partial class SequenceAligner
{
    /// <summary>
    /// The multiple alignment the T-Coffee 13.41 program (<c>t_coffee</c>, Notredame, Higgins &amp; Heringa 2000) produces for <paramref name="sequences"/>
    /// with its own primary library (<c>proba_pair</c> by default = <c>t_coffee seqs.fa</c>; <c>slow_pair</c>, <c>lalign_id_pair</c> and any combination,
    /// merged the way the program merges them), library relaxation, the k-tuple neighbour-joining guide tree and the sparse linked-pair-wise progressive
    /// stage (<see cref="MultipleAlignTCoffee(IReadOnlyList{string}, TCoffeeLibrary, TCoffeeTreeNode, bool)"/>). The rows equal
    /// <c>t_coffee -seq s.fa [-method m1,m2] -output fasta_aln</c> (any <c>-n_core</c>), see the Evidence file for the committed oracle sets.
    /// </summary>
    /// <remarks>
    /// The sequence type is decided as the program does (DNA when more than 98 % of the letters are <c>AGCTUNRYMKSWHBVD</c>, otherwise protein) and selects the
    /// matrices, gap penalties, pair-HMM and k-tuple alphabet. Letters are used as given (the program keeps the case of the input in the library; both cases score alike).
    /// RNA input (more than 98 % nucleotide letters and at least one <c>U</c>) uses the program's separate RNA pair-HMM (5 matrix types, 2 insert state pairs,
    /// <c>get_tot_prob</c> with the residues themselves instead of the 20-letter profile) for <see cref="TCoffeePairMethod.ProbaPair"/>.
    /// The reported <see cref="MultipleAlignmentResult.TotalScore"/> is the natural affine sum of pairs under
    /// <paramref name="scoring"/> (default EDNAFULL 10 / 0.5 for DNA, BLOSUM62 10 / 0.2 for proteins), not the program's objective.
    /// </remarks>
    /// <param name="sequences">The sequences (letters, none empty).</param>
    /// <param name="options">Options; <c>null</c> = the program's defaults (<c>proba_pair</c>, relaxation on).</param>
    /// <param name="scoring">Scoring of the reported total score only.</param>
    /// <returns>Aligned rows in input order, consensus, total score.</returns>
    public static MultipleAlignmentResult MultipleAlignTCoffeeProgram(IEnumerable<string> sequences, TCoffeeProgramOptions? options = null, MatrixScoring? scoring = null)
    {
        ArgumentNullException.ThrowIfNull(sequences);
        var list = sequences.ToList();
        if (list.Count == 0) return MultipleAlignmentResult.Empty;
        var rows = TCoffeeProgramRows(list, options);
        bool protein = (options ?? new TCoffeeProgramOptions()).Protein ?? TcIsProtein(list);
        var sp = scoring ?? (protein ? MatrixScoring.Emboss(SubstitutionMatrix.Blosum62, 10, 0.2) : MatrixScoring.Emboss(SubstitutionMatrix.EdnaFull, 10, 0.5));
        string consensus = BuildConsensus(rows.ToList(), rows.Max(r => r.Length), anySymbol: true);
        string[] scored = protein ? rows : rows.Select(r => r.ToUpperInvariant().Replace('U', 'T')).ToArray();
        double exact = rows.Length == 1 ? 0.0 : MatrixSumOfPairsScore(scored, sp, MsaGapModel.Affine);
        return new MultipleAlignmentResult(rows, consensus, (int)Math.Round(exact, MidpointRounding.AwayFromZero)) { ExactTotalScore = exact };
    }

    /// <summary>The rows of <see cref="MultipleAlignTCoffeeProgram(IEnumerable{string}, TCoffeeProgramOptions?, MatrixScoring?)"/> without the scoring (input order).</summary>
    internal static string[] TCoffeeProgramRows(IReadOnlyList<string> sequences, TCoffeeProgramOptions? options = null)
    {
        options ??= new TCoffeeProgramOptions();
        ArgumentNullException.ThrowIfNull(sequences);
        if (sequences.Count == 0) return Array.Empty<string>();
        foreach (var s in sequences)
            if (string.IsNullOrEmpty(s) || s.Any(c => !char.IsAsciiLetter(c))) throw new ArgumentException("Sequences must be non-empty strings of letters.", nameof(sequences));
        if (sequences.Count == 1) return new[] { sequences[0] };
        int k = sequences.Count;

        // -inorder aligned (the program's default): the sequences are sorted by name (byte order) before anything else is computed
        var perm = Enumerable.Range(0, k).ToArray();
        if (!options.InOrderInput && options.Names != null)
        {
            if (options.Names.Count != k) throw new ArgumentException("Names must have one entry per sequence.", nameof(options));
            perm = perm.OrderBy(i => options.Names[i], StringComparer.Ordinal).ToArray();
        }

        var seqs = perm.Select(i => sequences[i]).ToList();
        var kind = TcKindOf(seqs, options.Protein);
        bool matrixGroups = options.Methods.Any(m => m != TCoffeePairMethod.ProbaPair);
        var cl = TcBuildChecked(seqs, options.Methods, options.LalignTopAlignments, kind, options.Matrix, options.MethodParameters);
        var index = TcIndex.FromConstraintList(cl, seqs.Select(s => s.Length).ToArray());
        var tree = options.DistanceMatrixMode == TCoffeeDistanceMode.Ktup && options.TreeMode == TCoffeeTreeMode.Nj
            ? TCoffeeGuideTree(seqs, kind == TcSeqKind.Protein)
            : TcGuideTreeFor(seqs, kind, options, matrixGroups);

        float[]? weights = null;
        if (options.SequenceWeights != null)
        {
            if (options.SequenceWeights.Count != k) throw new ArgumentException("SequenceWeights must have one entry per sequence.", nameof(options));
            weights = perm.Select(i => options.SequenceWeights[i]).ToArray();
        }
        else if (options.TCoffeeSequenceWeights)
        {
            weights = TcTCoffeeSequenceWeights(TcDistanceMatrix(seqs, kind, options, matrixGroups).Sim, k);
        }

        // CL->gop / CL->gep: -gapopen / -gapext, or -cosmetic_penalty as the opening penalty when neither is given
        int gop = options.GapOpen, gep = options.GapExtend;
        if (gop == 0 && gep == 0) gop = options.CosmeticPenalty;
        int relax = options.RelaxLibIterations ?? (options.RelaxLibrary ? 1 : 0);
        var (sortedRows, _) = TcAlignDetailed(seqs, index, tree, relax, options.FilterLib, weights, new TcPenalties(gop, gep, options.NoMatch), options.DoNormalise);

        // -outorder aligned writes the rows in the program's internal order (the names sorted, see Names); the default (input) restores the input order
        if (options.OutOrderAligned) return sortedRows;
        var rows = new string[k];
        for (int i = 0; i < k; i++) rows[perm[i]] = sortedRows[i];
        return rows;
    }
}
