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
        IReadOnlyList<string> seqs, int x, int y, TCoffeePairMethod method, int lalignTop, bool protein)
    {
        if (method == TCoffeePairMethod.ProbaPair)
            return TcProbaPairEntries(seqs[x], seqs[y], protein).Select(e => (x, e.R1, y, e.R2, e.Weight)).ToList();

        var rows = new List<(string, string)>();
        switch (method)
        {
            case TCoffeePairMethod.SlowPair:
                rows.Add(TCoffeeSlowPairAlign(seqs[x], seqs[y], protein));
                break;
            case TCoffeePairMethod.LalignIdPair:
                {
                    var alns = TCoffeeLalignTopAlignments(seqs[x], seqs[y], lalignTop, protein);
                    if (alns.Count == 0)
                    {
                        // sim_pair_wise_lalign found nothing: the program then turns the padded (unaligned) pair into the library
                        int len = Math.Max(seqs[x].Length, seqs[y].Length);
                        rows.Add((seqs[x].PadRight(len, '-'), seqs[y].PadRight(len, '-')));
                    }

                    foreach (var a in alns) rows.Add((a.RowA, a.RowB));
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
        TcConstraintList start, IReadOnlyList<string> seqs, List<(int X, int Y)> pairs, TCoffeePairMethod method, int lalignTop, bool protein)
    {
        var methodList = start.Clone();
        int before = methodList.EntryCount;
        var child = methodList.Clone();
        var frozen = child.Freeze();
        foreach (var (x, y) in pairs)
            foreach (var (s1, r1, s2, r2, w) in TcPairJobEntries(seqs, x, y, method, lalignTop, protein))
                child.AddBoth(s1, r1, s2, r2, w, 1);
        methodList.Undump(child.Dump(frozen));
        return methodList.EntryCount > before ? methodList.Dump(null) : null;
    }

    /// <summary>
    /// The program's constraint list after the pair-wise methods <paramref name="methods"/> were read in this order
    /// (see the file header for the process structure that is reproduced). A method child starts from the list merged so far.
    /// </summary>
    internal static TcConstraintList TcBuildConstraintList(
        IReadOnlyList<string> seqs, IReadOnlyList<TCoffeePairMethod> methods, int lalignTop, bool protein)
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
            var dump = TcRunMethod(main, seqs, pairs, method, lalignTop, protein);
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
        var cl = TcBuildChecked(sequences, methods, lalignTopAlignments, protein);
        var list = new List<TCoffeeLibraryEntry>();
        for (int s1 = 0; s1 < sequences.Count; s1++)
            for (int r1 = 1; r1 <= cl.Length(s1); r1++)
                foreach (var e in cl.Entries(s1, r1))
                    if (s1 < e.Seq2 && e.Res2 > 0)
                        list.Add(new TCoffeeLibraryEntry(s1, r1 - 1, e.Seq2, e.Res2 - 1, e.Weight, e.Count));
        return list.OrderBy(e => e.SeqA).ThenBy(e => e.SeqB).ThenBy(e => e.PosA).ThenBy(e => e.PosB).ToList();
    }

    private static TcConstraintList TcBuildChecked(IReadOnlyList<string> sequences, IEnumerable<TCoffeePairMethod> methods, int lalignTop, bool? protein)
    {
        ArgumentNullException.ThrowIfNull(sequences);
        ArgumentNullException.ThrowIfNull(methods);
        foreach (var s in sequences)
            if (string.IsNullOrEmpty(s)) throw new ArgumentException("Empty sequences are not supported.", nameof(sequences));
        return TcBuildConstraintList(sequences, methods.ToList(), lalignTop, protein ?? TcIsProtein(sequences));
    }
}

/// <summary>Options of <see cref="SequenceAligner.MultipleAlignTCoffeeProgram"/> = the command line of the program they mirror.</summary>
public sealed record TCoffeeProgramOptions
{
    /// <summary><c>-method m1,m2</c>: the primary library methods in command-line order. Default = <c>proba_pair</c>, the program's own default of 13.41 (<c>t_coffee seqs.fa</c>).</summary>
    public IReadOnlyList<TCoffeePairMethod> Methods { get; init; } = new[] { TCoffeePairMethod.ProbaPair };

    /// <summary><c>-lalign_n_top</c>: alignments reported by <c>lalign_id_pair</c> per pair (default 10).</summary>
    public int LalignTopAlignments { get; init; } = 10;

    /// <summary><c>-relax_lib 1</c> (default) or <c>0</c>.</summary>
    public bool RelaxLibrary { get; init; } = true;

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
    /// Not available: RNA input (a <c>U</c>) for <see cref="TCoffeePairMethod.ProbaPair"/> (the program uses a different model and <c>get_tot_prob</c> for it).
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
        bool protein = options.Protein ?? TcIsProtein(sequences);
        if (options.Methods.Contains(TCoffeePairMethod.ProbaPair) && !protein && sequences.Any(s => s.Contains('U') || s.Contains('u')))
            throw new NotSupportedException("proba_pair of RNA (a sequence containing U) uses the program's separate RNA model and is not available; use DNA letters or another method.");
        var cl = TcBuildChecked(sequences, options.Methods, options.LalignTopAlignments, protein);
        var index = TcIndex.FromConstraintList(cl, sequences.Select(s => s.Length).ToArray());
        return TcAlignWithIndex(sequences, index, TCoffeeGuideTree(sequences, protein), options.RelaxLibrary);
    }
}
