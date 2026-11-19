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
    private static int TcSeqSim(string r1, string r2, bool singleLetterGroups = false)
    {
        int pos0 = 0, sim = 0;
        for (int a = 0; a < r1.Length; a++)
        {
            if (r1[a] == '-' || r2[a] == '-') continue;
            pos0++;
            // singleLetterGroups: the groups of make_group_aa("idmat") (26 groups of one letter: J, O, U and X are equal to themselves, T and U differ) instead of those of the matrix idmat
            if (singleLetterGroups ? char.ToUpperInvariant(r1[a]) == char.ToUpperInvariant(r2[a]) : TcSameGroup(r1[a], r2[a])) sim++;
        }

        return pos0 == 0 ? 0 : sim * 100 / pos0;
    }

    /// <summary><c>fix_seq_seq</c> for one sequence: index[r] = coordinate (1-based) in <paramref name="full"/> of the r-th (1-based) residue of <paramref name="segment"/>.</summary>
    private static int[] TcFixIndex(string segment, string full) => TcFixIndex(segment, full, out _);

    /// <summary>As <see cref="TcFixIndex(string, string)"/>; <paramref name="aligned"/> tells that the sequences differ and were reconciled by an alignment (whose <c>aln2sim(B, "idmat")</c> is a call of <c>is_in_same_group_aa</c>).</summary>
    private static int[] TcFixIndex(string segment, string full, out bool aligned)
    {
        aligned = false;
        var index = new int[segment.Length + 1];
        if (string.Equals(segment, full, StringComparison.OrdinalIgnoreCase))
        {
            for (int r = 1; r <= segment.Length; r++) index[r] = r;
            return index;
        }

        // align_two_sequences(..., "idmat" | "blosum62mt", -4, -1): the matrix follows the type of the sequences built by aln2seq, which is never set -> idmat
        aligned = true;
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

    /// <summary><c>is_in_set(letter, "subset_objOBJ-")</c> of <c>aln2constraint_list</c> for the weight mode "sim_subset_objOBJ-": the program takes the text after the FIRST '_' (so "subset_" belongs to the set), tested case-sensitively.</summary>
    private static bool TcInExonSubset(char c) => "subset_objOBJ-".Contains(c);

    /// <summary>The entries one pair job (x, y) adds, in the order of <c>aln2constraint_list</c> (alignments of the chain, columns left to right).</summary>
    private static List<(int S1, int R1, int S2, int R2, int W)> TcPairJobEntries(
        IReadOnlyList<string> seqs, int x, int y, TCoffeePairMethod method, int lalignTop, TcSeqKind kind, string? matrix, TCoffeeMethodParameters? mp, string[]? work = null, string?[]? templates = null, TcProcessGroups? groups = null, TcStructure? structure = null)
    {
        bool protein = kind == TcSeqKind.Protein;
        if (method is TCoffeePairMethod.PromoPair or TCoffeePairMethod.BlastrPair)
            return TcExtendedPairEntries(work ?? seqs.ToArray(), x, y, method, mp);
        if (method == TCoffeePairMethod.ProbaPair)
            return TcProbaPairEntries(seqs[x], seqs[y], kind, matrix).Select(e => (x, e.R1, y, e.R2, e.Weight)).ToList();
        bool structureMethod = method is TCoffeePairMethod.TMalignPair or TCoffeePairMethod.MustangPair or TCoffeePairMethod.SapPair or TCoffeePairMethod.DaliPair;
        var structureTool = method == TCoffeePairMethod.SapPair ? TCoffeePairMethod.TMalignPair : method; // SAP is not installed: the program runs TMalign
        (string, string)? structureRows = null;
        if (structureMethod)
        {
            // SEQ_TYPE P: both sequences need a usable _P_ template, otherwise the job is proba_pair's ("Method TMalign cannot be applied to [x vs y], proba_pair will be used instead")
            if (structure is null || !structure.Usable[x] || !structure.Usable[y])
                return TcProbaPairEntries(seqs[x], seqs[y], kind, matrix).Select(e => (x, e.R1, y, e.R2, e.Weight)).ToList();
            // dali_pair: pdb_pair with EXECUTABLE2 DaliLite, which is not installed: the first pair of two usable structures stops the program ("ERROR: UNSPECIFIED UNSPECIFIED", exit status 1)
            if (method == TCoffeePairMethod.DaliPair) throw new TCoffeeAbortedException("ERROR: UNSPECIFIED UNSPECIFIED -- dali_pair needs DaliLite (EXECUTABLE2 DaliLite), which is not installed");
            string? output = structure.Output(structureTool, x, y);
            if (output is not null) structureRows = structureTool == TCoffeePairMethod.TMalignPair ? TCoffeeTMalignAlignment(output) : TCoffeeMustangAlignment(output);
            if (structureRows is null) return new List<(int, int, int, int, int)>(); // the structure aligner failed: no entry
        }

        var rows = new List<(string, string)>();
        bool subset = false; // exon2 / exon3: only the columns with both letters in "subset_objOBJ-"
        switch (method)
        {
            case TCoffeePairMethod.SlowPair:
                rows.Add(TCoffeeSlowPairAlign(seqs[x], seqs[y], protein, matrix, mp));
                break;
            case TCoffeePairMethod.LalignIdPair:
            case TCoffeePairMethod.LalignRsSPair:
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
            case TCoffeePairMethod.BiphasicPair:
                rows.Add(TCoffeeBiphasicPairAlign(seqs[x], seqs[y], protein, matrix, mp));
                break;
            case TCoffeePairMethod.ExonPair:
            case TCoffeePairMethod.Exon2Pair:
            case TCoffeePairMethod.Exon3Pair:
                {
                    // SEQ_TYPE G: the pair job aligns the genomic templates (the exon methods have their own branch in method2pw_cl: dp mode myers_miller_pair_wise, the x rows of the matrix zeroed);
                    // MATRIX exon2mt, GOP -10 (exon, exon2) / 0 (exon3), GEP -1
                    var matrix2 = (TcResolveMatrix(mp?.Matrix) ?? TcMatrix.ByName("exon2mt")!).WithZeroX();
                    rows.Add(TcMyersMiller(templates![x]!, templates[y]!, matrix2, mp?.GapOpen ?? (method == TCoffeePairMethod.Exon3Pair ? 0 : -10), mp?.GapExtend ?? -1));
                    subset = method != TCoffeePairMethod.ExonPair;
                    break;
                }

            case TCoffeePairMethod.FastPair:
                {
                    var p = TcFastPairParams(protein, matrix, mp);
                    rows.Add(TcFastaPairRows(seqs[x], seqs[y], p.Matrix, p.Gop, p.Gep, protein ? 2 : 5, 0, 0, 0f, 1));
                    break;
                }

            case TCoffeePairMethod.ClustalwPair:
                {
                    // the external program decides the residue type of the pair file itself (-TYPE=auto: the SECOND sequence of the pair), not T-Coffee's type of the whole set:
                    // the pair YAGR / NSMYC consists of nucleotide letters only for T-Coffee (AGCTUNRYMKSWHBVD) but is protein for clustalw
                    var pair = new[] { seqs[x], seqs[y] };
                    var cw = GuessClustalWSequenceType(pair) == ClustalWSequenceType.Dna
                        ? WithClustalWLenientLetters(() => MultipleAlignClustalW(pair))
                        : MultipleAlignClustalW(pair, new ClustalWProteinOptions());
                    rows.Add((cw.AlignedSequences[0], cw.AlignedSequences[1]));
                    break;
                }

            case TCoffeePairMethod.TMalignPair:
            case TCoffeePairMethod.MustangPair:
            case TCoffeePairMethod.SapPair:
                rows.Add(structureRows!.Value);
                break;

            case TCoffeePairMethod.MafftPair:
                throw new NotSupportedException("mafft_pair runs MAFFT 7.505 on every pair; there is no MAFFT port in this repository yet (B09.md A2-19). 460 oracle sets of the program are committed (feature mafft_pair).");
            case TCoffeePairMethod.MusclePair:
                {
                    var mu = MultipleAlignMuscle(new[] { seqs[x], seqs[y] }, new MuscleOptions { MaxIterations = 8 });
                    rows.Add((mu.AlignedSequences[0], mu.AlignedSequences[1]));
                    break;
                }

            default:
                throw new NotSupportedException($"The {method} library is not available here.");
        }

        var entries = new List<(int, int, int, int, int)>();
        groups ??= new TcProcessGroups();
        foreach (var (rowA, rowB) in rows)
        {
            // aln2constraint_list_generic: fix_aln_seq_new (fix_seq_seq: a reconciling alignment calls aln2sim(B, "idmat") = is_in_same_group_aa with the mode "idmat") comes BEFORE the weights
            int[] fa = TcFixIndex(rowA.Replace("-", ""), seqs[x], out bool alignedA);
            int[] fb = TcFixIndex(rowB.Replace("-", ""), seqs[y], out bool alignedB);
            if (groups.SingleLetters == null) groups.SingleLetters = alignedA || alignedB; // the process-static groups are built by the FIRST call
            int weight = TcSeqSim(rowA, rowB, groups.SingleLetters.Value);
            if (weight == 0) weight = 1;
            weight *= 1000 / TcMaxId;
            int na = 0, nb = 0;
            for (int c = 0; c < rowA.Length; c++)
            {
                bool ga = rowA[c] == '-', gb = rowB[c] == '-';
                if (!ga) na++;
                if (!gb) nb++;
                if (!ga && !gb && (!subset || (TcInExonSubset(rowA[c]) && TcInExonSubset(rowB[c])))) entries.Add((x, fa[na], y, fb[nb], weight));
            }
        }

        return entries;
    }

    /// <summary>
    /// What the child process of one METHOD returns (<c>fork_read_n_constraint_list</c> child): its start list (the library merged so far) is cloned, the pair
    /// jobs run in a grand-child whose tail is merged back, and the whole list is returned when it grew.
    /// </summary>
    /// <summary>The residue groups of <c>is_in_same_group_aa</c> in the process that runs the pair jobs of one method: built once, by the first call (<c>null</c> = not built yet; true = the 26 single-letter groups of mode "idmat", false = the groups of the matrix idmat).</summary>
    private sealed class TcProcessGroups
    {
        public bool? SingleLetters;
    }

    private static List<TcEntry>? TcRunMethod(
        TcConstraintList start, IReadOnlyList<string> seqs, List<(int X, int Y)> pairs, TCoffeePairMethod method, int lalignTop, TcSeqKind kind, string? matrix, TCoffeeMethodParameters? mp, string?[]? templates = null, TcStructure? structure = null, IEnumerable<(int S1, int R1, int S2, int R2, int W)>? fixedEntries = null)
    {
        var methodList = start.Clone();
        int before = methodList.EntryCount;
        var child = methodList.Clone();
        var frozen = child.Freeze();
        // the jobs of one method run one after the other in one process on the program's shared sequence array (promo_pair / blastr_pair modify it)
        var work = seqs.ToArray();
        var groups = new TcProcessGroups();
        foreach (var (x, y) in pairs)
            foreach (var (s1, r1, s2, r2, w) in TcPairJobEntries(seqs, x, y, method, lalignTop, kind, matrix, mp, work, templates, groups, structure))
                child.AddBoth(s1, r1, s2, r2, w, 1);
        // a library file read by read_constraint_list_file (the single job of a multiple method): CONS 0
        if (fixedEntries != null)
            foreach (var (s1, r1, s2, r2, w) in fixedEntries)
                child.AddBoth(s1, r1, s2, r2, w, 0);
        methodList.Undump(child.Dump(frozen));
        return methodList.EntryCount > before ? methodList.Dump(null) : null;
    }

    /// <summary>
    /// The program's constraint list after the pair-wise methods <paramref name="methods"/> were read in this order
    /// (see the file header for the process structure that is reproduced). A method child starts from the list merged so far.
    /// </summary>
    internal static TcConstraintList TcBuildConstraintList(
        IReadOnlyList<string> seqs, IReadOnlyList<TCoffeePairMethod> methods, int lalignTop, TcSeqKind kind, string? matrix = null, IReadOnlyDictionary<TCoffeePairMethod, TCoffeeMethodParameters>? methodParameters = null, string?[]? templates = null, TcStructure? structure = null, TcBlastInput? blast = null)
    {
        var lengths = seqs.Select(s => s.Length).ToArray();
        var main = new TcConstraintList(lengths);
        var pairs = new List<(int X, int Y)>();
        var pairsWithSelf = new List<(int X, int Y)>(); // ALN_MODE s_pairwise: row-major (x, y >= x), the self pair first in every row
        for (int x = 0; x < seqs.Count; x++)
        {
            pairsWithSelf.Add((x, x));
            for (int y = x + 1; y < seqs.Count; y++)
            {
                pairs.Add((x, y));
                pairsWithSelf.Add((x, y));
            }
        }

        // a method named twice is read once
        foreach (var method in methods.Distinct())
        {
            if (pairs.Count > 0) TcRefuseMethod(method);
            var methodPairs = method == TCoffeePairMethod.LalignRsSPair ? pairsWithSelf : pairs;
            if (method is TCoffeePairMethod.BlastpMsa or TCoffeePairMethod.TblastxMsa)
            {
                // ALN_MODE multiple: ONE job for all the sequences (seq_msa / tc_generic_method.pl); the library comes from the HSPs of the search
                if (blast is null) throw new ArgumentException($"{(method == TCoffeePairMethod.BlastpMsa ? "blastp_msa" : "tblastx_msa")} needs the XML of the BLAST search (TCoffeeProgramOptions.BlastXml): the search itself is the external program.");
                var blastEntries = method == TCoffeePairMethod.BlastpMsa ? TcBlastMsaEntries(blast.Names, seqs, blast.Xml) : TcTblastxEntries(seqs, blast.Xml);
                var blastDump = TcRunMethod(main, seqs, new List<(int X, int Y)>(), method, lalignTop, kind, matrix, null, templates, structure, blastEntries);
                if (blastDump != null) main.Undump(blastDump);
                continue;
            }

            if (method is TCoffeePairMethod.ExonPair or TCoffeePairMethod.Exon2Pair or TCoffeePairMethod.Exon3Pair)
                methodPairs = pairs.Where(p => templates != null && templates[p.X] != null && templates[p.Y] != null).ToList(); // check_seq_type: SEQ_TYPE G needs a _G_ template on both sequences
            var dump = TcRunMethod(main, seqs, methodPairs, method, lalignTop, kind, matrix, methodParameters != null && methodParameters.TryGetValue(method, out var mp) ? mp : null, templates, structure);
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
    /// <param name="genomicTemplates">The <c>_G_</c> template of every sequence for <see cref="TCoffeePairMethod.ExonPair"/> / <see cref="TCoffeePairMethod.Exon2Pair"/> / <see cref="TCoffeePairMethod.Exon3Pair"/> (<c>null</c> entries = no template).</param>
    public static IReadOnlyList<TCoffeeLibraryEntry> TCoffeePrimaryLibraryEntries(
        IReadOnlyList<string> sequences, IEnumerable<TCoffeePairMethod> methods, int lalignTopAlignments = 10, bool? protein = null, IReadOnlyList<string?>? genomicTemplates = null)
    {
        var cl = TcBuildChecked(sequences, methods, lalignTopAlignments, TcKindOf(sequences, protein), templates: genomicTemplates?.ToArray());
        var list = new List<TCoffeeLibraryEntry>();
        for (int s1 = 0; s1 < sequences.Count; s1++)
            for (int r1 = 1; r1 <= cl.Length(s1); r1++)
                foreach (var e in cl.Entries(s1, r1))
                    if (s1 < e.Seq2 && e.Res2 > 0)
                        list.Add(new TCoffeeLibraryEntry(s1, r1 - 1, e.Seq2, e.Res2 - 1, e.Weight, e.Count));
        return list.OrderBy(e => e.SeqA).ThenBy(e => e.SeqB).ThenBy(e => e.PosA).ThenBy(e => e.PosB).ToList();
    }

    private static TcConstraintList TcBuildChecked(IReadOnlyList<string> sequences, IEnumerable<TCoffeePairMethod> methods, int lalignTop, TcSeqKind kind, string? matrix = null, IReadOnlyDictionary<TCoffeePairMethod, TCoffeeMethodParameters>? methodParameters = null, string?[]? templates = null, TcStructure? structure = null, TcBlastInput? blast = null)
    {
        ArgumentNullException.ThrowIfNull(sequences);
        ArgumentNullException.ThrowIfNull(methods);
        foreach (var s in sequences)
            if (string.IsNullOrEmpty(s)) throw new ArgumentException("Empty sequences are not supported.", nameof(sequences));
        return TcBuildConstraintList(sequences, methods.ToList(), lalignTop, kind, matrix, methodParameters, templates, structure, blast);
    }
}

/// <summary><c>-tree_mode</c> of the program.</summary>
public enum TCoffeeTreeMode
{
    /// <summary>Neighbour joining (ClustalW 1.x routine) with the trichotomy at the root turned into a binary tree; default.</summary>
    Nj = 0,

    /// <summary>UPGMA (<c>upgma_merge</c>) on the integer distances.</summary>
    Upgma = 1,

    /// <summary><c>cwph</c> (listed in the option's documentation): <c>make_tree</c> has no such mode; for three or more sequences the program core dumps (exit status 11); two sequences use a dummy tree and ignore the mode. Raises <see cref="TCoffeeAbortedException"/>.</summary>
    Cwph = 2,

    /// <summary><c>kmeans</c>: <c>seq2km_tree_old</c> builds a bootstrapped k-means tree ("Resolve Node" lines) and the program ends with exit status 0 WITHOUT writing an alignment (three or more sequences). Raises <see cref="TCoffeeAbortedException"/> with exit status 0.</summary>
    Kmeans = 3,

    /// <summary>Any other text: <c>make_tree</c> reports "unknown tree computation mode" and the program core dumps (exit status 11) for three or more sequences; the text is <see cref="TCoffeeProgramOptions.TreeModeText"/>.</summary>
    Other = 4,
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

    /// <summary><c>aln</c>: <c>aln2sim_mat</c> of the alignment that <c>cl2distance_matrix</c> never gets (<c>NOALN</c>): the program segfaults (exit status 11). Raises <see cref="TCoffeeAbortedException"/>.</summary>
    Aln = 4,

    /// <summary><c>ktup2</c>: <c>seq2distance_matrix</c> aligns with <c>very_fast_aln</c> and calls itself with the mode <c>aln</c>; the program stops with "ERROR: UNSPECIFIED UNSPECIFIED" (exit status 1). Raises <see cref="TCoffeeAbortedException"/>.</summary>
    Ktup2 = 5,

    /// <summary><c>slow</c>: every pair aligned with the relaxed library (<c>pair_wise</c> with the main dp mode): the identity of the alignment is the similarity, the DP score per column (rescaled by <c>-do_normalise</c>) the distance score.</summary>
    Slow = 6,

    /// <summary><c>cscore</c>: the consistency score of the pair alignment (<c>linked_pair_wise</c> score per column and per SCORE_K, float); falls back to <c>idscore</c> when the library is empty.</summary>
    Cscore = 7,

    /// <summary><c>geometric</c>: Euclidean distance of the vectors of FASTA-alignment (<c>fasta_pair_wise</c>, pam250mt) identities against the first three sequences.</summary>
    Geometric = 8,

    /// <summary><c>random</c>: <c>rand() % 100</c> as identity; the program's <c>rand</c> is seeded differently on every run (two runs of the same input give different weights), so the result cannot be reproduced; <see cref="NotSupportedException"/>.</summary>
    Random = 9,

    /// <summary><c>clustalw</c>: the DP score of the linked pair-wise alignment of the pair with the relaxed library as identity (the program has no mode of this name: it falls into the generic branch of <c>seq2distance_matrix</c>).</summary>
    Clustalw = 10,

    /// <summary><c>very_fast</c> (also <c>-quicktree</c>): <c>very_fast_gotoh_pair_wise</c> of every pair: the diagonals scored by the relaxed library entries, the 20 best with a shadow of 5, banded Gotoh DP, identity of the aligned rows.</summary>
    VeryFast = 12,

    /// <summary>Any other text: "Error: X is an unknown distance_matrix_mode [FATAL:T-COFFEE]", exit status 1 (<see cref="TCoffeeProgramOptions.DistanceModeText"/>). Raises <see cref="TCoffeeAbortedException"/>.</summary>
    Other = 11,
}

/// <summary>Options of <see cref="SequenceAligner.MultipleAlignTCoffeeProgram"/> = the command line of the program they mirror.</summary>
public sealed partial record TCoffeeProgramOptions
{
    /// <summary>
    /// Cancels the alignment (the same convention as <see cref="System.Threading.Tasks.ParallelOptions.CancellationToken"/>): the dynamic programs poll it at their
    /// row / iteration boundaries and the call throws <see cref="OperationCanceledException"/> without a partial result. Default: never cancelled.
    /// It does not take part in the semantics of the options (two option sets that differ only in the token align identically).
    /// </summary>
    public CancellationToken CancellationToken { get; init; }

    /// <summary><c>-method m1,m2</c>: the primary library methods in command-line order. Default = <c>proba_pair</c>, the program's own default of 13.41 (<c>t_coffee seqs.fa</c>).</summary>
    public IReadOnlyList<TCoffeePairMethod> Methods { get; init; } = new[] { TCoffeePairMethod.ProbaPair };

    /// <summary>
    /// <c>-template_file</c> with <c>_G_</c> templates: the genomic template sequence of every input sequence in input order (<c>null</c> = the sequence has none), the input of <see cref="TCoffeePairMethod.ExonPair"/>,
    /// <see cref="TCoffeePairMethod.Exon2Pair"/> and <see cref="TCoffeePairMethod.Exon3Pair"/>. The command-line form is <c>-template_file @t0,t1,...</c>.
    /// </summary>
    public IReadOnlyList<string?>? GenomicTemplates { get; init; }

    /// <summary>
    /// <c>-template_file</c> with <c>_P_</c> templates (A2-10): the sequence of the PDB structure of every input sequence in input order (<c>null</c> / empty = no structure), the input of
    /// <see cref="TCoffeePairMethod.TMalignPair"/> and <see cref="TCoffeePairMethod.MustangPair"/>. A structure is attached only when its sequence is similar enough to the sequence (identity &gt; 35 %, coverage &gt; 50 %);
    /// pairs without two structures are aligned by <c>proba_pair</c> as the program does. What the external structure aligner wrote for the pairs is <see cref="StructureAlignerOutputs"/>.
    /// </summary>
    public IReadOnlyList<string?>? StructureTemplates { get; init; }

    /// <summary>
    /// What the external BLAST search of <see cref="TCoffeePairMethod.BlastpMsa"/> wrote (A2-9): the XML (<c>blastall -p blastp -m7</c>, NCBI BLAST+ <c>-outfmt 5</c>) of the search of every sequence against a local database of all of them, one
    /// <c>Iteration</c> per sequence in the program's order. The hits are matched to the sequences by <see cref="Names"/> (omitted: q000, q001, ... by input index). The BLAST search is not reimplemented.
    /// </summary>
    public string? BlastXml { get; init; }

    /// <summary>
    /// What the external structure aligner of each method wrote for a pair of structures (the structure aligner is the external program; this library is T-Coffee's side): per method, keyed by the INPUT indices of the two sequences in the
    /// order the program calls the aligner (names sorted unless <see cref="InOrderInput"/>; first = the sequence that comes first), the standard output of <c>TMalign pdb1 pdb2</c> (its last four lines are what matters,
    /// see <see cref="SequenceAligner.TCoffeeTMalignAlignment"/>) or the content of the <c>results.afasta</c> of <c>mustang -i pdb1 pdb2 -F fasta</c> (<see cref="SequenceAligner.TCoffeeMustangAlignment"/>). A missing pair = the aligner failed.
    /// </summary>
    public IReadOnlyDictionary<TCoffeePairMethod, IReadOnlyDictionary<(int First, int Second), string>>? StructureAlignerOutputs { get; init; }

    /// <summary><c>-dp_mode</c>: the dynamic program of the progressive stage (default the sparse <c>linked_pair_wise</c>).</summary>
    public TCoffeeDpMode DpMode { get; init; } = TCoffeeDpMode.Linked;

    /// <summary>The text of <c>-dp_mode</c> when <see cref="DpMode"/> is <see cref="TCoffeeDpMode.Other"/>.</summary>
    public string? DpModeText { get; init; }

    /// <summary><c>-extend_mode</c> (default <c>very_fast_triplet</c>): how the library scores a pair of profile columns in the generic dynamic programs.</summary>
    public TCoffeeExtendMode ExtendMode { get; init; } = TCoffeeExtendMode.VeryFastTriplet;

    /// <summary>The text of <c>-extend_mode</c> when <see cref="ExtendMode"/> is <see cref="TCoffeeExtendMode.Other"/>.</summary>
    public string? ExtendModeText { get; init; }

    /// <summary><c>-max_n_pair</c> (default 10): the number of pairs of sequences compared for two profile columns by <c>very_fast_triplet</c> (0 = every pair).</summary>
    public int MaxNPair { get; init; } = 10;

    /// <summary><c>-diag_threshold</c> (default 0 = the maximum z-score / 25): the z-score above which a diagonal is used by <c>-dp_mode fasta_pair_wise</c>.</summary>
    public int DiagonalThreshold { get; init; }

    /// <summary><c>-weight no</c>: the sequence weighting stage (and with it the distance matrix of two sequences) is skipped.</summary>
    public bool WeightNo { get; init; }

    /// <summary><c>-tg_mode</c> (default 1): the terminal gap rule of the heuristic dynamic programs (0 = opening and extension, 1 = extension only, 2 = free); the sparse DP of <c>linked_pair_wise</c> has its own rule.</summary>
    public int TgMode { get; init; } = 1;

    /// <summary><c>-ndiag</c> (<c>fasta_step</c>, default 0): the number of best diagonals the FASTA heuristics use instead of the z-score threshold (<see cref="TCoffeePairMethod.FastPair"/>, <c>-dp_mode fasta_pair_wise</c>).</summary>
    public int DiagonalCount { get; init; }

    /// <summary><c>-ktuple</c> (default 1, DNA at least 2): the k-tuple size of the heuristic dynamic programming modes (<c>-dp_mode fasta_pair_wise</c>, <c>cfasta_pair_wise</c>).</summary>
    public int Ktuple { get; init; } = 1;

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

    /// <summary>The text of <c>-tree_mode</c> when <see cref="TreeMode"/> is <see cref="TCoffeeTreeMode.Other"/>.</summary>
    public string? TreeModeText { get; init; }

    /// <summary>The text of <c>-distance_matrix_mode</c> when <see cref="DistanceMatrixMode"/> is <see cref="TCoffeeDistanceMode.Other"/>.</summary>
    public string? DistanceModeText { get; init; }

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
        var list = ToListNoNull(sequences, nameof(sequences));
        if (list.Count == 0) return MultipleAlignmentResult.Empty;
        using var cancelScope = AlignmentCancellation.Enter(options?.CancellationToken ?? default);
        var rows = TCoffeeProgramRows(list, options);
        // ktup_pair_wise / precomputed_pair_wise / gotoh_pair_wise_sw write rows of different lengths (no alignment of the groups): there is no column to score or to take a consensus of
        if (rows.Select(r => r.Length).Distinct().Count() > 1) return new MultipleAlignmentResult(rows, string.Empty, 0);
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
        RequireProgramMembers(sequences, "t_coffee", nameof(sequences));
        foreach (var s in sequences)
            if (s.Any(c => !char.IsAsciiLetter(c))) throw new ArgumentException("Sequences must be strings of letters.", nameof(sequences));
        TcRefuseModes(options, sequences.Count);
        if (sequences.Count == 1) return new[] { sequences[0] };
        int k = sequences.Count;
        // -extend_mode rna0 = residue_pair_extended_list + slow_get_dp_cost (the same column score as slow_triplet)
        // rna1-4 / rna read the contact library `-rna_lib` (none given: R is empty and residue_pair_extended_list4rna returns the triplet score of the pair itself); rna4 would run RNAplfold (not installed): the
        // column score is then that of rna0
        if (options.ExtendMode == TCoffeeExtendMode.Other && options.ExtendModeText is "rna0" or "rna" or "rna1" or "rna2" or "rna3" or "rna4") options = options with { ExtendMode = TCoffeeExtendMode.SlowTriplet };

        // -inorder aligned (the program's default): the sequences are sorted by name (byte order) before anything else is computed
        var perm = Enumerable.Range(0, k).ToArray();
        if (!options.InOrderInput && options.Names != null)
        {
            if (options.Names.Count != k) throw new ArgumentException("Names must have one entry per sequence.", nameof(options));
            perm = perm.OrderBy(i => options.Names[i], StringComparer.Ordinal).ToArray();
        }

        var seqs = perm.Select(i => sequences[i]).ToList();
        var kind = TcKindOf(seqs, options.Protein);
        if (options.GenomicTemplates != null && options.GenomicTemplates.Count != k) throw new ArgumentException("GenomicTemplates must have one entry per sequence.", nameof(options));
        string?[]? templates = options.GenomicTemplates == null ? null : perm.Select(i => options.GenomicTemplates[i]).ToArray();
        TcStructure? structure = null;
        if (options.StructureTemplates != null)
        {
            if (options.StructureTemplates.Count != k) throw new ArgumentException("StructureTemplates must have one entry per sequence.", nameof(options));
            var outputs = options.StructureAlignerOutputs;
            structure = new TcStructure
            {
                Usable = perm.Select(i => !string.IsNullOrEmpty(options.StructureTemplates[i]) && TcPdbTemplateUsable(sequences[i], options.StructureTemplates[i]!)).ToArray(),
                Output = (method, x, y) => outputs != null && outputs.TryGetValue(method, out var byPair) && byPair.TryGetValue((perm[x], perm[y]), out var text) ? text : null,
            };
        }

        TcBlastInput? blast = null;
        if (options.BlastXml != null)
            blast = new TcBlastInput { Names = perm.Select(i => options.Names != null ? options.Names[i] : $"q{i:D3}").ToArray(), Xml = options.BlastXml };
        var cl = TcBuildChecked(seqs, options.Methods, options.LalignTopAlignments, kind, options.Matrix, options.MethodParameters, templates, structure, blast);
        // CL->S->nseq > 1 with an empty library: "You have not provided any method or enough Sequences [FATAL]" (t_coffee.c L4999-5010)
        if (cl.EntryCount == 0)
            throw new TCoffeeAbortedException("You have not provided any method or enough Sequences [FATAL] (no library entry: a method without suitable sequences produced no pair)");
        var index = TcIndex.FromConstraintList(cl, seqs.Select(s => s.Length).ToArray());

        float[]? explicitWeights = null;
        if (options.SequenceWeights != null)
        {
            if (options.SequenceWeights.Count != k) throw new ArgumentException("SequenceWeights must have one entry per sequence.", nameof(options));
            explicitWeights = perm.Select(i => options.SequenceWeights[i]).ToArray();
        }

        // CL->gop / CL->gep: -gapopen / -gapext, or -cosmetic_penalty as the opening penalty when neither is given
        int gop = options.GapOpen, gep = options.GapExtend;
        if (gop == 0 && gep == 0) gop = options.CosmeticPenalty;
        int relax = options.RelaxLibIterations ?? (options.RelaxLibrary ? 1 : 0);

        // the distance matrix of the program is computed after the library was relaxed (modes slow / cscore / clustalw / very_fast read it): guide tree and sequence weights come from that stage
        var dpSettings = new TcDpSettings(options.DpMode, options.ExtendMode, gop, gep, options.TgMode, kind == TcSeqKind.Protein ? options.Ktuple : Math.Max(2, options.Ktuple), options.DiagonalCount, options.DiagonalThreshold, options.MaxNPair, options.DoNormalise);
        (TCoffeeTreeNode, float[]?, int[,]?) LateStage(TcIndex relaxed, TcPenalties pen)
        {
            var state = new TcRelaxedState(relaxed, pen, options.DoNormalise, options.TgMode, dpSettings);
            var tree = options.DistanceMatrixMode == TCoffeeDistanceMode.Ktup && options.TreeMode == TCoffeeTreeMode.Nj
                ? TCoffeeGuideTree(seqs, kind == TcSeqKind.Protein)
                : TcGuideTreeFor(seqs, kind, options, state);

            // the program computes the distance matrix of two sequences as well (cl2distance_matrix in t_coffee.c); only its side effect matters here: the segmentation fault of 'fast' + an idscore sim mode
            if (k == 2 && !options.WeightNo) _ = TcDistanceMatrix(seqs, kind, options, state);

            // the distance matrix is the program's CL->DM: the generic dynamic programs sort the pairs of sequences of two profile columns by its identities
            int[,]? sim = null;
            float[]? weights = explicitWeights;
            if (options.DpMode != TCoffeeDpMode.Linked || (weights == null && options.TCoffeeSequenceWeights))
            {
                var dm = TcDistanceMatrix(seqs, kind, options, state);
                sim = dm.Sim;
                if (weights == null && options.TCoffeeSequenceWeights) weights = TcTCoffeeSequenceWeights(dm.Sim, k);
            }

            return (tree, weights, sim);
        }

        var (sortedRows, _) = TcAlignDetailed(seqs, index, null, relax, options.FilterLib, null, new TcPenalties(gop, gep, options.NoMatch), options.DoNormalise, LateStage, dpSettings);

        // -outorder aligned writes the rows in the program's internal order (the names sorted, see Names); the default (input) restores the input order
        if (options.OutOrderAligned) return sortedRows;
        var rows = new string[k];
        for (int i = 0; i < k; i++) rows[perm[i]] = sortedRows[i];
        return rows;
    }
}
