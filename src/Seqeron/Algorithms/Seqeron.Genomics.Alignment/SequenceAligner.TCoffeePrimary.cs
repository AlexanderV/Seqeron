
namespace Seqeron.Genomics.Alignment;

// =====================================================================================================================
// T-Coffee 13.41 PRIMARY LIBRARIES (audit A1-13b, F66): the pair-wise methods the program runs itself.
//
//   * slow_pair       (util_dp_mm_nw.c `myers_miller_pair_wise` / `diff`): global alignment by the Myers-Miller linear-space
//                     divide and conquer with affine gaps gap(k) = g + h k, g = gop x SCORE_K, h = gep x SCORE_K, terminal gaps pay
//                     h per residue only (the recursion starts with tb = te = 0), midpoint tie rule `(CC != DD && RR == SS)`.
//                     DNA: matrix `dna_idmat` (read_matrice maps the NAME `dna_idmat` to the array `idmat`: +10 diagonal, 0 elsewhere),
//                     GOP -10, GEP -1 (the method file of util_constraints_list.c); protein: blosum62mt, GOP = (int)(10 x average
//                     off-diagonal score of the 20 amino acids) = (int)(get_avg_matrix_mm x 10), GEP -1.
//   * lalign_id_pair  (util_dp_sim.c `sim_pair_wise_lalign`): Huang-Miller SIM, the K = `-lalign_n_top` (10) best NON-INTERSECTING local
//                     alignments, gaps q + r k with q = |gop| x SCORE_K, r = |gep| x SCORE_K (DNA: idmat, GOP -10, GEP -1; protein:
//                     blosum50mt, GOP -10, GEP -4).
//   * library         (util_constraints_list.c `aln2constraint_list`): every column of a pair-wise alignment that has a residue in both
//                     rows becomes the library pair with weight (NORM_F / MAXID) x get_seq_sim = 10 x percent identity of THAT alignment
//                     (integer division, 0 -> 1); identical pairs coming from several alignments / methods keep the maximum weight
//                     (`update_entry`: WE = MAX, CONS += 1).
// All of it is integer arithmetic, so the port is exact by construction; every library of the oracle sets is compared pair by pair
// (weights and CONS counts) with `t_coffee -lib_only -out_lib` (see docs/Evidence/ALIGN-MULTI-001-Evidence.md).
// =====================================================================================================================

/// <summary>Raised when the T-Coffee 13.41 program itself would abort (core dump) on the input.</summary>
/// <remarks>
/// SIM (<c>lalign_id_pair</c>, util_dp_sim.c) keeps the edit script of an alignment in an array <c>S</c> of <c>min(M, N) x 5 / 4 + 1</c> integers. With a small gap open penalty
/// (e.g. <c>-method lalign_id_pair@EP@GOP@-3</c>) an alignment can need more entries: the write runs past the array. Up to the slack of the heap chunk (glibc rounds a request of R bytes up to
/// ((R + 8 + 15) / 16 x 16) - 8 usable bytes) the overflow is harmless; beyond it the program core dumps (<c>COREDUMP</c>, exit status 1, no alignment). Measured on 25 sets with a very
/// small gap open (3 core dumps, all and only the sets with a pair beyond the chunk) and 240 random pairs.
/// </remarks>
public sealed class TCoffeeAbortedException : InvalidOperationException
{
    /// <summary>Creates the exception for a deterministic refusal / core dump of the program that is not SIM's script overflow (an unknown method, an unknown weight mode, a missing dispatcher branch, an empty library).</summary>
    /// <param name="reason">What the program reports (or, for a silent core dump, why it dies).</param>
    /// <param name="exitStatus">The exit status of the program (1 for every case reproduced).</param>
    /// <param name="heapDependent">The program reads memory that is not part of any object (undefined behaviour); whether it then dies depends on the heap contents, so a set may also run to the end (see <see cref="HeapDependent"/>).</param>
    public TCoffeeAbortedException(string reason, int exitStatus = 1, bool heapDependent = false)
        : base($"t_coffee 13.41 stops with exit status {exitStatus}: {reason}" + (heapDependent ? " (heap dependent: the program may also run through)" : string.Empty))
    {
        Reason = reason;
        ExitStatus = exitStatus;
        HeapDependent = heapDependent;
    }

    /// <summary>
    /// <c>false</c> = the program stops on this input every time (a refusal, a missing branch, SIM's script overflow, a wild pointer read in the first calls); <c>true</c> = undefined behaviour whose outcome depends on the
    /// heap contents: the committed oracle sets with such an input crash in most but not all runs (see <c>-distance_matrix_sim_mode idscore</c> of <c>-distance_matrix_mode fast</c>).
    /// </summary>
    public bool HeapDependent { get; }

    /// <summary>What the program reports; empty for the SIM script overflow.</summary>
    public string Reason { get; } = string.Empty;

    /// <summary>The exit status of the program.</summary>
    public int ExitStatus { get; } = 1;

    /// <summary>Creates the exception.</summary>
    /// <param name="scriptLength">Entries the edit script needed.</param>
    /// <param name="usableEntries">Entries the heap chunk of the script array can hold.</param>
    public TCoffeeAbortedException(int scriptLength, int usableEntries)
        : base($"t_coffee 13.41 core dumps: SIM's edit script needs {scriptLength} entries but its array holds only {usableEntries} (heap overflow); the program writes no alignment.")
    {
        ScriptLength = scriptLength;
        UsableEntries = usableEntries;
    }

    /// <summary>Entries the edit script needed.</summary>
    public int ScriptLength { get; }

    /// <summary>Entries the heap chunk of the script array can hold.</summary>
    public int UsableEntries { get; }
}

/// <summary>The pair-wise methods of the T-Coffee program that build the primary library.</summary>
public enum TCoffeePairMethod
{
    /// <summary><c>slow_pair</c>: Myers-Miller global alignment (the program's DNA default is <c>slow_pair</c> + <c>lalign_id_pair</c> in the paper; see <see cref="ProbaPair"/> for the 13.41 default).</summary>
    SlowPair = 0,

    /// <summary><c>lalign_id_pair</c>: Huang-Miller SIM, the <see cref="TCoffeeProgramOptions.LalignTopAlignments"/> best non-intersecting local alignments.</summary>
    LalignIdPair = 1,

    /// <summary><c>proba_pair</c>: ProbCons-style pair-HMM posterior library (the default method of T-Coffee 13.41).</summary>
    ProbaPair = 2,

    /// <summary>
    /// <c>clustalw_pair</c>: the program runs the external <c>clustalw -INFILE -OUTFILE -OUTORDER=INPUT -align</c> on every pair of sequences and reads the alignment; here the ClustalW 2.1 port of this
    /// repository (<see cref="SequenceAligner.MultipleAlignClustalW(IEnumerable{string}, ClustalWOptions?, ScoringMatrix?)"/> and its PROTEIN counterpart) is the aligner. Weights = 10 x percent identity as for the other methods.
    /// </summary>
    ClustalwPair = 3,

    /// <summary><c>biphasic_pair</c>: global alignment with two affine gap models at once (opening x 2 / extension / 2 and opening / 2 / extension x 2; <c>biphasic_pair_wise</c>, util_dp_suboptimal_nw.c).</summary>
    BiphasicPair = 4,

    /// <summary>
    /// <c>lalign_rs_s_pair</c> ("for Mocca"): <c>lalign_id_pair</c> with <c>ALN_MODE s_pairwise</c> - the pair jobs also contain every sequence against ITSELF. <c>fast_pair</c> always gives the two alignment rows the
    /// indices 0 and 1, so SIM sees two sequences (<c>nseq == 2</c>) and the main diagonal is an ordinary candidate: the library gets entries between residues of one sequence (the identity alignment first).
    /// </summary>
    LalignRsSPair = 5,

    /// <summary>
    /// <c>promo_pair</c>: <c>slow_pair</c> on the dinucleotide-extended alphabet (<c>EXTEND_SEQ 1</c>, matrix promoter_tf1, GOP -60, GEP -1, <c>WEIGHT extend_sim</c>); the shared sequence array is extended and
    /// unextended around every job, which is lossy for letters other than a c g t u n x.
    /// </summary>
    PromoPair = 6,

    /// <summary><c>blastr_pair</c>: as <see cref="PromoPair"/> with the matrix blosumR, GOP -20, GEP 0.</summary>
    BlastrPair = 7,

    /// <summary>
    /// <c>clean_slow_pair</c>: the method file asks for <c>WEIGHT clean</c>, which <c>seqpair2weight</c> does not know: every pair job dies with "Weight Mode clean is unknown [FATAL:T-COFFEE]" and a core dump (exit
    /// status 1) as soon as there are two sequences. Raises <see cref="TCoffeeAbortedException"/>.
    /// </summary>
    CleanSlowPair = 8,

    /// <summary>
    /// <c>hash_pair</c>: <c>method_file2... "EXECUTABLE hash_pair"</c> has no branch in the pair-job dispatcher (<c>THE FUNCTION hash_pair DOES NOT EXIST</c>, written to the suppressed local stderr): core dump, exit status 1,
    /// as soon as there are two sequences. Raises <see cref="TCoffeeAbortedException"/>.
    /// </summary>
    HashPair = 9,

    /// <summary><c>subop1_pair</c>: the dispatcher knows the function but the method list of the program has no such method: "Msubop1_pair not recognised [FATAL:T-COFFEE]", exit status 1. Raises <see cref="TCoffeeAbortedException"/>.</summary>
    Subop1Pair = 10,

    /// <summary><c>subop2_pair</c>: as <see cref="Subop1Pair"/>.</summary>
    Subop2Pair = 11,

    /// <summary><c>cfast_pair</c>: "Mcfast_pair not recognised [FATAL:T-COFFEE]", exit status 1. Raises <see cref="TCoffeeAbortedException"/>.</summary>
    CfastPair = 12,

    /// <summary>
    /// <c>fast_pair</c>: the FASTA heuristic of the program (<c>fasta_gotoh_pair_wise</c>, util_dp_fasta_nw.c): the DP is restricted to the diagonals with the best k-tuple scores (k-tuple 2 protein / 5 DNA,
    /// z-score threshold, or the <see cref="TCoffeeProgramOptions.DiagonalCount"/> best diagonals); DNA: matrix dna_idmat, GOP -50, GEP -1; protein: BLOSUM62 with GOP (int)(10 x mean mismatch), GEP -1; terminal gaps pay the extension only.
    /// </summary>
    FastPair = 13,

    /// <summary>
    /// <c>exon_pair</c> (<c>SEQ_TYPE G</c>): the Myers-Miller alignment (as <c>slow_pair</c>) of the GENOMIC templates (<c>_G_</c>, <see cref="TCoffeeProgramOptions.GenomicTemplates"/>) of the two sequences with the matrix exon2mt (BLOSUM62, the diagonal of o / b / j = 100, the x rows zeroed), GOP -10, GEP -1;
    /// the alignment of the templates is expressed in the residues of the sequences (<c>fix_seq_seq</c>). Only pairs of sequences that both have a template are aligned.
    /// </summary>
    ExonPair = 14,

    /// <summary><c>exon2_pair</c>: as <see cref="ExonPair"/>, but only the columns whose two template letters are in the set "subset_objOBJ-" (the program passes the text after the first '_' of the weight mode "sim_subset_objOBJ-", prefix included) enter the library.</summary>
    Exon2Pair = 15,

    /// <summary><c>exon3_pair</c>: as <see cref="Exon2Pair"/> with GOP 0.</summary>
    Exon3Pair = 16,

    /// <summary>
    /// <c>muscle_pair</c> (audit A2-11): the program runs the external <c>muscle -in x -out y -quiet -maxmb 0</c> on every pair of sequences and reads the alignment; here the MUSCLE 3.8.1551 port of this repository
    /// (<see cref="SequenceAligner.MultipleAlignMuscle(IEnumerable{string}, MuscleOptions?, MatrixScoring?)"/>, the program's default <c>-maxiters 8</c>) is the aligner. Weights = 10 x percent identity as for the other methods.
    /// The installed program finds MUSCLE only when the plugin directory of T-Coffee (<c>$HOME/.t_coffee/plugins/linux</c>) holds a MUSCLE 3.x binary named <c>muscle</c>.
    /// </summary>
    MusclePair = 17,

    /// <summary>
    /// <c>TMalign_pair</c> (audit A2-10): <c>pdb_pair</c> with the external structure aligner TM-align on the <c>_P_</c> templates of the two sequences (<see cref="TCoffeeProgramOptions.StructureTemplates"/>); the alignment it wrote
    /// (<see cref="TCoffeeProgramOptions.StructureAlignerOutputs"/>) becomes the library of the pair (weight 10 x identity). Pairs without two usable templates are aligned by <c>proba_pair</c>.
    /// </summary>
    TMalignPair = 18,

    /// <summary><c>mustang_pair</c> (audit A2-10): as <see cref="TMalignPair"/> with MUSTANG (<c>mustang -i pdb1 pdb2 -F fasta</c>, the file results.afasta).</summary>
    MustangPair = 19,

    /// <summary>
    /// <c>blastp_msa</c> (audit A2-9; ALN_MODE multiple): the sequences are searched against a local BLAST database made of themselves (<c>blastall -p blastp -m7</c>) and the high scoring pairs become the library (weight = percent
    /// identity of the HSP, truncated; <see cref="SequenceAligner.TCoffeeBlastMsaLibraryEntries"/>). The BLAST search is the external program: its XML is <see cref="TCoffeeProgramOptions.BlastXml"/>.
    /// </summary>
    BlastpMsa = 20,

    /// <summary>
    /// <c>sap_pair</c> (audit A2-12): when SAP is not installed - the case of the installed program - the method list of the program switches to TM-align ("SAP is not installed. TMalign will be used instead"), so this method behaves exactly as
    /// <see cref="TMalignPair"/> (<c>-mode 3dcoffee</c> adds it to the default <c>proba_pair</c>: 7 of 8 sets equal <c>-method proba_pair,sap_pair</c>, see <see cref="TCoffeeProgramOptions"/>). The real SAP (Taylor, Orengo 1989) is not reproduced: its glue (<c>sap_pair</c> of util_dp_drivers.c) is another parser (docs/Validation/review-2026-09/handoff/A2-12).
    /// </summary>
    SapPair = 21,

    /// <summary>
    /// <c>dali_pair</c> (audit A2-12): <c>pdb_pair</c> with <c>EXECUTABLE2 DaliLite</c>; DaliLite is not installed, so the program stops ("ERROR: UNSPECIFIED UNSPECIFIED", exit status 1) at the first pair of two sequences that both carry a usable
    /// structure (<see cref="TCoffeeProgramOptions.StructureTemplates"/>) and raises <see cref="TCoffeeAbortedException"/>; pairs without two structures are aligned by <c>proba_pair</c> as for every structure method.
    /// </summary>
    DaliPair = 22,

    /// <summary><c>fugue_pair</c> (audit A2-12): needs <c>fugueali</c> (FUGUE: licence-gated, no public source), which is not installed: "The Program fugueali Needed by T-COFFEE Could not be found [FATAL:T-COFFEE]", exit status 1. Raises <see cref="TCoffeeAbortedException"/>.</summary>
    FuguePair = 23,

    /// <summary>
    /// <c>mafft_pair</c> (audit A2-11): the program runs the external <c>mafft --anysymbol in &gt; out</c> on every pair (2 x 23 - the installed mafft 7.505, FFT-NS-1 progressive alignment of the two sequences). The earlier "rc 1" was not an
    /// impossibility: t_coffee points <c>MAFFT_BINARIES</c> at <c>$HOME/.t_coffee/plugins/linux</c>, which must hold the helper binaries of mafft (<c>/usr/lib/mafft/bin</c> and <c>/usr/lib/mafft/lib/mafft</c>) - 460 oracle sets are
    /// committed (feature <c>mafft_pair</c> of <c>tcoffee_program_options_oracle.tsv</c>). There is no MAFFT port in this repository yet (B09.md A2-19), so the method raises <see cref="NotSupportedException"/>.
    /// </summary>
    MafftPair = 24,
}

/// <summary>
/// Per-method overrides of the method-file parameters, the <c>-method slow_pair@EP@GOP@-30@GEP@-3@MATRIX@blosum50mt</c> syntax of the program (each <c>@EP@KEY@value</c> pair appends the line
/// <c>KEY value</c> to the method file, the last line wins). <c>null</c> = the method's own value.
/// </summary>
/// <param name="GapOpen"><c>GOP</c> (negative).</param>
/// <param name="GapExtend"><c>GEP</c>.</param>
/// <param name="Matrix"><c>MATRIX</c>: any matrix name of the program.</param>
public sealed record TCoffeeMethodParameters(int? GapOpen = null, int? GapExtend = null, string? Matrix = null);

/// <summary>Residue pair of a T-Coffee library: 0-based sequence indices and 0-based residue positions, weight and the number of times the pair was added (<c>CONS</c>).</summary>
public readonly record struct TCoffeeLibraryEntry(int SeqA, int PosA, int SeqB, int PosB, int Weight, int Count);

/// <summary>A local alignment reported by SIM (<see cref="SequenceAligner.TCoffeeLalignTopAlignments"/>).</summary>
/// <param name="RowA">Aligned row of the first sequence ('-' = gap).</param>
/// <param name="RowB">Aligned row of the second sequence.</param>
/// <param name="StartA">0-based position of the first aligned residue of the first sequence.</param>
/// <param name="StartB">0-based position of the first aligned residue of the second sequence.</param>
/// <param name="Score">SIM score (matrix x SCORE_K, gaps q + r k).</param>
public sealed record TCoffeeLocalAlignment(string RowA, string RowB, int StartA, int StartB, int Score);

public static partial class SequenceAligner
{
    private const int TcScoreK = 10;

    // ----- scoring matrices (matrices.h + read_matrice) -----------------------------------------------------------------

    /// <summary>
    /// A T-Coffee scoring matrix as <c>read_matrice</c> builds it: the triangular vector over the 23 letters "abcdefghiklmnpqrstvwxyz" copied to
    /// every letter case combination, with U taking the row/column of T (the "RNA correction"); letters outside the 23 (J, O) score 0
    /// against everything including themselves.
    /// </summary>
    internal sealed class TcMatrix
    {
        private const string Alphabet = "abcdefghiklmnpqrstvwxyz";
        private readonly int[,] _m = new int[256, 256];

        private TcMatrix(IReadOnlyList<int> triangle, bool exon2 = false)
        {
            int l = Alphabet.Length;
            var raw = new int[l, l];
            for (int a = 0; a < l; a++)
                for (int b = 0; b <= a; b++)
                    raw[a, b] = raw[b, a] = triangle[(a * a + a) / 2 + b];

            for (int b = 0; b < l; b++)
                for (int c = 0; c < l; c++)
                {
                    char x = Alphabet[b], y = Alphabet[c];
                    Set(char.ToUpperInvariant(x), char.ToUpperInvariant(y), raw[b, c]);
                    Set(char.ToUpperInvariant(x), y, raw[b, c]);
                    Set(x, char.ToUpperInvariant(y), raw[b, c]);
                    Set(x, y, raw[b, c]);
                }

            // exon2mt: the diagonal of the letters o, b and j is 100 (EXONVALUE_4_TCOFFEE is read with atoi of the variable's NAME, so the environment cannot change it: 100 unless the variable exists, then 0; it does not)
            if (exon2)
                foreach (char e in "obj")
                {
                    char up = char.ToUpperInvariant(e);
                    _m[up - 'A', up - 'A'] = _m[e - 'A', up - 'A'] = _m[up - 'A', e - 'A'] = _m[e - 'A', e - 'A'] = 100;
                }

            // "Correct for RNA: Add the U cost"
            for (int b = 0; b < l; b++)
            {
                char up = char.ToUpperInvariant(Alphabet[b]), lo = Alphabet[b];
                _m['U' - 'A', up - 'A'] = _m['T' - 'A', up - 'A'];
                _m[up - 'A', 'U' - 'A'] = _m[up - 'A', 'T' - 'A'];
                _m['U' - 'A', up - 'A'] = _m['t' - 'A', lo - 'A'];
                _m[lo - 'A', 'u' - 'A'] = _m[lo - 'A', 't' - 'A'];
            }

            _m['U' - 'A', 'U' - 'A'] = _m['T' - 'A', 'T' - 'A'];
            _m['u' - 'A', 'U' - 'A'] = _m['t' - 'A', 'T' - 'A'];
            _m['U' - 'A', 'u' - 'A'] = _m['T' - 'A', 't' - 'A'];
            _m['u' - 'A', 'u' - 'A'] = _m['t' - 'A', 't' - 'A'];
        }

        private TcMatrix(TcMatrix other) => Array.Copy(other._m, _m, other._m.Length);

        /// <summary>A copy whose rows and columns of the letter x (both cases) are 0 (<c>PW_CL->M['x' - 'A'][a] = PW_CL->M[a]['x' - 'A'] = 0</c> for a &lt; 60 in <c>method2pw_cl</c> of the exon methods).</summary>
        public TcMatrix WithZeroX()
        {
            var copy = new TcMatrix(this);
            for (int a = 0; a < 60; a++)
                foreach (char x in "xX")
                    copy._m[x - 'A', a] = copy._m[a, x - 'A'] = 0;
            return copy;
        }

        private void Set(char x, char y, int v) => _m[x - 'A', y - 'A'] = v;

        /// <summary><c>CL->M[a - 'A'][b - 'A']</c> (without the SCORE_K factor).</summary>
        public int Score(char a, char b) => _m[a - 'A', b - 'A'];

        private static TcMatrix? _idmat, _blosum62, _blosum50;
        private static readonly Dictionary<string, TcMatrix?> ByNameCache = new(StringComparer.Ordinal);

        /// <summary><c>read_matrice(name)</c>: any header matrix of <c>matrices.h</c> (see <see cref="TcMatrixCatalog"/>); <c>null</c> when the program does not know the name.</summary>
        public static TcMatrix? ByName(string name)
        {
            lock (ByNameCache)
            {
                if (ByNameCache.TryGetValue(name, out var cached)) return cached;
                var vector = TcMatrixCatalog.Find(name);
                var m = vector == null ? null : new TcMatrix(vector, name == "exon2mt");
                ByNameCache[name] = m;
                return m;
            }
        }

        /// <summary><c>idmat</c> = what <c>read_matrice("dna_idmat")</c> returns: +10 on the diagonal, 0 elsewhere.</summary>
        public static TcMatrix IdMat => _idmat ??= new TcMatrix(Enumerable.Range(0, 276).Select(i => IsDiagonal(i) ? 10 : 0).ToArray());

        /// <summary><c>blosum62mt</c> (matrices.h).</summary>
        public static TcMatrix Blosum62 => _blosum62 ??= new TcMatrix(TcMatrixTables.Blosum62Mt);

        /// <summary><c>blosum50mt</c> (matrices.h).</summary>
        public static TcMatrix Blosum50 => _blosum50 ??= new TcMatrix(TcMatrixTables.Blosum50Mt);

        private static bool IsDiagonal(int index)
        {
            // triangular vector: row a occupies [(a*a+a)/2, (a*a+a)/2 + a]; the diagonal is the last cell of the row
            int a = (int)((Math.Sqrt(8.0 * index + 1) - 1) / 2);
            return index == (a * a + a) / 2 + a;
        }

        /// <summary><c>get_avg_matrix_mm(M, AA_ALPHABET)</c>: float average of the off-diagonal scores of the first 20 letters "acdefghiklmnpqrstvwy".</summary>
        public float AverageMismatch()
        {
            const string aa = "acdefghiklmnpqrstvwy";
            float naa = 0, gop = 0;
            foreach (char a in aa)
                foreach (char b in aa)
                    if (a != b)
                    {
                        gop += Score(a, b);
                        naa++;
                    }

            return gop / naa;
        }
    }

    // ----- parameters of the program's pair-wise methods -----------------------------------------------------------------

    private readonly record struct TcPairParams(TcMatrix Matrix, int Gop, int Gep);

    /// <summary>
    /// The parameters of <c>slow_pair</c>: DNA / RNA take <c>MATRIX dna_idmat</c>, GOP -10, GEP -1 from the method file (so <c>-matrix</c> does not apply), PROTEIN has no MATRIX line
    /// and uses <c>-matrix</c> (<paramref name="matrixName"/>, default blosum62mt) with GOP <c>(int)(10 x get_avg_matrix_mm)</c> of that matrix.
    /// </summary>
    private static TcPairParams TcSlowPairParams(bool protein, string? matrixName = null, TCoffeeMethodParameters? mp = null)
    {
        if (!protein)
            return new TcPairParams(TcResolveMatrix(mp?.Matrix) ?? TcMatrix.IdMat, mp?.GapOpen ?? -10, mp?.GapExtend ?? -1);
        var matrix = TcResolveMatrix(mp?.Matrix) ?? TcResolveMatrix(matrixName) ?? TcMatrix.Blosum62;
        return new TcPairParams(matrix, mp?.GapOpen ?? (int)(matrix.AverageMismatch() * 10), mp?.GapExtend ?? -1);
    }

    /// <summary>The matrix of <c>-matrix NAME</c> (<c>null</c> / "default" = the method's own); an unknown name is an error as in the program.</summary>
    internal static TcMatrix? TcResolveMatrix(string? name)
    {
        if (string.IsNullOrEmpty(name) || name == "default") return null;
        return TcMatrix.ByName(name) ?? throw new ArgumentException($"-matrix {name} is not a matrix of the program (read_matrice).");
    }

    /// <summary>The parameters of <c>fast_pair</c>: DNA / RNA <c>MATRIX dna_idmat</c>, GOP -50, GEP -1 (method file); PROTEIN <c>-matrix</c> (default blosum62mt), GOP <c>(int)(10 x get_avg_matrix_mm)</c>, GEP -1.</summary>
    private static TcPairParams TcFastPairParams(bool protein, string? matrixName = null, TCoffeeMethodParameters? mp = null)
    {
        if (!protein)
            return new TcPairParams(TcResolveMatrix(mp?.Matrix) ?? TcMatrix.IdMat, mp?.GapOpen ?? -50, mp?.GapExtend ?? -1);
        var matrix = TcResolveMatrix(mp?.Matrix) ?? TcResolveMatrix(matrixName) ?? TcMatrix.Blosum62;
        return new TcPairParams(matrix, mp?.GapOpen ?? (int)(matrix.AverageMismatch() * 10), mp?.GapExtend ?? -1);
    }

    private static TcPairParams TcLalignParams(bool protein, TCoffeeMethodParameters? mp = null)
    {
        var d = protein ? new TcPairParams(TcMatrix.Blosum50, -10, -4) : new TcPairParams(TcMatrix.IdMat, -10, -1);
        return new TcPairParams(TcResolveMatrix(mp?.Matrix) ?? d.Matrix, mp?.GapOpen ?? d.Gop, mp?.GapExtend ?? d.Gep);
    }

    /// <summary>
    /// The program's sequence type decision (<c>get_string_type</c> on the concatenation of all sequences): DNA/RNA when more than 98 % of the residues are letters of
    /// <c>AGCTUNRYMKSWHBVD</c> (either case), protein otherwise. The type selects the matrices, gap penalties and k-tuple alphabet of the methods.
    /// </summary>
    internal static bool TcIsProtein(IEnumerable<string> sequences) => TcSeqType(sequences) == TcSeqKind.Protein;

    /// <summary>The full type decision: <c>RNA</c> when the nucleotide share is above 98 % and at least one residue is U (<c>is_rna</c>), DNA without a U, protein otherwise.</summary>
    internal static TcSeqKind TcSeqType(IEnumerable<string> sequences)
    {
        int dna = 0, tot = 0, rna = 0;
        foreach (var s in sequences)
            foreach (char c in s)
            {
                if ("-.#*~".Contains(c)) continue;
                tot++;
                if ("AGCTUNRYMKSWHBVDagctunrymkswhbvd".Contains(c)) dna++;
                if (c == 'U' || c == 'u') rna++;
            }

        if (tot == 0 || dna * 100 / tot <= 98) return TcSeqKind.Protein;
        return rna > 0 ? TcSeqKind.Rna : TcSeqKind.Dna;
    }

    /// <summary>The type to use: the caller's choice (<c>protein</c>: true = protein, false = nucleotide, nucleotide = DNA or RNA as the letters say) or the program's decision.</summary>
    internal static TcSeqKind TcKindOf(IEnumerable<string> sequences, bool? protein)
    {
        var decided = TcSeqType(sequences);
        if (protein == null) return decided;
        if (protein.Value) return TcSeqKind.Protein;
        return decided == TcSeqKind.Rna ? TcSeqKind.Rna : TcSeqKind.Dna;
    }

    // ----- Myers-Miller (slow_pair) ---------------------------------------------------------------------------------------

    /// <summary>
    /// <c>slow_pair</c>: the global alignment of two sequences computed by the Myers-Miller routine of the T-Coffee program
    /// (<c>myers_miller_pair_wise</c>, util_dp_mm_nw.c), with the program's DNA (<c>dna_idmat</c>, gap open -10, extend -1) or protein
    /// (BLOSUM62, gap open (int)(10 x mean mismatch score), extend -1) parameters.
    /// </summary>
    /// <param name="sequenceA">First sequence (letters).</param>
    /// <param name="sequenceB">Second sequence.</param>
    /// <param name="protein"><c>true</c> = the program's protein parameters; <c>false</c> = DNA/RNA.</param>
    /// <param name="matrix"><c>-matrix</c> (protein only; DNA / RNA always use dna_idmat): any matrix name of the program, default blosum62mt.</param>
    /// <param name="parameters">The <c>@EP@</c> overrides of the method (GOP, GEP, MATRIX; they also apply to DNA / RNA).</param>
    /// <returns>The two aligned rows.</returns>
    public static (string RowA, string RowB) TCoffeeSlowPairAlign(string sequenceA, string sequenceB, bool protein = false, string? matrix = null, TCoffeeMethodParameters? parameters = null)
    {
        ArgumentNullException.ThrowIfNull(sequenceA);
        ArgumentNullException.ThrowIfNull(sequenceB);
        var p = TcSlowPairParams(protein, matrix, parameters);
        return TcMyersMiller(sequenceA, sequenceB, p.Matrix, p.Gop, p.Gep);
    }

    private sealed class MyersMillerState
    {
        private readonly Func<int, int, int> _cost;
        public readonly int G, H;
        public readonly int[] CC, DD, RR, SS;
        public readonly List<int> Script = new();
        public int Last;

        /// <summary>The Myers-Miller state over a score function of the two 0-based positions (<c>get_dp_cost</c>; profile columns or residues).</summary>
        public MyersMillerState(int lengthA, int lengthB, Func<int, int, int> cost, int gop, int gep)
        {
            _cost = cost;
            G = gop * TcScoreK;
            H = gep * TcScoreK;
            int l = lengthA + lengthB + 1;
            CC = new int[l];
            DD = new int[l];
            RR = new int[l];
            SS = new int[l];
        }

        public MyersMillerState(string a, string b, TcMatrix matrix, int gop, int gep)
            : this(a.Length, b.Length, (i, j) => matrix.Score(a[i], b[j]) * TcScoreK, gop, gep)
        {
        }

        public int Cost(int i, int j) => _cost(i, j);

        private int Gap(int k) => k <= 0 ? 0 : G + H * k;

        private void Del(int k)
        {
            if (Last < 0)
            {
                Script[^1] -= k;
                Last = Script[^1];
            }
            else
            {
                Script.Add(-k);
                Last = -k;
            }
        }

        private void Ins(int k)
        {
            if (Last < 0)
            {
                int last = Last;
                Script[^1] = k;
                Script.Add(last);
            }
            else
            {
                Script.Add(k);
                Last = k;
            }
        }

        private void Rep()
        {
            Script.Add(0);
            Last = 0;
        }

        public int Diff(int s1, int M, int s2, int N, int tb, int te)
        {
            int g = G, h = H, m = g + h;
            int c, e, d, s, t;

            if (N <= 0)
            {
                if (M > 0) Del(M);
                return Gap(M);
            }

            if (M <= 1)
            {
                if (M <= 0)
                {
                    Ins(N);
                    return Gap(N);
                }

                if (tb > te) tb = te;
                int midc0 = (tb + h) + Gap(N);
                int midj0 = 0;
                for (int j = 1; j <= N; j++)
                {
                    c = Gap(j - 1) + Cost(s1, j - 1 + s2) + Gap(N - j);
                    if (c > midc0)
                    {
                        midc0 = c;
                        midj0 = j;
                    }
                }

                if (midj0 == 0)
                {
                    Del(1);
                    Ins(N);
                }
                else
                {
                    if (midj0 > 1) Ins(midj0 - 1);
                    Rep();
                    if (midj0 < N) Ins(N - midj0);
                }

                return midc0;
            }

            // Divide: find the optimum midpoint (midi, midj) of cost midc
            int midi = M / 2;
            CC[0] = 0;
            t = tb;
            for (int j = 1; j <= N; j++)
            {
                CC[j] = t = t + h;
                DD[j] = t + g;
            }

            t = tb;
            for (int i = 1; i <= midi; i++)
            {
                s = CC[0];
                CC[0] = c = t = t + h;
                e = t + g;
                for (int j = 1; j <= N; j++)
                {
                    if ((c = c + m) > (e = e + h)) e = c;
                    if ((c = CC[j] + m) > (d = DD[j] + h)) d = c;
                    c = s + Cost(i - 1 + s1, j - 1 + s2);
                    if (e > c) c = e;
                    if (d > c) c = d;
                    s = CC[j];
                    CC[j] = c;
                    DD[j] = d;
                }
            }

            DD[0] = CC[0];

            RR[N] = 0;
            t = te;
            for (int j = N - 1; j >= 0; j--)
            {
                RR[j] = t = t + h;
                SS[j] = t + g;
            }

            t = te;
            for (int i = M - 1; i >= midi; i--)
            {
                s = RR[N];
                RR[N] = c = t = t + h;
                e = t + g;
                for (int j = N - 1; j >= 0; j--)
                {
                    if ((c = c + m) > (e = e + h)) e = c;
                    if ((c = RR[j] + m) > (d = SS[j] + h)) d = c;
                    c = s + Cost(i + s1, j + s2);
                    if (e > c) c = e;
                    if (d > c) c = d;
                    s = RR[j];
                    RR[j] = c;
                    SS[j] = d;
                }
            }

            SS[N] = RR[N];
            int midc = CC[0] + RR[0];
            int midj = 0;
            int type = 1;
            for (int j = 0; j <= N; j++)
                if ((c = CC[j] + RR[j]) >= midc && (c > midc || (CC[j] != DD[j] && RR[j] == SS[j])))
                {
                    midc = c;
                    midj = j;
                }

            for (int j = N; j >= 0; j--)
                if ((c = DD[j] + SS[j] - g) > midc)
                {
                    midc = c;
                    midj = j;
                    type = 2;
                }

            // Conquer
            if (type == 1)
            {
                Diff(s1, midi, s2, midj, tb, G);
                Diff(s1 + midi, M - midi, s2 + midj, N - midj, G, te);
            }
            else
            {
                Diff(s1, midi - 1, s2, midj, tb, 0);
                Del(2);
                Diff(s1 + midi + 1, M - midi - 1, s2 + midj, N - midj, 0, te);
            }

            return midc;
        }
    }

    private static (string, string) TcMyersMiller(string a, string b, TcMatrix matrix, int gop, int gep)
    {
        if (a.Length == 0 || b.Length == 0) throw new ArgumentException("Empty sequences are not supported.");
        var st = new MyersMillerState(a, b, matrix, gop, gep);
        st.Diff(0, a.Length, 0, b.Length, 0, 0);
        return TcScriptToRows(a, b, 0, 0, a.Length, b.Length, st.Script);
    }

    /// <summary>Applies an edit script (0 = pair, k &gt; 0 = k residues of the second sequence against gaps, k &lt; 0 = -k residues of the first sequence against gaps).</summary>
    private static (string, string) TcScriptToRows(string a, string b, int startA, int startB, int lenA, int lenB, List<int> script)
    {
        var ra = new System.Text.StringBuilder();
        var rb = new System.Text.StringBuilder();
        int i = startA, j = startB, endA = startA + lenA, endB = startB + lenB;
        int idx = 0;
        while (!(i == endA && j == endB))
        {
            int op = script[idx++];
            if (op == 0)
            {
                ra.Append(a[i++]);
                rb.Append(b[j++]);
            }
            else if (op > 0)
            {
                for (int k = 0; k < op; k++)
                {
                    ra.Append('-');
                    rb.Append(b[j++]);
                }
            }
            else
            {
                for (int k = 0; k < -op; k++)
                {
                    ra.Append(a[i++]);
                    rb.Append('-');
                }
            }
        }

        return (ra.ToString(), rb.ToString());
    }

    // ----- SIM (lalign_id_pair) ------------------------------------------------------------------------------------------

    /// <summary>
    /// <c>lalign_id_pair</c>: the <paramref name="topAlignments"/> best NON-INTERSECTING local alignments of two sequences as the T-Coffee
    /// program's SIM (Huang &amp; Miller 1991, util_dp_sim.c <c>sim_pair_wise_lalign</c>) reports them, in report order (score order).
    /// DNA: <c>dna_idmat</c> (+10 / 0), gap open 10 + 1 per residue scaled x 10; protein: BLOSUM50, open 10 / extend 4.
    /// </summary>
    /// <param name="sequenceA">First sequence.</param>
    /// <param name="sequenceB">Second sequence (must be a different sequence of the set: the self-comparison mode is not used by the program's libraries).</param>
    /// <param name="topAlignments"><c>-lalign_n_top</c> (the program's default is 10).</param>
    /// <param name="protein"><c>true</c> = the program's protein parameters.</param>
    /// <param name="parameters">The <c>@EP@</c> overrides of the method (GOP, GEP, MATRIX).</param>
    /// <returns>The reported alignments (possibly none).</returns>
    public static IReadOnlyList<TCoffeeLocalAlignment> TCoffeeLalignTopAlignments(string sequenceA, string sequenceB, int topAlignments = 10, bool protein = false, TCoffeeMethodParameters? parameters = null)
    {
        ArgumentNullException.ThrowIfNull(sequenceA);
        ArgumentNullException.ThrowIfNull(sequenceB);
        ArgumentOutOfRangeException.ThrowIfLessThan(topAlignments, 1);
        var p = TcLalignParams(protein, parameters);
        return new TcSim(sequenceA, sequenceB, p.Matrix, p.Gop, p.Gep).Run(topAlignments);
    }

    internal static IReadOnlyList<TCoffeeLocalAlignment> TcSimRun(string a, string b, int gop, int gep, int k) =>
        new TcSim(a, b, TcMatrix.IdMat, gop, gep).Run(k);

    internal static (int MaxScript, int Capacity) TcSimScriptUse(string a, string b, int gop, int gep, int k)
    {
        var sim = new TcSim(a, b, TcMatrix.IdMat, gop, gep);
        sim.Run(k);
        return (sim.MaxScriptLength, sim.ScriptCapacity);
    }

    private sealed class TcSim
    {
        private sealed class Vertex
        {
            public int Score, StarI, StarJ, EndI, EndJ, Top, Bot, Left, Right;
        }

        private readonly string _a, _b;
        private readonly TcMatrix _matrix;
        private readonly int _m, _n;
        private readonly int _q, _r, _qr;

        private int[] _cc = null!, _dd = null!, _rr = null!, _ss = null!, _ee = null!, _ff = null!;
        private int[] _hh = null!, _ww = null!, _ii = null!, _jj = null!, _xx = null!, _yy = null!;
        private HashSet<int>[] _row = null!;
        private Vertex[] _list = null!;
        private Vertex? _low, _most;
        private int _numnode, _lmin, _flag;
        private int _m1, _mm, _n1, _nn, _rl, _cl;

        // diff_sim state
        private int[] _script = null!;
        private int _maxScript;

        /// <summary>The longest edit script of any reported alignment (the program's buffer <c>S</c> holds min(M, N) x 5 / 4 + 1 entries).</summary>
        public int MaxScriptLength => _maxScript;

        public int ScriptCapacity => Math.Min(_m, _n) * 5 / 4 + 1;

        /// <summary>Integers the glibc chunk of the script array holds: a request of R bytes gets ((R + 8 + 15) / 16 x 16) - 8 bytes.</summary>
        public int ScriptUsableEntries => (((ScriptCapacity * 4 + 8 + 15) / 16 * 16) - 8) / 4;
        private int _scriptLen, _last, _i, _j;

        public TcSim(string a, string b, TcMatrix matrix, int gop, int gep)
        {
            _a = a;
            _b = b;
            _matrix = matrix;
            _m = a.Length;
            _n = b.Length;
            _q = Math.Max(gop, -gop) * TcScoreK;
            _r = Math.Max(gep, -gep) * TcScoreK;
            _qr = _q + _r;
        }

        // TC_SCORE_2(x, y): x, y are 0-based residue indices
        private int Score(int x, int y) => _matrix.Score(_a[x], _b[y]) * TcScoreK;

        private bool Used(int ii, int jj) => _row[ii].Contains(jj);

        private static void Order(ref int ss1, ref int xx1, ref int yy1, int ss2, int xx2, int yy2)
        {
            if (ss1 < ss2)
            {
                ss1 = ss2;
                xx1 = xx2;
                yy1 = yy2;
            }
            else if (ss1 == ss2)
            {
                if (xx1 < xx2)
                {
                    xx1 = xx2;
                    yy1 = yy2;
                }
                else if (xx1 == xx2 && yy1 < yy2)
                {
                    yy1 = yy2;
                }
            }
        }

        public IReadOnlyList<TCoffeeLocalAlignment> Run(int topAlignments)
        {
            var result = new List<TCoffeeLocalAlignment>();
            int M = _m, N = _n;
            int K = topAlignments;

            _cc = new int[N + 1]; _dd = new int[N + 1]; _rr = new int[N + 1]; _ss = new int[N + 1]; _ee = new int[N + 1]; _ff = new int[N + 1];
            _hh = new int[M + 1]; _ww = new int[M + 1]; _ii = new int[M + 1]; _jj = new int[M + 1]; _xx = new int[M + 1]; _yy = new int[M + 1];
            _script = new int[M + N + 2];
            _row = new HashSet<int>[M + 1];
            for (int i = 0; i <= M; i++) _row[i] = new HashSet<int>();

            _list = new Vertex[K];
            for (int i = 0; i < K; i++) _list[i] = new Vertex();
            _numnode = _lmin = 0;
            _low = _most = null;

            BigPass(K);

            for (int count = K - 1; count >= 0; count--)
            {
                if (_numnode == 0) break;
                var cur = FindMax();
                int score = cur.Score;
                int stari = ++cur.StarI;
                int starj = ++cur.StarJ;
                int endi = cur.EndI, endj = cur.EndJ;
                _m1 = cur.Top;
                _mm = cur.Bot;
                _n1 = cur.Left;
                _nn = cur.Right;
                _rl = endi - stari + 1;
                _cl = endj - starj + 1;
                _i = stari - 1;
                _j = starj - 1;
                _scriptLen = 0;
                _last = 0;
                DiffSim(stari - 1, starj - 1, _rl, _cl, _q, _q);
                _maxScript = Math.Max(_maxScript, _scriptLen);
                if (_scriptLen > ScriptUsableEntries) throw new TCoffeeAbortedException(_scriptLen, ScriptUsableEntries);

                var script = _script.Take(_scriptLen).ToList();
                var (ra, rb) = TcScriptToRows(_a, _b, stari - 1, starj - 1, _rl, _cl, script);
                result.Add(new TCoffeeLocalAlignment(ra, rb, stari - 1, starj - 1, score));

                if (count > 0)
                {
                    _flag = 0;
                    Locate();
                    if (_flag != 0) SmallPass(count);
                }
            }

            return result;
        }

        // ----- list of the K best classes -----------------------------------------------------------------------------

        private int AddNode(int c, int ci, int cj, int i, int j, int K, int cost)
        {
            bool found = false;
            if (_most != null && _most.StarI == ci && _most.StarJ == cj)
                found = true;
            else
                for (int d = 0; d < _numnode; d++)
                {
                    _most = _list[d];
                    if (_most.StarI == ci && _most.StarJ == cj)
                    {
                        found = true;
                        break;
                    }
                }

            if (found)
            {
                if (_most!.Score < c)
                {
                    _most.Score = c;
                    _most.EndI = i;
                    _most.EndJ = j;
                }

                if (_most.Top > i) _most.Top = i;
                if (_most.Bot < i) _most.Bot = i;
                if (_most.Left > j) _most.Left = j;
                if (_most.Right < j) _most.Right = j;
            }
            else
            {
                _most = _numnode == K ? _low : _list[_numnode++];
                _most!.Score = c;
                _most.StarI = ci;
                _most.StarJ = cj;
                _most.EndI = i;
                _most.EndJ = j;
                _most.Top = _most.Bot = i;
                _most.Left = _most.Right = j;
            }

            if (_numnode == K)
            {
                if (_low == _most || _low == null)
                {
                    _low = _list[0];
                    for (int d = 1; d < _numnode; d++)
                        if (_list[d].Score < _low.Score)
                            _low = _list[d];
                }

                return _low!.Score;
            }

            return cost;
        }

        private Vertex FindMax()
        {
            int j = 0;
            for (int i = 1; i < _numnode; i++)
                if (_list[i].Score > _list[j].Score)
                    j = i;
            var cur = _list[j];
            if (j != --_numnode)
            {
                _list[j] = _list[_numnode];
                _list[_numnode] = cur;
            }

            _most = _list[0];
            if (_low == cur) _low = _list[0];
            return cur;
        }

        // ----- the forward / reverse / recomputation passes (nseq == 2) ----------------------------------------------------

        private void BigPass(int K)
        {
            int M = _m, N = _n;
            for (int j = 1; j <= N; j++)
            {
                _cc[j] = 0;
                _rr[j] = 0;
                _ee[j] = j;
                _dd[j] = -_q;
                _ss[j] = 0;
                _ff[j] = j;
            }

            for (int i = 1; i <= M; i++)
            {
                int c = 0;
                int f = -_q;
                int ci = i, fi = i;
                int p = 0;
                int pi = i - 1;
                int cj = 0, fj = 0, pj = 0;
                for (int j = 1; j <= N; j++)
                {
                    f -= _r;
                    c -= _qr;
                    Order(ref f, ref fi, ref fj, c, ci, cj);
                    c = _cc[j] - _qr;
                    ci = _rr[j];
                    cj = _ee[j];
                    int d = _dd[j] - _r;
                    int di = _ss[j];
                    int dj = _ff[j];
                    Order(ref d, ref di, ref dj, c, ci, cj);
                    c = 0;
                    if (!Used(i, j)) c = p + Score(i - 1, j - 1);
                    if (c <= 0)
                    {
                        c = 0;
                        ci = i;
                        cj = j;
                    }
                    else
                    {
                        ci = pi;
                        cj = pj;
                    }

                    Order(ref c, ref ci, ref cj, d, di, dj);
                    Order(ref c, ref ci, ref cj, f, fi, fj);
                    p = _cc[j];
                    _cc[j] = c;
                    pi = _rr[j];
                    pj = _ee[j];
                    _rr[j] = ci;
                    _ee[j] = cj;
                    _dd[j] = d;
                    _ss[j] = di;
                    _ff[j] = dj;
                    if (c > _lmin) _lmin = AddNode(c, ci, cj, i, j, K, _lmin);
                }
            }
        }

        /// <summary>
        /// <c>locate</c>: the reverse pass that determines the left and top boundaries of the area to recompute after an alignment was reported
        /// (nseq == 2 paths; the variables are function-wide as in the C source: the last cell's end points are tested after the loops).
        /// </summary>
        private void Locate()
        {
            int c = 0, f = 0, d = 0, p = 0, ci = 0, cj = 0, di = 0, dj = 0, fi = 0, fj = 0, pi = 0, pj = 0;
            bool cflag, rflag;

            for (int j = _nn; j >= _n1; j--)
            {
                _cc[j] = 0;
                _ee[j] = j;
                _dd[j] = -_q;
                _ff[j] = j;
                _rr[j] = _ss[j] = _mm + 1;
            }

            for (int i = _mm; i >= _m1; i--)
            {
                c = p = 0;
                f = -_q;
                ci = fi = i;
                pi = i + 1;
                cj = fj = pj = _nn + 1;
                for (int j = _nn; j >= _n1; j--)
                {
                    f -= _r;
                    c -= _qr;
                    Order(ref f, ref fi, ref fj, c, ci, cj);
                    c = _cc[j] - _qr;
                    ci = _rr[j];
                    cj = _ee[j];
                    d = _dd[j] - _r;
                    di = _ss[j];
                    dj = _ff[j];
                    Order(ref d, ref di, ref dj, c, ci, cj);
                    c = 0;
                    if (!Used(i, j)) c = p + Score(i - 1, j - 1);
                    if (c <= 0)
                    {
                        c = 0;
                        ci = i;
                        cj = j;
                    }
                    else
                    {
                        ci = pi;
                        cj = pj;
                    }

                    Order(ref c, ref ci, ref cj, d, di, dj);
                    Order(ref c, ref ci, ref cj, f, fi, fj);
                    p = _cc[j];
                    _cc[j] = c;
                    pi = _rr[j];
                    pj = _ee[j];
                    _rr[j] = ci;
                    _ee[j] = cj;
                    _dd[j] = d;
                    _ss[j] = di;
                    _ff[j] = dj;
                    if (c > _lmin) _flag = 1;
                }

                _hh[i] = _cc[_n1];
                _ii[i] = _rr[_n1];
                _jj[i] = _ee[_n1];
                _ww[i] = _dd[_n1];
                _xx[i] = _ss[_n1];
                _yy[i] = _ff[_n1];
            }

            _rl = _m1;
            _cl = _n1;
            while (true)
            {
                rflag = cflag = true;
                while ((rflag && _m1 > 1) || (cflag && _n1 > 1))
                {
                    if (rflag && _m1 > 1) // compute one row
                    {
                        rflag = false;
                        _m1--;
                        c = 0;
                        f = -_q;
                        ci = fi = _m1;
                        pi = _m1 + 1;
                        cj = fj = pj = _nn + 1;
                        for (int j = _nn; j >= _n1; j--)
                        {
                            f -= _r;
                            c -= _qr;
                            Order(ref f, ref fi, ref fj, c, ci, cj);
                            c = _cc[j] - _qr;
                            ci = _rr[j];
                            cj = _ee[j];
                            d = _dd[j] - _r;
                            di = _ss[j];
                            dj = _ff[j];
                            Order(ref d, ref di, ref dj, c, ci, cj);
                            c = 0;
                            if (!Used(_m1, j)) c = Score(_m1 - 1, j - 1);
                            if (c <= 0)
                            {
                                c = 0;
                                ci = _m1;
                                cj = j;
                            }
                            else
                            {
                                ci = pi;
                                cj = pj;
                            }

                            Order(ref c, ref ci, ref cj, d, di, dj);
                            Order(ref c, ref ci, ref cj, f, fi, fj);
                            _cc[j] = c;
                            pi = _rr[j];
                            pj = _ee[j];
                            _rr[j] = ci;
                            _ee[j] = cj;
                            _dd[j] = d;
                            _ss[j] = di;
                            _ff[j] = dj;
                            if (c > _lmin) _flag = 1;
                            if (!rflag && Crosses(ci, cj, di, dj, fi, fj)) rflag = true;
                        }

                        _hh[_m1] = _cc[_n1];
                        _ii[_m1] = _rr[_n1];
                        _jj[_m1] = _ee[_n1];
                        _ww[_m1] = _dd[_n1];
                        _xx[_m1] = _ss[_n1];
                        _yy[_m1] = _ff[_n1];
                        if (!cflag && Crosses(ci, cj, di, dj, fi, fj)) cflag = true;
                    }

                    if (cflag && _n1 > 1) // compute one column
                    {
                        cflag = false;
                        _n1--;
                        c = 0;
                        f = -_q;
                        cj = fj = _n1;
                        p = 0;
                        ci = fi = pi = _mm + 1;
                        pj = _n1 + 1;
                        int limit = _mm;
                        for (int i = limit; i >= _m1; i--)
                        {
                            f -= _r;
                            c -= _qr;
                            Order(ref f, ref fi, ref fj, c, ci, cj);
                            c = _hh[i] - _qr;
                            ci = _ii[i];
                            cj = _jj[i];
                            d = _ww[i] - _r;
                            di = _xx[i];
                            dj = _yy[i];
                            Order(ref d, ref di, ref dj, c, ci, cj);
                            c = 0;
                            if (!Used(i, _n1)) c = p + Score(i - 1, _n1 - 1);
                            if (c <= 0)
                            {
                                c = 0;
                                ci = i;
                                cj = _n1;
                            }
                            else
                            {
                                ci = pi;
                                cj = pj;
                            }

                            Order(ref c, ref ci, ref cj, d, di, dj);
                            Order(ref c, ref ci, ref cj, f, fi, fj);
                            p = _hh[i];
                            _hh[i] = c;
                            pi = _ii[i];
                            pj = _jj[i];
                            _ii[i] = ci;
                            _jj[i] = cj;
                            _ww[i] = d;
                            _xx[i] = di;
                            _yy[i] = dj;
                            if (c > _lmin) _flag = 1;
                            if (!cflag && Crosses(ci, cj, di, dj, fi, fj)) cflag = true;
                        }

                        _cc[_n1] = _hh[_m1];
                        _rr[_n1] = _ii[_m1];
                        _ee[_n1] = _jj[_m1];
                        _dd[_n1] = _ww[_m1];
                        _ss[_n1] = _xx[_m1];
                        _ff[_n1] = _yy[_m1];
                        if (!rflag && Crosses(ci, cj, di, dj, fi, fj)) rflag = true;
                    }
                }

                if ((_m1 == 1 && _n1 == 1) || NoCross()) break;
            }

            _m1--;
            _n1--;
        }

        private bool Crosses(int ci, int cj, int di, int dj, int fi, int fj) =>
            (ci > _rl && cj > _cl) || (di > _rl && dj > _cl) || (fi > _rl && fj > _cl);

        /// <summary><c>no_cross</c>: true when no node of LIST shares vertices with the area (otherwise the area is enlarged).</summary>
        private bool NoCross()
        {
            int i;
            for (i = 0; i < _numnode; i++)
            {
                var cur = _list[i];
                if (cur.StarI <= _mm && cur.StarJ <= _nn && cur.Bot >= _m1 - 1 && cur.Right >= _n1 - 1 && (cur.StarI < _rl || cur.StarJ < _cl))
                {
                    if (cur.StarI < _rl) _rl = cur.StarI;
                    if (cur.StarJ < _cl) _cl = cur.StarJ;
                    _flag = 1;
                    break;
                }
            }

            return i == _numnode;
        }

        /// <summary><c>small_pass</c>: recomputation of the area on a forward pass.</summary>
        private void SmallPass(int count)
        {
            for (int j = _n1 + 1; j <= _nn; j++)
            {
                _cc[j] = 0;
                _rr[j] = _m1;
                _ee[j] = j;
                _dd[j] = -_q;
                _ss[j] = _m1;
                _ff[j] = j;
            }

            for (int i = _m1 + 1; i <= _mm; i++)
            {
                int c = 0;
                int f = -_q;
                int ci = i, fi = i;
                int p = 0;
                int pi = i - 1;
                int cj = _n1, fj = _n1, pj = _n1;
                int limit = _n1 + 1;
                for (int j = limit; j <= _nn; j++)
                {
                    f -= _r;
                    c -= _qr;
                    Order(ref f, ref fi, ref fj, c, ci, cj);
                    c = _cc[j] - _qr;
                    ci = _rr[j];
                    cj = _ee[j];
                    int d = _dd[j] - _r;
                    int di = _ss[j];
                    int dj = _ff[j];
                    Order(ref d, ref di, ref dj, c, ci, cj);
                    c = 0;
                    if (!Used(i, j)) c = p + Score(i - 1, j - 1);
                    if (c <= 0)
                    {
                        c = 0;
                        ci = i;
                        cj = j;
                    }
                    else
                    {
                        ci = pi;
                        cj = pj;
                    }

                    Order(ref c, ref ci, ref cj, d, di, dj);
                    Order(ref c, ref ci, ref cj, f, fi, fj);
                    p = _cc[j];
                    _cc[j] = c;
                    pi = _rr[j];
                    pj = _ee[j];
                    _rr[j] = ci;
                    _ee[j] = cj;
                    _dd[j] = d;
                    _ss[j] = di;
                    _ff[j] = dj;
                    if (c > _lmin) _lmin = AddNode(c, ci, cj, i, j, count, _lmin);
                }
            }
        }

        // ----- diff_sim: the alignment of a class, marking its pairs as used --------------------------------------------------

        private int Gap(int k) => k <= 0 ? 0 : _q + _r * k;

        private void Del(int k)
        {
            _i += k;
            if (_last < 0)
            {
                _script[_scriptLen - 1] -= k;
                _last = _script[_scriptLen - 1];
            }
            else
            {
                _script[_scriptLen++] = -k;
                _last = -k;
            }
        }

        private void Ins(int k)
        {
            _j += k;
            if (_last < 0)
            {
                int last = _last;
                _script[_scriptLen - 1] = k;
                _script[_scriptLen++] = last;
            }
            else
            {
                _script[_scriptLen++] = k;
                _last = k;
            }
        }

        private void Rep()
        {
            _script[_scriptLen++] = 0;
            _last = 0;
        }

        private void DiffSim(int aOff, int bOff, int M, int N, int tb, int te)
        {
            int midi, midj, type, midc;
            int c, e, d, s, t;

            if (N <= 0)
            {
                if (M > 0) Del(M);
                return;
            }

            if (M <= 1)
            {
                if (M <= 0)
                {
                    Ins(N);
                    return;
                }

                if (tb > te) tb = te;
                midc = -(tb + _r + Gap(N));
                midj = 0;
                for (int j = 1; j <= N; j++)
                {
                    if (!Used(_i + 1, j + _j))
                    {
                        c = Score(aOff, bOff + j - 1) - (Gap(j - 1) + Gap(N - j));
                        if (c > midc)
                        {
                            midc = c;
                            midj = j;
                        }
                    }
                }

                if (midj == 0)
                {
                    Ins(N);
                    Del(1);
                }
                else
                {
                    if (midj > 1) Ins(midj - 1);
                    Rep();
                    _i++;
                    _j++;
                    _row[_i].Add(_j);
                    if (midj < N) Ins(N - midj);
                }

                return;
            }

            midi = M / 2;
            _cc[0] = 0;
            t = -_q;
            for (int j = 1; j <= N; j++)
            {
                _cc[j] = t = t - _r;
                _dd[j] = t - _q;
            }

            t = -tb;
            for (int i = 1; i <= midi; i++)
            {
                s = _cc[0];
                _cc[0] = c = t = t - _r;
                e = t - _q;
                for (int j = 1; j <= N; j++)
                {
                    if ((c = c - _qr) > (e = e - _r)) e = c;
                    if ((c = _cc[j] - _qr) > (d = _dd[j] - _r)) d = c;
                    if (!Used(i + _i, j + _j)) c = s + Score(aOff + i - 1, bOff + j - 1);
                    if (c < d) c = d;
                    if (c < e) c = e;
                    s = _cc[j];
                    _cc[j] = c;
                    _dd[j] = d;
                }
            }

            _dd[0] = _cc[0];

            _rr[N] = 0;
            t = -_q;
            for (int j = N - 1; j >= 0; j--)
            {
                _rr[j] = t = t - _r;
                _ss[j] = t - _q;
            }

            t = -te;
            for (int i = M - 1; i >= midi; i--)
            {
                s = _rr[N];
                _rr[N] = c = t = t - _r;
                e = t - _q;
                for (int j = N - 1; j >= 0; j--)
                {
                    if ((c = c - _qr) > (e = e - _r)) e = c;
                    if ((c = _rr[j] - _qr) > (d = _ss[j] - _r)) d = c;
                    if (!Used(i + 1 + _i, j + 1 + _j)) c = s + Score(aOff + i, bOff + j);
                    if (c < d) c = d;
                    if (c < e) c = e;
                    s = _rr[j];
                    _rr[j] = c;
                    _ss[j] = d;
                }
            }

            _ss[N] = _rr[N];

            midc = _cc[0] + _rr[0];
            midj = 0;
            type = 1;
            for (int j = 0; j <= N; j++)
                if ((c = _cc[j] + _rr[j]) >= midc && (c > midc || (_cc[j] != _dd[j] && _rr[j] == _ss[j])))
                {
                    midc = c;
                    midj = j;
                }

            for (int j = N; j >= 0; j--)
                if ((c = _dd[j] + _ss[j] + _q) > midc)
                {
                    midc = c;
                    midj = j;
                    type = 2;
                }

            if (type == 1)
            {
                DiffSim(aOff, bOff, midi, midj, tb, _q);
                DiffSim(aOff + midi, bOff + midj, M - midi, N - midj, _q, te);
            }
            else
            {
                DiffSim(aOff, bOff, midi - 1, midj, tb, 0);
                Del(2);
                DiffSim(aOff + midi + 1, bOff + midj, M - midi - 1, N - midj, 0, te);
            }
        }
    }
}
