
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

/// <summary>The pair-wise methods of the T-Coffee program that build the primary library.</summary>
public enum TCoffeePairMethod
{
    /// <summary><c>slow_pair</c>: Myers-Miller global alignment (the program's DNA default is <c>slow_pair</c> + <c>lalign_id_pair</c> in the paper; see <see cref="ProbaPair"/> for the 13.41 default).</summary>
    SlowPair = 0,

    /// <summary><c>lalign_id_pair</c>: Huang-Miller SIM, the <see cref="TCoffeeProgramOptions.LalignTopAlignments"/> best non-intersecting local alignments.</summary>
    LalignIdPair = 1,

    /// <summary><c>proba_pair</c>: ProbCons-style pair-HMM posterior library (the default method of T-Coffee 13.41).</summary>
    ProbaPair = 2,
}

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

        private TcMatrix(IReadOnlyList<int> triangle)
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

        private void Set(char x, char y, int v) => _m[x - 'A', y - 'A'] = v;

        /// <summary><c>CL->M[a - 'A'][b - 'A']</c> (without the SCORE_K factor).</summary>
        public int Score(char a, char b) => _m[a - 'A', b - 'A'];

        private static TcMatrix? _idmat, _blosum62, _blosum50;

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

    private static TcPairParams TcSlowPairParams(bool protein) =>
        protein
            ? new TcPairParams(TcMatrix.Blosum62, (int)(TcMatrix.Blosum62.AverageMismatch() * 10), -1)
            : new TcPairParams(TcMatrix.IdMat, -10, -1);

    private static TcPairParams TcLalignParams(bool protein) =>
        protein ? new TcPairParams(TcMatrix.Blosum50, -10, -4) : new TcPairParams(TcMatrix.IdMat, -10, -1);

    /// <summary>
    /// The program's sequence type decision (<c>get_string_type</c> on the concatenation of all sequences): DNA/RNA when more than 98 % of the residues are letters of
    /// <c>AGCTUNRYMKSWHBVD</c> (either case), protein otherwise. The type selects the matrices, gap penalties and k-tuple alphabet of the methods.
    /// </summary>
    internal static bool TcIsProtein(IEnumerable<string> sequences)
    {
        int dna = 0, tot = 0;
        foreach (var s in sequences)
            foreach (char c in s)
            {
                if ("-.#*~".Contains(c)) continue;
                tot++;
                if ("AGCTUNRYMKSWHBVDagctunrymkswhbvd".Contains(c)) dna++;
            }

        return tot == 0 || dna * 100 / tot <= 98;
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
    /// <returns>The two aligned rows.</returns>
    public static (string RowA, string RowB) TCoffeeSlowPairAlign(string sequenceA, string sequenceB, bool protein = false)
    {
        ArgumentNullException.ThrowIfNull(sequenceA);
        ArgumentNullException.ThrowIfNull(sequenceB);
        var p = TcSlowPairParams(protein);
        return TcMyersMiller(sequenceA, sequenceB, p.Matrix, p.Gop, p.Gep);
    }

    private sealed class MyersMillerState
    {
        public readonly string A, B;
        public readonly TcMatrix Matrix;
        public readonly int G, H;
        public readonly int[] CC, DD, RR, SS;
        public readonly List<int> Script = new();
        public int Last;

        public MyersMillerState(string a, string b, TcMatrix matrix, int gop, int gep)
        {
            A = a;
            B = b;
            Matrix = matrix;
            G = gop * TcScoreK;
            H = gep * TcScoreK;
            int l = a.Length + b.Length + 1;
            CC = new int[l];
            DD = new int[l];
            RR = new int[l];
            SS = new int[l];
        }

        public int Cost(int i, int j) => Matrix.Score(A[i], B[j]) * TcScoreK;

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

        public void Diff(int s1, int M, int s2, int N, int tb, int te)
        {
            int g = G, h = H, m = g + h;
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

                return;
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
    /// <returns>The reported alignments (possibly none).</returns>
    public static IReadOnlyList<TCoffeeLocalAlignment> TCoffeeLalignTopAlignments(string sequenceA, string sequenceB, int topAlignments = 10, bool protein = false)
    {
        ArgumentNullException.ThrowIfNull(sequenceA);
        ArgumentNullException.ThrowIfNull(sequenceB);
        ArgumentOutOfRangeException.ThrowIfLessThan(topAlignments, 1);
        var p = TcLalignParams(protein);
        return new TcSim(sequenceA, sequenceB, p.Matrix, p.Gop, p.Gep).Run(topAlignments);
    }

    internal static IReadOnlyList<TCoffeeLocalAlignment> TcSimRun(string a, string b, int gop, int gep, int k) =>
        new TcSim(a, b, TcMatrix.IdMat, gop, gep).Run(k);

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
