namespace Seqeron.Genomics.Alignment;

/// <summary>
/// Exact sum-of-pairs multiple alignment (<see cref="SequenceAligner.MultipleAlignExactWithStatistics(IEnumerable{DnaSequence}, ScoringMatrix?, MsaOptions?)"/>) with the size of the problem
/// that was actually solved.
/// </summary>
/// <param name="Alignment">The optimal alignment (rows in input order, majority consensus, optimal sum-of-pairs score).</param>
/// <param name="LowerBound">Sum-of-pairs score of the best heuristic alignment (star / progressive / iterative) that seeds the
/// Carrillo–Lipman bound (0 for fewer than three sequences or when pruning is off).</param>
/// <param name="LatticeCells">Number of cells of the full k-dimensional lattice, Π (length_i + 1).</param>
/// <param name="VisitedCells">Number of lattice cells the dynamic program evaluated (equals <paramref name="LatticeCells"/> without pruning).</param>
/// <param name="StateValues">Number of dynamic-programming values held (cells × gap-history states), the quantity limited by
/// <see cref="MsaOptions.ExactMaxStateCells"/>.</param>
public sealed record ExactMsaResult(
    MultipleAlignmentResult Alignment,
    int LowerBound,
    long LatticeCells,
    long VisitedCells,
    long StateValues)
{
    /// <summary>
    /// <see cref="LowerBound"/> is an integer in units of 1/<c>ScoreScale</c> of the real score: 1 for the scalar <see cref="ScoringMatrix"/>
    /// overloads, the gap-score denominator (2 for 0.5, 10 for 0.1, ...) for the <see cref="MatrixScoring"/> overloads.
    /// </summary>
    public int ScoreScale { get; init; } = 1;
}

public static partial class SequenceAligner
{
    #region Exact sum-of-pairs multiple alignment (k-dimensional dynamic program + Carrillo-Lipman bound)

    // Exact multiple sequence alignment under the sum-of-pairs (SP) objective (audit A1-11).
    //
    // OBJECTIVE. The one of the heuristic aligners (SumOfPairsScore): MsaGapModel.Linear = column sum of pairs (residue pair Match /
    // Mismatch, residue against gap GapExtend, gap against gap 0); MsaGapModel.Affine = the NATURAL affine sum of pairs of Altschul (1989,
    // J Theor Biol 138:297): every pair of rows is projected out, gap-gap columns are removed and every maximal run of columns of the
    // projection in which the same row has the gap costs GapOpen + length * GapExtend. The SP score of an alignment is the sum over all
    // pairs of the score of the pairwise projection, so k = 2 is exactly Needleman-Wunsch / Gotoh (delegated to GlobalAlign / GlobalAlignAffine).
    //
    // ALGORITHM. The optimum is the best path of the k-dimensional lattice (x_1..x_k), 0 <= x_i <= len_i, from the origin to the
    // far corner; a step is a non-empty subset S of the sequences that each consume one residue (column of the alignment): the
    // 2^k - 1 neighbours of Needleman-Wunsch generalised to k sequences (Altschul & Lipman 1989 "Trees, stars, and multiple biological sequence
    // alignment"; for k = 3 the classic 7-neighbour three-dimensional dynamic program, Murata, Richardson & Sussman 1985).
    // Linear gaps: the column cost depends on S only. Natural affine gaps (Altschul 1989): the cost of a gap column of a pair depends on
    // whether the pair's previous projected column was a gap of the same row, i.e. on the state (0 / X / Y) of every pair; the state of
    // the pair (a, b) is decided by which of the two sequences moved last, so the whole state vector is the weak order of the sequences by
    // the time of their last move (13 / 75 / 541 states for k = 3 / 4 / 5, enumerated from the origin by the transition table), and the
    // dynamic program keeps one value per (cell, weak order): the exact optimum, no heuristic profile approximation.
    //
    // CARRILLO-LIPMAN BOUND (Carrillo & Lipman 1988, SIAM J Appl Math 48:1073; used by the MSA program of Lipman, Altschul & Kececioglu
    // 1989, PNAS 86:4412): the SP score of ANY alignment is the sum of the scores of its pairwise projections, and the projection
    // of a path through the cell x onto the pair (a, b) is a pairwise alignment through the cell (x_a, x_b). Hence every alignment
    // through x scores at most  sum over pairs U_ab(x_a, x_b),  U_ab(i, j) = best pairwise alignment of a and b through (i, j) (forward
    // plus backward Gotoh program; with natural affine gaps the maximum over the gap state at the cell, the backward value taking the state
    // into account). A heuristic alignment gives a lower bound L of the optimum; a cell with sum U_ab < L cannot lie on an optimal
    // path and is never created. Additionally (the usual A*-style tightening, same pairwise tables) a state value v at x is dropped when
    // v + sum_ab G_ab(x_a, x_b | the state of the pair) < L, G = best pairwise suffix given the pair's gap state. Both tests only discard
    // paths that cannot reach L, so the optimum is exact; the optimum is found whenever L <= optimum (the heuristic alignment itself is
    // a path of the surviving lattice). The surviving cells are generated by pushing from the origin in lattice order.
    //
    // COMPLEXITY GUARD. Memory = surviving cells x states values (4 bytes each); the dynamic program throws an ArgumentException when
    // that exceeds MsaOptions.ExactMaxStateCells (default 50 000 000) - never an out-of-memory failure. Without pruning the lattice has
    // Π (len_i + 1) cells (NP-completeness of the general problem, Wang & Jiang 1994); the bound makes similar sequences cheap.

    private const int ExactMaxSequencesLinear = 8;
    private const int ExactMaxSequencesAffine = 6;
    private const int ExactNegInf = int.MinValue / 4;

    /// <summary>
    /// Exact sum-of-pairs-optimal multiple alignment: no other alignment of the sequences has a larger
    /// <see cref="SumOfPairsScore(IReadOnlyList{string}, ScoringMatrix?, MsaGapModel)"/> under the chosen gap model
    /// (<see cref="MsaOptions.GapModel"/>: linear = column sum of pairs, affine = Altschul's natural affine sum of pairs).
    /// k = 2 is <see cref="GlobalAlign(DnaSequence, DnaSequence, ScoringMatrix?)"/> / <see cref="GlobalAlignAffine(DnaSequence, DnaSequence, ScoringMatrix?)"/>;
    /// k ≥ 3 is the k-dimensional dynamic program with the Carrillo–Lipman (1988) bound, see <see cref="MultipleAlignExactWithStatistics(IEnumerable{DnaSequence}, ScoringMatrix?, MsaOptions?)"/>.
    /// </summary>
    /// <remarks>Optimal alignments are not unique; the returned one is deterministic. The general problem is NP-complete: the cost grows with
    /// the product of the sequence lengths unless the sequences are similar (the bound prunes the lattice); exceeding
    /// <see cref="MsaOptions.ExactMaxStateCells"/> throws <see cref="ArgumentException"/>. At most 8 sequences (linear) / 6 (affine).</remarks>
    /// <param name="sequences">Sequences to align.</param>
    /// <param name="scoring">Scoring matrix (default <see cref="SimpleDna"/>).</param>
    /// <param name="options">Options: <see cref="MsaOptions.GapModel"/>, <see cref="MsaOptions.ExactMaxStateCells"/>, <see cref="MsaOptions.ExactCarrilloLipman"/>.</param>
    /// <exception cref="ArgumentNullException"><paramref name="sequences"/> is null.</exception>
    /// <exception cref="ArgumentException">Too many sequences, or the dynamic program would exceed <see cref="MsaOptions.ExactMaxStateCells"/>.</exception>
    [System.Runtime.CompilerServices.OverloadResolutionPriority(1)]
    public static MultipleAlignmentResult MultipleAlignExact(
        IEnumerable<DnaSequence> sequences, ScoringMatrix? scoring = null, MsaOptions? options = null) =>
        MultipleAlignExactWithStatistics(sequences, scoring, options).Alignment;

    /// <summary>
    /// <see cref="MultipleAlignExact(IEnumerable{DnaSequence}, ScoringMatrix?, MsaOptions?)"/> together with the lower bound used and the number of lattice cells evaluated.
    /// </summary>
    [System.Runtime.CompilerServices.OverloadResolutionPriority(1)]
    public static ExactMsaResult MultipleAlignExactWithStatistics(
        IEnumerable<DnaSequence> sequences, ScoringMatrix? scoring = null, MsaOptions? options = null)
    {
        ArgumentNullException.ThrowIfNull(sequences);
        return ExactAlignCore(DnaRows(sequences), scoring ?? SimpleDna, options);
    }

    /// <summary>
    /// <see cref="MultipleAlignExact(IEnumerable{DnaSequence}, ScoringMatrix?, MsaOptions?)"/> for <b>protein</b> sequences (or any alphabet of the
    /// <see cref="SubstitutionMatrix"/>): no other alignment of the sequences has a larger sum-of-pairs score under the matrix and the (possibly
    /// fractional) gap scores — the same k-dimensional dynamic program with the Carrillo–Lipman bound as the DNA overload, with the substitution
    /// matrix replacing match/mismatch. At most 8 sequences (linear) / 6 (affine); k = 2 is <c>GlobalAlign</c> / <c>GlobalAlignAffine</c> with the same
    /// <see cref="MatrixScoring"/>. See <see cref="MultipleAlign(IEnumerable{ProteinSequence}, MatrixScoring, MsaOptions?)"/> for the input contract.
    /// </summary>
    public static MultipleAlignmentResult MultipleAlignExact(
        IEnumerable<ProteinSequence> sequences, MatrixScoring scoring, MsaOptions? options = null) =>
        MultipleAlignExactWithStatistics(sequences, scoring, options).Alignment;

    /// <summary><see cref="MultipleAlignExact(IEnumerable{ProteinSequence}, MatrixScoring, MsaOptions?)"/> for raw sequence strings.</summary>
    public static MultipleAlignmentResult MultipleAlignExact(
        IEnumerable<string> sequences, MatrixScoring scoring, MsaOptions? options = null) =>
        MultipleAlignExactWithStatistics(sequences, scoring, options).Alignment;

    /// <summary><see cref="MultipleAlignExactWithStatistics(IEnumerable{DnaSequence}, ScoringMatrix?, MsaOptions?)"/> for protein sequences.</summary>
    public static ExactMsaResult MultipleAlignExactWithStatistics(
        IEnumerable<ProteinSequence> sequences, MatrixScoring scoring, MsaOptions? options = null)
    {
        var (rows, ps) = PrepareMatrixMsa(sequences, scoring);
        return ScaleExactResult(ExactAlignCore(rows, ps, options), ps);
    }

    /// <summary><see cref="MultipleAlignExactWithStatistics(IEnumerable{DnaSequence}, ScoringMatrix?, MsaOptions?)"/> for raw sequence strings.</summary>
    public static ExactMsaResult MultipleAlignExactWithStatistics(
        IEnumerable<string> sequences, MatrixScoring scoring, MsaOptions? options = null)
    {
        var (rows, ps) = PrepareMatrixMsa(sequences, scoring);
        return ScaleExactResult(ExactAlignCore(rows, ps, options), ps);
    }

    private static ExactMsaResult ScaleExactResult(ExactMsaResult r, PairScoring ps) =>
        ps.Scale == 1 ? r : r with { Alignment = ScaleMsaResult(r.Alignment, ps), ScoreScale = ps.Scale };

    private static ExactMsaResult ExactAlignCore(string[] raw, PairScoring sc, MsaOptions? options)
    {
        using var cancelScope = AlignmentCancellation.Enter(options?.CancellationToken ?? default);
        var gapModel = ResolveGapModel(options);
        var trivial = TrivialMultipleAlignment(raw);
        if (trivial != null)
            return new ExactMsaResult(trivial, 0, 1, 1, 1);

        bool affine = gapModel == MsaGapModel.Affine;

        if (raw.Length == 2)
        {
            // k = 2: the canonical pairwise program (the lattice of the engine below with k = 2 is exactly its matrix).
            var r = affine ? AffineAlignCore(raw[0], raw[1], sc, AlignmentType.Global) : GlobalAlignCore(raw[0], raw[1], sc);
            var rows = new List<string> { r.AlignedSequence1, r.AlignedSequence2 };
            long cells = (long)(raw[0].Length + 1) * (raw[1].Length + 1);
            return new ExactMsaResult(
                new MultipleAlignmentResult(rows.ToArray(), BuildConsensus(rows, rows[0].Length, sc.IsMatrix), r.Score), 0, cells, cells, cells);
        }

        int maxK = affine ? ExactMaxSequencesAffine : ExactMaxSequencesLinear;
        if (raw.Length > maxK)
            throw new ArgumentException(
                $"Exact multiple alignment supports at most {maxK} sequences with {(affine ? "affine" : "linear")} gaps ({raw.Length} given): the dynamic program is exponential in the number of sequences.",
                nameof(raw));

        var opts = options ?? new MsaOptions();
        if (opts.ExactMaxStateCells < 1)
            throw new ArgumentOutOfRangeException(nameof(options), opts.ExactMaxStateCells, "ExactMaxStateCells must be at least 1.");

        int lowerBound = int.MinValue;
        if (opts.ExactCarrilloLipman)
        {
            // Lower bound: the best of the heuristic aligners of this class under the same objective.
            var heuristicOptions = new MsaOptions { GapModel = gapModel };
            foreach (var h in new[]
            {
                StarAlignCore(raw, sc, heuristicOptions),
                ProgressiveAlignCore(raw, sc, heuristicOptions),
                IterativeAlignCore(raw, sc, heuristicOptions, DefaultMaxRefinementIterations),
            })
                lowerBound = Math.Max(lowerBound, ComputeSumOfPairsScore(h.AlignedSequences, sc, gapModel));
        }

        var engine = new ExactMsaEngine(raw, sc, affine, opts.ExactCarrilloLipman, lowerBound, opts.ExactMaxStateCells);
        var aligned = engine.Solve(out int optimum);
        int check = ComputeSumOfPairsScore(aligned, sc, gapModel);
        if (check != optimum)
            throw new InvalidOperationException($"Exact multiple alignment is inconsistent: dynamic-program score {optimum} != sum of pairs {check} of its rows.");

        var alignment = new MultipleAlignmentResult(aligned.ToArray(), BuildConsensus(aligned, aligned[0].Length, sc.IsMatrix), optimum);
        return new ExactMsaResult(alignment, opts.ExactCarrilloLipman ? lowerBound : 0, engine.LatticeCells, engine.VisitedCells, engine.StateValues);
    }

    /// <summary>The k-dimensional dynamic program (see the region comment).</summary>
    private sealed class ExactMsaEngine
    {
        private readonly string[] _s;
        private readonly int _k;
        private readonly int _masks;            // 1 << k
        private readonly PairScoring _sc;
        private readonly bool _affine;
        private readonly bool _prune;
        private readonly long _lower;
        private readonly long _limit;
        private readonly int[] _len;
        private readonly long[] _stride;
        private readonly long[] _maskStride;
        private readonly long _lastCell;

        // pair (a < b) index tables
        private readonly int[] _pairA;
        private readonly int[] _pairB;
        private readonly int _pairs;
        private readonly int[][] _maskPairs;    // for every mask: the pairs with both members in it

        // gap-history states
        private int _states;
        private int[] _next = Array.Empty<int>();       // [state * _masks + mask] -> state
        private int[] _gapCost = Array.Empty<int>();    // [state * _masks + mask] -> summed gap cost of the step
        private byte[] _pairState = Array.Empty<byte>(); // [state * _pairs + pair] -> 0 / 1 (X: a moved last) / 2 (Y: b moved last)

        // pairwise bounds
        private int[][]? _u;                   // [pair][i * (Lb + 1) + j]  best pairwise alignment through (i, j)
        private int[][][]? _g;                 // [pair][state 0..2][...]    best pairwise suffix given the pair's state

        // dynamic-program storage
        private const int ChunkCells = 4096;
        private readonly Dictionary<long, int> _ids = new();
        private readonly List<int[]> _chunks = new();
        private int _count;

        public long LatticeCells { get; }
        public long VisitedCells => _count;
        public long StateValues => (long)_count * _states;

        public ExactMsaEngine(string[] seqs, PairScoring sc, bool affine, bool prune, long lowerBound, long limit)
        {
            _s = seqs;
            _k = seqs.Length;
            _masks = 1 << _k;
            _sc = sc;
            _affine = affine;
            _prune = prune;
            _lower = lowerBound;
            _limit = limit;
            _len = seqs.Select(x => x.Length).ToArray();

            _stride = new long[_k];
            long total = 1;
            try
            {
                checked
                {
                    for (int i = 0; i < _k; i++)
                    {
                        _stride[i] = total;
                        total *= _len[i] + 1;
                    }
                }
            }
            catch (OverflowException)
            {
                throw new ArgumentException("The sequences are too long for the exact multiple alignment lattice (the number of lattice cells overflows).", nameof(seqs));
            }

            LatticeCells = total;
            _lastCell = total - 1;
            _maskStride = new long[_masks];
            for (int m = 1; m < _masks; m++)
                for (int i = 0; i < _k; i++)
                    if ((m >> i & 1) != 0)
                        _maskStride[m] += _stride[i];

            var pa = new List<int>();
            var pb = new List<int>();
            for (int a = 0; a < _k; a++)
                for (int b = a + 1; b < _k; b++)
                {
                    pa.Add(a);
                    pb.Add(b);
                }

            _pairA = pa.ToArray();
            _pairB = pb.ToArray();
            _pairs = _pairA.Length;
            _maskPairs = new int[_masks][];
            for (int m = 0; m < _masks; m++)
            {
                var list = new List<int>();
                for (int p = 0; p < _pairs; p++)
                    if ((m >> _pairA[p] & 1) != 0 && (m >> _pairB[p] & 1) != 0)
                        list.Add(p);
                _maskPairs[m] = list.ToArray();
            }

            BuildStates();
            if (prune)
                BuildPairBounds();
        }

        #region gap-history states

        /// <summary>
        /// States = weak orders of the sequences by the time of their last move, all tied at the origin (every pair in state 0).
        /// A step S puts the moved sequences on top (tied); the pair (a, b) is in state X (1) when a moved after b, Y (2) when b moved
        /// after a, 0 when tied (last projected column residue / residue, or none yet). Linear gaps: one state.
        /// </summary>
        private void BuildStates()
        {
            if (!_affine)
            {
                _states = 1;
                _next = new int[_masks];
                _gapCost = new int[_masks];
                _pairState = new byte[_pairs];
                for (int m = 1; m < _masks; m++)
                {
                    int exactlyOne = 0;
                    for (int p = 0; p < _pairs; p++)
                        if ((m >> _pairA[p] & 1) != (m >> _pairB[p] & 1))
                            exactlyOne++;
                    _gapCost[m] = exactlyOne * _sc.GapExtend;
                }

                return;
            }

            var ids = new Dictionary<int, int>();
            var ranks = new List<int[]>();
            int Intern(int[] r)
            {
                int key = 0;
                for (int i = 0; i < _k; i++) key = key * 8 + r[i];
                if (!ids.TryGetValue(key, out int id))
                {
                    id = ranks.Count;
                    ids[key] = id;
                    ranks.Add(r);
                }

                return id;
            }

            Intern(new int[_k]);
            var next = new List<int[]>();
            var cost = new List<int[]>();
            for (int st = 0; st < ranks.Count; st++)   // ranks grows while we scan
            {
                var r = ranks[st];
                var nx = new int[_masks];
                var gc = new int[_masks];
                for (int m = 1; m < _masks; m++)
                {
                    // dense ranks of the sequences that do not move, then the moved ones on top
                    var distinct = new SortedSet<int>();
                    for (int i = 0; i < _k; i++)
                        if ((m >> i & 1) == 0) distinct.Add(r[i]);
                    var dense = new Dictionary<int, int>();
                    foreach (int v in distinct) dense[v] = dense.Count;
                    var nr = new int[_k];
                    for (int i = 0; i < _k; i++)
                        nr[i] = (m >> i & 1) != 0 ? distinct.Count : dense[r[i]];
                    nx[m] = Intern(nr);

                    int c = 0;
                    for (int p = 0; p < _pairs; p++)
                    {
                        bool inA = (m >> _pairA[p] & 1) != 0, inB = (m >> _pairB[p] & 1) != 0;
                        if (inA == inB) continue;
                        int kind = inA ? 1 : 2;
                        int cur = PairStateOf(r[_pairA[p]], r[_pairB[p]]);
                        c += _sc.GapExtend + (cur != kind ? _sc.GapOpen : 0);
                    }

                    gc[m] = c;
                }

                next.Add(nx);
                cost.Add(gc);
            }

            _states = ranks.Count;
            _next = new int[_states * _masks];
            _gapCost = new int[_states * _masks];
            _pairState = new byte[_states * _pairs];
            for (int st = 0; st < _states; st++)
            {
                Array.Copy(next[st], 0, _next, st * _masks, _masks);
                Array.Copy(cost[st], 0, _gapCost, st * _masks, _masks);
                for (int p = 0; p < _pairs; p++)
                {
                    int ra = ranks[st][_pairA[p]], rb = ranks[st][_pairB[p]];
                    _pairState[st * _pairs + p] = (byte)PairStateOf(ra, rb);
                }
            }
        }

        #endregion

        #region pairwise bounds (Carrillo-Lipman)

        /// <summary>
        /// For every pair the forward Gotoh matrices (M, X, Y) and the backward program G_s(i, j) = best score of a·[i..] / b·[j..] given
        /// that the previous column of the pair was in state s, so that U(i, j) = max_s F_s(i, j) + G_s(i, j) is the best pairwise
        /// alignment through (i, j). Linear gaps = the same program with GapOpen = 0.
        /// </summary>
        private void BuildPairBounds()
        {
            int open = _affine ? _sc.GapOpen : 0;
            int ext = _sc.GapExtend;
            int oe = open + ext;
            _u = new int[_pairs][];
            _g = new int[_pairs][][];
            for (int p = 0; p < _pairs; p++)
            {
                string a = _s[_pairA[p]], b = _s[_pairB[p]];
                int m = a.Length, n = b.Length, w = n + 1;
                var fm = new int[(m + 1) * w];
                var fx = new int[(m + 1) * w];
                var fy = new int[(m + 1) * w];
                for (int i = 0; i <= m; i++)
                    for (int j = 0; j <= n; j++)
                    {
                        int c = i * w + j;
                        if (i == 0 && j == 0) { fm[c] = 0; fx[c] = ExactNegInf; fy[c] = ExactNegInf; continue; }
                        fm[c] = i > 0 && j > 0
                            ? Max3(fm[c - w - 1], fx[c - w - 1], fy[c - w - 1]) + _sc.Sub(a[i - 1], b[j - 1])
                            : ExactNegInf;
                        fx[c] = i > 0 ? Max3(fm[c - w] + oe, fx[c - w] + ext, fy[c - w] + oe) : ExactNegInf;
                        fy[c] = j > 0 ? Max3(fm[c - 1] + oe, fy[c - 1] + ext, fx[c - 1] + oe) : ExactNegInf;
                        if (fm[c] < ExactNegInf) fm[c] = ExactNegInf;
                        if (fx[c] < ExactNegInf) fx[c] = ExactNegInf;
                        if (fy[c] < ExactNegInf) fy[c] = ExactNegInf;
                    }

                var g0 = new int[(m + 1) * w];
                var g1 = new int[(m + 1) * w];
                var g2 = new int[(m + 1) * w];
                var u = new int[(m + 1) * w];
                for (int i = m; i >= 0; i--)
                    for (int j = n; j >= 0; j--)
                    {
                        int c = i * w + j;
                        if (i == m && j == n) { g0[c] = g1[c] = g2[c] = 0; }
                        else
                        {
                            int diag = i < m && j < n ? _sc.Sub(a[i], b[j]) + g0[c + w + 1] : ExactNegInf;
                            int down = i < m ? g1[c + w] : ExactNegInf;    // a[i] against a gap
                            int right = j < n ? g2[c + 1] : ExactNegInf;   // gap against b[j]
                            g0[c] = Max3(diag, down == ExactNegInf ? ExactNegInf : down + oe, right == ExactNegInf ? ExactNegInf : right + oe);
                            g1[c] = Max3(diag, down == ExactNegInf ? ExactNegInf : down + ext, right == ExactNegInf ? ExactNegInf : right + oe);
                            g2[c] = Max3(diag, down == ExactNegInf ? ExactNegInf : down + oe, right == ExactNegInf ? ExactNegInf : right + ext);
                        }

                        u[c] = Math.Max(Math.Max(fm[c] + g0[c], fx[c] + g1[c]), fy[c] + g2[c]);
                    }

                _u[p] = u;
                _g[p] = _affine ? new[] { g0, g1, g2 } : new[] { g0, g0, g0 };
            }
        }

        /// <summary>State of a pair from the move recency ranks of its members: 0 tied, 1 (X) first moved last, 2 (Y) second moved last.</summary>
        private static int PairStateOf(int rankA, int rankB)
        {
            if (rankA == rankB)
                return 0;
            return rankA > rankB ? 1 : 2;
        }

        private static int Max3(int a, int b, int c) => Math.Max(Math.Max(a, b), c);

        #endregion

        #region cells

        private int IdOf(long cell, out bool created)
        {
            if (_ids.TryGetValue(cell, out int id))
            {
                created = false;
                return id;
            }

            if ((long)(_count + 1) * _states > _limit)
                throw new ArgumentException(
                    string.Create(System.Globalization.CultureInfo.InvariantCulture, $"Exact multiple alignment exceeds the complexity limit: more than {_limit:N0} dynamic-program values ") +
                    string.Create(System.Globalization.CultureInfo.InvariantCulture, $"({_count:N0} cells x {_states} gap states; lattice {LatticeCells:N0} cells). Raise MsaOptions.ExactMaxStateCells, use more similar / shorter sequences, ") +
                    "or use a heuristic aligner (MultipleAlign / MultipleAlignProgressive / MultipleAlignIterative).");
            id = _count++;
            if (id / ChunkCells >= _chunks.Count)
                _chunks.Add(new int[(long)ChunkCells * _states]);
            Array.Fill(_chunks[id / ChunkCells], ExactNegInf, id % ChunkCells * _states, _states);
            _ids[cell] = id;
            created = true;
            return id;
        }

        private (int[] Chunk, int Offset) Slot(int id) => (_chunks[id / ChunkCells], id % ChunkCells * _states);

        private void Decode(long cell, int[] x)
        {
            for (int i = _k - 1; i >= 0; i--)
            {
                x[i] = (int)(cell / _stride[i]);
                cell %= _stride[i];
            }
        }

        /// <summary>Sum of the substitution scores of the pairs inside the step, residues taken at the positions x.</summary>
        private int Substitution(int[] x, int mask)
        {
            int total = 0;
            foreach (int p in _maskPairs[mask])
                total += _sc.Sub(_s[_pairA[p]][x[_pairA[p]]], _s[_pairB[p]][x[_pairB[p]]]);
            return total;
        }

        private long UpperBoundSum(int[] y)
        {
            long total = 0;
            for (int p = 0; p < _pairs; p++)
                total += _u![p][y[_pairA[p]] * (_len[_pairB[p]] + 1) + y[_pairB[p]]];
            return total;
        }

        /// <summary>Best pairwise suffixes of every pair at y given the gap state of the history state.</summary>
        private long SuffixBound(int[] y, int state)
        {
            long total = 0;
            for (int p = 0; p < _pairs; p++)
                total += _g![p][_pairState[state * _pairs + p]][y[_pairA[p]] * (_len[_pairB[p]] + 1) + y[_pairB[p]]];
            return total;
        }

        #endregion

        #region solve

        public List<string> Solve(out int optimum)
        {
            var queue = new PriorityQueue<long, long>();
            int originId = IdOf(0, out _);
            {
                var (chunk, off) = Slot(originId);
                chunk[off] = 0;   // state 0 = all pairs tied
            }

            queue.Enqueue(0, 0);
            var x = new int[_k];
            var y = new int[_k];
            while (queue.TryDequeue(out long cell, out _))
            {
                if (cell == _lastCell)
                    break;
                Decode(cell, x);
                int id = _ids[cell];
                var (chunk, off) = Slot(id);
                for (int mask = 1; mask < _masks; mask++)
                {
                    bool inside = true;
                    for (int i = 0; i < _k && inside; i++)
                        if ((mask >> i & 1) != 0 && x[i] >= _len[i])
                            inside = false;
                    if (!inside)
                        continue;

                    for (int i = 0; i < _k; i++)
                        y[i] = x[i] + (mask >> i & 1);
                    if (_prune && UpperBoundSum(y) < _lower)
                        continue;

                    long target = cell + _maskStride[mask];
                    int subst = Substitution(x, mask);
                    int tid = -1;
                    for (int st = 0; st < _states; st++)
                    {
                        int v = chunk[off + st];
                        if (v <= ExactNegInf)
                            continue;
                        int ns = _next[st * _masks + mask];
                        int cand = v + _gapCost[st * _masks + mask] + subst;
                        if (_prune && cand + SuffixBound(y, ns) < _lower)
                            continue;
                        if (tid < 0)
                        {
                            tid = IdOf(target, out bool created);
                            if (created)
                                queue.Enqueue(target, target);
                            // the previous chunk reference stays valid: chunks are only appended
                        }

                        var (tc, to) = Slot(tid);
                        if (cand > tc[to + ns])
                            tc[to + ns] = cand;
                    }
                }
            }

            if (!_ids.TryGetValue(_lastCell, out int lastId))
                throw new InvalidOperationException("Exact multiple alignment: the far corner was not reached (internal error: the lower bound exceeds the optimum).");

            var (lc, lo) = Slot(lastId);
            int best = ExactNegInf, bestState = -1;
            for (int st = 0; st < _states; st++)
                if (lc[lo + st] > best) { best = lc[lo + st]; bestState = st; }
            if (bestState < 0)
                throw new InvalidOperationException("Exact multiple alignment: the far corner holds no value (internal error).");

            optimum = best;
            return Traceback(bestState, best);
        }

        /// <summary>
        /// Walks back from the far corner: at every (cell, state) the first predecessor (largest step first) whose stored value plus the
        /// step cost reproduces the value; the stored values of all surviving cells are kept, so no pointers are needed.
        /// </summary>
        private List<string> Traceback(int state, int value)
        {
            var steps = new List<int>();
            long cell = _lastCell;
            var x = new int[_k];
            var y = new int[_k];
            while (cell != 0)
            {
                Decode(cell, x);
                bool found = false;
                for (int mask = _masks - 1; mask >= 1 && !found; mask--)
                {
                    bool inside = true;
                    for (int i = 0; i < _k && inside; i++)
                        if ((mask >> i & 1) != 0 && x[i] < 1)
                            inside = false;
                    if (!inside)
                        continue;

                    long pcell = cell - _maskStride[mask];
                    if (!_ids.TryGetValue(pcell, out int pid))
                        continue;
                    for (int i = 0; i < _k; i++)
                        y[i] = x[i] - (mask >> i & 1);
                    int subst = Substitution(y, mask);
                    var (pc, po) = Slot(pid);
                    for (int ps = 0; ps < _states; ps++)
                    {
                        int pv = pc[po + ps];
                        if (pv <= ExactNegInf || _next[ps * _masks + mask] != state)
                            continue;
                        if (pv + _gapCost[ps * _masks + mask] + subst != value)
                            continue;
                        steps.Add(mask);
                        cell = pcell;
                        state = ps;
                        value = pv;
                        found = true;
                        break;
                    }
                }

                if (!found)
                    throw new InvalidOperationException("Exact multiple alignment: traceback failed (internal error).");
            }

            steps.Reverse();
            var rows = Enumerable.Range(0, _k).Select(_ => new System.Text.StringBuilder(steps.Count)).ToArray();
            var pos = new int[_k];
            foreach (int mask in steps)
                for (int i = 0; i < _k; i++)
                    if ((mask >> i & 1) != 0)
                        rows[i].Append(_s[i][pos[i]++]);
                    else
                        rows[i].Append('-');
            return rows.Select(r => r.ToString()).ToList();
        }

        #endregion
    }

    #endregion
}
