namespace Seqeron.Genomics.Alignment;

// =====================================================================================================================
// Co-optimal enumeration for the GLOBAL (Needleman-Wunsch) and ends-free (semi-global) problem (audit A1-7b).
//
// Biopython's PairwiseAligner.align(a, b) iterates every optimal alignment of the Needleman-Wunsch / Gotoh matrix; this file
// enumerates the same set (linear and affine gaps, every end-gap policy of EndGaps, scalar and substitution-matrix scoring)
// from the SAME filled matrices as the single-alignment aligners: LinearFillMatrix / AffineInitBorders / AffineFill /
// EndCellCandidates (the code that LinearAlignCore / AffineAlignCore run), only the traceback differs: a depth-first walk
// over ALL predecessors that reproduce the value of a cell instead of the first one.
//
// Definition of the set (exact, checked against exhaustive enumeration of every column sequence and Biopython 1.88):
//   an alignment = a gapped column sequence; its score = sum of substitution columns + gap costs (linear: GapExtend per gap
//   position; affine: GapOpen + k·GapExtend per maximal run of gap columns of the SAME row — a gap in row 1 directly followed
//   by a gap in row 2 are two runs, Gotoh/Flouri model of GlobalAlignAffine) where leading/trailing gap columns of a FREE end
//   (EndGaps) cost 0; the set = all column sequences of maximum score, each listed once (distinct gapped rows).
//
// Walk (backwards from an end cell to a start cell):
//   * End cells = EndCellCandidates (m, n; or the last row / column when trailing ends are free) whose value equals the optimum
//     (affine: any of M/X/Y with that value); the order is the one of LinearAlignCore / AffineAlignCore (last row ascending j,
//     then last column ascending i), so the first alignment equals GlobalAlign / GlobalAlignAffine / SemiGlobalAlign* (same
//     tie order diagonal > up > left, affine M > X > Y state with Y's own state before X).
//   * Linear: from (i, j) go to every predecessor that reproduces F(i, j): diagonal (a_i, b_j), up (a_i against a gap), left.
//     Affine: from (i, j, state) go to every predecessor state that reproduces the value (M: any state at (i-1, j-1) + s; X: M/X/Y at
//     (i-1, j) + o + e / e / o + e; Y symmetric) — a never-unreachable (−∞) state is not a predecessor.
//   * Start: the origin (0, 0), or, when the matching start is free, the first cell reached in column 0 (free seq-1 prefix,
//     EndGaps.FreeStart1) / row 0 (free seq-2 prefix, FreeStart2); the free flank is shown against gaps exactly as the single
//     aligners do. Walking on along a free border would re-list the same rows (and is not possible: the border gap states are −∞).
//   * Each alignment is listed once (distinct gapped rows): with zero gap costs a free trailing flank is also a run of zero-cost
//     gap columns of the path (the same rows end at several cells), so the walk keeps a set of the rows already listed. The first
//     element is always the single aligner's traceback (same end cell, same tie order). An empty side has the single all-gap alignment.
//   * Coordinates are the consumed intervals of the rows (the first column decides which side's leading flank is free, the last
//     column the trailing one), a function of the rows alone.
//
// The number of co-optimal alignments grows exponentially with ties (tandem repeats, long gap-placement ties), so the list is
// capped (maxAlignments); O(m·n) memory as the Waterman-Eggert declumper (TopNMaxCells).
// =====================================================================================================================
public static partial class SequenceAligner
{
    #region Co-optimal global / ends-free alignments (Biopython PairwiseAligner.align set)

    /// <summary>
    /// All co-optimal <b>global</b> (Needleman-Wunsch, linear gaps) alignments of <paramref name="sequence1"/> and
    /// <paramref name="sequence2"/> — the set that Biopython's <c>PairwiseAligner(mode='global').align(a, b)</c> iterates — capped at
    /// <paramref name="maxAlignments"/>.
    /// </summary>
    /// <param name="sequence1">First sequence.</param>
    /// <param name="sequence2">Second sequence.</param>
    /// <param name="scoring">Scoring (default <see cref="SimpleDna"/>); only <see cref="ScoringMatrix.GapExtend"/> is used.</param>
    /// <param name="maxAlignments">Cap on the number of alignments returned (≥ 1; the number of co-optimal alignments grows
    /// exponentially with tandem repeats and long gap-placement ties).</param>
    /// <returns>Distinct gapped-row pairs of the optimal score, in a deterministic order: the first element is
    /// <see cref="GlobalAlign(DnaSequence, DnaSequence, ScoringMatrix?)"/>'s result (diagonal &gt; up &gt; left at every branch of the
    /// depth-first walk from (m, n)); an empty side gives the single all-gap alignment.</returns>
    /// <remarks>
    /// The matrix is the one of <see cref="GlobalAlign(DnaSequence, DnaSequence, ScoringMatrix?)"/>; an alignment is a path of the
    /// traceback graph from (m, n) to the origin, and every path is listed (two paths = two column sequences = two different gapped
    /// rows). Needs the full O(m·n) matrix (<see cref="TopNMaxCells"/> cells, 4 bytes per cell; larger inputs throw
    /// <see cref="ArgumentException"/> — use <see cref="GlobalAlign(DnaSequence, DnaSequence, ScoringMatrix?)"/>, which falls back to
    /// linear space, when one optimal alignment is enough).
    /// </remarks>
    /// <exception cref="ArgumentNullException">A sequence or <paramref name="scoring"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="maxAlignments"/> &lt; 1.</exception>
    /// <exception cref="ArgumentException">The matrix would exceed <see cref="TopNMaxCells"/> cells.</exception>
    [System.Runtime.CompilerServices.OverloadResolutionPriority(1)]
    public static IReadOnlyList<AlignmentResult> GlobalAlignCoOptimal(
        DnaSequence sequence1, DnaSequence sequence2, ScoringMatrix? scoring = null, int maxAlignments = 1000)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        return GlobalCoOptimalCore(sequence1.Sequence, sequence2.Sequence, scoring ?? SimpleDna, affine: false, EndGaps.None, maxAlignments);
    }

    /// <summary><see cref="GlobalAlignCoOptimal(DnaSequence, DnaSequence, ScoringMatrix?, int)"/> on raw strings (upper-cased internally).</summary>
    [System.Runtime.CompilerServices.OverloadResolutionPriority(1)]
    public static IReadOnlyList<AlignmentResult> GlobalAlignCoOptimal(
        string sequence1, string sequence2, ScoringMatrix? scoring = null, int maxAlignments = 1000)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        return GlobalCoOptimalCore(sequence1.ToUpperInvariant(), sequence2.ToUpperInvariant(), scoring ?? SimpleDna, affine: false, EndGaps.None, maxAlignments);
    }

    /// <summary>
    /// All co-optimal global alignments under the <b>affine</b> gap model (a gap of length k scores
    /// <see cref="ScoringMatrix.GapOpen"/> + k·<see cref="ScoringMatrix.GapExtend"/>; Biopython <c>mode='global'</c> with
    /// <c>open_gap_score = GapOpen + GapExtend</c>, <c>extend_gap_score = GapExtend</c>); the first element is
    /// <see cref="GlobalAlignAffine(DnaSequence, DnaSequence, ScoringMatrix?)"/>'s result; see
    /// <see cref="GlobalAlignCoOptimal(DnaSequence, DnaSequence, ScoringMatrix?, int)"/>.
    /// </summary>
    [System.Runtime.CompilerServices.OverloadResolutionPriority(1)]
    public static IReadOnlyList<AlignmentResult> GlobalAlignAffineCoOptimal(
        DnaSequence sequence1, DnaSequence sequence2, ScoringMatrix? scoring = null, int maxAlignments = 1000)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        return GlobalCoOptimalCore(sequence1.Sequence, sequence2.Sequence, scoring ?? SimpleDna, affine: true, EndGaps.None, maxAlignments);
    }

    /// <summary><see cref="GlobalAlignAffineCoOptimal(DnaSequence, DnaSequence, ScoringMatrix?, int)"/> on raw strings.</summary>
    [System.Runtime.CompilerServices.OverloadResolutionPriority(1)]
    public static IReadOnlyList<AlignmentResult> GlobalAlignAffineCoOptimal(
        string sequence1, string sequence2, ScoringMatrix? scoring = null, int maxAlignments = 1000)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        return GlobalCoOptimalCore(sequence1.ToUpperInvariant(), sequence2.ToUpperInvariant(), scoring ?? SimpleDna, affine: true, EndGaps.None, maxAlignments);
    }

    /// <summary>All co-optimal global alignments with a substitution matrix and linear gap scores
    /// (<see cref="GlobalAlign(string, string, MatrixScoring)"/>; <see cref="AlignmentResult.ExactScore"/> = real score).</summary>
    public static IReadOnlyList<AlignmentResult> GlobalAlignCoOptimal(
        string sequence1, string sequence2, MatrixScoring scoring, int maxAlignments = 1000) =>
        GlobalCoOptimalWithMatrix(sequence1, sequence2, scoring, affine: false, EndGaps.None, maxAlignments);

    /// <summary><see cref="GlobalAlignCoOptimal(string, string, MatrixScoring, int)"/> for protein sequences.</summary>
    public static IReadOnlyList<AlignmentResult> GlobalAlignCoOptimal(
        ProteinSequence sequence1, ProteinSequence sequence2, MatrixScoring scoring, int maxAlignments = 1000)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        return GlobalCoOptimalWithMatrix(sequence1.Sequence, sequence2.Sequence, scoring, affine: false, EndGaps.None, maxAlignments);
    }

    /// <summary><see cref="GlobalAlignCoOptimal(string, string, MatrixScoring, int)"/> for DNA sequences.</summary>
    public static IReadOnlyList<AlignmentResult> GlobalAlignCoOptimal(
        DnaSequence sequence1, DnaSequence sequence2, MatrixScoring scoring, int maxAlignments = 1000)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        return GlobalCoOptimalWithMatrix(sequence1.Sequence, sequence2.Sequence, scoring, affine: false, EndGaps.None, maxAlignments);
    }

    /// <summary>All co-optimal affine-gap global alignments with a substitution matrix
    /// (<see cref="GlobalAlignAffine(string, string, MatrixScoring)"/>).</summary>
    public static IReadOnlyList<AlignmentResult> GlobalAlignAffineCoOptimal(
        string sequence1, string sequence2, MatrixScoring scoring, int maxAlignments = 1000) =>
        GlobalCoOptimalWithMatrix(sequence1, sequence2, scoring, affine: true, EndGaps.None, maxAlignments);

    /// <summary><see cref="GlobalAlignAffineCoOptimal(string, string, MatrixScoring, int)"/> for protein sequences.</summary>
    public static IReadOnlyList<AlignmentResult> GlobalAlignAffineCoOptimal(
        ProteinSequence sequence1, ProteinSequence sequence2, MatrixScoring scoring, int maxAlignments = 1000)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        return GlobalCoOptimalWithMatrix(sequence1.Sequence, sequence2.Sequence, scoring, affine: true, EndGaps.None, maxAlignments);
    }

    /// <summary><see cref="GlobalAlignAffineCoOptimal(string, string, MatrixScoring, int)"/> for DNA sequences.</summary>
    public static IReadOnlyList<AlignmentResult> GlobalAlignAffineCoOptimal(
        DnaSequence sequence1, DnaSequence sequence2, MatrixScoring scoring, int maxAlignments = 1000)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        return GlobalCoOptimalWithMatrix(sequence1.Sequence, sequence2.Sequence, scoring, affine: true, EndGaps.None, maxAlignments);
    }

    /// <summary>
    /// All co-optimal <b>ends-free</b> (semi-global, any <see cref="EndGaps"/> policy) alignments with linear gaps — Biopython
    /// <c>PairwiseAligner(mode='global')</c> with <c>left/right_insertion/deletion_score = 0</c> for the free ends; the first element
    /// is <see cref="SemiGlobalAlign(string, string, ScoringMatrix?, EndGaps)"/>'s result. <see cref="EndGaps.None"/> = the global
    /// problem (results tagged Global). A free flank is shown against '-' in the rows; coordinates are the consumed intervals
    /// (as of <see cref="SemiGlobalAlign(string, string, ScoringMatrix?, EndGaps)"/>) of each alignment.
    /// </summary>
    /// <remarks>An alignment that merely re-labels a free flank (same rows) is listed once; see
    /// <see cref="GlobalAlignCoOptimal(DnaSequence, DnaSequence, ScoringMatrix?, int)"/> for the model, cap and memory bound.</remarks>
    [System.Runtime.CompilerServices.OverloadResolutionPriority(1)]
    public static IReadOnlyList<AlignmentResult> SemiGlobalAlignCoOptimal(
        DnaSequence sequence1, DnaSequence sequence2, ScoringMatrix? scoring, EndGaps endGaps, int maxAlignments = 1000)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        return GlobalCoOptimalCore(sequence1.Sequence, sequence2.Sequence, scoring ?? SimpleDna, affine: false, endGaps, maxAlignments);
    }

    /// <summary><see cref="SemiGlobalAlignCoOptimal(DnaSequence, DnaSequence, ScoringMatrix?, EndGaps, int)"/> on raw strings.</summary>
    [System.Runtime.CompilerServices.OverloadResolutionPriority(1)]
    public static IReadOnlyList<AlignmentResult> SemiGlobalAlignCoOptimal(
        string sequence1, string sequence2, ScoringMatrix? scoring, EndGaps endGaps, int maxAlignments = 1000)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        return GlobalCoOptimalCore(sequence1.ToUpperInvariant(), sequence2.ToUpperInvariant(), scoring ?? SimpleDna, affine: false, endGaps, maxAlignments);
    }

    /// <summary>All co-optimal ends-free alignments with <b>affine</b> gaps (GapOpen + k·GapExtend; first element =
    /// <see cref="SemiGlobalAlignAffine(string, string, ScoringMatrix?, EndGaps)"/>).</summary>
    [System.Runtime.CompilerServices.OverloadResolutionPriority(1)]
    public static IReadOnlyList<AlignmentResult> SemiGlobalAlignAffineCoOptimal(
        DnaSequence sequence1, DnaSequence sequence2, ScoringMatrix? scoring, EndGaps endGaps, int maxAlignments = 1000)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        return GlobalCoOptimalCore(sequence1.Sequence, sequence2.Sequence, scoring ?? SimpleDna, affine: true, endGaps, maxAlignments);
    }

    /// <summary><see cref="SemiGlobalAlignAffineCoOptimal(DnaSequence, DnaSequence, ScoringMatrix?, EndGaps, int)"/> on raw strings.</summary>
    [System.Runtime.CompilerServices.OverloadResolutionPriority(1)]
    public static IReadOnlyList<AlignmentResult> SemiGlobalAlignAffineCoOptimal(
        string sequence1, string sequence2, ScoringMatrix? scoring, EndGaps endGaps, int maxAlignments = 1000)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        return GlobalCoOptimalCore(sequence1.ToUpperInvariant(), sequence2.ToUpperInvariant(), scoring ?? SimpleDna, affine: true, endGaps, maxAlignments);
    }

    /// <summary>All co-optimal ends-free alignments with a substitution matrix and linear gap scores.</summary>
    public static IReadOnlyList<AlignmentResult> SemiGlobalAlignCoOptimal(
        string sequence1, string sequence2, MatrixScoring scoring, EndGaps endGaps, int maxAlignments = 1000) =>
        GlobalCoOptimalWithMatrix(sequence1, sequence2, scoring, affine: false, endGaps, maxAlignments);

    /// <summary><see cref="SemiGlobalAlignCoOptimal(string, string, MatrixScoring, EndGaps, int)"/> for protein sequences.</summary>
    public static IReadOnlyList<AlignmentResult> SemiGlobalAlignCoOptimal(
        ProteinSequence sequence1, ProteinSequence sequence2, MatrixScoring scoring, EndGaps endGaps, int maxAlignments = 1000)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        return GlobalCoOptimalWithMatrix(sequence1.Sequence, sequence2.Sequence, scoring, affine: false, endGaps, maxAlignments);
    }

    /// <summary><see cref="SemiGlobalAlignCoOptimal(string, string, MatrixScoring, EndGaps, int)"/> for DNA sequences.</summary>
    public static IReadOnlyList<AlignmentResult> SemiGlobalAlignCoOptimal(
        DnaSequence sequence1, DnaSequence sequence2, MatrixScoring scoring, EndGaps endGaps, int maxAlignments = 1000)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        return GlobalCoOptimalWithMatrix(sequence1.Sequence, sequence2.Sequence, scoring, affine: false, endGaps, maxAlignments);
    }

    /// <summary>All co-optimal ends-free alignments with a substitution matrix and affine gap scores.</summary>
    public static IReadOnlyList<AlignmentResult> SemiGlobalAlignAffineCoOptimal(
        string sequence1, string sequence2, MatrixScoring scoring, EndGaps endGaps, int maxAlignments = 1000) =>
        GlobalCoOptimalWithMatrix(sequence1, sequence2, scoring, affine: true, endGaps, maxAlignments);

    /// <summary><see cref="SemiGlobalAlignAffineCoOptimal(string, string, MatrixScoring, EndGaps, int)"/> for protein sequences.</summary>
    public static IReadOnlyList<AlignmentResult> SemiGlobalAlignAffineCoOptimal(
        ProteinSequence sequence1, ProteinSequence sequence2, MatrixScoring scoring, EndGaps endGaps, int maxAlignments = 1000)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        return GlobalCoOptimalWithMatrix(sequence1.Sequence, sequence2.Sequence, scoring, affine: true, endGaps, maxAlignments);
    }

    /// <summary><see cref="SemiGlobalAlignAffineCoOptimal(string, string, MatrixScoring, EndGaps, int)"/> for DNA sequences.</summary>
    public static IReadOnlyList<AlignmentResult> SemiGlobalAlignAffineCoOptimal(
        DnaSequence sequence1, DnaSequence sequence2, MatrixScoring scoring, EndGaps endGaps, int maxAlignments = 1000)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        return GlobalCoOptimalWithMatrix(sequence1.Sequence, sequence2.Sequence, scoring, affine: true, endGaps, maxAlignments);
    }

    private static IReadOnlyList<AlignmentResult> GlobalCoOptimalWithMatrix(
        string sequence1, string sequence2, MatrixScoring scoring, bool affine, EndGaps ends, int maxAlignments)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        ArgumentNullException.ThrowIfNull(scoring);

        string a = sequence1.ToUpperInvariant();
        string b = sequence2.ToUpperInvariant();
        var ps = PairScoring.Create(scoring, a, b);
        var raw = GlobalCoOptimalCore(a, b, ps, affine, ends, maxAlignments);
        var list = new List<AlignmentResult>(raw.Count);
        foreach (var r in raw)
            list.Add(WithExactScore(r, ps));
        return list;
    }

    /// <summary>The result with the scaled integer score converted back to the real score (<see cref="AlignmentResult.ExactScore"/>).</summary>
    private static AlignmentResult WithExactScore(AlignmentResult r, PairScoring ps)
    {
        double exact = (double)r.Score / ps.Scale;
        return r with { Score = (int)Math.Round(exact, MidpointRounding.AwayFromZero), ExactScore = exact };
    }

    internal static IReadOnlyList<AlignmentResult> GlobalCoOptimalCore(
        string a, string b, PairScoring scoring, bool affine, EndGaps ends, int maxAlignments)
    {
        if (maxAlignments < 1)
            throw new ArgumentOutOfRangeException(nameof(maxAlignments), maxAlignments, "maxAlignments must be at least 1.");
        if ((ends & ~EndGaps.All) != 0)
            throw new ArgumentOutOfRangeException(nameof(ends), ends, "Unknown end-gap policy.");
        RequireTopNFits(a.Length, b.Length, nameof(a));
        scoring.EnsureRange(a.Length, b.Length);

        // An empty side has exactly one alignment (everything against gaps): the single aligner's result.
        if (a.Length == 0 || b.Length == 0)
            return new[] { affine ? AffineAlignCore(a, b, scoring, ResultAlignmentType(ends), ends) : LinearAlignCore(a, b, scoring, ends) };

        var results = new List<AlignmentResult>();
        new GlobalTraceEnumerator(a, b, scoring, affine, ends).Enumerate(maxAlignments, results);
        return results;
    }

    /// <summary>Type argument of <see cref="AffineAlignCore"/> for an end-gap policy (Global without free ends, else SemiGlobal).</summary>
    private static AlignmentType ResultAlignmentType(EndGaps ends) => ResultType(local: false, ends);

    /// <summary>
    /// The filled global / ends-free matrices (linear: F; affine: M, X, Y — the very matrices of <see cref="LinearAlignCore"/> /
    /// <see cref="AffineAlignCore"/>) and the depth-first enumeration of every optimal traceback path.
    /// </summary>
    private sealed class GlobalTraceEnumerator
    {
        private readonly string _a;
        private readonly string _b;
        private readonly PairScoring _sc;
        private readonly bool _affine;
        private readonly EndGaps _ends;
        private readonly int _m;
        private readonly int _n;
        private readonly bool _freeStart1;
        private readonly bool _freeStart2;
        private readonly bool _freeEnd1;
        private readonly bool _freeEnd2;
        private readonly int[,]? _f;
        private readonly int[,]? _mm;
        private readonly int[,]? _xx;
        private readonly int[,]? _yy;
        private readonly HashSet<(string, string)> _seen = new();

        public GlobalTraceEnumerator(string a, string b, PairScoring sc, bool affine, EndGaps ends)
        {
            _a = a;
            _b = b;
            _sc = sc;
            _affine = affine;
            _ends = ends;
            _m = a.Length;
            _n = b.Length;
            _freeStart1 = ends.HasFlag(EndGaps.FreeStart1);
            _freeStart2 = ends.HasFlag(EndGaps.FreeStart2);
            _freeEnd1 = ends.HasFlag(EndGaps.FreeEnd1);
            _freeEnd2 = ends.HasFlag(EndGaps.FreeEnd2);
            if (affine)
            {
                (_mm, _xx, _yy) = AffineNewMatrices(_m, _n);
                AffineInitBorders(_mm, _xx, _yy, sc, local: false, ends, StateMatch);
                AffineFill(_mm, _xx, _yy, a, b, sc, local: false);
            }
            else
            {
                _f = LinearFillMatrix(a, b, sc, ends);
            }
        }

        private int Value(int i, int j, byte state) =>
            !_affine ? _f![i, j] : state switch
            {
                StateMatch => _mm![i, j],
                StateGapInSeq2 => _xx![i, j],
                _ => _yy![i, j],
            };

        private int Best(int i, int j) => _affine ? Max3(_mm![i, j], _xx![i, j], _yy![i, j]) : _f![i, j];

        public void Enumerate(int cap, List<AlignmentResult> output)
        {
            // Optimum and the end cells that attain it, in the order of the single aligners (first = their end cell).
            int best = int.MinValue;
            var cells = new List<(int I, int J)>();
            foreach (var (ci, cj) in EndCellCandidates(_m, _n, _ends))
            {
                int v = Best(ci, cj);
                if (v > best) { best = v; cells.Clear(); cells.Add((ci, cj)); }
                else if (v == best && !cells.Contains((ci, cj))) cells.Add((ci, cj));
            }

            var steps = new List<byte>();   // columns from the end backwards: 0 = pair, 1 = a_i vs gap, 2 = gap vs b_j
            foreach (var (ei, ej) in cells)
            {
                if (output.Count >= cap)
                    return;
                if (!_affine)
                {
                    Walk(ei, ej, StateMatch, best, steps, cap, output);
                    continue;
                }

                // Affine: every state of the end cell that attains the optimum, M > X > Y.
                foreach (byte state in new[] { StateMatch, StateGapInSeq2, StateGapInSeq1 })
                {
                    if (output.Count >= cap)
                        return;
                    if (Value(ei, ej, state) != best)
                        continue;
                    Walk(ei, ej, state, best, steps, cap, output);
                }
            }
        }

        /// <summary>The start of the walk is reached: the origin, or the first cell in the free border (column 0 when the
        /// sequence-1 start is free, row 0 when the sequence-2 start is free).</summary>
        private bool IsStart(int i, int j) =>
            (i == 0 && (j == 0 || _freeStart2)) || (j == 0 && _freeStart1);

        // Iterative depth-first walk over the traceback graph from an end node to every start node. A frame (i, j, state, next)
        // resumes with candidate `next`; steps.Count = number of columns between the end cell and the node being processed.
        private void Walk(int endI, int endJ, byte endState, int best, List<byte> steps, int cap, List<AlignmentResult> output)
        {
            steps.Clear();
            var stack = new Stack<(int I, int J, byte State, int Next)>();
            stack.Push((endI, endJ, endState, 0));
            int visited = 0;
            while (stack.Count > 0 && output.Count < cap)
            {
                if ((++visited & 0x3FF) == 0)
                    AlignmentCancellation.ThrowIfRequested();
                var (i, j, state, next) = stack.Pop();
                if (IsStart(i, j))
                {
                    var built = Build(endI, endJ, i, j, best, steps);
                    if (_seen.Add((built.AlignedSequence1, built.AlignedSequence2)))
                        output.Add(built);
                    if (steps.Count > 0) steps.RemoveAt(steps.Count - 1);
                    continue;
                }

                bool descended = false;
                for (int cand = next; cand < 3 && !descended; cand++)
                {
                    if (!Predecessor(i, j, state, cand, out int pi, out int pj, out byte pstate, out byte kind))
                        continue;
                    stack.Push((i, j, state, cand + 1));
                    steps.Add(kind);
                    stack.Push((pi, pj, pstate, 0));
                    descended = true;
                }

                if (!descended && steps.Count > 0)
                    steps.RemoveAt(steps.Count - 1);   // node exhausted: leave the column that led to it
            }
        }

        /// <summary>
        /// Candidate predecessor number <paramref name="cand"/> (diagonal / M, up / X, left / Y — the fixed tie order of the single
        /// aligners; Y tries its own state before X) of the node (i, j, state) when it reproduces the node's value.
        /// </summary>
        private bool Predecessor(int i, int j, byte state, int cand, out int pi, out int pj, out byte pstate, out byte kind)
        {
            pi = i; pj = j; pstate = StateMatch; kind = 0;
            if (!_affine)
            {
                var f = _f!;
                int v = f[i, j];
                switch (cand)
                {
                    case 0:
                        if (i == 0 || j == 0 || v != f[i - 1, j - 1] + _sc.Sub(_a[i - 1], _b[j - 1])) return false;
                        pi = i - 1; pj = j - 1; kind = 0; return true;
                    case 1:
                        if (i == 0 || v != f[i - 1, j] + _sc.XExt(j, _n)) return false;
                        pi = i - 1; kind = 1; return true;
                    default:
                        if (j == 0 || v != f[i, j - 1] + _sc.YExt(i, _m)) return false;
                        pj = j - 1; kind = 2; return true;
                }
            }

            switch (state)
            {
                case StateMatch:
                {
                    if (i == 0 || j == 0) return false;
                    int prev = _mm![i, j] - _sc.Sub(_a[i - 1], _b[j - 1]);
                    pi = i - 1; pj = j - 1; kind = 0;
                    pstate = cand switch { 0 => StateMatch, 1 => StateGapInSeq2, _ => StateGapInSeq1 };
                    int pv = Value(pi, pj, pstate);
                    return pv > AffineNegInf && pv == prev;
                }
                case StateGapInSeq2:
                {
                    if (i == 0) return false;
                    int cur = _xx![i, j];
                    pi = i - 1; kind = 1;
                    pstate = cand switch { 0 => StateMatch, 1 => StateGapInSeq2, _ => StateGapInSeq1 };
                    int pv = Value(pi, pj, pstate);
                    int xe = _sc.XExt(j, _n);
                    return pv > AffineNegInf && pv + (pstate == StateGapInSeq2 ? xe : _sc.XOpen(j, _n) + xe) == cur;
                }
                default:
                {
                    if (j == 0) return false;
                    int cur = _yy![i, j];
                    pj = j - 1; kind = 2;
                    pstate = cand switch { 0 => StateMatch, 1 => StateGapInSeq1, _ => StateGapInSeq2 };
                    int pv = Value(pi, pj, pstate);
                    int ye = _sc.YExt(i, _m);
                    return pv > AffineNegInf && pv + (pstate == StateGapInSeq1 ? ye : _sc.YOpen(i, _m) + ye) == cur;
                }
            }
        }

        private AlignmentResult Build(int endI, int endJ, int startI, int startJ, int score, List<byte> steps)
        {
            var r1 = new List<char>(_m + _n);
            var r2 = new List<char>(_m + _n);

            // Free leading flank (emitted first): sequence-2 residues against gaps (row 0) or sequence-1 residues (column 0).
            for (int k = 0; k < startJ && startI == 0; k++) { r1.Add('-'); r2.Add(_b[k]); }
            for (int k = 0; k < startI && startJ == 0; k++) { r1.Add(_a[k]); r2.Add('-'); }

            // The path: the steps run from the end cell backwards.
            int i = endI, j = endJ;
            var pathCols = new (char, char)[steps.Count];
            for (int k = 0; k < steps.Count; k++)
            {
                switch (steps[k])
                {
                    case 0: pathCols[steps.Count - 1 - k] = (_a[i - 1], _b[j - 1]); i--; j--; break;
                    case 1: pathCols[steps.Count - 1 - k] = (_a[i - 1], '-'); i--; break;
                    default: pathCols[steps.Count - 1 - k] = ('-', _b[j - 1]); j--; break;
                }
            }
            foreach (var (x, y) in pathCols) { r1.Add(x); r2.Add(y); }

            // Free trailing flank.
            for (int k = endJ; k < _n; k++) { r1.Add('-'); r2.Add(_b[k]); }
            for (int k = endI; k < _m; k++) { r1.Add(_a[k]); r2.Add('-'); }

            string row1 = new string(r1.ToArray()), row2 = new string(r2.ToArray());
            var (s1, s2, e1, e2) = ConsumedInterval(row1, row2);
            return new AlignmentResult(
                AlignedSequence1: row1,
                AlignedSequence2: row2,
                Score: score,
                AlignmentType: ResultType(local: false, _ends),
                StartPosition1: s1,
                StartPosition2: s2,
                EndPosition1: e1,
                EndPosition2: e2);
        }

        /// <summary>Consumed intervals of a row pair: a leading run of gap columns of a free side is the free start flank, a
        /// trailing run the free end flank (the first / last column decides which side; FreeStart1 / FreeEnd1 are checked first).</summary>
        private (int Start1, int Start2, int End1, int End2) ConsumedInterval(string row1, string row2)
        {
            int len = row1.Length, lead = 0, trail = 0;
            int start1 = 0, start2 = 0, end1 = _m - 1, end2 = _n - 1;
            if (len > 0 && row2[0] == '-' && _freeStart1)
            {
                while (lead < len && row2[lead] == '-') lead++;
                start1 = lead;
            }
            else if (len > 0 && row1[0] == '-' && _freeStart2)
            {
                while (lead < len && row1[lead] == '-') lead++;
                start2 = lead;
            }

            if (len > 0 && row2[len - 1] == '-' && _freeEnd1)
            {
                while (trail < len - lead && row2[len - 1 - trail] == '-') trail++;
                end1 = _m - 1 - trail;
            }
            else if (len > 0 && row1[len - 1] == '-' && _freeEnd2)
            {
                while (trail < len - lead && row1[len - 1 - trail] == '-') trail++;
                end2 = _n - 1 - trail;
            }

            return (start1, start2, end1, end2);
        }
    }

    #endregion
}
