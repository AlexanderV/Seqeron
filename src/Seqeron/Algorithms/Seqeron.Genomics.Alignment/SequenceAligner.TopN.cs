namespace Seqeron.Genomics.Alignment;

// =====================================================================================================================
// Suboptimal local alignments (audit A1-7).
//
//   * Waterman-Eggert declumping (Waterman MS, Eggert M (1987) "A new algorithm for best subsequence alignments with
//     application to tRNA-rRNA comparisons", J Mol Biol 197:723-728; the algorithm of FASTA lalign, EMBOSS matcher
//     -alternatives (Huang-Miller SIM implements the same definition in linear space), Bioconductor/Biopython-free
//     "top-N Smith-Waterman"): report the best local alignment, forbid the aligned residue pairs it used, report the best
//     local alignment that uses none of the forbidden pairs, and so on. Exact definition used here (identical to the one
//     implemented by EMBOSS matcher, checked on 2 902 reported alignments, see docs/Evidence/ALIGN-GLOBAL-001-Evidence.md (A1-7)):
//
//       A pair is an ALIGNED PAIR (i, j) = a_i opposite b_j in a diagonal column (match or mismatch); a column that
//       aligns a residue to a gap is not a pair. Alignment k is a local alignment of maximum score among those that
//       contain none of the pairs of alignments 1..k-1 (so any two reported alignments share no aligned pair, although
//       their paths may cross through gap columns).
//
//     Implementation (the original recomputation scheme): the full DP matrices are kept; after an alignment is reported its
//     pairs are forbidden (the substitution/diagonal move into that cell is removed, gap moves are unchanged) and ONLY the
//     cells whose value can change are recomputed (a cell is recomputed when it is a new forbidden pair or one of its
//     three predecessors changed), row by row, until no cell of a row changes. The recurrences are the very same
//     LocalLinearCell / AffineCell that LocalAlignCore / AffineAlignCore fill the matrix with, and the traceback is the
//     one of LocalAlign / LocalAlignAffine (diagonal > up > left, M > X > Y; the alignment ends at the first maximal cell in
//     row-major order and stops at the first zero predecessor), so alignment 1 is exactly LocalAlign(...).
//
//   * Co-optimal enumeration (Biopython PairwiseAligner.align(...) iteration for mode='local'): all distinct local alignments
//     of the maximum score, i.e. every traceback path of the matrix from every maximal cell, in a deterministic order
//     (end cell in row-major order, then diagonal > up > left / M > X > Y at every branch), capped.
// =====================================================================================================================
public static partial class SequenceAligner
{
    #region Waterman-Eggert top-N local alignments and co-optimal enumeration

    /// <summary>
    /// Largest number of dynamic-programming cells <c>(m + 1)·(n + 1)</c> for which the Waterman-Eggert declumper and the co-optimal
    /// enumerator keep the matrices, as the original algorithm does (4 bytes per cell for the linear model, 12 for the affine
    /// one — 25 M cells = 100 MB / 300 MB). Above it the top-N alignments (<c>LocalAlignTopN</c> / <c>LocalAlignAffineTopN</c>, all
    /// overloads) are computed by the exact LINEAR-SPACE program (<see cref="LocalAlignTopNLinearSpace(string, string, ScoringMatrix?, int, int, long)"/>:
    /// the same alignments in O(m + n) memory per round); the co-optimal enumerator (<c>LocalAlignCoOptimal</c>) still keeps
    /// the matrices and throws <see cref="ArgumentException"/> above it.
    /// </summary>
    public const long TopNMaxCells = LinearSpaceCellThreshold;

    /// <summary>
    /// Waterman-Eggert (1987) non-intersecting local alignments, <b>linear</b> gap model: the best local alignment
    /// (identical to <see cref="LocalAlign(DnaSequence, DnaSequence, ScoringMatrix?)"/>), then the best one that shares no
    /// aligned residue pair with it, then the best one sharing no pair with either, … up to
    /// <paramref name="maxAlignments"/> alignments of score ≥ <paramref name="minScore"/> (and always &gt; 0).
    /// </summary>
    /// <param name="sequence1">First sequence.</param>
    /// <param name="sequence2">Second sequence.</param>
    /// <param name="scoring">Scoring (default <see cref="SimpleDna"/>); only <see cref="ScoringMatrix.GapExtend"/> is used.</param>
    /// <param name="maxAlignments">Maximum number of alignments to report (≥ 1).</param>
    /// <param name="minScore">Smallest accepted score (inclusive). Scores ≤ 0 are never reported.</param>
    /// <returns>The alignments in the order found = non-increasing score. Fewer than <paramref name="maxAlignments"/> are
    /// returned when no further positive alignment ≥ <paramref name="minScore"/> exists (empty list for an empty side or no
    /// positive-scoring pair).</returns>
    /// <remarks>
    /// <para>
    /// Definition (Waterman &amp; Eggert 1987; the same one EMBOSS <c>matcher -alternatives</c> implements): an
    /// <i>aligned pair</i> is a diagonal column (a_i opposite b_j); alignment k is a maximum-score local alignment among
    /// those that contain none of the aligned pairs of alignments 1..k-1. Columns that align a residue with a gap are not pairs,
    /// so two reported alignments may pass through the same matrix cell by a gap move but never share a matched/mismatched pair.
    /// </para>
    /// <para>
    /// Ties: the choice among co-optimal alignments is the fixed one of <see cref="LocalAlign(DnaSequence, DnaSequence, ScoringMatrix?)"/>
    /// (first maximal cell in row-major order, traceback diagonal &gt; up &gt; left, stop at the first zero cell), applied to the
    /// matrix with the pairs of the earlier alignments removed. Because a later alignment depends on which co-optimal alignment was chosen
    /// earlier, other tools (matcher, lalign) may report a different — equally valid — sequence of alignments on tie-rich inputs
    /// (the score of every alignment is always the true maximum of the remaining problem).
    /// </para>
    /// <para>Memory is O(m·n) up to <see cref="TopNMaxCells"/> cells, as in the original algorithm (matrices kept, only the changed cells are recomputed
    /// after each report). Above it the very same alignments are computed in O(m + n) memory by
    /// <see cref="LocalAlignTopNLinearSpace(string, string, ScoringMatrix?, int, int, long)"/> (one run of the exact linear-space local program per
    /// reported alignment, Huang-Miller SIM / EMBOSS matcher scheme; needs non-positive gap scores).</para>
    /// </remarks>
    /// <exception cref="ArgumentNullException">A sequence is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="maxAlignments"/> &lt; 1.</exception>
    /// <exception cref="ArgumentException">A gap score is positive or int.MinValue, or the integer score range is exceeded
    /// (see <see cref="ScoringMatrix"/>).</exception>
    [System.Runtime.CompilerServices.OverloadResolutionPriority(1)]
    public static IReadOnlyList<AlignmentResult> LocalAlignTopN(
        DnaSequence sequence1, DnaSequence sequence2, ScoringMatrix? scoring = null, int maxAlignments = 10, int minScore = 0)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        return LocalTopNCore(sequence1.Sequence, sequence2.Sequence, scoring ?? SimpleDna, affine: false, maxAlignments, minScore);
    }

    /// <summary><see cref="LocalAlignTopN(DnaSequence, DnaSequence, ScoringMatrix?, int, int)"/> on raw strings (upper-cased internally).</summary>
    [System.Runtime.CompilerServices.OverloadResolutionPriority(1)]
    public static IReadOnlyList<AlignmentResult> LocalAlignTopN(
        string sequence1, string sequence2, ScoringMatrix? scoring = null, int maxAlignments = 10, int minScore = 0)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        return LocalTopNCore(sequence1.ToUpperInvariant(), sequence2.ToUpperInvariant(), scoring ?? SimpleDna, affine: false, maxAlignments, minScore);
    }

    /// <summary>
    /// Waterman-Eggert non-intersecting local alignments with <b>affine</b> gaps (a gap of length k scores
    /// <see cref="ScoringMatrix.GapOpen"/> + k·<see cref="ScoringMatrix.GapExtend"/>): alignment 1 is
    /// <see cref="LocalAlignAffine(DnaSequence, DnaSequence, ScoringMatrix?)"/>; definition, ties and limits as
    /// <see cref="LocalAlignTopN(DnaSequence, DnaSequence, ScoringMatrix?, int, int)"/> (EMBOSS <c>matcher</c> with
    /// <c>-gapopen G -gapextend E</c> = GapOpen −(G − E), GapExtend −E).
    /// </summary>
    [System.Runtime.CompilerServices.OverloadResolutionPriority(1)]
    public static IReadOnlyList<AlignmentResult> LocalAlignAffineTopN(
        DnaSequence sequence1, DnaSequence sequence2, ScoringMatrix? scoring = null, int maxAlignments = 10, int minScore = 0)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        return LocalTopNCore(sequence1.Sequence, sequence2.Sequence, scoring ?? SimpleDna, affine: true, maxAlignments, minScore);
    }

    /// <summary><see cref="LocalAlignAffineTopN(DnaSequence, DnaSequence, ScoringMatrix?, int, int)"/> on raw strings.</summary>
    [System.Runtime.CompilerServices.OverloadResolutionPriority(1)]
    public static IReadOnlyList<AlignmentResult> LocalAlignAffineTopN(
        string sequence1, string sequence2, ScoringMatrix? scoring = null, int maxAlignments = 10, int minScore = 0)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        return LocalTopNCore(sequence1.ToUpperInvariant(), sequence2.ToUpperInvariant(), scoring ?? SimpleDna, affine: true, maxAlignments, minScore);
    }

    /// <summary>Waterman-Eggert top-N with a substitution matrix and linear gap scores
    /// (<see cref="LocalAlign(string, string, MatrixScoring)"/>); <paramref name="minScore"/> is a real score, each result carries
    /// <see cref="AlignmentResult.ExactScore"/>.</summary>
    public static IReadOnlyList<AlignmentResult> LocalAlignTopN(
        string sequence1, string sequence2, MatrixScoring scoring, int maxAlignments = 10, double minScore = 0) =>
        TopNWithMatrix(sequence1, sequence2, scoring, affine: false, maxAlignments, minScore, coOptimal: false);

    /// <summary><see cref="LocalAlignTopN(string, string, MatrixScoring, int, double)"/> for protein sequences.</summary>
    public static IReadOnlyList<AlignmentResult> LocalAlignTopN(
        ProteinSequence sequence1, ProteinSequence sequence2, MatrixScoring scoring, int maxAlignments = 10, double minScore = 0)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        return TopNWithMatrix(sequence1.Sequence, sequence2.Sequence, scoring, affine: false, maxAlignments, minScore, coOptimal: false);
    }

    /// <summary><see cref="LocalAlignTopN(string, string, MatrixScoring, int, double)"/> for DNA sequences.</summary>
    public static IReadOnlyList<AlignmentResult> LocalAlignTopN(
        DnaSequence sequence1, DnaSequence sequence2, MatrixScoring scoring, int maxAlignments = 10, double minScore = 0)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        return TopNWithMatrix(sequence1.Sequence, sequence2.Sequence, scoring, affine: false, maxAlignments, minScore, coOptimal: false);
    }

    /// <summary>Waterman-Eggert top-N with a substitution matrix and affine gap scores
    /// (<see cref="LocalAlignAffine(string, string, MatrixScoring)"/>; EMBOSS matcher with <see cref="MatrixScoring.Emboss"/> numbers).</summary>
    public static IReadOnlyList<AlignmentResult> LocalAlignAffineTopN(
        string sequence1, string sequence2, MatrixScoring scoring, int maxAlignments = 10, double minScore = 0) =>
        TopNWithMatrix(sequence1, sequence2, scoring, affine: true, maxAlignments, minScore, coOptimal: false);

    /// <summary><see cref="LocalAlignAffineTopN(string, string, MatrixScoring, int, double)"/> for protein sequences.</summary>
    public static IReadOnlyList<AlignmentResult> LocalAlignAffineTopN(
        ProteinSequence sequence1, ProteinSequence sequence2, MatrixScoring scoring, int maxAlignments = 10, double minScore = 0)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        return TopNWithMatrix(sequence1.Sequence, sequence2.Sequence, scoring, affine: true, maxAlignments, minScore, coOptimal: false);
    }

    /// <summary><see cref="LocalAlignAffineTopN(string, string, MatrixScoring, int, double)"/> for DNA sequences.</summary>
    public static IReadOnlyList<AlignmentResult> LocalAlignAffineTopN(
        DnaSequence sequence1, DnaSequence sequence2, MatrixScoring scoring, int maxAlignments = 10, double minScore = 0)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        return TopNWithMatrix(sequence1.Sequence, sequence2.Sequence, scoring, affine: true, maxAlignments, minScore, coOptimal: false);
    }

    /// <summary>
    /// All co-optimal local alignments (linear gap model): every distinct local alignment whose score equals the Smith-Waterman
    /// optimum — the set that Biopython's <c>PairwiseAligner(mode='local').align(a, b)</c> iterates — capped at
    /// <paramref name="maxAlignments"/>.
    /// </summary>
    /// <param name="sequence1">First sequence.</param>
    /// <param name="sequence2">Second sequence.</param>
    /// <param name="scoring">Scoring (default <see cref="SimpleDna"/>); only <see cref="ScoringMatrix.GapExtend"/> is used.</param>
    /// <param name="maxAlignments">Cap on the number of alignments returned (≥ 1; the number of co-optimal alignments can grow
    /// exponentially with tandem repeats).</param>
    /// <remarks>
    /// An alignment is the path of a traceback that starts (reading backwards) at a maximal cell of the matrix reached by a
    /// diagonal move (an alignment ends with an aligned pair) and follows every predecessor that reproduces the cell value, until
    /// the predecessor value is 0 (the local alignment then begins with the aligned pair just consumed). As in Biopython ("zero-score
    /// extensions to alignments are not allowed") a path never passes through another end cell of an optimal alignment, so an
    /// optimal alignment is not also reported prolonged by a tail of score 0 (e.g. a gap run followed by an equal number of
    /// compensating matches, or a 0-scoring mismatch). Two different paths are two different alignments (different gapped rows or coordinates).
    /// Order: end cell in row-major order, then diagonal &gt; up &gt; left at each branch (depth first), so the first element is
    /// <see cref="LocalAlign(DnaSequence, DnaSequence, ScoringMatrix?)"/>'s result. Empty list when no pair scores &gt; 0.
    /// Same memory bound as <see cref="LocalAlignTopN(DnaSequence, DnaSequence, ScoringMatrix?, int, int)"/>.
    /// </remarks>
    [System.Runtime.CompilerServices.OverloadResolutionPriority(1)]
    public static IReadOnlyList<AlignmentResult> LocalAlignCoOptimal(
        DnaSequence sequence1, DnaSequence sequence2, ScoringMatrix? scoring = null, int maxAlignments = 1000)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        return LocalCoOptimalCore(sequence1.Sequence, sequence2.Sequence, scoring ?? SimpleDna, affine: false, maxAlignments);
    }

    /// <summary><see cref="LocalAlignCoOptimal(DnaSequence, DnaSequence, ScoringMatrix?, int)"/> on raw strings.</summary>
    [System.Runtime.CompilerServices.OverloadResolutionPriority(1)]
    public static IReadOnlyList<AlignmentResult> LocalAlignCoOptimal(
        string sequence1, string sequence2, ScoringMatrix? scoring = null, int maxAlignments = 1000)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        return LocalCoOptimalCore(sequence1.ToUpperInvariant(), sequence2.ToUpperInvariant(), scoring ?? SimpleDna, affine: false, maxAlignments);
    }

    /// <summary>All co-optimal local alignments under the <b>affine</b> gap model (GapOpen + k·GapExtend; Biopython
    /// <c>mode='local'</c> with open = GapOpen + GapExtend, extend = GapExtend), capped; see
    /// <see cref="LocalAlignCoOptimal(DnaSequence, DnaSequence, ScoringMatrix?, int)"/>.</summary>
    [System.Runtime.CompilerServices.OverloadResolutionPriority(1)]
    public static IReadOnlyList<AlignmentResult> LocalAlignAffineCoOptimal(
        DnaSequence sequence1, DnaSequence sequence2, ScoringMatrix? scoring = null, int maxAlignments = 1000)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        return LocalCoOptimalCore(sequence1.Sequence, sequence2.Sequence, scoring ?? SimpleDna, affine: true, maxAlignments);
    }

    /// <summary><see cref="LocalAlignAffineCoOptimal(DnaSequence, DnaSequence, ScoringMatrix?, int)"/> on raw strings.</summary>
    [System.Runtime.CompilerServices.OverloadResolutionPriority(1)]
    public static IReadOnlyList<AlignmentResult> LocalAlignAffineCoOptimal(
        string sequence1, string sequence2, ScoringMatrix? scoring = null, int maxAlignments = 1000)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        return LocalCoOptimalCore(sequence1.ToUpperInvariant(), sequence2.ToUpperInvariant(), scoring ?? SimpleDna, affine: true, maxAlignments);
    }

    /// <summary>All co-optimal local alignments with a substitution matrix and linear gap scores
    /// (<see cref="AlignmentResult.ExactScore"/> = real score); see <see cref="LocalAlignCoOptimal(DnaSequence, DnaSequence, ScoringMatrix?, int)"/>.</summary>
    public static IReadOnlyList<AlignmentResult> LocalAlignCoOptimal(
        string sequence1, string sequence2, MatrixScoring scoring, int maxAlignments = 1000) =>
        TopNWithMatrix(sequence1, sequence2, scoring, affine: false, maxAlignments, 0, coOptimal: true);

    /// <summary><see cref="LocalAlignCoOptimal(string, string, MatrixScoring, int)"/> for protein sequences.</summary>
    public static IReadOnlyList<AlignmentResult> LocalAlignCoOptimal(
        ProteinSequence sequence1, ProteinSequence sequence2, MatrixScoring scoring, int maxAlignments = 1000)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        return TopNWithMatrix(sequence1.Sequence, sequence2.Sequence, scoring, affine: false, maxAlignments, 0, coOptimal: true);
    }

    /// <summary><see cref="LocalAlignCoOptimal(string, string, MatrixScoring, int)"/> for DNA sequences.</summary>
    public static IReadOnlyList<AlignmentResult> LocalAlignCoOptimal(
        DnaSequence sequence1, DnaSequence sequence2, MatrixScoring scoring, int maxAlignments = 1000)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        return TopNWithMatrix(sequence1.Sequence, sequence2.Sequence, scoring, affine: false, maxAlignments, 0, coOptimal: true);
    }

    /// <summary>All co-optimal affine-gap local alignments with a substitution matrix.</summary>
    public static IReadOnlyList<AlignmentResult> LocalAlignAffineCoOptimal(
        string sequence1, string sequence2, MatrixScoring scoring, int maxAlignments = 1000) =>
        TopNWithMatrix(sequence1, sequence2, scoring, affine: true, maxAlignments, 0, coOptimal: true);

    /// <summary><see cref="LocalAlignAffineCoOptimal(string, string, MatrixScoring, int)"/> for protein sequences.</summary>
    public static IReadOnlyList<AlignmentResult> LocalAlignAffineCoOptimal(
        ProteinSequence sequence1, ProteinSequence sequence2, MatrixScoring scoring, int maxAlignments = 1000)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        return TopNWithMatrix(sequence1.Sequence, sequence2.Sequence, scoring, affine: true, maxAlignments, 0, coOptimal: true);
    }

    /// <summary><see cref="LocalAlignAffineCoOptimal(string, string, MatrixScoring, int)"/> for DNA sequences.</summary>
    public static IReadOnlyList<AlignmentResult> LocalAlignAffineCoOptimal(
        DnaSequence sequence1, DnaSequence sequence2, MatrixScoring scoring, int maxAlignments = 1000)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        return TopNWithMatrix(sequence1.Sequence, sequence2.Sequence, scoring, affine: true, maxAlignments, 0, coOptimal: true);
    }

    /// <summary>
    /// The aligned residue pairs (1-based (i, j): a_i opposite b_j in a diagonal column) of a local alignment — the cells
    /// that Waterman-Eggert declumping forbids. Public so that callers and tests can state non-intersection exactly.
    /// </summary>
    public static IReadOnlyList<(int Row, int Column)> AlignedPairs(AlignmentResult alignment)
    {
        ArgumentNullException.ThrowIfNull(alignment);
        var pairs = new List<(int, int)>();
        if (alignment.AlignedSequence1.Length == 0 || alignment.StartPosition1 < 0)
            return pairs;
        int i = alignment.StartPosition1;
        int j = alignment.StartPosition2;
        for (int c = 0; c < alignment.AlignedSequence1.Length; c++)
        {
            bool gap1 = alignment.AlignedSequence1[c] == '-';
            bool gap2 = alignment.AlignedSequence2[c] == '-';
            if (!gap1) i++;
            if (!gap2) j++;
            if (!gap1 && !gap2) pairs.Add((i, j));
        }
        return pairs;
    }

    private static IReadOnlyList<AlignmentResult> TopNWithMatrix(
        string sequence1, string sequence2, MatrixScoring scoring, bool affine, int maxAlignments, double minScore, bool coOptimal)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        ArgumentNullException.ThrowIfNull(scoring);
        if (double.IsNaN(minScore))
            throw new ArgumentOutOfRangeException(nameof(minScore), "minScore must not be NaN.");

        string a = sequence1.ToUpperInvariant();
        string b = sequence2.ToUpperInvariant();
        var ps = PairScoring.Create(scoring, a, b);

        long minScaled = (long)Math.Ceiling((decimal)Math.Clamp(minScore, -1e12, 1e12) * ps.Scale);
        var raw = coOptimal
            ? LocalCoOptimalCore(a, b, ps, affine, maxAlignments)
            : LocalTopNCore(a, b, ps, affine, maxAlignments, minScaled);
        var list = new List<AlignmentResult>(raw.Count);
        foreach (var r in raw)
        {
            double exact = (double)r.Score / ps.Scale;
            list.Add(r with { Score = (int)Math.Round(exact, MidpointRounding.AwayFromZero), ExactScore = exact });
        }
        return list;
    }

    private static void RequireTopNFits(int m, int n, string paramName)
    {
        if ((long)(m + 1) * (n + 1) > TopNMaxCells)
            throw new ArgumentException(
                $"Waterman-Eggert declumping and co-optimal enumeration keep the O(m*n) matrices: ({m} + 1) x ({n} + 1) cells exceed the limit of {TopNMaxCells.ToString("N0", System.Globalization.CultureInfo.InvariantCulture)} cells.",
                paramName);
    }

    /// <param name="a">First sequence (upper-case).</param>
    /// <param name="b">Second sequence (upper-case).</param>
    /// <param name="scoring">Exact integer scoring model.</param>
    /// <param name="affine">Affine (Gotoh) instead of linear gaps.</param>
    /// <param name="maxAlignments">At most this many alignments (&gt;= 1).</param>
    /// <param name="minScaled">Smallest accepted score in the integer units of <paramref name="scoring"/>.</param>
    /// <param name="fullMatrixCells">Inputs of at most this many cells use the in-memory declumper, larger ones the linear-space
    /// program (tests lower it to force the linear path on small inputs).</param>
    /// <param name="baseCells">Base-case size of the linear-space program's divide and conquer.</param>
    internal static IReadOnlyList<AlignmentResult> LocalTopNCore(
        string a, string b, PairScoring scoring, bool affine, int maxAlignments, long minScaled,
        long fullMatrixCells = TopNMaxCells, long baseCells = LinearSpaceBaseCells)
    {
        if (maxAlignments < 1)
            throw new ArgumentOutOfRangeException(nameof(maxAlignments), maxAlignments, "maxAlignments must be at least 1.");
        var results = new List<AlignmentResult>();
        if (a.Length == 0 || b.Length == 0)
            return results;
        scoring.EnsureRange(a.Length, b.Length);
        if ((long)(a.Length + 1) * (b.Length + 1) > fullMatrixCells)
            return LocalTopNLinearSpaceCore(a, b, scoring, affine, maxAlignments, minScaled, baseCells);

        var declumper = new LocalDeclumper(a, b, scoring, affine);
        while (results.Count < maxAlignments)
        {
            var next = declumper.Next(minScaled);
            if (next is null)
                break;
            results.Add(next);
        }
        return results;
    }

    internal static IReadOnlyList<AlignmentResult> LocalCoOptimalCore(
        string a, string b, PairScoring scoring, bool affine, int maxAlignments)
    {
        if (maxAlignments < 1)
            throw new ArgumentOutOfRangeException(nameof(maxAlignments), maxAlignments, "maxAlignments must be at least 1.");
        var results = new List<AlignmentResult>();
        if (a.Length == 0 || b.Length == 0)
            return results;
        RequireTopNFits(a.Length, b.Length, nameof(a));
        scoring.EnsureRange(a.Length, b.Length);

        var matrices = new LocalDeclumper(a, b, scoring, affine);
        matrices.EnumerateCoOptimal(maxAlignments, results);
        return results;
    }

    /// <summary>
    /// The Waterman-Eggert state: the filled local-alignment matrices (linear: H; affine: M, X, Y), the set of forbidden
    /// aligned pairs and the best score of every row. Cells are filled and recomputed with the very same recurrences as
    /// <see cref="LocalAlignCore"/> / <see cref="AffineAlignCore"/>.
    /// </summary>
    private sealed class LocalDeclumper
    {
        private readonly string _a;
        private readonly string _b;
        private readonly PairScoring _sc;
        private readonly bool _affine;
        private readonly int _m;
        private readonly int _n;
        private readonly int[,]? _h;
        private readonly int[,]? _mm;
        private readonly int[,]? _xx;
        private readonly int[,]? _yy;
        private readonly HashSet<long> _forbidden = new();
        private readonly int[] _rowBest;
        private readonly int[] _rowBestCol;

        public LocalDeclumper(string a, string b, PairScoring sc, bool affine)
        {
            _a = a;
            _b = b;
            _sc = sc;
            _affine = affine;
            _m = a.Length;
            _n = b.Length;
            _rowBest = new int[_m + 1];
            _rowBestCol = new int[_m + 1];

            if (affine)
            {
                (_mm, _xx, _yy) = AffineNewMatrices(_m, _n);
                // The very borders of AffineAlignCore (local): every border state is unreachable.
                AffineInitBorders(_mm, _xx, _yy, sc, local: true, EndGaps.None, StateMatch);
            }
            else
            {
                _h = new int[_m + 1, _n + 1];
            }

            for (int i = 1; i <= _m; i++)
            {
                AlignmentCancellation.ThrowIfRequested();
                for (int j = 1; j <= _n; j++)
                    ComputeCell(i, j);
                ScanRow(i);
            }
        }

        private int[,] Score => _affine ? _mm! : _h!;

        /// <summary>Recomputes cell (i, j); returns true when any of its values changed.</summary>
        private bool ComputeCell(int i, int j)
        {
            bool forbidden = _forbidden.Count > 0 && _forbidden.Contains(PackCell(i, j));
            int sub = _sc.Sub(_a[i - 1], _b[j - 1]);
            if (!_affine)
            {
                int old = _h![i, j];
                int value = LocalLinearCell(_h, i, j, sub, _sc.GapExtend, forbidden);
                _h[i, j] = value;
                return value != old;
            }

            int om = _mm![i, j], ox = _xx![i, j], oy = _yy![i, j];
            AffineCell(_mm, _xx, _yy, i, j, sub, startFloor: 0, _sc.GapExtend, _sc.GapOpen + _sc.GapExtend, forbidden);
            return om != _mm[i, j] || ox != _xx[i, j] || oy != _yy[i, j];
        }

        /// <summary>First maximal cell of row i (smallest column on ties; 0 when the row has no positive cell).</summary>
        private void ScanRow(int i)
        {
            var score = Score;
            int best = 0, bestJ = 0;
            for (int j = 1; j <= _n; j++)
                if (score[i, j] > best) { best = score[i, j]; bestJ = j; }
            _rowBest[i] = best;
            _rowBestCol[i] = bestJ;
        }

        /// <summary>The next alignment (best remaining), or null when none scores &gt; 0 and ≥ minScaled.</summary>
        public AlignmentResult? Next(long minScaled)
        {
            int best = 0, bi = 0;
            for (int i = 1; i <= _m; i++)
                if (_rowBest[i] > best) { best = _rowBest[i]; bi = i; }
            if (best <= 0 || best < minScaled)
                return null;
            int bj = _rowBestCol[bi];

            var result = _affine
                ? AffineTraceback(_a, _b, _mm!, _xx!, _yy!, _sc, local: true, EndGaps.None, bi, bj, StateMatch, best)
                : TracebackLocal(_a, _b, _h!, bi, bj, _sc, _forbidden);

            var fresh = new List<(int Row, int Column)>(AlignedPairs(result));
            foreach (var (row, col) in fresh)
                _forbidden.Add(PackCell(row, col));
            Recompute(fresh);
            return result;
        }

        /// <summary>
        /// Waterman-Eggert recomputation: only cells that are new forbidden pairs, or whose diagonal / upper / left
        /// predecessor changed, are recomputed; the sweep stops at the first row with no change.
        /// </summary>
        private void Recompute(List<(int Row, int Column)> fresh)
        {
            if (fresh.Count == 0)
                return;
            var byRow = new Dictionary<int, List<int>>();
            int firstRow = int.MaxValue, lastRow = 0;
            foreach (var (row, col) in fresh)
            {
                if (!byRow.TryGetValue(row, out var list)) byRow[row] = list = new List<int>();
                list.Add(col);
                firstRow = Math.Min(firstRow, row);
                lastRow = Math.Max(lastRow, row);
            }

            var prev = new bool[_n + 2];   // changed flags of row i-1, indexed by column
            var cur = new bool[_n + 2];
            int prevFirst = _n + 1, prevLast = 0;

            for (int i = firstRow; i <= _m; i++)
            {
                AlignmentCancellation.ThrowIfRequested();
                byRow.TryGetValue(i, out var freshCols);
                if (freshCols is null && prevLast == 0)
                {
                    if (i > lastRow) break;
                    continue;
                }

                int jStart = prevLast == 0 ? _n + 1 : prevFirst;
                int limit = prevLast == 0 ? 0 : Math.Min(_n, prevLast + 1);
                HashSet<int>? freshSet = null;
                if (freshCols is not null)
                {
                    freshSet = new HashSet<int>(freshCols);
                    foreach (int c in freshCols) { jStart = Math.Min(jStart, c); limit = Math.Max(limit, c); }
                }

                Array.Clear(cur);
                int curFirst = _n + 1, curLast = 0;
                bool left = false;
                for (int j = jStart; j <= _n; j++)
                {
                    if (!left && j > limit)
                        break;
                    bool dirty = left || prev[j] || prev[j - 1] || (freshSet is not null && freshSet.Contains(j));
                    if (!dirty)
                    {
                        left = false;
                        continue;
                    }
                    bool changed = ComputeCell(i, j);
                    cur[j] = changed;
                    left = changed;
                    if (changed)
                    {
                        if (j < curFirst) curFirst = j;
                        curLast = j;
                    }
                }

                if (curLast != 0 || freshCols is not null)
                    ScanRow(i);

                (prev, cur) = (cur, prev);
                prevFirst = curFirst;
                prevLast = curLast;
            }
        }

        /// <summary>
        /// Depth-first enumeration of all traceback paths from every maximal cell (row-major), diagonal &gt; up &gt; left
        /// (M &gt; X &gt; Y), until the predecessor value is 0, never passing through an end cell of an optimal alignment; at most <paramref name="cap"/> alignments are added to <paramref name="output"/>.
        /// </summary>
        public void EnumerateCoOptimal(int cap, List<AlignmentResult> output)
        {
            int best = 0;
            for (int i = 1; i <= _m; i++) best = Math.Max(best, _rowBest[i]);
            if (best <= 0)
                return;
            _best = best;

            var steps = new List<byte>(); // columns from the end backwards: 0 = pair, 1 = a_i vs gap, 2 = gap vs b_j
            for (int ei = 1; ei <= _m && output.Count < cap; ei++)
            {
                AlignmentCancellation.ThrowIfRequested();
                if (_rowBest[ei] != best)
                    continue;
                for (int ej = 1; ej <= _n && output.Count < cap; ej++)
                {
                    if (!IsEndPoint(ei, ej))
                        continue;
                    steps.Clear();
                    Walk(ei, ej, StateMatch, ei, ej, steps, cap, output);
                }
            }
        }

        private int _best;

        /// <summary>
        /// A cell where an optimal alignment may end: it holds the optimum and is reached by a diagonal move (Biopython's
        /// ENDPOINT; an alignment cannot end with a gap column). A path may not pass through such a cell and continue:
        /// that would be a zero-score extension of a shorter optimal alignment (Biopython: "we don't allow zero-score
        /// extensions to alignments").
        /// </summary>
        private bool IsEndPoint(int i, int j)
        {
            if (i == 0 || j == 0 || Score[i, j] != _best)
                return false;
            return _affine || _h![i, j] == _h[i - 1, j - 1] + _sc.Sub(_a[i - 1], _b[j - 1]);
        }

        // Iterative DFS over the (cell, state) traceback graph; frames hold the next branch to try.
        private void Walk(int endI, int endJ, byte startState, int i0, int j0, List<byte> steps, int cap, List<AlignmentResult> output)
        {
            var stack = new Stack<(int I, int J, byte State, int Branch, int Depth)>();
            stack.Push((i0, j0, startState, 0, 0));
            int e = _sc.GapExtend;
            int oe = _sc.GapOpen + e;
            while (stack.Count > 0 && output.Count < cap)
            {
                var (i, j, state, branch, depth) = stack.Pop();
                steps.RemoveRange(depth, steps.Count - depth);

                // The successor candidates of this node in the fixed tie order; branch = number already explored.
                int tried = 0;
                bool pushed = false;
                for (int cand = 0; cand < 3 && !pushed; cand++)
                {
                    if (!TryPredecessor(i, j, state, cand, e, oe, out int pi, out int pj, out byte pstate, out byte kind, out bool starts))
                        continue;
                    if (!starts && pstate == StateMatch && IsEndPoint(pi, pj))
                        continue;   // would extend a shorter optimal alignment by a zero-score tail
                    if (tried++ < branch)
                        continue;
                    // Re-push this node with the next branch, then descend.
                    stack.Push((i, j, state, branch + 1, depth));
                    steps.Add(kind);
                    if (starts)
                    {
                        output.Add(BuildAlignment(endI, endJ, pi, pj, steps));
                        steps.RemoveAt(steps.Count - 1);
                    }
                    else
                    {
                        stack.Push((pi, pj, pstate, 0, depth + 1));
                    }
                    pushed = true;
                }
            }
        }

        /// <summary>
        /// Candidate predecessor number <paramref name="cand"/> (0 = diagonal / M, 1 = up / X, 2 = left / Y) of the node
        /// (i, j, state), if it reproduces the node's value. <paramref name="starts"/>: the predecessor has value 0, i.e. the
        /// alignment begins with the consumed column.
        /// </summary>
        private bool TryPredecessor(
            int i, int j, byte state, int cand, int e, int oe,
            out int pi, out int pj, out byte pstate, out byte kind, out bool starts)
        {
            pi = i; pj = j; pstate = StateMatch; kind = 0; starts = false;
            if (!_affine)
            {
                var h = _h!;
                int v = h[i, j];
                switch (cand)
                {
                    case 0:
                        if (_forbidden.Contains(PackCell(i, j)) || v != h[i - 1, j - 1] + _sc.Sub(_a[i - 1], _b[j - 1])) return false;
                        pi = i - 1; pj = j - 1; kind = 0; break;
                    case 1:
                        if (v != h[i - 1, j] + e) return false;
                        pi = i - 1; kind = 1; break;
                    default:
                        if (v != h[i, j - 1] + e) return false;
                        pj = j - 1; kind = 2; break;
                }
                starts = pi == 0 || pj == 0 || h[pi, pj] <= 0;
                return true;
            }

            var mm = _mm!; var xx = _xx!; var yy = _yy!;
            if (state == StateMatch)
            {
                int prev = mm[i, j] - _sc.Sub(_a[i - 1], _b[j - 1]);
                kind = 0; pi = i - 1; pj = j - 1;
                if (prev == 0)
                {
                    // The alignment begins with this pair: a single "start" candidate.
                    if (cand != 0) return false;
                    starts = true;
                    return true;
                }
                int pv = cand switch { 0 => mm[pi, pj], 1 => xx[pi, pj], _ => yy[pi, pj] };
                if (pv != prev || pi == 0 || pj == 0) return false;
                pstate = cand switch { 0 => StateMatch, 1 => StateGapInSeq2, _ => StateGapInSeq1 };
                return true;
            }
            if (state == StateGapInSeq2)
            {
                int cur = xx[i, j];
                pi = i - 1; kind = 1;
                if (pi == 0) return false;
                int pv = cand switch { 0 => mm[pi, j] + oe, 1 => xx[pi, j] + e, _ => yy[pi, j] + oe };
                if (pv != cur) return false;
                pstate = cand switch { 0 => StateMatch, 1 => StateGapInSeq2, _ => StateGapInSeq1 };
                return true;
            }
            else
            {
                int cur = yy[i, j];
                pj = j - 1; kind = 2;
                if (pj == 0) return false;
                int pv = cand switch { 0 => mm[i, pj] + oe, 1 => yy[i, pj] + e, _ => xx[i, pj] + oe };
                if (pv != cur) return false;
                pstate = cand switch { 0 => StateMatch, 1 => StateGapInSeq1, _ => StateGapInSeq2 };
                return true;
            }
        }

        private AlignmentResult BuildAlignment(int endI, int endJ, int startI, int startJ, List<byte> steps)
        {
            var r1 = new char[steps.Count];
            var r2 = new char[steps.Count];
            int i = endI, j = endJ;
            for (int k = 0; k < steps.Count; k++)
            {
                int col = steps.Count - 1 - k;   // steps are stored from the end backwards
                switch (steps[k])
                {
                    case 0: r1[col] = _a[i - 1]; r2[col] = _b[j - 1]; i--; j--; break;
                    case 1: r1[col] = _a[i - 1]; r2[col] = '-'; i--; break;
                    default: r1[col] = '-'; r2[col] = _b[j - 1]; j--; break;
                }
            }
            return new AlignmentResult(
                new string(r1), new string(r2), Score[endI, endJ], AlignmentType.Local,
                StartPosition1: startI, StartPosition2: startJ, EndPosition1: endI - 1, EndPosition2: endJ - 1);
        }
    }

    #endregion
}
