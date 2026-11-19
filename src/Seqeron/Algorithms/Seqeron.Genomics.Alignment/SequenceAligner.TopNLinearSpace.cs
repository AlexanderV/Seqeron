using System.Runtime.CompilerServices;

namespace Seqeron.Genomics.Alignment;

// =====================================================================================================================
// Waterman-Eggert top-N local alignments in LINEAR space (audit A1-7c).
//
// The in-memory declumper (SequenceAligner.TopN.cs) keeps the O(m·n) matrices of the 1987 algorithm. Huang-Miller SIM
// (EMBOSS matcher) obtains the very same alignments in O(m + n) space by RECOMPUTING: every round is one run of the
// exact linear-space local program (SequenceAligner.LinearSpace.cs) on the matrix in which the aligned residue pairs of
// the earlier rounds have no diagonal move. The forbidden pairs are threaded through
//   * the locate passes (LocateLinearLocal / LocateAffine: the diagonal into a forbidden cell is removed before the cell is
//     maximised, so the first row-major maximum, its traceback tie order and its first aligned pair are those of the
//     matrix the in-memory declumper holds after the same rounds),
//   * the midpoint passes of the interval solve (LinearMidpoint / AffineMidpoint) and the full-matrix base cases of the
//     divide and conquer (LinearAlignCore / AffineAlignCore + their fill and traceback), all addressed by ABSOLUTE
//     coordinates (a sub-problem knows its row/column offset).
// The interval solve is a global alignment of the located rectangle; removing diagonal edges changes neither the proof
// that the carried traceback crossing is the one of the full-matrix traceback nor the identity of rows and ties
// (docs/algorithms/Alignment/Linear_Space_Alignment.md, "Waterman-Eggert in linear space"). Memory: O(m + n) per round
// plus the forbidden pairs (the total length of the reported alignments); time: N + 1 rounds of the local program.
// =====================================================================================================================

/// <summary>
/// The aligned residue pairs (1-based absolute row/column of sequence 1 / sequence 2) that Waterman-Eggert declumping has
/// forbidden: stored per row as sorted column arrays so that a DP row can be scanned with a forward-moving cursor
/// (<see cref="ForbiddenRow"/>) at no cost for rows without forbidden pairs.
/// </summary>
internal sealed class ForbiddenPairs
{
    /// <summary>Diagonal candidate of a forbidden cell: far below every reachable value, never equal to a cell value.</summary>
    internal const int NoDiagonal = int.MinValue / 2;

    private readonly Dictionary<int, int[]> _rows = new();

    public int Count { get; private set; }

    /// <summary>Adds pairs (row, column), both 1-based; duplicates are ignored.</summary>
    public void AddRange(IEnumerable<(int Row, int Column)> pairs)
    {
        var added = new Dictionary<int, List<int>>();
        foreach (var (row, col) in pairs)
        {
            if (!added.TryGetValue(row, out var list)) added[row] = list = new List<int>();
            list.Add(col);
        }

        foreach (var (row, list) in added)
        {
            if (_rows.TryGetValue(row, out var existing))
                list.AddRange(existing);
            list.Sort();
            var merged = new List<int>(list.Count);
            foreach (int c in list)
                if (merged.Count == 0 || merged[^1] != c)
                    merged.Add(c);
            Count += merged.Count - (existing?.Length ?? 0);
            _rows[row] = merged.ToArray();
        }
    }

    /// <summary>The sorted forbidden columns of a row, or null.</summary>
    public int[]? Row(int row) => _rows.TryGetValue(row, out var cols) ? cols : null;

    public bool Contains(int row, int col) =>
        _rows.TryGetValue(row, out var cols) && Array.BinarySearch(cols, col) >= 0;
}

/// <summary>Forward-only cursor over the forbidden columns of one DP row (columns must be tested in increasing order).</summary>
internal struct ForbiddenRow
{
    private readonly int[]? _cols;
    private int _next;

    public ForbiddenRow(ForbiddenPairs? forbidden, int absoluteRow)
    {
        _cols = forbidden?.Row(absoluteRow);
        _next = 0;
    }

    /// <summary>True when the pair (row, <paramref name="absoluteColumn"/>) is forbidden.</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public bool Test(int absoluteColumn)
    {
        if (_cols is null)
            return false;
        while (_next < _cols.Length && _cols[_next] < absoluteColumn)
            _next++;
        return _next < _cols.Length && _cols[_next] == absoluteColumn;
    }
}

public static partial class SequenceAligner
{
    #region Waterman-Eggert in linear space

    /// <summary>
    /// Waterman-Eggert rounds in O(m + n) memory: round k runs the exact linear-space local program with the aligned
    /// pairs of rounds 1..k-1 forbidden. Results are identical (rows, score, coordinates, order, ties) to the in-memory
    /// declumper <see cref="LocalTopNCore"/> on the same input.
    /// </summary>
    internal static IReadOnlyList<AlignmentResult> LocalTopNLinearSpaceCore(
        string a, string b, PairScoring scoring, bool affine, int maxAlignments, long minScaled,
        long baseCells, CancellationToken ct = default)
    {
        if (maxAlignments < 1)
            throw new ArgumentOutOfRangeException(nameof(maxAlignments), maxAlignments, "maxAlignments must be at least 1.");
        var results = new List<AlignmentResult>();
        if (a.Length == 0 || b.Length == 0)
            return results;
        scoring.EnsureRange(a.Length, b.Length);

        var forbidden = new ForbiddenPairs();
        while (results.Count < maxAlignments)
        {
            AlignmentCancellation.Poll(ct);
            var next = LinearSpaceAlign(a, b, scoring, affine, local: true, EndGaps.None, baseCells, ct, forbidden);
            if (next.Score <= 0 || next.Score < minScaled)
                break;
            results.Add(next);
            forbidden.AddRange(AlignedPairs(next));
        }
        return results;
    }

    private static IReadOnlyList<AlignmentResult> TopNLinearSpaceCore(
        string a, string b, ScoringMatrix scoring, bool affine, int maxAlignments, int minScore, long maxFullMatrixCells) =>
        LocalTopNLinearSpaceCore(a, b, new PairScoring(scoring), affine, maxAlignments, minScore,
            Math.Max(1, maxFullMatrixCells));

    private static IReadOnlyList<AlignmentResult> TopNLinearSpaceWithMatrix(
        string sequence1, string sequence2, MatrixScoring scoring, bool affine, int maxAlignments, double minScore, long maxFullMatrixCells)
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
        var raw = LocalTopNLinearSpaceCore(a, b, ps, affine, maxAlignments, minScaled, Math.Max(1, maxFullMatrixCells));
        var list = new List<AlignmentResult>(raw.Count);
        foreach (var r in raw)
        {
            double exact = (double)r.Score / ps.Scale;
            list.Add(r with { Score = (int)Math.Round(exact, MidpointRounding.AwayFromZero), ExactScore = exact });
        }
        return list;
    }

    /// <summary>
    /// Waterman-Eggert top-N local alignments (linear gaps) in exact <b>linear space</b> (O(m + n) memory per round instead of
    /// the O(m·n) matrices): the very same alignments, in the same order, with the same ties as
    /// <see cref="LocalAlignTopN(string, string, ScoringMatrix?, int, int)"/>, which switches to this program by itself above
    /// <see cref="TopNMaxCells"/> cells; call this method to force it for any size. Cost: one run of the linear-space local
    /// program per reported alignment (plus the final empty round).
    /// </summary>
    /// <param name="sequence1">First sequence (upper-cased internally).</param>
    /// <param name="sequence2">Second sequence.</param>
    /// <param name="scoring">Scoring (default <see cref="SimpleDna"/>).</param>
    /// <param name="maxAlignments">At most this many alignments (&gt;= 1).</param>
    /// <param name="minScore">Smallest accepted score (alignments scoring &lt;= 0 are never reported).</param>
    /// <param name="maxFullMatrixCells">Sub-problems of at most this many cells are solved with the full-matrix core
    /// (clamped to <see cref="LinearSpaceCellThreshold"/>); smaller values recurse deeper.</param>
    /// <exception cref="ArgumentException">A gap score is positive (rejected by every scoring model), or the integer score range is exceeded.</exception>
    [System.Runtime.CompilerServices.OverloadResolutionPriority(1)]
    public static IReadOnlyList<AlignmentResult> LocalAlignTopNLinearSpace(
        string sequence1, string sequence2, ScoringMatrix? scoring = null, int maxAlignments = 10, int minScore = 0,
        long maxFullMatrixCells = LinearSpaceBaseCells)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        return TopNLinearSpaceCore(sequence1.ToUpperInvariant(), sequence2.ToUpperInvariant(), scoring ?? SimpleDna, affine: false, maxAlignments, minScore, maxFullMatrixCells);
    }

    /// <summary><see cref="LocalAlignTopNLinearSpace(string, string, ScoringMatrix?, int, int, long)"/> for DNA sequences.</summary>
    [System.Runtime.CompilerServices.OverloadResolutionPriority(1)]
    public static IReadOnlyList<AlignmentResult> LocalAlignTopNLinearSpace(
        DnaSequence sequence1, DnaSequence sequence2, ScoringMatrix? scoring = null, int maxAlignments = 10, int minScore = 0,
        long maxFullMatrixCells = LinearSpaceBaseCells)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        return TopNLinearSpaceCore(sequence1.Sequence, sequence2.Sequence, scoring ?? SimpleDna, affine: false, maxAlignments, minScore, maxFullMatrixCells);
    }

    /// <summary>Linear-space Waterman-Eggert top-N with a substitution matrix and (possibly fractional) linear gap scores;
    /// identical to <see cref="LocalAlignTopN(string, string, MatrixScoring, int, double)"/>.</summary>
    public static IReadOnlyList<AlignmentResult> LocalAlignTopNLinearSpace(
        string sequence1, string sequence2, MatrixScoring scoring, int maxAlignments = 10, double minScore = 0,
        long maxFullMatrixCells = LinearSpaceBaseCells) =>
        TopNLinearSpaceWithMatrix(sequence1, sequence2, scoring, affine: false, maxAlignments, minScore, maxFullMatrixCells);

    /// <summary><see cref="LocalAlignTopNLinearSpace(string, string, MatrixScoring, int, double, long)"/> for protein sequences.</summary>
    public static IReadOnlyList<AlignmentResult> LocalAlignTopNLinearSpace(
        ProteinSequence sequence1, ProteinSequence sequence2, MatrixScoring scoring, int maxAlignments = 10, double minScore = 0,
        long maxFullMatrixCells = LinearSpaceBaseCells)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        return TopNLinearSpaceWithMatrix(sequence1.Sequence, sequence2.Sequence, scoring, affine: false, maxAlignments, minScore, maxFullMatrixCells);
    }

    /// <summary><see cref="LocalAlignTopNLinearSpace(string, string, MatrixScoring, int, double, long)"/> for DNA sequences.</summary>
    public static IReadOnlyList<AlignmentResult> LocalAlignTopNLinearSpace(
        DnaSequence sequence1, DnaSequence sequence2, MatrixScoring scoring, int maxAlignments = 10, double minScore = 0,
        long maxFullMatrixCells = LinearSpaceBaseCells)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        return TopNLinearSpaceWithMatrix(sequence1.Sequence, sequence2.Sequence, scoring, affine: false, maxAlignments, minScore, maxFullMatrixCells);
    }

    /// <summary>Waterman-Eggert top-N local alignments with <b>affine</b> gaps in exact linear space; the same alignments as
    /// <see cref="LocalAlignAffineTopN(string, string, ScoringMatrix?, int, int)"/> (which switches to this program above
    /// <see cref="TopNMaxCells"/> cells); see <see cref="LocalAlignTopNLinearSpace(string, string, ScoringMatrix?, int, int, long)"/>.</summary>
    [System.Runtime.CompilerServices.OverloadResolutionPriority(1)]
    public static IReadOnlyList<AlignmentResult> LocalAlignAffineTopNLinearSpace(
        string sequence1, string sequence2, ScoringMatrix? scoring = null, int maxAlignments = 10, int minScore = 0,
        long maxFullMatrixCells = LinearSpaceBaseCells)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        return TopNLinearSpaceCore(sequence1.ToUpperInvariant(), sequence2.ToUpperInvariant(), scoring ?? SimpleDna, affine: true, maxAlignments, minScore, maxFullMatrixCells);
    }

    /// <summary><see cref="LocalAlignAffineTopNLinearSpace(string, string, ScoringMatrix?, int, int, long)"/> for DNA sequences.</summary>
    [System.Runtime.CompilerServices.OverloadResolutionPriority(1)]
    public static IReadOnlyList<AlignmentResult> LocalAlignAffineTopNLinearSpace(
        DnaSequence sequence1, DnaSequence sequence2, ScoringMatrix? scoring = null, int maxAlignments = 10, int minScore = 0,
        long maxFullMatrixCells = LinearSpaceBaseCells)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        return TopNLinearSpaceCore(sequence1.Sequence, sequence2.Sequence, scoring ?? SimpleDna, affine: true, maxAlignments, minScore, maxFullMatrixCells);
    }

    /// <summary>Linear-space affine Waterman-Eggert top-N with a substitution matrix; identical to
    /// <see cref="LocalAlignAffineTopN(string, string, MatrixScoring, int, double)"/>.</summary>
    public static IReadOnlyList<AlignmentResult> LocalAlignAffineTopNLinearSpace(
        string sequence1, string sequence2, MatrixScoring scoring, int maxAlignments = 10, double minScore = 0,
        long maxFullMatrixCells = LinearSpaceBaseCells) =>
        TopNLinearSpaceWithMatrix(sequence1, sequence2, scoring, affine: true, maxAlignments, minScore, maxFullMatrixCells);

    /// <summary><see cref="LocalAlignAffineTopNLinearSpace(string, string, MatrixScoring, int, double, long)"/> for protein sequences.</summary>
    public static IReadOnlyList<AlignmentResult> LocalAlignAffineTopNLinearSpace(
        ProteinSequence sequence1, ProteinSequence sequence2, MatrixScoring scoring, int maxAlignments = 10, double minScore = 0,
        long maxFullMatrixCells = LinearSpaceBaseCells)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        return TopNLinearSpaceWithMatrix(sequence1.Sequence, sequence2.Sequence, scoring, affine: true, maxAlignments, minScore, maxFullMatrixCells);
    }

    /// <summary><see cref="LocalAlignAffineTopNLinearSpace(string, string, MatrixScoring, int, double, long)"/> for DNA sequences.</summary>
    public static IReadOnlyList<AlignmentResult> LocalAlignAffineTopNLinearSpace(
        DnaSequence sequence1, DnaSequence sequence2, MatrixScoring scoring, int maxAlignments = 10, double minScore = 0,
        long maxFullMatrixCells = LinearSpaceBaseCells)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        return TopNLinearSpaceWithMatrix(sequence1.Sequence, sequence2.Sequence, scoring, affine: true, maxAlignments, minScore, maxFullMatrixCells);
    }

    #endregion
}
