using System.Runtime.CompilerServices;
using System.Text;

namespace Seqeron.Genomics.Alignment;


/// <summary>
/// Internal exact-integer scoring model shared by every dynamic program: either the scalar
/// <see cref="ScoringMatrix"/> (scale 1) or a <see cref="MatrixScoring"/> converted to integers in units of 1/Scale.
/// </summary>
internal readonly struct PairScoring
{
    public readonly int Match;
    public readonly int Mismatch;
    public readonly int GapOpen;
    public readonly int GapExtend;

    /// <summary>Every score of this model is an integer multiple of 1/Scale of the reported (real) score.</summary>
    public readonly int Scale;

    /// <summary>True when the matrix alphabet contains a residue that is not a nucleotide code (a protein matrix): selects ClustalW's PROTEIN pairwise stage.</summary>
    public readonly bool Protein;

    private readonly int[]? _table;
    private readonly sbyte[]? _index;
    private readonly int _dim;

    /// <summary>Largest |score| of one step (a substitution or the gap open + extend of one gap position) in the integer units of this model.</summary>
    private readonly long _maxAbs;

    /// <summary>
    /// The scalar match / mismatch model. The gap scores follow the contract of <see cref="MatrixScoring"/>: <b>finite scores at most 0</b>
    /// (a penalty; a positive gap score is rejected exactly like in the matrix model), and the whole model must fit the integer dynamic
    /// program: a path of length L may not reach <c>2^28</c> in absolute value (see <see cref="EnsureRange(long, long)"/>, called by every
    /// aligner with the real sequence lengths), so no score can silently wrap around <c>int</c>.
    /// </summary>
    /// <exception cref="ArgumentException">A gap score is positive or <c>int.MinValue</c>.</exception>
    public PairScoring(ScoringMatrix scoring)
    {
        ArgumentNullException.ThrowIfNull(scoring);
        if (scoring.GapOpen > 0 || scoring.GapOpen == int.MinValue)
            throw new ArgumentException($"{nameof(ScoringMatrix.GapOpen)} must be a finite score <= 0 (a penalty) and greater than int.MinValue, got {scoring.GapOpen}.", nameof(scoring));
        if (scoring.GapExtend > 0 || scoring.GapExtend == int.MinValue)
            throw new ArgumentException($"{nameof(ScoringMatrix.GapExtend)} must be a finite score <= 0 (a penalty) and greater than int.MinValue, got {scoring.GapExtend}.", nameof(scoring));
        _maxAbs = Math.Max(Math.Max(Math.Abs((long)scoring.Match), Math.Abs((long)scoring.Mismatch)),
            Math.Abs((long)scoring.GapOpen) + Math.Abs((long)scoring.GapExtend));
        Match = scoring.Match;
        Mismatch = scoring.Mismatch;
        GapOpen = scoring.GapOpen;
        GapExtend = scoring.GapExtend;
        Scale = 1;
        _table = null;
        _index = null;
        _dim = 0;
    }

    private PairScoring(int gapOpen, int gapExtend, int scale, int[] table, sbyte[] index, int dim, bool protein, long maxAbs)
    {
        _maxAbs = maxAbs;
        Protein = protein;
        Match = 0;
        Mismatch = 0;
        GapOpen = gapOpen;
        GapExtend = gapExtend;
        Scale = scale;
        _table = table;
        _index = index;
        _dim = dim;
    }

    public static implicit operator PairScoring(ScoringMatrix scoring) => new(scoring);

    /// <summary>True for a <see cref="MatrixScoring"/> model (any residue alphabet), false for the scalar match/mismatch model.</summary>
    public bool IsMatrix => _table is not null;

    /// <summary>Score of the aligned pair (a, b) (upper-case symbols already validated for matrix models).</summary>
    public int Sub(char a, char b)
    {
        if (_table is not null)
            return _table[_index![a] * _dim + _index[b]];
        return a == b ? Match : Mismatch;
    }

    /// <summary>Largest decimal precision of a gap score (1e-6, the resolution of Biopython's/EMBOSS's float comparisons).</summary>
    private const int MaxGapDecimals = 6;

    /// <summary>Headroom: |score| of any path must stay below this (the affine −∞ sentinel is int.MinValue / 4).</summary>
    private const long MaxPathScore = 1L << 28;

    /// <summary>
    /// Builds the exact integer model of <paramref name="scoring"/> for the (upper-case) sequences <paramref name="a"/>
    /// and <paramref name="b"/>, validating residues, gap scores and the integer range.
    /// </summary>
    public static PairScoring Create(MatrixScoring scoring, string a, string b)
    {
        var matrix = scoring.Matrix ?? throw new ArgumentException("MatrixScoring.Matrix must not be null.", nameof(scoring));
        ValidateResidues(matrix, a, "sequence1");
        ValidateResidues(matrix, b, "sequence2");
        return Build(scoring, matrix, (long)a.Length + b.Length + 2, 1,
            $"sequence lengths {a.Length} and {b.Length}");
    }

    /// <summary>
    /// Builds the exact integer model of <paramref name="scoring"/> for a SET of (upper-case) sequences, e.g. the multiple aligners: every
    /// sequence is validated against the matrix alphabet and the integer range is checked for the sum of pairs of all C(k, 2) pairs.
    /// </summary>
    public static PairScoring Create(MatrixScoring scoring, IReadOnlyList<string> sequences)
    {
        var matrix = scoring.Matrix ?? throw new ArgumentException("MatrixScoring.Matrix must not be null.", nameof(scoring));
        long longest = 0, second = 0;
        for (int i = 0; i < sequences.Count; i++)
        {
            ValidateResidues(matrix, sequences[i], $"sequence{i + 1}");
            long len = sequences[i].Length;
            if (len > longest) { second = longest; longest = len; }
            else if (len > second) second = len;
        }

        long k = sequences.Count;
        return Build(scoring, matrix, longest + second + 2, Math.Max(1, k * (k - 1) / 2),
            $"{sequences.Count} sequences, the two longest have {longest} and {second} residues");
    }

    private static PairScoring Build(MatrixScoring scoring, SubstitutionMatrix matrix, long pathLength, long pairs, string sizeText)
    {
        decimal open = ToGapDecimal(scoring.GapOpen, nameof(MatrixScoring.GapOpen));
        decimal extend = ToGapDecimal(scoring.GapExtend, nameof(MatrixScoring.GapExtend));

        // Smallest integer scale that makes both gap scores integral (0.5 -> 2, 0.1 -> 10, integers -> 1).
        int decimals = 0;
        while (decimals <= MaxGapDecimals && !(IsIntegral(open * Pow10(decimals)) && IsIntegral(extend * Pow10(decimals))))
            decimals++;
        if (decimals > MaxGapDecimals)
            throw new ArgumentException($"Gap scores are resolved to 1e-{MaxGapDecimals}: {scoring.GapOpen}, {scoring.GapExtend}.", nameof(scoring));

        long denominator = Pow10(decimals);
        long g = Gcd(denominator, Gcd(Math.Abs((long)(open * denominator)), Math.Abs((long)(extend * denominator))));
        long scale = denominator / g;
        long gapOpen = (long)(open * denominator) / g;
        long gapExtend = (long)(extend * denominator) / g;

        long maxAbs = Math.Max((long)matrix.MaxAbsScore * scale, Math.Abs(gapOpen) + Math.Abs(gapExtend));
        if (maxAbs * pathLength * pairs >= MaxPathScore)
            throw RangeError(scale, sizeText, nameof(scoring));

        int dim = matrix.Alphabet.Length;
        var index = new sbyte[128];
        Array.Fill(index, (sbyte)-1);
        for (int i = 0; i < dim; i++)
            index[matrix.Alphabet[i]] = (sbyte)i;

        var table = new int[dim * dim];
        for (int i = 0; i < dim; i++)
            for (int j = 0; j < dim; j++)
                table[i * dim + j] = checked((int)(matrix.ScoreAt(i, j) * scale));

        return new PairScoring((int)gapOpen, (int)gapExtend, (int)scale, table, index, dim, !IsNucleotideAlphabet(matrix.Alphabet), maxAbs);
    }

    private static ArgumentException RangeError(long scale, string sizeText, string paramName) =>
        new($"Score range too large for the integer dynamic program (scale {scale}, {sizeText}); use coarser gap scores or shorter sequences.", paramName);

    /// <summary>
    /// The one integer-range check of every dynamic program (scalar and matrix model alike): the score of any path through sequences of
    /// <paramref name="length1"/> and <paramref name="length2"/> residues, <c>max|step| · (length1 + length2 + 2)</c>, must stay below
    /// 2^28 (the affine −∞ sentinel is <c>int.MinValue / 4</c>), otherwise the DP would silently wrap around. The matrix model has been
    /// checked for the same bound (and the whole set of pairs) by <see cref="Create(MatrixScoring, string, string)"/>; this call repeats the pairwise
    /// part at the point where the DP is about to run, so the scalar <see cref="ScoringMatrix"/> is guarded identically.
    /// </summary>
    /// <exception cref="ArgumentException">"Score range too large for the integer dynamic program ...".</exception>
    internal void EnsureRange(long length1, long length2)
    {
        if (_maxAbs * (length1 + length2 + 2) >= MaxPathScore)
            throw RangeError(Scale, $"sequence lengths {length1} and {length2}", "scoring");
    }

    private static bool IsNucleotideAlphabet(string alphabet)
    {
        foreach (char c in alphabet)
            if ("ACGTUNRYSWKMBDHV".IndexOf(c) < 0)
                return false;
        return true;
    }

    internal static void ValidateResidues(SubstitutionMatrix matrix, string seq, string name)
    {
        for (int i = 0; i < seq.Length; i++)
            if (matrix.IndexOf(seq[i]) < 0)
                throw new ArgumentException(
                    $"Residue '{seq[i]}' at position {i} of {name} is not in substitution matrix {matrix.Name} ({matrix.Alphabet}).", name);
    }

    private static decimal ToGapDecimal(double value, string name)
    {
        if (double.IsNaN(value) || double.IsInfinity(value) || value > 0)
            throw new ArgumentException($"{name} must be a finite score <= 0 (a penalty), got {value}.", name);
        return (decimal)value;
    }

    private static bool IsIntegral(decimal value) => value == decimal.Truncate(value);

    private static long Pow10(int n)
    {
        long r = 1;
        for (int i = 0; i < n; i++) r *= 10;
        return r;
    }

    private static long Gcd(long x, long y)
    {
        while (y != 0) (x, y) = (y, x % y);
        return x;
    }
}

/// <summary>
/// Performs pairwise sequence alignment using dynamic programming algorithms.
/// Supports global (Needleman-Wunsch), local (Smith-Waterman), and semi-global alignment.
/// </summary>
public static partial class SequenceAligner
{
    #region Scoring Matrices

    /// <summary>
    /// Simple DNA scoring: +1 match, -1 mismatch.
    /// </summary>
    public static readonly ScoringMatrix SimpleDna = new(
        Match: 1,
        Mismatch: -1,
        GapOpen: -2,
        GapExtend: -1);

    /// <summary>
    /// NCBI BLAST+ <c>blastn</c>-task default DNA scoring: reward +2, penalty -3, gap existence 5,
    /// gap extension 2 (BLAST+ User Manual, blastn application options). BLAST charges a gap of
    /// length k as existence + k·extension, which is exactly the <see cref="ScoringMatrix"/> convention
    /// used by <see cref="GlobalAlignAffine(DnaSequence, DnaSequence, ScoringMatrix?)"/> and
    /// <see cref="LocalAlignAffine(DnaSequence, DnaSequence, ScoringMatrix?)"/> (GapOpen + k·GapExtend). The linear-gap aligners (<c>GlobalAlign</c>, <c>LocalAlign</c>,
    /// <c>SemiGlobalAlign</c>) use only <see cref="ScoringMatrix.GapExtend"/> (-2 per gap position).
    /// </summary>
    public static readonly ScoringMatrix BlastDna = new(
        Match: 2,
        Mismatch: -3,
        GapOpen: -5,
        GapExtend: -2);

    /// <summary>
    /// High identity DNA scoring for closely related sequences.
    /// </summary>
    public static readonly ScoringMatrix HighIdentityDna = new(
        Match: 5,
        Mismatch: -4,
        GapOpen: -10,
        GapExtend: -1);

    #endregion

    #region Global Alignment (Needleman-Wunsch)

    /// <summary>
    /// Performs global alignment using the Needleman-Wunsch algorithm.
    /// Aligns entire sequences end-to-end.
    /// </summary>
    /// <remarks>
    /// Linear gap model (Needleman &amp; Wunsch 1970): every gap position scores
    /// <see cref="ScoringMatrix.GapExtend"/>; <see cref="ScoringMatrix.GapOpen"/> is ignored.
    /// For affine gap costs (GapOpen + k·GapExtend) use
    /// <see cref="GlobalAlignAffine(DnaSequence, DnaSequence, ScoringMatrix?)"/>.
    /// Empty/null contract (identical for every <c>GlobalAlign</c>/<c>GlobalAlignAffine</c> overload, string or
    /// <see cref="DnaSequence"/>): null → <see cref="ArgumentNullException"/>; an empty side is not special-cased —
    /// the Needleman-Wunsch optimum is the all-gap alignment (<c>""</c> vs <c>ACGT</c> → <c>----</c>/<c>ACGT</c>, score
    /// 4·GapExtend, as Biopython <c>PairwiseAligner.score</c> and EMBOSS needle); two empty inputs score 0 with empty rows.
    /// </remarks>
    /// <param name="sequence1">First DNA sequence.</param>
    /// <param name="sequence2">Second DNA sequence.</param>
    /// <param name="scoring">Scoring matrix (default: SimpleDna).</param>
    /// <returns>Alignment result with aligned sequences and score.</returns>
    [System.Runtime.CompilerServices.OverloadResolutionPriority(1)]
    public static AlignmentResult GlobalAlign(
        DnaSequence sequence1,
        DnaSequence sequence2,
        ScoringMatrix? scoring = null)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);

        return GlobalAlignCore(sequence1.Sequence, sequence2.Sequence, scoring ?? SimpleDna);
    }

    /// <summary>
    /// Performs global alignment on raw sequence strings.
    /// </summary>
    [System.Runtime.CompilerServices.OverloadResolutionPriority(1)]
    public static AlignmentResult GlobalAlign(
        string sequence1,
        string sequence2,
        ScoringMatrix? scoring = null)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);

        return GlobalAlignCore(
            sequence1.ToUpperInvariant(),
            sequence2.ToUpperInvariant(),
            scoring ?? SimpleDna);
    }

    /// <summary>
    /// Performs global alignment with cancellation support.
    /// Recommended for aligning long sequences.
    /// </summary>
    /// <param name="sequence1">First sequence string.</param>
    /// <param name="sequence2">Second sequence string.</param>
    /// <param name="scoring">Scoring matrix (default: SimpleDna).</param>
    /// <param name="cancellationToken">Cancellation token for long-running operations.</param>
    /// <param name="progress">Optional progress reporter (0.0 to 1.0).</param>
    /// <returns>Alignment result with aligned sequences and score.</returns>
    /// <remarks>Same linear-gap Needleman-Wunsch dynamic program as
    /// <see cref="GlobalAlign(string, string, ScoringMatrix?)"/> (shared core), with periodic
    /// cancellation checks and progress reports.</remarks>
    public static AlignmentResult GlobalAlign(
        string sequence1,
        string sequence2,
        ScoringMatrix? scoring,
        CancellationToken cancellationToken,
        IProgress<double>? progress = null)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);

        return GlobalAlignCore(
            sequence1.ToUpperInvariant(),
            sequence2.ToUpperInvariant(),
            scoring ?? SimpleDna,
            cancellationToken,
            progress);
    }

    /// <summary>
    /// Performs global alignment on DNA sequences with cancellation support.
    /// </summary>
    public static AlignmentResult GlobalAlign(
        DnaSequence sequence1,
        DnaSequence sequence2,
        ScoringMatrix? scoring,
        CancellationToken cancellationToken,
        IProgress<double>? progress = null)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);

        return GlobalAlignCore(sequence1.Sequence, sequence2.Sequence, scoring ?? SimpleDna, cancellationToken, progress);
    }

    /// <summary>
    /// Linear-gap Needleman-Wunsch core shared by every <c>GlobalAlign</c> overload.
    /// Gap penalty d = <see cref="ScoringMatrix.GapExtend"/> per gap position;
    /// <see cref="ScoringMatrix.GapOpen"/> is not used (see <see cref="GlobalAlignAffine(string, string, ScoringMatrix?)"/>).
    /// </summary>
    private static AlignmentResult GlobalAlignCore(
        string seq1,
        string seq2,
        PairScoring scoring,
        CancellationToken cancellationToken = default,
        IProgress<double>? progress = null) =>
        LinearAlignCore(seq1, seq2, scoring, EndGaps.None, cancellationToken, progress);

    /// <summary>
    /// Linear-gap Needleman-Wunsch dynamic program shared by <c>GlobalAlign</c> and
    /// <c>SemiGlobalAlign</c> (every end-gap policy). d = <see cref="ScoringMatrix.GapExtend"/> per gap position.
    /// <list type="bullet">
    /// <item><b>Global</b> (<see cref="EndGaps.None"/>): F(i,0) = d·i, F(0,j) = d·j; traceback from (m, n).</item>
    /// <item><b>Ends-free</b>: F(i,0) = 0 when <see cref="EndGaps.FreeStart1"/> (a free prefix of sequence 1),
    /// F(0,j) = 0 when <see cref="EndGaps.FreeStart2"/>; the alignment ends at the first maximum of the last row
    /// (<see cref="EndGaps.FreeEnd2"/>, smallest j) and/or of the last column (<see cref="EndGaps.FreeEnd1"/>, smallest i;
    /// row candidates first when both are free). The classic fitting alignment is
    /// <see cref="EndGaps.Fitting"/>.</item>
    /// </list>
    /// Traceback ties are broken diagonal &gt; up &gt; left.
    /// </summary>
    internal static AlignmentResult LinearAlignCore(
        string seq1,
        string seq2,
        PairScoring scoring,
        EndGaps ends,
        CancellationToken cancellationToken = default,
        IProgress<double>? progress = null,
        bool forceFullMatrix = false,
        ForbiddenPairs? forbidden = null,
        int rowOffset = 0,
        int colOffset = 0)
    {
        int m = seq1.Length;
        int n = seq2.Length;
        scoring.EnsureRange(m, n);

        AlignmentCancellation.Poll(cancellationToken);

        // Above LinearSpaceCellThreshold cells the O(m·n) matrix is replaced by the exact linear-space program
        // (same rows, same score; see the Linear-space region).
        if (!forceFullMatrix && ExceedsFullMatrixLimit(m, n))
        {
            var big = LinearSpaceAlign(seq1, seq2, scoring, affine: false, local: false, ends, LinearSpaceBaseCells, cancellationToken);
            progress?.Report(1.0);
            return big;
        }

        var score = LinearFillMatrix(seq1, seq2, scoring, ends, cancellationToken, progress, forbidden, rowOffset, colOffset);

        AlignmentCancellation.Poll(cancellationToken);
        progress?.Report(0.75);

        // End cell: (m, n) for a global alignment; with free trailing residues the first (smallest j / smallest i)
        // maximum of the last row (sequence-2 suffix free) and/or last column (sequence-1 suffix free).
        int endI = m, endJ = n;
        if (ends.HasFlag(EndGaps.FreeEnd2) || ends.HasFlag(EndGaps.FreeEnd1))
        {
            int bestValue = int.MinValue;
            foreach (var (ci, cj) in EndCellCandidates(m, n, ends))
                if (score[ci, cj] > bestValue) { bestValue = score[ci, cj]; endI = ci; endJ = cj; }
        }

        var result = Traceback(seq1, seq2, score, endI, endJ, scoring, ends, forbidden, rowOffset, colOffset);

        progress?.Report(1.0);
        return result;
    }

    /// <summary>
    /// Border initialisation and fill of the linear-gap Needleman-Wunsch matrix (the one dynamic program of every linear global /
    /// ends-free alignment; also the matrix the co-optimal enumerator walks back):
    /// F(i,0) = d·i (0 when <see cref="EndGaps.FreeStart1"/>), F(0,j) = d·j (0 when <see cref="EndGaps.FreeStart2"/>),
    /// F(i,j) = max(F(i-1,j-1) + s(a_i,b_j), F(i-1,j) + d, F(i,j-1) + d), d = <see cref="ScoringMatrix.GapExtend"/>.
    /// Standard Needleman-Wunsch (<see href="https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm"/>); ends-free: a free
    /// prefix of seq1 / seq2 may precede the alignment.
    /// </summary>
    private static int[,] LinearFillMatrix(
        string seq1, string seq2, PairScoring scoring, EndGaps ends,
        CancellationToken cancellationToken = default, IProgress<double>? progress = null,
        ForbiddenPairs? forbidden = null, int rowOffset = 0, int colOffset = 0)
    {
        int m = seq1.Length;
        int n = seq2.Length;
        var score = new int[m + 1, n + 1];

        for (int i = 0; i <= m; i++)
            score[i, 0] = ends.HasFlag(EndGaps.FreeStart1) ? 0 : i * scoring.GapExtend;
        for (int j = 0; j <= n; j++)
            score[0, j] = ends.HasFlag(EndGaps.FreeStart2) ? 0 : j * scoring.GapExtend;

        // F(i,j) = max(F(i-1,j-1)+S(a_i,b_j), F(i-1,j)+d, F(i,j-1)+d)
        for (int i = 1; i <= m; i++)
        {
            if (i % 100 == 0)
            {
                AlignmentCancellation.Poll(cancellationToken);
                progress?.Report((double)i / m * 0.5); // First half is matrix fill
            }

            char c1 = seq1[i - 1];
            var forbiddenRow = new ForbiddenRow(forbidden, rowOffset + i);
            for (int j = 1; j <= n; j++)
            {
                int matchScore = scoring.Sub(c1, seq2[j - 1]);

                // A forbidden aligned pair (Waterman-Eggert declumping, linear-space top-N) has no diagonal move.
                int diag = forbiddenRow.Test(colOffset + j) ? ForbiddenPairs.NoDiagonal : score[i - 1, j - 1] + matchScore;
                int up = score[i - 1, j] + scoring.GapExtend;
                int left = score[i, j - 1] + scoring.GapExtend;

                score[i, j] = Math.Max(diag, Math.Max(up, left));
            }
        }

        return score;
    }

    /// <summary>
    /// The cells where a global / ends-free alignment may end, in the fixed tie order: (m, n) when no trailing end is free,
    /// otherwise the last row from j = 0 (<see cref="EndGaps.FreeEnd2"/>: trailing sequence-2 residues free) and then the last
    /// column from i = 0 (<see cref="EndGaps.FreeEnd1"/>); (m, n) is listed twice when both are free (callers keep the first).
    /// </summary>
    private static IEnumerable<(int I, int J)> EndCellCandidates(int m, int n, EndGaps ends)
    {
        bool freeEnd2 = ends.HasFlag(EndGaps.FreeEnd2);
        bool freeEnd1 = ends.HasFlag(EndGaps.FreeEnd1);
        if (!freeEnd2 && !freeEnd1)
        {
            yield return (m, n);
            yield break;
        }

        if (freeEnd2)
            for (int j = 0; j <= n; j++)
                yield return (m, j);
        if (freeEnd1)
            for (int i = 0; i <= m; i++)
                yield return (i, n);
    }

    #endregion

    #region Global Alignment with affine gaps (Gotoh)

    // Sentinel for "state unreachable" (-infinity). Kept far from int.MinValue so that adding a
    // bounded number of penalties can never overflow (the classic "±inf header overflows on the
    // first iteration" implementation error).
    private const int AffineNegInf = int.MinValue / 4;

    private const byte StateMatch = 0; // M: a_i aligned to b_j
    private const byte StateGapInSeq2 = 1; // X: a_i aligned to '-' (vertical move)
    private const byte StateGapInSeq1 = 2; // Y: '-' aligned to b_j (horizontal move)

    /// <summary>
    /// Performs optimal global alignment with <b>affine gap costs</b> (Gotoh 1982, three-state
    /// dynamic program with the corrected initialization of Flouri et al. 2015).
    /// </summary>
    /// <param name="sequence1">First DNA sequence.</param>
    /// <param name="sequence2">Second DNA sequence.</param>
    /// <param name="scoring">Scoring matrix (default: <see cref="SimpleDna"/>).</param>
    /// <remarks>
    /// <para>
    /// Gap-cost convention (Gotoh 1982 w(k) = v + u·k; NCBI BLAST "existence + extension"):
    /// a gap of length k scores <c>GapOpen + k·GapExtend</c>, i.e. <see cref="ScoringMatrix.GapOpen"/>
    /// does <b>not</b> include the first extension. Equivalent settings in other tools:
    /// Biopython <c>PairwiseAligner(open_gap_score = GapOpen + GapExtend, extend_gap_score = GapExtend)</c>;
    /// parasail <c>nw(open = -(GapOpen + GapExtend), extend = -GapExtend)</c>; BLAST+
    /// <c>-gapopen -GapOpen -gapextend -GapExtend</c>. With <c>GapOpen = 0</c> the result equals the
    /// linear-gap <see cref="GlobalAlign(DnaSequence, DnaSequence, ScoringMatrix?)"/>.
    /// </para>
    /// <para>
    /// Recurrences (M = ends in a substitution, X = ends in a gap in sequence 2, Y = ends in a gap in
    /// sequence 1; o = GapOpen, e = GapExtend):
    /// M(i,j) = max(M,X,Y)(i-1,j-1) + s(a_i,b_j);
    /// X(i,j) = max(M(i-1,j) + o + e, X(i-1,j) + e, Y(i-1,j) + o + e);
    /// Y(i,j) = max(M(i,j-1) + o + e, Y(i,j-1) + e, X(i,j-1) + o + e).
    /// Initialization: M(0,0) = 0; X(i,0) = o + i·e, Y(0,j) = o + j·e; all other border cells −∞
    /// — a leading/trailing end gap is charged exactly one opening (Flouri et al. 2015 correction of
    /// Gotoh's border initialization). The traceback follows the state that produced each value
    /// (switching matrices explicitly), so it always reproduces an alignment with the optimal score.
    /// Ties are broken M &gt; X &gt; Y.
    /// </para>
    /// <para>Sources: Gotoh O (1982) J Mol Biol 162:705-708; Flouri T, Kobert K, Rognes T, Stamatakis A
    /// (2015) bioRxiv 10.1101/031500; cross-checked against Biopython 1.88 PairwiseAligner and
    /// parasail 1.3.4 nw.</para>
    /// </remarks>
    [System.Runtime.CompilerServices.OverloadResolutionPriority(1)]
    public static AlignmentResult GlobalAlignAffine(
        DnaSequence sequence1,
        DnaSequence sequence2,
        ScoringMatrix? scoring = null)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);

        return GlobalAlignAffineCore(sequence1.Sequence, sequence2.Sequence, scoring ?? SimpleDna);
    }

    /// <summary>
    /// Affine-gap global alignment (Gotoh) on raw sequence strings (uppercased internally).
    /// See <see cref="GlobalAlignAffine(DnaSequence, DnaSequence, ScoringMatrix?)"/> for the model.
    /// </summary>
    [System.Runtime.CompilerServices.OverloadResolutionPriority(1)]
    public static AlignmentResult GlobalAlignAffine(
        string sequence1,
        string sequence2,
        ScoringMatrix? scoring = null)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);

        return GlobalAlignAffineCore(
            sequence1.ToUpperInvariant(),
            sequence2.ToUpperInvariant(),
            scoring ?? SimpleDna);
    }

    private static AlignmentResult GlobalAlignAffineCore(string seq1, string seq2, PairScoring scoring)
    {
        return AffineAlignCore(seq1, seq2, scoring, AlignmentType.Global);
    }

    /// <summary>
    /// Three-state (Gotoh) affine-gap dynamic program shared by <c>GlobalAlignAffine</c>,
    /// <c>LocalAlignAffine</c> and <c>SemiGlobalAlignAffine</c>. o = GapOpen, e = GapExtend; a gap of
    /// length k scores o + k·e.
    /// <list type="bullet">
    /// <item><b>Global</b>: Flouri et al. (2015) border initialization, traceback from the best of
    /// M/X/Y at (m, n) to the origin.</item>
    /// <item><b>Local</b> (Smith-Waterman-Gotoh): every border cell is −∞ and the substitution state
    /// may start a new alignment from score 0, M(i,j) = max(0, M, X, Y)(i-1,j-1) + s(a_i,b_j) — the
    /// zero floor of Smith &amp; Waterman (1981) applied to Gotoh's (1982) recurrences. The optimal
    /// local alignment ends at the first (row-major) cell holding the maximum of M, and the traceback
    /// stops as soon as the predecessor value is 0.</item>
    /// <item><b>SemiGlobal</b> (fitting, sequence 1 fitted into sequence 2): M(0,j) = 0 for every j
    /// (a free prefix of sequence 2 may precede the alignment), X(i,0) = o + i·e, all other border
    /// cells −∞; the alignment ends in the first (smallest j) column of the last row holding
    /// max(M, X, Y)(m, j) (free sequence-2 suffix). The traceback stops on reaching row 0 and the
    /// skipped sequence-2 prefix/suffix is emitted against gaps.</item>
    /// </list>
    /// Traceback ties are broken M &gt; X &gt; Y.
    /// </summary>
    internal static AlignmentResult AffineAlignCore(
        string seq1, string seq2, PairScoring scoring, AlignmentType type, EndGaps ends = EndGaps.None,
        int startState = StateMatch, int endState = -1, bool forceFullMatrix = false,
        ForbiddenPairs? forbidden = null, int rowOffset = 0, int colOffset = 0)
    {
        bool local = type == AlignmentType.Local;
        if (local) ends = EndGaps.None;
        scoring.EnsureRange(seq1.Length, seq2.Length);
        if (!forceFullMatrix && ExceedsFullMatrixLimit(seq1.Length, seq2.Length))
            return LinearSpaceAlign(seq1, seq2, scoring, affine: true, local, ends, LinearSpaceBaseCells, default);
        bool freeEnd1 = ends.HasFlag(EndGaps.FreeEnd1);
        bool freeEnd2 = ends.HasFlag(EndGaps.FreeEnd2);
        int m = seq1.Length;
        int n = seq2.Length;

        var (mm, xx, yy) = AffineNewMatrices(m, n);
        AffineInitBorders(mm, xx, yy, scoring, local, ends, startState);
        var (bestLocal, bestI, bestJ) = AffineFill(mm, xx, yy, seq1, seq2, scoring, local, forbidden, rowOffset, colOffset);

        int best, ci, cj;
        byte state;
        if (local)
        {
            if (bestLocal == 0)
                return EmptyLocalResult();
            best = bestLocal;
            ci = bestI;
            cj = bestJ;
            state = StateMatch;
        }
        else
        {
            // Global: end at (m, n). Ends-free: the first maximum (strict >) of max(M, X, Y) over the last row
            // (free trailing seq2 residues, smallest j) and/or the last column (free trailing seq1 residues, smallest i;
            // row candidates first). Optimal end state ties: M > X > Y.
            ci = m;
            cj = n;
            if (freeEnd2 || freeEnd1)
            {
                int bestValue = int.MinValue;
                foreach (var (ei, ej) in EndCellCandidates(m, n, ends))
                {
                    int v = Max3(mm[ei, ej], xx[ei, ej], yy[ei, ej]);
                    if (v > bestValue) { bestValue = v; ci = ei; cj = ej; }
                }
            }
            best = Max3(mm[ci, cj], xx[ci, cj], yy[ci, cj]);
            state = PickState(best, mm[ci, cj], StateMatch, xx[ci, cj], StateGapInSeq2, StateGapInSeq1);
            if (endState >= 0)
            {
                // Forced end state (the linear-space program's sub-problems end in the state of the split vertex).
                state = (byte)endState;
                best = endState switch
                {
                    StateMatch => mm[ci, cj],
                    StateGapInSeq2 => xx[ci, cj],
                    _ => yy[ci, cj],
                };
            }
        }

        return AffineTraceback(seq1, seq2, mm, xx, yy, scoring, local, ends, ci, cj, state, best);
    }

    private static (int[,] M, int[,] X, int[,] Y) AffineNewMatrices(int m, int n) =>
        (new int[m + 1, n + 1], new int[m + 1, n + 1], new int[m + 1, n + 1]);

    /// <summary>
    /// Border initialisation of the three Gotoh matrices (shared by every affine program: global, ends-free, local and the
    /// co-optimal enumerator). Origin vector: the path leaves (0,0) from the state it arrived in (M for a whole alignment; the
    /// linear-space program starts the sub-problems below the split row from the state of the split vertex, <paramref name="startState"/>).
    /// Global: Flouri et al. (2015) — X(i,0) = o + i·e, Y(0,j) = o + j·e (a leading end gap is charged exactly one opening), other
    /// border cells −∞; a free start (<see cref="EndGaps.FreeStart1"/> / <see cref="EndGaps.FreeStart2"/>) sets M of the matching
    /// border to 0 and the gap state of that border to −∞; local: every border cell is −∞.
    /// </summary>
    private static void AffineInitBorders(
        int[,] mm, int[,] xx, int[,] yy, PairScoring scoring, bool local, EndGaps ends, int startState)
    {
        int m = mm.GetLength(0) - 1;
        int n = mm.GetLength(1) - 1;
        int e = scoring.GapExtend;
        int oe = scoring.GapOpen + e;
        bool freeStart1 = ends.HasFlag(EndGaps.FreeStart1);
        bool freeStart2 = ends.HasFlag(EndGaps.FreeStart2);

        mm[0, 0] = !local && startState == StateMatch ? 0 : AffineNegInf;
        xx[0, 0] = !local && startState == StateGapInSeq2 ? 0 : AffineNegInf;
        yy[0, 0] = !local && startState == StateGapInSeq1 ? 0 : AffineNegInf;
        for (int i = 1; i <= m; i++)
        {
            mm[i, 0] = freeStart1 ? 0 : AffineNegInf;
            // X(i,0) = o + i·e from the default origin (a leading gap is charged exactly one opening, Flouri 2015).
            xx[i, 0] = local || freeStart1
                ? AffineNegInf
                : Math.Max(AffineNegInf, Max3(mm[i - 1, 0] + oe, xx[i - 1, 0] + e, yy[i - 1, 0] + oe));
            yy[i, 0] = AffineNegInf;
        }
        for (int j = 1; j <= n; j++)
        {
            mm[0, j] = freeStart2 ? 0 : AffineNegInf;
            xx[0, j] = AffineNegInf;
            yy[0, j] = local || freeStart2
                ? AffineNegInf
                : Math.Max(AffineNegInf, Max3(mm[0, j - 1] + oe, yy[0, j - 1] + e, xx[0, j - 1] + oe));
        }
    }

    /// <summary>
    /// Fills the interior of the three Gotoh matrices row by row with <see cref="AffineCell"/> (local: the zero floor of
    /// Smith-Waterman); returns the first (row-major) maximum of M for a local alignment (0, 0, 0 otherwise / when no positive cell).
    /// </summary>
    private static (int Best, int I, int J) AffineFill(
        int[,] mm, int[,] xx, int[,] yy, string seq1, string seq2, PairScoring scoring, bool local,
        ForbiddenPairs? forbidden = null, int rowOffset = 0, int colOffset = 0)
    {
        int m = seq1.Length;
        int n = seq2.Length;
        int e = scoring.GapExtend;
        int oe = scoring.GapOpen + e;
        int startFloor = local ? 0 : AffineNegInf;
        int bestLocal = 0, bestI = 0, bestJ = 0;
        for (int i = 1; i <= m; i++)
        {
            AlignmentCancellation.ThrowIfRequested();
            char c1 = seq1[i - 1];
            var forbiddenRow = new ForbiddenRow(forbidden, rowOffset + i);
            for (int j = 1; j <= n; j++)
            {
                AffineCell(mm, xx, yy, i, j, scoring.Sub(c1, seq2[j - 1]), startFloor, e, oe, forbiddenRow.Test(colOffset + j));

                if (local && mm[i, j] > bestLocal)
                {
                    bestLocal = mm[i, j];
                    bestI = i;
                    bestJ = j;
                }
            }
        }

        return (bestLocal, bestI, bestJ);
    }

    /// <summary>
    /// Gotoh traceback from state <paramref name="state"/> at cell (<paramref name="ci"/>, <paramref name="cj"/>) of the
    /// filled matrices; ties M &gt; X &gt; Y, a local alignment stops when the predecessor value is 0. Shared by
    /// <see cref="AffineAlignCore"/> and the Waterman-Eggert declumper.
    /// </summary>
    private static AlignmentResult AffineTraceback(
        string seq1, string seq2, int[,] mm, int[,] xx, int[,] yy, PairScoring scoring,
        bool local, EndGaps ends, int ci, int cj, byte state, int best)
    {
        int m = seq1.Length;
        int n = seq2.Length;
        int e = scoring.GapExtend;
        int oe = scoring.GapOpen + e;
        bool freeStart1 = ends.HasFlag(EndGaps.FreeStart1);
        bool freeStart2 = ends.HasFlag(EndGaps.FreeStart2);
        int endI = ci;
        int endJ = cj;
        var chars1 = new List<char>(m + n);
        var chars2 = new List<char>(m + n);
        if (!local)
        {
            // Free trailing overhang (emitted reversed; the whole buffer is reversed below). At most one side is non-empty.
            for (int k = n; k > endJ; k--)
            {
                chars1.Add('-');
                chars2.Add(seq2[k - 1]);
            }
            for (int k = m; k > endI; k--)
            {
                chars1.Add(seq1[k - 1]);
                chars2.Add('-');
            }
        }

        // Stop at the origin, in row 0 when the sequence-2 start is free, in column 0 when the sequence-1 start is free.
        while (!(ci == 0 && (cj == 0 || freeStart2)) && !(cj == 0 && freeStart1))
        {
            if (state == StateMatch)
            {
                int prev = mm[ci, cj] - (scoring.Sub(seq1[ci - 1], seq2[cj - 1]));
                chars1.Add(seq1[ci - 1]);
                chars2.Add(seq2[cj - 1]);
                ci--; cj--;
                if (local && prev == 0) break; // Smith-Waterman: the local alignment starts here
                state = PickState(prev, mm[ci, cj], StateMatch, xx[ci, cj], StateGapInSeq2, StateGapInSeq1);
            }
            else if (state == StateGapInSeq2)
            {
                int cur = xx[ci, cj];
                chars1.Add(seq1[ci - 1]);
                chars2.Add('-');
                ci--;
                if (ci == 0 && cj == 0) break; // leading gap reached the origin
                state = PickState(cur, mm[ci, cj] + oe, StateMatch, xx[ci, cj] + e, StateGapInSeq2, StateGapInSeq1);
            }
            else // StateGapInSeq1
            {
                int cur = yy[ci, cj];
                chars1.Add('-');
                chars2.Add(seq2[cj - 1]);
                cj--;
                if (ci == 0 && cj == 0) break;
                state = PickState(cur, mm[ci, cj] + oe, StateMatch, yy[ci, cj] + e, StateGapInSeq1, StateGapInSeq2);
            }
        }

        if (!local)
            AppendFreeLeadingOverhang(chars1, chars2, seq1, seq2, ci, cj);

        chars1.Reverse();
        chars2.Reverse();

        return new AlignmentResult(
            AlignedSequence1: new string(chars1.ToArray()),
            AlignedSequence2: new string(chars2.ToArray()),
            Score: best,
            AlignmentType: ResultType(local, ends),
            StartPosition1: ci,
            StartPosition2: cj,
            EndPosition1: endI - 1,
            EndPosition2: endJ - 1);
    }

    /// <summary>
    /// One cell of Gotoh's three-state recurrences (M, X, Y at (i,j)) with the unreachable states clamped to
    /// <see cref="AffineNegInf"/>. With <paramref name="forbidden"/> the aligned pair (a_i, b_j) may not be used
    /// (Waterman-Eggert declumping): the substitution state M(i,j) is unreachable, the gap states are unchanged.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    internal static void AffineCell(
        int[,] mm, int[,] xx, int[,] yy, int i, int j, int sub, int startFloor, int e, int oe, bool forbidden)
    {
        int m = forbidden
            ? AffineNegInf
            : Math.Max(startFloor, Max3(mm[i - 1, j - 1], xx[i - 1, j - 1], yy[i - 1, j - 1])) + sub;
        int x = Max3(mm[i - 1, j] + oe, xx[i - 1, j] + e, yy[i - 1, j] + oe);
        int y = Max3(mm[i, j - 1] + oe, yy[i, j - 1] + e, xx[i, j - 1] + oe);

        // Clamp unreachable states so repeated penalties never drift towards overflow.
        xx[i, j] = x < AffineNegInf ? AffineNegInf : x;
        yy[i, j] = y < AffineNegInf ? AffineNegInf : y;
        mm[i, j] = m < AffineNegInf ? AffineNegInf : m;
    }

    /// <summary>Tag of a result: Local, otherwise Global when no end gap is free, else SemiGlobal.</summary>
    private static AlignmentType ResultType(bool local, EndGaps ends)
    {
        if (local) return AlignmentType.Local;
        return ends == EndGaps.None ? AlignmentType.Global : AlignmentType.SemiGlobal;
    }

    private static int Max3(int a, int b, int c) => Math.Max(a, Math.Max(b, c));

    /// <summary>Returns the first predecessor state whose candidate value produced <paramref name="target"/>
    /// (first, then second, otherwise the remaining state) — the fixed traceback tie order.</summary>
    private static byte PickState(int target, int firstValue, byte firstState, int secondValue, byte secondState, byte otherState)
    {
        if (firstValue == target) return firstState;
        if (secondValue == target) return secondState;
        return otherState;
    }

    #endregion

    #region Local Alignment (Smith-Waterman)

    /// <summary>
    /// Performs local alignment using the Smith-Waterman algorithm.
    /// Finds the best local alignment between subsequences.
    /// </summary>
    /// <remarks>
    /// Linear gap model (Smith &amp; Waterman 1981 with W_k = k·d): every gap position scores
    /// <see cref="ScoringMatrix.GapExtend"/>; <see cref="ScoringMatrix.GapOpen"/> is ignored. For affine
    /// gap costs (GapOpen + k·GapExtend) use <see cref="LocalAlignAffine(DnaSequence, DnaSequence, ScoringMatrix?)"/>.
    /// H(i,j) = max(0, H(i-1,j-1) + s(a_i,b_j), H(i-1,j) + d, H(i,j-1) + d), zero first row/column.
    /// The reported alignment ends at the first maximal cell in row-major order (smallest end in
    /// sequence 1, then in sequence 2; parasail instead prefers the smallest end in sequence 2) and the
    /// traceback (ties: diagonal &gt; up &gt; left) stops at the first zero cell; the result is always one
    /// of the co-optimal alignments Biopython's local PairwiseAligner enumerates. Start/End positions
    /// are 0-based inclusive indices of the aligned substrings; with no positive-scoring pair the result
    /// has Score 0, empty strings, <see cref="AlignmentType.Local"/> and all positions −1 — this is also the
    /// result for an empty side or two empty inputs (every <c>LocalAlign</c>/<c>LocalAlignAffine</c> overload, string
    /// or <see cref="DnaSequence"/>); null → <see cref="ArgumentNullException"/>.
    /// </remarks>
    /// <param name="sequence1">First DNA sequence.</param>
    /// <param name="sequence2">Second DNA sequence.</param>
    /// <param name="scoring">Scoring matrix (default: SimpleDna).</param>
    /// <returns>Alignment result with aligned subsequences and score.</returns>
    [System.Runtime.CompilerServices.OverloadResolutionPriority(1)]
    public static AlignmentResult LocalAlign(
        DnaSequence sequence1,
        DnaSequence sequence2,
        ScoringMatrix? scoring = null)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);

        return LocalAlignCore(sequence1.Sequence, sequence2.Sequence, scoring ?? SimpleDna);
    }

    /// <summary>
    /// Performs local alignment on raw sequence strings.
    /// </summary>
    [System.Runtime.CompilerServices.OverloadResolutionPriority(1)]
    public static AlignmentResult LocalAlign(
        string sequence1,
        string sequence2,
        ScoringMatrix? scoring = null)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);

        return LocalAlignCore(
            sequence1.ToUpperInvariant(),
            sequence2.ToUpperInvariant(),
            scoring ?? SimpleDna);
    }

    internal static AlignmentResult LocalAlignCore(string seq1, string seq2, PairScoring scoring, bool forceFullMatrix = false)
    {
        int m = seq1.Length;
        int n = seq2.Length;
        scoring.EnsureRange(m, n);

        if (!forceFullMatrix && ExceedsFullMatrixLimit(m, n))
            return LinearSpaceAlign(seq1, seq2, scoring, affine: false, local: true, EndGaps.None, LinearSpaceBaseCells, default);

        var score = new int[m + 1, n + 1];
        int maxScore = 0;
        int maxI = 0, maxJ = 0;

        // Fill the matrix (with zero floor for local alignment)
        for (int i = 1; i <= m; i++)
        {
            AlignmentCancellation.ThrowIfRequested();
            for (int j = 1; j <= n; j++)
            {
                score[i, j] = LocalLinearCell(score, i, j, scoring.Sub(seq1[i - 1], seq2[j - 1]), scoring.GapExtend, forbidden: false);

                if (score[i, j] > maxScore)
                {
                    maxScore = score[i, j];
                    maxI = i;
                    maxJ = j;
                }
            }
        }

        // Traceback from max score
        return TracebackLocal(seq1, seq2, score, maxI, maxJ, scoring);
    }

    /// <summary>
    /// One cell of the linear-gap local recurrence H(i,j) = max(0, H(i-1,j-1) + s, H(i-1,j) + d, H(i,j-1) + d).
    /// With <paramref name="forbidden"/> the aligned pair (a_i, b_j) may not be used (Waterman-Eggert declumping,
    /// <see cref="LocalAlignTopN(string, string, ScoringMatrix?, int, int)"/>): the diagonal move is dropped, the gap moves stay.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    internal static int LocalLinearCell(int[,] score, int i, int j, int sub, int gap, bool forbidden)
    {
        int best = Math.Max(score[i - 1, j] + gap, score[i, j - 1] + gap);
        if (!forbidden)
            best = Math.Max(best, score[i - 1, j - 1] + sub);
        return Math.Max(0, best);
    }

    private static AlignmentResult TracebackLocal(
        string seq1, string seq2, int[,] score, int endI, int endJ, PairScoring scoring,
        HashSet<long>? forbidden = null)
    {
        var chars1 = new List<char>();
        var chars2 = new List<char>();

        int i = endI, j = endJ;
        int startI = endI, startJ = endJ;

        while (i > 0 && j > 0 && score[i, j] > 0)
        {
            startI = i;
            startJ = j;

            int matchScore = scoring.Sub(seq1[i - 1], seq2[j - 1]);

            if ((forbidden is null || !forbidden.Contains(PackCell(i, j))) && score[i, j] == score[i - 1, j - 1] + matchScore)
            {
                chars1.Add(seq1[i - 1]);
                chars2.Add(seq2[j - 1]);
                i--; j--;
            }
            else if (score[i, j] == score[i - 1, j] + scoring.GapExtend)
            {
                chars1.Add(seq1[i - 1]);
                chars2.Add('-');
                i--;
            }
            else
            {
                chars1.Add('-');
                chars2.Add(seq2[j - 1]);
                j--;
            }
        }

        chars1.Reverse();
        chars2.Reverse();

        return new AlignmentResult(
            AlignedSequence1: new string(chars1.ToArray()),
            AlignedSequence2: new string(chars2.ToArray()),
            Score: score[endI, endJ],
            AlignmentType: AlignmentType.Local,
            StartPosition1: startI - 1,
            StartPosition2: startJ - 1,
            EndPosition1: endI - 1,
            EndPosition2: endJ - 1);
    }

    #endregion

    #region Local Alignment with affine gaps (Smith-Waterman-Gotoh)

    /// <summary>
    /// Performs optimal local alignment with <b>affine gap costs</b> (Smith-Waterman-Gotoh: the
    /// Smith &amp; Waterman 1981 zero floor on Gotoh's 1982 three-state recurrences).
    /// </summary>
    /// <param name="sequence1">First DNA sequence.</param>
    /// <param name="sequence2">Second DNA sequence.</param>
    /// <param name="scoring">Scoring matrix (default: <see cref="SimpleDna"/>).</param>
    /// <remarks>
    /// <para>
    /// Gap-cost convention identical to <see cref="GlobalAlignAffine(DnaSequence, DnaSequence, ScoringMatrix?)"/>:
    /// a gap of length k scores <c>GapOpen + k·GapExtend</c>. Equivalent settings: Biopython
    /// <c>PairwiseAligner(mode='local', open_gap_score = GapOpen + GapExtend, extend_gap_score = GapExtend)</c>;
    /// parasail <c>sw(open = -(GapOpen + GapExtend), extend = -GapExtend)</c>; BLAST+
    /// <c>-gapopen -GapOpen -gapextend -GapExtend</c>. With <c>GapOpen = 0</c> the score equals the
    /// linear-gap <see cref="LocalAlign(DnaSequence, DnaSequence, ScoringMatrix?)"/>. Gap scores must be
    /// non-positive (GapOpen ≤ 0, GapExtend ≤ 0 and greater than int.MinValue, as in every standard scheme and in <see cref="MatrixScoring"/>); a positive
    /// gap score is an <see cref="ArgumentException"/>.
    /// </para>
    /// <para>
    /// Recurrences (o = GapOpen, e = GapExtend):
    /// M(i,j) = max(0, M(i-1,j-1), X(i-1,j-1), Y(i-1,j-1)) + s(a_i,b_j);
    /// X(i,j) = max(M(i-1,j) + o + e, X(i-1,j) + e, Y(i-1,j) + o + e);
    /// Y(i,j) = max(M(i,j-1) + o + e, Y(i,j-1) + e, X(i,j-1) + o + e);
    /// all border cells −∞. Score = max M(i,j). The reported alignment ends at the first maximal cell
    /// in row-major order (smallest end in sequence 1, then in sequence 2; parasail instead prefers the
    /// smallest end in sequence 2) and the traceback stops as soon as the predecessor value is 0;
    /// traceback ties are broken M &gt; X &gt; Y. The result is always one of the co-optimal alignments
    /// Biopython's local PairwiseAligner enumerates.
    /// </para>
    /// <para>
    /// Coordinates are 0-based inclusive indices into the inputs. When no positive-scoring pair exists,
    /// the result has <c>Score = 0</c>, empty aligned strings and all coordinates −1 (same as
    /// <see cref="LocalAlign(DnaSequence, DnaSequence, ScoringMatrix?)"/>).
    /// </para>
    /// <para>Sources: Smith TF, Waterman MS (1981) J Mol Biol 147:195-197; Gotoh O (1982) J Mol Biol
    /// 162:705-708; cross-checked against Biopython 1.88 PairwiseAligner (local) and parasail 1.3.4 sw_trace.</para>
    /// </remarks>
    [System.Runtime.CompilerServices.OverloadResolutionPriority(1)]
    public static AlignmentResult LocalAlignAffine(
        DnaSequence sequence1,
        DnaSequence sequence2,
        ScoringMatrix? scoring = null)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);

        return AffineAlignCore(sequence1.Sequence, sequence2.Sequence, scoring ?? SimpleDna, AlignmentType.Local);
    }

    /// <summary>
    /// Affine-gap local alignment (Smith-Waterman-Gotoh) on raw sequence strings (uppercased internally).
    /// See <see cref="LocalAlignAffine(DnaSequence, DnaSequence, ScoringMatrix?)"/> for the model.
    /// </summary>
    [System.Runtime.CompilerServices.OverloadResolutionPriority(1)]
    public static AlignmentResult LocalAlignAffine(
        string sequence1,
        string sequence2,
        ScoringMatrix? scoring = null)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);

        return AffineAlignCore(
            sequence1.ToUpperInvariant(),
            sequence2.ToUpperInvariant(),
            scoring ?? SimpleDna,
            AlignmentType.Local);
    }

    /// <summary>Local result with no positive-scoring region: empty strings, score 0, coordinates −1
    /// (the zero-size traceback endpoints, as returned by the linear <c>LocalAlign</c>).</summary>
    private static AlignmentResult EmptyLocalResult() =>
        new("", "", 0, AlignmentType.Local, -1, -1, -1, -1);

    #endregion

    #region Semi-Global Alignment

    /// <summary>
    /// Performs semi-global alignment of the <b>fitting</b> (query-in-reference, "glocal") kind:
    /// <paramref name="sequence1"/> is aligned end-to-end against the best-scoring substring of
    /// <paramref name="sequence2"/>; unaligned leading/trailing residues of sequence 2 are free.
    /// </summary>
    /// <param name="sequence1">Query (aligned in full; typically the shorter sequence).</param>
    /// <param name="sequence2">Reference (its unaligned prefix/suffix cost nothing).</param>
    /// <param name="scoring">Scoring matrix (default: SimpleDna).</param>
    /// <returns>Alignment result.</returns>
    /// <remarks>
    /// <para>
    /// Only end gaps <i>in sequence 1</i> (i.e. overhanging sequence-2 residues) are free; end gaps in
    /// sequence 2 (overhanging query residues) are charged like any other gap. This is Rosalind's
    /// "fitting alignment" (SIMS), Biopython <c>PairwiseAligner(mode='global')</c> with
    /// <c>end_insertion_score = 0</c> (formerly <c>target_end_gap_score</c>; target = sequence 1), and
    /// parasail <c>sg_dx</c> (s1 = sequence 1, s2 = sequence 2). Linear gap model: every gap position
    /// scores <see cref="ScoringMatrix.GapExtend"/>; <see cref="ScoringMatrix.GapOpen"/> is ignored —
    /// use <see cref="SemiGlobalAlignAffine(DnaSequence, DnaSequence, ScoringMatrix?)"/> for affine costs.
    /// </para>
    /// <para>
    /// F(0,j) = 0, F(i,0) = d·i, F(i,j) = max(F(i-1,j-1) + s(a_i,b_j), F(i-1,j) + d, F(i,j-1) + d);
    /// Score = max_j F(m,j). Among tied end columns the smallest j wins (with all-negative scores this
    /// can be j = 0, i.e. the query aligned entirely against gaps); traceback ties are broken
    /// diagonal &gt; up &gt; left. The aligned strings span both full inputs (the free sequence-2
    /// prefix/suffix is shown against '-'). Coordinates: sequence 1 is aligned in full (Start1 = 0, End1 = length − 1);
    /// Start2/End2 are the <b>fitted interval</b> of sequence 2 (0-based inclusive: first and last reference index consumed
    /// by the alignment; an empty interval has End2 = Start2 − 1) — Biopython <c>aligned</c>/<c>coordinates</c> of the fitted part,
    /// parasail <c>end_ref</c>. The interval is that of the returned alignment; co-optimal alignments may fit other intervals.
    /// </para>
    /// </remarks>
    [System.Runtime.CompilerServices.OverloadResolutionPriority(1)]
    public static AlignmentResult SemiGlobalAlign(
        DnaSequence sequence1,
        DnaSequence sequence2,
        ScoringMatrix? scoring = null)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);

        return LinearAlignCore(sequence1.Sequence, sequence2.Sequence, scoring ?? SimpleDna, EndGaps.Fitting);
    }

    /// <summary>
    /// Linear-gap semi-global (fitting) alignment on raw sequence strings (uppercased internally; null →
    /// <see cref="ArgumentNullException"/>, an empty side is not special-cased). See
    /// <see cref="SemiGlobalAlign(DnaSequence, DnaSequence, ScoringMatrix?)"/> for the model and coordinates.
    /// </summary>
    [System.Runtime.CompilerServices.OverloadResolutionPriority(1)]
    public static AlignmentResult SemiGlobalAlign(
        string sequence1,
        string sequence2,
        ScoringMatrix? scoring = null)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);

        return LinearAlignCore(
            sequence1.ToUpperInvariant(),
            sequence2.ToUpperInvariant(),
            scoring ?? SimpleDna,
            EndGaps.Fitting);
    }

    /// <summary>
    /// Performs semi-global (fitting) alignment with <b>affine gap costs</b>: Gotoh's (1982) three-state
    /// recurrences with free leading/trailing sequence-2 residues.
    /// </summary>
    /// <param name="sequence1">Query (aligned in full).</param>
    /// <param name="sequence2">Reference (its unaligned prefix/suffix cost nothing).</param>
    /// <param name="scoring">Scoring matrix (default: <see cref="SimpleDna"/>).</param>
    /// <remarks>
    /// <para>
    /// Same fitting variant as <see cref="SemiGlobalAlign(DnaSequence, DnaSequence, ScoringMatrix?)"/> and
    /// the same gap-cost convention as <see cref="GlobalAlignAffine(DnaSequence, DnaSequence, ScoringMatrix?)"/>:
    /// an internal gap (or a gap in sequence 2 at an end) of length k scores <c>GapOpen + k·GapExtend</c>.
    /// Equivalent settings: Biopython <c>PairwiseAligner(mode='global', open_gap_score = GapOpen + GapExtend,
    /// extend_gap_score = GapExtend, end_insertion_score = 0)</c>; parasail
    /// <c>sg_dx(open = -(GapOpen + GapExtend), extend = -GapExtend)</c>. With <c>GapOpen = 0</c> the score
    /// equals the linear <see cref="SemiGlobalAlign(DnaSequence, DnaSequence, ScoringMatrix?)"/>.
    /// </para>
    /// <para>
    /// Initialization: M(0,j) = 0 for all j, X(i,0) = GapOpen + i·GapExtend, other border cells −∞;
    /// Score = max_j max(M, X, Y)(m, j) (first maximal column wins; end-state ties M &gt; X &gt; Y;
    /// traceback ties M &gt; X &gt; Y). The aligned strings span both full inputs; Start1/End1 are 0 and length − 1,
    /// Start2/End2 the fitted interval of sequence 2 exactly as for
    /// <see cref="SemiGlobalAlign(DnaSequence, DnaSequence, ScoringMatrix?)"/>.
    /// </para>
    /// <para>Sources: Gotoh O (1982) J Mol Biol 162:705-708; Flouri et al. (2015) bioRxiv 10.1101/031500;
    /// cross-checked against Biopython 1.88 PairwiseAligner and parasail 1.3.4 sg_dx_trace.</para>
    /// </remarks>
    [System.Runtime.CompilerServices.OverloadResolutionPriority(1)]
    public static AlignmentResult SemiGlobalAlignAffine(
        DnaSequence sequence1,
        DnaSequence sequence2,
        ScoringMatrix? scoring = null)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);

        return AffineAlignCore(sequence1.Sequence, sequence2.Sequence, scoring ?? SimpleDna, AlignmentType.SemiGlobal, EndGaps.Fitting);
    }

    /// <summary>
    /// Affine-gap semi-global (fitting) alignment on raw sequence strings (uppercased internally).
    /// See <see cref="SemiGlobalAlignAffine(DnaSequence, DnaSequence, ScoringMatrix?)"/> for the model.
    /// </summary>
    [System.Runtime.CompilerServices.OverloadResolutionPriority(1)]
    public static AlignmentResult SemiGlobalAlignAffine(
        string sequence1,
        string sequence2,
        ScoringMatrix? scoring = null)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);

        return AffineAlignCore(
            sequence1.ToUpperInvariant(),
            sequence2.ToUpperInvariant(),
            scoring ?? SimpleDna,
            AlignmentType.SemiGlobal,
            EndGaps.Fitting);
    }

    #endregion

    #region Ends-free (generalised semi-global) alignment

    // One option on the existing aligners: the end-gap policy (EndGaps) of the four alignment ends. The very same linear /
    // affine dynamic programs run (LinearAlignCore / AffineAlignCore); the policy only changes the border initialisation (a
    // free start has value 0 on the corresponding first row / column) and the end cell (maximum of the last row and/or
    // column), so EndGaps.Fitting is bit-for-bit the classic SemiGlobalAlign and EndGaps.None is GlobalAlign.
    // Oracles: Biopython 1.88 PairwiseAligner left/right_insertion/deletion_score = 0, parasail 1.3.4 sg_* (see Semi_Global_Alignment.md §5.4).

    /// <summary>
    /// Ends-free (generalised semi-global) alignment with <b>linear</b> gaps: each of the four end gaps (leading/trailing
    /// residues of sequence 1/2) is free or penalised according to <paramref name="endGaps"/>.
    /// </summary>
    /// <param name="sequence1">First sequence.</param>
    /// <param name="sequence2">Second sequence.</param>
    /// <param name="scoring">Scoring matrix (null = <see cref="SimpleDna"/>); only <see cref="ScoringMatrix.GapExtend"/> is used (d per gap position).</param>
    /// <param name="endGaps">Which end gaps are free: <see cref="EndGaps.Fitting"/> (= the two-argument overload), the overlap
    /// (suffix-prefix) policies <see cref="EndGaps.OverlapSuffix1Prefix2"/> / <see cref="EndGaps.OverlapSuffix2Prefix1"/>,
    /// <see cref="EndGaps.All"/> (EMBOSS needle's default), any combination of the four flags, or <see cref="EndGaps.None"/> (global).</param>
    /// <remarks>
    /// <para>
    /// Model: F(i,0) = 0 if <see cref="EndGaps.FreeStart1"/> else d·i; F(0,j) = 0 if <see cref="EndGaps.FreeStart2"/> else d·j;
    /// the recurrence is Needleman-Wunsch's; Score = F(m,n) when no end is free, otherwise the maximum over the last row
    /// (<see cref="EndGaps.FreeEnd2"/>: trailing sequence-2 residues free) and/or last column (<see cref="EndGaps.FreeEnd1"/>:
    /// trailing sequence-1 residues free). Ties: the first maximal last-row column (smallest j), then the first maximal last-column
    /// row (smallest i; row candidates first); traceback diagonal &gt; up &gt; left. The unaligned free flanks are shown against '-', so
    /// the aligned strings always span both complete inputs.
    /// </para>
    /// <para>
    /// Coordinates (0-based inclusive; the path = everything between the free flanks): Start1/Start2 = first residue of sequence
    /// 1/2 consumed by the alignment path (the number of free leading residues), End1/End2 = last consumed residue (an empty
    /// interval has End = Start − 1). <see cref="AlignmentResult.AlignmentType"/> is <see cref="AlignmentType.SemiGlobal"/>, or
    /// <see cref="AlignmentType.Global"/> for <see cref="EndGaps.None"/>.
    /// </para>
    /// <para>Sources: Biopython 1.88 <c>PairwiseAligner</c> (<c>left/right_insertion/deletion_score</c>), parasail 1.3.4 <c>sg_*</c>, EMBOSS needle
    /// <c>-endweight</c>; Durbin et al. 1998 §2.3 (overlap/fitting/free end-gap variants), Gusfield 1997 §11.6.4.</para>
    /// </remarks>
    [System.Runtime.CompilerServices.OverloadResolutionPriority(1)]
    public static AlignmentResult SemiGlobalAlign(
        string sequence1,
        string sequence2,
        ScoringMatrix? scoring,
        EndGaps endGaps)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);

        return LinearAlignCore(
            sequence1.ToUpperInvariant(),
            sequence2.ToUpperInvariant(),
            scoring ?? SimpleDna,
            endGaps);
    }

    /// <summary><see cref="SemiGlobalAlign(string, string, ScoringMatrix?, EndGaps)"/> for DNA sequences.</summary>
    [System.Runtime.CompilerServices.OverloadResolutionPriority(1)]
    public static AlignmentResult SemiGlobalAlign(
        DnaSequence sequence1,
        DnaSequence sequence2,
        ScoringMatrix? scoring,
        EndGaps endGaps)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);

        return LinearAlignCore(sequence1.Sequence, sequence2.Sequence, scoring ?? SimpleDna, endGaps);
    }

    /// <summary>
    /// Ends-free alignment with <b>affine</b> gaps (Gotoh's three states): a gap of length k scores
    /// <see cref="ScoringMatrix.GapOpen"/> + k·<see cref="ScoringMatrix.GapExtend"/> unless it is a free end gap (then 0,
    /// whatever its length). See <see cref="SemiGlobalAlign(string, string, ScoringMatrix?, EndGaps)"/> for the policy, ties and coordinates.
    /// </summary>
    /// <remarks>
    /// Initialisation: M(i,0) = 0 if <see cref="EndGaps.FreeStart1"/> else X(i,0) = o + i·e; M(0,j) = 0 if
    /// <see cref="EndGaps.FreeStart2"/> else Y(0,j) = o + j·e (Flouri et al. 2015 for the penalised starts); the end cell is the
    /// first maximum of max(M, X, Y) over the free last row/column. Equivalent Biopython settings: <c>open_gap_score = GapOpen +
    /// GapExtend</c>, <c>extend_gap_score = GapExtend</c> and <c>left/right_insertion/deletion_score = 0</c> for the free ends
    /// (<see cref="EndGaps"/>); parasail <c>sg_*</c> with <c>open = -(GapOpen + GapExtend)</c>, <c>extend = -GapExtend</c>.
    /// </remarks>
    [System.Runtime.CompilerServices.OverloadResolutionPriority(1)]
    public static AlignmentResult SemiGlobalAlignAffine(
        string sequence1,
        string sequence2,
        ScoringMatrix? scoring,
        EndGaps endGaps)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);

        return AffineAlignCore(
            sequence1.ToUpperInvariant(),
            sequence2.ToUpperInvariant(),
            scoring ?? SimpleDna,
            AlignmentType.SemiGlobal,
            endGaps);
    }

    /// <summary><see cref="SemiGlobalAlignAffine(string, string, ScoringMatrix?, EndGaps)"/> for DNA sequences.</summary>
    [System.Runtime.CompilerServices.OverloadResolutionPriority(1)]
    public static AlignmentResult SemiGlobalAlignAffine(
        DnaSequence sequence1,
        DnaSequence sequence2,
        ScoringMatrix? scoring,
        EndGaps endGaps)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);

        return AffineAlignCore(sequence1.Sequence, sequence2.Sequence, scoring ?? SimpleDna, AlignmentType.SemiGlobal, endGaps);
    }

    /// <summary>Ends-free alignment with a substitution matrix and linear gap scores (EMBOSS needle's default
    /// <c>-endweight N</c> = <see cref="EndGaps.All"/>); see <see cref="SemiGlobalAlign(string, string, ScoringMatrix?, EndGaps)"/>.</summary>
    public static AlignmentResult SemiGlobalAlign(string sequence1, string sequence2, MatrixScoring scoring, EndGaps endGaps) =>
        AlignWithMatrix(sequence1, sequence2, scoring, AlignmentType.SemiGlobal, affine: false, endGaps);

    /// <summary><see cref="SemiGlobalAlign(string, string, MatrixScoring, EndGaps)"/> for protein sequences.</summary>
    public static AlignmentResult SemiGlobalAlign(ProteinSequence sequence1, ProteinSequence sequence2, MatrixScoring scoring, EndGaps endGaps) =>
        AlignWithMatrix(sequence1, sequence2, scoring, AlignmentType.SemiGlobal, affine: false, endGaps);

    /// <summary><see cref="SemiGlobalAlign(string, string, MatrixScoring, EndGaps)"/> for DNA sequences.</summary>
    public static AlignmentResult SemiGlobalAlign(DnaSequence sequence1, DnaSequence sequence2, MatrixScoring scoring, EndGaps endGaps) =>
        AlignWithMatrix(sequence1, sequence2, scoring, AlignmentType.SemiGlobal, affine: false, endGaps);

    /// <summary>Ends-free alignment with a substitution matrix and affine gap scores (needle default
    /// <c>-endweight N</c> with <see cref="MatrixScoring.Emboss"/> numbers = <see cref="EndGaps.All"/>); see
    /// <see cref="SemiGlobalAlignAffine(string, string, ScoringMatrix?, EndGaps)"/>.</summary>
    public static AlignmentResult SemiGlobalAlignAffine(string sequence1, string sequence2, MatrixScoring scoring, EndGaps endGaps) =>
        AlignWithMatrix(sequence1, sequence2, scoring, AlignmentType.SemiGlobal, affine: true, endGaps);

    /// <summary><see cref="SemiGlobalAlignAffine(string, string, MatrixScoring, EndGaps)"/> for protein sequences.</summary>
    public static AlignmentResult SemiGlobalAlignAffine(ProteinSequence sequence1, ProteinSequence sequence2, MatrixScoring scoring, EndGaps endGaps) =>
        AlignWithMatrix(sequence1, sequence2, scoring, AlignmentType.SemiGlobal, affine: true, endGaps);

    /// <summary><see cref="SemiGlobalAlignAffine(string, string, MatrixScoring, EndGaps)"/> for DNA sequences.</summary>
    public static AlignmentResult SemiGlobalAlignAffine(DnaSequence sequence1, DnaSequence sequence2, MatrixScoring scoring, EndGaps endGaps) =>
        AlignWithMatrix(sequence1, sequence2, scoring, AlignmentType.SemiGlobal, affine: true, endGaps);

    #endregion

    #region Substitution-matrix alignment (BLOSUM62 / PAM250 / EDNAFULL, fractional gap scores)

    // The matrix overloads below add residue-pair scoring (protein, IUPAC DNA) and real-valued gap scores to the
    // six aligners without a second dynamic program: a MatrixScoring is converted to the exact integer model
    // "every score = k / scale" (PairScoring) and the very same LinearAlignCore / LocalAlignCore / AffineAlignCore
    // run on it; the score is divided by the scale afterwards. Tie order (diagonal > up > left linear, M > X > Y
    // affine), coordinates and empty-input contracts are therefore identical to the integer overloads.
    //
    // Oracles (see docs/algorithms/Alignment/Global_Alignment_Needleman_Wunsch.md §5.4 and docs/Evidence/ALIGN-GLOBAL-001-Evidence.md): EMBOSS needle/water 6.6.0
    // (EBLOSUM62/EPAM250/EDNAFULL, gapopen G gapextend E == MatrixScoring.Emboss(matrix, G, E)), Biopython 1.88
    // PairwiseAligner(substitution_matrix=...) and parasail 1.3.4.

    /// <summary>
    /// Global (Needleman-Wunsch) alignment scored with a substitution matrix and <b>linear</b> gap scores
    /// (<see cref="MatrixScoring.GapExtend"/> per gap position; <see cref="MatrixScoring.GapOpen"/> unused).
    /// </summary>
    /// <remarks>Residue symbols are compared case-insensitively against the matrix alphabet; a symbol outside it throws
    /// <see cref="ArgumentException"/>. <see cref="AlignmentResult.ExactScore"/> carries the real-valued score.
    /// Empty/null contract and tie order as for <see cref="GlobalAlign(string, string, ScoringMatrix?)"/>.</remarks>
    /// <exception cref="ArgumentNullException">A sequence or <paramref name="scoring"/> is null.</exception>
    /// <exception cref="ArgumentException">Unknown residue symbol, gap score &gt; 0/non-finite/finer than 1e-6, or a
    /// score range that would overflow the integer dynamic program.</exception>
    public static AlignmentResult GlobalAlign(string sequence1, string sequence2, MatrixScoring scoring) =>
        AlignWithMatrix(sequence1, sequence2, scoring, AlignmentType.Global, affine: false);

    /// <summary><see cref="GlobalAlign(string, string, MatrixScoring)"/> for protein sequences.</summary>
    public static AlignmentResult GlobalAlign(ProteinSequence sequence1, ProteinSequence sequence2, MatrixScoring scoring) =>
        AlignWithMatrix(sequence1, sequence2, scoring, AlignmentType.Global, affine: false);

    /// <summary><see cref="GlobalAlign(string, string, MatrixScoring)"/> for DNA sequences (e.g. <see cref="SubstitutionMatrix.EdnaFull"/>).</summary>
    public static AlignmentResult GlobalAlign(DnaSequence sequence1, DnaSequence sequence2, MatrixScoring scoring) =>
        AlignWithMatrix(sequence1, sequence2, scoring, AlignmentType.Global, affine: false);

    /// <summary>
    /// Global alignment with a substitution matrix and <b>affine</b> gap scores: a gap of length k scores
    /// <see cref="MatrixScoring.GapOpen"/> + k·<see cref="MatrixScoring.GapExtend"/> (EMBOSS needle with
    /// <see cref="MatrixScoring.Emboss"/> numbers: gapopen is the cost of a gap of length 1).
    /// </summary>
    /// <remarks>Same model as <see cref="GlobalAlignAffine(string, string, ScoringMatrix?)"/> (Gotoh, Flouri border).</remarks>
    public static AlignmentResult GlobalAlignAffine(string sequence1, string sequence2, MatrixScoring scoring) =>
        AlignWithMatrix(sequence1, sequence2, scoring, AlignmentType.Global, affine: true);

    /// <summary><see cref="GlobalAlignAffine(string, string, MatrixScoring)"/> for protein sequences.</summary>
    public static AlignmentResult GlobalAlignAffine(ProteinSequence sequence1, ProteinSequence sequence2, MatrixScoring scoring) =>
        AlignWithMatrix(sequence1, sequence2, scoring, AlignmentType.Global, affine: true);

    /// <summary><see cref="GlobalAlignAffine(string, string, MatrixScoring)"/> for DNA sequences.</summary>
    public static AlignmentResult GlobalAlignAffine(DnaSequence sequence1, DnaSequence sequence2, MatrixScoring scoring) =>
        AlignWithMatrix(sequence1, sequence2, scoring, AlignmentType.Global, affine: true);

    /// <summary>Local (Smith-Waterman) alignment with a substitution matrix and linear gap scores
    /// (see <see cref="GlobalAlign(string, string, MatrixScoring)"/>; EMBOSS water with gapopen = gapextend).</summary>
    public static AlignmentResult LocalAlign(string sequence1, string sequence2, MatrixScoring scoring) =>
        AlignWithMatrix(sequence1, sequence2, scoring, AlignmentType.Local, affine: false);

    /// <summary><see cref="LocalAlign(string, string, MatrixScoring)"/> for protein sequences.</summary>
    public static AlignmentResult LocalAlign(ProteinSequence sequence1, ProteinSequence sequence2, MatrixScoring scoring) =>
        AlignWithMatrix(sequence1, sequence2, scoring, AlignmentType.Local, affine: false);

    /// <summary><see cref="LocalAlign(string, string, MatrixScoring)"/> for DNA sequences.</summary>
    public static AlignmentResult LocalAlign(DnaSequence sequence1, DnaSequence sequence2, MatrixScoring scoring) =>
        AlignWithMatrix(sequence1, sequence2, scoring, AlignmentType.Local, affine: false);

    /// <summary>Local (Smith-Waterman-Gotoh) alignment with a substitution matrix and affine gap scores
    /// (EMBOSS water with <see cref="MatrixScoring.Emboss"/> numbers).</summary>
    public static AlignmentResult LocalAlignAffine(string sequence1, string sequence2, MatrixScoring scoring) =>
        AlignWithMatrix(sequence1, sequence2, scoring, AlignmentType.Local, affine: true);

    /// <summary><see cref="LocalAlignAffine(string, string, MatrixScoring)"/> for protein sequences.</summary>
    public static AlignmentResult LocalAlignAffine(ProteinSequence sequence1, ProteinSequence sequence2, MatrixScoring scoring) =>
        AlignWithMatrix(sequence1, sequence2, scoring, AlignmentType.Local, affine: true);

    /// <summary><see cref="LocalAlignAffine(string, string, MatrixScoring)"/> for DNA sequences.</summary>
    public static AlignmentResult LocalAlignAffine(DnaSequence sequence1, DnaSequence sequence2, MatrixScoring scoring) =>
        AlignWithMatrix(sequence1, sequence2, scoring, AlignmentType.Local, affine: true);

    /// <summary>Semi-global (fitting: sequence 1 inside sequence 2, free sequence-2 flanks) alignment with a substitution
    /// matrix and linear gap scores; see <see cref="SemiGlobalAlign(string, string, ScoringMatrix?)"/> for the variant and coordinates.</summary>
    public static AlignmentResult SemiGlobalAlign(string sequence1, string sequence2, MatrixScoring scoring) =>
        AlignWithMatrix(sequence1, sequence2, scoring, AlignmentType.SemiGlobal, affine: false);

    /// <summary><see cref="SemiGlobalAlign(string, string, MatrixScoring)"/> for protein sequences.</summary>
    public static AlignmentResult SemiGlobalAlign(ProteinSequence sequence1, ProteinSequence sequence2, MatrixScoring scoring) =>
        AlignWithMatrix(sequence1, sequence2, scoring, AlignmentType.SemiGlobal, affine: false);

    /// <summary><see cref="SemiGlobalAlign(string, string, MatrixScoring)"/> for DNA sequences.</summary>
    public static AlignmentResult SemiGlobalAlign(DnaSequence sequence1, DnaSequence sequence2, MatrixScoring scoring) =>
        AlignWithMatrix(sequence1, sequence2, scoring, AlignmentType.SemiGlobal, affine: false);

    /// <summary>Semi-global (fitting) alignment with a substitution matrix and affine gap scores; see
    /// <see cref="SemiGlobalAlignAffine(string, string, ScoringMatrix?)"/>.</summary>
    public static AlignmentResult SemiGlobalAlignAffine(string sequence1, string sequence2, MatrixScoring scoring) =>
        AlignWithMatrix(sequence1, sequence2, scoring, AlignmentType.SemiGlobal, affine: true);

    /// <summary><see cref="SemiGlobalAlignAffine(string, string, MatrixScoring)"/> for protein sequences.</summary>
    public static AlignmentResult SemiGlobalAlignAffine(ProteinSequence sequence1, ProteinSequence sequence2, MatrixScoring scoring) =>
        AlignWithMatrix(sequence1, sequence2, scoring, AlignmentType.SemiGlobal, affine: true);

    /// <summary><see cref="SemiGlobalAlignAffine(string, string, MatrixScoring)"/> for DNA sequences.</summary>
    public static AlignmentResult SemiGlobalAlignAffine(DnaSequence sequence1, DnaSequence sequence2, MatrixScoring scoring) =>
        AlignWithMatrix(sequence1, sequence2, scoring, AlignmentType.SemiGlobal, affine: true);

    private static AlignmentResult AlignWithMatrix(
        ProteinSequence sequence1, ProteinSequence sequence2, MatrixScoring scoring, AlignmentType type, bool affine,
        EndGaps ends = EndGaps.Fitting, long linearSpaceBaseCells = 0)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        return AlignWithMatrix(sequence1.Sequence, sequence2.Sequence, scoring, type, affine, ends, linearSpaceBaseCells);
    }

    private static AlignmentResult AlignWithMatrix(
        DnaSequence sequence1, DnaSequence sequence2, MatrixScoring scoring, AlignmentType type, bool affine,
        EndGaps ends = EndGaps.Fitting, long linearSpaceBaseCells = 0)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        return AlignWithMatrix(sequence1.Sequence, sequence2.Sequence, scoring, type, affine, ends, linearSpaceBaseCells);
    }

    /// <param name="sequence1">First sequence.</param>
    /// <param name="sequence2">Second sequence.</param>
    /// <param name="scoring">Substitution matrix and gap scores.</param>
    /// <param name="type">Global, local or semi-global.</param>
    /// <param name="affine">Affine (Gotoh) instead of linear gaps.</param>
    /// <param name="ends">End-gap policy, used only by <see cref="AlignmentType.SemiGlobal"/> (Global = none, Local = n/a).</param>
    /// <param name="linearSpaceBaseCells">0: the usual full-matrix core (automatic linear-space switch above
    /// <see cref="LinearSpaceCellThreshold"/>); &gt; 0: always the linear-space program with that base-case size.</param>
    private static AlignmentResult AlignWithMatrix(
        string sequence1, string sequence2, MatrixScoring scoring, AlignmentType type, bool affine,
        EndGaps ends = EndGaps.Fitting, long linearSpaceBaseCells = 0)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        ArgumentNullException.ThrowIfNull(scoring);

        string a = sequence1.ToUpperInvariant();
        string b = sequence2.ToUpperInvariant();
        var ps = PairScoring.Create(scoring, a, b);

        EndGaps affineEnds = type == AlignmentType.SemiGlobal ? ends : EndGaps.None;
        AlignmentResult r;
        if (linearSpaceBaseCells > 0)
            r = LinearSpaceAlign(a, b, ps, affine, type == AlignmentType.Local,
                type == AlignmentType.SemiGlobal ? ends : EndGaps.None, linearSpaceBaseCells, default);
        else
            r = affine
                ? AffineAlignCore(a, b, ps, type, affineEnds)
                : type switch
                {
                    AlignmentType.Global => GlobalAlignCore(a, b, ps),
                    AlignmentType.Local => LocalAlignCore(a, b, ps),
                    _ => LinearAlignCore(a, b, ps, ends),
                };

        double exact = (double)r.Score / ps.Scale;
        return r with { Score = (int)Math.Round(exact, MidpointRounding.AwayFromZero), ExactScore = exact };
    }

    #endregion

    #region Traceback

    private static AlignmentResult Traceback(
        string seq1, string seq2, int[,] score, int i, int j,
        PairScoring scoring, EndGaps ends,
        ForbiddenPairs? forbidden = null, int rowOffset = 0, int colOffset = 0)
    {
        // Preserve the traceback origin for correct score retrieval: score[endI, endJ] is the optimum
        // (global: the bottom-right cell; ends-free: the chosen last-row / last-column cell).
        int endI = i, endJ = j;
        bool freeStart1 = ends.HasFlag(EndGaps.FreeStart1);
        bool freeStart2 = ends.HasFlag(EndGaps.FreeStart2);

        var chars1 = new List<char>();
        var chars2 = new List<char>();

        // Free trailing overhang (appended first, will be at the end after reverse; both buffers are reversed in their
        // entirety below, so the unmatched suffix must be appended in REVERSE order here, otherwise it would itself be
        // reversed in the output). At most one of the two suffixes is non-empty: the end cell lies in the last row or column.
        for (int k = seq2.Length; k > j; k--)
        {
            chars1.Add('-');
            chars2.Add(seq2[k - 1]);
        }
        for (int k = seq1.Length; k > i; k--)
        {
            chars1.Add(seq1[k - 1]);
            chars2.Add('-');
        }

        // The traceback stops at the origin or, when the matching start is free, on reaching row 0 (free sequence-2 prefix
        // seq2[0 .. j-1]) / column 0 (free sequence-1 prefix seq1[0 .. i-1]); the consumed intervals start at (i, j).
        while (!(i == 0 && (j == 0 || freeStart2)) && !(j == 0 && freeStart1))
        {
            if (i > 0 && j > 0)
            {
                int matchScore = scoring.Sub(seq1[i - 1], seq2[j - 1]);

                if (score[i, j] == score[i - 1, j - 1] + matchScore
                    && (forbidden is null || !forbidden.Contains(rowOffset + i, colOffset + j)))
                {
                    chars1.Add(seq1[i - 1]);
                    chars2.Add(seq2[j - 1]);
                    i--; j--;
                    continue;
                }
            }

            if (i > 0 && (j == 0 || score[i, j] == score[i - 1, j] + scoring.GapExtend))
            {
                chars1.Add(seq1[i - 1]);
                chars2.Add('-');
                i--;
            }
            else
            {
                chars1.Add('-');
                chars2.Add(seq2[j - 1]);
                j--;
            }
        }

        AppendFreeLeadingOverhang(chars1, chars2, seq1, seq2, i, j);

        chars1.Reverse();
        chars2.Reverse();

        return new AlignmentResult(
            AlignedSequence1: new string(chars1.ToArray()),
            AlignedSequence2: new string(chars2.ToArray()),
            Score: score[endI, endJ],
            AlignmentType: ResultType(local: false, ends),
            StartPosition1: i,
            StartPosition2: j,
            EndPosition1: endI - 1,
            EndPosition2: endJ - 1);
    }

    /// <summary>Appends (in reverse order, see the callers) the free leading overhang the traceback stopped in front of:
    /// sequence-2 residues <c>seq2[0 .. j-1]</c> against gaps (stop in row 0) or sequence-1 residues <c>seq1[0 .. i-1]</c>
    /// against gaps (stop in column 0); nothing when the path reached the origin.</summary>
    private static void AppendFreeLeadingOverhang(
        List<char> chars1, List<char> chars2, string seq1, string seq2, int i, int j)
    {
        for (int k = j; k > 0 && i == 0; k--)
        {
            chars1.Add('-');
            chars2.Add(seq2[k - 1]);
        }
        for (int k = i; k > 0 && j == 0; k--)
        {
            chars1.Add(seq1[k - 1]);
            chars2.Add('-');
        }
    }

    #endregion

    #region Alignment Statistics

    // EMBOSS needle/water srspair markup, as written by EMBOSS 6.6 (ajalign.c alignConsStats + alignSim,
    // verified against the needle 6.6.0 binary: L/I, A/S, K/R in BLOSUM62 = ':' ; A/T = 0 -> '.' ;
    // gap -> ' '):
    //   '|' = identical column, ':' = similar column (positive substitution score, not identical),
    //   '.' = aligned but neither identical nor similar (mismatch), ' ' = column with a gap or an identical
    //   pair whose own score is not positive (X/X in BLOSUM62, N/N in EDNAFULL).
    private const char IdentityMark = '|';
    private const char SimilarityMark = ':';
    private const char MismatchMark = '.';
    private const char GapMark = ' ';
    private const char GapChar = '-';

    /// <summary>
    /// Calculates alignment statistics (identity, similarity, gaps, gap openings) from a pairwise
    /// alignment, following the EMBOSS needle/water report (Rice, Longden &amp; Bleasby 2000) and the BLAST
    /// identities / positives / gaps definitions.
    /// </summary>
    /// <param name="alignment">The pairwise alignment to summarise (rows of equal length, '-' = gap).</param>
    /// <param name="scoring">
    /// Scalar scoring model deciding which aligned pairs are <i>similar</i>: an identical pair is
    /// similar iff <see cref="ScoringMatrix.Match"/> &gt; 0 and a non-identical pair iff
    /// <see cref="ScoringMatrix.Mismatch"/> &gt; 0 (EMBOSS: "positive comparisons"). For the DNA presets
    /// (Match &gt; 0, Mismatch &lt; 0) Similarity equals Identity. Defaults to <see cref="SimpleDna"/>.
    /// For residue-pair specific scores (BLOSUM62, EDNAFULL, ...) use
    /// <see cref="CalculateStatistics(AlignmentResult, SubstitutionMatrix)"/>.
    /// </param>
    /// <remarks>
    /// <para>
    /// Conventions (EMBOSS needle documentation, <c>ajalign.c</c> and the needle 6.6.0 binary; Biopython 1.88
    /// <c>Alignment.counts</c> agrees on every count):
    /// Length = number of columns including gap columns, and it is the denominator of Identity, Similarity and
    /// Gaps (not the shorter sequence, not the number of non-gap columns);
    /// a column is a gap column when <i>either</i> row has '-' ("one or more sequences with a gap"), so Matches +
    /// Mismatches + Gaps = Length; residues are compared case-insensitively (EMBOSS and the aligners of this
    /// class upper-case their input); GapOpens = number of maximal '-' runs per row (BLAST tabular
    /// <c>gapopen</c>, Biopython <c>open_gaps</c>), independent of any gap scoring.
    /// </para>
    /// <para>
    /// A column that is a gap in both rows cannot occur in a pairwise alignment; if present (a pair cut out
    /// of an MSA) it counts as a gap column exactly as EMBOSS counts it.
    /// </para>
    /// </remarks>
    /// <exception cref="ArgumentNullException"><paramref name="alignment"/> is null.</exception>
    /// <exception cref="ArgumentException">The two aligned rows differ in length.</exception>
    public static AlignmentStatistics CalculateStatistics(
        AlignmentResult alignment,
        ScoringMatrix? scoring = null)
    {
        return CalculateStatisticsCore(alignment, ScalarClassifier(scoring ?? SimpleDna));
    }

    /// <summary>
    /// Alignment statistics with a residue-pair <see cref="SubstitutionMatrix"/> deciding similarity: a
    /// non-gap column is similar (a BLAST "positive") iff its substitution score is &gt; 0, whether or not the
    /// residues are identical (so an identical pair with a non-positive score, e.g. X/X in BLOSUM62, is not
    /// similar). This is the EMBOSS needle/water rule (ajalign.c alignConsStats: score &gt; 0 and at least one of
    /// the two residues has a positive self score — true for every column of the BLOSUM62, PAM250 and EDNAFULL
    /// tables whose pair score is positive) and agrees with Biopython <c>counts(matrix).positives</c>. A pair
    /// with a symbol outside the matrix alphabet is not similar. See
    /// <see cref="CalculateStatistics(AlignmentResult, ScoringMatrix?)"/> for the other conventions.
    /// </summary>
    public static AlignmentStatistics CalculateStatistics(
        AlignmentResult alignment,
        SubstitutionMatrix matrix)
    {
        ArgumentNullException.ThrowIfNull(matrix);
        return CalculateStatisticsCore(alignment, MatrixClassifier(matrix));
    }

    // Column classifiers. Similar = the column is counted in Positives / Similarity; Mark = srspair markup char.

    // Scalar model: identical pair is similar iff Match > 0, non-identical iff Mismatch > 0.
    private static Func<char, char, (bool Similar, char Mark)> ScalarClassifier(ScoringMatrix score) =>
        (a, b) =>
        {
            bool identical = a == b;
            bool similar = identical ? score.Match > 0 : score.Mismatch > 0;
            if (identical)
                return (similar, similar ? IdentityMark : GapMark);
            return (similar, similar ? SimilarityMark : MismatchMark);
        };

    // Residue-pair model, reproducing EMBOSS 6.6.0 alignConsStats for two sequences (ajalign.c; weights 1,
    // plurality 51 % => threshold 1.02). With self(x) = [M(x,x) > 0] and pair = [M(x,y) > 0]:
    //   Similarity count : pair AND (self(x) OR self(y))   (== BLAST "positive" for every standard matrix);
    //   markup           : h = self(x) + pair  ->  h = 2: '|' (x = y) or ':' ; h = 1: '.' ; h = 0: ' '.
    // So an identical pair whose own score is not positive (X/X in BLOSUM62, N/N in EDNAFULL) is an identity but
    // is marked ' ', and a pair whose FIRST residue has a non-positive self score (an IUPAC ambiguity code or X in
    // the first row) is marked '.' / ' ' even when M(x,y) > 0 (an EMBOSS display quirk, reproduced verbatim:
    // needle: R/A -> '.', A/R -> ':', yet both count as similar in the header). A symbol outside the matrix
    // alphabet contributes 0 (EMBOSS "Bad character", code 0).
    private static Func<char, char, (bool Similar, char Mark)> MatrixClassifier(SubstitutionMatrix matrix) =>
        (a, b) =>
        {
            bool self1 = matrix.TryGetScore(a, a, out int s1) && s1 > 0;
            bool self2 = matrix.TryGetScore(b, b, out int s2) && s2 > 0;
            bool pair = matrix.TryGetScore(a, b, out int sp) && sp > 0;
            int himatch = (self1 ? 1 : 0) + (pair ? 1 : 0);
            char mark = GapMark;
            if (himatch == 1)
                mark = MismatchMark;
            else if (himatch == 2)
                mark = a == b ? IdentityMark : SimilarityMark;
            return (pair && (self1 || self2), mark);
        };

    private static AlignmentStatistics CalculateStatisticsCore(
        AlignmentResult alignment,
        Func<char, char, (bool Similar, char Mark)> classify)
    {
        ArgumentNullException.ThrowIfNull(alignment);
        string row1 = alignment.AlignedSequence1;
        string row2 = alignment.AlignedSequence2;
        if (row1.Length != row2.Length)
            throw new ArgumentException(
                $"Aligned sequences must have equal length (got {row1.Length} and {row2.Length}).", nameof(alignment));

        int alignmentLength = row1.Length;
        if (alignmentLength == 0)
            return AlignmentStatistics.Empty;

        int matches = 0, mismatches = 0, gaps = 0, positives = 0, gapOpens = 0, blastPositives = 0;
        bool inGap1 = false, inGap2 = false;

        for (int i = 0; i < alignmentLength; i++)
        {
            char c1 = row1[i];
            char c2 = row2[i];
            bool gap1 = c1 == GapChar;
            bool gap2 = c2 == GapChar;

            if (gap1 && !inGap1) gapOpens++;
            if (gap2 && !inGap2) gapOpens++;
            inGap1 = gap1;
            inGap2 = gap2;

            if (gap1 || gap2)
            {
                gaps++;
                continue;
            }

            char u1 = char.ToUpperInvariant(c1);
            char u2 = char.ToUpperInvariant(c2);
            bool similar = classify(u1, u2).Similar;
            if (u1 == u2) matches++; else mismatches++;
            if (similar) positives++;
            if (u1 == u2 || similar) blastPositives++;
        }

        // Length (gap columns included) is the denominator of all three percentages (EMBOSS needle).
        return new AlignmentStatistics(
            Matches: matches,
            Mismatches: mismatches,
            Gaps: gaps,
            AlignmentLength: alignmentLength,
            Identity: (double)matches / alignmentLength * 100,
            Similarity: (double)positives / alignmentLength * 100,
            GapPercent: (double)gaps / alignmentLength * 100,
            Positives: positives,
            GapOpens: gapOpens,
            BlastPositives: blastPositives);
    }

    /// <summary>
    /// Generates a visual three-line pairwise alignment using the EMBOSS needle/water srspair
    /// markup line: '|' identity, ':' similarity (positive score), '.' mismatch, ' ' gap column (also an
    /// identical pair whose own substitution score is not positive, e.g. X/X in BLOSUM62 — EMBOSS 6.6.0 behaviour).
    /// For the <see cref="SubstitutionMatrix"/> overload the rule of EMBOSS alignConsStats is reproduced exactly,
    /// including its quirk for a first-row residue with a non-positive self score (see the classifier comment).
    /// </summary>
    /// <param name="alignment">The pairwise alignment to render.</param>
    /// <param name="lineWidth">Residues per block (must be positive; default 60).</param>
    /// <param name="scoring">
    /// Scalar scoring model used to mark similar (positive-score) columns; defaults to
    /// <see cref="SimpleDna"/>. See <see cref="CalculateStatistics(AlignmentResult, ScoringMatrix?)"/>.
    /// </param>
    public static string FormatAlignment(
        AlignmentResult alignment,
        int lineWidth = 60,
        ScoringMatrix? scoring = null)
    {
        return FormatAlignmentCore(alignment, lineWidth, ScalarClassifier(scoring ?? SimpleDna));
    }

    /// <summary>
    /// Three-line alignment with similarity marks (':') taken from a residue-pair
    /// <see cref="SubstitutionMatrix"/> (positive score), see
    /// <see cref="CalculateStatistics(AlignmentResult, SubstitutionMatrix)"/>.
    /// </summary>
    public static string FormatAlignment(
        AlignmentResult alignment,
        int lineWidth,
        SubstitutionMatrix matrix)
    {
        ArgumentNullException.ThrowIfNull(matrix);
        return FormatAlignmentCore(alignment, lineWidth, MatrixClassifier(matrix));
    }

    private static string FormatAlignmentCore(
        AlignmentResult alignment,
        int lineWidth,
        Func<char, char, (bool Similar, char Mark)> classify)
    {
        ArgumentNullException.ThrowIfNull(alignment);
        ArgumentOutOfRangeException.ThrowIfNegativeOrZero(lineWidth);
        if (alignment.AlignedSequence1.Length != alignment.AlignedSequence2.Length)
            throw new ArgumentException(
                $"Aligned sequences must have equal length (got {alignment.AlignedSequence1.Length} and {alignment.AlignedSequence2.Length}).",
                nameof(alignment));

        if (string.IsNullOrEmpty(alignment.AlignedSequence1))
            return "";

        var sb = new StringBuilder();
        int length = alignment.AlignedSequence1.Length;

        for (int start = 0; start < length; start += lineWidth)
        {
            int end = Math.Min(start + lineWidth, length);

            // Emit '\n' explicitly (not AppendLine/Environment.NewLine) so the formatted alignment
            // is byte-identical across platforms — the EMBOSS-style block must not depend on OS CRLF.
            sb.Append(alignment.AlignedSequence1[start..end]).Append('\n');

            // Markup line (EMBOSS srspair).
            for (int i = start; i < end; i++)
            {
                char c1 = alignment.AlignedSequence1[i];
                char c2 = alignment.AlignedSequence2[i];

                if (c1 == GapChar || c2 == GapChar)
                {
                    sb.Append(GapMark);
                    continue;
                }

                sb.Append(classify(char.ToUpperInvariant(c1), char.ToUpperInvariant(c2)).Mark);
            }
            sb.Append('\n');

            sb.Append(alignment.AlignedSequence2[start..end]).Append('\n');
            sb.Append('\n');
        }

        return sb.ToString();
    }

    #endregion

    #region Multiple Sequence Alignment (Center Star)

    // Center-star multiple sequence alignment (Gusfield 1993, Bull Math Biol 55:141-154, DOI
    // 10.1007/BF02460299; Gusfield 1997 "Algorithms on Strings, Trees and Sequences" §14.6.1):
    //   1. All C(k,2) pairwise OPTIMAL global alignments give a pair score D(i,j).
    //   2. The center S_c is the sequence with the best total pairwise score, argmax_c Σ_{j≠c} D(c,j)
    //      (Gusfield: the string minimising Σ_j distance(S_c, S_j)).
    //   3. The optimal pairwise alignments of every other sequence against S_c are merged into one
    //      MSA by "once a gap, always a gap": the number of residues of S_j inserted between two
    //      consecutive center residues is made equal for all rows (the maximum over rows), and
    //      residues are never moved between columns.
    // The pairwise stage is the canonical linear-gap Needleman-Wunsch (GlobalAlignCore), so for two
    // sequences the MSA is exactly an optimal pairwise alignment. An approximate, faster variant that
    // replaces step 1-2 by a 4-mer-cosine center and the pairwise alignments by suffix-tree anchored
    // alignments (AnchorBasedAligner, MUMmer-style) is MultipleAlignAnchored.

    /// <summary>
    /// Performs <b>center-star</b> multiple sequence alignment (Gusfield 1993) with exact pairwise
    /// alignments: the center is the sequence with the highest total optimal pairwise
    /// Needleman-Wunsch score to all others, every other sequence is aligned to it with the canonical
    /// <see cref="GlobalAlign(DnaSequence, DnaSequence, ScoringMatrix?)"/> dynamic program, and the
    /// pairwise alignments are merged column-wise (once a gap, always a gap).
    /// </summary>
    /// <remarks>
    /// <para>
    /// Linear gap model (default): every gap position scores <see cref="ScoringMatrix.GapExtend"/>
    /// (<see cref="ScoringMatrix.GapOpen"/> is not used, as in <c>GlobalAlign</c>); affine gaps are the additive option
    /// <see cref="MultipleAlign(IEnumerable{DnaSequence}, ScoringMatrix?, MsaOptions?)"/>. For two sequences
    /// the result is an optimal pairwise global alignment. Ties for the center are broken by the
    /// lowest input index.
    /// </para>
    /// <para>
    /// Complexity O(k²·L²) time for the pairwise scores (parallelised over pairs), O(L) memory per
    /// concurrent pair for scoring and bounded (≤ 16 MB, exact Hirschberg above that) for alignment. For large, closely related inputs the approximate
    /// <see cref="MultipleAlignAnchored(IEnumerable{DnaSequence}, ScoringMatrix?)"/> is faster.
    /// </para>
    /// <para>
    /// <see cref="MultipleAlignmentResult.TotalScore"/> is the column-based sum-of-pairs score of the
    /// final MSA (match/mismatch from the matrix, residue-vs-gap = <c>GapExtend</c>, gap-vs-gap = 0).
    /// </para>
    /// </remarks>
    /// <param name="sequences">Collection of sequences to align.</param>
    /// <param name="scoring">Scoring matrix (default: SimpleDna).</param>
    /// <returns>Multiple alignment result.</returns>
    /// <exception cref="ArgumentNullException">When <paramref name="sequences"/> is null.</exception>
    public static MultipleAlignmentResult MultipleAlign(
        IEnumerable<DnaSequence> sequences,
        ScoringMatrix? scoring = null) =>
        MultipleAlign(sequences, scoring, new MsaOptions());

    /// <summary>
    /// <see cref="MultipleAlign(IEnumerable{DnaSequence}, ScoringMatrix?)"/> with additive <see cref="MsaOptions"/>.
    /// With <see cref="MsaGapModel.Affine"/> the pairwise stage is the exact affine (Gotoh) global aligner
    /// (<see cref="GlobalAlignAffineLinearSpace(string, string, ScoringMatrix?, long)"/>, identical to <c>GlobalAlignAffine</c>): the center maximises
    /// the total optimal affine pair score, every sequence is aligned to it with a gap of length k costing
    /// <see cref="ScoringMatrix.GapOpen"/> + k·<see cref="ScoringMatrix.GapExtend"/>, the alignments are merged by once-a-gap-always-a-gap and
    /// <see cref="MultipleAlignmentResult.TotalScore"/> is the natural affine sum of pairs (Altschul 1989, see <see cref="SumOfPairsScore"/>).
    /// For two sequences the result is an optimal pairwise affine alignment. For k ≥ 3 the star method is a heuristic for either gap model.
    /// </summary>
    /// <param name="sequences">Collection of sequences to align.</param>
    /// <param name="scoring">Scoring matrix (default: SimpleDna).</param>
    /// <param name="options">Options; <c>null</c> = defaults (linear gaps).</param>
    [System.Runtime.CompilerServices.OverloadResolutionPriority(1)]
    public static MultipleAlignmentResult MultipleAlign(
        IEnumerable<DnaSequence> sequences,
        ScoringMatrix? scoring,
        MsaOptions? options)
    {
        ArgumentNullException.ThrowIfNull(sequences);
        return StarAlignCore(DnaRows(sequences), scoring ?? SimpleDna, options);
    }

    /// <summary>
    /// <see cref="MultipleAlign(IEnumerable{DnaSequence}, ScoringMatrix?, MsaOptions?)"/> for <b>protein</b> sequences (or any alphabet of the
    /// <see cref="SubstitutionMatrix"/>): the same exact center-star (Gusfield 1993) — optimal pairwise alignments under the matrix with linear
    /// (<see cref="MatrixScoring.GapExtend"/> per position) or, with <see cref="MsaGapModel.Affine"/>, affine
    /// (<see cref="MatrixScoring.GapOpen"/> + k·<see cref="MatrixScoring.GapExtend"/>) gap scores, possibly fractional (EMBOSS-style 0.5);
    /// the matrix aligners share the dynamic programs, merge and objective of the DNA overload (no second implementation).
    /// </summary>
    /// <remarks>
    /// Residues are compared case-insensitively against the matrix alphabet (rows are returned upper-case); a residue outside it throws
    /// <see cref="ArgumentException"/>. <see cref="MultipleAlignmentResult.ExactTotalScore"/> is the real-valued sum of pairs
    /// (<see cref="MatrixSumOfPairsScore(IReadOnlyList{string}, MatrixScoring, MsaGapModel)"/>); the consensus is the column majority (ties: smaller
    /// residue letter, a residue is preferred to a gap). For two sequences the result is an optimal pairwise alignment (<c>GlobalAlign</c> /
    /// <c>GlobalAlignAffine</c> with the same <see cref="MatrixScoring"/>).
    /// </remarks>
    /// <param name="sequences">Protein sequences.</param>
    /// <param name="scoring">Substitution matrix (BLOSUM62, PAM250, ...) and gap scores.</param>
    /// <param name="options">Options; <c>null</c> = defaults (linear gaps).</param>
    public static MultipleAlignmentResult MultipleAlign(
        IEnumerable<ProteinSequence> sequences,
        MatrixScoring scoring,
        MsaOptions? options = null)
    {
        var (rows, ps) = PrepareMatrixMsa(sequences, scoring);
        return ScaleMsaResult(StarAlignCore(rows, ps, options), ps);
    }

    /// <summary><see cref="MultipleAlign(IEnumerable{ProteinSequence}, MatrixScoring, MsaOptions?)"/> for raw sequence strings (protein, or DNA with
    /// <see cref="SubstitutionMatrix.EdnaFull"/>; upper-cased internally).</summary>
    public static MultipleAlignmentResult MultipleAlign(
        IEnumerable<string> sequences,
        MatrixScoring scoring,
        MsaOptions? options = null)
    {
        var (rows, ps) = PrepareMatrixMsa(sequences, scoring);
        return ScaleMsaResult(StarAlignCore(rows, ps, options), ps);
    }

    private static MultipleAlignmentResult StarAlignCore(string[] rawSeqs, PairScoring effectiveScoring, MsaOptions? options)
    {
        using var cancelScope = AlignmentCancellation.Enter(options?.CancellationToken ?? default);
        var gapModel = ResolveGapModel(options);

        var trivial = TrivialMultipleAlignment(rawSeqs);
        if (trivial != null)
            return trivial;

        bool affine = gapModel == MsaGapModel.Affine;
        int dop = options?.MaxDegreeOfParallelism ?? -1;
        int centerIdx = SelectCenterByPairwiseScore(rawSeqs, affine
            ? (a, b) => LinearSpaceAlign(a, b, effectiveScoring, true, false, EndGaps.None, StarFullMatrixCells, AlignmentCancellation.Token).Score
            : (a, b) => LinearScoreOnly(a, b, effectiveScoring), dop);
        string centerStr = rawSeqs[centerIdx];

        var pairwise = AlignAllToCenter(rawSeqs.Length, centerIdx,
            i => LinearSpaceAlign(centerStr, rawSeqs[i], effectiveScoring, affine, false, EndGaps.None, StarFullMatrixCells, AlignmentCancellation.Token), dop);
        return FinishStarAlignment(MergeCenterAlignments(centerStr, pairwise, centerIdx), effectiveScoring, gapModel);
    }

    private static string[] DnaRows(IEnumerable<DnaSequence> sequences) => DnaStringsNoNull(sequences, nameof(sequences)).ToArray();

    /// <summary>Upper-cased rows and the exact integer model of a <see cref="MatrixScoring"/> for a whole set of sequences (validates residues, gaps, range).</summary>
    private static (string[] Rows, PairScoring Scoring) PrepareMatrixMsa(IEnumerable<string> sequences, MatrixScoring scoring)
    {
        ArgumentNullException.ThrowIfNull(sequences);
        ArgumentNullException.ThrowIfNull(scoring);
        var list = new List<string>();
        foreach (var s in sequences)
        {
            if (s is null)
                throw new ArgumentNullException(nameof(sequences), "A sequence is null.");
            list.Add(s.ToUpperInvariant());
        }

        var rows = list.ToArray();
        return (rows, PairScoring.Create(scoring, rows));
    }

    private static (string[] Rows, PairScoring Scoring) PrepareMatrixMsa(IEnumerable<ProteinSequence> sequences, MatrixScoring scoring)
    {
        ArgumentNullException.ThrowIfNull(sequences);
        var list = new List<string>();
        foreach (var p in sequences)
        {
            if (p is null)
                throw new ArgumentNullException(nameof(sequences), "A sequence is null.");
            list.Add(p.Sequence);
        }

        return PrepareMatrixMsa(list, scoring);
    }

    /// <summary>Converts the scaled integer sum of pairs of a matrix run into the real score (<see cref="MultipleAlignmentResult.ExactTotalScore"/>).</summary>
    private static MultipleAlignmentResult ScaleMsaResult(MultipleAlignmentResult r, PairScoring ps)
    {
        if (ps.Scale == 1)
            return r;
        double exact = (double)r.TotalScore / ps.Scale;
        return r with { TotalScore = (int)Math.Round(exact, MidpointRounding.AwayFromZero), ExactTotalScore = exact };
    }

    /// <summary>
    /// Approximate, faster center-star alignment: the center is chosen by 4-mer count cosine
    /// similarity (no pairwise DP) and every other sequence is aligned to it with the
    /// suffix-tree anchored aligner (<see cref="AnchorBasedAligner"/>: exact-match anchors, gaps
    /// between anchors filled by Needleman-Wunsch). Unlike
    /// <see cref="MultipleAlign(IEnumerable{DnaSequence}, ScoringMatrix?)"/> the pairwise alignments
    /// are <b>not guaranteed optimal</b> (an anchor chain constrains the path), so the resulting
    /// sum-of-pairs score can be lower.
    /// </summary>
    /// <remarks>
    /// For closely related sequences where anchors cover a fraction α of the length the DP work
    /// shrinks by about (1−α)². Same output contract and linear gap model as
    /// <see cref="MultipleAlign(IEnumerable{DnaSequence}, ScoringMatrix?)"/>.
    /// </remarks>
    /// <param name="sequences">Collection of sequences to align.</param>
    /// <param name="scoring">Scoring matrix (default: SimpleDna).</param>
    /// <returns>Multiple alignment result.</returns>
    /// <exception cref="ArgumentNullException">When <paramref name="sequences"/> is null.</exception>
    public static MultipleAlignmentResult MultipleAlignAnchored(
        IEnumerable<DnaSequence> sequences,
        ScoringMatrix? scoring = null) =>
        MultipleAlignAnchored(sequences, scoring, new MsaOptions());

    /// <summary>
    /// <see cref="MultipleAlignAnchored(IEnumerable{DnaSequence}, ScoringMatrix?)"/> with additive <see cref="MsaOptions"/>
    /// (<see cref="MsaGapModel.Affine"/>: the inter-anchor segments are aligned with the exact affine aligner, the reported
    /// <see cref="MultipleAlignmentResult.TotalScore"/> is the natural affine sum of pairs; still approximate because the anchors constrain the path).
    /// </summary>
    /// <param name="sequences">Collection of sequences to align.</param>
    /// <param name="scoring">Scoring matrix (default: SimpleDna).</param>
    /// <param name="options">Options; <c>null</c> = defaults (linear gaps).</param>
    [System.Runtime.CompilerServices.OverloadResolutionPriority(1)]
    public static MultipleAlignmentResult MultipleAlignAnchored(
        IEnumerable<DnaSequence> sequences,
        ScoringMatrix? scoring,
        MsaOptions? options)
    {
        ArgumentNullException.ThrowIfNull(sequences);
        return AnchoredStarAlignCore(DnaRows(sequences), scoring ?? SimpleDna, options);
    }

    /// <summary>
    /// <see cref="MultipleAlignAnchored(IEnumerable{DnaSequence}, ScoringMatrix?, MsaOptions?)"/> for <b>protein</b> sequences (or any alphabet of
    /// the <see cref="SubstitutionMatrix"/>): the exact-match anchors (maximal exact matches of at least 8 residues, any alphabet) are scored
    /// as Σ matrix(c, c), the inter-anchor segments are aligned with the matrix aligner and the center is chosen by the 4-mer cosine
    /// similarity of the residue strings; still approximate, because the anchor chain constrains the path.
    /// </summary>
    public static MultipleAlignmentResult MultipleAlignAnchored(
        IEnumerable<ProteinSequence> sequences,
        MatrixScoring scoring,
        MsaOptions? options = null)
    {
        var (rows, ps) = PrepareMatrixMsa(sequences, scoring);
        return ScaleMsaResult(AnchoredStarAlignCore(rows, ps, options), ps);
    }

    /// <summary><see cref="MultipleAlignAnchored(IEnumerable{ProteinSequence}, MatrixScoring, MsaOptions?)"/> for raw sequence strings.</summary>
    public static MultipleAlignmentResult MultipleAlignAnchored(
        IEnumerable<string> sequences,
        MatrixScoring scoring,
        MsaOptions? options = null)
    {
        var (rows, ps) = PrepareMatrixMsa(sequences, scoring);
        return ScaleMsaResult(AnchoredStarAlignCore(rows, ps, options), ps);
    }

    private static MultipleAlignmentResult AnchoredStarAlignCore(string[] rawSeqs, PairScoring effectiveScoring, MsaOptions? options)
    {
        using var cancelScope = AlignmentCancellation.Enter(options?.CancellationToken ?? default);
        var gapModel = ResolveGapModel(options);

        var trivial = TrivialMultipleAlignment(rawSeqs);
        if (trivial != null)
            return trivial;

        int centerIdx = SelectCenterSequence(rawSeqs);
        string centerStr = rawSeqs[centerIdx];
        var centerTree = SuffixTree.SuffixTree.Build(centerStr);

        var pairwise = AlignAllToCenter(rawSeqs.Length, centerIdx,
            i => AnchorBasedAligner.AlignWithAnchors(centerStr, centerTree, rawSeqs[i], effectiveScoring,
                affineGaps: gapModel == MsaGapModel.Affine), options?.MaxDegreeOfParallelism ?? -1);
        var merged = MergeCenterAlignments(centerStr, pairwise, centerIdx);

        // An empty input yields an empty anchored result: pad to the common length.
        int maxLen = merged.Max(s => s.Length);
        for (int i = 0; i < merged.Count; i++)
            if (merged[i].Length < maxLen)
                merged[i] = merged[i].PadRight(maxLen, '-');

        return FinishStarAlignment(merged, effectiveScoring, gapModel);
    }

    // Pairs above this many DP cells are aligned by the exact linear-space recursion (bounded memory
    // per concurrent pair: 4M cells = 16 MB); smaller pairs use the canonical full-matrix NW.
    private const long StarFullMatrixCells = 4_000_000;

    /// <summary>Empty / single-sequence results shared by every MSA entry point; null when k ≥ 2.</summary>
    private static MultipleAlignmentResult? TrivialMultipleAlignment(string[] seqList)
    {
        if (seqList.Length == 0)
            return MultipleAlignmentResult.Empty;

        if (seqList.Length == 1)
        {
            return new MultipleAlignmentResult(
                AlignedSequences: new[] { seqList[0] },
                Consensus: seqList[0],
                TotalScore: 0);
        }

        return null;
    }

    private static MultipleAlignmentResult FinishStarAlignment(List<string> merged, PairScoring scoring, MsaGapModel gapModel = MsaGapModel.Linear)
    {
        int maxLen = merged.Max(s => s.Length);
        return new MultipleAlignmentResult(
            AlignedSequences: merged.ToArray(),
            Consensus: BuildConsensus(merged, maxLen, scoring.IsMatrix),
            TotalScore: ComputeSumOfPairsScore(merged, scoring, gapModel));
    }

    /// <summary>
    /// Runs <paramref name="alignToCenter"/> for every non-center sequence (in parallel for k ≥ 4;
    /// each slot is independent and written once, so the result does not depend on scheduling).
    /// </summary>
    private static AlignmentResult[] AlignAllToCenter(int k, int centerIdx, Func<int, AlignmentResult> alignToCenter, int maxDegreeOfParallelism = -1)
    {
        var results = new AlignmentResult[k];

        if (k >= 4)
        {
            ParallelFor(0, k, i =>
            {
                if (i != centerIdx)
                    results[i] = alignToCenter(i);
            }, maxDegreeOfParallelism);
        }
        else
        {
            for (int i = 0; i < k; i++)
                if (i != centerIdx)
                    results[i] = alignToCenter(i);
        }

        return results;
    }

    /// <summary>
    /// Gusfield (1993) center: the sequence maximising the total optimal pair score to the other sequences. The pair score is
    /// <paramref name="optimalPairScore"/>: the last cell of the linear-space Needleman-Wunsch pass (identical to the canonical
    /// <see cref="GlobalAlignCore"/> score) or the exact affine (Gotoh) optimum. The C(k, 2) pair dynamic programs run in parallel
    /// (<paramref name="maxDegreeOfParallelism"/>; an exception of a pair is rethrown with its original type). Ties go to the lowest index.
    /// </summary>
    private static int SelectCenterByPairwiseScore(string[] seqs, Func<string, string, int> optimalPairScore, int maxDegreeOfParallelism = -1)
    {
        int k = seqs.Length;
        var pairScore = new int[k, k];

        var pairs = new List<(int I, int J)>(k * (k - 1) / 2);
        for (int i = 0; i < k; i++)
            for (int j = i + 1; j < k; j++)
                pairs.Add((i, j));

        ParallelFor(0, pairs.Count, p =>
        {
            // Optimal pair score (linear NW: last row of the NW matrix in O(n) memory, equals GlobalAlignCore(..).Score).
            var pair = pairs[p];
            int s = optimalPairScore(seqs[pair.I], seqs[pair.J]);
            pairScore[pair.I, pair.J] = s;
            pairScore[pair.J, pair.I] = s;
        }, maxDegreeOfParallelism);

        long bestTotal = long.MinValue;
        int bestIdx = 0;
        for (int c = 0; c < k; c++)
        {
            long total = 0;
            for (int j = 0; j < k; j++)
                if (j != c)
                    total += pairScore[c, j];

            if (total > bestTotal)
            {
                bestTotal = total;
                bestIdx = c;
            }
        }

        return bestIdx;
    }

    /// <summary>
    /// Selects the center for the approximate anchored star alignment: the sequence with the highest
    /// total 4-mer count cosine similarity to all others. K-mers are counted by the canonical Core
    /// counter <see cref="SequenceExtensions.CountKmersSpan"/> (the rows are upper-case residue strings of any alphabet). Always prefers a non-empty center (an empty center has no
    /// positions to carry the anchored reconciliation).
    /// </summary>
    private static int SelectCenterSequence(string[] sequences)
    {
        int k = sequences.Length;

        int firstNonEmpty = -1;
        for (int i = 0; i < k; i++)
        {
            if (sequences[i].Length > 0)
            {
                firstNonEmpty = i;
                break;
            }
        }

        // All sequences empty (or none at all): the trivial empty-coordinate MSA is correct.
        if (firstNonEmpty < 0) return 0;

        if (k <= 2) return firstNonEmpty;

        const int KmerLen = 4;
        var profiles = new Dictionary<string, int>[k];
        var norms = new double[k];

        for (int i = 0; i < k; i++)
        {
            profiles[i] = sequences[i].AsSpan().CountKmersSpan(KmerLen);
            double norm = 0;
            foreach (var count in profiles[i].Values)
                norm += (double)count * count;
            norms[i] = Math.Sqrt(norm);
        }

        var totalSimilarity = new double[k];

        for (int i = 0; i < k; i++)
        {
            for (int j = i + 1; j < k; j++)
            {
                // Iterate the smaller profile; integer products are exact in double.
                var (small, large) = profiles[i].Count <= profiles[j].Count
                    ? (profiles[i], profiles[j])
                    : (profiles[j], profiles[i]);

                double dot = 0;
                foreach (var (kmer, count) in small)
                    if (large.TryGetValue(kmer, out int other))
                        dot += (double)count * other;

                double sim = (norms[i] > 0 && norms[j] > 0)
                    ? dot / (norms[i] * norms[j])
                    : 0;

                totalSimilarity[i] += sim;
                totalSimilarity[j] += sim;
            }
        }

        int bestIdx = firstNonEmpty;
        for (int i = firstNonEmpty + 1; i < k; i++)
        {
            if (sequences[i].Length == 0)
                continue;
            if (totalSimilarity[i] > totalSimilarity[bestIdx])
                bestIdx = i;
        }

        return bestIdx;
    }

    /// <summary>
    /// Merges the pairwise center-vs-sequence alignments into one MSA ("once a gap, always a gap",
    /// Feng &amp; Doolittle 1987 / Gusfield 1993): for every center position the number of
    /// non-center residues inserted before it is the maximum over all pairwise alignments, and each
    /// row's own inserted residues are left-justified in that slot.
    /// </summary>
    private static List<string> MergeCenterAlignments(string centerStr, AlignmentResult[] pairwise, int centerIdx)
    {
        int k = pairwise.Length;

        // Center aligns to itself (identity).
        pairwise[centerIdx] = new AlignmentResult(
            centerStr, centerStr, 0, AlignmentType.Global, 0, 0,
            centerStr.Length - 1, centerStr.Length - 1);

        // gapsBefore[p] = number of columns to insert before center position p (p = Length: trailing).
        var gapsBefore = new int[centerStr.Length + 1];

        for (int i = 0; i < k; i++)
        {
            if (i == centerIdx) continue;

            string alignedCenter = pairwise[i].AlignedSequence1;
            int origPos = 0;
            int gapCount = 0;

            for (int j = 0; j < alignedCenter.Length; j++)
            {
                if (alignedCenter[j] == '-')
                {
                    gapCount++;
                }
                else
                {
                    gapsBefore[origPos] = Math.Max(gapsBefore[origPos], gapCount);
                    origPos++;
                    gapCount = 0;
                }
            }
            gapsBefore[centerStr.Length] = Math.Max(gapsBefore[centerStr.Length], gapCount);
        }

        // Re-project each pairwise alignment into the merged coordinate space.
        var merged = new List<string>(k);

        for (int i = 0; i < k; i++)
        {
            string alignedCenter = pairwise[i].AlignedSequence1;
            string alignedOther = pairwise[i].AlignedSequence2;

            var sb = new StringBuilder();
            int alignPos = 0;

            for (int p = 0; p <= centerStr.Length; p++)
            {
                int neededGaps = gapsBefore[p];

                // Residues of this row inserted before center position p.
                int thisGaps = 0;
                int tempAlignPos = alignPos;
                while (tempAlignPos < alignedCenter.Length && alignedCenter[tempAlignPos] == '-')
                {
                    thisGaps++;
                    tempAlignPos++;
                }

                for (int g = 0; g < thisGaps; g++)
                {
                    sb.Append(alignPos < alignedOther.Length ? alignedOther[alignPos] : '-');
                    alignPos++;
                }

                // Pad the slot to the maximum insertion width.
                for (int g = thisGaps; g < neededGaps; g++)
                    sb.Append('-');

                // The row's character at center position p.
                if (p < centerStr.Length && alignPos < alignedOther.Length)
                {
                    sb.Append(alignedOther[alignPos]);
                    alignPos++;
                }
            }

            merged.Add(sb.ToString());
        }

        return merged;
    }


    /// <summary>
    /// Builds a majority-vote consensus from aligned sequences.
    /// Per Wikipedia MSA: consensus derived from aligned columns using majority voting.
    /// Gap characters participate in the vote; on tie, residues are preferred to the gap and, among residues,
    /// the smaller letter (A &lt; C &lt; G &lt; T for DNA). DNA rows count only A, C, G, T and the gap; with
    /// <paramref name="anySymbol"/> (substitution-matrix runs: protein, IUPAC DNA) every character is a residue.
    /// </summary>
    private static string BuildConsensus(List<string> aligned, int length, bool anySymbol = false)
    {
        var consensus = new StringBuilder(length);
        var counts = new Dictionary<char, int>();

        for (int pos = 0; pos < length; pos++)
        {
            counts.Clear();
            if (!anySymbol)
            {
                counts['A'] = 0; counts['C'] = 0; counts['G'] = 0; counts['T'] = 0; counts['-'] = 0;
            }

            foreach (var seq in aligned)
            {
                if (pos >= seq.Length)
                    continue;
                char c = seq[pos];
                if (anySymbol)
                    counts[c] = counts.GetValueOrDefault(c) + 1;
                else if (counts.TryGetValue(c, out int value))
                    counts[c] = value + 1;
            }

            // Include all characters (including gaps) in majority vote.
            // On tie between gap and residue, prefer the residue; between residues, the smaller letter.
            if (counts.Count == 0)
            {
                consensus.Append('-');
                continue;
            }

            char mostCommon = counts.OrderByDescending(kv => kv.Value)
                                   .ThenBy(kv => kv.Key == '-' ? 1 : 0)
                                   .ThenBy(kv => kv.Key)
                                   .First().Key;

            consensus.Append(mostCommon);
        }

        return consensus.ToString();
    }

    /// <summary>
    /// Computes the column-based sum-of-pairs (SP) score for a multiple sequence alignment.
    /// Per Wikipedia MSA: "sum of all of the pairs of characters at each position in the alignment."
    /// Scores all C(k,2) sequence pairs using column-based scoring:
    /// match/mismatch from scoring matrix, gap-nucleotide penalty per position, gap-gap = 0.
    /// </summary>
    private static int ComputeSumOfPairsScore(IReadOnlyList<string> aligned, PairScoring scoring, MsaGapModel gapModel)
    {
        long total = SumOfPairsScoreLong(aligned, scoring, gapModel);
        if (total > int.MaxValue || total < int.MinValue)
            throw new ArgumentException(
                $"The sum-of-pairs score {total} of {aligned.Count} rows exceeds the 32-bit integer range of the reported score; use fewer / shorter sequences or smaller scores.",
                "scoring");
        return (int)total;
    }

    /// <summary>The sum-of-pairs score accumulated in 64 bits (the exact value; the reported <c>int</c> scores are range-checked from it).</summary>
    private static long SumOfPairsScoreLong(IReadOnlyList<string> aligned, PairScoring scoring, MsaGapModel gapModel) =>
        gapModel == MsaGapModel.Affine ? AffineSumOfPairsScoreCore(aligned, scoring) : LinearSumOfPairsScoreCore(aligned, scoring);

    private static long LinearSumOfPairsScoreCore(IReadOnlyList<string> aligned, PairScoring scoring)
    {
        if (aligned.Count < 2) return 0;

        int length = aligned[0].Length;
        long totalScore = 0;

        for (int pos = 0; pos < length; pos++)
        {
            for (int i = 0; i < aligned.Count; i++)
            {
                char ci = pos < aligned[i].Length ? aligned[i][pos] : '-';

                for (int j = i + 1; j < aligned.Count; j++)
                {
                    char cj = pos < aligned[j].Length ? aligned[j][pos] : '-';

                    // Gap-gap is neutral (standard SP convention) and simply contributes nothing.
                    if (ci != '-' || cj != '-')
                    {
                        if (ci == '-' || cj == '-')
                            totalScore += scoring.GapExtend; // Gap-residue penalty
                        else
                            totalScore += scoring.Sub(ci, cj);
                    }
                }
            }
        }

        return totalScore;
    }

    #endregion

    #region Multiple Sequence Alignment (Guide-Tree Progressive / Feng-Doolittle)

    // Progressive (guide-tree) multiple sequence alignment, the Feng-Doolittle /
    // Clustal-style method, added as a SECOND aligner alongside the star MSA above.
    // The star MSA (MultipleAlign) is unchanged.
    //
    // Algorithm (Feng & Doolittle 1987; Wikipedia "Multiple sequence alignment"):
    //   1. Compute all C(k,2) pairwise global (Needleman-Wunsch) alignments and derive a
    //      distance matrix. Distance = 1 - identities / min(len1, len2): the ClustalW pairwise
    //      distance (Thompson, Higgins & Gibson 1994; pairalign.c tmat = (100 - mm_score) / 100,
    //      mm_score = 100 * count / MIN(len1, len2)), see PairwiseIdentityDistance. (Feng-Doolittle
    //      converts pairwise alignment scores to distances; identity-based distance is the
    //      standard conversion used by Clustal-style tools.)
    //   2. Build a guide tree over that distance matrix — "an efficient clustering
    //      method such as neighbor-joining or unweighted pair group method with arithmetic
    //      mean (UPGMA)" (Wikipedia "Multiple sequence alignment"): UPGMA (the classic
    //      Feng-Doolittle choice, default) or neighbour-joining (ClustalW's default), both built by the
    //      canonical PhylogeneticAnalyzer builders (MsaOptions.GuideTree; SequenceAligner.GuideTree.cs).
    //   3. Progressively align along the guide tree from the tips: sequence-sequence, then
    //      sequence-profile, then profile-profile, using the Needleman-Wunsch recurrence over
    //      *columns* (a profile is the set of already-aligned rows). Column-vs-column score is
    //      the average match/mismatch over all cross-profile residue pairs (sum-of-pairs profile
    //      scoring). The "once a gap, always a gap" rule (Feng & Doolittle 1987: "Once a gap is
    //      introduced ... it is preserved within all subsequent fusions") is enforced: gaps
    //      already present in a profile are never removed; merging two profiles only inserts
    //      whole new gap columns, never edits existing columns.
    //
    // Sources retrieved this session:
    //   - Feng & Doolittle (1987) "Progressive sequence alignment as a prerequisite to correct
    //     phylogenetic trees", J Mol Evol 25:351-360 (PubMed 3118049) — progressive method,
    //     guide tree, "once a gap, always a gap".
    //   - Wikipedia "Multiple sequence alignment", § Progressive alignment construction
    //     (https://en.wikipedia.org/wiki/Multiple_sequence_alignment) — guide tree built by
    //     NJ/UPGMA; sequences added sequentially along the guide tree.
    //   - Wikipedia UPGMA / Sokal & Michener (1958) — cluster averaging formula.

    /// <summary>
    /// Performs guide-tree <b>progressive</b> multiple sequence alignment (the Feng-Doolittle /
    /// Clustal-style method), as an additive alternative to the star-alignment
    /// <see cref="MultipleAlign(IEnumerable{DnaSequence}, ScoringMatrix?)"/>.
    /// </summary>
    /// <remarks>
    /// <para>
    /// Pipeline (Feng &amp; Doolittle 1987; Wikipedia "Multiple sequence alignment"):
    /// </para>
    /// <list type="number">
    /// <item>All pairwise Needleman-Wunsch alignments give a distance matrix using the ClustalW distance
    /// d = 1 − (identical aligned residue pairs) / (length of the shorter sequence)
    /// (<c>pairalign.c</c>: <c>tmat = (100 − 100·count / MIN(len1, len2)) / 100</c>).</item>
    /// <item>A guide tree is built over that matrix: UPGMA (the classic Feng-Doolittle choice, default) or ClustalW's
    /// neighbour-joining tree (<see cref="MsaOptions.GuideTree"/>).</item>
    /// <item>Profiles are aligned progressively from the tips with the NW recurrence over
    /// columns and sum-of-pairs profile scoring, enforcing "once a gap, always a gap".</item>
    /// </list>
    /// <para>
    /// Linear gaps by default; <see cref="MsaGapModel.Affine"/> (<see cref="MsaOptions.GapModel"/>) switches to a Gotoh three-state
    /// profile DP with occupancy-scaled <see cref="ScoringMatrix.GapOpen"/>/<see cref="ScoringMatrix.GapExtend"/> and the natural affine
    /// sum-of-pairs <c>TotalScore</c> (see <see cref="SumOfPairsScore"/>); exact for two sequences, heuristic for more.
    /// </para>
    /// </remarks>
    /// <param name="sequences">Collection of sequences to align.</param>
    /// <param name="scoring">Scoring matrix (default: <see cref="SimpleDna"/>).</param>
    /// <returns>Multiple alignment result (aligned rows, majority consensus, sum-of-pairs score).</returns>
    /// <exception cref="ArgumentNullException">When <paramref name="sequences"/> is null.</exception>
    public static MultipleAlignmentResult MultipleAlignProgressive(
        IEnumerable<DnaSequence> sequences,
        ScoringMatrix? scoring = null) =>
        MultipleAlignProgressive(sequences, scoring, new MsaOptions());

    /// <summary>
    /// <see cref="MultipleAlignProgressive(IEnumerable{DnaSequence}, ScoringMatrix?)"/> with additive
    /// <see cref="MsaOptions"/> (guide-tree distance source: the repository's linear NW alignment or ClustalW's own
    /// pairwise stage, <see cref="MsaDistanceMode"/>).
    /// </summary>
    /// <param name="sequences">Collection of sequences to align.</param>
    /// <param name="scoring">Scoring matrix (default: <see cref="SimpleDna"/>).</param>
    /// <param name="options">Options; <c>null</c> = defaults.</param>
    [System.Runtime.CompilerServices.OverloadResolutionPriority(1)]
    public static MultipleAlignmentResult MultipleAlignProgressive(
        IEnumerable<DnaSequence> sequences,
        ScoringMatrix? scoring,
        MsaOptions? options)
    {
        ArgumentNullException.ThrowIfNull(sequences);
        return ProgressiveAlignCore(DnaRows(sequences), scoring ?? SimpleDna, options);
    }

    /// <summary>
    /// <see cref="MultipleAlignProgressive(IEnumerable{DnaSequence}, ScoringMatrix?, MsaOptions?)"/> for <b>protein</b> sequences (or any alphabet of the
    /// <see cref="SubstitutionMatrix"/>): the same Feng-Doolittle pipeline — identity distances of the optimal pairwise alignments under the
    /// matrix, the configured guide tree, and profile-profile alignment with the column score Σ matrix(a, b) averaged over the cross pairs
    /// (linear gaps) or the occupancy-scaled Gotoh program (<see cref="MsaGapModel.Affine"/>), real-valued (possibly fractional) gap scores.
    /// See <see cref="MultipleAlign(IEnumerable{ProteinSequence}, MatrixScoring, MsaOptions?)"/> for the input contract.
    /// </summary>
    public static MultipleAlignmentResult MultipleAlignProgressive(
        IEnumerable<ProteinSequence> sequences,
        MatrixScoring scoring,
        MsaOptions? options = null)
    {
        var (rows, ps) = PrepareMatrixMsa(sequences, scoring);
        return ScaleMsaResult(ProgressiveAlignCore(rows, ps, options), ps);
    }

    /// <summary><see cref="MultipleAlignProgressive(IEnumerable{ProteinSequence}, MatrixScoring, MsaOptions?)"/> for raw sequence strings.</summary>
    public static MultipleAlignmentResult MultipleAlignProgressive(
        IEnumerable<string> sequences,
        MatrixScoring scoring,
        MsaOptions? options = null)
    {
        var (rows, ps) = PrepareMatrixMsa(sequences, scoring);
        return ScaleMsaResult(ProgressiveAlignCore(rows, ps, options), ps);
    }

    private static MultipleAlignmentResult ProgressiveAlignCore(string[] rawSeqs, PairScoring effectiveScoring, MsaOptions? options)
    {
        using var cancelScope = AlignmentCancellation.Enter(options?.CancellationToken ?? default);
        options ??= new MsaOptions();
        var gapModel = ResolveGapModel(options);

        var trivial = TrivialMultipleAlignment(rawSeqs);
        if (trivial != null)
            return trivial;

        int k = rawSeqs.Length;

        // Step 1: pairwise alignments -> distance matrix (ClustalW: 1 - identities / min(len1, len2)).
        var distance = ComputeGuideDistances(rawSeqs, effectiveScoring, options.Distance, gapModel, options.GuideTree, options.QuickPairwise, options.ClustalWProtein);

        // Step 2: guide tree (UPGMA by default, neighbour-joining on request).
        ProgressiveGuideNode root = BuildProgressiveGuideTree(k, distance, options.GuideTree);

        // Step 3: progressive alignment along the tree (post-order: tips first).
        Profile aligned = AlignProfileSubtree(root, rawSeqs, effectiveScoring, gapModel);

        // Reproject rows back into input order (the tree visited leaves in tree order).
        var orderedRows = new string[k];
        for (int r = 0; r < aligned.Rows.Count; r++)
            orderedRows[aligned.SequenceIndices[r]] = aligned.Rows[r];

        var alignedList = orderedRows.ToList();
        int maxLen = alignedList.Count == 0 ? 0 : alignedList.Max(s => s.Length);
        string consensus = BuildConsensus(alignedList, maxLen, effectiveScoring.IsMatrix);
        int spScore = ComputeSumOfPairsScore(alignedList, effectiveScoring, gapModel);

        return new MultipleAlignmentResult(
            AlignedSequences: orderedRows,
            Consensus: consensus,
            TotalScore: spScore);
    }

    /// <summary>
    /// Pairwise distance for the guide tree, the ClustalW pairwise distance (Thompson, Higgins &amp;
    /// Gibson 1994): 1 − (identical aligned residue pairs of the optimal pairwise alignment) /
    /// (length of the shorter sequence). This is <c>tmat = (100 − mm_score) / 100</c> with
    /// <c>mm_score = 100·count / MIN(len1, len2)</c> in <c>pairalign.c</c> of ClustalW 1.83 (the
    /// reference alignment there uses affine gaps; this method uses the repository's linear
    /// Needleman-Wunsch). Two empty sequences (no comparable residue) are treated as distance 0, an
    /// empty versus a non-empty sequence as distance 1.
    /// </summary>
    private static double PairwiseIdentityDistance(string a, string b, PairScoring scoring, MsaGapModel gapModel = MsaGapModel.Linear)
    {
        if (a.Length == 0 && b.Length == 0)
            return 0.0;
        if (a.Length == 0 || b.Length == 0)
            return 1.0; // no shared residues possible

        var pa = gapModel == MsaGapModel.Affine ? AffineAlignCore(a, b, scoring, AlignmentType.Global) : GlobalAlignCore(a, b, scoring);

        int identical = 0;
        for (int i = 0; i < pa.AlignedSequence1.Length; i++)
        {
            char c1 = pa.AlignedSequence1[i];
            if (c1 != '-' && c1 == pa.AlignedSequence2[i])
                identical++;
        }

        return 1.0 - (double)identical / Math.Min(a.Length, b.Length);
    }

    /// <summary>
    /// Symmetric guide-tree distance matrix of the progressive / iterative / consistency aligners:
    /// <see cref="MsaDistanceMode.Linear"/> = the repository's own global alignment under <paramref name="scoring"/> (linear NW, or
    /// the affine Gotoh alignment when <paramref name="gapModel"/> is <see cref="MsaGapModel.Affine"/>; <see cref="PairwiseIdentityDistance"/>); <see cref="MsaDistanceMode.ClustalWAffine"/> = ClustalW's own pairwise stage
    /// (<see cref="ClustalWPairwiseAlign(string, string, double, double)"/>, independent of <paramref name="scoring"/>);
    /// <see cref="MsaDistanceMode.ClustalWQuick"/> = ClustalW's fast k-tuple stage (<see cref="ClustalWQuickPairwiseAlign(string, string, ClustalWQuickOptions?, int?)"/>).
    /// </summary>
    private static double[,] ComputeGuideDistances(string[] rawSeqs, PairScoring scoring, MsaDistanceMode mode, MsaGapModel gapModel, MsaGuideTree guideTree = MsaGuideTree.Upgma, ClustalWQuickOptions? quickOptions = null, ClustalWProteinPairwiseOptions? proteinOptions = null)
    {
        if (!Enum.IsDefined(mode))
            throw new ArgumentOutOfRangeException(nameof(mode), mode, "Unknown guide-tree distance mode.");

        // Protein matrices select ClustalW's PROTEIN pairwise stages (GONNET / BLOSUM30 / PAM350 / identity, gap 10 / 0.1, k-tuple 1), nucleotide matrices
        // and the scalar scoring the DNA stages.
        bool protein = scoring.Protein;
        var proteinParameters = proteinOptions ?? new ClustalWProteinPairwiseOptions();
        var quick = protein ? proteinParameters.Quick : quickOptions ?? new ClustalWQuickOptions();
        int longest = 0;
        foreach (var s in rawSeqs)
            longest = Math.Max(longest, s.Length);

        int k = rawSeqs.Length;
        var distance = new double[k, k];
        for (int i = 0; i < k; i++)
        {
            for (int j = i + 1; j < k; j++)
            {
                double d;
                if (mode == MsaDistanceMode.ClustalWQuick)
                    d = protein
                        ? ClustalWQuickPairwiseAlign(rawSeqs[i], rawSeqs[j], ClustalWSequenceType.Protein, quick, longest).Distance
                        : ClustalWQuickPairwiseDistance(rawSeqs[i], rawSeqs[j], quick, longest);
                else if (mode != MsaDistanceMode.ClustalWAffine)
                    d = PairwiseIdentityDistance(rawSeqs[i], rawSeqs[j], scoring, gapModel);
                else if (protein)
                    d = guideTree == MsaGuideTree.NeighborJoining
                        ? ClustalWProteinSinglePrecisionDistance(rawSeqs[i], rawSeqs[j], proteinParameters)
                        : ClustalWPairwiseAlign(rawSeqs[i], rawSeqs[j], proteinParameters).Distance;
                else if (guideTree == MsaGuideTree.NeighborJoining)
                    d = ClustalWSinglePrecisionDistance(rawSeqs[i], rawSeqs[j]);
                else
                    d = ClustalWPairwiseDistance(rawSeqs[i], rawSeqs[j]);
                distance[i, j] = d;
                distance[j, i] = d;
            }
        }

        return distance;
    }

    /// <summary>
    /// A profile: a block of aligned rows of equal length, each tagged with the input-sequence
    /// index it represents (so the final output can be reordered back to input order).
    /// </summary>
    private sealed class Profile
    {
        public List<string> Rows { get; } = new();
        public List<int> SequenceIndices { get; } = new();
        public int Width => Rows.Count == 0 ? 0 : Rows[0].Length;
    }

    /// <summary>
    /// Post-order traversal of the guide tree: align each subtree's profile, combining child
    /// profiles by profile-profile NW alignment. Leaves start as a single-row profile.
    /// </summary>
    private static Profile AlignProfileSubtree(
        ProgressiveGuideNode node, string[] rawSeqs, PairScoring scoring, MsaGapModel gapModel)
    {
        if (node.IsLeaf)
        {
            var leaf = new Profile();
            leaf.Rows.Add(rawSeqs[node.LeafIndex]);
            leaf.SequenceIndices.Add(node.LeafIndex);
            return leaf;
        }

        var left = AlignProfileSubtree(node.Left!, rawSeqs, scoring, gapModel);
        var right = AlignProfileSubtree(node.Right!, rawSeqs, scoring, gapModel);
        return gapModel == MsaGapModel.Affine
            ? AlignProfilesAffine(left, right, scoring)
            : AlignProfiles(left, right, scoring);
    }

    /// <summary>
    /// Largest profile-vs-profile dynamic-programming matrix (<c>(m + 1)·(n + 1)</c> doubles, 8 bytes per cell = 400 MB) of the
    /// progressive / iterative / consistency multiple aligners. The exact center-star <see cref="MultipleAlign(IEnumerable{DnaSequence}, ScoringMatrix?, MsaOptions?)"/> needs no such
    /// matrix (its pairwise stage is linear-space); a profile-profile DP has no published linear-space variant that keeps the
    /// column-pair scores of these methods, so longer inputs get a clear error instead of an <see cref="OutOfMemoryException"/>.
    /// </summary>
    public const long ProfileMatrixMaxCells = 50_000_000;

    private static void RequireProfileMatrixFits(int m, int n)
    {
        long cells = (long)(m + 1) * (n + 1);
        if (cells > ProfileMatrixMaxCells)
            throw new ArgumentException(
                $"The profile-profile dynamic program needs O(m·n) memory: {m} x {n} alignment columns = {cells:N0} cells exceeds the limit of {ProfileMatrixMaxCells:N0}. " +
                "Use MultipleAlign (exact center-star, linear-space pairwise stage) for long sequences.");
    }

    /// <summary>
    /// Aligns two profiles with the Needleman-Wunsch recurrence over columns. The column-vs-column
    /// substitution score is the average match/mismatch over all cross-profile residue pairs
    /// (sum-of-pairs profile scoring); a gap column scores GapExtend per residue averaged.
    /// Existing columns are never edited — gaps are inserted as whole new all-gap columns into one
    /// profile, enforcing "once a gap, always a gap" (Feng &amp; Doolittle 1987).
    /// </summary>
    private static Profile AlignProfiles(Profile p1, Profile p2, PairScoring scoring)
    {
        int m = p1.Width;
        int n = p2.Width;
        RequireProfileMatrixFits(m, n);

        // F(i,j): best score aligning first i columns of p1 with first j columns of p2.
        var f = new double[m + 1, n + 1];
        for (int i = 1; i <= m; i++)
            f[i, 0] = f[i - 1, 0] + ProfileColumnVsGap(p1, i - 1, p2.Rows.Count, scoring);
        for (int j = 1; j <= n; j++)
            f[0, j] = f[0, j - 1] + ProfileColumnVsGap(p2, j - 1, p1.Rows.Count, scoring);

        for (int i = 1; i <= m; i++)
        {
            AlignmentCancellation.ThrowIfRequested();
            for (int j = 1; j <= n; j++)
            {
                double diag = f[i - 1, j - 1] + ProfileColumnScore(p1, i - 1, p2, j - 1, scoring);
                double up = f[i - 1, j] + ProfileColumnVsGap(p1, i - 1, p2.Rows.Count, scoring);
                double left = f[i, j - 1] + ProfileColumnVsGap(p2, j - 1, p1.Rows.Count, scoring);
                f[i, j] = Math.Max(diag, Math.Max(up, left));
            }
        }

        // Traceback, building column lists for both profiles.
        var cols1 = new List<int>(); // column index in p1, or -1 for a gap column
        var cols2 = new List<int>();
        int ii = m, jj = n;
        while (ii > 0 || jj > 0)
        {
            if (ii > 0 && jj > 0)
            {
                double diag = f[ii - 1, jj - 1] + ProfileColumnScore(p1, ii - 1, p2, jj - 1, scoring);
                if (AreClose(f[ii, jj], diag))
                {
                    cols1.Add(ii - 1);
                    cols2.Add(jj - 1);
                    ii--; jj--;
                    continue;
                }
            }

            if (ii > 0)
            {
                double up = f[ii - 1, jj] + ProfileColumnVsGap(p1, ii - 1, p2.Rows.Count, scoring);
                if (jj == 0 || AreClose(f[ii, jj], up))
                {
                    cols1.Add(ii - 1);
                    cols2.Add(-1);
                    ii--;
                    continue;
                }
            }

            cols1.Add(-1);
            cols2.Add(jj - 1);
            jj--;
        }

        cols1.Reverse();
        cols2.Reverse();

        var merged = new Profile();
        merged.SequenceIndices.AddRange(p1.SequenceIndices);
        merged.SequenceIndices.AddRange(p2.SequenceIndices);
        AppendReprojectedRows(merged, p1, cols1);
        AppendReprojectedRows(merged, p2, cols2);
        return merged;
    }

    /// <summary>
    /// Re-emits each row of <paramref name="source"/> according to <paramref name="columnMap"/>:
    /// a non-negative entry copies that existing column (gaps in it preserved), a -1 inserts a
    /// fresh gap. Existing columns are copied verbatim ("once a gap, always a gap").
    /// </summary>
    private static void AppendReprojectedRows(Profile dest, Profile source, List<int> columnMap)
    {
        foreach (var row in source.Rows)
        {
            var sb = new StringBuilder(columnMap.Count);
            foreach (int col in columnMap)
                sb.Append(col < 0 ? '-' : row[col]);
            dest.Rows.Add(sb.ToString());
        }
    }

    /// <summary>Profile-profile column score: average match/mismatch/gap over all cross pairs.</summary>
    private static double ProfileColumnScore(Profile p1, int col1, Profile p2, int col2, PairScoring scoring)
    {
        double total = 0;
        int pairs = p1.Rows.Count * p2.Rows.Count;
        for (int r1 = 0; r1 < p1.Rows.Count; r1++)
        {
            char a = p1.Rows[r1][col1];
            for (int r2 = 0; r2 < p2.Rows.Count; r2++)
            {
                char b = p2.Rows[r2][col2];
                if (a == '-' && b == '-')
                    total += 0; // gap-gap neutral
                else if (a == '-' || b == '-')
                    total += scoring.GapExtend;
                else
                    total += scoring.Sub(a, b);
            }
        }
        return total / pairs;
    }

    /// <summary>
    /// Score for aligning one profile column against an all-gap column of the other profile of
    /// <paramref name="otherRowCount"/> rows: average gap penalty over all cross pairs.
    /// Residue-vs-gap costs GapExtend; gap-vs-gap is neutral.
    /// </summary>
    private static double ProfileColumnVsGap(Profile p, int col, int otherRowCount, PairScoring scoring)
    {
        double total = 0;
        int pairs = p.Rows.Count * otherRowCount;
        for (int r = 0; r < p.Rows.Count; r++)
        {
            char a = p.Rows[r][col];
            // The other side is a gap for all otherRowCount rows.
            if (a != '-')
                total += scoring.GapExtend * otherRowCount;
            // a == '-' : gap-gap, neutral.
        }
        return total / pairs;
    }

    private const double ProfileScoreEpsilon = 1e-9;

    private static bool AreClose(double x, double y) => Math.Abs(x - y) <= ProfileScoreEpsilon;

    #endregion

    #region Multiple Sequence Alignment (Iterative Refinement — MUSCLE tree-dependent restricted partitioning)

    // Iterative refinement of a progressive MSA, added as a THIRD aligner alongside the star
    // (MultipleAlign) and progressive (MultipleAlignProgressive) methods. Both of those remain
    // byte-for-byte unchanged. This removes the single-pass "once a gap, always a gap" limitation
    // of the progressive seed: an early gap-placement error CAN now be corrected, because the
    // alignment is repeatedly re-split and re-aligned and a strictly better arrangement replaces it.
    //
    // Algorithm — MUSCLE Stage 3, "tree-dependent restricted partitioning" (Edgar 2004):
    //   Quoting Edgar (2004), Nucleic Acids Res 32(5):1792-1797, §"Stage 3, Refinement":
    //   3.1 "An edge is chosen from TREE2 (edges are visited in order of decreasing distance from
    //        the root). TREE2 is divided into two subtrees by deleting the edge."
    //   3.2/3.3 "The profile of the multiple alignment in each subtree is computed. A new multiple
    //        alignment is produced by re-aligning the two profiles."
    //   3.4 "If the SP score is improved, the new alignment is kept, otherwise it is discarded."
    //   "Steps 3.1-3.4 are repeated until convergence or until a user-defined limit is reached."
    //
    //   The same scheme also realises the Barton-Sternberg (1987) idea of returning to an existing
    //   alignment and re-aligning a part of it against the rest; here the partition is the
    //   guide-tree edge rather than a single removed sequence, which is exactly Edgar's restricted
    //   partitioning. Re-aligning the two sub-profiles uses the EXISTING profile-profile NW
    //   (AlignProfiles) and is accepted only on a non-decreasing sum-of-pairs (SP) score, so the
    //   refined alignment is provably never worse than the progressive seed.
    //
    //   SP score: "the MSA program optimizes the sum of all of the pairs of characters at each
    //   position in the alignment (the so-called sum of pair score)" (Wikipedia "Multiple sequence
    //   alignment"). Computed by the existing ComputeSumOfPairsScore (column-based: match/mismatch
    //   from the matrix, residue-gap = GapExtend, gap-gap neutral).
    //
    // Determinism: edges are enumerated in a fixed order (decreasing distance from the root, i.e.
    // internal nodes nearest the leaves first, with a deterministic tie-break); profile splitting,
    // re-projection and the accept-on-improvement rule are deterministic. No RNG is used.
    //
    // Sources retrieved this session (2026-06-23):
    //   - Edgar RC (2004) "MUSCLE: multiple sequence alignment with high accuracy and high
    //     throughput", Nucleic Acids Res 32(5):1792-1797.
    //     https://academic.oup.com/nar/article/32/5/1792/2380623 — Stage 3 steps 3.1-3.4 quoted above.
    //   - Barton GJ, Sternberg MJ (1987) "A strategy for the rapid multiple alignment of protein
    //     sequences...", J Mol Biol 198(2):327-337. https://pubmed.ncbi.nlm.nih.gov/3430611/ —
    //     iterative refinement: re-align a part of the alignment to the rest, iterate to a final
    //     alignment.
    //   - Wikipedia "Multiple sequence alignment", §Iterative methods / sum-of-pairs score
    //     https://en.wikipedia.org/wiki/Multiple_sequence_alignment — SP-score definition; iterative
    //     methods "can return to previously calculated ... sub-MSAs ... optimizing a general
    //     objective function".

    // Default cap on full refinement passes over all edges (Edgar's "user-defined limit"). A pass
    // that makes no accepted change means convergence, so this is only an upper bound.
    private const int DefaultMaxRefinementIterations = 16;

    /// <summary>
    /// Performs <b>iterative refinement</b> of a progressive multiple sequence alignment using
    /// MUSCLE-style tree-dependent restricted partitioning (Edgar 2004, Stage 3), removing the
    /// single-pass "once a gap, always a gap" limitation of
    /// <see cref="MultipleAlignProgressive(IEnumerable{DnaSequence}, ScoringMatrix?)"/>.
    /// </summary>
    /// <remarks>
    /// <para>
    /// Pipeline (Edgar 2004, Nucleic Acids Res 32(5):1792-1797; Barton &amp; Sternberg 1987):
    /// </para>
    /// <list type="number">
    /// <item>Build the progressive (Feng-Doolittle; UPGMA by default, <see cref="MsaOptions.GuideTree"/>) seed alignment and keep its guide tree.</item>
    /// <item>Visit each internal guide-tree edge (leaves-first, deterministic order). Deleting the
    /// edge partitions the sequences into two groups.</item>
    /// <item>Project the current alignment onto each group, drop columns that became all-gap, and
    /// realign the two sub-profiles with the existing profile-profile Needleman-Wunsch.</item>
    /// <item>Accept the re-alignment only if the sum-of-pairs (SP) score does not decrease.</item>
    /// <item>Repeat full passes until no edge yields an improvement (convergence) or the iteration
    /// cap is reached.</item>
    /// </list>
    /// <para>
    /// The result is therefore guaranteed to have an SP score <b>≥</b> that of the progressive seed,
    /// and is deterministic (no RNG; fixed edge order). The star
    /// <see cref="MultipleAlign(IEnumerable{DnaSequence}, ScoringMatrix?)"/> and progressive
    /// aligners are unchanged.
    /// </para>
    /// </remarks>
    /// <param name="sequences">Collection of sequences to align.</param>
    /// <param name="scoring">Scoring matrix (default: <see cref="SimpleDna"/>).</param>
    /// <param name="maxIterations">
    /// Upper bound on full refinement passes (Edgar's "user-defined limit"); must be positive.
    /// Defaults to <see cref="DefaultMaxRefinementIterations"/>. Convergence usually stops earlier.
    /// </param>
    /// <returns>Multiple alignment result (aligned rows, majority consensus, sum-of-pairs score).</returns>
    /// <exception cref="ArgumentNullException">When <paramref name="sequences"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException">When <paramref name="maxIterations"/> &lt; 1.</exception>
    public static MultipleAlignmentResult MultipleAlignIterative(
        IEnumerable<DnaSequence> sequences,
        ScoringMatrix? scoring = null,
        int maxIterations = DefaultMaxRefinementIterations) =>
        MultipleAlignIterative(sequences, scoring, new MsaOptions(), maxIterations);

    /// <summary>
    /// <see cref="MultipleAlignIterative(IEnumerable{DnaSequence}, ScoringMatrix?, int)"/> with additive
    /// <see cref="MsaOptions"/> (guide-tree distance source, <see cref="MsaDistanceMode"/>).
    /// </summary>
    /// <param name="sequences">Collection of sequences to align.</param>
    /// <param name="scoring">Scoring matrix (default: <see cref="SimpleDna"/>).</param>
    /// <param name="options">Options; <c>null</c> = defaults.</param>
    /// <param name="maxIterations">Upper bound on full refinement passes (&gt;= 1).</param>
    [System.Runtime.CompilerServices.OverloadResolutionPriority(1)]
    public static MultipleAlignmentResult MultipleAlignIterative(
        IEnumerable<DnaSequence> sequences,
        ScoringMatrix? scoring,
        MsaOptions? options,
        int maxIterations = DefaultMaxRefinementIterations)
    {
        ArgumentNullException.ThrowIfNull(sequences);
        ArgumentOutOfRangeException.ThrowIfLessThan(maxIterations, 1);
        return IterativeAlignCore(DnaRows(sequences), scoring ?? SimpleDna, options, maxIterations);
    }

    /// <summary>
    /// <see cref="MultipleAlignIterative(IEnumerable{DnaSequence}, ScoringMatrix?, MsaOptions?, int)"/> for <b>protein</b> sequences (or any alphabet of
    /// the <see cref="SubstitutionMatrix"/>): the same MUSCLE-style tree-dependent restricted partitioning of the progressive seed, accepted on a
    /// strictly better sum-of-pairs score under the matrix and gap model. See
    /// <see cref="MultipleAlign(IEnumerable{ProteinSequence}, MatrixScoring, MsaOptions?)"/> for the input contract.
    /// </summary>
    public static MultipleAlignmentResult MultipleAlignIterative(
        IEnumerable<ProteinSequence> sequences,
        MatrixScoring scoring,
        MsaOptions? options = null,
        int maxIterations = DefaultMaxRefinementIterations)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(maxIterations, 1);
        var (rows, ps) = PrepareMatrixMsa(sequences, scoring);
        return ScaleMsaResult(IterativeAlignCore(rows, ps, options, maxIterations), ps);
    }

    /// <summary><see cref="MultipleAlignIterative(IEnumerable{ProteinSequence}, MatrixScoring, MsaOptions?, int)"/> for raw sequence strings.</summary>
    public static MultipleAlignmentResult MultipleAlignIterative(
        IEnumerable<string> sequences,
        MatrixScoring scoring,
        MsaOptions? options = null,
        int maxIterations = DefaultMaxRefinementIterations)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(maxIterations, 1);
        var (rows, ps) = PrepareMatrixMsa(sequences, scoring);
        return ScaleMsaResult(IterativeAlignCore(rows, ps, options, maxIterations), ps);
    }

    private static MultipleAlignmentResult IterativeAlignCore(string[] rawSeqs, PairScoring effectiveScoring, MsaOptions? options, int maxIterations)
    {
        using var cancelScope = AlignmentCancellation.Enter(options?.CancellationToken ?? default);
        options ??= new MsaOptions();
        var gapModel = ResolveGapModel(options);

        var trivial = TrivialMultipleAlignment(rawSeqs);
        if (trivial != null)
            return trivial;

        int k = rawSeqs.Length;

        // Step 1: progressive seed (same pipeline as MultipleAlignProgressive) + retain guide tree.
        var distance = ComputeGuideDistances(rawSeqs, effectiveScoring, options.Distance, gapModel, options.GuideTree, options.QuickPairwise, options.ClustalWProtein);

        ProgressiveGuideNode root = BuildProgressiveGuideTree(k, distance, options.GuideTree);
        Profile current = AlignProfileSubtree(root, rawSeqs, effectiveScoring, gapModel);

        // Step 2-5: iterative refinement over the guide-tree edges.
        current = RefineByTreePartitioning(current, root, effectiveScoring, maxIterations, gapModel);

        // Reproject rows back into input order.
        var orderedRows = new string[k];
        for (int r = 0; r < current.Rows.Count; r++)
            orderedRows[current.SequenceIndices[r]] = current.Rows[r];

        var alignedList = orderedRows.ToList();
        int maxLen = alignedList.Count == 0 ? 0 : alignedList.Max(s => s.Length);
        string consensus = BuildConsensus(alignedList, maxLen, effectiveScoring.IsMatrix);
        int spScore = ComputeSumOfPairsScore(alignedList, effectiveScoring, gapModel);

        return new MultipleAlignmentResult(
            AlignedSequences: orderedRows,
            Consensus: consensus,
            TotalScore: spScore);
    }

    /// <summary>
    /// Repeatedly splits the current alignment at each internal guide-tree edge into two
    /// sub-profiles, realigns them with profile-profile NW, and keeps the result only when the SP
    /// score does not decrease (Edgar 2004 Stage 3, steps 3.1-3.4). Iterates full passes until a
    /// pass accepts no change (convergence) or <paramref name="maxIterations"/> is reached.
    /// </summary>
    private static Profile RefineByTreePartitioning(
        Profile current, ProgressiveGuideNode root, PairScoring scoring, int maxIterations, MsaGapModel gapModel)
    {
        // 3.1 enumerate the partitions induced by deleting each internal edge of the guide tree.
        // An internal edge is the edge above each non-root internal node; deleting it isolates that
        // node's leaf set as one group and the remaining leaves as the other. Edges are ordered by
        // decreasing distance from the root (deepest internal nodes — nearest the leaves — first),
        // with a deterministic tie-break on the sorted leaf-set, matching Edgar's visiting order.
        var partitions = EnumerateEdgePartitions(root);
        if (partitions.Count == 0)
            return current; // 2 leaves: the single edge is the root split, no internal edge to refine.

        // The SP score the existing pipeline reports for this alignment (input-order independent).
        int CurrentSpScore(Profile p) => ComputeSumOfPairsScore(p.Rows, scoring, gapModel);

        int bestScore = CurrentSpScore(current);

        for (int iter = 0; iter < maxIterations; iter++)
        {
            bool improvedThisPass = false;

            foreach (var group in partitions)
            {
                // 3.1 split the current alignment into the two leaf groups.
                var groupSet = group;
                var sideA = SplitProfile(current, idx => groupSet.Contains(idx));
                var sideB = SplitProfile(current, idx => !groupSet.Contains(idx));
                if (sideA.Rows.Count == 0 || sideB.Rows.Count == 0)
                    continue;

                // 3.2/3.3 realign the two profiles with the existing profile-profile NW.
                Profile candidate = gapModel == MsaGapModel.Affine
                    ? AlignProfilesAffine(sideA, sideB, scoring)
                    : AlignProfiles(sideA, sideB, scoring);
                int candidateScore = CurrentSpScore(candidate);

                // 3.4 keep only if the SP score does not decrease. Strict ">" would also be valid;
                // ">=" lets an equal-score rearrangement settle, but acceptance still never lowers
                // the score, so the monotonic-non-decreasing guarantee holds. To keep convergence
                // well-defined we only treat a STRICT improvement as "made progress".
                if (candidateScore > bestScore)
                {
                    current = candidate;
                    bestScore = candidateScore;
                    improvedThisPass = true;
                }
            }

            if (!improvedThisPass)
                break; // convergence: a full pass changed nothing.
        }

        return current;
    }

    /// <summary>
    /// Returns, for each internal (non-root) node of the guide tree, the set of input-sequence
    /// indices in its subtree — i.e. one side of the partition obtained by deleting that node's
    /// parent edge. Ordered by decreasing distance from the root (deepest first) with a
    /// deterministic tie-break, per Edgar (2004) "edges are visited in order of decreasing distance
    /// from the root".
    /// </summary>
    private static List<HashSet<int>> EnumerateEdgePartitions(ProgressiveGuideNode root)
    {
        var result = new List<(int Depth, List<int> Sorted, HashSet<int> Set)>();

        void Visit(ProgressiveGuideNode node, int depth, bool isRoot)
        {
            if (node.IsLeaf)
                return;

            // The root edge is not an internal edge (deleting it gives the trivial whole-vs-empty
            // split handled by the seed); only record proper internal nodes.
            if (!isRoot)
            {
                var leaves = new List<int>();
                CollectLeaves(node, leaves);
                leaves.Sort();
                result.Add((depth, leaves, new HashSet<int>(leaves)));
            }

            Visit(node.Left!, depth + 1, false);
            Visit(node.Right!, depth + 1, false);
        }

        Visit(root, 0, true);

        // Decreasing distance from the root = larger depth first; tie-break deterministically on the
        // sorted leaf list so the iteration order is fully reproducible.
        result.Sort((a, b) =>
        {
            int byDepth = b.Depth.CompareTo(a.Depth);
            if (byDepth != 0) return byDepth;
            return CompareIntLists(a.Sorted, b.Sorted);
        });

        return result.Select(r => r.Set).ToList();
    }

    private static int CompareIntLists(List<int> a, List<int> b)
    {
        int n = Math.Min(a.Count, b.Count);
        for (int i = 0; i < n; i++)
        {
            int c = a[i].CompareTo(b[i]);
            if (c != 0) return c;
        }
        return a.Count.CompareTo(b.Count);
    }

    private static void CollectLeaves(ProgressiveGuideNode node, List<int> leaves)
    {
        if (node.IsLeaf)
        {
            leaves.Add(node.LeafIndex);
            return;
        }
        CollectLeaves(node.Left!, leaves);
        CollectLeaves(node.Right!, leaves);
    }

    /// <summary>
    /// Projects the current alignment onto the rows whose input-sequence index satisfies
    /// <paramref name="keep"/>, then drops every column that is all-gap within the selected rows.
    /// The result is a valid profile of the selected subset; gaps inside retained columns are
    /// preserved (no existing column is edited — "once a gap, always a gap" within the sub-profile).
    /// </summary>
    private static Profile SplitProfile(Profile current, Func<int, bool> keep)
    {
        var rowIdx = new List<int>();
        for (int r = 0; r < current.Rows.Count; r++)
            if (keep(current.SequenceIndices[r]))
                rowIdx.Add(r);

        var sub = new Profile();
        if (rowIdx.Count == 0)
            return sub;

        int width = current.Width;
        // Identify columns that are NOT all-gap within the selected rows.
        var keptCols = new List<int>(width);
        for (int c = 0; c < width; c++)
        {
            bool allGap = true;
            foreach (int r in rowIdx)
            {
                if (current.Rows[r][c] != '-')
                {
                    allGap = false;
                    break;
                }
            }
            if (!allGap)
                keptCols.Add(c);
        }

        foreach (int r in rowIdx)
        {
            var sb = new StringBuilder(keptCols.Count);
            foreach (int c in keptCols)
                sb.Append(current.Rows[r][c]);
            sub.Rows.Add(sb.ToString());
            sub.SequenceIndices.Add(current.SequenceIndices[r]);
        }

        return sub;
    }

    #endregion

    #region Multiple Sequence Alignment (Consistency-based — T-Coffee)

    // Consistency-based multiple sequence alignment, the T-Coffee method (Notredame, Higgins &
    // Heringa 2000), added as a FOURTH aligner alongside the star (MultipleAlign), progressive
    // (MultipleAlignProgressive) and iterative (MultipleAlignIterative) methods. All three of those
    // remain byte-for-byte unchanged. This optimises the T-Coffee CONSISTENCY objective — a distinct
    // objective class from the fixed-matrix sum-of-pairs (SP) score the other methods optimise.
    //
    // Algorithm — Notredame, Higgins & Heringa (2000) "T-Coffee: A novel method for fast and accurate
    // multiple sequence alignment", J Mol Biol 302:205-217 (DOI 10.1006/jmbi.2000.4042). Figure 1
    // pipeline: primary library -> extension -> extended library -> progressive alignment.
    //
    //   1. PRIMARY LIBRARY (p.207). For every pair of sequences (i,j) compute pairwise alignments —
    //      a GLOBAL (Needleman-Wunsch, the ClustalW library in the paper) and a LOCAL (Smith-Waterman,
    //      the Lalign library). Each aligned residue pair (residue at position p in Si aligned to the
    //      residue at position q in Sj) "receives a weight equal to percent identity within the
    //      pair-wise alignment it comes from". The global and local libraries are combined by
    //      "signal addition": "If any pair is duplicated between the two libraries, it is merged into
    //      a single entry that has a weight equal to the sum of the two weights." Pairs that never
    //      occur have weight 0.
    //
    //   2. LIBRARY EXTENSION — the triplet consistency transformation (pp.208-209). For each library
    //      pair (Si.p, Sj.q), and for every intermediate sequence Sk, if Si.p aligns to Sk.r (weight
    //      W1) and Sk.r aligns to Sj.q (weight W2), the triplet through Sk contributes min(W1, W2):
    //      "we associate that alignment with a weight equal to the minimum of W1 and W2". The extended
    //      weight is "the sum of all the weights gathered through the examination of all the triplets
    //      involving that pair" PLUS the direct primary weight — the paper's worked example: direct
    //      88, through-C min(77,100)=77, extended = 88 + 77 = 165. Uninformative triplets contribute 0,
    //      so extension never lowers a weight.
    //
    //   3. PROGRESSIVE ALIGNMENT (pp.209-210). A guide tree (the shared
    //      BuildProgressiveGuideTree: UPGMA by default; the paper uses NJ, available as MsaGuideTree.NeighborJoining —
    //      merge order only, same library and objective)
    //      directs a progressive alignment. Profiles are aligned by dynamic programming whose
    //      position-specific column score is the SUM of the EXTENDED-LIBRARY weights over all
    //      cross-profile residue pairs in the two columns (group-vs-group: "the average library scores
    //      in each column of existing alignment are taken" — here summed over the cross pairs, which
    //      is the same objective up to a per-merge constant). No fixed substitution matrix is used and
    //      "gap-opening penalties and gap-extension penalties [are] set to zero" (p.210): a gap column
    //      scores 0. "Once a gap is introduced ... it cannot be shifted later" is enforced as in the
    //      progressive aligner (whole all-gap columns inserted, existing columns never edited).
    //
    // Determinism: no RNG. Primary alignments, the guide tree, the extension sum and the DP tie-breaks
    // are all deterministic.
    //
    // Sources retrieved this session (2026-06-23):
    //   - Notredame C, Higgins DG, Heringa J (2000) J Mol Biol 302:205-217. Full text fetched/read:
    //     https://web.stanford.edu/class/gene211/pdfs/Notredame-Tcoffee.pdf (Figures 1-2, pp.207-211).
    //   - T-Coffee Technical Documentation (min() rule for triplet legs):
    //     https://tcoffee.readthedocs.io/en/latest/tcoffee_technical_documentation.html
    //   - Gotoh O (1982) J Mol Biol 162:705-708 (the DP used in the progressive phase).

    /// <summary>
    /// Performs <b>consistency-based</b> multiple sequence alignment — the T-Coffee method of
    /// Notredame, Higgins &amp; Heringa (2000) — optimising the T-Coffee consistency objective rather
    /// than the fixed-matrix sum-of-pairs score used by
    /// <see cref="MultipleAlign(IEnumerable{DnaSequence}, ScoringMatrix?)"/>,
    /// <see cref="MultipleAlignProgressive(IEnumerable{DnaSequence}, ScoringMatrix?)"/> and
    /// <see cref="MultipleAlignIterative(IEnumerable{DnaSequence}, ScoringMatrix?, int)"/>
    /// (all of which are unchanged).
    /// </summary>
    /// <remarks>
    /// <para>
    /// Pipeline (Notredame et al. 2000, J Mol Biol 302:205-217):
    /// </para>
    /// <list type="number">
    /// <item>Build a <b>primary library</b> from all pairwise global (Needleman-Wunsch) and local
    /// (Smith-Waterman) alignments; each aligned residue pair is weighted by the pairwise percent
    /// identity, and global+local weights are summed (signal addition).</item>
    /// <item><b>Extend</b> the library by the triplet consistency transformation: each residue pair
    /// accumulates, over every intermediate sequence, the minimum of the two connecting weights,
    /// added to the direct weight.</item>
    /// <item><b>Progressively align</b> along a guide tree (UPGMA by default, <see cref="MsaOptions.GuideTree"/>) with dynamic programming whose
    /// column score is the summed extended-library weight (no substitution matrix, zero gap
    /// penalty).</item>
    /// </list>
    /// </remarks>
    /// <param name="sequences">Collection of sequences to align.</param>
    /// <param name="scoring">
    /// Scoring matrix used only to compute the pairwise alignments that seed the library
    /// (default: <see cref="SimpleDna"/>). The progressive phase uses the library weights, not this
    /// matrix, as the substitution score.
    /// </param>
    /// <returns>Multiple alignment result (aligned rows, majority consensus, sum-of-pairs score).</returns>
    /// <exception cref="ArgumentNullException">When <paramref name="sequences"/> is null.</exception>
    public static MultipleAlignmentResult MultipleAlignConsistency(
        IEnumerable<DnaSequence> sequences,
        ScoringMatrix? scoring = null) =>
        MultipleAlignConsistency(sequences, scoring, new MsaOptions());

    /// <summary>
    /// <see cref="MultipleAlignConsistency(IEnumerable{DnaSequence}, ScoringMatrix?)"/> with additive
    /// <see cref="MsaOptions"/> (guide-tree distance source, <see cref="MsaDistanceMode"/>).
    /// </summary>
    /// <param name="sequences">Collection of sequences to align.</param>
    /// <param name="scoring">Scoring matrix seeding the pairwise library (default: <see cref="SimpleDna"/>).</param>
    /// <param name="options">Options; <c>null</c> = defaults.</param>
    [System.Runtime.CompilerServices.OverloadResolutionPriority(1)]
    public static MultipleAlignmentResult MultipleAlignConsistency(
        IEnumerable<DnaSequence> sequences,
        ScoringMatrix? scoring,
        MsaOptions? options)
    {
        ArgumentNullException.ThrowIfNull(sequences);
        return ConsistencyAlignCore(DnaRows(sequences), scoring ?? SimpleDna, options);
    }

    /// <summary>
    /// <see cref="MultipleAlignConsistency(IEnumerable{DnaSequence}, ScoringMatrix?, MsaOptions?)"/> for <b>protein</b> sequences (or any alphabet of the
    /// <see cref="SubstitutionMatrix"/>): the primary library is built from the global and local pairwise alignments under the matrix (linear or
    /// affine, possibly fractional gap scores), weighted by percent identity; the extension, guide tree and progressive stages are those of the
    /// DNA overload (including <see cref="MsaConsistencyEngine.TCoffeeLinked"/>). See
    /// <see cref="MultipleAlign(IEnumerable{ProteinSequence}, MatrixScoring, MsaOptions?)"/> for the input contract.
    /// </summary>
    public static MultipleAlignmentResult MultipleAlignConsistency(
        IEnumerable<ProteinSequence> sequences,
        MatrixScoring scoring,
        MsaOptions? options = null)
    {
        var (rows, ps) = PrepareMatrixMsa(sequences, scoring);
        return ScaleMsaResult(ConsistencyAlignCore(rows, ps, options), ps);
    }

    /// <summary><see cref="MultipleAlignConsistency(IEnumerable{ProteinSequence}, MatrixScoring, MsaOptions?)"/> for raw sequence strings.</summary>
    public static MultipleAlignmentResult MultipleAlignConsistency(
        IEnumerable<string> sequences,
        MatrixScoring scoring,
        MsaOptions? options = null)
    {
        var (rows, ps) = PrepareMatrixMsa(sequences, scoring);
        return ScaleMsaResult(ConsistencyAlignCore(rows, ps, options), ps);
    }

    private static MultipleAlignmentResult ConsistencyAlignCore(string[] rawSeqs, PairScoring effectiveScoring, MsaOptions? options)
    {
        using var cancelScope = AlignmentCancellation.Enter(options?.CancellationToken ?? default);
        options ??= new MsaOptions();
        var gapModel = ResolveGapModel(options);

        var trivial = TrivialMultipleAlignment(rawSeqs);
        if (trivial != null)
            return trivial;

        int k = rawSeqs.Length;

        // Steps 1-2: primary library (global + local, signal-added) then triplet extension.
        if (options.PrimaryLocalAlignments < 1)
            throw new ArgumentOutOfRangeException(nameof(options), options.PrimaryLocalAlignments, "PrimaryLocalAlignments must be at least 1.");
        if (!Enum.IsDefined(options.SequenceWeights))
            throw new ArgumentOutOfRangeException(nameof(options), "Unknown sequence weighting.");
        var primaryLibrary = BuildPrimaryLibrary(rawSeqs, effectiveScoring, gapModel, options.PrimaryLocalAlignments);
        float[]? seqWeights = options.SequenceWeights == MsaSequenceWeights.TCoffee && k > 1 && rawSeqs.All(r => r.Length > 0)
            ? TcTCoffeeSequenceWeights(TcDistanceMatrix(rawSeqs, TcKindOf(rawSeqs, null), new TCoffeeProgramOptions(), false).Sim, k)
            : null;
        var library = new ExtendedLibrary(primaryLibrary, seqWeights);

        // Step 3a: guide tree (UPGMA by default, neighbour-joining on request) over the ClustalW distances (1 - identities / min(len1, len2)), reusing the progressive pipeline.
        var distance = ComputeGuideDistances(rawSeqs, effectiveScoring, options.Distance, gapModel, options.GuideTree, options.QuickPairwise, options.ClustalWProtein);
        ProgressiveGuideNode root = BuildProgressiveGuideTree(k, distance, options.GuideTree);

        if (options.ConsistencyEngine == MsaConsistencyEngine.TCoffeeLinked)
            return BuildTCoffeeLinkedResult(rawSeqs, primaryLibrary, root, effectiveScoring, gapModel, seqWeights);

        // Step 3b: progressive alignment driven by the EXTENDED-LIBRARY weights.
        ConsistencyProfile aligned = AlignConsistencySubtree(root, rawSeqs, library);

        // Reproject rows back to input order.
        var orderedRows = new string[k];
        for (int r = 0; r < aligned.Rows.Count; r++)
            orderedRows[aligned.SequenceIndices[r]] = aligned.Rows[r];

        var alignedList = orderedRows.ToList();
        int maxLen = alignedList.Count == 0 ? 0 : alignedList.Max(s => s.Length);
        string consensus = BuildConsensus(alignedList, maxLen, effectiveScoring.IsMatrix);
        // TotalScore is reported in the same SP currency as the sibling aligners for comparability.
        int spScore = ComputeSumOfPairsScore(alignedList, effectiveScoring, gapModel);

        return new MultipleAlignmentResult(
            AlignedSequences: orderedRows,
            Consensus: consensus,
            TotalScore: spScore);
    }

    // ----- Library representation -------------------------------------------------------------

    // A residue-pair key (Si.posI, Sj.posJ) with i &lt; j, used to store library weights. Stored
    // canonically with the smaller sequence index first so lookups are symmetric.
    private readonly record struct ResiduePairKey(int SeqA, int PosA, int SeqB, int PosB);

    private static ResiduePairKey MakeKey(int s1, int p1, int s2, int p2) =>
        s1 < s2 || (s1 == s2 && p1 <= p2)
            ? new ResiduePairKey(s1, p1, s2, p2)
            : new ResiduePairKey(s2, p2, s1, p1);

    /// <summary>
    /// The T-Coffee <b>extended</b> library (Notredame et al. 2000, pp.207-209): the primary library
    /// (global + local pairwise alignments, weighted by percent identity and combined by signal
    /// addition) plus the triplet consistency transformation, evaluated for <b>any</b> residue pair
    /// on demand (the paper: "Extension will be carried out on each pair of residues of A and B";
    /// the reference implementation <c>residue_pair_extended_list</c> in <c>lib/dp_lib/evaluate.c</c>
    /// likewise scores every pair the dynamic program asks for, not only pairs of the primary library).
    /// </summary>
    private sealed class ExtendedLibrary
    {
        private readonly Dictionary<ResiduePairKey, double> _primary;
        private readonly Dictionary<(int Seq, int Pos), Dictionary<(int Seq, int Pos), double>> _neighbours;
        private readonly Dictionary<ResiduePairKey, double> _cache = new();
        private readonly float[]? _seqWeights;

        public ExtendedLibrary(Dictionary<ResiduePairKey, double> primary, float[]? seqWeights = null)
        {
            _primary = primary;
            _neighbours = BuildNeighbourIndex(primary);
            _seqWeights = seqWeights;
        }

        /// <summary>
        /// Inverts the primary library into a per-residue adjacency: (seq,pos) -&gt; {(otherSeq,otherPos):
        /// weight}. Both directions of each stored pair are emitted so the triplet extension can walk from
        /// either endpoint.
        /// </summary>
        private static Dictionary<(int Seq, int Pos), Dictionary<(int Seq, int Pos), double>>
            BuildNeighbourIndex(Dictionary<ResiduePairKey, double> primary)
        {
            var index = new Dictionary<(int, int), Dictionary<(int, int), double>>();

            void Add(int s1, int p1, int s2, int p2, double w)
            {
                if (!index.TryGetValue((s1, p1), out var inner))
                {
                    inner = new Dictionary<(int, int), double>();
                    index[(s1, p1)] = inner;
                }
                inner[(s2, p2)] = w;
            }

            foreach (var (key, w) in primary)
            {
                Add(key.SeqA, key.PosA, key.SeqB, key.PosB, w);
                Add(key.SeqB, key.PosB, key.SeqA, key.PosA, w);
            }

            return index;
        }

        public double Primary(ResiduePairKey key) => _primary.GetValueOrDefault(key);

        /// <summary>
        /// Extended weight of (Si.p, Sj.q) = direct primary weight + Σ over every residue Sk.r of a
        /// third sequence that is in the primary library with both Si.p (W1) and Sj.q (W2) of
        /// min(W1, W2) (worked example of the paper: 88 + min(77, 100) = 165). Zero when the pair has
        /// neither direct nor triplet support.
        /// </summary>
        public double Extended(ResiduePairKey key)
        {
            if (_cache.TryGetValue(key, out double cached))
                return cached;

            double total = _primary.GetValueOrDefault(key);

            if (_neighbours.TryGetValue((key.SeqA, key.PosA), out var partnersOfA) &&
                _neighbours.TryGetValue((key.SeqB, key.PosB), out var partnersOfB))
            {
                // Intersect on (Sk, r); walk the smaller adjacency list.
                var (walk, other) = partnersOfA.Count <= partnersOfB.Count
                    ? (partnersOfA, partnersOfB)
                    : (partnersOfB, partnersOfA);

                foreach (var (partner, w1) in walk)
                {
                    if (partner.Seq == key.SeqA || partner.Seq == key.SeqB) continue; // third sequences only
                    if (other.TryGetValue(partner, out double w2))
                        total += Math.Min(w1, w2);
                }
            }

            // T-Coffee -seq_weight (re_weight_constraint_list, after the extension): the extended weight is multiplied by the smaller of the two sequence weights
            if (_seqWeights is not null)
                total *= Math.Min(_seqWeights[key.SeqA], _seqWeights[key.SeqB]);

            _cache[key] = total;
            return total;
        }
    }

    /// <summary>Primary library: global (NW) + local (SW) pairwise alignments of every sequence pair, signal-added.</summary>
    private static Dictionary<ResiduePairKey, double> BuildPrimaryLibrary(
        string[] seqs, PairScoring scoring, MsaGapModel gapModel = MsaGapModel.Linear, int localAlignments = 1)
    {
        bool affine = gapModel == MsaGapModel.Affine;
        var primary = new Dictionary<ResiduePairKey, double>();
        int k = seqs.Length;

        for (int i = 0; i < k; i++)
        {
            for (int j = i + 1; j < k; j++)
            {
                AddAlignmentToPrimary(primary, i, j, seqs[i], seqs[j],
                    affine ? AffineAlignCore(seqs[i], seqs[j], scoring, AlignmentType.Global) : GlobalAlignCore(seqs[i], seqs[j], scoring));

                if (seqs[i].Length > 0 && seqs[j].Length > 0)
                {
                    if (localAlignments > 1)
                    {
                        // Lalign library of the paper = lalign_id_pair: the Waterman-Eggert top-N non-intersecting local alignments, each weighted by its own percent identity
                        foreach (var local in LocalTopNCore(seqs[i], seqs[j], scoring, affine, localAlignments, 0))
                            AddLocalAlignmentToPrimary(primary, i, j, seqs[i], seqs[j], local);
                    }
                    else
                    {
                        var local = affine ? AffineAlignCore(seqs[i], seqs[j], scoring, AlignmentType.Local) : LocalAlignCore(seqs[i], seqs[j], scoring);
                        AddLocalAlignmentToPrimary(primary, i, j, seqs[i], seqs[j], local);
                    }
                }
            }
        }

        return primary;
    }

    /// <summary>
    /// Adds a global pairwise alignment to the primary library. Every aligned residue pair (both
    /// non-gap) gets the alignment's percent identity weight; duplicate pairs are signal-added.
    /// </summary>
    private static void AddAlignmentToPrimary(
        Dictionary<ResiduePairKey, double> primary,
        int i, int j, string si, string sj, AlignmentResult pa)
    {
        double weight = PercentIdentity(pa.AlignedSequence1, pa.AlignedSequence2);
        if (weight <= 0) return;

        int posI = 0, posJ = 0;
        string a = pa.AlignedSequence1, b = pa.AlignedSequence2;
        for (int c = 0; c < a.Length; c++)
        {
            char ca = a[c], cb = b[c];
            if (ca != '-' && cb != '-')
            {
                AddPairWeight(primary, i, posI, j, posJ, weight);
                posI++; posJ++;
            }
            else if (ca != '-') posI++;
            else if (cb != '-') posJ++;
        }
    }

    /// <summary>
    /// Adds a local pairwise alignment to the primary library. The local alignment covers a segment
    /// of each sequence starting at <c>StartPosition1</c>/<c>StartPosition2</c>; matched residue pairs
    /// in that segment are weighted by the local segment's percent identity and signal-added.
    /// </summary>
    private static void AddLocalAlignmentToPrimary(
        Dictionary<ResiduePairKey, double> primary,
        int i, int j, string si, string sj, AlignmentResult la)
    {
        double weight = PercentIdentity(la.AlignedSequence1, la.AlignedSequence2);
        if (weight <= 0) return;

        int posI = la.StartPosition1, posJ = la.StartPosition2;
        string a = la.AlignedSequence1, b = la.AlignedSequence2;
        for (int c = 0; c < a.Length; c++)
        {
            char ca = a[c], cb = b[c];
            if (ca != '-' && cb != '-')
            {
                AddPairWeight(primary, i, posI, j, posJ, weight);
                posI++; posJ++;
            }
            else if (ca != '-') posI++;
            else if (cb != '-') posJ++;
        }
    }

    private static void AddPairWeight(
        Dictionary<ResiduePairKey, double> lib, int s1, int p1, int s2, int p2, double weight)
    {
        var key = MakeKey(s1, p1, s2, p2);
        lib[key] = lib.GetValueOrDefault(key) + weight;
    }

    /// <summary>
    /// Percent identity of a pairwise alignment = ⌊100 × identical columns / aligned columns⌋, where an
    /// aligned column has a residue in both rows (gap columns are not counted). This is the library
    /// weight of T-Coffee's default <c>-weight=sim</c> (<c>get_seq_sim</c> in
    /// <c>lib/util_lib/aln_convertion_util.c</c>, mode sim1: <c>(sim*MAXID)/pos0</c>, integer division;
    /// Notredame et al. 2000 p.207: "a weight equal to percent identity within the pair-wise
    /// alignment it comes from").
    /// </summary>
    private static double PercentIdentity(string a, string b)
    {
        int aligned = 0;
        int identical = 0;
        for (int c = 0; c < a.Length; c++)
        {
            if (a[c] == '-' || b[c] == '-')
                continue;
            aligned++;
            if (a[c] == b[c])
                identical++;
        }

        return aligned == 0 ? 0 : (int)(PercentIdentityScale * identical) / aligned;
    }

    // Weight is expressed as a percentage (0-100), matching the paper's integer "Prim Weight" values.
    private const double PercentIdentityScale = 100.0;

    // ----- Consistency profile + progressive DP ----------------------------------------------

    /// <summary>
    /// A profile for consistency alignment: like <see cref="Profile"/> but every cell also records the
    /// original residue position in its source sequence (or -1 for a gap), so the library weight for a
    /// residue pair can be looked up during the position-specific DP.
    /// </summary>
    private sealed class ConsistencyProfile
    {
        public List<string> Rows { get; } = new();
        public List<int> SequenceIndices { get; } = new();
        // PositionMaps[r][col] = original residue index in sequence SequenceIndices[r], or -1 for gap.
        public List<int[]> PositionMaps { get; } = new();
        public int Width => Rows.Count == 0 ? 0 : Rows[0].Length;
    }

    private static ConsistencyProfile AlignConsistencySubtree(
        ProgressiveGuideNode node, string[] rawSeqs, ExtendedLibrary library)
    {
        if (node.IsLeaf)
        {
            int idx = node.LeafIndex;
            string seq = rawSeqs[idx];
            var leaf = new ConsistencyProfile();
            leaf.Rows.Add(seq);
            leaf.SequenceIndices.Add(idx);
            var map = new int[seq.Length];
            for (int p = 0; p < seq.Length; p++) map[p] = p;
            leaf.PositionMaps.Add(map);
            return leaf;
        }

        var left = AlignConsistencySubtree(node.Left!, rawSeqs, library);
        var right = AlignConsistencySubtree(node.Right!, rawSeqs, library);
        return AlignConsistencyProfiles(left, right, library);
    }

    /// <summary>
    /// Aligns two consistency profiles with the Needleman-Wunsch recurrence over columns, scoring a
    /// column pair by the SUM of extended-library weights over all cross-profile residue pairs (gap
    /// columns score 0 — zero gap penalty, Notredame et al. 2000 p.210). Existing columns are never
    /// edited; gaps are inserted as whole all-gap columns ("once a gap, always a gap").
    /// </summary>
    private static ConsistencyProfile AlignConsistencyProfiles(
        ConsistencyProfile p1, ConsistencyProfile p2, ExtendedLibrary library)
    {
        int m = p1.Width;
        int n = p2.Width;
        RequireProfileMatrixFits(m, n);

        var f = new double[m + 1, n + 1];
        // First row/column: aligning a column against a gap scores 0 (zero gap penalty).
        // (f[i,0] = f[0,j] = 0 for all i,j.)

        for (int i = 1; i <= m; i++)
        {
            AlignmentCancellation.ThrowIfRequested();
            for (int j = 1; j <= n; j++)
            {
                double diag = f[i - 1, j - 1] + ColumnLibraryScore(p1, i - 1, p2, j - 1, library);
                double up = f[i - 1, j];     // gap in p2: 0
                double left = f[i, j - 1];   // gap in p1: 0
                f[i, j] = Math.Max(diag, Math.Max(up, left));
            }
        }

        var cols1 = new List<int>();
        var cols2 = new List<int>();
        int ii = m, jj = n;
        while (ii > 0 || jj > 0)
        {
            if (ii > 0 && jj > 0)
            {
                double diag = f[ii - 1, jj - 1] + ColumnLibraryScore(p1, ii - 1, p2, jj - 1, library);
                if (AreClose(f[ii, jj], diag))
                {
                    cols1.Add(ii - 1);
                    cols2.Add(jj - 1);
                    ii--; jj--;
                    continue;
                }
            }

            if (ii > 0 && (jj == 0 || AreClose(f[ii, jj], f[ii - 1, jj])))
            {
                cols1.Add(ii - 1);
                cols2.Add(-1);
                ii--;
                continue;
            }

            cols1.Add(-1);
            cols2.Add(jj - 1);
            jj--;
        }

        cols1.Reverse();
        cols2.Reverse();

        var merged = new ConsistencyProfile();
        merged.SequenceIndices.AddRange(p1.SequenceIndices);
        merged.SequenceIndices.AddRange(p2.SequenceIndices);
        AppendConsistencyRows(merged, p1, cols1);
        AppendConsistencyRows(merged, p2, cols2);
        return merged;
    }

    private static void AppendConsistencyRows(
        ConsistencyProfile dest, ConsistencyProfile source, List<int> columnMap)
    {
        for (int r = 0; r < source.Rows.Count; r++)
        {
            string row = source.Rows[r];
            int[] srcMap = source.PositionMaps[r];
            var sb = new StringBuilder(columnMap.Count);
            var newMap = new int[columnMap.Count];
            for (int c = 0; c < columnMap.Count; c++)
            {
                int col = columnMap[c];
                if (col < 0)
                {
                    sb.Append('-');
                    newMap[c] = -1;
                }
                else
                {
                    sb.Append(row[col]);
                    newMap[c] = srcMap[col];
                }
            }
            dest.Rows.Add(sb.ToString());
            dest.PositionMaps.Add(newMap);
        }
    }

    /// <summary>
    /// Column-vs-column score = sum over all cross-profile residue pairs of their extended-library
    /// weight (0 when a residue is a gap or the pair is absent from the library). This realises the
    /// position-specific library scoring of Notredame et al. (2000): the DP maximises the summed
    /// consistency weight rather than a fixed substitution matrix.
    /// </summary>
    private static double ColumnLibraryScore(
        ConsistencyProfile p1, int col1, ConsistencyProfile p2, int col2,
        ExtendedLibrary library)
    {
        double total = 0;
        for (int r1 = 0; r1 < p1.Rows.Count; r1++)
        {
            int posA = p1.PositionMaps[r1][col1];
            if (posA < 0) continue; // gap contributes 0
            int seqA = p1.SequenceIndices[r1];
            for (int r2 = 0; r2 < p2.Rows.Count; r2++)
            {
                int posB = p2.PositionMaps[r2][col2];
                if (posB < 0) continue;
                int seqB = p2.SequenceIndices[r2];
                total += library.Extended(MakeKey(seqA, posA, seqB, posB));
            }
        }
        return total;
    }

    /// <summary>
    /// Test/diagnostic accessor: builds the T-Coffee extended library for the given sequences and
    /// returns the extended weight of the residue pair (Si.posI, Sj.posJ), and that pair's primary
    /// (pre-extension) weight, so the consistency transformation can be verified directly
    /// (Notredame et al. 2000, p.209).
    /// </summary>
    internal static (double Primary, double Extended) GetLibraryWeights(
        string[] seqs, int seqI, int posI, int seqJ, int posJ, ScoringMatrix? scoring = null)
    {
        var effectiveScoring = scoring ?? SimpleDna;
        var library = new ExtendedLibrary(BuildPrimaryLibrary(seqs, effectiveScoring));
        var key = MakeKey(seqI, posI, seqJ, posJ);
        return (library.Primary(key), library.Extended(key));
    }

    /// <summary>
    /// Test/diagnostic accessor (audit A2-16): the primary library of <see cref="MultipleAlignConsistency(IEnumerable{string}, MatrixScoring, MsaOptions?)"/> under <paramref name="options"/>
    /// (<see cref="MsaOptions.PrimaryLocalAlignments"/>) as (seqA, posA, seqB, posB) -> weight, and the sequence weights of <see cref="MsaOptions.SequenceWeights"/> (null = none).
    /// </summary>
    internal static (Dictionary<(int SeqA, int PosA, int SeqB, int PosB), double> Primary, float[]? SeqWeights) GetPrimaryLibrary(string[] seqs, ScoringMatrix? scoring, MsaOptions options)
    {
        var lib = BuildPrimaryLibrary(seqs, scoring ?? SimpleDna, ResolveGapModel(options), options.PrimaryLocalAlignments);
        float[]? w = options.SequenceWeights == MsaSequenceWeights.TCoffee && seqs.Length > 1
            ? TcTCoffeeSequenceWeights(TcDistanceMatrix(seqs, TcKindOf(seqs, null), new TCoffeeProgramOptions(), false).Sim, seqs.Length)
            : null;
        return (lib.ToDictionary(kv => (kv.Key.SeqA, kv.Key.PosA, kv.Key.SeqB, kv.Key.PosB), kv => kv.Value), w);
    }

    /// <summary>Extended weight of a residue pair under <paramref name="options"/> (diagnostic; see <see cref="GetPrimaryLibrary"/>).</summary>
    internal static double GetExtendedWeight(string[] seqs, int seqI, int posI, int seqJ, int posJ, ScoringMatrix? scoring, MsaOptions options)
    {
        var (_, w) = GetPrimaryLibrary(seqs, scoring, options);
        var library = new ExtendedLibrary(BuildPrimaryLibrary(seqs, scoring ?? SimpleDna, ResolveGapModel(options), options.PrimaryLocalAlignments), w);
        return library.Extended(MakeKey(seqI, posI, seqJ, posJ));
    }

    #endregion
}

