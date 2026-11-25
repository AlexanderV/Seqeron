using System.Runtime.CompilerServices;
using System.Text;

namespace Seqeron.Genomics.Alignment;


/// <summary>
/// Internal exact-integer scoring model shared by every dynamic program: either the scalar
/// <see cref="ScoringMatrix"/> (scale 1) or a <see cref="MatrixScoring"/> converted to integers in units of 1/Scale.
/// </summary>
public readonly struct PairScoring
{
    public readonly int Match;
    public readonly int Mismatch;
    public readonly int GapOpen;
    public readonly int GapExtend;

    /// <summary>The direction- and position-specific gap scores (<see cref="GapScoreProfile"/>) in this model's integer units; null = every gap run scores <see cref="GapOpen"/> / <see cref="GapExtend"/>.</summary>
    public readonly ScaledGapProfile? Gaps;

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

    private PairScoring(int gapOpen, int gapExtend, int scale, int[]? table, sbyte[]? index, int dim, bool protein, long maxAbs,
        ScaledGapProfile? gaps = null, int match = 0, int mismatch = 0)
    {
        _maxAbs = maxAbs;
        Protein = protein;
        Match = match;
        Mismatch = mismatch;
        GapOpen = gapOpen;
        GapExtend = gapExtend;
        Scale = scale;
        _table = table;
        _index = index;
        _dim = dim;
        Gaps = gaps;
    }

    // ---- position-dependent gap scores (GapScoreProfile): a vertical run (gap in sequence 2) stays in one matrix column j, a horizontal run
    // (gap in sequence 1) in one row i, so "start run" = column / row 0 and "end run" = the last column / row. Without a profile: the uniform scores.

    /// <summary>Score per position of a gap in sequence 2 (vertical move) in column <paramref name="j"/> of a matrix with <paramref name="n"/> columns after column 0.</summary>
    public int XExt(int j, int n) => Gaps is null ? GapExtend : XGap(j, n).Extend;

    /// <summary>Open score of a gap in sequence 2 in column <paramref name="j"/> (see <see cref="XExt"/>).</summary>
    public int XOpen(int j, int n) => Gaps is null ? GapOpen : XGap(j, n).Open;

    /// <summary>Interior score per position of a gap in sequence 2 (local alignments never start or end in a gap, so only the interior scores apply).</summary>
    public int XInteriorExt => Gaps is null ? GapExtend : Gaps.Interior1.Extend;

    /// <summary>Interior open score of a gap in sequence 2.</summary>
    public int XInteriorOpen => Gaps is null ? GapOpen : Gaps.Interior1.Open;

    /// <summary>Interior score per position of a gap in sequence 1.</summary>
    public int YInteriorExt => Gaps is null ? GapExtend : Gaps.Interior2.Extend;

    /// <summary>Interior open score of a gap in sequence 1.</summary>
    public int YInteriorOpen => Gaps is null ? GapOpen : Gaps.Interior2.Open;

    /// <summary>Score per position of a gap in sequence 1 (horizontal move) in row <paramref name="i"/> of a matrix with <paramref name="m"/> rows after row 0.</summary>
    public int YExt(int i, int m) => Gaps is null ? GapExtend : YGap(i, m).Extend;

    /// <summary>Open score of a gap in sequence 1 in row <paramref name="i"/> (see <see cref="YExt"/>).</summary>
    public int YOpen(int i, int m) => Gaps is null ? GapOpen : YGap(i, m).Open;

    private ScaledGap XGap(int j, int n)
    {
        var g = Gaps!;
        if (_windowed)
        {
            j += _colOffset;
            n = _totalCols;
        }

        if (j == 0) return g.Start1;
        return j == n ? g.End1 : g.Interior1;
    }

    private ScaledGap YGap(int i, int m)
    {
        var g = Gaps!;
        if (_windowed)
        {
            i += _rowOffset;
            m = _totalRows;
        }

        if (i == 0) return g.Start2;
        return i == m ? g.End2 : g.Interior2;
    }

    // ---- sub-problems of the linear-space program (A4-8): a sub-rectangle [i0, i1) x [j0, j1) of an m x n matrix is a matrix of its own for the dynamic program, but its column / row 0 and last column /
    // row are start / end runs of the WHOLE alignment only where they lie on the border of the whole matrix. A window carries that position (a vertical run stays in one column, a horizontal run in one row).

    private readonly bool _windowed;
    private readonly int _rowOffset, _colOffset, _totalRows, _totalCols;

    private PairScoring(in PairScoring o, ScaledGapProfile? gaps, bool windowed, int rowOffset, int colOffset, int totalRows, int totalCols)
    {
        _maxAbs = o._maxAbs;
        Protein = o.Protein;
        Match = o.Match;
        Mismatch = o.Mismatch;
        GapOpen = o.GapOpen;
        GapExtend = o.GapExtend;
        Scale = o.Scale;
        _table = o._table;
        _index = o._index;
        _dim = o._dim;
        Gaps = gaps;
        _windowed = windowed;
        _rowOffset = rowOffset;
        _colOffset = colOffset;
        _totalRows = totalRows;
        _totalCols = totalCols;
    }

    /// <summary>
    /// This model for the sub-matrix whose row 0 / column 0 are the rows / columns <paramref name="rowOffset"/> / <paramref name="colOffset"/> of a matrix of <paramref name="totalRows"/> x <paramref name="totalCols"/> cells after row / column 0:
    /// <see cref="XExt"/> / <see cref="XOpen"/> / <see cref="YExt"/> / <see cref="YOpen"/> then resolve the start / interior / end member of the whole matrix (the index arguments are relative to the window). Without a profile: unchanged.
    /// </summary>
    public PairScoring Window(int rowOffset, int colOffset, int totalRows, int totalCols) =>
        Gaps is null ? this : new PairScoring(this, Gaps, true, rowOffset, colOffset, totalRows, totalCols);

    /// <summary>
    /// This model with the interior gap scores at every position (a local alignment never starts or ends in a gap: Biopython's local mode ignores the left / right scores); without a profile: unchanged.
    /// </summary>
    public PairScoring InteriorOnly() =>
        Gaps is null ? this : new PairScoring(this, new ScaledGapProfile(Gaps.Interior1, Gaps.Interior1, Gaps.Interior1, Gaps.Interior2, Gaps.Interior2, Gaps.Interior2), false, 0, 0, 0, 0);

    /// <summary>
    /// The scalar match / mismatch model with a <see cref="GapScoreProfile"/>: every score is scaled to the smallest integer unit that makes the
    /// match, mismatch and all gap scores integral (fractional gap scores down to 1e-6, as for <see cref="MatrixScoring"/>).
    /// </summary>
    public static PairScoring Create(ScoringMatrix scoring, GapScoreProfile profile)
    {
        ArgumentNullException.ThrowIfNull(scoring);
        ArgumentNullException.ThrowIfNull(profile);
        _ = new PairScoring(scoring); // validates the interior gap scores like the uniform model
        decimal open = scoring.GapOpen, extend = scoring.GapExtend;
        var values = ProfileValues(profile, open, extend);
        var (scale, ints) = ScaleOf(values, scoring);
        var gaps = ints.Gaps;
        long maxAbs = Math.Max(Math.Max(Math.Abs((long)scoring.Match), Math.Abs((long)scoring.Mismatch)) * scale, gaps.MaxAbsStep);
        return new PairScoring((int)ints.Open, (int)ints.Extend, (int)scale, null, null, 0, false, maxAbs, gaps,
            checked((int)(scoring.Match * scale)), checked((int)(scoring.Mismatch * scale)));
    }

    /// <summary>The resolved (open, extend) decimals of the interior scores and of the six profile members (null member = interior).</summary>
    private static decimal[] ProfileValues(GapScoreProfile profile, decimal open, decimal extend)
    {
        (decimal O, decimal E) R(GapScore? g, string name)
        {
            if (g is not { } v) return (open, extend);
            // A profile gap is the Biopython pair (first position open + extend, further positions extend): Open itself may be a bonus, but no gap position may score above 0.
            if (!double.IsFinite(v.Open) || !double.IsFinite(v.Extend) || v.Extend > 0 || v.Open + v.Extend > 0)
                throw new ArgumentException($"{name} must be finite with Extend <= 0 and Open + Extend <= 0 (a gap position is a penalty), got Open {v.Open}, Extend {v.Extend}.");
            return ((decimal)v.Open, (decimal)v.Extend);
        }
        var i1 = R(profile.Interior1, nameof(GapScoreProfile.Interior1));
        var s1 = R(profile.Start1, nameof(GapScoreProfile.Start1));
        var e1 = R(profile.End1, nameof(GapScoreProfile.End1));
        var i2 = R(profile.Interior2, nameof(GapScoreProfile.Interior2));
        var s2 = R(profile.Start2, nameof(GapScoreProfile.Start2));
        var e2 = R(profile.End2, nameof(GapScoreProfile.End2));
        return [open, extend, i1.O, i1.E, s1.O, s1.E, e1.O, e1.E, i2.O, i2.E, s2.O, s2.E, e2.O, e2.E];
    }

    private static (long Scale, (long Open, long Extend, ScaledGapProfile Gaps) Ints) ScaleOf(decimal[] values, ScoringMatrix? scalar)
    {
        int decimals = 0;
        while (decimals <= MaxGapDecimals && !values.All(v => IsIntegral(v * Pow10(decimals))))
            decimals++;
        if (decimals > MaxGapDecimals)
            throw new ArgumentException($"Gap scores are resolved to 1e-{MaxGapDecimals}: {string.Join(", ", values)}.");
        long denominator = Pow10(decimals);
        long g = denominator;
        foreach (var v in values) g = Gcd(g, Math.Abs((long)(v * denominator)));
        long scale = denominator / g;
        long Q(decimal v) => (long)(v * denominator) / g;
        ScaledGap G(int k) => new(checked((int)Q(values[k])), checked((int)Q(values[k + 1])));
        var gaps = new ScaledGapProfile(G(2), G(4), G(6), G(8), G(10), G(12));
        return (scale, (Q(values[0]), Q(values[1]), gaps));
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
    public static PairScoring Create(MatrixScoring scoring, string a, string b, GapScoreProfile? profile = null)
        => Create(scoring, a, b, profile, allowPositiveGapOpen: false);

    /// <summary>
    /// <see cref="Create(MatrixScoring, string, string, GapScoreProfile?)"/>; with <paramref name="allowPositiveGapOpen"/> a POSITIVE gap-open score is legal: the SIM program of EMBOSS <c>matcher</c> charges a gap of length k
    /// <c>q + r·k</c> with q = gapopen − gapextend, which is negative whenever <c>-gapopen</c> is below <c>-gapextend</c> (audit A6-3: <c>matcher -gapopen 2 -gapextend 4</c> is a legal run).
    /// </summary>
    public static PairScoring Create(MatrixScoring scoring, string a, string b, GapScoreProfile? profile, bool allowPositiveGapOpen)
    {
        var matrix = scoring.Matrix ?? throw new ArgumentException("MatrixScoring.Matrix must not be null.", nameof(scoring));
        ValidateResidues(matrix, a, "sequence1");
        ValidateResidues(matrix, b, "sequence2");
        return Build(scoring, matrix, (long)a.Length + b.Length + 2, 1,
            $"sequence lengths {a.Length} and {b.Length}", profile, allowPositiveGapOpen);
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

    private static PairScoring Build(MatrixScoring scoring, SubstitutionMatrix matrix, long pathLength, long pairs, string sizeText, GapScoreProfile? profile = null, bool allowPositiveGapOpen = false)
    {
        decimal open = ToGapDecimal(scoring.GapOpen, nameof(MatrixScoring.GapOpen), allowPositiveGapOpen);
        decimal extend = ToGapDecimal(scoring.GapExtend, nameof(MatrixScoring.GapExtend));

        // Smallest integer scale that makes every gap score integral (0.5 -> 2, 0.1 -> 10, integers -> 1); a profile adds its six members.
        var values = ProfileValues(profile ?? new GapScoreProfile(), open, extend);
        // R14-8: the cells of a fractional matrix (BENNER22, GONNET1992, JONES, ...) join the gap scores in the search for the smallest integer unit, so the program scores them exactly
        if (matrix.IsFractional)
        {
            var cells = new HashSet<decimal>();
            for (int i = 0; i < matrix.Alphabet.Length; i++)
            {
                for (int j = 0; j < matrix.Alphabet.Length; j++)
                    cells.Add(matrix.ValueAt(i, j));
            }

            values = [.. values, .. cells];
        }

        var (scale, ints) = ScaleOf(values, null);
        long gapOpen = ints.Open;
        long gapExtend = ints.Extend;
        ScaledGapProfile? scaledGaps = profile is null ? null : ints.Gaps;

        long maxAbs = Math.Max((long)Math.Ceiling(matrix.MaxAbsValue * scale), Math.Max(Math.Abs(gapOpen) + Math.Abs(gapExtend), scaledGaps?.MaxAbsStep ?? 0));
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
                table[i * dim + j] = checked((int)(matrix.ValueAt(i, j) * scale));

        return new PairScoring((int)gapOpen, (int)gapExtend, (int)scale, table, index, dim, !IsNucleotideAlphabet(matrix.Alphabet), maxAbs, scaledGaps);
    }

    private static ArgumentException RangeError(long scale, string sizeText, string paramName) =>
        new($"Score range too large for the integer dynamic program (scale {scale}, {sizeText}); use coarser gap scores or shorter sequences.", paramName);

    /// <summary>
    /// The one integer-range check of every dynamic program (scalar and matrix model alike): the score of any path through sequences of
    /// <paramref name="length1"/> and <paramref name="length2"/> residues, <c>max|step| · (length1 + length2 + 2)</c>, must stay below
    /// 2^28 (the affine −∞ sentinel is <c>int.MinValue / 4</c>), otherwise the DP would silently wrap around. The matrix model has been
    /// checked for the same bound (and the whole set of pairs) by <see cref="PairScoring.Create(MatrixScoring, string, string, GapScoreProfile?)"/>; this call repeats the pairwise
    /// part at the point where the DP is about to run, so the scalar <see cref="ScoringMatrix"/> is guarded identically.
    /// </summary>
    /// <exception cref="ArgumentException">"Score range too large for the integer dynamic program ...".</exception>
    public void EnsureRange(long length1, long length2)
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

    private static decimal ToGapDecimal(double value, string name, bool allowPositive = false)
    {
        if (double.IsNaN(value) || double.IsInfinity(value) || (value > 0 && !allowPositive))
            throw new ArgumentException($"{name} must be a finite score{(allowPositive ? "" : " <= 0 (a penalty)")}, got {value}.", name);
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

/// <summary>One gap score (open added once per run, score per position) in the integer units of a <see cref="PairScoring"/>.</summary>
public readonly record struct ScaledGap(int Open, int Extend);

/// <summary>The six direction / position gap scores of a <see cref="GapScoreProfile"/> in integer units (members already resolved against the interior score).</summary>
public sealed record ScaledGapProfile(ScaledGap Interior1, ScaledGap Start1, ScaledGap End1, ScaledGap Interior2, ScaledGap Start2, ScaledGap End2)
{
    /// <summary>Largest |score| of one gap position (open + extend) over the six members.</summary>
    public long MaxAbsStep => new[] { Interior1, Start1, End1, Interior2, Start2, End2 }.Max(g => Math.Abs(g.Open) + Math.Abs(g.Extend));
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
    [System.Runtime.CompilerServices.OverloadResolutionPriority(1)]
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
    [System.Runtime.CompilerServices.OverloadResolutionPriority(1)]
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
    public static AlignmentResult GlobalAlignCore(
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
    public static AlignmentResult LinearAlignCore(
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

        // Gap score per position of a vertical move in each column (up) / a horizontal move in each row (left): uniform d, or the position-dependent
        // scores of a GapScoreProfile (start run = column / row 0, end run = last column / row, otherwise interior).
        var upGap = new int[n + 1];
        for (int j = 0; j <= n; j++)
            upGap[j] = scoring.XExt(j, n);

        for (int i = 0; i <= m; i++)
            score[i, 0] = ends.HasFlag(EndGaps.FreeStart1) ? 0 : i * upGap[0];
        for (int j = 0; j <= n; j++)
            score[0, j] = ends.HasFlag(EndGaps.FreeStart2) ? 0 : j * scoring.YExt(0, m);

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
            int leftGap = scoring.YExt(i, m);
            for (int j = 1; j <= n; j++)
            {
                int matchScore = scoring.Sub(c1, seq2[j - 1]);

                // A forbidden aligned pair (Waterman-Eggert declumping, linear-space top-N) has no diagonal move.
                int diag = forbiddenRow.Test(colOffset + j) ? ForbiddenPairs.NoDiagonal : score[i - 1, j - 1] + matchScore;
                int up = score[i - 1, j] + upGap[j];
                int left = score[i, j - 1] + leftGap;

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
    public static AlignmentResult AffineAlignCore(
        string seq1, string seq2, PairScoring scoring, AlignmentType type, EndGaps ends = EndGaps.None,
        int startState = StateMatch, int endState = -1, bool forceFullMatrix = false,
        ForbiddenPairs? forbidden = null, int rowOffset = 0, int colOffset = 0)
    {
        bool local = type == AlignmentType.Local;
        if (local) ends = EndGaps.None;
        scoring.EnsureRange(seq1.Length, seq2.Length);
        if (!forceFullMatrix && ExceedsFullMatrixLimit(seq1.Length, seq2.Length))
        {
            return LinearSpaceAlign(seq1, seq2, scoring, affine: true, local, ends, LinearSpaceBaseCells, default);
        }
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

        return AffineTraceback(seq1, seq2, new FullAffineStates(mm, xx, yy), scoring, local, ends, ci, cj, state, best);
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
        // Column 0 (leading gap in sequence 2) and row 0 (leading gap in sequence 1) carry the start-run scores of a GapScoreProfile.
        int e = scoring.XExt(0, n);
        int oe = scoring.XOpen(0, n) + e;
        int ey = scoring.YExt(0, m);
        int oey = scoring.YOpen(0, m) + ey;
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
                : Math.Max(AffineNegInf, Max3(mm[0, j - 1] + oey, yy[0, j - 1] + ey, xx[0, j - 1] + oey));
        }
    }

    /// <summary>
    /// Fills the interior of the three Gotoh matrices row by row with <see cref="AffineCell(int[,], int[,], int[,], int, int, int, int, int, int, bool)"/> (local: the zero floor of
    /// Smith-Waterman); returns the first (row-major) maximum of M for a local alignment (0, 0, 0 otherwise / when no positive cell).
    /// </summary>
    private static (int Best, int I, int J) AffineFill(
        int[,] mm, int[,] xx, int[,] yy, string seq1, string seq2, PairScoring scoring, bool local,
        ForbiddenPairs? forbidden = null, int rowOffset = 0, int colOffset = 0)
    {
        int m = seq1.Length;
        int n = seq2.Length;
        int startFloor = local ? 0 : AffineNegInf;
        int bestLocal = 0, bestI = 0, bestJ = 0;
        // Gap scores of a vertical run per column and of a horizontal run per row: uniform, or the GapScoreProfile's (a local alignment never starts or ends in a gap: interior scores only).
        var xe = new int[n + 1];
        var xoe = new int[n + 1];
        for (int j = 0; j <= n; j++)
        {
            xe[j] = local ? scoring.XInteriorExt : scoring.XExt(j, n);
            xoe[j] = (local ? scoring.XInteriorOpen : scoring.XOpen(j, n)) + xe[j];
        }

        for (int i = 1; i <= m; i++)
        {
            AlignmentCancellation.ThrowIfRequested();
            char c1 = seq1[i - 1];
            var forbiddenRow = new ForbiddenRow(forbidden, rowOffset + i);
            int ye = local ? scoring.YInteriorExt : scoring.YExt(i, m);
            int yoe = (local ? scoring.YInteriorOpen : scoring.YOpen(i, m)) + ye;
            for (int j = 1; j <= n; j++)
            {
                AffineCell(mm, xx, yy, i, j, scoring.Sub(c1, seq2[j - 1]), startFloor, xe[j], xoe[j], ye, yoe, forbiddenRow.Test(colOffset + j));

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

    /// <summary>Read access to the three Gotoh matrices (M, X, Y) of the traceback: the full matrices of the aligners, or the diagonal band of the banded aligners (cells outside the band are unreachable).</summary>
    internal interface IAffineStates
    {
        /// <summary>M(i, j).</summary>
        int M(int i, int j);

        /// <summary>X(i, j).</summary>
        int X(int i, int j);

        /// <summary>Y(i, j).</summary>
        int Y(int i, int j);
    }

    /// <summary>The full matrices of the aligners.</summary>
    internal readonly struct FullAffineStates(int[,] mm, int[,] xx, int[,] yy) : IAffineStates
    {
        public int M(int i, int j) => mm[i, j];

        public int X(int i, int j) => xx[i, j];

        public int Y(int i, int j) => yy[i, j];
    }

    /// <summary>
    /// Gotoh traceback from state <paramref name="state"/> at cell (<paramref name="ci"/>, <paramref name="cj"/>) of the
    /// filled matrices; ties M &gt; X &gt; Y, a local alignment stops when the predecessor value is 0. Shared by
    /// <see cref="AffineAlignCore"/> and the Waterman-Eggert declumper.
    /// </summary>
    private static AlignmentResult AffineTraceback<TStates>(
        string seq1, string seq2, TStates st, PairScoring scoring,
        bool local, EndGaps ends, int ci, int cj, byte state, int best)
        where TStates : struct, IAffineStates
    {
        int m = seq1.Length;
        int n = seq2.Length;
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
                int prev = st.M(ci, cj) - (scoring.Sub(seq1[ci - 1], seq2[cj - 1]));
                chars1.Add(seq1[ci - 1]);
                chars2.Add(seq2[cj - 1]);
                ci--; cj--;
                if (local && prev == 0) break; // Smith-Waterman: the local alignment starts here
                state = PickState(prev, st.M(ci, cj), StateMatch, st.X(ci, cj), StateGapInSeq2, StateGapInSeq1);
            }
            else if (state == StateGapInSeq2)
            {
                int cur = st.X(ci, cj);
                int e = local ? scoring.XInteriorExt : scoring.XExt(cj, n);
                int oe = (local ? scoring.XInteriorOpen : scoring.XOpen(cj, n)) + e;
                chars1.Add(seq1[ci - 1]);
                chars2.Add('-');
                ci--;
                if (ci == 0 && cj == 0) break; // leading gap reached the origin
                state = PickState(cur, st.M(ci, cj) + oe, StateMatch, st.X(ci, cj) + e, StateGapInSeq2, StateGapInSeq1);
            }
            else // StateGapInSeq1
            {
                int cur = st.Y(ci, cj);
                int e = local ? scoring.YInteriorExt : scoring.YExt(ci, m);
                int oe = (local ? scoring.YInteriorOpen : scoring.YOpen(ci, m)) + e;
                chars1.Add('-');
                chars2.Add(seq2[cj - 1]);
                cj--;
                if (ci == 0 && cj == 0) break;
                state = PickState(cur, st.M(ci, cj) + oe, StateMatch, st.Y(ci, cj) + e, StateGapInSeq1, StateGapInSeq2);
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
        int[,] mm, int[,] xx, int[,] yy, int i, int j, int sub, int startFloor, int e, int oe, bool forbidden) =>
        AffineCell(mm, xx, yy, i, j, sub, startFloor, e, oe, e, oe, forbidden);

    /// <summary><see cref="AffineCell(int[,], int[,], int[,], int, int, int, int, int, int, bool)"/> with separate gap scores for a gap in sequence 2 (X: <paramref name="ex"/> / <paramref name="oex"/>) and in sequence 1 (Y).</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    internal static void AffineCell(
        int[,] mm, int[,] xx, int[,] yy, int i, int j, int sub, int startFloor, int ex, int oex, int ey, int oey, bool forbidden)
    {
        int m = forbidden
            ? AffineNegInf
            : Math.Max(startFloor, Max3(mm[i - 1, j - 1], xx[i - 1, j - 1], yy[i - 1, j - 1])) + sub;
        int x = Max3(mm[i - 1, j] + oex, xx[i - 1, j] + ex, yy[i - 1, j] + oex);
        int y = Max3(mm[i, j - 1] + oey, yy[i, j - 1] + ey, xx[i, j - 1] + oey);

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

    public static AlignmentResult LocalAlignCore(string seq1, string seq2, PairScoring scoring, bool forceFullMatrix = false)
    {
        int m = seq1.Length;
        int n = seq2.Length;
        scoring.EnsureRange(m, n);

        if (!forceFullMatrix && ExceedsFullMatrixLimit(m, n))
        {
            return LinearSpaceAlign(seq1, seq2, scoring, affine: false, local: true, EndGaps.None, LinearSpaceBaseCells, default);
        }

        var score = new int[m + 1, n + 1];
        int maxScore = 0;
        int maxI = 0, maxJ = 0;

        // Fill the matrix (with zero floor for local alignment)
        for (int i = 1; i <= m; i++)
        {
            AlignmentCancellation.ThrowIfRequested();
            for (int j = 1; j <= n; j++)
            {
                score[i, j] = LocalLinearCell(score, i, j, scoring.Sub(seq1[i - 1], seq2[j - 1]), scoring.XInteriorExt, scoring.YInteriorExt, forbidden: false);

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
    internal static int LocalLinearCell(int[,] score, int i, int j, int sub, int gap, bool forbidden) =>
        LocalLinearCell(score, i, j, sub, gap, gap, forbidden);

    /// <summary><see cref="LocalLinearCell(int[,], int, int, int, int, bool)"/> with a separate gap score for a gap in sequence 2 (up) and in sequence 1 (left).</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    internal static int LocalLinearCell(int[,] score, int i, int j, int sub, int gapUp, int gapLeft, bool forbidden)
    {
        int best = Math.Max(score[i - 1, j] + gapUp, score[i, j - 1] + gapLeft);
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
            else if (score[i, j] == score[i - 1, j] + scoring.XInteriorExt)
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

            if (i > 0 && (j == 0 || score[i, j] == score[i - 1, j] + scoring.XExt(j, seq2.Length)))
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
    public const char IdentityMark = '|';
    public const char SimilarityMark = ':';
    private const char MismatchMark = '.';
    public const char GapMark = ' ';
    public const char GapChar = '-';

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
    public static Func<char, char, (bool Similar, char Mark)> MatrixClassifier(SubstitutionMatrix matrix) =>
        (a, b) =>
        {
            bool self1 = matrix.IsPositive(a, a);
            bool self2 = matrix.IsPositive(b, b);
            bool pair = matrix.IsPositive(a, b);
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

    public static string[] DnaRows(IEnumerable<DnaSequence> sequences) => DnaStringsNoNull(sequences, nameof(sequences)).ToArray();

    /// <summary>Upper-cased rows and the exact integer model of a <see cref="MatrixScoring"/> for a whole set of sequences (validates residues, gaps, range).</summary>
    public static (string[] Rows, PairScoring Scoring) PrepareMatrixMsa(IEnumerable<string> sequences, MatrixScoring scoring)
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

    public static (string[] Rows, PairScoring Scoring) PrepareMatrixMsa(IEnumerable<ProteinSequence> sequences, MatrixScoring scoring)
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
    public static MultipleAlignmentResult ScaleMsaResult(MultipleAlignmentResult r, PairScoring ps)
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
    public static MultipleAlignmentResult? TrivialMultipleAlignment(string[] seqList)
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
    public static string BuildConsensus(List<string> aligned, int length, bool anySymbol = false)
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
    public static int ComputeSumOfPairsScore(IReadOnlyList<string> aligned, PairScoring scoring, MsaGapModel gapModel)
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



}

