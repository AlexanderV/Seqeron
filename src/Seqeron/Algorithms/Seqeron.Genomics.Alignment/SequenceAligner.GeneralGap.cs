using Seqeron.Genomics.Core;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Alignment;

/// <summary>
/// General gap-cost pairwise alignment (Waterman–Smith–Beyer 1976): the score of a gap of length k is an arbitrary
/// function <c>w(k)</c> — concave (logarithmic), convex (quadratic), capped, step-wise, … — instead of the linear/affine
/// models of the other aligners.
/// </summary>
public static partial class SequenceAligner
{
    #region General gap-cost alignment (Waterman–Smith–Beyer)

    // Model (Waterman, Smith & Beyer 1976, Bull Math Biol 38:1-12; Biopython PairwiseAligner with a callable gap_score, which
    // selects its "Waterman-Smith-Beyer global/local alignment algorithm"): every maximal run of k gap positions scores w(k),
    // runs of the two kinds (residue of sequence 1 against a gap / residue of sequence 2 against a gap) may be adjacent and are
    // charged separately, a gap run is never split (so a convex w is honoured, unlike the single-matrix textbook recurrence
    // H(i,j) = max(..., max_k H(i-k,j) + w(k), ...), which lets a gap be cut into cheaper pieces).
    //   M(i,j) = max(M, X, Y)(i-1,j-1) + s(a_i, b_j)                          (local: max(0, ...))
    //   X(i,j) = max_{1<=k<=i} max(M, Y)(i-k, j) + w(k)                         (a_{i-k+1..i} against a gap)
    //   Y(i,j) = max_{1<=k<=j} max(M, X)(i, j-k) + w(k)                         (b_{j-k+1..j} against a gap)
    // O(m·n·(m+n)) time, O(m·n) memory (back-pointers). The end-gap policy (EndGaps) is the same as for the affine aligners:
    // a free start has M(i,0) = 0 / M(0,j) = 0, a penalised start X(i,0) = w(i) / Y(0,j) = w(j), a free end is the maximum
    // over the last column / row.
    //
    // Fast paths: w is evaluated once for k = 1 .. max(m, n). When all those values lie on a line w(k) = w(1) + (k-1)·e with
    // e <= 0 and w(1) - e <= 0 (a linear or affine gap function with a representable gap score; the scalar scoring needs
    // integer values), the alignment is delegated to the O(m·n) linear (open == 0) or affine dynamic program — the same score; the
    // returned alignment is one of the co-optimal ones under that program's tie order, which can differ from the generic
    // program's (diagonal > up > left ties of the generic program: M > X > Y, then the smallest k).

    private sealed record GapSubstitution(Func<char, char, double> Score, ScoringMatrix? Scalar, SubstitutionMatrix? Matrix);

    private static GapSubstitution ScalarSubstitution(ScoringMatrix? scoring)
    {
        var sc = scoring ?? SimpleDna;
        return new GapSubstitution((a, b) => a == b ? sc.Match : sc.Mismatch, sc, null);
    }

    private static GapSubstitution MatrixSubstitution(SubstitutionMatrix matrix)
    {
        ArgumentNullException.ThrowIfNull(matrix);
        return new GapSubstitution((a, b) => matrix.GetScore(a, b), null, matrix);
    }

    // ----- global -----

    /// <summary>
    /// Global alignment (Needleman-Wunsch with Waterman-Smith-Beyer gaps) in which a gap of length k scores
    /// <paramref name="gapScore"/>(k) — an arbitrary function (concave, convex, capped, step, ...). Residues score
    /// <see cref="ScoringMatrix.Match"/> / <see cref="ScoringMatrix.Mismatch"/> (the gap fields of <paramref name="scoring"/> are not used).
    /// </summary>
    /// <param name="sequence1">First sequence.</param>
    /// <param name="sequence2">Second sequence.</param>
    /// <param name="gapScore">Score of a gap of length k &gt;= 1 (a score, negative = penalty, as <see cref="ScoringMatrix.GapExtend"/>).
    /// Must be a deterministic pure function; it is evaluated once for each k = 1 .. max(length1, length2). NaN and +∞ are rejected, −∞
    /// forbids a gap length (an alignment must still exist).</param>
    /// <param name="scoring">Match/mismatch scores (null = <see cref="SimpleDna"/>).</param>
    /// <remarks>
    /// <para>Every maximal gap run is charged once, w(k) for a run of k positions (a run of the other kind directly next to it is a
    /// separate run); end gaps are ordinary gaps. This is the model of Biopython 1.88 <c>PairwiseAligner</c> with a callable
    /// <c>gap_score</c> (Waterman-Smith-Beyer algorithm), cross-checked on random pairs with concave, convex, capped and step gap
    /// functions (see Evidence ALIGN-GLOBAL-001). Complexity O(m·n·(m+n)) time, O(m·n) memory; a linear or affine <paramref name="gapScore"/>
    /// is recognised and aligned with the O(m·n) dynamic programs of <see cref="GlobalAlign(string, string, ScoringMatrix?)"/> /
    /// <see cref="GlobalAlignAffine(string, string, ScoringMatrix?)"/> (identical score; a co-optimal alignment). Ties of the general
    /// program: diagonal M &gt; X &gt; Y, then the shortest gap run; empty/null contract as the other aligners (null throws, empty is the all-gap
    /// alignment of cost w(length)).</para>
    /// <para>Source: Waterman MS, Smith TF, Beyer WA (1976) Some biological sequence metrics. Adv Math 20:367-387 / Bull Math Biol 38:1-12;
    /// Gusfield 1997 §11.2 (general gap weights); Durbin et al. 1998 §2.4.</para>
    /// </remarks>
    /// <exception cref="ArgumentNullException">A sequence or <paramref name="gapScore"/> is null.</exception>
    /// <exception cref="ArgumentException"><paramref name="gapScore"/> returns NaN/+∞, or forbids every alignment.</exception>
    public static AlignmentResult GlobalAlignGeneralGap(
        string sequence1, string sequence2, Func<int, double> gapScore, ScoringMatrix? scoring = null) =>
        AlignGeneralGap(sequence1, sequence2, gapScore, ScalarSubstitution(scoring), AlignmentType.Global, EndGaps.None);

    /// <summary><see cref="GlobalAlignGeneralGap(string, string, Func{int, double}, ScoringMatrix?)"/> for DNA sequences.</summary>
    public static AlignmentResult GlobalAlignGeneralGap(
        DnaSequence sequence1, DnaSequence sequence2, Func<int, double> gapScore, ScoringMatrix? scoring = null)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        return AlignGeneralGap(sequence1.Sequence, sequence2.Sequence, gapScore, ScalarSubstitution(scoring), AlignmentType.Global, EndGaps.None);
    }

    /// <summary>General-gap global alignment scored with a <see cref="SubstitutionMatrix"/> (BLOSUM62, PAM250, EDNAFULL, ...); a residue outside
    /// its alphabet throws <see cref="ArgumentException"/>. <see cref="AlignmentResult.ExactScore"/> is the real-valued score.</summary>
    public static AlignmentResult GlobalAlignGeneralGap(
        string sequence1, string sequence2, Func<int, double> gapScore, SubstitutionMatrix matrix) =>
        AlignGeneralGap(sequence1, sequence2, gapScore, MatrixSubstitution(matrix), AlignmentType.Global, EndGaps.None);

    /// <summary><see cref="GlobalAlignGeneralGap(string, string, Func{int, double}, SubstitutionMatrix)"/> for protein sequences.</summary>
    public static AlignmentResult GlobalAlignGeneralGap(
        ProteinSequence sequence1, ProteinSequence sequence2, Func<int, double> gapScore, SubstitutionMatrix matrix)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        return AlignGeneralGap(sequence1.Sequence, sequence2.Sequence, gapScore, MatrixSubstitution(matrix), AlignmentType.Global, EndGaps.None);
    }

    /// <summary><see cref="GlobalAlignGeneralGap(string, string, Func{int, double}, SubstitutionMatrix)"/> for DNA sequences.</summary>
    public static AlignmentResult GlobalAlignGeneralGap(
        DnaSequence sequence1, DnaSequence sequence2, Func<int, double> gapScore, SubstitutionMatrix matrix)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        return AlignGeneralGap(sequence1.Sequence, sequence2.Sequence, gapScore, MatrixSubstitution(matrix), AlignmentType.Global, EndGaps.None);
    }

    // ----- local -----

    /// <summary>
    /// Local alignment (Smith-Waterman with Waterman-Smith-Beyer gaps): a gap of length k scores <paramref name="gapScore"/>(k).
    /// Model, ties and coordinates as <see cref="LocalAlignAffine(string, string, ScoringMatrix?)"/> (maximum at the first row-major cell,
    /// traceback stops when the predecessor is the zero floor; no positive-scoring pair → the tagged-Local empty result); see
    /// <see cref="GlobalAlignGeneralGap(string, string, Func{int, double}, ScoringMatrix?)"/> for the gap model, complexity and fast paths.
    /// Oracle: Biopython 1.88 <c>PairwiseAligner(mode='local', gap_score=callable)</c>.
    /// </summary>
    public static AlignmentResult LocalAlignGeneralGap(
        string sequence1, string sequence2, Func<int, double> gapScore, ScoringMatrix? scoring = null) =>
        AlignGeneralGap(sequence1, sequence2, gapScore, ScalarSubstitution(scoring), AlignmentType.Local, EndGaps.None);

    /// <summary><see cref="LocalAlignGeneralGap(string, string, Func{int, double}, ScoringMatrix?)"/> for DNA sequences.</summary>
    public static AlignmentResult LocalAlignGeneralGap(
        DnaSequence sequence1, DnaSequence sequence2, Func<int, double> gapScore, ScoringMatrix? scoring = null)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        return AlignGeneralGap(sequence1.Sequence, sequence2.Sequence, gapScore, ScalarSubstitution(scoring), AlignmentType.Local, EndGaps.None);
    }

    /// <summary>General-gap local alignment scored with a <see cref="SubstitutionMatrix"/>.</summary>
    public static AlignmentResult LocalAlignGeneralGap(
        string sequence1, string sequence2, Func<int, double> gapScore, SubstitutionMatrix matrix) =>
        AlignGeneralGap(sequence1, sequence2, gapScore, MatrixSubstitution(matrix), AlignmentType.Local, EndGaps.None);

    /// <summary><see cref="LocalAlignGeneralGap(string, string, Func{int, double}, SubstitutionMatrix)"/> for protein sequences.</summary>
    public static AlignmentResult LocalAlignGeneralGap(
        ProteinSequence sequence1, ProteinSequence sequence2, Func<int, double> gapScore, SubstitutionMatrix matrix)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        return AlignGeneralGap(sequence1.Sequence, sequence2.Sequence, gapScore, MatrixSubstitution(matrix), AlignmentType.Local, EndGaps.None);
    }

    /// <summary><see cref="LocalAlignGeneralGap(string, string, Func{int, double}, SubstitutionMatrix)"/> for DNA sequences.</summary>
    public static AlignmentResult LocalAlignGeneralGap(
        DnaSequence sequence1, DnaSequence sequence2, Func<int, double> gapScore, SubstitutionMatrix matrix)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        return AlignGeneralGap(sequence1.Sequence, sequence2.Sequence, gapScore, MatrixSubstitution(matrix), AlignmentType.Local, EndGaps.None);
    }

    // ----- semi-global / ends-free -----

    /// <summary>
    /// Ends-free (semi-global) alignment with Waterman-Smith-Beyer gaps: a gap of length k scores <paramref name="gapScore"/>(k) unless it is a
    /// free end gap (<paramref name="endGaps"/>, default <see cref="EndGaps.Fitting"/> = sequence 1 fitted into sequence 2). Policy, ties,
    /// coordinates: <see cref="SemiGlobalAlign(string, string, ScoringMatrix?, EndGaps)"/>; gap model:
    /// <see cref="GlobalAlignGeneralGap(string, string, Func{int, double}, ScoringMatrix?)"/>. (Biopython 1.88 ignores the end-gap properties when
    /// <c>gap_score</c> is a callable, so the end-gap handling is locked against an exhaustive enumeration and against the affine aligner on
    /// affine gap functions instead.)
    /// </summary>
    public static AlignmentResult SemiGlobalAlignGeneralGap(
        string sequence1, string sequence2, Func<int, double> gapScore, ScoringMatrix? scoring = null, EndGaps endGaps = EndGaps.Fitting) =>
        AlignGeneralGap(sequence1, sequence2, gapScore, ScalarSubstitution(scoring), AlignmentType.SemiGlobal, endGaps);

    /// <summary><see cref="SemiGlobalAlignGeneralGap(string, string, Func{int, double}, ScoringMatrix?, EndGaps)"/> for DNA sequences.</summary>
    public static AlignmentResult SemiGlobalAlignGeneralGap(
        DnaSequence sequence1, DnaSequence sequence2, Func<int, double> gapScore, ScoringMatrix? scoring = null, EndGaps endGaps = EndGaps.Fitting)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        return AlignGeneralGap(sequence1.Sequence, sequence2.Sequence, gapScore, ScalarSubstitution(scoring), AlignmentType.SemiGlobal, endGaps);
    }

    /// <summary>General-gap ends-free alignment scored with a <see cref="SubstitutionMatrix"/>.</summary>
    public static AlignmentResult SemiGlobalAlignGeneralGap(
        string sequence1, string sequence2, Func<int, double> gapScore, SubstitutionMatrix matrix, EndGaps endGaps = EndGaps.Fitting) =>
        AlignGeneralGap(sequence1, sequence2, gapScore, MatrixSubstitution(matrix), AlignmentType.SemiGlobal, endGaps);

    /// <summary><see cref="SemiGlobalAlignGeneralGap(string, string, Func{int, double}, SubstitutionMatrix, EndGaps)"/> for protein sequences.</summary>
    public static AlignmentResult SemiGlobalAlignGeneralGap(
        ProteinSequence sequence1, ProteinSequence sequence2, Func<int, double> gapScore, SubstitutionMatrix matrix, EndGaps endGaps = EndGaps.Fitting)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        return AlignGeneralGap(sequence1.Sequence, sequence2.Sequence, gapScore, MatrixSubstitution(matrix), AlignmentType.SemiGlobal, endGaps);
    }

    /// <summary><see cref="SemiGlobalAlignGeneralGap(string, string, Func{int, double}, SubstitutionMatrix, EndGaps)"/> for DNA sequences.</summary>
    public static AlignmentResult SemiGlobalAlignGeneralGap(
        DnaSequence sequence1, DnaSequence sequence2, Func<int, double> gapScore, SubstitutionMatrix matrix, EndGaps endGaps = EndGaps.Fitting)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        return AlignGeneralGap(sequence1.Sequence, sequence2.Sequence, gapScore, MatrixSubstitution(matrix), AlignmentType.SemiGlobal, endGaps);
    }

    // ----- shared driver -----

    private static AlignmentResult AlignGeneralGap(
        string sequence1, string sequence2, Func<int, double> gapScore, GapSubstitution substitution, AlignmentType type, EndGaps ends,
        bool allowFastPath = true)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        ArgumentNullException.ThrowIfNull(gapScore);

        string a = sequence1.ToUpperInvariant();
        string b = sequence2.ToUpperInvariant();
        if (type != AlignmentType.SemiGlobal) ends = EndGaps.None;
        if (substitution.Matrix is { } matrix)
        {
            PairScoring.ValidateResidues(matrix, a, "sequence1");
            PairScoring.ValidateResidues(matrix, b, "sequence2");
        }

        int maxLen = Math.Max(a.Length, b.Length);
        var w = new double[maxLen + 1];
        for (int k = 1; k <= maxLen; k++)
        {
            double v = gapScore(k);
            if (double.IsNaN(v) || double.IsPositiveInfinity(v))
                throw new ArgumentException($"gapScore({k}) = {v}: a gap score must be a real number or -infinity.", nameof(gapScore));
            w[k] = v;
        }

        if (allowFastPath && TryDelegateToAffine(a, b, w, substitution, type, ends, out var delegated))
            return delegated;

        return GeneralGapCore(a, b, w, substitution.Score, type, ends);
    }

    /// <summary>Internal entry that always runs the general program (the fast path is what the tests compare it with).</summary>
    internal static AlignmentResult AlignGeneralGapGeneric(
        string sequence1, string sequence2, Func<int, double> gapScore, ScoringMatrix scoring, AlignmentType type, EndGaps ends) =>
        AlignGeneralGap(sequence1, sequence2, gapScore, ScalarSubstitution(scoring), type, ends, allowFastPath: false);

    /// <summary>Same as <see cref="AlignGeneralGapGeneric"/> for a substitution matrix.</summary>
    internal static AlignmentResult AlignGeneralGapGeneric(
        string sequence1, string sequence2, Func<int, double> gapScore, SubstitutionMatrix matrix, AlignmentType type, EndGaps ends) =>
        AlignGeneralGap(sequence1, sequence2, gapScore, MatrixSubstitution(matrix), type, ends, allowFastPath: false);

    private const double GapFitTolerance = 1e-9;

    private static bool TryDelegateToAffine(
        string a, string b, double[] w, GapSubstitution substitution, AlignmentType type, EndGaps ends, out AlignmentResult result)
    {
        result = null!;
        int maxLen = w.Length - 1;
        if (maxLen < 1)
            return false;

        // w(k) = w(1) + (k-1)·e  <=>  gap of length k = open + k·e with open = w(1) - e
        double w1 = w[1];
        double e = maxLen >= 2 ? w[2] - w1 : 0;
        if (double.IsInfinity(w1) || double.IsInfinity(e))
            return false;
        for (int k = 1; k <= maxLen; k++)
        {
            double expected = w1 + (k - 1) * e;
            if (double.IsInfinity(w[k]) || Math.Abs(w[k] - expected) > GapFitTolerance * Math.Max(1.0, Math.Abs(expected)))
                return false;
        }

        // representable at the 1e-6 resolution of the integer dynamic programs, and a penalty (open <= 0, extend <= 0)
        double openRaw = w1 - e;
        double open = Math.Round(openRaw, 6);
        e = Math.Round(e, 6);
        if (Math.Abs(openRaw - open) > GapFitTolerance || (maxLen >= 2 && Math.Abs(w[2] - w1 - e) > GapFitTolerance) || e > 0 || open > 0)
            return false;

        if (substitution.Scalar is { } sc)
        {
            if (open != Math.Round(open) || e != Math.Round(e))
                return false; // the scalar scoring model is integral
            var scoring = new ScoringMatrix(sc.Match, sc.Mismatch, (int)open, (int)e);
            var ps = new PairScoring(scoring);
            AlignmentResult r = RunDelegated(a, b, ps, open == 0, type, ends);
            result = r with { ExactScore = r.Score };
            return true;
        }

        var matrix = substitution.Matrix!;
        PairScoring pm;
        try
        {
            pm = PairScoring.Create(new MatrixScoring(matrix, open, e), a, b);
        }
        catch (ArgumentException)
        {
            return false; // integer range / resolution: the general program handles it in doubles
        }
        AlignmentResult rm = RunDelegated(a, b, pm, open == 0, type, ends);
        double exact = (double)rm.Score / pm.Scale;
        result = rm with { Score = (int)Math.Round(exact, MidpointRounding.AwayFromZero), ExactScore = exact };
        return true;
    }

    private static AlignmentResult RunDelegated(string a, string b, PairScoring ps, bool linear, AlignmentType type, EndGaps ends)
    {
        if (type == AlignmentType.Local)
            return linear ? LocalAlignCore(a, b, ps) : AffineAlignCore(a, b, ps, AlignmentType.Local);
        if (linear)
            return LinearAlignCore(a, b, ps, ends);
        var affineType = type == AlignmentType.Global ? AlignmentType.Global : AlignmentType.SemiGlobal;
        return AffineAlignCore(a, b, ps, affineType, ends);
    }

    // back-pointer codes of the M state: predecessor state, or a fresh local start
    private const byte FromM = 0, FromX = 1, FromY = 2, FromStart = 3;

    /// <summary>
    /// The Waterman-Smith-Beyer dynamic program (three states, every gap run charged w(k) once), see the region comment. <paramref name="w"/>[k]
    /// = gap score of length k (index 0 unused).
    /// </summary>
    private static AlignmentResult GeneralGapCore(
        string seq1, string seq2, double[] w, Func<char, char, double> sub, AlignmentType type, EndGaps ends)
    {
        bool local = type == AlignmentType.Local;
        if (local) ends = EndGaps.None;
        bool freeStart1 = ends.HasFlag(EndGaps.FreeStart1);
        bool freeEnd1 = ends.HasFlag(EndGaps.FreeEnd1);
        bool freeStart2 = ends.HasFlag(EndGaps.FreeStart2);
        bool freeEnd2 = ends.HasFlag(EndGaps.FreeEnd2);
        int m = seq1.Length;
        int n = seq2.Length;
        const double NegInf = double.NegativeInfinity;

        var mm = new double[m + 1, n + 1];
        var xx = new double[m + 1, n + 1];
        var yy = new double[m + 1, n + 1];
        var pm = new byte[m + 1, n + 1];
        var px = new int[m + 1, n + 1]; // 2·k + (1 if the run followed a Y state, 0 if it followed M)
        var py = new int[m + 1, n + 1]; // 2·k + (1 if the run followed an X state, 0 if it followed M)
        for (int i = 0; i <= m; i++)
            for (int j = 0; j <= n; j++)
                mm[i, j] = xx[i, j] = yy[i, j] = NegInf;

        mm[0, 0] = local ? NegInf : 0;
        for (int i = 1; i <= m; i++)
        {
            if (freeStart1 && !local) mm[i, 0] = 0;
            else if (!local) { xx[i, 0] = w[i]; px[i, 0] = 2 * i; } // one leading run of i positions
        }
        for (int j = 1; j <= n; j++)
        {
            if (freeStart2 && !local) mm[0, j] = 0;
            else if (!local) { yy[0, j] = w[j]; py[0, j] = 2 * j; }
        }

        double bestLocal = 0;
        int bestI = 0, bestJ = 0;
        for (int i = 1; i <= m; i++)
        {
            for (int j = 1; j <= n; j++)
            {
                double s = sub(seq1[i - 1], seq2[j - 1]);

                // M: predecessor with the best value, ties M > X > Y; local: a non-positive predecessor is replaced by a fresh start
                double pv = mm[i - 1, j - 1];
                byte from = FromM;
                if (xx[i - 1, j - 1] > pv) { pv = xx[i - 1, j - 1]; from = FromX; }
                if (yy[i - 1, j - 1] > pv) { pv = yy[i - 1, j - 1]; from = FromY; }
                if (local && pv <= 0) { pv = 0; from = FromStart; }
                mm[i, j] = pv + s;
                pm[i, j] = from;

                // X: the run a_{i-k+1..i} against a gap, after M or Y at (i-k, j); smallest k first, M before Y
                double bestX = NegInf;
                int codeX = 0;
                for (int k = 1; k <= i; k++)
                {
                    double wk = w[k];
                    double c = mm[i - k, j] + wk;
                    if (c > bestX) { bestX = c; codeX = 2 * k; }
                    c = yy[i - k, j] + wk;
                    if (c > bestX) { bestX = c; codeX = 2 * k + 1; }
                }
                xx[i, j] = bestX;
                px[i, j] = codeX;

                // Y: the run b_{j-k+1..j} against a gap, after M or X at (i, j-k)
                double bestY = NegInf;
                int codeY = 0;
                for (int k = 1; k <= j; k++)
                {
                    double wk = w[k];
                    double c = mm[i, j - k] + wk;
                    if (c > bestY) { bestY = c; codeY = 2 * k; }
                    c = xx[i, j - k] + wk;
                    if (c > bestY) { bestY = c; codeY = 2 * k + 1; }
                }
                yy[i, j] = bestY;
                py[i, j] = codeY;

                if (local && mm[i, j] > bestLocal)
                {
                    bestLocal = mm[i, j];
                    bestI = i;
                    bestJ = j;
                }
            }
        }

        int ci, cj;
        byte state;
        double best;
        if (local)
        {
            if (bestLocal <= 0)
                return EmptyLocalResult();
            best = bestLocal;
            ci = bestI;
            cj = bestJ;
            state = FromM;
        }
        else
        {
            ci = m;
            cj = n;
            double bestValue = NegInf;
            bool any = false;
            if (freeEnd2)
                for (int j = 0; j <= n; j++)
                {
                    double v = Math.Max(mm[m, j], Math.Max(xx[m, j], yy[m, j]));
                    if (!any || v > bestValue) { bestValue = v; ci = m; cj = j; any = true; }
                }
            if (freeEnd1)
                for (int i = 0; i <= m; i++)
                {
                    double v = Math.Max(mm[i, n], Math.Max(xx[i, n], yy[i, n]));
                    if (!any || v > bestValue) { bestValue = v; ci = i; cj = n; any = true; }
                }
            best = Math.Max(mm[ci, cj], Math.Max(xx[ci, cj], yy[ci, cj]));
            if (double.IsNegativeInfinity(best))
                throw new ArgumentException("The gap score function forbids every alignment of these sequences.");
            state = FromY;
            if (mm[ci, cj] == best) state = FromM;
            else if (xx[ci, cj] == best) state = FromX;
        }

        int endI = ci, endJ = cj;
        var chars1 = new List<char>(m + n);
        var chars2 = new List<char>(m + n);
        if (!local)
        {
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

        while (!(ci == 0 && (cj == 0 || freeStart2)) && !(cj == 0 && freeStart1))
        {
            if (state == FromM)
            {
                byte from = pm[ci, cj];
                chars1.Add(seq1[ci - 1]);
                chars2.Add(seq2[cj - 1]);
                ci--; cj--;
                if (from == FromStart) break; // local alignment starts here
                state = from;
            }
            else if (state == FromX)
            {
                int code = px[ci, cj];
                int k = code >> 1;
                for (int t = 0; t < k; t++)
                {
                    chars1.Add(seq1[ci - 1 - t]);
                    chars2.Add('-');
                }
                ci -= k;
                state = (code & 1) == 0 ? FromM : FromY;
            }
            else
            {
                int code = py[ci, cj];
                int k = code >> 1;
                for (int t = 0; t < k; t++)
                {
                    chars1.Add('-');
                    chars2.Add(seq2[cj - 1 - t]);
                }
                cj -= k;
                state = (code & 1) == 0 ? FromM : FromX;
            }
        }

        if (!local)
            AppendFreeLeadingOverhang(chars1, chars2, seq1, seq2, ci, cj);

        chars1.Reverse();
        chars2.Reverse();

        return new AlignmentResult(
            AlignedSequence1: new string(chars1.ToArray()),
            AlignedSequence2: new string(chars2.ToArray()),
            Score: (int)Math.Round(best, MidpointRounding.AwayFromZero),
            AlignmentType: ResultType(local, ends),
            StartPosition1: ci,
            StartPosition2: cj,
            EndPosition1: endI - 1,
            EndPosition2: endJ - 1)
        {
            ExactScore = best,
        };
    }

    #endregion
}
