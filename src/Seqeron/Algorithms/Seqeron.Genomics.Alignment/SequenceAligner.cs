using System.Text;

namespace Seqeron.Genomics.Alignment;

/// <summary>
/// Performs pairwise sequence alignment using dynamic programming algorithms.
/// Supports global (Needleman-Wunsch), local (Smith-Waterman), and semi-global alignment.
/// </summary>
public static class SequenceAligner
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
    /// </remarks>
    /// <param name="sequence1">First DNA sequence.</param>
    /// <param name="sequence2">Second DNA sequence.</param>
    /// <param name="scoring">Scoring matrix (default: SimpleDna).</param>
    /// <returns>Alignment result with aligned sequences and score.</returns>
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
    public static AlignmentResult GlobalAlign(
        string sequence1,
        string sequence2,
        ScoringMatrix? scoring = null)
    {
        if (string.IsNullOrEmpty(sequence1) || string.IsNullOrEmpty(sequence2))
            return AlignmentResult.Empty;

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
        if (string.IsNullOrEmpty(sequence1) || string.IsNullOrEmpty(sequence2))
            return AlignmentResult.Empty;

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
        ScoringMatrix scoring,
        CancellationToken cancellationToken = default,
        IProgress<double>? progress = null) =>
        LinearAlignCore(seq1, seq2, scoring, fitting: false, cancellationToken, progress);

    /// <summary>
    /// Linear-gap Needleman-Wunsch dynamic program shared by <c>GlobalAlign</c> and
    /// <c>SemiGlobalAlign</c>. d = <see cref="ScoringMatrix.GapExtend"/> per gap position.
    /// <list type="bullet">
    /// <item><b>Global</b>: F(i,0) = d·i, F(0,j) = d·j; traceback from (m, n).</item>
    /// <item><b>Fitting</b> (semi-global, sequence 1 fitted into sequence 2): F(0,j) = 0 (leading
    /// sequence-2 residues free), F(i,0) = d·i; traceback from the first (smallest j) maximum of the
    /// last row, max_j F(m,j) (trailing sequence-2 residues free).</item>
    /// </list>
    /// Traceback ties are broken diagonal &gt; up &gt; left.
    /// </summary>
    private static AlignmentResult LinearAlignCore(
        string seq1,
        string seq2,
        ScoringMatrix scoring,
        bool fitting,
        CancellationToken cancellationToken = default,
        IProgress<double>? progress = null)
    {
        int m = seq1.Length;
        int n = seq2.Length;

        cancellationToken.ThrowIfCancellationRequested();

        var score = new int[m + 1, n + 1];

        // Initialize first row and column with linear gap penalty d = GapExtend.
        // Standard Needleman-Wunsch: F(i,0) = d*i, F(0,j) = d*j
        // Source: https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm
        // Fitting alignment: F(0,j) = 0 — a free prefix of seq2 may precede the alignment.
        for (int i = 0; i <= m; i++)
            score[i, 0] = i * scoring.GapExtend;
        for (int j = 0; j <= n; j++)
            score[0, j] = fitting ? 0 : j * scoring.GapExtend;

        // F(i,j) = max(F(i-1,j-1)+S(a_i,b_j), F(i-1,j)+d, F(i,j-1)+d)
        for (int i = 1; i <= m; i++)
        {
            if (i % 100 == 0)
            {
                cancellationToken.ThrowIfCancellationRequested();
                progress?.Report((double)i / m * 0.5); // First half is matrix fill
            }

            char c1 = seq1[i - 1];
            for (int j = 1; j <= n; j++)
            {
                int matchScore = c1 == seq2[j - 1] ? scoring.Match : scoring.Mismatch;

                int diag = score[i - 1, j - 1] + matchScore;
                int up = score[i - 1, j] + scoring.GapExtend;
                int left = score[i, j - 1] + scoring.GapExtend;

                score[i, j] = Math.Max(diag, Math.Max(up, left));
            }
        }

        cancellationToken.ThrowIfCancellationRequested();
        progress?.Report(0.75);

        int endJ = n;
        if (fitting)
        {
            // Free trailing seq2 residues: end in the first column holding max_j F(m,j).
            endJ = 0;
            for (int j = 1; j <= n; j++)
                if (score[m, j] > score[m, endJ])
                    endJ = j;
        }

        var result = Traceback(seq1, seq2, score, m, endJ, scoring,
            fitting ? AlignmentType.SemiGlobal : AlignmentType.Global);

        progress?.Report(1.0);
        return result;
    }

    /// <summary>
    /// Exact linear-space global alignment (Hirschberg 1975, Commun ACM 18:341-343; Myers &amp; Miller
    /// 1988) under the same linear-gap model as <see cref="GlobalAlign(DnaSequence, DnaSequence, ScoringMatrix?)"/>:
    /// the score equals the Needleman-Wunsch optimum, with O(m+n) working memory per recursion
    /// level instead of O(m·n). The recursion splits <paramref name="seq1"/> in the middle, finds the
    /// optimal split column of <paramref name="seq2"/> from a forward and a backward last-row score
    /// pass, and recurses; sub-problems of at most <paramref name="maxFullMatrixCells"/> cells are
    /// solved by the canonical <c>GlobalAlignCore</c>. Tie-breaking of co-optimal alignments may differ
    /// from <c>GlobalAlign</c> (only the score is guaranteed identical).
    /// </summary>
    internal static AlignmentResult GlobalAlignLinearSpace(
        string seq1, string seq2, ScoringMatrix scoring, long maxFullMatrixCells = 25_000_000)
    {
        ArgumentNullException.ThrowIfNull(seq1);
        ArgumentNullException.ThrowIfNull(seq2);
        ArgumentNullException.ThrowIfNull(scoring);

        var out1 = new StringBuilder(seq1.Length + seq2.Length);
        var out2 = new StringBuilder(seq1.Length + seq2.Length);
        int score = 0;
        HirschbergAlign(seq1, 0, seq1.Length, seq2, 0, seq2.Length, scoring, maxFullMatrixCells, out1, out2, ref score);

        return new AlignmentResult(out1.ToString(), out2.ToString(), score, AlignmentType.Global,
            0, 0, seq1.Length - 1, seq2.Length - 1);
    }

    private static void HirschbergAlign(
        string a, int aStart, int aEnd, string b, int bStart, int bEnd,
        ScoringMatrix scoring, long maxFullMatrixCells, StringBuilder out1, StringBuilder out2, ref int score)
    {
        int m = aEnd - aStart;
        int n = bEnd - bStart;

        if (m == 0 || n == 0)
        {
            // One side empty: the other is an all-gap run.
            if (m > 0)
            {
                out1.Append(a, aStart, m);
                out2.Append('-', m);
            }
            else if (n > 0)
            {
                out1.Append('-', n);
                out2.Append(b, bStart, n);
            }
            score += (m + n) * scoring.GapExtend;
            return;
        }

        if (m == 1 || (long)(m + 1) * (n + 1) <= maxFullMatrixCells)
        {
            var r = GlobalAlignCore(a.Substring(aStart, m), b.Substring(bStart, n), scoring);
            out1.Append(r.AlignedSequence1);
            out2.Append(r.AlignedSequence2);
            score += r.Score;
            return;
        }

        int mid = m / 2;
        int[] forward = LinearLastRow(a, aStart, aStart + mid, b, bStart, bEnd, scoring, reverse: false);
        int[] backward = LinearLastRow(a, aStart + mid, aEnd, b, bStart, bEnd, scoring, reverse: true);

        int split = 0;
        long bestSplit = long.MinValue;
        for (int j = 0; j <= n; j++)
        {
            long total = (long)forward[j] + backward[n - j];
            if (total > bestSplit)
            {
                bestSplit = total;
                split = j;
            }
        }

        HirschbergAlign(a, aStart, aStart + mid, b, bStart, bStart + split, scoring, maxFullMatrixCells, out1, out2, ref score);
        HirschbergAlign(a, aStart + mid, aEnd, b, bStart + split, bEnd, scoring, maxFullMatrixCells, out1, out2, ref score);
    }

    /// <summary>
    /// Last row of the linear-gap NW matrix of a[aStart..aEnd) against b[bStart..bEnd): entry j is the
    /// optimal score of aligning all of a with the first j characters of b (or, with
    /// <paramref name="reverse"/>, with the last j characters of b — both strings read backwards).
    /// </summary>
    private static int[] LinearLastRow(
        string a, int aStart, int aEnd, string b, int bStart, int bEnd, ScoringMatrix scoring, bool reverse)
    {
        int m = aEnd - aStart;
        int n = bEnd - bStart;
        var prev = new int[n + 1];
        var cur = new int[n + 1];

        for (int j = 0; j <= n; j++)
            prev[j] = j * scoring.GapExtend;

        for (int i = 1; i <= m; i++)
        {
            char ca = reverse ? a[aEnd - i] : a[aStart + i - 1];
            cur[0] = i * scoring.GapExtend;
            for (int j = 1; j <= n; j++)
            {
                char cb = reverse ? b[bEnd - j] : b[bStart + j - 1];
                int diag = prev[j - 1] + (ca == cb ? scoring.Match : scoring.Mismatch);
                int up = prev[j] + scoring.GapExtend;
                int left = cur[j - 1] + scoring.GapExtend;
                cur[j] = Math.Max(diag, Math.Max(up, left));
            }
            (prev, cur) = (cur, prev);
        }

        return prev;
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

    private static AlignmentResult GlobalAlignAffineCore(string seq1, string seq2, ScoringMatrix scoring)
    {
        if (seq1.Length == 0 && seq2.Length == 0)
            return AlignmentResult.Empty;

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
    private static AlignmentResult AffineAlignCore(string seq1, string seq2, ScoringMatrix scoring, AlignmentType type)
    {
        bool local = type == AlignmentType.Local;
        bool fitting = type == AlignmentType.SemiGlobal;
        int m = seq1.Length;
        int n = seq2.Length;

        int o = scoring.GapOpen;
        int e = scoring.GapExtend;
        int oe = o + e;
        int startFloor = local ? 0 : AffineNegInf;

        var mm = new int[m + 1, n + 1];
        var xx = new int[m + 1, n + 1];
        var yy = new int[m + 1, n + 1];

        mm[0, 0] = local ? AffineNegInf : 0;
        xx[0, 0] = AffineNegInf;
        yy[0, 0] = AffineNegInf;
        for (int i = 1; i <= m; i++)
        {
            mm[i, 0] = AffineNegInf;
            xx[i, 0] = local ? AffineNegInf : o + i * e;
            yy[i, 0] = AffineNegInf;
        }
        for (int j = 1; j <= n; j++)
        {
            mm[0, j] = fitting ? 0 : AffineNegInf;
            xx[0, j] = AffineNegInf;
            yy[0, j] = local || fitting ? AffineNegInf : o + j * e;
        }

        int bestLocal = 0, bestI = 0, bestJ = 0;
        for (int i = 1; i <= m; i++)
        {
            char c1 = seq1[i - 1];
            for (int j = 1; j <= n; j++)
            {
                int s = c1 == seq2[j - 1] ? scoring.Match : scoring.Mismatch;

                mm[i, j] = Math.Max(startFloor, Max3(mm[i - 1, j - 1], xx[i - 1, j - 1], yy[i - 1, j - 1])) + s;
                xx[i, j] = Max3(mm[i - 1, j] + oe, xx[i - 1, j] + e, yy[i - 1, j] + oe);
                yy[i, j] = Max3(mm[i, j - 1] + oe, yy[i, j - 1] + e, xx[i, j - 1] + oe);

                // Clamp unreachable states so repeated penalties never drift towards overflow.
                if (xx[i, j] < AffineNegInf) xx[i, j] = AffineNegInf;
                if (yy[i, j] < AffineNegInf) yy[i, j] = AffineNegInf;
                if (mm[i, j] < AffineNegInf) mm[i, j] = AffineNegInf;

                if (local && mm[i, j] > bestLocal)
                {
                    bestLocal = mm[i, j];
                    bestI = i;
                    bestJ = j;
                }
            }
        }

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
            // Global: end at (m, n). Fitting: end in the first column of the last row holding the
            // maximum (free trailing seq2 residues). Optimal end state ties: M > X > Y.
            ci = m;
            cj = n;
            if (fitting)
            {
                cj = 0;
                for (int j = 1; j <= n; j++)
                    if (Max3(mm[m, j], xx[m, j], yy[m, j]) > Max3(mm[m, cj], xx[m, cj], yy[m, cj]))
                        cj = j;
            }
            best = Max3(mm[ci, cj], xx[ci, cj], yy[ci, cj]);
            state = PickState(best, mm[ci, cj], StateMatch, xx[ci, cj], StateGapInSeq2, StateGapInSeq1);
        }

        int endJ = cj;
        var chars1 = new List<char>(m + n);
        var chars2 = new List<char>(m + n);
        if (fitting)
        {
            // Free trailing seq2 suffix (emitted reversed; the whole buffer is reversed below).
            for (int k = n; k > endJ; k--)
            {
                chars1.Add('-');
                chars2.Add(seq2[k - 1]);
            }
        }

        while (fitting ? ci > 0 : ci > 0 || cj > 0)
        {
            if (state == StateMatch)
            {
                int prev = mm[ci, cj] - (seq1[ci - 1] == seq2[cj - 1] ? scoring.Match : scoring.Mismatch);
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

        if (fitting)
        {
            // Free leading seq2 prefix: the alignment started at M(0, cj) = 0.
            for (int k = cj; k > 0; k--)
            {
                chars1.Add('-');
                chars2.Add(seq2[k - 1]);
            }
        }

        chars1.Reverse();
        chars2.Reverse();

        if (fitting)
            return new AlignmentResult(
                AlignedSequence1: new string(chars1.ToArray()),
                AlignedSequence2: new string(chars2.ToArray()),
                Score: best,
                AlignmentType: AlignmentType.SemiGlobal,
                StartPosition1: 0,
                StartPosition2: 0,
                EndPosition1: m - 1,
                EndPosition2: n - 1);

        return local
            ? new AlignmentResult(
                AlignedSequence1: new string(chars1.ToArray()),
                AlignedSequence2: new string(chars2.ToArray()),
                Score: best,
                AlignmentType: AlignmentType.Local,
                StartPosition1: ci,
                StartPosition2: cj,
                EndPosition1: bestI - 1,
                EndPosition2: bestJ - 1)
            : new AlignmentResult(
                AlignedSequence1: new string(chars1.ToArray()),
                AlignedSequence2: new string(chars2.ToArray()),
                Score: best,
                AlignmentType: AlignmentType.Global,
                StartPosition1: 0,
                StartPosition2: 0,
                EndPosition1: m - 1,
                EndPosition2: n - 1);
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
    /// has Score 0, empty strings and all positions −1.
    /// </remarks>
    /// <param name="sequence1">First DNA sequence.</param>
    /// <param name="sequence2">Second DNA sequence.</param>
    /// <param name="scoring">Scoring matrix (default: SimpleDna).</param>
    /// <returns>Alignment result with aligned subsequences and score.</returns>
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
    public static AlignmentResult LocalAlign(
        string sequence1,
        string sequence2,
        ScoringMatrix? scoring = null)
    {
        if (string.IsNullOrEmpty(sequence1) || string.IsNullOrEmpty(sequence2))
            return AlignmentResult.Empty;

        return LocalAlignCore(
            sequence1.ToUpperInvariant(),
            sequence2.ToUpperInvariant(),
            scoring ?? SimpleDna);
    }

    private static AlignmentResult LocalAlignCore(string seq1, string seq2, ScoringMatrix scoring)
    {
        int m = seq1.Length;
        int n = seq2.Length;

        var score = new int[m + 1, n + 1];
        int maxScore = 0;
        int maxI = 0, maxJ = 0;

        // Fill the matrix (with zero floor for local alignment)
        for (int i = 1; i <= m; i++)
        {
            for (int j = 1; j <= n; j++)
            {
                int matchScore = seq1[i - 1] == seq2[j - 1] ? scoring.Match : scoring.Mismatch;

                int diag = score[i - 1, j - 1] + matchScore;
                int up = score[i - 1, j] + scoring.GapExtend;
                int left = score[i, j - 1] + scoring.GapExtend;

                score[i, j] = Math.Max(0, Math.Max(diag, Math.Max(up, left)));

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

    private static AlignmentResult TracebackLocal(
        string seq1, string seq2, int[,] score, int endI, int endJ, ScoringMatrix scoring)
    {
        var chars1 = new List<char>();
        var chars2 = new List<char>();

        int i = endI, j = endJ;
        int startI = endI, startJ = endJ;

        while (i > 0 && j > 0 && score[i, j] > 0)
        {
            startI = i;
            startJ = j;

            int matchScore = seq1[i - 1] == seq2[j - 1] ? scoring.Match : scoring.Mismatch;

            if (score[i, j] == score[i - 1, j - 1] + matchScore)
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
    /// linear-gap <see cref="LocalAlign(DnaSequence, DnaSequence, ScoringMatrix?)"/>. Gap scores are
    /// expected to be non-positive (GapOpen ≤ 0, GapExtend ≤ 0), as in every standard scheme.
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
    /// prefix/suffix is shown against '-'), so Start/End positions are 0 and length − 1.
    /// </para>
    /// </remarks>
    public static AlignmentResult SemiGlobalAlign(
        DnaSequence sequence1,
        DnaSequence sequence2,
        ScoringMatrix? scoring = null)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);

        return LinearAlignCore(sequence1.Sequence, sequence2.Sequence, scoring ?? SimpleDna, fitting: true);
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
    /// traceback ties M &gt; X &gt; Y). The aligned strings span both full inputs, so Start/End positions
    /// are 0 and length − 1.
    /// </para>
    /// <para>Sources: Gotoh O (1982) J Mol Biol 162:705-708; Flouri et al. (2015) bioRxiv 10.1101/031500;
    /// cross-checked against Biopython 1.88 PairwiseAligner and parasail 1.3.4 sg_dx_trace.</para>
    /// </remarks>
    public static AlignmentResult SemiGlobalAlignAffine(
        DnaSequence sequence1,
        DnaSequence sequence2,
        ScoringMatrix? scoring = null)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);

        return AffineAlignCore(sequence1.Sequence, sequence2.Sequence, scoring ?? SimpleDna, AlignmentType.SemiGlobal);
    }

    /// <summary>
    /// Affine-gap semi-global (fitting) alignment on raw sequence strings (uppercased internally).
    /// See <see cref="SemiGlobalAlignAffine(DnaSequence, DnaSequence, ScoringMatrix?)"/> for the model.
    /// </summary>
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
            AlignmentType.SemiGlobal);
    }

    #endregion

    #region Traceback

    private static AlignmentResult Traceback(
        string seq1, string seq2, int[,] score, int i, int j,
        ScoringMatrix scoring, AlignmentType alignType)
    {
        // Preserve the traceback origin for correct score retrieval.
        // For semi-global (fitting): score[m, maxJ] is the optimal fitting score.
        // For global (NW): score[m, n] is the bottom-right cell.
        int endI = i, endJ = j;

        var chars1 = new List<char>();
        var chars2 = new List<char>();

        // Add trailing gaps for semi-global (appended first, will be at end after reverse).
        // chars2 is reversed in its entirety below, so the unmatched trailing reference
        // suffix seq2[j .. n-1] must be appended in REVERSE order here; otherwise the
        // suffix would itself be reversed in the final output and AlignedSequence2 would
        // no longer reproduce the reference.
        if (alignType == AlignmentType.SemiGlobal)
        {
            for (int k = seq2.Length; k > j; k--)
            {
                chars1.Add('-');
                chars2.Add(seq2[k - 1]);
            }
        }

        while (i > 0 || j > 0)
        {
            if (i > 0 && j > 0)
            {
                int matchScore = seq1[i - 1] == seq2[j - 1] ? scoring.Match : scoring.Mismatch;

                if (score[i, j] == score[i - 1, j - 1] + matchScore)
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

        chars1.Reverse();
        chars2.Reverse();

        string aligned1Str = new string(chars1.ToArray());
        string aligned2Str = new string(chars2.ToArray());

        return new AlignmentResult(
            AlignedSequence1: aligned1Str,
            AlignedSequence2: aligned2Str,
            Score: score[endI, endJ],
            AlignmentType: alignType,
            StartPosition1: 0,
            StartPosition2: 0,
            EndPosition1: seq1.Length - 1,
            EndPosition2: seq2.Length - 1);
    }

    #endregion

    #region Alignment Statistics

    // EMBOSS srspair markup legend (Rice, Longden & Bleasby 2000; EMBOSS AlignFormats):
    //   '|' = identical column, ':' = similar column (positive substitution score, not
    //   identical), ' ' = gap or non-positive (mismatch) column.
    //   Source: https://emboss.sourceforge.net/docs/themes/AlignFormats.html
    private const char IdentityMark = '|';
    private const char SimilarityMark = ':';
    private const char GapOrMismatchMark = ' ';

    /// <summary>
    /// Calculates alignment statistics (Identity, Similarity, Gaps) from an alignment
    /// result, following the EMBOSS needle/water convention.
    /// </summary>
    /// <param name="alignment">The pairwise alignment to summarise.</param>
    /// <param name="scoring">
    /// Scoring model used to decide whether a non-identical aligned pair counts as
    /// "similar". A column is similar when its substitution score is positive
    /// (Rice et al. 2000). For the DNA models used by this class (Match &gt; 0,
    /// Mismatch &lt; 0) no non-identical pair scores positively, so Similarity equals
    /// Identity. Defaults to <see cref="SimpleDna"/>.
    /// </param>
    /// <remarks>
    /// Identity = identical columns / alignment length × 100 (length includes gap
    /// columns); Similarity = (identical + similar) columns / length × 100;
    /// Gaps = gap columns / length × 100. Denominator is the full alignment length
    /// including gaps, per the EMBOSS needle definition.
    /// </remarks>
    public static AlignmentStatistics CalculateStatistics(
        AlignmentResult alignment,
        ScoringMatrix? scoring = null)
    {
        ArgumentNullException.ThrowIfNull(alignment);

        if (string.IsNullOrEmpty(alignment.AlignedSequence1))
            return AlignmentStatistics.Empty;

        var score = scoring ?? SimpleDna;

        int matches = 0, mismatches = 0, gaps = 0, similarSubstitutions = 0;
        int alignmentLength = alignment.AlignedSequence1.Length;

        for (int i = 0; i < alignmentLength; i++)
        {
            char c1 = alignment.AlignedSequence1[i];
            char c2 = alignment.AlignedSequence2[i];

            if (c1 == '-' || c2 == '-')
            {
                gaps++;
            }
            else if (c1 == c2)
            {
                matches++;
            }
            else
            {
                mismatches++;
                // A non-identical column is "similar" iff its substitution score is
                // positive (Rice, Longden & Bleasby 2000).
                if (score.Mismatch > 0)
                    similarSubstitutions++;
            }
        }

        // Length (including gaps) is the EMBOSS needle denominator for all three
        // percentages. Source: EMBOSS needle documentation (Rice et al. 2000).
        int similar = matches + similarSubstitutions;
        double identity = (double)matches / alignmentLength * 100;
        double similarity = (double)similar / alignmentLength * 100;
        double gapPercent = (double)gaps / alignmentLength * 100;

        return new AlignmentStatistics(
            Matches: matches,
            Mismatches: mismatches,
            Gaps: gaps,
            AlignmentLength: alignmentLength,
            Identity: identity,
            Similarity: similarity,
            GapPercent: gapPercent);
    }

    /// <summary>
    /// Generates a visual three-line pairwise alignment using the EMBOSS srspair
    /// markup legend ('|' identity, ':' similarity, ' ' gap/mismatch).
    /// </summary>
    /// <param name="alignment">The pairwise alignment to render.</param>
    /// <param name="lineWidth">Residues per block (must be positive; default 60).</param>
    /// <param name="scoring">
    /// Scoring model used to mark similar (positive-score) columns; defaults to
    /// <see cref="SimpleDna"/>. See <see cref="CalculateStatistics"/>.
    /// </param>
    public static string FormatAlignment(
        AlignmentResult alignment,
        int lineWidth = 60,
        ScoringMatrix? scoring = null)
    {
        ArgumentNullException.ThrowIfNull(alignment);
        ArgumentOutOfRangeException.ThrowIfNegativeOrZero(lineWidth);

        if (string.IsNullOrEmpty(alignment.AlignedSequence1))
            return "";

        var score = scoring ?? SimpleDna;
        var sb = new StringBuilder();
        int length = alignment.AlignedSequence1.Length;

        for (int start = 0; start < length; start += lineWidth)
        {
            int end = Math.Min(start + lineWidth, length);

            // Emit '\n' explicitly (not AppendLine/Environment.NewLine) so the formatted alignment
            // is byte-identical across platforms — the EMBOSS-style block must not depend on OS CRLF.
            sb.Append(alignment.AlignedSequence1[start..end]).Append('\n');

            // Markup line (EMBOSS srspair legend).
            for (int i = start; i < end; i++)
            {
                char c1 = alignment.AlignedSequence1[i];
                char c2 = alignment.AlignedSequence2[i];

                if (c1 == '-' || c2 == '-')
                    sb.Append(GapOrMismatchMark);
                else if (c1 == c2)
                    sb.Append(IdentityMark);
                else if (score.Mismatch > 0)
                    sb.Append(SimilarityMark);
                else
                    sb.Append(GapOrMismatchMark);
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
    /// Linear gap model: every gap position scores <see cref="ScoringMatrix.GapExtend"/>
    /// (<see cref="ScoringMatrix.GapOpen"/> is not used, as in <c>GlobalAlign</c>). For two sequences
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
        ScoringMatrix? scoring = null)
    {
        ArgumentNullException.ThrowIfNull(sequences);

        var seqList = sequences.ToList();
        var trivial = TrivialMultipleAlignment(seqList);
        if (trivial != null)
            return trivial;

        var effectiveScoring = scoring ?? SimpleDna;
        var rawSeqs = seqList.Select(s => s.Sequence).ToArray();

        int centerIdx = SelectCenterByPairwiseScore(rawSeqs, effectiveScoring);
        string centerStr = rawSeqs[centerIdx];

        var pairwise = AlignAllToCenter(rawSeqs.Length, centerIdx,
            i => GlobalAlignLinearSpace(centerStr, rawSeqs[i], effectiveScoring, StarFullMatrixCells));
        return FinishStarAlignment(MergeCenterAlignments(centerStr, pairwise, centerIdx), effectiveScoring);
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
        ScoringMatrix? scoring = null)
    {
        ArgumentNullException.ThrowIfNull(sequences);

        var seqList = sequences.ToList();
        var trivial = TrivialMultipleAlignment(seqList);
        if (trivial != null)
            return trivial;

        var effectiveScoring = scoring ?? SimpleDna;
        var rawSeqs = seqList.Select(s => s.Sequence).ToArray();

        int centerIdx = SelectCenterSequence(seqList);
        string centerStr = rawSeqs[centerIdx];
        var centerTree = SuffixTree.SuffixTree.Build(centerStr);

        var pairwise = AlignAllToCenter(rawSeqs.Length, centerIdx,
            i => AnchorBasedAligner.AlignWithAnchors(centerStr, centerTree, rawSeqs[i], effectiveScoring));
        var merged = MergeCenterAlignments(centerStr, pairwise, centerIdx);

        // An empty input yields an empty anchored result: pad to the common length.
        int maxLen = merged.Max(s => s.Length);
        for (int i = 0; i < merged.Count; i++)
            if (merged[i].Length < maxLen)
                merged[i] = merged[i].PadRight(maxLen, '-');

        return FinishStarAlignment(merged, effectiveScoring);
    }

    // Pairs above this many DP cells are aligned by the exact linear-space recursion (bounded memory
    // per concurrent pair: 4M cells = 16 MB); smaller pairs use the canonical full-matrix NW.
    private const long StarFullMatrixCells = 4_000_000;

    /// <summary>Empty / single-sequence results shared by every MSA entry point; null when k ≥ 2.</summary>
    private static MultipleAlignmentResult? TrivialMultipleAlignment(List<DnaSequence> seqList)
    {
        if (seqList.Count == 0)
            return MultipleAlignmentResult.Empty;

        if (seqList.Count == 1)
        {
            return new MultipleAlignmentResult(
                AlignedSequences: new[] { seqList[0].Sequence },
                Consensus: seqList[0].Sequence,
                TotalScore: 0);
        }

        return null;
    }

    private static MultipleAlignmentResult FinishStarAlignment(List<string> merged, ScoringMatrix scoring)
    {
        int maxLen = merged.Max(s => s.Length);
        return new MultipleAlignmentResult(
            AlignedSequences: merged.ToArray(),
            Consensus: BuildConsensus(merged, maxLen),
            TotalScore: ComputeSumOfPairsScore(merged, scoring));
    }

    /// <summary>
    /// Runs <paramref name="alignToCenter"/> for every non-center sequence (in parallel for k ≥ 4;
    /// each slot is independent and written once, so the result does not depend on scheduling).
    /// </summary>
    private static AlignmentResult[] AlignAllToCenter(int k, int centerIdx, Func<int, AlignmentResult> alignToCenter)
    {
        var results = new AlignmentResult[k];

        if (k >= 4)
        {
            Parallel.For(0, k, i =>
            {
                if (i != centerIdx)
                    results[i] = alignToCenter(i);
            });
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
    /// Classic star alignment with the FIRST sequence as center and exact pairwise alignments
    /// (benchmark baseline for <see cref="MultipleAlignAnchored(IEnumerable{DnaSequence}, ScoringMatrix?)"/>).
    /// </summary>
    internal static MultipleAlignmentResult MultipleAlignClassic(
        List<DnaSequence> seqList, ScoringMatrix scoring)
    {
        var raw = seqList.Select(s => s.Sequence).ToArray();
        var pairwise = AlignAllToCenter(raw.Length, 0, i => GlobalAlignCore(raw[0], raw[i], scoring));
        return FinishStarAlignment(MergeCenterAlignments(raw[0], pairwise, 0), scoring);
    }

    /// <summary>
    /// Gusfield (1993) center: the sequence maximising the total optimal Needleman-Wunsch score to the
    /// other sequences (pair scores = last cell of the linear-space NW pass, identical to the
    /// canonical <see cref="GlobalAlignCore"/> score).
    /// Ties go to the lowest index.
    /// </summary>
    private static int SelectCenterByPairwiseScore(string[] seqs, ScoringMatrix scoring)
    {
        int k = seqs.Length;
        var pairScore = new int[k, k];

        var pairs = new List<(int I, int J)>(k * (k - 1) / 2);
        for (int i = 0; i < k; i++)
            for (int j = i + 1; j < k; j++)
                pairs.Add((i, j));

        Parallel.ForEach(pairs, pair =>
        {
            // Optimal NW score in O(n) memory (last row of the NW matrix; equals GlobalAlignCore(..).Score).
            int s = LinearLastRow(seqs[pair.I], 0, seqs[pair.I].Length, seqs[pair.J], 0, seqs[pair.J].Length, scoring, reverse: false)[seqs[pair.J].Length];
            pairScore[pair.I, pair.J] = s;
            pairScore[pair.J, pair.I] = s;
        });

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
    /// counter <see cref="SequenceExtensions.CountKmersSpan"/> (DNA input is already upper-case ACGT, so
    /// every window is a valid 4-mer). Always prefers a non-empty center (an empty center has no
    /// positions to carry the anchored reconciliation).
    /// </summary>
    private static int SelectCenterSequence(List<DnaSequence> sequences)
    {
        int k = sequences.Count;

        int firstNonEmpty = -1;
        for (int i = 0; i < k; i++)
        {
            if (sequences[i].Sequence.Length > 0)
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
            profiles[i] = sequences[i].Sequence.AsSpan().CountKmersSpan(KmerLen);
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
            if (sequences[i].Sequence.Length == 0)
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
    /// Gap characters participate in the vote; on tie, nucleotides are preferred.
    /// </summary>
    private static string BuildConsensus(List<string> aligned, int length)
    {
        var consensus = new StringBuilder(length);

        for (int pos = 0; pos < length; pos++)
        {
            var counts = new Dictionary<char, int> { ['A'] = 0, ['C'] = 0, ['G'] = 0, ['T'] = 0, ['-'] = 0 };

            foreach (var seq in aligned)
            {
                if (pos < seq.Length && counts.TryGetValue(seq[pos], out int value))
                    counts[seq[pos]] = ++value;
            }

            // Include all characters (including gaps) in majority vote.
            // On tie between gap and nucleotide, prefer nucleotide.
            char mostCommon = counts.OrderByDescending(kv => kv.Value)
                                   .ThenBy(kv => kv.Key == '-' ? 1 : 0)
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
    private static int ComputeSumOfPairsScore(List<string> aligned, ScoringMatrix scoring)
    {
        if (aligned.Count < 2) return 0;

        int length = aligned[0].Length;
        int totalScore = 0;

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
                            totalScore += scoring.GapExtend; // Gap-nucleotide penalty
                        else if (ci == cj)
                            totalScore += scoring.Match;
                        else
                            totalScore += scoring.Mismatch;
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
    //      distance matrix. Distance = 1 - fractional identity, where fractional identity =
    //      identical columns / alignment length of the pairwise NW alignment. (Feng-Doolittle
    //      converts pairwise alignment scores to distances; identity-based distance is the
    //      standard, simplest such conversion used by Clustal-style tools.)
    //   2. Build a guide tree by UPGMA over that distance matrix — "an efficient clustering
    //      method such as neighbor-joining or unweighted pair group method with arithmetic
    //      mean (UPGMA)" (Wikipedia "Multiple sequence alignment"). UPGMA is the classic
    //      Feng-Doolittle choice. The UPGMA cluster-merge formula mirrors
    //      PhylogeneticAnalyzer.BuildUPGMA (proportional averaging, Sokal & Michener 1958);
    //      it is reimplemented here because the Phylogenetics project already depends on this
    //      Alignment project, so a reverse reference would be circular.
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
    /// <item>All pairwise Needleman-Wunsch alignments give a distance matrix using
    /// distance = 1 − fractional identity (identical columns / pairwise alignment length).</item>
    /// <item>A UPGMA guide tree is built over that matrix (the classic Feng-Doolittle choice).</item>
    /// <item>Profiles are aligned progressively from the tips with the NW recurrence over
    /// columns and sum-of-pairs profile scoring, enforcing "once a gap, always a gap".</item>
    /// </list>
    /// </remarks>
    /// <param name="sequences">Collection of sequences to align.</param>
    /// <param name="scoring">Scoring matrix (default: <see cref="SimpleDna"/>).</param>
    /// <returns>Multiple alignment result (aligned rows, majority consensus, sum-of-pairs score).</returns>
    /// <exception cref="ArgumentNullException">When <paramref name="sequences"/> is null.</exception>
    public static MultipleAlignmentResult MultipleAlignProgressive(
        IEnumerable<DnaSequence> sequences,
        ScoringMatrix? scoring = null)
    {
        ArgumentNullException.ThrowIfNull(sequences);

        var seqList = sequences.ToList();
        if (seqList.Count == 0)
            return MultipleAlignmentResult.Empty;

        if (seqList.Count == 1)
        {
            return new MultipleAlignmentResult(
                AlignedSequences: new[] { seqList[0].Sequence },
                Consensus: seqList[0].Sequence,
                TotalScore: 0);
        }

        var effectiveScoring = scoring ?? SimpleDna;
        int k = seqList.Count;
        var rawSeqs = seqList.Select(s => s.Sequence).ToArray();

        // Step 1: pairwise NW alignments -> distance matrix (1 - fractional identity).
        var distance = new double[k, k];
        for (int i = 0; i < k; i++)
        {
            for (int j = i + 1; j < k; j++)
            {
                double d = PairwiseIdentityDistance(rawSeqs[i], rawSeqs[j], effectiveScoring);
                distance[i, j] = d;
                distance[j, i] = d;
            }
        }

        // Step 2: UPGMA guide tree.
        ProgressiveGuideNode root = BuildProgressiveGuideTree(k, distance);

        // Step 3: progressive alignment along the tree (post-order: tips first).
        Profile aligned = AlignProfileSubtree(root, rawSeqs, effectiveScoring);

        // Reproject rows back into input order (the tree visited leaves in tree order).
        var orderedRows = new string[k];
        for (int r = 0; r < aligned.Rows.Count; r++)
            orderedRows[aligned.SequenceIndices[r]] = aligned.Rows[r];

        var alignedList = orderedRows.ToList();
        int maxLen = alignedList.Count == 0 ? 0 : alignedList.Max(s => s.Length);
        string consensus = BuildConsensus(alignedList, maxLen);
        int spScore = ComputeSumOfPairsScore(alignedList, effectiveScoring);

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
    private static double PairwiseIdentityDistance(string a, string b, ScoringMatrix scoring)
    {
        if (a.Length == 0 && b.Length == 0)
            return 0.0;
        if (a.Length == 0 || b.Length == 0)
            return 1.0; // no shared residues possible

        var pa = GlobalAlignCore(a, b, scoring);

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
    /// A node of the progressive-alignment guide tree. Leaves carry a single input-sequence
    /// index; internal nodes carry two children. (Distinct from PhylogeneticAnalyzer.PhyloNode,
    /// which lives in a project that already depends on this one.)
    /// </summary>
    private sealed class ProgressiveGuideNode
    {
        public int LeafIndex = -1; // >= 0 for leaves, -1 for internal nodes
        public ProgressiveGuideNode? Left;
        public ProgressiveGuideNode? Right;
        public bool IsLeaf => Left == null && Right == null;
    }

    /// <summary>
    /// Builds a UPGMA guide tree over the pairwise distance matrix. Uses proportional
    /// (cluster-size weighted) averaging, identical to PhylogeneticAnalyzer.BuildUPGMA
    /// (Sokal &amp; Michener 1958). Tie-breaking is deterministic: the lowest (i, j) index
    /// pair wins (distances equal within <see cref="ProfileScoreEpsilon"/> tie), so hand-derived test cases have a unique tree.
    /// </summary>
    private static ProgressiveGuideNode BuildProgressiveGuideTree(int k, double[,] distance)
    {
        var nodes = new Dictionary<int, ProgressiveGuideNode>();
        var sizes = new Dictionary<int, int>();
        var dist = new Dictionary<(int, int), double>();

        for (int i = 0; i < k; i++)
        {
            nodes[i] = new ProgressiveGuideNode { LeafIndex = i };
            sizes[i] = 1;
        }
        for (int i = 0; i < k; i++)
            for (int j = i + 1; j < k; j++)
            {
                dist[(i, j)] = distance[i, j];
                dist[(j, i)] = distance[i, j];
            }

        // Keep active clusters in a sorted list so tie-breaking is by ascending index.
        var active = Enumerable.Range(0, k).ToList();

        while (active.Count > 1)
        {
            double minDist = double.MaxValue;
            int minI = -1, minJ = -1;
            for (int ii = 0; ii < active.Count; ii++)
            {
                for (int jj = ii + 1; jj < active.Count; jj++)
                {
                    int a = active[ii];
                    int b = active[jj];
                    var key = a < b ? (a, b) : (b, a);
                    double d = dist.GetValueOrDefault(key, 0);
                    // Strictly smaller beyond rounding noise: mathematically equal averaged distances
                    // (e.g. 0.3 + 0.4 vs 0.7) must tie so that the lowest-index rule decides.
                    if (d < minDist - ProfileScoreEpsilon)
                    {
                        minDist = d;
                        minI = a;
                        minJ = b;
                    }
                }
            }

            int newSize = sizes[minI] + sizes[minJ];
            var newNode = new ProgressiveGuideNode
            {
                Left = nodes[minI],
                Right = nodes[minJ]
            };

            foreach (int c in active)
            {
                if (c == minI || c == minJ) continue;
                var keyIK = minI < c ? (minI, c) : (c, minI);
                var keyJK = minJ < c ? (minJ, c) : (c, minJ);
                double dIK = dist.GetValueOrDefault(keyIK, 0);
                double dJK = dist.GetValueOrDefault(keyJK, 0);
                double newDist = (dIK * sizes[minI] + dJK * sizes[minJ]) / newSize;
                var newKey = minI < c ? (minI, c) : (c, minI);
                dist[newKey] = newDist;
            }

            nodes[minI] = newNode;
            sizes[minI] = newSize;
            nodes.Remove(minJ);
            sizes.Remove(minJ);
            active.Remove(minJ);
            // active stays ascending: minJ removed, minI keeps its (smaller) slot.
        }

        return nodes[active[0]];
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
        ProgressiveGuideNode node, string[] rawSeqs, ScoringMatrix scoring)
    {
        if (node.IsLeaf)
        {
            var leaf = new Profile();
            leaf.Rows.Add(rawSeqs[node.LeafIndex]);
            leaf.SequenceIndices.Add(node.LeafIndex);
            return leaf;
        }

        var left = AlignProfileSubtree(node.Left!, rawSeqs, scoring);
        var right = AlignProfileSubtree(node.Right!, rawSeqs, scoring);
        return AlignProfiles(left, right, scoring);
    }

    /// <summary>
    /// Aligns two profiles with the Needleman-Wunsch recurrence over columns. The column-vs-column
    /// substitution score is the average match/mismatch over all cross-profile residue pairs
    /// (sum-of-pairs profile scoring); a gap column scores GapExtend per residue averaged.
    /// Existing columns are never edited — gaps are inserted as whole new all-gap columns into one
    /// profile, enforcing "once a gap, always a gap" (Feng &amp; Doolittle 1987).
    /// </summary>
    private static Profile AlignProfiles(Profile p1, Profile p2, ScoringMatrix scoring)
    {
        int m = p1.Width;
        int n = p2.Width;

        // F(i,j): best score aligning first i columns of p1 with first j columns of p2.
        var f = new double[m + 1, n + 1];
        for (int i = 1; i <= m; i++)
            f[i, 0] = f[i - 1, 0] + ProfileColumnVsGap(p1, i - 1, p2.Rows.Count, scoring);
        for (int j = 1; j <= n; j++)
            f[0, j] = f[0, j - 1] + ProfileColumnVsGap(p2, j - 1, p1.Rows.Count, scoring);

        for (int i = 1; i <= m; i++)
        {
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
    private static double ProfileColumnScore(Profile p1, int col1, Profile p2, int col2, ScoringMatrix scoring)
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
                    total += a == b ? scoring.Match : scoring.Mismatch;
            }
        }
        return total / pairs;
    }

    /// <summary>
    /// Score for aligning one profile column against an all-gap column of the other profile of
    /// <paramref name="otherRowCount"/> rows: average gap penalty over all cross pairs.
    /// Residue-vs-gap costs GapExtend; gap-vs-gap is neutral.
    /// </summary>
    private static double ProfileColumnVsGap(Profile p, int col, int otherRowCount, ScoringMatrix scoring)
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
    /// <item>Build the progressive (Feng-Doolittle / UPGMA) seed alignment and keep its guide tree.</item>
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
        int maxIterations = DefaultMaxRefinementIterations)
    {
        ArgumentNullException.ThrowIfNull(sequences);
        ArgumentOutOfRangeException.ThrowIfLessThan(maxIterations, 1);

        var seqList = sequences.ToList();
        if (seqList.Count == 0)
            return MultipleAlignmentResult.Empty;

        if (seqList.Count == 1)
        {
            return new MultipleAlignmentResult(
                AlignedSequences: new[] { seqList[0].Sequence },
                Consensus: seqList[0].Sequence,
                TotalScore: 0);
        }

        var effectiveScoring = scoring ?? SimpleDna;
        int k = seqList.Count;
        var rawSeqs = seqList.Select(s => s.Sequence).ToArray();

        // Step 1: progressive seed (same pipeline as MultipleAlignProgressive) + retain guide tree.
        var distance = new double[k, k];
        for (int i = 0; i < k; i++)
            for (int j = i + 1; j < k; j++)
            {
                double d = PairwiseIdentityDistance(rawSeqs[i], rawSeqs[j], effectiveScoring);
                distance[i, j] = d;
                distance[j, i] = d;
            }

        ProgressiveGuideNode root = BuildProgressiveGuideTree(k, distance);
        Profile current = AlignProfileSubtree(root, rawSeqs, effectiveScoring);

        // Step 2-5: iterative refinement over the guide-tree edges.
        current = RefineByTreePartitioning(current, root, effectiveScoring, maxIterations);

        // Reproject rows back into input order.
        var orderedRows = new string[k];
        for (int r = 0; r < current.Rows.Count; r++)
            orderedRows[current.SequenceIndices[r]] = current.Rows[r];

        var alignedList = orderedRows.ToList();
        int maxLen = alignedList.Count == 0 ? 0 : alignedList.Max(s => s.Length);
        string consensus = BuildConsensus(alignedList, maxLen);
        int spScore = ComputeSumOfPairsScore(alignedList, effectiveScoring);

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
        Profile current, ProgressiveGuideNode root, ScoringMatrix scoring, int maxIterations)
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
        int CurrentSpScore(Profile p) => ComputeSumOfPairsScore(p.Rows, scoring);

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
                Profile candidate = AlignProfiles(sideA, sideB, scoring);
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
    //   3. PROGRESSIVE ALIGNMENT (pp.209-210). A UPGMA guide tree (the existing
    //      BuildProgressiveGuideTree; the paper uses NJ — merge order only, same library and objective)
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
    /// <item><b>Progressively align</b> along a UPGMA guide tree with dynamic programming whose
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
        ScoringMatrix? scoring = null)
    {
        ArgumentNullException.ThrowIfNull(sequences);

        var seqList = sequences.ToList();
        if (seqList.Count == 0)
            return MultipleAlignmentResult.Empty;

        if (seqList.Count == 1)
        {
            return new MultipleAlignmentResult(
                AlignedSequences: new[] { seqList[0].Sequence },
                Consensus: seqList[0].Sequence,
                TotalScore: 0);
        }

        var effectiveScoring = scoring ?? SimpleDna;
        int k = seqList.Count;
        var rawSeqs = seqList.Select(s => s.Sequence).ToArray();

        // Steps 1-2: primary library (global + local, signal-added) then triplet extension.
        var library = new ExtendedLibrary(BuildPrimaryLibrary(rawSeqs, effectiveScoring));

        // Step 3a: UPGMA guide tree over (1 - fractional identity), reusing the progressive pipeline.
        var distance = new double[k, k];
        for (int i = 0; i < k; i++)
            for (int j = i + 1; j < k; j++)
            {
                double d = PairwiseIdentityDistance(rawSeqs[i], rawSeqs[j], effectiveScoring);
                distance[i, j] = d;
                distance[j, i] = d;
            }
        ProgressiveGuideNode root = BuildProgressiveGuideTree(k, distance);

        // Step 3b: progressive alignment driven by the EXTENDED-LIBRARY weights.
        ConsistencyProfile aligned = AlignConsistencySubtree(root, rawSeqs, library);

        // Reproject rows back to input order.
        var orderedRows = new string[k];
        for (int r = 0; r < aligned.Rows.Count; r++)
            orderedRows[aligned.SequenceIndices[r]] = aligned.Rows[r];

        var alignedList = orderedRows.ToList();
        int maxLen = alignedList.Count == 0 ? 0 : alignedList.Max(s => s.Length);
        string consensus = BuildConsensus(alignedList, maxLen);
        // TotalScore is reported in the same SP currency as the sibling aligners for comparability.
        int spScore = ComputeSumOfPairsScore(alignedList, effectiveScoring);

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

        public ExtendedLibrary(Dictionary<ResiduePairKey, double> primary)
        {
            _primary = primary;
            _neighbours = BuildNeighbourIndex(primary);
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

            _cache[key] = total;
            return total;
        }
    }

    /// <summary>Primary library: global (NW) + local (SW) pairwise alignments of every sequence pair, signal-added.</summary>
    private static Dictionary<ResiduePairKey, double> BuildPrimaryLibrary(string[] seqs, ScoringMatrix scoring)
    {
        var primary = new Dictionary<ResiduePairKey, double>();
        int k = seqs.Length;

        for (int i = 0; i < k; i++)
        {
            for (int j = i + 1; j < k; j++)
            {
                AddAlignmentToPrimary(primary, i, j, seqs[i], seqs[j],
                    GlobalAlignCore(seqs[i], seqs[j], scoring));

                if (seqs[i].Length > 0 && seqs[j].Length > 0)
                {
                    var local = LocalAlignCore(seqs[i], seqs[j], scoring);
                    AddLocalAlignmentToPrimary(primary, i, j, seqs[i], seqs[j], local);
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

        var f = new double[m + 1, n + 1];
        // First row/column: aligning a column against a gap scores 0 (zero gap penalty).
        // (f[i,0] = f[0,j] = 0 for all i,j.)

        for (int i = 1; i <= m; i++)
        {
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

    #endregion
}

