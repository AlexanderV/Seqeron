using System.Text;

namespace Seqeron.Genomics.Alignment;

// =====================================================================================================================
// Linear-space pairwise alignment (audit A1-6).
//
// The full-matrix programs (LinearAlignCore / AffineAlignCore / LocalAlignCore) need O(m·n) memory (4 bytes per cell for
// the linear model, 12 for the affine one: 20 kb × 20 kb = 4.8 GB). This file adds the exact linear-space counterparts:
//
//   * Hirschberg (1975, Commun ACM 18:341): the optimal path of a global alignment is found by divide and conquer on the
//     middle row of sequence 1; only one DP row is alive at a time, the path is rebuilt from the sub-problems.
//   * Myers & Miller (1988, CABIOS 4:11): the same idea under the affine gap model (Gotoh's three states). The split
//     vertex carries its STATE (a vertical/horizontal gap may run through the middle row, so the gap opening must be
//     charged once, not once per half); here a sub-problem is "global alignment from an arrival state s0 at the origin
//     to a forced state s1 at the end", which is exactly the type-1/type-2 midpoint bookkeeping of the paper.
//   * Huang, Hardison & Miller (1990, CABIOS 6:373; SIM, EMBOSS matcher) / Huang & Miller (1991): linear-space local
//     alignment = locate the optimal local interval with an O(n)-memory forward pass, then solve the interval as a
//     global problem. The same two-phase scheme gives the ends-free (semi-global) policies: locate the start cell and
//     the end cell, then solve the global rectangle between them.
//
// Midpoint search. Instead of the backward score pass of Hirschberg/Myers-Miller, the crossing of the middle row is
// obtained by CARRYING it: after the forward pass over the top half, the forward pass continues through the bottom half
// and every DP value (state) carries the column (and state) at which the traceback of the full-matrix program would
// leave through the middle row (payload[cell] = payload[predecessor chosen by the traceback]). The traceback of
// the full-matrix aligners is a deterministic function of the forward values (diagonal > up > left, resp.
// M > X > Y, first tight predecessor), so the split vertex is a vertex of THE path the full-matrix program returns, and
// - because every sub-problem again has the lexicographically preferred optimal path of the whole - the concatenation
// of the sub-alignments is row-for-row identical to the full-matrix result, ties included (proof sketch in
// docs/algorithms/Alignment/Linear_Space_Alignment.md). Score, rows and coordinates are therefore identical to
// GlobalAlign/GlobalAlignAffine/LocalAlign(Affine)/SemiGlobalAlign(Affine) and the automatic switch above
// LinearSpaceCellThreshold cells does not change any output.
// =====================================================================================================================
public static partial class SequenceAligner
{
    #region Linear-space alignment

    /// <summary>
    /// Number of dynamic-programming cells <c>(m + 1)·(n + 1)</c> above which every pairwise aligner of this class
    /// (<c>GlobalAlign</c>, <c>GlobalAlignAffine</c>, <c>LocalAlign</c>, <c>LocalAlignAffine</c>, <c>SemiGlobalAlign</c>,
    /// <c>SemiGlobalAlignAffine</c>, scalar and <see cref="MatrixScoring"/> overloads, the ends-free policies and the
    /// linear/affine special cases of the general-gap aligners) switches automatically from the full O(m·n) matrix to the
    /// exact linear-space program. 25 000 000 cells = 100 MB (linear model) / 300 MB (affine) — the same value as
    /// <c>AnchorBasedAligner</c>'s full-matrix limit. The result is identical either way (score, rows, coordinates);
    /// only the memory profile (O(m + n) instead of O(m·n)) and the running time (about 1.3-2× the full-matrix fill) change.
    /// </summary>
    public const long LinearSpaceCellThreshold = 25_000_000;

    /// <summary>
    /// Default size (cells) of the sub-problems that the linear-space program solves with the full-matrix core
    /// (2^22 = 4 194 304 cells = 17 MB linear / 50 MB affine). Smaller values recurse deeper (used by the tests to force the
    /// divide and conquer on tiny inputs); larger values are clamped to <see cref="LinearSpaceCellThreshold"/>.
    /// </summary>
    public const long LinearSpaceBaseCells = 1L << 22;

    private static bool ExceedsFullMatrixLimit(int m, int n) => (long)(m + 1) * (n + 1) > LinearSpaceCellThreshold;

    /// <summary>
    /// True when the automatic switch applies. The local program additionally needs non-positive gap scores (a local
    /// alignment never starts with a gap, which the interval-then-global scheme relies on); with a positive gap score
    /// (meaningless for alignment) the full matrix is used as before.
    /// </summary>
    private static bool LinearSpaceApplies(int m, int n, PairScoring scoring, bool local) =>
        ExceedsFullMatrixLimit(m, n) && (!local || (scoring.GapOpen <= 0 && scoring.GapExtend <= 0));

    private static long PackCell(int i, int j) => ((long)i << 32) | (uint)j;

    private static int CellRow(long packed) => (int)(packed >> 32);

    private static int CellCol(long packed) => (int)(packed & 0xFFFFFFFFL);

    /// <summary>Result of the locate phase: the start vertex, end cell/state and the optimal score of the aligned interval.</summary>
    private readonly record struct LinearSpacePath(int StartRow, int StartCol, int EndRow, int EndCol, int EndState, int Score);

    /// <summary>
    /// Exact linear-space alignment: same score, rows and coordinates as the full-matrix program of the same model
    /// (<paramref name="affine"/>: Gotoh three states, otherwise linear gaps; <paramref name="local"/>; end-gap policy
    /// <paramref name="ends"/>), O(m + n) working memory plus the O(m + n) output.
    /// </summary>
    /// <param name="baseCells">Sub-problems of at most this many cells (clamped to <see cref="LinearSpaceCellThreshold"/>)
    /// are solved by the full-matrix core.</param>
    private static AlignmentResult LinearSpaceAlign(
        string a, string b, PairScoring sc, bool affine, bool local, EndGaps ends, long baseCells,
        CancellationToken ct)
    {
        if (local)
        {
            if (sc.GapOpen > 0 || sc.GapExtend > 0)
                throw new ArgumentException("Linear-space local alignment needs non-positive gap scores (GapOpen <= 0 and GapExtend <= 0).", nameof(sc));
            ends = EndGaps.None;
        }
        baseCells = Math.Clamp(baseCells, 1, LinearSpaceCellThreshold);
        int m = a.Length;
        int n = b.Length;
        var type = ResultType(local, ends);
        ct.ThrowIfCancellationRequested();

        // ---- Phase 1: locate the aligned interval (start vertex, end cell, optimal score) ----
        LinearSpacePath path;
        if (local)
        {
            path = affine ? LocateAffine(a, b, sc, local: true, EndGaps.None, ct) : LocateLinearLocal(a, b, sc, ct);
            if (path.Score <= 0)
                return EmptyLocalResult();
        }
        else if (ends == EndGaps.None)
        {
            path = new LinearSpacePath(0, 0, m, n, -1, 0); // the whole matrix; the score comes from the solve
        }
        else
        {
            path = affine ? LocateAffine(a, b, sc, local: false, ends, ct) : LocateLinearEndsFree(a, b, sc, ends, ct);
        }

        // ---- Phase 2: global alignment of the rectangle [start, end] (Hirschberg / Myers-Miller) ----
        var out1 = new StringBuilder(m + n);
        var out2 = new StringBuilder(m + n);
        int startI = path.StartRow, startJ = path.StartCol;
        if (!local)
        {
            // Free leading overhang: the traceback stopped in row 0 / column 0 (nothing when it reached the origin).
            for (int k = 0; k < startJ && startI == 0; k++) { out1.Append('-'); out2.Append(b[k]); }
            for (int k = 0; k < startI && startJ == 0; k++) { out1.Append(a[k]); out2.Append('-'); }
        }
        else
        {
            // A local alignment starts with its first aligned pair (the traceback stops after consuming it).
            out1.Append(a[startI - 1]);
            out2.Append(b[startJ - 1]);
        }

        int score = affine
            ? AffineSolve(a, startI, path.EndRow, b, startJ, path.EndCol, sc, StateMatch, path.EndState, baseCells, out1, out2, ct)
            : LinearSolve(a, startI, path.EndRow, b, startJ, path.EndCol, sc, baseCells, out1, out2, ct);
        if (ends != EndGaps.None || local)
            score = path.Score;

        if (!local)
        {
            for (int k = path.EndCol; k < n; k++) { out1.Append('-'); out2.Append(b[k]); }
            for (int k = path.EndRow; k < m; k++) { out1.Append(a[k]); out2.Append('-'); }
        }

        return new AlignmentResult(
            out1.ToString(), out2.ToString(), score, type,
            startI - (local ? 1 : 0),
            startJ - (local ? 1 : 0),
            path.EndRow - 1,
            path.EndCol - 1);
    }

    // ------------------------------------------------------------------------------------------------------------
    // Linear gap model
    // ------------------------------------------------------------------------------------------------------------

    /// <summary>
    /// Global (Needleman-Wunsch, linear gaps) alignment of a[i0..i1) against b[j0..j1) by divide and conquer; appends the
    /// rows to the builders and returns the optimal score. Sub-problems up to <paramref name="baseCells"/> cells are solved
    /// by <see cref="LinearAlignCore"/>.
    /// </summary>
    private static int LinearSolve(
        string a, int i0, int i1, string b, int j0, int j1, PairScoring sc, long baseCells,
        StringBuilder out1, StringBuilder out2, CancellationToken ct)
    {
        int m = i1 - i0;
        int n = j1 - j0;
        if (m <= 1 || n <= 1 || (long)(m + 1) * (n + 1) <= baseCells)
        {
            var r = LinearAlignCore(a.Substring(i0, m), b.Substring(j0, n), sc, EndGaps.None, ct, null, forceFullMatrix: true);
            out1.Append(r.AlignedSequence1);
            out2.Append(r.AlignedSequence2);
            return r.Score;
        }

        int mid = m / 2;
        int cross = LinearMidpoint(a, i0, b, j0, m, n, mid, sc, ct, out int total);
        LinearSolve(a, i0, i0 + mid, b, j0, j0 + cross, sc, baseCells, out1, out2, ct);
        LinearSolve(a, i0 + mid, i1, b, j0 + cross, j1, sc, baseCells, out1, out2, ct);
        return total;
    }

    /// <summary>Optimal linear-gap Needleman-Wunsch score of a against b in O(n) memory (the last cell of the forward pass;
    /// equals <c>GlobalAlignCore(..).Score</c>).</summary>
    private static int LinearScoreOnly(string a, string b, PairScoring sc)
    {
        LinearMidpoint(a, 0, b, 0, a.Length, b.Length, a.Length, sc, default, out int total);
        return total;
    }

    /// <summary>
    /// Forward pass over all rows of the sub-problem; from row <paramref name="mid"/> on every cell carries the column at
    /// which the full-matrix traceback (diagonal &gt; up &gt; left) crosses row <paramref name="mid"/>. Returns that
    /// column for the end cell (m, n); <paramref name="total"/> is the optimal score F(m, n).
    /// </summary>
    private static int LinearMidpoint(
        string a, int i0, string b, int j0, int m, int n, int mid, PairScoring sc, CancellationToken ct, out int total)
    {
        int d = sc.GapExtend;
        var prev = new int[n + 1];
        var cur = new int[n + 1];
        var crossPrev = new int[n + 1];
        var crossCur = new int[n + 1];
        for (int j = 0; j <= n; j++)
            prev[j] = j * d;

        for (int r = 1; r <= m; r++)
        {
            if ((r & 63) == 0)
                ct.ThrowIfCancellationRequested();

            char ca = a[i0 + r - 1];
            bool carry = r > mid;
            cur[0] = r * d;
            if (carry)
                crossCur[0] = crossPrev[0];

            for (int j = 1; j <= n; j++)
            {
                int diag = prev[j - 1] + sc.Sub(ca, b[j0 + j - 1]);
                int up = prev[j] + d;
                int left = cur[j - 1] + d;
                int v = Math.Max(diag, Math.Max(up, left));
                cur[j] = v;
                if (carry)
                {
                    if (v == diag) crossCur[j] = crossPrev[j - 1];
                    else if (v == up) crossCur[j] = crossPrev[j];
                    else crossCur[j] = crossCur[j - 1];
                }
            }

            (prev, cur) = (cur, prev);
            (crossPrev, crossCur) = (crossCur, crossPrev);
            if (r == mid)
                for (int j = 0; j <= n; j++)
                    crossPrev[j] = j; // the traceback leaves row mid at the cell it is standing on
        }

        total = prev[n];
        return crossPrev[n];
    }

    /// <summary>
    /// Locate pass of the linear local program: forward pass with a rolling row; every cell with H &gt; 0 carries the cell
    /// of the FIRST aligned pair of its traceback (the traceback stops in front of the first predecessor with H = 0).
    /// The end cell is the first row-major maximum, as in <see cref="LocalAlignCore"/>.
    /// </summary>
    private static LinearSpacePath LocateLinearLocal(string a, string b, PairScoring sc, CancellationToken ct)
    {
        int m = a.Length;
        int n = b.Length;
        int d = sc.GapExtend;
        var hp = new int[n + 1];
        var hc = new int[n + 1];
        var sp = new long[n + 1];
        var scur = new long[n + 1];
        int best = 0, bestI = 0, bestJ = 0;
        long bestStart = 0;

        for (int i = 1; i <= m; i++)
        {
            if ((i & 63) == 0)
                ct.ThrowIfCancellationRequested();

            char ca = a[i - 1];
            hc[0] = 0;
            for (int j = 1; j <= n; j++)
            {
                int diag = hp[j - 1] + sc.Sub(ca, b[j - 1]);
                int up = hp[j] + d;
                int left = hc[j - 1] + d;
                int v = Math.Max(diag, Math.Max(up, left));
                if (v <= 0)
                {
                    hc[j] = 0;
                    continue;
                }

                hc[j] = v;
                int predValue;
                long predStart;
                if (v == diag) { predValue = hp[j - 1]; predStart = sp[j - 1]; }
                else if (v == up) { predValue = hp[j]; predStart = sp[j]; }
                else { predValue = hc[j - 1]; predStart = scur[j - 1]; }

                long start = predValue == 0 ? PackCell(i, j) : predStart;
                scur[j] = start;
                if (v > best)
                {
                    best = v;
                    bestI = i;
                    bestJ = j;
                    bestStart = start;
                }
            }

            (hp, hc) = (hc, hp);
            (sp, scur) = (scur, sp);
        }

        return new LinearSpacePath(CellRow(bestStart), CellCol(bestStart), bestI, bestJ, -1, best);
    }

    /// <summary>
    /// Locate pass of the linear ends-free program: forward pass with a rolling row; every cell carries the stop vertex of
    /// its traceback (row 0 when the sequence-2 start is free, column 0 when the sequence-1 start is free, otherwise the origin).
    /// The end cell follows <see cref="LinearAlignCore"/>: first (smallest j) maximum of the last row when
    /// <see cref="EndGaps.FreeEnd2"/>, then the first strictly larger maximum of the last column when <see cref="EndGaps.FreeEnd1"/>.
    /// </summary>
    private static LinearSpacePath LocateLinearEndsFree(string a, string b, PairScoring sc, EndGaps ends, CancellationToken ct)
    {
        bool freeStart1 = ends.HasFlag(EndGaps.FreeStart1);
        bool freeStart2 = ends.HasFlag(EndGaps.FreeStart2);
        bool freeEnd1 = ends.HasFlag(EndGaps.FreeEnd1);
        bool freeEnd2 = ends.HasFlag(EndGaps.FreeEnd2);
        int m = a.Length;
        int n = b.Length;
        int d = sc.GapExtend;
        var fp = new int[n + 1];
        var fc = new int[n + 1];
        var op = new long[n + 1];
        var oc = new long[n + 1];

        for (int j = 0; j <= n; j++)
        {
            fp[j] = freeStart2 ? 0 : j * d;
            op[j] = j == 0 || freeStart2 ? PackCell(0, j) : op[j - 1];
        }

        int colBest = int.MinValue, colRow = 0;
        long colOrigin = 0;
        if (freeEnd1)
        {
            colBest = fp[n];
            colOrigin = op[n];
        }

        for (int i = 1; i <= m; i++)
        {
            if ((i & 63) == 0)
                ct.ThrowIfCancellationRequested();

            char ca = a[i - 1];
            fc[0] = freeStart1 ? 0 : i * d;
            oc[0] = freeStart1 ? PackCell(i, 0) : op[0];
            for (int j = 1; j <= n; j++)
            {
                int diag = fp[j - 1] + sc.Sub(ca, b[j - 1]);
                int up = fp[j] + d;
                int left = fc[j - 1] + d;
                int v = Math.Max(diag, Math.Max(up, left));
                fc[j] = v;
                if (v == diag) oc[j] = op[j - 1];
                else if (v == up) oc[j] = op[j];
                else oc[j] = oc[j - 1];
            }

            (fp, fc) = (fc, fp);
            (op, oc) = (oc, op);
            if (freeEnd1 && fp[n] > colBest)
            {
                colBest = fp[n];
                colRow = i;
                colOrigin = op[n];
            }
        }

        int endI = m, endJ = n;
        int best = fp[n];
        long origin = op[n];
        if (freeEnd2 || freeEnd1)
        {
            best = int.MinValue;
            if (freeEnd2)
                for (int j = 0; j <= n; j++)
                    if (fp[j] > best) { best = fp[j]; endI = m; endJ = j; origin = op[j]; }
            if (freeEnd1 && colBest > best)
            {
                best = colBest;
                endI = colRow;
                endJ = n;
                origin = colOrigin;
            }
        }

        return new LinearSpacePath(CellRow(origin), CellCol(origin), endI, endJ, -1, best);
    }

    // ------------------------------------------------------------------------------------------------------------
    // Affine gap model (Gotoh three states, Myers-Miller state-carrying midpoint)
    // ------------------------------------------------------------------------------------------------------------

    /// <summary>
    /// Global affine alignment of a[i0..i1) against b[j0..j1) from arrival state <paramref name="startState"/> at the origin
    /// to <paramref name="endState"/> (-1: the best state at the end cell, ties M &gt; X &gt; Y) by divide and conquer;
    /// appends the rows and returns the score of the end state.
    /// </summary>
    private static int AffineSolve(
        string a, int i0, int i1, string b, int j0, int j1, PairScoring sc, int startState, int endState, long baseCells,
        StringBuilder out1, StringBuilder out2, CancellationToken ct)
    {
        int m = i1 - i0;
        int n = j1 - j0;
        if (m <= 1 || n <= 1 || (long)(m + 1) * (n + 1) <= baseCells)
        {
            var r = AffineAlignCore(a.Substring(i0, m), b.Substring(j0, n), sc, AlignmentType.Global, EndGaps.None,
                startState, endState, forceFullMatrix: true);
            out1.Append(r.AlignedSequence1);
            out2.Append(r.AlignedSequence2);
            return r.Score;
        }

        int mid = m / 2;
        AffineMidpoint(a, i0, b, j0, m, n, mid, sc, startState, endState, ct,
            out int crossCol, out int crossState, out int finalState, out int total);
        AffineSolve(a, i0, i0 + mid, b, j0, j0 + crossCol, sc, startState, crossState, baseCells, out1, out2, ct);
        AffineSolve(a, i0 + mid, i1, b, j0 + crossCol, j1, sc, crossState, finalState, baseCells, out1, out2, ct);
        return total;
    }

    private static int StatePayload(int state, int m, int x, int y) => state switch
    {
        StateMatch => m,
        StateGapInSeq2 => x,
        _ => y,
    };

    /// <summary>
    /// Forward Gotoh pass over all rows of the sub-problem with the origin vector of <paramref name="startState"/>. From row
    /// <paramref name="mid"/> on every (cell, state) carries "4·column + state" of the vertex at which the full-matrix
    /// traceback (M &gt; X &gt; Y, first tight predecessor) crosses row <paramref name="mid"/>.
    /// </summary>
    private static void AffineMidpoint(
        string a, int i0, string b, int j0, int m, int n, int mid, PairScoring sc, int startState, int endState,
        CancellationToken ct, out int crossCol, out int crossState, out int finalState, out int total)
    {
        int o = sc.GapOpen;
        int e = sc.GapExtend;
        int oe = o + e;
        var mp = new int[n + 1];
        var xp = new int[n + 1];
        var yp = new int[n + 1];
        var mc = new int[n + 1];
        var xc = new int[n + 1];
        var yc = new int[n + 1];
        var pmp = new int[n + 1];
        var pxp = new int[n + 1];
        var pyp = new int[n + 1];
        var pmc = new int[n + 1];
        var pxc = new int[n + 1];
        var pyc = new int[n + 1];

        mp[0] = startState == StateMatch ? 0 : AffineNegInf;
        xp[0] = startState == StateGapInSeq2 ? 0 : AffineNegInf;
        yp[0] = startState == StateGapInSeq1 ? 0 : AffineNegInf;
        for (int j = 1; j <= n; j++)
        {
            mp[j] = AffineNegInf;
            xp[j] = AffineNegInf;
            yp[j] = Math.Max(AffineNegInf, Max3(mp[j - 1] + oe, yp[j - 1] + e, xp[j - 1] + oe));
        }

        for (int r = 1; r <= m; r++)
        {
            if ((r & 63) == 0)
                ct.ThrowIfCancellationRequested();

            char ca = a[i0 + r - 1];
            bool carry = r > mid;

            mc[0] = AffineNegInf;
            yc[0] = AffineNegInf;
            xc[0] = Math.Max(AffineNegInf, Max3(mp[0] + oe, xp[0] + e, yp[0] + oe));
            if (carry)
            {
                pmc[0] = 0;
                pyc[0] = 0;
                int st = PickState(xc[0], mp[0] + oe, StateMatch, xp[0] + e, StateGapInSeq2, StateGapInSeq1);
                pxc[0] = StatePayload(st, pmp[0], pxp[0], pyp[0]);
            }

            for (int j = 1; j <= n; j++)
            {
                int s = sc.Sub(ca, b[j0 + j - 1]);

                int basis = Max3(mp[j - 1], xp[j - 1], yp[j - 1]);
                mc[j] = Math.Max(AffineNegInf, basis + s);
                int xv = Math.Max(AffineNegInf, Max3(mp[j] + oe, xp[j] + e, yp[j] + oe));
                xc[j] = xv;
                int yv = Math.Max(AffineNegInf, Max3(mc[j - 1] + oe, yc[j - 1] + e, xc[j - 1] + oe));
                yc[j] = yv;

                if (carry)
                {
                    int sm = PickState(basis, mp[j - 1], StateMatch, xp[j - 1], StateGapInSeq2, StateGapInSeq1);
                    pmc[j] = StatePayload(sm, pmp[j - 1], pxp[j - 1], pyp[j - 1]);
                    int sx = PickState(xv, mp[j] + oe, StateMatch, xp[j] + e, StateGapInSeq2, StateGapInSeq1);
                    pxc[j] = StatePayload(sx, pmp[j], pxp[j], pyp[j]);
                    int sy = PickState(yv, mc[j - 1] + oe, StateMatch, yc[j - 1] + e, StateGapInSeq1, StateGapInSeq2);
                    pyc[j] = StatePayload(sy, pmc[j - 1], pxc[j - 1], pyc[j - 1]);
                }
            }

            (mp, mc) = (mc, mp);
            (xp, xc) = (xc, xp);
            (yp, yc) = (yc, yp);
            (pmp, pmc) = (pmc, pmp);
            (pxp, pxc) = (pxc, pxp);
            (pyp, pyc) = (pyc, pyp);
            if (r == mid)
                for (int j = 0; j <= n; j++)
                {
                    pmp[j] = 4 * j + StateMatch;
                    pxp[j] = 4 * j + StateGapInSeq2;
                    pyp[j] = 4 * j + StateGapInSeq1;
                }
        }

        int bestEnd = Max3(mp[n], xp[n], yp[n]);
        finalState = endState >= 0 ? endState : PickState(bestEnd, mp[n], StateMatch, xp[n], StateGapInSeq2, StateGapInSeq1);
        total = StatePayload(finalState, mp[n], xp[n], yp[n]);
        int payload = StatePayload(finalState, pmp[n], pxp[n], pyp[n]);
        crossCol = payload >> 2;
        crossState = payload & 3;
    }

    /// <summary>
    /// Locate pass of the affine local and ends-free programs: Gotoh forward pass with rolling rows in which every
    /// (cell, state) carries the start vertex of its traceback: for the local program the cell of the first aligned pair
    /// (the M cell whose predecessor value is 0, where the traceback stops), otherwise the stop cell (row 0 when the
    /// sequence-2 start is free, column 0 when the sequence-1 start is free, the origin). End cell and state follow
    /// <see cref="AffineAlignCore"/> (local: first row-major maximum of M; ends-free: first maximum of max(M,X,Y)
    /// over the free last row, then the first strictly larger one over the free last column; state tie M &gt; X &gt; Y).
    /// </summary>
    private static LinearSpacePath LocateAffine(
        string a, string b, PairScoring sc, bool local, EndGaps ends, CancellationToken ct)
    {
        bool freeStart1 = ends.HasFlag(EndGaps.FreeStart1);
        bool freeStart2 = ends.HasFlag(EndGaps.FreeStart2);
        bool freeEnd1 = ends.HasFlag(EndGaps.FreeEnd1);
        bool freeEnd2 = ends.HasFlag(EndGaps.FreeEnd2);
        int m = a.Length;
        int n = b.Length;
        int o = sc.GapOpen;
        int e = sc.GapExtend;
        int oe = o + e;

        var mp = new int[n + 1];
        var xp = new int[n + 1];
        var yp = new int[n + 1];
        var mc = new int[n + 1];
        var xc = new int[n + 1];
        var yc = new int[n + 1];
        var pmp = new long[n + 1];
        var pxp = new long[n + 1];
        var pyp = new long[n + 1];
        var pmc = new long[n + 1];
        var pxc = new long[n + 1];
        var pyc = new long[n + 1];

        // Row 0 (borders exactly as AffineAlignCore).
        mp[0] = local ? AffineNegInf : 0;
        xp[0] = AffineNegInf;
        yp[0] = AffineNegInf;
        pmp[0] = PackCell(0, 0);
        for (int j = 1; j <= n; j++)
        {
            mp[j] = freeStart2 ? 0 : AffineNegInf;
            pmp[j] = PackCell(0, j);
            xp[j] = AffineNegInf;
            if (local || freeStart2)
            {
                yp[j] = AffineNegInf;
            }
            else
            {
                yp[j] = Math.Max(AffineNegInf, Max3(mp[j - 1] + oe, yp[j - 1] + e, xp[j - 1] + oe));
                int sy = PickState(yp[j], mp[j - 1] + oe, StateMatch, yp[j - 1] + e, StateGapInSeq1, StateGapInSeq2);
                pyp[j] = PayloadOf(sy, pmp[j - 1], pxp[j - 1], pyp[j - 1]);
            }
        }

        // End candidates.
        int bestLocal = 0, bestI = 0, bestJ = 0;
        long bestStart = 0;
        int colBest = int.MinValue, colRow = 0, colState = 0;
        long colStart = 0;
        if (freeEnd1)
        {
            int v = Max3(mp[n], xp[n], yp[n]);
            colBest = v;
            colState = PickState(v, mp[n], StateMatch, xp[n], StateGapInSeq2, StateGapInSeq1);
            colStart = PayloadOf(colState, pmp[n], pxp[n], pyp[n]);
        }

        for (int i = 1; i <= m; i++)
        {
            if ((i & 63) == 0)
                ct.ThrowIfCancellationRequested();

            char ca = a[i - 1];

            // Column 0.
            mc[0] = freeStart1 ? 0 : AffineNegInf;
            pmc[0] = PackCell(i, 0);
            yc[0] = AffineNegInf;
            if (local || freeStart1)
            {
                xc[0] = AffineNegInf;
            }
            else
            {
                xc[0] = Math.Max(AffineNegInf, Max3(mp[0] + oe, xp[0] + e, yp[0] + oe));
                int sx0 = PickState(xc[0], mp[0] + oe, StateMatch, xp[0] + e, StateGapInSeq2, StateGapInSeq1);
                pxc[0] = PayloadOf(sx0, pmp[0], pxp[0], pyp[0]);
            }

            for (int j = 1; j <= n; j++)
            {
                int s = sc.Sub(ca, b[j - 1]);

                int basis = Max3(mp[j - 1], xp[j - 1], yp[j - 1]);
                if (local && basis < 0)
                    basis = 0;
                int mv = Math.Max(AffineNegInf, basis + s);
                mc[j] = mv;
                int xv = Math.Max(AffineNegInf, Max3(mp[j] + oe, xp[j] + e, yp[j] + oe));
                xc[j] = xv;
                int yv = Math.Max(AffineNegInf, Max3(mc[j - 1] + oe, yc[j - 1] + e, xc[j - 1] + oe));
                yc[j] = yv;

                if (local && basis == 0)
                {
                    pmc[j] = PackCell(i, j); // first aligned pair: the traceback stops here
                }
                else
                {
                    int sm = PickState(basis, mp[j - 1], StateMatch, xp[j - 1], StateGapInSeq2, StateGapInSeq1);
                    pmc[j] = PayloadOf(sm, pmp[j - 1], pxp[j - 1], pyp[j - 1]);
                }

                int sx = PickState(xv, mp[j] + oe, StateMatch, xp[j] + e, StateGapInSeq2, StateGapInSeq1);
                pxc[j] = PayloadOf(sx, pmp[j], pxp[j], pyp[j]);
                int sy2 = PickState(yv, mc[j - 1] + oe, StateMatch, yc[j - 1] + e, StateGapInSeq1, StateGapInSeq2);
                pyc[j] = PayloadOf(sy2, pmc[j - 1], pxc[j - 1], pyc[j - 1]);

                if (local && mv > bestLocal)
                {
                    bestLocal = mv;
                    bestI = i;
                    bestJ = j;
                    bestStart = pmc[j];
                }
            }

            (mp, mc) = (mc, mp);
            (xp, xc) = (xc, xp);
            (yp, yc) = (yc, yp);
            (pmp, pmc) = (pmc, pmp);
            (pxp, pxc) = (pxc, pxp);
            (pyp, pyc) = (pyc, pyp);

            if (freeEnd1)
            {
                int v = Max3(mp[n], xp[n], yp[n]);
                if (v > colBest)
                {
                    colBest = v;
                    colRow = i;
                    colState = PickState(v, mp[n], StateMatch, xp[n], StateGapInSeq2, StateGapInSeq1);
                    colStart = PayloadOf(colState, pmp[n], pxp[n], pyp[n]);
                }
            }
        }

        if (local)
            return new LinearSpacePath(CellRow(bestStart), CellCol(bestStart), bestI, bestJ, StateMatch, bestLocal);

        int endI = m, endJ = n;
        int best = Max3(mp[n], xp[n], yp[n]);
        int endState = PickState(best, mp[n], StateMatch, xp[n], StateGapInSeq2, StateGapInSeq1);
        long start = PayloadOf(endState, pmp[n], pxp[n], pyp[n]);
        if (freeEnd2 || freeEnd1)
        {
            best = int.MinValue;
            if (freeEnd2)
                for (int j = 0; j <= n; j++)
                {
                    int v = Max3(mp[j], xp[j], yp[j]);
                    if (v > best)
                    {
                        best = v;
                        endI = m;
                        endJ = j;
                        endState = PickState(v, mp[j], StateMatch, xp[j], StateGapInSeq2, StateGapInSeq1);
                        start = PayloadOf(endState, pmp[j], pxp[j], pyp[j]);
                    }
                }
            if (freeEnd1 && colBest > best)
            {
                best = colBest;
                endI = colRow;
                endJ = n;
                endState = colState;
                start = colStart;
            }
        }

        return new LinearSpacePath(CellRow(start), CellCol(start), endI, endJ, endState, best);
    }

    private static long PayloadOf(int state, long m, long x, long y) => state switch
    {
        StateMatch => m,
        StateGapInSeq2 => x,
        _ => y,
    };


    // ---- public linear-space entry points ----------------------------------------------------------------------

    private static AlignmentResult ScalarLinearSpace(
        string sequence1, string sequence2, ScoringMatrix? scoring, bool affine, bool local, EndGaps ends, long baseCells)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        return LinearSpaceAlign(
            sequence1.ToUpperInvariant(), sequence2.ToUpperInvariant(), new PairScoring(scoring ?? SimpleDna),
            affine, local, ends, baseCells, default);
    }

    /// <summary>
    /// Exact linear-space global (Needleman-Wunsch, linear gaps) alignment (O(m + n) memory instead of O(m·n); see the region comment of this file and
    /// <c>docs/algorithms/Alignment/Linear_Space_Alignment.md</c>). Score, aligned rows, coordinates and tie behaviour are
    /// identical to <see cref="GlobalAlign(string, string, ScoringMatrix?)"/>, which switches to this program by itself above
    /// <see cref="LinearSpaceCellThreshold"/> cells; call this method to force it for any size.
    /// </summary>
    /// <param name="sequence1">First sequence (upper-cased internally).</param>
    /// <param name="sequence2">Second sequence.</param>
    /// <param name="scoring">Scoring (default <see cref="SimpleDna"/>).</param>
    /// <param name="maxFullMatrixCells">Sub-problems of at most this many cells are solved with the full-matrix core (clamped to
    /// <see cref="LinearSpaceCellThreshold"/>); smaller values recurse deeper.</param>
    public static AlignmentResult GlobalAlignLinearSpace(
        string sequence1, string sequence2, ScoringMatrix? scoring = null, long maxFullMatrixCells = LinearSpaceBaseCells) =>
        ScalarLinearSpace(sequence1, sequence2, scoring, affine: false, local: false, EndGaps.None, maxFullMatrixCells);

    /// <summary><see cref="GlobalAlignLinearSpace(string, string, ScoringMatrix?, long)"/> for DNA sequences.</summary>
    public static AlignmentResult GlobalAlignLinearSpace(
        DnaSequence sequence1, DnaSequence sequence2, ScoringMatrix? scoring = null, long maxFullMatrixCells = LinearSpaceBaseCells)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        return ScalarLinearSpace(sequence1.Sequence, sequence2.Sequence, scoring, affine: false, local: false, EndGaps.None, maxFullMatrixCells);
    }

    /// <summary>Linear-space global (Needleman-Wunsch, linear gaps) alignment with a substitution matrix and (possibly fractional) gap scores;
    /// identical to <see cref="GlobalAlign(string, string, MatrixScoring)"/>.</summary>
    public static AlignmentResult GlobalAlignLinearSpace(
        string sequence1, string sequence2, MatrixScoring scoring, long maxFullMatrixCells = LinearSpaceBaseCells) =>
        AlignWithMatrix(sequence1, sequence2, scoring, AlignmentType.Global, affine: false, linearSpaceBaseCells: Math.Max(1, maxFullMatrixCells));

    /// <summary><see cref="GlobalAlignLinearSpace(string, string, MatrixScoring, long)"/> for protein sequences.</summary>
    public static AlignmentResult GlobalAlignLinearSpace(
        ProteinSequence sequence1, ProteinSequence sequence2, MatrixScoring scoring, long maxFullMatrixCells = LinearSpaceBaseCells) =>
        AlignWithMatrix(sequence1, sequence2, scoring, AlignmentType.Global, affine: false, linearSpaceBaseCells: Math.Max(1, maxFullMatrixCells));

    /// <summary><see cref="GlobalAlignLinearSpace(string, string, MatrixScoring, long)"/> for DNA sequences.</summary>
    public static AlignmentResult GlobalAlignLinearSpace(
        DnaSequence sequence1, DnaSequence sequence2, MatrixScoring scoring, long maxFullMatrixCells = LinearSpaceBaseCells) =>
        AlignWithMatrix(sequence1, sequence2, scoring, AlignmentType.Global, affine: false, linearSpaceBaseCells: Math.Max(1, maxFullMatrixCells));

    /// <summary>
    /// Exact linear-space global affine-gap (Gotoh) alignment (O(m + n) memory instead of O(m·n); see the region comment of this file and
    /// <c>docs/algorithms/Alignment/Linear_Space_Alignment.md</c>). Score, aligned rows, coordinates and tie behaviour are
    /// identical to <see cref="GlobalAlignAffine(string, string, ScoringMatrix?)"/>, which switches to this program by itself above
    /// <see cref="LinearSpaceCellThreshold"/> cells; call this method to force it for any size.
    /// </summary>
    /// <param name="sequence1">First sequence (upper-cased internally).</param>
    /// <param name="sequence2">Second sequence.</param>
    /// <param name="scoring">Scoring (default <see cref="SimpleDna"/>).</param>
    /// <param name="maxFullMatrixCells">Sub-problems of at most this many cells are solved with the full-matrix core (clamped to
    /// <see cref="LinearSpaceCellThreshold"/>); smaller values recurse deeper.</param>
    public static AlignmentResult GlobalAlignAffineLinearSpace(
        string sequence1, string sequence2, ScoringMatrix? scoring = null, long maxFullMatrixCells = LinearSpaceBaseCells) =>
        ScalarLinearSpace(sequence1, sequence2, scoring, affine: true, local: false, EndGaps.None, maxFullMatrixCells);

    /// <summary><see cref="GlobalAlignAffineLinearSpace(string, string, ScoringMatrix?, long)"/> for DNA sequences.</summary>
    public static AlignmentResult GlobalAlignAffineLinearSpace(
        DnaSequence sequence1, DnaSequence sequence2, ScoringMatrix? scoring = null, long maxFullMatrixCells = LinearSpaceBaseCells)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        return ScalarLinearSpace(sequence1.Sequence, sequence2.Sequence, scoring, affine: true, local: false, EndGaps.None, maxFullMatrixCells);
    }

    /// <summary>Linear-space global affine-gap (Gotoh) alignment with a substitution matrix and (possibly fractional) gap scores;
    /// identical to <see cref="GlobalAlignAffine(string, string, MatrixScoring)"/>.</summary>
    public static AlignmentResult GlobalAlignAffineLinearSpace(
        string sequence1, string sequence2, MatrixScoring scoring, long maxFullMatrixCells = LinearSpaceBaseCells) =>
        AlignWithMatrix(sequence1, sequence2, scoring, AlignmentType.Global, affine: true, linearSpaceBaseCells: Math.Max(1, maxFullMatrixCells));

    /// <summary><see cref="GlobalAlignAffineLinearSpace(string, string, MatrixScoring, long)"/> for protein sequences.</summary>
    public static AlignmentResult GlobalAlignAffineLinearSpace(
        ProteinSequence sequence1, ProteinSequence sequence2, MatrixScoring scoring, long maxFullMatrixCells = LinearSpaceBaseCells) =>
        AlignWithMatrix(sequence1, sequence2, scoring, AlignmentType.Global, affine: true, linearSpaceBaseCells: Math.Max(1, maxFullMatrixCells));

    /// <summary><see cref="GlobalAlignAffineLinearSpace(string, string, MatrixScoring, long)"/> for DNA sequences.</summary>
    public static AlignmentResult GlobalAlignAffineLinearSpace(
        DnaSequence sequence1, DnaSequence sequence2, MatrixScoring scoring, long maxFullMatrixCells = LinearSpaceBaseCells) =>
        AlignWithMatrix(sequence1, sequence2, scoring, AlignmentType.Global, affine: true, linearSpaceBaseCells: Math.Max(1, maxFullMatrixCells));

    /// <summary>
    /// Exact linear-space local (Smith-Waterman, linear gaps) alignment (O(m + n) memory instead of O(m·n); see the region comment of this file and
    /// <c>docs/algorithms/Alignment/Linear_Space_Alignment.md</c>). Score, aligned rows, coordinates and tie behaviour are
    /// identical to <see cref="LocalAlign(string, string, ScoringMatrix?)"/>, which switches to this program by itself above
    /// <see cref="LinearSpaceCellThreshold"/> cells; call this method to force it for any size.
    /// </summary>
    /// <param name="sequence1">First sequence (upper-cased internally).</param>
    /// <param name="sequence2">Second sequence.</param>
    /// <param name="scoring">Scoring (default <see cref="SimpleDna"/>).</param>
    /// <param name="maxFullMatrixCells">Sub-problems of at most this many cells are solved with the full-matrix core (clamped to
    /// <see cref="LinearSpaceCellThreshold"/>); smaller values recurse deeper.</param>
    public static AlignmentResult LocalAlignLinearSpace(
        string sequence1, string sequence2, ScoringMatrix? scoring = null, long maxFullMatrixCells = LinearSpaceBaseCells) =>
        ScalarLinearSpace(sequence1, sequence2, scoring, affine: false, local: true, EndGaps.None, maxFullMatrixCells);

    /// <summary><see cref="LocalAlignLinearSpace(string, string, ScoringMatrix?, long)"/> for DNA sequences.</summary>
    public static AlignmentResult LocalAlignLinearSpace(
        DnaSequence sequence1, DnaSequence sequence2, ScoringMatrix? scoring = null, long maxFullMatrixCells = LinearSpaceBaseCells)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        return ScalarLinearSpace(sequence1.Sequence, sequence2.Sequence, scoring, affine: false, local: true, EndGaps.None, maxFullMatrixCells);
    }

    /// <summary>Linear-space local (Smith-Waterman, linear gaps) alignment with a substitution matrix and (possibly fractional) gap scores;
    /// identical to <see cref="LocalAlign(string, string, MatrixScoring)"/>.</summary>
    public static AlignmentResult LocalAlignLinearSpace(
        string sequence1, string sequence2, MatrixScoring scoring, long maxFullMatrixCells = LinearSpaceBaseCells) =>
        AlignWithMatrix(sequence1, sequence2, scoring, AlignmentType.Local, affine: false, linearSpaceBaseCells: Math.Max(1, maxFullMatrixCells));

    /// <summary><see cref="LocalAlignLinearSpace(string, string, MatrixScoring, long)"/> for protein sequences.</summary>
    public static AlignmentResult LocalAlignLinearSpace(
        ProteinSequence sequence1, ProteinSequence sequence2, MatrixScoring scoring, long maxFullMatrixCells = LinearSpaceBaseCells) =>
        AlignWithMatrix(sequence1, sequence2, scoring, AlignmentType.Local, affine: false, linearSpaceBaseCells: Math.Max(1, maxFullMatrixCells));

    /// <summary><see cref="LocalAlignLinearSpace(string, string, MatrixScoring, long)"/> for DNA sequences.</summary>
    public static AlignmentResult LocalAlignLinearSpace(
        DnaSequence sequence1, DnaSequence sequence2, MatrixScoring scoring, long maxFullMatrixCells = LinearSpaceBaseCells) =>
        AlignWithMatrix(sequence1, sequence2, scoring, AlignmentType.Local, affine: false, linearSpaceBaseCells: Math.Max(1, maxFullMatrixCells));

    /// <summary>
    /// Exact linear-space local affine-gap (Smith-Waterman-Gotoh) alignment (O(m + n) memory instead of O(m·n); see the region comment of this file and
    /// <c>docs/algorithms/Alignment/Linear_Space_Alignment.md</c>). Score, aligned rows, coordinates and tie behaviour are
    /// identical to <see cref="LocalAlignAffine(string, string, ScoringMatrix?)"/>, which switches to this program by itself above
    /// <see cref="LinearSpaceCellThreshold"/> cells; call this method to force it for any size.
    /// </summary>
    /// <param name="sequence1">First sequence (upper-cased internally).</param>
    /// <param name="sequence2">Second sequence.</param>
    /// <param name="scoring">Scoring (default <see cref="SimpleDna"/>).</param>
    /// <param name="maxFullMatrixCells">Sub-problems of at most this many cells are solved with the full-matrix core (clamped to
    /// <see cref="LinearSpaceCellThreshold"/>); smaller values recurse deeper.</param>
    public static AlignmentResult LocalAlignAffineLinearSpace(
        string sequence1, string sequence2, ScoringMatrix? scoring = null, long maxFullMatrixCells = LinearSpaceBaseCells) =>
        ScalarLinearSpace(sequence1, sequence2, scoring, affine: true, local: true, EndGaps.None, maxFullMatrixCells);

    /// <summary><see cref="LocalAlignAffineLinearSpace(string, string, ScoringMatrix?, long)"/> for DNA sequences.</summary>
    public static AlignmentResult LocalAlignAffineLinearSpace(
        DnaSequence sequence1, DnaSequence sequence2, ScoringMatrix? scoring = null, long maxFullMatrixCells = LinearSpaceBaseCells)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        return ScalarLinearSpace(sequence1.Sequence, sequence2.Sequence, scoring, affine: true, local: true, EndGaps.None, maxFullMatrixCells);
    }

    /// <summary>Linear-space local affine-gap (Smith-Waterman-Gotoh) alignment with a substitution matrix and (possibly fractional) gap scores;
    /// identical to <see cref="LocalAlignAffine(string, string, MatrixScoring)"/>.</summary>
    public static AlignmentResult LocalAlignAffineLinearSpace(
        string sequence1, string sequence2, MatrixScoring scoring, long maxFullMatrixCells = LinearSpaceBaseCells) =>
        AlignWithMatrix(sequence1, sequence2, scoring, AlignmentType.Local, affine: true, linearSpaceBaseCells: Math.Max(1, maxFullMatrixCells));

    /// <summary><see cref="LocalAlignAffineLinearSpace(string, string, MatrixScoring, long)"/> for protein sequences.</summary>
    public static AlignmentResult LocalAlignAffineLinearSpace(
        ProteinSequence sequence1, ProteinSequence sequence2, MatrixScoring scoring, long maxFullMatrixCells = LinearSpaceBaseCells) =>
        AlignWithMatrix(sequence1, sequence2, scoring, AlignmentType.Local, affine: true, linearSpaceBaseCells: Math.Max(1, maxFullMatrixCells));

    /// <summary><see cref="LocalAlignAffineLinearSpace(string, string, MatrixScoring, long)"/> for DNA sequences.</summary>
    public static AlignmentResult LocalAlignAffineLinearSpace(
        DnaSequence sequence1, DnaSequence sequence2, MatrixScoring scoring, long maxFullMatrixCells = LinearSpaceBaseCells) =>
        AlignWithMatrix(sequence1, sequence2, scoring, AlignmentType.Local, affine: true, linearSpaceBaseCells: Math.Max(1, maxFullMatrixCells));

    /// <summary>
    /// Exact linear-space semi-global (ends-free, default: fitting) linear-gap alignment (O(m + n) memory instead of O(m·n); see the region comment of this file and
    /// <c>docs/algorithms/Alignment/Linear_Space_Alignment.md</c>). Score, aligned rows, coordinates and tie behaviour are
    /// identical to <see cref="SemiGlobalAlign(string, string, ScoringMatrix?)"/>, which switches to this program by itself above
    /// <see cref="LinearSpaceCellThreshold"/> cells; call this method to force it for any size.
    /// </summary>
    /// <param name="sequence1">First sequence (upper-cased internally).</param>
    /// <param name="sequence2">Second sequence.</param>
    /// <param name="scoring">Scoring (default <see cref="SimpleDna"/>).</param>
    /// <param name="endGaps">Which end gaps are free (default <see cref="EndGaps.Fitting"/>); every one of the 16 policies is supported.</param>
    /// <param name="maxFullMatrixCells">Sub-problems of at most this many cells are solved with the full-matrix core (clamped to
    /// <see cref="LinearSpaceCellThreshold"/>); smaller values recurse deeper.</param>
    public static AlignmentResult SemiGlobalAlignLinearSpace(
        string sequence1, string sequence2, ScoringMatrix? scoring = null, EndGaps endGaps = EndGaps.Fitting, long maxFullMatrixCells = LinearSpaceBaseCells) =>
        ScalarLinearSpace(sequence1, sequence2, scoring, affine: false, local: false, endGaps, maxFullMatrixCells);

    /// <summary><see cref="SemiGlobalAlignLinearSpace(string, string, ScoringMatrix?, EndGaps, long)"/> for DNA sequences.</summary>
    public static AlignmentResult SemiGlobalAlignLinearSpace(
        DnaSequence sequence1, DnaSequence sequence2, ScoringMatrix? scoring = null, EndGaps endGaps = EndGaps.Fitting, long maxFullMatrixCells = LinearSpaceBaseCells)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        return ScalarLinearSpace(sequence1.Sequence, sequence2.Sequence, scoring, affine: false, local: false, endGaps, maxFullMatrixCells);
    }

    /// <summary>Linear-space semi-global (ends-free, default: fitting) linear-gap alignment with a substitution matrix and (possibly fractional) gap scores;
    /// identical to <see cref="SemiGlobalAlign(string, string, MatrixScoring, EndGaps)"/>.</summary>
    public static AlignmentResult SemiGlobalAlignLinearSpace(
        string sequence1, string sequence2, MatrixScoring scoring, EndGaps endGaps = EndGaps.Fitting, long maxFullMatrixCells = LinearSpaceBaseCells) =>
        AlignWithMatrix(sequence1, sequence2, scoring, AlignmentType.SemiGlobal, affine: false, endGaps, linearSpaceBaseCells: Math.Max(1, maxFullMatrixCells));

    /// <summary><see cref="SemiGlobalAlignLinearSpace(string, string, MatrixScoring, EndGaps, long)"/> for protein sequences.</summary>
    public static AlignmentResult SemiGlobalAlignLinearSpace(
        ProteinSequence sequence1, ProteinSequence sequence2, MatrixScoring scoring, EndGaps endGaps = EndGaps.Fitting, long maxFullMatrixCells = LinearSpaceBaseCells) =>
        AlignWithMatrix(sequence1, sequence2, scoring, AlignmentType.SemiGlobal, affine: false, endGaps, linearSpaceBaseCells: Math.Max(1, maxFullMatrixCells));

    /// <summary><see cref="SemiGlobalAlignLinearSpace(string, string, MatrixScoring, EndGaps, long)"/> for DNA sequences.</summary>
    public static AlignmentResult SemiGlobalAlignLinearSpace(
        DnaSequence sequence1, DnaSequence sequence2, MatrixScoring scoring, EndGaps endGaps = EndGaps.Fitting, long maxFullMatrixCells = LinearSpaceBaseCells) =>
        AlignWithMatrix(sequence1, sequence2, scoring, AlignmentType.SemiGlobal, affine: false, endGaps, linearSpaceBaseCells: Math.Max(1, maxFullMatrixCells));

    /// <summary>
    /// Exact linear-space semi-global (ends-free, default: fitting) affine-gap alignment (O(m + n) memory instead of O(m·n); see the region comment of this file and
    /// <c>docs/algorithms/Alignment/Linear_Space_Alignment.md</c>). Score, aligned rows, coordinates and tie behaviour are
    /// identical to <see cref="SemiGlobalAlignAffine(string, string, ScoringMatrix?)"/>, which switches to this program by itself above
    /// <see cref="LinearSpaceCellThreshold"/> cells; call this method to force it for any size.
    /// </summary>
    /// <param name="sequence1">First sequence (upper-cased internally).</param>
    /// <param name="sequence2">Second sequence.</param>
    /// <param name="scoring">Scoring (default <see cref="SimpleDna"/>).</param>
    /// <param name="endGaps">Which end gaps are free (default <see cref="EndGaps.Fitting"/>); every one of the 16 policies is supported.</param>
    /// <param name="maxFullMatrixCells">Sub-problems of at most this many cells are solved with the full-matrix core (clamped to
    /// <see cref="LinearSpaceCellThreshold"/>); smaller values recurse deeper.</param>
    public static AlignmentResult SemiGlobalAlignAffineLinearSpace(
        string sequence1, string sequence2, ScoringMatrix? scoring = null, EndGaps endGaps = EndGaps.Fitting, long maxFullMatrixCells = LinearSpaceBaseCells) =>
        ScalarLinearSpace(sequence1, sequence2, scoring, affine: true, local: false, endGaps, maxFullMatrixCells);

    /// <summary><see cref="SemiGlobalAlignAffineLinearSpace(string, string, ScoringMatrix?, EndGaps, long)"/> for DNA sequences.</summary>
    public static AlignmentResult SemiGlobalAlignAffineLinearSpace(
        DnaSequence sequence1, DnaSequence sequence2, ScoringMatrix? scoring = null, EndGaps endGaps = EndGaps.Fitting, long maxFullMatrixCells = LinearSpaceBaseCells)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        return ScalarLinearSpace(sequence1.Sequence, sequence2.Sequence, scoring, affine: true, local: false, endGaps, maxFullMatrixCells);
    }

    /// <summary>Linear-space semi-global (ends-free, default: fitting) affine-gap alignment with a substitution matrix and (possibly fractional) gap scores;
    /// identical to <see cref="SemiGlobalAlignAffine(string, string, MatrixScoring, EndGaps)"/>.</summary>
    public static AlignmentResult SemiGlobalAlignAffineLinearSpace(
        string sequence1, string sequence2, MatrixScoring scoring, EndGaps endGaps = EndGaps.Fitting, long maxFullMatrixCells = LinearSpaceBaseCells) =>
        AlignWithMatrix(sequence1, sequence2, scoring, AlignmentType.SemiGlobal, affine: true, endGaps, linearSpaceBaseCells: Math.Max(1, maxFullMatrixCells));

    /// <summary><see cref="SemiGlobalAlignAffineLinearSpace(string, string, MatrixScoring, EndGaps, long)"/> for protein sequences.</summary>
    public static AlignmentResult SemiGlobalAlignAffineLinearSpace(
        ProteinSequence sequence1, ProteinSequence sequence2, MatrixScoring scoring, EndGaps endGaps = EndGaps.Fitting, long maxFullMatrixCells = LinearSpaceBaseCells) =>
        AlignWithMatrix(sequence1, sequence2, scoring, AlignmentType.SemiGlobal, affine: true, endGaps, linearSpaceBaseCells: Math.Max(1, maxFullMatrixCells));

    /// <summary><see cref="SemiGlobalAlignAffineLinearSpace(string, string, MatrixScoring, EndGaps, long)"/> for DNA sequences.</summary>
    public static AlignmentResult SemiGlobalAlignAffineLinearSpace(
        DnaSequence sequence1, DnaSequence sequence2, MatrixScoring scoring, EndGaps endGaps = EndGaps.Fitting, long maxFullMatrixCells = LinearSpaceBaseCells) =>
        AlignWithMatrix(sequence1, sequence2, scoring, AlignmentType.SemiGlobal, affine: true, endGaps, linearSpaceBaseCells: Math.Max(1, maxFullMatrixCells));

    #endregion
}
