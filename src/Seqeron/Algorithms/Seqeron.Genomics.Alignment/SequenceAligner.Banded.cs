using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Alignment;

/// <summary>The settings of <see cref="SequenceAligner.BandedAlign(string, string, BandedAlignOptions)"/>.</summary>
public sealed record BandedAlignOptions
{
    /// <summary>Global (default), local or semi-global (<see cref="EndGaps"/>).</summary>
    public AlignmentType Type { get; init; } = AlignmentType.Global;

    /// <summary>Affine gaps (Gotoh) instead of linear gaps.</summary>
    public bool Affine { get; init; }

    /// <summary>The free end gaps of a semi-global alignment (the type is then <see cref="AlignmentType.SemiGlobal"/>).</summary>
    public EndGaps EndGaps { get; init; }

    /// <summary>Match / mismatch / gap scoring (default <see cref="SequenceAligner.SimpleDna"/>); ignored when <see cref="MatrixScoring"/> is set.</summary>
    public ScoringMatrix? Scoring { get; init; }

    /// <summary>A substitution matrix scoring (fractional gap scores allowed).</summary>
    public MatrixScoring? MatrixScoring { get; init; }

    /// <summary>The lowest diagonal j - i of the band (a cell (i, j) is inside when <c>DiagonalLow &lt;= j - i &lt;= DiagonalHigh</c>).</summary>
    public int DiagonalLow { get; init; }

    /// <summary>The highest diagonal of the band.</summary>
    public int DiagonalHigh { get; init; }

    /// <summary>
    /// X-drop pruning (BLAST gapped extension): a cell whose best state scores more than X below the best score seen so far is dropped. A heuristic: the result equals the banded optimum when X is at least the largest
    /// fall of the score along the optimal path (and always when null).
    /// </summary>
    public int? XDrop { get; init; }
}

/// <summary>The result of <see cref="SequenceAligner.GlobalAlignBandedExact(string, string, ScoringMatrix?, bool, int)"/>: the alignment, the band that proved it optimal and the work done.</summary>
/// <param name="Alignment">The alignment (the rows of the full dynamic program when <paramref name="Proven"/>).</param>
/// <param name="BandWidth">The half-width w of the final band (diagonals from min(0, n - m) - w to max(0, n - m) + w).</param>
/// <param name="Proven">True when every path that leaves the band scores strictly less than the banded optimum (Ukkonen's bound), or the band is the whole matrix.</param>
/// <param name="Rounds">Number of bandings computed (the band doubles each round).</param>
/// <param name="Cells">Total number of band cells filled over all rounds.</param>
public sealed record BandedAlignment(AlignmentResult Alignment, int BandWidth, bool Proven, int Rounds, long Cells);

public static partial class SequenceAligner
{
    private const long BandedMaxCells = 64_000_000;

    // the band of the three Gotoh matrices: row i holds the cells j = i + DiagonalLow .. i + DiagonalHigh
    private sealed class BandStates : IAffineStates
    {
        private readonly int _n;
        private readonly int _low;
        private readonly int _width;

        public BandStates(int m, int n, int low, int high)
        {
            _n = n;
            _low = low;
            _width = high - low + 1;
            Mm = new int[m + 1][];
            Xx = new int[m + 1][];
            Yy = new int[m + 1][];
            for (int i = 0; i <= m; i++)
            {
                Mm[i] = new int[_width];
                Xx[i] = new int[_width];
                Yy[i] = new int[_width];
                Array.Fill(Mm[i], AffineNegInf);
                Array.Fill(Xx[i], AffineNegInf);
                Array.Fill(Yy[i], AffineNegInf);
            }
        }

        public int[][] Mm { get; }

        public int[][] Xx { get; }

        public int[][] Yy { get; }

        public bool Has(int i, int j)
        {
            int k = j - i - _low;
            return j >= 0 && j <= _n && k >= 0 && k < _width;
        }

        public int Index(int i, int j) => j - i - _low;

        public int M(int i, int j) => Has(i, j) ? Mm[i][j - i - _low] : AffineNegInf;

        public int X(int i, int j) => Has(i, j) ? Xx[i][j - i - _low] : AffineNegInf;

        public int Y(int i, int j) => Has(i, j) ? Yy[i][j - i - _low] : AffineNegInf;
    }

    private readonly struct BandView(BandStates band) : IAffineStates
    {
        public int M(int i, int j) => band.M(i, j);

        public int X(int i, int j) => band.X(i, j);

        public int Y(int i, int j) => band.Y(i, j);
    }

    private static AlignmentResult BandedCore(string a, string b, PairScoring scoring, AlignmentType type, EndGaps ends, int diagonalLow, int diagonalHigh, int? xDrop, out long cells)
    {
        bool local = type == AlignmentType.Local;
        if (local)
            ends = EndGaps.None;
        else if (ends != EndGaps.None && type == AlignmentType.Global)
            type = AlignmentType.SemiGlobal;
        int m = a.Length, n = b.Length;
        scoring.EnsureRange(m, n);
        if (diagonalLow > diagonalHigh)
            throw new ArgumentException("DiagonalLow must not exceed DiagonalHigh.");
        if (!local && (diagonalLow > 0 || diagonalHigh < 0))
            throw new ArgumentException("A global or semi-global band must contain the diagonal 0 (the origin).");
        if (type == AlignmentType.Global && (n - m < diagonalLow || n - m > diagonalHigh))
            throw new ArgumentException($"The band [{diagonalLow}, {diagonalHigh}] does not contain the end cell (diagonal {n - m}).");
        // a band wider than the matrix is the matrix
        diagonalLow = Math.Max(diagonalLow, -m);
        diagonalHigh = Math.Min(diagonalHigh, n);
        long width = (long)diagonalHigh - diagonalLow + 1;
        if ((m + 1L) * width > BandedMaxCells)
            throw new ArgumentException($"The band needs {(m + 1L) * width} cells per matrix; the limit is {BandedMaxCells}: narrow the band.");
        bool freeStart1 = ends.HasFlag(EndGaps.FreeStart1), freeStart2 = ends.HasFlag(EndGaps.FreeStart2);
        bool freeEnd1 = ends.HasFlag(EndGaps.FreeEnd1), freeEnd2 = ends.HasFlag(EndGaps.FreeEnd2);
        var band = new BandStates(m, n, diagonalLow, diagonalHigh);
        cells = 0;

        // borders (AffineInitBorders restricted to the band)
        int e = scoring.XExt(0, n), oe = scoring.XOpen(0, n) + e;
        int ey = scoring.YExt(0, m), oey = scoring.YOpen(0, m) + ey;
        if (!local)
            band.Mm[0][band.Index(0, 0)] = 0;
        for (int i = 1; i <= m && band.Has(i, 0); i++)
        {
            int k = band.Index(i, 0);
            band.Mm[i][k] = freeStart1 ? 0 : AffineNegInf;
            band.Xx[i][k] = local || freeStart1 ? AffineNegInf : Math.Max(AffineNegInf, Max3(band.M(i - 1, 0) + oe, band.X(i - 1, 0) + e, band.Y(i - 1, 0) + oe));
        }

        for (int j = 1; j <= n && band.Has(0, j); j++)
        {
            int k = band.Index(0, j);
            band.Mm[0][k] = freeStart2 ? 0 : AffineNegInf;
            band.Yy[0][k] = local || freeStart2 ? AffineNegInf : Math.Max(AffineNegInf, Max3(band.M(0, j - 1) + oey, band.Y(0, j - 1) + ey, band.X(0, j - 1) + oey));
        }

        int startFloor = local ? 0 : AffineNegInf;
        int bestLocal = 0, bestI = 0, bestJ = 0;
        int bestSoFar = int.MinValue;
        for (int i = 1; i <= m; i++)
        {
            AlignmentCancellation.ThrowIfRequested();
            int jLo = Math.Max(1, i + diagonalLow), jHi = Math.Min(n, i + diagonalHigh);
            int ye = local ? scoring.YInteriorExt : scoring.YExt(i, m);
            int yoe = (local ? scoring.YInteriorOpen : scoring.YOpen(i, m)) + ye;
            char c1 = a[i - 1];
            for (int j = jLo; j <= jHi; j++)
            {
                int xe = local ? scoring.XInteriorExt : scoring.XExt(j, n);
                int xoe = (local ? scoring.XInteriorOpen : scoring.XOpen(j, n)) + xe;
                int sub = scoring.Sub(c1, b[j - 1]);
                int mv = Math.Max(startFloor, Max3(band.M(i - 1, j - 1), band.X(i - 1, j - 1), band.Y(i - 1, j - 1))) + sub;
                int xv = Max3(band.M(i - 1, j) + xoe, band.X(i - 1, j) + xe, band.Y(i - 1, j) + xoe);
                int yv = Max3(band.M(i, j - 1) + yoe, band.Y(i, j - 1) + ye, band.X(i, j - 1) + yoe);
                mv = Math.Max(mv, AffineNegInf);
                xv = Math.Max(xv, AffineNegInf);
                yv = Math.Max(yv, AffineNegInf);
                if (xDrop is { } drop)
                {
                    int top = Max3(mv, xv, yv);
                    if (bestSoFar != int.MinValue && top < bestSoFar - drop)
                    {
                        mv = xv = yv = AffineNegInf;
                    }
                    else if (top > bestSoFar)
                    {
                        bestSoFar = top;
                    }
                }

                int k = band.Index(i, j);
                band.Mm[i][k] = mv;
                band.Xx[i][k] = xv;
                band.Yy[i][k] = yv;
                cells++;
                if (local && mv > bestLocal)
                {
                    bestLocal = mv;
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
            ci = m;
            cj = n;
            if (freeEnd2 || freeEnd1)
            {
                int bestValue = int.MinValue;
                bool found = false;
                foreach (var (ei, ej) in EndCellCandidates(m, n, ends))
                {
                    if (!band.Has(ei, ej))
                        continue;
                    int v = Max3(band.M(ei, ej), band.X(ei, ej), band.Y(ei, ej));
                    if (v > bestValue)
                    {
                        bestValue = v;
                        ci = ei;
                        cj = ej;
                        found = true;
                    }
                }

                if (!found)
                    throw new ArgumentException("The band does not contain any end cell.");
            }

            best = Max3(band.M(ci, cj), band.X(ci, cj), band.Y(ci, cj));
            if (best <= AffineNegInf / 2)
                throw new ArgumentException("No alignment lies inside the band (the optimum needs cells outside it, or X-drop pruned every path).");
            state = PickState(best, band.M(ci, cj), StateMatch, band.X(ci, cj), StateGapInSeq2, StateGapInSeq1);
        }

        return AffineTraceback(a, b, new BandView(band), scoring, local, ends, ci, cj, state, best);
    }

    private static AlignmentResult BandedFinish(AlignmentResult r, PairScoring ps)
    {
        if (ps.Scale == 1)
            return r;
        double exact = (double)r.Score / ps.Scale;
        return r with { Score = (int)Math.Round(exact, MidpointRounding.AwayFromZero), ExactScore = exact };
    }

    private static PairScoring BandedScoring(BandedAlignOptions o, string a, string b)
    {
        if (o.MatrixScoring is { } ms)
            return PairScoring.Create(o.Affine ? ms : ms with { GapOpen = 0 }, a, b);
        var s = o.Scoring ?? SimpleDna;
        return o.Affine ? s : s with { GapOpen = 0 };
    }

    /// <summary>
    /// Banded alignment (Ukkonen 1985; Chao, Pearson &amp; Miller 1992): only the cells of the diagonals <see cref="BandedAlignOptions.DiagonalLow"/> .. <see cref="BandedAlignOptions.DiagonalHigh"/> are computed, O(m · band) time and memory
    /// instead of O(m · n), for global, local and semi-global alignment with linear or affine gaps and scalar or matrix scoring. The result is the alignment of the full dynamic program (same score, rows, coordinates, ties) whenever the
    /// full program's alignment lies inside the band; otherwise the best alignment inside the band (never better than the full optimum). <see cref="BandedAlignOptions.XDrop"/> adds BLAST-style pruning.
    /// </summary>
    /// <exception cref="ArgumentException">The band does not contain the origin / the end cell, is too wide for memory, or no alignment lies inside it.</exception>
    public static AlignmentResult BandedAlign(string sequence1, string sequence2, BandedAlignOptions options)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        ArgumentNullException.ThrowIfNull(options);
        string a = sequence1.ToUpperInvariant(), b = sequence2.ToUpperInvariant();
        var ps = BandedScoring(options, a, b);
        return BandedFinish(BandedCore(a, b, ps, options.Type, options.EndGaps, options.DiagonalLow, options.DiagonalHigh, options.XDrop, out _), ps);
    }

    private static BandedAlignOptions GlobalBand(string a, string b, ScoringMatrix? scoring, bool affine, int bandWidth)
    {
        if (bandWidth < 0)
            throw new ArgumentOutOfRangeException(nameof(bandWidth), bandWidth, "The band width must not be negative.");
        int d = b.Length - a.Length;
        return new BandedAlignOptions { Type = AlignmentType.Global, Affine = affine, Scoring = scoring, DiagonalLow = Math.Min(0, d) - bandWidth, DiagonalHigh = Math.Max(0, d) + bandWidth };
    }

    /// <summary>Banded global alignment with linear gaps: the diagonals from min(0, n - m) - <paramref name="bandWidth"/> to max(0, n - m) + <paramref name="bandWidth"/> (the rows of <see cref="GlobalAlign(string, string, ScoringMatrix?)"/> when its alignment lies inside).</summary>
    public static AlignmentResult GlobalAlignBanded(string sequence1, string sequence2, ScoringMatrix? scoring, int bandWidth)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        return BandedAlign(sequence1, sequence2, GlobalBand(sequence1, sequence2, scoring, false, bandWidth));
    }

    /// <summary>Banded global alignment with affine gaps (see <see cref="GlobalAlignBanded(string, string, ScoringMatrix?, int)"/>; the rows of <see cref="GlobalAlignAffine(string, string, ScoringMatrix?)"/> when inside).</summary>
    public static AlignmentResult GlobalAlignAffineBanded(string sequence1, string sequence2, ScoringMatrix? scoring, int bandWidth)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        return BandedAlign(sequence1, sequence2, GlobalBand(sequence1, sequence2, scoring, true, bandWidth));
    }

    /// <summary>Banded local alignment (Smith-Waterman) with linear gaps inside the diagonals <paramref name="diagonalLow"/> .. <paramref name="diagonalHigh"/> (j - i): the alignment around a seed diagonal.</summary>
    public static AlignmentResult LocalAlignBanded(string sequence1, string sequence2, ScoringMatrix? scoring, int diagonalLow, int diagonalHigh) =>
        BandedAlign(sequence1, sequence2, new BandedAlignOptions { Type = AlignmentType.Local, Scoring = scoring, DiagonalLow = diagonalLow, DiagonalHigh = diagonalHigh });

    /// <summary>Banded local alignment with affine gaps (see <see cref="LocalAlignBanded(string, string, ScoringMatrix?, int, int)"/>).</summary>
    public static AlignmentResult LocalAlignAffineBanded(string sequence1, string sequence2, ScoringMatrix? scoring, int diagonalLow, int diagonalHigh) =>
        BandedAlign(sequence1, sequence2, new BandedAlignOptions { Type = AlignmentType.Local, Affine = true, Scoring = scoring, DiagonalLow = diagonalLow, DiagonalHigh = diagonalHigh });

    /// <summary>Banded semi-global alignment with linear gaps: <paramref name="endGaps"/> free end gaps, the diagonals <paramref name="diagonalLow"/> .. <paramref name="diagonalHigh"/> (the band must contain the diagonal 0).</summary>
    public static AlignmentResult SemiGlobalAlignBanded(string sequence1, string sequence2, ScoringMatrix? scoring, EndGaps endGaps, int diagonalLow, int diagonalHigh) =>
        BandedAlign(sequence1, sequence2, new BandedAlignOptions { Type = AlignmentType.SemiGlobal, EndGaps = endGaps, Scoring = scoring, DiagonalLow = diagonalLow, DiagonalHigh = diagonalHigh });

    /// <summary>Banded semi-global alignment with affine gaps (see <see cref="SemiGlobalAlignBanded(string, string, ScoringMatrix?, EndGaps, int, int)"/>).</summary>
    public static AlignmentResult SemiGlobalAlignAffineBanded(string sequence1, string sequence2, ScoringMatrix? scoring, EndGaps endGaps, int diagonalLow, int diagonalHigh) =>
        BandedAlign(sequence1, sequence2, new BandedAlignOptions { Type = AlignmentType.SemiGlobal, Affine = true, EndGaps = endGaps, Scoring = scoring, DiagonalLow = diagonalLow, DiagonalHigh = diagonalHigh });

    /// <summary>
    /// Global alignment that is exact AND banded (Ukkonen's doubling): the band starts at <paramref name="initialBandWidth"/> and doubles until the banded score is strictly above the upper bound of every path that leaves it
    /// (such a path has at least |n - m| + 2(w + 1) gap columns, so it scores at most ((m + n - G) / 2) * maxPair + G * gapExtend with G that number), or the band is the whole matrix. The result then equals
    /// <see cref="GlobalAlign(string, string, ScoringMatrix?)"/> / <see cref="GlobalAlignAffine(string, string, ScoringMatrix?)"/> row for row, in O(s · min(m, n)) time for an alignment of s gap columns.
    /// </summary>
    [System.Runtime.CompilerServices.OverloadResolutionPriority(1)]
    public static BandedAlignment GlobalAlignBandedExact(string sequence1, string sequence2, ScoringMatrix? scoring = null, bool affine = false, int initialBandWidth = 16)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        return BandedExact(sequence1, sequence2, scoring, null, affine, initialBandWidth);
    }

    /// <summary><see cref="GlobalAlignBandedExact(string, string, ScoringMatrix?, bool, int)"/> with a substitution matrix.</summary>
    public static BandedAlignment GlobalAlignBandedExact(string sequence1, string sequence2, MatrixScoring scoring, bool affine = false, int initialBandWidth = 16)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        ArgumentNullException.ThrowIfNull(scoring);
        return BandedExact(sequence1, sequence2, null, scoring, affine, initialBandWidth);
    }

    private static BandedAlignment BandedExact(string sequence1, string sequence2, ScoringMatrix? scalar, MatrixScoring? matrix, bool affine, int initialBandWidth)
    {
        if (initialBandWidth < 1)
            throw new ArgumentOutOfRangeException(nameof(initialBandWidth), initialBandWidth, "The initial band width must be at least 1.");
        string a = sequence1.ToUpperInvariant(), b = sequence2.ToUpperInvariant();
        int m = a.Length, n = b.Length, d = n - m;
        var options = new BandedAlignOptions { Affine = affine, Scoring = scalar, MatrixScoring = matrix };
        var ps = BandedScoring(options, a, b);
        int maxPair = int.MinValue;
        foreach (char x in a.Distinct())
            foreach (char y in b.Distinct())
                maxPair = Math.Max(maxPair, ps.Sub(x, y));
        int ext = ps.GapExtend, open = ps.GapOpen;
        int full = Math.Max(m, n) + 1;
        int w = Math.Min(initialBandWidth, full);
        long cells = 0;
        int rounds = 0;
        while (true)
        {
            rounds++;
            int low = Math.Min(0, d) - w, high = Math.Max(0, d) + w;
            var r = BandedCore(a, b, ps, AlignmentType.Global, EndGaps.None, low, high, null, out long c);
            cells += c;
            bool whole = low <= -m && high >= n;
            bool proven = whole;
            if (!whole && a.Length > 0 && b.Length > 0 && maxPair > 0 && ext <= 0 && open <= 0)
            {
                long gmin = Math.Abs((long)d) + 2L * (w + 1);
                double bound = (m + n - gmin) / 2.0 * maxPair + (double)gmin * ext;
                proven = r.Score > bound;
            }

            if (proven || whole)
                return new BandedAlignment(BandedFinish(r, ps), w, proven, rounds, cells);
            w = (int)Math.Min((long)w * 2, full);
        }
    }
}
