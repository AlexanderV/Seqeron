using Seqeron.Genomics.Core;

namespace Seqeron.Genomics.Analysis;

/// <summary>
/// PWM scoring beyond the forward-strand threshold scan: per-window scores (Biopython
/// <c>PositionSpecificScoringMatrix.calculate</c>), the both-strand search (Biopython
/// <c>search(sequence, threshold, both=True)</c>), count-matrix PWMs (Biopython
/// <c>counts.normalize(pseudocounts).log_odds(background)</c>, JASPAR pseudocounts) and the
/// discretised score distribution with FPR/FNR/balanced/patser thresholds (Biopython
/// <c>Bio.motifs.thresholds.ScoreDistribution</c>, N. Dojer 2008).
/// </summary>
public static partial class MotifFinder
{
    #region PWM window scoring (shared kernel)

    /// <summary>
    /// Score of the window starting at <paramref name="start"/>: Σ_j W[s[start+j], j]; NaN when the
    /// window contains a symbol other than upper-case A/C/G/T. The single scoring kernel of
    /// <see cref="ScanWithPwm"/>, <see cref="ScanWithPwmBothStrands"/> and <see cref="CalculatePwmScores(string, PositionWeightMatrix)"/>.
    /// </summary>
    internal static double ScorePwmWindow(string upperSequence, int start, PositionWeightMatrix pwm)
    {
        double score = 0;
        double[,] matrix = pwm.Matrix;
        for (int j = 0; j < pwm.Length; j++)
        {
            int baseIndex = AcgtIndex(upperSequence[start + j]);
            if (baseIndex < 0)
                return double.NaN;
            score += matrix[baseIndex, j];
        }

        return score;
    }

    /// <summary>
    /// PWM score of every window of <paramref name="sequence"/> (Biopython
    /// <c>PositionSpecificScoringMatrix.calculate(sequence)</c>): element i is the forward-strand score
    /// of the window starting at i, for i = 0 … n − m. Windows that contain a symbol other than
    /// A/C/G/T (e.g. N) score NaN, exactly as Biopython. Case-insensitive.
    /// </summary>
    /// <remarks>
    /// Biopython returns float32 values and a bare scalar when n = m; this method always returns a
    /// double array (length n − m + 1, or 0 when the sequence is shorter than the motif — where
    /// Biopython raises).
    /// </remarks>
    /// <exception cref="ArgumentNullException">An argument is null.</exception>
    public static double[] CalculatePwmScores(string sequence, PositionWeightMatrix pwm)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentNullException.ThrowIfNull(pwm);

        string seq = sequence.ToUpperInvariant();
        int windows = seq.Length - pwm.Length + 1;
        if (windows <= 0)
            return Array.Empty<double>();

        var scores = new double[windows];
        for (int i = 0; i < windows; i++)
            scores[i] = ScorePwmWindow(seq, i, pwm);
        return scores;
    }

    /// <inheritdoc cref="CalculatePwmScores(string, PositionWeightMatrix)"/>
    public static double[] CalculatePwmScores(DnaSequence sequence, PositionWeightMatrix pwm)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        return CalculatePwmScores(sequence.Sequence, pwm);
    }

    #endregion

    #region Both-strand PWM search

    /// <summary>
    /// Scans both strands with a PWM — Biopython
    /// <c>PositionSpecificScoringMatrix.search(sequence, threshold, both=True)</c>: the plus strand is
    /// scored with <paramref name="pwm"/> and the minus strand with
    /// <see cref="PositionWeightMatrix.ReverseComplement"/> over the same forward windows; every hit
    /// with score ≥ <paramref name="threshold"/> is reported.
    /// </summary>
    /// <remarks>
    /// <para>Coordinates: <see cref="PwmStrandMatch.Position"/> is the 0-based forward-strand start of the
    /// window on both strands; <see cref="PwmStrandMatch.BiopythonPosition"/> is Biopython's
    /// coordinate (plus: the start; minus: start − n, i.e. negative).</para>
    /// <para>Order: ascending <see cref="PwmStrandMatch.Position"/>; at equal positions (e.g. a
    /// palindromic PWM, which hits both strands of the same window) the plus-strand hit comes
    /// first. Biopython sorts by the same key with NumPy's default (unstable) <c>argsort</c>, so its
    /// order of such ties is implementation-defined; the set of hits is identical.</para>
    /// <para>Windows with a non-ACGT symbol are skipped (NaN score). Biopython scores in float32;
    /// scores here are double.</para>
    /// </remarks>
    /// <param name="sequence">DNA sequence to scan.</param>
    /// <param name="pwm">Position Weight Matrix (plus-strand orientation).</param>
    /// <param name="threshold">Minimum score (inclusive).</param>
    /// <returns>Strand-annotated hits; the minus-strand <see cref="PwmStrandMatch.MatchedSequence"/> is
    /// the site read 5'→3' on the minus strand (reverse complement of the forward window).</returns>
    /// <exception cref="ArgumentNullException">An argument is null.</exception>
    public static IEnumerable<PwmStrandMatch> ScanWithPwmBothStrands(
        DnaSequence sequence,
        PositionWeightMatrix pwm,
        double threshold = 0.0)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentNullException.ThrowIfNull(pwm);
        return ScanWithPwmBothStrandsCore(sequence, pwm, threshold);
    }

    private static IEnumerable<PwmStrandMatch> ScanWithPwmBothStrandsCore(
        DnaSequence sequence, PositionWeightMatrix pwm, double threshold)
    {
        string seq = sequence.Sequence;
        int n = seq.Length;
        int m = pwm.Length;
        PositionWeightMatrix rc = pwm.ReverseComplement();

        for (int i = 0; i <= n - m; i++)
        {
            double plus = ScorePwmWindow(seq, i, pwm);
            if (plus >= threshold)
                yield return new PwmStrandMatch(i, i, '+', seq.Substring(i, m), pwm.Consensus, plus);

            double minus = ScorePwmWindow(seq, i, rc);
            if (minus >= threshold)
            {
                yield return new PwmStrandMatch(i, i - n, '-',
                    DnaSequence.GetReverseComplementString(seq.Substring(i, m)), pwm.Consensus, minus);
            }
        }
    }

    #endregion

    #region Count-matrix PWMs (JASPAR / Biopython)

    /// <summary>
    /// JASPAR pseudocounts of a count matrix — Biopython <c>Bio.motifs.jaspar.calculate_pseudocounts</c>
    /// (Wasserman &amp; Sandelin 2004): p[b] = √(N̄) · q[b], with N̄ the mean column total and q the
    /// background (uniform when null), normalised to sum 1.
    /// </summary>
    /// <param name="counts">4 × L count matrix, rows A, C, G, T (L ≥ 1).</param>
    /// <param name="background">Optional background (A, C, G, T); uniform when null.</param>
    /// <returns>Pseudocounts for A, C, G, T.</returns>
    public static double[] JasparPseudocounts(double[,] counts, IReadOnlyList<double>? background = null)
    {
        ValidateCountMatrix(counts);
        int length = counts.GetLength(1);
        if (length == 0)
            throw new ArgumentException("Count matrix must have at least one column.", nameof(counts));
        double[] bg = ResolveBackground(background);

        double total = 0;
        for (int i = 0; i < length; i++)
            for (int b = 0; b < PwmAlphabetSize; b++)
                total += counts[b, i];
        double sqrtInstances = Math.Sqrt(total / length);

        var pseudo = new double[PwmAlphabetSize];
        for (int b = 0; b < PwmAlphabetSize; b++)
            pseudo[b] = sqrtInstances * bg[b];
        return pseudo;
    }

    /// <summary>Normalised background (A, C, G, T); the shared uniform array when null (read-only).</summary>
    internal static double[] ResolveBackground(IReadOnlyList<double>? background)
        => background is null ? UniformBackground : NormalizeBackground(background);

    internal static void ValidateCountMatrix(double[,] counts)
    {
        ArgumentNullException.ThrowIfNull(counts);
        if (counts.GetLength(0) != PwmAlphabetSize)
            throw new ArgumentException("Count matrix must have 4 rows (A, C, G, T).", nameof(counts));
        foreach (double c in counts)
        {
            if (!double.IsFinite(c) || c < 0)
                throw new ArgumentOutOfRangeException(nameof(counts), c, "Counts must be finite and non-negative.");
        }
    }

    internal static PositionWeightMatrix PwmFromCounts(
        double[,] counts, IReadOnlyList<double> pseudocounts, IReadOnlyList<double>? background)
    {
        ValidateCountMatrix(counts);
        ArgumentNullException.ThrowIfNull(pseudocounts);
        if (pseudocounts.Count != PwmAlphabetSize)
            throw new ArgumentException("Exactly 4 pseudocounts (A, C, G, T) are required.", nameof(pseudocounts));

        var pseudo = new double[PwmAlphabetSize];
        double pseudoSum = 0;
        for (int b = 0; b < PwmAlphabetSize; b++)
        {
            double p = pseudocounts[b];
            if (!double.IsFinite(p) || p < 0)
                throw new ArgumentOutOfRangeException(nameof(pseudocounts), p,
                    "Pseudocounts must be finite and non-negative.");
            pseudo[b] = p;
            pseudoSum += p;
        }

        double[] bg = ResolveBackground(background);

        int length = counts.GetLength(1);
        for (int i = 0; i < length; i++)
        {
            double columnTotal = pseudoSum;
            for (int b = 0; b < PwmAlphabetSize; b++)
                columnTotal += counts[b, i];
            if (columnTotal <= 0)
                throw new ArgumentException(
                    $"Column {i} has zero counts and zero pseudocounts; its frequencies are undefined.", nameof(counts));
        }

        return new PositionWeightMatrix(LogOddsFromCounts(counts, pseudo, pseudoSum, bg), length);
    }

    #endregion
}

/// <summary>
/// A strand-annotated PWM hit (<see cref="MotifFinder.ScanWithPwmBothStrands"/>).
/// </summary>
/// <param name="Position">0-based forward-strand start of the scored window (both strands).</param>
/// <param name="BiopythonPosition">Biopython <c>search(both=True)</c> coordinate: <paramref name="Position"/> on '+', <paramref name="Position"/> − n on '−'.</param>
/// <param name="Strand">'+' or '-'.</param>
/// <param name="MatchedSequence">The site read 5'→3' on its own strand.</param>
/// <param name="Pattern">Consensus of the (plus-strand) PWM.</param>
/// <param name="Score">Log-odds score.</param>
public readonly record struct PwmStrandMatch(
    int Position,
    int BiopythonPosition,
    char Strand,
    string MatchedSequence,
    string Pattern,
    double Score);

public sealed partial class PositionWeightMatrix
{
    /// <summary>
    /// Builds a log-odds PWM from a 4 × L count (position frequency) matrix with a scalar pseudocount
    /// added to every cell — Biopython <c>counts.normalize(pseudocounts=p).log_odds(background)</c>:
    /// W[b,j] = log2( ((c[b,j] + p) / (Σ_b c[b,j] + 4p)) / q[b] ). Column totals may differ (as in JASPAR).
    /// </summary>
    /// <param name="counts">Count matrix, rows A, C, G, T; finite, non-negative.</param>
    /// <param name="pseudocount">Pseudocount per cell (finite, ≥ 0; 0 = Biopython default <c>None</c>).</param>
    /// <param name="background">Background (A, C, G, T), normalised to sum 1; uniform when null.</param>
    public static PositionWeightMatrix FromCounts(double[,] counts, double pseudocount = 0.0, IReadOnlyList<double>? background = null)
        => MotifFinder.PwmFromCounts(counts, new[] { pseudocount, pseudocount, pseudocount, pseudocount }, background);

    /// <summary>
    /// Builds a log-odds PWM from a count matrix with per-base pseudocounts (Biopython
    /// <c>normalize(pseudocounts={A:…, C:…, G:…, T:…})</c>), e.g. <see cref="MotifFinder.JasparPseudocounts"/>.
    /// </summary>
    /// <param name="counts">Count matrix, rows A, C, G, T; finite, non-negative.</param>
    /// <param name="pseudocounts">Pseudocounts for A, C, G, T (finite, ≥ 0).</param>
    /// <param name="background">Background (A, C, G, T), normalised to sum 1; uniform when null.</param>
    public static PositionWeightMatrix FromCounts(double[,] counts, IReadOnlyList<double> pseudocounts, IReadOnlyList<double>? background = null)
        => MotifFinder.PwmFromCounts(counts, pseudocounts, background);

    /// <summary>
    /// Expected score of a random background window — Biopython <c>PositionSpecificScoringMatrix.mean(background)</c>:
    /// Σ_j Σ_b q[b]·2^W[b,j]·W[b,j] (cells that are NaN or −∞ are skipped, as in Biopython).
    /// For a PWM built against background q this is the relative entropy (information content) of the motif.
    /// </summary>
    /// <param name="background">Background (A, C, G, T); uniform when null.</param>
    public double Mean(IReadOnlyList<double>? background = null)
    {
        double[] bg = ResolveBackground(background);
        double sx = 0.0;
        for (int i = 0; i < Length; i++)
        {
            for (int b = 0; b < 4; b++)
            {
                double w = Matrix[b, i];
                if (double.IsNaN(w) || double.IsNegativeInfinity(w))
                    continue;
                double p = bg[b] * Math.Pow(2, w);
                sx += p * w;
            }
        }

        return sx;
    }

    /// <summary>
    /// Standard deviation of the score of a random window — Biopython
    /// <c>PositionSpecificScoringMatrix.std(background)</c> (per-column variances summed; NaN / −∞ cells skipped).
    /// </summary>
    /// <param name="background">Background (A, C, G, T); uniform when null.</param>
    public double Std(IReadOnlyList<double>? background = null)
    {
        double[] bg = ResolveBackground(background);
        double variance = 0.0;
        for (int i = 0; i < Length; i++)
        {
            double sx = 0.0, sxx = 0.0;
            for (int b = 0; b < 4; b++)
            {
                double w = Matrix[b, i];
                if (double.IsNaN(w) || double.IsNegativeInfinity(w))
                    continue;
                double p = bg[b] * Math.Pow(2, w);
                sx += p * w;
                sxx += p * w * w;
            }
            sxx -= sx * sx;
            variance += sxx;
        }

        variance = Math.Max(variance, 0);
        return Math.Sqrt(variance);
    }

    /// <summary>
    /// Discretised score distribution of this PWM under the motif and background models — Biopython
    /// <c>pssm.distribution(background, precision)</c> (<c>Bio.motifs.thresholds.ScoreDistribution</c>).
    /// </summary>
    /// <param name="background">Background (A, C, G, T); uniform when null.</param>
    /// <param name="precision">Grid points per motif position (Biopython default 10³).</param>
    /// <exception cref="InvalidOperationException">The matrix has a non-finite cell (the grid is undefined).</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="precision"/> × length &lt; 2.</exception>
    public PwmScoreDistribution ScoreDistribution(IReadOnlyList<double>? background = null, int precision = 1000)
        => new(this, ResolveBackground(background), precision);

    private static double[] ResolveBackground(IReadOnlyList<double>? background)
        => MotifFinder.ResolveBackground(background);
}

/// <summary>
/// Approximate (grid-discretised) distribution of PWM scores under the motif model and the background
/// model, and the score thresholds derived from it — a line-by-line port of Biopython
/// <c>Bio.motifs.thresholds.ScoreDistribution</c> (N. Dojer 2008, adapted by B. Wilczynski), including
/// Python's float floor-division when mapping scores to grid indices, so thresholds are reproduced exactly.
/// </summary>
/// <remarks>
/// Grid: minScore = min(0, pssm.min), interval = max(0, pssm.max) − minScore, n = precision · L points,
/// step = interval / (n − 1). Both densities start as a point mass at score 0 and are convolved column by
/// column (motif probability of base b at column j = q[b]·2^W[b,j]; background probability = q[b]); shifts
/// that leave the grid are clamped to its ends.
/// </remarks>
public sealed class PwmScoreDistribution
{
    private readonly double[] _motifDensity;
    private readonly double[] _backgroundDensity;

    internal PwmScoreDistribution(PositionWeightMatrix pwm, double[] background, int precision)
    {
        foreach (double w in pwm.Matrix)
        {
            if (!double.IsFinite(w))
                throw new InvalidOperationException(
                    "The score distribution needs a finite matrix (use a positive pseudocount).");
        }
        if (precision < 1 || (long)precision * pwm.Length < 2)
            throw new ArgumentOutOfRangeException(nameof(precision), precision,
                "precision × motif length must be at least 2.");

        MinScore = Math.Min(0.0, pwm.MinScore);
        double interval = Math.Max(0.0, pwm.MaxScore) - MinScore;
        PointCount = checked(precision * pwm.Length);
        Step = interval / (PointCount - 1);
        MeanScore = pwm.Mean(background);

        var mo = new double[PointCount];
        var bgd = new double[PointCount];
        int origin = -IndexDiff(MinScore);
        mo[origin] = 1.0;
        bgd[origin] = 1.0;

        for (int position = 0; position < pwm.Length; position++)
        {
            var moNew = new double[PointCount];
            var bgNew = new double[PointCount];
            for (int b = 0; b < 4; b++)
            {
                double score = pwm.Matrix[b, position];
                double q = background[b];
                double moProb = Math.Pow(2, score) * q;
                int d = IndexDiff(score);
                for (int i = 0; i < PointCount; i++)
                {
                    int k = Math.Max(0, Math.Min(PointCount - 1, i + d));
                    moNew[k] += mo[i] * moProb;
                    bgNew[k] += bgd[i] * q;
                }
            }
            mo = moNew;
            bgd = bgNew;
        }

        _motifDensity = mo;
        _backgroundDensity = bgd;
    }

    /// <summary>Score of grid point 0: min(0, PWM minimum).</summary>
    public double MinScore { get; }

    /// <summary>Grid spacing.</summary>
    public double Step { get; }

    /// <summary>Number of grid points (precision · motif length).</summary>
    public int PointCount { get; }

    /// <summary>Expected background score <see cref="PositionWeightMatrix.Mean"/> (Biopython <c>ic</c>), used by <see cref="ThresholdPatser"/>.</summary>
    public double MeanScore { get; }

    /// <summary>Probability mass of each grid point under the motif model.</summary>
    public IReadOnlyList<double> MotifDensity => _motifDensity;

    /// <summary>Probability mass of each grid point under the background model.</summary>
    public IReadOnlyList<double> BackgroundDensity => _backgroundDensity;

    /// <summary>
    /// Score threshold with background false-positive rate ≈ <paramref name="fpr"/> (Biopython <c>threshold_fpr</c>):
    /// the lowest grid score whose upper-tail background mass reaches <paramref name="fpr"/>.
    /// </summary>
    /// <param name="fpr">Type-I error rate in [0, 1].</param>
    public double ThresholdFpr(double fpr)
    {
        ValidateRate(fpr, nameof(fpr));
        int i = PointCount;
        double prob = 0.0;
        while (prob < fpr && i > 0)
        {
            i--;
            prob += _backgroundDensity[i];
        }
        return MinScore + i * Step;
    }

    /// <summary>
    /// Score threshold with motif false-negative rate ≈ <paramref name="fnr"/> (Biopython <c>threshold_fnr</c>).
    /// </summary>
    /// <param name="fnr">Type-II error rate in [0, 1].</param>
    public double ThresholdFnr(double fnr)
    {
        ValidateRate(fnr, nameof(fnr));
        int i = -1;
        double prob = 0.0;
        while (prob < fnr && i < PointCount - 1)
        {
            i++;
            prob += _motifDensity[i];
        }
        return MinScore + i * Step;
    }

    /// <summary>
    /// Threshold where FNR ≈ FPR · <paramref name="rateProportion"/> (Biopython <c>threshold_balanced</c>).
    /// </summary>
    public double ThresholdBalanced(double rateProportion = 1.0) => ThresholdBalanced(rateProportion, out _);

    /// <summary>
    /// Threshold where FNR ≈ FPR · <paramref name="rateProportion"/>, also returning the FPR reached
    /// (Biopython <c>threshold_balanced(rate_proportion, return_rate=True)</c>).
    /// </summary>
    public double ThresholdBalanced(double rateProportion, out double falsePositiveRate)
    {
        if (!double.IsFinite(rateProportion) || rateProportion < 0)
            throw new ArgumentOutOfRangeException(nameof(rateProportion), rateProportion,
                "Rate proportion must be finite and non-negative.");
        int i = PointCount;
        double fpr = 0.0;
        double fnr = 1.0;
        while (fpr * rateProportion < fnr && i > 0)
        {
            i--;
            fpr += _backgroundDensity[i];
            fnr -= _motifDensity[i];
        }
        falsePositiveRate = fpr;
        return MinScore + i * Step;
    }

    /// <summary>
    /// patser-style threshold (Hertz &amp; Stormo 1999; Biopython <c>threshold_patser</c>):
    /// <see cref="ThresholdFpr"/>(2^−<see cref="MeanScore"/>), i.e. log2(FPR) = −information content.
    /// </summary>
    public double ThresholdPatser() => ThresholdFpr(Math.Pow(2, -MeanScore));

    private static void ValidateRate(double rate, string name)
    {
        if (!(rate >= 0 && rate <= 1))
            throw new ArgumentOutOfRangeException(name, rate, "Rate must be in [0, 1].");
    }

    /// <summary>Biopython <c>_index_diff(x)</c> = int((x + 0.5·step) // step) with Python float floor-division.</summary>
    private int IndexDiff(double x) => (int)PythonFloorDiv(x - 0.0 + 0.5 * Step, Step);

    /// <summary>CPython <c>float_floor_div</c> (Objects/floatobject.c, <c>_float_div_mod</c>).</summary>
    internal static double PythonFloorDiv(double vx, double wx)
    {
        double mod = vx % wx; // C fmod semantics
        double div = (vx - mod) / wx;
        if (mod != 0.0 && (wx < 0) != (mod < 0))
            div -= 1.0;

        if (div != 0.0)
        {
            double floorDiv = Math.Floor(div);
            if (div - floorDiv > 0.5)
                floorDiv += 1.0;
            return floorDiv;
        }

        return Math.CopySign(0.0, vx / wx);
    }
}
