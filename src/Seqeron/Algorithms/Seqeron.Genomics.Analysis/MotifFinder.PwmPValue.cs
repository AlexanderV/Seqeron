namespace Seqeron.Genomics.Analysis;

/// <summary>
/// Exact PWM score p-values and p-value thresholds under an i.i.d. background — the successive
/// refinement of integer-rounded matrices of Touzet &amp; Varré (2007), Algorithms Mol Biol 2:15
/// (TFM-Pvalue; reference C++ source: CRAN package TFMPvalue, <c>src/Matrix.cpp</c> / <c>src/TFMpvalue.cpp</c>).
/// </summary>
public static partial class MotifFinder
{
    /// <summary>
    /// Exact p-value of a PWM score: P(S ≥ <paramref name="score"/>) for a random word of length L drawn
    /// i.i.d. from <paramref name="background"/>, where S is the window score exactly as
    /// <see cref="CalculatePwmScores(string, PositionWeightMatrix)"/> computes it (left-to-right double sum
    /// of W[b, j]) — so the p-value of an observed window's score always counts that window.
    /// </summary>
    /// <remarks>
    /// <para><b>Method</b> (Touzet &amp; Varré 2007, TFM-Pvalue). For a scale g = 10, 100, 1000, … the matrix is
    /// rounded down, M[b,j] = ⌊g·W[b,j]⌋, so for every word g·S − M ∈ [0, E] with
    /// E = Σ_j max_b (g·W[b,j] − M[b,j]) (the paper's error bound). A dynamic programme over the integer
    /// partial scores (sparse maps, with the look-ahead pruning of TFM-Pvalue's <c>fastPvalue</c>: a prefix
    /// whose every completion is above / below the threshold is settled at once) gives the mass of words
    /// certainly above the threshold (M ≥ ⌈gα⌉ + 2) and of the undecided band (gα − E − 2 ≤ M &lt; ⌈gα⌉ + 2;
    /// the ±2 absorbs floating-point rounding). TFM-Pvalue refines g until the band is empty; here the band
    /// words are also enumerated exactly — a layered DP over (integer partial score, exact double partial window
    /// score) states, pruned by the exact sets of reachable integer suffix scores of the last columns
    /// (meet-in-the-middle) — and their window scores compared with α as soon as the band fits the budget. This
    /// also resolves scores that some word attains exactly (TFM-Pvalue does not converge there and returns 0 for
    /// the consensus score). The result is exact (only the summation of word probabilities is rounded).</para>
    /// <para>If the band can still not be resolved at the finest scale allowed (g·Σ_j max_b |W| ≤ 2⁴⁰) or the
    /// work budget is exhausted, the certified bounds are returned with <see cref="PwmPValueResult.IsExact"/>
    /// = false and <see cref="PwmPValueResult.PValue"/> = the upper (conservative) bound.</para>
    /// <para>Unlike <see cref="PwmScoreDistribution"/> (Biopython's fixed-grid approximation), no grid
    /// discretisation error remains.</para>
    /// </remarks>
    /// <param name="pwm">Position weight matrix with finite cells.</param>
    /// <param name="score">Score threshold α (finite).</param>
    /// <param name="background">Background (A, C, G, T), normalised to sum 1; uniform when null.</param>
    /// <returns>The p-value with its certified bounds.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="pwm"/> is null.</exception>
    /// <exception cref="ArgumentException">The matrix has a non-finite cell.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="score"/> is NaN or infinite.</exception>
    public static PwmPValueResult PwmScorePValue(
        PositionWeightMatrix pwm, double score, IReadOnlyList<double>? background = null)
    {
        if (!double.IsFinite(score))
            throw new ArgumentOutOfRangeException(nameof(score), score, "Score must be finite.");
        var engine = new PwmPValueEngine(pwm, ResolveBackground(background));
        return engine.ScoreToPValue(score);
    }

    /// <summary>
    /// Exact score threshold of a p-value (TFM-Pvalue <c>pv2sc</c>): the smallest score t attained by a word
    /// such that P(S ≥ t) ≤ <paramref name="pValue"/>, together with that exact p-value P(S ≥ t) — the
    /// largest achievable p-value not exceeding <paramref name="pValue"/>.
    /// </summary>
    /// <remarks>
    /// <para>Scores and the background model are as in <see cref="PwmScorePValue"/>. The returned
    /// <see cref="PwmPValueResult.Score"/> is the exact window score of a word (TFM-Pvalue returns the
    /// rounded value (α − offset)/g instead), so
    /// <see cref="PwmScorePValue"/>(pwm, t).PValue equals the returned p-value (same word set; only the order of the probability summation differs).</para>
    /// <para>Method: at scale g the integer score distribution (restricted, as in TFM-Pvalue's PvalueToScore loop,
    /// to the window certified by the previous scale) brackets t between the last bucket whose tail exceeds
    /// <paramref name="pValue"/> and the next bucket; the words of the bracket (± E) are enumerated with their
    /// exact scores as in <see cref="PwmScorePValue"/> and the tail is accumulated from the top. Refined until
    /// the bracket is resolved; otherwise the certified conservative threshold is returned with
    /// <see cref="PwmPValueResult.IsExact"/> = false.</para>
    /// <para>Edge cases: <paramref name="pValue"/> ≥ 1 gives the minimum word score with p-value 1; when even the
    /// best-scoring word(s) have P(S ≥ max) &gt; <paramref name="pValue"/> no word score qualifies and
    /// Score = +∞ with p-value 0.</para>
    /// </remarks>
    /// <param name="pwm">Position weight matrix with finite cells.</param>
    /// <param name="pValue">Target p-value in [0, 1].</param>
    /// <param name="background">Background (A, C, G, T), normalised to sum 1; uniform when null.</param>
    /// <returns>The threshold score and its p-value.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="pwm"/> is null.</exception>
    /// <exception cref="ArgumentException">The matrix has a non-finite cell.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="pValue"/> is not in [0, 1].</exception>
    public static PwmPValueResult PwmScoreThresholdForPValue(
        PositionWeightMatrix pwm, double pValue, IReadOnlyList<double>? background = null)
    {
        if (!(pValue >= 0 && pValue <= 1))
            throw new ArgumentOutOfRangeException(nameof(pValue), pValue, "P-value must be in [0, 1].");
        var engine = new PwmPValueEngine(pwm, ResolveBackground(background));
        return engine.PValueToScore(pValue);
    }

    /// <summary>Touzet–Varré integer-rounding engine shared by the two p-value entry points.</summary>
    private sealed class PwmPValueEngine
    {
        /// <summary>Rounding slack (integer units) absorbing floating-point error in g·W and in the window sums.</summary>
        private const long Pad = 2;

        /// <summary>Largest g·Σ_j max_b |W[b,j]| allowed (keeps every integer score and its rounding exact in a double).</summary>
        private const double MaxScaledRange = 1099511627776.0; // 2^40

        /// <summary>Maximum entries of one DP layer (integer or band enumeration) before the scale is refined or abandoned.</summary>
        private const int StateBudget = 1 << 21;
        private const int SuffixSetBudget = 1 << 20;

        private readonly double[,] _w;
        private readonly double[] _q;
        private readonly int _length;
        private readonly double _maxWord;   // exact (kernel-order) score of the best word
        private readonly double _minWord;   // exact (kernel-order) score of the worst word
        private readonly double _absRange;  // Σ_j max_b |W[b,j]|

        // Per-scale state.
        private long[,] _m = new long[4, 0];
        private long[] _sufMin = Array.Empty<long>();
        private long[] _sufMax = Array.Empty<long>();
        private double _errorMax;
        private long[]?[] _sufSet = Array.Empty<long[]?>();

        public PwmPValueEngine(PositionWeightMatrix pwm, double[] background)
        {
            ArgumentNullException.ThrowIfNull(pwm);
            foreach (double w in pwm.Matrix)
            {
                if (!double.IsFinite(w))
                    throw new ArgumentException(
                        "PWM p-values need a finite matrix (use a positive pseudocount).", nameof(pwm));
            }

            _w = pwm.Matrix;
            _q = background;
            _length = pwm.Length;

            double max = 0, min = 0, abs = 0;
            for (int j = 0; j < _length; j++)
            {
                double cMax = _w[0, j], cMin = _w[0, j], cAbs = Math.Abs(_w[0, j]);
                for (int b = 1; b < 4; b++)
                {
                    cMax = Math.Max(cMax, _w[b, j]);
                    cMin = Math.Min(cMin, _w[b, j]);
                    cAbs = Math.Max(cAbs, Math.Abs(_w[b, j]));
                }
                // Same left-to-right accumulation as ScorePwmWindow; IEEE addition is monotone, so these are the
                // largest / smallest window scores any word can produce.
                max += cMax;
                min += cMin;
                abs += cAbs;
            }
            _maxWord = max;
            _minWord = min;
            _absRange = abs;
        }

        public PwmPValueResult ScoreToPValue(double alpha)
        {
            if (alpha > _maxWord)
                return new PwmPValueResult(alpha, 0, 0, 0, true, 0);
            if (alpha <= _minWord)
                return new PwmPValueResult(alpha, 1, 1, 1, true, 0);

            double lower = 0, upper = 1, lastG = 0;
            foreach (double g in Scales())
            {
                SetScale(g);
                double ga = g * alpha;
                long lo = (long)Math.Floor(ga - _errorMax) - Pad;
                long inCut = (long)Math.Ceiling(ga) + Pad;
                if (!TryBandDistribution(lo, inCut, out double pIn, out Dictionary<long, (double P, double N)> band))
                    break;

                double bandMass = 0, bandWords = 0;
                foreach (var v in band.Values)
                {
                    bandMass += v.P;
                    bandWords += v.N;
                }
                lower = pIn;
                upper = Math.Min(1.0, pIn + bandMass);
                lastG = g;
                if (band.Count == 0)
                    return new PwmPValueResult(alpha, pIn, pIn, pIn, true, g);

                double hit = 0;
                bool done = EnumerateBand(lo, inCut, EnumerationBudget(bandWords), (s, p) => { if (s >= alpha) hit += p; });
                if (done)
                {
                    double pv = pIn + hit;
                    return new PwmPValueResult(alpha, pv, pv, pv, true, g);
                }
            }

            return new PwmPValueResult(alpha, upper, lower, upper, false, lastG);
        }

        public PwmPValueResult PValueToScore(double p)
        {
            if (p >= 1)
                return new PwmPValueResult(_minWord, 1, 1, 1, true, 0);

            PwmPValueResult top = ScoreToPValue(_maxWord);
            if (top.PValueLowerBound > p)
                return new PwmPValueResult(double.PositiveInfinity, 0, 0, 0, true, top.Granularity);

            // Certified bracket of the threshold t: lowScore < t ≤ highScore (refined at every scale, as the
            // [min, max] window of TFM-Pvalue's PvalueToScore loop).
            double lowScore = double.NegativeInfinity, highScore = double.PositiveInfinity;
            foreach (double g in Scales())
            {
                SetScale(g);
                long errCeil = (long)Math.Ceiling(_errorMax);
                long lo = double.IsNegativeInfinity(lowScore)
                    ? _sufMin[0]
                    : (long)Math.Floor(g * lowScore) - 2 * errCeil - 4 * Pad - 1;
                long hiCut = double.IsPositiveInfinity(highScore)
                    ? _sufMax[0] + 1
                    : (long)Math.Ceiling(g * highScore) + errCeil + 4 * Pad + 1;
                if (!TryBandDistribution(lo, hiCut, out double pAbove, out Dictionary<long, (double P, double N)> dist))
                    break;

                long[] keys = dist.Keys.ToArray();
                Array.Sort(keys);
                int n = keys.Length;
                var tail = new double[n];
                double acc = pAbove;
                for (int i = n - 1; i >= 0; i--)
                {
                    acc += dist[keys[i]].P;
                    tail[i] = acc;
                }

                int bi = -1;
                for (int i = n - 1; i >= 0; i--)
                {
                    if (tail[i] > p) { bi = i; break; }
                }
                if (bi < 0)
                    continue; // summation round-off: refine

                // b: last bucket whose tail exceeds p; a: the next one (tail ≤ p). Every word above the threshold
                // has M ≥ b − E − 2·pad; every word with M ≥ cIn scores at least (cIn − pad)/g.
                long bKey = keys[bi];
                bool hasA = bi + 1 < n;
                long cIn = hasA ? Math.Min(keys[bi + 1] + errCeil + 2 * Pad, hiCut) : hiCut;
                double pC = TailFrom(keys, tail, cIn, pAbove);
                double y = (cIn - Pad) / g;
                long bandLo = bKey - errCeil - 2 * Pad;

                lowScore = Math.Max(lowScore, (bKey - Pad) / g);
                if (hasA)
                    highScore = Math.Min(highScore, (keys[bi + 1] + _errorMax + Pad) / g);

                double bandWords = bandLo < lo ? double.PositiveInfinity : 0;
                foreach (long key in keys)
                {
                    if (key >= bandLo && key < cIn)
                        bandWords += dist[key].N;
                }

                var words = new List<(double Score, double P)>();
                if (!EnumerateBand(bandLo, cIn, EnumerationBudget(bandWords), (s, pr) => words.Add((s, pr))))
                    continue;

                words.Sort((u, v) => v.Score.CompareTo(u.Score));
                double cum = pC;
                double best = double.NaN, bestP = 0;
                bool infeasibleSeen = false;
                int k = 0;
                while (k < words.Count)
                {
                    double s = words[k].Score;
                    while (k < words.Count && words[k].Score == s)
                    {
                        cum += words[k].P;
                        k++;
                    }
                    if (s > y)
                        continue; // C ∩ U(s) unknown; P(U(s)) ≤ p anyway
                    if (cum <= p)
                    {
                        best = s;
                        bestP = cum;
                    }
                    else
                    {
                        infeasibleSeen = true;
                        break;
                    }
                }

                if (infeasibleSeen && !double.IsNaN(best))
                    return new PwmPValueResult(best, bestP, bestP, bestP, true, g);
            }

            double fallback = double.IsPositiveInfinity(highScore) ? _maxWord : Math.Min(highScore, _maxWord);
            PwmPValueResult conservative = ScoreToPValue(fallback);
            return conservative with { IsExact = false };
        }

        private static double TailFrom(long[] keys, double[] tail, long key, double above)
        {
            int idx = Array.BinarySearch(keys, key);
            if (idx < 0) idx = ~idx;
            return idx < keys.Length ? tail[idx] : above;
        }

        private IEnumerable<double> Scales()
        {
            for (double g = 10; g * Math.Max(_absRange, 1.0) <= MaxScaledRange; g *= 10)
                yield return g;
        }

        private void SetScale(double g)
        {
            _m = new long[4, _length];
            _errorMax = 0;
            for (int j = 0; j < _length; j++)
            {
                double colErr = 0;
                for (int b = 0; b < 4; b++)
                {
                    double x = _w[b, j] * g;
                    double f = Math.Floor(x);
                    _m[b, j] = (long)f;
                    colErr = Math.Max(colErr, x - f);
                }
                _errorMax += colErr;
            }

            _sufMin = new long[_length + 1];
            _sufMax = new long[_length + 1];
            for (int j = _length - 1; j >= 0; j--)
            {
                long cMin = _m[0, j], cMax = _m[0, j];
                for (int b = 1; b < 4; b++)
                {
                    cMin = Math.Min(cMin, _m[b, j]);
                    cMax = Math.Max(cMax, _m[b, j]);
                }
                _sufMin[j] = _sufMin[j + 1] + cMin;
                _sufMax[j] = _sufMax[j + 1] + cMax;
            }

            // Exact sets of reachable suffix scores for the last columns (meet-in-the-middle pruning of the band
            // enumeration); earlier layers fall back to the [sufMin, sufMax] bounds.
            _sufSet = new long[]?[_length + 1];
            _sufSet[_length] = new long[] { 0 };
            for (int j = _length - 1; j >= 0; j--)
            {
                long[] next = _sufSet[j + 1]!;
                if ((long)next.Length * 4 > SuffixSetBudget)
                    break;
                var all = new long[next.Length * 4];
                for (int b = 0; b < 4; b++)
                    for (int i = 0; i < next.Length; i++)
                        all[b * next.Length + i] = next[i] + _m[b, j];
                Array.Sort(all);
                int distinct = 0;
                for (int i = 0; i < all.Length; i++)
                {
                    if (distinct == 0 || all[i] != all[distinct - 1])
                        all[distinct++] = all[i];
                }
                Array.Resize(ref all, distinct);
                _sufSet[j] = all;
            }
        }

        /// <summary>
        /// State budget of a band enumeration: the full budget when the band holds at most that many words (every
        /// layer then fits, as each kept prefix ends in a distinct band word); otherwise a short attempt that still
        /// succeeds when ties merge the states (e.g. matrices with repeated values), before the scale is refined.
        /// </summary>
        private static int EnumerationBudget(double bandWords)
            => bandWords <= StateBudget ? StateBudget : StateBudget / 16;

        /// <summary>Whether a prefix with integer score <paramref name="key"/> (columns before j) can end in [lo, hiExclusive).</summary>
        private bool CanReach(int j, long key, long lo, long hiExclusive)
        {
            long[]? set = _sufSet[j];
            if (set is null)
                return key + _sufMax[j] >= lo && (hiExclusive == long.MaxValue || key + _sufMin[j] < hiExclusive);

            int idx = Array.BinarySearch(set, lo - key);
            if (idx < 0) idx = ~idx;
            return idx < set.Length && (hiExclusive == long.MaxValue || set[idx] < hiExclusive - key);
        }

        /// <summary>
        /// DP over integer partial scores keeping only prefixes that can still end in [lo, inCut); prefixes
        /// whose every completion reaches inCut are added to <paramref name="pIn"/>.
        /// </summary>
        private bool TryBandDistribution(long lo, long inCut, out double pIn, out Dictionary<long, (double P, double N)> band)
        {
            pIn = 0;
            var cur = new Dictionary<long, (double P, double N)> { [0] = (1.0, 1.0) };
            for (int j = 0; j < _length; j++)
            {
                var next = new Dictionary<long, (double P, double N)>();
                long sMin = _sufMin[j + 1], sMax = _sufMax[j + 1];
                foreach (var (key, v) in cur)
                {
                    for (int b = 0; b < 4; b++)
                    {
                        long nk = key + _m[b, j];
                        double np = v.P * _q[b];
                        if (nk + sMin >= inCut)
                            pIn += np;
                        else if (nk + sMax >= lo)
                        {
                            var old = next.GetValueOrDefault(nk);
                            next[nk] = (old.P + np, old.N + v.N);
                        }
                    }
                }
                if (next.Count > StateBudget)
                {
                    band = next;
                    return false;
                }
                cur = next;
            }

            band = cur;
            return true;
        }

        /// <summary>
        /// Exact enumeration of the words whose integer score lies in [lo, hiExclusive): a layered DP whose state is
        /// (integer partial score, exact double partial window score) — the window score is accumulated left to
        /// right exactly as <see cref="ScorePwmWindow"/>, and two prefixes with the same state have identical
        /// completions, so merging them is exact (branch-and-bound with memoisation). Reports each distinct final
        /// (score, probability) pair; false when the state budget is exceeded.
        /// </summary>
        private bool EnumerateBand(long lo, long hiExclusive, int budget, Action<double, double> visit)
        {
            var cur = new Dictionary<(long Key, double Score), double> { [(0L, 0.0)] = 1.0 };
            for (int j = 0; j < _length; j++)
            {
                var next = new Dictionary<(long Key, double Score), double>();
                foreach (var (state, pr) in cur)
                {
                    for (int b = 0; b < 4; b++)
                    {
                        long nk = state.Key + _m[b, j];
                        if (!CanReach(j + 1, nk, lo, hiExclusive))
                            continue;
                        var ns = (nk, state.Score + _w[b, j]);
                        next[ns] = next.GetValueOrDefault(ns) + pr * _q[b];
                    }
                }
                if (next.Count > budget)
                    return false;
                cur = next;
            }

            foreach (var (state, pr) in cur)
                visit(state.Score, pr);
            return true;
        }
    }
}

/// <summary>
/// Result of <see cref="MotifFinder.PwmScorePValue"/> / <see cref="MotifFinder.PwmScoreThresholdForPValue"/>.
/// </summary>
/// <param name="Score">Score threshold (the requested score, or the computed threshold).</param>
/// <param name="PValue">P(S ≥ <paramref name="Score"/>) under the background; the upper bound when not exact.</param>
/// <param name="PValueLowerBound">Certified lower bound of the p-value.</param>
/// <param name="PValueUpperBound">Certified upper bound of the p-value.</param>
/// <param name="IsExact">True when the p-value is exact (bounds coincide).</param>
/// <param name="Granularity">Scale g of the integer-rounded matrix that resolved the value (0 = decided without rounding).</param>
public readonly record struct PwmPValueResult(
    double Score,
    double PValue,
    double PValueLowerBound,
    double PValueUpperBound,
    bool IsExact,
    double Granularity);
