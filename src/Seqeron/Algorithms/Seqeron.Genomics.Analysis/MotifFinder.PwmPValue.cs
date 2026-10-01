namespace Seqeron.Genomics.Analysis;

/// <summary>
/// Exact PWM score p-values and p-value thresholds — the successive refinement of integer-rounded matrices of
/// Touzet &amp; Varré (2007), Algorithms Mol Biol 2:15 (TFM-Pvalue; reference C++ source: CRAN package TFMPvalue,
/// <c>src/Matrix.cpp</c> / <c>src/TFMpvalue.cpp</c>), under an i.i.d. background or an order-m Markov background
/// (state = integer partial score × last m letters, as the dinucleotide-background score distribution of MACRO-APE,
/// Vorontsov et al. 2013, <c>DiPWMScoresGenerator</c>; Markov sources: Boeva et al. 2007 AhoPro), for the DNA
/// <see cref="PositionWeightMatrix"/> and for <see cref="AlphabetPositionWeightMatrix"/> (K rows) through one engine.
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
        => PwmScorePValue(pwm, score, background, null);

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
        => PwmScoreThresholdForPValue(pwm, pValue, background, null);

    /// <summary>
    /// <see cref="PwmScorePValue(PositionWeightMatrix, double, IReadOnlyList{double}?)"/> with explicit search options
    /// (TFM-Pvalue initial granularity / maximal granularity / decrease factor, work budgets, exhaustive mode).
    /// </summary>
    /// <param name="pwm">Position weight matrix with finite cells.</param>
    /// <param name="score">Score threshold α (finite).</param>
    /// <param name="background">Background (A, C, G, T), normalised to sum 1; uniform when null.</param>
    /// <param name="options">Search options; null = <see cref="PwmPValueOptions.Default"/> (bit-identical to the 3-argument overload).</param>
    /// <exception cref="ArgumentException">Invalid options.</exception>
    public static PwmPValueResult PwmScorePValue(
        PositionWeightMatrix pwm, double score, IReadOnlyList<double>? background, PwmPValueOptions? options)
    {
        ValidatePValueScore(score);
        var engine = new PwmPValueEngine(pwm, ResolveBackground(background), options);
        return engine.ScoreToPValue(score);
    }

    /// <summary>
    /// <see cref="PwmScoreThresholdForPValue(PositionWeightMatrix, double, IReadOnlyList{double}?)"/> with explicit search options.
    /// </summary>
    /// <param name="pwm">Position weight matrix with finite cells.</param>
    /// <param name="pValue">Target p-value in [0, 1].</param>
    /// <param name="background">Background (A, C, G, T), normalised to sum 1; uniform when null.</param>
    /// <param name="options">Search options; null = <see cref="PwmPValueOptions.Default"/>.</param>
    public static PwmPValueResult PwmScoreThresholdForPValue(
        PositionWeightMatrix pwm, double pValue, IReadOnlyList<double>? background, PwmPValueOptions? options)
    {
        ValidatePValue(pValue);
        var engine = new PwmPValueEngine(pwm, ResolveBackground(background), options);
        return engine.PValueToScore(pValue);
    }

    /// <summary>
    /// Exact p-value of a PWM score under a Markov background: P(S ≥ <paramref name="score"/>) for a random word
    /// whose probability is the RSAT <c>segment_proba</c> of <paramref name="background"/>,
    /// P(w) = P(w[0..m−1]) · ∏_{c=m}^{L−1} P(w_c | w[c−m..c−1]) (for L &lt; m, the marginal of the prefix distribution).
    /// </summary>
    /// <remarks>
    /// <para>The Touzet–Varré refinement of <see cref="PwmScorePValue(PositionWeightMatrix, double, IReadOnlyList{double}?)"/>
    /// with the DP state extended by the last m letters (Markov context), exactly as MACRO-APE's dinucleotide-background
    /// score distribution keeps one score map per previous letter (<c>DiPWMScoresGenerator.recalc_score_hash</c>);
    /// the band enumeration keeps (integer partial score, context, exact partial window score). Exact as for the
    /// i.i.d. case (same <see cref="PwmPValueResult.IsExact"/> contract).</para>
    /// <para>Scores are the log-odds window scores of <see cref="CalculatePwmScores(string, PositionWeightMatrix)"/>
    /// (a PWM is scored against the background it was built with; only the null distribution of the words changes).</para>
    /// <para>Supported models: <see cref="OligoBackgroundModel.Equiprobable"/>, <see cref="OligoBackgroundModel.Bernoulli"/>
    /// (i.i.d.; identical to the <c>IReadOnlyList&lt;double&gt;</c> overload) and
    /// <see cref="OligoBackgroundModel.MarkovFromOligoFrequencies"/> (order m = table word length − 1 ≤ 10). The
    /// input-estimated models (<see cref="OligoBackgroundModel.BernoulliFromInput"/>, <see cref="OligoBackgroundModel.MarkovFromInput"/>)
    /// need sequences and RSAT's input Markov estimate is not a normalised word distribution (m-mer and (m+1)-mer
    /// frequencies come from different window sets), and <see cref="OligoBackgroundModel.Lexicon"/> is a maximal
    /// segmentation frequency, not a probability: these throw — build the chain with
    /// <see cref="OligoBackgroundModel.MarkovFromOligoFrequencies"/> from the (m+1)-mer counts instead.</para>
    /// </remarks>
    /// <param name="pwm">Position weight matrix with finite cells.</param>
    /// <param name="score">Score threshold α (finite).</param>
    /// <param name="background">The background model.</param>
    /// <param name="options">Search options; null = <see cref="PwmPValueOptions.Default"/>.</param>
    /// <exception cref="ArgumentNullException">An argument is null.</exception>
    /// <exception cref="ArgumentException">Non-finite matrix cell, unsupported model kind, or invalid options.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="score"/> is NaN or infinite.</exception>
    public static PwmPValueResult PwmMarkovScorePValue(
        PositionWeightMatrix pwm, double score, OligoBackgroundModel background, PwmPValueOptions? options = null)
    {
        ValidatePValueScore(score);
        return CreateMarkovEngine(pwm, background, options).ScoreToPValue(score);
    }

    /// <summary>
    /// Exact score threshold of a p-value under a Markov background (the <see cref="PwmScoreThresholdForPValue(PositionWeightMatrix, double, IReadOnlyList{double}?)"/>
    /// contract with the word distribution of <see cref="PwmMarkovScorePValue"/>).
    /// </summary>
    /// <param name="pwm">Position weight matrix with finite cells.</param>
    /// <param name="pValue">Target p-value in [0, 1].</param>
    /// <param name="background">The background model (see <see cref="PwmMarkovScorePValue"/>).</param>
    /// <param name="options">Search options; null = <see cref="PwmPValueOptions.Default"/>.</param>
    public static PwmPValueResult PwmMarkovScoreThresholdForPValue(
        PositionWeightMatrix pwm, double pValue, OligoBackgroundModel background, PwmPValueOptions? options = null)
    {
        ValidatePValue(pValue);
        return CreateMarkovEngine(pwm, background, options).PValueToScore(pValue);
    }

    /// <summary>
    /// Exact p-value of a score of a generic-alphabet PWM (protein, RNA, …): P(S ≥ <paramref name="score"/>) for a
    /// random word of length L drawn i.i.d. from <paramref name="background"/> over the K symbols, S being the window
    /// score of <see cref="CalculateAlphabetPwmScores"/> — the engine of
    /// <see cref="PwmScorePValue(PositionWeightMatrix, double, IReadOnlyList{double}?)"/> with K rows.
    /// </summary>
    /// <param name="pwm">Generic-alphabet PWM with finite cells (use a positive pseudocount).</param>
    /// <param name="score">Score threshold α (finite).</param>
    /// <param name="background">Background in alphabet order (finite, &gt; 0, normalised); uniform when null.</param>
    /// <param name="options">Search options; null = <see cref="PwmPValueOptions.Default"/>.</param>
    /// <exception cref="ArgumentNullException"><paramref name="pwm"/> is null.</exception>
    /// <exception cref="ArgumentException">The matrix has a non-finite cell, or invalid options.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="score"/> is NaN or infinite.</exception>
    public static PwmPValueResult AlphabetPwmScorePValue(
        AlphabetPositionWeightMatrix pwm, double score, IReadOnlyList<double>? background = null, PwmPValueOptions? options = null)
    {
        ArgumentNullException.ThrowIfNull(pwm);
        ValidatePValueScore(score);
        return CreateAlphabetEngine(pwm, background, options).ScoreToPValue(score);
    }

    /// <summary>
    /// Exact score threshold of a p-value for a generic-alphabet PWM (TFM-Pvalue <c>pv2sc</c>; contract of
    /// <see cref="PwmScoreThresholdForPValue(PositionWeightMatrix, double, IReadOnlyList{double}?)"/> with K rows).
    /// </summary>
    /// <param name="pwm">Generic-alphabet PWM with finite cells.</param>
    /// <param name="pValue">Target p-value in [0, 1].</param>
    /// <param name="background">Background in alphabet order; uniform when null.</param>
    /// <param name="options">Search options; null = <see cref="PwmPValueOptions.Default"/>.</param>
    public static PwmPValueResult AlphabetPwmScoreThresholdForPValue(
        AlphabetPositionWeightMatrix pwm, double pValue, IReadOnlyList<double>? background = null, PwmPValueOptions? options = null)
    {
        ArgumentNullException.ThrowIfNull(pwm);
        ValidatePValue(pValue);
        return CreateAlphabetEngine(pwm, background, options).PValueToScore(pValue);
    }

    private static void ValidatePValueScore(double score)
    {
        if (!double.IsFinite(score))
            throw new ArgumentOutOfRangeException(nameof(score), score, "Score must be finite.");
    }

    private static void ValidatePValue(double pValue)
    {
        if (!(pValue >= 0 && pValue <= 1))
            throw new ArgumentOutOfRangeException(nameof(pValue), pValue, "P-value must be in [0, 1].");
    }

    private static PwmPValueEngine CreateMarkovEngine(
        PositionWeightMatrix pwm, OligoBackgroundModel background, PwmPValueOptions? options)
    {
        ArgumentNullException.ThrowIfNull(pwm);
        ArgumentNullException.ThrowIfNull(background);
        PwmBackgroundChain chain = background.ToPValueChain(nameof(background));
        return chain.Order == 0
            ? new PwmPValueEngine(pwm, chain.Residues, options)
            : new PwmPValueEngine(pwm.Matrix, PwmAlphabetSize, pwm.Length, chain.Residues, chain, options, nameof(pwm));
    }

    private static PwmPValueEngine CreateAlphabetEngine(
        AlphabetPositionWeightMatrix pwm, IReadOnlyList<double>? background, PwmPValueOptions? options)
    {
        int k = pwm.Alphabet.Length;
        return new PwmPValueEngine(pwm.RawMatrix, k, pwm.Length,
            AlphabetPositionWeightMatrix.ResolveBackground(background, k), null, options, nameof(pwm));
    }

    /// <summary>
    /// Touzet–Varré integer-rounding engine shared by all p-value entry points: K matrix rows (4 for DNA, |alphabet| for
    /// <see cref="AlphabetPositionWeightMatrix"/>) and an i.i.d. or order-m Markov background. A DP state packs the integer
    /// partial score and the Markov context (last min(j, m) letters, base-K code) into one long:
    /// (key &lt;&lt; ctxBits) | ctx; with an i.i.d. background ctxBits = 0, so the state is the integer score itself and
    /// every operation is the one of the original 4-row i.i.d. engine.
    /// </summary>
    private sealed class PwmPValueEngine
    {
        /// <summary>Rounding slack (integer units) absorbing floating-point error in g·W and in the window sums.</summary>
        private const long Pad = 2;

        /// <summary>Largest g·Σ_j max_b |W[b,j]| allowed (keeps every integer score and its rounding exact in a double).</summary>
        private const double MaxScaledRange = 1099511627776.0; // 2^40

        private readonly double[,] _w;
        private readonly double[] _q;
        private readonly int _k;
        private readonly int _length;
        private readonly double _maxWord;   // exact (kernel-order) score of the best word
        private readonly double _minWord;   // exact (kernel-order) score of the worst word
        private readonly double _absRange;  // Σ_j max_b |W[b,j]|

        // Background chain (order 0 = i.i.d. _q).
        private readonly int _order;
        private readonly int _ctxCount;      // K^m
        private readonly int _ctxBits;
        private readonly long _ctxMask;
        private readonly double[][]? _prefixProb;  // [j] = P(first j letters), base-K code, j = 1 … m
        private readonly double[]? _transition;    // [ctx·K + b] = P(b | ctx)
        private readonly double[][]? _completion;  // [j][ctx] = mass of all completions of a column-j state (null: ≡ 1)
        private readonly double _totalMass = 1.0;  // Σ_w P(w) (< 1 only for a sub-stochastic Markov table)
        private readonly bool _sparse;             // some words have probability 0 (Markov zeros or a zero residue)

        // Options.
        private readonly int _stateBudget;
        private readonly int _suffixSetBudget;
        private readonly double _initialScale;
        private readonly double _scaleFactor;
        private readonly double _maxScale;
        private readonly bool _exhaustive;

        // Per-scale state.
        private double _scale;
        private long[,] _m = new long[4, 0];
        private long[] _sufMin = Array.Empty<long>();
        private long[] _sufMax = Array.Empty<long>();
        private double _errorMax;
        private long[]?[] _sufSet = Array.Empty<long[]?>();

        public PwmPValueEngine(PositionWeightMatrix pwm, double[] background, PwmPValueOptions? options = null)
            : this(RequirePwm(pwm).Matrix, PwmAlphabetSize, pwm.Length, background, null, options, nameof(pwm))
        {
        }

        private static PositionWeightMatrix RequirePwm(PositionWeightMatrix pwm)
        {
            ArgumentNullException.ThrowIfNull(pwm);
            return pwm;
        }

        public PwmPValueEngine(double[,] matrix, int k, int length, double[] background, PwmBackgroundChain? chain,
            PwmPValueOptions? options, string paramName)
        {
            foreach (double w in matrix)
            {
                if (!double.IsFinite(w))
                    throw new ArgumentException(
                        "PWM p-values need a finite matrix (use a positive pseudocount).", paramName);
            }

            options ??= PwmPValueOptions.Default;
            options.Validate(nameof(options));
            _stateBudget = options.MaxStates;
            _suffixSetBudget = options.MaxSuffixSet;
            _initialScale = 1.0 / options.InitialGranularity;
            _scaleFactor = options.DecreaseFactor;
            _maxScale = options.MaxGranularity is double gMax ? (1.0 / gMax) * (1 + 1e-9) : double.PositiveInfinity;
            _exhaustive = options.Exhaustive;

            _w = matrix;
            _q = background;
            _k = k;
            _length = length;
            _m = new long[k, 0];

            if (chain is { Order: > 0 })
            {
                _order = chain.Order;
                _ctxCount = chain.ContextCount;
                _prefixProb = chain.PrefixProbabilities;
                _transition = chain.Transitions;
                while ((1L << _ctxBits) < _ctxCount)
                    _ctxBits++;
                _ctxMask = (1L << _ctxBits) - 1;
                _completion = CompletionMasses(out _totalMass);
            }
            else
            {
                _ctxCount = 1;
            }

            double max = 0, min = 0, abs = 0;
            for (int j = 0; j < _length; j++)
            {
                double cMax = _w[0, j], cMin = _w[0, j], cAbs = Math.Abs(_w[0, j]);
                for (int b = 1; b < _k; b++)
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
            _sparse = _order > 0 || Array.IndexOf(_q, 0.0) >= 0;
            if (_sparse)
                (_maxWord, _minWord) = PossibleWordExtremes(paramName);
        }

        /// <summary>
        /// Largest / smallest window score over the words of positive probability under a Markov chain (a chain with
        /// zero transitions makes some words impossible): DP over (column, context) keeping the extreme left-to-right
        /// partial sums (IEEE addition is monotone, so the extreme prefix of a context gives the extreme completion).
        /// </summary>
        private (double Max, double Min) PossibleWordExtremes(string paramName)
        {
            var hi = new Dictionary<int, double> { [0] = 0.0 };
            var lo = new Dictionary<int, double> { [0] = 0.0 };
            for (int j = 0; j < _length; j++)
            {
                var nHi = new Dictionary<int, double>();
                var nLo = new Dictionary<int, double>();
                foreach (var (ctx, h) in hi)
                {
                    double l = lo[ctx];
                    for (int b = 0; b < _k; b++)
                    {
                        if (Step(1.0, j, ctx, b, out int nctx) == 0)
                            continue;
                        double sh = h + _w[b, j], sl = l + _w[b, j];
                        nHi[nctx] = nHi.TryGetValue(nctx, out double oh) ? Math.Max(oh, sh) : sh;
                        nLo[nctx] = nLo.TryGetValue(nctx, out double ol) ? Math.Min(ol, sl) : sl;
                    }
                }
                hi = nHi;
                lo = nLo;
            }

            if (hi.Count == 0)
                throw new ArgumentException(
                    $"The background gives every word of length {_length} probability 0.", paramName);
            return (hi.Values.Max(), lo.Values.Min());
        }

        /// <summary>
        /// Markov chains whose transition rows do not all sum to 1 (an RSAT table without pseudo-frequency where a
        /// context has no observed successor) give a word measure of total mass &lt; 1: the mass of all completions
        /// of a state, C[j][ctx], then multiplies the prefixes settled early by the look-ahead. Null when every
        /// C is 1 (to 1e-12), so a stochastic chain adds the prefix probability itself.
        /// </summary>
        private double[][]? CompletionMasses(out double totalMass)
        {
            var c = new double[_length + 1][];
            int Width(int j) => j < _order ? IntPow(_k, j) : _ctxCount;
            c[_length] = new double[Width(_length)];
            Array.Fill(c[_length], 1.0);
            bool allOne = true;
            for (int j = _length - 1; j >= 0; j--)
            {
                var col = new double[Width(j)];
                for (int ctx = 0; ctx < col.Length; ctx++)
                {
                    double sum = 0;
                    for (int b = 0; b < _k; b++)
                    {
                        int code = ctx * _k + b;
                        if (j < _order)
                            sum += _prefixProb![j + 1][code] * c[j + 1][code];
                        else
                            sum += _transition![code] * c[j + 1][code % _ctxCount];
                    }
                    if (j < _order)
                    {
                        double prefix = _prefixProb![j][ctx];
                        sum = prefix > 0 ? sum / prefix : 0;
                    }
                    col[ctx] = sum;
                    if (Math.Abs(sum - 1) > 1e-12 && (j >= _order || _prefixProb![j][ctx] > 0))
                        allOne = false;
                }
                c[j] = col;
            }

            totalMass = allOne ? 1.0 : c[0][0];
            return allOne ? null : c;
        }

        private static int IntPow(int k, int e)
        {
            int r = 1;
            for (int i = 0; i < e; i++) r *= k;
            return r;
        }

        public PwmPValueResult ScoreToPValue(double alpha)
        {
            if (alpha > _maxWord)
                return new PwmPValueResult(alpha, 0, 0, 0, true, 0);
            if (alpha <= _minWord)
                return new PwmPValueResult(alpha, _totalMass, _totalMass, _totalMass, true, 0);

            double lower = 0, upper = _totalMass, lastG = 0, firstG = 0;
            foreach (double g in Scales())
            {
                if (firstG == 0) firstG = g;
                if (TryScoreAtScale(alpha, g, _stateBudget, ref lower, ref upper, out bool overBudget, out PwmPValueResult result))
                    return result;
                if (overBudget)
                    break; // integer DP over budget at this scale (bounds of the previous scale kept)
                lastG = g;
            }

            if (_exhaustive)
            {
                // Exhaustive mode: resolve the band at the finest scale whose integer DP fitted, without a state limit.
                double g = firstG > 0 ? firstG : _initialScale;
                if (lastG > 0)
                    g = lastG;
                double lo2 = lower, up2 = upper;
                if (TryScoreAtScale(alpha, g, int.MaxValue, ref lo2, ref up2, out _, out PwmPValueResult exact))
                    return exact;
            }

            return new PwmPValueResult(alpha, upper, lower, upper, false, lastG);
        }

        /// <summary>
        /// One refinement step at scale g. True with the exact result when resolved; otherwise updates the certified
        /// bounds (or, when the integer DP exceeds the budget, leaves them and sets <paramref name="overBudget"/>).
        /// </summary>
        private bool TryScoreAtScale(double alpha, double g, int budget, ref double lower, ref double upper,
            out bool overBudget, out PwmPValueResult result)
        {
            result = default;
            overBudget = false;
            SetScale(g);
            double ga = g * alpha;
            long lo = (long)Math.Floor(ga - _errorMax) - Pad;
            long inCut = (long)Math.Ceiling(ga) + Pad;
            if (!TryBandDistribution(lo, inCut, budget, out double pIn, out Dictionary<long, (double P, double N)> band))
            {
                overBudget = true;
                return false;
            }

            double bandMass = 0, bandWords = 0;
            foreach (var v in band.Values)
            {
                bandMass += v.P;
                bandWords += v.N;
            }
            lower = pIn;
            upper = Math.Min(_totalMass, pIn + bandMass);
            if (band.Count == 0)
            {
                result = new PwmPValueResult(alpha, pIn, pIn, pIn, true, g);
                return true;
            }

            double hit = 0;
            int enumerationBudget = budget == int.MaxValue ? int.MaxValue : EnumerationBudget(bandWords);
            bool done = EnumerateBand(lo, inCut, enumerationBudget, (s, p) => { if (s >= alpha) hit += p; });
            if (done)
            {
                double pv = pIn + hit;
                result = new PwmPValueResult(alpha, pv, pv, pv, true, g);
                return true;
            }

            return false;
        }

        public PwmPValueResult PValueToScore(double p)
        {
            if (p >= _totalMass)
                return new PwmPValueResult(_minWord, _totalMass, _totalMass, _totalMass, true, 0);

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
                if (!TryBandDistribution(lo, hiCut, _stateBudget, out double pAbove, out Dictionary<long, (double P, double N)> packed))
                    break;
                Dictionary<long, (double P, double N)> dist = ByIntegerScore(packed);

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
                // Words with M ≥ cIn score ≥ y; when no word reaches cIn every band score is decided.
                double y = cIn > _sufMax[0] ? double.PositiveInfinity : (cIn - Pad) / g;
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
                double best = double.NaN, bestP = 0, infeasibleScore = double.NaN;
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
                        infeasibleScore = s;
                        break;
                    }
                }

                if (infeasibleSeen && !double.IsNaN(best))
                    return new PwmPValueResult(best, bestP, bestP, bestP, true, g);
                if (infeasibleSeen)
                {
                    // Every attained score ≤ y has an exact tail > p, so t > infeasibleScore and t is an attained score
                    // above y: P(S ≥ highScore) ≤ p does not bound t from above when no word scores in [highScore, t)
                    // (sparse score sets, e.g. Markov backgrounds with impossible words). Widen the window upwards.
                    lowScore = Math.Max(lowScore, infeasibleScore);
                    highScore = double.PositiveInfinity;
                }
            }

            if (_exhaustive)
                return ExhaustiveThreshold(p, lowScore, highScore);

            double fallback = double.IsPositiveInfinity(highScore) ? _maxWord : Math.Min(highScore, _maxWord);
            PwmPValueResult conservative = ScoreToPValue(fallback);
            return conservative with { IsExact = false };
        }

        /// <summary>
        /// Exhaustive inverse: at the last scale, the words of the certified bracket lowScore &lt; t ≤ highScore are
        /// enumerated without a state limit. With cIn = ⌈g·highScore⌉ + pad + 1 every word with integer score ≥ cIn
        /// scores above highScore (≥ t) and is counted in P(M ≥ cIn); every word scoring ≥ lowScore has
        /// M ≥ ⌊g·lowScore − E⌋ − pad. So for every enumerated score s ≤ y = (cIn − pad)/g the accumulated tail is
        /// the exact P(S ≥ s), t &lt; y, and P(S ≥ lowScore) &gt; p: the first score whose tail exceeds p ends the scan.
        /// </summary>
        private PwmPValueResult ExhaustiveThreshold(double p, double lowScore, double highScore)
        {
            if (_scale == 0)
                SetScale(_initialScale);
            double g = _scale;
            long lo = double.IsNegativeInfinity(lowScore)
                ? _sufMin[0]
                : (long)Math.Floor(g * lowScore - _errorMax) - Pad - 1;
            long cIn = _sufMax[0] + 1;
            if (!double.IsPositiveInfinity(highScore))
                cIn = Math.Min(cIn, (long)Math.Ceiling(g * highScore) + Pad + 1);
            double y = cIn > _sufMax[0] ? double.PositiveInfinity : (cIn - Pad) / g;

            TryBandDistribution(lo, cIn, int.MaxValue, out double pC, out _);
            var words = new List<(double Score, double P)>();
            EnumerateBand(lo, cIn, int.MaxValue, (s, pr) => words.Add((s, pr)));
            words.Sort((u, v) => v.Score.CompareTo(u.Score));

            double cum = pC, best = double.NaN, bestP = 0;
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
                    continue;
                if (cum > p)
                {
                    if (double.IsNaN(best) && !double.IsPositiveInfinity(highScore))
                        return ExhaustiveThreshold(p, s, double.PositiveInfinity); // t lies above y (see PValueToScore)
                    break;
                }
                best = s;
                bestP = cum;
            }

            return double.IsNaN(best)
                ? new PwmPValueResult(double.PositiveInfinity, 0, 0, 0, true, g)
                : new PwmPValueResult(best, bestP, bestP, bestP, true, g);
        }

        private static double TailFrom(long[] keys, double[] tail, long key, double above)
        {
            int idx = Array.BinarySearch(keys, key);
            if (idx < 0) idx = ~idx;
            return idx < keys.Length ? tail[idx] : above;
        }

        private IEnumerable<double> Scales()
        {
            for (double g = _initialScale; g * Math.Max(_absRange, 1.0) <= MaxScaledRange && g <= _maxScale; g *= _scaleFactor)
                yield return g;
        }

        private void SetScale(double g)
        {
            _scale = g;
            _m = new long[_k, _length];
            _errorMax = 0;
            for (int j = 0; j < _length; j++)
            {
                double colErr = 0;
                for (int b = 0; b < _k; b++)
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
                for (int b = 1; b < _k; b++)
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
                if ((long)next.Length * _k > _suffixSetBudget)
                    break;
                var all = new long[next.Length * _k];
                for (int b = 0; b < _k; b++)
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
        private int EnumerationBudget(double bandWords)
            => bandWords <= _stateBudget ? _stateBudget : _stateBudget / 16;

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
        /// Probability of a prefix in context <paramref name="ctx"/> (probability <paramref name="p"/>) extended by
        /// symbol <paramref name="b"/> at column <paramref name="j"/>, and the new context. i.i.d.: p·q[b]. Markov of
        /// order m: for j &lt; m the context is the whole prefix, so the probability is the prefix marginal
        /// P(w[0..j]); afterwards p·P(b | last m letters).
        /// </summary>
        private double Step(double p, int j, int ctx, int b, out int nextCtx)
        {
            if (_order == 0)
            {
                nextCtx = 0;
                return p * _q[b];
            }

            int code = ctx * _k + b;
            if (j < _order)
            {
                nextCtx = code;
                return _prefixProb![j + 1][code];
            }

            nextCtx = code % _ctxCount;
            return p * _transition![code];
        }

        private long Pack(long key, int ctx) => (key << _ctxBits) | (long)ctx;

        /// <summary>Sums packed (integer score, context) states by integer score (identity for an i.i.d. background).</summary>
        private Dictionary<long, (double P, double N)> ByIntegerScore(Dictionary<long, (double P, double N)> packed)
        {
            if (_ctxBits == 0)
                return packed;
            var byKey = new Dictionary<long, (double P, double N)>();
            foreach (var (state, v) in packed)
            {
                long key = state >> _ctxBits;
                var old = byKey.GetValueOrDefault(key);
                byKey[key] = (old.P + v.P, old.N + v.N);
            }
            return byKey;
        }

        /// <summary>
        /// DP over integer partial scores (× Markov context) keeping only prefixes that can still end in [lo, inCut);
        /// prefixes whose every completion reaches inCut are added to <paramref name="pIn"/>.
        /// </summary>
        private bool TryBandDistribution(long lo, long inCut, int budget, out double pIn, out Dictionary<long, (double P, double N)> band)
        {
            pIn = 0;
            var cur = new Dictionary<long, (double P, double N)> { [0] = (1.0, 1.0) };
            for (int j = 0; j < _length; j++)
            {
                var next = new Dictionary<long, (double P, double N)>();
                long sMin = _sufMin[j + 1], sMax = _sufMax[j + 1];
                foreach (var (state, v) in cur)
                {
                    long key = state >> _ctxBits;
                    int ctx = (int)(state & _ctxMask);
                    for (int b = 0; b < _k; b++)
                    {
                        long nk = key + _m[b, j];
                        double np = Step(v.P, j, ctx, b, out int nctx);
                        if (np == 0 && _sparse)
                            continue; // impossible under the Markov chain
                        if (nk + sMin >= inCut)
                            pIn += _completion is null ? np : np * _completion[j + 1][nctx];
                        else if (nk + sMax >= lo)
                        {
                            long ns = Pack(nk, nctx);
                            var old = next.GetValueOrDefault(ns);
                            next[ns] = (old.P + np, old.N + v.N);
                        }
                    }
                }
                if (next.Count > budget)
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
        /// (integer partial score × Markov context, exact double partial window score) — the window score is accumulated
        /// left to right exactly as <see cref="ScorePwmWindow"/>, and two prefixes with the same state have identical
        /// completions, so merging them is exact (branch-and-bound with memoisation). Reports each distinct final
        /// (score, probability) state; false when the state budget is exceeded.
        /// </summary>
        private bool EnumerateBand(long lo, long hiExclusive, int budget, Action<double, double> visit)
        {
            var cur = new Dictionary<(long Key, double Score), double> { [(0L, 0.0)] = 1.0 };
            for (int j = 0; j < _length; j++)
            {
                var next = new Dictionary<(long Key, double Score), double>();
                foreach (var (state, pr) in cur)
                {
                    long key = state.Key >> _ctxBits;
                    int ctx = (int)(state.Key & _ctxMask);
                    for (int b = 0; b < _k; b++)
                    {
                        long nk = key + _m[b, j];
                        if (!CanReach(j + 1, nk, lo, hiExclusive))
                            continue;
                        double np = Step(pr, j, ctx, b, out int nctx);
                        if (np == 0 && _sparse)
                            continue;
                        var ns = (Pack(nk, nctx), state.Score + _w[b, j]);
                        next[ns] = next.GetValueOrDefault(ns) + np;
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

/// <summary>
/// Search options of the exact PWM p-value engine (<see cref="MotifFinder.PwmScorePValue(PositionWeightMatrix, double, IReadOnlyList{double}?, PwmPValueOptions?)"/>
/// and the other p-value / threshold entry points) — the TFM-Pvalue driver parameters (Touzet &amp; Varré 2007; CRAN
/// TFMPvalue <c>src/TFMpvalue.cpp</c> <c>testScoreToPvalue</c> / <c>testPvalueToScore</c>: <c>initialGranularity</c>,
/// <c>maxGranularity</c>, <c>decrgr</c>) plus the work budgets of this implementation and an exhaustive mode.
/// </summary>
/// <remarks>
/// A granularity is the rounding step 1/g of the integer matrix ⌊g·W⌋ (TFM-Pvalue's convention);
/// <see cref="PwmPValueResult.Granularity"/> reports the scale g itself. Scales tried: g = 1/<see cref="InitialGranularity"/>,
/// then ×<see cref="DecreaseFactor"/>, while g ≤ 1/<see cref="MaxGranularity"/> and g·Σ_j max_b |W[b,j]| ≤ 2⁴⁰ (the
/// latter always applies: it keeps every integer score exact in a double). The defaults reproduce the 3-argument
/// overloads bit for bit.
/// </remarks>
public sealed record PwmPValueOptions
{
    /// <summary>Default options (initial granularity 0.1, decrease factor 10, budgets 2²¹ / 2²⁰, no exhaustive fallback).</summary>
    public static PwmPValueOptions Default { get; } = new();

    /// <summary>
    /// Exhaustive mode: the default refinement, then — when it does not resolve the value — one band enumeration
    /// without a state limit, so the result is always exact (<see cref="PwmPValueResult.IsExact"/> = true); time and
    /// memory are then bounded only by the number of distinct (partial score, context) states (≤ K^L words).
    /// </summary>
    public static PwmPValueOptions Exact { get; } = new() { Exhaustive = true };

    /// <summary>Rounding step of the first integer matrix (TFM-Pvalue <c>initialGranularity</c>, default 0.1 → g = 10).</summary>
    public double InitialGranularity { get; init; } = 0.1;

    /// <summary>
    /// Finest rounding step tried (TFM-Pvalue <c>maxGranularity</c>, 1e-9 / 1e-10 there); null (default) = limited only
    /// by g·Σ_j max_b |W[b,j]| ≤ 2⁴⁰.
    /// </summary>
    public double? MaxGranularity { get; init; }

    /// <summary>Factor by which the granularity decreases (the scale grows) between refinements (TFM-Pvalue <c>decrgr</c>, default 10).</summary>
    public double DecreaseFactor { get; init; } = 10;

    /// <summary>Maximum DP states per column (integer distribution and band enumeration) before a scale is abandoned (default 2²¹).</summary>
    public int MaxStates { get; init; } = 1 << 21;

    /// <summary>Maximum size of an exact reachable-suffix-score set used for pruning (default 2²⁰).</summary>
    public int MaxSuffixSet { get; init; } = 1 << 20;

    /// <summary>Resolve an unresolved band by unbounded enumeration (see <see cref="Exact"/>; default false).</summary>
    public bool Exhaustive { get; init; }

    internal void Validate(string paramName)
    {
        if (!(double.IsFinite(InitialGranularity) && InitialGranularity > 0))
            throw new ArgumentException("InitialGranularity must be finite and positive.", paramName);
        if (MaxGranularity is double m && !(double.IsFinite(m) && m > 0 && m <= InitialGranularity))
            throw new ArgumentException("MaxGranularity must be finite, positive and at most InitialGranularity.", paramName);
        if (!(double.IsFinite(DecreaseFactor) && DecreaseFactor > 1))
            throw new ArgumentException("DecreaseFactor must be finite and greater than 1.", paramName);
        if (MaxStates < 1)
            throw new ArgumentException("MaxStates must be at least 1.", paramName);
        if (MaxSuffixSet < 1)
            throw new ArgumentException("MaxSuffixSet must be at least 1.", paramName);
    }
}

/// <summary>
/// An order-m background chain over K symbols for the p-value engine: prefix marginals P(w[0..j)) (j = 1 … m) and
/// transitions P(b | last m letters); contexts are base-K codes, first letter most significant.
/// </summary>
internal sealed class PwmBackgroundChain
{
    public PwmBackgroundChain(int k, double[] residues)
    {
        K = k;
        Order = 0;
        Residues = residues;
        ContextCount = 1;
        PrefixProbabilities = Array.Empty<double[]>();
        Transitions = Array.Empty<double>();
    }

    /// <param name="k">Alphabet size.</param>
    /// <param name="order">m ≥ 1.</param>
    /// <param name="initial">P(x) for the K^m m-letter prefixes.</param>
    /// <param name="transitions">P(b | x) at [x·K + b].</param>
    public PwmBackgroundChain(int k, int order, double[] initial, double[] transitions)
    {
        K = k;
        Order = order;
        ContextCount = initial.Length;
        Transitions = transitions;
        Residues = Array.Empty<double>();
        PrefixProbabilities = new double[order + 1][];
        PrefixProbabilities[order] = initial;
        for (int j = order - 1; j >= 1; j--)
        {
            double[] longer = PrefixProbabilities[j + 1];
            var shorter = new double[longer.Length / k];
            for (int x = 0; x < shorter.Length; x++)
            {
                double sum = 0;
                for (int b = 0; b < k; b++)
                    sum += longer[x * k + b];
                shorter[x] = sum;
            }
            PrefixProbabilities[j] = shorter;
        }
        PrefixProbabilities[0] = new[] { 1.0 };
    }

    public int K { get; }
    public int Order { get; }
    public int ContextCount { get; }
    public double[] Residues { get; }
    public double[][] PrefixProbabilities { get; }
    public double[] Transitions { get; }
}
