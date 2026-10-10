namespace Seqeron.Genomics.Infrastructure
{
    /// <summary>
    /// Common statistical functions used across genomics analyzers.
    /// </summary>
    public static class StatisticsHelper
    {
        /// <summary>
        /// Calculates the cumulative distribution function of the standard normal distribution.
        /// </summary>
        public static double NormalCDF(double x)
        {
            return 0.5 * (1 + Erf(x / Math.Sqrt(2)));
        }

        /// <summary>
        /// Calculates the error function using Horner's method.
        /// </summary>
        public static double Erf(double x)
        {
            double a1 = 0.254829592, a2 = -0.284496736, a3 = 1.421413741;
            double a4 = -1.453152027, a5 = 1.061405429, p = 0.3275911;

            int sign = x < 0 ? -1 : 1;
            x = Math.Abs(x);

            double t = 1.0 / (1.0 + p * x);
            double y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * Math.Exp(-x * x);

            return sign * y;
        }

        /// <summary>
        /// Population variance σ² = Σ(xᵢ − μ)² / N (division by N, not the Bessel-corrected N − 1),
        /// for data that form the complete population. Returns 0 for an empty list.
        /// Matches <c>numpy.var(x)</c> (default <c>ddof=0</c>); e.g. {12, 13, 12, 14, 19} → 6.8.
        /// </summary>
        /// <param name="values">The complete population of values.</param>
        /// <exception cref="ArgumentNullException"><paramref name="values"/> is null.</exception>
        public static double PopulationVariance(IReadOnlyList<double> values)
        {
            ArgumentNullException.ThrowIfNull(values);
            if (values.Count == 0) return 0;

            double sum = 0;
            for (int i = 0; i < values.Count; i++) sum += values[i];
            double mean = sum / values.Count;

            double sumSq = 0;
            for (int i = 0; i < values.Count; i++)
            {
                double d = values[i] - mean;
                sumSq += d * d;
            }

            return sumSq / values.Count;
        }

        /// <summary>
        /// Sample median with R's <c>median.default</c> convention: the central order statistic for an odd count,
        /// the mean of the two central order statistics for an even count. Returns <see cref="double.NaN"/> when any
        /// value is NaN (R: <c>median(c(1, NaN))</c> is <c>NA</c>). The input is not mutated.
        /// Matches R <c>median</c> and <c>numpy.median</c>; e.g. {0.2, 0.4, 0.6, 0.8} → 0.5, {3, 1, 2} → 2.
        /// </summary>
        /// <param name="values">The values (at least one).</param>
        /// <exception cref="ArgumentNullException"><paramref name="values"/> is null.</exception>
        /// <exception cref="ArgumentException"><paramref name="values"/> is empty.</exception>
        public static double Median(IReadOnlyList<double> values)
        {
            ArgumentNullException.ThrowIfNull(values);
            int n = values.Count;
            if (n == 0)
            {
                throw new ArgumentException("The median of an empty set is undefined.", nameof(values));
            }

            var sorted = new double[n];
            for (int i = 0; i < n; i++)
            {
                double v = values[i];
                if (double.IsNaN(v)) return double.NaN;
                sorted[i] = v;
            }

            Array.Sort(sorted);
            int mid = n / 2;
            if (n % 2 == 1) return sorted[mid];

            double lo = sorted[mid - 1], hi = sorted[mid];
            double m = (lo + hi) / 2.0;
            // Overflow guard for two huge finite central values (R's long-double mean does not overflow).
            return double.IsInfinity(m) && double.IsFinite(lo) && double.IsFinite(hi) ? lo / 2.0 + hi / 2.0 : m;
        }

        /// <summary>
        /// Shannon diversity index H′ = −Σ pᵢ·ln pᵢ (natural logarithm) over class counts, pᵢ = cᵢ / Σc; zero counts
        /// contribute nothing (0·ln 0 = 0). Shannon (1948); identical to scikit-bio
        /// <c>skbio.diversity.alpha.shannon(counts, base=math.e)</c> and <c>scipy.stats.entropy(counts)</c>;
        /// e.g. {2, 2} → ln 2, {1, 1, 1, 1} → ln 4, {5} → 0.
        /// </summary>
        /// <param name="counts">Non-negative class counts with a positive total.</param>
        /// <exception cref="ArgumentNullException"><paramref name="counts"/> is null.</exception>
        /// <exception cref="ArgumentException">a count is negative, or the total is 0.</exception>
        public static double ShannonIndex(IReadOnlyList<int> counts)
        {
            ArgumentNullException.ThrowIfNull(counts);
            long total = 0;
            for (int i = 0; i < counts.Count; i++)
            {
                if (counts[i] < 0)
                {
                    throw new ArgumentException($"Counts must be non-negative; got {counts[i]} at index {i}.", nameof(counts));
                }

                total += counts[i];
            }

            if (total == 0)
            {
                throw new ArgumentException("The total count must be positive.", nameof(counts));
            }

            double h = 0.0;
            for (int i = 0; i < counts.Count; i++)
            {
                if (counts[i] == 0) continue;
                double p = (double)counts[i] / total;
                h -= p * Math.Log(p);
            }

            return h;
        }

        /// <summary>
        /// Shannon diversity index H′ = −Σ pᵢ·ln pᵢ (natural logarithm) over non-negative real weights, pᵢ = wᵢ / Σw;
        /// zero weights contribute nothing (0·ln 0 = 0). The real-valued counterpart of
        /// <see cref="ShannonIndex(IReadOnlyList{int})"/> (same accumulation), equal to <c>scipy.stats.entropy(weights)</c>
        /// (which normalises its input the same way); e.g. {0.6, 0.3, 0.1} → 0.8979457248567798, {0.5, 0.5} → ln 2.
        /// </summary>
        /// <param name="weights">Non-negative finite weights (proportions, fractions or abundances) with a positive total.</param>
        /// <exception cref="ArgumentNullException"><paramref name="weights"/> is null.</exception>
        /// <exception cref="ArgumentException">a weight is negative, NaN or infinite, or the total is 0.</exception>
        public static double ShannonIndexOfWeights(IReadOnlyList<double> weights)
        {
            ArgumentNullException.ThrowIfNull(weights);
            double total = 0.0;
            for (int i = 0; i < weights.Count; i++)
            {
                double w = weights[i];
                if (!double.IsFinite(w) || w < 0.0)
                {
                    throw new ArgumentException($"Weights must be finite and non-negative; got {w} at index {i}.", nameof(weights));
                }

                total += w;
            }

            if (total <= 0.0 || !double.IsFinite(total))
            {
                throw new ArgumentException("The total weight must be positive and finite.", nameof(weights));
            }

            double h = 0.0;
            for (int i = 0; i < weights.Count; i++)
            {
                if (weights[i] == 0.0) continue;
                double p = weights[i] / total;
                h -= p * Math.Log(p);
            }

            return h;
        }

        /// <summary>
        /// Binomial upper tail P(X ≥ <paramref name="successes"/>) for X ~ Binomial(<paramref name="trials"/>, <paramref name="p"/>)
        /// — the "sum of binomials" of RSAT <c>RSAT::stats::binomial_boe</c> / <c>sum_of_binomials</c> (van Helden et al. 1998)
        /// and <c>scipy.stats.binom.sf(successes − 1, trials, p)</c>; e.g. (6, 60, 1/256) → 1.4844942103072820e-07.
        /// Evaluated by <see cref="LogBinomialUpperTail"/> (log space), so values below the double range underflow to 0
        /// only at the very end; use <see cref="LogBinomialUpperTail"/> to keep them.
        /// </summary>
        /// <param name="successes">Observed number of successes (any integer; ≤ 0 gives 1, &gt; trials gives 0).</param>
        /// <param name="trials">Number of Bernoulli trials (≥ 0).</param>
        /// <param name="p">Success probability in [0, 1].</param>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="trials"/> &lt; 0 or <paramref name="p"/> outside [0, 1] / NaN.</exception>
        public static double BinomialUpperTail(long successes, long trials, double p)
        {
            if (!(p >= 0.0 && p <= 1.0))
                throw new ArgumentOutOfRangeException(nameof(p), p, "Probability must be in [0, 1].");
            return Math.Exp(LogBinomialUpperTail(successes, trials, p == 0.0 ? double.NegativeInfinity : Math.Log(p)));
        }

        /// <summary>
        /// Natural logarithm of the binomial upper tail ln P(X ≥ <paramref name="successes"/>), X ~ Binomial(n, p), with the
        /// success probability given as <paramref name="logP"/> = ln p (so p may be far below the double range, e.g. the
        /// expected frequency 4^−k of a long oligonucleotide). Equals <c>scipy.stats.binom.logsf(successes − 1, n, p)</c>.
        /// </summary>
        /// <remarks>
        /// Each probability mass term is computed with Loader's saddle-point expansion (C. Loader 2000, "Fast and accurate
        /// computation of binomial probabilities"; the <c>dbinom_raw</c>/<c>stirlerr</c>/<c>bd0</c> algorithm of R nmath),
        /// which is accurate to a few ulps for any n. The tail is summed from the first term outward with the exact
        /// ratio recurrence P(x+1)/P(x) = (n − x)/(x + 1) · p/q until the terms no longer change the sum (the same stopping
        /// rule as RSAT <c>sum_of_binomials</c>). When the requested tail contains the mean (successes ≤ n·p) the
        /// complementary lower tail is summed instead, P(X ≥ s) = 1 − P(X ≤ s − 1), which is then ≥ ~½ and has no
        /// cancellation problem.
        /// </remarks>
        /// <param name="successes">Observed number of successes.</param>
        /// <param name="trials">Number of trials (≥ 0).</param>
        /// <param name="logP">ln p, in [−∞, 0]; −∞ means p = 0.</param>
        /// <returns>ln P(X ≥ successes), in [−∞, 0].</returns>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="trials"/> &lt; 0 or <paramref name="logP"/> &gt; 0 / NaN.</exception>
        public static double LogBinomialUpperTail(long successes, long trials, double logP)
        {
            ArgumentOutOfRangeException.ThrowIfNegative(trials);
            if (double.IsNaN(logP) || logP > 0.0)
                throw new ArgumentOutOfRangeException(nameof(logP), logP, "ln p must be in [-inf, 0].");

            if (successes <= 0) return 0.0;
            if (successes > trials) return double.NegativeInfinity;
            if (double.IsNegativeInfinity(logP)) return double.NegativeInfinity; // p = 0: X ≡ 0
            if (logP == 0.0) return 0.0;                                         // p = 1: X ≡ n ≥ successes

            double p = Math.Exp(logP);
            double q = -ExpM1(logP);           // 1 − p without cancellation
            double logQ = Math.Log(q);
            double mean = trials * p;

            if (successes > mean)
            {
                // Right tail: terms decrease monotonically from x = successes.
                double logFirst = LogBinomialPmf(successes, trials, logP, p, logQ);
                double pOverQ = Math.Exp(logP - logQ);
                double sum = 1.0, term = 1.0;
                for (long x = successes; x < trials; x++)
                {
                    term *= (trials - x) / (x + 1.0) * pOverQ;
                    double next = sum + term;
                    if (next == sum) break;
                    sum = next;
                }
                return logFirst + Math.Log(sum);
            }
            else
            {
                // The tail contains the mean: 1 − P(X ≤ successes − 1), lower terms decrease from x = successes − 1.
                long start = successes - 1;
                double logFirst = LogBinomialPmf(start, trials, logP, p, logQ);
                double qOverP = Math.Exp(logQ - logP);   // ≤ n here, since n·p ≥ successes ≥ 1
                double sum = 1.0, term = 1.0;
                for (long x = start; x > 0; x--)
                {
                    term *= x / (trials - x + 1.0) * qOverP;
                    double next = sum + term;
                    if (next == sum) break;
                    sum = next;
                }
                double lower = Math.Exp(logFirst + Math.Log(sum));
                return lower >= 1.0 ? double.NegativeInfinity : Log1P(-lower);
            }
        }

        /// <summary>
        /// ln P(<paramref name="from"/> ≤ X ≤ <paramref name="to"/>) for X ~ Poisson(<paramref name="lambda"/>) — RSAT
        /// <c>RSAT::stats::sum_of_poisson($lambda, $from, $to)</c> (used by <c>oligo-analysis -calibN</c> with the right tail
        /// from occ to the number of positions); equals <c>log(scipy.stats.poisson.cdf(to, λ) − poisson.cdf(from − 1, λ))</c>
        /// without cancellation. Terms ln P(x) = −λ + x·ln λ − ln x! (exact log-factorial), summed outward from the end of the
        /// range nearest the mode with the ratio recurrence P(x + 1) = P(x)·λ/(x + 1) until the sum stops changing; a range
        /// that contains the mode is 1 − P(X &lt; from) − P(X &gt; to).
        /// </summary>
        /// <param name="from">Lower bound (values &lt; 0 are treated as 0).</param>
        /// <param name="to">Upper bound (inclusive).</param>
        /// <param name="lambda">Mean, finite and &gt; 0.</param>
        /// <returns>ln probability in [−∞, 0] (−∞ for an empty range).</returns>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="lambda"/> is not finite and positive.</exception>
        public static double LogPoissonRangeProbability(long from, long to, double lambda)
        {
            if (!(lambda > 0) || double.IsPositiveInfinity(lambda))
                throw new ArgumentOutOfRangeException(nameof(lambda), lambda, "The Poisson mean must be finite and positive.");
            double logLambda = Math.Log(lambda);
            long mode = Math.Max(0L, (long)Math.Ceiling(lambda - 1));
            return LogDiscreteRangeSum(Math.Max(0L, from), to, mode,
                x => -lambda + x * logLambda - LogFactorial(x),
                x => lambda / (x + 1.0));
        }

        /// <summary>
        /// ln P(<paramref name="from"/> ≤ X ≤ <paramref name="to"/>) for the negative binomial with the given mean m and
        /// variance v &gt; m, parameterised as RSAT <c>RSAT::stats::sum_of_negbin2</c>: p = v/m − 1, k = m/p, q = 1 + p,
        /// P(X = x) = C(k + x − 1, x)·p^x / q^(k + x) (= <c>scipy.stats.nbinom(k, 1/q)</c>). Terms follow RSAT's recurrence
        /// ln P(0) = −k·ln q, ln P(x) = ln P(x − 1) + ln p + ln(k + x − 1) − ln q − ln x; summed as in
        /// <see cref="LogPoissonRangeProbability"/>. RSAT rounds every term to 5 significant digits (<c>LogToEng</c>); the
        /// exact sum is returned here.
        /// </summary>
        /// <param name="from">Lower bound (values &lt; 0 are treated as 0).</param>
        /// <param name="to">Upper bound (inclusive).</param>
        /// <param name="mean">Mean m, finite and &gt; 0.</param>
        /// <param name="variance">Variance v, finite and &gt; m.</param>
        /// <exception cref="ArgumentOutOfRangeException">m ≤ 0, v ≤ m, or a non-finite parameter (RSAT <c>negbin2</c> checks).</exception>
        public static double LogNegativeBinomialRangeProbability(long from, long to, double mean, double variance)
        {
            if (!(mean > 0) || double.IsPositiveInfinity(mean))
                throw new ArgumentOutOfRangeException(nameof(mean), mean, "The mean must be finite and strictly positive.");
            if (!(variance > mean) || double.IsPositiveInfinity(variance))
                throw new ArgumentOutOfRangeException(nameof(variance), variance, "The variance must be finite and greater than the mean.");
            double p = variance / mean - 1;
            double k = mean / p;
            double q = 1 + p;
            double logP = Math.Log(p), logQ = Math.Log(q);
            double theta = p / q;
            double crossing = (k * theta - 1) / (1 - theta);
            long mode = crossing <= 0 ? 0 : (long)Math.Ceiling(crossing);
            return LogDiscreteRangeSum(Math.Max(0L, from), to, mode,
                x =>
                {
                    double log = -logQ * k;
                    for (long i = 1; i <= x; i++)
                        log += logP + Math.Log(k + i - 1) - logQ - Math.Log(i);
                    return log;
                },
                x => (k + x) / (x + 1.0) * theta);
        }

        // ln Σ_{x=from}^{to} P(x) for a unimodal pmf: non-decreasing below mode, non-increasing from mode on, with
        // ratioUp(x) = P(x + 1)/P(x). Each partial sum runs away from the mode so that the terms decrease.
        private static double LogDiscreteRangeSum(long from, long to, long mode, Func<long, double> logPmf, Func<long, double> ratioUp)
        {
            if (from > to) return double.NegativeInfinity;
            if (from == 0 && to == long.MaxValue) return 0.0;

            // Σ_{x=start}^{end} running upward (terms non-increasing for start ≥ mode).
            double Upward(long start, long end)
            {
                double first = logPmf(start), sum = 1.0, term = 1.0;
                for (long x = start; x < end; x++)
                {
                    term *= ratioUp(x);
                    double next = sum + term;
                    if (next == sum) break;
                    sum = next;
                }
                return first + Math.Log(sum);
            }

            // Σ_{x=end}^{start} running downward (terms non-increasing for start ≤ mode).
            double Downward(long start, long end)
            {
                double first = logPmf(start), sum = 1.0, term = 1.0;
                for (long x = start; x > end; x--)
                {
                    term /= ratioUp(x - 1);
                    double next = sum + term;
                    if (next == sum) break;
                    sum = next;
                }
                return first + Math.Log(sum);
            }

            if (from == 0)
            {
                // 1 − P(X > to), never above 1.
                if (to < mode) return Downward(to, 0);
                return Log1P(-Math.Exp(Upward(to + 1, long.MaxValue)));
            }
            if (from >= mode) return Upward(from, to);
            if (to <= mode) return Downward(to, from);

            double lower = from == 0 ? 0.0 : Math.Exp(Downward(from - 1, 0));
            double upper = to == long.MaxValue ? 0.0 : Math.Exp(Upward(to + 1, long.MaxValue));
            return Log1P(-(lower + upper));
        }

        /// <summary>
        /// ln(1 + x) without cancellation for small |x| (Goldberg 1991, "What every computer scientist should know about
        /// floating-point arithmetic", Theorem 4). <c>double.LogP1</c> evaluates ln(x + 1) directly and returns 0 for
        /// x = 1e−20; this returns 1e−20 (= numpy.log1p).
        /// </summary>
        public static double Log1P(double x)
        {
            double u = 1.0 + x;
            return u == 1.0 ? x : Math.Log(u) * x / (u - 1.0);
        }

        /// <summary>
        /// ln(eᵃ + eᵇ) without overflow (numpy <c>logaddexp</c>): max + ln(e^(a−max) + e^(b−max)); when either argument is
        /// −∞ (probability 0) the other is returned exactly, so ln 0 + ln 0 = −∞ (never NaN). As numpy's
        /// <c>npy_logaddexp</c> (numpy/_core/src/npymath/npy_math_internal.h.src), equal arguments give a + ln 2 — so
        /// logaddexp(+∞, +∞) = +∞ — and a +∞ argument gives +∞; NaN propagates.
        /// </summary>
        public static double LogAddExp(double a, double b)
        {
            if (a == b) return a + Math.Log(2.0); // handles infinities of the same sign (numpy); bit-identical for finite a
            if (double.IsNegativeInfinity(a)) return b;
            if (double.IsNegativeInfinity(b)) return a;
            double max = Math.Max(a, b);
            if (double.IsPositiveInfinity(max)) return max; // (+∞, finite): never ∞ − ∞ = NaN
            return max + Math.Log(Math.Exp(a - max) + Math.Exp(b - max));
        }

        /// <summary>
        /// eˣ − 1 without cancellation for small |x| (W. Kahan's expm1 identity). <c>double.ExpM1</c> evaluates eˣ − 1
        /// directly and returns 0 for x = 1e−20; this returns 1e−20 (= numpy.expm1).
        /// </summary>
        public static double ExpM1(double x)
        {
            double u = Math.Exp(x);
            if (u == 1.0) return x;
            double um1 = u - 1.0;
            if (um1 == -1.0 || double.IsPositiveInfinity(u)) return um1;
            return um1 * x / Math.Log(u);
        }

        /// <summary>
        /// Digamma function ψ(x) = d/dx ln Γ(x) for x &gt; 0: upward recurrence ψ(x) = ψ(x + 1) − 1/x until x ≥ 10, then
        /// the asymptotic expansion ψ(x) ~ ln x − 1/(2x) − Σ_{k≥1} B_{2k} / (2k·x^{2k}) (Abramowitz &amp; Stegun 6.3.18)
        /// through B₁₄ (next term &lt; 5e−17 at x = 10). Agrees with <c>scipy.special.digamma</c> / mpmath to ≈ 1e−15
        /// absolute; e.g. ψ(1) = −γ = −0.5772156649015329, ψ(½) = −γ − 2 ln 2.
        /// </summary>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="x"/> is not a finite positive number.</exception>
        public static double Digamma(double x)
        {
            if (!(x > 0) || double.IsPositiveInfinity(x))
                throw new ArgumentOutOfRangeException(nameof(x), x, "Digamma is implemented for finite x > 0.");

            double result = 0.0;
            while (x < 10.0)
            {
                result -= 1.0 / x;
                x += 1.0;
            }

            double inv2 = 1.0 / (x * x);
            // Σ B_{2k}/(2k) x^{-2k}, k = 1..7: 1/12, −1/120, 1/252, −1/240, 1/132, −691/32760, 1/12.
            double series = inv2 * (1.0 / 12 - inv2 * (1.0 / 120 - inv2 * (1.0 / 252 - inv2 * (1.0 / 240
                - inv2 * (1.0 / 132 - inv2 * (691.0 / 32760 - inv2 / 12))))));
            return result + Math.Log(x) - 0.5 / x - series;
        }

        /// <summary>
        /// Regularized incomplete beta function I_x(a, b) = B(x; a, b) / B(a, b) (R <c>pbeta(x, a, b)</c>,
        /// <c>scipy.special.betainc(a, b, x)</c>), for a, b &gt; 0 and x ∈ [0, 1]. Evaluated by the modified-Lentz continued
        /// fraction (Press et al., <i>Numerical Recipes</i> 3rd ed. §6.4, <c>betacf</c>) on the side
        /// x &lt; (a + 1)/(a + b + 2), else through I_x(a, b) = 1 − I_{1−x}(b, a); the prefactor xᵃ(1 − x)ᵇ/B(a, b) uses
        /// R's <c>lbeta</c> (nmath/lbeta.c: Stirling-correction form for large arguments), or TOMS 708 <c>brcomp</c>'s
        /// large-parameter form when a, b ≥ 8; the continued fraction runs in double-double arithmetic. Agrees with R 4.3.3
        /// <c>pbeta</c> to ≤ 9e−15 relative in 90 % and ≤ 1.5e−13 in all of a 300-case sweep (a, b ∈ [0.01, 3·10⁴];
        /// Evidence ONCO-ASCAT-001 §F39).
        /// </summary>
        /// <param name="x">Upper integration limit, in [0, 1].</param>
        /// <param name="a">Shape a &gt; 0 (finite).</param>
        /// <param name="b">Shape b &gt; 0 (finite).</param>
        /// <exception cref="ArgumentOutOfRangeException">x ∉ [0, 1] or a shape is not a finite positive number.</exception>
        public static double RegularizedIncompleteBeta(double x, double a, double b)
        {
            if (!(x >= 0.0 && x <= 1.0))
                throw new ArgumentOutOfRangeException(nameof(x), x, "x must lie in [0, 1].");
            ValidateBetaShape(a, nameof(a));
            ValidateBetaShape(b, nameof(b));
            return RegularizedIncompleteBetaCore(new DoubleDouble(x), DoubleDouble.Sum(1.0, -x), Math.Log(x), Log1P(-x), a, b);
        }

        /// <summary>
        /// Cumulative distribution function of Student's t distribution with <paramref name="degreesOfFreedom"/> degrees
        /// of freedom (R <c>pt(t, df)</c>, <c>scipy.stats.t.cdf</c>): with q = t²/ν, P(T ≤ t) = ½·I_{1/(1+q)}(ν/2, ½) for
        /// t &lt; 0 and 1 − ½·I_{1/(1+q)}(ν/2, ½) for t ≥ 0 (Abramowitz &amp; Stegun 26.7.1 / 26.5.27; R nmath/pt.c), with
        /// 1/(1+q) and q/(1+q) formed in double-double. Agrees with R 4.3.3 <c>pt</c> to ≤ 2e−14 relative for
        /// probabilities above 1e−30 and ≤ 2e−13 in the far tail (400-case sweep, ν ∈ [1, 3·10⁵]); for ν &gt; 4·10⁵ R switches
        /// to a normal approximation, this keeps the exact form.
        /// </summary>
        /// <param name="t">The quantile (±∞ allowed: 0 / 1).</param>
        /// <param name="degreesOfFreedom">ν &gt; 0 (finite; non-integer allowed).</param>
        /// <exception cref="ArgumentOutOfRangeException">ν is not a finite positive number, or t is NaN.</exception>
        public static double StudentTCdf(double t, double degreesOfFreedom)
        {
            if (double.IsNaN(t))
                throw new ArgumentOutOfRangeException(nameof(t), t, "t must not be NaN.");
            ValidateBetaShape(degreesOfFreedom, nameof(degreesOfFreedom));
            double halfTail = 0.5 * StudentTTwoSidedTail(t, degreesOfFreedom);
            return t < 0.0 ? halfTail : 1.0 - halfTail;
        }

        /// <summary>
        /// Two-sided one-sample Student t-test p-value of H₀: mean = <paramref name="mu"/>, exactly as R
        /// <c>t.test(x, mu = mu, alternative = "two.sided")$p.value</c> (stats/R/t.test.R): t = (x̄ − μ)/√(s²/n),
        /// ν = n − 1, p = 2·pt(−|t|, ν) = I_{ν/(ν+t²)}(ν/2, ½) (s² the n − 1 sample variance). Returns
        /// <see cref="double.NaN"/> in the two cases where R's <c>t.test</c> stops with an error: fewer than two values
        /// ("not enough 'x' observations") and "data are essentially constant" (√(s²/n) &lt; 10·ε·|x̄|, ε = 2⁻⁵²).
        /// E.g. x = {0.61, 0.64, 0.66, 0.59, 0.70}, μ = 0.6 ⇒ p = 0.10608282922154902 (R 4.3.3).
        /// </summary>
        /// <param name="values">The sample (finite values).</param>
        /// <param name="mu">Hypothesised mean (finite).</param>
        /// <exception cref="ArgumentNullException"><paramref name="values"/> is null.</exception>
        /// <exception cref="ArgumentException">a value or <paramref name="mu"/> is not finite.</exception>
        public static double OneSampleTTestPValue(IReadOnlyList<double> values, double mu)
        {
            ArgumentNullException.ThrowIfNull(values);
            if (!double.IsFinite(mu))
                throw new ArgumentException("mu must be finite.", nameof(mu));
            int n = values.Count;
            double sum = 0.0;
            for (int i = 0; i < n; i++)
            {
                if (!double.IsFinite(values[i]))
                    throw new ArgumentException("Every value must be finite.", nameof(values));
                sum += values[i];
            }

            if (n < 2)
                return double.NaN; // R: "not enough 'x' observations"

            // R mean.default: sum/n refined by the mean residual; var: two-pass Σ(x − x̄)²/(n − 1).
            double mean = sum / n;
            double residual = 0.0;
            for (int i = 0; i < n; i++) residual += values[i] - mean;
            mean += residual / n;
            double sumSq = 0.0;
            for (int i = 0; i < n; i++)
            {
                double d = values[i] - mean;
                sumSq += d * d;
            }

            double stdErr = Math.Sqrt(sumSq / (n - 1) / n);
            if (stdErr < 10.0 * MachineEpsilon * Math.Abs(mean))
                return double.NaN; // R: "data are essentially constant"

            double tStat = (mean - mu) / stdErr;
            return StudentTTwoSidedTail(tStat, n - 1);
        }

        /// <summary>
        /// Two-sided Welch (unequal-variance) two-sample t-test p-value of H₀: mean(x) = mean(y), exactly as R
        /// <c>t.test(x, y)$p.value</c> (stats/R/t.test.R, defaults <c>var.equal = FALSE</c>, <c>mu = 0</c>):
        /// sₓ = √(var(x)/nₓ), s_y = √(var(y)/n_y), s = √(sₓ² + s_y²), ν = s⁴/(sₓ⁴/(nₓ − 1) + s_y⁴/(n_y − 1)),
        /// t = (x̄ − ȳ)/s, p = 2·pt(−|t|, ν) = I_{ν/(ν+t²)}(ν/2, ½) (Welch 1947; Satterthwaite 1946). Returns
        /// <see cref="double.NaN"/> where R's <c>t.test</c> stops with an error: fewer than two values in either sample
        /// ("not enough observations") and "data are essentially constant" (s &lt; 10·ε·max(|x̄|, |ȳ|), ε = 2⁻⁵²). NaN
        /// values must be removed by the caller (R drops <c>NA</c> before testing). For ν &gt; 4·10⁵ R's <c>pt</c> switches to
        /// a normal approximation; this keeps the exact form (see <see cref="StudentTCdf"/>).
        /// E.g. x = {0.61, 0.64, 0.66, 0.59, 0.70}, y = {0.55, 0.58, 0.62, 0.57} ⇒ p = 0.042862029945299411 (R 4.3.3).
        /// </summary>
        /// <param name="x">First sample (finite values).</param>
        /// <param name="y">Second sample (finite values).</param>
        /// <exception cref="ArgumentNullException"><paramref name="x"/> or <paramref name="y"/> is null.</exception>
        /// <exception cref="ArgumentException">a value is not finite.</exception>
        public static double WelchTTestPValue(IReadOnlyList<double> x, IReadOnlyList<double> y)
        {
            ArgumentNullException.ThrowIfNull(x);
            ArgumentNullException.ThrowIfNull(y);
            (double mx, double vx) = RMeanAndVariance(x, nameof(x));
            (double my, double vy) = RMeanAndVariance(y, nameof(y));
            int nx = x.Count, ny = y.Count;
            if (nx < 2 || ny < 2)
                return double.NaN; // R: "not enough 'x' / 'y' observations"

            double stderrX = Math.Sqrt(vx / nx);
            double stderrY = Math.Sqrt(vy / ny);
            double stderr = Math.Sqrt(stderrX * stderrX + stderrY * stderrY);
            double df = Math.Pow(stderr, 4) / (Math.Pow(stderrX, 4) / (nx - 1) + Math.Pow(stderrY, 4) / (ny - 1));
            if (stderr < 10.0 * MachineEpsilon * Math.Max(Math.Abs(mx), Math.Abs(my)))
                return double.NaN; // R: "data are essentially constant"

            return StudentTTwoSidedTail((mx - my) / stderr, df);
        }

        // R mean.default (sum/n refined by the mean residual) and var (two-pass Σ(x − x̄)²/(n − 1)); NaN variance for n < 2.
        private static (double Mean, double Variance) RMeanAndVariance(IReadOnlyList<double> values, string name)
        {
            int n = values.Count;
            double sum = 0.0;
            for (int i = 0; i < n; i++)
            {
                if (!double.IsFinite(values[i]))
                    throw new ArgumentException("Every value must be finite.", name);
                sum += values[i];
            }

            if (n == 0)
                return (double.NaN, double.NaN);
            double mean = sum / n;
            double residual = 0.0;
            for (int i = 0; i < n; i++) residual += values[i] - mean;
            mean += residual / n;
            if (n < 2)
                return (mean, double.NaN);
            double sumSq = 0.0;
            for (int i = 0; i < n; i++)
            {
                double d = values[i] - mean;
                sumSq += d * d;
            }

            return (mean, sumSq / (n - 1));
        }

        private const double MachineEpsilon = 2.220446049250313e-16; // 2⁻⁵² (R .Machine$double.eps)

        // P(|T| ≥ |t|) = I_{ν/(ν+t²)}(ν/2, ½), with x = 1/(1+q), 1 − x = q/(1+q), q = t²/ν (R nmath/pt.c); x and 1 − x are
        // formed in double-double so the continued fraction sees the exact pair.
        private static double StudentTTwoSidedTail(double t, double df)
        {
            if (double.IsNaN(t))
                return double.NaN; // R: an all-zero sample tested against mu = 0 gives t = 0/0
            double r = t / Math.Sqrt(df);
            double q = r * r;
            if (double.IsPositiveInfinity(q))
                return 0.0;
            if (q == 0.0)
                return 1.0;
            DoubleDouble onePlusQ = DoubleDouble.Sum(1.0, q);
            DoubleDouble x = DoubleDouble.One / onePlusQ;
            DoubleDouble y = new DoubleDouble(q) / onePlusQ;
            return RegularizedIncompleteBetaCore(x, y, -Log1P(q), Math.Log(y.Hi), 0.5 * df, 0.5);
        }

        private static void ValidateBetaShape(double value, string name)
        {
            if (!(value > 0.0) || double.IsPositiveInfinity(value))
                throw new ArgumentOutOfRangeException(name, value, "Must be a finite positive number.");
        }

        // I_x(a, b) given x and y = 1 − x (double-double, so neither carries the other's cancellation) and ln x, ln y.
        // Continued fraction on the side where it converges fast (NR betacf), else the complement.
        private static double RegularizedIncompleteBetaCore(
            DoubleDouble x, DoubleDouble y, double logX, double logY, double a, double b)
        {
            if (x.Hi <= 0.0) return 0.0;
            if (y.Hi <= 0.0) return 1.0;
            double front = BetaPrefactor(x, y, logX, logY, a, b);
            if (x.Hi < (a + 1.0) / (a + b + 2.0))
                return front * BetaContinuedFraction(x, a, b) / a;
            // Symmetry I_x(a, b) = 1 − I_{1−x}(b, a): the shapes are swapped on purpose.
#pragma warning disable S2234
            return 1.0 - front * BetaContinuedFraction(y, b, a) / b;
#pragma warning restore S2234
        }

        // xᵃ·yᵇ / B(a, b). For a, b ≥ 8 the TOMS 708 brcomp form (Didonato & Morris 1992; R nmath/toms708.c), which
        // avoids the cancellation of a·ln x + b·ln y − ln B(a, b) between large terms; otherwise exp of that sum.
        private static double BetaPrefactor(DoubleDouble x, DoubleDouble y, double logX, double logY, double a, double b)
        {
            if (a < 8.0 || b < 8.0)
                return Math.Exp(a * logX + b * logY - LogBeta(a, b));

            double h, x0, y0, lambda;
            if (a > b)
            {
                h = b / a;
                x0 = 1.0 / (1.0 + h);
                y0 = h / (1.0 + h);
                lambda = (DoubleDouble.Sum(a, b) * y - new DoubleDouble(b)).ToDouble();
            }
            else
            {
                h = a / b;
                x0 = h / (1.0 + h);
                y0 = 1.0 / (1.0 + h);
                lambda = (new DoubleDouble(a) - DoubleDouble.Sum(a, b) * x).ToDouble();
            }

            double e = -lambda / a;
            double u = Math.Abs(e) > 0.6 ? e - (logX - Math.Log(x0)) : Rlog1(e);
            e = lambda / b;
            double v = Math.Abs(e) > 0.6 ? e - (logY - Math.Log(y0)) : Rlog1(e);
            double z = Math.Exp(-(a * u + b * v));
            double bcorr = LogGammaCorrection(a) + LogGammaCorrection(b) - LogGammaCorrection(a + b);
            return InvSqrt2Pi * Math.Sqrt(b * x0) * z * Math.Exp(-bcorr);
        }

        // rlog1(x) = x − ln(1 + x) for |x| ≤ 0.6 without cancellation: with r = x/(2 + x), ln(1 + x) = 2·atanh r, so
        // x − ln(1 + x) = r·x − 2(r³/3 + r⁵/5 + …) (r² ≤ 0.18, summed to machine precision).
        private static double Rlog1(double x)
        {
            double r = x / (2.0 + x);
            double r2 = r * r;
            double power = r * r2;
            double series = 0.0;
            for (int k = 3; k < 200; k += 2)
            {
                double term = power / k;
                series += term;
                if (Math.Abs(term) <= 1e-17 * Math.Abs(series))
                    break;
                power *= r2;
            }

            return r * x - 2.0 * series;
        }

        // Modified Lentz evaluation of the incomplete-beta continued fraction (NR 3rd ed. §6.4, betacf) in
        // double-double arithmetic: near x ≈ (a + 1)/(a + b + 2) with large shapes the recurrence 1 + aa·d cancels by
        // several orders of magnitude, which plain doubles turn into ~1e−11 relative error.
        private static double BetaContinuedFraction(DoubleDouble x, double a, double b)
        {
            const double Tiny = 1e-300;
            const int MaxIterations = 1_000_000;
            DoubleDouble one = DoubleDouble.One;
            DoubleDouble da = new(a), qab = DoubleDouble.Sum(a, b), qap = DoubleDouble.Sum(a, 1.0), qam = DoubleDouble.Sum(a, -1.0);
            DoubleDouble c = one;
            DoubleDouble d = one - qab * x / qap;
            if (Math.Abs(d.Hi) < Tiny) d = new DoubleDouble(Tiny);
            d = one / d;
            DoubleDouble h = d;
            for (int m = 1; m <= MaxIterations; m++)
            {
                double m2 = 2.0 * m;
                DoubleDouble aa = new DoubleDouble(m) * DoubleDouble.Sum(b, -m) * x / ((qam + m2) * (da + m2));
                d = one + aa * d;
                if (Math.Abs(d.Hi) < Tiny) d = new DoubleDouble(Tiny);
                c = one + aa / c;
                if (Math.Abs(c.Hi) < Tiny) c = new DoubleDouble(Tiny);
                d = one / d;
                h = h * d * c;
                aa = -(DoubleDouble.Sum(a, m) * (qab + m) * x / ((da + m2) * (qap + m2)));
                d = one + aa * d;
                if (Math.Abs(d.Hi) < Tiny) d = new DoubleDouble(Tiny);
                c = one + aa / c;
                if (Math.Abs(c.Hi) < Tiny) c = new DoubleDouble(Tiny);
                d = one / d;
                DoubleDouble del = d * c;
                h *= del;
                if (Math.Abs((del - one).Hi) <= 1e-17)
                    break;
            }

            return h.ToDouble();
        }

        // Double-double number Hi + Lo (|Lo| ≤ ½ ulp(Hi)), ~106-bit significand (Dekker 1971; Hida, Li &amp; Bailey 2001 QD).
        private readonly struct DoubleDouble
        {
            public DoubleDouble(double hi, double lo = 0.0)
            {
                Hi = hi;
                Lo = lo;
            }

            public double Hi { get; }

            public double Lo { get; }

            public static DoubleDouble One => new(1.0);

            public double ToDouble() => Hi + Lo;

            // Exact a + b.
            public static DoubleDouble Sum(double a, double b)
            {
                double s = a + b;
                double bb = s - a;
                return new DoubleDouble(s, (a - (s - bb)) + (b - bb));
            }

            private static DoubleDouble QuickSum(double a, double b)
            {
                double s = a + b;
                return new DoubleDouble(s, b - (s - a));
            }

            public static DoubleDouble operator -(DoubleDouble v) => new(-v.Hi, -v.Lo);

            public static DoubleDouble operator +(DoubleDouble l, DoubleDouble r)
            {
                DoubleDouble s = Sum(l.Hi, r.Hi);
                DoubleDouble t = Sum(l.Lo, r.Lo);
                DoubleDouble u = QuickSum(s.Hi, s.Lo + t.Hi);
                return QuickSum(u.Hi, u.Lo + t.Lo);
            }

            public static DoubleDouble operator +(DoubleDouble l, double r) => l + new DoubleDouble(r);

            public static DoubleDouble operator -(DoubleDouble l, DoubleDouble r) => l + -r;

            public static DoubleDouble operator *(DoubleDouble l, DoubleDouble r)
            {
                double p = l.Hi * r.Hi;
                double e = Math.FusedMultiplyAdd(l.Hi, r.Hi, -p);
                e += l.Hi * r.Lo + l.Lo * r.Hi;
                return QuickSum(p, e);
            }

            public static DoubleDouble operator /(DoubleDouble l, DoubleDouble r)
            {
                double q1 = l.Hi / r.Hi;
                DoubleDouble rem = l - r * new DoubleDouble(q1);
                double q2 = rem.Hi / r.Hi;
                rem -= r * new DoubleDouble(q2);
                double q3 = rem.Hi / r.Hi;
                return QuickSum(q1, q2) + new DoubleDouble(q3);
            }
        }

        // ln B(a, b) as R nmath/lbeta.c: Stirling-correction forms when an argument is ≥ 10, else ln Γ(p)Γ(q)/Γ(p+q).
        private static double LogBeta(double a, double b)
        {
            double p = Math.Min(a, b), q = Math.Max(a, b);
            if (p >= 10.0)
            {
                double corr = LogGammaCorrection(p) + LogGammaCorrection(q) - LogGammaCorrection(p + q);
                return Math.Log(q) * -0.5 + LogSqrt2Pi + corr + (p - 0.5) * Math.Log(p / (p + q)) + q * Log1P(-p / (p + q));
            }

            if (q >= 10.0)
            {
                double corr = LogGammaCorrection(q) - LogGammaCorrection(p + q);
                return LogGammaSmall(p) + corr + p - p * Math.Log(p + q) + (q - 0.5) * Log1P(-p / (p + q));
            }

            return LogGammaSmall(p) + LogGammaSmall(q) - LogGammaSmall(p + q);
        }

        // lgammacor(x) = ln Γ(x) − ((x − ½)·ln x − x + ln √(2π)) for x ≥ 10: Stirling series through B₁₄ (next term
        // < 3e−17 at x = 10).
        private static double LogGammaCorrection(double x)
        {
            double inv2 = 1.0 / (x * x);
            return (1.0 / 12 - inv2 * (1.0 / 360 - inv2 * (1.0 / 1260 - inv2 * (1.0 / 1680 - inv2 * (1.0 / 1188
                - inv2 * (691.0 / 360360 - inv2 / 156)))))) / x;
        }

        // ln Γ(x) for 0 < x < 20: shift to z = x + n ≥ 10, Stirling with lgammacor, minus ln Π_{k<n}(x + k).
        private static double LogGammaSmall(double x)
        {
            double product = 1.0;
            while (x < 10.0)
            {
                product *= x;
                x += 1.0;
            }

            return (x - 0.5) * Math.Log(x) - x + LogSqrt2Pi + LogGammaCorrection(x) - Math.Log(product);
        }

        // ln P(X = x), X ~ Binomial(n, p), 0 ≤ x ≤ n, 0 < p < 1 — Loader (2000) / R nmath dbinom_raw.
        private static double LogBinomialPmf(long x, long n, double logP, double p, double logQ)
        {
            double q = -ExpM1(logP);
            double np = n * p;
            if (!(np > 0.0) || !double.IsNormal(np) || !double.IsNormal(n * q))
            {
                // p (or q) is so small that n·p is not a normal double: exact log-factorial form, where the
                // cancellation of large terms is irrelevant because x·ln p dominates.
                return LogFactorial(n) - LogFactorial(x) - LogFactorial(n - x) + x * logP + (n - x) * logQ;
            }

            if (x == 0)
                return p < 0.1 ? -BinomialDeviance(n, n * q) - np : n * logQ;
            if (x == n)
                return q < 0.1 ? -BinomialDeviance(n, np) - n * q : n * logP;

            double lc = StirlingError(n) - StirlingError(x) - StirlingError(n - x)
                        - BinomialDeviance(x, np) - BinomialDeviance(n - x, n * q);
            double lf = Log2Pi + Math.Log(x) + Log1P(-(double)x / n);
            return lc - 0.5 * lf;
        }

        private const double Log2Pi = 1.8378770664093454835606594728112;      // ln(2π)
        private const double LogSqrt2Pi = 0.91893853320467274178032973640562;  // ln √(2π)

        // ln n! = stirlerr(n) + (n + ½)·ln n − n + ln √(2π) (definition of the Stirling error), exact for n = 0.
        private static double LogFactorial(long n)
            => n == 0 ? 0.0 : StirlingError(n) + (n + 0.5) * Math.Log(n) - n + LogSqrt2Pi;

        // stirlerr(n) = ln n! − (n + ½)·ln n + n − ln √(2π) for integer n ≥ 1 (R nmath stirlerr.c; table values for
        // n ≤ 15 re-derived with mpmath at 40 digits).
        private static readonly double[] StirlingErrorTable =
        {
            0.0,
            0.08106146679532725821967026, 0.04134069595540929409382208, 0.02767792568499833914878929,
            0.02079067210376509311152277, 0.01664469118982119216319487, 0.01387612882307074799874573,
            0.01189670994589177009505572, 0.01041126526197209649747857, 0.009255462182712732917728637,
            0.008330563433362871256469319, 0.007573675487951840794972024, 0.006942840107209529865664153,
            0.006408994188004207068439631, 0.005951370112758847735624416, 0.00555473355196280137103869,
        };

        private static double StirlingError(long n)
        {
            const double S0 = 1.0 / 12, S1 = 1.0 / 360, S2 = 1.0 / 1260, S3 = 1.0 / 1680, S4 = 1.0 / 1188;
            if (n <= 15) return StirlingErrorTable[n];
            double nd = n, nn = nd * nd;
            if (nd > 500) return (S0 - S1 / nn) / nd;
            if (nd > 80) return (S0 - (S1 - S2 / nn) / nn) / nd;
            if (nd > 35) return (S0 - (S1 - (S2 - S3 / nn) / nn) / nn) / nd;
            return (S0 - (S1 - (S2 - (S3 - S4 / nn) / nn) / nn) / nn) / nd;
        }

        // bd0(x, np) = x·ln(x/np) + np − x, evaluated by a series when x ≈ np (Loader 2000; R nmath bd0.c).
        private static double BinomialDeviance(double x, double np)
        {
            if (Math.Abs(x - np) < 0.1 * (x + np))
            {
                double v = (x - np) / (x + np);
                double s = (x - np) * v;
                if (Math.Abs(s) < 2.2250738585072014e-308) return s; // |s| < DBL_MIN
                double ej = 2 * x * v;
                v *= v;
                for (int j = 1; j < 1000; j++)
                {
                    ej *= v;
                    double s1 = s + ej / ((j << 1) + 1);
                    if (s1 == s) return s1;
                    s = s1;
                }
            }
            return x * Math.Log(x / np) + np - x;
        }
        /// <summary>
        /// R <c>sum(x)</c> on x86-64: the values are accumulated left to right in an x87 80-bit <c>long double</c>
        /// (64-bit significand, round-to-nearest-even after every addition; R <c>summary.c</c> <c>rsum</c>, also
        /// <c>colSums</c> and the <c>matprod = "internal"</c> dot product) and the total is rounded once more to double.
        /// The accumulator is emulated exactly with integer arithmetic, so the result is bit-identical to R's (± overflow
        /// to ±∞ as R). With <paramref name="skipNaN"/> (R <c>na.rm = TRUE</c>) NaN values are skipped, otherwise a NaN
        /// makes the sum NaN; +∞ and −∞ together give NaN.
        /// </summary>
        /// <param name="values">Values in summation order.</param>
        /// <param name="skipNaN">Skip NaN values (R <c>na.rm = TRUE</c>).</param>
        /// <returns>The sum rounded to double.</returns>
        public static double ExtendedPrecisionSum(ReadOnlySpan<double> values, bool skipNaN = false)
        {
            System.Numerics.BigInteger mantissa = System.Numerics.BigInteger.Zero; // value = mantissa · 2^exponent
            int exponent = 0;
            bool nan = false, positiveInfinity = false, negativeInfinity = false;
            foreach (double v in values)
            {
                if (double.IsNaN(v))
                {
                    if (!skipNaN) nan = true;
                    continue;
                }

                if (double.IsPositiveInfinity(v)) { positiveInfinity = true; continue; }
                if (double.IsNegativeInfinity(v)) { negativeInfinity = true; continue; }
                if (v == 0.0) continue;

                long bits = BitConverter.DoubleToInt64Bits(v);
                int biased = (int)((bits >> 52) & 0x7FF);
                long fraction = bits & 0xFFFFFFFFFFFFFL;
                long m = biased == 0 ? fraction : fraction | (1L << 52);
                int e = biased == 0 ? -1074 : biased - 1075;
                System.Numerics.BigInteger term = v < 0 ? -(System.Numerics.BigInteger)m : m;
                if (mantissa.IsZero)
                {
                    mantissa = term;
                    exponent = e;
                }
                else if (e >= exponent)
                {
                    mantissa += term << (e - exponent);
                }
                else
                {
                    mantissa = (mantissa << (exponent - e)) + term;
                    exponent = e;
                }

                (mantissa, exponent) = RoundToSignificantBits(mantissa, exponent, 64);
            }

            if (nan || (positiveInfinity && negativeInfinity)) return double.NaN;
            if (positiveInfinity) return double.PositiveInfinity;
            if (negativeInfinity) return double.NegativeInfinity;
            if (mantissa.IsZero) return 0.0;
            (mantissa, exponent) = RoundToSignificantBits(mantissa, exponent, 53);
            return Math.ScaleB((double)mantissa, exponent); // |mantissa| ≤ 2^53: exact; ScaleB overflows to ±∞ as R
        }

        // Rounds mantissa·2^exponent to `bits` significant bits, ties to even.
        private static (System.Numerics.BigInteger Mantissa, int Exponent) RoundToSignificantBits(
            System.Numerics.BigInteger mantissa, int exponent, int bits)
        {
            if (mantissa.IsZero) return (mantissa, 0);
            bool negative = mantissa.Sign < 0;
            System.Numerics.BigInteger magnitude = System.Numerics.BigInteger.Abs(mantissa);
            long length = magnitude.GetBitLength();
            if (length > bits)
            {
                int shift = (int)(length - bits);
                System.Numerics.BigInteger quotient = magnitude >> shift;
                System.Numerics.BigInteger remainder = magnitude - (quotient << shift);
                System.Numerics.BigInteger half = System.Numerics.BigInteger.One << (shift - 1);
                if (remainder > half || (remainder == half && !quotient.IsEven)) quotient += 1;
                magnitude = quotient;
                exponent += shift;
                if (magnitude.GetBitLength() > bits)
                {
                    magnitude >>= 1; // rounding carried to a power of two: exact
                    exponent += 1;
                }
            }

            // Drop trailing zero bits so the exponent stays bounded (does not change the value).
            while (!magnitude.IsZero && magnitude.IsEven)
            {
                int tz = (int)System.Numerics.BigInteger.TrailingZeroCount(magnitude);
                magnitude >>= tz;
                exponent += tz;
            }

            return (negative ? -magnitude : magnitude, exponent);
        }

        /// <summary>
        /// R <c>dbinom(x, size, prob, log = TRUE)</c> (R 4.3 nmath <c>dbinom.c</c> → <c>dbinom_raw(x, n, p, 1 − p)</c>:
        /// Loader's saddle-point form with <c>stirlerr</c> and <c>bd0</c>): ln P(X = x), X ~ Binomial(size, prob).
        /// NaN arguments propagate; prob ∉ [0, 1] or a negative / non-integer size give NaN; x &lt; 0, non-integer
        /// or &gt; size give −∞.
        /// </summary>
        /// <param name="x">Number of successes.</param>
        /// <param name="size">Number of trials.</param>
        /// <param name="prob">Success probability.</param>
        /// <returns>The log probability.</returns>
        public static double BinomialLogDensity(double x, double size, double prob)
        {
            if (double.IsNaN(x) || double.IsNaN(size) || double.IsNaN(prob)) return x + size + prob;
            if (prob < 0 || prob > 1 || size < 0 || Math.Abs(size - Math.Round(size)) > 1e-7 * Math.Max(1.0, Math.Abs(size)))
                return double.NaN;
            if (Math.Abs(x - Math.Round(x)) > 1e-7 * Math.Max(1.0, Math.Abs(x))) return double.NegativeInfinity; // R_D_nonint_check
            if (x < 0 || double.IsInfinity(x)) return double.NegativeInfinity;
            double n = Math.Round(size), k = Math.Round(x), p = prob, q = 1 - prob;

            if (p == 0) return k == 0 ? 0.0 : double.NegativeInfinity;
            if (q == 0) return k == n ? 0.0 : double.NegativeInfinity;
            if (k == 0)
            {
                if (n == 0) return 0.0;
                return p < 0.1 ? -BinomialDeviance(n, n * q) - n * p : n * Math.Log(q);
            }

            if (k == n) return q < 0.1 ? -BinomialDeviance(n, n * p) - n * q : n * Math.Log(p);
            if (k < 0 || k > n) return double.NegativeInfinity;
            double lc = StirlingError((long)n) - StirlingError((long)k) - StirlingError((long)(n - k))
                        - BinomialDeviance(k, n * p) - BinomialDeviance(n - k, n * q);
            double lf = Log2Pi + Math.Log(k) + Log1P(-k / n);
            return lc - 0.5 * lf;
        }

        /// <summary>
        /// Silverman's rule-of-thumb bandwidth for a Gaussian kernel density estimate, exactly as R
        /// <c>stats::bw.nrd0</c>: <c>h = 0.9 · lo · n^(−1/5)</c> with <c>lo = min(sd(x), IQR(x)/1.34)</c> (sample SD,
        /// n − 1 divisor; IQR from R's default type-7 quantiles); when that is 0, <c>lo</c> falls back to sd(x), then
        /// to |x₁| (the first value as supplied), then to 1 (R: <c>(lo &lt;- hi) || (lo &lt;- abs(x[1L])) || (lo &lt;- 1)</c>).
        /// Silverman (1986) <i>Density Estimation</i>, eq. 3.31. Agrees with R 4.3.3 <c>bw.nrd0</c> to ≈ 1e−15 relative.
        /// </summary>
        /// <param name="values">The sample (at least two finite values).</param>
        /// <exception cref="ArgumentNullException"><paramref name="values"/> is null.</exception>
        /// <exception cref="ArgumentException">fewer than two values, or a value is not finite.</exception>
        public static double BandwidthNrd0(IReadOnlyList<double> values)
        {
            ArgumentNullException.ThrowIfNull(values);
            int n = values.Count;
            if (n < 2)
            {
                throw new ArgumentException("bw.nrd0 needs at least 2 data points.", nameof(values));
            }

            double sum = 0.0;
            for (int i = 0; i < n; i++)
            {
                if (!double.IsFinite(values[i]))
                {
                    throw new ArgumentException("All values must be finite.", nameof(values));
                }

                sum += values[i];
            }

            double mean = sum / n;
            double ss = 0.0;
            for (int i = 0; i < n; i++)
            {
                double d = values[i] - mean;
                ss += d * d;
            }

            double hi = Math.Sqrt(ss / (n - 1));
            var sorted = new double[n];
            for (int i = 0; i < n; i++) sorted[i] = values[i];
            Array.Sort(sorted);
            double iqr = QuantileType7(sorted, 0.75) - QuantileType7(sorted, 0.25);

            double lo = Math.Min(hi, iqr / 1.34);
            if (lo == 0.0)
            {
                lo = hi;
                if (lo == 0.0) lo = Math.Abs(values[0]);
                if (lo == 0.0) lo = 1.0;
            }

            return 0.9 * lo * Math.Pow(n, -0.2);
        }

        /// <summary>
        /// R <c>quantile(x, probability, type = 7, names = FALSE)</c> (stats/R/quantile.R, Hyndman &amp; Fan 1996 type 7),
        /// with R's own 1-based arithmetic: <c>index = 1 + (n − 1)·p</c>, <c>lo = floor(index)</c>, <c>hi = ceiling(index)</c>,
        /// <c>qs = x[lo]</c>, and <c>(1 − h)·qs + h·x[hi]</c> (h = index − lo) only when <c>index &gt; lo</c> and
        /// <c>x[hi] ≠ qs</c>. E.g. x = {3, 1, 4, 1, 5}, p = 0.85 ⇒ 4.4 (R 4.3.3). Added for Battenberg's
        /// <c>filterMarkS4</c> (FIN-B24 F40).
        /// </summary>
        /// <param name="values">The sample (unsorted; at least one finite value).</param>
        /// <param name="probability">Probability p ∈ [0, 1].</param>
        /// <exception cref="ArgumentNullException"><paramref name="values"/> is null.</exception>
        /// <exception cref="ArgumentException">empty sample or a non-finite value.</exception>
        /// <exception cref="ArgumentOutOfRangeException">p outside [0, 1].</exception>
        public static double SampleQuantileType7(IReadOnlyList<double> values, double probability)
        {
            ArgumentNullException.ThrowIfNull(values);
            if (!(probability >= 0.0 && probability <= 1.0))
                throw new ArgumentOutOfRangeException(nameof(probability), probability, "p must lie in [0, 1].");
            int n = values.Count;
            if (n == 0)
                throw new ArgumentException("The sample must not be empty.", nameof(values));
            var x = new double[n];
            for (int i = 0; i < n; i++)
            {
                if (!double.IsFinite(values[i]))
                    throw new ArgumentException("Every value must be finite.", nameof(values));
                x[i] = values[i];
            }

            Array.Sort(x);
            double index = 1.0 + (n - 1) * probability;
            double lo = Math.Floor(index);
            double hi = Math.Ceiling(index);
            double qs = x[(int)lo - 1];
            double xHi = x[(int)hi - 1];
            if (index > lo && xHi != qs)
            {
                double h = index - lo;
                qs = ((1 - h) * qs) + (h * xHi);
            }

            return qs;
        }

        // R stats::quantile type 7 on an ascending-sorted sample: index = 1 + (n − 1)·p; interpolate
        // (1 − h)·x[lo] + h·x[hi] only when index > lo and x[hi] ≠ x[lo] (R quantile.default).
        private static double QuantileType7(double[] sorted, double probability)
        {
            int n = sorted.Length;
            double index = (n - 1) * probability; // 0-based
            int lo = (int)Math.Floor(index);
            int hi = (int)Math.Ceiling(index);
            double qs = sorted[lo];
            if (index > lo && sorted[hi] != qs)
            {
                double h = index - lo;
                qs = ((1 - h) * qs) + (h * sorted[hi]);
            }

            return qs;
        }

        /// <summary>
        /// Gaussian kernel density estimate exactly as R <c>stats::density.default(x, bw = "nrd0", adjust, kernel =
        /// "gaussian", n, cut)</c> with unit weights: bandwidth <c>bw = adjust · bw.nrd0(x)</c>; output grid
        /// <c>seq(min(x) − cut·bw, max(x) + cut·bw, length.out = n)</c>; the data are linearly binned
        /// (<c>C_BinDist</c>, weight 1/N) onto <c>N_g = max(n, 512)</c> (rounded up to a power of two above 512) points
        /// spanning <c>[from − 4·bw, to + 4·bw]</c>, convolved with the Gaussian kernel, clamped at 0 and linearly
        /// interpolated (<c>approx</c>) onto the output grid. R evaluates the circular convolution by FFT; here it is
        /// evaluated directly (the binned mass occupies only the first half of the zero-padded 2·N_g buffer, so the
        /// circular and linear sums coincide), which agrees with R to ≈ 1e−15 absolute.
        /// <para>
        /// Kernel lattice: R ≥ 4.4 (default <c>old.coords = FALSE</c>) evaluates the kernel at multiples of the bin
        /// spacing (up − lo)/(N_g − 1); R ≤ 4.3 (or <c>old.coords = TRUE</c>) used 2·(up − lo)/(2·N_g − 1), which
        /// rescales densities by ≈ 0.999 (R PR#18337). <paramref name="legacyCoordinates"/> selects the R ≤ 4.3 lattice.
        /// </para>
        /// </summary>
        /// <param name="values">The sample (at least two finite values).</param>
        /// <param name="adjust">Bandwidth multiplier (R <c>adjust</c>), &gt; 0.</param>
        /// <param name="points">Number of output grid points (R <c>n</c>, default 512), ≥ 2.</param>
        /// <param name="cut">Grid extension beyond the data in bandwidths (R <c>cut</c>, default 3), ≥ 0.</param>
        /// <param name="legacyCoordinates">true reproduces R ≤ 4.3 (<c>old.coords = TRUE</c>) values.</param>
        /// <exception cref="ArgumentNullException"><paramref name="values"/> is null.</exception>
        /// <exception cref="ArgumentException">fewer than two values or a non-finite value.</exception>
        /// <exception cref="ArgumentOutOfRangeException">invalid <paramref name="adjust"/>, <paramref name="points"/> or <paramref name="cut"/>.</exception>
        public static KernelDensityEstimate GaussianKernelDensity(
            IReadOnlyList<double> values,
            double adjust = 1.0,
            int points = 512,
            double cut = 3.0,
            bool legacyCoordinates = false)
        {
            ArgumentNullException.ThrowIfNull(values);
            if (!(adjust > 0.0) || double.IsPositiveInfinity(adjust))
                throw new ArgumentOutOfRangeException(nameof(adjust), adjust, "adjust must be a finite positive number.");
            if (points < 2)
                throw new ArgumentOutOfRangeException(nameof(points), points, "At least two grid points are required.");
            if (!(cut >= 0.0) || double.IsPositiveInfinity(cut))
                throw new ArgumentOutOfRangeException(nameof(cut), cut, "cut must be a finite non-negative number.");

            double bw = adjust * BandwidthNrd0(values);
            int nx = values.Count;
            double min = double.PositiveInfinity, max = double.NegativeInfinity;
            for (int i = 0; i < nx; i++)
            {
                min = Math.Min(min, values[i]);
                max = Math.Max(max, values[i]);
            }

            int ng = Math.Max(points, 512);
            if (ng > 512) ng = 1 << (int)Math.Ceiling(Math.Log2(ng));

            double from = min - (cut * bw);
            double to = max + (cut * bw);
            double lo = from - (4 * bw);
            double up = to + (4 * bw);

            // C_BinDist: linear binning of weight 1/N onto ng points over [lo, up].
            var binned = new double[ng];
            double weight = 1.0 / nx;
            double xdelta = (up - lo) / (ng - 1);
            int ixmax = ng - 2;
            for (int i = 0; i < nx; i++)
            {
                double xpos = (values[i] - lo) / xdelta;
                int ix = (int)Math.Floor(xpos);
                double fx = xpos - ix;
                if (ix >= 0 && ix <= ixmax)
                {
                    binned[ix] += (1 - fx) * weight;
                    binned[ix + 1] += fx * weight;
                }
                else if (ix == -1)
                {
                    binned[0] += fx * weight;
                }
                else if (ix == ixmax + 1)
                {
                    binned[ix] += (1 - fx) * weight;
                }
            }

            // Kernel ordinates kords = seq.int(0, L, length.out = 2·ng) with L = ((2ng − 1)/(ng − 1))·(up − lo)
            // (R ≥ 4.4) or 2·(up − lo) (old.coords); the lag-d kernel weight is dnorm(d·step, sd = bw).
            double span = legacyCoordinates ? 2 * (up - lo) : (2.0 * ng - 1) / (ng - 1) * (up - lo);
            double step = span / ((2 * ng) - 1);
            var kernel = new double[ng];
            for (int d = 0; d < ng; d++)
            {
                double z = d * step / bw;
                kernel[d] = InvSqrt2Pi * Math.Exp(-0.5 * z * z) / bw; // R dnorm4: M_1_SQRT_2PI·exp(−x²/2)/σ
            }

            var kords = new double[ng];
            for (int j = 0; j < ng; j++)
            {
                double acc = 0.0;
                for (int m = 0; m < ng; m++)
                {
                    if (binned[m] != 0.0) acc += binned[m] * kernel[Math.Abs(m - j)];
                }

                kords[j] = Math.Max(0.0, acc);
            }

            double[] xords = RSequence(lo, up, ng);
            double[] gridX = RSequence(from, to, points);
            var gridY = new double[points];
            for (int k = 0; k < points; k++)
            {
                gridY[k] = LinearInterpolate(xords, kords, gridX[k]);
            }

            return new KernelDensityEstimate(gridX, gridY, bw);
        }

        private const double InvSqrt2Pi = 0.398942280401432677939946059934; // R M_1_SQRT_2PI

        // R seq.int(from, to, length.out = n): from + i·((to − from)/(n − 1)), last element exactly `to`.
        private static double[] RSequence(double from, double to, int n)
        {
            var seq = new double[n];
            double by = (to - from) / (n - 1);
            seq[0] = from;
            for (int i = 1; i < n - 1; i++) seq[i] = from + (i * by);
            seq[n - 1] = to;
            return seq;
        }

        // R stats approx1 (method = "linear", rule = 1): bisection, exact knots returned verbatim.
        private static double LinearInterpolate(double[] x, double[] y, double v)
        {
            int i = 0, j = x.Length - 1;
            if (v < x[i] || v > x[j]) return double.NaN;
            while (i < j - 1)
            {
                int ij = (i + j) / 2;
                if (v < x[ij]) j = ij; else i = ij;
            }

            if (v == x[j]) return y[j];
            if (v == x[i]) return y[i];
            return y[i] + ((y[j] - y[i]) * ((v - x[i]) / (x[j] - x[i])));
        }

        /// <summary>
        /// Smooth-peak detection of R <c>peakPick::peakpick</c> (v0.11; Weber, Ramachandran &amp; Henikoff 2014,
        /// <i>Mol Cell</i> 53:819) on one series: (1) centred derivative <c>der[i] = (v[i+1] − v[i−1])/2</c>; a candidate is
        /// a point adjacent to a non-flat +→− derivative sign change with <c>|der| &lt; derivativeLimit</c>; (2) candidates
        /// within <paramref name="peakPositions"/> of either end are dropped, and a candidate survives only if
        /// <c>v[i] &gt; mean + peakMinSd · sd / √(2·npos + 1)</c> over its ±npos window; (3) while two candidates are
        /// ≤ <paramref name="neighborLimit"/> apart, the lower of the closest pair is removed (first on ties).
        /// Returns a boolean peak mask of the same length.
        /// </summary>
        /// <param name="series">The series (e.g. a density's y values).</param>
        /// <param name="neighborLimit">peakpick <c>neighlim</c> (≥ 0).</param>
        /// <param name="derivativeLimit">peakpick <c>deriv.lim</c> (default 0.04).</param>
        /// <param name="peakMinSd">peakpick <c>peak.min.sd</c> (default 0.5).</param>
        /// <param name="peakPositions">peakpick <c>peak.npos</c> (default 10).</param>
        /// <exception cref="ArgumentNullException"><paramref name="series"/> is null.</exception>
        /// <exception cref="ArgumentOutOfRangeException">negative <paramref name="neighborLimit"/> or <paramref name="peakPositions"/>.</exception>
        public static bool[] PeakPick(
            IReadOnlyList<double> series,
            int neighborLimit,
            double derivativeLimit = 0.04,
            double peakMinSd = 0.5,
            int peakPositions = 10)
        {
            ArgumentNullException.ThrowIfNull(series);
            if (neighborLimit < 0)
                throw new ArgumentOutOfRangeException(nameof(neighborLimit), neighborLimit, "neighlim must be non-negative.");
            if (peakPositions < 0)
                throw new ArgumentOutOfRangeException(nameof(peakPositions), peakPositions, "peak.npos must be non-negative.");

            int n = series.Count;
            var candidates = new bool[n];
            if (n < 3) return candidates;

            // der = rbind(NA, diff(mat, lag = 2)/2, NA); NaN plays R's NA.
            var der = new double[n];
            der[0] = double.NaN;
            der[n - 1] = double.NaN;
            for (int i = 1; i < n - 1; i++) der[i] = (series[i + 1] - series[i - 1]) / 2;

            // pos2neg[i] = der[i] ≥ 0 & der[i+1] ≤ 0 & !(der[i] == 0 & der[i+1] == 0), R three-valued logic.
            var pos2neg = new bool?[n - 1];
            for (int i = 0; i < n - 1; i++)
            {
                bool? a = Cmp(der[i], v => v >= 0), b = Cmp(der[i + 1], v => v <= 0);
                bool? bothZero = And(Cmp(der[i], v => v == 0), Cmp(der[i + 1], v => v == 0));
                pos2neg[i] = And(And(a, b), Not(bothZero));
            }

            for (int j = 0; j < n; j++)
            {
                bool? sign = Or(j < n - 1 ? pos2neg[j] : null, j > 0 ? pos2neg[j - 1] : null);
                bool? small = Cmp(der[j], v => Math.Abs(v) < derivativeLimit);
                candidates[j] = And(small, sign) == true;
            }

            // smallpeaks: drop the npos head/tail positions, then candidates not rising above mean + nsd·SEM.
            for (int j = 0; j < Math.Min(peakPositions, n); j++)
            {
                candidates[j] = false;
                candidates[n - 1 - j] = false;
            }

            var toDelete = new List<int>();
            for (int pos = 0; pos < n; pos++)
            {
                if (!candidates[pos]) continue;
                int count = (2 * peakPositions) + 1;
                double mean = 0.0;
                for (int k = pos - peakPositions; k <= pos + peakPositions; k++) mean += series[k];
                mean /= count;
                if (count == 1)
                {
                    continue; // R: sd() of one value is NA, ifelse(NA) is NA and which() drops it — kept
                }

                double ss = 0.0;
                for (int k = pos - peakPositions; k <= pos + peakPositions; k++)
                {
                    double d = series[k] - mean;
                    ss += d * d;
                }

                double sd = Math.Sqrt(ss / (count - 1));
                double limit = mean + (peakMinSd * sd / Math.Sqrt(count));
                if (!(series[pos] > limit)) toDelete.Add(pos);
            }

            foreach (int pos in toDelete) candidates[pos] = false;

            // keepmax: repeatedly drop the lower member of the closest pair while any pair is ≤ neighlim apart.
            while (true)
            {
                var positions = new List<int>();
                for (int k = 0; k < n; k++) if (candidates[k]) positions.Add(k);
                int best = -1, bestDist = int.MaxValue;
                for (int k = 0; k + 1 < positions.Count; k++)
                {
                    int dist = positions[k + 1] - positions[k];
                    if (dist < bestDist)
                    {
                        bestDist = dist;
                        best = k;
                    }
                }

                if (best < 0 || bestDist > neighborLimit) return candidates;
                int p1 = positions[best], p2 = positions[best + 1];
                candidates[series[p2] < series[p1] ? p2 : p1] = false;
            }
        }

        private static bool? Cmp(double v, Func<double, bool> predicate) => double.IsNaN(v) ? null : predicate(v);

        // R three-valued `&`: FALSE dominates, then NA.
        private static bool? And(bool? a, bool? b)
        {
            if (a == false || b == false) return false;
            if (a == true && b == true) return true;
            return null;
        }

        // R three-valued `|`: TRUE dominates, then NA.
        private static bool? Or(bool? a, bool? b)
        {
            if (a == true || b == true) return true;
            if (a == false && b == false) return false;
            return null;
        }

        private static bool? Not(bool? a) => a.HasValue ? !a.Value : null;
    }

    /// <summary>
    /// A kernel density estimate on an evenly spaced grid (R <c>density</c> object: <c>$x</c>, <c>$y</c>, <c>$bw</c>).
    /// </summary>
    /// <param name="X">Grid abscissae.</param>
    /// <param name="Y">Density ordinates at <paramref name="X"/>.</param>
    /// <param name="Bandwidth">The kernel standard deviation used (after <c>adjust</c>).</param>
    public sealed record KernelDensityEstimate(IReadOnlyList<double> X, IReadOnlyList<double> Y, double Bandwidth);
}
