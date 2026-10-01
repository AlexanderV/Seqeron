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
    }
}
