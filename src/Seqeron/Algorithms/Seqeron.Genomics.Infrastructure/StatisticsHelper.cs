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
    }
}
