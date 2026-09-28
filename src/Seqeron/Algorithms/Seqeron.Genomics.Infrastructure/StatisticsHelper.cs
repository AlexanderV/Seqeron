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
    }
}
