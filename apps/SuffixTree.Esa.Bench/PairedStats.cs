namespace SuffixTree.Esa.Bench;

/// <summary>
/// Statistics of the paired, interleaved A/B protocol (PROCESS.md amendment 2026-11-03, review MAJ-2). The statistic of a comparison is the median of the
/// per-pair ratios B/A (pairs are adjacent fresh-process runs, the order inside a pair alternates, so drift and position effects cancel). The
/// A/A calibration of the same session runs the same configuration on both sides; its ratios give the noise band. This file is also linked into
/// the test project, which pins the arithmetic on synthetic data.
/// </summary>
internal static class PairedStats
{
    /// <summary>The ratio is accepted as an improvement only if it is at most this value (3 % faster).</summary>
    internal const double AcceptRatio = 0.97;

    internal static double Median(IReadOnlyList<double> values)
    {
        if (values.Count == 0) throw new ArgumentException("no values", nameof(values));
        double[] sorted = [.. values.Order()];
        int n = sorted.Length;
        return n % 2 == 1 ? sorted[n / 2] : (sorted[(n / 2) - 1] + sorted[n / 2]) / 2;
    }

    /// <summary>Median absolute deviation from the median (not scaled).</summary>
    internal static double Mad(IReadOnlyList<double> values)
    {
        double m = Median(values);
        return Median([.. values.Select(v => Math.Abs(v - m))]);
    }

    /// <summary>The per-pair ratios B/A.</summary>
    internal static List<double> Ratios(IReadOnlyList<double> a, IReadOnlyList<double> b)
    {
        if (a.Count != b.Count) throw new ArgumentException("A and B must have the same number of pairs");
        return [.. a.Zip(b, (x, y) => y / x)];
    }

    /// <summary>
    /// The noise band <c>[1 - h, 1 + h]</c> from the A/A ratios <paramref name="aaRatios"/> of the same session:
    /// <c>h = |median - 1| + 2 * se</c>, where <c>se = 1.2533 * sigma / sqrt(N)</c> is the standard error of a median of N ratios and
    /// <c>sigma = 1.4826 * MAD</c> the robust spread of one ratio. The first term is the systematic offset of the setup (for example a
    /// persistent order effect), the second twice the sampling error of the statistic that is compared with the band.
    /// </summary>
    internal static (double Low, double High) NoiseBand(IReadOnlyList<double> aaRatios)
    {
        double bias = Math.Abs(Median(aaRatios) - 1);
        double sigma = 1.4826 * Mad(aaRatios);
        double se = 1.2533 * sigma / Math.Sqrt(aaRatios.Count);
        double h = bias + (2 * se);
        return (1 - h, 1 + h);
    }

    /// <summary>
    /// Distribution-free interval for the median from two order statistics, the widest pair with coverage of at least 95 % for the sample size
    /// (binomial, p = 1/2); for fewer than 6 values it is the full range.
    /// </summary>
    internal static (double Low, double High) MedianInterval(IReadOnlyList<double> values)
    {
        double[] sorted = [.. values.Order()];
        int n = sorted.Length;
        int k = 0; // 0-based index of the lower order statistic: the largest i with P(Bin(n, 1/2) <= i) <= 2.5 % (coverage 1 - 2 P >= 95 %)
        double cumulative = 0;
        for (int i = 0; i < n; i++)
        {
            cumulative += BinomialHalf(n, i);
            if (cumulative > 0.025) break;
            k = i;
        }
        return (sorted[k], sorted[n - 1 - k]);
    }

    /// <summary>
    /// The acceptance decision of PROCESS.md: the median B/A ratio is at most 0.97 (3 % faster), lies below the lower edge of the A/A noise band,
    /// and the 95 % interval of the median stays below 1.0 (the gain is not zero). <c>Borderline</c> marks a median within one percentage point of the
    /// bar: such a session is repeated with more pairs, whatever its verdict (a session of the S1a review series showed 2.7 % between identical binaries
    /// with 9 pairs).
    /// </summary>
    internal static (bool Accept, bool Borderline, string Reason) Decide(double medianRatio, (double Low, double High) band, (double Low, double High) interval)
    {
        bool fast = medianRatio <= AcceptRatio;
        bool outside = medianRatio < band.Low;
        bool significant = interval.High < 1.0;
        bool borderline = Math.Abs(medianRatio - AcceptRatio) <= 0.01;
        string text = $"median paired ratio {medianRatio:F4} (B/A), acceptance bar <= {AcceptRatio:F2}, A/A noise band [{band.Low:F4}, {band.High:F4}], 95 % interval of the median [{interval.Low:F4}, {interval.High:F4}]";
        string tail = borderline ? "; BORDERLINE (within 1 percentage point of the bar): repeat with more pairs" : string.Empty;
        if (fast && outside && significant) return (true, borderline, text + ": at least 3 % faster, outside the noise band, interval below 1" + tail);
        if (!fast) return (false, borderline, text + ": less than 3 % faster (or slower)" + tail);
        return (false, borderline, text + (outside ? ": the interval of the median reaches 1.0" : ": inside the noise band") + tail);
    }

    private static double BinomialHalf(int n, int k)
    {
        double logC = 0;
        for (int i = 1; i <= k; i++) logC += Math.Log((n - k + i) / (double)i);
        return Math.Exp(logC - (n * Math.Log(2)));
    }
}
