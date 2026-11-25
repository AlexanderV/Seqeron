// SPDX-License-Identifier: GPL-3.0-or-later
// Port of CNAqc 1.1.5 (caravagnalab/CNAqc, "License: GPL (>= 3) + file LICENSE"), BMix 1.0.0 (caravagnalab/BMix, "License: GPL (>= 3)") and the R stats kmeans (Hartigan-Wong, MacQueen).
// Part of Seqeron.Genomics.PurityPeakQc, licensed under the GNU General Public License, version 3 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
namespace Seqeron.Genomics.PurityPeakQc;

/// <summary>BMix binomial mixture fit with R kmeans (port of BMix 1.0.0 and R stats kmeans).</summary>
public static partial class PurityPeakAnalyzer
{
    /// <summary>One BMix Binomial-mixture fit (BMix <c>bmixfit_EM</c> result) for <c>K</c> components.</summary>
    /// <param name="Components">Number of Binomial components K.</param>
    /// <param name="Means">Component success probabilities <c>B.params</c> ("Bin 1".."Bin K").</param>
    /// <param name="MixingProportions">Mixing proportions <c>pi</c> after the E-step.</param>
    /// <param name="NegativeLogLikelihood">BMix <c>NLL</c>.</param>
    /// <param name="Bic">2·NLL + ln(N)·2K.</param>
    /// <param name="Icl">BIC + entropy of the responsibilities.</param>
    public sealed record BinomialMixtureFit(
        int Components,
        IReadOnlyList<double> Means,
        IReadOnlyList<double> MixingProportions,
        double NegativeLogLikelihood,
        double Bic,
        double Icl);

    /// <summary>
    /// Port of caravagnalab/BMix <c>bmixfit(data, K.Binomials = <paramref name="componentCounts"/>, K.BetaBinomials = 0,
    /// epsilon = 1e-8, samples = 2, score = "ICL")</c> — the Binomial mixture CNAqc fits to (NV, DP) per karyotype
    /// (<c>mixture_peak_detector</c>), reproducing R's random stream after <c>set.seed(<paramref name="seed"/>)</c>.
    /// For every K (two runs each, grid order K ascending, run 1 then 2): R <c>kmeans(NV/DP, K, nstart = 100)</c>
    /// (Hartigan–Wong, AS 136 <c>kmns.f</c>; MacQueen for K = 1; starts = <c>sample.int</c> over the distinct values),
    /// <c>sample()</c> of the centres, jitter <c>runif(K, −0.025, 0.025)</c> (redrawn until all means are in (0, 1)), mixing
    /// proportions = cluster shares; then BMix's EM loop, which stops after its first iteration because its convergence
    /// test compares the change with the initial NLL <c>.Machine$integer.max</c>: one E-step (log-sum-exp
    /// responsibilities, <c>NLL = −Σ log Σ_k π_k Binom(NV | DP, p_k)</c>) and one M-step (<c>p_k = Σ z·NV / Σ z·DP</c>,
    /// <c>π_k = Σ z / N</c>). Score: BIC = 2·NLL + ln N·2K, ICL = BIC + entropy; the run with the smallest ICL wins. A K with
    /// fewer distinct NV/DP values than K fails (R error, caught by <c>easypar</c>). R sums run in 80-bit long double
    /// (<see cref="StatisticsHelper.ExtendedPrecisionSum"/>); R's <c>%*%</c> uses the system BLAS, here the
    /// <c>matprod = "internal"</c> long-double dot product (≤ 1e−15 relative difference to OpenBLAS).
    /// </summary>
    /// <param name="successes">NV per mutation.</param>
    /// <param name="trials">DP per mutation (≥ 1, ≥ NV).</param>
    /// <param name="seed">R <c>set.seed</c> value.</param>
    /// <param name="componentCounts">K values (default CNAqc's 1:4).</param>
    /// <returns>The ICL-best fit and all fits in grid order (null where BMix fails).</returns>
    public static (BinomialMixtureFit Best, IReadOnlyList<BinomialMixtureFit?> Grid) FitBinomialMixture(
        IReadOnlyList<int> successes, IReadOnlyList<int> trials, int seed, IReadOnlyList<int>? componentCounts = null)
    {
        ArgumentNullException.ThrowIfNull(successes);
        ArgumentNullException.ThrowIfNull(trials);
        if (successes.Count != trials.Count || successes.Count == 0)
            throw new ArgumentException("successes and trials must be non-empty and of equal length.", nameof(trials));
        for (int i = 0; i < trials.Count; i++)
        {
            if (trials[i] < 1 || successes[i] < 0 || successes[i] > trials[i])
                throw new ArgumentOutOfRangeException(nameof(trials), "Each mutation needs 1 ≤ DP and 0 ≤ NV ≤ DP.");
        }

        componentCounts ??= new[] { 1, 2, 3, 4 };
        if (componentCounts.Count == 0 || componentCounts.Any(k => k < 1))
            throw new ArgumentOutOfRangeException(nameof(componentCounts), "Component counts must be ≥ 1.");
        return BMixFit(successes, trials, componentCounts, new RMersenneTwister(seed));
    }

    private static (BinomialMixtureFit Best, IReadOnlyList<BinomialMixtureFit?> Grid) BMixFit(
        IReadOnlyList<int> successes, IReadOnlyList<int> trials, IReadOnlyList<int> componentCounts, RMersenneTwister rng)
    {
        const int Samples = 2;
        var grid = new List<BinomialMixtureFit?>();
        foreach (int k in componentCounts)
        {
            for (int sample = 1; sample <= Samples; sample++) grid.Add(BMixRunner(successes, trials, k, rng));
        }

        BinomialMixtureFit? best = null;
        foreach (BinomialMixtureFit? fit in grid) // which.min(ICL): first minimum, NaN ignored
        {
            if (fit is null || double.IsNaN(fit.Icl)) continue;
            if (best is null || fit.Icl < best.Icl) best = fit;
        }

        return best is null
            ? throw new InvalidOperationException("All tasks returned error - cannot analyse this with BMix.")
            : (best, grid);
    }

    // BMix runner(): up to 14 attempts of bmixfit_EM; an error retries, the 15th "attempt" propagates it.
    private static BinomialMixtureFit? BMixRunner(IReadOnlyList<int> successes, IReadOnlyList<int> trials, int k, RMersenneTwister rng)
    {
        for (int attempt = 0; attempt < 14; attempt++)
        {
            BinomialMixtureFit? fit = BMixEm(successes, trials, k, rng);
            if (fit is not null) return fit;
        }

        return null;
    }

    // bmixfit_EM(data, K = c(k, 0), epsilon = 1e-8); null = R error.
    private static BinomialMixtureFit? BMixEm(IReadOnlyList<int> successes, IReadOnlyList<int> trials, int k, RMersenneTwister rng)
    {
        int n = successes.Count;
        var frequencies = new double[n];
        for (int i = 0; i < n; i++) frequencies[i] = (double)successes[i] / trials[i];

        (double[] centers, int[] cluster)? km = RKMeans(frequencies, k, 100, rng);
        if (km is null) return null;
        (double[] kmCenters, int[] kmCluster) = km.Value;

        // centres = sample(km$centers[, 1]): a permutation (K = 1: sample.int(1) still draws once).
        int[] order = k == 1 ? new[] { rng.UnifIndex(1) } : rng.SampleWithoutReplacement(k, k);
        var centres = order.Select(o => kmCenters[o]).ToArray();

        var means = new double[k];
        while (true)
        {
            for (int j = 0; j < k; j++) means[j] = centres[j] + RUnif(rng, -0.025, 0.025);
            if (means.All(b => b > 0 && b < 1)) break;
        }

        // pi = table(km$cluster) / nrow(data), reordered as the centres.
        var pi = new double[k];
        for (int j = 0; j < k; j++) pi[j] = (double)kmCluster.Count(c => c == order[j] + 1) / n;

        // E-step.
        var z = new double[n, k];
        var normaliser = new double[n];
        var row = new double[k];
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < k; j++)
                z[i, j] = StatisticsHelper.BinomialLogDensity(successes[i], trials[i], means[j]) + Math.Log(pi[j]);
            double offset = double.NegativeInfinity;
            for (int j = 0; j < k; j++) offset = Math.Max(offset, z[i, j]); // max(): NaN-propagating in R, finite here
            for (int j = 0; j < k; j++) row[j] = Math.Exp(z[i, j] - offset);
            double s = Math.Log(StatisticsHelper.ExtendedPrecisionSum(row)) + offset;
            normaliser[i] = double.IsFinite(s) ? s : offset;
        }

        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < k; j++) z[i, j] = Math.Exp(z[i, j] - normaliser[i]);
        }

        double nll = -StatisticsHelper.ExtendedPrecisionSum(normaliser);
        var column = new double[n];
        var products = new double[n];
        for (int j = 0; j < k; j++)
        {
            for (int i = 0; i < n; i++) column[i] = z[i, j];
            pi[j] = StatisticsHelper.ExtendedPrecisionSum(column) / n; // colSums (long double)

            // M-step: B[k] = z[, k] %*% NV / z[, k] %*% DP.
            for (int i = 0; i < n; i++) products[i] = column[i] * successes[i];
            double numerator = StatisticsHelper.ExtendedPrecisionSum(products);
            for (int i = 0; i < n; i++) products[i] = column[i] * trials[i];
            means[j] = numerator / StatisticsHelper.ExtendedPrecisionSum(products);
        }

        // The EM stops here: epsilon.conv = NLL − .Machine$integer.max < epsilon.
        int parameters = k + k;
        double bic = (2 * nll) + (Math.Log(n) * parameters);
        var entropyTerms = new double[n * k];
        for (int j = 0; j < k; j++)
        {
            for (int i = 0; i < n; i++) entropyTerms[(j * n) + i] = z[i, j] * Math.Log(z[i, j]); // column-major, 0·log 0 = NaN
        }

        double entropy = -StatisticsHelper.ExtendedPrecisionSum(entropyTerms, skipNaN: true);
        return new BinomialMixtureFit(k, means, pi, nll, bic, bic + entropy);
    }

    // R runif(1, a, b) (nmath/runif.c): a + (b − a)·u.
    private static double RUnif(RMersenneTwister rng, double a, double b)
    {
        double u;
        do
        {
            u = rng.UnifRand();
        }
        while (u <= 0 || u >= 1);
        return a + ((b - a) * u);
    }

    // R stats::kmeans(x, centers = k, nstart, algorithm = "Hartigan-Wong", iter.max = 10) for one-column x;
    // null = R error (fewer distinct values than k; kmns IFAULT 1 or 3). Cluster labels are 1-based.
    private static (double[] Centers, int[] Cluster)? RKMeans(double[] x, int k, int nstart, RMersenneTwister rng)
    {
        var distinct = new List<double>(); // unique(x): first occurrences
        var seen = new HashSet<double>();
        foreach (double v in x)
        {
            if (seen.Add(v)) distinct.Add(v);
        }

        if (distinct.Count < k) return null; // "more cluster centers than distinct data points."
        const int IterMax = 10;
        (double[] Centers, int[] Cluster, double[] Wss, int Fault) DoOne()
        {
            int[] pick = rng.SampleWithoutReplacement(distinct.Count, k);
            double[] start = pick.Select(i => distinct[i]).ToArray();
            return k == 1 ? KMeansMacQueen(x, start, IterMax) : KMeansHartiganWong(x, start, IterMax);
        }

        var best = DoOne();
        if (best.Fault == 1 || best.Fault == 3) return null;
        double bestWss = StatisticsHelper.ExtendedPrecisionSum(best.Wss);
        for (int s = 2; s <= nstart; s++)
        {
            var candidate = DoOne();
            if (candidate.Fault == 1 || candidate.Fault == 3) return null;
            double wss = StatisticsHelper.ExtendedPrecisionSum(candidate.Wss);
            if (wss < bestWss)
            {
                best = candidate;
                bestWss = wss;
            }
        }

        return (best.Centers, best.Cluster);
    }

    // R kmeans.c kmeans_MacQueen (one column).
    private static (double[] Centers, int[] Cluster, double[] Wss, int Fault) KMeansMacQueen(double[] x, double[] centers, int maxIter)
    {
        int n = x.Length, k = centers.Length;
        var cen = (double[])centers.Clone();
        var cl = new int[n];
        var nc = new int[k];
        for (int i = 0; i < n; i++)
        {
            double best = double.PositiveInfinity;
            int inew = 0;
            for (int j = 0; j < k; j++)
            {
                double t = x[i] - cen[j];
                double dd = 0.0 + (t * t);
                if (dd < best)
                {
                    best = dd;
                    inew = j + 1;
                }
            }

            cl[i] = inew;
        }

        Array.Clear(cen);
        for (int i = 0; i < n; i++)
        {
            nc[cl[i] - 1]++;
            cen[cl[i] - 1] += x[i];
        }

        for (int j = 0; j < k; j++) cen[j] /= nc[j];
        for (int iter = 0; iter < maxIter; iter++)
        {
            bool updated = false;
            for (int i = 0; i < n; i++)
            {
                double best = double.PositiveInfinity;
                int inew = 0;
                for (int j = 0; j < k; j++)
                {
                    double t = x[i] - cen[j];
                    double dd = 0.0 + (t * t);
                    if (dd < best)
                    {
                        best = dd;
                        inew = j;
                    }
                }

                int iold = cl[i] - 1;
                if (iold != inew)
                {
                    updated = true;
                    cl[i] = inew + 1;
                    nc[iold]--;
                    nc[inew]++;
                    cen[iold] += (cen[iold] - x[i]) / nc[iold];
                    cen[inew] += (x[i] - cen[inew]) / nc[inew];
                }
            }

            if (!updated) break;
        }

        var wss = new double[k];
        for (int i = 0; i < n; i++)
        {
            double t = x[i] - cen[cl[i] - 1];
            wss[cl[i] - 1] += t * t;
        }

        return (cen, cl, wss, 0);
    }

    // AS 136 (Hartigan & Wong 1979) as R stats kmns.f (one column, 1-based arrays as the Fortran).
    private static (double[] Centers, int[] Cluster, double[] Wss, int Fault) KMeansHartiganWong(double[] x, double[] centers, int iterMax)
    {
        int m = x.Length, k = centers.Length;
        double big = (double)1.0e30f; // DATA BIG /1.E30/: a REAL constant
        var a = new double[m + 1];
        for (int i = 0; i < m; i++) a[i + 1] = x[i];
        var c = new double[k + 1];
        for (int l = 0; l < k; l++) c[l + 1] = centers[l];
        var ic1 = new int[m + 1];
        var ic2 = new int[m + 1];
        var nc = new int[k + 1];
        var ncp = new int[k + 1];
        var itran = new int[k + 1];
        var live = new int[k + 1];
        var an1 = new double[k + 1];
        var an2 = new double[k + 1];
        var d = new double[m + 1];
        var wss = new double[k + 1];
        int maxQtr = (int)Math.Min(int.MaxValue, 50L * m);
        if (k <= 1 || k >= m) return (centers, new int[m], new double[k], 3);
        int fault = 0;

        for (int i = 1; i <= m; i++)
        {
            ic1[i] = 1;
            ic2[i] = 2;
            double da1 = a[i] - c[1], da2 = a[i] - c[2];
            double dt1 = 0.0 + (da1 * da1), dt2 = 0.0 + (da2 * da2);
            if (dt1 > dt2)
            {
                ic1[i] = 2;
                ic2[i] = 1;
                (dt1, dt2) = (dt2, dt1);
            }

            for (int l = 3; l <= k; l++)
            {
                double dc = a[i] - c[l];
                double db = 0.0 + (dc * dc);
                if (db >= dt2) continue;
                if (db >= dt1)
                {
                    dt2 = db;
                    ic2[i] = l;
                }
                else
                {
                    dt2 = dt1;
                    ic2[i] = ic1[i];
                    dt1 = db;
                    ic1[i] = l;
                }
            }
        }

        for (int l = 1; l <= k; l++)
        {
            nc[l] = 0;
            c[l] = 0.0;
        }

        for (int i = 1; i <= m; i++)
        {
            nc[ic1[i]]++;
            c[ic1[i]] += a[i];
        }

        for (int l = 1; l <= k; l++)
        {
            if (nc[l] == 0) return (centers, new int[m], new double[k], 1);
            double aa = nc[l];
            c[l] /= aa;
            an2[l] = aa / (aa + 1.0);
            an1[l] = big;
            if (aa > 1.0) an1[l] = aa / (aa - 1.0);
            itran[l] = 1;
            ncp[l] = -1;
        }

        int indx = 0;
        bool converged = false;
        for (int ij = 1; ij <= iterMax; ij++)
        {
            KmnsOptimalTransfer(a, m, c, k, ic1, ic2, nc, an1, an2, ncp, d, itran, live, ref indx, big);
            if (indx == m)
            {
                converged = true;
                break;
            }

            KmnsQuickTransfer(a, m, c, k, ic1, ic2, nc, an1, an2, ncp, d, itran, ref indx, ref maxQtr, big);
            if (maxQtr < 0)
            {
                fault = 4;
                converged = true;
                break;
            }

            if (k == 2)
            {
                converged = true;
                break;
            }

            for (int l = 1; l <= k; l++) ncp[l] = 0;
        }

        if (!converged) fault = 2;

        for (int l = 1; l <= k; l++)
        {
            wss[l] = 0.0;
            c[l] = 0.0;
        }

        for (int i = 1; i <= m; i++) c[ic1[i]] += a[i];
        for (int l = 1; l <= k; l++) c[l] /= nc[l];
        for (int i = 1; i <= m; i++)
        {
            double da = a[i] - c[ic1[i]];
            wss[ic1[i]] += da * da;
        }

        var cluster = new int[m];
        for (int i = 0; i < m; i++) cluster[i] = ic1[i + 1];
        return (c.Skip(1).ToArray(), cluster, wss.Skip(1).ToArray(), fault);
    }

    // AS 136.1 OPTRA: optimal-transfer stage.
    private static void KmnsOptimalTransfer(double[] a, int m, double[] c, int k, int[] ic1, int[] ic2, int[] nc,
        double[] an1, double[] an2, int[] ncp, double[] d, int[] itran, int[] live, ref int indx, double big)
    {
        for (int l = 1; l <= k; l++)
        {
            if (itran[l] == 1) live[l] = m + 1;
        }

        for (int i = 1; i <= m; i++)
        {
            indx++;
            int l1 = ic1[i], l2 = ic2[i], ll = l2;
            if (nc[l1] != 1)
            {
                if (ncp[l1] != 0)
                {
                    double df = a[i] - c[l1];
                    d[i] = (0.0 + (df * df)) * an1[l1];
                }

                double db = a[i] - c[l2];
                double r2 = (0.0 + (db * db)) * an2[l2];
                for (int l = 1; l <= k; l++)
                {
                    if ((i >= live[l1] && i >= live[l]) || l == l1 || l == ll) continue;
                    double rr = r2 / an2[l];
                    double dd = a[i] - c[l];
                    double dc = 0.0 + (dd * dd);
                    if (dc >= rr) continue;
                    r2 = dc * an2[l];
                    l2 = l;
                }

                if (r2 >= d[i])
                {
                    ic2[i] = l2;
                }
                else
                {
                    indx = 0;
                    live[l1] = m + i;
                    live[l2] = m + i;
                    ncp[l1] = i;
                    ncp[l2] = i;
                    KmnsTransfer(a[i], c, nc, an1, an2, l1, l2, big);
                    ic1[i] = l2;
                    ic2[i] = l1;
                }
            }

            if (indx == m) return;
        }

        for (int l = 1; l <= k; l++)
        {
            itran[l] = 0;
            live[l] -= m;
        }
    }

    // AS 136.2 QTRAN: quick-transfer stage.
    private static void KmnsQuickTransfer(double[] a, int m, double[] c, int k, int[] ic1, int[] ic2, int[] nc,
        double[] an1, double[] an2, int[] ncp, double[] d, int[] itran, ref int indx, ref int maxQtr, double big)
    {
        int icoun = 0, istep = 0;
        while (true)
        {
            for (int i = 1; i <= m; i++)
            {
                icoun++;
                istep++;
                if (istep >= maxQtr)
                {
                    maxQtr = -1;
                    return;
                }

                int l1 = ic1[i], l2 = ic2[i];
                if (nc[l1] != 1)
                {
                    if (istep <= ncp[l1])
                    {
                        double db = a[i] - c[l1];
                        d[i] = (0.0 + (db * db)) * an1[l1];
                    }

                    if (istep < ncp[l1] || istep < ncp[l2])
                    {
                        double r2 = d[i] / an2[l2];
                        double de = a[i] - c[l2];
                        double dd = 0.0 + (de * de);
                        if (dd < r2)
                        {
                            icoun = 0;
                            indx = 0;
                            itran[l1] = 1;
                            itran[l2] = 1;
                            ncp[l1] = istep + m;
                            ncp[l2] = istep + m;
                            KmnsTransfer(a[i], c, nc, an1, an2, l1, l2, big);
                            ic1[i] = l2;
                            ic2[i] = l1;
                        }
                    }
                }

                if (icoun == m) return;
            }
        }
    }

    // Move one point from cluster l1 to l2 (shared OPTRA/QTRAN update of centres, counts, AN1, AN2).
    private static void KmnsTransfer(double ai, double[] c, int[] nc, double[] an1, double[] an2, int l1, int l2, double big)
    {
        double al1 = nc[l1], alw = al1 - 1.0, al2 = nc[l2], alt = al2 + 1.0;
        c[l1] = ((c[l1] * al1) - ai) / alw;
        c[l2] = ((c[l2] * al2) + ai) / alt;
        nc[l1]--;
        nc[l2]++;
        an2[l1] = alw / al1;
        an1[l1] = big;
        if (alw > 1.0) an1[l1] = alw / (alw - 1.0);
        an1[l2] = alt / al2;
        an2[l2] = alt / (alt + 1.0);
    }
}
