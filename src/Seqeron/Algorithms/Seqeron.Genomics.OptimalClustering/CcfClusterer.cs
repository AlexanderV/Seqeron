// SPDX-License-Identifier: LGPL-3.0-or-later
// Port of the R package Ckmeans.1d.dp 4.3.6 (Wang & Song 2011), LGPL (>= 3).
// Part of Seqeron.Genomics.OptimalClustering, licensed under the GNU Lesser General Public License, version 3 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using System.Text;
using Seqeron.Genomics.Oncology;
using static Seqeron.Genomics.Oncology.OncologyAnalyzer;

namespace Seqeron.Genomics.OptimalClustering;

/// <summary>
/// Optimal one-dimensional k-means (Ckmeans.1d.dp 4.3.6, Wang &amp; Song 2011) of cancer cell fractions: the exact partition minimising the within-cluster sum of squares, with the BIC choice of the number of clusters.
/// Moved out of the MIT <c>OncologyAnalyzer</c> (B09 WP80): the port follows the source of the R package (LGPL-3.0-or-later).
/// </summary>
public static partial class CcfClusterer
{

    /// <summary>
    /// Upper bound on the number of cells (clusters × values) of the Ckmeans.1d.dp backtrack matrix J that
    /// <see cref="ClusterCcfValues"/> allocates (4 bytes per cell, i.e. at most 400 MB).
    /// </summary>
    private const long MaxCcfClusteringDpCells = 100_000_000;

    /// <summary>
    /// Clusters cancer cell fractions into <paramref name="clusterCount"/> clones/subclones by <b>optimal</b>
    /// one-dimensional k-means: the partition minimising the within-cluster sum of squares
    /// Σ_j Σ_{x∈S_j} (x − μ_j)² (the k-means objective, Lloyd 1982), found exactly by the dynamic program of
    /// Wang &amp; Song (2011), <i>The R Journal</i> 3(2):29–33 (R package Ckmeans.1d.dp). This is a line-by-line port
    /// of Ckmeans.1d.dp 4.3.6 <c>EWL2::fill_dp_matrix</c> (median-shifted prefix sums, R's default
    /// <c>method = "linear"</c> SMAWK row fill <c>fill_row_q_SMAWK</c>, whose tie-breaking selects the same split as
    /// R among equal-WCSS optima) and <c>backtrack</c>, i.e. it reproduces R <c>Ckmeans.1d.dp(x, k)</c> (unweighted,
    /// L2 criterion — the R defaults <c>y = 1</c>, <c>criterion "L2"</c>); unlike Lloyd iterations it cannot stop in a local
    /// optimum and needs no seeding, so the result is deterministic and independent of input order. As in
    /// Ckmeans.1d.dp, when the input has fewer distinct values than <paramref name="clusterCount"/> the number of
    /// clusters is reduced to the number of distinct values (every returned cluster is non-empty). The clonal
    /// cluster is the one with the highest centroid (Tarabichi et al. 2021, <i>Nat. Methods</i> 18:144–155: "the
    /// cluster with the highest CP can be deemed clonal").
    /// </summary>
    /// <param name="ccfValues">Cancer cell fractions to cluster (each finite).</param>
    /// <param name="clusterCount">Number of clusters k, in [1, count of values].</param>
    /// <returns>Ascending centroids (min(k, number of distinct values) of them, each the mean of a non-empty
    /// cluster), per-value cluster assignments (input order), and the clonal cluster index (the last centroid).</returns>
    /// <exception cref="ArgumentNullException"><paramref name="ccfValues"/> is null.</exception>
    /// <exception cref="ArgumentException">no values are supplied, or a value is NaN/infinite.</exception>
    /// <exception cref="ArgumentOutOfRangeException">clusterCount ∉ [1, count], or the dynamic-programming
    /// matrix (effective k × count) would exceed 10⁸ cells.</exception>
    public static CcfClustering ClusterCcfValues(IReadOnlyList<double> ccfValues, int clusterCount)
    {
        (double[] x, int[] order, int distinct) = CkmeansSortInput(ccfValues);
        int n = x.Length;

        if (clusterCount < 1 || clusterCount > n)
        {
            throw new ArgumentOutOfRangeException(
                nameof(clusterCount), clusterCount, $"Cluster count must be in [1, {n}].");
        }

        // Ckmeans.1d.dp: Kmax = min(k, number of unique values).
        int k = Math.Min(clusterCount, distinct);
        if (k > 1 && k < n)
        {
            CkmeansCheckDpSize(k, n, nameof(clusterCount), clusterCount);
        }

        return CkmeansBuildClustering(x, order, CkmeansClusterStarts(x, k));
    }

    /// <summary>
    /// Clusters cancer cell fractions by optimal one-dimensional k-means (Ckmeans.1d.dp, Wang &amp; Song 2011) with
    /// the number of clusters <b>chosen automatically</b> in [<paramref name="minClusters"/>,
    /// <paramref name="maxClusters"/>] by the Bayesian information criterion — R
    /// <c>Ckmeans.1d.dp(x, k = c(minClusters, maxClusters))</c> with its defaults (<c>method = "linear"</c>,
    /// <c>estimate.k = "BIC"</c>, unweighted). Port of Ckmeans.1d.dp 4.3.6: the R wrapper's bound adjustment
    /// (<c>cluster.1d.dp</c>: if the number of distinct values is below k_min both bounds become that number,
    /// otherwise k_max is capped at it), one dynamic program for K = k_max (<c>EWL2::fill_dp_matrix</c>, SMAWK), and
    /// <c>select_levels</c>: for each K, a Gaussian mixture with the K clusters' proportions λ = n_k/N, means and
    /// sample variances (a zero variance becomes d²/36 and a singleton's variance d², d = gap to the nearest
    /// neighbouring cluster's value), log-likelihood Σ_i log Σ_k λ_k φ(x_i; μ_k, σ²_k), and
    /// BIC = 2·logL − (3K − 1)·ln N; the first K with the largest BIC wins (strict <c>&gt;</c>, ties keep the
    /// smaller K). All values equal → one cluster.
    /// </summary>
    /// <param name="ccfValues">Cancer cell fractions to cluster (each finite).</param>
    /// <param name="minClusters">Smallest number of clusters considered (≥ 1; R <c>k[1]</c>).</param>
    /// <param name="maxClusters">Largest number of clusters considered (≥ <paramref name="minClusters"/>; may exceed
    /// the number of values — as in R it is capped at the number of distinct values).</param>
    /// <returns>The clustering for the selected number of clusters; the chosen k is
    /// <c>Centroids.Count</c>. Centroids ascend; the clonal cluster is the last.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="ccfValues"/> is null.</exception>
    /// <exception cref="ArgumentException">no values are supplied, or a value is NaN/infinite.</exception>
    /// <exception cref="ArgumentOutOfRangeException">minClusters &lt; 1, maxClusters &lt; minClusters, or the
    /// dynamic-programming matrix (effective k_max × count) would exceed 10⁸ cells.</exception>
    public static CcfClustering ClusterCcfValues(IReadOnlyList<double> ccfValues, int minClusters, int maxClusters) =>
        ClusterCcfValuesSelectingK(ccfValues, minClusters, maxClusters).Clustering;

    /// <summary>
    /// <see cref="ClusterCcfValues(IReadOnlyList{double}, int, int)"/> that also returns R's <c>$BIC</c> vector
    /// (one value per K in the adjusted [k_min, k_max]; all zero when every value is equal, as in Ckmeans.1d.dp).
    /// </summary>
    internal static (CcfClustering Clustering, double[] Bic) ClusterCcfValuesSelectingK(
        IReadOnlyList<double> ccfValues, int minClusters, int maxClusters)
    {
        (double[] x, int[] order, int distinct) = CkmeansSortInput(ccfValues);
        int n = x.Length;

        if (minClusters < 1)
        {
            throw new ArgumentOutOfRangeException(nameof(minClusters), minClusters, "Minimum cluster count must be ≥ 1.");
        }

        if (maxClusters < minClusters)
        {
            throw new ArgumentOutOfRangeException(
                nameof(maxClusters), maxClusters, $"Maximum cluster count must be ≥ minClusters ({minClusters}).");
        }

        // R cluster.1d.dp bound adjustment to the number of unique values.
        int kMin = minClusters;
        int kMax = maxClusters;
        if (distinct < kMin)
        {
            kMin = distinct;
            kMax = distinct;
        }
        else if (distinct < kMax)
        {
            kMax = distinct;
        }

        double[] bic = new double[kMax - kMin + 1];
        if (distinct == 1)
        {
            // kmeans_1d_dp: all elements equal → a single cluster; BIC is left at 0.
            return (CkmeansBuildClustering(x, order, new int[1]), bic);
        }

        if (kMax > 1)
        {
            CkmeansCheckDpSize(kMax, n, nameof(maxClusters), maxClusters);
        }

        int[][] j = CkmeansFillDpMatrix(x, kMax);
        int kOpt = CkmeansSelectLevels(x, j, kMin, kMax, bic);
        return (CkmeansBuildClustering(x, order, CkmeansBacktrackStarts(j, kOpt)), bic);
    }

    /// <summary>
    /// Validates the CCF values and returns them sorted ascending (stable; Ckmeans.1d.dp works on sorted data),
    /// the original index of each sorted value, and the number of distinct values.
    /// </summary>
    private static (double[] Sorted, int[] Order, int Distinct) CkmeansSortInput(IReadOnlyList<double> ccfValues)
    {
        ArgumentNullException.ThrowIfNull(ccfValues);

        int n = ccfValues.Count;
        if (n == 0)
        {
            throw new ArgumentException("At least one CCF value is required.", nameof(ccfValues));
        }

        for (int i = 0; i < n; i++)
        {
            if (double.IsNaN(ccfValues[i]) || double.IsInfinity(ccfValues[i]))
            {
                throw new ArgumentException($"CCF value must be finite; got {ccfValues[i]} at index {i}.", nameof(ccfValues));
            }
        }

        int[] order = Enumerable.Range(0, n).OrderBy(i => ccfValues[i]).ToArray();
        double[] x = new double[n];
        for (int s = 0; s < n; s++)
        {
            x[s] = ccfValues[order[s]];
        }

        int distinct = 1;
        for (int s = 1; s < n; s++)
        {
            if (x[s] != x[s - 1])
            {
                distinct++;
            }
        }

        return (x, order, distinct);
    }

    private static void CkmeansCheckDpSize(int k, int n, string paramName, int requested)
    {
        if ((long)k * n > MaxCcfClusteringDpCells)
        {
            throw new ArgumentOutOfRangeException(
                paramName, requested,
                $"Optimal 1-D k-means needs a {k} × {n} matrix, above the {MaxCcfClusteringDpCells}-cell limit.");
        }
    }

    /// <summary>
    /// Ckmeans.1d.dp <c>backtrack</c> output: centre = arithmetic mean of each contiguous sorted block, labels mapped
    /// back to input order; the clonal cluster (highest centroid) is the last block.
    /// </summary>
    private static CcfClustering CkmeansBuildClustering(double[] x, int[] order, int[] clusterStart)
    {
        int n = x.Length;
        int k = clusterStart.Length;
        double[] centroids = new double[k];
        int[] sortedCluster = new int[n];
        for (int q = 0; q < k; q++)
        {
            int left = clusterStart[q];
            int right = q + 1 < k ? clusterStart[q + 1] - 1 : n - 1;
            double sum = 0.0;
            for (int s = left; s <= right; s++)
            {
                sum += x[s];
                sortedCluster[s] = q;
            }

            centroids[q] = sum / (right - left + 1);
        }

        int[] assignments = new int[n];
        for (int s = 0; s < n; s++)
        {
            assignments[order[s]] = sortedCluster[s];
        }

        return new CcfClustering(centroids, assignments, k - 1);
    }

    /// <summary>
    /// Ckmeans.1d.dp 4.3.6 <c>select_levels</c> (<c>src/select_levels.cpp</c>, R default <c>estimate.k = "BIC"</c>):
    /// returns the K in [<paramref name="kMin"/>, <paramref name="kMax"/>] maximising the Gaussian-mixture BIC and
    /// writes each BIC to <paramref name="bic"/>[K − kMin].
    /// </summary>
    private static int CkmeansSelectLevels(double[] x, int[][] j, int kMin, int kMax, double[] bic)
    {
        int n = x.Length;
        if (kMin > kMax || n < 2)
        {
            return Math.Min(kMin, kMax);
        }

        int kOpt = kMin;
        double maxBic = 0.0;
        double[] lambda = new double[kMax];
        double[] mu = new double[kMax];
        double[] sigma2 = new double[kMax];
        double[] coeff = new double[kMax];

        for (int kk = kMin; kk <= kMax; kk++)
        {
            int[] starts = CkmeansBacktrackStarts(j, kk);
            int indexLeft = 0;
            for (int k = 0; k < kk; k++)
            {
                int size = (k + 1 < kk ? starts[k + 1] : n) - starts[k];
                lambda[k] = size / (double)n;
                int indexRight = indexLeft + size - 1;
                (mu[k], sigma2[k]) = CkmeansShiftedDataVariance(x, indexLeft, indexRight);

                if (sigma2[k] == 0 || size == 1)
                {
                    double dmin;
                    if (indexLeft > 0 && indexRight < n - 1)
                    {
                        dmin = Math.Min(x[indexLeft] - x[indexLeft - 1], x[indexRight + 1] - x[indexRight]);
                    }
                    else if (indexLeft > 0)
                    {
                        dmin = x[indexLeft] - x[indexLeft - 1];
                    }
                    else
                    {
                        dmin = x[indexRight + 1] - x[indexRight];
                    }

                    if (sigma2[k] == 0)
                    {
                        sigma2[k] = dmin * dmin / 4.0 / 9.0;
                    }

                    if (size == 1)
                    {
                        sigma2[k] = dmin * dmin;
                    }
                }

                coeff[k] = lambda[k] / Math.Sqrt(2.0 * Math.PI * sigma2[k]);
                indexLeft = indexRight + 1;
            }

            double logLikelihood = 0;
            for (int i = 0; i < n; i++)
            {
                double l = 0;
                for (int k = 0; k < kk; k++)
                {
                    l += coeff[k] * Math.Exp(-(x[i] - mu[k]) * (x[i] - mu[k]) / (2.0 * sigma2[k]));
                }

                logLikelihood += Math.Log(l);
            }

            double b = 2 * logLikelihood - (3 * kk - 1) * Math.Log(n);
            bic[kk - kMin] = b;

            if (kk == kMin)
            {
                maxBic = b;
                kOpt = kMin;
            }
            else if (b > maxBic)
            {
                maxBic = b;
                kOpt = kk;
            }
        }

        return kOpt;
    }

    /// <summary>
    /// Ckmeans.1d.dp <c>shifted_data_variance</c>: mean and unbiased sample variance of sorted x[left..right],
    /// computed around the block median for stability (variance 0 for a single value).
    /// </summary>
    private static (double Mean, double Variance) CkmeansShiftedDataVariance(double[] x, int left, int right)
    {
        double sum = 0.0;
        double sumsq = 0.0;
        double mean = 0.0;
        double variance = 0.0;
        int n = right - left + 1;

        if (right >= left)
        {
            double median = x[(left + right) / 2];
            for (int i = left; i <= right; ++i)
            {
                sum += x[i] - median;
                sumsq += (x[i] - median) * (x[i] - median);
            }

            mean = sum / n + median;
            if (n > 1)
            {
                variance = (sumsq - sum * sum / n) / (n - 1);
            }
        }

        return (mean, variance);
    }

    /// <summary>
    /// Ckmeans.1d.dp (Wang &amp; Song 2011) optimal 1-D k-means on sorted <paramref name="x"/> with
    /// <paramref name="k"/> ≤ number of distinct values: returns the first sorted index of each of the k clusters.
    /// Port of <c>EWL2::fill_dp_matrix</c> with R's default <c>method = "linear"</c> (SMAWK row fill) + <c>backtrack</c>.
    /// </summary>
    private static int[] CkmeansClusterStarts(double[] x, int k)
    {
        int n = x.Length;
        if (k == 1)
        {
            return new int[1];
        }

        if (k == n)
        {
            // Every value distinct and its own cluster: the unique zero-cost optimum.
            return Enumerable.Range(0, n).ToArray();
        }

        return CkmeansBacktrackStarts(CkmeansFillDpMatrix(x, k), k);
    }

    /// <summary>
    /// Ckmeans.1d.dp <c>EWL2::fill_dp_matrix(x, w, S, J, "linear")</c>: fills the K × N backtrack matrix J of
    /// optimal 1-D k-means on sorted <paramref name="x"/> for K = <paramref name="kMax"/>, with median-shifted
    /// prefix sums and the O(N) SMAWK row fill <c>fill_row_q_SMAWK</c> (R's default <c>method = "linear"</c>).
    /// Only two rows of S are kept; J keeps every row (needed by <c>backtrack</c> and <c>select_levels</c>).
    /// </summary>
    private static int[][] CkmeansFillDpMatrix(double[] x, int kMax)
    {
        int n = x.Length;

        // Median-shifted running sums for numerical stability (Ckmeans.1d.dp: shift = x[N/2]).
        double shift = x[n / 2];
        double[] sumX = new double[n];
        double[] sumXSq = new double[n];
        sumX[0] = x[0] - shift;
        sumXSq[0] = (x[0] - shift) * (x[0] - shift);

        double[] sPrev = new double[n];
        double[] sCur = new double[n];
        int[][] j = new int[kMax][];
        j[0] = new int[n];
        for (int i = 1; i < n; i++)
        {
            sumX[i] = sumX[i - 1] + x[i] - shift;
            sumXSq[i] = sumXSq[i - 1] + (x[i] - shift) * (x[i] - shift);
            sPrev[i] = CkmeansSsq(0, i, sumX, sumXSq);
        }

        for (int q = 1; q < kMax; q++)
        {
            j[q] = new int[n];
            int imin = q < kMax - 1 ? Math.Max(1, q) : n - 1;
            var row = new CkmeansSmawkRow(q, sPrev, sCur, j[q - 1], j[q], sumX, sumXSq);
            row.Fill(imin, n - 1);
            (sPrev, sCur) = (sCur, sPrev);
        }

        return j;
    }

    /// <summary>
    /// Ckmeans.1d.dp <c>backtrack</c>: walks J from the last value back through <paramref name="k"/> clusters and
    /// returns the first sorted index of each cluster.
    /// </summary>
    private static int[] CkmeansBacktrackStarts(int[][] j, int k)
    {
        int[] starts = new int[k];
        int right = j[0].Length - 1;
        for (int q = k - 1; q >= 0; q--)
        {
            int left = j[q][right];
            starts[q] = left;
            if (q > 0)
            {
                right = left - 1;
            }
        }

        return starts;
    }

    /// <summary>
    /// One row q of the Ckmeans.1d.dp dynamic program filled by <c>EWL2_fill_SMAWK.cpp</c> (Ckmeans.1d.dp 4.3.6):
    /// <c>fill_row_q_SMAWK</c>, <c>SMAWK</c>, <c>reduce_in_place</c>, <c>fill_even_positions</c> and
    /// <c>find_min_from_candidates</c>, ported line by line including their tie-breaking (strict <c>&lt;</c> in the
    /// column reduction, <c>&lt;=</c> — i.e. the largest optimal j — in the candidate scans), which decides the
    /// split chosen among equal-WCSS optima. S[q−1] / S[q] are <c>sPrev</c> / <c>sCur</c>.
    /// </summary>
    private readonly struct CkmeansSmawkRow(
        int q, double[] sPrev, double[] sCur, int[] jPrev, int[] jCur, double[] sumX, double[] sumXSq)
    {
        /// <summary><c>fill_row_q_SMAWK(imin, imax, q, …)</c>: candidate columns js = q..imax.</summary>
        public void Fill(int imin, int imax)
        {
            int[] js = new int[imax - q + 1];
            for (int r = 0; r < js.Length; r++)
            {
                js[r] = q + r;
            }

            Smawk(imin, imax, 1, js);
        }

        private void Smawk(int imin, int imax, int istep, int[] js)
        {
            if (imax - imin <= 0)
            {
                // Base case: only one element left.
                FindMinFromCandidates(imin, imax, istep, js);
                return;
            }

            int[] jsOdd = ReduceInPlace(imin, imax, istep, js);

            int istepx2 = istep << 1;
            int iminOdd = imin + istep;
            int imaxOdd = iminOdd + (imax - iminOdd) / istepx2 * istepx2;

            Smawk(iminOdd, imaxOdd, istepx2, jsOdd);

            FillEvenPositions(imin, imax, istep, js);
        }

        private int[] ReduceInPlace(int imin, int imax, int istep, int[] js)
        {
            int n = (imax - imin) / istep + 1;
            int[] jsRed = (int[])js.Clone();
            if (n >= js.Length)
            {
                return jsRed;
            }

            int left = -1;  // last favourable position / column
            int right = 0;  // current position / column
            int m = jsRed.Length;

            while (m > n)
            {
                int p = left + 1;
                int i = imin + p * istep;
                int j = jsRed[right];
                double sl = sPrev[j - 1] + CkmeansSsq(j, i, sumX, sumXSq);
                int jplus1 = jsRed[right + 1];
                double slplus1 = sPrev[jplus1 - 1] + CkmeansSsq(jplus1, i, sumX, sumXSq);

                if (sl < slplus1 && p < n - 1)
                {
                    jsRed[++left] = j;
                    right++;
                }
                else if (sl < slplus1 && p == n - 1)
                {
                    jsRed[++right] = j;  // delete column p+1
                    m--;
                }
                else
                {
                    if (p > 0)
                    {
                        jsRed[right] = jsRed[left--];
                    }
                    else
                    {
                        right++;  // delete column 0
                    }

                    m--;
                }
            }

            for (int r = left + 1; r < m; ++r)
            {
                jsRed[r] = jsRed[right++];
            }

            return jsRed[..m];
        }

        private void FillEvenPositions(int imin, int imax, int istep, int[] js)
        {
            int n = js.Length;
            int istepx2 = istep << 1;
            int jl = js[0];
            for (int i = imin, r = 0; i <= imax; i += istepx2)
            {
                while (js[r] < jl)
                {
                    r++;
                }

                sCur[i] = sPrev[js[r] - 1] + CkmeansSsq(js[r], i, sumX, sumXSq);
                jCur[i] = js[r];

                int jh = i + istep <= imax ? jCur[i + istep] : js[n - 1];
                int jmax = Math.Min(jh, i);
                double sjimin = CkmeansSsq(jmax, i, sumX, sumXSq);

                for (++r; r < n && js[r] <= jmax; r++)
                {
                    int jabs = js[r];
                    if (jabs > i)
                    {
                        break;
                    }

                    if (jabs < jPrev[i])
                    {
                        continue;
                    }

                    double s = CkmeansSsq(jabs, i, sumX, sumXSq);
                    double sj = sPrev[jabs - 1] + s;
                    if (sj <= sCur[i])
                    {
                        sCur[i] = sj;
                        jCur[i] = js[r];
                    }
                    else if (sPrev[jabs - 1] + sjimin > sCur[i])
                    {
                        break;
                    }
                }

                r--;
                jl = jh;
            }
        }

        private void FindMinFromCandidates(int imin, int imax, int istep, int[] js)
        {
            int rminPrev = 0;
            for (int i = imin; i <= imax; i += istep)
            {
                int rmin = rminPrev;
                sCur[i] = sPrev[js[rmin] - 1] + CkmeansSsq(js[rmin], i, sumX, sumXSq);
                jCur[i] = js[rmin];

                for (int r = rmin + 1; r < js.Length; ++r)
                {
                    int jAbs = js[r];
                    if (jAbs < jPrev[i])
                    {
                        continue;
                    }

                    if (jAbs > i)
                    {
                        break;
                    }

                    double sj = sPrev[jAbs - 1] + CkmeansSsq(jAbs, i, sumX, sumXSq);
                    if (sj <= sCur[i])
                    {
                        sCur[i] = sj;
                        jCur[i] = js[r];
                        rminPrev = r;
                    }
                }
            }
        }
    }

    /// <summary>Ckmeans.1d.dp <c>EWL2::ssq</c>: within-cluster sum of squares of sorted x[j..i] from prefix sums.</summary>
    private static double CkmeansSsq(int j, int i, double[] sumX, double[] sumXSq)
    {
        double sji;
        if (j >= i)
        {
            sji = 0.0;
        }
        else if (j > 0)
        {
            double muji = (sumX[i] - sumX[j - 1]) / (i - j + 1);
            sji = sumXSq[i] - sumXSq[j - 1] - (i - j + 1) * muji * muji;
        }
        else
        {
            sji = sumXSq[i] - sumX[i] * sumX[i] / (i + 1);
        }

        return sji < 0 ? 0 : sji;
    }
}
