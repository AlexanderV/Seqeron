namespace Seqeron.Genomics.Oncology;

public static partial class OncologyAnalyzer
{
    /// <summary>
    /// Optimal one-dimensional k-means of cancer cell fractions into exactly <paramref name="clusterCount"/> clusters (or as many as there are distinct values): the partition of the sorted values into contiguous
    /// blocks that minimises the within-cluster sum of squares, found by the dynamic program of Bellman (1973) / Fisher (1958) with the divide-and-conquer row fill (the optimal split point is monotone in the
    /// end index, Gries / Knuth). An independent implementation of the published principle: <c>AnalyzeHeterogeneity</c> uses it; the port of the R package Ckmeans.1d.dp (BIC choice of the number of
    /// clusters, the exact tie-breaking of R) is a separate assembly.
    /// </summary>
    private static CcfClustering ClusterCcfValuesOptimal(IReadOnlyList<double> ccfValues, int clusterCount)
    {
        ArgumentNullException.ThrowIfNull(ccfValues);
        int n = ccfValues.Count;
        if (n == 0)
            throw new ArgumentException("At least one CCF value is required.", nameof(ccfValues));
        for (int i = 0; i < n; i++)
        {
            if (double.IsNaN(ccfValues[i]) || double.IsInfinity(ccfValues[i]))
                throw new ArgumentException($"CCF value must be finite; got {ccfValues[i]} at index {i}.", nameof(ccfValues));
        }

        if (clusterCount < 1 || clusterCount > n)
            throw new ArgumentOutOfRangeException(nameof(clusterCount), clusterCount, $"Cluster count must be in [1, {n}].");

        int[] order = Enumerable.Range(0, n).OrderBy(i => ccfValues[i]).ToArray();
        double[] x = new double[n];
        for (int s = 0; s < n; s++)
            x[s] = ccfValues[order[s]];
        int distinct = 1;
        for (int s = 1; s < n; s++)
        {
            if (x[s] != x[s - 1])
                distinct++;
        }

        int k = Math.Min(clusterCount, distinct);
        int[] starts = k == 1 ? new int[1] : OptimalBlockStarts(x, k);

        var centroids = new double[k];
        var blockOfSorted = new int[n];
        for (int c = 0; c < k; c++)
        {
            int from = starts[c], to = c + 1 < k ? starts[c + 1] : n;
            double sum = 0;
            for (int s = from; s < to; s++)
            {
                sum += x[s];
                blockOfSorted[s] = c;
            }

            centroids[c] = sum / (to - from);
        }

        var assignments = new int[n];
        for (int s = 0; s < n; s++)
            assignments[order[s]] = blockOfSorted[s];
        return new CcfClustering(centroids, assignments, k - 1);
    }

    // first sorted index of each of the k blocks of the optimal partition of sorted x (every block non-empty)
    private static int[] OptimalBlockStarts(double[] x, int k)
    {
        int n = x.Length;
        double shift = x[n / 2];
        var prefix = new double[n + 1];
        var prefixSq = new double[n + 1];
        for (int i = 0; i < n; i++)
        {
            double d = x[i] - shift;
            prefix[i + 1] = prefix[i] + d;
            prefixSq[i + 1] = prefixSq[i] + (d * d);
        }

        // sum of squares about the mean of x[a..b] (inclusive)
        double Cost(int a, int b)
        {
            double sum = prefix[b + 1] - prefix[a];
            double sq = prefixSq[b + 1] - prefixSq[a];
            return Math.Max(0, sq - (sum * sum / (b - a + 1)));
        }

        var previous = new double[n];
        var current = new double[n];
        var split = new int[k][];
        for (int i = 0; i < n; i++)
            previous[i] = Cost(0, i);
        split[0] = new int[n];
        for (int row = 1; row < k; row++)
        {
            split[row] = new int[n];
            Array.Fill(current, double.PositiveInfinity);
            Fill(row, row, n - 1, row - 1, n - 1);
            (previous, current) = (current, previous);
        }

        var starts = new int[k];
        int end = n - 1;
        for (int row = k - 1; row >= 1; row--)
        {
            starts[row] = split[row][end];
            end = starts[row] - 1;
        }

        starts[0] = 0;
        return starts;

        // current[i] = min over j of previous[j - 1] + Cost(j, i), the last block being x[j..i]; the best j is non-decreasing in i
        void Fill(int row, int lo, int hi, int jLo, int jHi)
        {
            if (lo > hi)
                return;
            int mid = (lo + hi) / 2;
            double best = double.PositiveInfinity;
            int bestJ = jLo;
            for (int j = Math.Max(jLo, row); j <= Math.Min(mid, jHi); j++)
            {
                double c = previous[j - 1] + Cost(j, mid);
                if (double.IsPositiveInfinity(best) || c < best - (1e-12 * Math.Abs(best)))
                {
                    best = c;
                    bestJ = j;
                }
            }

            current[mid] = best;
            split[row][mid] = bestJ;
            Fill(row, lo, mid - 1, jLo, bestJ);
            Fill(row, mid + 1, hi, bestJ, jHi);
        }
    }
}
