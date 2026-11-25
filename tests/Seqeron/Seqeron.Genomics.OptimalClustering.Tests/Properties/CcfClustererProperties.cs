// SPDX-License-Identifier: LGPL-3.0-or-later
// Tests of the port of Ckmeans.1d.dp 4.3.6 (CcfClusterer), part of Seqeron.Genomics.OptimalClustering.Tests; licensed under the GNU Lesser General Public License, version 3 or later (see LICENSE in the production project).
using FsCheck;
using FsCheck.Fluent;
using Seqeron.Genomics.Oncology;
using static Seqeron.Genomics.Oncology.OncologyAnalyzer;

namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Property-based tests of the Ckmeans.1d.dp port: the returned partition attains the optimal within-cluster sum of squares
/// (independent dynamic program / brute force over the sorted values), centroids ascend and are the cluster means, the
/// BIC-selected k lies in the adjusted bounds and equals the fixed-k clustering for that k; constant input gives one cluster.
/// </summary>
[TestFixture]
[Category("Property")]
[Category("Oncology")]
public class CcfClustererProperties
{
    private const double CcfTolerance = 1e-9;

    private static Arbitrary<(double[] Values, int K, int Seed)> CcfArbitrary() =>
        (from n in Gen.Choose(1, 40)
         from raw in Gen.Choose(0, 1000).ArrayOf(n)
         from k in Gen.Choose(1, Math.Max(1, Math.Min(6, n)))
         from seed in Gen.Choose(0, int.MaxValue)
         select (raw.Select(v => v / 1000.0).ToArray(), k, seed)).ToArbitrary();

    private static double OptimalWcss(double[] values, int k)
    {
        double[] x = values.OrderBy(v => v).ToArray();
        int n = x.Length;
        double Cost(int i, int j) // inclusive [i, j]
        {
            double mean = 0;
            for (int t = i; t <= j; t++) mean += x[t];
            mean /= j - i + 1;
            double s = 0;
            for (int t = i; t <= j; t++) s += (x[t] - mean) * (x[t] - mean);
            return s;
        }

        var d = new double[k + 1, n + 1];
        for (int j = 1; j <= n; j++) d[0, j] = double.PositiveInfinity;
        for (int c = 1; c <= k; c++)
        {
            for (int j = 0; j <= n; j++)
            {
                if (j < c) { d[c, j] = j == 0 ? 0 : double.PositiveInfinity; continue; }
                double best = double.PositiveInfinity;
                for (int i = c; i <= j; i++) // last cluster = x[i-1 .. j-1]
                {
                    best = Math.Min(best, d[c - 1, i - 1] + Cost(i - 1, j - 1));
                }

                d[c, j] = best;
            }
        }

        return d[k, n];
    }

    private static double Wcss(double[] values, CcfClustering c) =>
        values.Select((v, i) => (v - c.Centroids[c.Assignments[i]]) * (v - c.Centroids[c.Assignments[i]])).Sum();

    /// <summary>
    /// F46 (SMAWK row fill): the returned clustering attains the brute-force optimal WCSS for k* = min(k, #distinct)
    /// clusters (only the choice among equal-WCSS optima may differ from the quadratic fill), centroids are ascending
    /// cluster means and the clonal index is the last centroid.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property ClusterCcfValues_AttainsBruteForceOptimalWcss()
    {
        return Prop.ForAll(CcfArbitrary(), t =>
        {
            var c = CcfClusterer.ClusterCcfValues(t.Values, t.K);
            int kEff = Math.Min(t.K, t.Values.Distinct().Count());
            double wcss = Wcss(t.Values, c);
            double optimum = OptimalWcss(t.Values, kEff);
            bool means = Enumerable.Range(0, c.Centroids.Count).All(j =>
                Math.Abs(c.Centroids[j] - t.Values.Where((_, i) => c.Assignments[i] == j).Average()) <= 1e-12);
            bool ascending = c.Centroids.Zip(c.Centroids.Skip(1), (a, b) => a < b).All(ok => ok);
            return (c.Centroids.Count == kEff && Math.Abs(wcss - optimum) <= 1e-9 * (1 + optimum) && means && ascending
                    && c.ClonalClusterIndex == c.Centroids.Count - 1)
                .Label($"k*={kEff}, centroids={c.Centroids.Count}, wcss={wcss}, optimum={optimum}");
        });
    }

    /// <summary>
    /// F47 (Ckmeans.1d.dp <c>select_levels</c> BIC): the selected k lies in [min(k_min, d), min(k_max, d)] (d = #distinct,
    /// R's bound adjustment) and the clustering is the fixed-k clustering for that k; all-equal input gives one cluster.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property ClusterCcfValues_BicSelection_WithinBounds_EqualsFixedK()
    {
        return Prop.ForAll(CcfArbitrary(), t =>
        {
            int kMin = Math.Max(1, t.K - 1), kMax = t.K + 2;
            int distinct = t.Values.Distinct().Count();
            var auto = CcfClusterer.ClusterCcfValues(t.Values, kMin, kMax);
            int k = auto.Centroids.Count;
            int lo = distinct < kMin ? distinct : kMin, hi = Math.Min(kMax, distinct);
            var fixedK = CcfClusterer.ClusterCcfValues(t.Values, k);
            var constant = CcfClusterer.ClusterCcfValues(Enumerable.Repeat(t.Values[0], t.Values.Length).ToArray(), kMin, kMax);
            return (k >= lo && k <= hi && auto.Assignments.SequenceEqual(fixedK.Assignments)
                    && auto.Centroids.Zip(fixedK.Centroids, (a, b) => Math.Abs(a - b) <= 1e-12).All(ok => ok)
                    && constant.Centroids.Count == 1)
                .Label($"k={k} ∉ [{lo}, {hi}] or differs from fixed-k");
        });
    }

    /// <summary>
    /// <c>ClusterCcfValues</c> produces min(k, distinct) ascending centroids, every cluster non-empty, one valid
    /// assignment per input value, and reports the highest-centroid cluster (last index) as clonal.
    /// (Tarabichi 2021 highest-CP-clonal; Ckmeans.1d.dp Kmax = min(k, unique))
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property ClusterCcfValues_AscendingCentroids_ValidAssignments_HighestIsClonal()
    {
        var arb = (from n in Gen.Choose(1, 12)
                   from values in Gen.Choose(0, 1000).Select(v => v / 1000.0).ArrayOf(n)
                   from k in Gen.Choose(1, n)
                   select (values, k)).ToArbitrary();

        return Prop.ForAll(arb, t =>
        {
            var clustering = CcfClusterer.ClusterCcfValues(t.values, t.k);

            // Ckmeans.1d.dp (Wang & Song 2011): Kmax = min(k, number of distinct values); no empty clusters.
            int expectedCount = Math.Min(t.k, t.values.Distinct().Count());
            bool centroidCount = clustering.Centroids.Count == expectedCount;
            bool ascending = true;
            for (int i = 1; i < clustering.Centroids.Count; i++)
            {
                ascending &= clustering.Centroids[i] >= clustering.Centroids[i - 1] - CcfTolerance;
            }

            bool assignmentsOk = clustering.Assignments.Count == t.values.Length
                && clustering.Assignments.All(a => a >= 0 && a < expectedCount)
                && clustering.Assignments.Distinct().Count() == expectedCount;
            bool clonalIsHighest = clustering.ClonalClusterIndex == expectedCount - 1;

            return (centroidCount && ascending && assignmentsOk && clonalIsHighest)
                .Label($"k={t.k}, centroids=[{string.Join(",", clustering.Centroids)}], clonalIndex={clustering.ClonalClusterIndex}");
        });
    }

    /// <summary>
    /// O (optimality, Wang &amp; Song 2011): <c>ClusterCcfValues</c> attains the minimum within-cluster sum of
    /// squares over <i>all</i> partitions of the sorted values into min(k, distinct) contiguous non-empty blocks
    /// (an optimal 1-D k-means partition is contiguous). The oracle enumerates every partition by brute force,
    /// independently of production.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property ClusterCcfValues_AttainsBruteForceMinimumWcss()
    {
        var arb = (from n in Gen.Choose(1, 9)
                   from values in Gen.Choose(0, 1000).Select(v => v / 1000.0).ArrayOf(n)
                   from k in Gen.Choose(1, n)
                   select (values, k)).ToArbitrary();

        return Prop.ForAll(arb, t =>
        {
            var clustering = CcfClusterer.ClusterCcfValues(t.values, t.k);
            double wcss = t.values.Select((v, i) => Math.Pow(v - clustering.Centroids[clustering.Assignments[i]], 2)).Sum();

            double[] sorted = t.values.OrderBy(v => v).ToArray();
            int kEff = Math.Min(t.k, sorted.Distinct().Count());
            double best = BruteForceMinWcss(sorted, 0, kEff);
            return (wcss <= best + 1e-12).Label($"WCSS {wcss} > brute-force optimum {best} (k={kEff})");
        });
    }

    private static double BruteForceMinWcss(double[] sorted, int start, int blocks)
    {
        if (blocks == 1)
        {
            return BlockSsq(sorted, start, sorted.Length - 1);
        }

        double best = double.PositiveInfinity;
        for (int end = start; end <= sorted.Length - blocks; end++)
        {
            best = Math.Min(best, BlockSsq(sorted, start, end) + BruteForceMinWcss(sorted, end + 1, blocks - 1));
        }

        return best;
    }

    private static double BlockSsq(double[] sorted, int from, int to)
    {
        double mean = 0.0;
        for (int i = from; i <= to; i++)
        {
            mean += sorted[i];
        }

        mean /= to - from + 1;
        double ssq = 0.0;
        for (int i = from; i <= to; i++)
        {
            ssq += (sorted[i] - mean) * (sorted[i] - mean);
        }

        return ssq;
    }


    /// <summary>
    /// Agreement of the two implementations: the MIT <c>AnalyzeHeterogeneity</c> (independent optimal one-dimensional k-means)
    /// reports the same subclone count and clone-size Shannon diversity as the Ckmeans.1d.dp port clustering of the same values.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property AnalyzeHeterogeneity_AgreesWithCkmeansPort()
    {
        var arb = (from n in Gen.Choose(1, 30)
                   from raw in Gen.Choose(0, 200).ArrayOf(n)
                   from k in Gen.Choose(1, n)
                   select (raw.Select(v => v / 200.0).ToArray(), k)).ToArbitrary();

        return Prop.ForAll(arb, t =>
        {
            var port = CcfClusterer.ClusterCcfValues(t.Item1, t.k);
            var sizes = new int[port.Centroids.Count];
            foreach (int label in port.Assignments)
                sizes[label]++;
            var mit = OncologyAnalyzer.AnalyzeHeterogeneity(t.Item1, t.Item1, t.k);
            return (mit.SubcloneCount == InferSubclones(port) && Math.Abs(mit.ShannonDiversity - StatisticsHelper.ShannonIndex(sizes)) <= 1e-12)
                .Label($"k={t.k}, MIT subclones {mit.SubcloneCount}, shannon {mit.ShannonDiversity}");
        });
    }
}
