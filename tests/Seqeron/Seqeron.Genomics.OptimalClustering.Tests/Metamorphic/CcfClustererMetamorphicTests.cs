// SPDX-License-Identifier: LGPL-3.0-or-later
// Tests of the port of Ckmeans.1d.dp 4.3.6 (CcfClusterer), part of Seqeron.Genomics.OptimalClustering.Tests; licensed under the GNU Lesser General Public License, version 3 or later (see LICENSE in the production project).
using FluentAssertions;
using Seqeron.Genomics.Oncology;

namespace Seqeron.Genomics.Tests.Metamorphic;

/// <summary>ONCO-CCF-001 (Ckmeans.1d.dp): the optimum is computed on the sorted values, so permuting the input gives the same centroids and the same (value, cluster) pairs.</summary>
[TestFixture]
[Category("Metamorphic")]
[Category("Oncology")]
public class CcfClustererMetamorphicTests
{
    [Test]
    [Description("INV (Ckmeans.1d.dp works on sorted data): permuting the CCFs gives identical centroids and identical (value, cluster) pairs.")]
    public void ClusterCcfValues_InputPermutation_IsInvariant()
    {
        for (int seed = 0; seed < 300; seed++)
        {
            var rng = new Random(seed);
            int n = rng.Next(1, 40);
            double[] values = Enumerable.Range(0, n).Select(_ => rng.Next(0, 60) / 50.0 * 0.8).ToArray();
            int k = rng.Next(1, n + 1);
            double[] permuted = values.OrderBy(_ => rng.Next()).ToArray();

            var a = CcfClusterer.ClusterCcfValues(values, k);
            var b = CcfClusterer.ClusterCcfValues(permuted, k);

            b.Centroids.Should().Equal(a.Centroids, "seed {0}", seed);
            b.ClonalClusterIndex.Should().Be(a.ClonalClusterIndex, "seed {0}", seed);
            var pairsA = values.Select((v, i) => (v, a.Assignments[i])).OrderBy(p => p.v).ThenBy(p => p.Item2);
            var pairsB = permuted.Select((v, i) => (v, b.Assignments[i])).OrderBy(p => p.v).ThenBy(p => p.Item2);
            pairsB.Should().Equal(pairsA, "each value keeps its cluster (seed {0})", seed);
        }
    }
}
