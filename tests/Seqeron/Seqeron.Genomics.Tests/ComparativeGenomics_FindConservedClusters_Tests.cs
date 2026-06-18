// COMPGEN-CLUSTER-001 — Conserved Gene Clusters (common intervals of permutations)
// Evidence: docs/Evidence/COMPGEN-CLUSTER-001-Evidence.md
// TestSpec: tests/TestSpecs/COMPGEN-CLUSTER-001.md
// Source: Bui-Xuan B-M, Habib M, Paul C (2013). MinMax-Profiles. arXiv:1304.5140 (Def. 1, Example 1).
//         Uno T, Yagiura M (2000). Algorithmica 26(2):290-309 (common-interval model).
//         Didier G, Schmidt T, Stoye J, Tsur D (2013). arXiv:1310.4290 (sequences with duplicates).
//         Heber S, Stoye J (2001). CPM, LNCS 2089:207-218 (k-permutation common intervals).

using NUnit.Framework;
using Seqeron.Genomics.Analysis;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Seqeron.Genomics.Tests;

[TestFixture]
public class ComparativeGenomics_FindConservedClusters_Tests
{
    #region Helpers

    // Builds a genome from an ordered list of ortholog-group labels. Each label becomes one gene
    // whose id is "<genomeId>_<index>"; the ortholog map sends that gene id to the group label.
    private static (IReadOnlyList<ComparativeGenomics.Gene> genome, Dictionary<string, string> map)
        Genome(string genomeId, params string[] groupLabels)
    {
        var genes = new List<ComparativeGenomics.Gene>(groupLabels.Length);
        var map = new Dictionary<string, string>();
        for (int i = 0; i < groupLabels.Length; i++)
        {
            string geneId = $"{genomeId}_{i}";
            genes.Add(new ComparativeGenomics.Gene(geneId, genomeId, i * 100, i * 100 + 50, '+'));
            // A null/empty label means "no ortholog group" (gene absent from the map).
            if (!string.IsNullOrEmpty(groupLabels[i]))
                map[geneId] = groupLabels[i];
        }
        return (genes, map);
    }

    // Renders the returned clusters as a set of canonical "sorted-joined" keys for exact comparison.
    private static HashSet<string> KeySet(IEnumerable<IReadOnlyList<string>> clusters)
        => clusters.Select(c => string.Join(",", c.OrderBy(x => x, StringComparer.Ordinal)))
                   .ToHashSet(StringComparer.Ordinal);

    private static string Key(params string[] labels)
        => string.Join(",", labels.OrderBy(x => x, StringComparer.Ordinal));

    #endregion

    #region FindConservedClusters

    // M1 — Golden vector from MinMax-Profiles Example 1: P1 = Id7 = (1..7), P2 = (7 2 1 3 6 4 5).
    // The paper lists the common intervals as {1,2},{1,2,3},{3,4,5,6},{4,5},{4,5,6},{1..6},{1..7}.
    // (Independently recomputed by brute force.) With minClusterSize=2 the result must be EXACTLY these.
    [Test]
    public void FindConservedClusters_GoldenVectorExample1_ReturnsExactCommonIntervals()
    {
        var (g1, m1) = Genome("A", "1", "2", "3", "4", "5", "6", "7");
        var (g2, m2) = Genome("B", "7", "2", "1", "3", "6", "4", "5");
        var map = m1.Concat(m2).ToDictionary(kv => kv.Key, kv => kv.Value);
        var genomes = new List<IReadOnlyList<ComparativeGenomics.Gene>> { g1, g2 };

        var clusters = ComparativeGenomics.FindConservedClusters(genomes, map, minClusterSize: 2).ToList();
        var actual = KeySet(clusters);

        var expected = new HashSet<string>(StringComparer.Ordinal)
        {
            Key("1", "2"), Key("4", "5"), Key("1", "2", "3"), Key("4", "5", "6"),
            Key("3", "4", "5", "6"), Key("1", "2", "3", "4", "5", "6"),
            Key("1", "2", "3", "4", "5", "6", "7"),
        };

        Assert.That(actual, Is.EquivalentTo(expected),
            "Common intervals must equal exactly the seven sets of MinMax-Profiles Example 1 (Bui-Xuan et al. 2013).");
    }

    // M2 — A set contiguous in only ONE genome is not a common interval. In the golden vector {2,3}
    // is contiguous in P1 but 2 (pos 2) and 3 (pos 4) are non-adjacent in P2, so it must be excluded.
    [Test]
    public void FindConservedClusters_SetSplitInSecondGenome_NotReturned()
    {
        var (g1, m1) = Genome("A", "1", "2", "3", "4", "5", "6", "7");
        var (g2, m2) = Genome("B", "7", "2", "1", "3", "6", "4", "5");
        var map = m1.Concat(m2).ToDictionary(kv => kv.Key, kv => kv.Value);
        var genomes = new List<IReadOnlyList<ComparativeGenomics.Gene>> { g1, g2 };

        var actual = KeySet(ComparativeGenomics.FindConservedClusters(genomes, map, minClusterSize: 2));

        Assert.That(actual, Does.Not.Contain(Key("2", "3")),
            "{2,3} is contiguous in P1 but split in P2 (positions 2 and 4), so it is not a common interval.");
    }

    // M3 — Sequences with a repeated ortholog-group label (paralog). A set is a common interval iff
    // SOME contiguous window in each genome has exactly that set (Didier et al. 2013). Here genome B
    // duplicates group "x"; the set {a,b,c} still occupies a contiguous window in both genomes.
    [Test]
    public void FindConservedClusters_RepeatedLabels_WindowStillMatches()
    {
        var (g1, m1) = Genome("A", "a", "b", "c", "x");
        var (g2, m2) = Genome("B", "x", "c", "b", "a", "x"); // "x" duplicated; a,b,c contiguous (reversed)
        var map = m1.Concat(m2).ToDictionary(kv => kv.Key, kv => kv.Value);
        var genomes = new List<IReadOnlyList<ComparativeGenomics.Gene>> { g1, g2 };

        var actual = KeySet(ComparativeGenomics.FindConservedClusters(genomes, map, minClusterSize: 3));

        // Exact common-interval set (size >= 3), independently brute-forced over the sequence model
        // (Didier et al. 2013, Set(T[i..j])): {a,b,c}, {a,b,c,x} and {b,c,x} each occupy a contiguous
        // window in both genomes despite the duplicated 'x'.
        var expected = new HashSet<string>(StringComparer.Ordinal)
        {
            Key("a", "b", "c"), Key("a", "b", "c", "x"), Key("b", "c", "x"),
        };
        Assert.That(actual, Is.EquivalentTo(expected),
            "The repeated label 'x' does not stop {a,b,c}, {a,b,c,x}, {b,c,x} from each being a common interval (some window in each genome has exactly that label set).");
    }

    // M4 — minClusterSize filters out smaller common intervals. With size>=4 on the golden vector,
    // only {3,4,5,6}, {1..6}, {1..7} qualify.
    [Test]
    public void FindConservedClusters_MinClusterSizeFour_ReturnsOnlyLargeIntervals()
    {
        var (g1, m1) = Genome("A", "1", "2", "3", "4", "5", "6", "7");
        var (g2, m2) = Genome("B", "7", "2", "1", "3", "6", "4", "5");
        var map = m1.Concat(m2).ToDictionary(kv => kv.Key, kv => kv.Value);
        var genomes = new List<IReadOnlyList<ComparativeGenomics.Gene>> { g1, g2 };

        var actual = KeySet(ComparativeGenomics.FindConservedClusters(genomes, map, minClusterSize: 4));

        var expected = new HashSet<string>(StringComparer.Ordinal)
        {
            Key("3", "4", "5", "6"), Key("1", "2", "3", "4", "5", "6"),
            Key("1", "2", "3", "4", "5", "6", "7"),
        };
        Assert.That(actual, Is.EquivalentTo(expected),
            "Only common intervals of size >= 4 survive the minClusterSize filter.");
    }

    // M5 — An interval is the set of ALL elements of a window: a foreign group between members breaks
    // the window. Genome B places "z" between a,b,c, so {a,b,c} is NOT a common interval there.
    [Test]
    public void FindConservedClusters_ForeignGroupInsideWindow_BreaksCluster()
    {
        var (g1, m1) = Genome("A", "a", "b", "c");
        var (g2, m2) = Genome("B", "a", "z", "b", "c"); // z separates a from b,c
        var map = m1.Concat(m2).ToDictionary(kv => kv.Key, kv => kv.Value);
        var genomes = new List<IReadOnlyList<ComparativeGenomics.Gene>> { g1, g2 };

        var actual = KeySet(ComparativeGenomics.FindConservedClusters(genomes, map, minClusterSize: 3));

        // An interval is the set of ALL window elements; 'z' inside the only {a,b,c}-spanning window in
        // genome B means no window has set exactly {a,b,c}. No size>=3 set is common here, so the result
        // is EXACTLY empty (independently brute-forced over the sequence model). Lock the full set, not
        // just the absence of {a,b,c}, so a wrong implementation that over-reports cannot pass.
        Assert.That(actual, Is.Empty,
            "The foreign group 'z' inside the window in genome B leaves no common interval of size >= 3.");
    }

    // S1 — A common interval is a family (K>=2) notion; a single genome yields no conserved clusters.
    [Test]
    public void FindConservedClusters_SingleGenome_ReturnsEmpty()
    {
        var (g1, m1) = Genome("A", "1", "2", "3", "4");
        var genomes = new List<IReadOnlyList<ComparativeGenomics.Gene>> { g1 };

        var clusters = ComparativeGenomics.FindConservedClusters(genomes, m1).ToList();

        Assert.That(clusters, Is.Empty,
            "Common intervals require >= 2 genomes (family definition); one genome gives no clusters.");
    }

    // S2 — Identical gene order across THREE genomes: every window of size >= minClusterSize is a
    // common interval (identity vs identity, all intervals common). For order 1..5, minClusterSize=3
    // the conserved sets are exactly the contiguous windows of size 3,4,5.
    [Test]
    public void FindConservedClusters_IdenticalOrderThreeGenomes_AllWindowsConserved()
    {
        var (g1, m1) = Genome("A", "1", "2", "3", "4", "5");
        var (g2, m2) = Genome("B", "1", "2", "3", "4", "5");
        var (g3, m3) = Genome("C", "1", "2", "3", "4", "5");
        var map = m1.Concat(m2).Concat(m3).ToDictionary(kv => kv.Key, kv => kv.Value);
        var genomes = new List<IReadOnlyList<ComparativeGenomics.Gene>> { g1, g2, g3 };

        var actual = KeySet(ComparativeGenomics.FindConservedClusters(genomes, map, minClusterSize: 3));

        var expected = new HashSet<string>(StringComparer.Ordinal)
        {
            Key("1", "2", "3"), Key("2", "3", "4"), Key("3", "4", "5"),
            Key("1", "2", "3", "4"), Key("2", "3", "4", "5"),
            Key("1", "2", "3", "4", "5"),
        };
        Assert.That(actual, Is.EquivalentTo(expected),
            "With identical order in all genomes every contiguous window of size >= 3 is a common interval.");
    }

    // S3 — A set conserved in genomes 1 and 2 but split in genome 3 must NOT be reported: the cluster
    // must be an interval of EVERY genome (Heber & Stoye 2001, k-permutation common intervals).
    [Test]
    public void FindConservedClusters_ConservedInTwoButSplitInThird_NotReturned()
    {
        var (g1, m1) = Genome("A", "1", "2", "3", "4");
        var (g2, m2) = Genome("B", "1", "2", "3", "4");
        var (g3, m3) = Genome("C", "1", "9", "2", "3", "4"); // foreign "9" splits {1,2}
        var map = m1.Concat(m2).Concat(m3).ToDictionary(kv => kv.Key, kv => kv.Value);
        var genomes = new List<IReadOnlyList<ComparativeGenomics.Gene>> { g1, g2, g3 };

        var actual = KeySet(ComparativeGenomics.FindConservedClusters(genomes, map, minClusterSize: 2));

        // Exact common-interval set (size >= 2), independently brute-forced: '9' in genome C splits
        // every set containing both 1 and 2, leaving exactly {2,3}, {3,4}, {2,3,4}. Lock the full set
        // so the test fails against an implementation that wrongly keeps {1,2} or drops {2,3,4}.
        var expected = new HashSet<string>(StringComparer.Ordinal)
        {
            Key("2", "3"), Key("3", "4"), Key("2", "3", "4"),
        };
        Assert.That(actual, Is.EquivalentTo(expected),
            "Only sets that stay contiguous in all three genomes are common: {2,3}, {3,4}, {2,3,4}; {1,2} is split by '9' in genome C.");
    }

    // C1 — Determinism: identical inputs produce an identical cluster sequence on repeated runs.
    [Test]
    public void FindConservedClusters_RepeatedRuns_AreDeterministic()
    {
        var (g1, m1) = Genome("A", "1", "2", "3", "4", "5", "6", "7");
        var (g2, m2) = Genome("B", "7", "2", "1", "3", "6", "4", "5");
        var map = m1.Concat(m2).ToDictionary(kv => kv.Key, kv => kv.Value);
        var genomes = new List<IReadOnlyList<ComparativeGenomics.Gene>> { g1, g2 };

        var run1 = ComparativeGenomics.FindConservedClusters(genomes, map, minClusterSize: 2)
            .Select(c => string.Join(",", c)).ToList();
        var run2 = ComparativeGenomics.FindConservedClusters(genomes, map, minClusterSize: 2)
            .Select(c => string.Join(",", c)).ToList();

        Assert.That(run2, Is.EqualTo(run1),
            "The result is order-independent and reproducible across runs (deterministic enumeration).");
    }

    // C2 — When the only common interval of size >= minClusterSize is the trivial whole set, no
    // smaller cluster is reported. For A=(1,2,3,4), B=(2,4,1,3) every size-3 subset is split in B,
    // so only the full set {1,2,3,4} qualifies (independently verified by brute force).
    [Test]
    public void FindConservedClusters_OnlyTrivialWholeSetShared_ReturnsOnlyWholeSet()
    {
        var (g1, m1) = Genome("A", "1", "2", "3", "4");
        var (g2, m2) = Genome("B", "2", "4", "1", "3");
        var map = m1.Concat(m2).ToDictionary(kv => kv.Key, kv => kv.Value);
        var genomes = new List<IReadOnlyList<ComparativeGenomics.Gene>> { g1, g2 };

        var actual = KeySet(ComparativeGenomics.FindConservedClusters(genomes, map, minClusterSize: 3));

        Assert.That(actual, Is.EquivalentTo(new[] { Key("1", "2", "3", "4") }),
            "Every size-3 subset is split in B; only the trivial whole-set common interval remains.");
    }

    // C3 — Null arguments throw ArgumentNullException (defensive contract).
    [Test]
    public void FindConservedClusters_NullArguments_Throw()
    {
        var (g1, m1) = Genome("A", "1", "2", "3");
        var genomes = new List<IReadOnlyList<ComparativeGenomics.Gene>> { g1 };

        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentNullException>(
                () => ComparativeGenomics.FindConservedClusters(null!, m1).ToList(),
                "Null genomes must throw ArgumentNullException.");
            Assert.Throws<ArgumentNullException>(
                () => ComparativeGenomics.FindConservedClusters(genomes, null!).ToList(),
                "Null orthologGroups must throw ArgumentNullException.");
        });
    }

    // P1 (property-based, INV-01) — For randomly generated genomes (fixed seed for determinism),
    // EVERY returned cluster must be a contiguous window (interval) in EVERY genome. This is the
    // defining property of a common interval (Bui-Xuan et al. 2013, Def. 1) and must hold for any
    // input, not just the worked example.
    [Test]
    public void FindConservedClusters_RandomGenomes_EveryClusterIsIntervalOfEveryGenome()
    {
        var rng = new Random(20260614); // fixed, documented seed → deterministic.
        const int genomeCount = 3;
        const int geneCount = 12;

        var orders = new List<string[]>();
        var map = new Dictionary<string, string>();
        var genomes = new List<IReadOnlyList<ComparativeGenomics.Gene>>();
        var labels = Enumerable.Range(1, geneCount).Select(i => i.ToString()).ToArray();

        for (int g = 0; g < genomeCount; g++)
        {
            var shuffled = labels.OrderBy(_ => rng.Next()).ToArray();
            orders.Add(shuffled);
            var (genome, m) = Genome($"G{g}", shuffled);
            genomes.Add(genome);
            foreach (var kv in m) map[kv.Key] = kv.Value;
        }

        var clusters = ComparativeGenomics.FindConservedClusters(genomes, map, minClusterSize: 2).ToList();

        // Independent contiguity check: a label set is an interval of an order iff the min and max
        // positions of its members span exactly |set| positions.
        static bool IsContiguous(string[] order, ISet<string> set)
        {
            var positions = new List<int>();
            for (int i = 0; i < order.Length; i++)
                if (set.Contains(order[i])) positions.Add(i);
            return positions.Count == set.Count &&
                   positions.Max() - positions.Min() == set.Count - 1;
        }

        Assert.Multiple(() =>
        {
            foreach (var cluster in clusters)
            {
                var set = cluster.ToHashSet(StringComparer.Ordinal);
                foreach (var order in orders)
                    Assert.That(IsContiguous(order, set), Is.True,
                        $"Returned cluster {{{string.Join(",", cluster)}}} must be contiguous in every genome (INV-01).");
            }
        });
    }

    #endregion
}
