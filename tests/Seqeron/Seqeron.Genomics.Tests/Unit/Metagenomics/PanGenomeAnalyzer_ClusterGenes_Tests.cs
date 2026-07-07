// PANGEN-CLUSTER-001 — Gene Clustering (homolog grouping by global sequence identity)
// Evidence: docs/Evidence/PANGEN-CLUSTER-001-Evidence.md
// TestSpec: tests/TestSpecs/PANGEN-CLUSTER-001.md
// Source: Li W, Godzik A (2006) Bioinformatics 22(13):1658 (CD-HIT); CD-HIT User's Guide / Algorithm wiki;
//         Page AJ et al. (2015) Bioinformatics 31(22):3691 (Roary).

namespace Seqeron.Genomics.Tests.Unit.Metagenomics;

[TestFixture]
public class PanGenomeAnalyzer_ClusterGenes_Tests
{
    private static Dictionary<string, IReadOnlyList<(string GeneId, string Sequence)>> Genomes(
        params (string Genome, (string GeneId, string Sequence)[] Genes)[] entries)
    {
        var dict = new Dictionary<string, IReadOnlyList<(string GeneId, string Sequence)>>();
        foreach (var (genome, genes) in entries)
            dict[genome] = genes.ToList();
        return dict;
    }

    // Evidence worked-value sequences (CD-HIT global identity = identical positions / shorter length).
    private const string S = "ATGCATGC";          // 8 bp
    private const string S1Sub = "ATGCATGG";       // 7/8 = 0.875 vs S
    private const string SLong = "ATGCATGCAAAA";   // 12 bp; vs S -> 8/8 over shorter(8) = 1.0
    private const string Diff = "CGTACGTA";        // 0/8 vs S

    #region ClusterGenes — identity-based grouping (CD-HIT global identity)

    // M1 — Identical sequence in 3 genomes -> one cluster spanning all 3 genomes; AverageIdentity 1.0.
    [Test]
    public void ClusterGenes_IdenticalSequences_SingleClusterAllGenomes()
    {
        var genomes = Genomes(
            ("g1", new[] { ("a", S) }),
            ("g2", new[] { ("b", S) }),
            ("g3", new[] { ("c", S) }));

        var clusters = PanGenomeAnalyzer.ClusterGenes(genomes, identityThreshold: 0.9).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(clusters, Has.Count.EqualTo(1),
                "Three identical sequences (identity 1.0) collapse into one cluster (CD-HIT greedy join).");
            Assert.That(clusters[0].GenomeCount, Is.EqualTo(3),
                "All 3 genomes contribute a member -> GenomeCount 3 (INV-05).");
            Assert.That(clusters[0].AverageIdentity, Is.EqualTo(1.0).Within(1e-10),
                "Identical members -> mean pairwise global identity 1.0 (CD-HIT -G 1).");
        });
    }

    // M2 — Three fully distinct sequences -> three singleton clusters.
    [Test]
    public void ClusterGenes_DisjointSequences_SeparateClusters()
    {
        var genomes = Genomes(
            ("g1", new[] { ("a", S), ("b", Diff), ("c", "GGGGGGGG") }));

        var clusters = PanGenomeAnalyzer.ClusterGenes(genomes, identityThreshold: 0.9).ToList();

        Assert.That(clusters, Has.Count.EqualTo(3),
            "Pairwise identity 0.0 between all three distinct sequences -> 3 singleton clusters (CD-HIT).");
    }

    // M3 — Substitution identity 0.875 < threshold 0.9 -> two clusters.
    [Test]
    public void ClusterGenes_SubstitutionBelowThreshold_SeparateClusters()
    {
        var genomes = Genomes(
            ("g1", new[] { ("a", S) }),
            ("g2", new[] { ("b", S1Sub) }));

        var clusters = PanGenomeAnalyzer.ClusterGenes(genomes, identityThreshold: 0.9).ToList();

        Assert.That(clusters, Has.Count.EqualTo(2),
            "Identity 7/8=0.875 is below the 0.9 cutoff -> the two genes do not cluster (CD-HIT -c).");
    }

    // M4 — Substitution identity exactly at threshold (0.875 vs threshold 0.875) -> inclusive '>=' groups them.
    [Test]
    public void ClusterGenes_SubstitutionAtThreshold_SingleCluster()
    {
        var genomes = Genomes(
            ("g1", new[] { ("a", S) }),
            ("g2", new[] { ("b", S1Sub) }));

        var clusters = PanGenomeAnalyzer.ClusterGenes(genomes, identityThreshold: 0.875).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(clusters, Has.Count.EqualTo(1),
                "Identity 0.875 meets the inclusive 0.875 cutoff -> single cluster (CD-HIT '>=' rule).");
            Assert.That(clusters[0].GenomeCount, Is.EqualTo(2),
                "Both genomes contribute -> GenomeCount 2.");
        });
    }

    // M5 — Length difference, full identity over the shorter length; representative is the longest member.
    [Test]
    public void ClusterGenes_LengthDifferenceFullIdentity_GroupsAndPicksLongestRepresentative()
    {
        // SLong (12) and S (8) share an 8/8 identical prefix -> identity over shorter(8) = 1.0.
        var genomes = Genomes(
            ("g1", new[] { ("short", S) }),
            ("g2", new[] { ("long", SLong) }));

        var clusters = PanGenomeAnalyzer.ClusterGenes(genomes, identityThreshold: 1.0).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(clusters, Has.Count.EqualTo(1),
                "Global identity over the shorter length is 1.0 -> one cluster even at threshold 1.0 (CD-HIT shorter-length denominator).");
            Assert.That(clusters[0].ConsensusSequence, Is.EqualTo(SLong),
                "CD-HIT sorts long->short; the longest member (SLong) is the representative/consensus (INV-05).");
            Assert.That(clusters[0].GeneIds[0], Is.EqualTo("long"),
                "The representative (longest) member is listed first.");
        });
    }

    // M6 — Greedy hand-derived dataset (Evidence): R,Q1,Q2,Q3 at threshold 0.8 -> exactly 2 clusters.
    [Test]
    public void ClusterGenes_GreedyHandDerivedDataset_TwoClusters()
    {
        // Flatten order R,Q1,Q2,Q3. Long->short puts Q2(12) first as representative.
        // Q2 rep; R(10/10=1.0) joins; Q1(9/10=0.9>=0.8) joins; Q3(0/10=0.0) -> new cluster.
        var genomes = Genomes(
            ("g1", new[]
            {
                ("R", "AAAAAAAAAA"),    // 10
                ("Q1", "AAAAAAAAAT"),   // 10, 9/10=0.9 vs R/Q2
                ("Q2", "AAAAAAAAAAAA"), // 12 -> representative
                ("Q3", "CCCCCCCCCC"),   // 10, 0.0
            }));

        var clusters = PanGenomeAnalyzer.ClusterGenes(genomes, identityThreshold: 0.8).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(clusters, Has.Count.EqualTo(2),
                "Greedy long->short: {Q2,R,Q1} group (>=0.8 to Q2) and {Q3} alone -> 2 clusters (CD-HIT Algorithm wiki).");
            var big = clusters.OrderByDescending(c => c.GeneIds.Count).First();
            Assert.That(big.GeneIds, Has.Count.EqualTo(3),
                "The Q2 cluster contains Q2, R and Q1.");
            Assert.That(big.ConsensusSequence, Is.EqualTo("AAAAAAAAAAAA"),
                "Q2 (longest, 12 bp) is the representative.");
        });
    }

    // M7 — Lowering the threshold merges near-identical sequences that a higher threshold separates.
    [Test]
    public void ClusterGenes_LoweredThreshold_MergesNearIdentical()
    {
        var genomes = Genomes(
            ("g1", new[] { ("a", S) }),       // ATGCATGC
            ("g2", new[] { ("b", S1Sub) }));  // ATGCATGG, identity 0.875

        var at09 = PanGenomeAnalyzer.ClusterGenes(genomes, identityThreshold: 0.9).Count();
        var at07 = PanGenomeAnalyzer.ClusterGenes(genomes, identityThreshold: 0.7).Count();

        Assert.Multiple(() =>
        {
            Assert.That(at09, Is.EqualTo(2),
                "0.875 < 0.9 -> two clusters at threshold 0.9.");
            Assert.That(at07, Is.EqualTo(1),
                "0.875 >= 0.7 -> one merged cluster at threshold 0.7 (CD-HIT -c controls granularity).");
        });
    }

    // M8 — AverageIdentity of a 2-member cluster equals the global identity 0.875 (Evidence worked value).
    [Test]
    public void ClusterGenes_TwoMembers_AverageIdentityMatchesGlobalIdentity()
    {
        var genomes = Genomes(
            ("g1", new[] { ("a", S) }),       // ATGCATGC
            ("g2", new[] { ("b", S1Sub) }));  // ATGCATGG -> 7/8 = 0.875

        var clusters = PanGenomeAnalyzer.ClusterGenes(genomes, identityThreshold: 0.8).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(clusters, Has.Count.EqualTo(1), "Both members cluster at threshold 0.8.");
            Assert.That(clusters[0].AverageIdentity, Is.EqualTo(0.875).Within(1e-10),
                "Single pair -> AverageIdentity = global identity = 7/8 = 0.875 (CD-HIT -G 1).");
        });
    }

    // M8b — AverageIdentity of a 3-member cluster is the MEAN of all 3 pairwise global
    // identities (not a single pair). Members Q2(AAAAAAAAAAAA,12), R(AAAAAAAAAA,10),
    // Q1(AAAAAAAAAT,10): id(Q2,R)=10/10=1.0, id(Q2,Q1)=9/10=0.9, id(R,Q1)=9/10=0.9
    // -> mean = (1.0+0.9+0.9)/3 = 2.8/3 = 0.9333... (CD-HIT -G 1 global identity, mean over
    // all C(3,2) pairs). Exercises the multi-pair averaging path (sum/pairs) with an exact
    // hand-computed value, not just the 2-member single-pair shortcut (M8).
    [Test]
    public void ClusterGenes_ThreeMembers_AverageIdentityIsMeanOfAllPairs()
    {
        var genomes = Genomes(
            ("g1", new[] { ("R", "AAAAAAAAAA") }),     // 10
            ("g2", new[] { ("Q1", "AAAAAAAAAT") }),    // 10, 9/10 = 0.9
            ("g3", new[] { ("Q2", "AAAAAAAAAAAA") })); // 12 -> representative

        var clusters = PanGenomeAnalyzer.ClusterGenes(genomes, identityThreshold: 0.8).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(clusters, Has.Count.EqualTo(1),
                "All three group at threshold 0.8 (>= 0.9 to the Q2 representative).");
            Assert.That(clusters[0].GeneIds, Has.Count.EqualTo(3), "Cluster holds R, Q1 and Q2.");
            Assert.That(clusters[0].AverageIdentity, Is.EqualTo(2.8 / 3.0).Within(1e-10),
                "Mean of the 3 pairwise global identities (1.0+0.9+0.9)/3 = 0.9333... (CD-HIT -G 1).");
        });
    }

    // M9 — Partition invariant: every input gene appears in exactly one cluster (INV-01/INV-02).
    [Test]
    public void ClusterGenes_MixedInput_ClustersPartitionAllGenes()
    {
        var genomes = Genomes(
            ("g1", new[] { ("a", S), ("b", Diff) }),
            ("g2", new[] { ("c", S), ("d", "GGGGGGGG") }),
            ("g3", new[] { ("e", S1Sub) }));

        var clusters = PanGenomeAnalyzer.ClusterGenes(genomes, identityThreshold: 0.9).ToList();

        var totalMembers = clusters.Sum(c => c.GeneIds.Count);
        var distinctMembers = clusters.SelectMany(c => c.GeneIds).Distinct().Count();

        Assert.Multiple(() =>
        {
            Assert.That(totalMembers, Is.EqualTo(5),
                "INV-02: cluster sizes sum to the 5 input genes.");
            Assert.That(distinctMembers, Is.EqualTo(5),
                "INV-01: each gene appears in exactly one cluster (no gene duplicated or dropped).");
        });
    }

    #endregion

    #region Edge cases — empty / null / no-genes / singleton

    // S1 — Empty genomes -> no clusters.
    [Test]
    public void ClusterGenes_EmptyGenomes_ReturnsEmpty()
    {
        var clusters = PanGenomeAnalyzer.ClusterGenes(
            new Dictionary<string, IReadOnlyList<(string GeneId, string Sequence)>>()).ToList();

        Assert.That(clusters, Is.Empty, "No genomes -> no gene clusters.");
    }

    // S2 — Null genomes -> empty (no throw), matching the sibling ConstructPanGenome contract.
    [Test]
    public void ClusterGenes_NullGenomes_ReturnsEmpty()
    {
        var clusters = PanGenomeAnalyzer.ClusterGenes(null!).ToList();

        Assert.That(clusters, Is.Empty, "Null input -> empty enumeration (no exception), per documented contract.");
    }

    // S3 — Genome present but with an empty gene list -> no clusters.
    [Test]
    public void ClusterGenes_GenomeWithNoGenes_ReturnsEmpty()
    {
        var genomes = Genomes(("g1", System.Array.Empty<(string, string)>()));

        var clusters = PanGenomeAnalyzer.ClusterGenes(genomes).ToList();

        Assert.That(clusters, Is.Empty, "A genome with no genes contributes nothing to cluster.");
    }

    // S4 — Single gene -> one singleton cluster with AverageIdentity 1.0 (INV-06).
    [Test]
    public void ClusterGenes_SingleGene_SingletonAverageIdentityOne()
    {
        var genomes = Genomes(("g1", new[] { ("a", S) }));

        var clusters = PanGenomeAnalyzer.ClusterGenes(genomes).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(clusters, Has.Count.EqualTo(1), "One gene -> one singleton cluster.");
            Assert.That(clusters[0].GenomeCount, Is.EqualTo(1), "Singleton spans one genome.");
            Assert.That(clusters[0].AverageIdentity, Is.EqualTo(1.0).Within(1e-10),
                "INV-06: a singleton cluster has AverageIdentity 1.0 (self-identity).");
        });
    }

    // S5 — CreatePresenceAbsenceMatrix delegate: one row per genome with correct present-gene counts.
    [Test]
    public void CreatePresenceAbsenceMatrix_FromClusters_RowPerGenome()
    {
        var genomes = Genomes(
            ("g1", new[] { ("a", S), ("b", Diff) }),
            ("g2", new[] { ("c", S) }));

        var clusters = PanGenomeAnalyzer.ClusterGenes(genomes, identityThreshold: 0.9).ToList();
        var matrix = PanGenomeAnalyzer.CreatePresenceAbsenceMatrix(genomes, clusters)
            .ToDictionary(r => r.GenomeId);

        Assert.Multiple(() =>
        {
            Assert.That(matrix, Has.Count.EqualTo(2), "One presence/absence row per genome.");
            // S clusters g1.a + g2.c (identity 1.0); Diff is g1-only.
            Assert.That(matrix["g1"].PresentGenes, Is.EqualTo(2),
                "g1 carries both clusters (the shared S cluster and the Diff cluster).");
            Assert.That(matrix["g2"].PresentGenes, Is.EqualTo(1),
                "g2 carries only the shared S cluster.");
        });
    }

    // S6 — Empty-sequence identity branches: two empty sequences are identical (1.0) and
    // cluster together; an empty and a non-empty sequence are disjoint (0.0).
    [Test]
    public void ClusterGenes_EmptyAndNonEmptySequences_IdentityBranchesHonoured()
    {
        var genomes = Genomes(
            ("g1", new[] { ("empty1", ""), ("real", S) }),
            ("g2", new[] { ("empty2", "") }));

        var clusters = PanGenomeAnalyzer.ClusterGenes(genomes, identityThreshold: 1.0).ToList();

        Assert.Multiple(() =>
        {
            // empty1 + empty2 (identity 1.0) cluster; S stands alone (0.0 vs empty).
            Assert.That(clusters, Has.Count.EqualTo(2),
                "Two empty sequences cluster (identity 1.0); the non-empty sequence is separate (0.0).");
            var emptyCluster = clusters.First(c => c.ConsensusSequence.Length == 0);
            Assert.That(emptyCluster.GenomeCount, Is.EqualTo(2),
                "Both empty-sequence genes (across 2 genomes) form one cluster (identity 1.0).");
        });
    }

    #endregion

    #region Determinism and property-based (O(g^2 . s) algorithm)

    // C1 — Determinism: identical input yields identical clustering across calls.
    [Test]
    public void ClusterGenes_CalledTwice_ProducesIdenticalClustering()
    {
        var genomes = Genomes(
            ("g1", new[] { ("a", S), ("b", Diff) }),
            ("g2", new[] { ("c", S), ("d", S1Sub) }));

        var c1 = PanGenomeAnalyzer.ClusterGenes(genomes, identityThreshold: 0.9).ToList();
        var c2 = PanGenomeAnalyzer.ClusterGenes(genomes, identityThreshold: 0.9).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(c2, Has.Count.EqualTo(c1.Count), "Cluster count is deterministic.");
            for (int i = 0; i < c1.Count; i++)
                Assert.That(c2[i].GeneIds, Is.EqualTo(c1[i].GeneIds),
                    $"Cluster {i} membership is deterministic (stable long->short order).");
        });
    }

    // C2 — Property (INV-03): all AverageIdentity values lie in [0,1] across structured inputs.
    [Test]
    public void ClusterGenes_VariedPairs_AverageIdentityInUnitInterval()
    {
        string[] pool = { S, S1Sub, SLong, Diff, "GGGGGGGG", "TTTTTTTT" };

        Assert.Multiple(() =>
        {
            for (double t = 0.5; t <= 1.0 + 1e-9; t += 0.1)
            {
                var genomes = new Dictionary<string, IReadOnlyList<(string GeneId, string Sequence)>>();
                for (int g = 0; g < pool.Length; g++)
                    genomes[$"g{g}"] = new[] { ($"x{g}", pool[g]) };

                foreach (var cluster in PanGenomeAnalyzer.ClusterGenes(genomes, identityThreshold: t))
                    Assert.That(cluster.AverageIdentity, Is.InRange(0.0, 1.0),
                        $"INV-03: AverageIdentity must lie in [0,1] (threshold={t:F1}).");
            }
        });
    }

    #endregion
}
