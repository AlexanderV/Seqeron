using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Metagenomics;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Fuzz tests for the PanGenome area — gene clustering into homolog families
/// (PANGEN-CLUSTER-001), the CD-HIT greedy incremental clustering model.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate, boundary and out-of-domain inputs to a unit and
/// asserts that the code NEVER fails in an undisciplined way: no hang or infinite
/// loop, no state corruption, no nonsense output (a gene assigned to zero or to
/// more than one cluster, a cluster whose members are below the identity threshold
/// to its representative, a non-deterministic clustering), and no *unhandled*
/// runtime exception (DivideByZero on a zero-length gene where identity divides by
/// the shorter length, IndexOutOfRange on a single gene, NullReferenceException on
/// a null gene list). Every input must resolve to EITHER a well-defined,
/// theory-correct partition, OR a *documented, intentional* degenerate result
/// (null genomes → empty enumeration, no throw). A raw runtime exception, a hang,
/// a malformed partition, or a non-deterministic result is a bug, not a passing
/// test. — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: PANGEN-CLUSTER-001 — Gene clustering (homolog grouping by identity)
/// Checklist: docs/checklists/03_FUZZING.md, row 190.
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — the degenerate clustering boundaries called out
///          in the checklist row: a SINGLE gene (one gene → one singleton cluster),
///          ALL-IDENTICAL genes (every gene collapses into one cluster), and the
///          IDENTITY EDGE (genes at exactly the threshold land per the inclusive
///          `>=` rule, just-below split off, just-at/above join).
/// — docs/checklists/03_FUZZING.md §Description (strategy codes:
///    "BE = Boundary Exploitation (граничні значення: 0, -1, MaxInt, empty)").
///
/// ───────────────────────────────────────────────────────────────────────────
/// The clustering contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// PanGenomeAnalyzer.ClusterGenes implements the CD-HIT greedy incremental
/// clustering model (Li &amp; Godzik 2006; Gene_Clustering.md §2.2): genes are
/// flattened across genomes, sorted longest→shortest (stable), the longest seeds
/// the first cluster as its representative, and each subsequent gene joins the
/// FIRST existing representative whose GLOBAL SEQUENCE IDENTITY meets the
/// threshold (inclusive `>=`), otherwise it becomes a new representative. Global
/// identity = identical residues in the ungapped positional alignment over
/// min(|s1|,|s2|) (CD-HIT -G 1 default); two empty sequences are identical (1.0),
/// one empty + one non-empty is 0.0 (§2.2 eqn, §3.3). The API entry under test is
///   PanGenomeAnalyzer.ClusterGenes(
///       IReadOnlyDictionary&lt;string, IReadOnlyList&lt;(string GeneId, string Sequence)&gt;&gt; genomes,
///       double identityThreshold = 0.9)
///   (src/Seqeron/Algorithms/Seqeron.Genomics.Metagenomics/PanGenomeAnalyzer.cs
///    lines 210–306; private CalculateSequenceIdentity lines 314–343).
///
/// THE DOCUMENTED INVARIANTS (Gene_Clustering.md §2.4):
///   • INV-01: every input gene is in EXACTLY one cluster (a partition).
///   • INV-02: Σ cluster sizes = number of input genes.
///   • INV-03: 0 ≤ identity ≤ 1; identical = 1.0, disjoint = 0.0.
///   • INV-04: members of a cluster have identity ≥ threshold to its representative.
///   • INV-05: the cluster representative is the LONGEST member.
///   • INV-06: a singleton cluster has AverageIdentity = 1.0.
/// The documented edge cases (§6.1):
///   • null genomes → empty enumeration (no throw); empty genomes → empty.
///   • genome with empty gene list → contributes nothing; null inner list skipped.
///   • both sequences empty → identity 1.0; one empty → 0.0.
///   • idThreshold = 1.0 → only exact-identity-over-shorter-length sequences cluster.
///   • singleton cluster → AverageIdentity 1.0 (INV-06).
///
/// The three BE checklist targets map to these documented behaviours:
///   • single gene    → ONE singleton cluster (INV-01/06): the smallest non-empty
///                      input. The hazard probed is a crash / IndexOutOfRange /
///                      DivideByZero on a 1-gene clustering; the gene seeds the only
///                      cluster, AverageIdentity = 1.0, and it is assigned exactly once.
///   • all-identical  → ONE cluster holding ALL genes (INV-01/02/04), deterministically,
///                      with NO double-assignment: each gene after the representative
///                      has identity 1.0 ≥ threshold so all collapse into the seed.
///   • identity edge  → genes at EXACTLY the threshold join the representative (the
///                      inclusive `>=` rule, §4.2 / INV-04: NO off-by-one); a gene a
///                      hair BELOW the threshold splits into its own cluster; a gene at
///                      or above joins. The DivideByZero hazard at the edge of the
///                      length axis (zero-length gene → denominator = shorter length)
///                      is probed by the empty-sequence cases (§3.3 defines 1.0 / 0.0).
/// A positive-sanity suite pins the documented golden vectors (§7.1): the two-group
/// example (a,c identical + b distinct → two clusters), the 0.9-edge walk-through
/// (R seeds, Q1=9/10=0.9 joins, Q3=0/10 splits → two clusters), all-identical → one
/// cluster, single gene → singleton, and idThreshold = 1.0 behaviour.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class PanGenomeClusterFuzzTests
{
    #region Helpers

    /// <summary>The genomes-dictionary type the API expects.</summary>
    private static IReadOnlyDictionary<string, IReadOnlyList<(string GeneId, string Sequence)>> Genomes(
        params (string genomeId, (string geneId, string seq)[] genes)[] entries)
    {
        var dict = new Dictionary<string, IReadOnlyList<(string GeneId, string Sequence)>>(StringComparer.Ordinal);
        foreach (var (genomeId, genes) in entries)
            dict[genomeId] = genes.Select(g => (g.geneId, g.seq)).ToList();
        return dict;
    }

    /// <summary>Single-genome convenience: one genome "g" holding the given (id, seq) genes.</summary>
    private static IReadOnlyDictionary<string, IReadOnlyList<(string GeneId, string Sequence)>> OneGenome(
        params (string geneId, string seq)[] genes) => Genomes(("g", genes));

    /// <summary>Every gene id flattened across all input genomes (the domain of the partition).</summary>
    private static List<string> AllGeneIds(
        IReadOnlyDictionary<string, IReadOnlyList<(string GeneId, string Sequence)>> genomes) =>
        genomes.Values.SelectMany(genes => genes.Select(g => g.GeneId)).ToList();

    /// <summary>Maps every gene id to its sequence for the independent identity re-check.</summary>
    private static Dictionary<string, string> SequenceOf(
        IReadOnlyDictionary<string, IReadOnlyList<(string GeneId, string Sequence)>> genomes)
    {
        var map = new Dictionary<string, string>(StringComparer.Ordinal);
        foreach (var genes in genomes.Values)
            foreach (var (geneId, seq) in genes)
                map[geneId] = seq;
        return map;
    }

    /// <summary>
    /// Independent re-implementation of CD-HIT global sequence identity for the test's
    /// INV-04 cross-check (identical residues over min length; both-empty = 1.0, one-empty
    /// = 0.0). Mirrors Gene_Clustering.md §2.2 without reusing the unit under test.
    /// </summary>
    private static double Identity(string a, string b)
    {
        a ??= string.Empty;
        b ??= string.Empty;
        if (a.Length == 0 && b.Length == 0) return 1.0;
        if (a.Length == 0 || b.Length == 0) return 0.0;
        int shorter = Math.Min(a.Length, b.Length);
        int identical = 0;
        for (int i = 0; i < shorter; i++)
            if (a[i] == b[i]) identical++;
        return (double)identical / shorter;
    }

    /// <summary>
    /// A well-formed clustering result satisfies the documented structural invariants
    /// over the INPUT gene set (Gene_Clustering.md §2.4):
    ///   • INV-01/02: the clusters PARTITION the input genes — every input gene appears
    ///     in EXACTLY one cluster (no gene assigned to 0 or to &gt;1 clusters), and the
    ///     union of cluster members equals the input gene set (Σ sizes = #genes).
    ///   • cluster ids are the documented "cluster_1", "cluster_2", … emission order.
    ///   • GenomeCount equals the number of DISTINCT genomes among a cluster's members.
    ///   • INV-04: every member has identity ≥ threshold to the cluster representative
    ///     (the FIRST member; INV-05 the longest by construction).
    ///   • INV-06: a singleton cluster has AverageIdentity = 1.0.
    /// This is the strong partition net re-derived independently of the unit.
    /// </summary>
    private static void AssertWellFormedPartition(
        IReadOnlyList<PanGenomeAnalyzer.GeneCluster> clusters,
        IReadOnlyDictionary<string, IReadOnlyList<(string GeneId, string Sequence)>> genomes,
        double threshold)
    {
        var inputGeneIds = AllGeneIds(genomes);
        var seqOf = SequenceOf(genomes);

        // INV-01: no gene in two clusters; every assigned gene is a real input gene.
        var assigned = new List<string>();
        foreach (var c in clusters)
        {
            c.GeneIds.Should().NotBeNull("a cluster's member list must never be null");
            c.GeneIds.Should().NotBeEmpty("a reported cluster must have at least one member");
            c.ClusterId.Should().StartWith("cluster_", "cluster ids follow the documented cluster_N scheme (§3.2)");
            assigned.AddRange(c.GeneIds);
        }
        assigned.Should().OnlyHaveUniqueItems("no gene may be assigned to more than one cluster (INV-01)");

        // INV-01/02: the assigned set EQUALS the input set (no gene dropped, none invented).
        assigned.Should().BeEquivalentTo(inputGeneIds,
            "the clusters must partition exactly the input genes: every gene in exactly one cluster (INV-01/02)");

        foreach (var c in clusters)
        {
            // GenomeCount = distinct contributing genomes.
            int distinctGenomes = c.GenomeIds.Distinct().Count();
            c.GenomeCount.Should().Be(distinctGenomes, "GenomeCount is the number of distinct genomes (§3.2)");
            c.GenomeCount.Should().BeGreaterThan(0, "an occupied cluster spans at least one genome");

            // INV-05: representative is the FIRST member; INV-04: each member ≥ threshold to it.
            string repSeq = seqOf[c.GeneIds[0]];
            foreach (var memberId in c.GeneIds)
            {
                double id = Identity(seqOf[memberId], repSeq);
                id.Should().BeGreaterThanOrEqualTo(threshold - 1e-12,
                    $"member {memberId} must have identity ≥ threshold {threshold} to the representative (INV-04)");
            }

            // INV-06: singletons report AverageIdentity 1.0.
            if (c.GeneIds.Count == 1)
                c.AverageIdentity.Should().BeApproximately(1.0, 1e-12, "a singleton is 100% identical to itself (INV-06)");

            // INV-03 shape: AverageIdentity stays in [0,1].
            c.AverageIdentity.Should().BeInRange(0.0, 1.0 + 1e-12, "average pairwise identity lies in [0,1] (INV-03)");
        }
    }

    /// <summary>Canonical signature of a clustering: each cluster's member ids sorted, then the
    /// set of those signatures — order- and id-stable so it can compare two runs / two orderings.</summary>
    private static HashSet<string> Signature(IEnumerable<PanGenomeAnalyzer.GeneCluster> clusters) =>
        clusters
            .Select(c => string.Join(",", c.GeneIds.OrderBy(x => x, StringComparer.Ordinal)))
            .ToHashSet(StringComparer.Ordinal);

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  PANGEN-CLUSTER-001 — Gene clustering (CD-HIT greedy) : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region PANGEN-CLUSTER-001 — Gene clustering (homolog grouping by identity)

    #region BE — Boundary: single gene (one gene → one singleton cluster)

    /// <summary>
    /// BE: a single gene is the smallest non-empty input — the floor of the gene-count
    /// axis. It must seed exactly ONE cluster of which it is the sole member, with
    /// AverageIdentity 1.0 (Gene_Clustering.md INV-01/06), assigned exactly once. The
    /// hazard probed is a crash / IndexOutOfRange / DivideByZero on a 1-gene clustering
    /// (no representatives loop iterations, an average over zero pairs); none must occur.
    /// </summary>
    [Test]
    [CancelAfter(20000)]
    public void ClusterGenes_SingleGene_OneSingletonClusterNoCrash()
    {
        var genomes = OneGenome(("a", "ATGCATGCAT"));

        List<PanGenomeAnalyzer.GeneCluster> clusters = null!;
        FluentActions.Invoking(() => clusters = PanGenomeAnalyzer.ClusterGenes(genomes, 0.9).ToList())
            .Should().NotThrow("a single gene must seed one cluster, never crash on a clustering-of-one");

        clusters.Should().HaveCount(1, "one gene → exactly one cluster (INV-01)");
        clusters[0].GeneIds.Should().ContainSingle().Which.Should().Be("a", "the lone gene is the only member");
        clusters[0].AverageIdentity.Should().BeApproximately(1.0, 1e-12, "a singleton is 100% self-identical (INV-06)");
        AssertWellFormedPartition(clusters, genomes, 0.9);
    }

    /// <summary>
    /// BE: a single gene with the EMPTY sequence is the corner where the identity
    /// denominator (the shorter length) would be ZERO. It must NOT trigger a
    /// DivideByZero: a singleton is never compared to a representative (it IS the
    /// representative), AverageIdentity is the defined 1.0, and the gene is clustered
    /// once. This pins the zero-length boundary of the single-gene case.
    /// </summary>
    [Test]
    [CancelAfter(20000)]
    public void ClusterGenes_SingleEmptySequenceGene_NoDivideByZero()
    {
        var genomes = OneGenome(("e", ""));

        List<PanGenomeAnalyzer.GeneCluster> clusters = null!;
        FluentActions.Invoking(() => clusters = PanGenomeAnalyzer.ClusterGenes(genomes, 0.9).ToList())
            .Should().NotThrow("a single zero-length gene must not divide by a zero shorter-length");

        clusters.Should().HaveCount(1, "the empty-sequence gene still seeds exactly one cluster");
        clusters[0].GeneIds.Should().ContainSingle().Which.Should().Be("e");
        clusters[0].AverageIdentity.Should().BeApproximately(1.0, 1e-12, "singleton self-identity is 1.0 (INV-06)");
    }

    #endregion

    #region BE — Boundary: all-identical (collapse into one cluster)

    /// <summary>
    /// BE: when every gene is identical, all genes have identity 1.0 ≥ any threshold to
    /// the representative, so they ALL collapse into ONE cluster (Gene_Clustering.md
    /// INV-01/02/04). We pin EXACTLY one cluster holding ALL gene ids, with no gene
    /// double-assigned (INV-01) and AverageIdentity 1.0 (all pairs identical). This is
    /// the upper boundary of the merge axis and the canonical double-assignment hazard.
    /// </summary>
    [Test]
    [CancelAfter(20000)]
    public void ClusterGenes_AllIdentical_CollapseIntoOneCluster()
    {
        const string seq = "ATGCATGCATGC";
        var genomes = Genomes(
            ("g1", new[] { ("a", seq), ("b", seq) }),
            ("g2", new[] { ("c", seq) }),
            ("g3", new[] { ("d", seq) }));

        var clusters = PanGenomeAnalyzer.ClusterGenes(genomes, 0.9).ToList();

        AssertWellFormedPartition(clusters, genomes, 0.9);
        clusters.Should().HaveCount(1, "all-identical genes collapse into exactly one cluster (INV-01)");
        clusters[0].GeneIds.Should().BeEquivalentTo(new[] { "a", "b", "c", "d" },
            "every identical gene joins the one cluster — none dropped, none double-assigned (INV-01/02)");
        clusters[0].GenomeCount.Should().Be(3, "the cluster spans the three distinct genomes g1,g2,g3");
        clusters[0].AverageIdentity.Should().BeApproximately(1.0, 1e-12, "all pairs are identical → mean identity 1.0");
    }

    /// <summary>
    /// BE: all-identical clustering must be DETERMINISTIC and order-robust. Repeated runs
    /// and a reordering of the genomes must produce the SAME single cluster with the same
    /// member set — the long→short sort is stable and the greedy join is first-match, so
    /// there is no hash-iteration nondeterminism (§5.2). Pins that the upper-boundary
    /// collapse is reproducible, not an accident of dictionary order.
    /// </summary>
    [Test]
    [CancelAfter(20000)]
    public void ClusterGenes_AllIdentical_DeterministicAcrossRunsAndOrder()
    {
        const string seq = "GGGGCCCCAAAA";
        var forward = Genomes(
            ("g1", new[] { ("a", seq), ("b", seq) }),
            ("g2", new[] { ("c", seq) }));
        var reordered = Genomes(
            ("g2", new[] { ("c", seq) }),
            ("g1", new[] { ("b", seq), ("a", seq) }));

        var first = Signature(PanGenomeAnalyzer.ClusterGenes(forward, 0.9));
        for (int i = 0; i < 5; i++)
            Signature(PanGenomeAnalyzer.ClusterGenes(forward, 0.9)).Should().BeEquivalentTo(first,
                "repeated runs on identical input give the identical partition (§5.2 determinism)");

        // All-identical genes collapse regardless of input ordering → same single cluster set.
        Signature(PanGenomeAnalyzer.ClusterGenes(reordered, 0.9)).Should().BeEquivalentTo(first,
            "all-identical genes collapse into one cluster irrespective of genome order");
    }

    #endregion

    #region BE — Boundary: identity edge (inclusive >= threshold, no off-by-one)

    /// <summary>
    /// BE: the identity EDGE — a gene whose identity to the representative is EXACTLY the
    /// threshold must JOIN it (inclusive `>=`, Gene_Clustering.md §4.2 / INV-04: no
    /// off-by-one). Representative R = 10×A; query Q differs in exactly one position
    /// (9/10 = 0.9) so identity = 0.9 EXACTLY. At threshold 0.9 they must share ONE
    /// cluster; this is the load-bearing "just-at-the-boundary joins" check.
    /// </summary>
    [Test]
    [CancelAfter(20000)]
    public void ClusterGenes_IdentityExactlyAtThreshold_JoinsInclusive()
    {
        // R = AAAAAAAAAA, Q = AAAAAAAAAT → 9/10 identical = 0.9 exactly.
        var genomes = OneGenome(("R", "AAAAAAAAAA"), ("Q", "AAAAAAAAAT"));

        var clusters = PanGenomeAnalyzer.ClusterGenes(genomes, 0.9).ToList();

        AssertWellFormedPartition(clusters, genomes, 0.9);
        clusters.Should().HaveCount(1, "identity 0.9 == threshold 0.9 → inclusive `>=` joins them into ONE cluster (no off-by-one, §4.2)");
        clusters[0].GeneIds.Should().BeEquivalentTo(new[] { "R", "Q" }, "the at-threshold gene joins the representative's cluster");
    }

    /// <summary>
    /// BE: the identity edge from the OTHER side — a gene a hair BELOW the threshold must
    /// SPLIT into its own cluster. Representative R = 10×A; query Q differs in two
    /// positions (8/10 = 0.8 &lt; 0.9). At threshold 0.9 they must form TWO clusters,
    /// pinning that the `>=` boundary is strict below it (no spurious merge).
    /// </summary>
    [Test]
    [CancelAfter(20000)]
    public void ClusterGenes_IdentityJustBelowThreshold_SplitsIntoTwoClusters()
    {
        // R = AAAAAAAAAA, Q = AAAAAAAATT → 8/10 identical = 0.8 < 0.9.
        var genomes = OneGenome(("R", "AAAAAAAAAA"), ("Q", "AAAAAAAATT"));

        var clusters = PanGenomeAnalyzer.ClusterGenes(genomes, 0.9).ToList();

        AssertWellFormedPartition(clusters, genomes, 0.9);
        clusters.Should().HaveCount(2, "identity 0.8 < threshold 0.9 → the gene is below the cutoff and starts its own cluster");
        clusters.SelectMany(c => c.GeneIds).Should().BeEquivalentTo(new[] { "R", "Q" }, "both genes are still partitioned exactly once");
    }

    /// <summary>
    /// BE: the identity edge at the EXTREME threshold idThreshold = 1.0 — only sequences
    /// that are 100% identical over the shorter length cluster (Gene_Clustering.md §6.1).
    /// A gene at 0.9 to the representative must now SPLIT (0.9 &lt; 1.0), while an exact
    /// duplicate joins. Pins the boundary at the top of the threshold axis: `>=` 1.0 is
    /// satisfied only by exact identity, never by a near-match.
    /// </summary>
    [Test]
    [CancelAfter(20000)]
    public void ClusterGenes_ThresholdOne_OnlyExactDuplicatesCluster()
    {
        // R and Dup are identical; Near is 9/10 = 0.9 to R.
        var genomes = OneGenome(
            ("R", "AAAAAAAAAA"),
            ("Dup", "AAAAAAAAAA"),
            ("Near", "AAAAAAAAAT"));

        var clusters = PanGenomeAnalyzer.ClusterGenes(genomes, 1.0).ToList();

        AssertWellFormedPartition(clusters, genomes, 1.0);
        clusters.Should().HaveCount(2, "at threshold 1.0 only the exact duplicate joins R; the 0.9 near-match splits off");
        var withR = clusters.Single(c => c.GeneIds.Contains("R"));
        withR.GeneIds.Should().BeEquivalentTo(new[] { "R", "Dup" }, "exact duplicates cluster at 1.0 (§6.1)");
        clusters.Single(c => c.GeneIds.Contains("Near")).GeneIds.Should().BeEquivalentTo(new[] { "Near" },
            "a 0.9 match does NOT meet `>= 1.0` → its own cluster (no off-by-one at the extreme threshold)");
    }

    /// <summary>
    /// BE: the identity edge at the FLOOR threshold idThreshold = 0.0 — every gene meets
    /// `identity &gt;= 0.0` against the first representative (identity is ≥ 0 by INV-03),
    /// so ALL genes collapse into a single cluster regardless of sequence. Pins the
    /// bottom of the threshold axis: a zero cutoff merges everything, deterministically
    /// and with no double-assignment.
    /// </summary>
    [Test]
    [CancelAfter(20000)]
    public void ClusterGenes_ThresholdZero_EverythingCollapsesToOneCluster()
    {
        var genomes = OneGenome(
            ("a", "AAAAAAAAAA"),
            ("b", "CCCCCCCCCC"),   // 0.0 identity to a, but 0.0 >= 0.0
            ("c", "GGGGGGGGGG"));

        var clusters = PanGenomeAnalyzer.ClusterGenes(genomes, 0.0).ToList();

        AssertWellFormedPartition(clusters, genomes, 0.0);
        clusters.Should().HaveCount(1, "threshold 0.0: identity ≥ 0 always holds → all genes join the first representative");
        clusters[0].GeneIds.Should().BeEquivalentTo(new[] { "a", "b", "c" }, "everything collapses, none double-assigned (INV-01)");
    }

    #endregion

    #region BE — Boundary: zero-length genes (identity denominator = shorter length)

    /// <summary>
    /// BE: a ZERO-LENGTH gene is the boundary of the length axis where the identity
    /// denominator (the shorter length) is 0. The documented convention (Gene_Clustering.md
    /// §2.2/§3.3) is: two empty sequences are identical (1.0), one empty + one non-empty
    /// is 0.0 — NEVER a DivideByZero. With multiple empty-sequence genes at any threshold
    /// ≤ 1.0 they must collapse together (1.0 ≥ threshold) and the non-empty gene splits
    /// off (0.0 to an empty representative). The hazard probed is DivideByZero on the
    /// min-length-0 denominator.
    /// </summary>
    [Test]
    [CancelAfter(20000)]
    public void ClusterGenes_EmptyAndNonEmptyGenes_NoDivideByZero()
    {
        var genomes = OneGenome(
            ("nonEmpty", "ATGCATGC"),  // longest → seeds cluster 1 (rep)
            ("empty1", ""),
            ("empty2", ""));

        List<PanGenomeAnalyzer.GeneCluster> clusters = null!;
        FluentActions.Invoking(() => clusters = PanGenomeAnalyzer.ClusterGenes(genomes, 0.9).ToList())
            .Should().NotThrow("identity over a zero shorter-length must use the 1.0/0.0 convention, never DivideByZero");

        AssertWellFormedPartition(clusters, genomes, 0.9);
        clusters.Should().HaveCount(2, "the two empty genes are mutually identical (1.0) and split from the non-empty gene (0.0)");
        var emptyCluster = clusters.Single(c => c.GeneIds.Contains("empty1"));
        emptyCluster.GeneIds.Should().BeEquivalentTo(new[] { "empty1", "empty2" },
            "two empty sequences are identical (1.0 ≥ 0.9) → one cluster (§3.3)");
        clusters.Single(c => c.GeneIds.Contains("nonEmpty")).GeneIds.Should().BeEquivalentTo(new[] { "nonEmpty" },
            "an empty vs non-empty pair is 0.0 identity → separate clusters (§3.3)");
    }

    #endregion

    #region BE — Boundary: degenerate / empty inputs (no throw)

    /// <summary>
    /// BE: null genomes is the documented degenerate input — it yields an EMPTY
    /// enumeration with NO throw (Gene_Clustering.md §3.3 / §6.1, matching the sibling
    /// ConstructPanGenome contract). Pins that the lazy iterator does not defer a null
    /// deref to enumeration time.
    /// </summary>
    [Test]
    public void ClusterGenes_NullGenomes_EmptyNoThrow()
    {
        List<PanGenomeAnalyzer.GeneCluster> clusters = null!;
        FluentActions.Invoking(() => clusters = PanGenomeAnalyzer.ClusterGenes(null!, 0.9).ToList())
            .Should().NotThrow("null genomes is a documented degenerate input → empty enumeration, never a throw (§3.3)");
        clusters.Should().BeEmpty("null genomes → empty (§6.1)");
    }

    /// <summary>
    /// BE: an empty genomes dictionary, and a genome whose gene list is empty or null,
    /// all contribute zero genes → empty result with no crash (§6.1: "genome with empty
    /// gene list → contributes nothing; null inner list skipped"). Pins the floor of the
    /// genome/gene-count axes.
    /// </summary>
    [Test]
    public void ClusterGenes_EmptyAndNullGeneLists_EmptyNoCrash()
    {
        PanGenomeAnalyzer.ClusterGenes(
            new Dictionary<string, IReadOnlyList<(string, string)>>(), 0.9)
            .Should().BeEmpty("empty genomes dictionary → empty (§6.1)");

        PanGenomeAnalyzer.ClusterGenes(
            Genomes(("g", Array.Empty<(string, string)>())), 0.9)
            .Should().BeEmpty("a genome with an empty gene list contributes nothing → empty (§6.1)");

        var withNullList = new Dictionary<string, IReadOnlyList<(string GeneId, string Sequence)>>
        {
            ["g"] = null!
        };
        List<PanGenomeAnalyzer.GeneCluster> clusters = null!;
        FluentActions.Invoking(() => clusters = PanGenomeAnalyzer.ClusterGenes(withNullList, 0.9).ToList())
            .Should().NotThrow("a null inner gene list must be skipped, never a null deref (§6.1)");
        clusters.Should().BeEmpty("the only genome's gene list is null → nothing to cluster");
    }

    #endregion

    #region Positive sanity — the documented metric on hand-built examples

    /// <summary>
    /// Positive sanity (Gene_Clustering.md §7.1 worked example): genomes g1 = {a:ATGCATGCAT,
    /// b:TTTTTTTTTT}, g2 = {c:ATGCATGCAT}. Gene c is identical to a, so {a,c} form one
    /// cluster (occupancy 2) and b is its own cluster (occupancy 1). Pins the EXACT
    /// two-cluster golden vector — the load-bearing correctness check distinguishing a
    /// theory-correct clusterer from an over- or under-merging one.
    /// </summary>
    [Test]
    [CancelAfter(20000)]
    public void ClusterGenes_DocExample_TwoClustersByIdentity()
    {
        var genomes = Genomes(
            ("g1", new[] { ("a", "ATGCATGCAT"), ("b", "TTTTTTTTTT") }),
            ("g2", new[] { ("c", "ATGCATGCAT") }));

        var clusters = PanGenomeAnalyzer.ClusterGenes(genomes, 0.9).ToList();

        AssertWellFormedPartition(clusters, genomes, 0.9);
        clusters.Should().HaveCount(2, "two distinct identity groups → two clusters (§7.1)");
        var abCluster = clusters.Single(c => c.GeneIds.Contains("a"));
        abCluster.GeneIds.Should().BeEquivalentTo(new[] { "a", "c" }, "a and c are identical → same cluster");
        abCluster.GenomeCount.Should().Be(2, "the {a,c} cluster spans both genomes (occupancy 2, §7.1)");
        clusters.Single(c => c.GeneIds.Contains("b")).GenomeCount.Should().Be(1, "b is unique to g1 (occupancy 1)");
    }

    /// <summary>
    /// Positive sanity (Gene_Clustering.md §7.1 numerical walk-through): R=10×A, Q1=9A+T
    /// (9/10 = 0.9), Q3=10×C (0/10 = 0.0). At threshold 0.9: R seeds cluster 1, Q1 joins
    /// (0.9 ≥ 0.9), Q3 starts cluster 2 → TWO clusters. Pins the documented edge example
    /// exactly: the at-threshold join AND the zero-identity split in one vector.
    /// </summary>
    [Test]
    [CancelAfter(20000)]
    public void ClusterGenes_DocNumericalWalkthrough_TwoClustersAtPoint9()
    {
        var genomes = OneGenome(
            ("R", "AAAAAAAAAA"),
            ("Q1", "AAAAAAAAAT"),   // 0.9 → joins R
            ("Q3", "CCCCCCCCCC"));  // 0.0 → new cluster

        var clusters = PanGenomeAnalyzer.ClusterGenes(genomes, 0.9).ToList();

        AssertWellFormedPartition(clusters, genomes, 0.9);
        clusters.Should().HaveCount(2, "R+Q1 cluster (0.9 ≥ 0.9) and Q3 alone (0.0) → two clusters (§7.1 walk-through)");
        clusters.Single(c => c.GeneIds.Contains("R")).GeneIds.Should().BeEquivalentTo(new[] { "R", "Q1" },
            "Q1 at exactly 0.9 joins the representative R");
        clusters.Single(c => c.GeneIds.Contains("Q3")).GeneIds.Should().BeEquivalentTo(new[] { "Q3" },
            "Q3 at 0.0 cannot join → its own cluster");
    }

    /// <summary>
    /// Positive sanity — INV-05: the cluster representative (consensus) is the LONGEST
    /// member, a consequence of the long→short processing order. A short gene and a long
    /// gene sharing identity ≥ threshold over the shorter length must cluster with the
    /// LONG one as the representative/consensus. Pins the documented representative rule.
    /// </summary>
    [Test]
    [CancelAfter(20000)]
    public void ClusterGenes_RepresentativeIsLongestMember()
    {
        // "AAAAAA" (len 6) is a prefix of "AAAAAAAAAA" (len 10): identity over shorter (6) = 6/6 = 1.0.
        var genomes = OneGenome(
            ("shortGene", "AAAAAA"),
            ("longGene", "AAAAAAAAAA"));

        var clusters = PanGenomeAnalyzer.ClusterGenes(genomes, 0.9).ToList();

        AssertWellFormedPartition(clusters, genomes, 0.9);
        clusters.Should().HaveCount(1, "the short gene is 100% identical over its length to the long gene → one cluster");
        clusters[0].ConsensusSequence.Should().Be("AAAAAAAAAA", "the longest member is the representative/consensus (INV-05)");
        clusters[0].GeneIds[0].Should().Be("longGene", "the representative (first member) is the longest gene (INV-05)");
    }

    #endregion

    #region Fuzz sweep — random inputs always partition without crashing

    /// <summary>
    /// Fuzz sweep — across many random multi-genome inputs with varied gene counts,
    /// sequence lengths (including zero-length), alphabets and thresholds (including the
    /// 0.0/1.0 extremes), ClusterGenes must ALWAYS terminate, never throw (no DivideByZero
    /// on empty genes, no IndexOutOfRange), and return a well-formed PARTITION: every input
    /// gene in exactly one cluster, every member ≥ threshold to its representative (INV-04),
    /// singletons at AverageIdentity 1.0. The partition + identity re-check is the strong
    /// independent net.
    /// </summary>
    [Test]
    [CancelAfter(60000)]
    public void ClusterGenes_RandomInputs_AlwaysWellFormedPartition()
    {
        const string alphabet = "ACGT";
        for (int seed = 0; seed < 80; seed++)
        {
            var rng = new Random(seed);
            int genomeCount = rng.Next(1, 5);          // includes the single-genome case
            double threshold = new[] { 0.0, 0.5, 0.8, 0.9, 1.0 }[rng.Next(5)]; // includes both extremes

            var entries = new List<(string, (string, string)[])>();
            int geneCounter = 0;
            for (int g = 0; g < genomeCount; g++)
            {
                int geneCount = rng.Next(0, 5);        // includes the empty gene list
                var genes = new (string, string)[geneCount];
                for (int j = 0; j < geneCount; j++)
                {
                    int len = rng.Next(0, 9);          // includes the zero-length gene
                    var chars = new char[len];
                    for (int k = 0; k < len; k++) chars[k] = alphabet[rng.Next(alphabet.Length)];
                    genes[j] = ($"gene{geneCounter++}", new string(chars));
                }
                entries.Add(($"genome{g}", genes));
            }
            var genomes = Genomes(entries.ToArray());

            List<PanGenomeAnalyzer.GeneCluster> clusters = null!;
            FluentActions.Invoking(() => clusters = PanGenomeAnalyzer.ClusterGenes(genomes, threshold).ToList())
                .Should().NotThrow($"random input must never crash (seed {seed}, threshold {threshold})");

            AssertWellFormedPartition(clusters, genomes, threshold);

            // Determinism: a second run produces the identical partition.
            Signature(clusters).Should().BeEquivalentTo(
                Signature(PanGenomeAnalyzer.ClusterGenes(genomes, threshold)),
                $"clustering is deterministic on repeated runs (seed {seed})");
        }
    }

    #endregion

    #endregion
}
