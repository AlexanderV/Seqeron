using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Metagenomics;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Metamorphic tests for the PanGenome area.
///
/// Each test encodes a metamorphic relation (MR) — a property relating the outputs of
/// multiple runs under an input transformation, with no hardcoded oracle. The relations
/// are derived from the ALGORITHM DEFINITION, not from observed output.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: PANGEN-CLUSTER-001 — gene clustering (CD-HIT-style) (PanGenome).
/// Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 190.
///
/// API under test (PanGenomeAnalyzer.ClusterGenes):
///   Greedy identity-threshold clustering: a gene joins a cluster when its identity to the
///   representative is ≥ the threshold, else starts a new cluster.
///
/// Relations (derived from the identity-threshold rule, NOT from output):
///   • MON  (lower identity ⇒ coarser clusters): a lower cutoff lets more genes merge, so the
///          number of clusters cannot increase (it gets coarser).
///   • INV  (gene order independent): the length-sorted greedy clustering produces the same cluster
///          membership regardless of input gene order (distinct-length genes make the sort unambiguous).
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Metamorphic")]
public class PanGenomeMetamorphicTests
{
    // geneA and geneB are 90% identical (1 mismatch in 10); geneC is unrelated.
    private static IReadOnlyDictionary<string, IReadOnlyList<(string, string)>> Genomes(
        IReadOnlyList<(string GeneId, string Sequence)> genes) =>
        new Dictionary<string, IReadOnlyList<(string, string)>> { ["g1"] = genes };

    private static readonly (string, string)[] ClusterGenesData =
    {
        ("geneA", "ACGTACGTAC"),
        ("geneB", "ACGTACGTAG"), // 9/10 identical to geneA
        ("geneC", "TTTTGGGGCC"), // unrelated
    };

    private static HashSet<string> ClusterMembershipKeys(
        IReadOnlyList<(string, string)> genes, double identity) =>
        PanGenomeAnalyzer.ClusterGenes(Genomes(genes), identity)
            .Select(c => string.Join(",", c.GeneIds.OrderBy(g => g)))
            .ToHashSet();

    #region PANGEN-CLUSTER-001 MON — lower identity yields coarser (fewer) clusters

    [Test]
    [Description("MON: a lower identity cutoff lets more genes merge into a representative, so the number of clusters cannot increase as the threshold drops.")]
    public void ClusterGenes_LowerIdentity_Coarser()
    {
        int strict = PanGenomeAnalyzer.ClusterGenes(Genomes(ClusterGenesData), 0.95).Count();
        int loose = PanGenomeAnalyzer.ClusterGenes(Genomes(ClusterGenesData), 0.85).Count();

        loose.Should().BeLessThanOrEqualTo(strict, because: "a lower identity cutoff can only merge clusters, never split them");
        loose.Should().BeLessThan(strict, because: "geneA and geneB (90% identical) merge at 0.85 but not at 0.95");
    }

    #endregion

    #region PANGEN-CLUSTER-001 INV — clustering is independent of gene order

    [Test]
    [Description("INV: the length-sorted greedy clustering produces the same cluster membership regardless of the input gene order.")]
    public void ClusterGenes_GeneOrder_Invariant()
    {
        var reordered = new[] { ClusterGenesData[2], ClusterGenesData[0], ClusterGenesData[1] };

        ClusterMembershipKeys(reordered, 0.85).Should().BeEquivalentTo(ClusterMembershipKeys(ClusterGenesData, 0.85),
            because: "clustering depends on the gene set and the identity rule, not the input order");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: PANGEN-CORE-001 — core genome (PanGenome).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 191.
    //
    // API under test (PanGenomeAnalyzer.ClusterGenes; core = clusters present in every genome):
    //   The core genome is the set of ortholog clusters occurring in all genomes.
    //
    // Relations (derived from the "present in all genomes" definition, NOT from output):
    //   • MON  (more genomes ⇒ ⊆ core): adding a genome can only drop clusters absent from it, so the
    //          core shrinks (antitone in the genome set).
    //   • INV  (genome order independent): the core depends on which clusters span all genomes, not on
    //          the genome iteration order.
    // ───────────────────────────────────────────────────────────────────────────

    // Distinct (pairwise < 95% identical) gene sequences; identical sequences across genomes are orthologs.
    private const string S1 = "AAAACCCCGGGG";
    private const string S2 = "TTTTACGTACGT";
    private const string S3 = "GGGGCCCCAAAA";
    private const string S4 = "CACACACAGTGT";
    private const string S5 = "TGTGTGTGCACA";
    private const string S6 = "ACGTACGTTTTT";

    private static IReadOnlyList<(string, string)> Genome(string id, params string[] seqs) =>
        seqs.Select((s, i) => ($"{id}_g{i}", s)).ToList();

    // Core = consensus sequences of clusters present in every genome.
    private static HashSet<string> CoreConsensus(IReadOnlyDictionary<string, IReadOnlyList<(string, string)>> genomes)
    {
        int n = genomes.Count;
        return PanGenomeAnalyzer.ClusterGenes(genomes, 0.95)
            .Where(c => c.GenomeCount == n)
            .Select(c => c.ConsensusSequence).ToHashSet();
    }

    #region PANGEN-CORE-001 MON — adding a genome shrinks the core

    [Test]
    [Description("MON: a core cluster must occur in every genome, so adding a genome can only remove clusters it lacks — the core of more genomes is a subset of the core of fewer.")]
    public void Core_MoreGenomes_Subset()
    {
        var twoGenomes = new Dictionary<string, IReadOnlyList<(string, string)>>
        {
            ["gA"] = Genome("gA", S1, S2, S3),
            ["gB"] = Genome("gB", S1, S2, S4),
        };
        var threeGenomes = new Dictionary<string, IReadOnlyList<(string, string)>>(twoGenomes)
        {
            ["gC"] = Genome("gC", S1, S5, S6),
        };

        var core2 = CoreConsensus(twoGenomes);
        var core3 = CoreConsensus(threeGenomes);

        core3.IsSubsetOf(core2).Should().BeTrue(because: "a cluster in all three genomes is in the first two");
        core2.Count.Should().BeGreaterThan(core3.Count, because: "S2 leaves the core once genome gC (which lacks it) is added");
    }

    #endregion

    #region PANGEN-CORE-001 INV — core is independent of genome order

    [Test]
    [Description("INV: the core depends on which clusters span all genomes, so building the genome map in a different order yields the same core.")]
    public void Core_GenomeOrder_Invariant()
    {
        var forward = new Dictionary<string, IReadOnlyList<(string, string)>>
        {
            ["gA"] = Genome("gA", S1, S2, S3),
            ["gB"] = Genome("gB", S1, S2, S4),
            ["gC"] = Genome("gC", S1, S2, S5),
        };
        var reordered = new Dictionary<string, IReadOnlyList<(string, string)>>
        {
            ["gC"] = Genome("gC", S1, S2, S5),
            ["gA"] = Genome("gA", S1, S2, S3),
            ["gB"] = Genome("gB", S1, S2, S4),
        };

        CoreConsensus(reordered).Should().BeEquivalentTo(CoreConsensus(forward),
            because: "the core is a function of the genome set, not its iteration order");
    }

    #endregion
}
