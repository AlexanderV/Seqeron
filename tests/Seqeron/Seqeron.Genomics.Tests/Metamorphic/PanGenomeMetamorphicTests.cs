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
}
