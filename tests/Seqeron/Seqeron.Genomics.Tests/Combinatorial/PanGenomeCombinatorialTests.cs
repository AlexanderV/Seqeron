namespace Seqeron.Genomics.Tests.Combinatorial;

using Seqeron.Genomics.Metagenomics;

/// <summary>
/// Combinatorial (pairwise / full-grid) tests for the PanGenome area (PanGenomeAnalyzer,
/// Seqeron.Genomics.Metagenomics).
///
/// See <see cref="CompositionCombinatorialTests"/> for the rationale of combinatorial testing.
/// Each grid cell carries a real business assertion; small grids use the exhaustive
/// <c>[Combinatorial]</c> product (a strict superset of pairwise).
/// — docs/checklists/09_COMBINATORIAL_TESTING.md §Description.
/// </summary>
[TestFixture]
[Category("Combinatorial")]
[Category("PanGenome")]
public class PanGenomeCombinatorialTests
{
    // Pairwise-distinct gene families (≤ 50% cross-identity): identical within a family across genomes,
    // so families never merge at identity thresholds ≥ 0.7 and members of a family always cluster.
    private static readonly string[] Families =
    {
        "ACACACACACAC", "GTGTGTGTGTGT", "AGAGAGAGAGAG", "CTCTCTCTCTCT",
    };

    private static IReadOnlyDictionary<string, IReadOnlyList<(string GeneId, string Sequence)>> Genomes(
        params (string GenomeId, int[] FamilyIndices)[] spec)
    {
        var d = new Dictionary<string, IReadOnlyList<(string, string)>>();
        foreach (var (genomeId, families) in spec)
            d[genomeId] = families.Select(f => ($"{genomeId}_F{f}", Families[f])).ToList<(string, string)>();
        return d;
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: PANGEN-CLUSTER-001 — Ortholog clustering (CD-HIT greedy) (PanGenome)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 190.
    // Spec: tests/TestSpecs/PANGEN-CLUSTER-001.md (canonical ClusterGenes). ADVANCED §10.
    // Dimensions: nGenes(3) × identity(3). Grid 3×3 = 9 (full, exhaustive ⊇ pairwise).
    //
    // Model (Li & Godzik 2006, CD-HIT): greedy incremental clustering — each gene joins the first
    // representative whose global identity (identical residues / shorter length) meets the threshold,
    // else starts a new cluster. Clustering is a partition of the genes.
    //
    // Axis mapping (documented): nGenes → the number of distinct gene families per genome (2/3/4);
    // identity → the CD-HIT identity threshold. Each family is identical across the 3 genomes and ≤ 50%
    // identical to the others. The combinatorial point: at every threshold the families stay separate
    // (one cluster each, three genomes per cluster) and the clustering partitions all genes.
    // ═══════════════════════════════════════════════════════════════════════

    [Test, Combinatorial]
    public void PanGenCluster_FamiliesFormOnePartitionedCluster_AcrossGeneCountAndIdentity(
        [Values(2, 3, 4)] int nGenes,
        [Values(0.7, 0.8, 0.9)] double identity)
    {
        int[] fams = Enumerable.Range(0, nGenes).ToArray();
        var genomes = Genomes(("G1", fams), ("G2", fams), ("G3", fams));

        var clusters = PanGenomeAnalyzer.ClusterGenes(genomes, identity).ToList();

        clusters.Should().HaveCount(nGenes, "each distinct family is one cluster");
        clusters.Should().OnlyContain(c => c.GenomeCount == 3 && c.GeneIds.Count == 3, "a family is shared by all 3 genomes");
        clusters.Should().OnlyContain(c => c.AverageIdentity >= -1e-9 && c.AverageIdentity <= 1.0 + 1e-9);

        // Partition: every gene appears in exactly one cluster.
        var allGenes = clusters.SelectMany(c => c.GeneIds).ToList();
        allGenes.Should().OnlyHaveUniqueItems("clusters partition the genes");
        allGenes.Should().HaveCount(3 * nGenes, "every gene is clustered exactly once");
    }

    /// <summary>
    /// Interaction witness — the identity threshold gates merging: two genes at 60% identity cluster
    /// together at threshold 0.5 but split at 0.7.
    /// </summary>
    [Test]
    public void PanGenCluster_IdentityThresholdGatesMerging()
    {
        // Two length-10 genes differing at 4 of 10 positions ⇒ 60% identity.
        var genomes = new Dictionary<string, IReadOnlyList<(string, string)>>
        {
            ["G1"] = new List<(string, string)> { ("g1", "AAAAAAAAAA") },
            ["G2"] = new List<(string, string)> { ("g2", "AAAAAATTTT") }, // 6/10 identical
        };

        PanGenomeAnalyzer.ClusterGenes(genomes, 0.5).Should().HaveCount(1, "60% identity ≥ 0.5 ⇒ one cluster");
        PanGenomeAnalyzer.ClusterGenes(genomes, 0.7).Should().HaveCount(2, "60% identity < 0.7 ⇒ two clusters");
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: PANGEN-CORE-001 — Core / accessory / unique partition (PanGenome)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 191.
    // Spec: tests/TestSpecs/PANGEN-CORE-001.md (canonical ConstructPanGenome / GetCoreGeneClusters).
    // ADVANCED §10.
    // Dimensions: nGenomes(3) × coreThreshold(3). Grid 3×3 = 9 (full, exhaustive ⊇ pairwise).
    //
    // Model (Page 2015 Roary; Tettelin 2005): a gene cluster is CORE when its occupancy fraction
    // occupancy/N ≥ coreFraction (a fractional test), UNIQUE when present in exactly one genome, else
    // ACCESSORY. The three categories partition the clusters.
    //
    // Engineered construct: a CORE family present in all N genomes, a SHELL family in exactly 2, and a
    // UNIQUE gene in 1. The combinatorial point: across genome count and core threshold, the core set
    // is exactly the clusters whose occupancy fraction clears the threshold, and the core/accessory/
    // unique counts partition all clusters — verified against an independent occupancy recount.
    // ═══════════════════════════════════════════════════════════════════════

    [Test, Combinatorial]
    public void PanGenCore_OccupancyPartition_AcrossGenomeCountAndThreshold(
        [Values(3, 4, 5)] int nGenomes,
        [Values(0.5, 0.75, 1.0)] double coreThreshold)
    {
        // Family 0 (CORE) in every genome; family 1 (SHELL) in genomes 0 and 1; family 2 (UNIQUE) in genome 0.
        var spec = new List<(string, int[])>();
        for (int g = 0; g < nGenomes; g++)
        {
            var fams = new List<int> { 0 };           // core in all
            if (g < 2) fams.Add(1);                   // shell in first two
            if (g == 0) fams.Add(2);                  // unique in the first
            spec.Add(($"G{g}", fams.ToArray()));
        }
        var genomes = Genomes(spec.ToArray());

        var clusters = PanGenomeAnalyzer.ClusterGenes(genomes, 0.9).ToList();
        clusters.Should().HaveCount(3, "CORE, SHELL and UNIQUE families");

        // Independent occupancy ground truth.
        int ExpectedCore() => clusters.Count(c => (double)c.GenomeCount / nGenomes >= coreThreshold - 1e-9);

        PanGenomeAnalyzer.GetCoreGeneClusters(clusters, nGenomes, coreThreshold).Select(c => c.ClusterId)
            .Should().BeEquivalentTo(clusters.Where(c => (double)c.GenomeCount / nGenomes >= coreThreshold - 1e-9).Select(c => c.ClusterId),
                "core = clusters whose occupancy fraction clears the threshold");

        var result = PanGenomeAnalyzer.ConstructPanGenome(genomes, 0.9, coreThreshold);
        result.Statistics.CoreGeneCount.Should().Be(ExpectedCore(), "core count matches the occupancy rule");
        (result.Statistics.CoreGeneCount + result.Statistics.AccessoryGeneCount + result.Statistics.UniqueGeneCount)
            .Should().Be(clusters.Count, "core/accessory/unique partition all clusters");
        result.Statistics.CoreFraction.Should().BeInRange(0.0, 1.0);
        result.Statistics.GenomeFluidity.Should().BeInRange(0.0, 1.0);
    }

    /// <summary>
    /// Interaction witness — a SHELL family present in 2 of N genomes is core at a low threshold but
    /// accessory at a strict one; the always-present CORE family is core regardless.
    /// </summary>
    [Test]
    public void PanGenCore_ShellFamilyCrossesThreshold()
    {
        // 3 genomes: family 0 in all, family 1 in two ⇒ shell occupancy 2/3 ≈ 0.667.
        var genomes = Genomes(("G0", new[] { 0, 1 }), ("G1", new[] { 0, 1 }), ("G2", new[] { 0 }));
        var clusters = PanGenomeAnalyzer.ClusterGenes(genomes, 0.9).ToList();

        PanGenomeAnalyzer.GetCoreGeneClusters(clusters, 3, 0.5).Should().HaveCount(2, "2/3 ≥ 0.5 ⇒ shell is core");
        PanGenomeAnalyzer.GetCoreGeneClusters(clusters, 3, 0.75).Should().HaveCount(1, "2/3 < 0.75 ⇒ only the all-present family is core");
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: PANGEN-MARKER-001 — Phylogenetic-marker selection (PanGenome)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 193.
    // Spec: tests/TestSpecs/PANGEN-MARKER-001.md (canonical SelectPhylogeneticMarkers /
    //       CountParsimonyInformativeSites). ADVANCED §10.
    // Dimensions: nGenomes(3) × nMarkers(3). Grid 3×3 = 9 (full, exhaustive ⊇ pairwise).
    //
    // Model (panX, Ding 2018; Roary): a phylogenetic marker is a SINGLE-COPY CORE cluster (present in
    // every genome with exactly one gene each) carrying ≥ 1 parsimony-informative site; markers are
    // ranked by descending PIS and capped at maxMarkers. A column is parsimony-informative when ≥ 2
    // states each occur in ≥ 2 sequences (Zvelebil 2008) — so ≥ 4 genomes are needed for any PIS.
    //
    // The combinatorial point (nGenomes ∈ {4,5,6} so PIS can be positive; nMarkers = maxMarkers cap):
    // exactly the single-copy core, informative clusters are returned, PIS-ranked and capped at
    // maxMarkers; conserved (PIS 0), paralogous and absent clusters are excluded.
    // ═══════════════════════════════════════════════════════════════════════

    private static PanGenomeAnalyzer.GeneCluster Cluster(string id, IReadOnlyList<string> geneIds, int genomeCount) =>
        new(id, geneIds, geneIds.Select((_, i) => $"genome{i}").ToList(), genomeCount, 1.0, geneIds.Count > 0 ? "X" : "");

    [Test, Combinatorial]
    public void PanGenMarker_SelectsSingleCopyInformativeCore_AcrossGenomeCountAndCap(
        [Values(4, 5, 6)] int nGenomes,
        [Values(1, 2, 3)] int maxMarkers)
    {
        var geneSeq = new Dictionary<string, string>();
        var allGenes = new List<(string, string)>();

        // Build a family's per-genome member sequences and register them; return the gene ids.
        List<string> Family(string fam, Func<int, string> seqOf, int memberCount)
        {
            var ids = new List<string>();
            for (int g = 0; g < memberCount; g++)
            {
                string id = $"{fam}_g{g}";
                ids.Add(id);
                allGenes.Add((id, seqOf(g)));
            }
            return ids;
        }

        int half = nGenomes / 2;
        // famA: two informative columns (PIS = 2). famB: one informative column (PIS = 1).
        var a = Family("A", g => (g < half ? "A" : "T").ToString()[0] + (g < half ? "C" : "G") + new string('A', 8), nGenomes);
        var b = Family("B", g => (g < half ? "A" : "T") + new string('C', 9), nGenomes);
        var conserved = Family("Cons", _ => new string('G', 10), nGenomes);          // PIS 0 → excluded
        var paralog = Family("Par", g => new string('A', 10), nGenomes + 1);          // N+1 genes → not single-copy
        var absent = Family("Abs", g => "AAAAA" + (g < half ? "A" : "T") + "AAAA", nGenomes - 1); // present in N-1 → not core

        var genomes = new Dictionary<string, IReadOnlyList<(string, string)>> { ["all"] = allGenes };

        var clusters = new[]
        {
            Cluster("famA", a, nGenomes),
            Cluster("famB", b, nGenomes),
            Cluster("famConserved", conserved, nGenomes),
            Cluster("famParalog", paralog, nGenomes),       // GeneIds.Count = N+1
            Cluster("famAbsent", absent, nGenomes - 1),     // GenomeCount = N-1
        };

        var markers = PanGenomeAnalyzer.SelectPhylogeneticMarkers(genomes, clusters, nGenomes, maxMarkers)
            .Select(c => c.ClusterId).ToList();

        var expected = new[] { "famA", "famB" }.Take(maxMarkers).ToList(); // PIS 2 ranks above PIS 1
        markers.Should().Equal(expected, "single-copy informative core, PIS-ranked, capped at maxMarkers");
    }

    /// <summary>
    /// Interaction witness — the parsimony-informative-site definition: a split-by-two column is
    /// informative, a singleton variant is not, and a monomorphic column is not.
    /// </summary>
    [Test]
    public void PanGenMarker_ParsimonyInformativeSiteDefinition()
    {
        // col informative: two states each in ≥ 2 sequences.
        PanGenomeAnalyzer.CountParsimonyInformativeSites(new[] { "AA", "AA", "TT", "TT" }).Should().Be(2);
        // singleton variant (T in one row) is not informative.
        PanGenomeAnalyzer.CountParsimonyInformativeSites(new[] { "AAAA", "AAAT", "AAAA", "AAAA" }).Should().Be(0);
        // two sequences can never be informative (a state cannot occur in ≥ 2 rows on both sides).
        PanGenomeAnalyzer.CountParsimonyInformativeSites(new[] { "AT", "TA" }).Should().Be(0);
        // monomorphic alignment.
        PanGenomeAnalyzer.CountParsimonyInformativeSites(new[] { "GGGG", "GGGG", "GGGG", "GGGG" }).Should().Be(0);
    }
}
