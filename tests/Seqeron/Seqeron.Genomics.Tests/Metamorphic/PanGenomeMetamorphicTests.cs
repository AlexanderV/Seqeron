namespace Seqeron.Genomics.Tests.Metamorphic;

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

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: PANGEN-HEAP-001 — Heaps' law fit (PanGenome).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 192.
    //
    // API under test (PanGenomeAnalyzer.FitHeapsLaw):
    //   Fits n(N) = K·N^(−alpha) to the permutation-averaged new-gene-discovery curve (micropan).
    //
    // Relations (derived from the Heaps model, NOT from output):
    //   • INV  (genome order independent): the seeded permutation-averaged fit is reproducible and
    //          its open/closed classification characterises the genome SET, not its iteration order.
    //   • MON  (more genomes ⇒ better-constrained, diminishing-returns model): the fitted curve has
    //          alpha ≥ 0, so the predicted number of new genes is non-increasing in N — the Heaps
    //          diminishing-returns property that improves with each added genome.
    // ───────────────────────────────────────────────────────────────────────────

    // Distinct DNA sequence per integer id (adjacent ids differ in one base ⇒ < 95% identity).
    private static string SeqOf(int id)
    {
        const string a = "ACGT";
        var chars = new char[12];
        for (int i = 0; i < 12; i++) chars[i] = a[(id >> (2 * i)) & 3];
        return new string(chars);
    }

    // Open pan-genome: a shared 2-gene core plus 3 genome-unique genes each.
    private static Dictionary<string, IReadOnlyList<(string, string)>> OpenPanGenome()
    {
        var g = new Dictionary<string, IReadOnlyList<(string, string)>>();
        for (int gi = 0; gi < 5; gi++)
        {
            var genes = new List<(string, string)>
            {
                ($"g{gi}_core1", SeqOf(1000)), // shared core
                ($"g{gi}_core2", SeqOf(1001)),
            };
            for (int u = 0; u < 3; u++)
                genes.Add(($"g{gi}_u{u}", SeqOf(gi * 10 + u + 1))); // genome-unique
            g[$"genome{gi}"] = genes;
        }
        return g;
    }

    #region PANGEN-HEAP-001 INV — the fit is reproducible and order-robust

    [Test]
    [Description("INV: the seeded permutation-averaged fit is reproducible on identical input, and its open/closed classification characterises the genome set regardless of iteration order.")]
    public void Heaps_OrderIndependent_AndReproducible()
    {
        var genomes = OpenPanGenome();

        var fit1 = PanGenomeAnalyzer.FitHeapsLaw(genomes, 0.95);
        var fit2 = PanGenomeAnalyzer.FitHeapsLaw(genomes, 0.95);
        fit2.Alpha.Should().Be(fit1.Alpha, because: "the seeded fit is reproducible on identical input");
        fit2.Intercept.Should().Be(fit1.Intercept, because: "the seeded fit is reproducible on identical input");

        var reordered = new Dictionary<string, IReadOnlyList<(string, string)>>(
            genomes.Reverse().ToDictionary(kv => kv.Key, kv => kv.Value));
        PanGenomeAnalyzer.FitHeapsLaw(reordered, 0.95).IsOpen.Should().Be(fit1.IsOpen,
            because: "the open/closed nature of the pan-genome is a property of the genome set, not its order");
    }

    #endregion

    #region PANGEN-HEAP-001 MON — the fitted model has diminishing returns

    [Test]
    [Description("MON: the fitted Heaps model n(N)=K·N^(−alpha) has alpha ≥ 0, so the predicted number of new genes is non-increasing as more genomes are added.")]
    public void Heaps_PredictedNewGenes_NonIncreasing()
    {
        var fit = PanGenomeAnalyzer.FitHeapsLaw(OpenPanGenome(), 0.95);

        fit.Alpha.Should().BeGreaterThanOrEqualTo(0.0, because: "Heaps' decay exponent is non-negative");
        double previous = double.MaxValue;
        foreach (int n in new[] { 2, 3, 5, 10, 20 })
        {
            double predicted = fit.PredictNewGenes(n);
            predicted.Should().BeLessThanOrEqualTo(previous, because: $"the new-gene curve is non-increasing, so N={n} predicts no more than the previous N");
            previous = predicted;
        }
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: PANGEN-MARKER-001 — phylogenetic marker selection (PanGenome).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 193.
    //
    // API under test (PanGenomeAnalyzer.SelectPhylogeneticMarkers):
    //   Keeps the single-copy core clusters that carry at least one parsimony-informative site.
    //
    // Relations (derived from the filtering definition, NOT from output):
    //   • SUB  (markers ⊆ core): markers are selected from the supplied core clusters, so every
    //          marker is a core cluster (conserved core clusters with no informative site are dropped).
    //   • INV  (genome order independent): the selection depends on cluster content, so reordering the
    //          genomes yields the same marker set.
    // ───────────────────────────────────────────────────────────────────────────

    // A variable single-copy core gene (one parsimony-informative site) and a conserved core gene.
    private static string MarkerVariableGene(int genomeIndex)
    {
        var chars = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".ToCharArray(); // length 40
        chars[20] = genomeIndex < 2 ? 'A' : 'C'; // site split 2/2 across 4 genomes ⇒ parsimony-informative
        return new string(chars);
    }
    private const string MarkerConservedGene = "TTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGG"; // identical everywhere

    private static Dictionary<string, IReadOnlyList<(string, string)>> MarkerGenomes()
    {
        var g = new Dictionary<string, IReadOnlyList<(string, string)>>();
        for (int gi = 0; gi < 4; gi++)
            g[$"genome{gi}"] = new[] { ($"g{gi}_A", MarkerVariableGene(gi)), ($"g{gi}_B", MarkerConservedGene) };
        return g;
    }

    private static HashSet<string> ClusterKeysByGenes(IEnumerable<PanGenomeAnalyzer.GeneCluster> clusters) =>
        clusters.Select(c => string.Join(",", c.GeneIds.OrderBy(x => x))).ToHashSet();

    #region PANGEN-MARKER-001 SUB — markers are a subset of the core

    [Test]
    [Description("SUB: markers are selected from the core clusters, so every marker is a core cluster; the conserved core gene (no informative site) is in the core but excluded from markers.")]
    public void Markers_AreSubsetOfCore()
    {
        var genomes = MarkerGenomes();
        var core = PanGenomeAnalyzer.ClusterGenes(genomes, 0.95).Where(c => c.GenomeCount == 4).ToList();
        var markers = PanGenomeAnalyzer.SelectPhylogeneticMarkers(genomes, core, 4).ToList();

        var coreKeys = ClusterKeysByGenes(core);
        var markerKeys = ClusterKeysByGenes(markers);

        markers.Should().NotBeEmpty(because: "the variable core gene is an informative single-copy marker");
        markerKeys.IsSubsetOf(coreKeys).Should().BeTrue(because: "markers are filtered from the supplied core clusters");
        markerKeys.Count.Should().BeLessThan(coreKeys.Count, because: "the conserved core gene carries no informative site and is excluded");
    }

    #endregion

    #region PANGEN-MARKER-001 INV — marker selection is independent of genome order

    [Test]
    [Description("INV: marker selection depends on cluster content, so reordering the genomes yields the same set of markers (by member genes).")]
    public void Markers_GenomeOrder_Invariant()
    {
        var genomes = MarkerGenomes();
        var reordered = new Dictionary<string, IReadOnlyList<(string, string)>>(
            genomes.Reverse().ToDictionary(kv => kv.Key, kv => kv.Value));

        HashSet<string> Markers(Dictionary<string, IReadOnlyList<(string, string)>> g)
        {
            var core = PanGenomeAnalyzer.ClusterGenes(g, 0.95).Where(c => c.GenomeCount == 4).ToList();
            return ClusterKeysByGenes(PanGenomeAnalyzer.SelectPhylogeneticMarkers(g, core, 4));
        }

        Markers(reordered).Should().BeEquivalentTo(Markers(genomes),
            because: "the marker set is a function of the genome content, not the genome iteration order");
    }

    #endregion
}
