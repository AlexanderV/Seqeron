using FsCheck;
using FsCheck.Fluent;
using FsCheck.NUnit;
using Seqeron.Genomics.Metagenomics;

namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Property-based tests for pan-genome analysis (PanGenomeAnalyzer): gene clustering, core/accessory
/// partition, Heaps' law openness, and phylogenetic marker selection.
///
/// Test Units: PANGEN-CLUSTER-001, PANGEN-CORE-001, PANGEN-HEAP-001, PANGEN-MARKER-001
/// </summary>
[TestFixture]
[Category("Property")]
[Category("PanGenome")]
public class PanGenomeProperties
{
    private static string RandDna(Random rng, int len)
    {
        const string bases = "ACGT";
        var c = new char[len];
        for (int i = 0; i < len; i++) c[i] = bases[rng.Next(4)];
        return new string(c);
    }

    #region PANGEN-CLUSTER-001: P: each gene in exactly one cluster; M: lower identity → fewer clusters; D: deterministic

    // ClusterGenes is CD-HIT greedy clustering (Li & Godzik 2006): a gene joins the first
    // representative whose identity meets the threshold, else starts a new cluster.

    /// <summary>Builds 2..3 genomes whose genes are drawn from a small pool of ortholog families.</summary>
    private static Arbitrary<IReadOnlyDictionary<string, IReadOnlyList<(string GeneId, string Sequence)>>> GenomesArbitrary() =>
        Gen.Choose(0, int.MaxValue).Select(seed =>
        {
            var rng = new Random(seed);
            int families = 2 + rng.Next(3);
            var pool = new string[families];
            for (int f = 0; f < families; f++) pool[f] = RandDna(rng, 20);

            int nGenomes = 2 + rng.Next(2);
            var genomes = new Dictionary<string, IReadOnlyList<(string, string)>>();
            int gid = 0;
            for (int g = 0; g < nGenomes; g++)
            {
                int m = 2 + rng.Next(2);
                var genes = new List<(string, string)>();
                for (int i = 0; i < m; i++)
                    genes.Add(($"gene{gid++}", pool[rng.Next(families)]));
                genomes[$"G{g}"] = genes;
            }
            return (IReadOnlyDictionary<string, IReadOnlyList<(string, string)>>)genomes;
        }).ToArbitrary();

    /// <summary>
    /// INV-1 (P): clustering is a partition — every gene appears in exactly one cluster, with no gene
    /// lost or duplicated.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Clustering_PartitionsAllGenes()
    {
        return Prop.ForAll(GenomesArbitrary(), genomes =>
        {
            var clusters = PanGenomeAnalyzer.ClusterGenes(genomes, 0.9).ToList();
            var allGenes = genomes.SelectMany(kv => kv.Value.Select(g => g.GeneId)).ToList();
            var clustered = clusters.SelectMany(c => c.GeneIds).ToList();
            bool ok = clustered.Count == allGenes.Count
                      && clustered.Distinct().Count() == clustered.Count
                      && clustered.ToHashSet().SetEquals(allGenes.ToHashSet());
            return ok.Label($"clustering not a partition: {clustered.Count} vs {allGenes.Count} genes");
        });
    }

    /// <summary>
    /// INV-2 (R): each cluster is non-empty, its genome count equals its distinct genomes, and its
    /// average identity is in [0,1].
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Clustering_ClustersAreWellFormed()
    {
        return Prop.ForAll(GenomesArbitrary(), genomes =>
        {
            var clusters = PanGenomeAnalyzer.ClusterGenes(genomes, 0.9).ToList();
            return clusters.All(c =>
                c.GeneIds.Count >= 1 &&
                c.GenomeCount == c.GenomeIds.Distinct().Count() &&
                c.AverageIdentity is >= 0.0 and <= 1.0)
                .Label("a cluster was malformed");
        });
    }

    /// <summary>
    /// INV-3 (M): a lower identity threshold never produces more clusters (it merges at least as much).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Clustering_LowerThreshold_FewerClusters()
    {
        return Prop.ForAll(GenomesArbitrary(), genomes =>
        {
            int loose = PanGenomeAnalyzer.ClusterGenes(genomes, 0.5).Count();
            int strict = PanGenomeAnalyzer.ClusterGenes(genomes, 0.95).Count();
            return (loose <= strict).Label($"lower threshold gave more clusters ({loose} > {strict})");
        });
    }

    /// <summary>
    /// INV-4 (M, explicit + D): two ~70%-identical genes split at high identity but merge at low
    /// identity; clustering is deterministic.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Clustering_IntermediateIdentity_MergesAtLowThreshold()
    {
        var genomes = new Dictionary<string, IReadOnlyList<(string, string)>>
        {
            ["G"] = new[]
            {
                ("a", "AAAAAAAAAAAAAAAAAAAA"),   // 20 A
                ("b", "AAAAAAAAAAAAAACCCCCC"),   // 14 A + 6 C → 70% identity to a
            }
        };
        int strict = PanGenomeAnalyzer.ClusterGenes(genomes, 0.9).Count();
        int loose = PanGenomeAnalyzer.ClusterGenes(genomes, 0.6).Count();
        int loose2 = PanGenomeAnalyzer.ClusterGenes(genomes, 0.6).Count();

        Assert.Multiple(() =>
        {
            Assert.That(strict, Is.EqualTo(2), "70%-identical genes are separate at 0.9");
            Assert.That(loose, Is.EqualTo(1), "70%-identical genes merge at 0.6");
            Assert.That(loose2, Is.EqualTo(loose), "deterministic");
        });
    }

    #endregion

    #region PANGEN-CORE-001: P: core ⊆ every genome; P: core+accessory+unique = pan; M: more genomes → ≤ core; D: deterministic

    // ConstructPanGenome partitions ortholog clusters into core (≥ coreFraction of genomes),
    // accessory, and unique (single-genome) sets (Tettelin et al. 2005; Roary core rule, Page 2015).

    /// <summary>
    /// INV-1 (P): the core, accessory and unique cluster sets are disjoint and together exhaust the
    /// pan-genome (the total cluster count), with counts matching the statistics.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property PanGenome_CategoriesPartitionThePan()
    {
        return Prop.ForAll(GenomesArbitrary(), genomes =>
        {
            var r = PanGenomeAnalyzer.ConstructPanGenome(genomes, 0.9, 0.99);
            var all = r.CoreGenes.Concat(r.AccessoryGenes).Concat(r.UniqueGenes).ToList();
            bool disjoint = all.Distinct().Count() == all.Count;
            bool exhausts = all.Count == r.Statistics.TotalGenes;
            bool countsMatch = r.CoreGenes.Count == r.Statistics.CoreGeneCount
                               && r.AccessoryGenes.Count == r.Statistics.AccessoryGeneCount
                               && r.UniqueGenes.Count == r.Statistics.UniqueGeneCount;
            return (disjoint && exhausts && countsMatch).Label("core/accessory/unique do not partition the pan-genome");
        });
    }

    /// <summary>
    /// INV-2 (P): with coreFraction = 1.0 every core cluster is present in every genome.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property PanGenome_CoreIsInEveryGenome()
    {
        return Prop.ForAll(GenomesArbitrary(), genomes =>
        {
            var clusterOccupancy = PanGenomeAnalyzer.ClusterGenes(genomes, 0.9)
                .ToDictionary(c => c.ClusterId, c => c.GenomeCount);
            var r = PanGenomeAnalyzer.ConstructPanGenome(genomes, 0.9, 1.0);
            return r.CoreGenes.All(id => clusterOccupancy[id] == genomes.Count)
                .Label("a core gene was not present in every genome");
        });
    }

    /// <summary>
    /// INV-3 (M): adding a genome of entirely novel genes cannot increase the (strict) core size.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property PanGenome_MoreGenomes_CoreDoesNotGrow()
    {
        return Prop.ForAll(GenomesArbitrary(), genomes =>
        {
            int baseCore = PanGenomeAnalyzer.ConstructPanGenome(genomes, 0.9, 1.0).Statistics.CoreGeneCount;

            var extended = new Dictionary<string, IReadOnlyList<(string, string)>>(
                genomes.ToDictionary(kv => kv.Key, kv => kv.Value));
            extended["GX"] = new[] { ("novel_gene", "TTTTGGGGCCCCAAAATTTT") }; // a fresh family
            int extCore = PanGenomeAnalyzer.ConstructPanGenome(extended, 0.9, 1.0).Statistics.CoreGeneCount;

            return (extCore <= baseCore).Label($"core grew when a genome was added ({extCore} > {baseCore})");
        });
    }

    /// <summary>
    /// INV-4 (D): Pan-genome construction is deterministic.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property PanGenome_IsDeterministic()
    {
        return Prop.ForAll(GenomesArbitrary(), genomes =>
        {
            var a = PanGenomeAnalyzer.ConstructPanGenome(genomes, 0.9, 0.99).Statistics;
            var b = PanGenomeAnalyzer.ConstructPanGenome(genomes, 0.9, 0.99).Statistics;
            return (a == b).Label("ConstructPanGenome must be deterministic");
        });
    }

    #endregion

    #region PANGEN-HEAP-001: R: α ∈ [0,2]; P: open ⟺ α < 1; predictor = K·N^(−α)

    // FitHeapsLaw fits n(N) = K·N^(−α) to the new-gene-discovery curve (Tettelin et al. 2008;
    // micropan heaps()). The pan-genome is OPEN when α < 1, CLOSED when α ≥ 1.

    private static Dictionary<string, IReadOnlyList<(string, string)>> BuildGenomes(params string[][] familiesPerGenome)
    {
        var rng = new Random(12345);
        var d = new Dictionary<string, IReadOnlyList<(string, string)>>();
        int gid = 0;
        for (int g = 0; g < familiesPerGenome.Length; g++)
        {
            var genes = familiesPerGenome[g].Select(fam => ($"gene{gid++}", fam)).ToList();
            d[$"G{g}"] = genes;
        }
        return d;
    }

    /// <summary>
    /// INV-1 (R + P): the fit has α ∈ [0,2], intercept ∈ [0,10000], the open flag equals α &lt; 1, and
    /// the predictor reproduces K·N^(−α) and is non-negative.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Heaps_FitIsWellFormed()
    {
        return Prop.ForAll(GenomesArbitrary(), genomes =>
        {
            var fit = PanGenomeAnalyzer.FitHeapsLaw(genomes);
            bool ranges = fit.Alpha is >= 0.0 and <= 2.0 && fit.Intercept is >= 0.0 and <= 10000.0;
            bool openFlag = fit.IsOpen == (fit.Alpha < 1.0);
            double pred = fit.PredictNewGenes(2);
            bool predictor = Math.Abs(pred - fit.Intercept * Math.Pow(2, -fit.Alpha)) < 1e-6 && pred >= 0;
            return (ranges && openFlag && predictor).Label($"α={fit.Alpha}, open={fit.IsOpen}, K={fit.Intercept}");
        });
    }

    /// <summary>
    /// INV-2 (P, positive controls): genomes with entirely unique gene content give an OPEN pan-genome
    /// (α &lt; 1); genomes with identical gene content give a CLOSED one (α ≥ 1).
    /// </summary>
    [Test]
    [Category("Property")]
    public void Heaps_OpenVsClosed_AreClassified()
    {
        var shared = new[] { "AAAACCCCGGGGTTTTACGT", "TTTTGGGGCCCCAAAATGCA", "ACACACACGTGTGTGTACGT" };
        // Closed: every genome carries the same gene families.
        var closed = BuildGenomes(shared, shared, shared, shared);
        // Open: every genome carries its own unique families.
        var open = BuildGenomes(
            new[] { "AAAAAAAAAAAAAAAAAAAA", "CCCCCCCCCCCCCCCCCCCC" },
            new[] { "GGGGGGGGGGGGGGGGGGGG", "TTTTTTTTTTTTTTTTTTTT" },
            new[] { "ACGTACGTACGTACGTACGT", "TGCATGCATGCATGCATGCA" },
            new[] { "AGAGAGAGAGAGAGAGAGAG", "TCTCTCTCTCTCTCTCTCTC" });

        Assert.Multiple(() =>
        {
            Assert.That(PanGenomeAnalyzer.FitHeapsLaw(closed).IsOpen, Is.False, "identical gene content → closed");
            Assert.That(PanGenomeAnalyzer.FitHeapsLaw(open).IsOpen, Is.True, "all-unique gene content → open");
        });
    }

    #endregion

    #region PANGEN-MARKER-001: P: markers ⊆ single-copy core clusters; R: marker count ≤ requested; D: deterministic

    // SelectPhylogeneticMarkers keeps single-copy core clusters (present once in every genome) that
    // carry at least one parsimony-informative site, capped at maxMarkers (panX/Roary marker rule).

    /// <summary>
    /// INV-1 (P + R): every marker is one of the candidate clusters, is single-copy core (present once
    /// per genome), and the number of markers never exceeds the requested cap.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Markers_AreSingleCopyCore_AndCapped()
    {
        return Prop.ForAll(GenomesArbitrary(), genomes =>
        {
            const int maxMarkers = 2;
            int total = genomes.Count;
            var clusters = PanGenomeAnalyzer.ClusterGenes(genomes, 0.9).ToList();
            var clusterIds = clusters.Select(c => c.ClusterId).ToHashSet();
            var markers = PanGenomeAnalyzer.SelectPhylogeneticMarkers(genomes, clusters, total, maxMarkers).ToList();

            bool ok = markers.Count <= maxMarkers
                      && markers.All(m => clusterIds.Contains(m.ClusterId)
                                          && m.GenomeCount == total
                                          && m.GeneIds.Count == total);
            return ok.Label($"markers not single-copy core or exceeded cap ({markers.Count})");
        });
    }

    /// <summary>
    /// INV-2 (P, positive control): a single-copy core family with a parsimony-informative column is
    /// selected as a marker; markers are a subset of the core clusters.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Markers_VariableSingleCopyCore_IsSelected()
    {
        // Four genomes; each has one near-identical marker gene. Position 0 reads A,A,C,C across the
        // four genomes → a parsimony-informative site; the rest are conserved (clusters at 0.9).
        var genomes = new Dictionary<string, IReadOnlyList<(string, string)>>
        {
            ["G0"] = new[] { ("m0", "AAAAAAAAAAAAAAAAAAAA") },
            ["G1"] = new[] { ("m1", "AAAAAAAAAAAAAAAAAAAA") },
            ["G2"] = new[] { ("m2", "CAAAAAAAAAAAAAAAAAAA") },
            ["G3"] = new[] { ("m3", "CAAAAAAAAAAAAAAAAAAA") },
        };
        var clusters = PanGenomeAnalyzer.ClusterGenes(genomes, 0.9).ToList();
        var markers = PanGenomeAnalyzer.SelectPhylogeneticMarkers(genomes, clusters, 4, 100).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(markers, Has.Count.EqualTo(1), "the variable single-copy core family is a marker");
            Assert.That(markers[0].GenomeCount, Is.EqualTo(4));
            Assert.That(markers[0].GeneIds, Has.Count.EqualTo(4));
            Assert.That(clusters.Select(c => c.ClusterId), Does.Contain(markers[0].ClusterId), "marker ⊆ clusters");
        });
    }

    /// <summary>
    /// INV-3 (D): Marker selection is deterministic.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Markers_AreDeterministic()
    {
        return Prop.ForAll(GenomesArbitrary(), genomes =>
        {
            var clusters = PanGenomeAnalyzer.ClusterGenes(genomes, 0.9).ToList();
            var a = PanGenomeAnalyzer.SelectPhylogeneticMarkers(genomes, clusters, genomes.Count, 100).Select(m => m.ClusterId).ToList();
            var b = PanGenomeAnalyzer.SelectPhylogeneticMarkers(genomes, clusters, genomes.Count, 100).Select(m => m.ClusterId).ToList();
            return a.SequenceEqual(b).Label("SelectPhylogeneticMarkers must be deterministic");
        });
    }

    #endregion
}
