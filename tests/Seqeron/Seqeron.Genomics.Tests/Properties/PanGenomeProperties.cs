using FsCheck;
using FsCheck.Fluent;
using FsCheck.NUnit;
using Seqeron.Genomics.Metagenomics;

namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Property-based tests for pan-genome analysis (PanGenomeAnalyzer): gene clustering, core/accessory
/// partition, Heaps' law openness, and phylogenetic marker selection.
///
/// Test Units: PANGEN-CLUSTER-001
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
}
