using NUnit.Framework;
using Seqeron.Genomics;
using System.Collections.Generic;
using System.Linq;

namespace Seqeron.Genomics.Tests;

[TestFixture]
public class PanGenomeAnalyzerTests
{
    #region Pan-Genome Construction Tests

    [Test]
    public void ConstructPanGenome_ThreeGenomes_IdentifiesCoreAccessoryUnique()
    {
        var genomes = new Dictionary<string, IReadOnlyList<(string GeneId, string Sequence)>>
        {
            ["genome1"] = new List<(string, string)>
            {
                ("geneA", "ATGCGATCGATCGATCG"),
                ("geneB", "GCTAGCTAGCTAGCTAG"),
                ("geneC", "TACGTACGTACGTACGT")
            },
            ["genome2"] = new List<(string, string)>
            {
                ("geneA", "ATGCGATCGATCGATCG"),
                ("geneB", "GCTAGCTAGCTAGCTAG"),
                ("geneD", "AAAAAAAAAAAAAAAAAA")
            },
            ["genome3"] = new List<(string, string)>
            {
                ("geneA", "ATGCGATCGATCGATCG"),
                ("geneB", "GCTAGCTAGCTAGCTAG"),
                ("geneE", "TTTTTTTTTTTTTTTTTT")
            }
        };

        var result = PanGenomeAnalyzer.ConstructPanGenome(genomes, coreFraction: 0.99);

        Assert.That(result.Statistics.TotalGenomes, Is.EqualTo(3));
        Assert.That(result.Statistics.CoreGeneCount, Is.GreaterThan(0));
        Assert.That(result.Statistics.TotalGenes, Is.GreaterThan(0));
    }

    [Test]
    public void ConstructPanGenome_EmptyInput_ReturnsEmptyResult()
    {
        var genomes = new Dictionary<string, IReadOnlyList<(string GeneId, string Sequence)>>();
        var result = PanGenomeAnalyzer.ConstructPanGenome(genomes);

        Assert.That(result.CoreGenes.Count, Is.EqualTo(0));
        Assert.That(result.Statistics.TotalGenomes, Is.EqualTo(0));
    }

    [Test]
    public void ConstructPanGenome_SingleGenome_AllUnique()
    {
        var genomes = new Dictionary<string, IReadOnlyList<(string GeneId, string Sequence)>>
        {
            ["genome1"] = new List<(string, string)>
            {
                ("gene1", "ATGCGATCG"),
                ("gene2", "GCTAGCTAG")
            }
        };

        var result = PanGenomeAnalyzer.ConstructPanGenome(genomes);

        Assert.That(result.Statistics.TotalGenomes, Is.EqualTo(1));
        Assert.That(result.Statistics.TotalGenes, Is.GreaterThan(0));
    }

    [Test]
    public void ConstructPanGenome_CalculatesCoreFraction()
    {
        var genomes = new Dictionary<string, IReadOnlyList<(string GeneId, string Sequence)>>
        {
            ["g1"] = new List<(string, string)> { ("a", "ATGC"), ("b", "GCTA") },
            ["g2"] = new List<(string, string)> { ("a", "ATGC"), ("c", "TTTT") }
        };

        var result = PanGenomeAnalyzer.ConstructPanGenome(genomes, coreFraction: 0.5);

        Assert.That(result.Statistics.CoreFraction, Is.GreaterThanOrEqualTo(0).And.LessThanOrEqualTo(1));
    }

    [Test]
    public void ConstructPanGenome_IdentifiesAccessoryGenes()
    {
        var genomes = new Dictionary<string, IReadOnlyList<(string GeneId, string Sequence)>>
        {
            ["g1"] = new List<(string, string)> { ("core", "ATGC"), ("accessory", "GCTA") },
            ["g2"] = new List<(string, string)> { ("core", "ATGC"), ("accessory", "GCTA") },
            ["g3"] = new List<(string, string)> { ("core", "ATGC"), ("unique3", "AAAA") }
        };

        var result = PanGenomeAnalyzer.ConstructPanGenome(genomes, coreFraction: 0.99);

        Assert.That(result.Statistics.AccessoryGeneCount + result.Statistics.UniqueGeneCount, Is.GreaterThan(0));
    }

    #endregion

    #region Gene Clustering Tests

    [Test]
    public void ClusterGenes_SimilarSequences_ClustersTogether()
    {
        var genomes = new Dictionary<string, IReadOnlyList<(string GeneId, string Sequence)>>
        {
            ["g1"] = new List<(string, string)> { ("gene1_g1", "ATGCGATCGATCGATCGATCGATCGATCG") },
            ["g2"] = new List<(string, string)> { ("gene1_g2", "ATGCGATCGATCGATCGATCGATCGATCG") } // Identical
        };

        var clusters = PanGenomeAnalyzer.ClusterGenes(genomes, identityThreshold: 0.9).ToList();

        // Identical sequences should cluster together
        Assert.That(clusters.Any(c => c.GeneIds.Count == 2), Is.True);
    }

    [Test]
    public void ClusterGenes_DifferentSequences_SeparateClusters()
    {
        var genomes = new Dictionary<string, IReadOnlyList<(string GeneId, string Sequence)>>
        {
            ["g1"] = new List<(string, string)> { ("gene1", "AAAAAAAAAAAAAAAAAAAAAA") },
            ["g2"] = new List<(string, string)> { ("gene2", "TTTTTTTTTTTTTTTTTTTTTT") }
        };

        var clusters = PanGenomeAnalyzer.ClusterGenes(genomes, identityThreshold: 0.9).ToList();

        Assert.That(clusters.Count, Is.EqualTo(2));
    }

    [Test]
    public void ClusterGenes_EmptyGenomes_ReturnsEmpty()
    {
        var genomes = new Dictionary<string, IReadOnlyList<(string GeneId, string Sequence)>>();
        var clusters = PanGenomeAnalyzer.ClusterGenes(genomes).ToList();

        Assert.That(clusters, Is.Empty);
    }

    [Test]
    public void ClusterGenes_CalculatesAverageIdentity()
    {
        // Two identical genes from different genomes MUST cluster together
        var genomes = new Dictionary<string, IReadOnlyList<(string GeneId, string Sequence)>>
        {
            ["g1"] = new List<(string, string)> { ("gene1", "ATGCGATCGATCGATCGATCGATCGATCG") },
            ["g2"] = new List<(string, string)> { ("gene2", "ATGCGATCGATCGATCGATCGATCGATCG") }
        };

        var clusters = PanGenomeAnalyzer.ClusterGenes(genomes).ToList();

        // Two identical 30bp genes MUST be clustered together
        Assert.That(clusters.Any(c => c.GeneIds.Count > 1), Is.True,
            "Two identical genes MUST form a multi-gene cluster");

        var multiGeneCluster = clusters.First(c => c.GeneIds.Count > 1);
        Assert.That(multiGeneCluster.AverageIdentity, Is.GreaterThan(0.9),
            "Identical sequences should have >90% average identity");
    }

    [Test]
    public void ClusterGenes_RecordsGenomeCount()
    {
        var genomes = new Dictionary<string, IReadOnlyList<(string GeneId, string Sequence)>>
        {
            ["g1"] = new List<(string, string)> { ("gene1", "ATGCGATCGATCGATCGATCGATCGATCG") },
            ["g2"] = new List<(string, string)> { ("gene2", "ATGCGATCGATCGATCGATCGATCGATCG") },
            ["g3"] = new List<(string, string)> { ("gene3", "ATGCGATCGATCGATCGATCGATCGATCG") }
        };

        var clusters = PanGenomeAnalyzer.ClusterGenes(genomes).ToList();

        // Should have at least one cluster with 3 genomes
        Assert.That(clusters.Any(c => c.GenomeCount == 3), Is.True);
    }

    #endregion

    #region Presence/Absence Matrix Tests

    [Test]
    public void CreatePresenceAbsenceMatrix_CreatesCorrectMatrix()
    {
        var genomes = new Dictionary<string, IReadOnlyList<(string GeneId, string Sequence)>>
        {
            ["g1"] = new List<(string, string)> { ("gene1", "ATGC"), ("gene2", "GCTA") },
            ["g2"] = new List<(string, string)> { ("gene1", "ATGC") }
        };

        var clusters = PanGenomeAnalyzer.ClusterGenes(genomes).ToList();
        var matrix = PanGenomeAnalyzer.CreatePresenceAbsenceMatrix(genomes, clusters).ToList();

        Assert.That(matrix.Count, Is.EqualTo(2));
        Assert.That(matrix.First(r => r.GenomeId == "g1").PresentGenes, Is.EqualTo(2));
        Assert.That(matrix.First(r => r.GenomeId == "g2").PresentGenes, Is.EqualTo(1));
    }

    [Test]
    public void CreatePresenceAbsenceMatrix_RecordsTotalGenes()
    {
        var genomes = new Dictionary<string, IReadOnlyList<(string GeneId, string Sequence)>>
        {
            ["g1"] = new List<(string, string)> { ("a", "AAAA"), ("b", "TTTT"), ("c", "GGGG") },
            ["g2"] = new List<(string, string)> { ("a", "AAAA") }
        };

        var clusters = PanGenomeAnalyzer.ClusterGenes(genomes).ToList();
        var matrix = PanGenomeAnalyzer.CreatePresenceAbsenceMatrix(genomes, clusters).ToList();

        Assert.That(matrix[0].TotalGenes, Is.EqualTo(clusters.Count));
    }

    #endregion

    #region Heaps Law Tests

    [Test]
    public void FitHeapsLaw_WithMultipleGenomes_ReturnsValidFit()
    {
        var genomes = new Dictionary<string, IReadOnlyList<(string GeneId, string Sequence)>>
        {
            ["g1"] = new List<(string, string)> { ("a1", "ATGC"), ("b1", "GCTA"), ("c1", "TTTT") },
            ["g2"] = new List<(string, string)> { ("a2", "ATGC"), ("b2", "GCTA"), ("d2", "AAAA") },
            ["g3"] = new List<(string, string)> { ("a3", "ATGC"), ("b3", "GCTA"), ("e3", "GGGG") },
            ["g4"] = new List<(string, string)> { ("a4", "ATGC"), ("f4", "CCCC"), ("g4", "TATA") }
        };

        var fit = PanGenomeAnalyzer.FitHeapsLaw(genomes, permutations: 5);

        Assert.That(fit.K, Is.GreaterThan(0));
        Assert.That(fit.PredictPanGenomeSize, Is.Not.Null);
    }

    [Test]
    public void FitHeapsLaw_PredictorWorks()
    {
        var genomes = new Dictionary<string, IReadOnlyList<(string GeneId, string Sequence)>>
        {
            ["g1"] = new List<(string, string)> { ("a", "ATGC"), ("b", "GCTA") },
            ["g2"] = new List<(string, string)> { ("c", "TTTT"), ("d", "AAAA") },
            ["g3"] = new List<(string, string)> { ("e", "GGGG"), ("f", "CCCC") }
        };

        var fit = PanGenomeAnalyzer.FitHeapsLaw(genomes, permutations: 3);

        double predicted5 = fit.PredictPanGenomeSize(5);
        double predicted10 = fit.PredictPanGenomeSize(10);

        Assert.That(predicted10, Is.GreaterThanOrEqualTo(predicted5));
    }

    [Test]
    public void FitHeapsLaw_TooFewGenomes_ReturnsEmpty()
    {
        var genomes = new Dictionary<string, IReadOnlyList<(string GeneId, string Sequence)>>
        {
            ["g1"] = new List<(string, string)> { ("a", "ATGC") }
        };

        var fit = PanGenomeAnalyzer.FitHeapsLaw(genomes);

        Assert.That(fit.K, Is.EqualTo(0));
    }

    #endregion

    #region Core Genome Tests

    [Test]
    public void GetCoreGeneClusters_FiltersByThreshold()
    {
        var clusters = new List<PanGenomeAnalyzer.GeneCluster>
        {
            new("c1", new[] { "g1" }, new[] { "genome1", "genome2", "genome3" }, 3, 0.95, "ATGC"),
            new("c2", new[] { "g2" }, new[] { "genome1", "genome2" }, 2, 0.95, "GCTA"),
            new("c3", new[] { "g3" }, new[] { "genome1" }, 1, 1.0, "TTTT")
        };

        // threshold=1.0 means all genomes required, so only c1 (3/3) qualifies
        var core = PanGenomeAnalyzer.GetCoreGeneClusters(clusters, totalGenomes: 3, threshold: 1.0).ToList();

        Assert.That(core.Count, Is.EqualTo(1));
        Assert.That(core[0].ClusterId, Is.EqualTo("c1"));
    }

    [Test]
    public void CreateCoreGenomeAlignment_ConcatenatesGenes()
    {
        var genomes = new Dictionary<string, IReadOnlyList<(string GeneId, string Sequence)>>
        {
            ["g1"] = new List<(string, string)> { ("gene1", "ATGC"), ("gene2", "GCTA") }
        };

        var clusters = new List<PanGenomeAnalyzer.GeneCluster>
        {
            new("c1", new[] { "gene1" }, new[] { "g1" }, 1, 1.0, "ATGC"),
            new("c2", new[] { "gene2" }, new[] { "g1" }, 1, 1.0, "GCTA")
        };

        string alignment = PanGenomeAnalyzer.CreateCoreGenomeAlignment(genomes, clusters, "g1");

        Assert.That(alignment, Is.EqualTo("ATGCGCTA"));
    }

    [Test]
    public void CreateCoreGenomeAlignment_NonexistentGenome_ReturnsEmpty()
    {
        var genomes = new Dictionary<string, IReadOnlyList<(string GeneId, string Sequence)>>
        {
            ["g1"] = new List<(string, string)> { ("gene1", "ATGC") }
        };

        var clusters = new List<PanGenomeAnalyzer.GeneCluster>
        {
            new("c1", new[] { "gene1" }, new[] { "g1" }, 1, 1.0, "ATGC")
        };

        string alignment = PanGenomeAnalyzer.CreateCoreGenomeAlignment(genomes, clusters, "nonexistent");

        Assert.That(alignment, Is.Empty);
    }

    #endregion

    #region Accessory Genome Tests

    [Test]
    public void AnalyzeAccessoryGenes_FiltersCorrectly()
    {
        var clusters = new List<PanGenomeAnalyzer.GeneCluster>
        {
            new("core", new[] { "g1" }, new[] { "genome1", "genome2", "genome3" }, 3, 0.95, "ATGC"),
            new("accessory", new[] { "g2" }, new[] { "genome1", "genome2" }, 2, 0.95, "GCTA"),
            new("unique", new[] { "g3" }, new[] { "genome1" }, 1, 1.0, "TTTT")
        };

        var accessory = PanGenomeAnalyzer.AnalyzeAccessoryGenes(clusters, totalGenomes: 3).ToList();

        Assert.That(accessory.Count, Is.EqualTo(1));
        Assert.That(accessory[0].ClusterId, Is.EqualTo("accessory"));
    }

    [Test]
    public void AnalyzeAccessoryGenes_CalculatesFrequency()
    {
        var clusters = new List<PanGenomeAnalyzer.GeneCluster>
        {
            new("acc", new[] { "g1" }, new[] { "g1", "g2" }, 2, 0.95, "ATGC")
        };

        var accessory = PanGenomeAnalyzer.AnalyzeAccessoryGenes(clusters, totalGenomes: 4).ToList();

        Assert.That(accessory[0].Frequency, Is.EqualTo(0.5).Within(0.01));
    }

    [Test]
    public void FindGenomeSpecificGenes_FindsUnique()
    {
        var genomes = new Dictionary<string, IReadOnlyList<(string GeneId, string Sequence)>>
        {
            ["g1"] = new List<(string, string)> { ("unique1", "AAAA") },
            ["g2"] = new List<(string, string)> { ("unique2", "TTTT") }
        };

        var clusters = new List<PanGenomeAnalyzer.GeneCluster>
        {
            new("c1", new[] { "unique1" }, new[] { "g1" }, 1, 1.0, "AAAA"),
            new("c2", new[] { "unique2" }, new[] { "g2" }, 1, 1.0, "TTTT")
        };

        var unique = PanGenomeAnalyzer.FindGenomeSpecificGenes(genomes, clusters).ToList();

        Assert.That(unique.Count, Is.EqualTo(2));
    }

    #endregion

    #region Phylogenetic Marker Tests

    [Test]
    public void SelectPhylogeneticMarkers_FiltersAndLimits()
    {
        var clusters = new List<PanGenomeAnalyzer.GeneCluster>
        {
            new("c1", new[] { "g1" }, new[] { "g1", "g2", "g3" }, 3, 0.85, "ATGCATGCATGC"),
            new("c2", new[] { "g2" }, new[] { "g1", "g2", "g3" }, 3, 0.95, "ATGC"),
            new("c3", new[] { "g3" }, new[] { "g1", "g2", "g3" }, 3, 0.99, "ATGCAT"),
            new("c4", new[] { "g4" }, new[] { "g1", "g2", "g3" }, 3, 0.65, "ATGCATGCATGCAT") // Too divergent
        };

        var markers = PanGenomeAnalyzer.SelectPhylogeneticMarkers(
            clusters, maxMarkers: 2, minIdentity: 0.7, maxIdentity: 0.99).ToList();

        Assert.That(markers.Count, Is.LessThanOrEqualTo(2));
        Assert.That(markers.All(m => m.AverageIdentity >= 0.7 && m.AverageIdentity <= 0.99), Is.True);
    }

    [Test]
    public void SelectPhylogeneticMarkers_PrefersLongerSequences()
    {
        var clusters = new List<PanGenomeAnalyzer.GeneCluster>
        {
            new("short", new[] { "g1" }, new[] { "g1" }, 1, 0.9, "ATGC"),
            new("long", new[] { "g2" }, new[] { "g1" }, 1, 0.9, "ATGCATGCATGCATGC")
        };

        var markers = PanGenomeAnalyzer.SelectPhylogeneticMarkers(clusters, maxMarkers: 1).ToList();

        Assert.That(markers, Has.Count.EqualTo(1), "Should return exactly 1 marker");
        Assert.That(markers[0].ClusterId, Is.EqualTo("long"),
            "Should prefer longer sequence (16bp) over shorter (4bp)");
    }

    #endregion

    #region Pan-Genome Type Tests

    [Test]
    public void ConstructPanGenome_DeterminesPanGenomeType()
    {
        var genomes = new Dictionary<string, IReadOnlyList<(string GeneId, string Sequence)>>
        {
            ["g1"] = new List<(string, string)> { ("core", "ATGC"), ("u1", "AAAA") },
            ["g2"] = new List<(string, string)> { ("core", "ATGC"), ("u2", "TTTT") },
            ["g3"] = new List<(string, string)> { ("core", "ATGC"), ("u3", "GGGG") }
        };

        var result = PanGenomeAnalyzer.ConstructPanGenome(genomes);

        Assert.That(result.Statistics.Type, Is.EqualTo(PanGenomeAnalyzer.PanGenomeType.Open)
            .Or.EqualTo(PanGenomeAnalyzer.PanGenomeType.Closed));
    }

    [Test]
    public void ConstructPanGenome_CalculatesGenomeFluidity()
    {
        var genomes = new Dictionary<string, IReadOnlyList<(string GeneId, string Sequence)>>
        {
            ["g1"] = new List<(string, string)> { ("a", "ATGC"), ("b", "GCTA") },
            ["g2"] = new List<(string, string)> { ("a", "ATGC"), ("c", "TTTT") }
        };

        var result = PanGenomeAnalyzer.ConstructPanGenome(genomes);

        Assert.That(result.Statistics.GenomeFluidity, Is.GreaterThanOrEqualTo(0));
    }

    #endregion
}
