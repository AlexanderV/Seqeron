using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using Seqeron.Genomics.Metagenomics;
using static Seqeron.Genomics.Metagenomics.PanGenomeAnalyzer;

namespace Seqeron.Genomics.Tests.Mutation;

/// <summary>
/// PANGEN-* mutation killers: exact classification tests for the Roary-style core/accessory/unique
/// partition (occupancy/N >= coreFraction; Page et al. 2015) and the CD-HIT greedy clustering with
/// within-cluster average identity (Li &amp; Godzik 2006).
/// </summary>
[TestFixture]
public class PanGenomeAnalyzerMutationTests
{
    private const double Tol = 1e-9;

    private static IReadOnlyDictionary<string, IReadOnlyList<(string GeneId, string Sequence)>> Genomes(
        params (string Genome, (string GeneId, string Sequence)[] Genes)[] entries)
    {
        var d = new Dictionary<string, IReadOnlyList<(string, string)>>();
        foreach (var (g, genes) in entries)
            d[g] = genes.Select(x => (x.GeneId, x.Sequence)).ToList();
        return d;
    }

    [Test]
    public void ConstructPanGenome_GeneInAllGenomes_IsCoreAtInclusiveBoundary()
    {
        // One ortholog cluster present in both genomes; coreFraction = 1.0 ⇒ 2/2 = 1.0 ≥ 1.0 ⇒ core.
        var g = Genomes(
            ("g1", new[] { ("a", "ATGAAAAAA") }),
            ("g2", new[] { ("b", "ATGAAAAAA") }));
        var r = ConstructPanGenome(g, identityThreshold: 0.9, coreFraction: 1.0);

        Assert.That(r.Statistics.CoreGeneCount, Is.EqualTo(1)); // kills occupancy '>' (vs '>=') mutant
        Assert.That(r.Statistics.AccessoryGeneCount, Is.EqualTo(0));
        Assert.That(r.Statistics.UniqueGeneCount, Is.EqualTo(0));
        Assert.That(r.Statistics.CoreFraction, Is.EqualTo(1.0).Within(Tol));
    }

    [Test]
    public void ConstructPanGenome_AccessoryAndUniquePartition()
    {
        // geneA in 2 of 3 genomes (66.7% < 90% ⇒ accessory, not unique); geneB in 1 ⇒ unique.
        var g = Genomes(
            ("g1", new[] { ("a", "ATGAAAAAA") }),
            ("g2", new[] { ("b", "ATGAAAAAA") }),
            ("g3", new[] { ("c", "CCCGGGTTT") }));
        var r = ConstructPanGenome(g, identityThreshold: 0.9, coreFraction: 0.9);

        Assert.That(r.Statistics.CoreGeneCount, Is.EqualTo(0));
        Assert.That(r.Statistics.AccessoryGeneCount, Is.EqualTo(1));
        Assert.That(r.Statistics.UniqueGeneCount, Is.EqualTo(1));
    }

    [Test]
    public void ClusterGenes_TwoSimilarGenes_AverageIdentityIsExact()
    {
        // AAAA vs AAAT: global identity 3/4 = 0.75 ≥ threshold 0.7 ⇒ one 2-member cluster;
        // the within-cluster average over its single pair = 0.75.
        var g = Genomes(
            ("g1", new[] { ("x", "AAAA") }),
            ("g2", new[] { ("y", "AAAT") }));
        var clusters = ClusterGenes(g, identityThreshold: 0.7).ToList();

        Assert.That(clusters, Has.Count.EqualTo(1));
        Assert.That(clusters[0].GenomeCount, Is.EqualTo(2));
        Assert.That(clusters[0].AverageIdentity, Is.EqualTo(0.75).Within(Tol)); // kills sum/pairs & loop bound
    }

    [Test]
    public void ClusterGenes_DissimilarGenes_AreSeparateSingletons()
    {
        // 0% identity < threshold ⇒ two singleton clusters, each self-identity 1.0.
        var g = Genomes(
            ("g1", new[] { ("x", "AAAA") }),
            ("g2", new[] { ("y", "CCCC") }));
        var clusters = ClusterGenes(g, identityThreshold: 0.9).ToList();

        Assert.That(clusters, Has.Count.EqualTo(2));
        Assert.That(clusters.All(c => c.GenomeCount == 1), Is.True);
        Assert.That(clusters.All(c => Math.Abs(c.AverageIdentity - 1.0) < Tol), Is.True);
    }
}
