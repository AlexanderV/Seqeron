using System.Linq;
using NUnit.Framework;
using Seqeron.Genomics.Metagenomics;
using Seqeron.Mcp.Metagenomics.Tools;

namespace Seqeron.Mcp.Metagenomics.Tests;

// Wraps PanGenomeAnalyzer.AnalyzeAccessoryGenes: accessory clusters are those with
// 1 < GenomeCount < totalGenomes; Frequency = GenomeCount / totalGenomes.
// Reference values hand-derived from the algorithm contract (Tettelin 2005; Page 2015).
[TestFixture]
public class AccessoryGenesTests
{
    private static PanGenomeAnalyzer.GeneCluster Cluster(string id, int genomeCount) =>
        new(
            ClusterId: id,
            GeneIds: new[] { "gene_" + id },
            GenomeIds: Enumerable.Range(1, genomeCount).Select(i => $"genome{i}").ToArray(),
            GenomeCount: genomeCount,
            AverageIdentity: 1.0,
            ConsensusSequence: "ATGC");

    [Test]
    public void AccessoryGenes_Schema_ValidatesCorrectly()
    {
        // A well-formed call does not throw.
        Assert.DoesNotThrow(() =>
            MetagenomicsTools.AccessoryGenes(
                new[] { Cluster("c1", 2) },
                totalGenomes: 3));

        // Empty cluster list is valid input (no accessory genes) — no throw.
        Assert.DoesNotThrow(() =>
            MetagenomicsTools.AccessoryGenes(
                System.Array.Empty<PanGenomeAnalyzer.GeneCluster>(),
                totalGenomes: 3));
    }

    [Test]
    public void AccessoryGenes_Binding_InvokesSuccessfully()
    {
        // c1: 3/3 (core, excluded), c2: 2/3 (accessory), c3: 1/3 (unique, excluded).
        var clusters = new[]
        {
            Cluster("c1", 3),
            Cluster("c2", 2),
            Cluster("c3", 1),
        };

        var result = MetagenomicsTools.AccessoryGenes(clusters, totalGenomes: 3);

        Assert.Multiple(() =>
        {
            Assert.That(result.Items, Has.Count.EqualTo(1),
                "Only the 2-of-3 cluster is accessory (1 < GenomeCount < totalGenomes).");
            Assert.That(result.Items[0].ClusterId, Is.EqualTo("c2"));
            Assert.That(result.Items[0].Frequency, Is.EqualTo(2.0 / 3.0).Within(1e-12),
                "Frequency = GenomeCount / totalGenomes = 2/3.");
            Assert.That(result.Items[0].GenomesWithGene, Has.Count.EqualTo(2),
                "The accessory cluster spans 2 genomes.");
        });
    }
}
