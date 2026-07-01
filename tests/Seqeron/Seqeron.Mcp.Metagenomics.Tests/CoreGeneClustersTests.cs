using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using Seqeron.Genomics.Metagenomics;
using Seqeron.Mcp.Metagenomics.Tools;

namespace Seqeron.Mcp.Metagenomics.Tests;

// Wraps PanGenomeAnalyzer.GetCoreGeneClusters (Roary fractional core rule, Page 2015):
// a cluster is core when occupancy / totalGenomes >= threshold.
[TestFixture]
public class CoreGeneClustersTests
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
    public void CoreGeneClusters_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => MetagenomicsTools.CoreGeneClusters(
            new[] { Cluster("c1", 3) }, totalGenomes: 3, threshold: 1.0));

        // Empty cluster list is valid input (no core), not an error.
        Assert.DoesNotThrow(() => MetagenomicsTools.CoreGeneClusters(
            System.Array.Empty<PanGenomeAnalyzer.GeneCluster>(), totalGenomes: 5, threshold: 0.99));
    }

    [Test]
    public void CoreGeneClusters_Binding_InvokesSuccessfully()
    {
        // Occupancies {3,2,1} over 3 genomes at threshold 1.0 -> only the 3/3 cluster is core.
        var result = MetagenomicsTools.CoreGeneClusters(
            new[] { Cluster("c1", 3), Cluster("c2", 2), Cluster("c3", 1) },
            totalGenomes: 3, threshold: 1.0);

        Assert.Multiple(() =>
        {
            Assert.That(result.Items, Has.Count.EqualTo(1),
                "threshold 1.0 over 3 genomes requires occupancy >= 3.");
            Assert.That(result.Items[0].ClusterId, Is.EqualTo("c1"));
        });

        // Fractional boundary at 99% over 100 genomes: 99/100 is core, 98/100 is not.
        var boundary = MetagenomicsTools.CoreGeneClusters(
            new[] { Cluster("c99", 99), Cluster("c98", 98) },
            totalGenomes: 100, threshold: 0.99);

        Assert.Multiple(() =>
        {
            Assert.That(boundary.Items, Has.Count.EqualTo(1),
                "Only occupancy >= 99 (>= 99%) is core (Page 2015).");
            Assert.That(boundary.Items[0].ClusterId, Is.EqualTo("c99"));
        });
    }
}
