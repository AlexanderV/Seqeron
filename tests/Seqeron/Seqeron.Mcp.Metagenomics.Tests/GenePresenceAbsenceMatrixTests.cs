using System.Linq;
using NUnit.Framework;
using Seqeron.Genomics.Metagenomics;
using Seqeron.Mcp.Metagenomics.Tools;

namespace Seqeron.Mcp.Metagenomics.Tests;

// Wraps PanGenomeAnalyzer.CreatePresenceAbsenceMatrix: one row per genome; a cluster is
// present in a genome when the genome owns any of the cluster's gene ids. Reference values
// hand-derived from the algorithm contract (mirrors the Genomics ClusterGenes S5 case).
[TestFixture]
public class GenePresenceAbsenceMatrixTests
{
    private static GenomeInput Genome(string id, params string[] geneIds)
        => new(id, geneIds.Select(g => new GeneInput(g, "ATGC")).ToList());

    private static GenomeInput[] Genomes() => new[]
    {
        Genome("g1", "a", "b"),
        Genome("g2", "c"),
    };

    // Cluster c1 groups genes a (g1) and c (g2); cluster c2 holds only b (g1).
    private static PanGenomeAnalyzer.GeneCluster[] Clusters() => new[]
    {
        new PanGenomeAnalyzer.GeneCluster("c1", new[] { "a", "c" }, new[] { "g1", "g2" }, 2, 1.0, "ATGC"),
        new PanGenomeAnalyzer.GeneCluster("c2", new[] { "b" }, new[] { "g1" }, 1, 1.0, "ATGC"),
    };

    [Test]
    public void GenePresenceAbsenceMatrix_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => MetagenomicsTools.GenePresenceAbsenceMatrix(Genomes(), Clusters()));

        // No clusters -> defined (rows with zero present genes), not an error.
        Assert.DoesNotThrow(() => MetagenomicsTools.GenePresenceAbsenceMatrix(
            Genomes(), System.Array.Empty<PanGenomeAnalyzer.GeneCluster>()));
    }

    [Test]
    public void GenePresenceAbsenceMatrix_Binding_InvokesSuccessfully()
    {
        var result = MetagenomicsTools.GenePresenceAbsenceMatrix(Genomes(), Clusters());

        var g1 = result.Items.Single(r => r.GenomeId == "g1");
        var g2 = result.Items.Single(r => r.GenomeId == "g2");

        Assert.Multiple(() =>
        {
            Assert.That(result.Items, Has.Count.EqualTo(2), "One row per genome.");
            Assert.That(g1.TotalGenes, Is.EqualTo(2), "Two clusters (columns).");
            Assert.That(g1.PresentGenes, Is.EqualTo(2),
                "g1 owns gene a (cluster c1) and gene b (cluster c2).");
            Assert.That(g2.PresentGenes, Is.EqualTo(1),
                "g2 owns gene c (cluster c1) only.");
            // Per-cluster flags for g2.
            Assert.That(g2.GenePresence.Single(p => p.ClusterId == "c1").Present, Is.True);
            Assert.That(g2.GenePresence.Single(p => p.ClusterId == "c2").Present, Is.False);
        });
    }
}
