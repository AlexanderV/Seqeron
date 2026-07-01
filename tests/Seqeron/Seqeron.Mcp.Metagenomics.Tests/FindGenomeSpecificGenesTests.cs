using System.Linq;
using NUnit.Framework;
using Seqeron.Genomics.Metagenomics;
using Seqeron.Mcp.Metagenomics.Tools;

namespace Seqeron.Mcp.Metagenomics.Tests;

// Wraps PanGenomeAnalyzer.FindGenomeSpecificGenes: per genome, the singleton (GenomeCount==1)
// cluster ids that only that genome owns. Reference values from
// Seqeron.Genomics.Tests PanGenomeAnalyzerTests.
[TestFixture]
public class FindGenomeSpecificGenesTests
{
    private static GenomeInput Genome(string id, params (string GeneId, string Seq)[] genes)
        => new(id, genes.Select(g => new GeneInput(g.GeneId, g.Seq)).ToList());

    private static GenomeInput[] Genomes() => new[]
    {
        Genome("g1", ("unique1", "AAAA")),
        Genome("g2", ("unique2", "TTTT")),
    };

    private static PanGenomeAnalyzer.GeneCluster[] Clusters() => new[]
    {
        new PanGenomeAnalyzer.GeneCluster("c1", new[] { "unique1" }, new[] { "g1" }, 1, 1.0, "AAAA"),
        new PanGenomeAnalyzer.GeneCluster("c2", new[] { "unique2" }, new[] { "g2" }, 1, 1.0, "TTTT"),
    };

    [Test]
    public void FindGenomeSpecificGenes_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => MetagenomicsTools.FindGenomeSpecificGenes(Genomes(), Clusters()));

        // No clusters -> defined (no genome-specific genes), not an error.
        Assert.DoesNotThrow(() => MetagenomicsTools.FindGenomeSpecificGenes(
            Genomes(), System.Array.Empty<PanGenomeAnalyzer.GeneCluster>()));
    }

    [Test]
    public void FindGenomeSpecificGenes_Binding_InvokesSuccessfully()
    {
        var result = MetagenomicsTools.FindGenomeSpecificGenes(Genomes(), Clusters());

        Assert.Multiple(() =>
        {
            Assert.That(result.Items, Has.Count.EqualTo(2),
                "Each genome owns exactly one singleton cluster.");
            var g1 = result.Items.Single(i => i.GenomeId == "g1");
            var g2 = result.Items.Single(i => i.GenomeId == "g2");
            Assert.That(g1.UniqueGeneIds, Is.EqualTo(new[] { "c1" }));
            Assert.That(g2.UniqueGeneIds, Is.EqualTo(new[] { "c2" }));
        });
    }
}
