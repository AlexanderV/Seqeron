using NUnit.Framework;
using Seqeron.Genomics.Metagenomics;
using Seqeron.Mcp.Metagenomics.Tools;

namespace Seqeron.Mcp.Metagenomics.Tests;

// Wraps PanGenomeAnalyzer.CreateCoreGenomeAlignment.
// Reference values from Seqeron.Genomics.Tests PanGenomeAnalyzerTests.
[TestFixture]
public class CoreGenomeAlignmentTests
{
    private static GenomeInput Genome(string id, params (string GeneId, string Seq)[] genes)
        => new(id, genes.Select(g => new GeneInput(g.GeneId, g.Seq)).ToList());

    private static PanGenomeAnalyzer.GeneCluster Cluster(string id, string geneId) =>
        new(id, new[] { geneId }, new[] { "g1" }, 1, 1.0, geneId);

    private static GenomeInput[] Genomes() => new[]
    {
        Genome("g1", ("gene1", "ATGC"), ("gene2", "GCTA")),
    };

    private static PanGenomeAnalyzer.GeneCluster[] Clusters() => new[]
    {
        Cluster("c1", "gene1"),
        Cluster("c2", "gene2"),
    };

    [Test]
    public void CoreGenomeAlignment_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() =>
            MetagenomicsTools.CoreGenomeAlignment(Genomes(), Clusters(), "g1"));

        // A genome id absent from the input is defined -> empty alignment, not an error.
        Assert.DoesNotThrow(() =>
            MetagenomicsTools.CoreGenomeAlignment(Genomes(), Clusters(), "nonexistent"));
    }

    [Test]
    public void CoreGenomeAlignment_Binding_InvokesSuccessfully()
    {
        // Concatenate g1's gene1 ("ATGC") then gene2 ("GCTA") -> "ATGCGCTA".
        var aligned = MetagenomicsTools.CoreGenomeAlignment(Genomes(), Clusters(), "g1");
        Assert.That(aligned.Result, Is.EqualTo("ATGCGCTA"),
            "Core alignment = concatenation of each core cluster's member sequence in g1.");

        // Nonexistent genome -> empty string.
        var missing = MetagenomicsTools.CoreGenomeAlignment(Genomes(), Clusters(), "nonexistent");
        Assert.That(missing.Result, Is.Empty);
    }
}
