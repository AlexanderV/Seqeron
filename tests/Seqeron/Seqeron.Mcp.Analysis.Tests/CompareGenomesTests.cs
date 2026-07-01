using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>compare_genomes</c> MCP tool.
/// Vectors and expected core/dispensable partition taken from COMPGEN-COMPARE-001 /
/// the algorithm's own unit tests (ComparativeGenomics_CompareGenomes_Tests,
/// Tettelin 2005 core/dispensable), NOT the wrapper's output.
/// </summary>
[TestFixture]
public class CompareGenomesTests
{
    // Five distinct >=60-nt shared sequences and mutually-dissimilar unique sequences,
    // copied verbatim from the algorithm's unit-test fixture.
    private const string Shared0 = "ATGGCAAAGCTTGATCCGTACGGGTTAACCGGATCAGGTTCAAAGCTTGATCCGTACGGG";
    private const string Unique1 = "CCCCCCGGGGGGTTTTTTAAAAAACCCCCCGGGGGGTTTTTTAAAAAACCCCCCGGGGGG";
    private const string Unique2 = "TTGGAACCTTGGAACCTTGGAACCTTGGAACCTTGGAACCTTGGAACCTTGGAACCTTGG";
    private const string Unique3 = "GATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATT";
    private const string Unique4 = "ACGTAAACCCGGGTTTACGTAAACCCGGGTTTACGTAAACCCGGGTTTACGTAAACCCGG";

    private static GeneInput[] GenomeOf(string prefix, string genomeId, params string[] seqs)
    {
        var genes = new GeneInput[seqs.Length];
        for (int i = 0; i < seqs.Length; i++)
            genes[i] = new GeneInput($"{prefix}{i}", genomeId, i * 100, i * 100 + 60, '+', seqs[i]);
        return genes;
    }

    [Test]
    public void CompareGenomes_Schema_ValidatesCorrectly()
    {
        var g1 = GenomeOf("a", "G1", Shared0, Unique1);
        var g2 = GenomeOf("c", "G2", Shared0, Unique2);
        Assert.DoesNotThrow(() => AnalysisTools.CompareGenomes(g1, g2));
        Assert.Throws<ArgumentException>(() => AnalysisTools.CompareGenomes(Array.Empty<GeneInput>(), g2));
        Assert.Throws<ArgumentException>(() => AnalysisTools.CompareGenomes(g1, null!));
    }

    [Test]
    public void CompareGenomes_OneSharedOneUnique_PartitionsCoreAndSpecific()
    {
        var g1 = GenomeOf("a", "G1", Shared0, Unique1);
        var g2 = GenomeOf("c", "G2", Shared0, Unique2);

        var r = AnalysisTools.CompareGenomes(g1, g2);
        Assert.Multiple(() =>
        {
            Assert.That(r.ConservedGenes, Is.EqualTo(1));
            Assert.That(r.Orthologs, Has.Length.EqualTo(1));
            Assert.That(r.GenomeSpecificGenes1, Is.EqualTo(1));
            Assert.That(r.GenomeSpecificGenes2, Is.EqualTo(1));
        });
    }

    [Test]
    public void CompareGenomes_DisjointContent_NoCoreAllSpecific()
    {
        var g1 = GenomeOf("a", "G1", Unique1, Unique3);
        var g2 = GenomeOf("c", "G2", Unique2, Unique4);

        var r = AnalysisTools.CompareGenomes(g1, g2);
        Assert.Multiple(() =>
        {
            Assert.That(r.ConservedGenes, Is.EqualTo(0));
            Assert.That(r.Orthologs, Is.Empty);
            Assert.That(r.GenomeSpecificGenes1, Is.EqualTo(2));
            Assert.That(r.GenomeSpecificGenes2, Is.EqualTo(2));
        });
    }
}
