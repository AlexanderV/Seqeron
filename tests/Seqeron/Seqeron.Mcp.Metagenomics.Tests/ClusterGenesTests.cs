using NUnit.Framework;
using Seqeron.Mcp.Metagenomics.Tools;

namespace Seqeron.Mcp.Metagenomics.Tests;

// Wraps PanGenomeAnalyzer.ClusterGenes (CD-HIT greedy global-identity clustering).
// Reference values from Seqeron.Genomics.Tests PanGenomeAnalyzer_ClusterGenes_Tests
// (Li & Godzik 2006, CD-HIT).
[TestFixture]
public class ClusterGenesTests
{
    private const string S = "ATGCATGC";        // 8 bp
    private const string S1Sub = "ATGCATGG";     // 7/8 = 0.875 vs S
    private const string SLong = "ATGCATGCAAAA"; // 12 bp; 8/8 over shorter(8) = 1.0 vs S

    private static GenomeInput Genome(string id, params (string GeneId, string Seq)[] genes)
        => new(id, genes.Select(g => new GeneInput(g.GeneId, g.Seq)).ToList());

    [Test]
    public void ClusterGenes_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => MetagenomicsTools.ClusterGenes(
            new[] { Genome("g1", ("a", S)) }, identityThreshold: 0.9));

        // Empty genome set is valid input (no clusters), not an error.
        Assert.DoesNotThrow(() => MetagenomicsTools.ClusterGenes(
            System.Array.Empty<GenomeInput>(), identityThreshold: 0.9));
    }

    [Test]
    public void ClusterGenes_Binding_InvokesSuccessfully()
    {
        // Two members with 7/8 = 0.875 identity; threshold 0.8 -> single cluster,
        // AverageIdentity = 0.875 (single pair).
        var twoMember = MetagenomicsTools.ClusterGenes(
            new[] { Genome("g1", ("a", S)), Genome("g2", ("b", S1Sub)) },
            identityThreshold: 0.8);

        Assert.Multiple(() =>
        {
            Assert.That(twoMember.Items, Has.Count.EqualTo(1));
            Assert.That(twoMember.Items[0].GenomeCount, Is.EqualTo(2));
            Assert.That(twoMember.Items[0].AverageIdentity, Is.EqualTo(0.875).Within(1e-10),
                "Single pair -> AverageIdentity = global identity = 7/8 = 0.875 (CD-HIT -G 1).");
        });

        // Length difference, full identity over the shorter length; the longest member
        // (SLong) is the representative/consensus and is listed first.
        var lenDiff = MetagenomicsTools.ClusterGenes(
            new[] { Genome("g1", ("short", S)), Genome("g2", ("long", SLong)) },
            identityThreshold: 1.0);

        Assert.Multiple(() =>
        {
            Assert.That(lenDiff.Items, Has.Count.EqualTo(1),
                "Global identity over the shorter length is 1.0 -> one cluster even at threshold 1.0.");
            Assert.That(lenDiff.Items[0].ConsensusSequence, Is.EqualTo(SLong),
                "CD-HIT sorts long->short; the longest member is the representative/consensus.");
            Assert.That(lenDiff.Items[0].GeneIds[0], Is.EqualTo("long"),
                "The representative (longest) member is listed first.");
        });
    }
}
