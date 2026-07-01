using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>find_conserved_clusters</c> MCP tool.
/// Expected results derived from the common-interval model in
/// ComparativeGenomics.FindConservedClusters: a set of ortholog-group labels that is a
/// contiguous interval in every genome. NOT the wrapper's output.
/// </summary>
[TestFixture]
public class FindConservedClustersTests
{
    private static GeneInput G(string id, int start) => new(id, "genome", start, start + 10, '+');

    [Test]
    public void FindConservedClusters_Schema_ValidatesCorrectly()
    {
        var g = new[] { new[] { G("A", 0), G("B", 20) }, new[] { G("A", 0), G("B", 20) } };
        var groups = new Dictionary<string, string> { ["A"] = "g1", ["B"] = "g2" };

        Assert.DoesNotThrow(() => AnalysisTools.FindConservedClusters(g, groups, 2, 2));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindConservedClusters(null!, groups));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindConservedClusters(Array.Empty<GeneInput[]>(), groups));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindConservedClusters(g, null!));
    }

    [Test]
    public void FindConservedClusters_IdenticalOrder_ReturnsAllIntervals()
    {
        // Two genomes, same order A,B,C -> groups g1,g2,g3. With minClusterSize 2 the common
        // intervals are {g1,g2}, {g2,g3} (size 2) and {g1,g2,g3} (size 3),
        // sorted by size then lexicographically.
        var genome = new[] { G("A", 0), G("B", 20), G("C", 40) };
        var genomes = new[] { genome, genome };
        var groups = new Dictionary<string, string> { ["A"] = "g1", ["B"] = "g2", ["C"] = "g3" };

        var clusters = AnalysisTools.FindConservedClusters(genomes, groups, 2, 2).Clusters;
        Assert.That(clusters, Has.Length.EqualTo(3));
        Assert.Multiple(() =>
        {
            Assert.That(clusters[0], Is.EqualTo(new[] { "g1", "g2" }));
            Assert.That(clusters[1], Is.EqualTo(new[] { "g2", "g3" }));
            Assert.That(clusters[2], Is.EqualTo(new[] { "g1", "g2", "g3" }));
        });
    }

    [Test]
    public void FindConservedClusters_SingleGenome_ReturnsEmpty()
    {
        // A common interval needs >= 2 genomes; a single genome yields no clusters.
        var genomes = new[] { new[] { G("A", 0), G("B", 20) } };
        var groups = new Dictionary<string, string> { ["A"] = "g1", ["B"] = "g2" };
        var clusters = AnalysisTools.FindConservedClusters(genomes, groups, 2, 2).Clusters;
        Assert.That(clusters, Is.Empty);
    }
}
