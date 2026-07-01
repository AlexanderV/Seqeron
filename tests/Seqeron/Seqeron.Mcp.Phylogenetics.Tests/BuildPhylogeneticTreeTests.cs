using NUnit.Framework;
using Seqeron.Mcp.Phylogenetics.Tools;

namespace Seqeron.Mcp.Phylogenetics.Tests;

[TestFixture]
public class BuildPhylogeneticTreeTests
{
    [Test]
    public void BuildPhylogeneticTree_Schema_ValidatesCorrectly()
    {
        var ok = new Dictionary<string, string> { ["A"] = "ACGT", ["B"] = "TCGT" };
        Assert.DoesNotThrow(() => PhylogeneticsTools.BuildPhylogeneticTree(ok));

        // Null → ArgumentException.
        Assert.Throws<ArgumentException>(() => PhylogeneticsTools.BuildPhylogeneticTree(null!));

        // Single sequence → ArgumentException (a tree needs ≥2 taxa).
        Assert.Throws<ArgumentException>(() =>
            PhylogeneticsTools.BuildPhylogeneticTree(new Dictionary<string, string> { ["A"] = "ACGT" }));

        // Unknown methods → ArgumentException.
        Assert.Throws<ArgumentException>(() => PhylogeneticsTools.BuildPhylogeneticTree(ok, distanceMethod: "bogus"));
        Assert.Throws<ArgumentException>(() => PhylogeneticsTools.BuildPhylogeneticTree(ok, treeMethod: "bogus"));
    }

    [Test]
    public void BuildPhylogeneticTree_Binding_InvokesSuccessfully()
    {
        // A and B identical; C differs. UPGMA. Expected (Seqeron.Genomics.Tests PHYLO-TREE-001 M12):
        // distance[0,1] = 0 (identical sequences). Taxa preserved; Newick well-formed.
        var sequences = new Dictionary<string, string>
        {
            ["A"] = "ACGTACGT",
            ["B"] = "ACGTACGT",
            ["C"] = "TTTTTTTT",
        };

        var result = PhylogeneticsTools.BuildPhylogeneticTree(
            sequences, distanceMethod: "JukesCantor", treeMethod: "UPGMA");

        Assert.Multiple(() =>
        {
            Assert.That(result.Method, Is.EqualTo("UPGMA"));
            Assert.That(result.Taxa, Is.EquivalentTo(new[] { "A", "B", "C" }));
            Assert.That(result.Newick, Does.EndWith(";"));
            Assert.That(result.Newick, Does.Contain("A"));
            Assert.That(result.Newick, Does.Contain("B"));
            Assert.That(result.Newick, Does.Contain("C"));

            // 3×3 matrix, diagonal zero, symmetric.
            Assert.That(result.DistanceMatrix.Length, Is.EqualTo(3));
            Assert.That(result.DistanceMatrix[0].Length, Is.EqualTo(3));
            Assert.That(result.DistanceMatrix[0][0], Is.EqualTo(0.0).Within(1e-10));
            Assert.That(result.DistanceMatrix[1][1], Is.EqualTo(0.0).Within(1e-10));
            Assert.That(result.DistanceMatrix[0][1], Is.EqualTo(result.DistanceMatrix[1][0]).Within(1e-10));

            // A,B identical → distance 0 (PHYLO-TREE-001 M12).
            Assert.That(result.DistanceMatrix[0][1], Is.EqualTo(0.0).Within(1e-10));
            // A,C differ at every comparable site → JC saturates to +Infinity.
            Assert.That(double.IsPositiveInfinity(result.DistanceMatrix[0][2]), Is.True);
        });

        // NeighborJoining produces Method="NeighborJoining" (PHYLO-TREE-001 M04).
        var nj = PhylogeneticsTools.BuildPhylogeneticTree(sequences, treeMethod: "NeighborJoining");
        Assert.That(nj.Method, Is.EqualTo("NeighborJoining"));
    }
}
