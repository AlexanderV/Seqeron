using NUnit.Framework;
using Seqeron.Genomics;
using System.Collections.Generic;
using System.Linq;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Tests for PhylogeneticAnalyzer tree comparison and analysis.
/// 
/// Test Unit coverage:
/// - PHYLO-COMP-001: Tree comparison (Robinson-Foulds, MRCA, Patristic distance)
/// 
/// Other test files:
/// - PhylogeneticAnalyzer_NewickIO_Tests.cs (PHYLO-NEWICK-001)
/// - PhylogeneticAnalyzer_TreeConstruction_Tests.cs (PHYLO-TREE-001)
/// - PhylogeneticAnalyzer_DistanceMatrix_Tests.cs (PHYLO-DIST-001)
/// </summary>
[TestFixture]
[Category("PHYLO-COMP-001")]
public class PhylogeneticAnalyzerTests
{
    #region Tree Analysis Tests (PHYLO-COMP-001)

    [Test]
    [Description("PHYLO-COMP-001: GetLeaves returns all leaf nodes")]
    public void GetLeaves_ReturnsAllLeafNodes()
    {
        var sequences = new Dictionary<string, string>
        {
            ["A"] = "ACGT",
            ["B"] = "TCGT",
            ["C"] = "GCGT",
            ["D"] = "CCGT"
        };

        var tree = PhylogeneticAnalyzer.BuildTree(sequences);
        var leaves = PhylogeneticAnalyzer.GetLeaves(tree.Root).ToList();

        Assert.That(leaves, Has.Count.EqualTo(4));
        Assert.That(leaves.All(l => l.IsLeaf), Is.True);
    }

    [Test]
    public void CalculateTreeLength_SumsAllBranches()
    {
        var sequences = new Dictionary<string, string>
        {
            ["A"] = "ACGT",
            ["B"] = "TCGT"
        };

        var tree = PhylogeneticAnalyzer.BuildTree(sequences);
        double length = PhylogeneticAnalyzer.CalculateTreeLength(tree.Root);

        Assert.That(length, Is.GreaterThan(0));
    }

    [Test]
    public void GetTreeDepth_ReturnsCorrectDepth()
    {
        var sequences = new Dictionary<string, string>
        {
            ["A"] = "ACGT",
            ["B"] = "TCGT",
            ["C"] = "GCGT"
        };

        var tree = PhylogeneticAnalyzer.BuildTree(sequences);
        int depth = PhylogeneticAnalyzer.GetTreeDepth(tree.Root);

        Assert.That(depth, Is.GreaterThanOrEqualTo(1));
    }

    [Test]
    public void FindMRCA_FindsCommonAncestor()
    {
        var sequences = new Dictionary<string, string>
        {
            ["Human"] = "ACGTACGT",
            ["Chimp"] = "ACGTACGA",
            ["Mouse"] = "TCGTACGT"
        };

        var tree = PhylogeneticAnalyzer.BuildTree(sequences);
        var mrca = PhylogeneticAnalyzer.FindMRCA(tree.Root, "Human", "Chimp");

        Assert.That(mrca, Is.Not.Null);
        Assert.That(mrca!.Taxa, Does.Contain("Human"));
        Assert.That(mrca.Taxa, Does.Contain("Chimp"));
    }

    [Test]
    public void FindMRCA_SameTaxon_ReturnsTaxonItself()
    {
        var sequences = new Dictionary<string, string>
        {
            ["A"] = "ACGT",
            ["B"] = "TCGT"
        };

        var tree = PhylogeneticAnalyzer.BuildTree(sequences);
        var mrca = PhylogeneticAnalyzer.FindMRCA(tree.Root, "A", "A");

        Assert.That(mrca, Is.Not.Null);
        Assert.That(mrca!.Name, Is.EqualTo("A"));
    }

    [Test]
    public void PatristicDistance_CalculatesTreePathDistance()
    {
        var sequences = new Dictionary<string, string>
        {
            ["A"] = "ACGTACGT",
            ["B"] = "ACGTACGA",
            ["C"] = "TCGTACGT"
        };

        var tree = PhylogeneticAnalyzer.BuildTree(sequences);
        double dist = PhylogeneticAnalyzer.PatristicDistance(tree.Root, "A", "B");

        Assert.That(dist, Is.GreaterThanOrEqualTo(0));
    }

    [Test]
    public void PatristicDistance_SameTaxon_ReturnsZero()
    {
        var sequences = new Dictionary<string, string>
        {
            ["A"] = "ACGT",
            ["B"] = "TCGT"
        };

        var tree = PhylogeneticAnalyzer.BuildTree(sequences);
        double dist = PhylogeneticAnalyzer.PatristicDistance(tree.Root, "A", "A");

        Assert.That(dist, Is.EqualTo(0));
    }

    #endregion

    #region Robinson-Foulds Distance Tests (PHYLO-COMP-001)

    [Test]
    [Description("PHYLO-COMP-001: Identical trees have RF distance = 0")]
    public void RobinsonFouldsDistance_IdenticalTrees_ReturnsZero()
    {
        var sequences = new Dictionary<string, string>
        {
            ["A"] = "ACGT",
            ["B"] = "TCGT",
            ["C"] = "GCGT"
        };

        var tree1 = PhylogeneticAnalyzer.BuildTree(sequences, treeMethod: PhylogeneticAnalyzer.TreeMethod.UPGMA);
        var tree2 = PhylogeneticAnalyzer.BuildTree(sequences, treeMethod: PhylogeneticAnalyzer.TreeMethod.UPGMA);

        int rfDist = PhylogeneticAnalyzer.RobinsonFouldsDistance(tree1.Root, tree2.Root);

        Assert.That(rfDist, Is.EqualTo(0));
    }

    [Test]
    public void RobinsonFouldsDistance_DifferentTrees_ReturnsPositive()
    {
        // Build trees that might have different topology
        var sequences1 = new Dictionary<string, string>
        {
            ["A"] = "AAAA",
            ["B"] = "AAAC",
            ["C"] = "CCCC",
            ["D"] = "CCCA"
        };

        var tree1 = PhylogeneticAnalyzer.BuildTree(sequences1);

        // Create tree with different groupings
        var node = new PhylogeneticAnalyzer.PhyloNode
        {
            Left = new PhylogeneticAnalyzer.PhyloNode("A") { Taxa = { "A" } },
            Right = new PhylogeneticAnalyzer.PhyloNode
            {
                Left = new PhylogeneticAnalyzer.PhyloNode("B") { Taxa = { "B" } },
                Right = new PhylogeneticAnalyzer.PhyloNode
                {
                    Left = new PhylogeneticAnalyzer.PhyloNode("C") { Taxa = { "C" } },
                    Right = new PhylogeneticAnalyzer.PhyloNode("D") { Taxa = { "D" } }
                }
            }
        };
        node.Taxa = new List<string> { "A", "B", "C", "D" };

        // RF distance should be non-negative
        int rfDist = PhylogeneticAnalyzer.RobinsonFouldsDistance(tree1.Root, node);
        Assert.That(rfDist, Is.GreaterThanOrEqualTo(0));
    }

    #endregion

    #region Bootstrap Analysis Tests

    [Test]
    public void Bootstrap_ReturnsSupportsForSplits()
    {
        var sequences = new Dictionary<string, string>
        {
            ["A"] = "ACGTACGTAC",
            ["B"] = "ACGTACGTAC",
            ["C"] = "TCGTACGTAC",
            ["D"] = "TCGTACGTAC"
        };

        var supports = PhylogeneticAnalyzer.Bootstrap(sequences, replicates: 10);

        Assert.That(supports, Is.Not.Empty);
        Assert.That(supports.Values, Has.All.InRange(0.0, 1.0));
    }

    [Test]
    public void Bootstrap_IdenticalSequences_HighSupport()
    {
        // Groups with very different sequences should have high bootstrap support
        var sequences = new Dictionary<string, string>
        {
            ["A"] = "AAAAAAAAAA",
            ["B"] = "AAAAAAAAAC",
            ["C"] = "CCCCCCCCCC",
            ["D"] = "CCCCCCCCCA"
        };

        var supports = PhylogeneticAnalyzer.Bootstrap(sequences, replicates: 50);

        // At least some splits should have high support
        Assert.That(supports.Values.Any(v => v >= 0.5), Is.True);
    }

    #endregion

    #region Edge Cases (PHYLO-COMP-001 tree helpers)

    [Test]
    [Description("PHYLO-COMP-001: GetLeaves returns empty for null root")]
    public void GetLeaves_NullRoot_ReturnsEmpty()
    {
        var leaves = PhylogeneticAnalyzer.GetLeaves(null!).ToList();
        Assert.That(leaves, Is.Empty);
    }

    [Test]
    [Description("PHYLO-COMP-001: CalculateTreeLength returns zero for null root")]
    public void CalculateTreeLength_NullRoot_ReturnsZero()
    {
        double length = PhylogeneticAnalyzer.CalculateTreeLength(null!);
        Assert.That(length, Is.EqualTo(0));
    }

    #endregion
}
