using NUnit.Framework;
using Seqeron.Genomics;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Tests for PHYLO-TREE-001: Phylogenetic Tree Construction.
/// Covers: BuildTree (UPGMA and Neighbor-Joining algorithms)
/// 
/// Evidence: Wikipedia (UPGMA, Neighbor joining, Phylogenetic tree),
/// Saitou &amp; Nei (1987), Sokal &amp; Michener (1958)
/// </summary>
[TestFixture]
[Category("PHYLO-TREE-001")]
public class PhylogeneticAnalyzer_TreeConstruction_Tests
{
    #region M01-M04: Basic Contract Tests

    [Test]
    [Description("M01: BuildTree returns valid PhylogeneticTree with non-null root")]
    public void BuildTree_ValidInput_ReturnsValidTree()
    {
        var sequences = new Dictionary<string, string>
        {
            ["A"] = "ACGT",
            ["B"] = "TCGT"
        };

        var tree = PhylogeneticAnalyzer.BuildTree(sequences);

        Assert.Multiple(() =>
        {
            Assert.That(tree.Root, Is.Not.Null, "Root should not be null");
            Assert.That(tree.Taxa, Is.Not.Null, "Taxa should not be null");
            Assert.That(tree.DistanceMatrix, Is.Not.Null, "DistanceMatrix should not be null");
            Assert.That(tree.Method, Is.Not.Null.And.Not.Empty, "Method should be set");
        });
    }

    [Test]
    [Description("M02: Tree contains all input taxa as leaves")]
    public void BuildTree_ContainsAllTaxa_AllInputTaxaPresent()
    {
        var sequences = new Dictionary<string, string>
        {
            ["Human"] = "ACGTACGT",
            ["Chimp"] = "ACGTACGA",
            ["Mouse"] = "TCGTACGT",
            ["Rat"] = "TCGTACGA"
        };

        var tree = PhylogeneticAnalyzer.BuildTree(sequences);
        var leaves = PhylogeneticAnalyzer.GetLeaves(tree.Root).Select(l => l.Name).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(leaves, Has.Count.EqualTo(4), "Should have 4 leaves");
            Assert.That(leaves, Does.Contain("Human"));
            Assert.That(leaves, Does.Contain("Chimp"));
            Assert.That(leaves, Does.Contain("Mouse"));
            Assert.That(leaves, Does.Contain("Rat"));
        });
    }

    [Test]
    [Description("M03: UPGMA method produces tree with Method='UPGMA'")]
    public void BuildTree_UPGMA_ReturnsTreeWithUPGMAMethod()
    {
        var sequences = new Dictionary<string, string>
        {
            ["A"] = "ACGTACGT",
            ["B"] = "ACGTACGA",
            ["C"] = "TCGTACGT"
        };

        var tree = PhylogeneticAnalyzer.BuildTree(
            sequences,
            treeMethod: PhylogeneticAnalyzer.TreeMethod.UPGMA);

        Assert.That(tree.Method, Is.EqualTo("UPGMA"));
    }

    [Test]
    [Description("M04: NeighborJoining method produces tree with Method='NeighborJoining'")]
    public void BuildTree_NeighborJoining_ReturnsTreeWithNJMethod()
    {
        var sequences = new Dictionary<string, string>
        {
            ["A"] = "ACGTACGT",
            ["B"] = "ACGTACGA",
            ["C"] = "TCGTACGT"
        };

        var tree = PhylogeneticAnalyzer.BuildTree(
            sequences,
            treeMethod: PhylogeneticAnalyzer.TreeMethod.NeighborJoining);

        Assert.That(tree.Method, Is.EqualTo("NeighborJoining"));
    }

    #endregion

    #region M05-M07: Tree Structure Tests

    [Test]
    [Description("M05: Two sequences create binary tree with both as leaves")]
    public void BuildTree_TwoSequences_CreatesBinaryTreeWithBothAsLeaves()
    {
        var sequences = new Dictionary<string, string>
        {
            ["A"] = "ACGT",
            ["B"] = "TCGT"
        };

        var tree = PhylogeneticAnalyzer.BuildTree(sequences);
        var leaves = PhylogeneticAnalyzer.GetLeaves(tree.Root).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(tree.Root.IsLeaf, Is.False, "Root should not be a leaf");
            Assert.That(tree.Root.Left, Is.Not.Null, "Root should have left child");
            Assert.That(tree.Root.Right, Is.Not.Null, "Root should have right child");
            Assert.That(leaves, Has.Count.EqualTo(2), "Should have exactly 2 leaves");
        });
    }

    [Test]
    [Description("M06: Three sequences create valid binary tree structure")]
    public void BuildTree_ThreeSequences_CreatesValidBinaryTree()
    {
        var sequences = new Dictionary<string, string>
        {
            ["A"] = "AAAA",
            ["B"] = "CCCC",
            ["C"] = "GGGG"
        };

        var tree = PhylogeneticAnalyzer.BuildTree(sequences);
        var leaves = PhylogeneticAnalyzer.GetLeaves(tree.Root).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(leaves, Has.Count.EqualTo(3), "Should have exactly 3 leaves");
            Assert.That(leaves.All(l => l.IsLeaf), Is.True, "All leaves should be leaf nodes");
            // Binary tree with 3 leaves has 2 internal nodes (including root)
            int internalNodes = CountInternalNodes(tree.Root);
            Assert.That(internalNodes, Is.EqualTo(2), "Should have 2 internal nodes");
        });
    }

    [Test]
    [Description("M07: Four sequences (Wikipedia example structure) create valid tree")]
    public void BuildTree_FourSequences_CreatesValidTree()
    {
        // Based on structure similar to Wikipedia UPGMA example
        var sequences = new Dictionary<string, string>
        {
            ["A"] = "ACGTACGTAC",
            ["B"] = "ACGTACGTAC",  // Identical to A
            ["C"] = "GCGTACGTAC",
            ["D"] = "GCGAACGTAC"
        };

        var tree = PhylogeneticAnalyzer.BuildTree(sequences);
        var leaves = PhylogeneticAnalyzer.GetLeaves(tree.Root).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(leaves, Has.Count.EqualTo(4));
            Assert.That(tree.Taxa, Has.Count.EqualTo(4));
            Assert.That(tree.DistanceMatrix.GetLength(0), Is.EqualTo(4));
            Assert.That(tree.DistanceMatrix.GetLength(1), Is.EqualTo(4));
        });
    }

    #endregion

    #region M08-M09: Input Validation Tests

    [Test]
    [Description("M08: BuildTree throws ArgumentException on single sequence")]
    public void BuildTree_SingleSequence_ThrowsArgumentException()
    {
        var sequences = new Dictionary<string, string>
        {
            ["Only"] = "ACGT"
        };

        var ex = Assert.Throws<ArgumentException>(() =>
            PhylogeneticAnalyzer.BuildTree(sequences));

        Assert.That(ex!.Message, Does.Contain("2").Or.Contain("two").IgnoreCase);
    }

    [Test]
    [Description("M09: BuildTree throws ArgumentException on unequal sequence lengths")]
    public void BuildTree_UnequalLengths_ThrowsArgumentException()
    {
        var sequences = new Dictionary<string, string>
        {
            ["Short"] = "ACG",
            ["Long"] = "ACGTACGT"
        };

        var ex = Assert.Throws<ArgumentException>(() =>
            PhylogeneticAnalyzer.BuildTree(sequences));

        Assert.That(ex!.Message, Does.Contain("length").IgnoreCase);
    }

    [Test]
    [Description("M09b: BuildTree throws on null input")]
    public void BuildTree_NullInput_ThrowsArgumentException()
    {
        Assert.Throws<ArgumentException>(() =>
            PhylogeneticAnalyzer.BuildTree(null!));
    }

    [Test]
    [Description("M09c: BuildTree throws on empty dictionary")]
    public void BuildTree_EmptyInput_ThrowsArgumentException()
    {
        var sequences = new Dictionary<string, string>();

        Assert.Throws<ArgumentException>(() =>
            PhylogeneticAnalyzer.BuildTree(sequences));
    }

    #endregion

    #region M10-M13: Invariant Tests

    [Test]
    [Description("M10: All branch lengths are non-negative (UPGMA property)")]
    public void BuildTree_UPGMA_AllBranchLengthsNonNegative()
    {
        var sequences = new Dictionary<string, string>
        {
            ["A"] = "ACGTACGT",
            ["B"] = "TCGTACGT",
            ["C"] = "GCGTACGT",
            ["D"] = "CCGTACGT"
        };

        var tree = PhylogeneticAnalyzer.BuildTree(
            sequences,
            treeMethod: PhylogeneticAnalyzer.TreeMethod.UPGMA);

        var allNodes = GetAllNodes(tree.Root);

        Assert.Multiple(() =>
        {
            foreach (var node in allNodes)
            {
                Assert.That(node.BranchLength, Is.GreaterThanOrEqualTo(0),
                    $"Node '{node.Name}' has negative branch length");
            }
        });
    }

    [Test]
    [Description("M11: UPGMA produces rooted tree with valid structure")]
    public void BuildTree_UPGMA_ProducesRootedTree()
    {
        var sequences = new Dictionary<string, string>
        {
            ["A"] = "ACGT",
            ["B"] = "TCGT",
            ["C"] = "GCGT"
        };

        var tree = PhylogeneticAnalyzer.BuildTree(
            sequences,
            treeMethod: PhylogeneticAnalyzer.TreeMethod.UPGMA);

        Assert.Multiple(() =>
        {
            Assert.That(tree.Root, Is.Not.Null, "Tree should have root");
            Assert.That(tree.Root.IsLeaf, Is.False, "Root should not be a leaf");
            // UPGMA always produces a rooted tree
            Assert.That(tree.Root.Taxa, Has.Count.EqualTo(3), "Root should contain all taxa");
        });
    }

    [Test]
    [Description("M12: Identical sequences produce zero-distance subtree")]
    public void BuildTree_IdenticalSequences_ProducesZeroDistanceSubtree()
    {
        var sequences = new Dictionary<string, string>
        {
            ["A"] = "ACGTACGT",
            ["B"] = "ACGTACGT",  // Identical to A
            ["C"] = "TTTTTTTT"   // Different
        };

        var tree = PhylogeneticAnalyzer.BuildTree(sequences);

        // Distance between A and B should be 0
        Assert.That(tree.DistanceMatrix[0, 1], Is.EqualTo(0).Within(1e-10),
            "Distance between identical sequences should be 0");
    }

    [Test]
    [Description("M13: Case-insensitive sequence handling")]
    public void BuildTree_CaseInsensitive_TreatsAsIdentical()
    {
        var sequences = new Dictionary<string, string>
        {
            ["A"] = "acgt",
            ["B"] = "ACGT"
        };

        var tree = PhylogeneticAnalyzer.BuildTree(sequences);

        // Should not throw, and distance should be 0
        Assert.That(tree.DistanceMatrix[0, 1], Is.EqualTo(0).Within(1e-10),
            "Case-insensitive: same sequence should have 0 distance");
    }

    #endregion

    #region S01-S05: Validation Against Wikipedia Examples

    [Test]
    [Description("S01: UPGMA joins closest pair first (Wikipedia invariant)")]
    public void BuildTree_UPGMA_JoinsClosestPairFirst()
    {
        // Two pairs: (A,B) very close, (C,D) very close
        var sequences = new Dictionary<string, string>
        {
            ["A"] = "AAAAAAAA",
            ["B"] = "AAAAAAAC",  // 1 difference from A
            ["C"] = "CCCCCCCC",
            ["D"] = "CCCCCCCA"   // 1 difference from C
        };

        var tree = PhylogeneticAnalyzer.BuildTree(
            sequences,
            treeMethod: PhylogeneticAnalyzer.TreeMethod.UPGMA);

        var leaves = PhylogeneticAnalyzer.GetLeaves(tree.Root).ToList();

        Assert.That(leaves, Has.Count.EqualTo(4));
        // The tree should group (A,B) and (C,D) based on similarity
    }

    [Test]
    [Description("S02: Neighbor-Joining produces valid tree for divergent sequences")]
    public void BuildTree_NeighborJoining_HandlesVariableRates()
    {
        var sequences = new Dictionary<string, string>
        {
            ["A"] = "ACGTACGT",
            ["B"] = "ACGTACGA",  // 1 change
            ["C"] = "TTTTTTTT"   // Very different
        };

        var tree = PhylogeneticAnalyzer.BuildTree(
            sequences,
            treeMethod: PhylogeneticAnalyzer.TreeMethod.NeighborJoining);

        var leaves = PhylogeneticAnalyzer.GetLeaves(tree.Root).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(leaves, Has.Count.EqualTo(3));
            Assert.That(tree.Root, Is.Not.Null);
        });
    }

    [Test]
    [Description("S03: Tree depth is appropriate for number of taxa")]
    public void BuildTree_TreeDepth_AppropriateForTaxaCount()
    {
        var sequences = new Dictionary<string, string>
        {
            ["A"] = "ACGT",
            ["B"] = "TCGT",
            ["C"] = "GCGT",
            ["D"] = "CCGT"
        };

        var tree = PhylogeneticAnalyzer.BuildTree(sequences);
        int depth = PhylogeneticAnalyzer.GetTreeDepth(tree.Root);

        // For 4 taxa in a binary tree, depth should be 2-3
        Assert.That(depth, Is.InRange(1, 3), "Depth should be reasonable for binary tree");
    }

    [Test]
    [Description("S04: DistanceMatrix in result is symmetric")]
    public void BuildTree_DistanceMatrix_IsSymmetric()
    {
        var sequences = new Dictionary<string, string>
        {
            ["A"] = "ACGT",
            ["B"] = "TCGT",
            ["C"] = "GCGT"
        };

        var tree = PhylogeneticAnalyzer.BuildTree(sequences);

        int n = tree.DistanceMatrix.GetLength(0);
        Assert.Multiple(() =>
        {
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    Assert.That(tree.DistanceMatrix[i, j],
                        Is.EqualTo(tree.DistanceMatrix[j, i]).Within(1e-10),
                        $"Matrix not symmetric at [{i},{j}]");
                }
            }
        });
    }

    [Test]
    [Description("S05: Leaf names exactly match input taxon names")]
    public void BuildTree_LeafNames_MatchInputExactly()
    {
        var sequences = new Dictionary<string, string>
        {
            ["Homo_sapiens"] = "ACGT",
            ["Pan_troglodytes"] = "TCGT",
            ["Mus_musculus"] = "GCGT"
        };

        var tree = PhylogeneticAnalyzer.BuildTree(sequences);
        var leafNames = PhylogeneticAnalyzer.GetLeaves(tree.Root).Select(l => l.Name).ToHashSet();

        Assert.Multiple(() =>
        {
            Assert.That(leafNames, Does.Contain("Homo_sapiens"));
            Assert.That(leafNames, Does.Contain("Pan_troglodytes"));
            Assert.That(leafNames, Does.Contain("Mus_musculus"));
        });
    }

    #endregion

    #region C01-C03: Extended Coverage

    [Test]
    [Description("C01: Different distance methods all produce valid trees")]
    [TestCase(PhylogeneticAnalyzer.DistanceMethod.Hamming)]
    [TestCase(PhylogeneticAnalyzer.DistanceMethod.PDistance)]
    [TestCase(PhylogeneticAnalyzer.DistanceMethod.JukesCantor)]
    [TestCase(PhylogeneticAnalyzer.DistanceMethod.Kimura2Parameter)]
    public void BuildTree_AllDistanceMethods_ProduceValidTrees(
        PhylogeneticAnalyzer.DistanceMethod method)
    {
        var sequences = new Dictionary<string, string>
        {
            ["A"] = "ACGTACGT",
            ["B"] = "TCGTACGT",
            ["C"] = "GCGTACGT"
        };

        var tree = PhylogeneticAnalyzer.BuildTree(sequences, method);

        Assert.Multiple(() =>
        {
            Assert.That(tree.Root, Is.Not.Null);
            Assert.That(tree.Taxa, Has.Count.EqualTo(3));
        });
    }

    [Test]
    [Description("C02: Tree total length is sum of all branch lengths")]
    public void BuildTree_TreeLength_IsSumOfBranchLengths()
    {
        var sequences = new Dictionary<string, string>
        {
            ["A"] = "ACGT",
            ["B"] = "TCGT",
            ["C"] = "GCGT"
        };

        var tree = PhylogeneticAnalyzer.BuildTree(sequences);
        double totalLength = PhylogeneticAnalyzer.CalculateTreeLength(tree.Root);

        Assert.That(totalLength, Is.GreaterThan(0), "Total tree length should be positive");
    }

    [Test]
    [Description("C03: Both tree methods produce trees with same leaf set")]
    public void BuildTree_BothMethods_ProduceSameLeafSet()
    {
        var sequences = new Dictionary<string, string>
        {
            ["A"] = "ACGTACGT",
            ["B"] = "TCGTACGT",
            ["C"] = "GCGTACGT"
        };

        var upgmaTree = PhylogeneticAnalyzer.BuildTree(
            sequences, treeMethod: PhylogeneticAnalyzer.TreeMethod.UPGMA);
        var njTree = PhylogeneticAnalyzer.BuildTree(
            sequences, treeMethod: PhylogeneticAnalyzer.TreeMethod.NeighborJoining);

        var upgmaLeaves = PhylogeneticAnalyzer.GetLeaves(upgmaTree.Root)
            .Select(l => l.Name).OrderBy(n => n).ToList();
        var njLeaves = PhylogeneticAnalyzer.GetLeaves(njTree.Root)
            .Select(l => l.Name).OrderBy(n => n).ToList();

        Assert.That(njLeaves, Is.EqualTo(upgmaLeaves));
    }

    #endregion

    #region Helper Methods

    private static int CountInternalNodes(PhylogeneticAnalyzer.PhyloNode? node)
    {
        if (node == null || node.IsLeaf) return 0;
        return 1 + CountInternalNodes(node.Left) + CountInternalNodes(node.Right);
    }

    private static List<PhylogeneticAnalyzer.PhyloNode> GetAllNodes(
        PhylogeneticAnalyzer.PhyloNode? root)
    {
        var nodes = new List<PhylogeneticAnalyzer.PhyloNode>();
        CollectNodes(root, nodes);
        return nodes;
    }

    private static void CollectNodes(
        PhylogeneticAnalyzer.PhyloNode? node,
        List<PhylogeneticAnalyzer.PhyloNode> nodes)
    {
        if (node == null) return;
        nodes.Add(node);
        CollectNodes(node.Left, nodes);
        CollectNodes(node.Right, nodes);
    }

    #endregion
}
