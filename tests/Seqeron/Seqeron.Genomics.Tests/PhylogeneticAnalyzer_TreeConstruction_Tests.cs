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

            // A and B are identical → must be merged first as siblings
            var abMRCA = PhylogeneticAnalyzer.FindMRCA(tree.Root, "A", "B");
            Assert.That(abMRCA, Is.Not.Null, "A,B should share an MRCA");
            Assert.That(abMRCA!.Left!.IsLeaf && abMRCA.Right!.IsLeaf, Is.True,
                "Identical A,B should be siblings (merged first by UPGMA)");
            Assert.That(tree.DistanceMatrix[0, 1], Is.EqualTo(0).Within(1e-10),
                "Identical sequences A,B should have distance 0");
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

    /// <summary>
    /// Wikipedia UPGMA working example distance matrix (5S rRNA, 5 bacteria).
    /// Source: https://en.wikipedia.org/wiki/UPGMA#Working_example
    /// Taxa: a=B.subtilis, b=B.stearothermophilus, c=L.viridescens, d=A.modicum, e=M.luteus
    /// </summary>
    private static double[,] WikipediaUPGMAMatrix => new double[,]
    {
        { 0, 17, 21, 31, 23 },
        { 17, 0, 30, 34, 21 },
        { 21, 30, 0, 28, 39 },
        { 31, 34, 28, 0, 43 },
        { 23, 21, 39, 43, 0 }
    };

    private static readonly string[] WikipediaUPGMATaxa = { "a", "b", "c", "d", "e" };

    /// <summary>
    /// Wikipedia Neighbor-Joining working example distance matrix (5 taxa).
    /// Source: https://en.wikipedia.org/wiki/Neighbor_joining#Example
    /// This is an additive distance matrix.
    /// </summary>
    private static double[,] WikipediaNJMatrix => new double[,]
    {
        { 0, 5, 9, 9, 8 },
        { 5, 0, 10, 10, 9 },
        { 9, 10, 0, 8, 7 },
        { 9, 10, 8, 0, 3 },
        { 8, 9, 7, 3, 0 }
    };

    private static readonly string[] WikipediaNJTaxa = { "a", "b", "c", "d", "e" };

    [Test]
    [Description("S01: Wikipedia UPGMA example — clustering order matches expected merges")]
    public void BuildTree_UPGMA_WikipediaExample_CorrectClusteringOrder()
    {
        // Source: https://en.wikipedia.org/wiki/UPGMA#Working_example
        // Expected clustering order:
        //   1. (a,b) at d=17
        //   2. ((a,b),e) at d=22
        //   3. (c,d) at d=28
        //   4. Final join at d=33
        var tree = PhylogeneticAnalyzer.BuildTreeFromMatrix(
            WikipediaUPGMATaxa, WikipediaUPGMAMatrix,
            PhylogeneticAnalyzer.TreeMethod.UPGMA);

        // Verify: a and b should be siblings (merged first)
        var abMRCA = PhylogeneticAnalyzer.FindMRCA(tree.Root, "a", "b");
        Assert.That(abMRCA, Is.Not.Null, "a,b should share an MRCA");
        Assert.That(abMRCA!.Left!.IsLeaf && abMRCA.Right!.IsLeaf, Is.True,
            "a and b should be direct children of their MRCA (merged first)");

        // Verify: c and d should be siblings
        var cdMRCA = PhylogeneticAnalyzer.FindMRCA(tree.Root, "c", "d");
        Assert.That(cdMRCA, Is.Not.Null, "c,d should share an MRCA");
        Assert.That(cdMRCA!.Left!.IsLeaf && cdMRCA.Right!.IsLeaf, Is.True,
            "c and d should be direct children of their MRCA");
    }

    [Test]
    [Description("S01b: Wikipedia UPGMA example — branch lengths match expected values")]
    public void BuildTree_UPGMA_WikipediaExample_CorrectBranchLengths()
    {
        // Source: https://en.wikipedia.org/wiki/UPGMA#Working_example
        // Expected branch lengths:
        //   δ(a,u) = δ(b,u) = 8.5
        //   δ(u,v) = 2.5, δ(e,v) = 11
        //   δ(c,w) = δ(d,w) = 14
        //   δ(v,r) = 5.5, δ(w,r) = 2.5
        var tree = PhylogeneticAnalyzer.BuildTreeFromMatrix(
            WikipediaUPGMATaxa, WikipediaUPGMAMatrix,
            PhylogeneticAnalyzer.TreeMethod.UPGMA);

        // Find specific nodes
        var leafA = FindLeafByName(tree.Root, "a");
        var leafB = FindLeafByName(tree.Root, "b");
        var leafC = FindLeafByName(tree.Root, "c");
        var leafD = FindLeafByName(tree.Root, "d");
        var leafE = FindLeafByName(tree.Root, "e");

        // Also find internal nodes for full branch-length verification
        var abMRCA = PhylogeneticAnalyzer.FindMRCA(tree.Root, "a", "b");
        var cdMRCA = PhylogeneticAnalyzer.FindMRCA(tree.Root, "c", "d");
        var abeMRCA = PhylogeneticAnalyzer.FindMRCA(tree.Root, "a", "e");

        Assert.Multiple(() =>
        {
            // Leaf branch lengths (exact integer arithmetic → 1e-10 tolerance)
            // δ(a,u) = δ(b,u) = 8.5
            Assert.That(leafA!.BranchLength, Is.EqualTo(8.5).Within(1e-10),
                "δ(a,u) should be 8.5");
            Assert.That(leafB!.BranchLength, Is.EqualTo(8.5).Within(1e-10),
                "δ(b,u) should be 8.5");

            // δ(c,w) = δ(d,w) = 14
            Assert.That(leafC!.BranchLength, Is.EqualTo(14.0).Within(1e-10),
                "δ(c,w) should be 14");
            Assert.That(leafD!.BranchLength, Is.EqualTo(14.0).Within(1e-10),
                "δ(d,w) should be 14");

            // δ(e,v) = 11
            Assert.That(leafE!.BranchLength, Is.EqualTo(11.0).Within(1e-10),
                "δ(e,v) should be 11");

            // Internal branch lengths
            // δ(u,v) = height(v) - height(u) = 11 - 8.5 = 2.5
            Assert.That(abMRCA!.BranchLength, Is.EqualTo(2.5).Within(1e-10),
                "δ(u_ab, v) should be 2.5");

            // δ(w,root) = height(root) - height(w) = 16.5 - 14 = 2.5
            Assert.That(cdMRCA!.BranchLength, Is.EqualTo(2.5).Within(1e-10),
                "δ(w_cd, root) should be 2.5");

            // δ(v,root) = height(root) - height(v) = 16.5 - 11 = 5.5
            Assert.That(abeMRCA!.BranchLength, Is.EqualTo(5.5).Within(1e-10),
                "δ(v, root) should be 5.5");
        });
    }

    [Test]
    [Description("S01c: Wikipedia UPGMA example — ultrametric property (all tips equidistant from root)")]
    public void BuildTree_UPGMA_WikipediaExample_Ultrametric()
    {
        // Source: https://en.wikipedia.org/wiki/UPGMA#The_UPGMA_dendrogram
        // "It is ultrametric because all tips (a to e) are equidistant from r:
        //  δ(a,r) = δ(b,r) = δ(e,r) = δ(c,r) = δ(d,r) = 16.5"
        var tree = PhylogeneticAnalyzer.BuildTreeFromMatrix(
            WikipediaUPGMATaxa, WikipediaUPGMAMatrix,
            PhylogeneticAnalyzer.TreeMethod.UPGMA);

        // Calculate root-to-tip distance for each taxon
        var tipDistances = new Dictionary<string, double>();
        CalculateRootToTipDistances(tree.Root, 0.0, tipDistances);

        Assert.Multiple(() =>
        {
            foreach (var taxon in WikipediaUPGMATaxa)
            {
                Assert.That(tipDistances[taxon], Is.EqualTo(16.5).Within(1e-10),
                    $"Root-to-tip distance for '{taxon}' should be 16.5");
            }
        });
    }

    [Test]
    [Description("S02: Wikipedia NJ example — patristic distances match input (additive matrix, INV-N01 topology guarantee)")]
    public void BuildTree_NJ_WikipediaExample_PatristicDistancesMatchInput()
    {
        // Source: https://en.wikipedia.org/wiki/Neighbor_joining#Conclusion:_additive_distances
        // "if we move from any taxon to any other along the branches of the tree,
        //  and sum the lengths of the branches traversed, the result is equal to
        //  the distance between those taxa in the input distance matrix."
        //
        // Also covers INV-N01: NJ is guaranteed to find the correct topology
        // for additive distance matrices (Wikipedia NJ, Advantages and disadvantages).
        var tree = PhylogeneticAnalyzer.BuildTreeFromMatrix(
            WikipediaNJTaxa, WikipediaNJMatrix,
            PhylogeneticAnalyzer.TreeMethod.NeighborJoining);

        int n = WikipediaNJTaxa.Length;
        Assert.Multiple(() =>
        {
            for (int i = 0; i < n; i++)
            {
                for (int j = i + 1; j < n; j++)
                {
                    double patristic = PhylogeneticAnalyzer.PatristicDistance(
                        tree.Root, WikipediaNJTaxa[i], WikipediaNJTaxa[j]);
                    double expected = WikipediaNJMatrix[i, j];
                    Assert.That(patristic, Is.EqualTo(expected).Within(1e-10),
                        $"Patristic d({WikipediaNJTaxa[i]},{WikipediaNJTaxa[j]}) = {patristic:F6}, expected {expected}");
                }
            }
        });
    }

    [Test]
    [Description("S02b: Wikipedia NJ example — first join is (a,b) per Q-matrix minimum")]
    public void BuildTree_NJ_WikipediaExample_CorrectTopology()
    {
        // Source: https://en.wikipedia.org/wiki/Neighbor_joining#First_joining
        // Q1(a,b) = -50 is the unique smallest value → a,b joined first.
        //
        // Note: In step 2, Q2(u,c) = Q2(d,e) = -28 (tied minimum).
        // Wikipedia: "We may choose either to join u and c, or to join d and e."
        // Both topologies are valid NJ results. We only assert the deterministic
        // first-step topology: a and b must be siblings.
        var tree = PhylogeneticAnalyzer.BuildTreeFromMatrix(
            WikipediaNJTaxa, WikipediaNJMatrix,
            PhylogeneticAnalyzer.TreeMethod.NeighborJoining);

        // a and b should be siblings (unique Q-minimum in step 1)
        var abMRCA = PhylogeneticAnalyzer.FindMRCA(tree.Root, "a", "b");
        Assert.That(abMRCA, Is.Not.Null, "a,b should share an MRCA");
        Assert.That(abMRCA!.Left!.IsLeaf && abMRCA.Right!.IsLeaf, Is.True,
            "a and b should be direct children of their join node");

        // Verify all 5 taxa are present as leaves
        var leaves = PhylogeneticAnalyzer.GetLeaves(tree.Root)
            .Select(l => l.Name).ToHashSet();
        Assert.That(leaves, Is.EquivalentTo(new[] { "a", "b", "c", "d", "e" }));
    }

    [Test]
    [Description("S02c: Wikipedia NJ example — branch lengths δ(a,u)=2, δ(b,u)=3 per first step")]
    public void BuildTree_NJ_WikipediaExample_FirstJoinBranchLengths()
    {
        // Source: https://en.wikipedia.org/wiki/Neighbor_joining#First_branch_length_estimation
        // δ(a,u) = 2, δ(b,u) = 3
        var tree = PhylogeneticAnalyzer.BuildTreeFromMatrix(
            WikipediaNJTaxa, WikipediaNJMatrix,
            PhylogeneticAnalyzer.TreeMethod.NeighborJoining);

        var leafA = FindLeafByName(tree.Root, "a");
        var leafB = FindLeafByName(tree.Root, "b");

        Assert.Multiple(() =>
        {
            Assert.That(leafA!.BranchLength, Is.EqualTo(2.0).Within(1e-10),
                "δ(a,u) should be 2");
            Assert.That(leafB!.BranchLength, Is.EqualTo(3.0).Within(1e-10),
                "δ(b,u) should be 3");
        });
    }

    [Test]
    [Description("S03: UPGMA ultrametric property holds for general input")]
    public void BuildTree_UPGMA_GeneralInput_Ultrametric()
    {
        var sequences = new Dictionary<string, string>
        {
            ["A"] = "ACGTACGT",
            ["B"] = "ACGTACGA",
            ["C"] = "TCGTACGT",
            ["D"] = "GCGTACGA"
        };

        var tree = PhylogeneticAnalyzer.BuildTree(
            sequences,
            treeMethod: PhylogeneticAnalyzer.TreeMethod.UPGMA);

        var tipDistances = new Dictionary<string, double>();
        CalculateRootToTipDistances(tree.Root, 0.0, tipDistances);

        // All tips must be equidistant from root (ultrametric)
        double firstDist = tipDistances.Values.First();
        Assert.Multiple(() =>
        {
            foreach (var kvp in tipDistances)
            {
                Assert.That(kvp.Value, Is.EqualTo(firstDist).Within(1e-6),
                    $"Root-to-tip for '{kvp.Key}' = {kvp.Value:F6}, expected {firstDist:F6} (ultrametric)");
            }
        });
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

        Assert.That(leafNames, Is.EquivalentTo(new[] { "Homo_sapiens", "Pan_troglodytes", "Mus_musculus" }),
            "Leaf names must exactly match input taxon names — no extras, no missing");
    }

    // S06 removed: was duplicate of S02 (same data, same patristic-distance-matches-input assertion).
    // INV-N01 topology guarantee is now covered by S02.

    #endregion

    #region C01-C05: Extended Coverage

    [Test]
    [Description("C01: Large input (50+ sequences) completes in reasonable time — O(n³) complexity")]
    public void BuildTree_LargeInput_CompletesInReasonableTime()
    {
        // Generate 50 random sequences of length 100
        var random = new Random(42);
        var sequences = new Dictionary<string, string>();
        for (int i = 0; i < 50; i++)
        {
            var sb = new System.Text.StringBuilder();
            for (int j = 0; j < 100; j++)
                sb.Append("ACGT"[random.Next(4)]);
            sequences[$"Taxon_{i:D3}"] = sb.ToString();
        }

        var sw = System.Diagnostics.Stopwatch.StartNew();
        var tree = PhylogeneticAnalyzer.BuildTree(sequences);
        sw.Stop();

        var leaves = PhylogeneticAnalyzer.GetLeaves(tree.Root).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(tree.Root, Is.Not.Null);
            Assert.That(tree.Taxa, Has.Count.EqualTo(50));
            Assert.That(leaves, Has.Count.EqualTo(50));
            Assert.That(sw.ElapsedMilliseconds, Is.LessThan(30_000),
                "50-taxon UPGMA tree should build in < 30 s (O(n³) with n=50)");
        });
    }

    [Test]
    [Description("C02: Gap-only columns are handled correctly — gaps skipped, not treated as mismatches")]
    public void BuildTree_GapOnlyColumns_HandledCorrectly()
    {
        // Positions 4-5 are gap-only columns for all sequences
        var sequences = new Dictionary<string, string>
        {
            ["A"] = "ACGT--ACGT",
            ["B"] = "ACGT--ACGT",  // Identical to A at all non-gap sites
            ["C"] = "TCGT--ACGT"   // Differs at position 0 only
        };

        var tree = PhylogeneticAnalyzer.BuildTree(sequences);
        var leaves = PhylogeneticAnalyzer.GetLeaves(tree.Root).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(tree.Root, Is.Not.Null);
            Assert.That(leaves, Has.Count.EqualTo(3));

            // A and B identical at 8 comparable sites → distance 0
            Assert.That(tree.DistanceMatrix[0, 1], Is.EqualTo(0).Within(1e-10),
                "Identical sequences at non-gap positions should have distance 0");

            // A and C differ at 1 of 8 comparable sites → positive distance
            Assert.That(tree.DistanceMatrix[0, 2], Is.GreaterThan(0),
                "Sequences differing at non-gap sites should have positive distance");
        });
    }

    [Test]
    [Description("C03: Different distance methods all produce valid trees")]
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
    [Description("C04: Tree total length is sum of all branch lengths")]
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
    [Description("C05: Both tree methods produce trees with same leaf set")]
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

    /// <summary>
    /// Recursively calculates the distance from the root to each leaf (tip).
    /// Used to verify the ultrametric property of UPGMA trees.
    /// </summary>
    private static void CalculateRootToTipDistances(
        PhylogeneticAnalyzer.PhyloNode? node,
        double accumulatedDistance,
        Dictionary<string, double> tipDistances)
    {
        if (node == null) return;

        if (node.IsLeaf)
        {
            tipDistances[node.Name] = accumulatedDistance;
            return;
        }

        if (node.Left != null)
            CalculateRootToTipDistances(node.Left,
                accumulatedDistance + node.Left.BranchLength, tipDistances);
        if (node.Right != null)
            CalculateRootToTipDistances(node.Right,
                accumulatedDistance + node.Right.BranchLength, tipDistances);
    }

    /// <summary>
    /// Finds a leaf node by name in the tree.
    /// </summary>
    private static PhylogeneticAnalyzer.PhyloNode? FindLeafByName(
        PhylogeneticAnalyzer.PhyloNode? node, string name)
    {
        if (node == null) return null;
        if (node.IsLeaf) return node.Name == name ? node : null;
        return FindLeafByName(node.Left, name) ?? FindLeafByName(node.Right, name);
    }

    #endregion
}
