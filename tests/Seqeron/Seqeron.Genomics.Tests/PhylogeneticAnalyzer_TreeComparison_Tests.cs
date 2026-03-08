using NUnit.Framework;
using Seqeron.Genomics;
using System.Collections.Generic;
using System.Linq;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Canonical tests for PhylogeneticAnalyzer tree comparison methods.
/// 
/// Test Unit: PHYLO-COMP-001
/// Methods:
///   - RobinsonFouldsDistance (topology comparison)
///   - FindMRCA (most recent common ancestor)
///   - PatristicDistance (tree path distance)
/// 
/// Evidence: Wikipedia (Robinson-Foulds metric, MRCA, Phylogenetic tree),
///           Robinson &amp; Foulds (1981)
/// </summary>
[TestFixture]
[Category("PHYLO-COMP-001")]
[Description("Tree Comparison: Robinson-Foulds, MRCA, Patristic Distance")]
public class PhylogeneticAnalyzer_TreeComparison_Tests
{
    #region Test Data Helpers

    /// <summary>
    /// Creates a simple tree: ((A,B),(C,D))
    /// </summary>
    private static PhylogeneticAnalyzer.PhyloNode CreateFourTaxaTree_ABCD()
    {
        var a = new PhylogeneticAnalyzer.PhyloNode("A") { BranchLength = 0.5 };
        a.Taxa = new List<string> { "A" };

        var b = new PhylogeneticAnalyzer.PhyloNode("B") { BranchLength = 0.5 };
        b.Taxa = new List<string> { "B" };

        var c = new PhylogeneticAnalyzer.PhyloNode("C") { BranchLength = 1.0 };
        c.Taxa = new List<string> { "C" };

        var d = new PhylogeneticAnalyzer.PhyloNode("D") { BranchLength = 1.0 };
        d.Taxa = new List<string> { "D" };

        var ab = new PhylogeneticAnalyzer.PhyloNode { Left = a, Right = b, BranchLength = 1.5 };
        ab.Taxa = new List<string> { "A", "B" };

        var cd = new PhylogeneticAnalyzer.PhyloNode { Left = c, Right = d, BranchLength = 2.0 };
        cd.Taxa = new List<string> { "C", "D" };

        var root = new PhylogeneticAnalyzer.PhyloNode { Left = ab, Right = cd };
        root.Taxa = new List<string> { "A", "B", "C", "D" };

        return root;
    }

    /// <summary>
    /// Creates a tree with different topology: ((A,C),(B,D))
    /// </summary>
    private static PhylogeneticAnalyzer.PhyloNode CreateFourTaxaTree_ACBD()
    {
        var a = new PhylogeneticAnalyzer.PhyloNode("A") { BranchLength = 0.5 };
        a.Taxa = new List<string> { "A" };

        var c = new PhylogeneticAnalyzer.PhyloNode("C") { BranchLength = 0.5 };
        c.Taxa = new List<string> { "C" };

        var b = new PhylogeneticAnalyzer.PhyloNode("B") { BranchLength = 1.0 };
        b.Taxa = new List<string> { "B" };

        var d = new PhylogeneticAnalyzer.PhyloNode("D") { BranchLength = 1.0 };
        d.Taxa = new List<string> { "D" };

        var ac = new PhylogeneticAnalyzer.PhyloNode { Left = a, Right = c, BranchLength = 1.5 };
        ac.Taxa = new List<string> { "A", "C" };

        var bd = new PhylogeneticAnalyzer.PhyloNode { Left = b, Right = d, BranchLength = 2.0 };
        bd.Taxa = new List<string> { "B", "D" };

        var root = new PhylogeneticAnalyzer.PhyloNode { Left = ac, Right = bd };
        root.Taxa = new List<string> { "A", "B", "C", "D" };

        return root;
    }

    /// <summary>
    /// Creates a simple two-taxa tree: (A,B)
    /// </summary>
    private static PhylogeneticAnalyzer.PhyloNode CreateTwoTaxaTree()
    {
        var a = new PhylogeneticAnalyzer.PhyloNode("A") { BranchLength = 1.0 };
        a.Taxa = new List<string> { "A" };

        var b = new PhylogeneticAnalyzer.PhyloNode("B") { BranchLength = 2.0 };
        b.Taxa = new List<string> { "B" };

        var root = new PhylogeneticAnalyzer.PhyloNode { Left = a, Right = b };
        root.Taxa = new List<string> { "A", "B" };

        return root;
    }

    /// <summary>
    /// Creates a three-taxa tree: ((A,B),C)
    /// </summary>
    private static PhylogeneticAnalyzer.PhyloNode CreateThreeTaxaTree()
    {
        var a = new PhylogeneticAnalyzer.PhyloNode("A") { BranchLength = 0.5 };
        a.Taxa = new List<string> { "A" };

        var b = new PhylogeneticAnalyzer.PhyloNode("B") { BranchLength = 0.5 };
        b.Taxa = new List<string> { "B" };

        var c = new PhylogeneticAnalyzer.PhyloNode("C") { BranchLength = 1.5 };
        c.Taxa = new List<string> { "C" };

        var ab = new PhylogeneticAnalyzer.PhyloNode { Left = a, Right = b, BranchLength = 1.0 };
        ab.Taxa = new List<string> { "A", "B" };

        var root = new PhylogeneticAnalyzer.PhyloNode { Left = ab, Right = c };
        root.Taxa = new List<string> { "A", "B", "C" };

        return root;
    }

    /// <summary>
    /// Creates a three-taxa tree with different topology: ((A,C),B)
    /// </summary>
    private static PhylogeneticAnalyzer.PhyloNode CreateThreeTaxaTree_ACB()
    {
        var a = new PhylogeneticAnalyzer.PhyloNode("A") { BranchLength = 0.5 };
        a.Taxa = new List<string> { "A" };

        var c = new PhylogeneticAnalyzer.PhyloNode("C") { BranchLength = 0.5 };
        c.Taxa = new List<string> { "C" };

        var b = new PhylogeneticAnalyzer.PhyloNode("B") { BranchLength = 1.5 };
        b.Taxa = new List<string> { "B" };

        var ac = new PhylogeneticAnalyzer.PhyloNode { Left = a, Right = c, BranchLength = 1.0 };
        ac.Taxa = new List<string> { "A", "C" };

        var root = new PhylogeneticAnalyzer.PhyloNode { Left = ac, Right = b };
        root.Taxa = new List<string> { "A", "B", "C" };

        return root;
    }

    #endregion

    #region Robinson-Foulds Distance Tests

    [Test]
    [Description("RF-M01: Identical trees have Robinson-Foulds distance = 0 (Wikipedia: RF metric)")]
    public void RobinsonFouldsDistance_IdenticalTrees_ReturnsZero()
    {
        // Arrange
        var tree1 = CreateFourTaxaTree_ABCD();
        var tree2 = CreateFourTaxaTree_ABCD();

        // Act
        int rfDistance = PhylogeneticAnalyzer.RobinsonFouldsDistance(tree1, tree2);

        // Assert
        Assert.That(rfDistance, Is.EqualTo(0), "Identical trees must have RF distance of 0");
    }

    [Test]
    [Description("RF-M02: RF distance is symmetric (Robinson & Foulds 1981: metric property)")]
    public void RobinsonFouldsDistance_IsSymmetric()
    {
        // Arrange
        var tree1 = CreateFourTaxaTree_ABCD();
        var tree2 = CreateFourTaxaTree_ACBD();

        // Act
        int rf12 = PhylogeneticAnalyzer.RobinsonFouldsDistance(tree1, tree2);
        int rf21 = PhylogeneticAnalyzer.RobinsonFouldsDistance(tree2, tree1);

        // Assert
        Assert.That(rf12, Is.EqualTo(rf21), "RF distance must be symmetric: RF(T1,T2) = RF(T2,T1)");
    }

    [Test]
    [Description("RF-M03: RF distance is non-negative — verified via exact values " +
                 "(Robinson & Foulds 1981: metric property)")]
    public void RobinsonFouldsDistance_IsNonNegative()
    {
        // Verify non-negativity through exact values on diverse tree pairs
        Assert.Multiple(() =>
        {
            // Two-taxa identical: only one rooted topology → RF = 0
            Assert.That(PhylogeneticAnalyzer.RobinsonFouldsDistance(
                CreateTwoTaxaTree(), CreateTwoTaxaTree()), Is.EqualTo(0),
                "Two-taxa identical trees: RF = 0");

            // Three-taxa identical: RF = 0
            Assert.That(PhylogeneticAnalyzer.RobinsonFouldsDistance(
                CreateThreeTaxaTree(), CreateThreeTaxaTree()), Is.EqualTo(0),
                "Three-taxa identical trees: RF = 0");

            // Three-taxa different: RF = 2 (max for n=3)
            Assert.That(PhylogeneticAnalyzer.RobinsonFouldsDistance(
                CreateThreeTaxaTree(), CreateThreeTaxaTree_ACB()), Is.EqualTo(2),
                "Three-taxa different topology: RF = 2");
        });
    }

    [Test]
    [Description("RF-M04: Different topologies have positive RF distance — exact value " +
                 "(Wikipedia: symmetric difference)")]
    public void RobinsonFouldsDistance_DifferentTopologies_ReturnsPositive()
    {
        // Arrange — three-taxa trees with different groupings:
        // Tree1: ((A,B),C) — clade {A,B} → canonical split "C"
        // Tree2: ((A,C),B) — clade {A,C} → canonical split "B"
        // Symmetric diff = |{"C"} - {"B"}| + |{"B"} - {"C"}| = 2
        var tree1 = CreateThreeTaxaTree();
        var tree2 = CreateThreeTaxaTree_ACB();

        // Act
        int rfDistance = PhylogeneticAnalyzer.RobinsonFouldsDistance(tree1, tree2);

        // Assert — exact value rather than just > 0
        Assert.That(rfDistance, Is.EqualTo(2),
            "((A,B),C) vs ((A,C),B): one non-trivial clade differs in each → RF = 2");
    }

    [Test]
    [Description("RF-M06: Four-taxa trees with maximally different topologies have RF = 4 " +
                 "(rooted RF: max = 2(n-2) = 4 for n=4; Wikipedia RF metric + dummy-leaf approach)")]
    public void RobinsonFouldsDistance_FourTaxaMaxDifference_ReturnsExact4()
    {
        // Arrange
        // Tree1: ((A,B),(C,D)) — clades: {A,B}, {C,D}
        // Tree2: ((A,C),(B,D)) — clades: {A,C}, {B,D}
        // No shared non-trivial clades → RF = 4 = 2(n-2) = max
        var tree1 = CreateFourTaxaTree_ABCD();
        var tree2 = CreateFourTaxaTree_ACBD();

        // Act
        int rfDistance = PhylogeneticAnalyzer.RobinsonFouldsDistance(tree1, tree2);

        // Assert — exact value from clade symmetric difference
        Assert.That(rfDistance, Is.EqualTo(4),
            "((A,B),(C,D)) vs ((A,C),(B,D)): no shared clades → RF = 4");
    }

    [Test]
    [Description("RF-M05: RF distance is always even — verified with exact values on multiple tree sizes " +
                 "(Wikipedia: symmetric difference of two sets)")]
    public void RobinsonFouldsDistance_IsEven()
    {
        // Verify evenness through exact values on diverse tree pairs
        Assert.Multiple(() =>
        {
            // Four-taxa: RF = 4 (even ✓)
            int rf4 = PhylogeneticAnalyzer.RobinsonFouldsDistance(
                CreateFourTaxaTree_ABCD(), CreateFourTaxaTree_ACBD());
            Assert.That(rf4, Is.EqualTo(4), "Four-taxa different: RF = 4");
            Assert.That(rf4 % 2, Is.EqualTo(0), "RF = 4 is even");

            // Three-taxa: RF = 2 (even ✓)
            int rf3 = PhylogeneticAnalyzer.RobinsonFouldsDistance(
                CreateThreeTaxaTree(), CreateThreeTaxaTree_ACB());
            Assert.That(rf3, Is.EqualTo(2), "Three-taxa different: RF = 2");
            Assert.That(rf3 % 2, Is.EqualTo(0), "RF = 2 is even");

            // Identical: RF = 0 (even ✓)
            int rf0 = PhylogeneticAnalyzer.RobinsonFouldsDistance(
                CreateFourTaxaTree_ABCD(), CreateFourTaxaTree_ABCD());
            Assert.That(rf0, Is.EqualTo(0), "Identical trees: RF = 0");
            Assert.That(rf0 % 2, Is.EqualTo(0), "RF = 0 is even");
        });
    }

    [Test]
    [Description("RF-S01: Three taxa trees with same topology have RF = 0")]
    public void RobinsonFouldsDistance_ThreeTaxaSameTopology_ReturnsZero()
    {
        // Arrange
        var tree1 = CreateThreeTaxaTree();
        var tree2 = CreateThreeTaxaTree();

        // Act
        int rfDistance = PhylogeneticAnalyzer.RobinsonFouldsDistance(tree1, tree2);

        // Assert
        Assert.That(rfDistance, Is.EqualTo(0));
    }

    [Test]
    [Description("RF-S02: Three-taxa trees with different topologies have RF = 2 " +
                 "(rooted RF: max = 2(n-2) = 2 for n=3; each tree has 1 non-trivial clade, none shared)")]
    public void RobinsonFouldsDistance_ThreeTaxaDifferentTopology_ReturnsExact2()
    {
        // Arrange
        // Tree1: ((A,B),C) — clade: {A,B}
        // Tree2: ((A,C),B) — clade: {A,C}
        // Symmetric diff = {A,B} + {A,C} = 2
        var tree1 = CreateThreeTaxaTree();
        var tree2 = CreateThreeTaxaTree_ACB();

        // Act
        int rfDistance = PhylogeneticAnalyzer.RobinsonFouldsDistance(tree1, tree2);

        // Assert
        Assert.That(rfDistance, Is.EqualTo(2),
            "((A,B),C) vs ((A,C),B): one clade differs in each → RF = 2");
    }

    [Test]
    [Description("RF-S03: RF distance is bounded by 2(n-2) for rooted binary trees " +
                 "(Wikipedia: max partitions per tree = n-2 non-trivial clades)")]
    public void RobinsonFouldsDistance_BoundedByMaximum()
    {
        // For n taxa in a rooted binary tree, max RF = 2(n-2)
        var tree1 = CreateFourTaxaTree_ABCD();
        var tree2 = CreateFourTaxaTree_ACBD();

        int n = 4;
        int maxRF = 2 * (n - 2); // = 4

        int rfDistance = PhylogeneticAnalyzer.RobinsonFouldsDistance(tree1, tree2);

        Assert.That(rfDistance, Is.LessThanOrEqualTo(maxRF),
            $"RF distance must be ≤ 2(n-2) = {maxRF} for rooted binary trees with {n} taxa");
    }

    [Test]
    [Description("RF integration: Tree built from sequences compared to itself returns 0")]
    public void RobinsonFouldsDistance_SequenceBuiltTree_SelfComparison_ReturnsZero()
    {
        // Arrange
        var sequences = new Dictionary<string, string>
        {
            ["A"] = "ACGT",
            ["B"] = "TCGT",
            ["C"] = "GCGT"
        };
        var tree = PhylogeneticAnalyzer.BuildTree(sequences, treeMethod: PhylogeneticAnalyzer.TreeMethod.UPGMA);

        // Act
        int rfDistance = PhylogeneticAnalyzer.RobinsonFouldsDistance(tree.Root, tree.Root);

        // Assert
        Assert.That(rfDistance, Is.EqualTo(0));
    }

    #endregion

    #region FindMRCA Tests

    [Test]
    [Description("MRCA-M01: Same taxon queried twice returns the taxon node itself (Wikipedia: MRCA definition)")]
    public void FindMRCA_SameTaxon_ReturnsTaxonItself()
    {
        // Arrange
        var tree = CreateFourTaxaTree_ABCD();

        // Act
        var mrca = PhylogeneticAnalyzer.FindMRCA(tree, "A", "A");

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(mrca, Is.Not.Null, "MRCA should not be null for existing taxon");
            Assert.That(mrca!.Name, Is.EqualTo("A"), "MRCA(A,A) should return node A");
            Assert.That(mrca.IsLeaf, Is.True, "A is a leaf node");
        });
    }

    [Test]
    [Description("MRCA-M02: Sibling taxa return their parent node (Wikipedia: MRCA is deepest common ancestor)")]
    public void FindMRCA_SiblingTaxa_ReturnsParent()
    {
        // Arrange - In ((A,B),(C,D)), A and B are siblings
        var tree = CreateFourTaxaTree_ABCD();

        // Act
        var mrca = PhylogeneticAnalyzer.FindMRCA(tree, "A", "B");

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(mrca, Is.Not.Null);
            Assert.That(mrca!.Taxa, Does.Contain("A"), "MRCA must contain taxon A");
            Assert.That(mrca.Taxa, Does.Contain("B"), "MRCA must contain taxon B");
            Assert.That(mrca.Taxa.Count, Is.EqualTo(2), "MRCA of siblings should only contain the two siblings");
        });
    }

    [Test]
    [Description("MRCA-M03: Distant taxa return deeper common ancestor (Wikipedia: MRCA definition)")]
    public void FindMRCA_DistantTaxa_ReturnsDeepAncestor()
    {
        // Arrange - In ((A,B),(C,D)), A and C require going to root
        var tree = CreateFourTaxaTree_ABCD();

        // Act
        var mrca = PhylogeneticAnalyzer.FindMRCA(tree, "A", "C");

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(mrca, Is.Not.Null);
            Assert.That(mrca!.Taxa, Does.Contain("A"));
            Assert.That(mrca.Taxa, Does.Contain("C"));
            Assert.That(mrca.Taxa.Count, Is.EqualTo(4), "MRCA of A and C should be root (all 4 taxa)");
        });
    }

    [Test]
    [Description("MRCA-M04: MRCA contains both queried taxa in its subtree (Wikipedia: MRCA properties)")]
    public void FindMRCA_ContainsBothTaxa()
    {
        // Arrange
        var tree = CreateThreeTaxaTree();

        // Act
        var mrca = PhylogeneticAnalyzer.FindMRCA(tree, "A", "B");

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(mrca, Is.Not.Null);
            Assert.That(mrca!.Taxa, Does.Contain("A"));
            Assert.That(mrca.Taxa, Does.Contain("B"));
        });
    }

    [Test]
    [Description("MRCA-M05: Null root returns null (Edge case)")]
    public void FindMRCA_NullRoot_ReturnsNull()
    {
        // Act
        var mrca = PhylogeneticAnalyzer.FindMRCA(null!, "A", "B");

        // Assert
        Assert.That(mrca, Is.Null);
    }

    [Test]
    [Description("MRCA invariant: MRCA is symmetric — MRCA(x,y) = MRCA(y,x)")]
    public void FindMRCA_IsSymmetric()
    {
        // Arrange
        var tree = CreateFourTaxaTree_ABCD();

        // Act
        var mrcaAB = PhylogeneticAnalyzer.FindMRCA(tree, "A", "B");
        var mrcaBA = PhylogeneticAnalyzer.FindMRCA(tree, "B", "A");

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(mrcaAB, Is.Not.Null);
            Assert.That(mrcaBA, Is.Not.Null);
            Assert.That(mrcaAB!.Taxa, Is.EquivalentTo(mrcaBA!.Taxa));
        });
    }

    [Test]
    [Description("MRCA with sequence-built tree finds common ancestor")]
    public void FindMRCA_SequenceBuiltTree_FindsCommonAncestor()
    {
        // Arrange
        var sequences = new Dictionary<string, string>
        {
            ["Human"] = "ACGTACGT",
            ["Chimp"] = "ACGTACGA",
            ["Mouse"] = "TCGTACGT"
        };
        var tree = PhylogeneticAnalyzer.BuildTree(sequences);

        // Act
        var mrca = PhylogeneticAnalyzer.FindMRCA(tree.Root, "Human", "Chimp");

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(mrca, Is.Not.Null);
            Assert.That(mrca!.Taxa, Does.Contain("Human"));
            Assert.That(mrca.Taxa, Does.Contain("Chimp"));
        });
    }

    [Test]
    [Description("MRCA-S01: MRCA of taxa from different subtrees is root (Property verification)")]
    public void FindMRCA_CrossCladeTaxa_ReturnsRoot()
    {
        // Arrange — In ((A,B),(C,D)), any pair spanning both clades has root as MRCA
        var tree = CreateFourTaxaTree_ABCD();

        // Act & Assert — all cross-clade pairs must return root (4 taxa)
        Assert.Multiple(() =>
        {
            var mrcaAC = PhylogeneticAnalyzer.FindMRCA(tree, "A", "C");
            var mrcaAD = PhylogeneticAnalyzer.FindMRCA(tree, "A", "D");
            var mrcaBD = PhylogeneticAnalyzer.FindMRCA(tree, "B", "D");
            var mrcaBC = PhylogeneticAnalyzer.FindMRCA(tree, "B", "C");

            Assert.That(mrcaAC, Is.Not.Null);
            Assert.That(mrcaAC!.Taxa.Count, Is.EqualTo(4), "MRCA(A,C) is root");
            Assert.That(mrcaAD!.Taxa.Count, Is.EqualTo(4), "MRCA(A,D) is root");
            Assert.That(mrcaBD!.Taxa.Count, Is.EqualTo(4), "MRCA(B,D) is root");
            Assert.That(mrcaBC!.Taxa.Count, Is.EqualTo(4), "MRCA(B,C) is root");
        });
    }

    [Test]
    [Description("MRCA edge case: Non-existent taxon returns null " +
                 "(Evidence doc: 'Taxon not in tree: Returns null')")]
    public void FindMRCA_NonExistentTaxon_ReturnsNull()
    {
        // Arrange
        var tree = CreateFourTaxaTree_ABCD();

        // Act & Assert
        Assert.Multiple(() =>
        {
            // One existing + one non-existent
            Assert.That(PhylogeneticAnalyzer.FindMRCA(tree, "A", "NonExistent"), Is.Null,
                "MRCA with one non-existent taxon must return null");

            // Both non-existent
            Assert.That(PhylogeneticAnalyzer.FindMRCA(tree, "X", "Y"), Is.Null,
                "MRCA with both non-existent taxa must return null");
        });
    }

    #endregion

    #region Patristic Distance Tests

    [Test]
    [Description("PD-M01: Same taxon has patristic distance 0 (Definition: no path to traverse)")]
    public void PatristicDistance_SameTaxon_ReturnsZero()
    {
        // Arrange
        var tree = CreateFourTaxaTree_ABCD();

        // Act
        double distance = PhylogeneticAnalyzer.PatristicDistance(tree, "A", "A");

        // Assert
        Assert.That(distance, Is.EqualTo(0.0), "Distance from taxon to itself must be 0");
    }

    [Test]
    [Description("PD-M02: Sibling distance equals sum of their branch lengths (Definition)")]
    public void PatristicDistance_Siblings_ReturnsSumOfBranchLengths()
    {
        // Arrange - A and B are siblings with BL=0.5 each
        var tree = CreateTwoTaxaTree();

        // Act
        double distance = PhylogeneticAnalyzer.PatristicDistance(tree, "A", "B");

        // Assert
        // A (BL=1.0) + B (BL=2.0) = 3.0
        Assert.That(distance, Is.EqualTo(3.0).Within(0.0001),
            "Sibling distance should be sum of branch lengths");
    }

    [Test]
    [Description("PD-M03: Patristic distance is symmetric (Metric property)")]
    public void PatristicDistance_IsSymmetric()
    {
        // Arrange
        var tree = CreateFourTaxaTree_ABCD();

        // Act
        double distAB = PhylogeneticAnalyzer.PatristicDistance(tree, "A", "B");
        double distBA = PhylogeneticAnalyzer.PatristicDistance(tree, "B", "A");

        // Assert
        Assert.That(distAB, Is.EqualTo(distBA), "Patristic distance must be symmetric");
    }

    [Test]
    [Description("PD-M04: Non-existent taxon returns NaN (Edge case)")]
    public void PatristicDistance_NonExistentTaxon_ReturnsNaN()
    {
        // Arrange
        var tree = CreateFourTaxaTree_ABCD();

        // Act
        double distance = PhylogeneticAnalyzer.PatristicDistance(tree, "A", "NonExistent");

        // Assert
        Assert.That(double.IsNaN(distance), Is.True,
            "Distance with non-existent taxon should return NaN");
    }

    [Test]
    [Description("PD-M05: Patristic distance is non-negative — verified via exact values on three-taxa tree " +
                 "(Metric property)")]
    public void PatristicDistance_IsNonNegative()
    {
        // Arrange — three-taxa tree: ((A(0.5),B(0.5))(1.0),C(1.5))
        // Uses different tree from PD-S01 (which uses four-taxa tree)
        var tree = CreateThreeTaxaTree();

        // Act & Assert — exact values prove non-negativity
        Assert.Multiple(() =>
        {
            // A→B: 0.5 + 0.5 = 1.0 (siblings under AB)
            Assert.That(PhylogeneticAnalyzer.PatristicDistance(tree, "A", "B"),
                Is.EqualTo(1.0).Within(1e-10), "PD(A,B) = 0.5 + 0.5 = 1.0");

            // A→C: 0.5 + 1.0 + 1.5 = 3.0 (path through root)
            Assert.That(PhylogeneticAnalyzer.PatristicDistance(tree, "A", "C"),
                Is.EqualTo(3.0).Within(1e-10), "PD(A,C) = 0.5 + 1.0 + 1.5 = 3.0");

            // B→C: 0.5 + 1.0 + 1.5 = 3.0 (path through root)
            Assert.That(PhylogeneticAnalyzer.PatristicDistance(tree, "B", "C"),
                Is.EqualTo(3.0).Within(1e-10), "PD(B,C) = 0.5 + 1.0 + 1.5 = 3.0");
        });
    }

    [Test]
    [Description("PD-S01: All documented patristic distances for four-taxa tree " +
                 "(TestSpec tree: A-B=1.0, A-C=5.0, C-D=2.0)")]
    public void PatristicDistance_FourTaxaTree_AllDocumentedValues()
    {
        // Arrange — four-taxa tree with known branch lengths:
        //     root (BL=0)
        //    /          \
        //   AB (BL=1.5)  CD (BL=2.0)
        //  / \          / \
        // A   B        C   D
        // (0.5)(0.5)  (1.0)(1.0)
        var tree = CreateFourTaxaTree_ABCD();

        // Act & Assert — exact values from test spec
        Assert.Multiple(() =>
        {
            // A→B: 0.5 + 0.5 = 1.0 (siblings under AB)
            Assert.That(PhylogeneticAnalyzer.PatristicDistance(tree, "A", "B"),
                Is.EqualTo(1.0).Within(1e-10), "PD(A,B) = 0.5 + 0.5 = 1.0");

            // A→C: 0.5 + 1.5 + 2.0 + 1.0 = 5.0 (path through root)
            Assert.That(PhylogeneticAnalyzer.PatristicDistance(tree, "A", "C"),
                Is.EqualTo(5.0).Within(1e-10), "PD(A,C) = 0.5 + 1.5 + 2.0 + 1.0 = 5.0");

            // C→D: 1.0 + 1.0 = 2.0 (siblings under CD)
            Assert.That(PhylogeneticAnalyzer.PatristicDistance(tree, "C", "D"),
                Is.EqualTo(2.0).Within(1e-10), "PD(C,D) = 1.0 + 1.0 = 2.0");
        });
    }

    [Test]
    [Description("PD-S02: Patristic distance satisfies triangle inequality " +
                 "(Wikipedia LCA: d(v,w) is a metric on tree; metric axiom)")]
    public void PatristicDistance_SatisfiesTriangleInequality()
    {
        // Arrange
        var tree = CreateFourTaxaTree_ABCD();

        // Act
        double distAB = PhylogeneticAnalyzer.PatristicDistance(tree, "A", "B");
        double distBC = PhylogeneticAnalyzer.PatristicDistance(tree, "B", "C");
        double distAC = PhylogeneticAnalyzer.PatristicDistance(tree, "A", "C");
        double distAD = PhylogeneticAnalyzer.PatristicDistance(tree, "A", "D");
        double distBD = PhylogeneticAnalyzer.PatristicDistance(tree, "B", "D");
        double distCD = PhylogeneticAnalyzer.PatristicDistance(tree, "C", "D");

        // Assert — triangle inequality: PD(x,z) ≤ PD(x,y) + PD(y,z)
        Assert.Multiple(() =>
        {
            Assert.That(distAC, Is.LessThanOrEqualTo(distAB + distBC + 1e-10),
                "PD(A,C) ≤ PD(A,B) + PD(B,C)");
            Assert.That(distAD, Is.LessThanOrEqualTo(distAB + distBD + 1e-10),
                "PD(A,D) ≤ PD(A,B) + PD(B,D)");
            Assert.That(distBD, Is.LessThanOrEqualTo(distBC + distCD + 1e-10),
                "PD(B,D) ≤ PD(B,C) + PD(C,D)");
            Assert.That(distAC, Is.LessThanOrEqualTo(distAD + distCD + 1e-10),
                "PD(A,C) ≤ PD(A,D) + PD(C,D)");
        });
    }

    [Test]
    [Description("Patristic distance with sequence-built tree")]
    public void PatristicDistance_SequenceBuiltTree_CalculatesDistance()
    {
        // Arrange
        var sequences = new Dictionary<string, string>
        {
            ["A"] = "ACGTACGT",
            ["B"] = "ACGTACGA",
            ["C"] = "TCGTACGT"
        };
        var tree = PhylogeneticAnalyzer.BuildTree(sequences);

        // Act
        double distance = PhylogeneticAnalyzer.PatristicDistance(tree.Root, "A", "B");

        // Assert
        Assert.That(distance, Is.GreaterThanOrEqualTo(0));
    }

    #endregion

    #region Helper Methods (Smoke Tests)

    [Test]
    [Description("GetLeaves returns all leaf nodes")]
    public void GetLeaves_ReturnsAllLeafNodes()
    {
        // Arrange
        var tree = CreateFourTaxaTree_ABCD();

        // Act
        var leaves = PhylogeneticAnalyzer.GetLeaves(tree).ToList();

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(leaves, Has.Count.EqualTo(4));
            Assert.That(leaves.All(l => l.IsLeaf), Is.True);
            Assert.That(leaves.Select(l => l.Name), Is.EquivalentTo(new[] { "A", "B", "C", "D" }));
        });
    }

    [Test]
    [Description("GetLeaves with null root returns empty")]
    public void GetLeaves_NullRoot_ReturnsEmpty()
    {
        // Act
        var leaves = PhylogeneticAnalyzer.GetLeaves(null!).ToList();

        // Assert
        Assert.That(leaves, Is.Empty);
    }

    [Test]
    [Description("CalculateTreeLength sums all branch lengths")]
    public void CalculateTreeLength_SumsAllBranches()
    {
        // Arrange
        var tree = CreateTwoTaxaTree();
        // A (BL=1.0) + B (BL=2.0) = 3.0

        // Act
        double length = PhylogeneticAnalyzer.CalculateTreeLength(tree);

        // Assert
        Assert.That(length, Is.EqualTo(3.0).Within(0.0001));
    }

    [Test]
    [Description("CalculateTreeLength with null root returns zero")]
    public void CalculateTreeLength_NullRoot_ReturnsZero()
    {
        // Act
        double length = PhylogeneticAnalyzer.CalculateTreeLength(null!);

        // Assert
        Assert.That(length, Is.EqualTo(0));
    }

    [Test]
    [Description("GetTreeDepth returns correct depth")]
    public void GetTreeDepth_ReturnsCorrectDepth()
    {
        // Arrange - ((A,B),(C,D)) has depth 2
        var tree = CreateFourTaxaTree_ABCD();

        // Act
        int depth = PhylogeneticAnalyzer.GetTreeDepth(tree);

        // Assert
        Assert.That(depth, Is.EqualTo(2));
    }

    #endregion

    #region Bootstrap Tests (Different Scope - Kept for Coverage)

    [Test]
    [Description("Bootstrap returns supports in valid range")]
    public void Bootstrap_ReturnsSupportsInValidRange()
    {
        // Arrange
        var sequences = new Dictionary<string, string>
        {
            ["A"] = "ACGTACGTAC",
            ["B"] = "ACGTACGTAC",
            ["C"] = "TCGTACGTAC",
            ["D"] = "TCGTACGTAC"
        };

        // Act
        var supports = PhylogeneticAnalyzer.Bootstrap(sequences, replicates: 10);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(supports, Is.Not.Empty);
            Assert.That(supports.Values, Has.All.InRange(0.0, 1.0));
        });
    }

    [Test]
    [Description("Bootstrap with distinct groups shows high support")]
    public void Bootstrap_DistinctGroups_HighSupport()
    {
        // Arrange - Groups with very different sequences should have high bootstrap support
        var sequences = new Dictionary<string, string>
        {
            ["A"] = "AAAAAAAAAA",
            ["B"] = "AAAAAAAAAC",
            ["C"] = "CCCCCCCCCC",
            ["D"] = "CCCCCCCCCA"
        };

        // Act
        var supports = PhylogeneticAnalyzer.Bootstrap(sequences, replicates: 50);

        // Assert - At least some splits should have high support
        Assert.That(supports.Values.Any(v => v >= 0.5), Is.True);
    }

    #endregion
}
