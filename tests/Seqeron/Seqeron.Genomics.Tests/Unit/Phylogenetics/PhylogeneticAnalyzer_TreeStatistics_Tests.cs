// PHYLO-STATS-001 — Tree Statistics (leaves, total tree length, tree height/depth)
// Evidence: docs/Evidence/PHYLO-STATS-001-Evidence.md
// TestSpec: tests/TestSpecs/PHYLO-STATS-001.md
// Source: Wikipedia "Tree (graph theory)" / "Tree (abstract data type)" (leaf, height definitions);
//         Biopython Bio.Phylo.BaseTree (get_terminals, total_branch_length);
//         DendroPy Tree.length() ("sum of edge lengths of self");
//         Rzhetsky A, Nei M (1992) via Wikipedia "Minimum evolution" (total branch length).

using NUnit.Framework;
using Seqeron.Genomics.Phylogenetics;
using System.Linq;

namespace Seqeron.Genomics.Tests.Unit.Phylogenetics;

/// <summary>
/// Tests for PHYLO-STATS-001: tree statistics over a rooted binary phylogenetic tree.
/// Covers <see cref="PhylogeneticAnalyzer.GetLeaves"/>, <see cref="PhylogeneticAnalyzer.CalculateTreeLength"/>,
/// and <see cref="PhylogeneticAnalyzer.GetTreeDepth"/>. All expected values are derived from the
/// cited tree-statistics definitions (Evidence §Test Datasets), not from running the implementation.
/// </summary>
[TestFixture]
[Category("PHYLO-STATS-001")]
public class PhylogeneticAnalyzer_TreeStatistics_Tests
{
    // Balanced 4-taxon tree ((A:1,B:1):1,(C:1,D:1):1): Σ branch lengths = 6, height = 2 edges.
    private static PhylogeneticAnalyzer.PhyloNode BalancedTree()
        => PhylogeneticAnalyzer.ParseNewick("((A:1,B:1):1,(C:1,D:1):1);");

    // Caterpillar (ladder) tree (A:1,(B:1,(C:1,D:1):0.5):0.5): Σ branch lengths = 5, height = 3 edges.
    private static PhylogeneticAnalyzer.PhyloNode CaterpillarTree()
        => PhylogeneticAnalyzer.ParseNewick("(A:1,(B:1,(C:1,D:1):0.5):0.5);");

    #region GetLeaves (M1, M2, M3, C1)

    // M1 — A leaf is "a vertex with no children" (Tree (graph theory)); Biopython get_terminals
    // returns "all of this tree's terminal (leaf) nodes". The balanced tree has exactly A,B,C,D.
    [Test]
    [Description("M1: GetLeaves returns exactly the terminal nodes (no children) of a balanced tree.")]
    public void GetLeaves_BalancedTree_ReturnsAllTerminalNodesInOrder()
    {
        var root = BalancedTree();

        var leaves = PhylogeneticAnalyzer.GetLeaves(root).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(leaves.Select(l => l.Name), Is.EqualTo(new[] { "A", "B", "C", "D" }),
                "Leaf set of ((A,B),(C,D)) must be exactly {A,B,C,D} in pre-order (Evidence Dataset 1).");
            Assert.That(leaves.All(l => l.IsLeaf), Is.True,
                "Every node returned by GetLeaves must have no children (INV-01: leaf definition).");
        });
    }

    // M2 — count of terminal nodes equals N (Biopython count_terminals); N=4 for the balanced tree.
    [Test]
    [Description("M2: leaf count of a 4-taxon tree is exactly 4 (count_terminals semantics).")]
    public void GetLeaves_FourLeafTree_CountEqualsFour()
    {
        var root = BalancedTree();

        int count = PhylogeneticAnalyzer.GetLeaves(root).Count();

        Assert.That(count, Is.EqualTo(4),
            "A fully-bifurcating 4-taxon tree has exactly 4 leaves (INV-02).");
    }

    // M3 — a single node that is both root and leaf is itself the only leaf (graph-theory leaf def).
    [Test]
    [Description("M3: GetLeaves on a lone leaf returns that node itself.")]
    public void GetLeaves_SingleLeaf_ReturnsThatNode()
    {
        var leaf = new PhylogeneticAnalyzer.PhyloNode("X");

        var leaves = PhylogeneticAnalyzer.GetLeaves(leaf).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(leaves, Has.Count.EqualTo(1), "A lone node is a single leaf.");
            Assert.That(leaves[0].Name, Is.EqualTo("X"), "The returned leaf must be the node itself.");
        });
    }

    // C1 — leaves are produced in left-to-right pre-order traversal order.
    [Test]
    [Description("C1: GetLeaves yields leaves in left-to-right pre-order.")]
    public void GetLeaves_BalancedTree_PreOrderTraversalOrder()
    {
        var root = BalancedTree();

        var order = PhylogeneticAnalyzer.GetLeaves(root).Select(l => l.Name).ToArray();

        Assert.That(order, Is.EqualTo(new[] { "A", "B", "C", "D" }),
            "Pre-order traversal of ((A,B),(C,D)) visits A,B,C,D in that order.");
    }

    #endregion

    #region CalculateTreeLength (M4, M5, M6, S1)

    // M4 — total tree length = "sum of edge lengths" (DendroPy) / "sum of all branch lengths"
    // (Biopython). Balanced tree: 4 leaf edges (1 each) + 2 internal edges (1 each) = 6.
    [Test]
    [Description("M4: total tree length of ((A:1,B:1):1,(C:1,D:1):1) is the sum of all branch lengths = 6.")]
    public void CalculateTreeLength_BalancedTree_SumsAllBranchLengths()
    {
        var root = BalancedTree();

        double length = PhylogeneticAnalyzer.CalculateTreeLength(root);

        Assert.That(length, Is.EqualTo(6.0).Within(1e-10),
            "Σ branch lengths = 1+1+1+1 (leaves) + 1+1 (internal) = 6 (INV-03; DendroPy/Biopython).");
    }

    // M5 — caterpillar tree: 4 leaf edges (1 each) + 2 internal edges (0.5 each) = 5.
    [Test]
    [Description("M5: total tree length of (A:1,(B:1,(C:1,D:1):0.5):0.5) = 5.")]
    public void CalculateTreeLength_CaterpillarTree_SumsAllBranchLengths()
    {
        var root = CaterpillarTree();

        double length = PhylogeneticAnalyzer.CalculateTreeLength(root);

        Assert.That(length, Is.EqualTo(5.0).Within(1e-10),
            "Σ branch lengths = 1+1+1+1 (leaves) + 0.5+0.5 (internal) = 5 (INV-03).");
    }

    // M6 — the length sum includes the branch length of the root node itself.
    [Test]
    [Description("M6: CalculateTreeLength includes the root node's own branch length.")]
    public void CalculateTreeLength_IncludesRootBranchLength()
    {
        // Subtree (C:1,D:1):0.5 — as a standalone tree its length must include the 0.5 root edge.
        var root = PhylogeneticAnalyzer.ParseNewick("(C:1,D:1):0.5;");

        double length = PhylogeneticAnalyzer.CalculateTreeLength(root);

        Assert.That(length, Is.EqualTo(2.5).Within(1e-10),
            "Σ = 1 (C) + 1 (D) + 0.5 (root edge) = 2.5; the root's own branch length is counted (INV-03).");
    }

    // S1 — DendroPy: edges with no length defined are treated as length 0. Default BranchLength is 0.
    [Test]
    [Description("S1: a tree whose branch lengths are all default (0) has total length 0.")]
    public void CalculateTreeLength_DefaultBranchLengths_ReturnsZero()
    {
        var root = new PhylogeneticAnalyzer.PhyloNode
        {
            Left = new PhylogeneticAnalyzer.PhyloNode("A"),
            Right = new PhylogeneticAnalyzer.PhyloNode("B")
        };

        double length = PhylogeneticAnalyzer.CalculateTreeLength(root);

        Assert.That(length, Is.EqualTo(0.0).Within(1e-10),
            "Unset branch lengths default to 0 and sum to 0 (DendroPy: missing length -> 0).");
    }

    #endregion

    #region GetTreeDepth (M7, M8, M9, S2)

    // M7 — height = longest root-to-leaf path in edges. Balanced tree: root->internal->leaf = 2.
    [Test]
    [Description("M7: tree height of ((A,B),(C,D)) is 2 edges to the deepest leaf.")]
    public void GetTreeDepth_BalancedTree_ReturnsEdgeCountToDeepestLeaf()
    {
        var root = BalancedTree();

        int height = PhylogeneticAnalyzer.GetTreeDepth(root);

        Assert.That(height, Is.EqualTo(2),
            "Longest root->leaf path is root->(internal)->leaf = 2 edges (INV-05; height definition).");
    }

    // M8 — caterpillar tree (A,(B,(C,D))): deepest leaves C,D at root->.->.->leaf = 3 edges.
    [Test]
    [Description("M8: tree height of (A,(B,(C,D))) is 3 edges.")]
    public void GetTreeDepth_CaterpillarTree_ReturnsThree()
    {
        var root = CaterpillarTree();

        int height = PhylogeneticAnalyzer.GetTreeDepth(root);

        Assert.That(height, Is.EqualTo(3),
            "Deepest leaves C/D are 3 edges below the root (INV-05).");
    }

    // M9 — a single node (both root and leaf) has height 0.
    [Test]
    [Description("M9: a lone leaf node has height 0.")]
    public void GetTreeDepth_SingleLeaf_ReturnsZero()
    {
        var leaf = new PhylogeneticAnalyzer.PhyloNode("X");

        int height = PhylogeneticAnalyzer.GetTreeDepth(leaf);

        Assert.That(height, Is.EqualTo(0),
            "A tree with a single node has height 0 (Tree (graph theory) / ADT).");
    }

    // S2 — a two-leaf tree (A,B) has one edge from root to each leaf: height 1.
    [Test]
    [Description("S2: tree height of (A,B) is 1 edge.")]
    public void GetTreeDepth_TwoLeafTree_ReturnsOne()
    {
        var root = PhylogeneticAnalyzer.ParseNewick("(A:1,B:1);");

        int height = PhylogeneticAnalyzer.GetTreeDepth(root);

        Assert.That(height, Is.EqualTo(1),
            "Root to either leaf is a single edge (INV-05).");
    }

    #endregion

    #region Null / empty tree (M10)

    // M10 — null root maps to the empty-tree convention: no leaves, length 0, height -1
    // ("Conventionally, an empty tree ... has depth and height -1" — Tree (graph theory)).
    [Test]
    [Description("M10: null root yields no leaves, length 0, and height -1 (empty-tree convention).")]
    public void TreeStatistics_NullRoot_ReturnEmptyZeroAndMinusOne()
    {
        Assert.Multiple(() =>
        {
            Assert.That(PhylogeneticAnalyzer.GetLeaves(null!), Is.Empty,
                "An empty tree has no leaves (INV-06).");
            Assert.That(PhylogeneticAnalyzer.CalculateTreeLength(null!), Is.EqualTo(0.0).Within(1e-10),
                "An empty tree has total length 0 (no edges) (INV-06).");
            Assert.That(PhylogeneticAnalyzer.GetTreeDepth(null!), Is.EqualTo(-1),
                "An empty tree has height -1 by the cited convention (INV-06).");
        });
    }

    #endregion
}
