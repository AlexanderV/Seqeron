using FsCheck;
using FsCheck.Fluent;
using FsCheck.NUnit;

namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Property-based tests for phylogenetic analysis.
/// Verifies distance matrix and Newick I/O invariants.
///
/// Test Units: PHYLO-DIST-001, PHYLO-NEWICK-001, PHYLO-TREE-001, PHYLO-COMP-001
/// </summary>
[TestFixture]
[Category("Property")]
[Category("Phylogenetics")]
public class PhylogeneticProperties
{
    private static Arbitrary<string> DnaArbitrary(int len = 30) =>
        Gen.Elements('A', 'C', 'G', 'T')
            .ArrayOf()
            .Where(a => a.Length >= len)
            .Select(a => new string(a, 0, len))
            .ToArbitrary();

    #region PHYLO-DIST-001: S: d(a,b)=d(b,a); I: d(x,x)=0; R: d ≥ 0; triangle inequality

    /// <summary>
    /// INV-1: Pairwise distance is symmetric: d(A,B) == d(B,A).
    /// Evidence: All supported distance methods count differences between positions,
    /// which is symmetric by definition.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property PairwiseDistance_IsSymmetric()
    {
        return Prop.ForAll(DnaArbitrary(20), DnaArbitrary(20), (s1, s2) =>
        {
            double dAB = PhylogeneticAnalyzer.CalculatePairwiseDistance(s1, s2,
                PhylogeneticAnalyzer.DistanceMethod.PDistance);
            double dBA = PhylogeneticAnalyzer.CalculatePairwiseDistance(s2, s1,
                PhylogeneticAnalyzer.DistanceMethod.PDistance);
            return (Math.Abs(dAB - dBA) < 1e-10)
                .Label($"d(A,B)={dAB} ≠ d(B,A)={dBA}");
        });
    }

    /// <summary>
    /// INV-2: Self-distance is zero: d(X,X) == 0.
    /// Evidence: A sequence compared to itself has zero differences.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property PairwiseDistance_SelfIsZero()
    {
        return Prop.ForAll(DnaArbitrary(20), seq =>
        {
            double d = PhylogeneticAnalyzer.CalculatePairwiseDistance(seq, seq,
                PhylogeneticAnalyzer.DistanceMethod.PDistance);
            return (d == 0.0)
                .Label($"d(X,X)={d}, expected 0");
        });
    }

    /// <summary>
    /// INV-3: Pairwise distance is non-negative: d ≥ 0.
    /// Evidence: Distance is a count of differences divided by total sites, always ≥ 0.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property PairwiseDistance_NonNegative()
    {
        return Prop.ForAll(DnaArbitrary(20), DnaArbitrary(20), (s1, s2) =>
        {
            double d = PhylogeneticAnalyzer.CalculatePairwiseDistance(s1, s2,
                PhylogeneticAnalyzer.DistanceMethod.PDistance);
            return (d >= 0.0)
                .Label($"Distance={d} must be ≥ 0");
        });
    }

    /// <summary>
    /// INV-4: Triangle inequality holds: d(A,C) ≤ d(A,B) + d(B,C).
    /// Evidence: p-distance (proportion of differences) satisfies metric axioms
    /// because character mismatches at each position are independent Bernoulli events.
    /// </summary>
    [Test]
    [Category("Property")]
    public void PairwiseDistance_TriangleInequality()
    {
        string a = "ACGTACGTACGTACGTACGT";
        string b = "ACGTACGTAAGTACGTACGT";
        string c = "ACGTACATACGTACATACGT";

        double dAB = PhylogeneticAnalyzer.CalculatePairwiseDistance(a, b,
            PhylogeneticAnalyzer.DistanceMethod.PDistance);
        double dBC = PhylogeneticAnalyzer.CalculatePairwiseDistance(b, c,
            PhylogeneticAnalyzer.DistanceMethod.PDistance);
        double dAC = PhylogeneticAnalyzer.CalculatePairwiseDistance(a, c,
            PhylogeneticAnalyzer.DistanceMethod.PDistance);

        Assert.That(dAC, Is.LessThanOrEqualTo(dAB + dBC + 1e-10),
            $"Triangle inequality violated: d(A,C)={dAC} > d(A,B)+d(B,C)={dAB + dBC}");
    }

    /// <summary>
    /// INV-5: Distance matrix is symmetric: M[i,j] == M[j,i].
    /// Evidence: CalculateDistanceMatrix explicitly sets matrix[i,j] = matrix[j,i].
    /// </summary>
    [Test]
    [Category("Property")]
    public void DistanceMatrix_IsSymmetric()
    {
        var seqs = new List<string>
        {
            "ACGTACGTACGTACGTACGT",
            "ACGTACGTAAGTACGTACGT",
            "ACGTACATACGTACATACGT"
        };
        var matrix = PhylogeneticAnalyzer.CalculateDistanceMatrix(seqs);

        for (int i = 0; i < seqs.Count; i++)
            for (int j = i + 1; j < seqs.Count; j++)
                Assert.That(matrix[i, j], Is.EqualTo(matrix[j, i]).Within(1e-10),
                    $"Matrix[{i},{j}]={matrix[i, j]} ≠ Matrix[{j},{i}]={matrix[j, i]}");
    }

    /// <summary>
    /// INV-6: Distance matrix diagonal is zero: M[i,i] == 0.
    /// Evidence: Distance from a sequence to itself is always zero.
    /// </summary>
    [Test]
    [Category("Property")]
    public void DistanceMatrix_DiagonalIsZero()
    {
        var seqs = new List<string>
        {
            "ACGTACGTACGTACGTACGT",
            "ACGTACGTAAGTACGTACGT",
            "ACGTACATACGTACATACGT"
        };
        var matrix = PhylogeneticAnalyzer.CalculateDistanceMatrix(seqs);

        for (int i = 0; i < seqs.Count; i++)
            Assert.That(matrix[i, i], Is.EqualTo(0.0),
                $"Matrix[{i},{i}]={matrix[i, i]}, expected 0");
    }

    /// <summary>
    /// INV-7: Pairwise distance calculation is deterministic.
    /// Evidence: CalculatePairwiseDistance is a pure function.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property PairwiseDistance_IsDeterministic()
    {
        return Prop.ForAll(DnaArbitrary(20), DnaArbitrary(20), (s1, s2) =>
        {
            double d1 = PhylogeneticAnalyzer.CalculatePairwiseDistance(s1, s2,
                PhylogeneticAnalyzer.DistanceMethod.JukesCantor);
            double d2 = PhylogeneticAnalyzer.CalculatePairwiseDistance(s1, s2,
                PhylogeneticAnalyzer.DistanceMethod.JukesCantor);
            return (d1 == d2)
                .Label("CalculatePairwiseDistance must be deterministic");
        });
    }

    #endregion

    #region PHYLO-NEWICK-001 / PHYLO-TREE-001 / PHYLO-COMP-001

    /// <summary>
    /// Newick round-trip: ToNewick → ParseNewick preserves leaf names.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Newick_RoundTrip_PreservesLeafNames()
    {
        var seqs = new Dictionary<string, string>
        {
            ["TaxonA"] = "ACGTACGTACGTACGTACGT",
            ["TaxonB"] = "ACGTACGTAAGTACGTACGT",
            ["TaxonC"] = "ACGTACATACGTACATACGT"
        };
        var tree = PhylogeneticAnalyzer.BuildTree(seqs);
        string newick = PhylogeneticAnalyzer.ToNewick(tree.Root);
        var parsed = PhylogeneticAnalyzer.ParseNewick(newick);
        var leaves = PhylogeneticAnalyzer.GetLeaves(parsed).Select(n => n.Name).OrderBy(n => n).ToList();
        var expected = seqs.Keys.OrderBy(k => k).ToList();

        Assert.That(leaves, Is.EqualTo(expected));
    }

    /// <summary>
    /// Tree length is non-negative.
    /// </summary>
    [Test]
    [Category("Property")]
    public void TreeLength_IsNonNegative()
    {
        var seqs = new Dictionary<string, string>
        {
            ["A"] = "ACGTACGTACGTACGTACGT",
            ["B"] = "ACGTACGTAAGTACGTACGT",
            ["C"] = "ACGTACATACGTACATACGT"
        };
        var tree = PhylogeneticAnalyzer.BuildTree(seqs);
        double length = PhylogeneticAnalyzer.CalculateTreeLength(tree.Root);
        Assert.That(length, Is.GreaterThanOrEqualTo(0.0));
    }

    /// <summary>
    /// Robinson-Foulds distance of a tree to itself is zero.
    /// </summary>
    [Test]
    [Category("Property")]
    public void RobinsonFoulds_SelfDistance_IsZero()
    {
        var seqs = new Dictionary<string, string>
        {
            ["A"] = "ACGTACGTACGTACGTACGT",
            ["B"] = "ACGTACGTAAGTACGTACGT",
            ["C"] = "ACGTACATACGTACATACGT",
            ["D"] = "TCGTACATACGTACATACGT"
        };
        var tree = PhylogeneticAnalyzer.BuildTree(seqs);
        int rf = PhylogeneticAnalyzer.RobinsonFouldsDistance(tree.Root, tree.Root);
        Assert.That(rf, Is.EqualTo(0));
    }

    /// <summary>
    /// Robinson-Foulds distance is symmetric.
    /// </summary>
    [Test]
    [Category("Property")]
    public void RobinsonFoulds_IsSymmetric()
    {
        var seqs = new Dictionary<string, string>
        {
            ["A"] = "ACGTACGTACGTACGTACGT",
            ["B"] = "ACGTACGTAAGTACGTACGT",
            ["C"] = "ACGTACATACGTACATACGT",
            ["D"] = "TCGTACATACGTACATACGT"
        };
        var tree1 = PhylogeneticAnalyzer.BuildTree(seqs, treeMethod: PhylogeneticAnalyzer.TreeMethod.UPGMA);
        var tree2 = PhylogeneticAnalyzer.BuildTree(seqs, treeMethod: PhylogeneticAnalyzer.TreeMethod.NeighborJoining);

        int rf12 = PhylogeneticAnalyzer.RobinsonFouldsDistance(tree1.Root, tree2.Root);
        int rf21 = PhylogeneticAnalyzer.RobinsonFouldsDistance(tree2.Root, tree1.Root);
        Assert.That(rf12, Is.EqualTo(rf21));
    }

    #endregion

    #region PHYLO-TREE-001: R: N leaves = N input sequences; P: tree connected; R: branch lengths ≥ 0

    /// <summary>
    /// INV-T1: Leaf count equals input sequence count for both UPGMA and NJ.
    /// Evidence: BuildTree constructs a binary tree where each input sequence becomes a leaf.
    /// Source: Felsenstein (2004) Inferring Phylogenies, ch. 11.
    /// </summary>
    [Test]
    [Category("Property")]
    public void BuildTree_LeafCount_EqualsInputCount()
    {
        var seqs = new Dictionary<string, string>
        {
            ["A"] = "ACGTACGTACGTACGTACGT",
            ["B"] = "ACGTACGTAAGTACGTACGT",
            ["C"] = "ACGTACATACGTACATACGT",
            ["D"] = "TCGTACATACGTACATACGT",
            ["E"] = "ACGAACGTACGTACGTACGA"
        };

        foreach (var method in new[] { PhylogeneticAnalyzer.TreeMethod.UPGMA, PhylogeneticAnalyzer.TreeMethod.NeighborJoining })
        {
            var tree = PhylogeneticAnalyzer.BuildTree(seqs, treeMethod: method);
            var leaves = PhylogeneticAnalyzer.GetLeaves(tree.Root).ToList();
            Assert.That(leaves.Count, Is.EqualTo(seqs.Count),
                $"Method {method}: expected {seqs.Count} leaves, got {leaves.Count}");
        }
    }

    /// <summary>
    /// INV-T2: All input taxa are reachable as leaves from the root (tree is connected).
    /// Evidence: Binary tree construction algorithms produce a single connected tree.
    /// </summary>
    [Test]
    [Category("Property")]
    public void BuildTree_AllLeaves_ReachableFromRoot()
    {
        var seqs = new Dictionary<string, string>
        {
            ["TaxA"] = "ACGTACGTACGTACGTACGT",
            ["TaxB"] = "ACGTACGTAAGTACGTACGT",
            ["TaxC"] = "ACGTACATACGTACATACGT",
            ["TaxD"] = "TCGTACATACGTACATACGT"
        };

        foreach (var method in new[] { PhylogeneticAnalyzer.TreeMethod.UPGMA, PhylogeneticAnalyzer.TreeMethod.NeighborJoining })
        {
            var tree = PhylogeneticAnalyzer.BuildTree(seqs, treeMethod: method);
            var leafNames = PhylogeneticAnalyzer.GetLeaves(tree.Root)
                .Select(l => l.Name)
                .OrderBy(n => n)
                .ToList();
            var expected = seqs.Keys.OrderBy(k => k).ToList();

            Assert.That(leafNames, Is.EqualTo(expected),
                $"Method {method}: not all input taxa reachable from root");
        }
    }

    /// <summary>
    /// INV-T3: UPGMA tree has all branch lengths ≥ 0 (ultrametric property).
    /// Evidence: UPGMA assigns branch length = newHeight - childHeight, always ≥ 0
    /// because the minimum distance pair is merged first.
    /// Source: Sokal &amp; Michener (1958); Wikipedia UPGMA.
    /// </summary>
    [Test]
    [Category("Property")]
    public void BuildTree_UPGMA_AllBranchLengths_NonNegative()
    {
        var seqs = new Dictionary<string, string>
        {
            ["A"] = "ACGTACGTACGTACGTACGT",
            ["B"] = "ACGTACGTAAGTACGTACGT",
            ["C"] = "ACGTACATACGTACATACGT",
            ["D"] = "TCGTACATACGTACATACGT",
            ["E"] = "ACGAACGTACGTACGTACGA"
        };
        var tree = PhylogeneticAnalyzer.BuildTree(seqs, treeMethod: PhylogeneticAnalyzer.TreeMethod.UPGMA);
        AssertAllBranchLengthsNonNegative(tree.Root);
    }

    private static void AssertAllBranchLengthsNonNegative(PhylogeneticAnalyzer.PhyloNode? node)
    {
        if (node == null) return;
        Assert.That(node.BranchLength, Is.GreaterThanOrEqualTo(0.0),
            $"Node '{node.Name}' has negative branch length {node.BranchLength}");
        AssertAllBranchLengthsNonNegative(node.Left);
        AssertAllBranchLengthsNonNegative(node.Right);
    }

    /// <summary>
    /// INV-T4: Tree construction is deterministic (same input → same topology).
    /// Evidence: BuildTree is a pure function without random state.
    /// </summary>
    [Test]
    [Category("Property")]
    public void BuildTree_IsDeterministic()
    {
        var seqs = new Dictionary<string, string>
        {
            ["A"] = "ACGTACGTACGTACGTACGT",
            ["B"] = "ACGTACGTAAGTACGTACGT",
            ["C"] = "ACGTACATACGTACATACGT"
        };

        var newick1 = PhylogeneticAnalyzer.ToNewick(PhylogeneticAnalyzer.BuildTree(seqs).Root);
        var newick2 = PhylogeneticAnalyzer.ToNewick(PhylogeneticAnalyzer.BuildTree(seqs).Root);
        Assert.That(newick1, Is.EqualTo(newick2));
    }

    #endregion

    #region PHYLO-NEWICK-001 (extended): RT: parse(serialize(tree))=tree; P: leaf labels preserved; D: deterministic

    /// <summary>
    /// INV-N1: Newick round-trip preserves branch lengths within serialization precision.
    /// Evidence: ToNewick writes F4 format; ParseNewick reads them back.
    /// Source: Olsen (1990) Newick specification.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Newick_RoundTrip_PreservesBranchLengths()
    {
        var seqs = new Dictionary<string, string>
        {
            ["A"] = "ACGTACGTACGTACGTACGT",
            ["B"] = "ACGTACGTAAGTACGTACGT",
            ["C"] = "ACGTACATACGTACATACGT"
        };
        var tree = PhylogeneticAnalyzer.BuildTree(seqs);
        string newick = PhylogeneticAnalyzer.ToNewick(tree.Root, includeBranchLengths: true);
        var parsed = PhylogeneticAnalyzer.ParseNewick(newick);

        var originalLeaves = PhylogeneticAnalyzer.GetLeaves(tree.Root)
            .OrderBy(l => l.Name).ToList();
        var parsedLeaves = PhylogeneticAnalyzer.GetLeaves(parsed)
            .OrderBy(l => l.Name).ToList();

        Assert.That(parsedLeaves.Count, Is.EqualTo(originalLeaves.Count));
        for (int i = 0; i < originalLeaves.Count; i++)
        {
            Assert.That(parsedLeaves[i].BranchLength,
                Is.EqualTo(originalLeaves[i].BranchLength).Within(0.001),
                $"Branch length mismatch for leaf '{originalLeaves[i].Name}'");
        }
    }

    /// <summary>
    /// INV-N2: ToNewick serialization is deterministic.
    /// Evidence: ToNewick is a pure recursive traversal.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Newick_ToNewick_IsDeterministic()
    {
        var seqs = new Dictionary<string, string>
        {
            ["X"] = "ACGTACGTACGTACGTACGT",
            ["Y"] = "ACGTACGTAAGTACGTACGT",
            ["Z"] = "ACGTACATACGTACATACGT",
            ["W"] = "TCGTACATACGTACATACGT"
        };
        var tree = PhylogeneticAnalyzer.BuildTree(seqs);
        string n1 = PhylogeneticAnalyzer.ToNewick(tree.Root);
        string n2 = PhylogeneticAnalyzer.ToNewick(tree.Root);
        Assert.That(n1, Is.EqualTo(n2));
    }

    #endregion

    #region PHYLO-COMP-001 (extended): S: RF(a,b)=RF(b,a); I: RF(t,t)=0; R: RF ≥ 0; R: RF ≤ 2(n-2)

    /// <summary>
    /// INV-C1: Robinson-Foulds distance is non-negative (count of set-symmetric difference).
    /// Evidence: RF = |Clades(T1) △ Clades(T2)| ≥ 0 by definition.
    /// Source: Robinson &amp; Foulds (1981) "Comparison of phylogenetic trees".
    /// </summary>
    [Test]
    [Category("Property")]
    public void RobinsonFoulds_IsNonNegative_Explicit()
    {
        var seqs = new Dictionary<string, string>
        {
            ["A"] = "ACGTACGTACGTACGTACGT",
            ["B"] = "ACGTACGTAAGTACGTACGT",
            ["C"] = "ACGTACATACGTACATACGT",
            ["D"] = "TCGTACATACGTACATACGT"
        };
        var tree1 = PhylogeneticAnalyzer.BuildTree(seqs, treeMethod: PhylogeneticAnalyzer.TreeMethod.UPGMA);
        var tree2 = PhylogeneticAnalyzer.BuildTree(seqs, treeMethod: PhylogeneticAnalyzer.TreeMethod.NeighborJoining);

        int rf = PhylogeneticAnalyzer.RobinsonFouldsDistance(tree1.Root, tree2.Root);
        Assert.That(rf, Is.GreaterThanOrEqualTo(0),
            $"RF distance must be ≥ 0, got {rf}");
    }

    /// <summary>
    /// INV-C2: RF distance upper bound: RF ≤ 2(n-2) for rooted binary trees with n leaves.
    /// Evidence: A rooted binary tree with n leaves has at most n-2 non-trivial clades
    /// (excluding single-leaf and root clades). The symmetric difference ≤ 2(n-2).
    /// Source: Robinson &amp; Foulds (1981); Felsenstein (2004) ch. 30.
    /// </summary>
    [Test]
    [Category("Property")]
    public void RobinsonFoulds_UpperBound()
    {
        var seqs = new Dictionary<string, string>
        {
            ["A"] = "ACGTACGTACGTACGTACGT",
            ["B"] = "ACGTACGTAAGTACGTACGT",
            ["C"] = "ACGTACATACGTACATACGT",
            ["D"] = "TCGTACATACGTACATACGT",
            ["E"] = "ACGAACGTACGTACGTACGA"
        };
        int n = seqs.Count;
        var tree1 = PhylogeneticAnalyzer.BuildTree(seqs, treeMethod: PhylogeneticAnalyzer.TreeMethod.UPGMA);
        var tree2 = PhylogeneticAnalyzer.BuildTree(seqs, treeMethod: PhylogeneticAnalyzer.TreeMethod.NeighborJoining);

        int rf = PhylogeneticAnalyzer.RobinsonFouldsDistance(tree1.Root, tree2.Root);
        int upperBound = 2 * (n - 2);
        Assert.That(rf, Is.LessThanOrEqualTo(upperBound),
            $"RF={rf} exceeds upper bound 2(n-2)={upperBound} for n={n} leaves");
    }

    #endregion
}
