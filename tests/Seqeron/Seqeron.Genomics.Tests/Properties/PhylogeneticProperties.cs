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
}
