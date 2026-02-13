using FsCheck;
using FsCheck.Fluent;
using FsCheck.NUnit;

namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Property-based tests for phylogenetic analysis.
/// Verifies distance matrix and Newick I/O invariants.
///
/// Test Units: PHYLO-DIST-001, PHYLO-NEWICK-001, PHYLO-COMP-001 (Property Extensions)
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

    /// <summary>
    /// Distance matrix is symmetric: d[i,j] == d[j,i].
    /// </summary>
    [Test]
    [Category("Property")]
    public void DistanceMatrix_IsSymmetric()
    {
        var seqs = new[] { "ACGTACGTACGTACGTACGT", "ACGTACGTAAGTACGTACGT", "ACGTACGTACGTACATACGT" };
        var matrix = PhylogeneticAnalyzer.CalculateDistanceMatrix(seqs);

        for (int i = 0; i < seqs.Length; i++)
            for (int j = 0; j < seqs.Length; j++)
                Assert.That(matrix[i, j], Is.EqualTo(matrix[j, i]).Within(0.0001),
                    $"d[{i},{j}] != d[{j},{i}]");
    }

    /// <summary>
    /// Diagonal of distance matrix is zero: d[i,i] == 0.
    /// </summary>
    [Test]
    [Category("Property")]
    public void DistanceMatrix_DiagonalIsZero()
    {
        var seqs = new[] { "ACGTACGTACGTACGTACGT", "ACGTACGTAAGTACGTACGT", "ACGTACGTACGTACATACGT" };
        var matrix = PhylogeneticAnalyzer.CalculateDistanceMatrix(seqs);

        for (int i = 0; i < seqs.Length; i++)
            Assert.That(matrix[i, i], Is.EqualTo(0.0).Within(0.0001), $"d[{i},{i}] must be 0");
    }

    /// <summary>
    /// All off-diagonal distances are non-negative.
    /// </summary>
    [Test]
    [Category("Property")]
    public void DistanceMatrix_AllNonNegative()
    {
        var seqs = new[] { "ACGTACGTACGTACGTACGT", "TCGTACGTAAGTACGTACGT", "ACGTACGTACGTACATACGT" };
        var matrix = PhylogeneticAnalyzer.CalculateDistanceMatrix(seqs);

        for (int i = 0; i < seqs.Length; i++)
            for (int j = 0; j < seqs.Length; j++)
                Assert.That(matrix[i, j], Is.GreaterThanOrEqualTo(0.0),
                    $"d[{i},{j}] must be ≥ 0");
    }

    /// <summary>
    /// Identical sequences have distance 0.
    /// </summary>
    [Test]
    [Category("Property")]
    public void IdenticalSequences_HaveZeroDistance()
    {
        string seq = "ACGTACGTACGTACGTACGT";
        double d = PhylogeneticAnalyzer.CalculatePairwiseDistance(seq, seq);
        Assert.That(d, Is.EqualTo(0.0).Within(0.0001));
    }

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
}
