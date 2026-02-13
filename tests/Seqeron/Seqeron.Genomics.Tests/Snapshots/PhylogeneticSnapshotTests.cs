namespace Seqeron.Genomics.Tests.Snapshots;

/// <summary>
/// Snapshot tests for phylogenetic analysis.
/// Verifies tree construction, distance matrix, and Newick output stability.
///
/// Test Units: PHYLO-TREE-001, PHYLO-NEWICK-001, PHYLO-DIST-001
/// </summary>
[TestFixture]
[Category("Snapshot")]
[Category("Phylogenetics")]
public class PhylogeneticSnapshotTests
{
    private readonly Dictionary<string, string> _testSequences = new()
    {
        ["Human"] = "ACGTACGTACGTACGTACGT",
        ["Chimp"] = "ACGTACGTAAGTACGTACGT",
        ["Gorilla"] = "ACGTACATACGTACATACGT",
        ["Mouse"] = "TCATACATACGTACACACCG"
    };

    [Test]
    public Task BuildTree_UPGMA_MatchesSnapshot()
    {
        var tree = PhylogeneticAnalyzer.BuildTree(_testSequences, treeMethod: PhylogeneticAnalyzer.TreeMethod.UPGMA);
        string newick = PhylogeneticAnalyzer.ToNewick(tree.Root);
        var leaves = PhylogeneticAnalyzer.GetLeaves(tree.Root)
            .Select(l => new { l.Name, BranchLength = Math.Round(l.BranchLength, 4) })
            .OrderBy(l => l.Name)
            .ToList();
        double length = Math.Round(PhylogeneticAnalyzer.CalculateTreeLength(tree.Root), 4);

        return Verify(new
        {
            Method = tree.Method,
            Newick = newick,
            Leaves = leaves,
            TotalTreeLength = length,
            TaxaCount = tree.Taxa.Count
        });
    }

    [Test]
    public Task BuildTree_NeighborJoining_MatchesSnapshot()
    {
        var tree = PhylogeneticAnalyzer.BuildTree(_testSequences, treeMethod: PhylogeneticAnalyzer.TreeMethod.NeighborJoining);
        string newick = PhylogeneticAnalyzer.ToNewick(tree.Root);
        var leaves = PhylogeneticAnalyzer.GetLeaves(tree.Root)
            .Select(l => new { l.Name, BranchLength = Math.Round(l.BranchLength, 4) })
            .OrderBy(l => l.Name)
            .ToList();

        return Verify(new
        {
            Method = tree.Method,
            Newick = newick,
            Leaves = leaves,
            TaxaCount = tree.Taxa.Count
        });
    }

    [Test]
    public Task DistanceMatrix_JukesCantor_MatchesSnapshot()
    {
        var seqList = _testSequences.Values.ToList();
        var matrix = PhylogeneticAnalyzer.CalculateDistanceMatrix(seqList);
        var names = _testSequences.Keys.ToList();

        var distances = new Dictionary<string, double>();
        for (int i = 0; i < names.Count; i++)
            for (int j = i + 1; j < names.Count; j++)
                distances[$"{names[i]}-{names[j]}"] = Math.Round(matrix[i, j], 4);

        return Verify(new { Distances = distances });
    }
}
