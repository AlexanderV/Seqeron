using System.ComponentModel;
using ModelContextProtocol.Server;
using Seqeron.Genomics.Core;
using Seqeron.Genomics.Phylogenetics;
using Seqeron.Mcp.Phylogenetics.Models;

namespace Seqeron.Mcp.Phylogenetics.Tools;

[McpServerToolType]
public class PhylogeneticsTools
{
    #region PhylogeneticAnalyzer

    /// <summary>Build a phylogenetic tree from aligned sequences.</summary>
    [McpServerTool(Name = "build_phylogenetic_tree", Title = "Phylogenetics — Build Tree", ReadOnly = true)]
    [Description("Build a phylogenetic tree from a set of named, pre-aligned sequences. Computes a pairwise distance matrix using the chosen substitution model (PDistance | JukesCantor | Kimura2Parameter | Hamming) and constructs the tree using UPGMA or NeighborJoining. Returns the Newick string, taxa, distance matrix, and method name.")]
    public static PhyloTreeBuildResult BuildPhylogeneticTree(
        [Description("Aligned sequences keyed by taxon name. All values must be the same length.")] Dictionary<string, string> sequences,
        [Description("Distance method: PDistance | JukesCantor | Kimura2Parameter | Hamming. Default: JukesCantor.")] string? distanceMethod = null,
        [Description("Tree construction method: UPGMA | NeighborJoining. Default: UPGMA.")] string? treeMethod = null)
    {
        var dm = ParseDistanceMethod(distanceMethod);
        var tm = ParseTreeMethod(treeMethod);
        var tree = global::Seqeron.Genomics.Phylogenetics.PhylogeneticAnalyzer.BuildTree(sequences, dm, tm);
        var newick = global::Seqeron.Genomics.Phylogenetics.PhylogeneticAnalyzer.ToNewick(tree.Root);
        return new PhyloTreeBuildResult(
            newick,
            tree.Taxa.ToList(),
            ToJagged(tree.DistanceMatrix),
            tree.Method);
    }

    /// <summary>Build a phylogenetic tree from a precomputed distance matrix.</summary>
    [McpServerTool(Name = "build_tree_from_matrix", Title = "Phylogenetics — Build Tree From Matrix", ReadOnly = true)]
    [Description("Build a phylogenetic tree directly from a precomputed symmetric distance matrix (UPGMA or NeighborJoining). Useful for verifying tree-construction algorithms against reference matrices.")]
    public static PhyloTreeBuildResult BuildTreeFromMatrix(
        [Description("Taxon names in the same order as the matrix rows/columns.")] string[] taxa,
        [Description("Symmetric square distance matrix; size must equal the length of taxa.")] double[][] distanceMatrix,
        [Description("Tree construction method: UPGMA | NeighborJoining. Default: UPGMA.")] string? treeMethod = null)
    {
        var tm = ParseTreeMethod(treeMethod);
        var matrix = ToMultiDim(distanceMatrix, taxa?.Length ?? 0);
        var tree = global::Seqeron.Genomics.Phylogenetics.PhylogeneticAnalyzer.BuildTreeFromMatrix(taxa!, matrix, tm);
        var newick = global::Seqeron.Genomics.Phylogenetics.PhylogeneticAnalyzer.ToNewick(tree.Root);
        return new PhyloTreeBuildResult(
            newick,
            tree.Taxa.ToList(),
            ToJagged(tree.DistanceMatrix),
            tree.Method);
    }

    /// <summary>Compute the pairwise distance matrix for aligned sequences.</summary>
    [McpServerTool(Name = "distance_matrix", Title = "Phylogenetics — Distance Matrix", ReadOnly = true)]
    [Description("Compute the symmetric pairwise distance matrix for a list of aligned sequences using the chosen substitution model (PDistance | JukesCantor | Kimura2Parameter | Hamming). Diagonal is zero.")]
    public static DistanceMatrixResult DistanceMatrix(
        [Description("Aligned sequences (equal length).")] string[] alignedSequences,
        [Description("Distance method: PDistance | JukesCantor | Kimura2Parameter | Hamming. Default: JukesCantor.")] string? method = null)
    {
        var dm = ParseDistanceMethod(method);
        var matrix = global::Seqeron.Genomics.Phylogenetics.PhylogeneticAnalyzer.CalculateDistanceMatrix(alignedSequences, dm);
        return new DistanceMatrixResult(ToJagged(matrix));
    }

    /// <summary>Calculate evolutionary distance between two aligned sequences.</summary>
    [McpServerTool(Name = "pairwise_distance", Title = "Phylogenetics — Pairwise Distance", ReadOnly = true)]
    [Description("Calculate the evolutionary distance between two aligned sequences under a chosen substitution model: p-distance, Hamming, Jukes-Cantor (JC69), or Kimura 2-parameter (K80). JC69/K2P may return +Infinity at saturation.")]
    public static PairwiseDistanceResult PairwiseDistance(
        [Description("First aligned sequence.")] string seq1,
        [Description("Second aligned sequence (same length as seq1).")] string seq2,
        [Description("Distance method: PDistance | JukesCantor | Kimura2Parameter | Hamming. Default: JukesCantor.")] string? method = null)
    {
        var dm = ParseDistanceMethod(method);
        var d = global::Seqeron.Genomics.Phylogenetics.PhylogeneticAnalyzer.CalculatePairwiseDistance(seq1, seq2, dm);
        return new PairwiseDistanceResult(d);
    }

    /// <summary>Serialize a tree (provided as Newick) back to Newick format.</summary>
    [McpServerTool(Name = "to_newick", Title = "Phylogenetics — To Newick", ReadOnly = true)]
    [Description("Serialize a phylogenetic tree to Newick format. The tool accepts a Newick string, parses it internally, and re-emits canonical Newick. Round-trips with parse_newick for trees whose labels are valid unquoted Newick labels (Olsen 1990). Branch lengths are formatted with F4 invariant culture.")]
    public static ToNewickResult ToNewick(
        [Description("Input tree in Newick format.")] string newick,
        [Description("Whether to include branch lengths. Default: true.")] bool includeBranchLengths = true)
    {
        var node = global::Seqeron.Genomics.Phylogenetics.PhylogeneticAnalyzer.ParseNewick(newick);
        var serialized = global::Seqeron.Genomics.Phylogenetics.PhylogeneticAnalyzer.ToNewick(node, includeBranchLengths);
        return new ToNewickResult(serialized);
    }

    /// <summary>Parse a Newick tree and report a structural summary.</summary>
    [McpServerTool(Name = "parse_newick", Title = "Phylogenetics — Parse Newick", ReadOnly = true)]
    [Description("Parse a Newick-format tree string and report a summary: canonical re-serialization (round-trip), taxa list, leaf count, depth, and total branch length. Trailing ';' is stripped; an optional root branch length is supported. Parser is lenient with malformed numbers.")]
    public static ParseNewickResult ParseNewick(
        [Description("Tree in Newick format.")] string newick)
    {
        var node = global::Seqeron.Genomics.Phylogenetics.PhylogeneticAnalyzer.ParseNewick(newick);
        var canonical = global::Seqeron.Genomics.Phylogenetics.PhylogeneticAnalyzer.ToNewick(node);
        var leaves = global::Seqeron.Genomics.Phylogenetics.PhylogeneticAnalyzer.GetLeaves(node).ToList();
        var depth = global::Seqeron.Genomics.Phylogenetics.PhylogeneticAnalyzer.GetTreeDepth(node);
        var totalLen = global::Seqeron.Genomics.Phylogenetics.PhylogeneticAnalyzer.CalculateTreeLength(node);
        var taxa = leaves.Select(l => l.Name).ToList();
        return new ParseNewickResult(canonical, taxa, leaves.Count, depth, totalLen);
    }

    /// <summary>Enumerate the leaf (taxon) nodes of a tree.</summary>
    [McpServerTool(Name = "tree_leaves", Title = "Phylogenetics — Tree Leaves", ReadOnly = true)]
    [Description("Enumerate the leaf (taxon) nodes of a tree, returning each leaf's name and branch length.")]
    public static TreeLeavesResult TreeLeaves(
        [Description("Tree in Newick format.")] string newick)
    {
        var node = global::Seqeron.Genomics.Phylogenetics.PhylogeneticAnalyzer.ParseNewick(newick);
        var items = global::Seqeron.Genomics.Phylogenetics.PhylogeneticAnalyzer
            .GetLeaves(node)
            .Select(l => new TreeLeafItem(l.Name, l.BranchLength))
            .ToList();
        return new TreeLeavesResult(items);
    }

    /// <summary>Sum of all branch lengths in a tree.</summary>
    [McpServerTool(Name = "tree_length", Title = "Phylogenetics — Tree Length", ReadOnly = true)]
    [Description("Sum of all branch lengths in a tree. Negative branch lengths (possible from Neighbor-Joining output) are summed as-is.")]
    public static TreeLengthResult TreeLength(
        [Description("Tree in Newick format.")] string newick)
    {
        var node = global::Seqeron.Genomics.Phylogenetics.PhylogeneticAnalyzer.ParseNewick(newick);
        var len = global::Seqeron.Genomics.Phylogenetics.PhylogeneticAnalyzer.CalculateTreeLength(node);
        return new TreeLengthResult(len);
    }

    /// <summary>Maximum number of internal-node edges from root to any leaf.</summary>
    [McpServerTool(Name = "tree_depth", Title = "Phylogenetics — Tree Depth", ReadOnly = true)]
    [Description("Maximum number of internal-node edges from root to any leaf. A single-leaf tree has depth 0.")]
    public static TreeDepthResult TreeDepth(
        [Description("Tree in Newick format.")] string newick)
    {
        var node = global::Seqeron.Genomics.Phylogenetics.PhylogeneticAnalyzer.ParseNewick(newick);
        var d = global::Seqeron.Genomics.Phylogenetics.PhylogeneticAnalyzer.GetTreeDepth(node);
        return new TreeDepthResult(d);
    }

    /// <summary>Robinson-Foulds distance between two trees.</summary>
    [McpServerTool(Name = "robinson_foulds_distance", Title = "Phylogenetics — Robinson-Foulds Distance", ReadOnly = true)]
    [Description("Robinson-Foulds (symmetric clade-difference) distance between two rooted trees over the same taxon set. Compares non-trivial clades only; trivial clades (single leaf, full taxon set) are excluded.")]
    public static RobinsonFouldsDistanceResult RobinsonFouldsDistance(
        [Description("First tree in Newick format.")] string tree1,
        [Description("Second tree in Newick format.")] string tree2)
    {
        var n1 = global::Seqeron.Genomics.Phylogenetics.PhylogeneticAnalyzer.ParseNewick(tree1);
        var n2 = global::Seqeron.Genomics.Phylogenetics.PhylogeneticAnalyzer.ParseNewick(tree2);
        var d = global::Seqeron.Genomics.Phylogenetics.PhylogeneticAnalyzer.RobinsonFouldsDistance(n1, n2);
        return new RobinsonFouldsDistanceResult(d);
    }

    /// <summary>Most Recent Common Ancestor of two taxa in a rooted tree.</summary>
    [McpServerTool(Name = "mrca", Title = "Phylogenetics — Most Recent Common Ancestor", ReadOnly = true)]
    [Description("Most Recent Common Ancestor (MRCA) of two taxa in a rooted tree. Returns found=false (and empty fields) when the root is null or either taxon is missing. Self-MRCA is the leaf when taxon1 == taxon2 and that taxon exists.")]
    public static MrcaResult Mrca(
        [Description("Tree in Newick format.")] string newick,
        [Description("First taxon name.")] string taxon1,
        [Description("Second taxon name.")] string taxon2)
    {
        var root = global::Seqeron.Genomics.Phylogenetics.PhylogeneticAnalyzer.ParseNewick(newick);
        var mrca = global::Seqeron.Genomics.Phylogenetics.PhylogeneticAnalyzer.FindMRCA(root, taxon1, taxon2);
        if (mrca is null)
            return new MrcaResult(false, string.Empty, Array.Empty<string>());
        return new MrcaResult(true, mrca.Name, mrca.Taxa.ToList());
    }

    /// <summary>Patristic distance between two taxa in a rooted tree.</summary>
    [McpServerTool(Name = "patristic_distance", Title = "Phylogenetics — Patristic Distance", ReadOnly = true)]
    [Description("Sum of branch lengths along the unique path between two taxa in a rooted tree (via MRCA). Returns NaN when MRCA cannot be found (i.e., one or both taxa missing). Identical taxa returns 0.")]
    public static PatristicDistanceResult PatristicDistance(
        [Description("Tree in Newick format.")] string newick,
        [Description("First taxon name.")] string taxon1,
        [Description("Second taxon name.")] string taxon2)
    {
        var root = global::Seqeron.Genomics.Phylogenetics.PhylogeneticAnalyzer.ParseNewick(newick);
        var d = global::Seqeron.Genomics.Phylogenetics.PhylogeneticAnalyzer.PatristicDistance(root, taxon1, taxon2);
        return new PatristicDistanceResult(d);
    }

    /// <summary>Non-parametric bootstrap support for clades of the reference tree.</summary>
    [McpServerTool(Name = "bootstrap_support", Title = "Phylogenetics — Bootstrap Support", ReadOnly = true)]
    [Description("Non-parametric bootstrap support for clades of a reference tree. Resamples alignment columns with replacement and reports the proportion of replicates in which each non-trivial reference clade reappears. The internal RNG is seeded with constant 42, so results are deterministic for identical inputs. Cost is O(replicates * n^2 * L).")]
    public static BootstrapSupportResult BootstrapSupport(
        [Description("Aligned sequences keyed by taxon name (equal length).")] Dictionary<string, string> sequences,
        [Description("Number of bootstrap replicates. Default: 100.")] int replicates = 100,
        [Description("Distance method: PDistance | JukesCantor | Kimura2Parameter | Hamming. Default: JukesCantor.")] string? distanceMethod = null,
        [Description("Tree construction method: UPGMA | NeighborJoining. Default: UPGMA.")] string? treeMethod = null)
    {
        if (sequences is null)
            throw new ArgumentException("Sequences cannot be null.", nameof(sequences));
        if (sequences.Count < 2)
            throw new ArgumentException("At least 2 sequences required.", nameof(sequences));
        if (replicates < 1)
            throw new ArgumentException("At least 1 replicate required.", nameof(replicates));

        var dm = ParseDistanceMethod(distanceMethod);
        var tm = ParseTreeMethod(treeMethod);
        var support = global::Seqeron.Genomics.Phylogenetics.PhylogeneticAnalyzer.Bootstrap(sequences, replicates, dm, tm);
        var items = support.Select(kvp => new BootstrapSupportItem(kvp.Key, kvp.Value)).ToList();
        return new BootstrapSupportResult(items);
    }

    // ----------------------------------------------------------------------
    // Internal helpers (not MCP tools).
    // ----------------------------------------------------------------------

    private static global::Seqeron.Genomics.Phylogenetics.PhylogeneticAnalyzer.DistanceMethod ParseDistanceMethod(string? s)
    {
        var k = (s ?? string.Empty).Trim().ToLowerInvariant().Replace("-", "").Replace("_", "");
        return k switch
        {
            "" or "jukescantor" or "jc" or "jc69" => global::Seqeron.Genomics.Phylogenetics.PhylogeneticAnalyzer.DistanceMethod.JukesCantor,
            "pdistance" or "p" => global::Seqeron.Genomics.Phylogenetics.PhylogeneticAnalyzer.DistanceMethod.PDistance,
            "kimura2parameter" or "kimura" or "k2p" or "k80" => global::Seqeron.Genomics.Phylogenetics.PhylogeneticAnalyzer.DistanceMethod.Kimura2Parameter,
            "hamming" => global::Seqeron.Genomics.Phylogenetics.PhylogeneticAnalyzer.DistanceMethod.Hamming,
            _ => throw new ArgumentException($"Unknown distance method: '{s}'. Expected: PDistance | JukesCantor | Kimura2Parameter | Hamming.")
        };
    }

    private static global::Seqeron.Genomics.Phylogenetics.PhylogeneticAnalyzer.TreeMethod ParseTreeMethod(string? s)
    {
        var k = (s ?? string.Empty).Trim().ToLowerInvariant().Replace("-", "").Replace("_", "");
        return k switch
        {
            "" or "upgma" => global::Seqeron.Genomics.Phylogenetics.PhylogeneticAnalyzer.TreeMethod.UPGMA,
            "neighborjoining" or "nj" or "neighbourjoining" => global::Seqeron.Genomics.Phylogenetics.PhylogeneticAnalyzer.TreeMethod.NeighborJoining,
            _ => throw new ArgumentException($"Unknown tree method: '{s}'. Expected: UPGMA | NeighborJoining.")
        };
    }

    private static double[][] ToJagged(double[,] m)
    {
        int rows = m.GetLength(0);
        int cols = m.GetLength(1);
        var result = new double[rows][];
        for (int i = 0; i < rows; i++)
        {
            var row = new double[cols];
            for (int j = 0; j < cols; j++) row[j] = m[i, j];
            result[i] = row;
        }
        return result;
    }

    private static double[,] ToMultiDim(double[][] m, int expectedSize)
    {
        if (m is null)
            throw new ArgumentException("Distance matrix is required.", nameof(m));
        int rows = m.Length;
        if (rows != expectedSize)
            throw new ArgumentException($"Distance matrix has {rows} rows but expected {expectedSize} (taxa.Length).");
        var result = new double[rows, rows];
        for (int i = 0; i < rows; i++)
        {
            var row = m[i];
            if (row is null || row.Length != rows)
                throw new ArgumentException($"Distance matrix row {i} must have length {rows}.");
            for (int j = 0; j < rows; j++) result[i, j] = row[j];
        }
        return result;
    }

    #endregion
}
