namespace Seqeron.Mcp.Phylogenetics.Models;

// ================================
// Phylogenetics Results
// ================================

/// <summary>Result of build_phylogenetic_tree and build_tree_from_matrix.</summary>
public record PhyloTreeBuildResult(
    string Newick,
    IReadOnlyList<string> Taxa,
    double[][] DistanceMatrix,
    string Method);

/// <summary>Result of distance_matrix.</summary>
public record DistanceMatrixResult(double[][] Matrix);

/// <summary>Result of pairwise_distance.</summary>
public record PairwiseDistanceResult(double Distance);

/// <summary>Result of to_newick.</summary>
public record ToNewickResult(string Newick);

/// <summary>Result of parse_newick.</summary>
public record ParseNewickResult(
    string Newick,
    IReadOnlyList<string> Taxa,
    int LeafCount,
    int Depth,
    double TotalLength);

/// <summary>Single leaf entry returned by tree_leaves.</summary>
public record TreeLeafItem(string Name, double BranchLength);

/// <summary>Result of tree_leaves.</summary>
public record TreeLeavesResult(IReadOnlyList<TreeLeafItem> Leaves);

/// <summary>Result of tree_length.</summary>
public record TreeLengthResult(double Length);

/// <summary>Result of tree_depth.</summary>
public record TreeDepthResult(int Depth);

/// <summary>Result of robinson_foulds_distance.</summary>
public record RobinsonFouldsDistanceResult(int Distance);

/// <summary>Result of mrca.</summary>
public record MrcaResult(bool Found, string Name, IReadOnlyList<string> Taxa);

/// <summary>Result of patristic_distance.</summary>
public record PatristicDistanceResult(double Distance);

/// <summary>Single bootstrap support entry returned by bootstrap_support.</summary>
public record BootstrapSupportItem(string Clade, double Support);

/// <summary>Result of bootstrap_support.</summary>
public record BootstrapSupportResult(IReadOnlyList<BootstrapSupportItem> Support);
