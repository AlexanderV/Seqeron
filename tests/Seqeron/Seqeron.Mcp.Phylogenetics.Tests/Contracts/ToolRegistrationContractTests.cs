using Seqeron.Mcp.Phylogenetics.Tools;
using Seqeron.Mcp.Tests.Contracts;
using NUnit.Framework;

namespace Seqeron.Mcp.Phylogenetics.Tests;

/// <summary>
/// Pins the Seqeron.Mcp.Phylogenetics server's public MCP tool surface. Assertions live in
/// <see cref="ToolRegistrationContractTestsBase"/>; this fixture supplies the tool class and the
/// documented tool-name set.
/// </summary>
[TestFixture]
public sealed class ToolRegistrationContractTests : ToolRegistrationContractTestsBase
{
    protected override Type ToolsType => typeof(PhylogeneticsTools);

    protected override string[] ExpectedToolNames { get; } =
    {
        "bootstrap_support",
        "build_phylogenetic_tree",
        "build_tree_from_matrix",
        "distance_matrix",
        "mrca",
        "pairwise_distance",
        "parse_newick",
        "patristic_distance",
        "robinson_foulds_distance",
        "to_newick",
        "tree_depth",
        "tree_leaves",
        "tree_length",
    };
}
