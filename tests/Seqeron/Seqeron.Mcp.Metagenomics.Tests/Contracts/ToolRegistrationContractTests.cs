using Seqeron.Mcp.Metagenomics.Tools;
using Seqeron.Mcp.Tests.Contracts;
using NUnit.Framework;

namespace Seqeron.Mcp.Metagenomics.Tests;

/// <summary>
/// Pins the Seqeron.Mcp.Metagenomics server's public MCP tool surface. Assertions live in
/// <see cref="ToolRegistrationContractTestsBase"/>; this fixture supplies the tool class and the
/// documented tool-name set.
/// </summary>
[TestFixture]
public sealed class ToolRegistrationContractTests : ToolRegistrationContractTestsBase
{
    protected override Type ToolsType => typeof(MetagenomicsTools);

    protected override string[] ExpectedToolNames { get; } =
    {
        "accessory_genes",
        "alpha_diversity",
        "beta_diversity",
        "bin_contigs",
        "build_kmer_database",
        "classify_reads",
        "cluster_genes",
        "construct_pangenome",
        "core_gene_clusters",
        "core_genome_alignment",
        "differential_abundance",
        "find_genome_specific_genes",
        "find_resistance_genes",
        "fit_heaps_law",
        "functional_diversity",
        "gene_presence_absence_matrix",
        "predict_functions",
        "select_phylogenetic_markers",
        "taxonomic_profile",
    };
}
