using Seqeron.Mcp.Population.Tools;
using Seqeron.Mcp.Tests.Contracts;
using NUnit.Framework;

namespace Seqeron.Mcp.Population.Tests;

/// <summary>
/// Pins the Seqeron.Mcp.Population server's public MCP tool surface. Assertions live in
/// <see cref="ToolRegistrationContractTestsBase"/>; this fixture supplies the tool class and the
/// documented tool-name set.
/// </summary>
[TestFixture]
public sealed class ToolRegistrationContractTests : ToolRegistrationContractTestsBase
{
    protected override Type ToolsType => typeof(PopulationTools);

    protected override string[] ExpectedToolNames { get; } =
    {
        "allele_frequencies",
        "diversity_statistics",
        "estimate_ancestry",
        "f_statistics",
        "filter_variants_by_maf",
        "fst",
        "haplotype_blocks",
        "hardy_weinberg_test",
        "inbreeding_from_roh",
        "integrated_haplotype_score",
        "linkage_disequilibrium",
        "minor_allele_frequency",
        "nucleotide_diversity",
        "pairwise_fst",
        "runs_of_homozygosity",
        "scan_selection_signals",
        "tajimas_d",
        "wattersons_theta",
    };
}
