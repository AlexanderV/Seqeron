using Seqeron.Mcp.Chromosome.Tools;
using Seqeron.Mcp.Tests.Contracts;
using NUnit.Framework;

namespace Seqeron.Mcp.Chromosome.Tests;

/// <summary>
/// Pins the Seqeron.Mcp.Chromosome server's public MCP tool surface. Assertions live in
/// <see cref="ToolRegistrationContractTestsBase"/>; this fixture supplies the tool class and the
/// documented tool-name set.
/// </summary>
[TestFixture]
public sealed class ToolRegistrationContractTests : ToolRegistrationContractTestsBase
{
    protected override Type ToolsType => typeof(ChromosomeTools);

    protected override string[] ExpectedToolNames { get; } =
    {
        "analyze_centromere",
        "analyze_karyotype",
        "analyze_scaffolds",
        "analyze_telomeres",
        "arm_ratio",
        "assembly_statistics",
        "assess_completeness",
        "au_n",
        "classify_chromosome_by_arm_ratio",
        "compare_assemblies",
        "detect_aneuploidy",
        "detect_ploidy",
        "detect_rearrangements",
        "estimate_cell_divisions_from_telomere_length",
        "estimate_completeness_from_kmers",
        "estimate_telomere_length_from_ts_ratio",
        "extract_contigs",
        "find_gaps",
        "find_heterochromatin_regions",
        "find_repetitive_regions",
        "find_suspicious_regions",
        "find_syntenic_blocks_assemblies",
        "find_synteny_blocks",
        "find_tandem_repeats",
        "gap_distribution",
        "identify_whole_chromosome_aneuploidy",
        "length_distribution",
        "local_quality",
        "nx_curve",
        "nx_statistics",
        "predict_g_bands",
        "repeat_content",
    };
}
