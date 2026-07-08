using Seqeron.Mcp.Alignment.Tools;
using Seqeron.Mcp.Tests.Contracts;
using NUnit.Framework;

namespace Seqeron.Mcp.Alignment.Tests;

/// <summary>
/// Pins the Seqeron.Mcp.Alignment server's public MCP tool surface. Assertions live in
/// <see cref="ToolRegistrationContractTestsBase"/>; this fixture supplies the tool class and the
/// documented tool-name set.
/// </summary>
[TestFixture]
public sealed class ToolRegistrationContractTests : ToolRegistrationContractTestsBase
{
    protected override Type ToolsType => typeof(AlignmentTools);

    protected override string[] ExpectedToolNames { get; } =
    {
        "alignment_statistics",
        "assemble_de_bruijn",
        "assemble_olc",
        "assembly_stats",
        "calculate_coverage",
        "compute_consensus",
        "error_correct_reads",
        "find_all_overlaps",
        "find_best_match",
        "find_overlap",
        "find_with_edits",
        "find_with_mismatches",
        "format_alignment",
        "frequent_kmers_with_mismatches",
        "global_align",
        "local_align",
        "merge_contigs",
        "multiple_align",
        "quality_trim_reads",
        "scaffold_contigs",
        "semi_global_align",
        "sequence_identity",
    };
}
