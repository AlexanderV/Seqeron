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
        "align_gap_profile",
        "alignment_counts",
        "alignment_operation",
        "alignment_statistics",
        "assemble_de_bruijn",
        "assemble_olc",
        "assembly_stats",
        "banded_align",
        "calculate_coverage",
        "compute_consensus",
        "damerau_alignment",
        "damerau_levenshtein_distance",
        "dotmatcher",
        "dual_affine_align",
        "edit_alignment",
        "error_correct_reads",
        "find_all_overlaps",
        "find_best_match",
        "find_edit_end_positions",
        "find_overlap",
        "find_with_edits",
        "find_with_mismatches",
        "format_alignment",
        "frequent_kmers_with_mismatches",
        "frequent_kmers_with_mismatches_and_revcomp",
        "general_gap_align",
        "global_align",
        "global_align_affine",
        "global_align_co_optimal",
        "ksw2_extend",
        "local_align",
        "local_align_affine",
        "local_align_co_optimal",
        "local_align_top_n",
        "local_align_zscore",
        "mafft_command",
        "merge_contigs",
        "merger",
        "minimap2_chains",
        "minimap2_map",
        "multiple_align",
        "multiple_align_anchored",
        "multiple_align_exact",
        "multiple_align_mafft",
        "multiple_align_muscle",
        "muscle_msa_command",
        "needle_align",
        "needleall_align",
        "pairwise_aligner",
        "codon_aligner",
        "multiple_seq_alignment",
        "parasail_score",
        "quality_trim_reads",
        "scaffold_contigs",
        "semi_global_align",
        "semi_global_align_affine",
        "sequence_identity",
        "spoa_command",
        "stretcher_align",
        "sum_of_pairs_score",
        "supermatcher_align",
        "water_align",
        "write_alignment",
        "write_coordinate_alignment",
        "parse_coordinate_alignment",
        "write_big_alignment",
        "parse_big_alignment",
    };
}
