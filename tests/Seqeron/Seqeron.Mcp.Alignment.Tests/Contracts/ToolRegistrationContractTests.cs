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
        "clustalw_command",
        "clustalw_pairwise_align",
        "clustalw_profile_align",
        "clustalw_quick_pairwise_align",
        "clustalw_sequences_to_profile",
        "compute_consensus",
        "damerau_alignment",
        "damerau_levenshtein_distance",
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
        "kalign_command",
        "local_align",
        "local_align_affine",
        "local_align_co_optimal",
        "local_align_top_n",
        "merge_contigs",
        "multiple_align",
        "multiple_align_anchored",
        "multiple_align_clustalw",
        "multiple_align_consistency",
        "multiple_align_exact",
        "multiple_align_iterative",
        "multiple_align_kalign",
        "multiple_align_mafft",
        "multiple_align_muscle",
        "multiple_align_progressive",
        "multiple_align_tcoffee",
        "muscle_msa_command",
        "quality_trim_reads",
        "scaffold_contigs",
        "semi_global_align",
        "semi_global_align_affine",
        "sequence_identity",
        "stretcher_align",
        "sum_of_pairs_score",
        "tcoffee_guide_tree",
        "tcoffee_ktup_similarity",
        "tcoffee_lalign_top_alignments",
        "tcoffee_primary_library",
        "tcoffee_slow_pair_align",
        "write_alignment",
    };
}
