// SPDX-License-Identifier: GPL-3.0-or-later
// Wrappers of the copyleft multiple-alignment assemblies (Kalign 3.3.5, GPL-3.0-or-later and others).
// Part of Seqeron.Mcp.MultipleAlignment, licensed under the GNU General Public License, version 3 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using Seqeron.Mcp.MultipleAlignment.Tools;
using Seqeron.Mcp.Tests.Contracts;
using NUnit.Framework;

namespace Seqeron.Mcp.MultipleAlignment.Tests;

/// <summary>Pins the Seqeron.Mcp.MultipleAlignment server's public MCP tool surface (assertions in <see cref="ToolRegistrationContractTestsBase"/>).</summary>
[TestFixture]
public sealed class ToolRegistrationContractTests : ToolRegistrationContractTestsBase
{
    protected override Type ToolsType => typeof(MultipleAlignmentTools);

    protected override string[] ExpectedToolNames { get; } =
    {
        "clustalw_command",
        "clustalw_guide_tree",
        "clustalw_pairwise_align",
        "clustalw_profile_align",
        "clustalw_quick_pairwise_align",
        "clustalw_sequences_to_profile",
        "kalign_command",
        "multiple_align_clustalw",
        "multiple_align_consistency",
        "multiple_align_iterative",
        "multiple_align_kalign",
        "multiple_align_poa",
        "multiple_align_probcons",
        "multiple_align_progressive",
        "multiple_align_tcoffee",
        "poa_command",
        "probcons_command",
        "tcoffee_blast_msa_library",
        "tcoffee_evaluate",
        "tcoffee_guide_tree",
        "tcoffee_ktup_similarity",
        "tcoffee_lalign_top_alignments",
        "tcoffee_primary_library",
        "tcoffee_slow_pair_align",
        "tcoffee_tblastx_library",
    };
}
