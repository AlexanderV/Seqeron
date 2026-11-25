// SPDX-License-Identifier: GPL-2.0-or-later
// Tests of the MCP tools of the port of EMBOSS 6.6.0 (the GNU General Public License, version 2 or later).
// Part of Seqeron.Mcp.EmbossPairwiseAlignment, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using Seqeron.Mcp.EmbossPairwiseAlignment.Tools;
using Seqeron.Mcp.Tests.Contracts;
using NUnit.Framework;

namespace Seqeron.Mcp.EmbossPairwiseAlignment.Tests;

/// <summary>Pins the Seqeron.Mcp.EmbossPairwiseAlignment server's public MCP tool surface (assertions in <see cref="ToolRegistrationContractTestsBase"/>).</summary>
[TestFixture]
public sealed class ToolRegistrationContractTests : ToolRegistrationContractTestsBase
{
    protected override Type ToolsType => typeof(EmbossPairwiseAlignmentTools);

    protected override string[] ExpectedToolNames { get; } =
    {
        "alignment_distance_matrix",
        "diffseq",
        "dotmatcher",
        "est2genome",
        "matcher_align",
        "megamerger",
        "merger",
        "multiple_alignment_statistics",
        "needle_align",
        "needleall_align",
        "plotcon",
        "showalign",
        "stretcher_align",
        "supermatcher_align",
        "water_align",
        "wordfinder",
    };
}
