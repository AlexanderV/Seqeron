// SPDX-License-Identifier: GPL-2.0-or-later
// Tests of the MCP tools of the port of EMBOSS 6.6.0 (the GNU General Public License, version 2 or later).
// Part of Seqeron.Mcp.EmbossMotifs, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using Seqeron.Mcp.EmbossMotifs.Tools;
using Seqeron.Mcp.Tests.Contracts;
using NUnit.Framework;

namespace Seqeron.Mcp.EmbossMotifs.Tests;

/// <summary>Pins the Seqeron.Mcp.EmbossMotifs server's public MCP tool surface (assertions in <see cref="ToolRegistrationContractTestsBase"/>).</summary>
[TestFixture]
public sealed class ToolRegistrationContractTests : ToolRegistrationContractTestsBase
{
    protected override Type ToolsType => typeof(EmbossMotifsTools);

    protected override string[] ExpectedToolNames { get; } =
    {
        "consambig",
        "emboss_report_file",
        "find_known_motifs",
        "find_motif",
        "find_pattern_hits",
        "fuzzpro",
        "fuzztran",
        "generate_emboss_consensus",
        "marscan",
        "patmatdb",
        "patmatmotifs",
        "prosextract",
    };
}
