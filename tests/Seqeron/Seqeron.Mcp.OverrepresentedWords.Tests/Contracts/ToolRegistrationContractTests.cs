// SPDX-License-Identifier: AGPL-3.0-only
// Tests of the MCP server Seqeron.Mcp.OverrepresentedWords (RSAT oligo-analysis / dyad-analysis); licensed under the GNU Affero General Public License, version 3 (see LICENSE in the server project).
using Seqeron.Mcp.OverrepresentedWords.Tools;
using Seqeron.Mcp.Tests.Contracts;
using NUnit.Framework;

namespace Seqeron.Mcp.OverrepresentedWords.Tests;

/// <summary>Pins the Seqeron.Mcp.OverrepresentedWords server's public MCP tool surface (assertions in <see cref="ToolRegistrationContractTestsBase"/>).</summary>
[TestFixture]
public sealed class ToolRegistrationContractTests : ToolRegistrationContractTestsBase
{
    protected override Type ToolsType => typeof(OverrepresentedWordsTools);

    protected override string[] ExpectedToolNames { get; } =
    {
        "dyad_analysis",
        "oligo_analysis",
        "shared_motifs_significance",
    };
}
