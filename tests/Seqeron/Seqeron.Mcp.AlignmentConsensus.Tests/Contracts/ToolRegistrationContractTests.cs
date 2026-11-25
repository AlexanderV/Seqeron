// SPDX-License-Identifier: GPL-3.0-only
// Wrapper of the DECIPHER ConsensusSequence port (GPL-3.0-only).
// Part of Seqeron.Mcp.AlignmentConsensus, licensed under the GNU General Public License, version 3 (only) (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using Seqeron.Genomics.AlignmentConsensus;
using Seqeron.Mcp.AlignmentConsensus.Tools;
using Seqeron.Mcp.Tests.Contracts;
using NUnit.Framework;

namespace Seqeron.Mcp.AlignmentConsensus.Tests;

/// <summary>Pins the Seqeron.Mcp.AlignmentConsensus server's public MCP tool surface (assertions in <see cref="ToolRegistrationContractTestsBase"/>).</summary>
[TestFixture]
public sealed class ToolRegistrationContractTests : ToolRegistrationContractTestsBase
{
    protected override Type ToolsType => typeof(AlignmentConsensusTools);

    protected override string[] ExpectedToolNames { get; } =
    {
        "generate_decipher_consensus",
    };
}
