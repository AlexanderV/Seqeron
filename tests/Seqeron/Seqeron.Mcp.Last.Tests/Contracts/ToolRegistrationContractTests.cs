// SPDX-License-Identifier: GPL-3.0-or-later
// Wrapper of the LAST 1256 port (GPL-3.0-or-later).
// Part of Seqeron.Mcp.Last, licensed under the GNU General Public License, version 3 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using Seqeron.Genomics.Last;
using Seqeron.Mcp.Last.Tools;
using Seqeron.Mcp.Tests.Contracts;
using NUnit.Framework;

namespace Seqeron.Mcp.Last.Tests;

/// <summary>Pins the Seqeron.Mcp.Last server's public MCP tool surface (assertions in <see cref="ToolRegistrationContractTestsBase"/>).</summary>
[TestFixture]
public sealed class ToolRegistrationContractTests : ToolRegistrationContractTestsBase
{
    protected override Type ToolsType => typeof(LastTools);

    protected override string[] ExpectedToolNames { get; } =
    {
        "mafft_command_last",
    };
}
