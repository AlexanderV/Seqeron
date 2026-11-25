// SPDX-License-Identifier: GPL-2.0-or-later
// Tests of the MCP tools of the port of EMBOSS 6.6.0 (the GNU General Public License, version 2 or later).
// Part of Seqeron.Mcp.EmbossFeatureTools, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using Seqeron.Mcp.EmbossFeatureTools.Tools;
using Seqeron.Mcp.Tests.Contracts;
using NUnit.Framework;

namespace Seqeron.Mcp.EmbossFeatureTools.Tests;

/// <summary>Pins the Seqeron.Mcp.EmbossFeatureTools server's public MCP tool surface (assertions in <see cref="ToolRegistrationContractTestsBase"/>).</summary>
[TestFixture]
public sealed class ToolRegistrationContractTests : ToolRegistrationContractTestsBase
{
    protected override Type ToolsType => typeof(EmbossFeatureToolsTools);

    protected override string[] ExpectedToolNames { get; } =
    {
        "coderet",
        "extractfeat",
        "featcopy",
        "featmerge",
        "featreport",
        "feattext",
        "maskfeat",
        "showfeat",
        "twofeat",
    };
}
