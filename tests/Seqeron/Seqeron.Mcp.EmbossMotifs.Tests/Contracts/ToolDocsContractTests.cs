// SPDX-License-Identifier: GPL-2.0-or-later
// Tests of the MCP tools of the port of EMBOSS 6.6.0 (the GNU General Public License, version 2 or later).
// Part of Seqeron.Mcp.EmbossMotifs, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using NUnit.Framework;
using Seqeron.Mcp.EmbossMotifs.Tools;
using Seqeron.Mcp.Tests.Contracts;

namespace Seqeron.Mcp.EmbossMotifs.Tests;

/// <summary>The documents of every EmbossMotifs tool (<c>docs/mcp/tools/embossmotifs</c>) equal the tool metadata of the code.</summary>
[TestFixture]
public sealed class ToolDocsContractTests : ToolDocsContractTestsBase
{
    protected override Type ToolsType => typeof(EmbossMotifsTools);

    protected override string DocsServer => "embossmotifs";
}
