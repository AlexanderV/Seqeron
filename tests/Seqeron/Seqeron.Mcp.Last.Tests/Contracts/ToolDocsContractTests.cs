// SPDX-License-Identifier: GPL-3.0-or-later
// Wrapper of the LAST 1256 port (GPL-3.0-or-later).
// Part of Seqeron.Mcp.Last, licensed under the GNU General Public License, version 3 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using Seqeron.Genomics.Last;
using NUnit.Framework;
using Seqeron.Mcp.Last.Tools;
using Seqeron.Mcp.Tests.Contracts;

namespace Seqeron.Mcp.Last.Tests;

/// <summary>The documents of every Last tool (<c>docs/mcp/tools/last</c>) equal the tool metadata of the code.</summary>
[TestFixture]
public sealed class ToolDocsContractTests : ToolDocsContractTestsBase
{
    protected override Type ToolsType => typeof(LastTools);

    protected override string DocsServer => "last";
}
