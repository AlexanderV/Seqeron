// SPDX-License-Identifier: GPL-2.0-or-later
// Part of Seqeron.Mcp.OligoThermodynamics (Primer3-derived primer / probe design); licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the server project).
using NUnit.Framework;
using Seqeron.Mcp.OligoThermodynamics.Tools;
using Seqeron.Mcp.Tests.Contracts;

namespace Seqeron.Mcp.OligoThermodynamics.Tests;

/// <summary>The documents of every OligoThermodynamics tool (<c>docs/mcp/tools/oligothermodynamics</c>) equal the tool metadata of the code.</summary>
[TestFixture]
public sealed class ToolDocsContractTests : ToolDocsContractTestsBase
{
    protected override Type ToolsType => typeof(OligoThermodynamicsTools);

    protected override string DocsServer => "oligothermodynamics";
}
