// SPDX-License-Identifier: GPL-3.0-or-later
// Wrappers of the copyleft multiple-alignment assemblies (Kalign 3.3.5, GPL-3.0-or-later and others).
// Part of Seqeron.Mcp.MultipleAlignment, licensed under the GNU General Public License, version 3 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using NUnit.Framework;
using Seqeron.Mcp.MultipleAlignment.Tools;
using Seqeron.Mcp.Tests.Contracts;

namespace Seqeron.Mcp.MultipleAlignment.Tests;

/// <summary>The documents of every MultipleAlignment tool (<c>docs/mcp/tools/multiplealignment</c>) equal the tool metadata of the code.</summary>
[TestFixture]
public sealed class ToolDocsContractTests : ToolDocsContractTestsBase
{
    protected override Type ToolsType => typeof(MultipleAlignmentTools);

    protected override string DocsServer => "multiplealignment";
}
