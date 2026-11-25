// SPDX-License-Identifier: AGPL-3.0-only
// Tests of the MCP server Seqeron.Mcp.OverrepresentedWords (RSAT oligo-analysis / dyad-analysis); licensed under the GNU Affero General Public License, version 3 (see LICENSE in the server project).
using NUnit.Framework;
using Seqeron.Mcp.OverrepresentedWords.Tools;
using Seqeron.Mcp.Tests.Contracts;

namespace Seqeron.Mcp.OverrepresentedWords.Tests;

/// <summary>The documents of every OverrepresentedWords tool (<c>docs/mcp/tools/overrepresentedwords</c>) equal the tool metadata of the code.</summary>
[TestFixture]
public sealed class ToolDocsContractTests : ToolDocsContractTestsBase
{
    protected override Type ToolsType => typeof(OverrepresentedWordsTools);

    protected override string DocsServer => "overrepresentedwords";
}
