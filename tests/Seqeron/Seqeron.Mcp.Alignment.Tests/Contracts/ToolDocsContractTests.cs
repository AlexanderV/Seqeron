using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;
using Seqeron.Mcp.Tests.Contracts;

namespace Seqeron.Mcp.Alignment.Tests;

/// <summary>Audit A5-5 (B09 F170): the documents of every Alignment tool (<c>docs/mcp/tools/alignment</c>) equal the tool metadata of the code (description, parameters, types, defaults, parameter descriptions, the Input Schema table).</summary>
[TestFixture]
public sealed class ToolDocsContractTests : ToolDocsContractTestsBase
{
    protected override Type ToolsType => typeof(AlignmentTools);

    protected override string DocsServer => "alignment";
}
