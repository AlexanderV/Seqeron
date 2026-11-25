// SPDX-License-Identifier: GPL-3.0-only
// Wrapper of the DECIPHER ConsensusSequence port (GPL-3.0-only).
// Part of Seqeron.Mcp.AlignmentConsensus, licensed under the GNU General Public License, version 3 (only) (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using Seqeron.Genomics.AlignmentConsensus;
using NUnit.Framework;
using Seqeron.Mcp.AlignmentConsensus.Tools;
using Seqeron.Mcp.Tests.Contracts;

namespace Seqeron.Mcp.AlignmentConsensus.Tests;

/// <summary>The documents of every AlignmentConsensus tool (<c>docs/mcp/tools/alignmentconsensus</c>) equal the tool metadata of the code.</summary>
[TestFixture]
public sealed class ToolDocsContractTests : ToolDocsContractTestsBase
{
    protected override Type ToolsType => typeof(AlignmentConsensusTools);

    protected override string DocsServer => "alignmentconsensus";
}
