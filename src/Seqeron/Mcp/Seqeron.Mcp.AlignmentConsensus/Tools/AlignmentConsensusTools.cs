// SPDX-License-Identifier: GPL-3.0-only
// This file is part of Seqeron.Mcp.AlignmentConsensus, licensed under the terms stated in LICENSE (this directory).
using ModelContextProtocol.Server;

namespace Seqeron.Mcp.AlignmentConsensus.Tools;

/// <summary>MCP tools: alignment consensus (DECIPHER ConsensusSequence port)</summary>
[McpServerToolType]
public partial class AlignmentConsensusTools
{
    // Utility holder for static MCP tools; never instantiated (S1118).
    private AlignmentConsensusTools() { }
}
