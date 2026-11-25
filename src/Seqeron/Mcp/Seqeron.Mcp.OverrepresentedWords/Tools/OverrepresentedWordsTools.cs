// SPDX-License-Identifier: AGPL-3.0-only
// This file is part of Seqeron.Mcp.OverrepresentedWords, licensed under the terms stated in LICENSE (this directory).
using ModelContextProtocol.Server;

namespace Seqeron.Mcp.OverrepresentedWords.Tools;

/// <summary>MCP tools: RSAT oligo-analysis, dyad-analysis and shared-motif significance (word over-representation statistics)</summary>
[McpServerToolType]
public partial class OverrepresentedWordsTools
{
    // Utility holder for static MCP tools; never instantiated (S1118).
    private OverrepresentedWordsTools() { }
}
