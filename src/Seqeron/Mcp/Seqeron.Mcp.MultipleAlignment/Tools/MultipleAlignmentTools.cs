// SPDX-License-Identifier: GPL-3.0-or-later
// This file is part of Seqeron.Mcp.MultipleAlignment, licensed under the terms stated in LICENSE (this directory).
using ModelContextProtocol.Server;

namespace Seqeron.Mcp.MultipleAlignment.Tools;

/// <summary>MCP tools: multiple sequence alignment programs (Kalign, POA, ClustalW, T-Coffee ports)</summary>
[McpServerToolType]
public partial class MultipleAlignmentTools
{
    // Utility holder for static MCP tools; never instantiated (S1118).
    private MultipleAlignmentTools() { }
}
