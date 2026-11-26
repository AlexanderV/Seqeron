// SPDX-License-Identifier: GPL-3.0-or-later
// This file is part of Seqeron.Mcp.Last, licensed under the terms stated in LICENSE (this directory).
using ModelContextProtocol.Server;

namespace Seqeron.Mcp.Last.Tools;

/// <summary>MCP tools: the script mafft with the LAST programs (lastdb, lastal) behind it - separately licensed (GPL-3.0-or-later).</summary>
[McpServerToolType]
public partial class LastTools
{
    // Utility holder for static MCP tools; never instantiated (S1118).
    private LastTools() { }
}
