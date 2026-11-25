// SPDX-License-Identifier: GPL-2.0-or-later
// This file is part of Seqeron.Mcp.EmbossRepeats, licensed under the terms stated in LICENSE (this directory).
using ModelContextProtocol.Server;

namespace Seqeron.Mcp.EmbossRepeats.Tools;

/// <summary>MCP tools: the EMBOSS 6.6.0 word-match and tandem-repeat programs (wordmatch, seqmatchall, etandem, equicktandem).</summary>
[McpServerToolType]
public partial class EmbossRepeatsTools
{
    // Utility holder for static MCP tools; never instantiated (S1118).
    private EmbossRepeatsTools() { }
}
