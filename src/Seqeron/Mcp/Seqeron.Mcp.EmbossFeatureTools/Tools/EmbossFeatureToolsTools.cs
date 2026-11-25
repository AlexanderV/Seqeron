// SPDX-License-Identifier: GPL-2.0-or-later
// This file is part of Seqeron.Mcp.EmbossFeatureTools, licensed under the terms stated in LICENSE (this directory).
using ModelContextProtocol.Server;

namespace Seqeron.Mcp.EmbossFeatureTools.Tools;

/// <summary>MCP tools: the EMBOSS 6.6.0 feature programs (showfeat, extractfeat, twofeat, featcopy, ...).</summary>
[McpServerToolType]
public partial class EmbossFeatureToolsTools
{
    // Utility holder for static MCP tools; never instantiated (S1118).
    private EmbossFeatureToolsTools() { }
}
