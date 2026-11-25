// SPDX-License-Identifier: GPL-2.0-or-later
// This file is part of Seqeron.Mcp.EmbossOpenReadingFrames, licensed under the terms stated in LICENSE (this directory).
using ModelContextProtocol.Server;

namespace Seqeron.Mcp.EmbossOpenReadingFrames.Tools;

/// <summary>MCP tools: the EMBOSS 6.6.0 open reading frame programs (getorf, plotorf, tcode, checktrans).</summary>
[McpServerToolType]
public partial class EmbossOpenReadingFramesTools
{
    // Utility holder for static MCP tools; never instantiated (S1118).
    private EmbossOpenReadingFramesTools() { }
}
