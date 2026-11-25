// SPDX-License-Identifier: GPL-2.0-or-later
// This file is part of Seqeron.Mcp.EmbossSequenceFormats, licensed under the terms stated in LICENSE (this directory).
using ModelContextProtocol.Server;

namespace Seqeron.Mcp.EmbossSequenceFormats.Tools;

/// <summary>MCP tools: the EMBOSS 6.6.0 sequence programs (seqret, union, revseq, transeq, extractseq, infoseq, ...).</summary>
[McpServerToolType]
public partial class EmbossSequenceFormatsTools
{
    // Utility holder for static MCP tools; never instantiated (S1118).
    private EmbossSequenceFormatsTools() { }
}
