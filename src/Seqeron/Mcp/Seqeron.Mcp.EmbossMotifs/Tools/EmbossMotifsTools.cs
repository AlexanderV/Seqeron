// SPDX-License-Identifier: GPL-2.0-or-later
// This file is part of Seqeron.Mcp.EmbossMotifs, licensed under the terms stated in LICENSE (this directory).
using ModelContextProtocol.Server;

namespace Seqeron.Mcp.EmbossMotifs.Tools;

/// <summary>MCP tools: the EMBOSS 6.6.0 pattern and motif programs (fuzzpro, fuzztran, patmatmotifs, patmatdb, prosextract, marscan, consambig).</summary>
[McpServerToolType]
public partial class EmbossMotifsTools
{
    // Utility holder for static MCP tools; never instantiated (S1118).
    private EmbossMotifsTools() { }
}
