// SPDX-License-Identifier: GPL-2.0-or-later
// This file is part of Seqeron.Mcp.EmbossPairwiseAlignment, licensed under the terms stated in LICENSE (this directory).
using ModelContextProtocol.Server;

namespace Seqeron.Mcp.EmbossPairwiseAlignment.Tools;

/// <summary>MCP tools: the EMBOSS 6.6.0 pairwise alignment programs (needle, water, stretcher, supermatcher, merger, est2genome, diffseq, ...).</summary>
[McpServerToolType]
public partial class EmbossPairwiseAlignmentTools
{
    // Utility holder for static MCP tools; never instantiated (S1118).
    private EmbossPairwiseAlignmentTools() { }
}
