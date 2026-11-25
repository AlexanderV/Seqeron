// SPDX-License-Identifier: GPL-2.0-or-later
// This file is part of Seqeron.Mcp.OligoThermodynamics, licensed under the terms stated in LICENSE (this directory).
using ModelContextProtocol.Server;

namespace Seqeron.Mcp.OligoThermodynamics.Tools;

/// <summary>MCP tools: Primer3-style primer and probe design and ntthal thermodynamic structure analysis (hairpins, dimers, beacons)</summary>
[McpServerToolType]
public partial class OligoThermodynamicsTools
{
    // Utility holder for static MCP tools; never instantiated (S1118).
    private OligoThermodynamicsTools() { }
}
