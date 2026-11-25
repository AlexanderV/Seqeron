// SPDX-License-Identifier: GPL-2.0-or-later
// Tests of the MCP tools of the copyleft MCP servers (see LICENSE of the server).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using NUnit.Framework;
using Seqeron.Mcp.EmbossSequenceFormats.Tools;

namespace Seqeron.Mcp.EmbossPairwiseAlignment.Tests;

/// <summary>The EMBOSS servers install the EMBOSS alignment formats into the format parameter of the core alignment tools at start-up; the tests do the same once.</summary>
[SetUpFixture]
public sealed class AlignmentFormatsSetUp
{
    [OneTimeSetUp]
    public void InstallEmbossAlignmentFormats() => EmbossSequenceFormatsTools.InstallAlignmentFormats();
}
