// SPDX-License-Identifier: GPL-2.0-or-later
// Tests of the MCP tools of the port of EMBOSS 6.6.0 (the GNU General Public License, version 2 or later).
// Part of Seqeron.Mcp.EmbossSequenceFormats, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using Seqeron.Mcp.EmbossSequenceFormats.Tools;
using Seqeron.Mcp.Tests.Contracts;
using NUnit.Framework;

namespace Seqeron.Mcp.EmbossSequenceFormats.Tests;

/// <summary>Pins the Seqeron.Mcp.EmbossSequenceFormats server's public MCP tool surface (assertions in <see cref="ToolRegistrationContractTestsBase"/>).</summary>
[TestFixture]
public sealed class ToolRegistrationContractTests : ToolRegistrationContractTestsBase
{
    protected override Type ToolsType => typeof(EmbossSequenceFormatsTools);

    protected override string[] ExpectedToolNames { get; } =
    {
        "cutseq",
        "dbtell",
        "degapseq",
        "descseq",
        "drfinddata",
        "entret",
        "extractseq",
        "infoseq",
        "listor",
        "makenucseq",
        "makeprotseq",
        "maskambignuc",
        "maskambigprot",
        "maskseq",
        "msbar",
        "newseq",
        "noreturn",
        "nospace",
        "notab",
        "notseq",
        "nthseq",
        "nthseqset",
        "pasteseq",
        "prettyseq",
        "revseq",
        "seqcount",
        "seqret",
        "seqretsplit",
        "seqxref",
        "showdb",
        "showorf",
        "shuffleseq",
        "skipseq",
        "splitter",
        "transeq",
        "trimest",
        "trimseq",
        "trimspace",
        "union",
        "wordcount",
        "write_alignment",
        "yank",
    };
}
