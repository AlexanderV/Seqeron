// SPDX-License-Identifier: GPL-2.0-or-later
// Tests of the MCP tools of the port of EMBOSS 6.6.0 (the GNU General Public License, version 2 or later).
// Part of Seqeron.Mcp.EmbossPairwiseAlignment, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using System.Text.Json.Nodes;
using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.EmbossPairwiseAlignment.Tests;

/// <summary>
/// B09 WP67 (F308, R14-13): the MCP tools <c>alignment_summary_info</c> (Bio.Align.AlignInfo.SummaryInfo of Biopython 1.79) and <c>msa_conservation</c> (Capra and Singh 2007) and the extended columns of
/// <c>multiple_alignment_statistics</c>. Expected values: the test suite of Biopython 1.79 (<c>Tests/test_AlignInfo.py</c>) and an independent implementation of the formulas of the paper.
/// </summary>
[TestFixture]
public class AlignInfoToolsTests
{

    [Test]
    public void MultipleAlignmentStatistics_ColumnsCarryTheShannonEntropyAndTheInformationContent()
    {
        var s = EmbossPairwiseAlignmentTools.MultipleAlignmentStatistics(["AAC-", "AAG-", "AAC-", "ACC-"]);
        Assert.Multiple(() =>
        {
            Assert.That(s.Columns[0].ShannonEntropy, Is.EqualTo(0.0));
            Assert.That(s.Columns[1].ShannonEntropy, Is.EqualTo(0.8112781244591328).Within(1e-15), "A x3, C x1: -(0.75 log2 0.75 + 0.25 log2 0.25)");
            Assert.That(s.Columns[0].InformationContent, Is.EqualTo(2.0).Within(1e-15), "a conserved nucleotide column against 0.25: log2(4) bits");
        });
    }}
