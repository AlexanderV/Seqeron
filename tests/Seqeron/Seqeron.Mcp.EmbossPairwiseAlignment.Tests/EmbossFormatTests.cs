// SPDX-License-Identifier: GPL-2.0-or-later
// Tests of the MCP tools of the port of EMBOSS 6.6.0 (the GNU General Public License, version 2 or later).
// Part of Seqeron.Mcp.EmbossPairwiseAlignment, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;
using Seqeron.Mcp.EmbossPairwiseAlignment.Tools;

namespace Seqeron.Mcp.EmbossPairwiseAlignment.Tests;

/// <summary>stretcher_align / write_alignment: EMBOSS 6.6 stretcher and the -aformat output formats (audit A2-27 / A2-28, B09 F120 / F121). The tools delegate to EmbossPairwisePrograms.StretcherAlign and EmbossAlignmentWriter; the oracle comparisons live in the library tests.</summary>
[TestFixture]
public class EmbossFormatTests
{

    [Test]
    public void StretcherAlign_EqualsTheLibrary_AndTheDocumentedEmbossOutput()
    {
        var tool = EmbossPairwiseAlignmentTools.StretcherAlign(A, B);
        var lib = global::Seqeron.Genomics.Emboss.PairwiseAlignment.EmbossPairwisePrograms.StretcherAlign(A, B);
        Assert.Multiple(() =>
        {
            Assert.That(tool.AlignedSequence1, Is.EqualTo(lib.AlignedSequence1));
            Assert.That(tool.AlignedSequence2, Is.EqualTo(lib.AlignedSequence2));
            Assert.That(tool.Score, Is.EqualTo(140), "stretcher a.fa b.fa: Score: 140");
            Assert.That(tool.AlignedSequence1, Is.EqualTo("ACGTACGTAGCTAGCTAGCATCGA-TCGATCGAT"));
        });
    }

    [Test]
    public void StretcherAlign_Validation()
    {
        Assert.Throws<ArgumentException>(() => EmbossPairwiseAlignmentTools.StretcherAlign(null!, A));
        Assert.Throws<ArgumentException>(() => EmbossPairwiseAlignmentTools.StretcherAlign(A, B, matrix: "NOSUCH"));
        Assert.Throws<ArgumentOutOfRangeException>(() => EmbossPairwiseAlignmentTools.StretcherAlign(A, B, gapOpen: -1));
    }

    private const string A = "ACGTACGTAGCTAGCTAGCATCGATCGATCGAT";

    private const string B = "ACGTACGTAGCTAGATAGCATCGAATCGATCGAT";
}
