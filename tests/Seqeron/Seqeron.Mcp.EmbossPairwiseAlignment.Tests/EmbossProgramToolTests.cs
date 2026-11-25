// SPDX-License-Identifier: GPL-2.0-or-later
// Tests of the MCP tools of the port of EMBOSS 6.6.0 (the GNU General Public License, version 2 or later).
// Part of Seqeron.Mcp.EmbossPairwiseAlignment, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;
using Seqeron.Mcp.EmbossPairwiseAlignment.Tools;

namespace Seqeron.Mcp.EmbossPairwiseAlignment.Tests;

/// <summary>supermatcher_align / needleall_align / ksw2_extend and the `format` parameter of the pairwise and multiple-alignment tools (audit A2-27b / A2-28c / A2-29b, B09 F124 - F126). The tools delegate to
/// EmbossPairwisePrograms.Supermatcher / NeedleAll, Ksw2Aligner and AlignmentFormats / EmbossAlignmentWriter; the oracle comparisons live in the library tests.</summary>
[TestFixture]
public class EmbossProgramToolTests
{

    [Test]
    public void SupermatcherAlign_IdenticalSequences_OneWordSeededAlignment()
    {
        var r = EmbossPairwiseAlignmentTools.SupermatcherAlign([Q], [Q]);
        Assert.Multiple(() =>
        {
            Assert.That(r.Count, Is.EqualTo(1));
            Assert.That(r.Alignments[0].AlignedSequence1, Is.EqualTo(Q));
            Assert.That(r.Alignments[0].AlignedSequence2, Is.EqualTo(Q));
            Assert.That(r.Alignments[0].Score, Is.EqualTo(140.0), "28 matches x EDNAFULL 5");
            Assert.That(r.Alignments[0].Name1, Is.EqualTo("query1"));
            Assert.That(r.ErrorText, Is.Empty);
        });
    }

    [Test]
    public void SupermatcherAlign_NoSharedWord_ReportsTheErrorFileLine()
    {
        var r = EmbossPairwiseAlignmentTools.SupermatcherAlign(["AAAAAAAAAAAAAAAAAAAA"], ["CCCCCCCCCCCCCCCCCCCC"], ["q"], ["t"]);
        Assert.That(r.Count, Is.EqualTo(0));
        Assert.That(r.ErrorText, Is.EqualTo("No wordmatch start points for q vs t. No alignment\n"));
    }

    [Test]
    public void SupermatcherAlign_Validation()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => EmbossPairwiseAlignmentTools.SupermatcherAlign([Q], [Q], wordLength: 2));
        Assert.Throws<ArgumentOutOfRangeException>(() => EmbossPairwiseAlignmentTools.SupermatcherAlign([Q], [Q], width: 0));
        Assert.Throws<ArgumentException>(() => EmbossPairwiseAlignmentTools.SupermatcherAlign([], [Q]));
        Assert.Throws<ArgumentException>(() => EmbossPairwiseAlignmentTools.SupermatcherAlign([Q], [Q], matrix: "NOSUCH"));
    }

    [Test]
    public void NeedleallAlign_PairsInTheProgramOrder_AndTheErrorLine()
    {
        var r = EmbossPairwiseAlignmentTools.NeedleallAlign(["ACGT", "AAAA"], ["ACGT", "CCCC"]);
        Assert.Multiple(() =>
        {
            Assert.That(r.Alignments.Select(a => a.Name1 + "/" + a.Name2), Is.EqualTo(new[] { "a1/b1", "a2/b1", "a1/b2" }));
            Assert.That(r.Alignments[0].Score, Is.EqualTo(20.0));
            Assert.That(r.ErrorText, Does.Contain("is less than minimum score(0.0) for sequences a2 vs b2"));
        });
    }

    private const string Q = "ACGTACGTACGGTTAACCGGTTAAACCC";
}
