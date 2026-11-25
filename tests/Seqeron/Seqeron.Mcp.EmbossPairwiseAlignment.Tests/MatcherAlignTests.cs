// SPDX-License-Identifier: GPL-2.0-or-later
// Tests of the MCP tools of the port of EMBOSS 6.6.0 (the GNU General Public License, version 2 or later).
// Part of Seqeron.Mcp.EmbossPairwiseAlignment, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;

namespace Seqeron.Mcp.EmbossPairwiseAlignment.Tests;

[TestFixture]
public class MatcherAlignTests
{
    [Test]
    public void TopN_TieOrderSim_ReproducesEmbossMatcher_WhereWatermanEggertOrderDiffers()
    {
        // topn_oracle.tsv line 46 (B09 F145, audit A2-28d): matcher -alternatives 6 -gapopen 8 -gapextend 2 on EDNAFULL reports 55,38,21,16,13,11.
        const string a = "GCTCTAAATTTTGTAGCTCTAATT", b = "TTCGCTCTAAATTTC";
        var sim = EmbossPairwiseAlignmentTools.MatcherAlign(a, b, maxAlignments: 6, affine: true, match: 5, mismatch: -4, gapOpen: -6, gapExtend: -2, tieOrder: "Sim");
        Assert.That(sim.Alignments.Select(x => x.Score), Is.EqualTo(new[] { 55, 38, 21, 16, 13, 11 }));
        var we = AlignmentTools.LocalAlignTopN(a, b, maxAlignments: 6, affine: true, match: 5, mismatch: -4, gapOpen: -6, gapExtend: -2);
        Assert.That(we.Alignments.Select(x => x.Score), Is.Not.EqualTo(sim.Alignments.Select(x => x.Score)));
        Assert.Throws<ArgumentException>(() => EmbossPairwiseAlignmentTools.MatcherAlign(a, b, tieOrder: "nonsense"));
    }
}
