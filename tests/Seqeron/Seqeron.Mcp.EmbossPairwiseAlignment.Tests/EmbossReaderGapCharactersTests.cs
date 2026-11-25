// SPDX-License-Identifier: GPL-2.0-or-later
// Tests of the MCP tools of the port of EMBOSS 6.6.0 (the GNU General Public License, version 2 or later).
// Part of Seqeron.Mcp.EmbossPairwiseAlignment, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;
using Seqeron.Mcp.EmbossPairwiseAlignment.Tools;

namespace Seqeron.Mcp.EmbossPairwiseAlignment.Tests;

/// <summary>A4-3 (B09 F153): the MCP EMBOSS tools read plain sequences like the programs: gap characters, digits, blanks and every other non-letter except <c>*</c> are dropped first; lengths and offsets of the formatted output count the cleaned sequence.</summary>
[TestFixture]
public class EmbossReaderGapCharactersTests
{
    [Test]
    public void NeedleAlign_OfGappedText_EqualsTheCleanedText_Score58AndLength11()
    {
        // needle 6.6.0: `AC-GTTAGC.TT*` vs `ACGTTAGCTT` prints Length 11 and Score 58.0
        var r = EmbossPairwiseAlignmentTools.NeedleAlign("AC-GTTAGC.TT*", "ACGTTAGCTT");
        Assert.That(r.ExactScore, Is.EqualTo(58.0));
        Assert.That(r.AlignedSequence1, Is.EqualTo("ACGTTAGCTT*"));
        Assert.That(r.AlignedSequence2, Is.EqualTo("ACGTTAGCTT-"));
        Assert.That(r.EndPosition1, Is.EqualTo(10));      // 0-based inclusive end (R14-1); the program's Length is 11
        var formatted = EmbossPairwiseAlignmentTools.NeedleAlign("AC-GTTAGC.TT*", "ACGTTAGCTT", format: "pair").Formatted;
        Assert.That(formatted, Does.Contain("# Length: 11").And.Contain("# Score: 58.0"));
        Assert.That(EmbossPairwiseAlignmentTools.WaterAlign("AC-GTTAGC.TT*", "ACGTTAGCTT").ExactScore, Is.EqualTo(58.0));
    }

    [Test]
    public void DotmatcherAndMerger_DropGapCharacters()
    {
        Assert.That(EmbossPairwiseAlignmentTools.Dotmatcher("ACGTACGT-AACCGGTT", "ACGT.ACGTAACCGGTT", 6, 20).Segments,
            Is.EqualTo(EmbossPairwiseAlignmentTools.Dotmatcher("ACGTACGTAACCGGTT", "ACGTACGTAACCGGTT", 6, 20).Segments));
        var gapped = EmbossPairwiseAlignmentTools.Merger("acgtacgtac-gtacgtt1ggcc", "gtacgtt ggccaa.ttacg");
        var clean = EmbossPairwiseAlignmentTools.Merger("acgtacgtacgtacgttggcc", "gtacgttggccaattacg");
        Assert.That(gapped.MergedSequence, Is.EqualTo(clean.MergedSequence));
    }
}
