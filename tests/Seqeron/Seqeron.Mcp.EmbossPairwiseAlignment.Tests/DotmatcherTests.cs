// SPDX-License-Identifier: GPL-2.0-or-later
// Tests of the MCP tools of the port of EMBOSS 6.6.0 (the GNU General Public License, version 2 or later).
// Part of Seqeron.Mcp.EmbossPairwiseAlignment, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;
using Seqeron.Mcp.EmbossPairwiseAlignment.Tools;

namespace Seqeron.Mcp.EmbossPairwiseAlignment.Tests;

/// <summary>MCP <c>dotmatcher</c> (B09 A2-28d, F152): EMBOSS dotmatcher 6.6.0 plot segments; expected values are the lines <c>-graph data</c> prints (x = the second sequence, y = the first).</summary>
[TestFixture]
public class DotmatcherTests
{
    [Test]
    public void Dotmatcher_ReportsTheProgramsSegmentsInDrawingOrder()
    {
        // dotmatcher a b -windowsize 6 -threshold 20 -graph data: "1 1 8 8", "1 5 6 10", "5 11 12 18" ... (x1 y1 x2 y2, x on b, y on a)
        var r = EmbossPairwiseAlignmentTools.Dotmatcher("ACGTACGTAACCGGTTACGTAAACCCGGGTTT", "ACGTACGGTAACCGGTTACGTAAACCCGGGTT", 6, 20);
        Assert.That(r.Segments, Has.Length.EqualTo(17));
        Assert.That(r.Segments[0], Is.EqualTo(new DotMatcherSegmentDto(1, 1, 8, 8)));
        Assert.That(r.Segments[1], Is.EqualTo(new DotMatcherSegmentDto(5, 1, 10, 6)));   // plot (x1 y1 x2 y2) = 1 5 6 10: x on the second sequence
        Assert.That(r.Segments[2], Is.EqualTo(new DotMatcherSegmentDto(11, 5, 18, 12)));
        Assert.That(EmbossPairwiseAlignmentTools.Dotmatcher("ACGTACGTAACCGGTT", "ACGTACGTAACCGGTT", 6, 20, "EDNAFULL").Segments, Is.Not.Empty);
        Assert.That(EmbossPairwiseAlignmentTools.Dotmatcher("ACGT", "ACGTACGT").Segments, Is.Empty);
        Assert.Throws<ArgumentException>(() => EmbossPairwiseAlignmentTools.Dotmatcher("ACGTACGTACGT", "ACGTACGTACGT", matrix: "NOPE"));
    }
}
