// SPDX-License-Identifier: GPL-2.0-or-later
// Tests of the MCP tools of the port of EMBOSS 6.6.0 (the GNU General Public License, version 2 or later).
// Part of Seqeron.Mcp.EmbossMotifs, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.EmbossMotifs.Tests;

/// <summary>
/// Tests for the <c>find_motif</c> MCP tool.
/// Expected positions from the exact-match / overlapping-occurrence definition
/// (suffix-tree order is unspecified), NOT the wrapper output.
/// </summary>
[TestFixture]
public class FindMotifTests
{
    [Test]
    public void FindMotif_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => EmbossMotifsTools.FindMotif("ATGATG", "ATG"));
        Assert.Throws<ArgumentException>(() => EmbossMotifsTools.FindMotif("", "ATG"));
        Assert.Throws<ArgumentException>(() => EmbossMotifsTools.FindMotif(null!, "ATG"));
        Assert.Throws<ArgumentException>(() => EmbossMotifsTools.FindMotif("ATGATG", ""));
        Assert.Throws<ArgumentException>(() => EmbossMotifsTools.FindMotif("ATGATG", null!));
        Assert.Throws<ArgumentException>(() => EmbossMotifsTools.FindMotif("XYZ123", "ATG"));
    }

    [Test]
    public void FindMotif_Binding_InvokesSuccessfully()
    {
        // ATG occurs at 0 and 3 (order from suffix tree is unspecified).
        var atg = EmbossMotifsTools.FindMotif("ATGATG", "ATG").Positions;
        Assert.That(atg, Is.EquivalentTo(new[] { 0, 3 }));

        // Overlapping AA in AAAA -> {0,1,2}.
        var aa = EmbossMotifsTools.FindMotif("AAAA", "AA").Positions;
        Assert.That(aa, Is.EquivalentTo(new[] { 0, 1, 2 }));

        // Absent motif -> empty.
        var none = EmbossMotifsTools.FindMotif("ATGATG", "CCC").Positions;
        Assert.That(none, Is.Empty);
    }

    // B05 request (B09 F18): positions are ascending, equal to the MotifFinder-backed find_exact_motif tool.
    [Test]
    public void FindMotif_Positions_AreAscending_AndEqualFindExactMotif()
    {
        foreach (var (seq, motif) in new[] { ("ACACACACAC", "AC"), ("GATATATGCATATACTT", "ATAT"), ("AAAAAA", "AA"), ("CATGCATGCATG", "ATG") })
        {
            var positions = EmbossMotifsTools.FindMotif(seq, motif).Positions;
            Assert.That(positions, Is.Ordered.Ascending, seq);
            Assert.That(positions, Is.EqualTo(AnalysisTools.FindExactMotif(seq, motif).Positions), seq);
        }

        // Rosalind SUBS sample (1-based 2 4 10).
        Assert.That(EmbossMotifsTools.FindMotif("GATATATGCATATACTT", "ATAT").Positions, Is.EqualTo(new[] { 1, 3, 9 }));
    }
}
