// SPDX-License-Identifier: GPL-2.0-or-later
// Tests of the MCP tools of the port of EMBOSS 6.6.0 (the GNU General Public License, version 2 or later).
// Part of Seqeron.Mcp.EmbossMotifs, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.EmbossMotifs.Tests;

/// <summary>
/// Tests for the <c>find_known_motifs</c> MCP tool.
/// Expected values from the exact set-matching definition (ascending positions,
/// misses omitted, overlapping reported), NOT the wrapper output.
/// </summary>
[TestFixture]
public class FindKnownMotifsTests
{
    [Test]
    public void FindKnownMotifs_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => EmbossMotifsTools.FindKnownMotifs("ATGATGCC", new[] { "ATG" }));
        Assert.Throws<ArgumentException>(() => EmbossMotifsTools.FindKnownMotifs("", new[] { "ATG" }));
        Assert.Throws<ArgumentException>(() => EmbossMotifsTools.FindKnownMotifs(null!, new[] { "ATG" }));
        Assert.Throws<ArgumentException>(() => EmbossMotifsTools.FindKnownMotifs("XYZ", new[] { "ATG" }));
        Assert.Throws<ArgumentException>(() => EmbossMotifsTools.FindKnownMotifs("ATGATGCC", null!));
    }

    [Test]
    public void FindKnownMotifs_Binding_InvokesSuccessfully()
    {
        // "ATGATGCC": ATG at 0,3; CC at 6; GG absent (omitted).
        var m = EmbossMotifsTools.FindKnownMotifs("ATGATGCC", new[] { "ATG", "CC", "GG" }).Matches;
        Assert.Multiple(() =>
        {
            Assert.That(m.ContainsKey("ATG"), Is.True);
            Assert.That(m["ATG"], Is.EqualTo(new[] { 0, 3 }));
            Assert.That(m.ContainsKey("CC"), Is.True);
            Assert.That(m["CC"], Is.EqualTo(new[] { 6 }));
            Assert.That(m.ContainsKey("GG"), Is.False);
        });

        // Overlapping AA in AAAA -> 0,1,2.
        var overlap = EmbossMotifsTools.FindKnownMotifs("AAAA", new[] { "AA" }).Matches;
        Assert.That(overlap["AA"], Is.EqualTo(new[] { 0, 1, 2 }));
    }

    // EMBOSS fuzznuc 6.6.0 (-pmismatch 0 -complement Y) on GAATTCAAAGAATTCTTTATATAAAGTTCAAGAAC (B09 F19):
    // GAATTC '+' 1-6, 10-15 and '-' 1-6, 10-15; TATAWAW '+' 18-24, '-' 17-23 (forward-strand 1-based coordinates).
    [Test]
    public void FindKnownMotifs_BothStrandsAndIupac_MatchFuzznuc()
    {
        const string seq = "GAATTCAAAGAATTCTTTATATAAAGTTCAAGAAC";

        var eco = EmbossMotifsTools.FindKnownMotifs(seq, new[] { "GAATTC" }, bothStrands: true);
        Assert.Multiple(() =>
        {
            Assert.That(eco.Matches["GAATTC"], Is.EqualTo(new[] { 0, 9 }));
            Assert.That(eco.ReverseMatches!["GAATTC"], Is.EqualTo(new[] { 0, 9 }));
        });

        var tata = EmbossMotifsTools.FindKnownMotifs(seq, new[] { "TATAWAW" }, bothStrands: true, iupac: true);
        Assert.Multiple(() =>
        {
            Assert.That(tata.Matches["TATAWAW"], Is.EqualTo(new[] { 17 }));
            Assert.That(tata.ReverseMatches!["TATAWAW"], Is.EqualTo(new[] { 16 }));
        });

        // Defaults keep the exact forward-only contract: IUPAC symbols are literal, no reverse hits.
        var plain = EmbossMotifsTools.FindKnownMotifs(seq, new[] { "TATAWAW", "GAATTC" });
        Assert.Multiple(() =>
        {
            Assert.That(plain.Matches.ContainsKey("TATAWAW"), Is.False);
            Assert.That(plain.Matches["GAATTC"], Is.EqualTo(new[] { 0, 9 }));
            Assert.That(plain.ReverseMatches, Is.Empty);
        });
        Assert.Throws<ArgumentException>(() => EmbossMotifsTools.FindKnownMotifs(seq, new[] { "AXG" }, iupac: true));
    }

    // EMBOSS fuzznuc 6.6.0 (-pmismatch k -complement Y) on GAATTCAAGATTCGGATTTTGAATTTT (B09 A1-20), 0-based forward-strand starts:
    // GAATTC k=1: + 0 (0 mm), 20 (1 mm); - 0 (0 mm), 20 (1 mm). RAATTC k=2 (IUPAC): + 0:0 7:1 13:2 14:2 20:1 21:2; - 0:0 7:2 13:1 14:1 20:0 21:2.
    [Test]
    public void FindKnownMotifs_Mismatches_MatchFuzznuc_WithDistances()
    {
        const string seq = "GAATTCAAGATTCGGATTTTGAATTTT";

        var eco = EmbossMotifsTools.FindKnownMotifs(seq, new[] { "GAATTC" }, bothStrands: true, mismatches: 1);
        Assert.Multiple(() =>
        {
            Assert.That(eco.Matches["GAATTC"], Is.EqualTo(new[] { 0, 20 }));
            Assert.That(eco.MatchMismatches!["GAATTC"], Is.EqualTo(new[] { 0, 1 }));
            Assert.That(eco.ReverseMatches!["GAATTC"], Is.EqualTo(new[] { 0, 20 }));
            Assert.That(eco.ReverseMismatches!["GAATTC"], Is.EqualTo(new[] { 0, 1 }));
        });

        var r = EmbossMotifsTools.FindKnownMotifs(seq, new[] { "RAATTC" }, bothStrands: true, iupac: true, mismatches: 2);
        Assert.Multiple(() =>
        {
            Assert.That(r.Matches["RAATTC"], Is.EqualTo(new[] { 0, 7, 13, 14, 20, 21 }));
            Assert.That(r.MatchMismatches!["RAATTC"], Is.EqualTo(new[] { 0, 1, 2, 2, 1, 2 }));
            Assert.That(r.ReverseMatches!["RAATTC"], Is.EqualTo(new[] { 0, 7, 13, 14, 20, 21 }));
            Assert.That(r.ReverseMismatches!["RAATTC"], Is.EqualTo(new[] { 0, 2, 1, 1, 0, 2 }));
        });
    }

    [Test]
    public void FindKnownMotifs_Mismatches_DefaultKeepsExactContract_AndRejectsNegative()
    {
        var exact = EmbossMotifsTools.FindKnownMotifs("GAATTCAAGATTCGGATTTTGAATTTT", new[] { "GAATTC" });
        Assert.Multiple(() =>
        {
            Assert.That(exact.Matches["GAATTC"], Is.EqualTo(new[] { 0 }));
            Assert.That(exact.MatchMismatches, Is.Null, "distance arrays only with mismatches > 0");
            Assert.That(exact.ReverseMismatches, Is.Null);
        });
        Assert.Throws<ArgumentOutOfRangeException>(() => EmbossMotifsTools.FindKnownMotifs("ACGT", new[] { "AC" }, mismatches: -1));
    }
}
