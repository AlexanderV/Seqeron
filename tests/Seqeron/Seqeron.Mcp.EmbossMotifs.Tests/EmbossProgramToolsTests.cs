// SPDX-License-Identifier: GPL-2.0-or-later
// Tests of the MCP tools of the port of EMBOSS 6.6.0 (the GNU General Public License, version 2 or later).
// Part of Seqeron.Mcp.EmbossMotifs, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.EmbossMotifs.Tests;

/// <summary>
/// B09 A3-8 / A3-12 / A3-13: MCP tools of the EMBOSS program ports. Expected values are those printed by the installed EMBOSS 6.6.0 programs (also locked by the oracle tests of
/// the library).
/// </summary>
[TestFixture]
public class EmbossProgramToolsTests
{

    [Test]
    public void Consambig_UnionOfTheColumnCodes_LowerCaseAtGaps()
    {
        // consambig 6.6.0: AC-GTN / AC-GTA / ARGGT- -> AVgGTn; D|N = B, I|L = J for proteins; an all-gap protein column is the NUL character (B09 A3-13d, F147)
        Assert.That(EmbossMotifsTools.Consambig(["AC-GTN", "AC-GTA", "ARGGT-"]).Consensus, Is.EqualTo("AVgGTn"));
        Assert.That(EmbossMotifsTools.Consambig(["DI", "NL"], "protein").Consensus, Is.EqualTo("BJ"));
        Assert.That(EmbossMotifsTools.Consambig(["MK-LA", "MK-LD"]).Consensus, Is.EqualTo("MK\0LX"));
        Assert.Throws<ArgumentException>(() => EmbossMotifsTools.Consambig(["ACGT"]));
        Assert.Throws<ArgumentException>(() => EmbossMotifsTools.Consambig(["ACGT", "ACGT"], "rna"));
    }

    [Test]
    public void Fuzzpro_WithMismatches_AsTheProgram()
    {
        // fuzzpro 6.6.0 -pmismatch 1 'H-H-K' on MKHHKLLHHKAAHHKMMHHRPPQQ (A5-8c2, F195): the last window HHR has one mismatch (score 2)
        var r = EmbossMotifsTools.Fuzzpro("MKHHKLLHHKAAHHKMMHHRPPQQ", "H-H-K", mismatch: 1);
        Assert.That(r.Hits.Select(h => (h.Position, h.Length, h.Mismatches)), Is.EqualTo(new[] { (2, 3, 0), (7, 3, 0), (12, 3, 0), (17, 3, 1) }));
        Assert.Throws<ArgumentException>(() => EmbossMotifsTools.Fuzzpro("MKHHK", "[HK]-[HK]", mismatch: 2));      // type 6: the mismatches must be below the real pattern length
    }

    [Test]
    public void Fuzzpro_ReportsEveryWindowOfTheProteinPattern()
    {
        // fuzzpro -pattern "[AG]-x(2)-[IV]" on MKVLAAGIVGLLLAQPSTAEEKV prints the 1-based windows 5-8 and 6-9
        var r = EmbossMotifsTools.Fuzzpro("MKVLAAGIVGLLLAQPSTAEEKV", "[AG]-x(2)-[IV]");
        Assert.That(r.Pattern, Is.EqualTo("[AG]x(2)[IV]"));
        Assert.That(r.Hits.Select(h => (h.Position, h.Length)).ToArray(), Is.EqualTo(new[] { (4, 4), (5, 4) }));
        Assert.Throws<ArgumentException>(() => EmbossMotifsTools.Fuzzpro("MKV", "[BA]"));
        Assert.Throws<ArgumentException>(() => EmbossMotifsTools.Fuzzpro("", "A"));
    }

    [Test]
    public void Fuzztran_SearchesTheTranslatedFrames()
    {
        // fuzztran -pattern KV on ATGAAAGTTCTGGCT (frame 1 = MKVLA) prints the feature 4 9 with *start 2 *end 3 *translation KV
        var r = EmbossMotifsTools.Fuzztran("ATGAAAGTTCTGGCT", "KV");
        var h = r.Hits.Single();
        Assert.That((h.Frame, h.ProteinStart, h.ProteinEnd, h.Translation, h.Position, h.Length, h.EmbossStart, h.EmbossEnd), Is.EqualTo((1, 2, 3, "KV", 3, 6, (int?)4, (int?)9)));
        // reverse frames carry the forward-strand location and the numbers the program prints for them (S - |frame| - 3 start + 4 / + 2, A5-7)
        var rev = EmbossMotifsTools.Fuzztran("GGATTA", "S", "-1").Hits.Single();
        Assert.That((rev.Frame, rev.Translation, rev.Position, rev.Length, rev.EmbossStart, rev.EmbossEnd), Is.EqualTo((-1, "S", 0, 3, (int?)(-2), (int?)(-4))));        // protein position 2: 1 - 1 - 3 * 2 + 4 and 1 - 1 - 3 * 2 + 2
        Assert.Throws<ArgumentException>(() => EmbossMotifsTools.Fuzztran("ACGT", "A", "7"));
        Assert.Throws<ArgumentException>(() => EmbossMotifsTools.Fuzztran("ACGT", "[X]"));
    }}
