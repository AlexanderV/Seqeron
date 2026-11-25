// SPDX-License-Identifier: GPL-2.0-or-later
// Tests of the MCP tools of the port of EMBOSS 6.6.0 (the GNU General Public License, version 2 or later).
// Part of Seqeron.Mcp.EmbossMotifs, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.EmbossMotifs.Tests;

/// <summary>
/// B09 A2-30 .. A2-36: the MCP parameters / tools added with the both-strand repeats, canonical tandem rotation, fuzznuc pattern language, common-region modes and strands,
/// k-mismatch common region and IUPAC input. Expected values are those locked by the reference-checked library tests (MUMmer 3.23 repeat-match / mummer, EMBOSS fuzznuc 6.6.0)
/// or derived by hand from the definitions.
/// </summary>
[TestFixture]
public class StrandPatternAndRepeatToolsTests
{

    [Test]
    public void FindPatternHits_SequenceAmbiguity_ReadsAmbiguityCodesAsFuzznucDoes()
    {
        // fuzznuc 6.6.0 on ARGRYN: pattern R -> bases 1-4 (+) and 5 (-); the default reading matches only the plain A and G (A3-12a, F146)
        var r = EmbossMotifsTools.FindPatternHits("ARGRYN", "R", bothStrands: true, sequenceAmbiguity: true);
        Assert.That(r.Hits.Select(h => (h.Position, h.Strand)), Is.EquivalentTo(new[] { (0, "forward"), (1, "forward"), (2, "forward"), (3, "forward"), (4, "reverse") }));
        var plain = EmbossMotifsTools.FindPatternHits("ARGRYN", "R");
        Assert.That(plain.Hits.Select(h => h.Position), Is.EqualTo(new[] { 0, 2 }));
    }

    [Test]
    public void FindPatternHits_FuzznucPatterns()
    {
        // fuzznuc -pattern "A(1,3)C" on AAAACAAAAT -> 1-based 4..5, 3..5, 2..5
        var r = EmbossMotifsTools.FindPatternHits("AAAACAAAAT", "A(1,3)C");
        Assert.That(r.Hits.Select(h => (h.Position, h.Length, h.Strand, h.Mismatches)), Is.EquivalentTo(new[] { (3, 2, "forward", 0), (2, 3, "forward", 0), (1, 4, "forward", 0) }));
        Assert.That(r.Patterns["pattern1"], Is.EqualTo("A(1,3)C"));
        // fuzznuc -complement Y: "<ACG" on AAACGT -> '-' 4..6
        var anchored = EmbossMotifsTools.FindPatternHits("AAACGT", "<ACG", bothStrands: true);
        Assert.That(anchored.Hits.Select(h => (h.Position, h.Length, h.Strand)), Is.EqualTo(new[] { (3, 3, "reverse") }));
        // several named patterns from a pattern file
        var many = EmbossMotifsTools.FindPatternHits("ACGTACGTTTGAACC", patternFile: ">p1\nACGT\n>p2\nTTGA\n");
        Assert.That(many.Hits.Select(h => (h.Pattern, h.Position)), Is.EquivalentTo(new[] { ("p1", 0), ("p1", 4), ("p2", 8) }));
        Assert.That(many.Patterns.Keys, Is.EqualTo(new[] { "p1", "p2" }));
        // mismatches, N in the sequence is a mismatch
        Assert.That(EmbossMotifsTools.FindPatternHits("ACNTAAAT", "ACNT", maxMismatches: 1).Hits.Select(h => h.Position), Is.EquivalentTo(new[] { 0, 4 }));
        Assert.Throws<ArgumentException>(() => EmbossMotifsTools.FindPatternHits("ACGT"));
        Assert.Throws<ArgumentException>(() => EmbossMotifsTools.FindPatternHits("ACGT", "AC", "AC"));
        Assert.Throws<ArgumentException>(() => EmbossMotifsTools.FindPatternHits("ACGT", "A(2,1)"));
        Assert.Throws<ArgumentOutOfRangeException>(() => EmbossMotifsTools.FindPatternHits("ACGT", "AC", maxMismatches: -1));
        Assert.Throws<ArgumentException>(() => EmbossMotifsTools.FindPatternHits("", "AC"));
    }

    [Test]
    public void FindMotif_And_KnownMotifs_AcceptIupacSequences()
    {
        Assert.That(EmbossMotifsTools.FindMotif("ACANACAGA", "AGA").Positions, Is.EqualTo(new[] { 6 }));
        Assert.That(EmbossMotifsTools.FindMotif("ACANAC", "ANA").Positions, Is.Empty);
        Assert.That(EmbossMotifsTools.FindKnownMotifs("ACNACGACG", new[] { "ACG" }).Matches["ACG"], Is.EqualTo(new[] { 3, 6 }));
        Assert.Throws<ArgumentException>(() => EmbossMotifsTools.FindMotif("ACGU", "A"));
    }}
