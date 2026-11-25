// SPDX-License-Identifier: GPL-2.0-or-later
// Tests of the MCP tools of the port of EMBOSS 6.6.0 (the GNU General Public License, version 2 or later).
// Part of Seqeron.Mcp.EmbossRepeats, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.EmbossRepeats.Tests;

/// <summary>
/// B09 A3-8 / A3-12 / A3-13: MCP tools of the EMBOSS program ports. Expected values are those printed by the installed EMBOSS 6.6.0 programs (also locked by the oracle tests of
/// the library).
/// </summary>
[TestFixture]
public class EmbossProgramToolsTests
{

    [Test]
    public void Wordmatch_FeatureOutputFormats_AsTheProgram()
    {
        // wordmatch 6.6.0 -offormat4 embl -offormat5 genbank on the same pair (A5-8c3d, F225): the feature files in the EMBL feature table and the GenBank FEATURES layout
        var r = EmbossRepeatsTools.Wordmatch(">p\nACGTACGTAC\n", ">g\nTTACGTACGTACTT\n", 6, featureQualifiers: "-offormat4 embl -offormat5 genbank", featureFile1: "pa.out", featureFile2: "gb.out");
        Assert.Multiple(() =>
        {
            Assert.That(r.FeatureFile1, Is.EqualTo("FH   Key             Location/Qualifiers\nFH\nFT   misc_feature    4..10\nFT                   /note=\"g\"\nFT   misc_feature    2..10\nFT                   /note=\"g\"\n"));
            Assert.That(r.FeatureFile2, Is.EqualTo("FEATURES             Location/Qualifiers\n     misc_feature    2..8\n                     /note=\"p\"\n     misc_feature    4..12\n                     /note=\"p\"\n"));
            Assert.That((r.FeatureFileName1, r.FeatureFileName2), Is.EqualTo(("pa.out", "gb.out")));
        });
    }

    [Test]
    public void Wordmatch_SequenceQualifiers_AsTheProgram()
    {
        // wordmatch 6.6.0 -sbegin1 3 -send1 18 -sreverse2 -sbegin2 4 -send2 40 (A5-8c3e, F239): the set is cut to residues 3..18, the scanned sequences are reversed (their -sbegin2 / -send2 are ignored by the program)
        var r = EmbossRepeatsTools.Wordmatch(">p1 desc one\nACGTACGTTAGCATCGATCG\n>p2\nTTAGCATCGGGATTACGAT\n", ">g1 genome\nGGACGTACGTTAGCATCGATCGAAATTAGCATCGGGATTACGATTT\n>g2\nCCCTTAGCATCGATCGATCGATGGG\n", 5,
            featureQualifiers: "-sbegin1 3 -send1 18 -sreverse2 -sbegin2 4 -send2 40");
        Assert.Multiple(() =>
        {
            Assert.That(r.AllMatches, Is.EqualTo(7));
            Assert.That(r.Alignments, Does.Contain("     5 p1              +       12..16       g1              -       28..24\n     6 p1              +       11..16       g1              -       32..27\n     6 p1              +        1..6        g1              -       44..39\n"));
            Assert.That(r.FeatureFile1, Does.Contain("##sequence-region p1 1 16\n"));
        });
    }

    [Test]
    public void Wordmatch_RabinKarpMatches_AsTheProgram()
    {
        // wordmatch 6.6.0 on p = ACGTACGTAC against g = TTACGTACGTACTT, -wordsize 6 (A2-28d, F198): the two words TACGTA / CGTACG are extended to the right, the match of the second starts inside the first
        var r = EmbossRepeatsTools.Wordmatch(">p\nACGTACGTAC\n", ">g\nTTACGTACGTACTT\n", 6);
        Assert.That(r.Matches, Is.EqualTo(new[] { new WordmatchHitItem("p", 3, "g", 1, 7), new WordmatchHitItem("p", 1, "g", 3, 9) }));
        Assert.That((r.PatternCount, r.AllMatches, r.SumLength), Is.EqualTo((4, 2L, 16L)));
        Assert.That(r.Alignments, Does.Contain("     7 p               +        4..10       g               +        2..8\n     9 p               +        2..10       g               +        4..12\n"));
        Assert.That(r.Log, Does.Contain("Average match length: 8.00").And.Contain("TACGTA :            1             1              7.00"));
        Assert.Throws<ArgumentException>(() => EmbossRepeatsTools.Wordmatch("", ">g\nAC\n"));
    }

    [Test]
    public void FindTandemRepeatsEmboss_Etandem_ReportsScoreIdentityAndConsensus()
    {
        // etandem -minrepeat 2 -maxrepeat 8 -threshold 10: "33 61 96 3 12 100.0 cag" and "26 137 172 6 6 94.4 acgtta"
        var r = EmbossRepeatsTools.FindTandemRepeatsEmboss(TandemSeq, "etandem", 2, 8, 10);
        Assert.That(r.Program, Is.EqualTo("etandem"));
        Assert.That(r.Repeats.Select(h => (h.Position + 1, h.Position + h.Length, h.Score, h.RepeatSize, h.RepeatCount, h.Identity, h.Consensus)),
            Is.EqualTo(new[] { (61, 96, 33, 3, 12, (double?)100.0, "cag"), (137, 172, 26, 6, 6, (double?)94.4, "acgtta") }));
    }

    [Test]
    public void FindTandemRepeatsEmboss_Equicktandem_ReportsDiagonalScores()
    {
        // equicktandem -maxrepeat 6 -threshold 10 on (ACG)x8: "21 1 24 3 8"
        var r = EmbossRepeatsTools.FindTandemRepeatsEmboss("ACGACGACGACGACGACGACGACG", "EquickTandem", maxRepeat: 6, threshold: 10);
        var h = r.Repeats.Single();
        Assert.That((h.Position, h.Length, h.Score, h.RepeatSize, h.RepeatCount, h.Identity, h.Consensus), Is.EqualTo((0, 24, 21, 3, 8, (double?)null, (string?)null)));
    }

    [Test]
    public void FindTandemRepeatsEmboss_RejectsBadArguments()
    {
        Assert.Throws<ArgumentException>(() => EmbossRepeatsTools.FindTandemRepeatsEmboss("ACGT", "tandem"));
        Assert.Throws<ArgumentException>(() => EmbossRepeatsTools.FindTandemRepeatsEmboss("", "etandem"));
        Assert.Throws<ArgumentOutOfRangeException>(() => EmbossRepeatsTools.FindTandemRepeatsEmboss("ACGT", "etandem", 1));
        Assert.Throws<ArgumentOutOfRangeException>(() => EmbossRepeatsTools.FindTandemRepeatsEmboss("ACGT", "equicktandem", threshold: 0));
    }

    [Test]
    public void SeqMatchAll_ReportsMaximalExactMatchesOfEveryPair()
    {
        // seqmatchall -wordsize 6 on these three sequences: "10 a 1..10 b 3..12", "6 a 1..6 b 7..12", "6 a 5..10 b 3..8" (1-based); the third is shorter than the word size
        var r = EmbossRepeatsTools.SeqMatchAll(["ACGTACGTAC", "GGACGTACGTACTT", "TTTT"], 6);
        Assert.That(r.Matches.Select(m => (m.Sequence1, m.Sequence2, m.Position1, m.Position2, m.Length)).ToArray(),
            Is.EqualTo(new[] { (0, 1, 0, 2, 10), (0, 1, 0, 6, 6), (0, 1, 4, 2, 6) }));
        Assert.Throws<ArgumentException>(() => EmbossRepeatsTools.SeqMatchAll([], 4));
        Assert.Throws<ArgumentException>(() => EmbossRepeatsTools.SeqMatchAll(["ACGT", ""], 4));
        Assert.Throws<ArgumentOutOfRangeException>(() => EmbossRepeatsTools.SeqMatchAll(["ACGT", "ACGT"], 0));
    }

    /// <summary>A4-3 (B09 F153): the EMBOSS tools read the sequence like the programs (gap characters and other non-letters dropped first; positions count the cleaned sequence).</summary>
    [Test]
    public void EmbossTools_DropGapCharactersFromThePlainSequence_LikeThePrograms()
    {
        // fuzznuc -pattern ACGT on `AC-GT.AACGT` (cleaned ACGTAACGT): hits 1-4 and 6-9
        var hits = EmbossMotifsTools.FindPatternHits("AC-GT.AACGT", "ACGT").Hits.Select(h => (h.Position, h.Length)).ToArray();
        Assert.That(hits, Is.EqualTo(new[] { (0, 4), (5, 4) }));
        // fuzztran -pattern KV on `ATG-AAA GTT.CTGGCT` is the clean-text hit
        var h = EmbossMotifsTools.Fuzztran("ATG-AAA GTT.CTGGCT", "KV").Hits.Single();
        Assert.That((h.EmbossStart, h.EmbossEnd), Is.EqualTo(((int?)4, (int?)9)));
        Assert.That(EmbossRepeatsTools.FindTandemRepeatsEmboss("ACACAC-ACACAC.ACAC", "etandem", 2, 4, 5).Repeats.Select(r => (r.Position, r.Length)).ToArray(),
            Is.EqualTo(EmbossRepeatsTools.FindTandemRepeatsEmboss("ACACACACACACACAC", "etandem", 2, 4, 5).Repeats.Select(r => (r.Position, r.Length)).ToArray()));
        // consambig / showalign: `? # +` are the unknown residue, digits and blanks vanish
        Assert.That(EmbossMotifsTools.Consambig(["ACGT#CGT", "ACGT?CGT"]).Consensus, Is.EqualTo(EmbossMotifsTools.Consambig(["ACGTNCGT", "ACGTNCGT"]).Consensus));
        Assert.That(EmbossMotifsTools.Consambig(["AC1GT CGT", "ACG_TACGT"]).Consensus, Is.EqualTo(EmbossMotifsTools.Consambig(["ACGTCGT", "ACGTACGT"]).Consensus));
        Assert.That(EmbossRepeatsTools.SeqMatchAll(["ACGT-ACGTA", "ACGT.ACGTA"], 4).Matches, Is.EqualTo(EmbossRepeatsTools.SeqMatchAll(["ACGTACGTA", "ACGTACGTA"], 4).Matches));
    }

    private const string TandemSeq = "CAGATTTTCATATTATGCAGAAAATCTACTTCGCCTGATACGAGTCGGTTATCTTCGGATCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGACTGTATAGTCCCACCTGGTGATCCTATGCTTGTGAGTACACGTTAACGTTAACGTTAACGTTAACGTTAacgtNN";
}
