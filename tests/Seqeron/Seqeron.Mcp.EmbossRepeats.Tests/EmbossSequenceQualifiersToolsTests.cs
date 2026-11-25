// SPDX-License-Identifier: GPL-2.0-or-later
// Tests of the MCP tools of the port of EMBOSS 6.6.0 (the GNU General Public License, version 2 or later).
// Part of Seqeron.Mcp.EmbossRepeats, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.EmbossRepeats.Tests;

/// <summary>B09 A5-7: the <c>sequenceQualifiers</c> parameter of fuzznuc / fuzzpro / fuzztran / etandem / equicktandem / seqmatchall / distmat. Values are what the installed EMBOSS 6.6.0 prints.</summary>
[TestFixture]
public class EmbossSequenceQualifiersToolsTests
{

    [Test]
    public void Etandem_IgnoresTheRange_Equicktandem_UsesIt()
    {
        const string seq = "GGTACCACACACACACACACGTTGCATGCAATGCAATGCAATGCAATGCATTAGC";
        // etandem et.fa -minrepeat 2 -maxrepeat 6 -threshold 8 -sbegin 10 -send 40 prints the same repeats as without
        var whole = EmbossRepeatsTools.FindTandemRepeatsEmboss(seq, "etandem", 2, 6, 8);
        var cut = EmbossRepeatsTools.FindTandemRepeatsEmboss(seq, "etandem", 2, 6, 8, sequenceQualifiers: "-sbegin 10 -send 40");
        Assert.That(cut.Repeats.Select(r => (r.Position, r.Length)).ToArray(), Is.EqualTo(whole.Repeats.Select(r => (r.Position, r.Length)).ToArray()));
        // equicktandem -maxrepeat 6 -threshold 8 -sbegin 10 -send 40: the first line of the program is  9  10  20  2  5
        var quick = EmbossRepeatsTools.FindTandemRepeatsEmboss(seq, "equicktandem", 2, 6, 8, sequenceQualifiers: "-sbegin 10 -send 40");
        // the whole output of the program: -origfile of equicktandem et.fa -maxrepeat 6 -threshold 8 -sbegin 10 -send 40 is "9 10 20 2 5" and "10 26 40 5 3" (score, start, end, size, count); without the range "13 6 20 2 7" and "20 26 50 5 5"
        Assert.That(quick.Repeats.Select(r => (r.Score, r.Position + 1, r.Position + r.Length, r.RepeatSize, r.RepeatCount)).ToArray(), Is.EqualTo(new[] { (9, 10, 20, 2, 5), (10, 26, 40, 5, 3) }));
        var quickWhole = EmbossRepeatsTools.FindTandemRepeatsEmboss(seq, "equicktandem", 2, 6, 8);
        Assert.That(quickWhole.Repeats.Select(r => (r.Score, r.Position + 1, r.Position + r.Length, r.RepeatSize, r.RepeatCount)).ToArray(), Is.EqualTo(new[] { (13, 6, 20, 2, 7), (20, 26, 50, 5, 5) }));
    }

    [Test]
    public void Seqmatchall_Reverse_NumbersFromTheEnd()
    {
        var r = EmbossRepeatsTools.SeqMatchAll(["TCCACTAGCGAGAGGTACCCTGATAAGGGGGAGGCAGCGTAGGTTA", "TCTCTTAGCGGGAACTAGCCTGAGCTGATCGGGGGTCTAGCCGCCTAGGGTTAGG"], 3, sequenceQualifiers: "-sreverse");
        // the complete match list of `seqmatchall -sequence s.fa -wordsize 3 -sreverse -aformat match` (installed 6.6.0): 42 matches, every one on the reverse strand of both sequences, numbered from the end (start > end)
        var expected = new (int Start1, int End1, int Start2, int End2, int Length)[] { (3, 1, 51, 49, 3), (5, 1, 7, 3, 5), (5, 3, 22, 20, 3), (7, 4, 4, 1, 4), (7, 4, 10, 7, 4), (7, 5, 18, 16, 3), (7, 5, 40, 38, 3), (7, 5, 50, 48, 3), (11, 8, 49, 46, 4), (11, 9, 17, 15, 3), (11, 9, 33, 31, 3), (11, 9, 39, 37, 3), (15, 13, 3, 1, 3), (15, 13, 9, 7, 3), (16, 14, 34, 32, 3), (18, 15, 46, 43, 4), (18, 16, 8, 6, 3), (18, 16, 25, 23, 3), (19, 16, 25, 22, 4), (19, 17, 8, 6, 3), (19, 17, 46, 44, 3), (20, 16, 25, 21, 5), (20, 17, 24, 21, 4), (20, 18, 23, 21, 3), (20, 18, 46, 44, 3), (21, 18, 9, 6, 4), (21, 19, 3, 1, 3), (27, 23, 31, 27, 5), (28, 24, 37, 33, 5), (28, 26, 12, 10, 3), (33, 31, 7, 5, 3), (33, 31, 22, 20, 3), (34, 32, 3, 1, 3), (34, 32, 9, 7, 3), (35, 33, 34, 32, 3), (37, 35, 34, 32, 3), (40, 38, 33, 31, 3), (41, 37, 50, 46, 5), (41, 39, 4, 2, 3), (42, 38, 19, 15, 5), (42, 39, 11, 8, 4), (43, 38, 42, 37, 6) };
        var actual = r.ReportedMatches!.OrderBy(m => m.Start1).ThenBy(m => m.End1).ThenBy(m => m.Start2).ThenBy(m => m.End2).ThenBy(m => m.Length).Select(m => (m.Start1, m.End1, m.Start2, m.End2, m.Length)).ToArray();
        Assert.That(actual, Is.EqualTo(expected));
        Assert.That(r.ReportedMatches!.All(m => m.Reversed), Is.True);
    }}
