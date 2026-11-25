// SPDX-License-Identifier: GPL-2.0-or-later
// Tests of the MCP tools of the port of EMBOSS 6.6.0 (the GNU General Public License, version 2 or later).
// Part of Seqeron.Mcp.EmbossMotifs, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.EmbossMotifs.Tests;

/// <summary>B09 A5-7: the <c>sequenceQualifiers</c> parameter of fuzznuc / fuzzpro / fuzztran / etandem / equicktandem / seqmatchall / distmat. Values are what the installed EMBOSS 6.6.0 prints.</summary>
[TestFixture]
public class EmbossSequenceQualifiersToolsTests
{
    private const string S = "ACGTACGTTAGCCGATTAGGCATCGATCGGATCGAAACCCTTTGGG";

    [Test]
    public void FindPatternHits_Reverse_ReportsTheMinusStrandInOriginalCoordinates()
    {
        // fuzznuc s.fa -pattern CGATC -sreverse -rformat excel: s 30 34 - and s 25 29 - (original numbering); without -sreverse: 24 28 +
        var r = EmbossMotifsTools.FindPatternHits(S, "CGATC", sequenceQualifiers: "-sreverse");
        Assert.That(r.Hits.Select(h => (h.Position + 1, h.Position + h.Length, h.Strand)).OrderBy(h => h.Item1).ToArray(),
            Is.EqualTo(new[] { (25, 29, "reverse"), (30, 34, "reverse") }));
        var plain = EmbossMotifsTools.FindPatternHits(S, "CGATC");
        Assert.That(plain.Hits.Select(h => (h.Position + 1, h.Strand)).ToArray(), Is.EqualTo(new[] { (24, "forward") }));
    }

    [Test]
    public void FindPatternHits_Range_AddsTheBeginOffset()
    {
        var r = EmbossMotifsTools.FindPatternHits(S, "CGATC", sequenceQualifiers: "-sbegin 5 -send 40");
        Assert.That(r.Hits.Single().Position, Is.EqualTo(23));
    }

    [Test]
    public void FindPatternHits_NumberedEqualsAndAbbreviatedForms_AreResolvedLikeTheProgram()
    {
        // A6-2: fuzznuc s.fa -pattern CGATC -sbegin1 5 -send1 40 == -sbegin 5 -send 40; also -sbegin1=5, the unique abbreviations -sbeg1 / -sen1 / --SEND=40
        var expected = EmbossMotifsTools.FindPatternHits(S, "CGATC", sequenceQualifiers: "-sbegin 5 -send 40").Hits.Select(h => (h.Position, h.Length, h.Strand)).ToArray();
        Assert.That(expected, Is.EqualTo(new[] { (23, 5, "forward") }));
        foreach (string text in new[] { "-sbegin1 5 -send1 40", "-sbegin1=5 -send1=40", "-sbeg1 5 -sen1 40", "-SBEGIN=5 --send 40", "-sbegin1 5 -send 40" })
            Assert.That(EmbossMotifsTools.FindPatternHits(S, "CGATC", sequenceQualifiers: text).Hits.Select(h => (h.Position, h.Length, h.Strand)).ToArray(), Is.EqualTo(expected), text);
        // a range that starts after the hit finds nothing; -se is ambiguous between -send and -sequence in fuzznuc, -sbegin2 names a sequence it does not have
        Assert.That(EmbossMotifsTools.FindPatternHits(S, "CGATC", sequenceQualifiers: "-sbeg1 30").Hits, Is.Empty);
        Assert.That(Assert.Throws<ArgumentException>(() => EmbossMotifsTools.FindPatternHits(S, "CGATC", sequenceQualifiers: "-se 40"))!.Message, Does.Contain("ambiguous qualifier 'se' (send,sequence)"));
        Assert.Throws<ArgumentException>(() => EmbossMotifsTools.FindPatternHits(S, "CGATC", sequenceQualifiers: "-sbegin2 5"));
        Assert.Throws<ArgumentException>(() => EmbossMotifsTools.FindPatternHits(S, "CGATC", sequenceQualifiers: "-pattern ACG"));
    }

    [Test]
    public void Fuzzpro_Range_CountsFromTheOriginalStart()
    {
        var plain = EmbossMotifsTools.Fuzzpro("MKVLAAGIVGLLLAQ", "K-V");
        var ranged = EmbossMotifsTools.Fuzzpro("MKVLAAGIVGLLLAQ", "K-V", sequenceQualifiers: "-sbegin 2");
        Assert.That(ranged.Hits.Single().Position, Is.EqualTo(plain.Hits.Single().Position));
        Assert.Throws<ArgumentException>(() => EmbossMotifsTools.Fuzzpro("MKVLAAGIVGLLLAQ", "K-V", sequenceQualifiers: "-snucleotide"));
    }

    [Test]
    public void ReportFormat_ReturnsTheProgramsReportText()
    {
        // B09 A5-8: fuzznuc a -pattern ACGTA -rformat gff / excel
        var gff = EmbossMotifsTools.FindPatternHits("ACGTAGGG", "ACGTA", reportFormat: "gff", sequenceName: "a");
        Assert.That(gff.Report, Does.StartWith("##gff-version 3\n##sequence-region a 1 8\n#!Date "));
        Assert.That(gff.Report, Does.EndWith("a\tfuzznuc\tnucleotide_motif\t1\t5\t5\t+\t.\tID=a.1;note=*pat pattern:ACGTA\n"));
        var excel = EmbossMotifsTools.FindPatternHits("ACGTAGGG", "ACGTA", reportFormat: "excel", sequenceName: "a");
        Assert.That(excel.Report, Is.EqualTo("SeqName\tStart\tEnd\tScore\tStrand\tPattern\tMismatch\na\t1\t5\t5\t+\tpattern:ACGTA\t.\n"));
        Assert.That(EmbossMotifsTools.FindPatternHits("ACGTAGGG", "ACGTA").Report, Is.Null);
        var table = EmbossMotifsTools.FindPatternHits("ACGTAGGG", "ACGTA", reportFormat: "seqtable", sequenceName: "a");          // A5-8b: every -rformat of the program
        Assert.That(table.Report, Does.Contain("  Start     End  Strand Pattern       Mismatch Sequence\n      1       5       + pattern:ACGTA        . ACGTA\n"));
        Assert.Throws<ArgumentException>(() => EmbossMotifsTools.FindPatternHits("ACGTAGGG", "ACGTA", reportFormat: "nosuch"));
        var pro = EmbossMotifsTools.Fuzzpro("ACDKLEFGHIKLMNPQRSTVWY", "K-L", reportFormat: "excel", sequenceName: "p");
        Assert.That(pro.Report, Is.EqualTo("SeqName\tStart\tEnd\tScore\tStrand\tPattern\tMismatch\np\t4\t5\t2\t+\tpattern:K-L\t.\np\t11\t12\t2\t+\tpattern:K-L\t.\n"));
    }
}
