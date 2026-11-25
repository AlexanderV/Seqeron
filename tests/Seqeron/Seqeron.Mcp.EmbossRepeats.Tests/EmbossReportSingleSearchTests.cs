// SPDX-License-Identifier: GPL-2.0-or-later
// Tests of the MCP tools of the port of EMBOSS 6.6.0 (the GNU General Public License, version 2 or later).
// Part of Seqeron.Mcp.EmbossRepeats, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using System.Reflection;
using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.EmbossRepeats.Tests;

/// <summary>
/// B09 A6-10: the report mode of <c>find_tandem_repeats_emboss</c>, <c>fuzzpro</c>, <c>fuzztran</c> and <c>find_pattern_hits</c> runs the search ONCE (it used to run the plain tool and then the library search again for the report). The scans poll the
/// ambient cancellation seam once per window / unit: with <c>reportFormat</c> the same number of polls must happen as without it (the seam is internal to the alignment assembly, reached by reflection).
/// </summary>
[TestFixture]
public class EmbossReportSingleSearchTests
{

    [TestCase("etandem")]
    [TestCase("equicktandem")]
    public void TandemRepeats_Report_RunsTheSearchOnce(string program)
    {
        int plain = Polls(() => EmbossRepeatsTools.FindTandemRepeatsEmboss(Tandem, program, 2, 6, 8));
        int report = Polls(() => EmbossRepeatsTools.FindTandemRepeatsEmboss(Tandem, program, 2, 6, 8, reportFormat: "excel", sequenceName: "t"));
        Assert.That(plain, Is.GreaterThan(0));
        Assert.That(report, Is.EqualTo(plain));
        var r = EmbossRepeatsTools.FindTandemRepeatsEmboss(Tandem, program, 2, 6, 8, reportFormat: "excel", sequenceName: "t");
        Assert.That(r.Repeats.Length, Is.EqualTo(EmbossRepeatsTools.FindTandemRepeatsEmboss(Tandem, program, 2, 6, 8).Repeats.Length));
        Assert.That(r.Report, Does.StartWith("SeqName\tStart\tEnd\tScore\tStrand\t"));
    }

    [Test]
    public void FindPatternHits_Report_RunsTheSearchOnce_AndTheFourToolsShareTheFormatCheck()
    {
        int plain = Polls(() => EmbossMotifsTools.FindPatternHits(Dna, "ATG"));
        int report = Polls(() => EmbossMotifsTools.FindPatternHits(Dna, "ATG", reportFormat: "gff", sequenceName: "d"));
        Assert.That(report, Is.EqualTo(plain));
        // one helper: a format the program refuses for the sequence type (embl / genbank for proteins, pir / swiss for nucleotides, diffseq: 7 tags) and an unknown name are refused identically by the four tools
        Assert.Throws<ArgumentException>(() => EmbossMotifsTools.Fuzzpro(Protein, "K-V", reportFormat: "genbank"));
        Assert.Throws<ArgumentException>(() => EmbossMotifsTools.FindPatternHits(Dna, "ATG", reportFormat: "pir"));
        Assert.Throws<ArgumentException>(() => EmbossMotifsTools.Fuzztran(Dna, "M", reportFormat: "swiss"));
        Assert.Throws<ArgumentException>(() => EmbossRepeatsTools.FindTandemRepeatsEmboss(Tandem, reportFormat: "diffseq"));
        foreach (var call in new Action[]
        {
            () => EmbossMotifsTools.FindPatternHits(Dna, "ATG", reportFormat: "nosuch"),
            () => EmbossMotifsTools.Fuzzpro(Protein, "K-V", reportFormat: "nosuch"),
            () => EmbossMotifsTools.Fuzztran(Dna, "M", reportFormat: "nosuch"),
            () => EmbossRepeatsTools.FindTandemRepeatsEmboss(Tandem, reportFormat: "nosuch"),
        })
        {
            Assert.That(Assert.Throws<ArgumentException>(() => call())!.ParamName, Is.EqualTo("reportFormat"));
        }

        // the report of a ranged search is written (A5-8b): the header names the range
        Assert.That(EmbossMotifsTools.FindPatternHits(Dna, "ATG", reportFormat: "table", sequenceQualifiers: "-sbegin 2", sequenceName: "d").Report, Does.Contain("# Sequence: d     from: 2   to: "));
        Assert.That(EmbossMotifsTools.Fuzzpro(Protein, "K-V", reportFormat: "-rformat simple -rdesshow", sequenceQualifiers: "-sbegin 2", sequenceName: "p").Report, Does.Contain("# Report_format: simple"));
    }

    private const string Tandem = "GGTACCACACACACACACACGTTGCATGCAATGCAATGCAATGCAATGCATTAGC";

    private static int Polls(Action work)
    {
        var type = typeof(global::Seqeron.Genomics.Alignment.SequenceAligner).Assembly.GetType("Seqeron.Genomics.Alignment.AlignmentCancellation")!;
        int polls = 0;
        using var scope = (IDisposable)type.GetMethod("OnPoll", BindingFlags.Public | BindingFlags.NonPublic | BindingFlags.Static)!.Invoke(null, [(Action)(() => polls++)])!;
        work();
        return polls;
    }

    private const string Dna = "ATGGCCTTAGGATCCATGAAATAGGCCTTAGGATCCATGAAATAGGATCCATGAAATAGGCCTTAGGATCC";

    private const string Protein = "MKVLAAGIVGLLLAQKVLMKVLAAGIVGLLLAQKVLMKVLAAGIVGLLLAQKVL";
}
