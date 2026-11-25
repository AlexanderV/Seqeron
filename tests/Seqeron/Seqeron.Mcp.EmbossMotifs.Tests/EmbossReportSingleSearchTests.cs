// SPDX-License-Identifier: GPL-2.0-or-later
// Tests of the MCP tools of the port of EMBOSS 6.6.0 (the GNU General Public License, version 2 or later).
// Part of Seqeron.Mcp.EmbossMotifs, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using System.Reflection;
using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.EmbossMotifs.Tests;

/// <summary>
/// B09 A6-10: the report mode of <c>find_tandem_repeats_emboss</c>, <c>fuzzpro</c>, <c>fuzztran</c> and <c>find_pattern_hits</c> runs the search ONCE (it used to run the plain tool and then the library search again for the report). The scans poll the
/// ambient cancellation seam once per window / unit: with <c>reportFormat</c> the same number of polls must happen as without it (the seam is internal to the alignment assembly, reached by reflection).
/// </summary>
[TestFixture]
public class EmbossReportSingleSearchTests
{
    private static int Polls(Action work)
    {
        var type = typeof(global::Seqeron.Genomics.Alignment.SequenceAligner).Assembly.GetType("Seqeron.Genomics.Alignment.AlignmentCancellation")!;
        int polls = 0;
        using var scope = (IDisposable)type.GetMethod("OnPoll", BindingFlags.Public | BindingFlags.NonPublic | BindingFlags.Static)!.Invoke(null, [(Action)(() => polls++)])!;
        work();
        return polls;
    }

    private const string Protein = "MKVLAAGIVGLLLAQKVLMKVLAAGIVGLLLAQKVLMKVLAAGIVGLLLAQKVL";
    private const string Dna = "ATGGCCTTAGGATCCATGAAATAGGCCTTAGGATCCATGAAATAGGATCCATGAAATAGGCCTTAGGATCC";
    private const string Tandem = "GGTACCACACACACACACACGTTGCATGCAATGCAATGCAATGCAATGCATTAGC";

    [Test]
    public void Fuzzpro_Report_RunsTheSearchOnce()
    {
        int plain = Polls(() => EmbossMotifsTools.Fuzzpro(Protein, "K-V-L"));
        int report = Polls(() => EmbossMotifsTools.Fuzzpro(Protein, "K-V-L", reportFormat: "gff", sequenceName: "p"));
        Assert.That(plain, Is.GreaterThan(0), "the search polls");
        Assert.That(report, Is.EqualTo(plain));
        var r = EmbossMotifsTools.Fuzzpro(Protein, "K-V-L", reportFormat: "excel", sequenceName: "p");
        Assert.That(r.Hits, Has.Length.EqualTo(EmbossMotifsTools.Fuzzpro(Protein, "K-V-L").Hits.Length));
        Assert.That(r.Report, Does.StartWith("SeqName\tStart\tEnd\tScore\tStrand\tPattern\tMismatch\np\t"));
    }

    [Test]
    public void Fuzztran_Report_RunsTheSearchOnce()
    {
        int plain = Polls(() => EmbossMotifsTools.Fuzztran(Dna, "M-A", "6"));
        int report = Polls(() => EmbossMotifsTools.Fuzztran(Dna, "M-A", "6", reportFormat: "gff", sequenceName: "d"));
        Assert.That(plain, Is.GreaterThan(0));
        Assert.That(report, Is.EqualTo(plain));
        var r = EmbossMotifsTools.Fuzztran(Dna, "M-A", "6", reportFormat: "gff", sequenceName: "d");
        Assert.That(r.Hits, Has.Length.EqualTo(EmbossMotifsTools.Fuzztran(Dna, "M-A", "6").Hits.Length));
        Assert.That(r.Report, Does.StartWith("##gff-version 3\n##sequence-region d 1 "));
    }
}
