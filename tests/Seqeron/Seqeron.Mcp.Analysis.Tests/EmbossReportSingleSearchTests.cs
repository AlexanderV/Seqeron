using System.Reflection;
using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

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
        int plain = Polls(() => AnalysisTools.Fuzzpro(Protein, "K-V-L"));
        int report = Polls(() => AnalysisTools.Fuzzpro(Protein, "K-V-L", reportFormat: "gff", sequenceName: "p"));
        Assert.That(plain, Is.GreaterThan(0), "the search polls");
        Assert.That(report, Is.EqualTo(plain));
        var r = AnalysisTools.Fuzzpro(Protein, "K-V-L", reportFormat: "excel", sequenceName: "p");
        Assert.That(r.Hits, Has.Length.EqualTo(AnalysisTools.Fuzzpro(Protein, "K-V-L").Hits.Length));
        Assert.That(r.Report, Does.StartWith("SeqName\tStart\tEnd\tScore\tStrand\tPattern\tMismatch\np\t"));
    }

    [Test]
    public void Fuzztran_Report_RunsTheSearchOnce()
    {
        int plain = Polls(() => AnalysisTools.Fuzztran(Dna, "M-A", "6"));
        int report = Polls(() => AnalysisTools.Fuzztran(Dna, "M-A", "6", reportFormat: "gff", sequenceName: "d"));
        Assert.That(plain, Is.GreaterThan(0));
        Assert.That(report, Is.EqualTo(plain));
        var r = AnalysisTools.Fuzztran(Dna, "M-A", "6", reportFormat: "gff", sequenceName: "d");
        Assert.That(r.Hits, Has.Length.EqualTo(AnalysisTools.Fuzztran(Dna, "M-A", "6").Hits.Length));
        Assert.That(r.Report, Does.StartWith("##gff-version 3\n##sequence-region d 1 "));
    }

    [TestCase("etandem")]
    [TestCase("equicktandem")]
    public void TandemRepeats_Report_RunsTheSearchOnce(string program)
    {
        int plain = Polls(() => AnalysisTools.FindTandemRepeatsEmboss(Tandem, program, 2, 6, 8));
        int report = Polls(() => AnalysisTools.FindTandemRepeatsEmboss(Tandem, program, 2, 6, 8, reportFormat: "excel", sequenceName: "t"));
        Assert.That(plain, Is.GreaterThan(0));
        Assert.That(report, Is.EqualTo(plain));
        var r = AnalysisTools.FindTandemRepeatsEmboss(Tandem, program, 2, 6, 8, reportFormat: "excel", sequenceName: "t");
        Assert.That(r.Repeats.Length, Is.EqualTo(AnalysisTools.FindTandemRepeatsEmboss(Tandem, program, 2, 6, 8).Repeats.Length));
        Assert.That(r.Report, Does.StartWith("SeqName\tStart\tEnd\tScore\tStrand\t"));
    }

    [Test]
    public void FindPatternHits_Report_RunsTheSearchOnce_AndTheFourToolsShareTheFormatCheck()
    {
        int plain = Polls(() => AnalysisTools.FindPatternHits(Dna, "ATG"));
        int report = Polls(() => AnalysisTools.FindPatternHits(Dna, "ATG", reportFormat: "gff", sequenceName: "d"));
        Assert.That(report, Is.EqualTo(plain));
        // one helper: a format the program refuses for the sequence type (embl / genbank for proteins, pir / swiss for nucleotides, diffseq: 7 tags) and an unknown name are refused identically by the four tools
        Assert.Throws<ArgumentException>(() => AnalysisTools.Fuzzpro(Protein, "K-V", reportFormat: "genbank"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindPatternHits(Dna, "ATG", reportFormat: "pir"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.Fuzztran(Dna, "M", reportFormat: "swiss"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindTandemRepeatsEmboss(Tandem, reportFormat: "diffseq"));
        foreach (var call in new Action[]
        {
            () => AnalysisTools.FindPatternHits(Dna, "ATG", reportFormat: "nosuch"),
            () => AnalysisTools.Fuzzpro(Protein, "K-V", reportFormat: "nosuch"),
            () => AnalysisTools.Fuzztran(Dna, "M", reportFormat: "nosuch"),
            () => AnalysisTools.FindTandemRepeatsEmboss(Tandem, reportFormat: "nosuch"),
        })
        {
            Assert.That(Assert.Throws<ArgumentException>(() => call())!.ParamName, Is.EqualTo("reportFormat"));
        }

        // the report of a ranged search is written (A5-8b): the header names the range
        Assert.That(AnalysisTools.FindPatternHits(Dna, "ATG", reportFormat: "table", sequenceQualifiers: "-sbegin 2", sequenceName: "d").Report, Does.Contain("# Sequence: d     from: 2   to: "));
        Assert.That(AnalysisTools.Fuzzpro(Protein, "K-V", reportFormat: "-rformat simple -rdesshow", sequenceQualifiers: "-sbegin 2", sequenceName: "p").Report, Does.Contain("# Report_format: simple"));
    }
}
