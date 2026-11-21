using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>B09 A5-7: the <c>sequenceQualifiers</c> parameter of fuzznuc / fuzzpro / fuzztran / etandem / equicktandem / seqmatchall / distmat. Values are what the installed EMBOSS 6.6.0 prints.</summary>
[TestFixture]
public class EmbossSequenceQualifiersToolsTests
{
    private const string S = "ACGTACGTTAGCCGATTAGGCATCGATCGGATCGAAACCCTTTGGG";

    [Test]
    public void FindPatternHits_Reverse_ReportsTheMinusStrandInOriginalCoordinates()
    {
        // fuzznuc s.fa -pattern CGATC -sreverse -rformat excel: s 30 34 - and s 25 29 - (original numbering); without -sreverse: 24 28 +
        var r = AnalysisTools.FindPatternHits(S, "CGATC", sequenceQualifiers: "-sreverse");
        Assert.That(r.Hits.Select(h => (h.Position + 1, h.Position + h.Length, h.Strand)).OrderBy(h => h.Item1).ToArray(),
            Is.EqualTo(new[] { (25, 29, "reverse"), (30, 34, "reverse") }));
        var plain = AnalysisTools.FindPatternHits(S, "CGATC");
        Assert.That(plain.Hits.Select(h => (h.Position + 1, h.Strand)).ToArray(), Is.EqualTo(new[] { (24, "forward") }));
    }

    [Test]
    public void FindPatternHits_Range_AddsTheBeginOffset()
    {
        var r = AnalysisTools.FindPatternHits(S, "CGATC", sequenceQualifiers: "-sbegin 5 -send 40");
        Assert.That(r.Hits.Single().Position, Is.EqualTo(23));
    }

    [Test]
    public void FindPatternHits_NumberedEqualsAndAbbreviatedForms_AreResolvedLikeTheProgram()
    {
        // A6-2: fuzznuc s.fa -pattern CGATC -sbegin1 5 -send1 40 == -sbegin 5 -send 40; also -sbegin1=5, the unique abbreviations -sbeg1 / -sen1 / --SEND=40
        var expected = AnalysisTools.FindPatternHits(S, "CGATC", sequenceQualifiers: "-sbegin 5 -send 40").Hits.Select(h => (h.Position, h.Length, h.Strand)).ToArray();
        Assert.That(expected, Is.EqualTo(new[] { (23, 5, "forward") }));
        foreach (string text in new[] { "-sbegin1 5 -send1 40", "-sbegin1=5 -send1=40", "-sbeg1 5 -sen1 40", "-SBEGIN=5 --send 40", "-sbegin1 5 -send 40" })
            Assert.That(AnalysisTools.FindPatternHits(S, "CGATC", sequenceQualifiers: text).Hits.Select(h => (h.Position, h.Length, h.Strand)).ToArray(), Is.EqualTo(expected), text);
        // a range that starts after the hit finds nothing; -se is ambiguous between -send and -sequence in fuzznuc, -sbegin2 names a sequence it does not have
        Assert.That(AnalysisTools.FindPatternHits(S, "CGATC", sequenceQualifiers: "-sbeg1 30").Hits, Is.Empty);
        Assert.That(Assert.Throws<ArgumentException>(() => AnalysisTools.FindPatternHits(S, "CGATC", sequenceQualifiers: "-se 40"))!.Message, Does.Contain("ambiguous qualifier 'se' (send,sequence)"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindPatternHits(S, "CGATC", sequenceQualifiers: "-sbegin2 5"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindPatternHits(S, "CGATC", sequenceQualifiers: "-pattern ACG"));
    }

    [Test]
    public void Fuzzpro_Range_CountsFromTheOriginalStart()
    {
        var plain = AnalysisTools.Fuzzpro("MKVLAAGIVGLLLAQ", "K-V");
        var ranged = AnalysisTools.Fuzzpro("MKVLAAGIVGLLLAQ", "K-V", sequenceQualifiers: "-sbegin 2");
        Assert.That(ranged.Hits.Single().Position, Is.EqualTo(plain.Hits.Single().Position));
        Assert.Throws<ArgumentException>(() => AnalysisTools.Fuzzpro("MKVLAAGIVGLLLAQ", "K-V", sequenceQualifiers: "-snucleotide"));
    }

    [Test]
    public void Etandem_IgnoresTheRange_Equicktandem_UsesIt()
    {
        const string seq = "GGTACCACACACACACACACGTTGCATGCAATGCAATGCAATGCAATGCATTAGC";
        // etandem et.fa -minrepeat 2 -maxrepeat 6 -threshold 8 -sbegin 10 -send 40 prints the same repeats as without
        var whole = AnalysisTools.FindTandemRepeatsEmboss(seq, "etandem", 2, 6, 8);
        var cut = AnalysisTools.FindTandemRepeatsEmboss(seq, "etandem", 2, 6, 8, sequenceQualifiers: "-sbegin 10 -send 40");
        Assert.That(cut.Repeats.Select(r => (r.Position, r.Length)).ToArray(), Is.EqualTo(whole.Repeats.Select(r => (r.Position, r.Length)).ToArray()));
        // equicktandem -maxrepeat 6 -threshold 8 -sbegin 10 -send 40: the first line of the program is  9  10  20  2  5
        var quick = AnalysisTools.FindTandemRepeatsEmboss(seq, "equicktandem", 2, 6, 8, sequenceQualifiers: "-sbegin 10 -send 40");
        // the whole output of the program: -origfile of equicktandem et.fa -maxrepeat 6 -threshold 8 -sbegin 10 -send 40 is "9 10 20 2 5" and "10 26 40 5 3" (score, start, end, size, count); without the range "13 6 20 2 7" and "20 26 50 5 5"
        Assert.That(quick.Repeats.Select(r => (r.Score, r.Position + 1, r.Position + r.Length, r.RepeatSize, r.RepeatCount)).ToArray(), Is.EqualTo(new[] { (9, 10, 20, 2, 5), (10, 26, 40, 5, 3) }));
        var quickWhole = AnalysisTools.FindTandemRepeatsEmboss(seq, "equicktandem", 2, 6, 8);
        Assert.That(quickWhole.Repeats.Select(r => (r.Score, r.Position + 1, r.Position + r.Length, r.RepeatSize, r.RepeatCount)).ToArray(), Is.EqualTo(new[] { (13, 6, 20, 2, 7), (20, 26, 50, 5, 5) }));
    }

    [Test]
    public void Seqmatchall_Reverse_NumbersFromTheEnd()
    {
        var r = AnalysisTools.SeqMatchAll(["TCCACTAGCGAGAGGTACCCTGATAAGGGGGAGGCAGCGTAGGTTA", "TCTCTTAGCGGGAACTAGCCTGAGCTGATCGGGGGTCTAGCCGCCTAGGGTTAGG"], 3, sequenceQualifiers: "-sreverse");
        // the complete match list of `seqmatchall -sequence s.fa -wordsize 3 -sreverse -aformat match` (installed 6.6.0): 42 matches, every one on the reverse strand of both sequences, numbered from the end (start > end)
        var expected = new (int Start1, int End1, int Start2, int End2, int Length)[] { (3, 1, 51, 49, 3), (5, 1, 7, 3, 5), (5, 3, 22, 20, 3), (7, 4, 4, 1, 4), (7, 4, 10, 7, 4), (7, 5, 18, 16, 3), (7, 5, 40, 38, 3), (7, 5, 50, 48, 3), (11, 8, 49, 46, 4), (11, 9, 17, 15, 3), (11, 9, 33, 31, 3), (11, 9, 39, 37, 3), (15, 13, 3, 1, 3), (15, 13, 9, 7, 3), (16, 14, 34, 32, 3), (18, 15, 46, 43, 4), (18, 16, 8, 6, 3), (18, 16, 25, 23, 3), (19, 16, 25, 22, 4), (19, 17, 8, 6, 3), (19, 17, 46, 44, 3), (20, 16, 25, 21, 5), (20, 17, 24, 21, 4), (20, 18, 23, 21, 3), (20, 18, 46, 44, 3), (21, 18, 9, 6, 4), (21, 19, 3, 1, 3), (27, 23, 31, 27, 5), (28, 24, 37, 33, 5), (28, 26, 12, 10, 3), (33, 31, 7, 5, 3), (33, 31, 22, 20, 3), (34, 32, 3, 1, 3), (34, 32, 9, 7, 3), (35, 33, 34, 32, 3), (37, 35, 34, 32, 3), (40, 38, 33, 31, 3), (41, 37, 50, 46, 5), (41, 39, 4, 2, 3), (42, 38, 19, 15, 5), (42, 39, 11, 8, 4), (43, 38, 42, 37, 6) };
        var actual = r.ReportedMatches!.OrderBy(m => m.Start1).ThenBy(m => m.End1).ThenBy(m => m.Start2).ThenBy(m => m.End2).ThenBy(m => m.Length).Select(m => (m.Start1, m.End1, m.Start2, m.End2, m.Length)).ToArray();
        Assert.That(actual, Is.EqualTo(expected));
        Assert.That(r.ReportedMatches!.All(m => m.Reversed), Is.True);
    }

    [Test]
    public void ReportFormat_ReturnsTheProgramsReportText()
    {
        // B09 A5-8: fuzznuc a -pattern ACGTA -rformat gff / excel
        var gff = AnalysisTools.FindPatternHits("ACGTAGGG", "ACGTA", reportFormat: "gff", sequenceName: "a");
        Assert.That(gff.Report, Does.StartWith("##gff-version 3\n##sequence-region a 1 8\n#!Date "));
        Assert.That(gff.Report, Does.EndWith("a\tfuzznuc\tnucleotide_motif\t1\t5\t5\t+\t.\tID=a.1;note=*pat pattern:ACGTA\n"));
        var excel = AnalysisTools.FindPatternHits("ACGTAGGG", "ACGTA", reportFormat: "excel", sequenceName: "a");
        Assert.That(excel.Report, Is.EqualTo("SeqName\tStart\tEnd\tScore\tStrand\tPattern\tMismatch\na\t1\t5\t5\t+\tpattern:ACGTA\t.\n"));
        Assert.That(AnalysisTools.FindPatternHits("ACGTAGGG", "ACGTA").Report, Is.Null);
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindPatternHits("ACGTAGGG", "ACGTA", reportFormat: "seqtable"));
        var pro = AnalysisTools.Fuzzpro("ACDKLEFGHIKLMNPQRSTVWY", "K-L", reportFormat: "excel", sequenceName: "p");
        Assert.That(pro.Report, Is.EqualTo("SeqName\tStart\tEnd\tScore\tStrand\tPattern\tMismatch\np\t4\t5\t2\t+\tpattern:K-L\t.\np\t11\t12\t2\t+\tpattern:K-L\t.\n"));
    }

    [Test]
    public void Megamerger_MergesTheTrimmedSequences()
    {
        // megamerger a.fa b.fa -wordsize 6 -sbegin1 37 (installed 6.6.0, oracle row of emboss_seqopts_oracle.tsv)
        const string a = "CAGACCGATGCGCGCCGTCTTCGCCGCTTAGCTTGACGACAGTACCCGGCGGGAATGCTGGCTGGCT";
        const string b = "AGACCGATGCGCGCCGTCTTCGCCGCTTAGCTTGACGACAGTACCCGGCGGGAATGCTGGCTGGCTACACCAGTCGGTCTACAACTTGCTGCGAA";
        var r = AnalysisTools.Megamerger(a, b, wordSize: 6, name1: "a", name2: "b", sequenceQualifiers: "-sbegin1 37");
        Assert.That(r.MergedSequence, Is.EqualTo("agaccgatgcgcgccgtcttcgccgcttagcttgacgacagtacccggcgggaatgctggctggctacaccagtcggtctacaacttgctgcgaa"));
        Assert.That(r.Report, Does.Contain("Using b 1-35 as the initial sequence"));
        Assert.That(r.Report, Does.Contain("Matching region a 1-31 : b 36-66"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.Megamerger(a, b, wordSize: 6, sequenceQualifiers: "-sprotein1"));
    }
}
