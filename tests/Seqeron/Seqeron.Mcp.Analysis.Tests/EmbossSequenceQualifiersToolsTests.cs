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
        Assert.That(quick.Repeats.Any(r => r.Score == 9 && r.Position + 1 == 10 && r.Position + r.Length == 20), Is.True);
    }

    [Test]
    public void Seqmatchall_Reverse_NumbersFromTheEnd()
    {
        var r = AnalysisTools.SeqMatchAll(["TCCACTAGCGAGAGGTACCCTGATAAGGGGGAGGCAGCGTAGGTTA", "TCTCTTAGCGGGAACTAGCCTGAGCTGATCGGGGGTCTAGCCGCCTAGGGTTAGG"], 3, sequenceQualifiers: "-sreverse");
        Assert.That(r.ReportedMatches!.Any(m => m.Start1 == 3 && m.End1 == 1 && m.Start2 == 51 && m.End2 == 49 && m.Length == 3 && m.Reversed), Is.True);
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
