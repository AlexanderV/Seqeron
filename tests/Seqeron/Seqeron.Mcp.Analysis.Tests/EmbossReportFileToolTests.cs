using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>emboss_report_file (B09 A5-8b, F184): the report FILE of the five EMBOSS pattern / repeat programs through the tool; the byte-for-byte oracle is EmbossReportFormats_Oracle_Tests in Seqeron.Genomics.Tests.</summary>
[TestFixture]
public class EmbossReportFileToolTests
{
    private static readonly string[] Entries = ["[-sequence] ex.fa", "-pattern ACGT", "-rformat table", "-auto", "-outfile ex.rep"];

    [Test]
    public void Fuzznuc_Table_TwoSequences_IsTheProgramsFile()
    {
        // fuzznuc ex.fa -pattern ACGT -rformat table (installed 6.6.0, run date faked to 2023-11-14 22:13:20 UTC)
        var r = AnalysisTools.EmbossReportFile("fuzznuc", ">seq1 a test\nACGTACGTAAACGT\n>seq2\nTTACGTGG\n", "ACGT", reportQualifiers: "-rformat table", commandLine: Entries, runDate: "2023-11-14T22:13:20Z", fastaFile: "ex.fa", reportFile: "ex.rep");
        Assert.That(r.OriginalFile, Is.Null);
        Assert.That(r.Report, Does.StartWith("########################################\n# Program: fuzznuc\n# Rundate: Tue 14 Nov 2023 22:13:20\n# Commandline: fuzznuc\n#    [-sequence] ex.fa\n#    -pattern ACGT\n"));
        Assert.That(r.Report, Does.Contain("  Start     End  Strand   Score Pattern      Mismatch\n      1       4       +       4 pattern:ACGT        .\n      5       8       +       4 pattern:ACGT        .\n     11      14       +       4 pattern:ACGT        .\n"));
        Assert.That(r.Report, Does.EndWith("# Total_sequences: 2\n# Total_length: 22\n# Reported_sequences: 2\n# Reported_hitcount: 4\n#---------------------------------------\n"));
    }

    [Test]
    public void Etandem_ReturnsTheOriginalFile_AndAFormatTheProgramRefusesIsAnArgumentException()
    {
        var r = AnalysisTools.EmbossReportFile("etandem", ">t\nACGACGACGACGTTTTTTTTTTGGGGGGGGGGGCACACACACACACAGTAGTAGTAGAAAA\n", minRepeat: 2, maxRepeat: 5, threshold: 8, reportQualifiers: "-rformat excel");
        Assert.That(r.Report, Does.StartWith("SeqName\tStart\tEnd\tScore\tStrand\tSize\tCount\tIdentity\tConsensus\nt\t34\t47\t12\t+\t2\t7\t100.0\tca\n"));
        Assert.That(r.OriginalFile, Is.EqualTo("    12         34         47  2   7 100.0 ca\n     9          1         12  3   4 100.0 acg\n"));
        var ex = Assert.Throws<ArgumentException>(() => AnalysisTools.EmbossReportFile("fuzzpro", ">p\nMKVLAAGIV\n", "L-A", reportQualifiers: "-rformat embl"));
        Assert.That(ex!.Message, Does.StartWith("Report format 'embl' not supported for protein sequences"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.EmbossReportFile("nosuch", ">p\nMKVLAAGIV\n", "L-A"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.EmbossReportFile("fuzznuc", ">p\nACGT\n"));
    }

    [Test]
    public void ReportOptions_AndSequenceQualifiers_AreTheProgramsOwn()
    {
        var r = AnalysisTools.EmbossReportFile("fuzzpro", ">p1 first one\nMKVLAAGIVGLLLAQ\n>p2\nMKKVLAGGIV\n", "L-A", reportQualifiers: "-rformat table -rdesshow -rmaxseq 1", runDate: "2023-11-14T22:13:20Z");
        Assert.That(r.Report, Does.Contain("# Max_hits_sequence: 1\n"));
        Assert.That(r.Report, Does.Contain("# Description: first one\n"));
        Assert.That(r.Report, Does.Contain("# HitLimit: 1/"));
        var ranged = AnalysisTools.EmbossReportFile("fuzznuc", ">s\nGGGAAATTTCCCGGGGAAAAAAAAAAAACCC\n", "GGG", sequenceQualifiers: "-sbegin 5 -send 25", reportQualifiers: "-rformat listfile");
        Assert.That(ranged.Report, Does.Contain("# Sequence: fasta::sequence.fasta:s[5:25]     from: 5   to: 25\n"));
    }
}

/// <summary>mash_screen (B09 A5-10b, F185): the tool over MashCommands.Screen; the byte-for-byte oracle is MashScreen_Oracle_Tests in Seqeron.Genomics.Tests.</summary>
[TestFixture]
public class MashScreenToolTests
{
    [Test]
    public void MashScreen_QueryContainedInReads_IsTheProgramsLine()
    {
        var r = AnalysisTools.MashScreen([">g\nACGTTGCATGCCGATAGCTAGGATCCATGCATGCAAGT\n"], [">r1\nACGTTGCATGCCGATAGCTAGGATCC\n>r2\nGATAGCTAGGATCCATGCATGCAAGT\n"], ["g.fa"], k: 8, sketchSize: 100);
        Assert.That(r.Text, Is.EqualTo("1\t30/30\t1\t6.59538e-101\tg.fa\tg \n"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.MashScreen([], [">r\nACGT\n"]));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.MashScreen([">g\nACGTTGCATGCC\n"], [">r\nACGTTGCATGCC\n"], k: 8, seed: -1));
    }
}

/// <summary>mash_dist / mash_triangle read-set mode and -C (B09 A5-10b, F186); the byte-for-byte oracle is MashReads_Oracle_Tests in Seqeron.Genomics.Tests.</summary>
[TestFixture]
public class MashReadsToolTests
{
    private const string A = ">a first desc\nACGTTGCATGCCGATAGCTAGGATCCATGCATGCAAGTAGCTAGCATCGATCGATCGACTAGCTAGCATGC\n>b second\nTTGACGTACGTAGCTAGCTAGGGATATCGATCGATCGGGGATATATAGCGCGATCGAG\n";

    [Test]
    public void ShowComments_AppendsTheComment_AndReadsModeValidatesItsOptions()
    {
        // mash dist -k 8 -s 50 -C c1.fa c2.fa (installed 2.3)
        var r = AnalysisTools.MashDist([A], [A], ["c1.fa"], ["c2.fa"], k: 8, sketchSize: 50, showComments: true);
        Assert.That(r.Text, Is.EqualTo("c1.fa:[2 seqs] a first desc [...]\tc2.fa:[2 seqs] a first desc [...]\t0\t4.29169e-151\t50/50\n"));
        var t = AnalysisTools.MashTriangle([A, A], ["c1.fa", "c2.fa"], k: 8, sketchSize: 50, showComments: true);
        Assert.That(t.Text, Does.StartWith("\t2\n[2 seqs] a first desc [...]\n"));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.MashDist([A], [A], k: 8, reads: true, minCopies: 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.MashDist([A], [A], k: 8, reads: true, targetCoverage: -1));
        Assert.Throws<ArgumentException>(() => AnalysisTools.MashDist([A], [A], k: 8, reads: true, individual: true));
    }

    [Test]
    public void FuzznucPatternFile_IsOnePatternPerRowOfTheHeader()
    {
        // fuzznuc -pattern @file on s1 ACGTACGTTTGACGTAAACCCGGGTTTACGTACGATCAAACCCGGGTTT with a FASTA style file with <mismatch=1> (A5-8c): measured with the installed 6.6.0
        var r = AnalysisTools.EmbossReportFile("fuzznuc", ">s1 d1\nACGTACGTTTGACGTAAACCCGGGTTTACGTACGATCAAACCCGGGTTT\n", patternFile: ">pa<mismatch=1>\nCCCGGG\n>pb desc x\nACGT\n>pc\nTTT\n", complement: true, reportQualifiers: "-rformat seqtable");
        Assert.That(r.Report, Does.Contain("# Pattern_name Mismatch Pattern\n# pa                  1 CCCGGG\n# pb desc x           0 ACGT\n# pc                  0 TTT\n"));
        Assert.That(r.Report, Does.Contain("     19      24       + pa:CCCGGG             . CCCGGG\n").And.Contain("     16      18       - pc:TTT                . TTT\n").And.Contain("# HitCount: 17"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.EmbossReportFile("fuzznuc", ">a\nACGT\n", "ACGT", patternFile: "ACGT\n"));
    }
}
