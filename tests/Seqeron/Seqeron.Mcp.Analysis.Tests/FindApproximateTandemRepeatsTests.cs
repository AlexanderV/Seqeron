using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>find_approximate_tandem_repeats</c> MCP tool. Expected rows = compiled TRF 4.10.0
/// (`trf craft.fa 2 7 7 80 10 50 500 -h -d -ngs`, 1-based in TRF), locked in RepeatFinder_TrfParameters_Tests
/// (Evidence REP-APPROX-001 §WP6). NOT the wrapper's output.
/// </summary>
[TestFixture]
public class FindApproximateTandemRepeatsTests
{
    private const string U1 =
        "CTATCCTAACCCGACCCTAGGAGCGGTTGGCGTGTATGCCGTGAATTTTCTCATTTCCGCTCATTGGTCATTGGTCANTGGTCATTGGTCATTGGTCATTNGTCATTGGTCATTGGTCATTGGTAGACATAATCGTTCTGCCTATATCTGGACAACATCCCGGCGACTTAGGCGACCCACAGA";
    private const string U3 =
        "CACACACACACACACACGCACACACACACCGACACACACACACACACACACACACACACAAGATCGTCAACGACCATCGTTCTTCTGCTTTAAGTGTGAGTTCTCTCTTACAAGCTAATACACCCCAGCGTTCTCCGTAC";

    [Test]
    public void FindApproximateTandemRepeats_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.FindApproximateTandemRepeats(U3));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindApproximateTandemRepeats(""));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindApproximateTandemRepeats(null!));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.FindApproximateTandemRepeats(U3, matchProbability: 70));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.FindApproximateTandemRepeats(U3, maxPeriod: 2001));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.FindApproximateTandemRepeats(U3, minPeriod: 0));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindApproximateTandemRepeats(U3, mismatchPenalty: 3, examineUpToMaxPeriodOnly: true));
    }

    [Test]
    public void FindApproximateTandemRepeats_Binding_RecommendedParameters_ReproduceTrfRow()
    {
        // TRF U1 row: 61 124 7 9.1 7 92 0 110 14 14 26 42 1.83 TCATTGG
        var r = AnalysisTools.FindApproximateTandemRepeats(U1).Items.Single();
        Assert.Multiple(() =>
        {
            Assert.That((r.Start + 1, r.Start + r.SpanLength, r.Period, r.ConsensusSize), Is.EqualTo((61, 124, 7, 7)));
            Assert.That(r.CopyNumber, Is.EqualTo(9.1).Within(0.05 + 1e-9));
            Assert.That(((int)(r.PercentMatches + 1e-9), (int)(r.PercentIndels + 1e-9), r.AlignmentScore), Is.EqualTo((92, 0, 110)));
            Assert.That(((int)(r.PercentA + 1e-9), (int)(r.PercentC + 1e-9), (int)(r.PercentG + 1e-9), (int)(r.PercentT + 1e-9)),
                Is.EqualTo((14, 14, 26, 42)));
            Assert.That(r.EntropyTrf, Is.EqualTo(1.829258111162015).Within(1e-12));
            Assert.That(r.Consensus, Is.EqualTo("TCATTGG"));
            Assert.That(r.AlignedSequence, Is.Not.Null);
        });
    }

    [Test]
    public void FindApproximateTandemRepeats_CustomWeights_ReproduceTrfRow()
    {
        // `trf 2 3 5 80 10 40 200` → U3 1 60 2 30.0 2 89 0 105 CA
        var r = AnalysisTools.FindApproximateTandemRepeats(U3, maxPeriod: 200, minScore: 40, mismatchPenalty: 3, indelPenalty: 5).Items.Single();
        Assert.That((r.Start + 1, r.Start + r.SpanLength, r.Period, (int)(r.PercentMatches + 1e-9), r.AlignmentScore, r.Consensus),
            Is.EqualTo((1, 60, 2, 89, 105, "CA")));
    }

    [Test]
    public void FindApproximateTandemRepeats_NoRedundancyElimination_ReportsTrfMultiples()
    {
        // `trf ... -r`: U1 reports periods 7, 14, 21 over 61..124 (all score 110).
        var items = AnalysisTools.FindApproximateTandemRepeats(U1, eliminateRedundancy: false).Items;
        Assert.That(items.Select(r => (r.Start + 1, r.Start + r.SpanLength, r.Period, r.AlignmentScore)),
            Is.EqualTo(new[] { (61, 124, 7, 110), (61, 124, 14, 110), (61, 124, 21, 110) }));
    }

    // RepeatFinder_TrfDetection_Tests D2: TRF reports nothing on this sequence; its apparent-size test rejects a clustered
    // period-28 candidate (148..209, score 92) that a table of zeros (no apparent-size rejection) admits.
    private const string Spread28 =
        "CAGTCACGGGCTCTGGATCCAGCAGCAGTGCAGCATGTTGGTACCCTATCCCCATACGACACTGTTTGGCGCTGTTGGTTTATGCACGAGTCGTTACTAT"
        + "ATAAAGACCTCGAAGTGCCAGAATTCATCTTTGACCTCAGCGCGTTCGTACTCCGATCGGAACCGCCCGTTCACTGTACTCCGATCGGAACCGCCCCGAT"
        + "ATGTACTCCATTAATCGTCCCTTTGAATTCGGAGATACGCGTGACGGACGTATCGCGTCTCCATTCTTAGCCGACTCCACGACCTCCTTAATGGTTAATC"
        + "AACATAAGAATATTCCCAGGAG";

    private static string Zeros() => string.Join(",", Enumerable.Repeat(0, 2001));

    // TRF 4.10.0 `trf U1.fa 2 7 7 80 10 50 500 -d -h` (.dat) / `-ngs -h` (stdout) / without -h (.1.html).
    [Test]
    public void FindApproximateTandemRepeats_Formats_ReproduceTrfOutput()
    {
        const string row = "61 124 7 9.1 7 92 0 110 14 14 26 42 1.83 TCATTGG TCATTGGTCATTGGTCANTGGTCATTGGTCATTGGTCATTNGTCATTGGTCATTGGTCATTGGT";
        var dat = AnalysisTools.FindApproximateTandemRepeats(U1, format: "dat", sequenceName: "U1");
        var ngs = AnalysisTools.FindApproximateTandemRepeats(U1, format: "ngs", sequenceName: "U1");
        var html = AnalysisTools.FindApproximateTandemRepeats(U1, format: "html", sequenceName: "U1.fa");
        var json = AnalysisTools.FindApproximateTandemRepeats(U1);
        Assert.Multiple(() =>
        {
            Assert.That(dat.Formatted, Is.EqualTo("Tandem Repeats Finder Program written by:\n\nGary Benson\nProgram in Bioinformatics\nBoston University\nVersion 4.10.0\n"
                + "\n\nSequence: U1\n\n\n\nParameters: 2 7 7 80 10 50 500\n\n\n" + row + "\n"));
            Assert.That(ngs.Formatted, Is.EqualTo("@U1\n" + row
                + " CCGACCCTAGGAGCGGTTGGCGTGTATGCCGTGAATTTTCTCATTTCCGC AGACATAATCGTTCTGCCTATATCTGGACAACATCCCGGCGACTTAGGCG\n"));
            Assert.That(html.HtmlPages!.Select(p => p.FileName), Is.EqualTo(new[] { "U1.fa.2.7.7.80.10.50.500.1.html" }));
            Assert.That(html.HtmlPages![0].Html, Does.Contain("<A HREF=\"U1.fa.2.7.7.80.10.50.500.1.txt.html#61--124,7,9.1,7,1\">61--124</A>"));
            Assert.That(html.HtmlPages![0].Html, Does.Contain("<TD><CENTER>9.1</CENTER></TD><TD><CENTER>7</CENTER></TD><TD><CENTER>92</CENTER></TD>"));
            Assert.That(json.Formatted, Is.Null);
            Assert.That(json.HtmlPages, Is.Null);
            Assert.That(json.Items.Single().OutputIndex, Is.EqualTo(1));
            Assert.That((json.Items.Single().CopyMatches + json.Items.Single().CopyMismatches + json.Items.Single().CopyIndels) > 0);
            Assert.Throws<ArgumentException>(() => AnalysisTools.FindApproximateTandemRepeats(U1, format: "xml"));
        });
    }

    // TRF 4.10.0 `trf U1.fa 2 7 7 80 10 50 500` → U1.fa.2.7.7.80.10.50.500.1.txt.html (the page the table links to); the MCP
    // tool uses sequenceName as both the description and the file prefix, so only the "Sequence:" line differs.
    [Test]
    public void FindApproximateTandemRepeats_Html_ReturnsTrfAlignmentPages()
    {
        const string trf = "<HTML><HEAD><TITLE>U1.fa.2.7.7.80.10.50.500.txt.html</TITLE></HEAD><BODY bgcolor=\"#FBF8BC\"><PRE>\n"
            + "Tandem Repeats Finder Program written by:\n"
            + "\n"
            + "                 Gary Benson\n"
            + "      Program in Bioinformatics\n"
            + "          Boston University\n"
            + "\n"
            + "Version 4.10.0\n"
            + "\n"
            + "Sequence: U1\n"
            + "\n"
            + "Parameters: 2 7 7 80 10 50 500\n"
            + "\n"
            + "Pmatch=0.80,Pindel=0.10\n"
            + "tuple sizes 0,4,5,7\n"
            + "tuple distances 0, 29, 159, 200\n"
            + "\n"
            + "Length: 183\n"
            + "ACGTcount: A:0.20, C:0.25, G:0.23, T:0.31\n"
            + "\n"
            + "Warning! 2 characters in sequence are not A, C, G, or T\n"
            + "\n"
            + "\n"
            + "Found at i:71 original size:7 final size:7\n"
            + "\n"
            + "<A NAME=\"61--124,7,9.1,7,1\"></A><A HREF=\"http://tandem.bu.edu/trf/trf.definitions.html#alignment\" target =\"explanation\">Alignment explanation</A><BR><BR>\n"
            + "    Indices: 61--124  Score: 110\n"
            + "    Period size: 7  Copynumber: 9.1  Consensus size: 7\n"
            + "\n"
            + "         51 TCATTTCCGC\n"
            + "\n"
            + "                   \n"
            + "         61 TCATTGG\n"
            + "          1 TCATTGG\n"
            + "\n"
            + "                   \n"
            + "         68 TCATTGG\n"
            + "          1 TCATTGG\n"
            + "\n"
            + "               *   \n"
            + "         75 TCANTGG\n"
            + "          1 TCATTGG\n"
            + "\n"
            + "                   \n"
            + "         82 TCATTGG\n"
            + "          1 TCATTGG\n"
            + "\n"
            + "                   \n"
            + "         89 TCATTGG\n"
            + "          1 TCATTGG\n"
            + "\n"
            + "                 * \n"
            + "         96 TCATTNG\n"
            + "          1 TCATTGG\n"
            + "\n"
            + "                   \n"
            + "        103 TCATTGG\n"
            + "          1 TCATTGG\n"
            + "\n"
            + "                   \n"
            + "        110 TCATTGG\n"
            + "          1 TCATTGG\n"
            + "\n"
            + "                   \n"
            + "        117 TCATTGG\n"
            + "          1 TCATTGG\n"
            + "\n"
            + "             \n"
            + "        124 T\n"
            + "          1 T\n"
            + "\n"
            + "        125 AGACATAATC\n"
            + "\n"
            + "\n"
            + "Statistics\n"
            + "Matches: 53,  Mismatches: 4, Indels: 0\n"
            + "        0.93            0.07        0.00\n"
            + "\n"
            + "Matches are distributed among these distances:\n"
            + "   7   53  1.00\n"
            + "\n"
            + "ACGTcount: A:0.14, C:0.14, G:0.27, T:0.42\n"
            + "\n"
            + "\n"
            + "Consensus pattern (7 bp):   \n"
            + "TCATTGG\n"
            + "\n"
            + "Done.\n"
            + "</PRE></BODY></HTML>\n";
        var html = AnalysisTools.FindApproximateTandemRepeats(U1, format: "html", sequenceName: "U1.fa");
        var none = AnalysisTools.FindApproximateTandemRepeats(U1, maxPeriod: 3, format: "html", sequenceName: "U1.fa");
        Assert.Multiple(() =>
        {
            Assert.That(html.AlignmentPages!.Select(p => p.FileName), Is.EqualTo(new[] { "U1.fa.2.7.7.80.10.50.500.1.txt.html" }));
            Assert.That(html.AlignmentPages![0].Html, Is.EqualTo(trf.Replace("Sequence: U1\n", "Sequence: U1.fa\n")));
            var item = html.Items.Single();
            Assert.That((item.DetectionPosition + 1, item.DetectionDistance, item.OutputCount), Is.EqualTo((71, 7, 3)));
            Assert.That(none.Items, Is.Empty);
            Assert.That(none.AlignmentPages!.Single().Html, Does.EndWith("T:0.31\n\nWarning! 2 characters in sequence are not A, C, G, or T\n\n\nDone.\n</PRE></BODY></HTML>\n"),
                "three alignments were reported and dropped (period > 3): TRF keeps the blank line before them");
            Assert.That(AnalysisTools.FindApproximateTandemRepeats(U1).AlignmentPages, Is.Null);
        });
    }

    [Test]
    public void FindApproximateTandemRepeats_ApparentSizeTable_ChangesDetection()
    {
        // y = 0 everywhere, given directly or as TRF waiting times w = max(d,20) - 1.
        string waits = string.Join(",", Enumerable.Range(0, 2001).Select(d => d == 0 ? 0 : Math.Max(d, 20) - 1));
        Assert.Multiple(() =>
        {
            Assert.That(AnalysisTools.FindApproximateTandemRepeats(Spread28).Items, Is.Empty);
            foreach (var items in new[]
            {
                AnalysisTools.FindApproximateTandemRepeats(Spread28, apparentSizeTable: Zeros()).Items,
                AnalysisTools.FindApproximateTandemRepeats(Spread28, apparentSizeTable: waits, apparentSizeTableKind: "trfWaitingTimes").Items,
            })
                Assert.That(items.Select(r => (r.Start + 1, r.Start + r.SpanLength, r.Period, r.AlignmentScore)), Is.EqualTo(new[] { (148, 209, 28, 92) }));
            Assert.Throws<ArgumentException>(() => AnalysisTools.FindApproximateTandemRepeats(Spread28, apparentSizeTable: "1,2,3"));
            Assert.Throws<ArgumentException>(() => AnalysisTools.FindApproximateTandemRepeats(Spread28, apparentSizeTable: "0,x"));
            Assert.Throws<ArgumentException>(() => AnalysisTools.FindApproximateTandemRepeats(Spread28, apparentSizeTable: Zeros(), apparentSizeTableKind: "other"));
            Assert.Throws<ArgumentException>(() => AnalysisTools.FindApproximateTandemRepeats(Spread28, apparentSizeTable: Zeros(), examineUpToMaxPeriodOnly: true));
        });
    }
}
