using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;

namespace Seqeron.Mcp.Alignment.Tests;

/// <summary>stretcher_align / write_alignment: EMBOSS 6.6 stretcher and the -aformat output formats (audit A2-27 / A2-28, B09 F120 / F121). The tools delegate to SequenceAligner.StretcherAlign and EmbossAlignmentWriter; the oracle comparisons live in the library tests.</summary>
[TestFixture]
public class EmbossFormatTests
{
    private const string A = "ACGTACGTAGCTAGCTAGCATCGATCGATCGAT";
    private const string B = "ACGTACGTAGCTAGATAGCATCGAATCGATCGAT";

    [Test]
    public void StretcherAlign_EqualsTheLibrary_AndTheDocumentedEmbossOutput()
    {
        var tool = AlignmentTools.StretcherAlign(A, B);
        var lib = Seqeron.Genomics.Alignment.SequenceAligner.StretcherAlign(A, B);
        Assert.Multiple(() =>
        {
            Assert.That(tool.AlignedSequence1, Is.EqualTo(lib.AlignedSequence1));
            Assert.That(tool.AlignedSequence2, Is.EqualTo(lib.AlignedSequence2));
            Assert.That(tool.Score, Is.EqualTo(140), "stretcher a.fa b.fa: Score: 140");
            Assert.That(tool.AlignedSequence1, Is.EqualTo("ACGTACGTAGCTAGCTAGCATCGA-TCGATCGAT"));
        });
    }

    [Test]
    public void StretcherAlign_Validation()
    {
        Assert.Throws<ArgumentException>(() => AlignmentTools.StretcherAlign(null!, A));
        Assert.Throws<ArgumentException>(() => AlignmentTools.StretcherAlign(A, B, matrix: "NOSUCH"));
        Assert.Throws<ArgumentOutOfRangeException>(() => AlignmentTools.StretcherAlign(A, B, gapOpen: -1));
    }

    [Test]
    public void WriteAlignment_PairFormat_IsTheNeedleLayout()
    {
        var r = AlignmentTools.WriteAlignment(["ACGTACGTAGCTAGCTAGCATCG-ATCGATCGAT", B], "pair", ["a", "b"], score: 146, matrix: "EDNAFULL", gapOpen: 10, gapExtend: 0.5,
            runDate: "2026-11-19T21:43:53", commandLine: ["-asequence a.fa", "-bsequence b.fa", "-aformat3 pair", "-auto", "-stdout"]);
        Assert.Multiple(() =>
        {
            Assert.That(r.Format, Is.EqualTo("pair"));
            Assert.That(r.Dialect, Is.EqualTo("emboss"));
            Assert.That(r.Text, Does.StartWith("########################################\n# Program: needle\n# Rundate: Thu 19 Nov 2026 21:43:53\n# Commandline: needle\n#    -asequence a.fa\n"));
            Assert.That(r.Text, Does.Contain("# Identity:      32/34 (94.1%)\n"));
            Assert.That(r.Text, Does.Contain("# Score: 146.0\n"));
            Assert.That(r.Text, Does.Contain("a                  1 ACGTACGTAGCTAGCTAGCATCG-ATCGATCGAT     33\n"));
        });
    }

    [Test]
    public void WriteAlignment_SequenceSetFormats_TakeAnyNumberOfRows()
    {
        var rows = new[] { "ACGT-ACGTA", "AC-TTACG-A", "ACGTTAC--A" };
        var fasta = AlignmentTools.WriteAlignment(rows, "fasta", ["s1", "s2", "s3"], descriptions: ["desc one", "desc one", "desc one"]);
        Assert.That(fasta.Text, Is.EqualTo(">s1 desc one\nACGT-ACGTA\n>s2 desc one\nAC-TTACG-A\n>s3 desc one\nACGTTAC--A\n"));
        var clustal = AlignmentTools.WriteAlignment(rows, "clustal", ["s1", "s2", "s3"]);
        Assert.That(clustal.Text, Does.StartWith("CLUSTAL W (1.83) multiple sequence alignment\n\n\ns1              ACGT-ACGTA\n"));
        var phylip = AlignmentTools.WriteAlignment(rows, "phylip", ["s1", "s2", "s3"]);
        Assert.That(phylip.Text, Is.EqualTo(" 3 10\ns1        ACGT-ACGTA\ns2        AC-TTACG-A\ns3        ACGTTAC--A\n"));
    }

    [Test]
    public void WriteAlignment_Validation()
    {
        Assert.Throws<ArgumentException>(() => AlignmentTools.WriteAlignment(["ACGT", "ACGT"], "bam"));
        Assert.Throws<ArgumentException>(() => AlignmentTools.WriteAlignment([], "fasta"));
        Assert.Throws<ArgumentException>(() => AlignmentTools.WriteAlignment(["ACGT", "ACGT", "ACGT"], "pair"), "pairwise formats take exactly two rows");
        Assert.Throws<ArgumentNullException>(() => AlignmentTools.WriteAlignment(null!, "fasta"));
        Assert.Throws<ArgumentException>(() => AlignmentTools.WriteAlignment(["ACGT"], "score"), "one row: score needs two");
    }
}
