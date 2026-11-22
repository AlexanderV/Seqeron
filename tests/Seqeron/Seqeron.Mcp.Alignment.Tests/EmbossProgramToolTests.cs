using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;

namespace Seqeron.Mcp.Alignment.Tests;

/// <summary>supermatcher_align / needleall_align / ksw2_extend and the `format` parameter of the pairwise and multiple-alignment tools (audit A2-27b / A2-28c / A2-29b, B09 F124 - F126). The tools delegate to
/// SequenceAligner.Supermatcher / NeedleAll, Ksw2Aligner and AlignmentFormats / EmbossAlignmentWriter; the oracle comparisons live in the library tests.</summary>
[TestFixture]
public class EmbossProgramToolTests
{
    private const string Q = "ACGTACGTACGGTTAACCGGTTAAACCC";

    [Test]
    public void SupermatcherAlign_IdenticalSequences_OneWordSeededAlignment()
    {
        var r = AlignmentTools.SupermatcherAlign([Q], [Q]);
        Assert.Multiple(() =>
        {
            Assert.That(r.Count, Is.EqualTo(1));
            Assert.That(r.Alignments[0].AlignedSequence1, Is.EqualTo(Q));
            Assert.That(r.Alignments[0].AlignedSequence2, Is.EqualTo(Q));
            Assert.That(r.Alignments[0].Score, Is.EqualTo(140.0), "28 matches x EDNAFULL 5");
            Assert.That(r.Alignments[0].Name1, Is.EqualTo("query1"));
            Assert.That(r.ErrorText, Is.Empty);
        });
    }

    [Test]
    public void SupermatcherAlign_NoSharedWord_ReportsTheErrorFileLine()
    {
        var r = AlignmentTools.SupermatcherAlign(["AAAAAAAAAAAAAAAAAAAA"], ["CCCCCCCCCCCCCCCCCCCC"], ["q"], ["t"]);
        Assert.That(r.Count, Is.EqualTo(0));
        Assert.That(r.ErrorText, Is.EqualTo("No wordmatch start points for q vs t. No alignment\n"));
    }

    [Test]
    public void SupermatcherAlign_Validation()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => AlignmentTools.SupermatcherAlign([Q], [Q], wordLength: 2));
        Assert.Throws<ArgumentOutOfRangeException>(() => AlignmentTools.SupermatcherAlign([Q], [Q], width: 0));
        Assert.Throws<ArgumentException>(() => AlignmentTools.SupermatcherAlign([], [Q]));
        Assert.Throws<ArgumentException>(() => AlignmentTools.SupermatcherAlign([Q], [Q], matrix: "NOSUCH"));
    }

    [Test]
    public void NeedleallAlign_PairsInTheProgramOrder_AndTheErrorLine()
    {
        var r = AlignmentTools.NeedleallAlign(["ACGT", "AAAA"], ["ACGT", "CCCC"]);
        Assert.Multiple(() =>
        {
            Assert.That(r.Alignments.Select(a => a.Name1 + "/" + a.Name2), Is.EqualTo(new[] { "a1/b1", "a2/b1", "a1/b2" }));
            Assert.That(r.Alignments[0].Score, Is.EqualTo(20.0));
            Assert.That(r.ErrorText, Does.Contain("is less than minimum score(0.0) for sequences a2 vs b2"));
        });
    }

    [Test]
    public void KSw2Extend_IdenticalSequences_AndTheZDropStop()
    {
        var r = AlignmentTools.Ksw2Extend("ACGTACGT", "ACGTACGT");
        Assert.Multiple(() =>
        {
            Assert.That(r.Max, Is.EqualTo(16));
            Assert.That(r.Score, Is.EqualTo(16));
            Assert.That(r.Cigar, Is.EqualTo("8M"));
            Assert.That(r.ZDropped, Is.False);
        });
        string common = "ACGTTGCAACGTAGCTAGCTAGGATCCATCGA";
        var z = AlignmentTools.Ksw2Extend(common + new string('A', 60), common + new string('C', 60), extensionOnly: true, zdrop: 10);
        Assert.Multiple(() =>
        {
            Assert.That(z.ZDropped, Is.True);
            Assert.That(z.Max, Is.EqualTo(common.Length * 2));
            Assert.That(z.Cigar, Is.EqualTo(common.Length + "M"));
        });
    }

    [Test]
    public void KSw2Extend_Validation()
    {
        Assert.Throws<ArgumentException>(() => AlignmentTools.Ksw2Extend("ACGT", "ACGT", variant: "nope"));
        Assert.Throws<ArgumentOutOfRangeException>(() => AlignmentTools.Ksw2Extend("ACGT", "ACGT", gapExtend: 0));
        Assert.Throws<ArgumentException>(() => AlignmentTools.Ksw2Extend(null!, "ACGT"));
        var empty = AlignmentTools.Ksw2Extend("", "ACGT");
        Assert.That((empty.Max, empty.Cigar), Is.EqualTo((0, "*")), "an empty query is the ksw2 no-op: nothing is aligned");
    }

    [Test]
    public void FormatParameter_PairwiseTools_WriteTheCanonicalFormats()
    {
        var plain = AlignmentTools.GlobalAlign("ACGT", "ACGT");
        Assert.That(plain.Formatted, Is.Null);
        var fasta = AlignmentTools.GlobalAlign("ACGT", "AC-T".Replace("-", "G"), format: "fasta", names: ["x", "y"]);
        Assert.That(fasta.Formatted, Is.EqualTo(">x\nACGT\n>y\nACGT\n"));
        var clustal = AlignmentTools.GlobalAlign("ACGT", "ACGT", format: "clustal");
        Assert.That(clustal.Formatted, Does.StartWith("Biopython 1.88 multiple sequence alignment\n\n\nseq1"));
        var pair = AlignmentTools.LocalAlignAffine("ACGTACGT", "ACGTACGT", format: "pair", names: ["a", "b"]);
        Assert.That(pair.Formatted, Does.StartWith("########################################\n# Program: water\n"), "a local result is a water report (R14-1)");
        Assert.That(AlignmentTools.GlobalAlignAffine("ACGT", "ACGT", format: "pair").Formatted, Does.StartWith("########################################\n# Program: needle\n"));
        var emboss = AlignmentTools.GlobalAlign("ACGT", "ACGT", format: "emboss:fasta", names: ["a", "b"]);
        Assert.That(emboss.Formatted, Is.EqualTo(AlignmentTools.WriteAlignment(["ACGT", "ACGT"], "fasta", ["a", "b"]).Text));
        Assert.Throws<ArgumentException>(() => AlignmentTools.GlobalAlign("ACGT", "ACGT", format: "bogus"));
    }

    [Test]
    public void FormatParameter_MultipleAlignmentTools_WriteTheCanonicalFormats()
    {
        var r = AlignmentTools.MultipleAlign(["ACGT", "ACGA", "ACGT"], format: "phylip");
        Assert.That(r.Formatted, Does.StartWith("3 4\nseq1      "));
        var s = AlignmentTools.MultipleAlign(["ACGT", "ACGA", "ACGT"], format: "stockholm", names: ["a", "b", "c"]);
        Assert.That(s.Formatted, Does.StartWith("# STOCKHOLM 1.0\n#=GF SQ   3\n"));
        var m = AlignmentTools.MultipleAlign(["ACGT", "ACGA", "ACGT"], format: "msf");
        Assert.That(m.Formatted, Does.Contain("MSF: 4"));
        Assert.That(AlignmentTools.MultipleAlign(["ACGT", "ACGA"]).Formatted, Is.Null);
    }

    [Test]
    public void WriteAlignment_BiopythonDialect_WritesTheBioAlignFiles()
    {
        var r = AlignmentTools.WriteAlignment(["AC-T", "ACGT"], "fasta", ["s1", "s2"], descriptions: ["first", ""], dialect: "biopython");
        Assert.That(r.Dialect, Is.EqualTo("biopython"));
        Assert.That(r.Text, Is.EqualTo(">s1 first\nAC-T\n>s2\nACGT\n"));
        var nexus = AlignmentTools.WriteAlignment(["AC-T", "ACGT"], "nexus", ["s1", "s2"], dialect: "biopython");
        Assert.That(nexus.Text, Is.EqualTo("#NEXUS\nbegin data;\ndimensions ntax=2 nchar=4;\nformat datatype=dna missing=? gap=-;\nmatrix\ns1 AC-T\ns2 ACGT\n;\nend;\n"));
        Assert.Throws<ArgumentException>(() => AlignmentTools.WriteAlignment(["AC-T"], "pair", dialect: "biopython"));
        Assert.Throws<ArgumentException>(() => AlignmentTools.WriteAlignment(["AC-T", "ACGT"], "fasta", dialect: "other"));
    }
}
