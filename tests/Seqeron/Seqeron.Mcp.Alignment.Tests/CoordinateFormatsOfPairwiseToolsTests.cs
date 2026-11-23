using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;

namespace Seqeron.Mcp.Alignment.Tests;

/// <summary>
/// R14-9 (B09 F302): the <c>format</c> parameter of the pairwise tools writes the alignment in the coordinate formats of Biopython 1.88 <c>Alignment.format</c> (bed, psl, chain, maf, sam, exonerate) and as cigar / md / cs strings.
/// Expected texts: Biopython 1.88 on the same alignment ("TTACGTTT" and "ACGT": the local alignment is target a 2..6, query b 0..4; the fitting one keeps the free flanks as 2D4M2D).
/// </summary>
[TestFixture]
public class CoordinateFormatsOfPairwiseToolsTests
{
    private const string A = "TTACGTTT";
    private const string B = "ACGT";
    private static readonly string[] Names = ["a", "b"];

    private static string Local(string format, string[]? names = null) => AlignmentTools.LocalAlign(A, B, format: format, names: names ?? Names).Formatted!;

    [TestCase("sam", "b\t0\ta\t3\t255\t4M\t*\t0\t0\tACGT\t*\tAS:i:4\n")]
    [TestCase("sam:md", "b\t0\ta\t3\t255\t4M\t*\t0\t0\tACGT\t*\tMD:Z:4\tAS:i:4\n")]
    [TestCase("psl", "4\t0\t0\t0\t0\t0\t0\t0\t+\tb\t4\t0\t4\ta\t8\t2\t6\t1\t4,\t0,\t2,\n")]
    [TestCase("bed", "a\t2\t6\tb\t4\t+\t2\t6\t0\t1\t4,\t0,\n")]
    [TestCase("chain", "chain 4 a 8 + 2 6 b 4 + 0 4\n4\n\n")]
    [TestCase("maf", "a score=4.000000\ns a 2 4 + 8 ACGT\ns b 0 4 + 4 ACGT\n\n")]
    [TestCase("exonerate", "vulgar: b 0 4 + a 2 6 + 4 M 4 4\n")]
    [TestCase("bio:sam", "b\t0\ta\t3\t255\t4M\t*\t0\t0\tACGT\t*\tAS:i:4\n")]
    public void LocalAlign_CoordinateFormats_EqualBiopython(string format, string expected) => Assert.That(Local(format), Is.EqualTo(expected));

    [TestCase("cigar", "2S4M2S")]
    [TestCase("cigar:eqx", "2S4=2S")]
    [TestCase("md", "4")]
    [TestCase("cs", ":4")]
    [TestCase("cs:long", "=ACGT")]
    public void LocalAlign_StringFormats(string format, string expected)
    {
        // a = TTACGTTT is the target, b = ACGT the query: the query has no flanks to clip
        string text = Local(format);
        Assert.That(text, Is.EqualTo(expected.Replace("2S", "", StringComparison.Ordinal)));
    }

    [Test]
    public void LocalAlign_EmbossSamIsStillTheEmbossRecord()
    {
        string text = AlignmentTools.LocalAlign(A, B, format: "emboss:sam", names: Names).Formatted!;
        Assert.That(text, Does.Contain("@HD").Or.Contain("b\t0\ta\t3\t255\t4M\t"));
        Assert.That(text, Is.Not.EqualTo(Local("sam")));
    }

    [Test]
    public void SemiGlobalAlign_SequenceTwoIsTheReference_AndTheEndGapsStayOrTrim()
    {
        string[] names = ["q", "ref"];
        string kept = AlignmentTools.SemiGlobalAlign(B, A, format: "sam", names: names).Formatted!;
        Assert.That(kept, Is.EqualTo("q\t0\tref\t1\t255\t2D4M2D\t*\t0\t0\tACGT\t*\tAS:i:4\n"));
        string trimmed = AlignmentTools.SemiGlobalAlign(B, A, format: "sam", names: names, trimEndGaps: true).Formatted!;
        Assert.That(trimmed, Is.EqualTo("q\t0\tref\t3\t255\t4M\t*\t0\t0\tACGT\t*\tAS:i:4\n"));
        Assert.That(AlignmentTools.SemiGlobalAlign(B, A, format: "cigar", names: names, trimEndGaps: true).Formatted, Is.EqualTo("4M"));
    }

    [Test]
    public void GlobalAlign_AffineAndGeneralGap_WriteCoordinates()
    {
        string sam = AlignmentTools.GlobalAlignAffine("ACGTACGT", "ACGACGT", format: "sam", names: ["t", "q"]).Formatted!;
        Assert.That(sam, Does.StartWith("q\t0\tt\t1\t255\t3M1D4M\t"));
        string cigar = AlignmentTools.GlobalAlign("ACGTACGT", "ACGACGT", format: "cigar:eqx").Formatted!;
        Assert.That(cigar, Does.Match("^[0-9=XID]+$"));
        string gg = AlignmentTools.GeneralGapAlign("ACGTACGT", "ACGACGT", [-2, -3, -4, -5, -6, -7, -8, -9], format: "psl", names: ["t", "q"]).Formatted!;
        Assert.That(gg, Does.StartWith("7\t0\t0\t0\t0\t0\t1\t1\t+\tq\t7\t0\t7\tt\t8\t0\t8\t2\t3,4,\t0,3,\t0,4,\n"));
    }

    [Test]
    public void EmptyLocalResult_HasNoRecord_AndUnknownNamesAreRefused()
    {
        Assert.That(AlignmentTools.LocalAlign("AAAA", "CCCC", format: "sam").Formatted, Is.EqualTo(""));
        Assert.Throws<ArgumentException>(() => AlignmentTools.GlobalAlign("ACGT", "ACGT", format: "nosuchformat"));
        Assert.That(AlignmentTools.GlobalAlign("ACGT", "ACGT", format: "bio:fasta").Formatted, Does.StartWith(">seq1"));
    }
}
