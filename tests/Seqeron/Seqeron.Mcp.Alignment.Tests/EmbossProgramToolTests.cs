using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;

namespace Seqeron.Mcp.Alignment.Tests;

/// <summary>supermatcher_align / needleall_align / ksw2_extend and the `format` parameter of the pairwise and multiple-alignment tools (audit A2-27b / A2-28c / A2-29b, B09 F124 - F126). The tools delegate to
/// EmbossPairwisePrograms.Supermatcher / NeedleAll, Ksw2Aligner and AlignmentFormats / EmbossAlignmentWriter; the oracle comparisons live in the library tests.</summary>
[TestFixture]
public class EmbossProgramToolTests
{
    private const string Q = "ACGTACGTACGGTTAACCGGTTAAACCC";

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
}
