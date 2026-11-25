using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;

namespace Seqeron.Mcp.Alignment.Tests;

/// <summary>stretcher_align / write_alignment: EMBOSS 6.6 stretcher and the -aformat output formats (audit A2-27 / A2-28, B09 F120 / F121). The tools delegate to EmbossPairwisePrograms.StretcherAlign and EmbossAlignmentWriter; the oracle comparisons live in the library tests.</summary>
[TestFixture]
public class EmbossFormatTests
{
    private const string A = "ACGTACGTAGCTAGCTAGCATCGATCGATCGAT";
    private const string B = "ACGTACGTAGCTAGATAGCATCGAATCGATCGAT";

    [Test]
    public void DualAffineAlign_EqualsTheLibrary_AndValidates()
    {
        var tool = AlignmentTools.DualAffineAlign("ACGTACGTAGCTAGCTAGCATCG", "ACGTACGTAGCTAGATAGCATCG", "local");
        var lib = Seqeron.Genomics.Alignment.SequenceAligner.LocalAlignDualAffine("ACGTACGTAGCTAGCTAGCATCG", "ACGTACGTAGCTAGATAGCATCG", Seqeron.Genomics.Alignment.DualAffineScoring.Minimap2);
        Assert.That(tool.Score, Is.EqualTo(lib.Score));
        Assert.That(tool.AlignmentType, Is.EqualTo("Local"));
        Assert.Throws<ArgumentException>(() => AlignmentTools.DualAffineAlign("ACGT", "ACGT", "semi"));
        Assert.Throws<ArgumentException>(() => AlignmentTools.DualAffineAlign("ACGT", "ACGT", gapOpen1: 3));
        Assert.That(AlignmentTools.DualAffineAlign("ACGTACGT", "ACGTACGT", band: 0).Score, Is.EqualTo(16));
    }
}
