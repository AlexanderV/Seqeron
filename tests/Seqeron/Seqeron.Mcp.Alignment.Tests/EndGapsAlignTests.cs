using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;

namespace Seqeron.Mcp.Alignment.Tests;

/// <summary>
/// MCP semi-global tools with the additive <c>endGaps</c> parameter (B09 F44, audit A1-3): the 16 end-gap policies delegate to
/// <c>SequenceAligner.SemiGlobalAlign(Affine)(…, EndGaps)</c>. Expected rows/scores: Biopython 1.88 PairwiseAligner
/// (left_deletion_score = right_insertion_score = 0 etc.) on a unique-optimum overlap example.
/// </summary>
[TestFixture]
public class EndGapsAlignTests
{
    [Test]
    public void DefaultEndGaps_IsTheClassicFitting()
    {
        var def = AlignmentTools.SemiGlobalAlign("GATTACA", "CCGATTTACAGG", 2, -3, -2, -2);
        var fit = AlignmentTools.SemiGlobalAlign("GATTACA", "CCGATTTACAGG", 2, -3, -2, -2, endGaps: "Fitting");
        Assert.Multiple(() =>
        {
            Assert.That(fit, Is.EqualTo(def));
            Assert.That((def.StartPosition2, def.EndPosition2), Is.EqualTo((2, 9)));
            Assert.That(def.AlignmentType, Is.EqualTo("SemiGlobal"));
        });
    }

    [Test]
    public void Overlap_Linear_SuffixOfSequence1_PrefixOfSequence2()
    {
        var r = AlignmentTools.SemiGlobalAlign("GGGGACGTAC", "ACGTACTTTT", 1, -1, 0, -1, endGaps: "OverlapSuffix1Prefix2");
        Assert.Multiple(() =>
        {
            Assert.That(r.AlignedSequence1, Is.EqualTo("GGGGACGTAC----"));
            Assert.That(r.AlignedSequence2, Is.EqualTo("----ACGTACTTTT"));
            Assert.That(r.Score, Is.EqualTo(6));
            Assert.That((r.StartPosition1, r.EndPosition1, r.StartPosition2, r.EndPosition2), Is.EqualTo((4, 9, 0, 5)));
        });
    }

    [Test]
    public void FlagList_EqualsThePreset_AndOtherOrientationFindsTheSameOverlap()
    {
        var preset = AlignmentTools.SemiGlobalAlignAffine("ACGTACTTTT", "GGGGACGTAC", 1, -1, 0, -1, endGaps: "overlapsuffix2prefix1");
        var flags = AlignmentTools.SemiGlobalAlignAffine("ACGTACTTTT", "GGGGACGTAC", 1, -1, 0, -1, endGaps: "FreeEnd1, FreeStart2");
        Assert.Multiple(() =>
        {
            Assert.That(flags, Is.EqualTo(preset));
            Assert.That(preset.AlignedSequence1, Is.EqualTo("----ACGTACTTTT"));
            Assert.That(preset.AlignedSequence2, Is.EqualTo("GGGGACGTAC----"));
            Assert.That(preset.Score, Is.EqualTo(6));
        });
    }

    [Test]
    public void None_IsGlobal_TaggedGlobal()
    {
        var none = AlignmentTools.SemiGlobalAlignAffine("ACGTT", "ACGGTT", 2, -3, -3, -1, endGaps: "None");
        var global = AlignmentTools.GlobalAlignAffine("ACGTT", "ACGGTT", 2, -3, -3, -1);
        Assert.That(none, Is.EqualTo(global));
        Assert.That(none.AlignmentType, Is.EqualTo("Global"));
    }

    [Test]
    public void All_WithMatrix_IsNeedleDefaultEndweight_AndProteinInputIsAccepted()
    {
        // all four end gaps free: a missing N-terminal / C-terminal piece costs nothing
        var r = AlignmentTools.SemiGlobalAlignAffine("WWWWMKVLAAGIVGLLLA", "MKVLAAGIVGLLLAPPPP", gapOpen: -9.5, gapExtend: -0.5, matrix: "BLOSUM62", endGaps: "All");
        Assert.Multiple(() =>
        {
            Assert.That(r.AlignedSequence1, Is.EqualTo("WWWWMKVLAAGIVGLLLA----"));
            Assert.That(r.AlignedSequence2, Is.EqualTo("----MKVLAAGIVGLLLAPPPP"));
            Assert.That(r.ExactScore, Is.GreaterThan(60));
        });
    }

    [TestCase("Sideways")]
    [TestCase("7")]
    [TestCase("FreeStart1,Nope")]
    public void UnknownEndGaps_IsRejected(string value)
    {
        Assert.Throws<ArgumentException>(() => AlignmentTools.SemiGlobalAlign("ACGT", "ACGT", endGaps: value));
        Assert.Throws<ArgumentException>(() => AlignmentTools.SemiGlobalAlignAffine("ACGT", "ACGT", endGaps: value));
    }
}
