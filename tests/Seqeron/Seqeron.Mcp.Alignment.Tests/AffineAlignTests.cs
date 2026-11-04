using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;

namespace Seqeron.Mcp.Alignment.Tests;

/// <summary>
/// MCP wrappers global_align_affine / local_align_affine / semi_global_align_affine delegate to
/// SequenceAligner.*Affine (gap of length k scores gapOpen + k·gapExtend). Expected values are the
/// Biopython 1.88 PairwiseAligner (open_gap_score = gapOpen + gapExtend, extend_gap_score = gapExtend)
/// = parasail 1.3.4 reference scores recorded in the SequenceAligner unit tests (B09 F1-F3).
/// </summary>
[TestFixture]
public class AffineAlignTests
{
    [Test]
    public void GlobalAlignAffine_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AlignmentTools.GlobalAlignAffine("ACGT", "ACGT"));
        Assert.Throws<ArgumentException>(() => AlignmentTools.GlobalAlignAffine("", "ACGT"));
        Assert.Throws<ArgumentException>(() => AlignmentTools.GlobalAlignAffine("ACGT", null!));
    }

    [Test]
    public void GlobalAlignAffine_ChargesGapOpen_UnlikeLinearTool()
    {
        // TTTTTTTTTTACGT / ACGT under blastn scoring (2/-3/-5/-2): 4 matches (+8) and one 10-gap
        // (-5 + 10*-2 = -25) = -17 (Biopython = parasail). The linear tool ignores gapOpen: 8 - 20 = -12.
        var affine = AlignmentTools.GlobalAlignAffine("TTTTTTTTTTACGT", "ACGT", 2, -3, -5, -2);
        var linear = AlignmentTools.GlobalAlign("TTTTTTTTTTACGT", "ACGT", 2, -3, -5, -2);

        Assert.Multiple(() =>
        {
            Assert.That(affine.Score, Is.EqualTo(-17));
            Assert.That(affine.AlignmentType, Is.EqualTo("Global"));
            Assert.That(affine.AlignedSequence1, Is.EqualTo("TTTTTTTTTTACGT"));
            Assert.That(affine.AlignedSequence2, Is.EqualTo("----------ACGT"));
            Assert.That(linear.Score, Is.EqualTo(-12));
        });
    }

    [Test]
    public void LocalAlignAffine_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AlignmentTools.LocalAlignAffine("ACGT", "ACGT"));
        Assert.Throws<ArgumentException>(() => AlignmentTools.LocalAlignAffine("ACGT", ""));
    }

    [Test]
    public void LocalAlignAffine_ReferenceValue_MatchesBiopythonAndParasail()
    {
        // GGGACGTTTACGTCCC / TTACGTACGTAA, blastn scoring -> 12 at seq1[7..12], seq2[0..5].
        var r = AlignmentTools.LocalAlignAffine("GGGACGTTTACGTCCC", "TTACGTACGTAA", 2, -3, -5, -2);

        Assert.Multiple(() =>
        {
            Assert.That(r.Score, Is.EqualTo(12));
            Assert.That(r.AlignmentType, Is.EqualTo("Local"));
            Assert.That((r.StartPosition1, r.EndPosition1, r.StartPosition2, r.EndPosition2), Is.EqualTo((7, 12, 0, 5)));
        });
    }

    [Test]
    public void SemiGlobalAlignAffine_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AlignmentTools.SemiGlobalAlignAffine("ACGT", "ACGTAC"));
        Assert.Throws<ArgumentException>(() => AlignmentTools.SemiGlobalAlignAffine("ACXT", "ACGT"));
        Assert.Throws<ArgumentException>(() => AlignmentTools.SemiGlobalAlignAffine("ACGT", "AC-T"));
    }

    [Test]
    public void SemiGlobalAlignAffine_ReferenceValue_MatchesBiopythonAndParasail()
    {
        // GATTACA fitted into CCGATTTACAGG, blastn scoring -> 7, --GA-TTACA-- / CCGATTTACAGG.
        var r = AlignmentTools.SemiGlobalAlignAffine("GATTACA", "CCGATTTACAGG", 2, -3, -5, -2);

        Assert.Multiple(() =>
        {
            Assert.That(r.Score, Is.EqualTo(7));
            Assert.That(r.AlignmentType, Is.EqualTo("SemiGlobal"));
            Assert.That(r.AlignedSequence1, Is.EqualTo("--GA-TTACA--"));
            Assert.That(r.AlignedSequence2, Is.EqualTo("CCGATTTACAGG"));
        });
    }

    [Test]
    public void LinearTools_IgnoreGapOpen_AsDocumented()
    {
        // The linear tools must give identical results for any gapOpen (documented in their descriptions).
        var a = AlignmentTools.GlobalAlign("ACGTTTTACGT", "ACGTACGT", 2, -3, -5, -2);
        var b = AlignmentTools.GlobalAlign("ACGTTTTACGT", "ACGTACGT", 2, -3, 0, -2);
        Assert.That(b, Is.EqualTo(a));
    }
}
