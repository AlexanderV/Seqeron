using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;

namespace Seqeron.Mcp.Alignment.Tests;

/// <summary>R14-12 (B09 F305): the <c>banded_align</c> tool equals the full tools when the optimum lies in the band and proves the optimum in the exact mode.</summary>
[TestFixture]
public class BandedAlignToolTests
{
    private const string A = "ACGTACGTTAGCTAGCTAGGATCCGATCGATTACGATCGATCGA";
    private const string B = "ACGTACGTAGCTAGCTAGGATCGATCGATTACGGATCGATCGA";

    [Test]
    public void GlobalExact_EqualsGlobalAlign_WithAProof()
    {
        var banded = AlignmentTools.BandedAlignPair(A, B, mode: "global", match: 2, mismatch: -3, gapExtend: -2);
        var full = AlignmentTools.GlobalAlign(A, B, match: 2, mismatch: -3, gapExtend: -2);
        Assert.That(banded.Proven, Is.True);
        Assert.That(banded.Alignment.AlignedSequence1, Is.EqualTo(full.AlignedSequence1));
        Assert.That(banded.Alignment.AlignedSequence2, Is.EqualTo(full.AlignedSequence2));
        Assert.That(banded.Alignment.Score, Is.EqualTo(full.Score));
        Assert.That(banded.Cells, Is.GreaterThan(0));
    }

    [Test]
    public void GlobalAffineBand_LocalAndSemiGlobal_EqualTheFullTools()
    {
        var g = AlignmentTools.BandedAlignPair(A, B, affine: true, match: 2, mismatch: -3, gapOpen: -4, gapExtend: -2, bandWidth: 8);
        var fullG = AlignmentTools.GlobalAlignAffine(A, B, match: 2, mismatch: -3, gapOpen: -4, gapExtend: -2);
        Assert.That(g.Alignment.AlignedSequence1, Is.EqualTo(fullG.AlignedSequence1));
        Assert.That(g.Alignment.Score, Is.EqualTo(fullG.Score));
        Assert.That(g.Proven, Is.Null);
        var l = AlignmentTools.BandedAlignPair(A, B, mode: "local", match: 2, mismatch: -3, gapExtend: -2, diagonalLow: -10, diagonalHigh: 10);
        var fullL = AlignmentTools.LocalAlign(A, B, match: 2, mismatch: -3, gapExtend: -2);
        Assert.That(l.Alignment.Score, Is.EqualTo(fullL.Score));
        var s = AlignmentTools.BandedAlignPair("ACGTAC", "TTTACGTACTTT", mode: "semi_global", diagonalLow: -2, diagonalHigh: 8);
        var fullS = AlignmentTools.SemiGlobalAlign("ACGTAC", "TTTACGTACTTT");
        Assert.That(s.Alignment.AlignedSequence2, Is.EqualTo(fullS.AlignedSequence2));
    }

    [Test]
    public void Refusals()
    {
        Assert.Throws<ArgumentException>(() => AlignmentTools.BandedAlignPair(A, B, mode: "local"));
        Assert.Throws<ArgumentException>(() => AlignmentTools.BandedAlignPair(A, B, mode: "banded"));
        Assert.Throws<ArgumentException>(() => AlignmentTools.BandedAlignPair(A, B, exact: false));
        Assert.Throws<ArgumentOutOfRangeException>(() => AlignmentTools.BandedAlignPair(A, B, bandWidth: -1));
    }
}
