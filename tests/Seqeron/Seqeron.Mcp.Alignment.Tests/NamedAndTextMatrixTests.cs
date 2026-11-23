using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;

namespace Seqeron.Mcp.Alignment.Tests;

/// <summary>
/// R14-8 (B09 F301): the pairwise tools take every named matrix of Biopython 1.88 (fractional ones included, scored exactly) and the TEXT of a matrix in the <c>matrix</c> parameter. Expected values:
/// Biopython 1.88 <c>PairwiseAligner</c> (GONNET1992, gap open -2.5 / extend -0.5 = open_gap_score -3, extend_gap_score -0.5: global 39.4, local 42.4; BLASTP, MKTAUYIAK / MKTAAYIAK, 11 / 1: 39.0).
/// </summary>
[TestFixture]
public class NamedAndTextMatrixTests
{
    [Test]
    public void FractionalNamedMatrix_GlobalAndLocal_EqualBiopython()
    {
        var g = AlignmentTools.GlobalAlignAffine("ACDEFGHIK", "ACDFGHIKL", gapOpen: -2.5, gapExtend: -0.5, matrix: "GONNET1992");
        Assert.That(g.ExactScore, Is.EqualTo(39.4).Within(1e-9));
        var l = AlignmentTools.LocalAlignAffine("ACDEFGHIK", "ACDFGHIKL", gapOpen: -2.5, gapExtend: -0.5, matrix: "gonnet1992");
        Assert.That(l.ExactScore, Is.EqualTo(42.4).Within(1e-9));
    }

    [Test]
    public void Blastp_ScoresSelenocysteineLikeBlastp()
    {
        var r = AlignmentTools.GlobalAlignAffine("MKTAUYIAK", "MKTAAYIAK", gapOpen: -11, gapExtend: -1, matrix: "BLASTP");
        Assert.That(r.Score, Is.EqualTo(39));
        Assert.Throws<ArgumentException>(() => AlignmentTools.GlobalAlignAffine("MKTAUYIAK", "MKTAAYIAK", gapOpen: -11, gapExtend: -1, matrix: "BLOSUM62"), "BLOSUM62 has no U");
    }

    [Test]
    public void MatrixText_EqualsTheNamedMatrix()
    {
        const string text = "# my DNA matrix\n   A  C  G  T\nA  5 -4 -4 -4\nC -4  5 -4 -4\nG -4 -4  5 -4\nT -4 -4 -4  5\n";
        var viaText = AlignmentTools.GlobalAlignAffine("ACGTACGT", "ACGACGT", gapOpen: -10, gapExtend: -0.5, matrix: text);
        // the same cells as EDNAFULL on A / C / G / T
        var viaName = AlignmentTools.GlobalAlignAffine("ACGTACGT", "ACGACGT", gapOpen: -10, gapExtend: -0.5, matrix: "EDNAFULL");
        Assert.That(viaText.ExactScore, Is.EqualTo(viaName.ExactScore));
        Assert.That(viaText.AlignedSequence1, Is.EqualTo(viaName.AlignedSequence1));
        var fractional = AlignmentTools.LocalAlign("AC", "AC", matrix: "  A  C\nA 1.5 -1\nC -1 2.25\n");
        Assert.That(fractional.ExactScore, Is.EqualTo(3.75));
    }

    [Test]
    public void BadMatrixText_IsRefused()
    {
        Assert.Throws<ArgumentException>(() => AlignmentTools.GlobalAlign("AC", "AC", matrix: "A C\nA 1 2\n"));
        Assert.Throws<ArgumentException>(() => AlignmentTools.GlobalAlign("AC", "AC", matrix: "NOSUCH"));
    }
}
