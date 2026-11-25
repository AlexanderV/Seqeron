using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;

namespace Seqeron.Mcp.Alignment.Tests;

/// <summary>
/// B09 A2-36: the MCP tools of the public stage ports (ClustalW pairwise stages, T-Coffee pair methods / primary library / guide tree / k-tuple similarity, sum of pairs).
/// Every expected value is copied from the reference-locked library tests (clustalw 2.1 / t_coffee 13.41 hand runs and oracles) or computed by hand.
/// </summary>
[TestFixture]
public class PortStageToolsTests
{
    private const string HbA = "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHF";
    private const string HbB = "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKV";
    private const string HbC = "MGLSDGEWQLVLNVWGKVEADIPGHGQEVLIRLFKGHPETLEKFDKFKHLKTEAEMKASEDLKKHG";

    [Test]
    public void SumOfPairsScore_LinearAffineAndMatrix()
    {
        // 3 identical rows x 2 columns: 3 pairs x 2 columns x match 1 = 6 (linear); BLOSUM62 A/A = 4 -> 24 (SequenceAligner_ScoringContract_Tests)
        string[] rows = ["AA", "AA", "AA"];
        Assert.That(AlignmentTools.SumOfPairsScore(rows).Score, Is.EqualTo(6));
        Assert.That(AlignmentTools.SumOfPairsScore(rows, matrix: "BLOSUM62", gapOpen: -1, gapExtend: -1).Score, Is.EqualTo(24.0));
        // AC- / AC- / ACA with match 2, mismatch -1, gap -2 per position (linear): pairs (0,1) 2+2+0(gap-gap) = 4; (0,2) 2+2-2 = 2; (1,2) 2; total 8
        Assert.That(AlignmentTools.SumOfPairsScore(["AC-", "AC-", "ACA"], match: 2, mismatch: -1, gapOpen: -3, gapExtend: -2).Score, Is.EqualTo(8));
        // affine: pair (0,2): the single gap column costs gapOpen + 1 x gapExtend = -3 + -2 = -5 -> 4 + (4 - 5) + (4 - 5) = 2
        var affine = AlignmentTools.SumOfPairsScore(["AC-", "AC-", "ACA"], match: 2, mismatch: -1, gapOpen: -3, gapExtend: -2, gapModel: "affine");
        Assert.That(affine.Score, Is.EqualTo(2));
        Assert.That(affine.GapModel, Is.EqualTo("Affine"));
        Assert.Throws<ArgumentException>(() => AlignmentTools.SumOfPairsScore(["AA", "A"]));
        Assert.Throws<ArgumentException>(() => AlignmentTools.SumOfPairsScore(rows, gapOpen: -0.5));
        Assert.Throws<ArgumentException>(() => AlignmentTools.SumOfPairsScore(rows, gapModel: "x"));
        Assert.Throws<ArgumentException>(() => AlignmentTools.SumOfPairsScore([]));
    }
}
