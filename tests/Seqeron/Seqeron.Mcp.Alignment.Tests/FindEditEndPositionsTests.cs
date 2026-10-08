using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;

namespace Seqeron.Mcp.Alignment.Tests;

/// <summary>
/// find_edit_end_positions — delegation to <c>ApproximateMatcher.FindEditEndPositions</c>.
/// Expected values: Navarro 2001 survey/surgery and edlib HW cases locked in
/// ApproximateMatcher_EditDistance_Tests.
/// </summary>
[TestFixture]
public class FindEditEndPositionsTests
{
    private static (int, int)[] Items(EditEndPositionsResult r) => r.Items.Select(i => (i.EndPosition, i.Distance)).ToArray();

    [Test]
    public void FindEditEndPositions_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AlignmentTools.FindEditEndPositions("ACGTACGT", "ACG", 1));
        Assert.Throws<ArgumentException>(() => AlignmentTools.FindEditEndPositions("", "ACG", 1));
        Assert.Throws<ArgumentException>(() => AlignmentTools.FindEditEndPositions("ACGT", null!, 1));
        Assert.Throws<ArgumentOutOfRangeException>(() => AlignmentTools.FindEditEndPositions("ACGT", "ACG", -1));
    }

    [Test]
    public void FindEditEndPositions_NavarroSurveySurgery_ReturnsSellersEnds()
    {
        Assert.That(Items(AlignmentTools.FindEditEndPositions("surgery", "survey", 2)),
            Is.EqualTo(new[] { (4, 2), (5, 2), (6, 2) }));
    }

    [Test]
    public void FindEditEndPositions_TtacInGattaca_MatchesEdlib_AndEqualsLibrary()
    {
        var r = AlignmentTools.FindEditEndPositions("GATTACAGATTTACA", "TTAC", 1);
        Assert.That(Items(r), Is.EqualTo(new[] { (4, 1), (5, 0), (6, 1), (12, 1), (13, 0), (14, 1) }));
        Assert.That(Items(r), Is.EqualTo(global::Seqeron.Genomics.Alignment.ApproximateMatcher
            .FindEditEndPositions("GATTACAGATTTACA", "TTAC", 1).ToArray()));
    }

    // Weighted Sellers (B05 F27): Python brute force min_i rapidfuzz Levenshtein.distance(p, t[i..j], weights).
    [Test]
    public void FindEditEndPositions_WeightedCosts_BruteForceValues()
    {
        var r = AlignmentTools.FindEditEndPositions("TTACGTAAGGCTACG", "ACGT", 2, insertionCost: 1, deletionCost: 3, substitutionCost: 1);
        Assert.Multiple(() =>
        {
            Assert.That(r.Items.Select(i => (i.EndPosition, i.Distance)),
                Is.EqualTo(new[] { (5, 0), (6, 1), (7, 2), (9, 2), (10, 2), (11, 2) }));
            Assert.Throws<ArgumentOutOfRangeException>(() => AlignmentTools.FindEditEndPositions("ACGT", "AC", 1, insertionCost: -1));
        });
    }
}
