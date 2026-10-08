using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>find_inverted_repeats_scored</c> MCP tool. Expected values = EMBOSS 6.6.0 einverted output
/// locked in RepeatFinder_InvertedRepeatOptions_Tests (Evidence REP-INV-001), NOT the wrapper's output.
/// </summary>
[TestFixture]
public class FindInvertedRepeatsScoredTests
{
    private const string PlantedImperfect =
        "CCCAACCCATGCGTACGTTAGCCTAGGATCCATTTTTTTTTGGATACTAGGCAACGTACGCATGGGAAGGG";

    [Test]
    public void FindInvertedRepeatsScored_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.FindInvertedRepeatsScored(PlantedImperfect));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindInvertedRepeatsScored(""));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindInvertedRepeatsScored(null!));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.FindInvertedRepeatsScored(PlantedImperfect, gapPenalty: -1));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.FindInvertedRepeatsScored(PlantedImperfect, mismatchScore: 1));
    }

    [Test]
    public void FindInvertedRepeatsScored_Binding_EinvertedWorkedExample()
    {
        // einverted defaults: "Score 60: 28/31 (90%) matches, 1 gaps", left 1..32, right 71..41 (1-based).
        var items = AnalysisTools.FindInvertedRepeatsScored(PlantedImperfect).Items;
        Assert.That(items, Has.Length.EqualTo(1));
        var r = items[0];
        Assert.Multiple(() =>
        {
            Assert.That((r.LeftArmStart, r.LeftArmEnd, r.RightArmStart, r.RightArmEnd), Is.EqualTo((0, 31, 40, 70)));
            Assert.That((r.Score, r.Matches, r.Mismatches, r.Gaps), Is.EqualTo((60, 28, 3, 1)));
            Assert.That(r.LeftArmAlignment, Is.EqualTo("CCCAACCCATGCGTACGTTAGCCTAGGATCCA"));
            Assert.That(r.MatchLine, Is.EqualTo("|||  |||||||||||||| |||||| |||||"));
            Assert.That(r.RightArmAlignment, Is.EqualTo("GGGAAGGGTACGCATGCAA-CGGATCATAGGT"));
            Assert.That((r.LeftArmLength, r.RightArmLength, r.LoopLength), Is.EqualTo((32, 31, 8)));
            Assert.That((int)r.PercentMatches, Is.EqualTo(90));
        });
    }

    [Test]
    public void FindInvertedRepeatsScored_CustomGap_MatchesEinverted()
    {
        // einverted -gap 8 -threshold 30: Score 64 (= 84 - 12 - 8).
        var r = AnalysisTools.FindInvertedRepeatsScored(PlantedImperfect, gapPenalty: 8, threshold: 30).Items.Single();
        Assert.That((r.Score, r.Matches, r.Mismatches, r.Gaps, r.LeftArmStart, r.RightArmEnd), Is.EqualTo((64, 28, 3, 1, 0, 70)));
    }
}
