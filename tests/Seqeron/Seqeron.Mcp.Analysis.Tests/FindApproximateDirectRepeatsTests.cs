using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>find_approximate_direct_repeats</c> MCP tool. Expected values = Vmatch 2.3.1 -h k (-allmax),
/// locked in RepeatFinder_RepeatVariants_Tests. NOT the wrapper's output.
/// </summary>
[TestFixture]
public class FindApproximateDirectRepeatsTests
{
    [Test]
    public void FindApproximateDirectRepeats_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.FindApproximateDirectRepeats("GATTACAGATTACA", 6));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindApproximateDirectRepeats(""));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindApproximateDirectRepeats(null!));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.FindApproximateDirectRepeats("ACGT", 1, 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.FindApproximateDirectRepeats("ACGT", 5, -1));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.FindApproximateDirectRepeats("ACGT", 5, 5));
    }

    [Test]
    public void FindApproximateDirectRepeats_Binding_TwoMismatches_MatchesVmatch()
    {
        // vmatch -l 10 -h 2: (0, 21, 17) with 2 mismatches.
        var items = AnalysisTools.FindApproximateDirectRepeats(
            "ACGTTGCAAGCTTACGGGGGGACGATGCAAGCATACGG", 10, 2, int.MaxValue, int.MinValue).Items;
        Assert.That(items, Has.Length.EqualTo(1));
        Assert.Multiple(() =>
        {
            Assert.That((items[0].FirstPosition, items[0].SecondPosition, items[0].Length, items[0].Mismatches), Is.EqualTo((0, 21, 17, 2)));
            Assert.That(items[0].FirstCopy, Is.EqualTo("ACGTTGCAAGCTTACGG"));
            Assert.That(items[0].SecondCopy, Is.EqualTo("ACGATGCAAGCATACGG"));
            Assert.That(items[0].Spacing, Is.EqualTo(4));
        });
    }

    [Test]
    public void FindApproximateDirectRepeats_ExcludeContained_EqualsVmatchAllmax()
    {
        var items = AnalysisTools.FindApproximateDirectRepeats(
            "AAAAAAAACGTTGCAACGTAAAA", 5, 1, int.MaxValue, int.MinValue, excludeContained: true).Items;
        Assert.That(string.Join(";", items.Select(x => $"{x.FirstPosition},{x.SecondPosition},{x.Length},{x.Mismatches}")),
            Is.EqualTo("0,1,8,1;0,18,5,1;1,18,5,1;2,18,5,1;3,18,5,1;5,13,6,1;6,14,6,1"));
    }
}
