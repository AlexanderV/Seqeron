using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// R14-17 (B09 F298): the <c>find_longest_common_region_multi</c> MCP tool (multi-sequence longest common substring with minSupport, Gusfield 1997 sec. 7.6, Rosalind LCSM).
/// Expected values: an independent brute force over every substring of every sequence (distinct longest strings in ordinal order, first occurrence per sequence, -1 = absent), not the wrapper's output.
/// </summary>
[TestFixture]
public class FindLongestCommonRegionMultiTests
{
    [Test]
    public void RosalindLcsm_Sample_ReturnsEveryLongestString_WithFirstOccurrences()
    {
        var r = AnalysisTools.FindLongestCommonRegionMulti(["GATTACA", "TAGACCA", "ATACA"]);
        Assert.That(r.Length, Is.EqualTo(2));
        Assert.That(r.Items.Select(i => i.Sequence), Is.EqualTo(new[] { "AC", "CA", "TA" }));
        Assert.That(r.Items[0].Positions, Is.EqualTo(new[] { 4, 3, 2 }));
        Assert.That(r.Items[1].Positions, Is.EqualTo(new[] { 5, 5, 3 }));
        Assert.That(r.Items[2].Positions, Is.EqualTo(new[] { 3, 0, 1 }));
        Assert.That(r.Items.All(i => i.Support == 3 && i.Length == 2), Is.True);
    }

    [Test]
    public void MinSupport_TwoOfThree_AllowsALongerString_AndMarksTheAbsentSequence()
    {
        // GATTACA / TAGACCA / ATACA: the longest string in at least 2 sequences is longer than the one in all 3
        var r = AnalysisTools.FindLongestCommonRegionMulti(["ACGTACGTAA", "CCGTACGGTT", "TTTCGTACGA"], minSupport: 2);
        Assert.That(r.Items.Select(i => (i.Sequence, i.Positions.ToArray())), Is.EqualTo(new[] { ("CGTACG", new[] { 1, 1, 3 }) }));
        var partial = AnalysisTools.FindLongestCommonRegionMulti(["AAAACCCC", "GGGGTTTT", "AAAACGCG"], minSupport: 2);
        Assert.That(partial.Items.Single().Sequence, Is.EqualTo("AAAAC"));
        Assert.That(partial.Items.Single().Positions, Is.EqualTo(new[] { 0, -1, 0 }));
        Assert.That(partial.Items.Single().Support, Is.EqualTo(2));
    }

    [Test]
    public void NRuns_AndLowerCase_EndARegion()
    {
        var r = AnalysisTools.FindLongestCommonRegionMulti(["ACGTNACGT", "ggnacgtaa", "TACGTNN"]);
        Assert.That(r.Items.Select(i => (i.Sequence, i.Positions.ToArray())), Is.EqualTo(new[] { ("ACGT", new[] { 0, 3, 1 }) }));
    }

    [Test]
    public void NoCommonString_GivesNoItems()
    {
        var r = AnalysisTools.FindLongestCommonRegionMulti(["AAAA", "CCCC", "GGGG"]);
        Assert.That(r.Items, Is.Empty);
        Assert.That(r.Length, Is.EqualTo(0));
    }

    [Test]
    public void Validation()
    {
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindLongestCommonRegionMulti(["ACGT"]));
        Assert.Throws<ArgumentNullException>(() => AnalysisTools.FindLongestCommonRegionMulti(null!));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindLongestCommonRegionMulti(["ACGT", ""]));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindLongestCommonRegionMulti(["ACGT", "AC1T"]));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.FindLongestCommonRegionMulti(["ACGT", "ACGA", "ACGG"], minSupport: 1));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.FindLongestCommonRegionMulti(["ACGT", "ACGA", "ACGG"], minSupport: 4));
        using var cts = new CancellationTokenSource();
        cts.Cancel();
        Assert.Throws<OperationCanceledException>(() => AnalysisTools.FindLongestCommonRegionMulti(["ACGT", "ACGA"], cancellationToken: cts.Token));
    }
}
