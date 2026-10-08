using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>find_reverse_complement_repeats</c> MCP tool. Expected values = mummer4 repeat-match (reverse
/// lines) / Vmatch -p, locked in RepeatFinder_RepeatVariants_Tests. NOT the wrapper's output.
/// </summary>
[TestFixture]
public class FindReverseComplementRepeatsTests
{
    [Test]
    public void FindReverseComplementRepeats_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.FindReverseComplementRepeats("GAATTCAAAAAGAATTC", 6));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindReverseComplementRepeats(""));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindReverseComplementRepeats(null!));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.FindReverseComplementRepeats("ACGT", 1));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.FindReverseComplementRepeats("ACGT", 5, 4));
    }

    [Test]
    public void FindReverseComplementRepeats_Binding_MatchesRepeatMatchAndVmatch()
    {
        // repeat-match: 8 19r 12 / 7 12r 6 / 16 19r 4; Vmatch -p -l 3: (7,7,12) (15,15,4) (6,6,6).
        var items = AnalysisTools.FindReverseComplementRepeats("AAAAAAAACGTTGCAACGTAAAA", 3, int.MaxValue, int.MinValue).Items;
        Assert.That(items.Select(r => (r.FirstPosition, r.SecondPosition, r.Length)),
            Is.EqualTo(new[] { (6, 6, 6), (7, 7, 12), (15, 15, 4) }));
        Assert.That(items[1].RepeatSequence, Is.EqualTo("ACGTTGCAACGT"));
        Assert.That(items[1].Spacing, Is.EqualTo(-12));
    }

    [Test]
    public void FindReverseComplementRepeats_Defaults_SeparatedCopiesOnly()
    {
        var items = AnalysisTools.FindReverseComplementRepeats("TTGCATGCAAAAAATTTTTTTGCATGCAA").Items;
        Assert.That(items.Select(r => (r.FirstPosition, r.SecondPosition, r.Length, r.Spacing)),
            Is.EqualTo(new[] { (0, 15, 14, 1), (8, 14, 5, 1), (9, 16, 5, 2) }));
        Assert.That(items[0].SecondSequence, Is.EqualTo("TTTTTTGCATGCAA"));
    }
}
