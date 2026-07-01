using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>find_direct_repeats</c> MCP tool.
/// Expected results derived from RepeatFinder.FindDirectRepeats: identical substrings appearing
/// twice with spacing = secondPos - firstPos - length. NOT the wrapper's output.
/// </summary>
[TestFixture]
public class FindDirectRepeatsTests
{
    [Test]
    public void FindDirectRepeats_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.FindDirectRepeats("ATGCGGATGC", 4, 4, 1));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindDirectRepeats("", 4, 4, 1));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindDirectRepeats(null!, 4, 4, 1));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.FindDirectRepeats("ATGC", 1, 4, 1)); // minLength < 2
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.FindDirectRepeats("ATGC", 4, 3, 1)); // maxLength < minLength
    }

    [Test]
    public void FindDirectRepeats_Binding_InvokesSuccessfully()
    {
        // "ATGCGGATGC": the only length-4 repeat is "ATGC" at positions 0 and 6.
        // spacing = 6 - 0 - 4 = 2.
        var items = AnalysisTools.FindDirectRepeats("ATGCGGATGC", 4, 4, 1).Items;
        Assert.That(items, Has.Length.EqualTo(1));
        Assert.Multiple(() =>
        {
            Assert.That(items[0].FirstPosition, Is.EqualTo(0));
            Assert.That(items[0].SecondPosition, Is.EqualTo(6));
            Assert.That(items[0].RepeatSequence, Is.EqualTo("ATGC"));
            Assert.That(items[0].Length, Is.EqualTo(4));
            Assert.That(items[0].Spacing, Is.EqualTo(2));
        });
    }
}
