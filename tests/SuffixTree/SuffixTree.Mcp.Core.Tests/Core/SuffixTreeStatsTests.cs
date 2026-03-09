using NUnit.Framework;
using SuffixTree.Mcp.Core.Tools;

namespace SuffixTree.Mcp.Core.Tests;

[TestFixture]
[Category("McpCore")]
public class SuffixTreeStatsTests
{
    [Test]
    public void SuffixTreeStats_InvalidArguments_ThrowArgumentException()
    {
        Assert.Throws<ArgumentException>(() => SuffixTreeCoreTools.SuffixTreeStats(""));
        Assert.Throws<ArgumentException>(() => SuffixTreeCoreTools.SuffixTreeStats(null!));
    }

    [TestCase("banana")]
    [TestCase("abcdef")]
    [TestCase("aaaaa")]
    [TestCase("mississippi")]
    public void SuffixTreeStats_MatchesFundamentalInvariants(string text)
    {
        var stats = SuffixTreeCoreTools.SuffixTreeStats(text);

        Assert.That(stats.TextLength, Is.EqualTo(text.Length));
        Assert.That(stats.LeafCount, Is.EqualTo(text.Length));
        Assert.That(stats.MaxDepth, Is.EqualTo(text.Length));
        Assert.That(stats.NodeCount, Is.GreaterThanOrEqualTo(text.Length + 1));
        Assert.That(stats.NodeCount, Is.LessThanOrEqualTo(2 * text.Length + 1));
    }

    [Test]
    public void SuffixTreeStats_IsConsistentWithOtherTools()
    {
        const string text = "banana";
        var stats = SuffixTreeCoreTools.SuffixTreeStats(text);
        var lrs = SuffixTreeCoreTools.SuffixTreeLrs(text);

        Assert.That(stats.MaxDepth, Is.GreaterThanOrEqualTo(lrs.Length));
        Assert.That(SuffixTreeCoreTools.SuffixTreeCount(text, "").Count, Is.EqualTo(stats.LeafCount));
    }
}

