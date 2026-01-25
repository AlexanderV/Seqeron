using NUnit.Framework;
using SuffixTree.Mcp.Core.Tools;

namespace SuffixTree.Mcp.Core.Tests;

[TestFixture]
public class SuffixTreeStatsTests
{
    [Test]
    public void SuffixTreeStats_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => SuffixTreeTools.SuffixTreeStats("banana"));
        Assert.Throws<ArgumentException>(() => SuffixTreeTools.SuffixTreeStats(""));
        Assert.Throws<ArgumentException>(() => SuffixTreeTools.SuffixTreeStats(null!));
    }

    [Test]
    public void SuffixTreeStats_Binding_InvokesSuccessfully()
    {
        var result = SuffixTreeTools.SuffixTreeStats("banana");
        Assert.That(result, Is.Not.Null);
        Assert.That(result.TextLength, Is.EqualTo(6));
        Assert.That(result.NodeCount, Is.GreaterThan(0));
        Assert.That(result.LeafCount, Is.GreaterThan(0));
        Assert.That(result.MaxDepth, Is.GreaterThan(0));
    }
}
