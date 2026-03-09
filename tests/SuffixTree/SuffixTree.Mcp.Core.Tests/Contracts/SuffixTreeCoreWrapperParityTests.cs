using System.Linq;
using NUnit.Framework;
using SuffixTree.Mcp.Core.Tools;

namespace SuffixTree.Mcp.Core.Tests;

[TestFixture]
[Category("McpCore")]
public class SuffixTreeCoreWrapperParityTests
{
    [Test]
    public void CoreTools_ProduceConsistentResults()
    {
        var contains = SuffixTreeCoreTools.SuffixTreeContains("banana", "ana");
        Assert.That(contains.Found, Is.True);

        var count = SuffixTreeCoreTools.SuffixTreeCount("banana", "ana");
        Assert.That(count.Count, Is.EqualTo(2));

        var findAll = SuffixTreeCoreTools.SuffixTreeFindAll("banana", "ana");
        Assert.That(findAll.Positions.OrderBy(x => x).ToArray(), Is.EqualTo(new[] { 1, 3 }));

        var lrs = SuffixTreeCoreTools.SuffixTreeLrs("banana");
        Assert.That(lrs.Substring, Is.EqualTo("ana"));
        Assert.That(lrs.Length, Is.EqualTo(3));

        var lcs = SuffixTreeCoreTools.SuffixTreeLcs("banana", "panama");
        Assert.That(lcs.Substring, Is.EqualTo("ana"));
        Assert.That(lcs.Length, Is.EqualTo(3));

        var stats = SuffixTreeCoreTools.SuffixTreeStats("banana");
        Assert.That(stats.TextLength, Is.EqualTo(6));
        Assert.That(stats.LeafCount, Is.GreaterThan(0));
    }

    [Test]
    public void CoreTools_ThrowArgumentErrors()
    {
        Assert.Throws<ArgumentException>(() => SuffixTreeCoreTools.SuffixTreeContains("", "a"));
        Assert.Throws<ArgumentException>(() => SuffixTreeCoreTools.SuffixTreeContains("abc", null!));
        Assert.Throws<ArgumentException>(() => SuffixTreeCoreTools.SuffixTreeCount("", "a"));
        Assert.Throws<ArgumentException>(() => SuffixTreeCoreTools.SuffixTreeFindAll("", "a"));
        Assert.Throws<ArgumentException>(() => SuffixTreeCoreTools.SuffixTreeLrs(""));
        Assert.Throws<ArgumentException>(() => SuffixTreeCoreTools.SuffixTreeLcs("", "abc"));
        Assert.Throws<ArgumentException>(() => SuffixTreeCoreTools.SuffixTreeStats(""));
    }
}
