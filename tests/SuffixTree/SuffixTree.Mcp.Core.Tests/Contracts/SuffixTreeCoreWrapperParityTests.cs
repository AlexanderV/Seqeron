using System.Linq;
using NUnit.Framework;
using SuffixTree.Mcp.Core.Tools;

namespace SuffixTree.Mcp.Core.Tests;

[TestFixture]
[Category("McpCore")]
public class SuffixTreeCoreWrapperParityTests
{
    [Test]
    public void CompatibilityWrappers_MatchCoreToolImplementations()
    {
        var containsWrapper = SuffixTreeTools.SuffixTreeContains("banana", "ana");
        var containsCore = SuffixTreeCoreTools.SuffixTreeContains("banana", "ana");
        Assert.That(containsWrapper.Found, Is.EqualTo(containsCore.Found));

        var countWrapper = SuffixTreeTools.SuffixTreeCount("banana", "ana");
        var countCore = SuffixTreeCoreTools.SuffixTreeCount("banana", "ana");
        Assert.That(countWrapper.Count, Is.EqualTo(countCore.Count));

        var findAllWrapper = SuffixTreeTools.SuffixTreeFindAll("banana", "ana");
        var findAllCore = SuffixTreeCoreTools.SuffixTreeFindAll("banana", "ana");
        Assert.That(findAllWrapper.Positions.OrderBy(x => x).ToArray(),
            Is.EqualTo(findAllCore.Positions.OrderBy(x => x).ToArray()));

        var lrsWrapper = SuffixTreeTools.SuffixTreeLrs("banana");
        var lrsCore = SuffixTreeCoreTools.SuffixTreeLrs("banana");
        Assert.That((lrsWrapper.Substring, lrsWrapper.Length), Is.EqualTo((lrsCore.Substring, lrsCore.Length)));

        var lcsWrapper = SuffixTreeTools.SuffixTreeLcs("banana", "panama");
        var lcsCore = SuffixTreeCoreTools.SuffixTreeLcs("banana", "panama");
        Assert.That((lcsWrapper.Substring, lcsWrapper.Length), Is.EqualTo((lcsCore.Substring, lcsCore.Length)));

        var statsWrapper = SuffixTreeTools.SuffixTreeStats("banana");
        var statsCore = SuffixTreeCoreTools.SuffixTreeStats("banana");
        Assert.That((statsWrapper.NodeCount, statsWrapper.LeafCount, statsWrapper.MaxDepth, statsWrapper.TextLength),
            Is.EqualTo((statsCore.NodeCount, statsCore.LeafCount, statsCore.MaxDepth, statsCore.TextLength)));
    }

    [Test]
    public void CompatibilityWrappers_AndCoreTools_ThrowEquivalentArgumentErrors()
    {
        AssertBothThrowArgumentException(
            () => SuffixTreeTools.SuffixTreeContains("", "a"),
            () => SuffixTreeCoreTools.SuffixTreeContains("", "a"));
        AssertBothThrowArgumentException(
            () => SuffixTreeTools.SuffixTreeContains("abc", null!),
            () => SuffixTreeCoreTools.SuffixTreeContains("abc", null!));

        AssertBothThrowArgumentException(
            () => SuffixTreeTools.SuffixTreeCount("", "a"),
            () => SuffixTreeCoreTools.SuffixTreeCount("", "a"));
        AssertBothThrowArgumentException(
            () => SuffixTreeTools.SuffixTreeFindAll("", "a"),
            () => SuffixTreeCoreTools.SuffixTreeFindAll("", "a"));
        AssertBothThrowArgumentException(
            () => SuffixTreeTools.SuffixTreeLrs(""),
            () => SuffixTreeCoreTools.SuffixTreeLrs(""));
        AssertBothThrowArgumentException(
            () => SuffixTreeTools.SuffixTreeLcs("", "abc"),
            () => SuffixTreeCoreTools.SuffixTreeLcs("", "abc"));
        AssertBothThrowArgumentException(
            () => SuffixTreeTools.SuffixTreeStats(""),
            () => SuffixTreeCoreTools.SuffixTreeStats(""));
    }

    private static void AssertBothThrowArgumentException(TestDelegate wrapperCall, TestDelegate coreCall)
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentException>(wrapperCall);
            Assert.Throws<ArgumentException>(coreCall);
        });
    }
}
