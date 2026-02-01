using NUnit.Framework;
using SuffixTree.Mcp.Core.Tools;

namespace SuffixTree.Mcp.Core.Tests;

[TestFixture]
public class FindLongestCommonRegionTests
{
    [Test]
    public void FindLongestCommonRegion_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => SuffixTreeTools.FindLongestCommonRegion("ATGCATGC", "CATGCAT"));
        Assert.Throws<ArgumentException>(() => SuffixTreeTools.FindLongestCommonRegion("", "ATGC"));
        Assert.Throws<ArgumentException>(() => SuffixTreeTools.FindLongestCommonRegion(null!, "ATGC"));
        Assert.Throws<ArgumentException>(() => SuffixTreeTools.FindLongestCommonRegion("ATGC", ""));
    }

    [Test]
    public void FindLongestCommonRegion_Binding_InvokesSuccessfully()
    {
        var result = SuffixTreeTools.FindLongestCommonRegion("ATGCATGC", "CATGCAT");
        Assert.That(result, Is.Not.Null);
        Assert.That(result.Length, Is.GreaterThan(0));
        Assert.That(result.Position1, Is.GreaterThanOrEqualTo(0));
        Assert.That(result.Position2, Is.GreaterThanOrEqualTo(0));

        var noCommon = SuffixTreeTools.FindLongestCommonRegion("AAAA", "TTTT");
        Assert.That(noCommon.Length, Is.EqualTo(0));
    }
}
