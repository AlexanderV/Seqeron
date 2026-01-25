using NUnit.Framework;
using SuffixTree.Mcp.Core.Tools;

namespace SuffixTree.Mcp.Core.Tests;

[TestFixture]
public class HammingDistanceTests
{
    [Test]
    public void HammingDistance_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => SuffixTreeTools.HammingDistance("ATGC", "ATGC"));
        Assert.Throws<ArgumentException>(() => SuffixTreeTools.HammingDistance("", "ATGC"));
        Assert.Throws<ArgumentException>(() => SuffixTreeTools.HammingDistance(null!, "ATGC"));
        Assert.Throws<ArgumentException>(() => SuffixTreeTools.HammingDistance("ATGC", ""));
        Assert.Throws<ArgumentException>(() => SuffixTreeTools.HammingDistance("ATGC", "AT")); // Length mismatch
    }

    [Test]
    public void HammingDistance_Binding_InvokesSuccessfully()
    {
        var identical = SuffixTreeTools.HammingDistance("ATGC", "ATGC");
        Assert.That(identical.Distance, Is.EqualTo(0));

        var oneDiff = SuffixTreeTools.HammingDistance("ATGC", "ATGG");
        Assert.That(oneDiff.Distance, Is.EqualTo(1));

        var allDiff = SuffixTreeTools.HammingDistance("AAAA", "TTTT");
        Assert.That(allDiff.Distance, Is.EqualTo(4));
    }
}
