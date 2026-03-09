using NUnit.Framework;
using SuffixTree.Mcp.Core.Tools;

namespace SuffixTree.Mcp.Core.Tests;

[TestFixture]
[Category("McpCore")]
public class HammingDistanceTests
{
    [Test]
    public void HammingDistance_InvalidArguments_ThrowArgumentException()
    {
        Assert.Throws<ArgumentException>(() => SuffixTreeGenomicsTools.HammingDistance("", "ATGC"));
        Assert.Throws<ArgumentException>(() => SuffixTreeGenomicsTools.HammingDistance(null!, "ATGC"));
        Assert.Throws<ArgumentException>(() => SuffixTreeGenomicsTools.HammingDistance("ATGC", ""));
        Assert.Throws<ArgumentException>(() => SuffixTreeGenomicsTools.HammingDistance("ATGC", "AT"));
    }

    [TestCase("ATGC", "ATGC", 0)]
    [TestCase("ATGC", "ATGG", 1)]
    [TestCase("AAAA", "TTTT", 4)]
    [TestCase("atgc", "ATGG", 1)] // case-insensitive path in ApproximateMatcher
    public void HammingDistance_ReturnsExpectedDistance(string sequence1, string sequence2, int expected)
    {
        var result = SuffixTreeGenomicsTools.HammingDistance(sequence1, sequence2);
        Assert.That(result.Distance, Is.EqualTo(expected));
    }
}

