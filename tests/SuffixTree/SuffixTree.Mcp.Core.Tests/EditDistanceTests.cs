using NUnit.Framework;
using SuffixTree.Mcp.Core.Tools;

namespace SuffixTree.Mcp.Core.Tests;

[TestFixture]
[Category("McpCore")]
public class EditDistanceTests
{
    [Test]
    public void EditDistance_InvalidArguments_ThrowArgumentException()
    {
        Assert.Throws<ArgumentException>(() => SuffixTreeTools.EditDistance("", "ATGC"));
        Assert.Throws<ArgumentException>(() => SuffixTreeTools.EditDistance(null!, "ATGC"));
        Assert.Throws<ArgumentException>(() => SuffixTreeTools.EditDistance("ATGC", ""));
    }

    [TestCase("ATGC", "ATGC", 0)]
    [TestCase("ATGC", "ATGG", 1)]
    [TestCase("ATGC", "ATG", 1)]
    [TestCase("kitten", "sitting", 3)]
    [TestCase("atgc", "ATGG", 1)] // case-insensitive in ApproximateMatcher.EditDistance
    public void EditDistance_ReturnsExpectedDistance(string sequence1, string sequence2, int expected)
    {
        var result = SuffixTreeTools.EditDistance(sequence1, sequence2);
        Assert.That(result.Distance, Is.EqualTo(expected));
    }
}

