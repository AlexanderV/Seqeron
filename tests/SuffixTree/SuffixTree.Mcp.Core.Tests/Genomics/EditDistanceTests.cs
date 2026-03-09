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
        Assert.Throws<ArgumentException>(() => SuffixTreeGenomicsTools.EditDistance("", "ATGC"));
        Assert.Throws<ArgumentException>(() => SuffixTreeGenomicsTools.EditDistance(null!, "ATGC"));
        Assert.Throws<ArgumentException>(() => SuffixTreeGenomicsTools.EditDistance("ATGC", ""));
    }

    [TestCase("ATGC", "ATGC", 0)]
    [TestCase("ATGC", "ATGG", 1)]
    [TestCase("ATGC", "ATG", 1)]
    [TestCase("kitten", "sitting", 3)]
    [TestCase("atgc", "ATGG", 4)] // ApproximateMatcher.EditDistance is case-sensitive
    public void EditDistance_ReturnsExpectedDistance(string sequence1, string sequence2, int expected)
    {
        var result = SuffixTreeGenomicsTools.EditDistance(sequence1, sequence2);
        Assert.That(result.Distance, Is.EqualTo(expected));
    }
}

