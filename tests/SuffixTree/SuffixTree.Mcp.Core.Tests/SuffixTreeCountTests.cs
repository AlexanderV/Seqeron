using NUnit.Framework;
using SuffixTree.Mcp.Core.Tools;

namespace SuffixTree.Mcp.Core.Tests;

[TestFixture]
[Category("McpCore")]
public class SuffixTreeCountTests
{
    [Test]
    public void SuffixTreeCount_InvalidArguments_ThrowArgumentException()
    {
        Assert.Throws<ArgumentException>(() => SuffixTreeTools.SuffixTreeCount("", "pattern"));
        Assert.Throws<ArgumentException>(() => SuffixTreeTools.SuffixTreeCount(null!, "pattern"));
        Assert.Throws<ArgumentException>(() => SuffixTreeTools.SuffixTreeCount("text", null!));
    }

    [TestCase("banana", "ana", 2)]
    [TestCase("banana", "b", 1)]
    [TestCase("banana", "xyz", 0)]
    [TestCase("aaaaa", "aa", 4)]
    [TestCase("abc", "", 3)]
    [TestCase("AbCd", "abcd", 0)]
    public void SuffixTreeCount_ReturnsExpectedCount(string text, string pattern, int expectedCount)
    {
        var result = SuffixTreeTools.SuffixTreeCount(text, pattern);
        Assert.That(result.Count, Is.EqualTo(expectedCount));
    }

    [Test]
    public void SuffixTreeCount_MatchesFindAllLength()
    {
        const string text = "mississippi";
        string[] patterns = { "", "i", "issi", "ss", "sip", "xyz" };

        foreach (string pattern in patterns)
        {
            int count = SuffixTreeTools.SuffixTreeCount(text, pattern).Count;
            int[] positions = SuffixTreeTools.SuffixTreeFindAll(text, pattern).Positions;
            Assert.That(count, Is.EqualTo(positions.Length), $"pattern={pattern}");
        }
    }
}

