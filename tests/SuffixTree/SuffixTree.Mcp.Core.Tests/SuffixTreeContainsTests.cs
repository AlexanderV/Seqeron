using NUnit.Framework;
using SuffixTree.Mcp.Core.Tools;

namespace SuffixTree.Mcp.Core.Tests;

[TestFixture]
[Category("McpCore")]
public class SuffixTreeContainsTests
{
    [Test]
    public void SuffixTreeContains_InvalidArguments_ThrowArgumentException()
    {
        Assert.Throws<ArgumentException>(() => SuffixTreeTools.SuffixTreeContains("", "pattern"));
        Assert.Throws<ArgumentException>(() => SuffixTreeTools.SuffixTreeContains(null!, "pattern"));
        Assert.Throws<ArgumentException>(() => SuffixTreeTools.SuffixTreeContains("text", null!));
    }

    [TestCase("banana", "ana", true)]
    [TestCase("banana", "xyz", false)]
    [TestCase("banana", "", true)]
    [TestCase("aaaaa", "aa", true)]
    [TestCase("AbCd", "abcd", false)]
    [TestCase("AbCd", "AbC", true)]
    public void SuffixTreeContains_ReturnsExpectedResult(string text, string pattern, bool expected)
    {
        var result = SuffixTreeTools.SuffixTreeContains(text, pattern);
        Assert.That(result.Found, Is.EqualTo(expected));
    }

    [Test]
    public void SuffixTreeContains_IsConsistentWithCountAndFindAllTools()
    {
        const string text = "banana";
        string[] patterns = { "", "a", "ana", "n", "xyz", "banana", "bananax" };

        foreach (string pattern in patterns)
        {
            bool found = SuffixTreeTools.SuffixTreeContains(text, pattern).Found;
            int count = SuffixTreeTools.SuffixTreeCount(text, pattern).Count;
            int[] positions = SuffixTreeTools.SuffixTreeFindAll(text, pattern).Positions;

            Assert.That(found, Is.EqualTo(count > 0));
            Assert.That(found, Is.EqualTo(positions.Length > 0));
        }
    }
}

