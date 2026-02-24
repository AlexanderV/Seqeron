using System.Linq;
using NUnit.Framework;
using SuffixTree.Mcp.Core.Tools;

namespace SuffixTree.Mcp.Core.Tests;

[TestFixture]
[Category("McpCore")]
public class SuffixTreeFindAllTests
{
    [Test]
    public void SuffixTreeFindAll_InvalidArguments_ThrowArgumentException()
    {
        Assert.Throws<ArgumentException>(() => SuffixTreeTools.SuffixTreeFindAll("", "pattern"));
        Assert.Throws<ArgumentException>(() => SuffixTreeTools.SuffixTreeFindAll(null!, "pattern"));
        Assert.Throws<ArgumentException>(() => SuffixTreeTools.SuffixTreeFindAll("text", null!));
    }

    [Test]
    public void SuffixTreeFindAll_ReturnsExactPositions()
    {
        Assert.That(SuffixTreeTools.SuffixTreeFindAll("banana", "ana").Positions.OrderBy(x => x).ToArray(),
            Is.EqualTo(new[] { 1, 3 }));
        Assert.That(SuffixTreeTools.SuffixTreeFindAll("aaaaa", "aa").Positions.OrderBy(x => x).ToArray(),
            Is.EqualTo(new[] { 0, 1, 2, 3 }));
        Assert.That(SuffixTreeTools.SuffixTreeFindAll("abc", "").Positions.OrderBy(x => x).ToArray(),
            Is.EqualTo(new[] { 0, 1, 2 }));
        Assert.That(SuffixTreeTools.SuffixTreeFindAll("abc", "xyz").Positions, Is.Empty);
    }

    [Test]
    public void SuffixTreeFindAll_MatchesCountTool()
    {
        const string text = "abracadabra";
        string[] patterns = { "", "a", "abra", "cad", "ra", "xyz" };

        foreach (string pattern in patterns)
        {
            int[] positions = SuffixTreeTools.SuffixTreeFindAll(text, pattern).Positions;
            int count = SuffixTreeTools.SuffixTreeCount(text, pattern).Count;
            Assert.That(positions.Length, Is.EqualTo(count), $"pattern={pattern}");
        }
    }
}

