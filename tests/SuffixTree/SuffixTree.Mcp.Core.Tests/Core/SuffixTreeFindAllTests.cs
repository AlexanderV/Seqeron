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
        Assert.Throws<ArgumentException>(() => SuffixTreeCoreTools.SuffixTreeFindAll("", "pattern"));
        Assert.Throws<ArgumentException>(() => SuffixTreeCoreTools.SuffixTreeFindAll(null!, "pattern"));
        Assert.Throws<ArgumentException>(() => SuffixTreeCoreTools.SuffixTreeFindAll("text", null!));
    }

    [Test]
    public void SuffixTreeFindAll_ReturnsExactPositions()
    {
        Assert.That(SuffixTreeCoreTools.SuffixTreeFindAll("banana", "ana").Positions.OrderBy(x => x).ToArray(),
            Is.EqualTo(new[] { 1, 3 }));
        Assert.That(SuffixTreeCoreTools.SuffixTreeFindAll("aaaaa", "aa").Positions.OrderBy(x => x).ToArray(),
            Is.EqualTo(new[] { 0, 1, 2, 3 }));
        Assert.That(SuffixTreeCoreTools.SuffixTreeFindAll("abc", "").Positions.OrderBy(x => x).ToArray(),
            Is.EqualTo(new[] { 0, 1, 2 }));
        Assert.That(SuffixTreeCoreTools.SuffixTreeFindAll("abc", "xyz").Positions, Is.Empty);
    }
}

