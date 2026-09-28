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

    // Tool output is reported in ascending position order (deterministic), as in the
    // Rosalind SUBS problem (https://rosalind.info/problems/subs/): s = GATATATGCATATACTT,
    // t = ATAT -> "2 4 10" (1-based) = {1, 3, 9} 0-based. Cross-checked with Python
    // re.finditer("(?=ATAT)") -> [1, 3, 9] and "mississippi"/"i" -> [1, 4, 7, 10].
    [Test]
    public void SuffixTreeFindAll_ReturnsPositionsInAscendingOrder()
    {
        Assert.That(SuffixTreeCoreTools.SuffixTreeFindAll("GATATATGCATATACTT", "ATAT").Positions,
            Is.EqualTo(new[] { 1, 3, 9 }));
        Assert.That(SuffixTreeCoreTools.SuffixTreeFindAll("mississippi", "i").Positions,
            Is.EqualTo(new[] { 1, 4, 7, 10 }));
        Assert.That(SuffixTreeCoreTools.SuffixTreeFindAll("abracadabra", "a").Positions,
            Is.EqualTo(new[] { 0, 3, 5, 7, 10 }));
    }
}

