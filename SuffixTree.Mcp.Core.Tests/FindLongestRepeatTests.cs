using NUnit.Framework;
using SuffixTree.Mcp.Core.Tools;

namespace SuffixTree.Mcp.Core.Tests;

[TestFixture]
public class FindLongestRepeatTests
{
    [Test]
    public void FindLongestRepeat_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => SuffixTreeTools.FindLongestRepeat("ATGATGATG"));
        Assert.Throws<ArgumentException>(() => SuffixTreeTools.FindLongestRepeat(""));
        Assert.Throws<ArgumentException>(() => SuffixTreeTools.FindLongestRepeat(null!));
    }

    [Test]
    public void FindLongestRepeat_Binding_InvokesSuccessfully()
    {
        var result = SuffixTreeTools.FindLongestRepeat("ATGATGATG");
        Assert.That(result, Is.Not.Null);
        Assert.That(result.Length, Is.GreaterThan(0));
        Assert.That(result.Positions, Is.Not.Empty);

        var noRepeat = SuffixTreeTools.FindLongestRepeat("ACGT");
        Assert.That(noRepeat, Is.Not.Null);
        Assert.That(noRepeat.Length, Is.EqualTo(0));
    }
}
