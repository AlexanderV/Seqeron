using NUnit.Framework;
using SuffixTree.Mcp.Core.Tools;

namespace SuffixTree.Mcp.Core.Tests;

[TestFixture]
public class SuffixTreeLrsTests
{
    [Test]
    public void SuffixTreeLrs_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => SuffixTreeTools.SuffixTreeLrs("banana"));
        Assert.Throws<ArgumentException>(() => SuffixTreeTools.SuffixTreeLrs(""));
        Assert.Throws<ArgumentException>(() => SuffixTreeTools.SuffixTreeLrs(null!));
    }

    [Test]
    public void SuffixTreeLrs_Binding_InvokesSuccessfully()
    {
        var result = SuffixTreeTools.SuffixTreeLrs("banana");
        Assert.That(result, Is.Not.Null);
        Assert.That(result.Substring, Is.EqualTo("ana"));
        Assert.That(result.Length, Is.EqualTo(3));

        var noRepeat = SuffixTreeTools.SuffixTreeLrs("abcdef");
        Assert.That(noRepeat.Substring, Is.Empty);
        Assert.That(noRepeat.Length, Is.EqualTo(0));
    }
}
