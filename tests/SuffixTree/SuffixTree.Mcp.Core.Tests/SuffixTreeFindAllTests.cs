using NUnit.Framework;
using SuffixTree.Mcp.Core.Tools;

namespace SuffixTree.Mcp.Core.Tests;

[TestFixture]
public class SuffixTreeFindAllTests
{
    [Test]
    public void SuffixTreeFindAll_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => SuffixTreeTools.SuffixTreeFindAll("banana", "ana"));
        Assert.Throws<ArgumentException>(() => SuffixTreeTools.SuffixTreeFindAll("", "pattern"));
        Assert.Throws<ArgumentException>(() => SuffixTreeTools.SuffixTreeFindAll(null!, "pattern"));
        Assert.Throws<ArgumentException>(() => SuffixTreeTools.SuffixTreeFindAll("text", null!));
    }

    [Test]
    public void SuffixTreeFindAll_Binding_InvokesSuccessfully()
    {
        var result = SuffixTreeTools.SuffixTreeFindAll("banana", "ana");
        Assert.That(result, Is.Not.Null);
        Assert.That(result.Positions, Is.Not.Null);
        Assert.That(result.Positions, Has.Length.EqualTo(2));
        Assert.That(result.Positions, Does.Contain(1));
        Assert.That(result.Positions, Does.Contain(3));

        var notFound = SuffixTreeTools.SuffixTreeFindAll("banana", "xyz");
        Assert.That(notFound.Positions, Is.Empty);
    }
}
