using NUnit.Framework;
using SuffixTree.Mcp.Core.Tools;

namespace SuffixTree.Mcp.Core.Tests;

[TestFixture]
public class SuffixTreeLcsTests
{
    [Test]
    public void SuffixTreeLcs_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => SuffixTreeTools.SuffixTreeLcs("banana", "bandana"));
        Assert.Throws<ArgumentException>(() => SuffixTreeTools.SuffixTreeLcs("", "text"));
        Assert.Throws<ArgumentException>(() => SuffixTreeTools.SuffixTreeLcs(null!, "text"));
        Assert.Throws<ArgumentException>(() => SuffixTreeTools.SuffixTreeLcs("text", ""));
        Assert.Throws<ArgumentException>(() => SuffixTreeTools.SuffixTreeLcs("text", null!));
    }

    [Test]
    public void SuffixTreeLcs_Binding_InvokesSuccessfully()
    {
        var result = SuffixTreeTools.SuffixTreeLcs("banana", "panama");
        Assert.That(result, Is.Not.Null);
        Assert.That(result.Substring, Is.EqualTo("ana"));
        Assert.That(result.Length, Is.EqualTo(3));

        var noCommon = SuffixTreeTools.SuffixTreeLcs("abc", "xyz");
        Assert.That(noCommon.Substring, Is.Empty);
        Assert.That(noCommon.Length, Is.EqualTo(0));
    }
}
