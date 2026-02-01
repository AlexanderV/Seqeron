using NUnit.Framework;
using SuffixTree.Mcp.Core.Tools;

namespace SuffixTree.Mcp.Core.Tests;

[TestFixture]
public class SuffixTreeCountTests
{
    /// <summary>
    /// Schema test: validates input parameters are handled correctly.
    /// </summary>
    [Test]
    public void SuffixTreeCount_Schema_ValidatesCorrectly()
    {
        // Valid inputs should not throw
        Assert.DoesNotThrow(() => SuffixTreeTools.SuffixTreeCount("banana", "ana"));

        // Empty text should throw
        Assert.Throws<ArgumentException>(() => SuffixTreeTools.SuffixTreeCount("", "pattern"));

        // Null text should throw
        Assert.Throws<ArgumentException>(() => SuffixTreeTools.SuffixTreeCount(null!, "pattern"));

        // Null pattern should throw
        Assert.Throws<ArgumentException>(() => SuffixTreeTools.SuffixTreeCount("text", null!));
    }

    /// <summary>
    /// Binding test: validates the tool invokes successfully and returns correct structure.
    /// </summary>
    [Test]
    public void SuffixTreeCount_Binding_InvokesSuccessfully()
    {
        // Invoke the tool - "ana" appears twice in "banana" (positions 1 and 3)
        var result = SuffixTreeTools.SuffixTreeCount("banana", "ana");

        // Verify result is not null and has expected structure
        Assert.That(result, Is.Not.Null);
        Assert.That(result.Count, Is.EqualTo(2));

        // Verify zero count case
        var notFound = SuffixTreeTools.SuffixTreeCount("banana", "xyz");
        Assert.That(notFound, Is.Not.Null);
        Assert.That(notFound.Count, Is.EqualTo(0));

        // Verify single occurrence
        var single = SuffixTreeTools.SuffixTreeCount("banana", "b");
        Assert.That(single.Count, Is.EqualTo(1));
    }
}
