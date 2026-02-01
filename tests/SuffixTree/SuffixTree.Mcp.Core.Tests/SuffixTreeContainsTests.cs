using NUnit.Framework;
using SuffixTree.Mcp.Core.Tools;

namespace SuffixTree.Mcp.Core.Tests;

[TestFixture]
public class SuffixTreeContainsTests
{
    /// <summary>
    /// Schema test: validates input parameters are handled correctly.
    /// </summary>
    [Test]
    public void SuffixTreeContains_Schema_ValidatesCorrectly()
    {
        // Valid inputs should not throw
        Assert.DoesNotThrow(() => SuffixTreeTools.SuffixTreeContains("banana", "ana"));

        // Empty text should throw
        Assert.Throws<ArgumentException>(() => SuffixTreeTools.SuffixTreeContains("", "pattern"));

        // Null text should throw
        Assert.Throws<ArgumentException>(() => SuffixTreeTools.SuffixTreeContains(null!, "pattern"));

        // Null pattern should throw
        Assert.Throws<ArgumentException>(() => SuffixTreeTools.SuffixTreeContains("text", null!));
    }

    /// <summary>
    /// Binding test: validates the tool invokes successfully and returns correct structure.
    /// </summary>
    [Test]
    public void SuffixTreeContains_Binding_InvokesSuccessfully()
    {
        // Invoke the tool
        var result = SuffixTreeTools.SuffixTreeContains("banana", "ana");

        // Verify result is not null and has expected structure
        Assert.That(result, Is.Not.Null);
        Assert.That(result.Found, Is.True);

        // Verify negative case
        var notFound = SuffixTreeTools.SuffixTreeContains("banana", "xyz");
        Assert.That(notFound, Is.Not.Null);
        Assert.That(notFound.Found, Is.False);

        // Verify empty pattern (should return true - empty string is substring of any string)
        var emptyPattern = SuffixTreeTools.SuffixTreeContains("banana", "");
        Assert.That(emptyPattern.Found, Is.True);
    }
}
