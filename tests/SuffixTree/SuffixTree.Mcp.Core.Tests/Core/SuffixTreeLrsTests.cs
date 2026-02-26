using NUnit.Framework;
using SuffixTree.Mcp.Core.Tools;

namespace SuffixTree.Mcp.Core.Tests;

[TestFixture]
[Category("McpCore")]
public class SuffixTreeLrsTests
{
    [Test]
    public void SuffixTreeLrs_InvalidArguments_ThrowArgumentException()
    {
        Assert.Throws<ArgumentException>(() => SuffixTreeTools.SuffixTreeLrs(""));
        Assert.Throws<ArgumentException>(() => SuffixTreeTools.SuffixTreeLrs(null!));
    }

    [Test]
    public void SuffixTreeLrs_ReturnsExpectedForRepresentativeInputs()
    {
        var banana = SuffixTreeTools.SuffixTreeLrs("banana");
        Assert.That(banana.Substring, Is.EqualTo("ana"));
        Assert.That(banana.Length, Is.EqualTo(3));

        var repeated = SuffixTreeTools.SuffixTreeLrs("aaaaa");
        Assert.That(repeated.Substring, Is.EqualTo("aaaa"));
        Assert.That(repeated.Length, Is.EqualTo(4));

        var none = SuffixTreeTools.SuffixTreeLrs("abcdef");
        Assert.That(none.Substring, Is.Empty);
        Assert.That(none.Length, Is.EqualTo(0));
    }

    [Test]
    public void SuffixTreeLrs_ResultSatisfiesCoreContract()
    {
        const string text = "mississippi";
        var result = SuffixTreeTools.SuffixTreeLrs(text);

        Assert.That(result.Length, Is.EqualTo(result.Substring.Length));
        if (result.Length > 0)
        {
            Assert.That(text.Contains(result.Substring, StringComparison.Ordinal), Is.True);
            Assert.That(CountOccurrences(text, result.Substring), Is.GreaterThanOrEqualTo(2));
        }
    }

    private static int CountOccurrences(string text, string pattern)
    {
        if (pattern.Length == 0)
            return text.Length;

        int count = 0;
        for (int i = 0; i <= text.Length - pattern.Length; i++)
        {
            if (text.AsSpan(i, pattern.Length).SequenceEqual(pattern.AsSpan()))
                count++;
        }

        return count;
    }
}

