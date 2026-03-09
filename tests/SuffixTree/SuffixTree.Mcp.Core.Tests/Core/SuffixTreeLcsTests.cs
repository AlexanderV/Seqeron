using NUnit.Framework;
using SuffixTree.Mcp.Core.Tools;

namespace SuffixTree.Mcp.Core.Tests;

[TestFixture]
[Category("McpCore")]
public class SuffixTreeLcsTests
{
    [Test]
    public void SuffixTreeLcs_InvalidArguments_ThrowArgumentException()
    {
        Assert.Throws<ArgumentException>(() => SuffixTreeCoreTools.SuffixTreeLcs("", "text"));
        Assert.Throws<ArgumentException>(() => SuffixTreeCoreTools.SuffixTreeLcs(null!, "text"));
        Assert.Throws<ArgumentException>(() => SuffixTreeCoreTools.SuffixTreeLcs("text", ""));
        Assert.Throws<ArgumentException>(() => SuffixTreeCoreTools.SuffixTreeLcs("text", null!));
    }

    [Test]
    public void SuffixTreeLcs_ReturnsExpectedForRepresentativeInputs()
    {
        var banana = SuffixTreeCoreTools.SuffixTreeLcs("banana", "panama");
        Assert.That(banana.Substring, Is.EqualTo("ana"));
        Assert.That(banana.Length, Is.EqualTo(3));

        var cde = SuffixTreeCoreTools.SuffixTreeLcs("abcdef", "xyzcdeww");
        Assert.That(cde.Substring, Is.EqualTo("cde"));
        Assert.That(cde.Length, Is.EqualTo(3));

        var none = SuffixTreeCoreTools.SuffixTreeLcs("abc", "xyz");
        Assert.That(none.Substring, Is.Empty);
        Assert.That(none.Length, Is.EqualTo(0));
    }

    [Test]
    public void SuffixTreeLcs_ResultSatisfiesCoreContract()
    {
        const string text1 = "abcabxabcd";
        const string text2 = "zzabcabczz";

        var result = SuffixTreeCoreTools.SuffixTreeLcs(text1, text2);
        Assert.That(result.Length, Is.EqualTo(result.Substring.Length));

        if (result.Length > 0)
        {
            Assert.That(text1.Contains(result.Substring, StringComparison.Ordinal), Is.True);
            Assert.That(text2.Contains(result.Substring, StringComparison.Ordinal), Is.True);
        }

        int expectedLength = LongestCommonSubstringLength(text1, text2);
        Assert.That(result.Length, Is.EqualTo(expectedLength));
    }

    private static int LongestCommonSubstringLength(string a, string b)
    {
        if (a.Length == 0 || b.Length == 0)
            return 0;

        var dp = new int[b.Length + 1];
        int best = 0;

        for (int i = 1; i <= a.Length; i++)
        {
            for (int j = b.Length; j >= 1; j--)
            {
                if (a[i - 1] == b[j - 1])
                {
                    dp[j] = dp[j - 1] + 1;
                    if (dp[j] > best)
                        best = dp[j];
                }
                else
                {
                    dp[j] = 0;
                }
            }
        }

        return best;
    }
}

