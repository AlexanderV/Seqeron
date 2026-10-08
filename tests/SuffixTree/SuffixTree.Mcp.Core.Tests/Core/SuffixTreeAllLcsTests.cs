using NUnit.Framework;
using SuffixTree.Mcp.Core.Tools;

namespace SuffixTree.Mcp.Core.Tests;

/// <summary>
/// suffix_tree_all_lcs — delegation to <c>ISuffixTreeAnalysis.FindAllDistinctLongestCommonSubstrings</c>.
/// Expected values are the brute-force-locked cases of SuffixTree.Tests AllDistinctLongestCommonSubstringsTests.
/// </summary>
[TestFixture]
[Category("McpCore")]
public class SuffixTreeAllLcsTests
{
    [Test]
    public void SuffixTreeAllLcs_InvalidArguments_Throw()
    {
        Assert.Throws<ArgumentException>(() => SuffixTreeCoreTools.SuffixTreeAllLcs("", "AC"));
        Assert.Throws<ArgumentException>(() => SuffixTreeCoreTools.SuffixTreeAllLcs(null!, "AC"));
        Assert.Throws<ArgumentException>(() => SuffixTreeCoreTools.SuffixTreeAllLcs("AC", ""));
        Assert.Throws<ArgumentException>(() => SuffixTreeCoreTools.SuffixTreeAllLcs("AC", null!));
    }

    [Test]
    public void SuffixTreeAllLcs_Ties_AllReported()
    {
        var r = SuffixTreeCoreTools.SuffixTreeAllLcs("abcxyzabc", "xyzqabcpxyz");
        Assert.Multiple(() =>
        {
            Assert.That(r.Length, Is.EqualTo(3));
            Assert.That(r.Substrings.Select(s => s.Substring), Is.EqualTo(new[] { "xyz", "abc" }));
            Assert.That(r.Substrings.Select(s => s.PositionsInText1), Is.EqualTo(new[] { new[] { 3 }, new[] { 0, 6 } }));
            Assert.That(r.Substrings.Select(s => s.PositionsInText2), Is.EqualTo(new[] { new[] { 0, 8 }, new[] { 4 } }));
        });
    }

    [Test]
    public void SuffixTreeAllLcs_NothingShared_Empty()
    {
        var r = SuffixTreeCoreTools.SuffixTreeAllLcs("AAAA", "CCC");
        Assert.That(r.Substrings, Is.Empty);
        Assert.That(r.Length, Is.Zero);
    }

    [Test]
    public void SuffixTreeAllLcs_EqualsLibrary()
    {
        var r = SuffixTreeCoreTools.SuffixTreeAllLcs("GATTACA", "TAGACCA");
        var library = global::SuffixTree.SuffixTree.Build("GATTACA").FindAllDistinctLongestCommonSubstrings("TAGACCA");
        Assert.That(r.Substrings.Select(s => s.Substring), Is.EqualTo(library.Select(l => l.Substring)));
        Assert.That(r.Substrings.Select(s => s.Substring), Is.EqualTo(new[] { "TA", "GA", "AC", "CA" }));
    }
}
