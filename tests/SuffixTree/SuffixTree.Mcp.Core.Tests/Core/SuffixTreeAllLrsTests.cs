using NUnit.Framework;
using SuffixTree.Mcp.Core.Tools;

namespace SuffixTree.Mcp.Core.Tests;

/// <summary>
/// suffix_tree_all_lrs — delegation to <c>ISuffixTree.FindAllLongestRepeatedSubstrings</c>.
/// Expected values are the brute-force-locked cases of SuffixTree.Tests AllLongestRepeatedSubstringsTests.
/// </summary>
[TestFixture]
[Category("McpCore")]
public class SuffixTreeAllLrsTests
{
    [Test]
    public void SuffixTreeAllLrs_InvalidArguments_Throw()
    {
        Assert.Throws<ArgumentException>(() => SuffixTreeCoreTools.SuffixTreeAllLrs(""));
        Assert.Throws<ArgumentException>(() => SuffixTreeCoreTools.SuffixTreeAllLrs(null!));
    }

    [Test]
    public void SuffixTreeAllLrs_Ties_AllReportedInFirstOccurrenceOrder()
    {
        var r = SuffixTreeCoreTools.SuffixTreeAllLrs("abcxbcaxcab");
        Assert.Multiple(() =>
        {
            Assert.That(r.Length, Is.EqualTo(2));
            Assert.That(r.Substrings.Select(s => s.Substring), Is.EqualTo(new[] { "ab", "bc", "ca" }));
            Assert.That(r.Substrings.Select(s => s.Positions), Is.EqualTo(new[] { new[] { 0, 9 }, new[] { 1, 4 }, new[] { 5, 8 } }));
        });
    }

    [TestCase("banana", "ana", new[] { 1, 3 })]
    [TestCase("mississippi", "issi", new[] { 1, 4 })]
    [TestCase("aaaa", "aaa", new[] { 0, 1 })]
    public void SuffixTreeAllLrs_SingleRepeat_MatchesLockedCases(string text, string substring, int[] positions)
    {
        var r = SuffixTreeCoreTools.SuffixTreeAllLrs(text);
        Assert.That(r.Substrings, Has.Length.EqualTo(1));
        Assert.That(r.Substrings[0].Substring, Is.EqualTo(substring));
        Assert.That(r.Substrings[0].Positions, Is.EqualTo(positions));
        Assert.That(r.Length, Is.EqualTo(substring.Length));
    }

    [Test]
    public void SuffixTreeAllLrs_NoRepeat_Empty_AndEqualsLibrary()
    {
        var none = SuffixTreeCoreTools.SuffixTreeAllLrs("abcd");
        Assert.That(none.Substrings, Is.Empty);
        Assert.That(none.Length, Is.EqualTo(0));

        const string text = "xyzzyxxz";
        var library = global::SuffixTree.SuffixTree.Build(text).FindAllLongestRepeatedSubstrings();
        var r = SuffixTreeCoreTools.SuffixTreeAllLrs(text);
        Assert.That(r.Substrings.Select(s => s.Substring), Is.EqualTo(library.Select(l => l.Substring)));
        Assert.That(r.Substrings.Select(s => s.Positions), Is.EqualTo(library.Select(l => l.Positions.ToArray())));
    }
}
