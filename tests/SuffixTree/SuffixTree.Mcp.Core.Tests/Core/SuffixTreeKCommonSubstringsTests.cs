using NUnit.Framework;
using SuffixTree.Mcp.Core.Tools;

namespace SuffixTree.Mcp.Core.Tests;

/// <summary>
/// suffix_tree_k_common_substrings — delegation to <c>SuffixTree.FindLongestCommonSubstrings</c> /
/// <c>LongestCommonSubstringLengthsBySupport</c>. Rosalind LCSM sample (GATTACA, TAGACCA, ATACA → AC/CA/TA)
/// and brute-force values locked in SuffixTree.Tests KCommonSubstringTests.
/// </summary>
[TestFixture]
[Category("McpCore")]
public class SuffixTreeKCommonSubstringsTests
{
    private static readonly string[] Lcsm = { "GATTACA", "TAGACCA", "ATACA" };

    [Test]
    public void SuffixTreeKCommonSubstrings_InvalidArguments_Throw()
    {
        Assert.Throws<ArgumentException>(() => SuffixTreeCoreTools.SuffixTreeKCommonSubstrings(null!));
        Assert.Throws<ArgumentException>(() => SuffixTreeCoreTools.SuffixTreeKCommonSubstrings(Array.Empty<string>()));
        Assert.Throws<ArgumentException>(() => SuffixTreeCoreTools.SuffixTreeKCommonSubstrings(new[] { "A", null! }));
        Assert.Throws<ArgumentOutOfRangeException>(() => SuffixTreeCoreTools.SuffixTreeKCommonSubstrings(Lcsm, 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => SuffixTreeCoreTools.SuffixTreeKCommonSubstrings(Lcsm, 4));
    }

    [Test]
    public void SuffixTreeKCommonSubstrings_RosalindLcsmSample()
    {
        var r = SuffixTreeCoreTools.SuffixTreeKCommonSubstrings(Lcsm);
        Assert.Multiple(() =>
        {
            Assert.That(r.Substrings, Is.EqualTo(new[] { "AC", "CA", "TA" }));
            Assert.That(r.Length, Is.EqualTo(2));
            Assert.That(r.MinSupport, Is.EqualTo(3));
            Assert.That(r.LengthsBySupport, Is.EqualTo(new[] { 7, 4, 2 }));
        });
    }

    [Test]
    public void SuffixTreeKCommonSubstrings_MinSupportTwo()
    {
        var r = SuffixTreeCoreTools.SuffixTreeKCommonSubstrings(Lcsm, 2);
        Assert.That(r.Substrings, Is.EqualTo(new[] { "TACA" }));
        Assert.That(r.Length, Is.EqualTo(4));
        Assert.That(r.Substrings, Is.EqualTo(global::SuffixTree.SuffixTree.FindLongestCommonSubstrings(Lcsm, 2)));
    }

    [Test]
    public void SuffixTreeKCommonSubstrings_NothingShared_Empty()
    {
        var r = SuffixTreeCoreTools.SuffixTreeKCommonSubstrings(new[] { "AAA", "CCC" });
        Assert.That(r.Substrings, Is.Empty);
        Assert.That(r.Length, Is.Zero);
    }
}
