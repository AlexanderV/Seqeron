using NUnit.Framework;
using Seqeron.Genomics.Analysis;
using Seqeron.Genomics.Core;
using SuffixTree.Mcp.Core.Tools;

namespace SuffixTree.Mcp.Core.Tests;

[TestFixture]
[Category("McpCore")]
public class FindLongestCommonRegionTests
{
    [Test]
    public void FindLongestCommonRegion_InvalidArguments_ThrowArgumentException()
    {
        Assert.Throws<ArgumentException>(() => SuffixTreeTools.FindLongestCommonRegion("", "ATGC"));
        Assert.Throws<ArgumentException>(() => SuffixTreeTools.FindLongestCommonRegion(null!, "ATGC"));
        Assert.Throws<ArgumentException>(() => SuffixTreeTools.FindLongestCommonRegion("ATGC", ""));
        Assert.Throws<ArgumentException>(() => SuffixTreeTools.FindLongestCommonRegion("ATGC", "ATGN")); // invalid DNA symbol
    }

    [Test]
    public void FindLongestCommonRegion_ReturnsExpectedForRepresentativeInputs()
    {
        var result = SuffixTreeTools.FindLongestCommonRegion("AAAACCC", "GGAAAAT");
        Assert.That(result.Region, Is.EqualTo("AAAA"));
        Assert.That(result.Length, Is.EqualTo(4));
        Assert.That(result.Position1, Is.EqualTo(0));
        Assert.That(result.Position2, Is.EqualTo(2));

        var noCommon = SuffixTreeTools.FindLongestCommonRegion("AAAA", "TTTT");
        Assert.That(noCommon.Region, Is.Empty);
        Assert.That(noCommon.Length, Is.EqualTo(0));
        Assert.That(noCommon.Position1, Is.EqualTo(-1));
        Assert.That(noCommon.Position2, Is.EqualTo(-1));
    }

    [Test]
    public void FindLongestCommonRegion_MatchesGenomicAnalyzer()
    {
        var cases = new (string S1, string S2)[]
        {
            ("ATGCATGC", "CATGCAT"),
            ("AAAACCC", "GGAAAAT"),
            ("ACGT", "ACGT"),
            ("AAAA", "TTTT")
        };

        foreach (var (s1, s2) in cases)
        {
            var expected = GenomicAnalyzer.FindLongestCommonRegion(new DnaSequence(s1), new DnaSequence(s2));
            var actual = SuffixTreeTools.FindLongestCommonRegion(s1, s2);

            Assert.That(actual.Region, Is.EqualTo(expected.Sequence), $"s1={s1}, s2={s2}: region");
            Assert.That(actual.Length, Is.EqualTo(expected.Length), $"s1={s1}, s2={s2}: length");
            Assert.That(actual.Position1, Is.EqualTo(expected.PositionInFirst), $"s1={s1}, s2={s2}: pos1");
            Assert.That(actual.Position2, Is.EqualTo(expected.PositionInSecond), $"s1={s1}, s2={s2}: pos2");

            if (actual.Length > 0)
            {
                Assert.That(s1.Substring(actual.Position1, actual.Length), Is.EqualTo(actual.Region));
                Assert.That(s2.Substring(actual.Position2, actual.Length), Is.EqualTo(actual.Region));
            }
        }
    }
}

