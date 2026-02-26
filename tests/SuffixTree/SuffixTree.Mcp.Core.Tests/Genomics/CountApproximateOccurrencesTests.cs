using NUnit.Framework;
using Seqeron.Genomics.Alignment;
using SuffixTree.Mcp.Core.Tools;

namespace SuffixTree.Mcp.Core.Tests;

[TestFixture]
[Category("McpCore")]
public class CountApproximateOccurrencesTests
{
    [Test]
    public void CountApproximateOccurrences_InvalidArguments_ThrowArgumentException()
    {
        Assert.Throws<ArgumentException>(() => SuffixTreeTools.CountApproximateOccurrences("", "ATGC", 1));
        Assert.Throws<ArgumentException>(() => SuffixTreeTools.CountApproximateOccurrences(null!, "ATGC", 1));
        Assert.Throws<ArgumentException>(() => SuffixTreeTools.CountApproximateOccurrences("ATGC", "", 1));
        Assert.Throws<ArgumentException>(() => SuffixTreeTools.CountApproximateOccurrences("ATGC", "AT", -1));
    }

    [TestCase("ATGATGATG", "ATG", 0, 3)]
    [TestCase("ATGCTG", "ATG", 1, 2)]
    [TestCase("AAAA", "GGG", 0, 0)]
    [TestCase("AACAA", "AAA", 1, 3)]
    public void CountApproximateOccurrences_ReturnsExpectedCount(
        string sequence,
        string pattern,
        int maxMismatches,
        int expectedCount)
    {
        var result = SuffixTreeTools.CountApproximateOccurrences(sequence, pattern, maxMismatches);
        Assert.That(result.Count, Is.EqualTo(expectedCount));
    }

    [Test]
    public void CountApproximateOccurrences_MatchesCoreAlgorithm()
    {
        var testCases = new (string Sequence, string Pattern, int MaxMismatches)[]
        {
            ("ATGATGATG", "ATG", 0),
            ("ATGCTG", "ATG", 1),
            ("AACAA", "AAA", 1),
            ("ACGTACGT", "CGT", 0),
            ("ACGTACGT", "CGA", 1)
        };

        foreach (var (sequence, pattern, maxMismatches) in testCases)
        {
            int expected = ApproximateMatcher.CountApproximateOccurrences(sequence, pattern, maxMismatches);
            int actual = SuffixTreeTools.CountApproximateOccurrences(sequence, pattern, maxMismatches).Count;
            Assert.That(actual, Is.EqualTo(expected),
                $"sequence={sequence}, pattern={pattern}, k={maxMismatches}");
        }
    }
}

