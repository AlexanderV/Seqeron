using NUnit.Framework;
using SuffixTree.Mcp.Core.Tools;

namespace SuffixTree.Mcp.Core.Tests;

[TestFixture]
[Category("McpCore")]
public class EditDistanceTests
{
    [Test]
    public void EditDistance_InvalidArguments_ThrowArgumentException()
    {
        Assert.Throws<ArgumentException>(() => SuffixTreeGenomicsTools.EditDistance("", "ATGC"));
        Assert.Throws<ArgumentException>(() => SuffixTreeGenomicsTools.EditDistance(null!, "ATGC"));
        Assert.Throws<ArgumentException>(() => SuffixTreeGenomicsTools.EditDistance("ATGC", ""));
    }

    [TestCase("ATGC", "ATGC", 0)]
    [TestCase("ATGC", "ATGG", 1)]
    [TestCase("ATGC", "ATG", 1)]
    [TestCase("kitten", "sitting", 3)]
    [TestCase("atgc", "ATGG", 4)] // ApproximateMatcher.EditDistance is case-sensitive
    public void EditDistance_ReturnsExpectedDistance(string sequence1, string sequence2, int expected)
    {
        var result = SuffixTreeGenomicsTools.EditDistance(sequence1, sequence2);
        Assert.That(result.Distance, Is.EqualTo(expected));
    }

    // Weighted costs (B05 F27): rapidfuzz 3.14.6 Levenshtein.distance(s1, s2, weights=(ins, del, sub)).
    [TestCase("kitten", "sitting", 1, 1, 2, 5)]
    [TestCase("kitten", "sitting", 2, 1, 1, 4)]
    [TestCase("ACGT", "AGT", 1, 5, 1, 5)]
    [TestCase("AGT", "ACGT", 1, 5, 1, 1)]
    [TestCase("ACGTACGT", "ACGTTACG", 1000, 2, 999, 1002)]
    public void EditDistance_WeightedCosts_RapidfuzzValues(string s1, string s2, int ins, int del, int sub, int expected)
    {
        Assert.That(SuffixTreeGenomicsTools.EditDistance(s1, s2, ins, del, sub).Distance, Is.EqualTo(expected));
    }

    [Test]
    public void EditDistance_NegativeCost_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => SuffixTreeGenomicsTools.EditDistance("A", "C", insertionCost: -1));
        Assert.Throws<ArgumentOutOfRangeException>(() => SuffixTreeGenomicsTools.EditDistance("A", "C", deletionCost: -1));
        Assert.Throws<ArgumentOutOfRangeException>(() => SuffixTreeGenomicsTools.EditDistance("A", "C", substitutionCost: -1));
    }
}
