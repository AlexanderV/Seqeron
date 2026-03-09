using NUnit.Framework;
using Seqeron.Genomics.Analysis;
using Seqeron.Genomics.Core;
using SuffixTree.Mcp.Core.Tools;

namespace SuffixTree.Mcp.Core.Tests;

[TestFixture]
[Category("McpCore")]
public class CalculateSimilarityTests
{
    [Test]
    public void CalculateSimilarity_InvalidArguments_ThrowArgumentException()
    {
        Assert.Throws<ArgumentException>(() => SuffixTreeGenomicsTools.CalculateSimilarity("", "ATGC"));
        Assert.Throws<ArgumentException>(() => SuffixTreeGenomicsTools.CalculateSimilarity(null!, "ATGC"));
        Assert.Throws<ArgumentException>(() => SuffixTreeGenomicsTools.CalculateSimilarity("ATGC", ""));
        Assert.Throws<ArgumentException>(() => SuffixTreeGenomicsTools.CalculateSimilarity("ATGN", "ATGC")); // invalid DNA symbol
    }

    [TestCase("ATAT", "ATAT", 2, 100.0)]
    [TestCase("ATAT", "CGCG", 2, 0.0)]
    [TestCase("ATGC", "ATGA", 2, 50.0)]
    public void CalculateSimilarity_ReturnsExpectedPercentage(
        string sequence1,
        string sequence2,
        int kmerSize,
        double expected)
    {
        var result = SuffixTreeGenomicsTools.CalculateSimilarity(sequence1, sequence2, kmerSize);
        Assert.That(result.Similarity, Is.EqualTo(expected).Within(1e-9));
    }

    [Test]
    public void CalculateSimilarity_MatchesGenomicAnalyzer()
    {
        var cases = new (string S1, string S2, int K)[]
        {
            ("ATGCATGC", "ATGCATGC", 5),
            ("ATGCATGC", "AAAATTTT", 3),
            ("ACGTACGT", "ACGTTCGT", 4),
            ("ATATAT", "TATATA", 2)
        };

        foreach (var (s1, s2, k) in cases)
        {
            double expected = GenomicAnalyzer.CalculateSimilarity(new DnaSequence(s1), new DnaSequence(s2), k);
            double actual = SuffixTreeGenomicsTools.CalculateSimilarity(s1, s2, k).Similarity;
            Assert.That(actual, Is.EqualTo(expected).Within(1e-9),
                $"s1={s1}, s2={s2}, k={k}");
        }
    }
}

