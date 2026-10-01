using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>lempel_ziv_complexity</c> MCP tool. Expected = antropy 0.2.2 (_lz_complexity /
/// lziv_complexity(normalize=True)) and the Lempel &amp; Ziv 1976 worked example, locked in
/// SequenceComplexity_EstimateCompressionRatio_Tests. NOT the wrapper's output.
/// </summary>
[TestFixture]
public class LempelZivComplexityTests
{
    [Test]
    public void LempelZivComplexity_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.LempelZivComplexity("ACGT"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.LempelZivComplexity(""));
        Assert.Throws<ArgumentException>(() => AnalysisTools.LempelZivComplexity(null!));
    }

    [TestCase("0001101001000101", 6)]
    [TestCase("1001111011000010", 6)]
    [TestCase("0000000000000000", 2)]
    public void LempelZivComplexity_Binding_RawMatchesLz76(string s, int expected)
    {
        Assert.That(AnalysisTools.LempelZivComplexity(s).Complexity, Is.EqualTo(expected));
    }

    [TestCase("HELLO WORLD! HELLO WORLD! HELLO WORLD! HELLO WORLD!", 11, 0.38596001132145313)]
    [TestCase("CAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAG", 4, 0.2484555351907228)]
    public void LempelZivComplexity_RawAndNormalized_MatchAntropy(string s, int raw, double normalized)
    {
        var r = AnalysisTools.LempelZivComplexity(s);
        Assert.That(r.Complexity, Is.EqualTo(raw));
        Assert.That(r.Normalized, Is.EqualTo(normalized).Within(1e-12));
        Assert.That(r.Normalized, Is.EqualTo(AnalysisTools.CompressionRatio(s).Ratio), "compression_ratio = normalized LZ76");
    }
}
