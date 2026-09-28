using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>compression_ratio</c> MCP tool.
/// Expected values taken from the algorithm's own unit tests
/// (SequenceComplexity_EstimateCompressionRatio_Tests, normalized Lempel-Ziv),
/// NOT the wrapper's output.
/// </summary>
[TestFixture]
public class CompressionRatioTests
{
    [Test]
    public void CompressionRatio_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.CompressionRatio("ACGTACGTACGTACGT"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.CompressionRatio(""));
        Assert.Throws<ArgumentException>(() => AnalysisTools.CompressionRatio(null!));
    }

    [Test]
    public void CompressionRatio_Binding_InvokesSuccessfully()
    {
        Assert.Multiple(() =>
        {
            // Classic LZ76 string (antropy doctest): c = 6 -> normalized LZ = 6/(16/log2 16) = 1.5.
            Assert.That(AnalysisTools.CompressionRatio("1001111011000010").Ratio,
                Is.EqualTo(1.5).Within(1e-10));
            // ACGT x4 (repetitive): A/C/G/T/ACGTACGTACGT, c = 5 -> 5/(16/log4 16) = 0.625 (antropy).
            Assert.That(AnalysisTools.CompressionRatio("ACGTACGTACGTACGT").Ratio,
                Is.EqualTo(0.625).Within(1e-10));
        });
    }
}
