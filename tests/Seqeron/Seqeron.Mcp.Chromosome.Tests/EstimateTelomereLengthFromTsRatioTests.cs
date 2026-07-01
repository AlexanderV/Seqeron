using NUnit.Framework;
using Seqeron.Mcp.Chromosome.Tools;

namespace Seqeron.Mcp.Chromosome.Tests;

/// <summary>
/// Tests for <c>estimate_telomere_length_from_ts_ratio</c>.
/// ChromosomeAnalyzer.EstimateTelomereLengthFromTSRatio = referenceLength * tsRatio / referenceRatio.
/// </summary>
[TestFixture]
public class EstimateTelomereLengthFromTsRatioTests
{
    [Test]
    public void EstimateTelomereLengthFromTsRatio_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => ChromosomeTools.EstimateTelomereLengthFromTsRatio(1.5));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.EstimateTelomereLengthFromTsRatio(1.0, referenceRatio: 0));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.EstimateTelomereLengthFromTsRatio(1.0, referenceRatio: -1));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.EstimateTelomereLengthFromTsRatio(1.0, referenceLength: 0));
    }

    [Test]
    public void EstimateTelomereLengthFromTsRatio_Binding_InvokesSuccessfully()
    {
        // 7000 * 2 / 1 = 14000
        Assert.That(ChromosomeTools.EstimateTelomereLengthFromTsRatio(2.0).TelomereLength, Is.EqualTo(14000.0).Within(1e-6));
        // 7000 * 1 / 1 = 7000 (reference)
        Assert.That(ChromosomeTools.EstimateTelomereLengthFromTsRatio(1.0).TelomereLength, Is.EqualTo(7000.0).Within(1e-6));
        // 10000 * 0.5 / 2 = 2500
        Assert.That(ChromosomeTools.EstimateTelomereLengthFromTsRatio(0.5, referenceRatio: 2.0, referenceLength: 10000).TelomereLength,
            Is.EqualTo(2500.0).Within(1e-6));
    }
}
