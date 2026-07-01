using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>at_skew</c> MCP tool.
/// Expected values taken from the algorithm's own unit tests
/// (GcSkewCalculator_CalculateAtSkew_Tests), NOT from the wrapper's output.
/// </summary>
[TestFixture]
public class AtSkewTests
{
    [Test]
    public void AtSkew_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.AtSkew("AAAT"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.AtSkew(""));
        Assert.Throws<ArgumentException>(() => AnalysisTools.AtSkew(null!));
    }

    [Test]
    public void AtSkew_Binding_InvokesSuccessfully()
    {
        Assert.Multiple(() =>
        {
            // (A-T)/(A+T): AAAA -> +1, TTTT -> -1, ATAT -> 0.
            Assert.That(AnalysisTools.AtSkew("AAAA").AtSkew, Is.EqualTo(1.0).Within(1e-10));
            Assert.That(AnalysisTools.AtSkew("TTTT").AtSkew, Is.EqualTo(-1.0).Within(1e-10));
            Assert.That(AnalysisTools.AtSkew("ATAT").AtSkew, Is.EqualTo(0.0).Within(1e-10));
            // AAAT -> (3-1)/4 = 0.5.
            Assert.That(AnalysisTools.AtSkew("AAAT").AtSkew, Is.EqualTo(0.5).Within(1e-10));
        });
    }

    [Test]
    public void AtSkew_IgnoresNonAtBases()
    {
        // Only A/T counted; G/C ignored: A=3, T=1 -> 0.5.
        Assert.That(AnalysisTools.AtSkew("AAATGGGCCC").AtSkew, Is.EqualTo(0.5).Within(1e-10));
        // No A/T at all -> defined as 0.
        Assert.That(AnalysisTools.AtSkew("GGCC").AtSkew, Is.EqualTo(0.0).Within(1e-10));
    }
}
