using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>gc_skew</c> MCP tool.
/// Expected values from the definition GC-skew = (G - C) / (G + C), NOT the wrapper output.
/// </summary>
[TestFixture]
public class GcSkewTests
{
    [Test]
    public void GcSkew_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.GcSkew("GGGG"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.GcSkew(""));
        Assert.Throws<ArgumentException>(() => AnalysisTools.GcSkew(null!));
    }

    [Test]
    public void GcSkew_Binding_InvokesSuccessfully()
    {
        Assert.Multiple(() =>
        {
            // (4-0)/4 = 1.
            Assert.That(AnalysisTools.GcSkew("GGGG").GcSkew, Is.EqualTo(1.0).Within(1e-12));
            // (0-4)/4 = -1.
            Assert.That(AnalysisTools.GcSkew("CCCC").GcSkew, Is.EqualTo(-1.0).Within(1e-12));
            // (4-4)/8 = 0.
            Assert.That(AnalysisTools.GcSkew("GGGGCCCC").GcSkew, Is.EqualTo(0.0).Within(1e-12));
            // No G/C -> 0.
            Assert.That(AnalysisTools.GcSkew("ATAT").GcSkew, Is.EqualTo(0.0).Within(1e-12));
        });
    }
}
