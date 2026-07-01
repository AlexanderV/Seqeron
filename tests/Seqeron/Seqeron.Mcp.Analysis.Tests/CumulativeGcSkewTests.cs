using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>cumulative_gc_skew</c> MCP tool.
/// Expected values computed by hand from GcSkewCalculator.CalculateCumulativeGcSkew:
/// per non-overlapping window (step == windowSize) skew = (G-C)/(G+C), cumulative is the
/// running sum, position = windowStart + windowSize/2. NOT from the wrapper's output.
/// </summary>
[TestFixture]
public class CumulativeGcSkewTests
{
    [Test]
    public void CumulativeGcSkew_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.CumulativeGcSkew("GGGGCCCC", 4));
        Assert.Throws<ArgumentException>(() => AnalysisTools.CumulativeGcSkew("", 4));
        Assert.Throws<ArgumentException>(() => AnalysisTools.CumulativeGcSkew(null!, 4));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.CumulativeGcSkew("GGGGCCCC", 0));
    }

    [Test]
    public void CumulativeGcSkew_Binding_InvokesSuccessfully()
    {
        // "GGGGCCCC", windowSize 4 (non-overlapping):
        //   window0 "GGGG": skew = (4-0)/4 = +1.0, cumulative = 1.0, position = 0 + 2 = 2
        //   window1 "CCCC": skew = (0-4)/4 = -1.0, cumulative = 0.0, position = 4 + 2 = 6
        var items = AnalysisTools.CumulativeGcSkew("GGGGCCCC", 4).Items;
        Assert.That(items, Has.Length.EqualTo(2));
        Assert.Multiple(() =>
        {
            Assert.That(items[0].Position, Is.EqualTo(2));
            Assert.That(items[0].GcSkew, Is.EqualTo(1.0).Within(1e-12));
            Assert.That(items[0].CumulativeGcSkew, Is.EqualTo(1.0).Within(1e-12));

            Assert.That(items[1].Position, Is.EqualTo(6));
            Assert.That(items[1].GcSkew, Is.EqualTo(-1.0).Within(1e-12));
            Assert.That(items[1].CumulativeGcSkew, Is.EqualTo(0.0).Within(1e-12));
        });
    }
}
