using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>windowed_gc_skew</c> MCP tool.
/// Expected values from per-window (G-C)/(G+C) with Position = i + windowSize/2 and
/// window bounds [i, i+windowSize-1], NOT the wrapper output.
/// </summary>
[TestFixture]
public class WindowedGcSkewTests
{
    [Test]
    public void WindowedGcSkew_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.WindowedGcSkew("GGGGCCCC", 4, 4));
        Assert.Throws<ArgumentException>(() => AnalysisTools.WindowedGcSkew("", 4, 4));
        Assert.Throws<ArgumentException>(() => AnalysisTools.WindowedGcSkew(null!, 4, 4));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.WindowedGcSkew("GGGGCCCC", 0, 4));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.WindowedGcSkew("GGGGCCCC", 4, 0));
    }

    [Test]
    public void WindowedGcSkew_Binding_InvokesSuccessfully()
    {
        // window 4 step 4 -> GGGG (skew 1, center 2) then CCCC (skew -1, center 6).
        var pts = AnalysisTools.WindowedGcSkew("GGGGCCCC", 4, 4).Items;
        Assert.Multiple(() =>
        {
            Assert.That(pts, Has.Length.EqualTo(2));
            Assert.That(pts[0].GcSkew, Is.EqualTo(1.0).Within(1e-12));
            Assert.That(pts[0].Position, Is.EqualTo(2));
            Assert.That(pts[0].WindowStart, Is.EqualTo(0));
            Assert.That(pts[0].WindowEnd, Is.EqualTo(3));
            Assert.That(pts[1].GcSkew, Is.EqualTo(-1.0).Within(1e-12));
            Assert.That(pts[1].Position, Is.EqualTo(6));
            Assert.That(pts[1].WindowStart, Is.EqualTo(4));
            Assert.That(pts[1].WindowEnd, Is.EqualTo(7));
        });
    }
}
