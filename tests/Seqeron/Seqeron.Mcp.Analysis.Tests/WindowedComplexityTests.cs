using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>windowed_complexity</c> MCP tool.
/// Expected values from SequenceComplexity's own unit test
/// (SequenceComplexity_CalculateWindowedComplexity_Tests, "ACGTACGTAAAAAAAAACGTACGT",
/// w=8 s=8: 3 windows, uniform window entropy 2.0, homopolymer window entropy 0),
/// NOT the wrapper output.
/// </summary>
[TestFixture]
public class WindowedComplexityTests
{
    private const string Profiled24 = "ACGTACGTAAAAAAAAACGTACGT";

    [Test]
    public void WindowedComplexity_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.WindowedComplexity(Profiled24, 8, 8));
        Assert.Throws<ArgumentException>(() => AnalysisTools.WindowedComplexity("", 8, 8));
        Assert.Throws<ArgumentException>(() => AnalysisTools.WindowedComplexity(null!, 8, 8));
        Assert.Throws<ArgumentException>(() => AnalysisTools.WindowedComplexity("XYZ", 8, 8));
    }

    [Test]
    public void WindowedComplexity_Binding_InvokesSuccessfully()
    {
        var pts = AnalysisTools.WindowedComplexity(Profiled24, 8, 8).Items;
        Assert.Multiple(() =>
        {
            Assert.That(pts, Has.Length.EqualTo(3));
            Assert.That(pts.Select(p => p.WindowStart), Is.EqualTo(new[] { 0, 8, 16 }));
            Assert.That(pts.Select(p => p.WindowEnd), Is.EqualTo(new[] { 7, 15, 23 }));
            Assert.That(pts.Select(p => p.Position), Is.EqualTo(new[] { 4, 12, 20 }));
            // Window 0 "ACGTACGT" is uniform -> Shannon entropy log2(4) = 2.0.
            Assert.That(pts[0].ShannonEntropy, Is.EqualTo(2.0).Within(1e-10));
            // Window 1 "AAAAAAAA" is a homopolymer -> Shannon entropy 0.
            Assert.That(pts[1].ShannonEntropy, Is.EqualTo(0.0).Within(1e-10));
        });
    }
}
