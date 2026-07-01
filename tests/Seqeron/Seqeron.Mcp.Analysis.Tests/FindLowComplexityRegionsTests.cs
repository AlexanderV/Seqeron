using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>find_low_complexity_regions</c> MCP tool.
/// Expected values from SequenceComplexity's own unit test
/// (SequenceComplexityTests.FindLowComplexityRegions_FindsPolyARegion: ATGCx20 + A64 +
/// ATGCx20, w=20 thr=0.5 -> region 79..146, minEntropy 0), NOT the wrapper output.
/// </summary>
[TestFixture]
public class FindLowComplexityRegionsTests
{
    private static string PolyAFlanked() =>
        string.Concat(Enumerable.Repeat("ATGC", 20)) + new string('A', 64) + string.Concat(Enumerable.Repeat("ATGC", 20));

    [Test]
    public void FindLowComplexityRegions_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.FindLowComplexityRegions(PolyAFlanked(), 20, 0.5));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindLowComplexityRegions("", 20, 0.5));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindLowComplexityRegions(null!, 20, 0.5));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindLowComplexityRegions("XYZ", 20, 0.5));
    }

    [Test]
    public void FindLowComplexityRegions_Binding_InvokesSuccessfully()
    {
        // The internal poly-A tract is the single low-complexity region at 79..146.
        var regions = AnalysisTools.FindLowComplexityRegions(PolyAFlanked(), 20, 0.5).Items;
        Assert.Multiple(() =>
        {
            Assert.That(regions, Has.Length.EqualTo(1));
            Assert.That(regions[0].Start, Is.EqualTo(79));
            Assert.That(regions[0].End, Is.EqualTo(146));
            Assert.That(regions[0].MinEntropy, Is.EqualTo(0.0).Within(1e-10));
        });

        // A pure ATGC repeat has uniformly high entropy -> no low-complexity region.
        var none = AnalysisTools.FindLowComplexityRegions(string.Concat(Enumerable.Repeat("ATGC", 20)), 20, 0.5).Items;
        Assert.That(none, Is.Empty);
    }
}
