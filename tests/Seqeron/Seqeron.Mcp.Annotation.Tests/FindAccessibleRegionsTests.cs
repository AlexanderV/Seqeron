using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class FindAccessibleRegionsTests
{
    // Positions 0,10,...,200 all at signal 0.9 (gaps of 10 < default maxGap 50).
    private static List<AccessibilitySignalDto> HighSignal() =>
        Enumerable.Range(0, 21).Select(i => new AccessibilitySignalDto(i * 10, 0.9)).ToList();

    [Test]
    public void FindAccessibleRegions_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.FindAccessibleRegions(HighSignal()));
        Assert.Throws<ArgumentNullException>(() => AnnotationTools.FindAccessibleRegions(null!));
    }

    [Test]
    public void FindAccessibleRegions_Binding_InvokesSuccessfully()
    {
        // EpigeneticsAnalyzer.FindAccessibleRegions: defaults threshold=0.5, minWidth=100, maxGap=50.
        // All 21 positions (0..200) exceed threshold with gaps of 10 -> one region 0..200,
        // width 200 >= 100, max signal 0.9 -> PeakType "Strong" (>0.8).
        var result = AnnotationTools.FindAccessibleRegions(HighSignal());

        Assert.That(result.Regions, Has.Count.EqualTo(1));
        Assert.Multiple(() =>
        {
            Assert.That(result.Regions[0].Start, Is.EqualTo(0));
            Assert.That(result.Regions[0].End, Is.EqualTo(200));
            Assert.That(result.Regions[0].AccessibilityScore, Is.EqualTo(0.9).Within(1e-9));
            Assert.That(result.Regions[0].PeakType, Is.EqualTo("Strong"));
            Assert.That(result.Regions[0].NearbyGenes, Is.Empty);
        });
    }

    [Test]
    public void FindAccessibleRegions_BelowThreshold_ReturnsEmpty()
    {
        // All signal below the 0.5 threshold -> no region opens.
        var lowSignal = Enumerable.Range(0, 21).Select(i => new AccessibilitySignalDto(i * 10, 0.1)).ToList();
        Assert.That(AnnotationTools.FindAccessibleRegions(lowSignal).Regions, Is.Empty);
    }

    [Test]
    public void FindAccessibleRegions_TooNarrow_ReturnsEmpty()
    {
        // A 50 bp accessible stretch is below the default minWidth of 100.
        var narrow = Enumerable.Range(0, 6).Select(i => new AccessibilitySignalDto(i * 10, 0.9)).ToList();
        Assert.That(AnnotationTools.FindAccessibleRegions(narrow).Regions, Is.Empty);
    }
}
