using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class DetectDifferentialSplicingTests
{
    private static List<SplicingComparisonInputDto> Data() =>
        new()
        {
            new SplicingComparisonInputDto("G1", 100, 200, 0.2, 0.8),  // deltaPSI +0.6 -> IncreasedInclusion
            new SplicingComparisonInputDto("G2", 300, 400, 0.50, 0.52), // deltaPSI 0.02 -> below threshold
            new SplicingComparisonInputDto("G3", 500, 600, 0.8, 0.3),  // deltaPSI -0.5 -> IncreasedSkipping
        };

    [Test]
    public void DetectDifferentialSplicing_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.DetectDifferentialSplicing(Data()));
        Assert.Throws<ArgumentNullException>(() => AnnotationTools.DetectDifferentialSplicing(null!));
    }

    [Test]
    public void DetectDifferentialSplicing_Binding_InvokesSuccessfully()
    {
        // TranscriptomeAnalyzer.DetectDifferentialSplicing: deltaPSI = psi2 - psi1; reported when
        // |deltaPSI| >= threshold (0.1). Positive -> IncreasedInclusion, negative -> IncreasedSkipping.
        var result = AnnotationTools.DetectDifferentialSplicing(Data());

        Assert.That(result.Events, Has.Count.EqualTo(2));
        Assert.Multiple(() =>
        {
            var g1 = result.Events.Single(e => e.GeneId == "G1");
            Assert.That(g1.EventType, Is.EqualTo("IncreasedInclusion"));
            Assert.That(g1.DeltaPsi, Is.EqualTo(0.6).Within(1e-9));
            Assert.That(g1.InclusionLevel, Is.EqualTo(0.8).Within(1e-9));

            var g3 = result.Events.Single(e => e.GeneId == "G3");
            Assert.That(g3.EventType, Is.EqualTo("IncreasedSkipping"));
            Assert.That(g3.DeltaPsi, Is.EqualTo(-0.5).Within(1e-9));

            Assert.That(result.Events.Any(e => e.GeneId == "G2"), Is.False);
        });
    }

    [Test]
    public void DetectDifferentialSplicing_CustomThreshold()
    {
        // A stricter threshold of 0.55 keeps only G1 (0.6).
        var result = AnnotationTools.DetectDifferentialSplicing(Data(), deltaPsiThreshold: 0.55);
        Assert.That(result.Events.Select(e => e.GeneId), Is.EquivalentTo(new[] { "G1" }));
    }
}
