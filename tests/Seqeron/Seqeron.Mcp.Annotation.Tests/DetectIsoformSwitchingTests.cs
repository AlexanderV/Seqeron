using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class DetectIsoformSwitchingTests
{
    private static TranscriptIsoformDto Iso(string tx, string gene) =>
        new(tx, gene, 1000, 5, 0, true, new List<ExonRangeDto>());

    private static IsoformExpressionInputDto In(string tx, string gene, double e1, double e2) =>
        new(Iso(tx, gene), e1, e2);

    // GENE1: TX1 dominant in condition 1, TX2 dominant in condition 2 -> a switch.
    private static List<IsoformExpressionInputDto> Data() =>
        new()
        {
            In("TX1", "GENE1", 90, 10),
            In("TX2", "GENE1", 10, 90),
        };

    [Test]
    public void DetectIsoformSwitching_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.DetectIsoformSwitching(Data()));
        Assert.Throws<ArgumentNullException>(() => AnnotationTools.DetectIsoformSwitching(null!));
    }

    [Test]
    public void DetectIsoformSwitching_Binding_InvokesSuccessfully()
    {
        // TranscriptomeAnalyzer.DetectIsoformSwitching: usage = expr/total per condition; a switch needs
        // one isoform up by > threshold and another down by > threshold.
        // TX1 usage 0.9->0.1 (delta -0.8), TX2 0.1->0.9 (delta +0.8). switchScore = 0.8 + 0.8 = 1.6.
        // TranscriptId1 = decreased (TX1), TranscriptId2 = increased (TX2).
        var result = AnnotationTools.DetectIsoformSwitching(Data());

        Assert.That(result.Switches, Has.Count.EqualTo(1));
        Assert.Multiple(() =>
        {
            Assert.That(result.Switches[0].GeneId, Is.EqualTo("GENE1"));
            Assert.That(result.Switches[0].TranscriptId1, Is.EqualTo("TX1"));
            Assert.That(result.Switches[0].TranscriptId2, Is.EqualTo("TX2"));
            Assert.That(result.Switches[0].SwitchScore, Is.EqualTo(1.6).Within(1e-9));
        });
    }

    [Test]
    public void DetectIsoformSwitching_NoSwitch_ReturnsEmpty()
    {
        // Stable usage (no isoform crosses the threshold) -> no switch.
        var stable = new List<IsoformExpressionInputDto>
        {
            In("TX1", "GENE1", 60, 58),
            In("TX2", "GENE1", 40, 42),
        };
        Assert.That(AnnotationTools.DetectIsoformSwitching(stable).Switches, Is.Empty);
    }
}
