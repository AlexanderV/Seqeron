using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class FindSkippedExonEventsTests
{
    private static List<SkippedExonInputDto> ExonData() =>
        new()
        {
            new SkippedExonInputDto("G1", 100, 200, 80, 20),  // PSI 0.8
            new SkippedExonInputDto("G2", 300, 400, 0, 0),    // no reads -> skipped
        };

    [Test]
    public void FindSkippedExonEvents_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.FindSkippedExonEvents(ExonData()));
        Assert.Throws<ArgumentNullException>(() => AnnotationTools.FindSkippedExonEvents(null!));
    }

    [Test]
    public void FindSkippedExonEvents_Binding_InvokesSuccessfully()
    {
        // TranscriptomeAnalyzer.FindSkippedExonEvents: PSI = inclusion / (inclusion + skipping).
        // G1: 80/(80+20) = 0.8, EventType "SkippedExon". G2 has zero total reads and is skipped.
        var result = AnnotationTools.FindSkippedExonEvents(ExonData());

        Assert.That(result.Events, Has.Count.EqualTo(1));
        Assert.Multiple(() =>
        {
            Assert.That(result.Events[0].GeneId, Is.EqualTo("G1"));
            Assert.That(result.Events[0].EventType, Is.EqualTo("SkippedExon"));
            Assert.That(result.Events[0].Start, Is.EqualTo(100));
            Assert.That(result.Events[0].End, Is.EqualTo(200));
            Assert.That(result.Events[0].InclusionLevel, Is.EqualTo(0.8).Within(1e-9));
            Assert.That(result.Events[0].DeltaPsi, Is.EqualTo(0.0));
        });
    }

    [Test]
    public void FindSkippedExonEvents_NoReads_ReturnsEmpty()
    {
        var data = new List<SkippedExonInputDto> { new SkippedExonInputDto("G", 1, 2, 0, 0) };
        Assert.That(AnnotationTools.FindSkippedExonEvents(data).Events, Is.Empty);
    }
}
