using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class DetectAlternativeSplicingTests
{
    // Two perfect GU donors (CAGGUAAGU) within the same 50-nt group -> Alt5SS event.
    private const string TwoDonors = "CAGGUAAGUCAGGUAAGU";

    [Test]
    public void DetectAlternativeSplicing_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.DetectAlternativeSplicing(TwoDonors, minScore: 0.2));
        Assert.Throws<ArgumentException>(() => AnnotationTools.DetectAlternativeSplicing(""));
        Assert.Throws<ArgumentException>(() => AnnotationTools.DetectAlternativeSplicing(null!));
    }

    [Test]
    public void DetectAlternativeSplicing_Binding_InvokesSuccessfully()
    {
        // SpliceSitePredictor.DetectAlternativeSplicing: >1 donor in a 50-nt window -> Alt5SS event,
        // reported at the first donor's position (the GU of the first CAGGUAAGU is at index 3).
        var result = AnnotationTools.DetectAlternativeSplicing(TwoDonors, minScore: 0.2);

        Assert.That(result.Events, Is.Not.Empty);
        var alt5 = result.Events.Where(e => e.Type == "Alt5SS").ToList();
        Assert.That(alt5, Is.Not.Empty);
        Assert.That(alt5[0].Position, Is.EqualTo(3));
    }

    [Test]
    public void DetectAlternativeSplicing_EmptySequence_ReturnsEmpty()
    {
        // The underlying algorithm returns empty for empty input; the wrapper guards empty first.
        Assert.Throws<ArgumentException>(() => AnnotationTools.DetectAlternativeSplicing(""));
        // A sequence with a single donor and no grouping yields no Alt5SS.
        var single = AnnotationTools.DetectAlternativeSplicing("CAGGUAAGU", minScore: 0.2);
        Assert.That(single.Events.Any(e => e.Type == "Alt5SS"), Is.False);
    }
}
