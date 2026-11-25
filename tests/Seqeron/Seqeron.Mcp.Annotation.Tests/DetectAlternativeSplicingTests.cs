using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class DetectAlternativeSplicingTests
{
    // Two GU-AG introns (CAGGUAAGU + 60 A + PPT + CAGG) separated by a 50-nt exon: the spanning candidate
    // intron makes the middle exon an SE (ExonSkipping) event (SUPPA / rMATS definition).
    private const string Unit = "CAGGUAAGU" + "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" + "UUUUUUUUUUUUUUUCAGG";
    private const string TwoDonors = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC" + Unit +
        "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC" + Unit + "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC";

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
        // ExonSkipping is reported at the donor (GU) of the first flanking intron (index 40 + 3).
        var result = AnnotationTools.DetectAlternativeSplicing(TwoDonors, minScore: 0.3);

        var se = result.Events.Where(e => e.Type == "ExonSkipping").ToList();
        Assert.That(se, Is.Not.Empty);
        Assert.That(se.Select(e => e.Position), Does.Contain(43));
    }

    [Test]
    public void DetectAlternativeSplicing_EmptySequence_ReturnsEmpty()
    {
        // The underlying algorithm returns empty for empty input; the wrapper guards empty first.
        Assert.Throws<ArgumentException>(() => AnnotationTools.DetectAlternativeSplicing(""));
        // A sequence without a candidate intron yields no event.
        var single = AnnotationTools.DetectAlternativeSplicing("CAGGUAAGU", minScore: 0.2);
        Assert.That(single.Events, Is.Empty);
    }
}
