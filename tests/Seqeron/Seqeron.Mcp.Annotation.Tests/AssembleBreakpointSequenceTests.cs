using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class AssembleBreakpointSequenceTests
{
    // Fixtures mirror StructuralVariantAnalyzerTests.AssembleBreakpointSequence_MultipleSplits_ReturnsLongest.
    private static List<SplitReadDto> Splits() => new()
    {
        new("r1", "chr1", 1000, 5000, 20, "ACGTACGTACGTACGTACGT"),
        new("r2", "chr1", 1005, 5005, 50, "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC"),
        new("r3", "chr1", 1010, 5010, 30, "ACGTACGTACGTACGTACGTACGTACGTAC")
    };

    [Test]
    public void AssembleBreakpointSequence_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.AssembleBreakpointSequence(Splits()));
        Assert.DoesNotThrow(() => AnnotationTools.AssembleBreakpointSequence(new List<SplitReadDto>()));
        Assert.Throws<ArgumentNullException>(() => AnnotationTools.AssembleBreakpointSequence(null!));
    }

    [Test]
    public void AssembleBreakpointSequence_MultipleSplits_ReturnsLongestClippedSequence()
    {
        // Assembly returns the ClippedSequence of the read with the largest ClipLength (r2, 50 nt).
        var result = AnnotationTools.AssembleBreakpointSequence(Splits());

        Assert.That(result.Sequence, Is.Not.Null);
        Assert.Multiple(() =>
        {
            Assert.That(result.Sequence!.Length, Is.EqualTo(50));
            Assert.That(result.Sequence, Is.EqualTo("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC"));
        });
    }

    [Test]
    public void AssembleBreakpointSequence_EmptyInput_ReturnsNull()
    {
        var result = AnnotationTools.AssembleBreakpointSequence(new List<SplitReadDto>());
        Assert.That(result.Sequence, Is.Null);
    }
}
