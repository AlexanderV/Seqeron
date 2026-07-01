using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class ClusterSplitReadsTests
{
    private static SplitReadDto Read(string id, int primary) =>
        new(id, "chr1", primary, 2000, 30, "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");

    // Two split reads with primary positions within clusterDistance (10) -> one breakpoint.
    private static List<SplitReadDto> Reads() =>
        new() { Read("r1", 1000), Read("r2", 1004) };

    [Test]
    public void ClusterSplitReads_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.ClusterSplitReads(Reads()));
        Assert.Throws<ArgumentNullException>(() => AnnotationTools.ClusterSplitReads(null!));
    }

    [Test]
    public void ClusterSplitReads_Binding_InvokesSuccessfully()
    {
        // StructuralVariantAnalyzer.ClusterSplitReads: same-chromosome primaries within clusterDistance
        // and >= minSupport (2) form a breakpoint at the average primary/supplementary positions.
        // avg(1000, 1004) = 1002; supp avg = 2000; quality = min(2*15, 100) = 30.
        var result = AnnotationTools.ClusterSplitReads(Reads());

        Assert.That(result.Breakpoints, Has.Count.EqualTo(1));
        Assert.Multiple(() =>
        {
            Assert.That(result.Breakpoints[0].Position1, Is.EqualTo(1002));
            Assert.That(result.Breakpoints[0].Position2, Is.EqualTo(2000));
            Assert.That(result.Breakpoints[0].SupportingReads, Is.EqualTo(2));
            Assert.That(result.Breakpoints[0].Quality, Is.EqualTo(30.0).Within(1e-9));
        });
    }

    [Test]
    public void ClusterSplitReads_BelowMinSupport_ReturnsEmpty()
    {
        // A single split read is below the default minSupport of 2.
        var result = AnnotationTools.ClusterSplitReads(new List<SplitReadDto> { Read("r1", 1000) });
        Assert.That(result.Breakpoints, Is.Empty);
    }
}
