using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class FindDmrsTests
{
    private static MethylationSiteDto Cpg(int pos, double level, int cov) =>
        new(pos, "CpG", "CGT", level, cov);

    // Sample 1 unmethylated, sample 2 fully methylated at positions 10/20/30 (one 1000-bp window).
    private static List<MethylationSiteDto> Control() =>
        new() { Cpg(10, 0.0, 10), Cpg(20, 0.0, 10), Cpg(30, 0.0, 10) };

    private static List<MethylationSiteDto> Treatment() =>
        new() { Cpg(10, 1.0, 10), Cpg(20, 1.0, 10), Cpg(30, 1.0, 10) };

    [Test]
    public void FindDmrs_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.FindDmrs(Control(), Treatment()));
        Assert.Throws<ArgumentNullException>(() => AnnotationTools.FindDmrs(null!, Treatment()));
        Assert.Throws<ArgumentNullException>(() => AnnotationTools.FindDmrs(Control(), null!));
        Assert.Throws<ArgumentException>(
            () => AnnotationTools.FindDmrs(new List<MethylationSiteDto>(), new List<MethylationSiteDto>()));
    }

    [Test]
    public void FindDmrs_Binding_InvokesSuccessfully()
    {
        // EpigeneticsAnalyzer.FindDMRs (methylKit tiling): positions 10/20/30 fall in one 1000-bp
        // window (>= minCpGCount=3). meanDiff = mean(level2 - level1) = 1.0 > 0.25 -> Hypermethylated.
        // Region spans start=10 (first position) to end=30 (last position), CpGCount=3.
        var result = AnnotationTools.FindDmrs(Control(), Treatment());

        Assert.That(result.Regions, Has.Count.EqualTo(1));
        Assert.Multiple(() =>
        {
            Assert.That(result.Regions[0].Start, Is.EqualTo(10));
            Assert.That(result.Regions[0].End, Is.EqualTo(30));
            Assert.That(result.Regions[0].MeanDifference, Is.EqualTo(1.0).Within(1e-9));
            Assert.That(result.Regions[0].CpGCount, Is.EqualTo(3));
            Assert.That(result.Regions[0].Annotation, Is.EqualTo("Hypermethylated"));
            // Pooled 2x2 table (0/30 vs 30/0) is maximally significant.
            Assert.That(result.Regions[0].PValue, Is.LessThan(0.01));
        });
    }

    [Test]
    public void FindDmrs_BelowThreshold_NotReported()
    {
        // Identical samples -> meanDiff = 0 <= minDifference -> no DMR.
        var result = AnnotationTools.FindDmrs(Control(), Control());
        Assert.That(result.Regions, Is.Empty);
    }

    [Test]
    public void FindDmrs_HypomethylatedDirection()
    {
        // Control methylated, treatment unmethylated -> meanDiff = -1.0 -> Hypomethylated.
        var result = AnnotationTools.FindDmrs(Treatment(), Control());
        Assert.That(result.Regions, Has.Count.EqualTo(1));
        Assert.That(result.Regions[0].MeanDifference, Is.EqualTo(-1.0).Within(1e-9));
        Assert.That(result.Regions[0].Annotation, Is.EqualTo("Hypomethylated"));
    }
}
