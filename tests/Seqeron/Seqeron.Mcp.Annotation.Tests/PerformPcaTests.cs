using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class PerformPcaTests
{
    // Two samples, 4 genes; every gene varies, so all are selected (topGenes >= geneCount).
    private static List<SampleExpressionInputDto> Samples() =>
        new()
        {
            new SampleExpressionInputDto("A", new List<double> { 4, 4, 4, 4 }),
            new SampleExpressionInputDto("B", new List<double> { 0, 0, 0, 0 }),
        };

    [Test]
    public void PerformPca_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.PerformPca(Samples()));
        Assert.Throws<ArgumentNullException>(() => AnnotationTools.PerformPca(null!));
        Assert.Throws<ArgumentException>(() => AnnotationTools.PerformPca(new List<SampleExpressionInputDto>()));
        // Unequal expression lengths are rejected.
        var jagged = new List<SampleExpressionInputDto>
        {
            new SampleExpressionInputDto("A", new List<double> { 1, 2, 3 }),
            new SampleExpressionInputDto("B", new List<double> { 1, 2 }),
        };
        Assert.Throws<ArgumentException>(() => AnnotationTools.PerformPca(jagged));
    }

    [Test]
    public void PerformPca_Binding_InvokesSuccessfully()
    {
        // TranscriptomeAnalyzer.PerformPCA: selects top-variance genes, then approximates
        // PC1 = sum of the first half of selected values, PC2 = sum of the second half.
        // All 4 genes selected; A -> PC1 = 4+4 = 8, PC2 = 8; B -> 0, 0.
        var result = AnnotationTools.PerformPca(Samples());

        Assert.That(result.Points, Has.Count.EqualTo(2));
        Assert.Multiple(() =>
        {
            var a = result.Points.Single(p => p.SampleId == "A");
            Assert.That(a.Pc1, Is.EqualTo(8.0).Within(1e-9));
            Assert.That(a.Pc2, Is.EqualTo(8.0).Within(1e-9));

            var b = result.Points.Single(p => p.SampleId == "B");
            Assert.That(b.Pc1, Is.EqualTo(0.0).Within(1e-9));
            Assert.That(b.Pc2, Is.EqualTo(0.0).Within(1e-9));
        });
    }

    [Test]
    public void PerformPca_SingleSample_ProjectsToOrigin()
    {
        // Fewer than 2 samples -> PC scores are 0.
        var one = new List<SampleExpressionInputDto>
        {
            new SampleExpressionInputDto("A", new List<double> { 1, 2, 3 }),
        };
        var result = AnnotationTools.PerformPca(one);
        Assert.Multiple(() =>
        {
            Assert.That(result.Points, Has.Count.EqualTo(1));
            Assert.That(result.Points[0].Pc1, Is.EqualTo(0.0));
            Assert.That(result.Points[0].Pc2, Is.EqualTo(0.0));
        });
    }
}
