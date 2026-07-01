using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class RnaseqQualityMetricsTests
{
    private static List<double> GeneCounts() => new() { 10, 0, 5, 0, 3 };

    [Test]
    public void RnaseqQualityMetrics_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.RnaseqQualityMetrics(1000, 900, 720, 45, GeneCounts()));
        Assert.Throws<ArgumentNullException>(() => AnnotationTools.RnaseqQualityMetrics(1000, 900, 720, 45, null!));
        Assert.Throws<ArgumentException>(() => AnnotationTools.RnaseqQualityMetrics(-1, 900, 720, 45, GeneCounts()));
        Assert.Throws<ArgumentException>(() => AnnotationTools.RnaseqQualityMetrics(1000, -1, 720, 45, GeneCounts()));
    }

    [Test]
    public void RnaseqQualityMetrics_Binding_InvokesSuccessfully()
    {
        // TranscriptomeAnalyzer.CalculateQualityMetrics: mappingRate = mapped/total, exonicRate =
        // exonic/mapped, rRnaRate = rRNA/mapped, detectedGenes = count(gene counts > 0).
        var result = AnnotationTools.RnaseqQualityMetrics(1000, 900, 720, 45, GeneCounts());

        Assert.Multiple(() =>
        {
            Assert.That(result.MappingRate, Is.EqualTo(0.9).Within(1e-9));
            Assert.That(result.ExonicRate, Is.EqualTo(0.8).Within(1e-9));   // 720/900
            Assert.That(result.RRnaRate, Is.EqualTo(0.05).Within(1e-9));    // 45/900
            Assert.That(result.DetectedGenes, Is.EqualTo(3));               // 10, 5, 3 > 0
        });
    }

    [Test]
    public void RnaseqQualityMetrics_ZeroTotals_YieldZeroRates()
    {
        var result = AnnotationTools.RnaseqQualityMetrics(0, 0, 0, 0, new List<double>());
        Assert.Multiple(() =>
        {
            Assert.That(result.MappingRate, Is.EqualTo(0));
            Assert.That(result.ExonicRate, Is.EqualTo(0));
            Assert.That(result.RRnaRate, Is.EqualTo(0));
            Assert.That(result.DetectedGenes, Is.EqualTo(0));
        });
    }
}
