using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class QuantileNormalizeTests
{
    private static List<IReadOnlyList<double>> Samples() =>
        new()
        {
            new List<double> { 2.0, 3.0, 4.0 },
            new List<double> { 8.0, 6.0, 5.0 },
        };

    [Test]
    public void QuantileNormalize_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.QuantileNormalize(Samples()));
        Assert.Throws<ArgumentNullException>(() => AnnotationTools.QuantileNormalize(null!));
        Assert.Throws<ArgumentException>(() => AnnotationTools.QuantileNormalize(new List<IReadOnlyList<double>>()));
        // Unequal lengths are rejected.
        var jagged = new List<IReadOnlyList<double>>
        {
            new List<double> { 1, 2, 3 },
            new List<double> { 1, 2 },
        };
        Assert.Throws<ArgumentException>(() => AnnotationTools.QuantileNormalize(jagged));
    }

    [Test]
    public void QuantileNormalize_Binding_InvokesSuccessfully()
    {
        // TranscriptomeAnalyzer.QuantileNormalize (Bolstad 2003): rank means over sorted columns
        // reassigned by within-column rank. rankMeans = [3.5, 4.5, 6.0].
        var result = AnnotationTools.QuantileNormalize(Samples());

        Assert.That(result.Normalized, Has.Count.EqualTo(2));
        Assert.Multiple(() =>
        {
            Assert.That(result.Normalized[0], Is.EqualTo(new[] { 3.5, 4.5, 6.0 }).Within(1e-9));
            Assert.That(result.Normalized[1], Is.EqualTo(new[] { 6.0, 4.5, 3.5 }).Within(1e-9));
        });
    }
}
