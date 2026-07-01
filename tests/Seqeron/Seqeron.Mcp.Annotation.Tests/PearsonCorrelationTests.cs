using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class PearsonCorrelationTests
{
    [Test]
    public void PearsonCorrelation_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.PearsonCorrelation(
            new List<double> { 1, 2, 3 }, new List<double> { 1, 2, 3 }));
        Assert.Throws<ArgumentNullException>(() => AnnotationTools.PearsonCorrelation(null!, new List<double> { 1 }));
        Assert.Throws<ArgumentNullException>(() => AnnotationTools.PearsonCorrelation(new List<double> { 1 }, null!));
        Assert.Throws<ArgumentException>(() => AnnotationTools.PearsonCorrelation(new List<double>(), new List<double> { 1 }));
        // Unequal lengths are rejected.
        Assert.Throws<ArgumentException>(() => AnnotationTools.PearsonCorrelation(
            new List<double> { 1, 2, 3 }, new List<double> { 1, 2 }));
    }

    [Test]
    public void PearsonCorrelation_Binding_InvokesSuccessfully()
    {
        // Mirrors Seqeron.Genomics.Tests: perfectly positive/negative linear relationships.
        Assert.Multiple(() =>
        {
            var pos = AnnotationTools.PearsonCorrelation(
                new List<double> { 1, 2, 3, 4, 5 }, new List<double> { 2, 4, 6, 8, 10 });
            Assert.That(pos.Correlation, Is.EqualTo(1.0).Within(1e-9));

            var neg = AnnotationTools.PearsonCorrelation(
                new List<double> { 1, 2, 3, 4, 5 }, new List<double> { 5, 4, 3, 2, 1 });
            Assert.That(neg.Correlation, Is.EqualTo(-1.0).Within(1e-9));
        });
    }
}
