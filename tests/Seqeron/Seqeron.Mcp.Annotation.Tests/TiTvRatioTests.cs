using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class TiTvRatioTests
{
    private static VariantDto Snp(int pos, string refA, string altA) => new(pos, refA, altA, "SNP", pos);

    [Test]
    public void TiTvRatio_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.TiTvRatio(new[] { Snp(0, "A", "G") }));
        Assert.Throws<ArgumentNullException>(() => AnnotationTools.TiTvRatio(null!));
    }

    [Test]
    public void TiTvRatio_OneTransitionOneTransversion_ReturnsOne()
    {
        // Mirrors VariantCaller_CallVariants_Tests.CalculateTiTvRatio_OneTransitionOneTransversion_ReturnsOne.
        var result = AnnotationTools.TiTvRatio(new[] { Snp(0, "A", "G"), Snp(1, "A", "C") });
        Assert.That(result.Ratio, Is.EqualTo(1.0).Within(1e-10));
    }

    [Test]
    public void TiTvRatio_TwoTransitionsOneTransversion_ReturnsTwo()
    {
        // Mirrors CalculateTiTvRatio_TwoTransitionsOneTransversion_ReturnsTwo.
        var result = AnnotationTools.TiTvRatio(new[]
        {
            Snp(0, "A", "G"), Snp(1, "C", "T"), Snp(2, "A", "C")
        });
        Assert.That(result.Ratio, Is.EqualTo(2.0).Within(1e-10));
    }

    [Test]
    public void TiTvRatio_NoTransversions_ReturnsZero()
    {
        // Mirrors CalculateTiTvRatio_NoTransversions_ReturnsZero: no transversion -> 0.
        var result = AnnotationTools.TiTvRatio(new[] { Snp(0, "A", "G") });
        Assert.That(result.Ratio, Is.EqualTo(0.0).Within(1e-10));
    }
}
