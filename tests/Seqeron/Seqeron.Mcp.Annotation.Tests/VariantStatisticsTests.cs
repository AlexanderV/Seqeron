using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class VariantStatisticsTests
{
    [Test]
    public void VariantStatistics_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.VariantStatistics("ATGC", "ATTC"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.VariantStatistics("", "ATTC"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.VariantStatistics("ATGC", ""));
    }

    [Test]
    public void VariantStatistics_SingleSnp_ReportsExactCounts()
    {
        // ATGC -> ATTC: one SNP (G>T, a transversion). Density = 1/4 * 1000 = 250.
        var s = AnnotationTools.VariantStatistics("ATGC", "ATTC");

        Assert.Multiple(() =>
        {
            Assert.That(s.TotalVariants, Is.EqualTo(1));
            Assert.That(s.Snps, Is.EqualTo(1));
            Assert.That(s.Insertions, Is.EqualTo(0));
            Assert.That(s.Deletions, Is.EqualTo(0));
            Assert.That(s.TiTvRatio, Is.EqualTo(0.0).Within(1e-10)); // 0 transitions / 1 transversion
            Assert.That(s.VariantDensity, Is.EqualTo(250.0).Within(1e-10));
            Assert.That(s.ReferenceLength, Is.EqualTo(4));
            Assert.That(s.QueryLength, Is.EqualTo(4));
        });
    }

    [Test]
    public void VariantStatistics_IdenticalSequences_ZeroVariants()
    {
        var s = AnnotationTools.VariantStatistics("ATGCATGC", "ATGCATGC");
        Assert.Multiple(() =>
        {
            Assert.That(s.TotalVariants, Is.EqualTo(0));
            Assert.That(s.VariantDensity, Is.EqualTo(0.0));
            Assert.That(s.ReferenceLength, Is.EqualTo(8));
        });
    }
}
