using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class GenotypeSvTests
{
    private static StructuralVariantDto Sv() =>
        new("sv1", "chr1", 1000, 2000, "Deletion", 1000, 40, 5, null);

    [Test]
    public void GenotypeSv_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.GenotypeSv(Sv(), 10, 10, 20));
        Assert.Throws<ArgumentNullException>(() => AnnotationTools.GenotypeSv(null!, 10, 10, 20));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnnotationTools.GenotypeSv(Sv(), -1, 10, 20));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnnotationTools.GenotypeSv(Sv(), 10, -1, 20));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnnotationTools.GenotypeSv(Sv(), 10, 10, -1));
    }

    [Test]
    public void GenotypeSv_Binding_InvokesSuccessfully()
    {
        // StructuralVariantAnalyzer.GenotypeSV: alt fraction determines genotype.
        Assert.Multiple(() =>
        {
            // Heterozygous: altFraction 0.5 in [0.3,0.7] -> "0/1", quality (ref+alt)*2 = 40.
            var het = AnnotationTools.GenotypeSv(Sv(), 10, 10, 20);
            Assert.That(het.Genotype, Is.EqualTo("0/1"));
            Assert.That(het.Quality, Is.EqualTo(40.0).Within(1e-9));

            // Homozygous reference: altFraction 0 < 0.1 -> "0/0", quality ref*3 = 60.
            var homRef = AnnotationTools.GenotypeSv(Sv(), 20, 0, 20);
            Assert.That(homRef.Genotype, Is.EqualTo("0/0"));
            Assert.That(homRef.Quality, Is.EqualTo(60.0).Within(1e-9));

            // Homozygous alternate: altFraction 20/20 = 1.0 > 0.9 -> "1/1", quality alt*3 capped at 99.
            var homAlt = AnnotationTools.GenotypeSv(Sv(), 0, 20, 20);
            Assert.That(homAlt.Genotype, Is.EqualTo("1/1"));
            Assert.That(homAlt.Quality, Is.EqualTo(60.0).Within(1e-9));

            // No reads -> missing genotype.
            var missing = AnnotationTools.GenotypeSv(Sv(), 0, 0, 0);
            Assert.That(missing.Genotype, Is.EqualTo("./."));
            Assert.That(missing.Quality, Is.EqualTo(0));
        });
    }

    [Test]
    public void GenotypeSv_QualityCappedAt99()
    {
        // Homozygous alt with alt*3 = 150 -> capped at 99.
        var result = AnnotationTools.GenotypeSv(Sv(), 0, 50, 50);
        Assert.That(result.Genotype, Is.EqualTo("1/1"));
        Assert.That(result.Quality, Is.EqualTo(99.0).Within(1e-9));
    }
}
