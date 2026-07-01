using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class FilterSvsTests
{
    private static StructuralVariantDto Sv(string id, int length, double q, int support) =>
        new(id, "chr1", 1000, 1000 + length, "Deletion", length, q, support, null);

    private static List<StructuralVariantDto> Variants() =>
        new()
        {
            Sv("pass",     500, 30, 5),   // quality 30>=20, support 5>=2, len 500 in [50,1e8] -> pass
            Sv("lowqual",  500, 10, 5),   // quality 10 < 20 -> filtered
            Sv("lowsupp",  500, 30, 1),   // support 1 < 2 -> filtered
            Sv("tooshort",  10, 30, 5),   // length 10 < 50 -> filtered
        };

    [Test]
    public void FilterSvs_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.FilterSvs(Variants()));
        Assert.Throws<ArgumentNullException>(() => AnnotationTools.FilterSvs(null!));
    }

    [Test]
    public void FilterSvs_Binding_InvokesSuccessfully()
    {
        // StructuralVariantAnalyzer.FilterSVs keeps SVs with quality>=minQuality (20), support>=minSupport (2),
        // and length within [minLength (50), maxLength (1e8)].
        var result = AnnotationTools.FilterSvs(Variants());
        Assert.That(result.Variants.Select(v => v.Id), Is.EquivalentTo(new[] { "pass" }));
    }

    [Test]
    public void FilterSvs_CustomThresholds()
    {
        // Relaxing minQuality/minSupport/minLength lets more variants through.
        var result = AnnotationTools.FilterSvs(Variants(), minQuality: 0, minSupport: 1, minLength: 1);
        Assert.That(result.Variants.Select(v => v.Id),
            Is.EquivalentTo(new[] { "pass", "lowqual", "lowsupp", "tooshort" }));
    }
}
