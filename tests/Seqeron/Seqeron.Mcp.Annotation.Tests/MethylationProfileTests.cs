using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class MethylationProfileTests
{
    private static List<MethylationSiteDto> Sites() =>
        new()
        {
            new MethylationSiteDto(10, "CpG", "CGT", 1.0, 10),
            new MethylationSiteDto(20, "CpG", "CGA", 0.0, 10),
            new MethylationSiteDto(30, "CHG", "CAG", 0.5, 4),
        };

    [Test]
    public void MethylationProfile_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.MethylationProfile(Sites()));
        Assert.Throws<ArgumentNullException>(() => AnnotationTools.MethylationProfile(null!));
        // Unknown methylation type is rejected by the DTO -> domain mapping.
        var bad = new List<MethylationSiteDto> { new MethylationSiteDto(1, "XYZ", "C", 0.5, 1) };
        Assert.Throws<ArgumentException>(() => AnnotationTools.MethylationProfile(bad));
    }

    [Test]
    public void MethylationProfile_Binding_InvokesSuccessfully()
    {
        // EpigeneticsAnalyzer.GenerateMethylationProfile: per-context weighted level (Schultz 2012)
        //   = Σ(level·coverage) / Σ(coverage).
        // Global: (1.0*10 + 0.0*10 + 0.5*4)/24 = 12/24 = 0.5.
        // CpG: (1.0*10 + 0.0*10)/20 = 0.5.  CHG: (0.5*4)/4 = 0.5.  CHH: no sites -> 0.
        // TotalCpG = 2; MethylatedCpG (level >= 0.5) = 1 (the fully methylated site).
        var result = AnnotationTools.MethylationProfile(Sites());

        Assert.Multiple(() =>
        {
            Assert.That(result.GlobalMethylation, Is.EqualTo(0.5).Within(1e-9));
            Assert.That(result.CpGMethylation, Is.EqualTo(0.5).Within(1e-9));
            Assert.That(result.CHGMethylation, Is.EqualTo(0.5).Within(1e-9));
            Assert.That(result.CHHMethylation, Is.EqualTo(0.0).Within(1e-9));
            Assert.That(result.TotalCpGSites, Is.EqualTo(2));
            Assert.That(result.MethylatedCpGSites, Is.EqualTo(1));
            Assert.That(result.MethylationByPosition, Has.Count.EqualTo(3));
            Assert.That(result.MethylationByPosition[0].Position, Is.EqualTo(10));
            Assert.That(result.MethylationByPosition[0].Level, Is.EqualTo(1.0).Within(1e-9));
        });
    }

    [Test]
    public void MethylationProfile_EmptyInput_ReturnsZeroProfile()
    {
        var result = AnnotationTools.MethylationProfile(new List<MethylationSiteDto>());
        Assert.Multiple(() =>
        {
            Assert.That(result.GlobalMethylation, Is.EqualTo(0));
            Assert.That(result.TotalCpGSites, Is.EqualTo(0));
            Assert.That(result.MethylationByPosition, Is.Empty);
        });
    }
}
