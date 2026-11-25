using NUnit.Framework;
using Seqeron.Genomics.Annotation;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class CountTargetScanOrf8merSitesTests
{
    // targetscan_count_8mers.pl on the kit-style worked example (let-7 seed GAGGUAG): 1 site; second ORF: 3 sites.
    private static readonly string WorkedExampleOrf =
        "ATG" + string.Concat(Enumerable.Repeat("GCTGCA", 200)) + "--ctacctca.." + "CUACCUCG"
        + string.Concat(Enumerable.Repeat("GGCAAT", 46)) + "GCTAA";
    private static readonly string SecondOrf =
        "atg" + string.Concat(Enumerable.Repeat("gcugca", 300)) + "cuaccuca" + new string('-', 10)
        + string.Concat(Enumerable.Repeat("ggcaat", 50)) + "CTACctcA" + "ccc.CUACCUCA.gg"
        + string.Concat(Enumerable.Repeat("ACGT", 100)) + "taa";

    [Test]
    public void CountTargetScanOrf8merSites_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.CountTargetScanOrf8merSites("ACGU", "GAGGUAG"));
        Assert.DoesNotThrow(() => AnnotationTools.CountTargetScanOrf8merSites("", "GAGGUAG"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.CountTargetScanOrf8merSites(null!, "GAGGUAG"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.CountTargetScanOrf8merSites("ACGU", ""));
        Assert.Throws<ArgumentException>(() => AnnotationTools.CountTargetScanOrf8merSites("ACGU", null!));
        Assert.Throws<ArgumentException>(() => AnnotationTools.CountTargetScanOrf8merSites("ACGU", "GAGGUA"));
    }

    [Test]
    public void CountTargetScanOrf8merSites_Binding_EqualsCountScript()
    {
        Assert.Multiple(() =>
        {
            Assert.That(AnnotationTools.CountTargetScanOrf8merSites(WorkedExampleOrf, "GAGGUAG").Count, Is.EqualTo(1));
            Assert.That(AnnotationTools.CountTargetScanOrf8merSites(SecondOrf, "GAGGUAG").Count, Is.EqualTo(3));
            Assert.That(AnnotationTools.CountTargetScanOrf8merSites(SecondOrf, "GAGGUAG").Count,
                Is.EqualTo(MiRnaAnalyzer.CountTargetScanOrf8merSites(SecondOrf, "GAGGUAG")));
            // non-overlapping, case-insensitive (perl /gi): AUUUUUUAUUUUUUA ... = 3, not 4
            Assert.That(AnnotationTools.CountTargetScanOrf8merSites("GGAUUUUUUAUUUUUUAccauuuuuuaATTTTTTAG", "AAAAAAU").Count, Is.EqualTo(3));
            Assert.That(AnnotationTools.CountTargetScanOrf8merSites("", "GAGGUAG").Count, Is.EqualTo(0));
        });
    }
}
