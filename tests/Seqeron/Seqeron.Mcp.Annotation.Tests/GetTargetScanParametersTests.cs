using NUnit.Framework;
using Seqeron.Genomics.Annotation;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class GetTargetScanParametersTests
{
    // Hand-copied source rows (MiRnaAnalyzer tests): TA_SPS_by_seed_region.txt GAGGUAG: SPS 8mer/7mer-m8 -9.250, 7mer-1a/6mer -6.720, TA 3.393;
    // 7mer-m8_PCT_parameters.txt let-7 row.
    private static readonly PctSigmoidParametersDto Let7M8 =
        new(-0.221367319203409, 1.24081845475817, 1.17370692612852, 3.41806979609985);

    [Test]
    public void GetTargetScanParameters_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.GetTargetScanParameters("GAGGUAG"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.GetTargetScanParameters(""));
        Assert.Throws<ArgumentException>(() => AnnotationTools.GetTargetScanParameters(null!));
        Assert.Throws<ArgumentException>(() => AnnotationTools.GetTargetScanParameters("GAGGUAG", "Bogus"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.GetTargetScanParameters("GAGGUAG", "2"), "numeric enum");
    }

    [Test]
    public void GetTargetScanParameters_Binding_Let7_SourceRows()
    {
        var m8 = AnnotationTools.GetTargetScanParameters("GAGGUAG", "Seed7merM8");
        var a1 = AnnotationTools.GetTargetScanParameters("gaggtag", "Seed7merA1");
        var none = AnnotationTools.GetTargetScanParameters("GAGGUAG");
        Assert.Multiple(() =>
        {
            Assert.That(m8.IsPctSeedFamily, Is.True);
            Assert.That(m8.PctParameters, Is.EqualTo(Let7M8));
            Assert.That(m8.HasTaSps, Is.True);
            Assert.That(m8.Ta, Is.EqualTo(3.393));
            Assert.That(m8.Sps, Is.EqualTo(-9.25));
            Assert.That(a1.PctParameters, Is.Not.Null);
            Assert.That(a1.Sps, Is.EqualTo(-6.72), "7mer-A1 SPS column; lower case / T as U");
            Assert.That(none.IsPctSeedFamily, Is.True);
            Assert.That(none.PctParameters, Is.Null);
            Assert.That(none.Ta, Is.EqualTo(3.393));
            Assert.That(none.Sps, Is.Null);
        });
        Assert.That(MiRnaAnalyzer.TryGetTargetScanPctParameters("GAGGUAG", MiRnaAnalyzer.TargetSiteType.Seed7merM8, out var p), Is.True);
        Assert.That(m8.PctParameters, Is.EqualTo(new PctSigmoidParametersDto(p.B0, p.B1, p.B2, p.B3)));
    }

    [Test]
    public void GetTargetScanParameters_TaSpsTable_FirstAndLastRows()
    {
        var first = AnnotationTools.GetTargetScanParameters("AAAAAAA", "Seed8mer");
        var last = AnnotationTools.GetTargetScanParameters("UUUUUUU", "Seed6mer");
        Assert.Multiple(() =>
        {
            Assert.That((first.Ta, first.Sps), Is.EqualTo((4.371, -0.59)));
            Assert.That((last.Ta, last.Sps), Is.EqualTo((4.075, 0.34)));
        });
    }

    [Test]
    public void GetTargetScanParameters_Misses()
    {
        var nonFamily = AnnotationTools.GetTargetScanParameters("ACACACU", "Seed8mer");
        var bad = AnnotationTools.GetTargetScanParameters("GAGGNAG", "Seed8mer");
        var sixMer = AnnotationTools.GetTargetScanParameters("GAGGUAG", "Seed6mer");
        var offset = AnnotationTools.GetTargetScanParameters("GAGGUAG", "Offset6mer");
        Assert.Multiple(() =>
        {
            Assert.That(nonFamily.IsPctSeedFamily, Is.False);
            Assert.That(nonFamily.PctParameters, Is.Null);
            Assert.That(nonFamily.HasTaSps, Is.True, "any ACGU 7-mer has a TA/SPS row");
            Assert.That(bad.HasTaSps, Is.False);
            Assert.That(bad.Ta, Is.Null);
            Assert.That(bad.Sps, Is.Null);
            Assert.That(sixMer.PctParameters, Is.Null, "no 6mer PCT row");
            Assert.That(sixMer.Sps, Is.EqualTo(-6.72));
            Assert.That(offset.PctParameters, Is.Null);
            Assert.That(offset.Sps, Is.Null, "no SPS column for Offset6mer");
            Assert.That(offset.Ta, Is.EqualTo(3.393));
        });
    }
}
