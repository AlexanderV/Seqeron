using NUnit.Framework;
using Seqeron.Genomics.Annotation;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class ComputeTa3UtrTests
{
    private const string Let7a = "UGAGGUAGUAGGUUGUAUAGUU";

    [Test]
    public void ComputeTa3Utr_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.ComputeTa3Utr(Let7a, new[] { "ACGU" }));
        Assert.DoesNotThrow(() => AnnotationTools.ComputeTa3Utr(Let7a, Array.Empty<string>()));
        Assert.Throws<ArgumentException>(() => AnnotationTools.ComputeTa3Utr("", new[] { "ACGU" }));
        Assert.Throws<ArgumentException>(() => AnnotationTools.ComputeTa3Utr(null!, new[] { "ACGU" }));
        Assert.Throws<ArgumentException>(() => AnnotationTools.ComputeTa3Utr("UGAGGUA", new[] { "ACGU" }));
        Assert.Throws<ArgumentException>(() => AnnotationTools.ComputeTa3Utr(Let7a, null!));
    }

    [Test]
    public void ComputeTa3Utr_Binding_CountsEightmerSevenmersNotSixmers_Log10()
    {
        // Seed GAGGUAG (nt 2-8) → 7mer-m8 CUACCUC, 7mer-A1 UACCUCA, 8mer CUACCUCA, bare 6mer UACCUC (Garcia 2011 /
        // TargetScan: 8mer + 7mer-m8 + 7mer-A1 are counted, a bare 6mer is not). Four sites → TA = log10(4).
        var utrs = new[]
        {
            "GGGCUACCUCAGGG",                 // 8mer
            "GGGCUACCUCGGGGGUACCUCAGG",       // 7mer-m8 + 7mer-A1
            "GGGUACCUCGGG",                   // bare 6mer: not counted
            "gggctacctcaggg",                 // DNA / lower case: second 8mer
        };
        var r = AnnotationTools.ComputeTa3Utr(Let7a, utrs);
        Assert.Multiple(() =>
        {
            Assert.That(r.SiteCount, Is.EqualTo(4));
            Assert.That(r.Ta, Is.EqualTo(Math.Log10(4)).Within(1e-15));
            Assert.That(r.UtrCount, Is.EqualTo(4));
        });
    }

    [Test]
    public void ComputeTa3Utr_NoSites_TaIsZero()
    {
        var r = AnnotationTools.ComputeTa3Utr(Let7a, new[] { "GGGGGGGGGG", "" });
        Assert.That(r.SiteCount, Is.EqualTo(0));
        Assert.That(r.Ta, Is.EqualTo(0.0));
    }

    [Test]
    public void ComputeTa3Utr_DelegatesToLibrary()
    {
        var utrs = new[] { "GGGCUACCUCAGGGCUACCUCAGG", "CUACCUCACUACCUCAAA", "AUGCUACCUAGGG" };
        var miRna = MiRnaAnalyzer.CreateMiRna("x", Let7a);
        var r = AnnotationTools.ComputeTa3Utr(Let7a, utrs);
        Assert.That(r.SiteCount, Is.EqualTo(MiRnaAnalyzer.CountSeedSites3Utr(miRna, utrs)));
        Assert.That(r.Ta, Is.EqualTo(MiRnaAnalyzer.ComputeTa3Utr(miRna, utrs)));
    }
}
