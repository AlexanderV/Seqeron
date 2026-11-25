using NUnit.Framework;
using Seqeron.Genomics.Annotation;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class ScoreTargetSiteGrimsonTests
{
    private const string Let7a = "UGAGGUAGUAGGUUGUAUAGUU";

    [Test]
    public void ScoreTargetSiteGrimson_Schema_ValidatesCorrectly()
    {
        string utr = new string('G', 100) + "CUACCUC" + "A" + new string('G', 100);
        Assert.DoesNotThrow(() => AnnotationTools.ScoreTargetSiteGrimson(utr, Let7a, 100, 107, "Seed8mer"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.ScoreTargetSiteGrimson("", Let7a, 0, 1, "Seed8mer"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.ScoreTargetSiteGrimson(null!, Let7a, 0, 1, "Seed8mer"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.ScoreTargetSiteGrimson(utr, "", 100, 107, "Seed8mer"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.ScoreTargetSiteGrimson(utr, null!, 100, 107, "Seed8mer"));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnnotationTools.ScoreTargetSiteGrimson(utr, Let7a, -1, 107, "Seed8mer"));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnnotationTools.ScoreTargetSiteGrimson(utr, Let7a, 100, utr.Length, "Seed8mer"));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnnotationTools.ScoreTargetSiteGrimson(utr, Let7a, 107, 100, "Seed8mer"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.ScoreTargetSiteGrimson(utr, Let7a, 100, 107, "Bogus"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.ScoreTargetSiteGrimson(utr, Let7a, 100, 107, "99"));
    }

    [Test]
    public void ScoreTargetSiteGrimson_Binding_8merPolyG_HandDerived()
    {
        // Mirrors Seqeron.Genomics.Tests WorkedExample_8merPolyG_HandDerived (TargetScan 5.0 regression constants).
        string utr = new string('G', 100) + "CUACCUC" + "A" + new string('G', 100);
        var r = AnnotationTools.ScoreTargetSiteGrimson(utr, Let7a, 100, 107, "Seed8mer");

        Assert.That(r.Scored, Is.True);
        double h30 = 0, down = 0;
        for (int i = 1; i <= 30; i++) h30 += 1.0 / i;
        for (int i = 0; i < 30; i++) down += 1.0 / (i + 2);
        double fraction = 0.5 / (h30 + down);
        var g = r.Score!;
        Assert.Multiple(() =>
        {
            Assert.That(g.Type, Is.EqualTo("Seed8mer"));
            Assert.That(g.SiteTypeContribution, Is.EqualTo(-0.310));
            Assert.That(g.LocalAuFraction, Is.EqualTo(fraction).Within(1e-12));
            Assert.That(g.LocalAuContribution, Is.EqualTo(Math.Round(-0.64 * fraction + 0.365, 3)).Within(1e-12));
            Assert.That(g.PositionContribution, Is.EqualTo(Math.Round(100 * 0.000172 - 0.07, 3)).Within(1e-12));
            Assert.That(g.ContextScore!.Value, Is.EqualTo(
                g.SiteTypeContribution + g.LocalAuContribution + g.PositionContribution + g.ThreePrimePairingContribution).Within(1e-12));
        });
    }

    [Test]
    public void ScoreTargetSiteGrimson_NearStopCodon_ContextScoreNull_ContributionsKept()
    {
        string utr = "GGGGGGG" + "CUACCUC" + "A" + new string('G', 60);   // 1-based start 8 < 15 -> too_close
        var r = AnnotationTools.ScoreTargetSiteGrimson(utr, Let7a, 7, 14, "seed8mer");
        Assert.That(r.Scored, Is.True);
        Assert.That(r.Score!.ContextScore, Is.Null);
        Assert.That(r.Score.SiteTypeContribution, Is.EqualTo(-0.310));
    }

    [TestCase("Seed6mer")]
    [TestCase("Offset6mer")]
    [TestCase("Centered")]
    public void ScoreTargetSiteGrimson_TypesWithoutCoefficients_NotScored(string type)
    {
        const string utr = "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGCUACCUGGGGGGGGGGGGGGGGGGGGGGGGGG";
        var r = AnnotationTools.ScoreTargetSiteGrimson(utr, Let7a, 40, 45, type);
        Assert.That(r.Scored, Is.False);
        Assert.That(r.Score, Is.Null);
    }

    [TestCase("Seed8mer", 100, 107)]
    [TestCase("Seed7merM8", 100, 106)]
    [TestCase("Seed7merA1", 100, 106)]
    public void ScoreTargetSiteGrimson_DelegatesToLibrary(string type, int start, int end)
    {
        string utr = new string('A', 60) + new string('G', 40) + "CUACCUC" + "A" + new string('U', 70);
        var direct = MiRnaAnalyzer.ScoreTargetSiteGrimson(utr, Let7a, start, end, Enum.Parse<MiRnaAnalyzer.TargetSiteType>(type))!.Value;
        var g = AnnotationTools.ScoreTargetSiteGrimson(utr, Let7a, start, end, type).Score!;
        Assert.Multiple(() =>
        {
            Assert.That(g.SiteTypeContribution, Is.EqualTo(direct.SiteTypeContribution));
            Assert.That(g.ThreePrimePairingContribution, Is.EqualTo(direct.ThreePrimePairingContribution));
            Assert.That(g.LocalAuContribution, Is.EqualTo(direct.LocalAuContribution));
            Assert.That(g.PositionContribution, Is.EqualTo(direct.PositionContribution));
            Assert.That(g.ThreePrimePairingScore, Is.EqualTo(direct.ThreePrimePairingScore));
            Assert.That(g.LocalAuFraction, Is.EqualTo(direct.LocalAuFraction));
            Assert.That(g.ContextScore, Is.EqualTo(direct.ContextScore));
        });
    }
}
