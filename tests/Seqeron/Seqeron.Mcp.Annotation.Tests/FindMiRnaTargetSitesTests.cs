using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class FindMiRnaTargetSitesTests
{
    // let-7a; seed pos 2-8 = GAGGUAG; seed reverse complement = CUACCUC.
    // Mirrors Seqeron.Genomics.Tests MiRnaAnalyzer_TargetPrediction_Tests (Bartel 2009 site types).
    private const string Let7aSeedRC = "CUACCUC";
    private static MiRnaDto Let7a() => AnnotationTools.CreateMiRna("let-7a", "UGAGGUAGUAGGUUGUAUAGUU").MiRna;

    [Test]
    public void FindMiRnaTargetSites_Schema_ValidatesCorrectly()
    {
        var mrna = "GGGGG" + Let7aSeedRC + "A" + "GGGGG";
        Assert.DoesNotThrow(() => AnnotationTools.FindMiRnaTargetSites(mrna, Let7a()));
        Assert.Throws<ArgumentException>(() => AnnotationTools.FindMiRnaTargetSites("", Let7a()));
        var emptyMiRna = new MiRnaDto("x", "", "", 1, 7);
        Assert.Throws<ArgumentException>(() => AnnotationTools.FindMiRnaTargetSites(mrna, emptyMiRna));
    }

    [Test]
    public void FindMiRnaTargetSites_8merSite_DetectedAndScoredHighest()
    {
        // 8mer: full seedRC + trailing A (Bartel 2009). Single highest-scoring site.
        string mrna = "GGGGG" + Let7aSeedRC + "A" + "GGGGG";
        var result = AnnotationTools.FindMiRnaTargetSites(mrna, Let7a(), minScore: 0.1);

        Assert.That(result.Sites, Has.Count.EqualTo(1));
        Assert.Multiple(() =>
        {
            Assert.That(result.Sites[0].Type, Is.EqualTo("Seed8mer"));
            Assert.That(result.Sites[0].SeedMatchLength, Is.EqualTo(8));
            Assert.That(result.Sites[0].Score, Is.GreaterThanOrEqualTo(0.9));
        });
    }

    [Test]
    public void FindMiRnaTargetSites_7merM8Site_Detected()
    {
        // 7mer-m8: full seedRC, no trailing A (trailing G).
        string mrna = "GGGGG" + Let7aSeedRC + "G" + "GGGGG";
        var result = AnnotationTools.FindMiRnaTargetSites(mrna, Let7a(), minScore: 0.1);

        Assert.That(result.Sites, Has.Count.EqualTo(1));
        Assert.Multiple(() =>
        {
            Assert.That(result.Sites[0].Type, Is.EqualTo("Seed7merM8"));
            Assert.That(result.Sites[0].SeedMatchLength, Is.EqualTo(7));
        });
    }

    [Test]
    public void FindMiRnaTargetSites_NoSeedMatch_ReturnsEmpty()
    {
        // Poly-A mRNA has no seed match for let-7a.
        var result = AnnotationTools.FindMiRnaTargetSites("AAAAAAAAAAAAAAAAAAAA", Let7a());
        Assert.That(result.Sites, Is.Empty);
    }
}
