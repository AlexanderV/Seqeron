using NUnit.Framework;
using Seqeron.Genomics.Annotation;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class FindCenteredSitesTests
{
    private const string Let7aSeq = "UGAGGUAGUAGGUUGUAUAGUU";
    private static readonly MiRnaDto Let7a = ToDto(MiRnaAnalyzer.CreateMiRna("let-7a", Let7aSeq));

    private static MiRnaDto ToDto(MiRnaAnalyzer.MiRna m) => new(m.Name, m.Sequence, m.SeedSequence, m.SeedStart, m.SeedEnd);

    private static string Rc(string s) => MiRnaAnalyzer.GetReverseComplement(s);

    [Test]
    public void FindCenteredSites_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.FindCenteredSites("AAAAAAAAAAAAAAAAAAAA", Let7a));
        Assert.Throws<ArgumentException>(() => AnnotationTools.FindCenteredSites("", Let7a));
        Assert.Throws<ArgumentException>(() => AnnotationTools.FindCenteredSites(null!, Let7a));
        Assert.Throws<ArgumentException>(() => AnnotationTools.FindCenteredSites("AAAA", null!));
        Assert.Throws<ArgumentException>(() => AnnotationTools.FindCenteredSites("AAAA", Let7a with { Sequence = "" }));
    }

    [Test]
    public void FindCenteredSites_Binding_12ntRun_Nt4to15()
    {
        // Mirrors Seqeron.Genomics.Tests FindCenteredSites_12ntRun_Nt4to15_Found (Shin 2010: 11-12 contiguous pairs at nt 4-15).
        string utr = "GGGGGGGGGG" + Rc(Let7aSeq.Substring(3, 12)) + "GGGGGGGGGG";
        var site = AnnotationTools.FindCenteredSites(utr, Let7a).Sites.Single();
        Assert.Multiple(() =>
        {
            Assert.That(site.Type, Is.EqualTo("Centered"));
            Assert.That(site.Start, Is.EqualTo(10));
            Assert.That(site.End, Is.EqualTo(21));
            Assert.That(site.SeedMatchLength, Is.EqualTo(12));
            Assert.That(site.Score, Is.EqualTo(0.0));
        });
    }

    [Test]
    public void FindCenteredSites_10ntRun_AndCanonicalSeedSite_Empty()
    {
        Assert.That(AnnotationTools.FindCenteredSites("AAAAAAAAAA" + Rc(Let7aSeq.Substring(4, 10)) + "AAAAAAAAAA", Let7a).Sites, Is.Empty);
        Assert.That(AnnotationTools.FindCenteredSites("AAAAAAAAAA" + Rc(Let7aSeq.Substring(1, 14)) + "AAAAAAAAAA", Let7a).Sites, Is.Empty,
            "perfect seed pairing = canonical seed site, not centered");
    }

    [Test]
    public void FindCenteredSites_DelegatesToLibrary()
    {
        string utr = "AAAAAAAAAA" + Rc(Let7aSeq.Substring(4, 11)) + "AAAAAAAAAA" + Rc(Let7aSeq.Substring(3, 12)) + "GGG";
        var direct = MiRnaAnalyzer.FindCenteredSites(utr, MiRnaAnalyzer.CreateMiRna("let-7a", Let7aSeq)).ToList();
        var wrapped = AnnotationTools.FindCenteredSites(utr, Let7a).Sites;
        Assert.That(direct, Is.Not.Empty);
        Assert.That(wrapped.Select(s => (s.Start, s.End, s.SeedMatchLength, s.TargetSequence, s.Alignment, s.FreeEnergy)),
            Is.EqualTo(direct.Select(s => (s.Start, s.End, s.SeedMatchLength, s.TargetSequence, s.Alignment, s.FreeEnergy))));
    }
}
