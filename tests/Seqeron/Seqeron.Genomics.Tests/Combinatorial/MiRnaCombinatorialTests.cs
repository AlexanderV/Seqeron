namespace Seqeron.Genomics.Tests.Combinatorial;

/// <summary>
/// Combinatorial (pairwise / full-grid) tests for the MiRNA area.
///
/// See <see cref="CompositionCombinatorialTests"/> for the rationale of
/// combinatorial testing. Each grid cell carries a real business assertion;
/// small grids use the exhaustive <c>[Combinatorial]</c> product.
/// — docs/checklists/09_COMBINATORIAL_TESTING.md §Description.
/// </summary>
[TestFixture]
[Category("Combinatorial")]
[Category("MiRNA")]
public class MiRnaCombinatorialTests
{
    // miRNA whose seed (positions 2–8) is CGUACGU; the 6-mer target core is then CGUACG.
    private const string MiRnaSeq = "ACGUACGUACGUACGUACGUAC";
    private const int SiteOffset = 30;

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: MIRNA-TARGET-001 — miRNA target-site prediction (MiRNA)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 75.
    // Spec: tests/TestSpecs/MIRNA-TARGET-001.md (canonical FindTargetSites + context analysis).
    // Dimensions: seedType(3) × utrLen(3) × scoringMethod(2). Grid 3×3×2 = 18.
    //
    // Model (Bartel 2009; TargetScan): canonical miRNA sites are seed-complementary matches in the
    // 3′UTR — 8mer (seed match 2–8 + A1), 7mer-m8 (2–8) and 6mer (2–7). FindTargetSites classifies
    // each. The scoring axis covers seed-match scoring vs context augmentation (AnalyzeTargetContext).
    //
    // The combinatorial point: site type, UTR length and scoring method interact — the planted
    // canonical site is detected and classified correctly at every UTR length, with a positive
    // seed score and a valid context score.
    // ═══════════════════════════════════════════════════════════════════════

    public enum SeedType { Site8mer, Site7merM8, Site6mer }
    public enum MiRnaScoring { SeedMatch, ContextAugmented }

    private static string PlantSite(SeedType type) => type switch
    {
        // seedRC = ACGUACG (RC of seed CGUACGU); core = CGUACG.
        SeedType.Site8mer => "ACGUACGA",   // pos8 match + A1 ⇒ 8mer
        SeedType.Site7merM8 => "ACGUACGC", // pos8 match, no A1 ⇒ 7mer-m8
        _ => "UCGUACGU",                   // core only, no pos8 / no A1 ⇒ 6mer
    };

    private static (MiRnaAnalyzer.TargetSiteType Type, int Len) Expected(SeedType type) => type switch
    {
        SeedType.Site8mer => (MiRnaAnalyzer.TargetSiteType.Seed8mer, 8),
        SeedType.Site7merM8 => (MiRnaAnalyzer.TargetSiteType.Seed7merM8, 7),
        _ => (MiRnaAnalyzer.TargetSiteType.Seed6mer, 6),
    };

    [Test, Combinatorial]
    public void MiRnaTarget_DetectsAndClassifiesCanonicalSite(
        [Values(SeedType.Site8mer, SeedType.Site7merM8, SeedType.Site6mer)] SeedType seedType,
        [Values(60, 120, 240)] int utrLen,
        [Values(MiRnaScoring.SeedMatch, MiRnaScoring.ContextAugmented)] MiRnaScoring scoring)
    {
        var miRna = MiRnaAnalyzer.CreateMiRna("miR-x", MiRnaSeq);
        string site = PlantSite(seedType);
        string utr = new string('U', SiteOffset) + site + new string('U', utrLen - SiteOffset - site.Length);

        var sites = MiRnaAnalyzer.FindTargetSites(utr, miRna, minScore: 0.0).ToList();
        var (expectedType, expectedLen) = Expected(seedType);

        sites.Should().Contain(s => s.Type == expectedType, "the planted canonical site is classified correctly");
        var found = sites.First(s => s.Type == expectedType);

        if (scoring == MiRnaScoring.SeedMatch)
        {
            found.Score.Should().BeGreaterThan(0, "a canonical seed match scores positively");
            found.SeedMatchLength.Should().Be(expectedLen);
        }
        else
        {
            var ctx = MiRnaAnalyzer.AnalyzeTargetContext(utr, found.Start, found.End);
            ctx.AuContent.Should().BeInRange(0.0, 1.0);
            ctx.ContextScore.Should().BeInRange(0.0, 1.0);
        }
    }

    /// <summary>
    /// Interaction witness: the three canonical site types are ordered by stringency — an 8mer
    /// scores at least as high as a 7mer-m8, which scores at least as high as a 6mer.
    /// </summary>
    [Test]
    public void MiRnaTarget_SiteHierarchy_ScoresOrdered()
    {
        var miRna = MiRnaAnalyzer.CreateMiRna("miR-x", MiRnaSeq);
        double Score(SeedType t)
        {
            string utr = new string('U', SiteOffset) + PlantSite(t) + new string('U', 80);
            var (type, _) = Expected(t);
            return MiRnaAnalyzer.FindTargetSites(utr, miRna, 0.0).First(s => s.Type == type).Score;
        }

        Score(SeedType.Site8mer).Should().BeGreaterThanOrEqualTo(Score(SeedType.Site7merM8));
        Score(SeedType.Site7merM8).Should().BeGreaterThanOrEqualTo(Score(SeedType.Site6mer));
    }

    /// <summary>
    /// Interaction witness: a UTR with no seed-complementary core yields no canonical target site.
    /// </summary>
    [Test]
    public void MiRnaTarget_NoSeedMatch_NoSite()
    {
        var miRna = MiRnaAnalyzer.CreateMiRna("miR-x", MiRnaSeq);
        MiRnaAnalyzer.FindTargetSites(new string('U', 100), miRna, 0.0)
            .Should().BeEmpty("a poly-U UTR has no seed-complementary site");
    }
}
