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

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: MIRNA-PRECURSOR-001 — Pre-miRNA hairpin detection (MiRNA)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 76.
    // Spec: tests/TestSpecs/MIRNA-PRECURSOR-001.md (canonical FindPreMiRnaHairpins).
    // Dimensions: precLen(3) × minStem(3) × maxLoop(3). Grid 3×3×3 = 27.
    //
    // Model (Bartel 2004): a pre-miRNA is a ≥55-nt hairpin with a long (~≥18 bp) stem and a small
    // (3–25 nt) terminal loop. The stem/loop thresholds are FIXED internal constants here, so the
    // minStem/maxLoop axes vary the PLANTED hairpin's stem/loop against them; precLen is the
    // search length window (min/maxHairpinLength).
    //
    // The combinatorial point: planted stem length, loop size and the length window jointly
    // determine detection — a hairpin is found exactly when stem ≥ 18, loop ∈ [3,25], and its
    // length is ≥ 55 and within the search window.
    // ═══════════════════════════════════════════════════════════════════════

    private static readonly int[] StemValues = { 16, 25, 32 };
    private static readonly int[] LoopValues = { 6, 20, 30 };
    private static readonly (int Min, int Max)[] LenWindows = { (40, 200), (40, 70), (40, 50) };

    [Test, Combinatorial]
    public void MiRnaPrecursor_DetectsHairpinAgainstThresholds(
        [Values(0, 1, 2)] int precLenIdx,
        [Values(0, 1, 2)] int minStemIdx,
        [Values(0, 1, 2)] int maxLoopIdx)
    {
        int stem = StemValues[minStemIdx];
        int loop = LoopValues[maxLoopIdx];
        var (minH, maxH) = LenWindows[precLenIdx];

        string hairpin = new string('G', stem) + new string('A', loop) + new string('C', stem);
        int innerLen = 2 * stem + loop;
        string seq = new string('U', 10) + hairpin + new string('U', 10);

        var found = MiRnaAnalyzer.FindPreMiRnaHairpins(seq, minH, maxH, matureLength: 22).ToList();

        bool expected = stem >= 18 && loop is >= 3 and <= 25 && innerLen >= 55 && innerLen >= minH && innerLen <= maxH;
        found.Any(p => p.Sequence == hairpin)
            .Should().Be(expected, "the planted hairpin is found iff stem/loop/length all qualify");

        foreach (var p in found)
        {
            p.Sequence.Length.Should().Be(p.End - p.Start + 1, "coordinates match the sequence span");
            p.Structure.Length.Should().Be(p.Sequence.Length, "dot-bracket spans the hairpin");
            p.MatureSequence.Should().NotBeEmpty("a mature arm is extracted");
        }
    }

    /// <summary>
    /// Interaction witness: each requirement independently rejects a hairpin — a short stem, an
    /// oversized loop, or a length below 55 nt all prevent detection.
    /// </summary>
    [Test]
    public void MiRnaPrecursor_EachRequirement_GatesDetection()
    {
        bool Found(int stem, int loop)
        {
            string hp = new string('G', stem) + new string('A', loop) + new string('C', stem);
            return MiRnaAnalyzer.FindPreMiRnaHairpins(new string('U', 10) + hp + new string('U', 10), 40, 200, 22)
                .Any(p => p.Sequence == hp);
        }

        Found(25, 6).Should().BeTrue("a 25-bp stem, 6-nt loop, 56-nt hairpin qualifies");
        Found(16, 6).Should().BeFalse("stem < 18 bp is rejected");
        Found(25, 30).Should().BeFalse("loop > 25 nt is rejected");
        Found(20, 6).Should().BeFalse("a 46-nt hairpin is below the 55-nt floor");
    }
}
