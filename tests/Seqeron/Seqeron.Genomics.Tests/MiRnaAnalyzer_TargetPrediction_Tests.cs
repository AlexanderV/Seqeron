using NUnit.Framework;
using Seqeron.Genomics.Annotation;
using static Seqeron.Genomics.Annotation.MiRnaAnalyzer;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// MIRNA-TARGET-001: Target Site Prediction
/// Canonical test file for FindTargetSites and target scoring.
/// Evidence: Bartel (2009), Lewis et al. (2005), Grimson et al. (2007), TargetScan.
/// </summary>
[TestFixture]
public class MiRnaAnalyzer_TargetPrediction_Tests
{
    #region Reference Data

    // hsa-let-7a-5p: UGAGGUAGUAGGUUGUAUAGUU
    // Seed (pos 2-8): GAGGUAG
    // Seed RC: CUACCUC
    private static readonly MiRna Let7a = CreateMiRna("let-7a", "UGAGGUAGUAGGUUGUAUAGUU");
    private const string Let7aSeedRC = "CUACCUC";

    // hsa-miR-21-5p: UAGCUUAUCAGACUGAUGUUGA
    // Seed (pos 2-8): AGCUUAU
    // Seed RC: AUAAGCU
    private static readonly MiRna MiR21 = CreateMiRna("miR-21", "UAGCUUAUCAGACUGAUGUUGA");
    private const string MiR21SeedRC = "AUAAGCU";

    #endregion

    #region M-001: 8mer site detected — Bartel (2009)

    [Test]
    public void FindTargetSites_8merSite_DetectedAndScoredHighest()
    {
        // 8mer: full seedRC (CUACCUC) + trailing A → CUACCUCA
        // Evidence: Bartel (2009) — 8mer = positions 2-8 match + A opposite position 1
        // Layout on mRNA (5'→3'): [seedRC] [A] where seedRC = RC of miRNA pos 2-8
        string mrna = "GGGGG" + Let7aSeedRC + "A" + "GGGGG";

        var sites = FindTargetSites(mrna, Let7a, minScore: 0.1).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(sites, Has.Count.EqualTo(1));
            Assert.That(sites[0].Type, Is.EqualTo(TargetSiteType.Seed8mer));
            Assert.That(sites[0].SeedMatchLength, Is.EqualTo(8));
            Assert.That(sites[0].Score, Is.GreaterThanOrEqualTo(0.9));
        });
    }

    #endregion

    #region M-002: 7mer-m8 site detected — Bartel (2009)

    [Test]
    public void FindTargetSites_7merM8Site_Detected()
    {
        // 7mer-m8: full seed RC (CUACCUC) but no trailing A → trailing G
        // Evidence: Bartel (2009) — 7mer-m8 = positions 2-8 match, no A at pos 1
        string mrna = "GGGGG" + Let7aSeedRC + "G" + "GGGGG";

        var sites = FindTargetSites(mrna, Let7a, minScore: 0.1).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(sites, Has.Count.EqualTo(1));
            Assert.That(sites[0].Type, Is.EqualTo(TargetSiteType.Seed7merM8));
            Assert.That(sites[0].SeedMatchLength, Is.EqualTo(7));
        });
    }

    #endregion

    #region M-003: 7mer-A1 site detected — Bartel (2009)

    [Test]
    public void FindTargetSites_7merA1Site_Detected()
    {
        // 7mer-A1: 6mer core (seedRC[1:7] = UACCUC, RC of pos 2-7) + trailing A
        // No upstream seedRC[0]='C' → prevents upgrade to 8mer
        // Evidence: Bartel (2009) — 7mer-A1 = positions 2-7 + A opposite position 1
        string sixmerCore = Let7aSeedRC[1..]; // UACCUC (RC of miRNA positions 2-7)
        string mrna = "GGGGG" + sixmerCore + "A" + "GGGGG";

        var sites = FindTargetSites(mrna, Let7a, minScore: 0.1).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(sites, Has.Count.EqualTo(1));
            Assert.That(sites[0].Type, Is.EqualTo(TargetSiteType.Seed7merA1));
            Assert.That(sites[0].SeedMatchLength, Is.EqualTo(7));
        });
    }

    #endregion

    #region M-004: 6mer site detected — Bartel (2009)

    [Test]
    public void FindTargetSites_6merSite_Detected()
    {
        // 6mer: 6mer core = seedRC[1:7] = UACCUC (RC of pos 2-7), no trailing A, no upstream seedRC[0]
        // Evidence: Bartel (2009) — 6mer = positions 2-7 match only
        string sixmerCore = Let7aSeedRC[1..]; // UACCUC (RC of miRNA positions 2-7)
        string mrna = "GGGGG" + sixmerCore + "G" + "GGGGG";

        var sites = FindTargetSites(mrna, Let7a, minScore: 0.01).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(sites, Has.Count.EqualTo(1));
            Assert.That(sites[0].Type, Is.EqualTo(TargetSiteType.Seed6mer));
            Assert.That(sites[0].SeedMatchLength, Is.EqualTo(6));
        });
    }

    #endregion

    #region M-005: Score monotonicity — Grimson (2007)

    [Test]
    public void FindTargetSites_ScoreMonotonicity_8mer_GT_7merM8_GT_7merA1_GT_6mer()
    {
        // Evidence: Grimson (2007) — efficacy hierarchy: 8mer > 7mer-m8 > 7mer-A1 > 6mer
        // Grimson weights: 8mer=0.310, 7mer-m8=0.161, 7mer-A1=0.099
        string sixmerCore = Let7aSeedRC[1..]; // UACCUC (RC of positions 2-7)

        string mrna8mer = "GGGGG" + Let7aSeedRC + "A" + "GGGGG";
        string mrna7merM8 = "GGGGG" + Let7aSeedRC + "G" + "GGGGG";
        string mrna7merA1 = "GGGGG" + sixmerCore + "A" + "GGGGG";
        string mrna6mer = "GGGGG" + sixmerCore + "G" + "GGGGG";

        double score8mer = FindTargetSites(mrna8mer, Let7a, 0.01).First().Score;
        double score7merM8 = FindTargetSites(mrna7merM8, Let7a, 0.01).First().Score;
        double score7merA1 = FindTargetSites(mrna7merA1, Let7a, 0.01).First().Score;
        double score6mer = FindTargetSites(mrna6mer, Let7a, 0.01).First().Score;

        Assert.Multiple(() =>
        {
            Assert.That(score8mer, Is.GreaterThan(score7merM8), "8mer > 7mer-m8");
            Assert.That(score7merM8, Is.GreaterThan(score7merA1), "7mer-m8 > 7mer-A1");
            Assert.That(score7merA1, Is.GreaterThan(score6mer), "7mer-A1 > 6mer");
        });
    }

    #endregion

    #region M-006: Empty/null inputs — Defensive contract

    [Test]
    public void FindTargetSites_EmptyMrna_ReturnsEmpty()
    {
        Assert.That(FindTargetSites("", Let7a).ToList(), Is.Empty);
    }

    [Test]
    public void FindTargetSites_EmptyMiRnaSequence_ReturnsEmpty()
    {
        var emptyMiRna = new MiRna();
        Assert.That(FindTargetSites("AUGCAUGCAUGC", emptyMiRna).ToList(), Is.Empty);
    }

    #endregion

    #region M-007: No match returns empty — Trivial correctness

    [Test]
    public void FindTargetSites_NoSeedMatch_ReturnsEmpty()
    {
        // All G's — no seed RC present
        string mrna = new string('G', 30);
        Assert.That(FindTargetSites(mrna, Let7a, 0.1).ToList(), Is.Empty);
    }

    #endregion

    #region M-008: Multiple sites found — Bartel (2009)

    [Test]
    public void FindTargetSites_MultipleSites_AllFound()
    {
        // Evidence: each seed match site functions independently — Bartel (2009)
        // Three 7mer-m8 sites (full seedRC + non-A) separated by GGG spacers
        string site = Let7aSeedRC + "G"; // 7mer-m8 site (full seedRC + non-A)
        string mrna = "GGG" + site + "GGG" + site + "GGG" + site + "GGG";

        var sites = FindTargetSites(mrna, Let7a, minScore: 0.01).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(sites, Has.Count.EqualTo(3));
            Assert.That(sites.All(s => s.Type == TargetSiteType.Seed7merM8), Is.True,
                "All three sites should be classified as 7mer-m8");
        });
    }

    #endregion

    #region M-009: Score range [0, 1] — Implementation contract

    [TestCase("GGGGG" + "CUACCUC" + "A" + "GGGGG")] // 8mer for let-7a
    [TestCase("GGGGG" + "CUACCUC" + "G" + "GGGGG")] // 7mer-m8 for let-7a
    [TestCase("GGGGG" + "UACCUC" + "A" + "GGGGG")]  // 7mer-A1 for let-7a
    [TestCase("GGGGG" + "UACCUC" + "G" + "GGGGG")]  // 6mer for let-7a
    [TestCase("GGGGG" + "CUACCU" + "G" + "GGGGG")]  // offset 6mer for let-7a
    public void FindTargetSites_AllSites_ScoreInRange(string mrna)
    {
        var sites = FindTargetSites(mrna, Let7a, minScore: 0.0).ToList();

        foreach (var site in sites)
        {
            Assert.That(site.Score, Is.InRange(0.0, 1.0),
                $"Site type {site.Type} score out of range");
        }
    }

    #endregion

    #region M-010: minScore filtering — API contract

    [Test]
    public void FindTargetSites_HighMinScore_FiltersLowScoringSites()
    {
        // 6mer (Grimson base score ~0.15) should be filtered with high threshold
        string sixmerCore = Let7aSeedRC[1..]; // UACCUC (RC of positions 2-7)
        string mrna = "GGGGG" + sixmerCore + "G" + "GGGGG";

        var withLowThreshold = FindTargetSites(mrna, Let7a, minScore: 0.01).ToList();
        var withHighThreshold = FindTargetSites(mrna, Let7a, minScore: 0.99).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(withLowThreshold, Is.Not.Empty, "Low threshold should find 6mer");
            Assert.That(withHighThreshold, Is.Empty, "High threshold should filter 6mer");
        });
    }

    #endregion

    #region M-011: AlignMiRnaToTarget — perfect complement — Watson-Crick

    [Test]
    public void AlignMiRnaToTarget_PerfectComplement_AllMatches()
    {
        // A pairs with U in antiparallel orientation
        var duplex = AlignMiRnaToTarget("AAAA", "UUUU");

        Assert.Multiple(() =>
        {
            Assert.That(duplex.Matches, Is.EqualTo(4));
            Assert.That(duplex.Mismatches, Is.EqualTo(0));
            Assert.That(duplex.GUWobbles, Is.EqualTo(0));
        });
    }

    #endregion

    #region M-012: AlignMiRnaToTarget — G:U wobble — Crick (1966)

    [Test]
    public void AlignMiRnaToTarget_GUWobblePairs_Detected()
    {
        // G pairs with U as wobble
        var duplex = AlignMiRnaToTarget("GGGG", "UUUU");

        Assert.Multiple(() =>
        {
            Assert.That(duplex.GUWobbles, Is.EqualTo(4));
            Assert.That(duplex.Matches, Is.EqualTo(0));
        });
    }

    #endregion

    #region M-013: AlignMiRnaToTarget — mismatches — Trivial

    [Test]
    public void AlignMiRnaToTarget_SameBases_AllMismatches()
    {
        // A does not pair with A
        var duplex = AlignMiRnaToTarget("AAAA", "AAAA");

        Assert.That(duplex.Mismatches, Is.EqualTo(4));
    }

    #endregion

    #region M-014: AlignMiRnaToTarget — empty input — Defensive

    [Test]
    public void AlignMiRnaToTarget_EmptyMiRna_ReturnsEmptyDuplex()
    {
        var duplex = AlignMiRnaToTarget("", "AAAA");

        Assert.Multiple(() =>
        {
            Assert.That(duplex.Matches, Is.EqualTo(0));
            Assert.That(duplex.MiRnaSequence, Is.Empty);
        });
    }

    [Test]
    public void AlignMiRnaToTarget_EmptyTarget_ReturnsEmptyDuplex()
    {
        var duplex = AlignMiRnaToTarget("AAAA", "");

        Assert.Multiple(() =>
        {
            Assert.That(duplex.Matches, Is.EqualTo(0));
            Assert.That(duplex.TargetSequence, Is.Empty);
        });
    }

    #endregion

    #region M-015: Free energy negative for paired duplex — Thermodynamics

    [Test]
    public void AlignMiRnaToTarget_WellPairedDuplex_NegativeFreeEnergy()
    {
        // Evidence: thermodynamic principle — stable duplexes have negative ΔG
        string mirna = "UGAGGUAGUAGGUUGUAUAGUU";
        string target = GetReverseComplement(mirna);

        var duplex = AlignMiRnaToTarget(mirna, target);

        Assert.That(duplex.FreeEnergy, Is.LessThan(0),
            "Stable duplex should have negative free energy");
    }

    #endregion

    #region M-016: Target site includes alignment string — API contract

    [Test]
    public void FindTargetSites_FoundSite_HasNonEmptyAlignment()
    {
        string mrna = "GGGGG" + Let7aSeedRC + "A" + "GGGGG";

        var sites = FindTargetSites(mrna, Let7a, minScore: 0.1).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(sites, Has.Count.EqualTo(1));
            Assert.That(sites[0].Alignment, Is.Not.Empty);
        });
    }

    #endregion

    #region M-017: DNA input (T) handled as RNA (U) — Implementation design

    [Test]
    public void FindTargetSites_DnaInput_ConvertedToRnaAndMatched()
    {
        // Same site but with T instead of U
        string mrnaRna = "GGGGG" + Let7aSeedRC + "A" + "GGGGG";
        string mrnaDna = mrnaRna.Replace('U', 'T');

        var sitesRna = FindTargetSites(mrnaRna, Let7a, minScore: 0.1).ToList();
        var sitesDna = FindTargetSites(mrnaDna, Let7a, minScore: 0.1).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(sitesDna, Has.Count.EqualTo(sitesRna.Count));
            Assert.That(sitesDna[0].Type, Is.EqualTo(sitesRna[0].Type));
            Assert.That(sitesDna[0].Score, Is.EqualTo(sitesRna[0].Score).Within(0.001));
        });
    }

    #endregion

    #region S-001: Offset 6mer detected

    [Test]
    public void FindTargetSites_Offset6merSite_Detected()
    {
        // Offset 6mer: match to positions 3-8 of miRNA = seedRC[0..6] = CUACCU
        // Must NOT have seedRC[6] match at position +6 (would make it part of full seedRC)
        string offset6 = Let7aSeedRC[..6]; // CUACCU (RC of miRNA positions 3-8)
        string mrna = "GGGGG" + offset6 + "G" + "GGGGG";

        var sites = FindTargetSites(mrna, Let7a, minScore: 0.01).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(sites, Has.Count.EqualTo(1));
            Assert.That(sites[0].Type, Is.EqualTo(TargetSiteType.Offset6mer));
            Assert.That(sites[0].SeedMatchLength, Is.EqualTo(6));
        });
    }

    #endregion

    #region S-002: TargetSite record fields populated

    [Test]
    public void FindTargetSites_FoundSite_AllFieldsPopulated()
    {
        string mrna = "GGGGG" + Let7aSeedRC + "A" + "GGGGG";

        var sites = FindTargetSites(mrna, Let7a, minScore: 0.1).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(sites, Has.Count.EqualTo(1));
            var site = sites[0];
            Assert.That(site.Start, Is.GreaterThanOrEqualTo(0));
            Assert.That(site.End, Is.GreaterThan(site.Start));
            Assert.That(site.TargetSequence, Is.Not.Empty);
            Assert.That(site.MiRnaName, Is.EqualTo("let-7a"));
            Assert.That(site.Score, Is.InRange(0.0, 1.0));
            Assert.That(site.FreeEnergy, Is.Not.EqualTo(0.0));
            Assert.That(site.Alignment, Is.Not.Empty);
        });
    }

    #endregion

    #region S-003: Real miRNA integration test

    [Test]
    public void FindTargetSites_Let7a_RealTargetSequence_FindsSite()
    {
        // Integration: Real let-7a against a constructed 3'UTR-like sequence
        // mRNA contains CUACCUCA at positions 14-21 = 8mer site for let-7a
        string mrna = "AUGGCUAAAGCUUUCUACCUCAGCUUAACCC";

        var sites = FindTargetSites(mrna, Let7a, minScore: 0.1).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(sites, Has.Count.EqualTo(1));
            Assert.That(sites[0].Type, Is.EqualTo(TargetSiteType.Seed8mer));
        });
    }

    [Test]
    public void FindTargetSites_MiR21_8merTarget_FindsSite()
    {
        // hsa-miR-21 seed RC = AUAAGCU; 8mer = AUAAGCUA
        string mrna = "GGGGG" + MiR21SeedRC + "A" + "GGGGG";

        var sites = FindTargetSites(mrna, MiR21, minScore: 0.1).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(sites, Has.Count.EqualTo(1));
            Assert.That(sites[0].Type, Is.EqualTo(TargetSiteType.Seed8mer));
        });
    }

    #endregion
}
