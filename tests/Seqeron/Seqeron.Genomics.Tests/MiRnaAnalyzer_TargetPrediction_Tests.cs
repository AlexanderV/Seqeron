using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using Seqeron.Genomics.Annotation;
using Seqeron.Genomics.Phylogenetics;
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
            // FreeEnergy is the Turner 2004 nearest-neighbor stacking sum (MIRNA-PAIR-001);
            // it is 0.0 when a sparse duplex has no consecutive paired stacks, so assert it is
            // a finite value rather than non-zero (the old != 0.0 check assumed the removed
            // invented per-position mismatch penalty).
            Assert.That(double.IsFinite(site.FreeEnergy), Is.True);
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

    // ─────────────────────────────────────────────────────────────────────────────────────
    // TargetScan context++ scoring (Agarwal et al. 2015) — opt-in ScoreTargetSiteContextPlusPlus
    //
    // Source (coefficients, retrieved verbatim this session):
    //   Agarwal_2015_parameters.txt — TargetScan distribution
    //   https://raw.githubusercontent.com/nsoranzo/targetscan/main/Agarwal_2015_parameters.txt
    // Source (feature computation/scaling, retrieved verbatim this session):
    //   targetscan_70_context_scores.pl getAgarwalContribution / getLocalAU_contribution /
    //   get_sRNA1_8_contributions / getSite8_contribution / get3primePairingContribution /
    //   getMinDist_weighted_contribution / get_len3UTR_weighted_contribution /
    //   getOffset6mer_weighted_contribution + getOffset6merSites
    //   https://raw.githubusercontent.com/nsoranzo/targetscan/main/targetscan_70_context_scores.pl
    // Peer-reviewed model: Agarwal V et al. (2015) eLife 4:e05005, doi:10.7554/eLife.05005.
    // 3P_score raw values below were cross-checked by running the perl on the same UTR/miRNA.
    //
    // Each expected value below is hand-derived from those retrieved coefficients × feature
    // values (see per-test derivation comments), NOT copied from the implementation output.
    // ─────────────────────────────────────────────────────────────────────────────────────

    #region CTX-001: 8mer context++ — intercept + scaled local-AU + sRNA8 indicator

    [Test]
    public void ScoreTargetSiteContextPlusPlus_8merLet7a_GcFlanks_MatchesHandDerivedScore()
    {
        // let-7a-5p UGAGGUAGUAGGUUGUAUAGUU: nt1=U (⇒ sRNA1=0), nt8=G (⇒ sRNA8G).
        // 8mer layout: GGGGG + CUACCUCA + GGGGG → site Start=5,End=12 (len 18); flanks all G ⇒ local-AU fraction=0.
        // Now-computed local UTR features (all derived independently below; 3P raw=0 verified vs perl):
        //   3P_score : raw=0 (no 3' complementarity) ⇒ scaled=(0-1)/(3.5-1)=-0.4 ; coeff(8mer)=-0.040 ⇒ +0.016
        //   Min_dist : perlStart=6,distTo5=5 ; perlEnd=13,distTo3=18-13=5 ; min=5 ; log10(5)=0.698970004336
        //              scaled=(0.698970004336-1.415)/(3.113-1.415) ; coeff=0.118 ⇒ -0.049759446106213065
        //   Len_3UTR : log10(18)=1.255272505103 ; scaled=(…-2.392)/(3.637-2.392) ; coeff=0.310 ⇒ -0.2830405810586145
        //   Off6m    : pattern "CUACCU" (first 6 of revcomp of seed GAGGUAG) occurs 1× ; coeff(8mer)=-0.020 ⇒ -0.020
        //   Intercept(8mer)=-0.589 ; Local_AU=-0.254×((0-0.308)/(0.814-0.308))=+0.154608695652174 ; sRNA8G(8mer)=+0.015
        //   ⇒ CS_partial = -0.7561913315126536
        string mrna = "GGGGG" + Let7aSeedRC + "A" + "GGGGG";
        var site = FindTargetSites(mrna, Let7a, minScore: 0.0).Single();

        var ctx = ScoreTargetSiteContextPlusPlus(mrna, Let7a, site);

        Assert.Multiple(() =>
        {
            Assert.That(site.Type, Is.EqualTo(TargetSiteType.Seed8mer), "fixture must be an 8mer");
            Assert.That(ctx.Intercept, Is.EqualTo(-0.589).Within(1e-9), "8mer intercept (Agarwal params)");
            Assert.That(ctx.LocalAuContribution, Is.EqualTo(0.154608695652174).Within(1e-9),
                "Local_AU = -0.254×((0-0.308)/(0.814-0.308)) with all-G flanks");
            Assert.That(ctx.SRna1Contribution, Is.EqualTo(0.0).Within(1e-9), "nt1=U ⇒ sRNA1 not scored");
            Assert.That(ctx.SRna8Contribution, Is.EqualTo(0.015).Within(1e-9), "nt8=G ⇒ sRNA8G(8mer)=0.015");
            Assert.That(ctx.Site8Contribution, Is.EqualTo(0.0).Within(1e-9), "Site8 undefined for 8mer");
            Assert.That(ctx.ThreePrimePairingContribution, Is.EqualTo(0.016).Within(1e-9),
                "3P raw=0 ⇒ -0.040×((0-1)/(3.5-1)) = +0.016");
            Assert.That(ctx.MinDistContribution, Is.EqualTo(-0.049759446106213065).Within(1e-9),
                "Min_dist = 0.118×((log10(5)-1.415)/(3.113-1.415))");
            Assert.That(ctx.Len3UtrContribution, Is.EqualTo(-0.2830405810586145).Within(1e-9),
                "Len_3UTR = 0.310×((log10(18)-2.392)/(3.637-2.392))");
            Assert.That(ctx.Off6mContribution, Is.EqualTo(-0.020).Within(1e-9),
                "Off6m = -0.020×1 (one offset-6mer 'CUACCU' in the UTR)");
            Assert.That(ctx.ContextScorePartial, Is.EqualTo(-0.7561913315126536).Within(1e-9),
                "8mer partial context++ = sum of realised contributions");
        });
    }

    #endregion

    #region CTX-002: 7mer-m8 context++ — scaled local-AU + sRNA8G(7mer-m8)

    [Test]
    public void ScoreTargetSiteContextPlusPlus_7merM8Let7a_MatchesHandDerivedScore()
    {
        // 7mer-m8 layout: GGGGG + CUACCUC + GGGG → site Start=5,End=11; all-G flanks ⇒ fraction=0.
        // CS_partial = Intercept(7mer-m8) + Local_AU + sRNA8G(7mer-m8)
        //   Intercept(7mer-m8)       = -0.224
        //   Local_AU = -0.177 × ((0 - 0.277)/(0.782 - 0.277)) = +0.097087128712871
        //   sRNA8G(7mer-m8)          = -0.008  ;  sRNA1=0, Site8=0 (7mer-m8)
        //   ⇒ CS_partial            = -0.134912871287129
        string mrna = "GGGGG" + Let7aSeedRC + "GGGG";
        var site = FindTargetSites(mrna, Let7a, minScore: 0.0).Single();

        var ctx = ScoreTargetSiteContextPlusPlus(mrna, Let7a, site);

        Assert.Multiple(() =>
        {
            Assert.That(site.Type, Is.EqualTo(TargetSiteType.Seed7merM8), "fixture must be a 7mer-m8");
            Assert.That(ctx.Intercept, Is.EqualTo(-0.224).Within(1e-9), "7mer-m8 intercept");
            Assert.That(ctx.LocalAuContribution, Is.EqualTo(0.097087128712871).Within(1e-9),
                "Local_AU = -0.177×((0-0.277)/(0.782-0.277))");
            Assert.That(ctx.SRna8Contribution, Is.EqualTo(-0.008).Within(1e-9), "nt8=G ⇒ sRNA8G(7mer-m8)=-0.008");
            Assert.That(ctx.Site8Contribution, Is.EqualTo(0.0).Within(1e-9), "Site8 undefined for 7mer-m8");
            // 3P raw=0 ⇒ -0.055×((0-1)/2.5)=+0.022 ; Off6m 'CUACCU' once ⇒ -0.011 ; len=16, site Start=5/End=11.
            Assert.That(ctx.ThreePrimePairingContribution, Is.EqualTo(0.022).Within(1e-9),
                "3P raw=0 ⇒ -0.055×((0-1)/(3.5-1)) = +0.022");
            Assert.That(ctx.MinDistContribution, Is.EqualTo(-0.031015975380457392).Within(1e-9),
                "Min_dist(7mer-m8) = 0.056×((log10(4)-1.491)/(3.096-1.491)) [distTo3=16-12=4]");
            Assert.That(ctx.Len3UtrContribution, Is.EqualTo(-0.15385698397262643).Within(1e-9),
                "Len_3UTR(7mer-m8) = 0.154×((log10(16)-2.409)/(3.615-2.409))");
            Assert.That(ctx.Off6mContribution, Is.EqualTo(-0.011).Within(1e-9), "Off6m = -0.011×1");
            Assert.That(ctx.ContextScorePartial, Is.EqualTo(-0.3087858306402126).Within(1e-9),
                "7mer-m8 partial context++");
        });
    }

    #endregion

    #region CTX-003: 7mer-A1 context++ — Site8 indicator path

    [Test]
    public void ScoreTargetSiteContextPlusPlus_7merA1Let7a_Site8G_MatchesHandDerivedScore()
    {
        // 7mer-A1 layout: GGGGG + UACCUC + A + GGGG → site Start=5,End=11; Site8 base = mrna[4]='G'.
        // CS_partial = Intercept(7mer-A1) + Local_AU + sRNA8G(7mer-A1) + Site8G(7mer-A1)
        //   Intercept(7mer-A1)       = -0.195
        //   Local_AU = -0.075 × ((0 - 0.342)/(0.801 - 0.342)) = +0.055882352941176
        //   sRNA8G(7mer-A1)          = -0.017
        //   Site8G(7mer-A1)          = +0.015  ;  sRNA1=0 (nt1=U)
        //   ⇒ CS_partial            = -0.141117647058824
        string mrna = "GGGGG" + "UACCUC" + "A" + "GGGG";
        var site = FindTargetSites(mrna, Let7a, minScore: 0.0).Single();

        var ctx = ScoreTargetSiteContextPlusPlus(mrna, Let7a, site);

        Assert.Multiple(() =>
        {
            Assert.That(site.Type, Is.EqualTo(TargetSiteType.Seed7merA1), "fixture must be a 7mer-A1");
            Assert.That(ctx.Intercept, Is.EqualTo(-0.195).Within(1e-9), "7mer-A1 intercept");
            Assert.That(ctx.LocalAuContribution, Is.EqualTo(0.055882352941176).Within(1e-9),
                "Local_AU = -0.075×((0-0.342)/(0.801-0.342))");
            Assert.That(ctx.SRna8Contribution, Is.EqualTo(-0.017).Within(1e-9), "nt8=G ⇒ sRNA8G(7mer-A1)=-0.017");
            Assert.That(ctx.Site8Contribution, Is.EqualTo(0.015).Within(1e-9),
                "Site8 base 'G' ⇒ Site8G(7mer-A1)=0.015 (only defined for 7mer-A1/6mer)");
            // 3P raw=0 ⇒ -0.060×((0-1)/2.5)=+0.024 ; Off6m: 'CUACCU' absent (UACCUC core only) ⇒ 0 ; len=16, Start=5/End=11.
            Assert.That(ctx.ThreePrimePairingContribution, Is.EqualTo(0.024).Within(1e-9),
                "3P raw=0 ⇒ -0.060×((0-1)/(3.5-1)) = +0.024");
            Assert.That(ctx.MinDistContribution, Is.EqualTo(-0.022124733327545488).Within(1e-9),
                "Min_dist(7mer-A1) = 0.045×((log10(4)-1.431)/(3.117-1.431))");
            Assert.That(ctx.Len3UtrContribution, Is.EqualTo(-0.12813929518273268).Within(1e-9),
                "Len_3UTR(7mer-A1) = 0.129×((log10(16)-2.413)/(3.630-2.413))");
            Assert.That(ctx.Off6mContribution, Is.EqualTo(0.0).Within(1e-9),
                "no 'CUACCU' offset-6mer in this UTR ⇒ Off6m=0");
            Assert.That(ctx.ContextScorePartial, Is.EqualTo(-0.2673816755691017).Within(1e-9),
                "7mer-A1 partial context++");
        });
    }

    #endregion

    #region CTX-004: 6mer context++ — non-zero local-AU fraction + Site8C

    [Test]
    public void ScoreTargetSiteContextPlusPlus_6merMiR21_MixedFlanks_MatchesHandDerivedScore()
    {
        // miR-21 seedRC AUAAGCU, 6mer core = UAAGCU. nt1=U (sRNA1=0), nt8=U (sRNA8=0).
        // Layout: GGGGC + UAAGCU + UGAGG → site Start=5,End=10; Site8 base = mrna[4]='C'.
        // Local-AU (6mer weights: up 1/(i+2), down 1/(i+1)): only downstream U(i0),A(i2) are A/U
        //   fraction = (1 + 1/3) / (Σ up + Σ down) = 0.357142857142857
        // CS_partial = Intercept(6mer) + Local_AU + Site8C(6mer)
        //   Intercept(6mer)          = -0.079
        //   Local_AU = -0.040 × ((0.357142857142857 - 0.295)/(0.772 - 0.295)) = -0.005211141060198
        //   Site8C(6mer)             = +0.015  ;  sRNA1=0, sRNA8=0
        //   ⇒ CS_partial            = -0.069211141060198
        string mrna = "GGGGC" + "UAAGCU" + "UGAGG";
        var site = FindTargetSites(mrna, MiR21, minScore: 0.0).Single();

        var ctx = ScoreTargetSiteContextPlusPlus(mrna, MiR21, site);

        Assert.Multiple(() =>
        {
            Assert.That(site.Type, Is.EqualTo(TargetSiteType.Seed6mer), "fixture must be a 6mer");
            Assert.That(ctx.Intercept, Is.EqualTo(-0.079).Within(1e-9), "6mer intercept");
            Assert.That(ctx.LocalAuContribution, Is.EqualTo(-0.005211141060198).Within(1e-9),
                "Local_AU with downstream A/U fraction 0.357142857142857, scaled ×(-0.040)");
            Assert.That(ctx.SRna1Contribution, Is.EqualTo(0.0).Within(1e-9), "nt1=U");
            Assert.That(ctx.SRna8Contribution, Is.EqualTo(0.0).Within(1e-9), "nt8=U ⇒ sRNA8 not scored");
            Assert.That(ctx.Site8Contribution, Is.EqualTo(0.015).Within(1e-9),
                "Site8 base 'C' ⇒ Site8C(6mer)=0.015");
            // 3P raw=0 ⇒ -0.024×((0-1)/2.5)=+0.0096 ; miR-21 seed AGCUUAU ⇒ Off6m pattern absent ⇒ 0 ; len=16, Start=5/End=10.
            Assert.That(ctx.ThreePrimePairingContribution, Is.EqualTo(0.0096).Within(1e-9),
                "3P raw=0 ⇒ -0.024×((0-1)/(3.5-1)) = +0.0096");
            Assert.That(ctx.MinDistContribution, Is.EqualTo(-0.017194033053347654).Within(1e-9),
                "Min_dist(6mer) = 0.036×((log10(5)-1.477)/(3.106-1.477)) [distTo3=16-11=5]");
            Assert.That(ctx.Len3UtrContribution, Is.EqualTo(-0.04447703767941017).Within(1e-9),
                "Len_3UTR(6mer) = 0.045×((log10(16)-2.405)/(3.620-2.405))");
            Assert.That(ctx.Off6mContribution, Is.EqualTo(0.0).Within(1e-9), "no miR-21 offset-6mer in this UTR");
            Assert.That(ctx.ContextScorePartial, Is.EqualTo(-0.12128221179295548).Within(1e-9),
                "6mer partial context++");
        });
    }

    #endregion

    #region CTX-005: sRNA1 non-U indicator branch — coefficient lookup

    [Test]
    public void ScoreTargetSiteContextPlusPlus_SRna1G_8mer_AddsSRna1GCoefficient()
    {
        // The sRNA1 indicator is only scored when miRNA nt1 ≠ U (perl: sRNA1_nt ne "U").
        // Construct a miRNA with nt1=G, nt8=G so both indicator branches fire on an 8mer site.
        // sRNA1G(8mer)=+0.060, sRNA8G(8mer)=+0.015.  miRNA: G GAGGUAG ... (seed pos2-8 = GAGGUAG = let-7a seed)
        var miR = CreateMiRna("synthetic-G1", "GGAGGUAGUAGGUUGUAUAGUU");
        string mrna = "GGGGG" + GetReverseComplement(miR.SeedSequence) + "A" + "GGGGG";
        var site = FindTargetSites(mrna, miR, minScore: 0.0).Single();

        var ctx = ScoreTargetSiteContextPlusPlus(mrna, miR, site);

        Assert.Multiple(() =>
        {
            Assert.That(site.Type, Is.EqualTo(TargetSiteType.Seed8mer), "fixture must be an 8mer");
            Assert.That(ctx.SRna1Contribution, Is.EqualTo(0.060).Within(1e-9),
                "nt1=G ⇒ sRNA1G(8mer)=0.060 (Agarwal params)");
            Assert.That(ctx.SRna8Contribution, Is.EqualTo(0.015).Within(1e-9),
                "nt8=G ⇒ sRNA8G(8mer)=0.015");
        });
    }

    #endregion

    #region CTX-006: invalid site type rejected — contract

    [Test]
    public void ScoreTargetSiteContextPlusPlus_NonSeedSiteType_Throws()
    {
        // context++ models are defined only for the four canonical seed-match types (Agarwal 2015).
        var offsetSite = new TargetSite(
            Start: 0, End: 5, TargetSequence: "AAAAAA", MiRnaName: "x",
            Type: TargetSiteType.Offset6mer, SeedMatchLength: 6, Score: 0, FreeEnergy: 0, Alignment: "");

        Assert.Throws<ArgumentException>(
            () => ScoreTargetSiteContextPlusPlus("AAAAAAAAAA", Let7a, offsetSite),
            "Offset6mer / Supplementary / Centered have no fitted context++ model");
    }

    #endregion

    #region CTX-007: residual features reported when no optional inputs supplied — honest residual

    [Test]
    public void ScoreTargetSiteContextPlusPlus_NoOptionalInputs_ReportsResidualFeatures()
    {
        // With no caller-supplied inputs on this SHORT layout, the residual is: SA (the 14-nt
        // accessibility window does not fit — windowStart0 = Start+7-13 = -1 < 0), PCT (always
        // blocked), and SPS, TA_3UTR, Len_ORF, ORF8m (data-blocked unless supplied). 3P_score,
        // Min_dist, Len_3UTR and Off6m are now COMPUTED, so they must NOT appear as omitted.
        string mrna = "GGGGG" + Let7aSeedRC + "A" + "GGGGG";
        var site = FindTargetSites(mrna, Let7a, minScore: 0.0).Single();

        var ctx = ScoreTargetSiteContextPlusPlus(mrna, Let7a, site);

        Assert.Multiple(() =>
        {
            Assert.That(ctx.OmittedFeatures, Has.Some.Contains("SA"),
                "SA residual here because the 14-nt window does not fit this short UTR");
            Assert.That(ctx.SaContribution, Is.EqualTo(0.0), "out-of-fit SA contributes 0");
            Assert.That(ctx.OmittedFeatures, Has.Some.Contains("PCT"), "PCT stays residual (needs alignment)");
            Assert.That(ctx.OmittedFeatures, Has.Some.Contains("SPS"), "SPS omitted unless supplied");
            Assert.That(ctx.OmittedFeatures, Has.Some.Contains("TA_3UTR"), "TA omitted unless supplied");
            Assert.That(ctx.OmittedFeatures, Has.Some.Contains("Len_ORF"), "Len_ORF omitted unless supplied");
            Assert.That(ctx.OmittedFeatures, Has.Some.Contains("ORF8m"), "ORF8m omitted unless supplied");
            Assert.That(ctx.OmittedFeatures, Has.None.Contains("3P_score"),
                "3' supplementary pairing is now computed, not omitted");
            Assert.That(ctx.OmittedFeatures, Has.None.Contains("Min_dist"), "Min_dist is now computed");
            Assert.That(ctx.OmittedFeatures, Has.None.Contains("Len_3UTR"), "Len_3UTR is now computed");
            Assert.That(ctx.OmittedFeatures, Has.None.Contains("Off6m"), "Off6m is now computed");
        });
    }

    #endregion

    #region CTX-008: 3P_score with real 3' supplementary pairing — DP raw score (vs targetscan perl)

    [Test]
    public void ScoreTargetSiteContextPlusPlus_3PrimeSupplementaryPairing_8mer_MatchesScaledRawScore()
    {
        // UTR carries a 3'-supplementary complement (ACAACCUA…) upstream of the let-7a seed match.
        // get3primePairingContribution (faithful port) yields raw 3P score = 6 for this site, which
        // was reproduced exactly by running targetscan_70_context_scores.pl on the same UTR/miRNA.
        // scaled = (6 - 1)/(3.5 - 1) = 2.0 ; coeff(8mer)=-0.040 ⇒ 3P contribution = -0.080.
        string mrna = "AAAAAACAACCUAACUACCUCAGGG";
        var site = FindTargetSites(mrna, Let7a, minScore: 0.0).First(s => s.Type == TargetSiteType.Seed8mer);

        var ctx = ScoreTargetSiteContextPlusPlus(mrna, Let7a, site);

        Assert.That(ctx.ThreePrimePairingContribution, Is.EqualTo(-0.080).Within(1e-9),
            "raw 3P=6 (perl-verified) ⇒ -0.040×((6-1)/(3.5-1)) = -0.080");
    }

    #endregion

    #region CTX-009: 3P_score across site types — raw scores reproduce the perl reference

    [TestCase("GGGACAACCUAGGGGGCUACCUCAGGG", TargetSiteType.Seed8mer, -0.040, 4.5)]   // raw 4.5
    [TestCase("AAAACAACCUAAUACCUCAGGGG", TargetSiteType.Seed7merA1, -0.060, 6.0)]     // raw 6
    public void ScoreTargetSiteContextPlusPlus_3PrimeRawScore_PerlReference(
        string mrna, TargetSiteType expectedType, double coeff3P, double rawScore)
    {
        // Raw 3P score for each fixture was computed by running the TargetScan reference perl
        // (get3primePairingContribution) on the identical UTR/miRNA/site coordinates.
        // Expected contribution = coeff(siteType,3P_score) × ((raw - 1)/(3.5 - 1)).
        double expected = coeff3P * ((rawScore - 1.0) / (3.5 - 1.0));
        var site = FindTargetSites(mrna, Let7a, minScore: 0.0).First(s => s.Type == expectedType);

        var ctx = ScoreTargetSiteContextPlusPlus(mrna, Let7a, site);

        Assert.Multiple(() =>
        {
            Assert.That(site.Type, Is.EqualTo(expectedType), "fixture must be the expected site type");
            Assert.That(ctx.ThreePrimePairingContribution, Is.EqualTo(expected).Within(1e-9),
                $"3P = {coeff3P}×((raw {rawScore} - 1)/2.5)");
        });
    }

    #endregion

    #region CTX-010: Off6m counts every offset-6mer occurrence (used raw, no scaling)

    [Test]
    public void ScoreTargetSiteContextPlusPlus_TwoOffset6mers_CountedRaw()
    {
        // Offset-6mer pattern = first 6 nt of revcomp(let-7a seed GAGGUAG) = "CUACCU".
        // This UTR contains it twice (a bare CUACCU and the CUACCU inside the 8mer CUACCUCA).
        // Off6m is used RAW (not min-max scaled): contribution = coeff(8mer) × count = -0.020 × 2.
        string mrna = "GGGCUACCUGGGGGCUACCUCAGGG";
        var site = FindTargetSites(mrna, Let7a, minScore: 0.0).First(s => s.Type == TargetSiteType.Seed8mer);

        var ctx = ScoreTargetSiteContextPlusPlus(mrna, Let7a, site);

        Assert.That(ctx.Off6mContribution, Is.EqualTo(-0.040).Within(1e-9),
            "two 'CUACCU' offset-6mers ⇒ -0.020 × 2 = -0.040");
    }

    #endregion

    #region CTX-011: caller-supplied SPS / TA / Len_ORF / ORF8m computed and removed from residual

    [Test]
    public void ScoreTargetSiteContextPlusPlus_SuppliedInputs_8mer_MatchHandDerivedAndDropFromResidual()
    {
        // Faithful Agarwal contributions for caller-supplied values on an 8mer site:
        //   SPS=-8.0   : scaled=(-8.0-(-11.13))/(-5.52-(-11.13)) ; coeff=0.210 ⇒ +0.11716577540106952
        //   TA=3.5     : scaled=(3.5-3.113)/(3.865-3.113)         ; coeff=0.222 ⇒ +0.11424734042553189
        //   Len_ORF=1000: log10=3 ; scaled=(3-2.788)/(3.753-2.788); coeff=0.205 ⇒ +0.04503626943005184
        //   ORF8m=2    : used raw ; coeff=-0.118 × 2              ⇒ -0.236
        string mrna = "AAAAAACAACCUAACUACCUCAGGG";
        var site = FindTargetSites(mrna, Let7a, minScore: 0.0).First(s => s.Type == TargetSiteType.Seed8mer);
        var inputs = new ContextPlusPlusInputs(Sps: -8.0, Ta: 3.5, OrfLength: 1000, Orf8mCount: 2);

        var ctx = ScoreTargetSiteContextPlusPlus(mrna, Let7a, site, inputs);

        Assert.Multiple(() =>
        {
            Assert.That(ctx.SpsContribution, Is.EqualTo(0.11716577540106952).Within(1e-9),
                "SPS = 0.210×((-8.0+11.13)/(-5.52+11.13))");
            Assert.That(ctx.TaContribution, Is.EqualTo(0.11424734042553189).Within(1e-9),
                "TA = 0.222×((3.5-3.113)/(3.865-3.113))");
            Assert.That(ctx.LenOrfContribution, Is.EqualTo(0.04503626943005184).Within(1e-9),
                "Len_ORF = 0.205×((log10(1000)-2.788)/(3.753-2.788))");
            Assert.That(ctx.Orf8mContribution, Is.EqualTo(-0.236).Within(1e-9),
                "ORF8m raw = -0.118×2");
            Assert.That(ctx.OmittedFeatures, Has.None.Contains("SPS"), "supplied SPS leaves residual");
            Assert.That(ctx.OmittedFeatures, Has.None.Contains("TA_3UTR"), "supplied TA leaves residual");
            Assert.That(ctx.OmittedFeatures, Has.None.Contains("Len_ORF"), "supplied Len_ORF leaves residual");
            Assert.That(ctx.OmittedFeatures, Has.None.Contains("ORF8m"), "supplied ORF8m leaves residual");
            // SA is now computed from the Turner-2004 McCaskill partition function (the 14-nt
            // window fits this UTR: windowStart0 = Start+7-13 = 8 ≥ 0), so it is NO LONGER residual.
            // Only PCT (multi-species conservation) remains residual once all data inputs are supplied.
            Assert.That(ctx.SaContribution, Is.Not.EqualTo(0.0),
                "SA is computed (window fits) — no longer an honest residual");
            Assert.That(ctx.OmittedFeatures, Has.None.Contains("SA"), "SA computed ⇒ not in residual");
            Assert.That(ctx.OmittedFeatures, Has.Some.Contains("PCT"), "PCT still residual");
        });
    }

    #endregion

    #region CTX-SA: structural accessibility wired from the Turner-2004 McCaskill partition function

    // CTX-SA-001 — SA contribution equals the verbatim getSA_contribution / getAgarwalContribution
    // arithmetic: coeff(SA,8mer) × (log10(plfold) - min)/(max - min), where plfold is the 14-nt
    // window unpaired probability from the Turner-2004 McCaskill partition function. The expected
    // value is recomputed INDEPENDENTLY from CalculateRegionUnpairedProbability + the verbatim
    // Agarwal_2015_parameters.txt SA row (8mer: coeff -0.115, min -4.356, max -0.661), so it would
    // fail a wrong coefficient, a wrong window, or an MFE-as-accessibility implementation.
    [Test]
    public void ScoreTargetSiteContextPlusPlus_SA_8mer_MatchesHandDerivedAccessibility()
    {
        // 48-nt UTR: 20-nt structured 5' flank + 8mer let-7a site (CUACCUCA) + 20-nt 3' flank.
        string mrna = "GGGGCCCCGGGGCCCCGGGG" + "CUACCUCA" + "GGGGCCCCGGGGCCCCGGGG";
        var site = FindTargetSites(mrna, Let7a, minScore: 0.0).First(s => s.Type == TargetSiteType.Seed8mer);

        // Independently reproduce the SA local-context computation (getSA_contribution semantics):
        // 8mer ⇒ no utrStart decrement; row read = utrStart+7 ⇒ windowEnd0 = Start+7; window L=14.
        string seq = mrna.ToUpperInvariant().Replace('T', 'U');
        const int W = 80, U = 14, RowOff = 7;
        int windowEnd0 = (site.Start + 1) + RowOff - 1;     // 1-based utrStart+7, back to 0-based
        int windowStart0 = windowEnd0 - U + 1;
        int contextStart = Math.Max(0, windowEnd0 - (W - U) / 2 - U + 1);
        contextStart = Math.Min(contextStart, windowStart0);
        int contextEnd = Math.Min(seq.Length - 1, contextStart + W - 1);
        contextStart = Math.Max(0, contextEnd - W + 1);
        string context = seq.Substring(contextStart, contextEnd - contextStart + 1);
        int localWindowEnd = windowEnd0 - contextStart;
        double plfold = Seqeron.Genomics.Analysis.RnaSecondaryStructure
            .CalculateRegionUnpairedProbability(context, localWindowEnd, U);
        double log10 = Math.Log10(plfold);
        // Verbatim SA row for 8mer (Agarwal_2015_parameters.txt): coeff -0.115, min -4.356, max -0.661.
        double expectedSa = -0.115 * ((log10 - (-4.356)) / ((-0.661) - (-4.356)));

        var ctx = ScoreTargetSiteContextPlusPlus(mrna, Let7a, site);

        Assert.Multiple(() =>
        {
            Assert.That(plfold, Is.GreaterThan(0.0).And.LessThanOrEqualTo(1.0),
                "the 14-nt window accessibility is a probability in (0,1]");
            Assert.That(ctx.SaContribution, Is.EqualTo(expectedSa).Within(1e-12),
                "SA = coeff × (log10(plfold) - min)/(max - min) with the verbatim 8mer SA parameters");
            Assert.That(ctx.SaContribution, Is.Not.EqualTo(0.0),
                "SA is computed (the window fits), not an honest residual");
            Assert.That(ctx.OmittedFeatures, Has.None.Contains("SA"),
                "computed SA is not reported as omitted");
            // SA must be part of the partial sum.
            Assert.That(ctx.ContextScorePartial, Is.EqualTo(
                ctx.Intercept + ctx.LocalAuContribution + ctx.SRna1Contribution + ctx.SRna8Contribution
                + ctx.Site8Contribution + ctx.SaContribution + ctx.ThreePrimePairingContribution
                + ctx.MinDistContribution + ctx.Len3UtrContribution + ctx.Off6mContribution
                + ctx.SpsContribution + ctx.TaContribution + ctx.LenOrfContribution
                + ctx.Orf8mContribution + ctx.PctContribution).Within(1e-12),
                "ContextScorePartial includes the SA contribution");
        });
    }

    #endregion

    #region CTX-PCT: Friedman 2009 branch-length score (Bls) + PCT wired into context++

    // Worked phylogenetic tree used across the Bls cases (Newick, explicit branch lengths):
    //   ((A:1.0,B:2.0):0.5,(C:1.5,D:3.0):4.0);
    // Internal node (A,B) connects to the root by an edge of length 0.5; (C,D) by 4.0.
    private const string WorkedTreeNewick = "((A:1.0,B:2.0):0.5,(C:1.5,D:3.0):4.0);";

    // CTX-PCT-001 — Bls = total branch length of the minimal subtree connecting the species in
    // which the site is conserved (Friedman et al. 2009 Genome Res 19:92, Methods). Hand-derived
    // for several conserved-species subsets on the worked tree:
    //   {A,B}      : A(1.0)+B(2.0) — the (A,B)→root edge is NOT counted (both A and B are below it,
    //                so no conserved species lies outside) ⇒ 3.0
    //   {A,C}      : A(1.0)+ (A,B)→root(0.5) + (C,D)→root(4.0) + C(1.5)               ⇒ 7.0
    //   {A,B,C,D}  : every leaf edge (1+2+1.5+3) + both internal edges (0.5+4.0)       ⇒ 12.0
    //   {A}        : a single species — no connecting subtree                          ⇒ 0.0
    [Test]
    [TestCase(new[] { "A", "B" }, 3.0)]
    [TestCase(new[] { "A", "C" }, 7.0)]
    [TestCase(new[] { "A", "B", "C", "D" }, 12.0)]
    [TestCase(new[] { "A" }, 0.0)]
    public void ComputeBranchLengthScore_WorkedTree_MatchesHandDerivedBls(string[] species, double expectedBls)
    {
        PhylogeneticAnalyzer.PhyloNode tree = PhylogeneticAnalyzer.ParseNewick(WorkedTreeNewick);

        double bls = ComputeBranchLengthScore(tree, species);

        Assert.That(bls, Is.EqualTo(expectedBls).Within(1e-9),
            $"Bls({string.Join(",", species)}) = total branch length of the minimal connecting subtree");
    }

    // CTX-PCT-002 — empty / null conserved-species set ⇒ Bls = 0 (single or no species ⇒ no subtree).
    [Test]
    public void ComputeBranchLengthScore_NoSpecies_IsZero()
    {
        PhylogeneticAnalyzer.PhyloNode tree = PhylogeneticAnalyzer.ParseNewick(WorkedTreeNewick);

        Assert.Multiple(() =>
        {
            Assert.That(ComputeBranchLengthScore(tree, Array.Empty<string>()), Is.EqualTo(0.0),
                "no conserved species ⇒ Bls 0");
            Assert.Throws<ArgumentNullException>(() => ComputeBranchLengthScore(tree, null!),
                "null species set rejected");
            Assert.Throws<ArgumentNullException>(() => ComputeBranchLengthScore(null!, new[] { "A" }),
                "null tree rejected");
        });
    }

    // CTX-PCT-003 — PCT(Bls) = B0 + B1/(1 + e^(−B2·Bls + B3)), truncated at 0
    // (targetscan_70_BL_PCT.pl, calculatePCTthisBL). Hand-derived with the simple parameters
    // (B0=0, B1=1, B2=1, B3=0):  PCT(3.0) = 1/(1+e^-3) = 0.952574126822433.
    [Test]
    public void PctFromBranchLength_SimpleSigmoid_MatchesHandDerivedValue()
    {
        var p = new PctSigmoidParameters(B0: 0.0, B1: 1.0, B2: 1.0, B3: 0.0);

        double pct = PctFromBranchLength(3.0, p);

        Assert.That(pct, Is.EqualTo(0.952574126822433).Within(1e-12),
            "PCT(3.0) = 0 + 1/(1+e^(−1·3+0)) per the published logistic relationship");
    }

    // CTX-PCT-004 — negative PCT values are truncated to 0 (perl: if ($pct < 0) { $pct = "0.0"; }).
    [Test]
    public void PctFromBranchLength_NegativeRaw_TruncatedToZero()
    {
        // B0=-0.5, B1=0.3, B3=5 with Bls=0 ⇒ raw = -0.5 + 0.3/(1+e^5) ≈ -0.498 < 0 ⇒ 0.
        var p = new PctSigmoidParameters(B0: -0.5, B1: 0.3, B2: 1.0, B3: 5.0);

        double pct = PctFromBranchLength(0.0, p);

        Assert.That(pct, Is.EqualTo(0.0),
            "a negative logistic output is truncated to 0 as in the TargetScan reference");
    }

    // CTX-PCT-005 — the PCT contribution = coeff(PCT) × (PCT − min)/(max − min) enters context++
    // and PCT leaves OmittedFeatures when a Conservation input is supplied. For an 8mer with the
    // worked tree and {A,B} (Bls=3.0), sigmoid (0,1,1,0): PCT=0.952574126822433.
    //   8mer PCT row (Agarwal_2015_parameters.txt): coeff -0.103, min 0, max 0.816.
    //   contribution = -0.103 × (0.952574126822433 / 0.816) = -0.120239136106263.
    [Test]
    public void ScoreTargetSiteContextPlusPlus_ConservationSupplied_8mer_PctEntersScoreAndDropsResidual()
    {
        string mrna = "AAAAAACAACCUAACUACCUCAGGG";
        var site = FindTargetSites(mrna, Let7a, minScore: 0.0).First(s => s.Type == TargetSiteType.Seed8mer);

        PhylogeneticAnalyzer.PhyloNode tree = PhylogeneticAnalyzer.ParseNewick(WorkedTreeNewick);
        var conservation = new PctConservation(
            tree,
            new[] { "A", "B" },
            new PctSigmoidParameters(B0: 0.0, B1: 1.0, B2: 1.0, B3: 0.0));
        var inputs = new ContextPlusPlusInputs(Conservation: conservation);

        var ctx = ScoreTargetSiteContextPlusPlus(mrna, Let7a, site, inputs);

        Assert.Multiple(() =>
        {
            Assert.That(ctx.BranchLengthScore, Is.EqualTo(3.0).Within(1e-9),
                "Bls({A,B}) on the worked tree = 3.0");
            Assert.That(ctx.Pct, Is.EqualTo(0.952574126822433).Within(1e-12),
                "PCT = 1/(1+e^-3) from the supplied sigmoid");
            Assert.That(ctx.PctContribution, Is.EqualTo(-0.120239136106263).Within(1e-9),
                "PCT contribution = -0.103 × (PCT/0.816) with the verbatim 8mer PCT parameters");
            Assert.That(ctx.OmittedFeatures, Has.None.Contains("PCT"),
                "supplied conservation ⇒ PCT no longer an honest residual");
            Assert.That(ctx.ContextScorePartial, Is.EqualTo(
                ctx.Intercept + ctx.LocalAuContribution + ctx.SRna1Contribution + ctx.SRna8Contribution
                + ctx.Site8Contribution + ctx.SaContribution + ctx.ThreePrimePairingContribution
                + ctx.MinDistContribution + ctx.Len3UtrContribution + ctx.Off6mContribution
                + ctx.SpsContribution + ctx.TaContribution + ctx.LenOrfContribution
                + ctx.Orf8mContribution + ctx.PctContribution).Within(1e-12),
                "ContextScorePartial includes the PCT contribution");
        });
    }

    // CTX-PCT-006 — per-site-type PCT parameters: a 7mer-m8 site with {A,C} (Bls=7.0) and the same
    // sigmoid uses the 7mer-m8 PCT row (coeff -0.048, min 0, max 0.364):
    //   PCT(7.0) = 1/(1+e^-7) = 0.999088948805599
    //   contribution = -0.048 × (0.999088948805599 / 0.364) = -0.131747993249090.
    [Test]
    public void ScoreTargetSiteContextPlusPlus_Conservation_7merM8_UsesSiteTypePctParameters()
    {
        // 7mer-m8 = seed match to miRNA pos 2-8 WITHOUT a trailing A (no 'A' opposite pos 1).
        // let-7a seedRC = CUACCUC ; place it with a non-A base 3' of the site to force 7mer-m8.
        string mrna = "GGGGG" + Let7aSeedRC + "G" + "GGGGG";
        var site = FindTargetSites(mrna, Let7a, minScore: 0.0).First(s => s.Type == TargetSiteType.Seed7merM8);

        PhylogeneticAnalyzer.PhyloNode tree = PhylogeneticAnalyzer.ParseNewick(WorkedTreeNewick);
        var conservation = new PctConservation(
            tree,
            new[] { "A", "C" },
            new PctSigmoidParameters(B0: 0.0, B1: 1.0, B2: 1.0, B3: 0.0));

        var ctx = ScoreTargetSiteContextPlusPlus(mrna, Let7a, site,
            new ContextPlusPlusInputs(Conservation: conservation));

        Assert.Multiple(() =>
        {
            Assert.That(ctx.BranchLengthScore, Is.EqualTo(7.0).Within(1e-9), "Bls({A,C}) = 7.0");
            Assert.That(ctx.Pct, Is.EqualTo(0.999088948805599).Within(1e-12), "PCT = 1/(1+e^-7)");
            Assert.That(ctx.PctContribution, Is.EqualTo(-0.131747993249090).Within(1e-9),
                "7mer-m8 PCT contribution uses coeff -0.048, max 0.364");
        });
    }

    #endregion
}
