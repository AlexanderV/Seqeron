// MIRNA-CONTEXT-001 — parity of ScoreTargetSiteContextPlusPlus with the TargetScan 7.0 reference code
// (targetscan_70_context_scores.pl + Agarwal_2015_parameters.txt + TA_SPS_by_seed_region.txt,
// raw.githubusercontent.com/nsoranzo/targetscan/main, run 2026-11-25 with perl 5.38.2).
//
// Reference run (B13 F15-F18). Sites from targetscan_70.pl on the 181-nt UTR below with family let-7
// (seed GAGGUAG, species 9606): 8mer-1a 31-38, 6mer 82-87, 7mer-1a 125-131. Context-score inputs:
// miRNA hsa-let-7a UGAGGUAGUAGGUUGUAUAGUU; ORF length 1500; ORF 8mer count 1; AIR file "UTR2 1 181
// 100"; PCT NA. TA/SPS come from TA_SPS_by_seed_region.txt row GAGGUAG (SPS -9.250 for 8mer/7mer-m8,
// -6.720 for 7mer-1a/6mer, TA 3.393). RNAplfold -L 40 -W 80 -u 20 was run as ViennaRNA 2.8.0
// RNA.pfl_fold_up(seq, 20, 80, 40) written in _lunp format; the l=14 values the perl reads (rows 38, 89,
// 131) are passed here as PlfoldUnpairedProbability. The perl output columns (each contribution
// printed with 3 decimals, then the context++ score) are locked below. B13 F33/F34 (WP-4/5): the same rows with
// no caller TA/SPS (bundled table) and with only an ORF sequence (targetscan_count_8mers.pl → 1500 nt, 1 8mer);
// whole-kit checks against TargetScan7_context_scores.zip (targetscan_70_context_scores_output.txt, ORF files).

using static Seqeron.Genomics.Annotation.MiRnaAnalyzer;

namespace Seqeron.Genomics.Tests.Unit.Annotation;

[TestFixture]
[Category("MIRNA-CONTEXT-001")]
public class MiRnaAnalyzer_ContextPlusPlusTargetScanParity_Tests
{
    private const string Utr =
        "GCUAAAGACAAUUACAUAACAUACACGUCACUACCUCAGCACGAAACUUGUUGGCCCAGUGUGAAUCGCUUAAGGGUUGCGUACCUCGAAGUAAG" +
        "UGUGAUGCAUACGCCUUUACUUGCUGUGGUACCUCAUCCACCCCAUCGGACUGGCAUUUUUAUUACACUCAGAAACAGAACUCGGG";

    private static readonly MiRna Let7a = CreateMiRna("hsa-let-7a", "UGAGGUAGUAGGUUGUAUAGUU");

    private static TargetSite Site(TargetSiteType type) =>
        FindTargetSites(Utr, Let7a, minScore: 0.0).Single(s => s.Type == type);

    private static ContextPlusPlusInputs PerlInputs(TargetSiteType type, double plfold) => new(
        Sps: type is TargetSiteType.Seed8mer or TargetSiteType.Seed7merM8 ? -9.25 : -6.72,
        Ta: 3.393,
        OrfLength: 1500,
        Orf8mCount: 1,
        PlfoldUnpairedProbability: plfold);

    // Columns of the perl output: intercept, 3P, local AU, Min_dist, sRNA1 (A+C+G), sRNA8 (A+C+G),
    // Site8 (A+C+G), Len_3UTR, SA, Len_ORF, ORF8m, Off6m, TA, SPS, context++ score.
    private static IEnumerable<TestCaseData> PerlRows()
    {
        yield return new TestCaseData(TargetSiteType.Seed8mer, 31, 0.009772092,
            new[] { -0.589, -0.048, -0.115, 0.004, 0.0, 0.015, 0.0, -0.033, -0.073, 0.082, -0.118, -0.020, 0.083, 0.070 },
            -0.742).SetName("TargetScan70Perl_8mer_31_38");
        yield return new TestCaseData(TargetSiteType.Seed6mer, 82, 0.008715454,
            new[] { -0.079, -0.005, -0.013, 0.010, 0.0, 0.006, 0.012, -0.005, -0.019, 0.012, -0.060, -0.010, 0.021, 0.012 },
            -0.118).SetName("TargetScan70Perl_6mer_82_87");
        yield return new TestCaseData(TargetSiteType.Seed7merA1, 125, 0.0005159016,
            new[] { -0.195, -0.036, -0.015, 0.007, 0.0, -0.017, 0.015, -0.016, -0.020, 0.027, -0.058, -0.020, 0.039, 0.032 },
            -0.257).SetName("TargetScan70Perl_7merA1_125_131");
    }

    [TestCaseSource(nameof(PerlRows))]
    public void ScoreTargetSiteContextPlusPlus_SuppliedPlfold_ReproducesTargetScanPerlOutput(
        TargetSiteType type, int utrStart1, double plfold, double[] perlColumns, double perlScore)
    {
        var site = Site(type);
        Assert.That(site.Start + 1, Is.EqualTo(utrStart1), "same 1-based UTR_start as targetscan_70.pl");

        var c = ScoreTargetSiteContextPlusPlus(Utr, Let7a, site, PerlInputs(type, plfold));

        double[] ours =
        {
            c.Intercept, c.ThreePrimePairingContribution, c.LocalAuContribution, c.MinDistContribution,
            c.SRna1Contribution, c.SRna8Contribution, c.Site8Contribution, c.Len3UtrContribution,
            c.SaContribution, c.LenOrfContribution, c.Orf8mContribution, c.Off6mContribution,
            c.TaContribution, c.SpsContribution
        };
        string[] names = { "intercept", "3P", "local AU", "Min_dist", "sRNA1", "sRNA8", "Site8", "Len_3UTR",
            "SA", "Len_ORF", "ORF8m", "Off6m", "TA", "SPS" };
        Assert.Multiple(() =>
        {
            for (int k = 0; k < ours.Length; k++)
                Assert.That(MiRnaAnalyzer.RoundPerl3(ours[k]), Is.EqualTo(perlColumns[k]).Within(1e-12),
                    $"{names[k]} contribution (perl prints {perlColumns[k]:F3})");
            Assert.That(c.ContextScore, Is.EqualTo(perlScore).Within(1e-12), "context++ score column");
            Assert.That(c.OmittedFeatures, Has.None.StartsWith("SA"), "supplied plfold ⇒ SA included");
        });
    }

    // The reported score sums the 3-dp contributions: for the 7mer-A1 site the unrounded sum is
    // -0.2580 but TargetScan prints -0.257 (= Σ of the 3-dp columns above).
    [Test]
    public void ContextScore_SumsRoundedContributions_NotTheRoundedPartial()
    {
        var c = ScoreTargetSiteContextPlusPlus(Utr, Let7a, Site(TargetSiteType.Seed7merA1),
            PerlInputs(TargetSiteType.Seed7merA1, 0.0005159016));
        Assert.That(MiRnaAnalyzer.RoundPerl3(c.ContextScorePartial), Is.EqualTo(-0.258).Within(1e-12));
        Assert.That(c.ContextScore, Is.EqualTo(-0.257).Within(1e-12));
    }

    // $maxContextScore{4} = 0 for a 6mer. With SPS at its max (-3.33), TA 4.4, ORF 10 000 nt and plfold
    // 0.5 the 6mer contributions sum above 0: -0.079 -0.005 -0.013 +0.010 +0.006 +0.012 -0.005 -0.030
    // +0.037 +0 -0.010 +0.096 +0.035 = +0.054 → TargetScan reports the ceiling 0.
    [Test]
    public void ContextScore_AboveSiteTypeCeiling_IsCapped()
    {
        var c = ScoreTargetSiteContextPlusPlus(Utr, Let7a, Site(TargetSiteType.Seed6mer),
            new ContextPlusPlusInputs(Sps: -3.33, Ta: 4.4, OrfLength: 10000, Orf8mCount: 0, PlfoldUnpairedProbability: 0.5));
        Assert.That(c.ContextScorePartial, Is.EqualTo(0.054).Within(5e-4), "uncapped sum");
        Assert.That(c.ContextScore, Is.EqualTo(0.0), "6mer ceiling $maxContextScore{4} = 0");
    }

    [TestCase(TargetSiteType.Seed8mer, -0.03)]
    [TestCase(TargetSiteType.Seed7merM8, -0.02)]
    [TestCase(TargetSiteType.Seed7merA1, -0.01)]
    [TestCase(TargetSiteType.Seed6mer, 0.0)]
    public void ContextScore_NeverAboveTargetScanCeiling(TargetSiteType type, double ceiling)
    {
        // A site of each type far from both ends; inputs pushed to their most positive values.
        string core = GetReverseComplement(Let7a.SeedSequence); // CUACCUC
        string site = type switch
        {
            TargetSiteType.Seed8mer => core + "A",
            TargetSiteType.Seed7merM8 => core + "G",
            TargetSiteType.Seed7merA1 => "G" + core[1..] + "A",
            _ => "G" + core[1..] + "G"
        };
        string utr = new string('G', 40) + site + new string('G', 40);
        var s = FindTargetSites(utr, Let7a, 0.0).Single(x => x.Type == type);
        var c = ScoreTargetSiteContextPlusPlus(utr, Let7a, s,
            new ContextPlusPlusInputs(Sps: 0, Ta: 6, OrfLength: 1e6, Orf8mCount: 0, PlfoldUnpairedProbability: 1.0));
        Assert.That(c.ContextScore, Is.LessThanOrEqualTo(ceiling));
    }

    // readTargets: "if ( $utrStart < $MIN_DIST_TO_CDS )" (= 15) → "too_close", no score.
    [TestCase(13, false)]
    [TestCase(14, true)]
    public void ContextScore_SiteWithin15NtOfCds_IsNull(int start0, bool scored)
    {
        string utr = new string('G', start0) + "CUACCUCA" + new string('G', 60);
        var s = FindTargetSites(utr, Let7a, 0.0).Single(x => x.Type == TargetSiteType.Seed8mer);
        Assert.That(s.Start, Is.EqualTo(start0));
        var c = ScoreTargetSiteContextPlusPlus(utr, Let7a, s);
        Assert.That(c.ContextScore.HasValue, Is.EqualTo(scored),
            $"1-based utrStart {start0 + 1} {(scored ? ">=" : "<")} 15");
    }

    // perl printf("%.3f") — the exact binary value, ties to even.
    [TestCase(0.0125, 0.013)]   // binary 0.01250000000000000069 → up
    [TestCase(0.0625, 0.062)]   // exact tie → even
    [TestCase(0.1875, 0.188)]   // exact tie → even
    [TestCase(-0.0625, -0.062)]
    [TestCase(0.0005, 0.001)]
    [TestCase(-0.0004999, 0.0)]
    public void RoundPerl3_MatchesPerlSprintf(double x, double perl)
        => Assert.That(MiRnaAnalyzer.RoundPerl3(x), Is.EqualTo(perl).Within(1e-15));

    // getLocalAU_contribution rounds the fraction to 3 dp before scaling. 8mer with 2 A/U upstream
    // at i = 0,1 of 30 upstream G-free positions and no downstream bases:
    // fraction = (1 + 1/2) / H30 = 1.5 / 3.9949871309203906 = 0.375470… → 0.375;
    // contribution = -0.254 × (0.375 - 0.308) / (0.814 - 0.308) = -0.033632411067193676.
    [Test]
    public void LocalAuContribution_FractionRoundedTo3Decimals()
    {
        string mrna = new string('G', 28) + "AU" + "CUACCUCA";
        double expected = -0.254 * (0.375 - 0.308) / (0.814 - 0.308);
        Assert.That(MiRnaAnalyzer.LocalAuContribution(mrna, 30, 37, TargetSiteType.Seed8mer),
            Is.EqualTo(expected).Within(1e-12));
    }

    // getOffset6merSites counts `while ($utr =~ /$siteToFind/gi)` — non-overlapping. miRNA nt 3-8 =
    // GUGUGU → offset-6mer target ACACAC; perl counts 2 in ACACACACACAC and 1 in ACACACACACGG
    // (perl -e checked), where an overlapping scan would give 4 and 3.
    [TestCase("ACACACACACAC", 2)]
    [TestCase("ACACACACACGG", 1)]
    public void Off6mContribution_CountsNonOverlappingLikePerlRegex(string mrna, int perlCount)
        => Assert.That(MiRnaAnalyzer.Off6mContribution(mrna, "UAGUGUGUAAAAAAAAAAAAAA", TargetSiteType.Seed8mer),
            Is.EqualTo(-0.020 * perlCount).Within(1e-12));

    // getSA_contribution: "if (! $plfold || $plfold eq "NA") { return 0; }".
    [Test]
    public void SaContribution_ZeroPlfold_ContributesZeroAndIsOmitted()
    {
        double sa = MiRnaAnalyzer.SaContribution(Utr, 30, TargetSiteType.Seed8mer, 0.0, out bool included);
        Assert.That(sa, Is.EqualTo(0.0));
        Assert.That(included, Is.False);
    }

    // `grep -A7 -P '^$utrStart\t' _lunp | tail -1`: a 7mer-m8 flush with the UTR 3' end (utrStart = n-6)
    // has only 6 rows after it, so the perl reads the LAST row (n): the 14 nt ending at the UTR end.
    [Test]
    public void SaContribution_SiteNear3PrimeEnd_ReadsLastRowLikeGrepTail()
    {
        string utr = Utr + "CUACCUC"; // 7mer-m8 occupying the last 7 nt
        int n = utr.Length;
        var s = FindTargetSites(utr, Let7a, 0.0).Single(x => x.Start == n - 7);
        Assert.That(s.Type, Is.EqualTo(TargetSiteType.Seed7merM8));

        double sa = MiRnaAnalyzer.SaContribution(utr, s.Start, s.Type, out bool included);
        Assert.That(included, Is.True, "row utrStart+7 = n+1 does not exist → tail -1 gives row n");

        // A 7mer-m8 starting at n-8 has row (n-8)+1+7 = n: the same window [n-14 .. n-1].
        double pAtLastRow = MiRnaAnalyzer.SaContribution(utr, n - 8, s.Type, out _);
        Assert.That(sa, Is.EqualTo(pAtLastRow).Within(1e-12));
    }

    // ── B13 F33: TA / SPS from the bundled TargetScan TA_SPS_by_seed_region.txt ─────────────────────

    // Every row of the kit's own output (TargetScan7_context_scores.zip, targetscan_70_context_scores_output.txt,
    // 675 rows → 158 distinct (miRNA, site type, TA, SPS)): the bundled table + the Agarwal scaling give the
    // printed TA and SPS contributions.
    [Test]
    public void KitOutput_AllRows_TaAndSpsContributions_FromBundledTable_EqualTargetScan()
    {
        Assert.That(MiRnaContextTargetScanKitData.TaSpsRows, Has.Length.EqualTo(158));
        Assert.Multiple(() =>
        {
            foreach (var (mature, type, perlTa, perlSps) in MiRnaContextTargetScanKitData.TaSpsRows)
            {
                string seed = mature.Substring(1, 7);
                Assert.That(TryGetTargetScanTaSps(seed, type, out double ta, out double sps), Is.True, seed);
                Assert.That(MiRnaAnalyzer.RoundPerl3(MiRnaAnalyzer.TaContribution(ta, type)), Is.EqualTo(perlTa).Within(1e-12),
                    $"TA {mature} {type}");
                Assert.That(MiRnaAnalyzer.RoundPerl3(MiRnaAnalyzer.SpsContribution(sps, type)), Is.EqualTo(perlSps).Within(1e-12),
                    $"SPS {mature} {type}");
            }
        });
    }

    // Rows copied by hand from TA_SPS_by_seed_region.txt (seed, SPS 8mer/7mer-m8, SPS 7mer-1a/6mer, TA): first,
    // last and three inner rows, so the base-4 indexing and the column order are both checked.
    [TestCase("AAAAAAA", -0.590, 0.340, 4.371)]
    [TestCase("AAAAAAU", -0.760, 0.340, 4.301)]
    [TestCase("ACGUACG", -8.230, -5.870, 2.362)]
    [TestCase("CCCCCCC", -15.470, -12.210, 3.298)]
    [TestCase("GAGGUAG", -9.250, -6.720, 3.393)]
    [TestCase("UGUGCUU", -7.900, -6.970, 3.523)]
    [TestCase("UUUUUUU", -0.590, 0.340, 4.075)]
    public void BundledTaSpsTable_SpotRows_EqualSourceFile(string seed, double sps1, double sps2, double ta)
    {
        Assert.Multiple(() =>
        {
            foreach (var (type, sps) in new[] { (TargetSiteType.Seed8mer, sps1), (TargetSiteType.Seed7merM8, sps1),
                         (TargetSiteType.Seed7merA1, sps2), (TargetSiteType.Seed6mer, sps2) })
            {
                Assert.That(TryGetTargetScanTaSps(seed, type, out double gotTa, out double gotSps), Is.True);
                Assert.That(gotTa, Is.EqualTo(ta), $"TA {type}");
                Assert.That(gotSps, Is.EqualTo(sps), $"SPS {type}");
            }
        });
    }

    [Test]
    public void TargetScanTaSpsLookup_KeyNormalisation_AndMisses()
    {
        Assert.That(TargetScanTaSpsSeedRegionCount, Is.EqualTo(16384));
        Assert.That(TryGetTargetScanTaSps("gaggtag", TargetSiteType.Seed8mer, out double ta, out double sps), Is.True,
            "T read as U, lower case upper-cased");
        Assert.That((ta, sps), Is.EqualTo((3.393, -9.25)));
        Assert.That(TryGetTargetScanTaSps("GAGGNAG", TargetSiteType.Seed8mer, out _, out _), Is.False, "N has no row");
        Assert.That(TryGetTargetScanTaSps("GAGGUA", TargetSiteType.Seed8mer, out _, out _), Is.False, "6 nt");
        Assert.That(TryGetTargetScanTaSps("GAGGUAG", TargetSiteType.Offset6mer, out _, out _), Is.False,
            "no SPS column for a non-canonical type");
        Assert.Throws<ArgumentNullException>(() => TryGetTargetScanTaSps(null!, TargetSiteType.Seed8mer, out _, out _));
    }

    // The worked example without any caller TA / SPS: the TA and SPS columns and the context++ score equal the perl
    // run (which always reads the table), i.e. the PerlRows above with Sps/Ta left null.
    [TestCaseSource(nameof(PerlRows))]
    public void ScoreTargetSiteContextPlusPlus_NoCallerTaSps_BundledRow_ReproducesTargetScanPerlOutput(
        TargetSiteType type, int utrStart1, double plfold, double[] perlColumns, double perlScore)
    {
        var c = ScoreTargetSiteContextPlusPlus(Utr, Let7a, Site(type),
            new ContextPlusPlusInputs(OrfLength: 1500, Orf8mCount: 1, PlfoldUnpairedProbability: plfold));
        Assert.Multiple(() =>
        {
            Assert.That(MiRnaAnalyzer.RoundPerl3(c.TaContribution), Is.EqualTo(perlColumns[12]).Within(1e-12), "TA");
            Assert.That(MiRnaAnalyzer.RoundPerl3(c.SpsContribution), Is.EqualTo(perlColumns[13]).Within(1e-12), "SPS");
            Assert.That(c.ContextScore, Is.EqualTo(perlScore).Within(1e-12), "context++ score column");
            Assert.That(c.OmittedFeatures, Has.None.StartsWith("SPS").And.None.StartsWith("TA_3UTR"));
        });
    }

    // A caller value still overrides the table, per feature: TA 4.4 (8mer) → 0.222 × (4.4 − 3.113)/(3.865 − 3.113),
    // while SPS keeps the bundled row.
    [Test]
    public void ScoreTargetSiteContextPlusPlus_CallerTa_OverridesBundledRow_SpsStaysBundled()
    {
        var site = Site(TargetSiteType.Seed8mer);
        var bundled = ScoreTargetSiteContextPlusPlus(Utr, Let7a, site);
        var withTa = ScoreTargetSiteContextPlusPlus(Utr, Let7a, site, new ContextPlusPlusInputs(Ta: 4.4));
        Assert.Multiple(() =>
        {
            Assert.That(withTa.TaContribution, Is.EqualTo(0.222 * (4.4 - 3.113) / (3.865 - 3.113)).Within(1e-12));
            Assert.That(withTa.SpsContribution, Is.EqualTo(bundled.SpsContribution));
            Assert.That(bundled.SpsContribution, Is.EqualTo(0.210 * (-9.25 + 11.13) / (-5.52 + 11.13)).Within(1e-12));
        });
    }

    // ── B13 F34: Len_ORF / ORF8m from an ORF sequence, as targetscan_count_8mers.pl ──────────────────

    // The kit's ORF rows (aligned, with gaps, lower case): our length and 8mer counts equal the kit's
    // ORF_Sequences_sample.lengths.txt and ORF_8mer_counts_sample.txt (the script's output, reproduced by
    // re-running it) for all 14 family seeds of miR_Family_info_sample.txt.
    [Test]
    public void KitOrfs_LengthAndEightmerCounts_EqualTargetScanCountScript()
    {
        Assert.Multiple(() =>
        {
            foreach (var (gene, species, orf, length, counts) in MiRnaContextTargetScanKitData.Orfs)
            {
                Assert.That(MiRnaAnalyzer.TargetScanOrfLength(orf), Is.EqualTo(length), $"{gene} {species} length");
                foreach (var (seed, count) in counts)
                    Assert.That(CountTargetScanOrf8merSites(orf, seed), Is.EqualTo(count), $"{gene} {species} {seed}");
            }
        });
    }

    // Every kit output row of those genes: Len_ORF and ORF8m contributions from the ORF sequence equal the
    // printed columns (91 distinct rows; non-zero ORF8m for miR-1 / miR-10 / miR-604 families).
    [Test]
    public void KitOutput_OrfRows_LenOrfAndOrf8mContributions_FromOrfSequence_EqualTargetScan()
    {
        var orfs = MiRnaContextTargetScanKitData.Orfs.ToDictionary(o => (o.Gene, o.Species), o => o.Orf);
        Assert.That(MiRnaContextTargetScanKitData.OrfRows, Has.Length.EqualTo(91));
        Assert.That(MiRnaContextTargetScanKitData.OrfRows.Count(r => r.Orf8m != 0), Is.GreaterThanOrEqualTo(4));
        Assert.Multiple(() =>
        {
            foreach (var (gene, species, mature, type, perlLenOrf, perlOrf8m) in MiRnaContextTargetScanKitData.OrfRows)
            {
                string orf = orfs[(gene, species)];
                double lenOrf = MiRnaAnalyzer.LenOrfContribution(MiRnaAnalyzer.TargetScanOrfLength(orf), type);
                double orf8m = MiRnaAnalyzer.Orf8mContribution(CountTargetScanOrf8merSites(orf, mature.Substring(1, 7)), type);
                Assert.That(MiRnaAnalyzer.RoundPerl3(lenOrf), Is.EqualTo(perlLenOrf).Within(1e-12), $"Len_ORF {gene} {species} {mature} {type}");
                Assert.That(MiRnaAnalyzer.RoundPerl3(orf8m), Is.EqualTo(perlOrf8m).Within(1e-12), $"ORF8m {gene} {species} {mature} {type}");
            }
        });
    }

    // Synthetic ORF for the worked example (gaps, lower case, DNA alphabet, one 8mer CUACCUCA and one 7mer
    // CUACCUCG that must not count). targetscan_count_8mers.pl on it: length 1500, let-7 count 1 — the inputs of
    // the perl run above — and targetscan_70_context_scores.pl on those files prints the same rows. So with only
    // the ORF sequence + RNAplfold value supplied (no TA / SPS / ORF numbers) every column equals the perl.
    private static readonly string WorkedExampleOrf =
        "ATG" + string.Concat(Enumerable.Repeat("GCTGCA", 200)) + "--ctacctca.." + "CUACCUCG"
        + string.Concat(Enumerable.Repeat("GGCAAT", 46)) + "GCTAA";

    [TestCaseSource(nameof(PerlRows))]
    public void ScoreTargetSiteContextPlusPlus_OrfSequenceOnly_ReproducesTargetScanPerlOutput(
        TargetSiteType type, int utrStart1, double plfold, double[] perlColumns, double perlScore)
    {
        var c = ScoreTargetSiteContextPlusPlus(Utr, Let7a, Site(type),
            new ContextPlusPlusInputs(OrfSequence: WorkedExampleOrf, PlfoldUnpairedProbability: plfold));
        double[] ours =
        {
            c.Intercept, c.ThreePrimePairingContribution, c.LocalAuContribution, c.MinDistContribution,
            c.SRna1Contribution, c.SRna8Contribution, c.Site8Contribution, c.Len3UtrContribution,
            c.SaContribution, c.LenOrfContribution, c.Orf8mContribution, c.Off6mContribution,
            c.TaContribution, c.SpsContribution
        };
        Assert.Multiple(() =>
        {
            Assert.That(MiRnaAnalyzer.TargetScanOrfLength(WorkedExampleOrf), Is.EqualTo(1500), "count script length");
            Assert.That(CountTargetScanOrf8merSites(WorkedExampleOrf, "GAGGUAG"), Is.EqualTo(1), "count script 8mers");
            for (int k = 0; k < ours.Length; k++)
                Assert.That(MiRnaAnalyzer.RoundPerl3(ours[k]), Is.EqualTo(perlColumns[k]).Within(1e-12), $"column {k}");
            Assert.That(c.ContextScore, Is.EqualTo(perlScore).Within(1e-12), "context++ score column");
            Assert.That(c.OmittedFeatures, Is.EquivalentTo(new[] { c.OmittedFeatures.Single(f => f.StartsWith("PCT")) }),
                "only PCT (no conservation input) is residual");
        });
    }

    // A second ORF (lower case, 3 let-7 8mers incl. mixed case / T and gap-split context, length 2535 after gap
    // removal): count script → 2535 / 3; perl context++ run on those files: ORF length, ORF 8mer and score columns.
    private static readonly string SecondOrf =
        "atg" + string.Concat(Enumerable.Repeat("gcugca", 300)) + "cuaccuca" + new string('-', 10)
        + string.Concat(Enumerable.Repeat("ggcaat", 50)) + "CTACctcA" + "ccc.CUACCUCA.gg"
        + string.Concat(Enumerable.Repeat("ACGT", 100)) + "taa";

    [TestCase(TargetSiteType.Seed8mer, 0.009772092, 0.131, -0.354, -0.929)]
    [TestCase(TargetSiteType.Seed6mer, 0.008715454, 0.019, -0.180, -0.231)]
    [TestCase(TargetSiteType.Seed7merA1, 0.0005159016, 0.042, -0.174, -0.358)]
    public void ScoreTargetSiteContextPlusPlus_SecondOrf_EqualsCountScriptThenPerl(
        TargetSiteType type, double plfold, double perlLenOrf, double perlOrf8m, double perlScore)
    {
        var c = ScoreTargetSiteContextPlusPlus(Utr, Let7a, Site(type),
            new ContextPlusPlusInputs(OrfSequence: SecondOrf, PlfoldUnpairedProbability: plfold));
        Assert.Multiple(() =>
        {
            Assert.That(MiRnaAnalyzer.TargetScanOrfLength(SecondOrf), Is.EqualTo(2535));
            Assert.That(CountTargetScanOrf8merSites(SecondOrf, "GAGGUAG"), Is.EqualTo(3));
            Assert.That(MiRnaAnalyzer.RoundPerl3(c.LenOrfContribution), Is.EqualTo(perlLenOrf).Within(1e-12));
            Assert.That(MiRnaAnalyzer.RoundPerl3(c.Orf8mContribution), Is.EqualTo(perlOrf8m).Within(1e-12));
            Assert.That(c.ContextScore, Is.EqualTo(perlScore).Within(1e-12));
        });
    }

    // `while ($sequence =~ /$siteToFind/gi)` is non-overlapping: seed AAAAAAU → site AUUUUUUA; in
    // GG AUUUUUUAUUUUUUA cc auuuuuua ATTTTTTAG the script prints 3 (an overlapping scan would give 4).
    [Test]
    public void CountTargetScanOrf8merSites_NonOverlappingCaseInsensitive_LikeCountScript()
    {
        Assert.That(CountTargetScanOrf8merSites("GGAUUUUUUAUUUUUUAccauuuuuuaATTTTTTAG", "AAAAAAU"), Is.EqualTo(3));
        Assert.That(MiRnaAnalyzer.TargetScanOrfLength("GGAUUUUUUAUUUUUUAccauuuuuuaATTTTTTAG"), Is.EqualTo(36));
        Assert.Throws<ArgumentException>(() => CountTargetScanOrf8merSites("ACGU", "GAGGUA"));
    }

    // Caller numbers win over the ORF sequence (per feature).
    [Test]
    public void ScoreTargetSiteContextPlusPlus_CallerOrfNumbers_OverrideOrfSequence()
    {
        var site = Site(TargetSiteType.Seed8mer);
        var c = ScoreTargetSiteContextPlusPlus(Utr, Let7a, site,
            new ContextPlusPlusInputs(OrfLength: 1500, OrfSequence: SecondOrf));
        Assert.Multiple(() =>
        {
            Assert.That(c.LenOrfContribution, Is.EqualTo(MiRnaAnalyzer.LenOrfContribution(1500, TargetSiteType.Seed8mer)));
            Assert.That(c.Orf8mContribution, Is.EqualTo(-0.118 * 3).Within(1e-12), "count from the sequence");
        });
    }
}
