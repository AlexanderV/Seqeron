// MIRNA-TARGET-001 / WP-6 — Grimson et al. (2007) context score as TargetScan 5.0 computes it
// (targetscan_50_context_scores.pl, run on the kit sample files; numbers in MiRnaGrimsonKitData),
// FindTargetSites.Score / TargetSite.ContextScore, AnalyzeTargetContext, centered sites (Shin et al. 2010)
// and the offset-6mer suppression independence from minScore.

using static Seqeron.Genomics.Annotation.MiRnaAnalyzer;

namespace Seqeron.Genomics.Tests.Unit.Annotation;

[TestFixture]
[Category("MIRNA-TARGET-001")]
public class MiRnaAnalyzer_GrimsonContext_Tests
{
    private static readonly MiRna Let7a = CreateMiRna("let-7a", "UGAGGUAGUAGGUUGUAUAGUU");

    private static TargetSiteType ToType(string t) => t switch
    {
        "8mer" => TargetSiteType.Seed8mer,
        "m8" => TargetSiteType.Seed7merM8,
        _ => TargetSiteType.Seed7merA1
    };

    private static (string Utr, string Mature) Fetch(string gene, string species, string mirna) =>
        (MiRnaGrimsonKitData.Utrs.Single(u => u.Gene == gene && u.Species == species).Utr,
         MiRnaGrimsonKitData.MiRnas.Single(m => m.Id == mirna).Mature);

    // TargetScan 5.0 coordinates exclude the A1 of an 8mer; TargetSite coordinates include it.
    private static (int Start, int End) LibraryCoords(string type, int utrStart, int utrEnd) =>
        (utrStart - 1, utrEnd - 1 + (type == "8mer" ? 1 : 0));

    [Test]
    public void KitRows_AllContributionsAndContextScore_EqualTargetScan50()
    {
        Assert.That(MiRnaGrimsonKitData.Rows.Length, Is.EqualTo(400));
        foreach (var r in MiRnaGrimsonKitData.Rows)
        {
            var (utr, mature) = Fetch(r.Gene, r.Species, r.MiRna);
            var (s, e) = LibraryCoords(r.Type, r.UtrStart, r.UtrEnd);
            var g = ScoreTargetSiteGrimson(utr, mature, s, e, ToType(r.Type))!.Value;
            string id = $"{r.Gene}/{r.Species}/{r.MiRna}/{r.Type}@{r.UtrStart}";
            Assert.That(g.ThreePrimePairingContribution, Is.EqualTo(r.ThreeP).Within(1e-9), id + " 3'");
            Assert.That(g.LocalAuContribution, Is.EqualTo(r.Au).Within(1e-9), id + " AU");
            Assert.That(g.PositionContribution, Is.EqualTo(r.Pos).Within(1e-9), id + " pos");
            Assert.That(g.ContextScore!.Value, Is.EqualTo(r.Context).Within(1e-9), id + " context");
        }
    }

    [Test]
    public void KitRows_SiteTypeContributions_Are_Minus0310_0161_0099()
    {
        var r8 = MiRnaGrimsonKitData.Rows.First(r => r.Type == "8mer");
        var (utr, mature) = Fetch(r8.Gene, r8.Species, r8.MiRna);
        var (s, e) = LibraryCoords("8mer", r8.UtrStart, r8.UtrEnd);
        Assert.That(ScoreTargetSiteGrimson(utr, mature, s, e, TargetSiteType.Seed8mer)!.Value.SiteTypeContribution, Is.EqualTo(-0.310));
        var rm = MiRnaGrimsonKitData.Rows.First(r => r.Type == "m8");
        (utr, mature) = Fetch(rm.Gene, rm.Species, rm.MiRna);
        (s, e) = LibraryCoords("m8", rm.UtrStart, rm.UtrEnd);
        Assert.That(ScoreTargetSiteGrimson(utr, mature, s, e, TargetSiteType.Seed7merM8)!.Value.SiteTypeContribution, Is.EqualTo(-0.161));
        var ra = MiRnaGrimsonKitData.Rows.First(r => r.Type == "1a");
        (utr, mature) = Fetch(ra.Gene, ra.Species, ra.MiRna);
        (s, e) = LibraryCoords("1a", ra.UtrStart, ra.UtrEnd);
        Assert.That(ScoreTargetSiteGrimson(utr, mature, s, e, TargetSiteType.Seed7merA1)!.Value.SiteTypeContribution, Is.EqualTo(-0.099));
    }

    // The one kit row with UTR start < 32 (FNDC3A 9615 cfa-miR-101, 7mer-A1 at 21): TargetScan 5.0 passes
    // substr(utr, 21-32 = -11, 30), a negative offset that reads from the END of the UTR (printed local AU
    // -0.095, context -0.2634); TargetScan fixed that in 5.3 by clipping at the UTR start. Independent python
    // derivation of the clipped value: 19 upstream nt (weights 1/(i+2)) + 30 downstream (1/(i+2)) →
    // fraction 0.9127687245, AU = -0.42*0.9127687245+0.236 = -0.147363 → -0.147; context = -0.099-0.147-0.031-0.0384.
    [Test]
    public void KitRow_NearUtrStart_ClipsInsteadOfNegativeSubstr()
    {
        var r = MiRnaGrimsonKitData.SubstrDefectRow;
        var (utr, mature) = Fetch(r.Gene, r.Species, r.MiRna);
        var (s, e) = LibraryCoords(r.Type, r.UtrStart, r.UtrEnd);
        var g = ScoreTargetSiteGrimson(utr, mature, s, e, TargetSiteType.Seed7merA1)!.Value;
        Assert.That(g.LocalAuFraction, Is.EqualTo(0.9127687245327373).Within(1e-12));
        Assert.That(g.LocalAuContribution, Is.EqualTo(-0.147).Within(1e-12));
        Assert.That(g.PositionContribution, Is.EqualTo(r.Pos).Within(1e-12));
        Assert.That(g.ThreePrimePairingContribution, Is.EqualTo(r.ThreeP).Within(1e-9));
        Assert.That(g.ContextScore!.Value, Is.EqualTo(-0.099 - 0.147 - 0.031 + r.ThreeP).Within(1e-9));
    }

    // Worked example (hand derivation, TargetScan constants): let-7a 8mer in poly-G flanks, 100 nt each side.
    //   utr = G*100 + "CUACCUC" + "A" + G*100   (8mer = RC of nt2-8 + A1; TargetSite coordinates 100..107)
    //   local AU: flanks all G except A1 which is inside the site → 0 (A1 of the 8mer is part of the 8-nt site here,
    //   but TargetScan 5.0 coordinates end 1 nt earlier, so A1 is the first downstream nt, weight 1/2, an A):
    //   up Σ 1/(i+1), i<30 = H30 = 3.994987; down Σ 1/(i+2), i<30 = H31−1 = 3.0...; fraction = (1/2)/(H30 + H31 − 1).
    [Test]
    public void WorkedExample_8merPolyG_HandDerived()
    {
        string utr = new string('G', 100) + "CUACCUC" + "A" + new string('G', 100);
        var g = ScoreTargetSiteGrimson(utr, Let7a.Sequence, 100, 107, TargetSiteType.Seed8mer)!.Value;

        double h30 = 0, down = 0;
        for (int i = 1; i <= 30; i++) h30 += 1.0 / i;
        for (int i = 0; i < 30; i++) down += 1.0 / (i + 2);
        double fraction = 0.5 / (h30 + down);
        Assert.That(g.LocalAuFraction, Is.EqualTo(fraction).Within(1e-12));
        Assert.That(g.LocalAuContribution, Is.EqualTo(Math.Round(-0.64 * fraction + 0.365, 3)).Within(1e-12));
        // distance to nearest end: min(utrStart-1 = 100, len - utrEnd(5.0) = 208 - 107 = 101) = 100
        Assert.That(g.PositionContribution, Is.EqualTo(Math.Round(100 * 0.000172 - 0.07, 3)).Within(1e-12));
        Assert.That(g.SiteTypeContribution, Is.EqualTo(-0.310));
        Assert.That(g.ContextScore!.Value,
            Is.EqualTo(g.SiteTypeContribution + g.LocalAuContribution + g.PositionContribution + g.ThreePrimePairingContribution).Within(1e-12));
    }

    [Test]
    public void SiteWithin15NtOfStopCodon_TooClose_ContextScoreNull_ContributionsKept()
    {
        string utr = "GGGGGGG" + "CUACCUC" + "A" + new string('G', 60);   // 1-based start 8 < 15
        var g = ScoreTargetSiteGrimson(utr, Let7a.Sequence, 7, 14, TargetSiteType.Seed8mer)!.Value;
        Assert.That(g.ContextScore, Is.Null);
        Assert.That(g.SiteTypeContribution, Is.EqualTo(-0.310));
    }

    [Test]
    public void UnsupportedTypesAndBadInput_ReturnNull()
    {
        const string utr = "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGCUACCUGGGGGGGGGGGGGGGGGGGGGGGGGG";
        foreach (var t in new[] { TargetSiteType.Seed6mer, TargetSiteType.Offset6mer, TargetSiteType.Centered })
            Assert.That(ScoreTargetSiteGrimson(utr, Let7a.Sequence, 40, 45, t), Is.Null, t.ToString());
        Assert.That(ScoreTargetSiteGrimson("", Let7a.Sequence, 0, 1, TargetSiteType.Seed8mer), Is.Null);
        Assert.That(ScoreTargetSiteGrimson(utr, Let7a.Sequence, 40, 500, TargetSiteType.Seed8mer), Is.Null);
        Assert.That(ScoreTargetSiteGrimson(utr, Let7a.Sequence, 45, 40, TargetSiteType.Seed8mer), Is.Null);
    }

    [Test]
    public void FindTargetSites_Score_IsNegatedContextScore_ClampedAtZero_UnscoredTypesZero()
    {
        string utr = new string('A', 60) + "CUACCUC" + "A" + new string('A', 60) + "UACCUC" + new string('G', 60);
        var sites = FindTargetSites(utr, Let7a, minScore: 0.0).ToList();
        var s8 = sites.Single(s => s.Type == TargetSiteType.Seed8mer);
        var g = ScoreTargetSiteGrimson(utr, Let7a.Sequence, s8.Start, s8.End, TargetSiteType.Seed8mer)!.Value;
        Assert.That(s8.ContextScore, Is.EqualTo(g.ContextScore));
        Assert.That(s8.Score, Is.EqualTo(Math.Max(0.0, -g.ContextScore!.Value)).Within(1e-12));
        Assert.That(s8.Score, Is.GreaterThan(0.3), "an AU-rich 8mer is predicted strongly repressed");

        var s6 = sites.Single(s => s.Type == TargetSiteType.Seed6mer);
        Assert.That(s6.Score, Is.EqualTo(0.0));
        Assert.That(s6.ContextScore, Is.Null, "TargetScan 5.0 has no 6mer coefficients");
    }

    [Test]
    public void FindTargetSites_MinScore_FiltersOnPredictedRepression()
    {
        string utr = new string('A', 60) + "CUACCUC" + "A" + new string('A', 60);
        double score = FindTargetSites(utr, Let7a, 0.0).Single(s => s.Type == TargetSiteType.Seed8mer).Score;
        Assert.That(FindTargetSites(utr, Let7a, score).Any(s => s.Type == TargetSiteType.Seed8mer), Is.True, "threshold is inclusive");
        Assert.That(FindTargetSites(utr, Let7a, score + 1e-9).Any(s => s.Type == TargetSiteType.Seed8mer), Is.False);
    }

    // AnalyzeTargetContext: ContextScore = round3(slope_AU·AuContent+intercept_AU) + round3(slope_pos·d+intercept_pos).
    [Test]
    public void AnalyzeTargetContext_ContextScore_IsAuPlusPositionContribution_OfTargetScan50()
    {
        string utr = new string('A', 100) + "CUACCUC" + new string('U', 100);      // 7mer-m8 at 100..106
        var (au, nearStart, nearEnd, ctx) = AnalyzeTargetContext(utr, 100, 106);
        Assert.That(au, Is.EqualTo(1.0).Within(1e-12));
        Assert.That(nearStart, Is.False);
        Assert.That(nearEnd, Is.False);
        // m8: -0.5*1 + 0.269 = -0.231 ; distance min(100, 207-107 = 100) -> 100*0.000091-0.037 = -0.0279 -> -0.028
        Assert.That(ctx, Is.EqualTo(-0.231 + -0.028).Within(1e-12));
    }

    [Test]
    public void AnalyzeTargetContext_8mer_UsesEightmerCoefficients()
    {
        string utr = new string('G', 100) + "CUACCUC" + "A" + new string('G', 100);   // 8mer incl. A1 at 100..107
        var (au, _, _, ctx) = AnalyzeTargetContext(utr, 100, 107, siteType: TargetSiteType.Seed8mer);
        var g = ScoreTargetSiteGrimson(utr, Let7a.Sequence, 100, 107, TargetSiteType.Seed8mer)!.Value;
        Assert.That(au, Is.EqualTo(g.LocalAuFraction).Within(1e-12));
        Assert.That(ctx, Is.EqualTo(g.LocalAuContribution + g.PositionContribution).Within(1e-12));
    }

    [Test]
    public void AnalyzeTargetContext_NearStart_ScoreZero_AndUnsupportedTypeScoreZero()
    {
        string utr = new string('A', 10) + "CUACCUC" + new string('U', 100);
        Assert.That(AnalyzeTargetContext(utr, 10, 16).ContextScore, Is.EqualTo(0.0));
        Assert.That(AnalyzeTargetContext(utr, 10, 16).NearStart, Is.True);
        string utr2 = new string('A', 100) + "CUACCU" + new string('U', 100);
        Assert.That(AnalyzeTargetContext(utr2, 100, 105, siteType: TargetSiteType.Seed6mer).ContextScore, Is.EqualTo(0.0));
    }

    [Test]
    public void AnalyzeTargetContext_PositionCap_At1500()
    {
        string utr = new string('G', 2000) + "CUACCUC" + new string('G', 2000);
        var (_, _, _, ctx) = AnalyzeTargetContext(utr, 2000, 2006);
        var au = AnalyzeTargetContext(utr, 2000, 2006).AuContent;
        // distance capped at 1500: 1500*0.000091-0.037 = 0.0995 -> 0.100 (IEEE, as perl on x86-64)
        Assert.That(ctx, Is.EqualTo(Math.Round(-0.5 * au + 0.269, 3) + 0.100).Within(1e-12));
    }

    // ── centered sites (Shin et al. 2010) ─────────────────────────────────────────────────────────────

    private static string Rc(string rna) => GetReverseComplement(rna);

    // Hand construction: let-7a nt 4..15 = "GGUAGUAGGUUG" (nt1 U,2 G,3 A,4 G,5 G,6 U,7 A,8 G,9 U,10 A,11 G,12 G,13 U,14 U,15 G).
    [Test]
    public void FindCenteredSites_12ntRun_Nt4to15_Found()
    {
        string pair415 = Rc(Let7a.Sequence.Substring(3, 12));
        string utr = "GGGGGGGGGG" + pair415 + "GGGGGGGGGG";
        var site = FindCenteredSites(utr, Let7a).Single();
        Assert.That(site.Type, Is.EqualTo(TargetSiteType.Centered));
        Assert.That(site.Start, Is.EqualTo(10));
        Assert.That(site.End, Is.EqualTo(21));
        Assert.That(site.SeedMatchLength, Is.EqualTo(12));
        Assert.That(site.Score, Is.EqualTo(0.0));
        Assert.That(site.ContextScore, Is.Null);
    }

    [Test]
    public void FindCenteredSites_11ntRuns_Nt4to14_And_Nt5to15_Found()
    {
        // nt 4-14 paired, nt 15 mismatched (G opposite the miRNA's G needs not-C): replace its partner by G.
        string p414 = Rc(Let7a.Sequence.Substring(3, 11));        // nt 4..14 → 11 nt
        string u1 = "AAAAAAAAAA" + p414 + "G" + "AAAAAAAAAA";        // nt 15's partner position (upstream of p414) is 'A' ≠ C
        var s1 = FindCenteredSites(u1, Let7a).Single();
        Assert.That(s1.SeedMatchLength, Is.EqualTo(11));
        Assert.That(s1.End - s1.Start + 1, Is.EqualTo(11));

        string p515 = Rc(Let7a.Sequence.Substring(4, 11));        // nt 5..15 → 11 nt
        string u2 = "AAAAAAAAAA" + p515 + "AAAAAAAAAA";
        var s2 = FindCenteredSites(u2, Let7a).Single();
        Assert.That(s2.SeedMatchLength, Is.EqualTo(11));
    }

    [Test]
    public void FindCenteredSites_10ntRun_NotACenteredSite()
    {
        string p = Rc(Let7a.Sequence.Substring(4, 10));            // nt 5..14 only
        Assert.That(FindCenteredSites("AAAAAAAAAA" + p + "AAAAAAAAAA", Let7a), Is.Empty);
    }

    [Test]
    public void FindCenteredSites_PerfectSeedPresent_IsASeedSiteNotCentered()
    {
        // full complement of nt 2-15 contains the seed: canonical site, not centered.
        string full = Rc(Let7a.Sequence.Substring(1, 14));
        Assert.That(FindCenteredSites("AAAAAAAAAA" + full + "AAAAAAAAAA", Let7a), Is.Empty);
    }

    [Test]
    public void FindCenteredSites_WobbleDoesNotCount_AndShortInputs()
    {
        string pair = Rc(Let7a.Sequence.Substring(3, 12));   // "CAACCUACUACC": target idx i+15-k pairs miRNA nt k
        // turn the C opposite nt 12 (G) into U: a G:U wobble in the middle destroys both the 4-14 and the 5-15 run
        char[] c = pair.ToCharArray();
        Assert.That(c[3], Is.EqualTo('C'));
        c[3] = 'U';
        Assert.That(FindCenteredSites("AAAAAAAAAA" + new string(c) + "AAAAAAAAAA", Let7a), Is.Empty);
        Assert.That(FindCenteredSites("", Let7a), Is.Empty);
        Assert.That(FindCenteredSites("ACGU", Let7a), Is.Empty);
        Assert.That(FindCenteredSites(pair, CreateMiRna("short", "UGAGGUAGUAGGUU")), Is.Empty, "miRNA shorter than nt 15");
    }

    [Test]
    public void FindTargetSites_NeverEmits_Centered()
    {
        string utr = new string('A', 40) + Rc(Let7a.Sequence.Substring(3, 12)) + new string('A', 40) + "CUACCUCA" + new string('A', 40);
        Assert.That(FindTargetSites(utr, Let7a, double.NegativeInfinity)
            .Any(s => s.Type == TargetSiteType.Centered), Is.False);
    }

    // ── A17: offset-6mer suppression independent of minScore ───────────────────────────────────────────

    [Test]
    public void OffsetSixmerSuppression_DoesNotDependOnMinScore()
    {
        // "UACCUC" core + the offset pattern "CUACCU" overlapping its last base: 6mer at 0..5, offset 6mer at 5..10.
        string utr = "UACCUCUACCU" + new string('G', 30);
        var all = FindTargetSites(utr, Let7a, double.NegativeInfinity).ToList();
        Assert.That(all.Any(s => s.Type == TargetSiteType.Offset6mer), Is.False, "offset overlapping a 6mer-core site is suppressed");
        foreach (double min in new[] { -1.0, 0.0, 1e-9, 0.1, 5.0 })
            Assert.That(FindTargetSites(utr, Let7a, min).Any(s => s.Type == TargetSiteType.Offset6mer), Is.False, $"minScore {min}");
    }
}
