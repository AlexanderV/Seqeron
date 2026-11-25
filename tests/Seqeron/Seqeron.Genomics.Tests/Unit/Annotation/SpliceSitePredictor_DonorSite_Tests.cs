using static Seqeron.Genomics.Annotation.SpliceSitePredictor;

namespace Seqeron.Genomics.Tests.Unit.Annotation;

/// <summary>
/// SPLICE-DONOR-001: Canonical test file for Donor (5') Splice Site Detection.
/// Tests <see cref="SpliceSitePredictor.FindDonorSites(string, double, bool)"/> (Shapiro &amp;
/// Senapathy 1987 score), its MaxEntScan mode and U12-type donors, and
/// <see cref="SpliceSitePredictor.ScoreDonorMaxEnt"/>.
/// Evidence: Shapiro &amp; Senapathy (1987) NAR 15:7155 (Table 1 percentages, Methods formula,
/// Table 5 published scores), Sharp &amp; Burge (1997), Burge et al. (1999), Yeo &amp; Burge (2004).
/// </summary>
[TestFixture]
public class SpliceSitePredictor_DonorSite_Tests
{
    // S&S (1987) score = (t - mint) / (maxt - mint) over the 8-nt subsequence -2..+6.
    // Primate Table 1: mint = 13+4+0+0+2+8+5+15 = 47, maxt = 58+78+100+100+57+71+84+47 = 595.
    private const double PrimateRange = 595 - 47;

    #region M1: Canonical GT Consensus Detected — Shapiro & Senapathy (1987)

    [Test]
    public void FindDonorSites_CanonicalGT_ConsensusMotif_ProducesDonorSite()
    {
        // CAG|GUAAGU: every scored position -2..+6 (AG|GUAAGU) carries its highest primate
        // percentage (58, 78, 100, 100, 57, 71, 84, 47) ⇒ t = maxt ⇒ score 1.0.
        // Confidence = (1.0 - 0.5) / (1.0 - 0.5) = 1.0
        const string sequence = "CAGGUAAGU";
        var sites = FindDonorSites(sequence, minScore: 0.3).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(sites, Has.Count.EqualTo(1), "Exactly one GU in the sequence");
            Assert.That(sites[0].Type, Is.EqualTo(SpliceSiteType.Donor));
            Assert.That(sites[0].Score, Is.EqualTo(1.0).Within(1e-12));
            Assert.That(sites[0].Confidence, Is.EqualTo(1.0).Within(1e-12));
            Assert.That(sites[0].Position, Is.EqualTo(3), "Position = index of the G of GU (intron +1)");
            Assert.That(sites[0].Motif, Is.EqualTo("CAGGUAAGU"), "Motif = -3..+6 context, RNA alphabet");
        });
    }

    #endregion

    #region M1b: Published S&S scores — Shapiro & Senapathy (1987) Table 5

    // Table 5 (p. 7168) lists the 5' ss score of real introns, each scored with the subtable of
    // its GENBANK category; the formula reproduces 39 of its 42 values exactly at 1 dp (the
    // other three: primate GC/GTGAGT 75.1 vs 75.4 and CG/GTGAGT 88.4 vs 88.5 — Table 1
    // percentages are rounded; rodent TG/GTACAG 68.3 vs 60.9 — a misprint, its rank >100 is
    // below intron 4's 67.2 at rank 71). Methods worked example: primate AAGTGAGT → 84.3.
    [TestCase("AAGTGAGT", ShapiroSenapathyCategory.Primate, 84.3)]
    [TestCase("GGGTACCA", ShapiroSenapathyCategory.Primate, 60.8)]
    [TestCase("GGGTAAGG", ShapiroSenapathyCategory.Primate, 87.6)]
    [TestCase("AGGTAAGG", ShapiroSenapathyCategory.Primate, 95.4)]
    [TestCase("AAGTAAGT", ShapiroSenapathyCategory.Primate, 87.6)]
    [TestCase("TGGTAAGT", ShapiroSenapathyCategory.Rodent, 92.8)]
    [TestCase("AGGTGAGA", ShapiroSenapathyCategory.Rodent, 89.5)]
    [TestCase("AGGTACAC", ShapiroSenapathyCategory.Rodent, 67.2)]
    [TestCase("GGGTGAGT", ShapiroSenapathyCategory.OtherMammal, 94.8)]
    [TestCase("AGGTTCGC", ShapiroSenapathyCategory.OtherMammal, 75.1)]
    [TestCase("GGGTAGAG", ShapiroSenapathyCategory.OtherMammal, 63.6)]
    [TestCase("AGGTAAGT", ShapiroSenapathyCategory.OtherMammal, 98.3)]
    [TestCase("AGGTACCG", ShapiroSenapathyCategory.OtherVertebrate, 71.5)]
    [TestCase("TTGTAAGT", ShapiroSenapathyCategory.OtherVertebrate, 80.4)]
    [TestCase("AGGTAAAA", ShapiroSenapathyCategory.OtherVertebrate, 80.8)]
    [TestCase("AGGTGAGT", ShapiroSenapathyCategory.Invertebrate, 91.2)]
    [TestCase("ACGTAAGT", ShapiroSenapathyCategory.Invertebrate, 87.8)]
    [TestCase("TTGTGAGT", ShapiroSenapathyCategory.Invertebrate, 71.6)]
    [TestCase("CTGTAAGT", ShapiroSenapathyCategory.Invertebrate, 80.1)]
    [TestCase("AGGTCCGC", ShapiroSenapathyCategory.Plant, 76.5)]
    [TestCase("GGGTATGT", ShapiroSenapathyCategory.Plant, 92.8)]
    [TestCase("TGGTAAGT", ShapiroSenapathyCategory.Plant, 96.7)]
    [TestCase("AGGTACAG", ShapiroSenapathyCategory.Plant, 73.8)]
    public void FindDonorSites_ShapiroSenapathy_ReproducesPublishedTable5Scores(
        string subsequence, ShapiroSenapathyCategory category, double published)
    {
        // The 8-mer is the scored subsequence -2..+6, so its GT sits at index 2.
        var site = FindDonorSites(subsequence, 0.0, false, DonorScoringModel.ShapiroSenapathy, category).Single();

        Assert.Multiple(() =>
        {
            Assert.That(site.Position, Is.EqualTo(2));
            Assert.That(Math.Round(site.Score * 100, 1, MidpointRounding.AwayFromZero), Is.EqualTo(published));
        });
    }

    [Test]
    public void FindDonorSites_ShapiroSenapathy_MethodsWorkedExample_ExactFraction()
    {
        // Methods (p. 7158): primate AAGTGAGT, t = 58+10+100+100+39+71+84+47 = 509 ⇒ (509-47)/548.
        var site = FindDonorSites("AAGTGAGT", minScore: 0.0).Single();
        Assert.That(site.Score, Is.EqualTo(462 / PrimateRange).Within(1e-15));
    }

    [Test]
    public void FindDonorSites_DefaultOverload_IsPrimateShapiroSenapathy()
    {
        const string sequence = "CCAGGUGAGUACGGUAAGCAUUGUAUGG";
        var legacy = FindDonorSites(sequence, 0.0, true).ToList();
        var explicitPrimate = FindDonorSites(sequence, 0.0, true,
            DonorScoringModel.ShapiroSenapathy, ShapiroSenapathyCategory.Primate).ToList();
        Assert.That(legacy, Is.EqualTo(explicitPrimate));
    }

    #endregion

    #region M2-M4: Negative controls

    [Test]
    public void FindDonorSites_NoGU_ReturnsEmpty()
    {
        const string sequence = "AAAAACCCCC";
        var sites = FindDonorSites(sequence, minScore: 0.3).ToList();
        Assert.That(sites, Is.Empty, "No GU dinucleotide → no donor sites");
    }

    [Test]
    public void FindDonorSites_EmptyString_ReturnsEmpty()
    {
        var sites = FindDonorSites("", minScore: 0.3).ToList();
        Assert.That(sites, Is.Empty, "Empty input yields no donor sites");
    }

    [Test]
    public void FindDonorSites_SequenceShorterThanScoredWindow_ReturnsEmpty()
    {
        // The S&S subsequence -2..+6 needs 8 nt; "GUAA" and "AGGUAAG" (7 nt) cannot be scored.
        Assert.Multiple(() =>
        {
            Assert.That(FindDonorSites("GUAA", minScore: 0.0), Is.Empty);
            Assert.That(FindDonorSites("AGGUAAG", minScore: 0.0), Is.Empty);
            Assert.That(FindDonorSites("AGGUAAGU", minScore: 0.0).Single().Position, Is.EqualTo(2));
        });
    }

    #endregion

    #region M5: Strong vs Weak Context — Shapiro & Senapathy (1987)

    [Test]
    public void FindDonorSites_StrongConsensus_ScoresHigherThanWeakContext()
    {
        // Strong: CAG|GUAAGU ⇒ 1.0.
        // Weak:   UUU|GUAAUU, -2..+6 = UU|GUAAUU ⇒ t = 15+8+100+100+57+71+5+47 = 403 ⇒ 356/548.
        var strongSites = FindDonorSites("CAGGUAAGU", minScore: 0.0).ToList();
        var weakSites = FindDonorSites("UUUGUAAUU", minScore: 0.0).ToList();

        Assert.That(strongSites, Has.Count.EqualTo(1));
        Assert.That(weakSites, Has.Count.EqualTo(1));
        Assert.Multiple(() =>
        {
            Assert.That(strongSites[0].Score, Is.EqualTo(1.0).Within(1e-12));
            Assert.That(weakSites[0].Score, Is.EqualTo(356 / PrimateRange).Within(1e-12));
            Assert.That(weakSites[0].Confidence, Is.EqualTo((356 / PrimateRange - 0.5) / 0.5).Within(1e-12));
            Assert.That(strongSites[0].Score, Is.GreaterThan(weakSites[0].Score));
        });
    }

    #endregion

    #region M6-M7: GC Non-Canonical Donor — Burge et al. (1999)

    [Test]
    public void FindDonorSites_GC_NonCanonical_DetectedWhenEnabled()
    {
        // CAG|GCAAGU: +2 C carries 0 % (U = 100 %) ⇒ t = 595 - 100 = 495 ⇒ 448/548.
        var sites = FindDonorSites("CAGGCAAGU", minScore: 0.0, includeNonCanonical: true).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(sites, Has.Count.EqualTo(1), "GC donor found; no GU/AU ⇒ no U12 candidate");
            Assert.That(sites[0].Type, Is.EqualTo(SpliceSiteType.Donor));
            Assert.That(sites[0].Score, Is.EqualTo(448 / PrimateRange).Within(1e-12));
            Assert.That(sites[0].Position, Is.EqualTo(3));
        });
    }

    [Test]
    public void FindDonorSites_GC_NotDetected_WhenCanonicalOnly()
    {
        var sites = FindDonorSites("CAGGCAACC", minScore: 0.0, includeNonCanonical: false).ToList();
        Assert.That(sites, Is.Empty, "GC donors are reported only with includeNonCanonical = true");
    }

    [Test]
    public void FindDonorSites_GC_Donor_ScoresLowerThanEquivalentGT()
    {
        var gt = FindDonorSites("CAGGUAAGU", minScore: 0.0).Single();
        var gc = FindDonorSites("CAGGCAAGU", minScore: 0.0, includeNonCanonical: true).Single();
        Assert.Multiple(() =>
        {
            Assert.That(gt.Score, Is.EqualTo(1.0).Within(1e-12));
            Assert.That(gc.Score, Is.EqualTo(448 / PrimateRange).Within(1e-12));
            Assert.That(gc.Score, Is.LessThan(gt.Score));
        });
    }

    #endregion

    #region M8-M9: Alphabet and case

    [Test]
    public void FindDonorSites_DNA_T_ProducesSameResultsAsRNA_U()
    {
        var dnaSites = FindDonorSites("CAGGTAAGT", minScore: 0.3).ToList();
        var rnaSites = FindDonorSites("CAGGUAAGU", minScore: 0.3).ToList();
        Assert.That(dnaSites, Is.EqualTo(rnaSites), "T is read as U; Motif is reported in RNA alphabet");
    }

    [Test]
    public void FindDonorSites_LowercaseInput_FindsSite()
    {
        var sites = FindDonorSites("cagguaagu", minScore: 0.3).ToList();
        Assert.That(sites, Is.EqualTo(FindDonorSites("CAGGUAAGU", minScore: 0.3).ToList()));
    }

    #endregion

    #region M10: Multiple GU Sites All Detected — Full Sequence Scan

    [Test]
    public void FindDonorSites_MultipleSites_AllDetected()
    {
        // GU at 3 (CAG|GUAAGU ⇒ 1.0), at 7 (-2..+6 = AA|GUUUUU ⇒ t = 58+10+200+2+9+5+47 = 331
        // ⇒ 284/548) and at 26 (⇒ 1.0; its +6 is the last base).
        const string sequence = "CAGGUAAGUUUUUUUUUUUUUUUCAGGUAAGU";
        var sites = FindDonorSites(sequence, minScore: 0.3).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(sites.Select(s => s.Position), Is.EqualTo(new[] { 3, 7, 26 }));
            Assert.That(sites, Has.All.Property(nameof(SpliceSite.Type)).EqualTo(SpliceSiteType.Donor));
            Assert.That(sites[0].Score, Is.EqualTo(1.0).Within(1e-12));
            Assert.That(sites[1].Score, Is.EqualTo(284 / PrimateRange).Within(1e-12));
            Assert.That(sites[2].Score, Is.EqualTo(1.0).Within(1e-12));
        });
    }

    #endregion

    #region S1-S5: Invariants

    [Test]
    public void FindDonorSites_AllScoresAndConfidences_InZeroOneRange()
    {
        const string sequence = "CAGGUAAGUAAAGUAAGUCCCCAAAGGUAAGUAUAUCCUUGU";
        var sites = FindDonorSites(sequence, minScore: 0.0, includeNonCanonical: true).ToList();

        Assert.That(sites, Is.Not.Empty);
        foreach (var site in sites)
        {
            Assert.That(site.Score, Is.InRange(0.0, 1.0));
            Assert.That(site.Confidence, Is.InRange(0.0, 1.0));
            Assert.That(site.Motif, Is.Not.Null.And.Not.Empty);
        }
    }

    [Test]
    public void FindDonorSites_HigherMinScore_IsSubsetOfLowerThreshold()
    {
        const string sequence = "CAGGUAAGUAAAGUAAGUCCCCAAAGGUAAGU";
        var low = FindDonorSites(sequence, minScore: 0.2).ToList();
        var high = FindDonorSites(sequence, minScore: 0.8).ToList();
        Assert.That(high, Is.SubsetOf(low));
        Assert.That(high.All(s => s.Score >= 0.8), Is.True);
    }

    [Test]
    public void FindDonorSites_NullInput_ReturnsEmpty()
    {
        Assert.That(FindDonorSites(null!, minScore: 0.3), Is.Empty);
    }

    [Test]
    public void FindDonorSites_UnknownModelOrCategory_Throws()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentOutOfRangeException>(() =>
                FindDonorSites("CAGGUAAGU", 0.0, false, (DonorScoringModel)99));
            Assert.Throws<ArgumentOutOfRangeException>(() =>
                FindDonorSites("CAGGUAAGU", 0.0, false, DonorScoringModel.ShapiroSenapathy, (ShapiroSenapathyCategory)99));
        });
    }

    #endregion

    #region C1: U12-type donors — /RTATCCTY (Sharp & Burge 1997; Burge et al. 1999)

    [Test]
    public void FindDonorSites_U12_AT_Donor_Detected_WhenNonCanonical()
    {
        // AT-AC U12 intron: /AUAUCCUU matches all 8 positions of /RUAUCCUY ⇒ 1.0.
        var sites = FindDonorSites("CAGAUAUCCUUA", minScore: 0.0, includeNonCanonical: true).ToList();

        Assert.That(sites, Has.Count.EqualTo(1));
        Assert.Multiple(() =>
        {
            Assert.That(sites[0].Type, Is.EqualTo(SpliceSiteType.U12Donor));
            Assert.That(sites[0].Position, Is.EqualTo(3));
            Assert.That(sites[0].Score, Is.EqualTo(1.0).Within(1e-12));
            Assert.That(sites[0].Motif, Is.EqualTo("CAGAUAUCCUU"), "U12 motif = -3..+8");
        });
        Assert.That(FindDonorSites("CAGAUAUCCUUA", minScore: 0.0), Is.Empty,
            "U12 donors are reported only with includeNonCanonical = true");
    }

    [Test]
    public void FindDonorSites_U12_GTAG_Donor_ReportedAsBothU2AndU12Candidates()
    {
        // GT-AG U12 introns (~70 % of U12 introns) start /GUAUCCUY; the terminal dinucleotide
        // does not decide the spliceosome (Dietrich et al. 1997), so the GU is scored both as a
        // U2 donor (S&S: AG|GUAUCC ⇒ t = 58+78+200+57+9+6+15 = 423 ⇒ 376/548) and as a U12 donor (8/8).
        var sites = FindDonorSites("CAGGUAUCCUUA", minScore: 0.0, includeNonCanonical: true).ToList();

        Assert.That(sites.Select(s => (s.Position, s.Type)),
            Is.EquivalentTo(new[] { (3, SpliceSiteType.Donor), (3, SpliceSiteType.U12Donor) }));
        Assert.Multiple(() =>
        {
            Assert.That(sites.Single(s => s.Type == SpliceSiteType.Donor).Score,
                Is.EqualTo(376 / PrimateRange).Within(1e-12));
            Assert.That(sites.Single(s => s.Type == SpliceSiteType.U12Donor).Score, Is.EqualTo(1.0).Within(1e-12));
        });
    }

    [Test]
    public void FindDonorSites_U12_PartialMatch_IsFractionOfEightPositions()
    {
        // /AUAACCAA: positions +4 (A≠U), +7 (A≠U), +8 (A∉Y) mismatch ⇒ 5/8.
        var site = FindDonorSites("CCCAUAACCAA", minScore: 0.0, includeNonCanonical: true)
            .Single(s => s.Type == SpliceSiteType.U12Donor);
        Assert.That(site.Score, Is.EqualTo(5.0 / 8.0).Within(1e-12));
    }

    #endregion

    #region C4: Non-ACGT characters

    [Test]
    public void FindDonorSites_NonAcgtOutsideScoredWindow_DoesNotAffectScore()
    {
        // S&S scores -2..+6 only; N at -3 is outside the subsequence ⇒ still 1.0.
        var site = FindDonorSites("NAGGUAAGU", minScore: 0.0).Single();
        Assert.Multiple(() =>
        {
            Assert.That(site.Position, Is.EqualTo(3));
            Assert.That(site.Score, Is.EqualTo(1.0).Within(1e-12));
        });
    }

    [Test]
    public void FindDonorSites_NonAcgtInsideScoredWindow_CandidateNotScored()
    {
        // N at +5 has no Table 1 percentage ⇒ the subsequence cannot be scored ⇒ no site.
        Assert.That(FindDonorSites("CAGGUAANU", minScore: 0.0), Is.Empty);
        Assert.That(FindDonorSites("NNNNNNNNN", minScore: 0.0), Is.Empty);
    }

    #endregion

    #region MaxEntScan mode — Yeo & Burge (2004) score5 over GU/GC candidates

    [Test]
    public void FindDonorSites_MaxEntScanMode_ScoresWindowWithScore5()
    {
        // maxentpy score5('CAGGTAAGT') = 10.858313101356437 bits; confidence 2^s/(1+2^s).
        var site = FindDonorSites("CAGGUAAGU", -100.0, false, DonorScoringModel.MaxEntScan).Single();
        Assert.Multiple(() =>
        {
            Assert.That(site.Position, Is.EqualTo(3));
            Assert.That(site.Type, Is.EqualTo(SpliceSiteType.Donor));
            Assert.That(site.Score, Is.EqualTo(10.858313101356437).Within(1e-12));
            Assert.That(site.Score, Is.EqualTo(ScoreDonorMaxEnt("CAGGUAAGU")));
            Assert.That(site.Confidence, Is.EqualTo(0.9994616209198077).Within(1e-12));
        });
    }

    [Test]
    public void FindDonorSites_MaxEntScanMode_GcDonorAndThresholdInBits()
    {
        // maxentpy score5('CAGGCAAGT') = 3.104161570996358 (GC donor, includeNonCanonical).
        Assert.Multiple(() =>
        {
            Assert.That(FindDonorSites("CAGGCAAGU", -100.0, true, DonorScoringModel.MaxEntScan).Single().Score,
                Is.EqualTo(3.104161570996358).Within(1e-12));
            Assert.That(FindDonorSites("CAGGCAAGU", -100.0, false, DonorScoringModel.MaxEntScan), Is.Empty);
            Assert.That(FindDonorSites("CAGGCAAGU", 3.2, true, DonorScoringModel.MaxEntScan), Is.Empty,
                "minScore is a threshold in bits");
            Assert.That(FindDonorSites("CAGAUAUCCUUA", -100.0, true, DonorScoringModel.MaxEntScan)
                .Any(s => s.Type == SpliceSiteType.U12Donor), Is.False, "score5 is a U2 GT model: no U12 candidates");
        });
    }

    [Test]
    public void FindDonorSites_MaxEntScanMode_WindowMustFitAndBeAcgu()
    {
        // The 9-nt window -3..+6 needs 3 exonic nt: "AGGUAAGU" (GU at 2) and a window with N are skipped.
        Assert.Multiple(() =>
        {
            Assert.That(FindDonorSites("AGGUAAGUAA", -100.0, false, DonorScoringModel.MaxEntScan), Is.Empty);
            Assert.That(FindDonorSites("CAGGUAANU", -100.0, false, DonorScoringModel.MaxEntScan), Is.Empty);
            Assert.That(FindDonorSites("AAGGUAAGUC", -100.0, false, DonorScoringModel.MaxEntScan).Single().Score,
                Is.EqualTo(10.998354285522845).Within(1e-12));
        });
    }

    #endregion

    #region ME1–ME9: MaxEntScan score5ss — Yeo & Burge (2004)

    // Expected values are the documented MaxEntScan score5 worked examples, taken from the
    // maxentpy score5 docstring (kepbod/maxentpy, MIT-licensed port retrieved 2026-06-25):
    //   score5('cagGTAAGT') -> 10.86
    //   score5('gagGTAAGT') -> 11.08
    //   score5('taaATAAGT') -> -0.12
    // These round-to-2dp targets were reproduced exactly this session (10.858313 / 11.078494 /
    // -0.116791) with the embedded table + the published factorisation. Asserted to 2 dp
    // (the documented precision) plus the full-precision value, so a wrong table or
    // factorisation fails the 10.86 cross-check rather than silently passing.

    // ME1 — canonical documented reference: score5('cagGTAAGT') == 10.86 (bits).
    [Test]
    public void ScoreDonorMaxEnt_CanonicalReferenceWindow_Returns10Point86()
    {
        // 9 nt: 3 exon + 6 intron; invariant GT at 0-based positions 3-4.
        const string window = "cagGTAAGT";

        double score = ScoreDonorMaxEnt(window);

        Assert.That(System.Math.Round(score, 2), Is.EqualTo(10.86),
            "MaxEntScan score5 of the canonical documented example must be 10.86 bits "
            + "(Yeo & Burge 2004 / maxentpy docstring)");
    }

    // ME2 — full-precision value behind the 10.86 reference (guards the factorisation/table).
    [Test]
    public void ScoreDonorMaxEnt_CanonicalReferenceWindow_MatchesFullPrecision()
    {
        const string window = "cagGTAAGT";

        double score = ScoreDonorMaxEnt(window);

        Assert.That(score, Is.EqualTo(10.858313).Within(1e-6),
            "Full-precision MaxEntScan score5 for the canonical window, reproduced from the "
            + "embedded table + published factorisation");
    }

    // ME3 — strong 5' site second documented reference: 11.08 bits.
    [Test]
    public void ScoreDonorMaxEnt_StrongerSiteReferenceWindow_Returns11Point08()
    {
        const string window = "gagGTAAGT";

        double score = ScoreDonorMaxEnt(window);

        Assert.Multiple(() =>
        {
            Assert.That(System.Math.Round(score, 2), Is.EqualTo(11.08),
                "Second documented MaxEntScan score5 example must be 11.08 bits");
            Assert.That(score, Is.EqualTo(11.078494).Within(1e-6),
                "Full-precision value behind the 11.08 reference");
        });
    }

    // ME4 — weak 5' site third documented reference: -0.12 bits (non-GT, negative score).
    [Test]
    public void ScoreDonorMaxEnt_WeakSiteReferenceWindow_ReturnsMinus0Point12()
    {
        const string window = "taaATAAGT";

        double score = ScoreDonorMaxEnt(window);

        Assert.Multiple(() =>
        {
            Assert.That(System.Math.Round(score, 2), Is.EqualTo(-0.12),
                "Third documented MaxEntScan score5 example must be -0.12 bits");
            Assert.That(score, Is.EqualTo(-0.116791).Within(1e-6),
                "Full-precision value behind the -0.12 reference");
        });
    }

    // ME5 — a strong site ranks above a weak site (ordering of the documented examples).
    [Test]
    public void ScoreDonorMaxEnt_StrongSite_RanksAboveWeakSite()
    {
        double strong = ScoreDonorMaxEnt("cagGTAAGT"); // 10.86
        double weak = ScoreDonorMaxEnt("taaATAAGT");   // -0.12

        Assert.That(strong, Is.GreaterThan(weak),
            "A strong 5' splice site must score higher than a weak one (10.86 > -0.12)");
    }

    // ME6 — DNA (T) and RNA (U) windows are equivalent (T==U in the rest key and GT model).
    [Test]
    public void ScoreDonorMaxEnt_DnaAndRnaWindows_ProduceIdenticalScores()
    {
        const string dna = "cagGTAAGT";
        string rna = dna.Replace('t', 'u').Replace('T', 'U');

        double dnaScore = ScoreDonorMaxEnt(dna);
        double rnaScore = ScoreDonorMaxEnt(rna);

        Assert.That(rnaScore, Is.EqualTo(dnaScore).Within(1e-12),
            "U-form (RNA) window must score identically to the T-form (DNA) window");
    }

    // ME7 — case-insensitive: upper-case window scores identically.
    [Test]
    public void ScoreDonorMaxEnt_UpperCaseWindow_ScoresIdenticallyToLowerCase()
    {
        const string lower = "cagGTAAGT";
        string upper = lower.ToUpperInvariant();

        Assert.That(ScoreDonorMaxEnt(upper), Is.EqualTo(ScoreDonorMaxEnt(lower)).Within(1e-12),
            "Scoring must be case-insensitive");
    }

    // ME8 — null window throws ArgumentNullException.
    [Test]
    public void ScoreDonorMaxEnt_NullWindow_Throws()
    {
        Assert.Throws<System.ArgumentNullException>(() => ScoreDonorMaxEnt(null!),
            "A null window is invalid input");
    }

    // ME9 — wrong length and invalid alphabet throw ArgumentException.
    [Test]
    public void ScoreDonorMaxEnt_InvalidWindow_Throws()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<System.ArgumentException>(() => ScoreDonorMaxEnt("cagGTAAG"),
                "An 8-nt window (wrong length) is invalid");
            Assert.Throws<System.ArgumentException>(() => ScoreDonorMaxEnt("cagGTAAGTT"),
                "A 10-nt window (wrong length) is invalid");
            Assert.Throws<System.ArgumentException>(() => ScoreDonorMaxEnt("cagGTAANT"),
                "A non-A/C/G/T(/U) character (N) is invalid");
        });
    }

    // ME10 — reference sweep (SPLICE-MAXENT5-001, B11 F19). Expected values are maxentpy score5
    // (kepbod/maxentpy maxent.py + score5_matrix.txt) full-precision outputs; the CPAN port of
    // score5.pl (Bio-RNA-SpliceSites-Scoring-MaxEntScan 0.05, me2x5 + splice5sequences) prints the
    // same 2-decimal values. The C# score equals maxentpy bit-for-bit on all 4^9 = 262144 ACGT 9-mers.
    // Includes non-GT dinucleotides (GC, AT, CC, ...), which MaxEntScan scores rather than rejects.
    [TestCase("CAGGTAAGT", 10.858313101356437, 10.86)]
    [TestCase("GGGGTAAGT", 9.652449523571205, 9.65)]
    [TestCase("CACGTAAGT", 10.155128271120882, 10.16)]
    [TestCase("AAAGTAAAA", 1.1774343078815603, 1.18)]
    [TestCase("CAGGCAAGT", 3.104161570996358, 3.10)]
    [TestCase("CAGATAAGT", 2.676285865160975, 2.68)]
    [TestCase("AACATCTGG", -12.126838926748352, -12.13)]
    [TestCase("GGACCATGG", -15.483798810300666, -15.48)]
    [TestCase("TGCTTAATG", -9.751934925348701, -9.75)]
    [TestCase("CTTCCGGAG", -31.261128275661, -31.26)]
    [TestCase("ACAGTGTTT", -8.544635130848441, -8.54)]
    [TestCase("TAGTACTCT", -20.020155712973068, -20.02)]
    [TestCase("CCTATTCCC", -30.870167664881997, -30.87)]
    [TestCase("TTCCGGCGT", -12.089519108744033, -12.09)]
    [TestCase("AAAAAAAAA", -15.188009382392348, -15.19)]
    public void ScoreDonorMaxEnt_ReferenceSweep_MatchesMaxEntScan(string window, double maxentpy, double score5Pl)
    {
        double score = ScoreDonorMaxEnt(window);
        Assert.Multiple(() =>
        {
            Assert.That(score, Is.EqualTo(maxentpy).Within(1e-12), "maxentpy score5 (full precision)");
            Assert.That(System.Math.Round(score, 2, System.MidpointRounding.AwayFromZero), Is.EqualTo(score5Pl),
                "score5.pl prints %.2f");
        });
    }

    #endregion
}
