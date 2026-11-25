using static Seqeron.Genomics.Annotation.SpliceSitePredictor;

namespace Seqeron.Genomics.Tests.Unit.Annotation;

/// <summary>
/// Canonical test file for SPLICE-ACCEPTOR-001: Acceptor Site Detection.
/// Tests <see cref="SpliceSitePredictor.FindAcceptorSites(string, double, bool)"/> (Shapiro &amp;
/// Senapathy 1987 3' ss score, MaxEntScan mode, U12 branch-site acceptors) and the branch-point
/// detectors.
///
/// Evidence sources:
///   - Shapiro &amp; Senapathy (1987) Nucleic Acids Res 15(17):7155–7174
///   - Burge, Tuschl &amp; Sharp (1999) The RNA World, 2nd ed., pp. 525–560
///   - Yeo &amp; Burge (2004) J Comput Biol 11(2–3):377–394
///   - Patel &amp; Steitz (2003) Nat Rev Mol Cell Biol 4(12):960–970
///   - Gao, Masuda, Matsuura &amp; Ohno (2008) Nucleic Acids Res 36(7):2257–2267
///     (branch-point consensus yUnAy; DOI 10.1093/nar/gkn073)
///   - Mercer et al. (2015) Genome Res 25(2):290–303 (branch-point distribution)
/// </summary>
[TestFixture]
public class SpliceSitePredictor_AcceptorSite_Tests
{
    #region M1: Canonical AG — Shapiro & Senapathy (1987) 3' ss score

    [Test]
    public void FindAcceptorSites_CanonicalAG_WithStrongPPT_FindsAcceptorSite()
    {
        // U×16 + CAG|GG: AG at index 17–18. S&S window -14..+1 = indices 5..19: ten U at -14..-5
        // (best 8 of 10 = the 8 highest primate U percentages, = h1, since row -8 has C 43 > U 42),
        // -4 = U (unscored), C A G | G = the row maxima 74+100+100+49 = h2 → score exactly 1.0.
        string sequence = "UUUUUUUUUUUUUUUUCAGGG";

        var sites = FindAcceptorSites(sequence, minScore: 0.1).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(sites, Has.Count.EqualTo(1), "Exactly one AG");
            Assert.That(sites[0].Type, Is.EqualTo(SpliceSiteType.Acceptor));
            Assert.That(sites[0].Position, Is.EqualTo(18), "Position = A_index + 1 = 18 (G of AG)");
            Assert.That(sites[0].Score, Is.EqualTo(1.0), "t1 = h1 and t2 = h2 → S&S score 100/100");
            Assert.That(sites[0].Confidence, Is.EqualTo(1.0));
            Assert.That(sites[0].Motif, Is.EqualTo("UUUUUUUUUUUUUUUUCAGGG"),
                "Motif = -20..+3 context, truncated at the sequence start (16 + 5 nt)");
        });
    }

    // Shapiro & Senapathy (1987) Table 5: published 3' ss scores of real introns (0-100), each
    // scored with its GENBANK category subtable of Table 1. Sequence = intron -14..-1 / exon +1.
    // 40 of the 41 published values are reproduced at 1 dp (rabbit MHC intron 3, 86.9 vs 86.8,
    // differs by table rounding and is listed separately below).
    [TestCase("CTCCTCGCTCCCAGG", ShapiroSenapathyCategory.Primate, 91.4)]
    [TestCase("CGGGGGCGGGCCAGG", ShapiroSenapathyCategory.Primate, 61.7)]
    [TestCase("CTCTTTCCCGTCAGA", ShapiroSenapathyCategory.Primate, 91.0)]
    [TestCase("CCCCCTTTTCCCAGA", ShapiroSenapathyCategory.Primate, 87.6)]
    [TestCase("TTTTCTTCCCACAGA", ShapiroSenapathyCategory.Primate, 93.8)]
    [TestCase("TTATTCTACTCCAGG", ShapiroSenapathyCategory.Primate, 95.8)]
    [TestCase("TATTCCCTTCCCAGG", ShapiroSenapathyCategory.Rodent, 94.3)]
    [TestCase("TTCGTCCTCCCCAGG", ShapiroSenapathyCategory.Rodent, 92.4)]
    [TestCase("GTCTCAATACCTAGG", ShapiroSenapathyCategory.Rodent, 72.5)]
    [TestCase("TCTTGTCCCCACAGG", ShapiroSenapathyCategory.Rodent, 95.9)]
    [TestCase("CTGTATCCCTGCAGC", ShapiroSenapathyCategory.Rodent, 87.6)]
    [TestCase("ACCTTTCTCTGCAGG", ShapiroSenapathyCategory.Rodent, 95.6)]
    [TestCase("CCTGTCGCCCGCAGG", ShapiroSenapathyCategory.OtherMammal, 87.3)]
    [TestCase("CGGGGCGGGGCCAGG", ShapiroSenapathyCategory.OtherMammal, 59.6)]
    [TestCase("CCTTCCCTTCTCAGA", ShapiroSenapathyCategory.OtherMammal, 86.2)]
    [TestCase("TTTTCTTCTCACAGA", ShapiroSenapathyCategory.OtherMammal, 92.7)]
    [TestCase("TTTGTTCACCCCAGG", ShapiroSenapathyCategory.OtherMammal, 96.0)]
    [TestCase("GGCTCTGCCTCCAGG", ShapiroSenapathyCategory.OtherVertebrate, 84.9)]
    [TestCase("TCCTGCTCTTGCAGC", ShapiroSenapathyCategory.OtherVertebrate, 86.7)]
    [TestCase("TCCCTTCTCCTCAGG", ShapiroSenapathyCategory.OtherVertebrate, 88.0)]
    [TestCase("TGTACCTCCTGTAGT", ShapiroSenapathyCategory.OtherVertebrate, 75.8)]
    [TestCase("CTGTTTTGGTGCAGG", ShapiroSenapathyCategory.OtherVertebrate, 91.0)]
    [TestCase("GATTTCCTTGATAGG", ShapiroSenapathyCategory.OtherVertebrate, 82.5)]
    [TestCase("TATTGTAATTTCAGG", ShapiroSenapathyCategory.Invertebrate, 90.0)]
    [TestCase("GGTTTTCAATGTAGT", ShapiroSenapathyCategory.Invertebrate, 70.6)]
    [TestCase("GTTTCCTATTTCAGA", ShapiroSenapathyCategory.Invertebrate, 94.5)]
    [TestCase("TATAAATTTTGCAGA", ShapiroSenapathyCategory.Invertebrate, 92.9)]
    [TestCase("TTCTTATTTTTCAGC", ShapiroSenapathyCategory.Invertebrate, 93.7)]
    [TestCase("TGTGAATCTTTCAGG", ShapiroSenapathyCategory.Invertebrate, 91.3)]
    [TestCase("TTCATGGTTTTCAGG", ShapiroSenapathyCategory.Invertebrate, 90.8)]
    [TestCase("CTTATTAATTTCAGG", ShapiroSenapathyCategory.Invertebrate, 91.5)]
    [TestCase("CTGGACCCGTGCAGC", ShapiroSenapathyCategory.Plant, 61.9)]
    [TestCase("CTTATCTGTCTCAGG", ShapiroSenapathyCategory.Plant, 82.8)]
    [TestCase("CTTGATTTTGCCAGT", ShapiroSenapathyCategory.Plant, 79.1)]
    [TestCase("CTTTCTCTGTTTAGG", ShapiroSenapathyCategory.Plant, 84.8)]
    [TestCase("TGCGCATGGTTAAGG", ShapiroSenapathyCategory.Plant, 66.5)]
    [TestCase("GTCCCATTTTTCAGC", ShapiroSenapathyCategory.Plant, 78.3)]
    [TestCase("CCTTTATGGTCTAGG", ShapiroSenapathyCategory.Plant, 79.9)]
    [TestCase("CCCTCGTGATCCAGG", ShapiroSenapathyCategory.Plant, 74.7)]
    [TestCase("GCAATTCTGCACAGG", ShapiroSenapathyCategory.Plant, 70.9)]
    public void FindAcceptorSites_ShapiroSenapathy_ReproducesPublishedTable5Scores(
        string site, ShapiroSenapathyCategory category, double published)
    {
        var sites = FindAcceptorSites(site, 0.0, false, AcceptorScoringModel.ShapiroSenapathy, category)
            .Where(s => s.Position == 13).ToList();

        Assert.That(sites, Has.Count.EqualTo(1));
        Assert.That(Math.Round(sites[0].Score * 100, 1, MidpointRounding.AwayFromZero), Is.EqualTo(published),
            $"S&S (1987) Table 5 published 3' ss score for {site} ({category})");
    }

    [Test]
    public void FindAcceptorSites_ShapiroSenapathy_Table5RoundingOutlier_ExactFraction()
    {
        // Rabbit MHC intron 3 (other mammal), published 86.9. Table 1 (rounded percentages) gives
        // t1 = 381, l1 = 27, h1 = 442; t2 = 79+100+100+13 = 292, l2 = 6, h2 = 330:
        // ((381-27)/415 + (292-6)/324)/2 = 0.86786… (86.8): a 0.1 table-rounding difference.
        var site = FindAcceptorSites("AATCTTTCCCTCAGA", 0.0, false,
            AcceptorScoringModel.ShapiroSenapathy, ShapiroSenapathyCategory.OtherMammal).Single();

        Assert.That(site.Score, Is.EqualTo((354.0 / 415 + 286.0 / 324) / 2).Within(1e-15));
    }

    [Test]
    public void FindAcceptorSites_ShapiroSenapathy_AgInPyrimidineRegion_Eliminated()
    {
        // S&S (1987) Methods: "any sequence containing an AG in those 10 locations [is]
        // eliminated from consideration". AG at -9/-8 (indices 8-9) → no site; AA instead → site.
        Assert.Multiple(() =>
        {
            Assert.That(FindAcceptorSites("UUUUUUUUAGUUUUUUCAGGG", 0.0).ToList(), Is.Empty);
            Assert.That(FindAcceptorSites("UUUUUUUUAAUUUUUUCAGGG", 0.0).Single().Score,
                Is.EqualTo(0.978125).Within(1e-15));
        });
    }

    [Test]
    public void FindAcceptorSites_ShapiroSenapathy_TwoPurinesScoreAsWellAsNone()
    {
        // S&S (1987) Methods: the best 8 of the 10 pyrimidine locations are used, "thus a pyrimidine
        // sequence containing two purines can score as well as one containing no purines". Purines
        // at -8/-7 (primate U 42, 42 — the two lowest U rows) leave t1 unchanged → 1.0.
        var site = FindAcceptorSites("UUUUUUUUUUUAAUUUCAGGG", 0.0).Single();

        Assert.That(site.Score, Is.EqualTo(1.0));
    }

    [Test]
    public void FindAcceptorSites_MaxEntScanModel_DelegatesToScore3()
    {
        // MaxEntScan score3 window -20..+3 = the whole 23-nt input (reference 2.89).
        var site = FindAcceptorSites("TTCCAAACGAACTTTTGTAGGGA", -100.0, false,
            AcceptorScoringModel.MaxEntScan).Single();

        Assert.Multiple(() =>
        {
            Assert.That(site.Position, Is.EqualTo(19));
            Assert.That(site.Score, Is.EqualTo(ScoreAcceptorMaxEnt("TTCCAAACGAACTTTTGTAGGGA")));
            Assert.That(site.Score, Is.EqualTo(2.8867730651152104));
            Assert.That(site.Confidence, Is.EqualTo(1.0 / (1.0 + Math.Pow(2.0, -site.Score))));
        });
    }

    [Test]
    public void FindAcceptorSites_MaxEntScanModel_PartialWindowSkipped_MinScoreInBits()
    {
        // 22 nt: no full 23-nt window → no site; threshold in bits filters the 2.89 site.
        Assert.Multiple(() =>
        {
            Assert.That(FindAcceptorSites("TCCAAACGAACTTTTGTAGGGA", -100.0, false,
                AcceptorScoringModel.MaxEntScan).ToList(), Is.Empty);
            Assert.That(FindAcceptorSites("TTCCAAACGAACTTTTGTAGGGA", 3.0, false,
                AcceptorScoringModel.MaxEntScan).ToList(), Is.Empty);
        });
    }

    [Test]
    public void FindAcceptorSites_UnknownModelOrCategory_Throws()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentOutOfRangeException>(() =>
                FindAcceptorSites("UUUUUUUUUUUUUUUUCAGGG", 0.5, false, (AcceptorScoringModel)7));
            Assert.Throws<ArgumentOutOfRangeException>(() =>
                FindAcceptorSites("UUUUUUUUUUUUUUUUCAGGG", 0.5, false,
                    AcceptorScoringModel.ShapiroSenapathy, (ShapiroSenapathyCategory)42));
        });
    }

    #endregion

    #region M2: No AG Returns Empty — Trivially Correct

    [Test]
    public void FindAcceptorSites_NoAGDinucleotide_ReturnsEmpty()
    {
        // All pyrimidines — no AG dinucleotide anywhere
        string sequence = "UUUUUUUUUUUUUUUUUUUUU";

        var sites = FindAcceptorSites(sequence, minScore: 0.1).ToList();

        Assert.That(sites, Is.Empty,
            "Sequence without AG dinucleotide should produce no acceptor sites");
    }

    #endregion

    #region M3: Empty/Null Input Returns Empty — Guard Behavior

    [Test]
    public void FindAcceptorSites_EmptyInput_ReturnsEmpty()
    {
        var empty = FindAcceptorSites("", minScore: 0.1).ToList();
        var nullResult = FindAcceptorSites(null!, minScore: 0.1).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(empty, Is.Empty,
                "Empty string should return empty (guard condition)");
            Assert.That(nullResult, Is.Empty,
                "Null input should return empty (guard condition)");
        });
    }

    #endregion

    #region M4: Window Shorter Than -14..+1 Returns Empty — full-window scoring

    [Test]
    public void FindAcceptorSites_ShortSequence_ReturnsEmpty()
    {
        // 14 nt: the AG has only 11 upstream nt, the S&S window needs -14..+1 (15 nt).
        // 15 nt: exactly the window → scored (1.0).
        Assert.Multiple(() =>
        {
            Assert.That(FindAcceptorSites("UUUUUUUUUUCAGG", minScore: 0.0).ToList(), Is.Empty,
                "Window -14..+1 runs off the sequence → no site");
            Assert.That(FindAcceptorSites("UUUUUUUUUUUCAGG", minScore: 0.0).Single().Position, Is.EqualTo(13),
                "Exactly the 15-nt window → scored");
            Assert.That(FindAcceptorSites("UUUUUUUUUUUUUUUUCAG", minScore: 0.0).ToList(), Is.Empty,
                "No exon nucleotide +1 → no site");
        });
    }

    #endregion

    #region M5: Strong PPT > purine region — Shapiro & Senapathy (1987)

    [Test]
    public void FindAcceptorSites_StrongPPT_ScoresHigherThanWeakContext()
    {
        // Pyrimidine region all U (1.0) vs all A: t1 = best 8 primate A percentages
        // 10+10+9+9+9+7+7+7 = 68; l1 = 8 lowest row minima = 57; h1 = 377; CAGG part maximal →
        // ((68-57)/(377-57) + 1)/2 = 0.5171875.
        var strong = FindAcceptorSites("UUUUUUUUUUUUUUUUCAGGG", minScore: 0.0).Single();
        var weak = FindAcceptorSites("AAAAAAAAAAAAAAAACAGGG", minScore: 0.0).Single();

        Assert.Multiple(() =>
        {
            Assert.That(strong.Score, Is.EqualTo(1.0));
            Assert.That(weak.Score, Is.EqualTo(0.5171875).Within(1e-15));
            Assert.That(strong.Score, Is.GreaterThan(weak.Score));
        });
    }

    #endregion

    #region M6: Score Range [0, 1] — Normalization (INV-1)

    [Test]
    public void FindAcceptorSites_AllScores_InZeroOneRange()
    {
        string sequence = "UUUUUUUUUUUUUUUUCAGGG";

        var sites = FindAcceptorSites(sequence, minScore: 0.0).ToList();

        Assert.That(sites, Is.Not.Empty, "Should find at least one site with minScore=0");

        foreach (var site in sites)
        {
            Assert.That(site.Score, Is.InRange(0.0, 1.0),
                $"Score for site at position {site.Position} must be in [0, 1] — normalization invariant");
        }
    }

    #endregion

    #region M7: Confidence Range [0, 1] — CalculateConfidence (INV-2)

    [Test]
    public void FindAcceptorSites_AllConfidences_InZeroOneRange()
    {
        string sequence = "UUUUUUUUUUUUUUUUCAGGG";

        var sites = FindAcceptorSites(sequence, minScore: 0.0).ToList();

        Assert.That(sites, Is.Not.Empty, "Should find at least one site with minScore=0");

        foreach (var site in sites)
        {
            Assert.That(site.Confidence, Is.InRange(0.0, 1.0),
                $"Confidence for site at position {site.Position} must be in [0, 1]");
        }
    }

    #endregion

    #region M8: DNA T Equivalence — Implementation T→U Conversion

    [Test]
    public void FindAcceptorSites_DnaInput_ProducesSameResults()
    {
        // RNA input
        string rna = "UUUUUUUUUUUUUUUUCAGGG";
        // DNA equivalent (T instead of U)
        string dna = "TTTTTTTTTTTTTTTTTCAGGG";

        var rnaSites = FindAcceptorSites(rna, minScore: 0.1).ToList();
        var dnaSites = FindAcceptorSites(dna, minScore: 0.1).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(dnaSites.Count, Is.EqualTo(rnaSites.Count),
                "DNA (T) and RNA (U) inputs should produce the same number of acceptor sites");

            if (rnaSites.Count > 0 && dnaSites.Count > 0)
            {
                Assert.That(dnaSites[0].Score, Is.EqualTo(rnaSites[0].Score).Within(1e-10),
                    "Scores should be identical for DNA vs RNA input");
            }
        });
    }

    #endregion

    #region M9: Case Insensitivity — Implementation ToUpperInvariant

    [Test]
    public void FindAcceptorSites_LowercaseInput_ProducesSameResults()
    {
        string upper = "UUUUUUUUUUUUUUUUCAGGG";
        string lower = "uuuuuuuuuuuuuuuucaggg";

        var upperSites = FindAcceptorSites(upper, minScore: 0.1).ToList();
        var lowerSites = FindAcceptorSites(lower, minScore: 0.1).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(lowerSites.Count, Is.EqualTo(upperSites.Count),
                "Lowercase input should produce same number of sites as uppercase — case insensitivity");
            if (upperSites.Count > 0)
            {
                Assert.That(lowerSites[0].Score, Is.EqualTo(upperSites[0].Score).Within(1e-10),
                    "Scores must be identical for lowercase vs uppercase — ToUpperInvariant normalization");
                Assert.That(lowerSites[0].Position, Is.EqualTo(upperSites[0].Position),
                    "Positions must match for lowercase vs uppercase");
            }
        });
    }

    #endregion

    #region M10: Multiple AG Sites — Scanning Algorithm

    [Test]
    public void FindAcceptorSites_MultipleAGSites_FindsAll()
    {
        // Two strong acceptors 18 nt apart: AG at 17 and 35; both windows are AG-free → both 1.0.
        string twoSites = "UUUUUUUUUUUUUUUUCAGGUUUUUUUUUUUUUUCAGG";
        // AG at 21 has the AG at 17-18 inside its -14..-5 region → eliminated (S&S 1987).
        string shadowed = "UUUUUUUUUUUUUUUUCAGCUAGGG";

        var sites = FindAcceptorSites(twoSites, minScore: 0.0).ToList();
        var shadow = FindAcceptorSites(shadowed, minScore: 0.0).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(sites.Select(s => s.Position), Is.EqualTo(new[] { 18, 36 }));
            Assert.That(sites.Select(s => s.Score), Is.EqualTo(new[] { 1.0, 1.0 }));
            Assert.That(shadow.Select(s => s.Position), Is.EqualTo(new[] { 18 }),
                "Second AG eliminated: AG inside its pyrimidine region");
            Assert.That(shadow[0].Score, Is.EqualTo(0.9423076923076923).Within(1e-15),
                "+1 = C (13) instead of G (49): t2 = 287, l2 = 11, h2 = 323 → (1 + 276/312)/2");
        });
    }

    #endregion

    #region S1: U12 acceptor — U12 branch site UUCCUUAAC (Brock et al. 2008; Levine & Durbin 2001)

    [Test]
    public void FindAcceptorSites_U12AcceptorYCCAC_DetectedWhenNonCanonicalEnabled()
    {
        // Consensus U12 branch site UUCCUUAAC, branch A (second A) at index 8; AT-AC terminus with
        // the terminal C at index 16 → branch-to-acceptor distance 8 (window 8-21) → 9/9 = 1.0.
        string sequence = "GUUCCUUAACGGGGGACGGG";

        var u12Sites = FindAcceptorSites(sequence, minScore: 0.1, includeNonCanonical: true)
            .Where(s => s.Type == SpliceSiteType.U12Acceptor).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(u12Sites, Has.Count.EqualTo(1));
            Assert.That(u12Sites[0].Position, Is.EqualTo(16), "Position = index of the terminal C of AC");
            Assert.That(u12Sites[0].Score, Is.EqualTo(1.0));
        });
    }

    [TestCase(4, 0.4444444444444444)]  // distance 7 → consensus out of window
    [TestCase(5, 1.0)]                 // distance 8 (window start)
    [TestCase(18, 1.0)]                // distance 21 (window end)
    [TestCase(19, 0.4444444444444444)] // distance 22 → out of window
    public void FindAcceptorSites_U12_BranchSiteDistanceWindow_8To21(int spacer, double expected)
    {
        string sequence = "G" + "UUCCUUAAC" + new string('G', spacer) + "AC" + "GGG";
        int last = 1 + 9 + spacer + 1;

        var site = FindAcceptorSites(sequence, 0.0, true)
            .Single(s => s.Type == SpliceSiteType.U12Acceptor && s.Position == last);

        Assert.That(site.Score, Is.EqualTo(expected).Within(1e-15));
    }

    [Test]
    public void FindAcceptorSites_U12_GtAgTerminus_ReportedAsBothU2AndU12Candidates()
    {
        // GT-AG U12 introns end in AG (Dietrich et al. 1997): the AG is scored by S&S as a U2
        // acceptor and by the branch site as a U12 candidate.
        var sites = FindAcceptorSites("GUUCCUUAACGGGGGGGGGGAGGGG", 0.0, true).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(sites.Select(s => (s.Position, s.Type)), Is.EqualTo(new[]
            {
                (21, SpliceSiteType.Acceptor), (21, SpliceSiteType.U12Acceptor)
            }));
            Assert.That(sites[1].Score, Is.EqualTo(1.0));
        });
    }

    [Test]
    public void FindAcceptorSites_U12_PartialBranchSite_ExactFraction()
    {
        // UUCCAUAAC: one mismatch (-3 relative to branch A) → 8/9.
        var site = FindAcceptorSites("GUUCCAUAACGGGGGGGGGGACGGG", 0.0, true)
            .Single(s => s.Type == SpliceSiteType.U12Acceptor && s.Position == 21);

        Assert.That(site.Score, Is.EqualTo(8.0 / 9).Within(1e-15));
    }

    #endregion

    #region S2: U12 AC Excluded When Canonical Only — Default Parameter

    [Test]
    public void FindAcceptorSites_U12AcceptorYCCAC_ExcludedByDefault()
    {
        // YCCAC dinucleotide but no AG — canonical-only mode should find no U12 sites
        string sequence = "UUUUUUUUUUUUUUCCACGGG";

        var sites = FindAcceptorSites(sequence, minScore: 0.1).ToList();

        Assert.That(sites.Where(s => s.Type == SpliceSiteType.U12Acceptor).ToList(), Is.Empty,
            "U12 YCCAC acceptor should not be detected with includeNonCanonical=false (default)");
    }

    #endregion

    #region S3: Position After AG — INV-5

    [Test]
    public void FindAcceptorSites_ReturnsPositionAfterAG()
    {
        // U×15 + A(15) + G(16) + G(17) + G(18) + G(19) + G(20) = 21 nt
        // AG at indices 15 (A) and 16 (G)
        string sequence = "UUUUUUUUUUUUUUUAGGGGG";

        var sites = FindAcceptorSites(sequence, minScore: 0.0).ToList();

        Assert.That(sites, Has.Count.EqualTo(1),
            "Exactly one AG in scannable range");

        // Position = i + 1 where i is the index of A in AG
        // So Position = 15 + 1 = 16 (index of G in AG — last intronic nucleotide)
        Assert.That(sites[0].Position, Is.EqualTo(16),
            "Position = A_index + 1 = 16 (G of AG, last intronic nucleotide) — INV-5");
    }

    #endregion

    #region S4: Motif Non-Empty — INV-6

    [Test]
    public void FindAcceptorSites_MotifIsNonEmpty()
    {
        string sequence = "UUUUUUUUUUUUUUUUCAGGG";

        var sites = FindAcceptorSites(sequence, minScore: 0.1).ToList();

        Assert.That(sites, Is.Not.Empty, "Should find at least one site");

        foreach (var site in sites)
        {
            Assert.That(site.Motif, Is.Not.Null.And.Not.Empty,
                $"Motif at position {site.Position} must be non-empty — context around AG");
        }
    }

    #endregion

    #region S5: Threshold Filtering — INV-8

    [Test]
    public void FindAcceptorSites_HigherThreshold_ProducesSubset()
    {
        string sequence = "UUUUUUUUUUUUUUUUCAGGG";

        var lowThreshold = FindAcceptorSites(sequence, minScore: 0.1).ToList();
        var highThreshold = FindAcceptorSites(sequence, minScore: 0.8).ToList();

        Assert.That(highThreshold.Count, Is.LessThanOrEqualTo(lowThreshold.Count),
            "Higher minScore should produce equal or fewer results than lower minScore — INV-8");
    }

    #endregion

    #region C1: AG without a full -14..+1 window not detected

    [Test]
    public void FindAcceptorSites_AGBeforeScanStart_NotDetected()
    {
        // AG at index 5: the S&S window would start at -7 → not scored.
        string sequence = "UUUUUAGUUUUUUUUUUUUUU";

        var sites = FindAcceptorSites(sequence, minScore: 0.0).ToList();

        Assert.That(sites, Is.Empty,
            "AG with fewer than 12 upstream nt has no -14..+1 window → not scored");
    }

    #endregion

    #region C2: Pyrimidine region best 8 of 10 — Shapiro & Senapathy (1987)

    [Test]
    public void FindAcceptorSites_PPTContribution_VerifiedByHelper()
    {
        // Independent recomputation of the S&S formula from the primate Table 1 rows for
        // UUUUUUUUAAUUUUUUCAGGG (AG at 17; -14..-5 = UUUAAUUUUU): t1 = best 8 of
        // {45,45,51,7,10,44,42,42,48,46} = 363; l1 = 8 lowest of the row minima
        // {9,9,7,7,7,10,7,8,6,6} = 57; h1 = 8 highest of the row maxima
        // {45,45,51,51,47,44,43,42,48,46} = 377; CAGG part maximal (t2 = h2).
        int[] minima = { 9, 9, 7, 7, 7, 10, 7, 8, 6, 6 };
        int[] maxima = { 45, 45, 51, 51, 47, 44, 43, 42, 48, 46 };
        int[] observed = { 45, 45, 51, 7, 10, 44, 42, 42, 48, 46 };
        int t1 = observed.OrderByDescending(v => v).Take(8).Sum();
        int l1 = minima.OrderBy(v => v).Take(8).Sum();
        int h1 = maxima.OrderByDescending(v => v).Take(8).Sum();
        double expected = ((double)(t1 - l1) / (h1 - l1) + 1.0) / 2;

        var site = FindAcceptorSites("UUUUUUUUAAUUUUUUCAGGG", minScore: 0.0).Single();

        Assert.Multiple(() =>
        {
            Assert.That((t1, l1, h1), Is.EqualTo((363, 57, 377)));
            Assert.That(site.Score, Is.EqualTo(expected).Within(1e-15));
            Assert.That(site.Score, Is.EqualTo(0.978125).Within(1e-15));
        });
    }

    #endregion

    #region C3: Non-ACGT characters — full-window scoring

    [Test]
    public void FindAcceptorSites_NonACGTCharacters_ToleratedAndAGStillFound()
    {
        // A percentage table has no N row: an N inside -14..+1 → candidate not scored; N outside
        // the window (indices 0-4 for the AG at 17) does not affect the score. No throw.
        Assert.Multiple(() =>
        {
            Assert.That(FindAcceptorSites("NNNNNNNNNNNNNNNNCAGGG", 0.0).ToList(), Is.Empty);
            Assert.That(FindAcceptorSites("UUUUUNUUUUUUUUUUCAGGG", 0.0).ToList(), Is.Empty);
            var site = FindAcceptorSites("NNNNNUUUUUUUUUUUCAGGG", 0.0).Single();
            Assert.That(site.Position, Is.EqualTo(18));
            Assert.That(site.Score, Is.EqualTo(1.0));
        });
    }

    #endregion

    #region BP: Branch-Point Detection — Gao et al. (2008) yUnAy consensus

    // BP1 — Canonical branch point at a known position and distance.
    // Consensus yUnAy (positions -3..+1); branch A at position 0.
    // Sequence: GG(pad) + C U U A C (branchA index 5) + U×22 (PPT) + A G + GGG.
    //   index 2=C(-3 y), 3=U(-2 U), 4=U(-1 n), 5=A(0 branch), 6=C(+1 y) → perfect yUnAy.
    //   AG: A at index 29, G at index 30 → FindAcceptorSites Position = 30.
    //   distance = agEnd(30) − branchA(5) = 25 nt (within 18..40 window).
    // Conservation-weighted score = (0.790+0.746+0.923+0.751)/3.210 = 1.0 (perfect).
    [Test]
    public void FindAcceptorBranchPoint_CanonicalYUnAy_DetectsPositionMotifAndScore()
    {
        string sequence = "UU" + "CUUAC" + new string('U', 22) + "AG" + "GGG";

        var bp = FindAcceptorBranchPoint(sequence, acceptorAgPosition: 30);

        Assert.Multiple(() =>
        {
            Assert.That(bp.Found, Is.True,
                "A canonical yUnAy branch point 25 nt upstream of the AG must be detected — Gao et al. (2008)");
            Assert.That(bp.BranchPointPosition, Is.EqualTo(5),
                "Branch-point adenosine is at index 5 (motif position 0)");
            Assert.That(bp.DistanceFromAg, Is.EqualTo(25),
                "Branch point is 25 nt upstream of the AG (within the 18–40 nt window) — Gao et al. (2008)");
            Assert.That(bp.Motif, Is.EqualTo("CUUAC"),
                "Reported motif is the 5-nt yUnAy window (positions -3..+1)");
            Assert.That(bp.Score, Is.EqualTo(1.0).Within(1e-10),
                "Perfect yUnAy (y@-3, U@-2, A@0, y@+1) → conservation-weighted score = 1.0");
            Assert.That(bp.PolypyrimidineTractFraction, Is.EqualTo(1.0).Within(1e-10),
                "Tract between branch point and AG is all U → pyrimidine fraction = 1.0");
        });
    }

    // BP2 — Integration with FindAcceptorSites: the AG position fed to the
    // branch-point detector is the one the acceptor scanner reports.
    [Test]
    public void FindAcceptorBranchPoint_UsingAcceptorSitePosition_LocatesBranchPoint()
    {
        string sequence = "UU" + "CUUAC" + new string('U', 22) + "AG" + "GGG";

        var acceptor = FindAcceptorSites(sequence, minScore: 0.0).Single();
        var bp = FindAcceptorBranchPoint(sequence, acceptor.Position);

        Assert.Multiple(() =>
        {
            Assert.That(acceptor.Position, Is.EqualTo(30),
                "Acceptor AG terminal G is at index 30");
            Assert.That(bp.Found, Is.True,
                "Branch point found upstream of the detected acceptor");
            Assert.That(bp.BranchPointPosition, Is.EqualTo(5),
                "Branch-point adenosine at index 5");
            Assert.That(bp.DistanceFromAg, Is.EqualTo(25),
                "25 nt between branch point and the acceptor AG");
        });
    }

    // BP3 — No branch adenosine in window → not found at a threshold requiring the A.
    // All-pyrimidine intron: every candidate lacks the conserved A (position 0),
    // so the best score is (0.790+0.746+0+0.751)/3.210 = 0.712461, below 0.8.
    [Test]
    public void FindAcceptorBranchPoint_NoBranchAdenosine_NotFoundAtAdenosineThreshold()
    {
        string sequence = "UU" + "CUUUC" + new string('U', 22) + "AG" + "GGG";

        var bp = FindAcceptorBranchPoint(sequence, acceptorAgPosition: 30, minScore: 0.8);

        Assert.That(bp.Found, Is.False,
            "Without the conserved branch adenosine, the best score (0.712) is below 0.8 — Gao et al. (2008)");
    }

    // BP4 — Degenerate branch point with a purine at position -3 scores below a
    // perfect yUnAy but is still found. Motif AUUAC: A@-3 (purine, loses the 0.790
    // pyrimidine term), U@-2, A@0, C@+1 → (0+0.746+0.923+0.751)/3.210 = 0.753894.
    // Purine (G) tract isolates the single branch A so no other candidate competes.
    [Test]
    public void FindAcceptorBranchPoint_PurineAtMinus3_ScoresBelowPerfectButFound()
    {
        string sequence = "GG" + "AUUAC" + new string('G', 22) + "AG" + "GGG";

        var bp = FindAcceptorBranchPoint(sequence, acceptorAgPosition: 30);

        Assert.Multiple(() =>
        {
            Assert.That(bp.Found, Is.True,
                "Branch adenosine present → candidate qualifies above default 0.5 threshold");
            Assert.That(bp.BranchPointPosition, Is.EqualTo(5),
                "Branch-point adenosine at index 5");
            Assert.That(bp.Motif, Is.EqualTo("AUUAC"),
                "Reported motif spans positions -3..+1");
            Assert.That(bp.Score, Is.EqualTo(0.753894).Within(1e-6),
                "Purine at -3 forfeits the 0.790 pyrimidine term → score = 0.753894 — Gao et al. (2008) frequencies");
            Assert.That(bp.PolypyrimidineTractFraction, Is.EqualTo(0.0).Within(1e-10),
                "Purine (G) tract between branch point and AG → pyrimidine fraction = 0.0");
        });
    }

    // BP5 — Branch point exactly at the near window edge (distance 18) is found;
    // one nt closer (distance 17) is outside the 18–40 nt window → not found.
    [Test]
    public void FindAcceptorBranchPoint_NearWindowEdge_18Found_17NotFound()
    {
        // distance = agEnd − branchA. branchA = 5 in both; vary the G fill so the AG lands
        // at distance 18 vs 17. Purine (G) fill so only the real branch A scores.
        // dist 18: GG + CUUAC + G×15 + AG + GGG  → G of AG index 23, distance 23−5 = 18.
        string seq18 = "GG" + "CUUAC" + new string('G', 15) + "AG" + "GGG";
        // dist 17: GG + CUUAC + G×14 + AG + GGG  → G of AG index 22, distance 22−5 = 17.
        string seq17 = "GG" + "CUUAC" + new string('G', 14) + "AG" + "GGG";

        var bp18 = FindAcceptorBranchPoint(seq18, acceptorAgPosition: 23);
        var bp17 = FindAcceptorBranchPoint(seq17, acceptorAgPosition: 22);

        Assert.Multiple(() =>
        {
            Assert.That(bp18.Found, Is.True,
                "Branch point 18 nt upstream is at the inclusive near edge of the window → found");
            Assert.That(bp18.DistanceFromAg, Is.EqualTo(18),
                "Distance at the near edge is 18 nt");
            Assert.That(bp17.Found, Is.False,
                "Branch point 17 nt upstream is closer than the 18 nt minimum → not found");
        });
    }

    // BP6 — Branch point at the far window edge (distance 40) is found; one nt
    // farther (distance 41) is outside the window → not found.
    [Test]
    public void FindAcceptorBranchPoint_FarWindowEdge_40Found_41NotFound()
    {
        // dist 40: GG + CUUAC + G×37 + AG + GGG → G of AG index 45, distance 45−5 = 40.
        string seq40 = "GG" + "CUUAC" + new string('G', 37) + "AG" + "GGG";
        // dist 41: GG + CUUAC + G×38 + AG + GGG → G of AG index 46, distance 46−5 = 41.
        string seq41 = "GG" + "CUUAC" + new string('G', 38) + "AG" + "GGG";

        var bp40 = FindAcceptorBranchPoint(seq40, acceptorAgPosition: 45);
        var bp41 = FindAcceptorBranchPoint(seq41, acceptorAgPosition: 46);

        Assert.Multiple(() =>
        {
            Assert.That(bp40.Found, Is.True,
                "Branch point 40 nt upstream is at the inclusive far edge of the window → found");
            Assert.That(bp40.DistanceFromAg, Is.EqualTo(40),
                "Distance at the far edge is 40 nt");
            Assert.That(bp41.Found, Is.False,
                "Branch point 41 nt upstream is farther than the 40 nt maximum → not found");
        });
    }

    // BP7 — Guards: empty/null sequence and out-of-range AG position → not found.
    [Test]
    public void FindAcceptorBranchPoint_InvalidInput_ReturnsNotFound()
    {
        var nullResult = FindAcceptorBranchPoint(null!, acceptorAgPosition: 30);
        var emptyResult = FindAcceptorBranchPoint("", acceptorAgPosition: 30);
        var oob = FindAcceptorBranchPoint("CUUAC" + new string('U', 30), acceptorAgPosition: 999);
        var zeroPos = FindAcceptorBranchPoint("CUUAC" + new string('U', 30), acceptorAgPosition: 0);

        Assert.Multiple(() =>
        {
            Assert.That(nullResult.Found, Is.False, "Null sequence → not found");
            Assert.That(nullResult.BranchPointPosition, Is.EqualTo(-1), "Not-found position sentinel is -1");
            Assert.That(emptyResult.Found, Is.False, "Empty sequence → not found");
            Assert.That(oob.Found, Is.False, "AG position past sequence end → not found");
            Assert.That(zeroPos.Found, Is.False, "AG position 0 leaves no upstream window → not found");
        });
    }

    // BP9 — Multiple qualifying candidates in the window: the HIGHEST-scoring one wins,
    // not the nearest. A weaker degenerate branch point (AUUAC, A@-3 → 0.753894) sits
    // 20 nt upstream (scanned first, nearest); a perfect yUnAy (CUUAC → 1.0) sits 30 nt
    // upstream (scanned later, farther). The detector must report the perfect far one.
    //   Layout: G-pad, CUUAC at index 7..11 (branchA 10, dist 30),
    //           AUUAC at index 17..21 (branchA 20, dist 20), G-fill, AG at 39..40.
    [Test]
    public void FindAcceptorBranchPoint_MultipleCandidates_SelectsHighestScoringNotNearest()
    {
        string sequence =
            "GGGGGGG" + "CUUAC" + "GGGGG" + "AUUAC" + new string('G', 17) + "AG" + "GGG";

        var bp = FindAcceptorBranchPoint(sequence, acceptorAgPosition: 40);

        Assert.Multiple(() =>
        {
            Assert.That(bp.Found, Is.True, "Both candidates qualify above the default 0.5 threshold");
            Assert.That(bp.BranchPointPosition, Is.EqualTo(10),
                "The perfect yUnAy at index 10 (30 nt upstream) wins over the weaker one at index 20 (20 nt)");
            Assert.That(bp.DistanceFromAg, Is.EqualTo(30),
                "Reported distance is that of the highest-scoring candidate (30 nt), not the nearest");
            Assert.That(bp.Motif, Is.EqualTo("CUUAC"), "Best candidate motif is the perfect yUnAy");
            Assert.That(bp.Score, Is.EqualTo(1.0).Within(1e-10),
                "Best candidate is the perfect-consensus branch point → score 1.0");
        });
    }

    // BP8 — DNA (T) input is treated identically to RNA (U) input.
    [Test]
    public void FindAcceptorBranchPoint_DnaInput_MatchesRnaResult()
    {
        string rna = "UU" + "CUUAC" + new string('U', 22) + "AG" + "GGG";
        string dna = "TT" + "CTTAC" + new string('T', 22) + "AG" + "GGG";

        var bpRna = FindAcceptorBranchPoint(rna, acceptorAgPosition: 30);
        var bpDna = FindAcceptorBranchPoint(dna, acceptorAgPosition: 30);

        Assert.Multiple(() =>
        {
            Assert.That(bpDna.Found, Is.EqualTo(bpRna.Found), "DNA and RNA inputs both detect the branch point");
            Assert.That(bpDna.BranchPointPosition, Is.EqualTo(bpRna.BranchPointPosition), "Same branch-point position");
            Assert.That(bpDna.Score, Is.EqualTo(bpRna.Score).Within(1e-10), "Same score (T↔U normalization)");
            Assert.That(bpDna.Motif, Is.EqualTo(bpRna.Motif), "Motif reported in RNA alphabet for both");
        });
    }

    #endregion

    #region FBP: FindBranchPoints — human yUnAy (Gao 2008) / yeast UACUAAC

    [Test]
    public void FindBranchPoints_HumanYUnAy_DefaultModel_PositionIsBranchA()
    {
        // Gao et al. (2008) yUnAy: CUAAC at 7..11 (branch A index 10) = perfect → 1.0; UUCAG
        // (A at 23): y,U,n,A + G@+1 → (0.790+0.746+0.923)/3.210.
        const string seq = "GGGGGUACUAACGGGGGUUUUUCAGG";
        var sites = FindBranchPoints(seq, minScore: 0.75).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(sites.Select(s => (s.Position, s.Motif)), Is.EqualTo(new[] { (10, "CUAAC"), (23, "UUCAG") }));
            Assert.That(sites[0].Score, Is.EqualTo(1.0));
            Assert.That(sites[1].Score, Is.EqualTo((0.790 + 0.746 + 0.923) / 3.210).Within(1e-15));
            Assert.That(sites.All(s => s.Type == SpliceSiteType.Branch), Is.True);
        });
    }

    [Test]
    public void FindBranchPoints_YeastModel_UACUAAC_MatchFraction()
    {
        // S. cerevisiae UACUAAC (branch A = 6th nt): exact match → 7/7 at A index 10; UACUGAC → 6/7.
        var perfect = FindBranchPoints("GGGGGUACUAACGGGG", 0, -1, 0.9, BranchPointModel.YeastUACUAAC).Single();
        var oneOff = FindBranchPoints("GGGGGUACUGACGGGG", 0, -1, 0.8, BranchPointModel.YeastUACUAAC).Single();

        Assert.Multiple(() =>
        {
            Assert.That((perfect.Position, perfect.Motif, perfect.Score), Is.EqualTo((10, "UACUAAC", 1.0)));
            Assert.That(oneOff.Position, Is.EqualTo(10));
            Assert.That(oneOff.Score, Is.EqualTo(6.0 / 7).Within(1e-15));
        });
    }

    [Test]
    public void FindBranchPoints_NonAcguInMotif_Skipped_UnknownModelThrows()
    {
        Assert.Multiple(() =>
        {
            Assert.That(FindBranchPoints("CUNAC", minScore: 0.0).ToList(), Is.Empty,
                "a placement with N is not scored");
            Assert.That(FindBranchPoints("CUUAC", minScore: 0.0).Single().Score, Is.EqualTo(1.0));
            Assert.Throws<ArgumentOutOfRangeException>(() =>
                FindBranchPoints("CUUAC", 0, -1, 0.5, (BranchPointModel)9));
        });
    }

    [Test]
    public void FindAcceptorBranchPoint_NonAcguInMotif_PlacementSkipped()
    {
        // Same as BP1 but the yUnAy at distance 25 contains an N → that placement is not scored;
        // the best remaining placement in 18..40 is reported instead (or none at 0.8).
        string sequence = "UU" + "CNUAC" + new string('U', 22) + "AG" + "GGG";

        var bp = FindAcceptorBranchPoint(sequence, acceptorAgPosition: 30, minScore: 0.8);

        Assert.That(bp.Found, Is.False);
    }

    #endregion

    #region ME1–ME9: MaxEntScan score3ss — Yeo & Burge (2004)

    // Expected values are the documented MaxEntScan score3 worked examples, taken from the
    // maxentpy score3 docstring (kepbod/maxentpy, MIT-licensed port retrieved 2026-06-24):
    //   score3('ttccaaacgaacttttgtAGgga') -> 2.89
    //   score3('tgtctttttctgtgtggcAGtgg') -> 8.19
    //   score3('ttctctcttcagacttatAGcaa') -> -0.08
    // These round-to-2dp targets were reproduced exactly this session (2.886773 / 8.190965 /
    // -0.080278) with the embedded tables + the published factorisation. Asserted to 2 dp
    // (the documented precision) plus the full-precision value, so a wrong table or
    // factorisation fails the 2.89 cross-check rather than silently passing.

    // ME1 — canonical documented reference: score3('...AGgga') == 2.89 (bits).
    [Test]
    public void ScoreAcceptorMaxEnt_CanonicalReferenceWindow_Returns2Point89()
    {
        // 23 nt: 20 intron + 3 exon; AG at 0-based positions 18-19.
        string window = "ttccaaacgaacttttgtAGgga";

        double score = ScoreAcceptorMaxEnt(window);

        Assert.That(System.Math.Round(score, 2), Is.EqualTo(2.89),
            "MaxEntScan score3 of the canonical documented example must be 2.89 bits "
            + "(Yeo & Burge 2004 / maxentpy docstring)");
    }

    // ME2 — full-precision value behind the 2.89 reference (guards the factorisation/tables).
    [Test]
    public void ScoreAcceptorMaxEnt_CanonicalReferenceWindow_MatchesFullPrecision()
    {
        string window = "ttccaaacgaacttttgtAGgga";

        double score = ScoreAcceptorMaxEnt(window);

        Assert.That(score, Is.EqualTo(2.886773).Within(1e-6),
            "Full-precision MaxEntScan score3 for the canonical window, reproduced from the "
            + "embedded tables + published factorisation");
    }

    // ME3 — strong 3' site second documented reference: 8.19 bits.
    [Test]
    public void ScoreAcceptorMaxEnt_StrongSiteReferenceWindow_Returns8Point19()
    {
        string window = "tgtctttttctgtgtggcAGtgg";

        double score = ScoreAcceptorMaxEnt(window);

        Assert.Multiple(() =>
        {
            Assert.That(System.Math.Round(score, 2), Is.EqualTo(8.19),
                "Second documented MaxEntScan score3 example must be 8.19 bits");
            Assert.That(score, Is.EqualTo(8.190965).Within(1e-6),
                "Full-precision value behind the 8.19 reference");
        });
    }

    // ME4 — weak 3' site third documented reference: -0.08 bits (negative score).
    [Test]
    public void ScoreAcceptorMaxEnt_WeakSiteReferenceWindow_ReturnsMinus0Point08()
    {
        string window = "ttctctcttcagacttatAGcaa";

        double score = ScoreAcceptorMaxEnt(window);

        Assert.Multiple(() =>
        {
            Assert.That(System.Math.Round(score, 2), Is.EqualTo(-0.08),
                "Third documented MaxEntScan score3 example must be -0.08 bits");
            Assert.That(score, Is.EqualTo(-0.080278).Within(1e-6),
                "Full-precision value behind the -0.08 reference");
        });
    }

    // ME5 — a strong site ranks above a weak site (ordering of the documented examples).
    [Test]
    public void ScoreAcceptorMaxEnt_StrongSite_RanksAboveWeakSite()
    {
        double strong = ScoreAcceptorMaxEnt("tgtctttttctgtgtggcAGtgg"); // 8.19
        double weak = ScoreAcceptorMaxEnt("ttctctcttcagacttatAGcaa");   // -0.08

        Assert.That(strong, Is.GreaterThan(weak),
            "A strong 3' splice site must score higher than a weak one (8.19 > -0.08)");
    }

    // ME6 — DNA (T) and RNA (U) windows are equivalent (T==U in the hash and AG model).
    [Test]
    public void ScoreAcceptorMaxEnt_DnaAndRnaWindows_ProduceIdenticalScores()
    {
        string dna = "ttccaaacgaacttttgtAGgga";
        string rna = dna.Replace('t', 'u').Replace('T', 'U');

        double dnaScore = ScoreAcceptorMaxEnt(dna);
        double rnaScore = ScoreAcceptorMaxEnt(rna);

        Assert.That(rnaScore, Is.EqualTo(dnaScore).Within(1e-12),
            "U-form (RNA) window must score identically to the T-form (DNA) window");
    }

    // ME7 — case-insensitive: upper-case window scores identically.
    [Test]
    public void ScoreAcceptorMaxEnt_UpperCaseWindow_ScoresIdenticallyToLowerCase()
    {
        string lower = "ttccaaacgaacttttgtAGgga";
        string upper = lower.ToUpperInvariant();

        Assert.That(ScoreAcceptorMaxEnt(upper), Is.EqualTo(ScoreAcceptorMaxEnt(lower)).Within(1e-12),
            "Scoring must be case-insensitive");
    }

    // ME8 — null window throws ArgumentNullException.
    [Test]
    public void ScoreAcceptorMaxEnt_NullWindow_Throws()
    {
        Assert.Throws<System.ArgumentNullException>(() => ScoreAcceptorMaxEnt(null!),
            "A null window is invalid input");
    }

    // ME9 — wrong length and invalid alphabet throw ArgumentException.
    [Test]
    public void ScoreAcceptorMaxEnt_InvalidWindow_Throws()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<System.ArgumentException>(() => ScoreAcceptorMaxEnt("ACGT"),
                "A window not exactly 23 nt is invalid");
            Assert.Throws<System.ArgumentException>(() => ScoreAcceptorMaxEnt("ttccaaacgaacttttgtAGggaa"),
                "A 24-nt window is invalid");
            Assert.Throws<System.ArgumentException>(() => ScoreAcceptorMaxEnt("ttccaaacgaacttttgtAGggN"),
                "A non-A/C/G/T(/U) character (N) is invalid");
        });
    }

    // ME10 — the model SCORES (does not reject) a window whose 18-19 dinucleotide is not AG;
    // the consensus AG term heavily penalises it, yielding a strongly negative score. This
    // matches maxentpy score3, which performs no AG validation. Expected values reproduced
    // from the independent oracle (maxentpy score3 over the embedded tables, 2026-06-25):
    //   score3('ttccaaacgaacttttgtCCgga') -> -13.220039
    //   score3('ttccaaacgaacttttgtTTgga') -> -14.078362
    // The canonical AG window (same flanks) scores +2.886773, so AG vastly outscores non-AG.
    [Test]
    public void ScoreAcceptorMaxEnt_NonAgDinucleotide_ScoredNotRejected_StronglyNegative()
    {
        const double canonicalAg = 2.886773; // ME2 reference (…AG…)

        double cc = ScoreAcceptorMaxEnt("ttccaaacgaacttttgtCCgga");
        double tt = ScoreAcceptorMaxEnt("ttccaaacgaacttttgtTTgga");

        Assert.Multiple(() =>
        {
            Assert.That(cc, Is.EqualTo(-13.220039).Within(1e-6),
                "Non-AG (CC) window is scored, not rejected — maxentpy score3 has no AG check");
            Assert.That(tt, Is.EqualTo(-14.078362).Within(1e-6),
                "Non-AG (TT) window is scored, not rejected — maxentpy score3 has no AG check");
            Assert.That(cc, Is.LessThan(canonicalAg),
                "A non-AG acceptor dinucleotide scores far below the canonical AG site");
            Assert.That(tt, Is.LessThan(canonicalAg),
                "A non-AG acceptor dinucleotide scores far below the canonical AG site");
        });
    }

    // ME10 — reference sweep (SPLICE-MAXENT3-001, B11 F21/F22). Expected values are maxentpy score3
    // (kepbod/maxentpy maxent.py + score3_matrix.txt) full-precision outputs (Python repr, so the C#
    // literal is the same double); the CPAN port of score3.pl (Bio-RNA-SpliceSites-Scoring-MaxEntScan
    // 0.05, original me2x3acc1-9 tables) prints the 2-decimal values. The C# score equals maxentpy
    // bit-for-bit on 200,025 windows (100,000 random + 100,000 AG 23-mers + 16 dinucleotides + examples);
    // with Math.Log2 instead of ln/ln2, 62,349 of them (incl. the 8.19 docstring example) were 1 ulp off.
    // Includes non-AG dinucleotides (AC, GG, TA, ...), which MaxEntScan scores rather than rejects.
    [TestCase("ttccaaacgaacttttgtAGgga", 2.8867730651152104, 2.89)]
    [TestCase("tgtctttttctgtgtggcAGtgg", 8.190965487162021, 8.19)]
    [TestCase("ttctctcttcagacttatAGcaa", -0.08027774948780869, -0.08)]
    [TestCase("CTTTCCCGGTCTCCCCGCAGGGC", 13.974478220613745, 13.97)]
    [TestCase("CCCCCCCCCCCCCCCCCCAGCCC", 10.350554250918666, 10.35)]
    [TestCase("AGGCACGGTCCGTCCGGTAGCGC", 2.3223895715101035, 2.32)]
    [TestCase("TCAATGATCTCCGGGTTAAGTTA", 2.2213509402965332, 2.22)]
    [TestCase("TTATCCAACTTGCAAGCAAGGGC", -3.3326429559313904, -3.33)]
    [TestCase("CGTATAGCCCTCGTGATGAGGAT", -6.014271156784488, -6.01)]
    [TestCase("ttccaaacgaacttttgtACGGA", -5.177714829032689, -5.18)]
    [TestCase("ttccaaacgaacttttgtGGGGA", -5.068088471471686, -5.07)]
    [TestCase("GGGGAGTTTATCCTTCTACGTAA", -9.482356641325355, -9.48)]
    [TestCase("TGCTTAAGGAGCTCACCTAAACC", -15.042687466794712, -15.04)]
    [TestCase("TTTTTTTTTTTTTTTTTTTTTTT", -5.737465487038207, -5.74)]
    [TestCase("AAAAAAAAAAAAAAAAAAAAAAA", -16.26511154573316, -16.27)]
    [TestCase("TCAAGCAGAGGAGAAATGCTTGA", -61.106105250879395, -61.11)]
    public void ScoreAcceptorMaxEnt_ReferenceSweep_MatchesMaxEntScan(string window, double maxentpy, double score3Pl)
    {
        double score = ScoreAcceptorMaxEnt(window);
        Assert.Multiple(() =>
        {
            Assert.That(score, Is.EqualTo(maxentpy), "maxentpy score3 (bit-identical double)");
            Assert.That(System.Math.Round(score, 2, System.MidpointRounding.AwayFromZero), Is.EqualTo(score3Pl),
                "score3.pl prints %.2f");
        });
    }

    #endregion
}
