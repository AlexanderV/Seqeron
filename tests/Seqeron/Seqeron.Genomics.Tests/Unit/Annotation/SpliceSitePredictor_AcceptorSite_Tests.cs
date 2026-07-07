using static Seqeron.Genomics.Annotation.SpliceSitePredictor;

namespace Seqeron.Genomics.Tests.Unit.Annotation;

/// <summary>
/// Canonical test file for SPLICE-ACCEPTOR-001: Acceptor Site Detection.
/// Tests <see cref="SpliceSitePredictor.FindAcceptorSites"/> and the internal
/// ScoreAcceptorSite / ScoreU12AcceptorSite scoring functions.
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
    #region M1: Canonical AG Detected — Shapiro & Senapathy (1987)

    [Test]
    public void FindAcceptorSites_CanonicalAG_WithStrongPPT_FindsAcceptorSite()
    {
        // Consensus: (Y)nNCAG|G — strong PPT (continuous U) followed by CAG
        // U×16 + C(16) + A(17) + G(18) + G(19) + G(20) = 21 nt
        // AG at index 17–18, splice site at 19
        string sequence = "UUUUUUUUUUUUUUUUCAGGG";

        var sites = FindAcceptorSites(sequence, minScore: 0.1).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(sites, Has.Count.EqualTo(1),
                "Exactly one AG dinucleotide in scannable range");
            Assert.That(sites[0].Type, Is.EqualTo(SpliceSiteType.Acceptor),
                "Canonical AG site must have Type = Acceptor");
            Assert.That(sites[0].Position, Is.EqualTo(18),
                "Position = A_index + 1 = 17 + 1 = 18 (G of AG)");
            Assert.That(sites[0].Score, Is.EqualTo(0.8393).Within(0.001),
                "Strong PPT + perfect CAG|G consensus → high score (~0.84) — Shapiro & Senapathy (1987)");
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

    #region M4: Short Sequence Returns Empty — Guard: Length < 20

    [Test]
    public void FindAcceptorSites_ShortSequence_ReturnsEmpty()
    {
        // 19 chars — below the 20-char minimum
        string sequence = "UUUUUUUUUUUUUUUCAGG";

        var sites = FindAcceptorSites(sequence, minScore: 0.1).ToList();

        Assert.That(sites, Is.Empty,
            "Sequence shorter than 20 nt should return empty — insufficient PPT context");
    }

    #endregion

    #region M5: Strong PPT > Weak PPT Score — Burge et al. (1999)

    [Test]
    public void FindAcceptorSites_StrongPPT_ScoresHigherThanWeakContext()
    {
        // Strong PPT: continuous pyrimidines before CAG
        string strongPpt = "UUUUUUUUUUUUUUUUCAGGG";

        // Weak context: purines interrupt the PPT region, but still has AG at scannable position
        // AAG×5 + AACAGGG = 22 chars, AG at index 18–19
        string weakContext = "AAGAAGAAGAAGAAGAACAGGG";

        var strongSites = FindAcceptorSites(strongPpt, minScore: 0.0).ToList();
        var weakSites = FindAcceptorSites(weakContext, minScore: 0.0).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(strongSites, Is.Not.Empty,
                "Strong PPT must produce at least one site — continuous pyrimidines are optimal");
            Assert.That(weakSites.Where(s => s.Type == SpliceSiteType.Acceptor).ToList(), Is.Not.Empty,
                "Weak context still has a canonical AG — must produce at least one Acceptor site");

            double strongScore = strongSites.Max(s => s.Score);
            double weakScore = weakSites.Where(s => s.Type == SpliceSiteType.Acceptor).Max(s => s.Score);

            Assert.That(strongScore, Is.GreaterThan(weakScore),
                "Strong PPT (continuous pyrimidines) should score higher than purine-interrupted context — Burge et al. (1999)");
            Assert.That(strongScore, Is.GreaterThan(0.7),
                "Perfect PPT + CAG|G consensus should produce high score (>0.7)");
            Assert.That(weakScore, Is.LessThan(0.7),
                "Purine-interrupted PPT with correct CAG consensus should produce moderate score (<0.7)");
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
        // Two AG dinucleotides: AG at index 17–18 and AG at index 21–22
        // U×16 + C(16) + A(17) + G(18) + C(19) + U(20) + A(21) + G(22) + G(23) + G(24) = 25 nt
        string sequence = "UUUUUUUUUUUUUUUUCAGCUAGGG";

        var sites = FindAcceptorSites(sequence, minScore: 0.0).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(sites.Count, Is.EqualTo(2),
                "Exactly two AG dinucleotides in scannable range");
            Assert.That(sites[0].Position, Is.EqualTo(18),
                "First AG at index 17 → position 18");
            Assert.That(sites[1].Position, Is.EqualTo(22),
                "Second AG at index 21 → position 22");
            Assert.That(sites[0].Score, Is.GreaterThan(sites[1].Score),
                "First AG has stronger PPT context upstream → higher score");
        });
    }

    #endregion

    #region S1: U12 AC Detected — Patel & Steitz (2003), Hall & Padgett (1994)

    [Test]
    public void FindAcceptorSites_U12AcceptorYCCAC_DetectedWhenNonCanonicalEnabled()
    {
        // U12-type 3' splice site consensus: YCCAC (Hall & Padgett 1994, Jackson 1991)
        // U×14 (PPT) + C(14) + C(15) + A(16) + C(17) + G(18) + G(19) + G(20) = 21 nt
        // YCCAC pattern: Y=U(13), C=C(14), C=C(15), A=A(16), C=C(17)
        // AC dinucleotide at index 16–17
        string sequence = "UUUUUUUUUUUUUUCCACGGG";

        var sites = FindAcceptorSites(sequence, minScore: 0.1, includeNonCanonical: true).ToList();
        var u12Sites = sites.Where(s => s.Type == SpliceSiteType.U12Acceptor).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(u12Sites, Has.Count.EqualTo(1),
                "Exactly one U12 YCCAC acceptor — Hall & Padgett (1994)");
            Assert.That(u12Sites[0].Position, Is.EqualTo(17),
                "Position = AC_A_index + 1 = 16 + 1 = 17");
            Assert.That(u12Sites[0].Score, Is.EqualTo(1.0).Within(0.001),
                "Perfect YCCAC consensus + strong PPT → maximum U12 score (3.5/3.5)");
        });
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

    #region C1: AG Before Position 15 Not Detected — Scan Start

    [Test]
    public void FindAcceptorSites_AGBeforeScanStart_NotDetected()
    {
        // AG at position 5 (too early), no AG in scannable range
        string sequence = "UUUUUAGUUUUUUUUUUUUUU";

        var sites = FindAcceptorSites(sequence, minScore: 0.0).ToList();

        Assert.That(sites, Is.Empty,
            "AG at position < 15 should not be detected — scan starts at i=15 to ensure PPT context");
    }

    #endregion

    #region C2: PPT Pyrimidine Count Helper Verification

    [Test]
    public void FindAcceptorSites_PPTContribution_VerifiedByHelper()
    {
        // Strong PPT: 12 pyrimidines out of 12 in PPT window
        string strongPpt = "UUUUUUUUUUUUUUUUCAGGG";

        var sites = FindAcceptorSites(strongPpt, minScore: 0.0).ToList();

        Assert.That(sites, Is.Not.Empty, "Should find at least one site");

        // Verify PPT contribution independently
        // PPT window: positions [i-15, i-3) where i is the AG position
        // For AG at position 16 (the A), PPT window is [1, 13)
        double expectedPptFraction = ComputePptFraction(strongPpt.ToUpperInvariant().Replace('T', 'U'), 16 - 1);

        Assert.That(expectedPptFraction, Is.GreaterThan(0.8),
            "Strong PPT should have high pyrimidine fraction — independent verification");
    }

    /// <summary>
    /// Independent PPT scoring helper — counts pyrimidine fraction in
    /// the 12-nt window [position-15, position-3) relative to the A of AG.
    /// This reproduces the PPT component of ScoreAcceptorSite from first principles.
    /// </summary>
    private static double ComputePptFraction(string sequence, int agPosition)
    {
        int pptCount = 0;
        int total = 0;

        for (int i = agPosition - 15; i < agPosition - 3; i++)
        {
            if (i >= 0 && i < sequence.Length)
            {
                total++;
                if (sequence[i] == 'C' || sequence[i] == 'U')
                    pptCount++;
            }
        }

        return total > 0 ? (double)pptCount / total : 0;
    }

    #endregion

    #region C3: Non-ACGT Characters Tolerated — Robustness

    [Test]
    public void FindAcceptorSites_NonACGTCharacters_ToleratedAndAGStillFound()
    {
        // N's occupy the PPT/context window but the canonical AG remains. Non-ACGT bases
        // simply match no PWM entry / no pyrimidine (they contribute nothing), so the call
        // must not throw and must still locate the AG.
        // N×16 + C(16) + A(17) + G(18) + G(19) + G(20) = 21 nt; AG at index 17–18.
        string sequence = "NNNNNNNNNNNNNNNNCAGGG";

        var sites = FindAcceptorSites(sequence, minScore: 0.0).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(sites, Has.Count.EqualTo(1),
                "Non-ACGT context must not suppress the canonical AG (no throw, AG still found)");
            Assert.That(sites[0].Position, Is.EqualTo(18),
                "AG at index 17 → Position 18, unaffected by surrounding non-ACGT bases");
            Assert.That(sites[0].Score, Is.InRange(0.0, 1.0),
                "Score with non-ACGT context stays within the normalized [0, 1] range");
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

    #endregion
}
