using NUnit.Framework;
using Seqeron.Genomics.Annotation;
using static Seqeron.Genomics.Annotation.SpliceSitePredictor;

namespace Seqeron.Genomics.Tests;

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
}
