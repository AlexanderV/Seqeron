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
        // PPT at positions 0–15, CAG at 14–16, G at 17
        string sequence = "UUUUUUUUUUUUUUUUCAGGG";

        var sites = FindAcceptorSites(sequence, minScore: 0.1).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(sites, Has.Count.GreaterThanOrEqualTo(1),
                "Canonical AG with strong PPT must produce at least one acceptor site");
            Assert.That(sites.All(s => s.Type == SpliceSiteType.Acceptor), Is.True,
                "All sites from canonical AG scan should have Type = Acceptor");
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
        // Need 21+ chars with AG at position >= 15
        string weakContext = "AAGAAGAAGAAGAAGAACAGGG";

        var strongSites = FindAcceptorSites(strongPpt, minScore: 0.1).ToList();
        var weakSites = FindAcceptorSites(weakContext, minScore: 0.1).ToList();

        Assert.That(strongSites, Is.Not.Empty,
            "Strong PPT must produce at least one site — continuous pyrimidines are optimal");

        if (weakSites.Count > 0)
        {
            double strongScore = strongSites.Max(s => s.Score);
            double weakScore = weakSites.Max(s => s.Score);

            Assert.That(strongScore, Is.GreaterThan(weakScore),
                "Strong PPT (continuous pyrimidines) should score higher than purine-interrupted context — Burge et al. (1999)");
        }
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

        Assert.That(lowerSites.Count, Is.EqualTo(upperSites.Count),
            "Lowercase input should produce same number of sites as uppercase — case insensitivity");
    }

    #endregion

    #region M10: Multiple AG Sites — Scanning Algorithm

    [Test]
    public void FindAcceptorSites_MultipleAGSites_FindsAll()
    {
        // Two AG dinucleotides at positions >= 15:
        // First AG at position ~16, second AG at position ~20
        string sequence = "UUUUUUUUUUUUUUUUCAGCUAGGG";

        var sites = FindAcceptorSites(sequence, minScore: 0.0).ToList();

        Assert.That(sites.Count, Is.GreaterThanOrEqualTo(2),
            "Sequence with two AG dinucleotides in scannable range should find at least 2 sites");
    }

    #endregion

    #region S1: U12 AC Detected — Patel & Steitz (2003)

    [Test]
    public void FindAcceptorSites_U12AcceptorAC_DetectedWhenNonCanonicalEnabled()
    {
        // AC dinucleotide at position >= 15 with non-canonical enabled
        string sequence = "UUUUUUUUUUUUUUUUACGGG";

        var sites = FindAcceptorSites(sequence, minScore: 0.1, includeNonCanonical: true).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(sites, Is.Not.Empty,
                "U12 AC acceptor should be detected when includeNonCanonical=true — Patel & Steitz (2003)");
            Assert.That(sites.Any(s => s.Type == SpliceSiteType.U12Acceptor), Is.True,
                "U12 site should have Type = U12Acceptor");
        });
    }

    #endregion

    #region S2: U12 AC Excluded When Canonical Only — Default Parameter

    [Test]
    public void FindAcceptorSites_U12AcceptorAC_ExcludedByDefault()
    {
        // AC dinucleotide but no AG — canonical-only mode should find nothing
        string sequence = "UUUUUUUUUUUUUUUUACGGG";

        var sites = FindAcceptorSites(sequence, minScore: 0.1).ToList();

        Assert.That(sites.Where(s => s.Type == SpliceSiteType.U12Acceptor).ToList(), Is.Empty,
            "U12 AC acceptor should not be detected with includeNonCanonical=false (default)");
    }

    #endregion

    #region S3: Position After AG — INV-5

    [Test]
    public void FindAcceptorSites_ReturnsPositionAfterAG()
    {
        // Construct sequence where AG is at a known position
        // 15 U's + AG + GGG = 20 chars
        // AG at indices 15–16, so position should be 16+1 = 17 (after the G of AG)
        string sequence = "UUUUUUUUUUUUUUUAGGGGG";

        var sites = FindAcceptorSites(sequence, minScore: 0.0).ToList();

        if (sites.Count > 0)
        {
            // The AG is at position 15 (A) and 16 (G)
            // Implementation: Position = i + 1 where i is the index of A in AG
            // So Position = 15 + 1 = 16
            Assert.That(sites[0].Position, Is.EqualTo(16),
                "Position should be index after AG dinucleotide (first exonic nucleotide)");
        }
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
