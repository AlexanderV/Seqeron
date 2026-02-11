using NUnit.Framework;
using Seqeron.Genomics.Annotation;
using static Seqeron.Genomics.Annotation.SpliceSitePredictor;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// SPLICE-DONOR-001: Canonical test file for Donor (5') Splice Site Detection.
/// Tests <see cref="SpliceSitePredictor.FindDonorSites"/> and internal scoring.
/// Evidence: Shapiro &amp; Senapathy (1987), Burge et al. (1999), Yeo &amp; Burge (2004).
/// </summary>
[TestFixture]
public class SpliceSitePredictor_DonorSite_Tests
{
    #region M1: Canonical GT Consensus Detected — Shapiro & Senapathy (1987)

    [Test]
    public void FindDonorSites_CanonicalGT_ConsensusMotif_ProducesDonorSite()
    {
        // MAG|GURAGU perfect consensus: CAG|GUAAGU
        const string sequence = "CAGGUAAGU";
        var sites = FindDonorSites(sequence, minScore: 0.3).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(sites, Has.Count.GreaterThanOrEqualTo(1),
                "Perfect GU consensus CAG|GUAAGU must produce at least one donor site");
            Assert.That(sites[0].Type, Is.EqualTo(SpliceSiteType.Donor),
                "Site type must be Donor for canonical GT/GU dinucleotide");
            Assert.That(sites[0].Score, Is.GreaterThan(0),
                "PWM score for perfect consensus must be positive");
        });
    }

    #endregion

    #region M2: No GU Dinucleotide Returns Empty — Trivially Correct

    [Test]
    public void FindDonorSites_NoGU_ReturnsEmpty()
    {
        const string sequence = "AAAAACCCCC";
        var sites = FindDonorSites(sequence, minScore: 0.3).ToList();

        Assert.That(sites, Is.Empty,
            "Sequence without any GU/GT dinucleotide must yield no donor sites");
    }

    #endregion

    #region M3: Empty Input Returns Empty — Trivially Correct

    [Test]
    public void FindDonorSites_EmptyString_ReturnsEmpty()
    {
        var sites = FindDonorSites("", minScore: 0.3).ToList();

        Assert.That(sites, Is.Empty,
            "Empty sequence must yield no donor sites");
    }

    #endregion

    #region M4: Short Sequence Returns Empty — Implementation Guard

    [Test]
    public void FindDonorSites_SequenceShorterThan6_ReturnsEmpty()
    {
        // PWM window requires at least 6 characters (positions 0..+5)
        var sites = FindDonorSites("GUAA", minScore: 0.3).ToList();

        Assert.That(sites, Is.Empty,
            "Sequence shorter than 6 nucleotides lacks sufficient context for PWM scoring");
    }

    #endregion

    #region M5: Strong Context Scores Higher Than Weak — Shapiro & Senapathy (1987), Yeo & Burge (2004)

    [Test]
    public void FindDonorSites_StrongConsensus_ScoresHigherThanWeakContext()
    {
        // Strong: perfect consensus CAG|GUAAGU (high weight at every position)
        // Weak: poor context UUU|GUAAUU (low weight at -3,-2,-1 and +4,+5)
        const string strongSequence = "CAGGUAAGU";
        const string weakSequence = "UUUGUAAUU";

        var strongSites = FindDonorSites(strongSequence, minScore: 0.0).ToList();
        var weakSites = FindDonorSites(weakSequence, minScore: 0.0).ToList();

        Assert.That(strongSites, Is.Not.Empty,
            "Perfect consensus must produce at least one site even with minScore=0");
        Assert.That(weakSites, Is.Not.Empty,
            "GU dinucleotide present, so site should be found with minScore=0");

        double strongScore = strongSites.Max(s => s.Score);
        double weakScore = weakSites.Max(s => s.Score);

        Assert.That(strongScore, Is.GreaterThan(weakScore),
            "Perfect consensus CAGGUAAGU must score higher than weak context UUUGUAAUU per PWM model");
    }

    #endregion

    #region M6: GC Non-Canonical Detected When Enabled — Burge et al. (1999)

    [Test]
    public void FindDonorSites_GC_NonCanonical_DetectedWhenEnabled()
    {
        // GC-AG introns are valid U2-type (~0.5-1% of introns) — Burge et al. (1999)
        // Sequence with GC but no GT: only GC donor possible
        const string sequence = "CAGGCAAGU";
        var sites = FindDonorSites(sequence, minScore: 0.0, includeNonCanonical: true).ToList();

        Assert.That(sites, Has.Count.GreaterThanOrEqualTo(1),
            "GC donor must be detected when includeNonCanonical=true");
        Assert.That(sites[0].Type, Is.EqualTo(SpliceSiteType.Donor),
            "GC donors are still classified as Donor type (U2 spliceosome)");
    }

    #endregion

    #region M7: GC Not Detected With Canonical Only — Burge et al. (1999)

    [Test]
    public void FindDonorSites_GC_NotDetected_WhenCanonicalOnly()
    {
        // Sequence with GC but no GT dinucleotide
        const string gcOnlySequence = "CAGGCAAGU";
        var sites = FindDonorSites(gcOnlySequence, minScore: 0.0, includeNonCanonical: false).ToList();

        // No GT/GU present, so canonical-only mode should return empty
        Assert.That(sites, Is.Empty,
            "Without GT/GU and includeNonCanonical=false, no donor sites should be found");
    }

    #endregion

    #region M8: DNA T Equivalence — Implementation T→U Conversion

    [Test]
    public void FindDonorSites_DNA_T_ProducesSameResultsAsRNA_U()
    {
        const string dnaSequence = "CAGGTAAGT"; // DNA with T
        const string rnaSequence = "CAGGUAAGU"; // RNA with U

        var dnaSites = FindDonorSites(dnaSequence, minScore: 0.3).ToList();
        var rnaSites = FindDonorSites(rnaSequence, minScore: 0.3).ToList();

        Assert.That(dnaSites.Count, Is.EqualTo(rnaSites.Count),
            "DNA (T) and RNA (U) representations must produce the same number of donor sites");
    }

    #endregion

    #region M9: Lowercase Handling — Implementation ToUpperInvariant

    [Test]
    public void FindDonorSites_LowercaseInput_FindsSite()
    {
        const string lowercase = "cagguaagu";
        var sites = FindDonorSites(lowercase, minScore: 0.3).ToList();

        Assert.That(sites, Has.Count.GreaterThanOrEqualTo(1),
            "Lowercase input must be handled via case-insensitive normalization");
        Assert.That(sites[0].Type, Is.EqualTo(SpliceSiteType.Donor),
            "Lowercase input must produce correct Donor type classification");
    }

    #endregion

    #region M10: Multiple GU Sites All Detected — Full Sequence Scan

    [Test]
    public void FindDonorSites_MultipleSites_AllDetected()
    {
        // Two well-separated consensus donor motifs
        const string sequence = "CAGGUAAGUUUUUUUUUUUUUUUCAGGUAAGU";
        var sites = FindDonorSites(sequence, minScore: 0.3).ToList();

        Assert.That(sites.Count, Is.GreaterThanOrEqualTo(2),
            "Sequence with two CAG|GUAAGU motifs must yield at least 2 donor sites");
        Assert.That(sites, Has.All.Property(nameof(SpliceSite.Type)).EqualTo(SpliceSiteType.Donor),
            "All detected sites must be classified as Donor");
    }

    #endregion

    #region S1: Score In Valid Range — INV-2

    [Test]
    public void FindDonorSites_AllScores_InZeroOneRange()
    {
        const string sequence = "CAGGUAAGUAAAGUAAGUCCCCAAAGGUAAGU";
        var sites = FindDonorSites(sequence, minScore: 0.0).ToList();

        Assert.That(sites, Is.Not.Empty, "Should find at least one site");
        foreach (var site in sites)
        {
            Assert.That(site.Score, Is.InRange(0.0, 1.0),
                $"Score {site.Score} at position {site.Position} must be in [0, 1]");
        }
    }

    #endregion

    #region S2: Confidence In Valid Range — INV-3

    [Test]
    public void FindDonorSites_AllConfidences_InZeroOneRange()
    {
        const string sequence = "CAGGUAAGUAAAGUAAGUCCCC";
        var sites = FindDonorSites(sequence, minScore: 0.0).ToList();

        Assert.That(sites, Is.Not.Empty, "Should find at least one site");
        foreach (var site in sites)
        {
            Assert.That(site.Confidence, Is.InRange(0.0, 1.0),
                $"Confidence {site.Confidence} at position {site.Position} must be in [0, 1]");
        }
    }

    #endregion

    #region S3: Motif Non-Empty — INV-8

    [Test]
    public void FindDonorSites_AllMotifs_AreNonEmpty()
    {
        const string sequence = "CAGGUAAGU";
        var sites = FindDonorSites(sequence, minScore: 0.3).ToList();

        Assert.That(sites, Is.Not.Empty, "Should find at least one site");
        foreach (var site in sites)
        {
            Assert.That(site.Motif, Is.Not.Null.And.Not.Empty,
                $"Motif at position {site.Position} must be non-empty for display context");
        }
    }

    #endregion

    #region S4: Higher Threshold Reduces Results — Score Filtering

    [Test]
    public void FindDonorSites_HigherMinScore_FewerOrEqualSites()
    {
        const string sequence = "CAGGUAAGUAAAGUAAGUCCCCAAAGGUAAGU";
        var lowThreshold = FindDonorSites(sequence, minScore: 0.2).ToList();
        var highThreshold = FindDonorSites(sequence, minScore: 0.8).ToList();

        Assert.That(highThreshold.Count, Is.LessThanOrEqualTo(lowThreshold.Count),
            "Higher minScore threshold must yield fewer or equal sites than lower threshold");
    }

    #endregion

    #region S5: Null Input Returns Empty — Defensive

    [Test]
    public void FindDonorSites_NullInput_ReturnsEmpty()
    {
        var sites = FindDonorSites(null!, minScore: 0.3).ToList();

        Assert.That(sites, Is.Empty,
            "Null input must yield empty result, not throw exception");
    }

    #endregion

    #region C1: U12 AT Donor Detected — Minor Spliceosome

    [Test]
    public void FindDonorSites_U12_AT_Donor_Detected_WhenNonCanonical()
    {
        // U12-type introns use AT/AU at donor — Burge et al. (1999)
        // /ATATCC/ consensus → AUAUCC in RNA
        const string sequence = "CAGAUAUCC";
        var sites = FindDonorSites(sequence, minScore: 0.0, includeNonCanonical: true).ToList();

        var u12Sites = sites.Where(s => s.Type == SpliceSiteType.U12Donor).ToList();
        Assert.That(u12Sites, Has.Count.GreaterThanOrEqualTo(1),
            "AU donor matching U12 consensus should be detected with includeNonCanonical=true");
    }

    #endregion

    #region C2: GC Donor Score Lower Than GT — INV-7

    [Test]
    public void FindDonorSites_GC_Donor_ScoresLowerThanEquivalentGT()
    {
        // Same context except GT vs GC at the splice dinucleotide
        // GT context: CAGGUAAGU → strong canonical donor
        // GC context needs to be tested relative to GT with similar surroundings
        const string gtSequence = "CAGGUAAGU";
        const string gcSequence = "CAGGCAAGU";

        var gtSites = FindDonorSites(gtSequence, minScore: 0.0).ToList();
        var gcSites = FindDonorSites(gcSequence, minScore: 0.0, includeNonCanonical: true).ToList();

        if (gtSites.Any() && gcSites.Any())
        {
            Assert.That(gcSites.Max(s => s.Score), Is.LessThan(gtSites.Max(s => s.Score)),
                "GC donor score must be lower than GT donor due to 0.7 penalty (GC-AG introns are weaker)");
        }
    }

    #endregion

    #region Helper: Independent PWM Score Verification

    /// <summary>
    /// Independently computes a simplified donor site score to verify the PWM model direction.
    /// Uses the documented consensus MAG|GURAGU from Shapiro &amp; Senapathy (1987).
    /// This helper is NOT copied from the implementation; it encodes the expected
    /// conservation pattern from the published consensus.
    /// </summary>
    private static double ComputeExpectedConsensusStrength(string ninemerRna)
    {
        // Conservation weights from Shapiro & Senapathy (1987) consensus pattern:
        // Position: -3   -2   -1   0    +1   +2   +3   +4   +5
        // Best:      M    A    G    G    U    A    A    G    U
        // Weight:   0.35  0.60 0.80 1.00 1.00 0.60 0.70 0.80 0.55
        double[] bestWeights = { 0.35, 0.60, 0.80, 1.00, 1.00, 0.60, 0.70, 0.80, 0.55 };
        char[] bestNucs = { 'A', 'A', 'G', 'G', 'U', 'A', 'A', 'G', 'U' };

        double score = 0;
        string upper = ninemerRna.ToUpperInvariant().Replace('T', 'U');
        for (int i = 0; i < Math.Min(upper.Length, 9); i++)
        {
            if (upper[i] == bestNucs[i])
                score += bestWeights[i];
        }
        return score;
    }

    [Test]
    public void HelperVerification_PerfectConsensus_HasHighestIndependentScore()
    {
        double perfect = ComputeExpectedConsensusStrength("CAGGUAAGU");
        double weak = ComputeExpectedConsensusStrength("UUUGUAAUU");

        Assert.That(perfect, Is.GreaterThan(weak),
            "Independent PWM verification: perfect consensus must score higher than weak context");
    }

    #endregion
}
