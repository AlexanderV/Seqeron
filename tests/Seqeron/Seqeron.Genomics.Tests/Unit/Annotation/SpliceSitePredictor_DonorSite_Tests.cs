using static Seqeron.Genomics.Annotation.SpliceSitePredictor;

namespace Seqeron.Genomics.Tests.Unit.Annotation;

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
        // Theory: all 9 positions match → score = 9/9 = 1.0
        // Confidence = (1.0 - 0.5) / (1.0 - 0.5) = 1.0
        const string sequence = "CAGGUAAGU";
        var sites = FindDonorSites(sequence, minScore: 0.3).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(sites, Has.Count.EqualTo(1),
                "Perfect GU consensus CAG|GUAAGU must produce exactly one donor site");
            Assert.That(sites[0].Type, Is.EqualTo(SpliceSiteType.Donor),
                "Site type must be Donor for canonical GT/GU dinucleotide");
            Assert.That(sites[0].Score, Is.EqualTo(1.0).Within(1e-10),
                "All 9 IUPAC positions match MAG|GURAGU: C∈M, A=A, G=G, G=G, U=U, A∈R, A=A, G=G, U=U → 9/9");
            Assert.That(sites[0].Confidence, Is.EqualTo(1.0).Within(1e-10),
                "Confidence = (1.0 − 0.5) / (1.0 − 0.5) = 1.0");
            Assert.That(sites[0].Position, Is.EqualTo(3),
                "GU dinucleotide starts at index 3 in CAGGUAAGU");
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
        // Strong: CAG|GUAAGU → all 9 positions match MAG|GURAGU → 9/9 = 1.0
        // Weak: UUU|GUAAUU → U∉M, U≠A, U≠G, G=G, U=U, A∈R, A=A, U≠G, U=U → 5/9
        const string strongSequence = "CAGGUAAGU";
        const string weakSequence = "UUUGUAAUU";

        var strongSites = FindDonorSites(strongSequence, minScore: 0.0).ToList();
        var weakSites = FindDonorSites(weakSequence, minScore: 0.0).ToList();

        Assert.That(strongSites, Has.Count.EqualTo(1));
        Assert.That(weakSites, Has.Count.EqualTo(1));

        Assert.Multiple(() =>
        {
            Assert.That(strongSites[0].Score, Is.EqualTo(1.0).Within(1e-10),
                "Perfect consensus: 9/9 = 1.0");
            Assert.That(weakSites[0].Score, Is.EqualTo(5.0 / 9).Within(1e-10),
                "Weak context: 5 of 9 match (G,U at +0,+1; A∈R at +2; A at +3; U at +5) → 5/9");
            Assert.That(strongSites[0].Score, Is.GreaterThan(weakSites[0].Score),
                "Perfect consensus must score higher than weak context per PWM model");
        });
    }

    #endregion

    #region M6: GC Non-Canonical Detected When Enabled — Burge et al. (1999)

    [Test]
    public void FindDonorSites_GC_NonCanonical_DetectedWhenEnabled()
    {
        // GC-AG introns are valid U2-type (~0.5-1% of introns) — Burge et al. (1999)
        // GC at position +1 mismatches invariant U consensus → 8/9 positions match
        // Theory: C∈M, A=A, G=G, G=G, C≠U, A∈R, A=A, G=G, U=U → 8/9
        const string sequence = "CAGGCAAGU";
        var sites = FindDonorSites(sequence, minScore: 0.0, includeNonCanonical: true).ToList();

        // GC at position 3 is the only donor-like dinucleotide in scan window
        Assert.Multiple(() =>
        {
            Assert.That(sites, Has.Count.EqualTo(1),
                "Exactly one GC donor in scan window");
            Assert.That(sites[0].Type, Is.EqualTo(SpliceSiteType.Donor),
                "GC donors are classified as Donor type (U2 spliceosome)");
            Assert.That(sites[0].Score, Is.EqualTo(8.0 / 9).Within(1e-10),
                "GC donor: position +1 (C) mismatches invariant U → 8/9");
            Assert.That(sites[0].Position, Is.EqualTo(3),
                "GC dinucleotide is at index 3");
        });
    }

    #endregion

    #region M7: GC Not Detected With Canonical Only — Burge et al. (1999)

    [Test]
    public void FindDonorSites_GC_NotDetected_WhenCanonicalOnly()
    {
        // Sequence with GC donor context but no GT/GU dinucleotide anywhere
        const string gcOnlySequence = "CAGGCAACC";
        var sites = FindDonorSites(gcOnlySequence, minScore: 0.0, includeNonCanonical: false).ToList();

        // No GT/GU present anywhere in sequence
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

        Assert.Multiple(() =>
        {
            Assert.That(dnaSites.Count, Is.EqualTo(rnaSites.Count),
                "DNA (T) and RNA (U) must produce the same number of donor sites");
            Assert.That(dnaSites[0].Score, Is.EqualTo(rnaSites[0].Score).Within(1e-10),
                "Scores must be identical after T→U normalization");
            Assert.That(dnaSites[0].Score, Is.EqualTo(1.0).Within(1e-10),
                "Both DNA and RNA perfect consensus score 9/9 = 1.0");
            Assert.That(dnaSites[0].Position, Is.EqualTo(rnaSites[0].Position),
                "Site positions must be identical");
        });
    }

    #endregion

    #region M9: Lowercase Handling — Implementation ToUpperInvariant

    [Test]
    public void FindDonorSites_LowercaseInput_FindsSite()
    {
        const string lowercase = "cagguaagu";
        var sites = FindDonorSites(lowercase, minScore: 0.3).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(sites, Has.Count.EqualTo(1),
                "Lowercase input must be handled via case-insensitive normalization");
            Assert.That(sites[0].Type, Is.EqualTo(SpliceSiteType.Donor),
                "Lowercase input must produce correct Donor type classification");
            Assert.That(sites[0].Score, Is.EqualTo(1.0).Within(1e-10),
                "Lowercase perfect consensus must score identically to uppercase: 9/9 = 1.0");
        });
    }

    #endregion

    #region M10: Multiple GU Sites All Detected — Full Sequence Scan

    [Test]
    public void FindDonorSites_MultipleSites_AllDetected()
    {
        // Two consensus motifs (pos 3, 26) plus one cryptic GU (pos 7) in scan window
        // Position 3: CAGGUAAGU → 9/9 = 1.0
        // Position 7: context UAAGUUUUU → 4/9 (G,U at +0,+1; A at -2; U at +5)
        // Position 26: CAGGUAAGU → 9/9 = 1.0
        const string sequence = "CAGGUAAGUUUUUUUUUUUUUUUCAGGUAAGU";
        var sites = FindDonorSites(sequence, minScore: 0.3).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(sites.Count, Is.EqualTo(3),
                "Three GU dinucleotides in scan window: positions 3, 7, 26");
            Assert.That(sites, Has.All.Property(nameof(SpliceSite.Type)).EqualTo(SpliceSiteType.Donor),
                "All detected sites must be classified as Donor");
            Assert.That(sites.Select(s => s.Position).ToArray(),
                Is.EqualTo(new[] { 3, 7, 26 }),
                "Site positions must match all GU locations within scan range");
            Assert.That(sites[0].Score, Is.EqualTo(1.0).Within(1e-10),
                "First consensus site (pos 3): 9/9");
            Assert.That(sites[1].Score, Is.EqualTo(4.0 / 9).Within(1e-10),
                "Cryptic site (pos 7): 4/9");
            Assert.That(sites[2].Score, Is.EqualTo(1.0).Within(1e-10),
                "Second consensus site (pos 26): 9/9");
        });
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
        // GT: CAGGUAAGU → 9/9 = 1.0 (all match)
        // GC: CAGGCAAGU → 8/9 (position +1: C mismatches invariant U)
        const string gtSequence = "CAGGUAAGU";
        const string gcSequence = "CAGGCAAGU";

        var gtSites = FindDonorSites(gtSequence, minScore: 0.0).ToList();
        var gcSites = FindDonorSites(gcSequence, minScore: 0.0, includeNonCanonical: true).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(gtSites, Has.Count.EqualTo(1), "GT sequence must yield one site");
            Assert.That(gcSites, Has.Count.EqualTo(1), "GC sequence must yield one site");
            Assert.That(gtSites[0].Score, Is.EqualTo(1.0).Within(1e-10),
                "GT donor: 9/9 = 1.0");
            Assert.That(gcSites[0].Score, Is.EqualTo(8.0 / 9).Within(1e-10),
                "GC donor: 8/9 (C≠U at position +1)");
            Assert.That(gcSites[0].Score, Is.LessThan(gtSites[0].Score),
                "GC donor must score lower than GT donor — Burge et al. (1999)");
        });
    }

    #endregion

    #region C3: Donor At Sequence Start — Partial-Window Boundary

    [Test]
    public void FindDonorSites_DonorAtSequenceStart_ScoresOverInBoundsPositionsOnly()
    {
        // GU at index 0: upstream offsets -3, -2, -1 fall before the sequence and are
        // out of bounds, so ScoreDonorSite normalizes over only the 6 in-bounds positions
        // (offsets 0..+5). GUAAGU... has G,U,A,A,G,U at offsets 0..+5 →
        // 0:G✓ +1:U✓ +2:A∈R✓ +3:A✓ +4:G✓ +5:U✓ → 6/6 = 1.0.
        // The score is a fraction over POSITIONS SCORED, so a partial window of perfect
        // matches still reaches 1.0 (hand-computed, not a code echo).
        const string sequence = "GUAAGUAAA";
        var sites = FindDonorSites(sequence, minScore: 0.0).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(sites, Has.Count.EqualTo(1),
                "Single GU at the sequence start must be detected");
            Assert.That(sites[0].Position, Is.EqualTo(0),
                "Donor at the very start reports Position 0 (G of GU)");
            Assert.That(sites[0].Score, Is.EqualTo(1.0).Within(1e-10),
                "Partial window (6 in-bounds positions, offsets 0..+5) all match → 6/6 = 1.0");
        });
    }

    #endregion

    #region C4: Non-ACGT Characters — Skipped, No Throw

    [Test]
    public void FindDonorSites_NonAcgtCharacterAtScoredPosition_IsSkippedNotCounted()
    {
        // A non-A/C/G/U base (N) at a scored offset is not a PWM key, so ScoreDonorSite
        // skips it WITHOUT incrementing the divisor (graceful degradation, no throw).
        // NAGGUAAGU: GU at index 3; offset -3 → index 0 = 'N' is skipped, leaving 8
        // scored positions, all matching → 8/8 = 1.0. Hand-computed from the consensus.
        const string sequence = "NAGGUAAGU";
        var sites = FindDonorSites(sequence, minScore: 0.0).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(sites, Has.Count.EqualTo(1),
                "GU is still located despite an upstream non-ACGT base");
            Assert.That(sites[0].Position, Is.EqualTo(3),
                "GU dinucleotide is at index 3");
            Assert.That(sites[0].Score, Is.EqualTo(1.0).Within(1e-10),
                "N at offset -3 is skipped (not counted); remaining 8 positions match → 8/8 = 1.0");
        });
    }

    [Test]
    public void FindDonorSites_AllNonAcgt_ReturnsEmpty()
    {
        // A sequence with no G/U at all (and non-ACGT noise) must yield no donor sites
        // and must not throw.
        const string sequence = "NNNNNNNNN";
        Assert.DoesNotThrow(() => FindDonorSites(sequence, minScore: 0.0).ToList());
        Assert.That(FindDonorSites(sequence, minScore: 0.0).ToList(), Is.Empty,
            "No GU dinucleotide present → no donor sites");
    }

    #endregion

    #region Helper: Independent PWM Score Verification

    /// <summary>
    /// Independently computes a simplified donor site score to verify the PWM model direction.
    /// Uses the IUPAC consensus MAG|GURAGU from Shapiro &amp; Senapathy (1987),
    /// Mount (1982), Burge et al. (1999).
    /// This helper is NOT copied from the implementation; it encodes the expected
    /// conservation pattern from the published consensus.
    /// </summary>
    private static double ComputeExpectedConsensusStrength(string ninemerRna)
    {
        // IUPAC consensus: MAG|GURAGU
        // M = A/C, R = A/G. Match = 1, no match = 0.
        char[][] allowed =
        {
            new[] { 'A', 'C' }, // -3: M
            new[] { 'A' },       // -2: A
            new[] { 'G' },       // -1: G
            new[] { 'G' },       //  0: G (invariant)
            new[] { 'U' },       //  1: U (invariant)
            new[] { 'A', 'G' }, //  2: R
            new[] { 'A' },       //  3: A
            new[] { 'G' },       //  4: G
            new[] { 'U' },       //  5: U
        };

        double matches = 0;
        string upper = ninemerRna.ToUpperInvariant().Replace('T', 'U');
        for (int i = 0; i < Math.Min(upper.Length, 9); i++)
        {
            if (allowed[i].Contains(upper[i]))
                matches++;
        }
        return matches;
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

    #endregion
}
